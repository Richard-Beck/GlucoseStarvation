#!/usr/bin/env python3
"""
Build segmentation artifacts and object-level embedding feature table from an annotation stack.

Outputs:
1) Object label stack TIFF (0 background, 1..N object ids per frame)
2) Object-level feature table with:
   - metadata (frame/object/area/centroid/border)
   - label: alive vs not_alive (object contains >=1 annotation point)
   - frame ResNet18 embedding
   - masked-object ResNet18 embedding
   - log(object_area)

Notes:
- Designed for GPU execution (Cellpose batched inference + batched embedding forward passes).
- Annotation points can come from any channel; channel is ignored for labeling.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import tifffile
import torch
import torch.nn as nn
from cellpose.models import CellposeModel
from scipy import ndimage
from torchvision import models
from torchvision.models import ResNet18_Weights


EPS = 1e-8


@dataclass
class ObjMeta:
    frame: int  # 1-based
    object_id: int
    area_px: int
    centroid_x: float
    centroid_y: float
    bbox_x0: int
    bbox_y0: int
    bbox_x1: int
    bbox_y1: int
    touches_border: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Segment annotation stack, derive alive/not_alive object labels, and compute embeddings."
    )
    parser.add_argument("--input_stack_tif", required=True, help="Input 90-image annotation stack (expected TCYX).")
    parser.add_argument(
        "--annotation_results_csv",
        required=True,
        help="ImageJ Results CSV with columns including X,Y,Frame (channel ignored).",
    )
    parser.add_argument(
        "--manifest_csv",
        default=None,
        help="Optional stack manifest CSV to append condition metadata by frame index.",
    )
    parser.add_argument(
        "--output_labels_tif",
        required=True,
        help="Output object labels stack TIFF (shape T,Y,X; 0 background, 1..N objects/frame).",
    )
    parser.add_argument(
        "--output_feature_table",
        required=True,
        help="Output feature table path (.parquet preferred; .csv also supported).",
    )
    parser.add_argument(
        "--output_meta_json",
        required=True,
        help="Output run metadata JSON.",
    )
    parser.add_argument(
        "--cellpose_batch_size",
        type=int,
        default=32,
        help="Batch size for Cellpose eval over crops.",
    )
    parser.add_argument(
        "--embed_batch_size",
        type=int,
        default=512,
        help="Batch size for ResNet18 embedding forwards.",
    )
    parser.add_argument(
        "--point_snap_radius",
        type=int,
        default=3,
        help="If point is on background, search this radius for nearest labeled pixel.",
    )
    parser.add_argument(
        "--exclude_border_objects",
        action="store_true",
        default=True,
        help="Exclude objects touching crop boundary from final training table.",
    )
    parser.add_argument(
        "--include_border_objects",
        action="store_true",
        help="Override and include border objects in final table.",
    )
    parser.add_argument(
        "--resnet_weights",
        default="imagenet",
        help="ResNet18 weights: 'imagenet', 'none', or path to .pth state_dict.",
    )
    return parser.parse_args()


def robust01(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=np.float32)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.zeros_like(arr, dtype=np.float32)
    vals = arr[finite]
    lo = float(np.percentile(vals, 1))
    hi = float(np.percentile(vals, 99))
    if hi <= lo:
        lo = float(np.min(vals))
        hi = float(np.max(vals))
    if hi <= lo:
        return np.zeros_like(arr, dtype=np.float32)
    out = (arr - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0)


def load_stack_tcxy(path: str) -> np.ndarray:
    arr = tifffile.imread(path)
    if arr.ndim != 4:
        raise ValueError(f"Expected 4D stack (TCYX), got shape {arr.shape}")
    # Heuristic: we expect channel dimension to be 3.
    if arr.shape[1] == 3:
        return arr
    if arr.shape[0] == 3:
        # likely CTYX -> transpose to TCYX
        return np.transpose(arr, (1, 0, 2, 3))
    raise ValueError(f"Could not interpret stack channel axis for shape {arr.shape}")


def make_cellpose_inputs(stack_tcxy: np.ndarray) -> List[np.ndarray]:
    # Input channel order expected from prior scripts: [phase, dead, alive]
    # Build Cellpose input with 2 channels [phase, fluor_sum] and channels=[2,1]
    frames = []
    t, c, h, w = stack_tcxy.shape
    if c != 3:
        raise ValueError(f"Expected 3 channels (phase/dead/alive), got {c}")
    for i in range(t):
        phase = stack_tcxy[i, 0].astype(np.float32)
        dead = stack_tcxy[i, 1].astype(np.float32)
        alive = stack_tcxy[i, 2].astype(np.float32)

        phase01 = robust01(phase)
        fluor01 = robust01(dead + alive)
        inp = np.stack([phase01, fluor01], axis=0).astype(np.float32, copy=False)
        frames.append(inp)
    return frames


def run_cellpose_batched(inputs: List[np.ndarray], batch_size: int, use_gpu: bool) -> List[np.ndarray]:
    model = CellposeModel(gpu=use_gpu)
    masks, _, _ = model.eval(
        inputs,
        diameter=None,
        channels=[2, 1],
        batch_size=batch_size,
    )
    # Cellpose may return list; normalize to list[np.ndarray]
    out = []
    for m in masks:
        out.append(np.asarray(m))
    return out


def save_labels_stack(labels: List[np.ndarray], out_path: str) -> str:
    max_id = max(int(m.max()) for m in labels) if labels else 0
    dtype = np.uint16 if max_id <= np.iinfo(np.uint16).max else np.uint32
    stack = np.stack([m.astype(dtype, copy=False) for m in labels], axis=0)  # TYX
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    tifffile.imwrite(out_path, stack, imagej=False)
    return str(dtype)


def extract_object_metadata_for_frame(labels: np.ndarray, frame_1based: int) -> List[ObjMeta]:
    max_id = int(labels.max())
    if max_id == 0:
        return []

    y_idx, x_idx = np.indices(labels.shape)
    flat = labels.ravel()
    areas = np.bincount(flat, minlength=max_id + 1)
    sum_x = np.bincount(flat, weights=x_idx.ravel(), minlength=max_id + 1)
    sum_y = np.bincount(flat, weights=y_idx.ravel(), minlength=max_id + 1)

    # bbox via scipy slices (index corresponds to object id starting at 1)
    slices = ndimage.find_objects(labels)

    h, w = labels.shape
    border_ids = set(np.unique(np.concatenate([
        labels[0, :],
        labels[h - 1, :],
        labels[:, 0],
        labels[:, w - 1],
    ])))
    border_ids.discard(0)

    out: List[ObjMeta] = []
    for oid in range(1, max_id + 1):
        area = int(areas[oid])
        if area <= 0:
            continue
        sl = slices[oid - 1]
        if sl is None:
            continue
        ys, xs = sl
        x0, x1 = int(xs.start), int(xs.stop)
        y0, y1 = int(ys.start), int(ys.stop)
        cx = float(sum_x[oid] / max(area, 1))
        cy = float(sum_y[oid] / max(area, 1))
        out.append(
            ObjMeta(
                frame=frame_1based,
                object_id=oid,
                area_px=area,
                centroid_x=cx,
                centroid_y=cy,
                bbox_x0=x0,
                bbox_y0=y0,
                bbox_x1=x1,
                bbox_y1=y1,
                touches_border=(oid in border_ids),
            )
        )
    return out


def load_points_by_frame(results_csv: str) -> Dict[int, List[Tuple[float, float]]]:
    points: Dict[int, List[Tuple[float, float]]] = {}
    with open(results_csv, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        needed = {"X", "Y", "Frame"}
        missing = needed - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Annotation CSV missing required columns: {sorted(missing)}")
        for row in reader:
            if row.get("Frame", "").strip() == "":
                continue
            fr = int(float(row["Frame"]))
            x = float(row["X"])
            y = float(row["Y"])
            points.setdefault(fr, []).append((x, y))
    return points


def point_to_object_id(
    labels: np.ndarray,
    x: float,
    y: float,
    snap_radius: int,
) -> int:
    h, w = labels.shape
    xi = int(np.floor(x))
    yi = int(np.floor(y))
    xi = min(max(xi, 0), w - 1)
    yi = min(max(yi, 0), h - 1)

    oid = int(labels[yi, xi])
    if oid != 0 or snap_radius <= 0:
        return oid

    x0 = max(0, xi - snap_radius)
    x1 = min(w, xi + snap_radius + 1)
    y0 = max(0, yi - snap_radius)
    y1 = min(h, yi + snap_radius + 1)
    patch = labels[y0:y1, x0:x1]
    ys, xs = np.where(patch > 0)
    if ys.size == 0:
        return 0
    # nearest labeled pixel
    dx = (xs + x0 - x)
    dy = (ys + y0 - y)
    dist2 = dx * dx + dy * dy
    j = int(np.argmin(dist2))
    return int(patch[ys[j], xs[j]])


def make_resnet18_embedder(weights_spec: str, device: torch.device) -> nn.Module:
    if weights_spec == "imagenet":
        try:
            weights = ResNet18_Weights.DEFAULT
            net = models.resnet18(weights=weights)
        except Exception as exc:
            raise RuntimeError(
                "Failed to load ImageNet weights. Use --resnet_weights none or provide local .pth"
            ) from exc
    elif weights_spec == "none":
        net = models.resnet18(weights=None)
    else:
        net = models.resnet18(weights=None)
        sd = torch.load(weights_spec, map_location="cpu")
        net.load_state_dict(sd)

    # penultimate 512-d embedding
    embedder = nn.Sequential(*list(net.children())[:-1])
    embedder.eval().to(device)
    return embedder


def imagenet_normalize(batch: torch.Tensor) -> torch.Tensor:
    # batch in [0,1], shape N,3,H,W
    mean = torch.tensor([0.485, 0.456, 0.406], dtype=batch.dtype, device=batch.device).view(1, 3, 1, 1)
    std = torch.tensor([0.229, 0.224, 0.225], dtype=batch.dtype, device=batch.device).view(1, 3, 1, 1)
    return (batch - mean) / std


def resize_batch(batch: torch.Tensor, out_hw: int = 224) -> torch.Tensor:
    return torch.nn.functional.interpolate(batch, size=(out_hw, out_hw), mode="bilinear", align_corners=False)


def run_embed_batches(embedder: nn.Module, imgs: np.ndarray, batch_size: int, device: torch.device) -> np.ndarray:
    # imgs shape N,H,W,3 float32 in [0,1]
    n = imgs.shape[0]
    out = np.zeros((n, 512), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, n, batch_size):
            x = torch.from_numpy(imgs[i : i + batch_size]).permute(0, 3, 1, 2).to(device=device, dtype=torch.float32)
            x = resize_batch(x, 224)
            x = imagenet_normalize(x)
            z = embedder(x).flatten(1)  # N,512
            out[i : i + x.shape[0], :] = z.detach().cpu().numpy().astype(np.float32, copy=False)
    return out


def run_embed_patch_batches(
    embedder: nn.Module,
    imgs: List[np.ndarray],
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    # imgs is a list of H,W,3 float32 arrays in [0,1] with potentially varying H,W.
    n = len(imgs)
    out = np.zeros((n, 512), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, n, batch_size):
            batch_tensors = []
            for img in imgs[i : i + batch_size]:
                x = torch.from_numpy(img).permute(2, 0, 1).unsqueeze(0).to(device=device, dtype=torch.float32)
                x = resize_batch(x, 224)
                batch_tensors.append(x)
            x = torch.cat(batch_tensors, dim=0)
            x = imagenet_normalize(x)
            z = embedder(x).flatten(1)  # N,512
            out[i : i + x.shape[0], :] = z.detach().cpu().numpy().astype(np.float32, copy=False)
    return out


def frame_rgb_for_embedding(stack_tcxy: np.ndarray) -> np.ndarray:
    # Build pseudo-RGB as [dead, alive, phase] then robust normalize per channel per frame.
    t, c, h, w = stack_tcxy.shape
    assert c == 3
    out = np.zeros((t, h, w, 3), dtype=np.float32)
    for i in range(t):
        phase = robust01(stack_tcxy[i, 0].astype(np.float32))
        dead = robust01(stack_tcxy[i, 1].astype(np.float32))
        alive = robust01(stack_tcxy[i, 2].astype(np.float32))
        out[i, :, :, 0] = dead
        out[i, :, :, 1] = alive
        out[i, :, :, 2] = phase
    return out


def object_rgb_patch_masked(stack_tcxy: np.ndarray, labels_by_frame: List[np.ndarray], obj: ObjMeta) -> np.ndarray:
    f = obj.frame - 1
    labels = labels_by_frame[f]
    y0, y1 = obj.bbox_y0, obj.bbox_y1
    x0, x1 = obj.bbox_x0, obj.bbox_x1

    phase = stack_tcxy[f, 0, y0:y1, x0:x1].astype(np.float32)
    dead = stack_tcxy[f, 1, y0:y1, x0:x1].astype(np.float32)
    alive = stack_tcxy[f, 2, y0:y1, x0:x1].astype(np.float32)

    mask = (labels[y0:y1, x0:x1] == obj.object_id)

    # normalize channels inside bbox, then mask background to 0
    phase_n = robust01(phase)
    dead_n = robust01(dead)
    alive_n = robust01(alive)

    rgb = np.zeros((phase.shape[0], phase.shape[1], 3), dtype=np.float32)
    rgb[:, :, 0] = dead_n
    rgb[:, :, 1] = alive_n
    rgb[:, :, 2] = phase_n
    rgb[~mask] = 0.0
    return rgb


def maybe_merge_manifest(df: pd.DataFrame, manifest_csv: Optional[str]) -> pd.DataFrame:
    if manifest_csv is None:
        return df
    m = pd.read_csv(manifest_csv)
    if "stack_index_1based" not in m.columns:
        return df
    m2 = m.rename(columns={"stack_index_1based": "frame"})
    keep = [
        c
        for c in [
            "frame",
            "image_key",
            "cellLine",
            "ploidy",
            "G0",
            "glucose_bin",
            "time_bin",
            "hours",
            "base_key",
            "target_x0",
            "target_y0",
            "target_x1",
            "target_y1",
        ]
        if c in m2.columns
    ]
    return df.merge(m2[keep], on="frame", how="left")


def main() -> None:
    args = parse_args()
    exclude_border = args.exclude_border_objects and (not args.include_border_objects)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    stack = load_stack_tcxy(args.input_stack_tif)
    t, c, h, w = stack.shape
    if c != 3:
        raise ValueError(f"Expected 3-channel stack, got shape {stack.shape}")

    cpose_inputs = make_cellpose_inputs(stack)
    labels_by_frame = run_cellpose_batched(cpose_inputs, args.cellpose_batch_size, use_gpu=(device.type == "cuda"))

    # Save segmentation artifact first
    labels_dtype = save_labels_stack(labels_by_frame, args.output_labels_tif)

    # Object metadata
    objs: List[ObjMeta] = []
    for fi, lbl in enumerate(labels_by_frame, start=1):
        objs.extend(extract_object_metadata_for_frame(lbl, frame_1based=fi))

    # Map points -> alive object IDs
    pts_by_frame = load_points_by_frame(args.annotation_results_csv)
    alive_ids_by_frame: Dict[int, set] = {}
    unmatched_points = 0
    for fr, pts in pts_by_frame.items():
        if fr < 1 or fr > len(labels_by_frame):
            continue
        lbl = labels_by_frame[fr - 1]
        aset = alive_ids_by_frame.setdefault(fr, set())
        for x, y in pts:
            oid = point_to_object_id(lbl, x=x, y=y, snap_radius=args.point_snap_radius)
            if oid > 0:
                aset.add(oid)
            else:
                unmatched_points += 1

    # Filter objects and assign labels
    rows = []
    kept_objs: List[ObjMeta] = []
    for o in objs:
        if exclude_border and o.touches_border:
            continue
        alive = (o.object_id in alive_ids_by_frame.get(o.frame, set()))
        rows.append(
            {
                "frame": o.frame,
                "object_id": o.object_id,
                "area_px": o.area_px,
                "log_area_px": float(math.log(max(o.area_px, 1))),
                "centroid_x": o.centroid_x,
                "centroid_y": o.centroid_y,
                "bbox_x0": o.bbox_x0,
                "bbox_y0": o.bbox_y0,
                "bbox_x1": o.bbox_x1,
                "bbox_y1": o.bbox_y1,
                "touches_border": bool(o.touches_border),
                "label": "alive" if alive else "not_alive",
                "label_int": 1 if alive else 0,
            }
        )
        kept_objs.append(o)

    df = pd.DataFrame(rows)

    # Frame embeddings (once/frame)
    embedder = make_resnet18_embedder(args.resnet_weights, device)
    frame_rgb = frame_rgb_for_embedding(stack)
    frame_emb = run_embed_batches(embedder, frame_rgb, batch_size=args.embed_batch_size, device=device)

    frame_emb_cols = [f"frame_emb_{i:03d}" for i in range(frame_emb.shape[1])]
    frame_emb_df = pd.DataFrame(frame_emb, columns=frame_emb_cols)
    frame_emb_df.insert(0, "frame", np.arange(1, t + 1, dtype=np.int32))

    # Object embeddings (masked object patches)
    obj_emb = np.zeros((len(kept_objs), 512), dtype=np.float32)
    with torch.no_grad():
        batch_imgs: List[np.ndarray] = []
        batch_idx: List[int] = []
        for i, obj in enumerate(kept_objs):
            patch = object_rgb_patch_masked(stack, labels_by_frame, obj)
            batch_imgs.append(patch)
            batch_idx.append(i)
            if len(batch_imgs) >= args.embed_batch_size:
                z = run_embed_patch_batches(embedder, batch_imgs, batch_size=args.embed_batch_size, device=device)
                obj_emb[np.array(batch_idx, dtype=np.int64), :] = z
                batch_imgs = []
                batch_idx = []
        if batch_imgs:
            z = run_embed_patch_batches(embedder, batch_imgs, batch_size=args.embed_batch_size, device=device)
            obj_emb[np.array(batch_idx, dtype=np.int64), :] = z

    obj_emb_cols = [f"obj_emb_{i:03d}" for i in range(obj_emb.shape[1])]
    obj_emb_df = pd.DataFrame(obj_emb, columns=obj_emb_cols)

    out_df = pd.concat([df.reset_index(drop=True), obj_emb_df], axis=1)
    out_df = out_df.merge(frame_emb_df, on="frame", how="left")
    out_df = maybe_merge_manifest(out_df, args.manifest_csv)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_feature_table)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_meta_json)), exist_ok=True)

    if args.output_feature_table.lower().endswith(".parquet"):
        out_df.to_parquet(args.output_feature_table, index=False)
    elif args.output_feature_table.lower().endswith(".csv"):
        out_df.to_csv(args.output_feature_table, index=False)
    else:
        raise ValueError("output_feature_table must end with .parquet or .csv")

    meta = {
        "input_stack_tif": os.path.abspath(args.input_stack_tif),
        "annotation_results_csv": os.path.abspath(args.annotation_results_csv),
        "manifest_csv": os.path.abspath(args.manifest_csv) if args.manifest_csv else None,
        "output_labels_tif": os.path.abspath(args.output_labels_tif),
        "output_feature_table": os.path.abspath(args.output_feature_table),
        "device": str(device),
        "cellpose_batch_size": int(args.cellpose_batch_size),
        "embed_batch_size": int(args.embed_batch_size),
        "point_snap_radius": int(args.point_snap_radius),
        "exclude_border_objects": bool(exclude_border),
        "resnet_weights": args.resnet_weights,
        "stack_shape_TCXY": [int(t), int(c), int(h), int(w)],
        "labels_tif_dtype": labels_dtype,
        "n_objects_total": int(len(objs)),
        "n_objects_kept": int(len(kept_objs)),
        "n_points_total": int(sum(len(v) for v in pts_by_frame.values())),
        "n_points_unmatched": int(unmatched_points),
        "label_counts": {
            "alive": int((out_df["label_int"] == 1).sum()),
            "not_alive": int((out_df["label_int"] == 0).sum()),
        },
        "feature_dims": {
            "obj_embedding": 512,
            "frame_embedding": 512,
            "total_numeric_embedding_plus_log_area": 1025,
        },
    }
    with open(args.output_meta_json, "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote labels stack: {os.path.abspath(args.output_labels_tif)}")
    print(f"Wrote feature table: {os.path.abspath(args.output_feature_table)}")
    print(f"Wrote run metadata: {os.path.abspath(args.output_meta_json)}")


if __name__ == "__main__":
    main()
