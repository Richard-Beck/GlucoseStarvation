#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import torch
from cellpose.models import CellposeModel
from scipy import ndimage
from skimage import filters, measure, morphology, segmentation


VALID_EXTENSIONS = (".tif", ".tiff")
CHANNEL_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$",
    re.IGNORECASE,
)
MAX_IMAGES_HARD_LIMIT = 1000
EPS = 1e-8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a bounded WP3 CPSAM cell-mask and nuclear-size pilot."
    )
    parser.add_argument("--manifest_csv", required=True, help="WP3 pilot manifest CSV.")
    parser.add_argument("--raw_data_dir", default="all_raw", help="Directory containing raw channel TIFF symlinks.")
    parser.add_argument(
        "--output_root",
        default="data/image_processing_runs/wp3_nuclear_size_pilot",
        help="Root directory under which a timestamped run folder is created.",
    )
    parser.add_argument("--run_name", default=None, help="Optional run folder name. Defaults to run_<timestamp>.")
    parser.add_argument("--max_images", type=int, default=1000, help="Extra safety cap on processed manifest rows.")
    parser.add_argument("--cellpose_batch_size", type=int, default=8, help="Cellpose/CPSAM internal batch size.")
    parser.add_argument("--min_cell_area_px", type=int, default=50, help="Objects smaller than this remain in masks but are marked below size gate.")
    parser.add_argument("--min_nucleus_area_px", type=int, default=10, help="Minimum nuclear connected-component area.")
    parser.add_argument("--qc_count", type=int, default=12, help="Number of representative QC PNG panels to write.")
    parser.add_argument("--save_masks", action="store_true", default=True, help="Save per-image cell and nuclear mask TIFFs.")
    parser.add_argument("--no_save_masks", action="store_false", dest="save_masks", help="Do not save per-image masks.")
    parser.add_argument("--allow_cpu", action="store_true", help="Allow running without CUDA. Intended only for debugging.")
    parser.add_argument("--dry_run", action="store_true", help="Validate inputs and write run metadata without segmentation.")
    return parser.parse_args()


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def make_run_dir(output_root: str, run_name: Optional[str]) -> str:
    name = run_name or f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_dir = os.path.join(output_root, name)
    os.makedirs(run_dir, exist_ok=False)
    return run_dir


def shlex_quote(text: str) -> str:
    if text == "":
        return "''"
    if all(ch.isalnum() or ch in "-_./=:" for ch in text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def read_manifest(path: str, max_images: int) -> list[dict]:
    if max_images < 1:
        raise ValueError("--max_images must be >= 1.")
    if max_images > MAX_IMAGES_HARD_LIMIT:
        raise ValueError("--max_images must be <= 1000 for the bounded WP3 pilot.")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_key", "cellLine", "ploidy", "G0", "glucose_bin", "time_bin", "hours"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")
        rows = list(reader)
    if len(rows) > max_images:
        print(f"Manifest has {len(rows)} rows; processing first {max_images} due to --max_images.", flush=True)
        rows = rows[:max_images]
    return rows


def parse_channel_filename(filename: str) -> Optional[tuple[str, str]]:
    match = CHANNEL_RE.match(filename)
    if match is None:
        return None
    channel = match.group("channel").lower()
    key = f"{match.group('prefix')}_{match.group('suffix')}"
    return key, channel


def build_raw_group_map(raw_data_dir: str, target_keys: set[str]) -> dict[str, dict[str, list[str]]]:
    groups: dict[str, dict[str, list[str]]] = {
        key: {"alive": [], "dead": [], "phase": []} for key in target_keys
    }
    with os.scandir(raw_data_dir) as it:
        for entry in it:
            if not entry.name.lower().endswith(VALID_EXTENSIONS):
                continue
            parsed = parse_channel_filename(entry.name)
            if parsed is None:
                continue
            key, channel = parsed
            if key not in groups:
                continue
            if channel.startswith("alive"):
                groups[key]["alive"].append(entry.path)
            elif channel == "dead":
                groups[key]["dead"].append(entry.path)
            elif channel == "phase":
                groups[key]["phase"].append(entry.path)
    return groups


def sum_channel(paths: list[str]) -> Optional[np.ndarray]:
    if not paths:
        return None
    imgs = [tifffile.imread(path).astype(np.float32) for path in sorted(paths)]
    ref_shape = imgs[0].shape
    for img in imgs[1:]:
        if img.shape != ref_shape:
            raise ValueError("Channel images do not share shape.")
    return np.sum(imgs, axis=0, dtype=np.float32)


def robust01(img: np.ndarray) -> np.ndarray:
    arr = np.asarray(img, dtype=np.float32)
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
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0).astype(np.float32, copy=False)


def load_channels(channel_paths: dict[str, list[str]]) -> dict[str, np.ndarray]:
    alive_raw = sum_channel(channel_paths["alive"])
    dead_raw = sum_channel(channel_paths["dead"])
    phase_raw = sum_channel(channel_paths["phase"])
    ref = next((x for x in (alive_raw, dead_raw, phase_raw) if x is not None), None)
    if ref is None:
        raise ValueError("No alive/dead/phase files found for image_key.")
    if alive_raw is None:
        alive_raw = np.zeros_like(ref, dtype=np.float32)
    if dead_raw is None:
        dead_raw = np.zeros_like(ref, dtype=np.float32)
    if phase_raw is None:
        phase_raw = np.zeros_like(ref, dtype=np.float32)
    return {
        "alive_raw": alive_raw,
        "dead_raw": dead_raw,
        "phase_raw": phase_raw,
        "alive01": robust01(alive_raw),
        "dead01": robust01(dead_raw),
        "phase01": robust01(phase_raw),
        "nuclear01": robust01(alive_raw + dead_raw),
    }


def segment_cell_masks(channels: dict[str, np.ndarray], model: CellposeModel, batch_size: int) -> np.ndarray:
    cpose_input = np.stack(
        [channels["phase01"], robust01(channels["alive_raw"] + channels["dead_raw"])],
        axis=0,
    )
    masks, _, _ = model.eval(cpose_input, diameter=None, channels=[2, 1], batch_size=batch_size)
    return np.asarray(masks)


def nuclear_mask_for_cell(
    nuclear_img: np.ndarray,
    cell_mask: np.ndarray,
    min_nucleus_area_px: int,
) -> np.ndarray:
    values = nuclear_img[cell_mask]
    if values.size < min_nucleus_area_px:
        return np.zeros_like(cell_mask, dtype=bool)
    if float(np.percentile(values, 95) - np.percentile(values, 5)) < 0.03:
        return np.zeros_like(cell_mask, dtype=bool)

    smooth = ndimage.gaussian_filter(nuclear_img, sigma=1.0)
    smooth_values = smooth[cell_mask]
    try:
        threshold = float(filters.threshold_otsu(smooth_values))
    except ValueError:
        threshold = float(np.percentile(smooth_values, 75))
    threshold = max(threshold, float(np.percentile(smooth_values, 60)))

    candidate = (smooth >= threshold) & cell_mask
    candidate = morphology.remove_small_objects(candidate, min_size=min_nucleus_area_px)
    candidate = morphology.binary_closing(candidate, morphology.disk(1))
    candidate = ndimage.binary_fill_holes(candidate)
    candidate = morphology.remove_small_objects(candidate, min_size=min_nucleus_area_px)
    return np.asarray(candidate, dtype=bool)


def measure_objects(
    cell_masks: np.ndarray,
    nuclear_img: np.ndarray,
    alive_img: np.ndarray,
    dead_img: np.ndarray,
    min_cell_area_px: int,
    min_nucleus_area_px: int,
) -> tuple[list[dict], np.ndarray]:
    object_rows: list[dict] = []
    nuclear_labels = np.zeros_like(cell_masks, dtype=cell_masks.dtype)
    props = measure.regionprops(cell_masks)
    for prop in props:
        oid = int(prop.label)
        cell_area = int(prop.area)
        y0, x0, y1, x1 = [int(v) for v in prop.bbox]
        local_cell_mask = prop.image.astype(bool, copy=False)
        local_nuclear_img = nuclear_img[y0:y1, x0:x1]
        local_alive_img = alive_img[y0:y1, x0:x1]
        local_dead_img = dead_img[y0:y1, x0:x1]

        nuc_mask = nuclear_mask_for_cell(local_nuclear_img, local_cell_mask, min_nucleus_area_px)
        nuc_area = int(np.sum(nuc_mask))
        if nuc_area > 0:
            nuclear_view = nuclear_labels[y0:y1, x0:x1]
            nuclear_view[nuc_mask] = oid
            nuc_labels, n_components = ndimage.label(nuc_mask)
            component_areas = np.bincount(nuc_labels.ravel())[1:]
            largest_nuc_area = int(component_areas.max()) if component_areas.size else 0
        else:
            n_components = 0
            largest_nuc_area = 0

        cy, cx = prop.centroid
        row = {
            "object_id": oid,
            "cell_area_px": cell_area,
            "cell_area_pass": int(cell_area >= min_cell_area_px),
            "centroid_x": f"{cx:.3f}",
            "centroid_y": f"{cy:.3f}",
            "bbox_x0": x0,
            "bbox_y0": y0,
            "bbox_x1": x1,
            "bbox_y1": y1,
            "nuclear_area_px": nuc_area,
            "largest_nuclear_area_px": largest_nuc_area,
            "nuclear_component_count": int(n_components),
            "nuclear_to_cell_area_ratio": f"{(nuc_area / max(cell_area, 1)):.8f}",
            "nuclear_mean_intensity": f"{float(np.mean(local_nuclear_img[nuc_mask])):.8f}" if nuc_area > 0 else "",
            "nuclear_integrated_intensity": f"{float(np.sum(local_nuclear_img[nuc_mask])):.8f}" if nuc_area > 0 else "",
            "cell_mean_nuclear_channel": f"{float(np.mean(local_nuclear_img[local_cell_mask])):.8f}",
            "cell_mean_alive_channel": f"{float(np.mean(local_alive_img[local_cell_mask])):.8f}",
            "cell_mean_dead_channel": f"{float(np.mean(local_dead_img[local_cell_mask])):.8f}",
        }
        object_rows.append(row)
    return object_rows, nuclear_labels


def overlay_boundaries(base: np.ndarray, cell_masks: np.ndarray, nuclear_masks: np.ndarray) -> np.ndarray:
    rgb = np.stack([base, base, base], axis=-1)
    cell_boundary = segmentation.find_boundaries(cell_masks, mode="inner")
    nuclear_boundary = segmentation.find_boundaries(nuclear_masks, mode="inner")
    rgb[cell_boundary, :] = np.array([0.0, 0.85, 1.0])
    rgb[nuclear_boundary, :] = np.array([1.0, 0.15, 0.15])
    return rgb


def write_qc_panel(
    path: str,
    manifest_row: dict,
    channels: dict[str, np.ndarray],
    cell_masks: np.ndarray,
    nuclear_masks: np.ndarray,
) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(14, 4), constrained_layout=True)
    axes[0].imshow(channels["phase01"], cmap="gray")
    axes[0].set_title("phase")
    axes[1].imshow(channels["nuclear01"], cmap="magma")
    axes[1].set_title("nuclear stain")
    axes[2].imshow(overlay_boundaries(channels["phase01"], cell_masks, nuclear_masks))
    axes[2].set_title("CPSAM cells + nuclei")
    axes[3].imshow(cell_masks, cmap="viridis", interpolation="nearest")
    axes[3].contour(nuclear_masks > 0, levels=[0.5], colors="red", linewidths=0.4)
    axes[3].set_title("label masks")
    for ax in axes:
        ax.set_axis_off()
    title = (
        f"{manifest_row['cellLine']} {manifest_row['ploidy']} "
        f"G0={manifest_row['G0']} {manifest_row['time_bin']} h={manifest_row['hours']}"
    )
    fig.suptitle(title, fontsize=10)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_run_provenance(args: argparse.Namespace, run_dir: str, manifest_rows: list[dict]) -> None:
    with open(os.path.join(run_dir, "invocation.txt"), "w", encoding="utf-8") as handle:
        handle.write(" ".join(shlex_quote(x) for x in sys.argv) + "\n")
    with open(os.path.join(run_dir, "run_args.json"), "w", encoding="utf-8") as handle:
        json.dump(vars(args), handle, indent=2, sort_keys=True)
        handle.write("\n")
    shutil.copy2(os.path.abspath(__file__), os.path.join(run_dir, os.path.basename(__file__)))
    shutil.copy2(os.path.abspath(args.manifest_csv), os.path.join(run_dir, "pilot_manifest.csv"))
    meta = {
        "created_at": datetime.now().isoformat(),
        "hostname": os.uname().nodename,
        "python_executable": sys.executable,
        "torch_version": getattr(torch, "__version__", None),
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "cellpose_backend": "CellposeModel / Cellpose v4 CPSAM",
        "n_manifest_rows": len(manifest_rows),
    }
    with open(os.path.join(run_dir, "run_meta.json"), "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")


def summarize_image(manifest_row: dict, object_rows: list[dict], status: str, error: str = "") -> dict:
    n_cells = len(object_rows)
    n_cells_with_nucleus = sum(1 for row in object_rows if int(row["nuclear_area_px"]) > 0)
    area_values = [int(row["cell_area_px"]) for row in object_rows]
    nuclear_values = [int(row["nuclear_area_px"]) for row in object_rows if int(row["nuclear_area_px"]) > 0]
    out = {
        "image_key": manifest_row["image_key"],
        "cellLine": manifest_row["cellLine"],
        "ploidy": manifest_row["ploidy"],
        "G0": manifest_row["G0"],
        "glucose_bin": manifest_row["glucose_bin"],
        "time_bin": manifest_row["time_bin"],
        "hours": manifest_row["hours"],
        "status": status,
        "error": error,
        "n_cells": n_cells,
        "n_cells_with_nucleus": n_cells_with_nucleus,
        "nuclear_success_fraction": f"{(n_cells_with_nucleus / max(n_cells, 1)):.8f}",
        "median_cell_area_px": f"{float(np.median(area_values)):.3f}" if area_values else "",
        "median_nuclear_area_px": f"{float(np.median(nuclear_values)):.3f}" if nuclear_values else "",
    }
    return out


def main() -> None:
    args = parse_args()
    manifest_rows = read_manifest(args.manifest_csv, args.max_images)
    run_dir = make_run_dir(args.output_root, args.run_name)
    for subdir in ("masks/cell", "masks/nuclear", "qc"):
        os.makedirs(os.path.join(run_dir, subdir), exist_ok=True)
    write_run_provenance(args, run_dir, manifest_rows)

    if not manifest_rows:
        raise RuntimeError("Manifest contains no rows after applying max_images cap.")
    if not args.allow_cpu and not args.dry_run and not torch.cuda.is_available():
        raise RuntimeError("CUDA is not available. Submit this script through the WP3 SLURM GPU wrapper or pass --allow_cpu only for debugging.")

    print(f"Run directory: {run_dir}", flush=True)
    print(f"Processing {len(manifest_rows)} manifest rows.", flush=True)
    print(f"CUDA available: {torch.cuda.is_available()}", flush=True)

    target_keys = {row["image_key"] for row in manifest_rows}
    raw_groups = build_raw_group_map(args.raw_data_dir, target_keys)
    missing = [key for key, groups in raw_groups.items() if not any(groups.values())]
    if missing:
        raise RuntimeError(f"No raw channel files found for {len(missing)} manifest image_key values. First missing: {missing[0]}")
    if args.dry_run:
        with open(os.path.join(run_dir, "summary.json"), "w", encoding="utf-8") as handle:
            json.dump({"dry_run": True, "n_manifest_rows": len(manifest_rows), "n_missing_raw_groups": len(missing)}, handle, indent=2)
            handle.write("\n")
        print("Dry run complete.", flush=True)
        return

    model = CellposeModel(gpu=True)
    object_csv = os.path.join(run_dir, "wp3_nuclear_object_features.csv")
    image_csv = os.path.join(run_dir, "wp3_nuclear_image_qc.csv")
    object_headers = [
        "image_key",
        "cellLine",
        "ploidy",
        "G0",
        "glucose_bin",
        "time_bin",
        "hours",
        "object_id",
        "cell_area_px",
        "cell_area_pass",
        "centroid_x",
        "centroid_y",
        "bbox_x0",
        "bbox_y0",
        "bbox_x1",
        "bbox_y1",
        "nuclear_area_px",
        "largest_nuclear_area_px",
        "nuclear_component_count",
        "nuclear_to_cell_area_ratio",
        "nuclear_mean_intensity",
        "nuclear_integrated_intensity",
        "cell_mean_nuclear_channel",
        "cell_mean_alive_channel",
        "cell_mean_dead_channel",
    ]
    image_headers = [
        "image_key",
        "cellLine",
        "ploidy",
        "G0",
        "glucose_bin",
        "time_bin",
        "hours",
        "status",
        "error",
        "n_cells",
        "n_cells_with_nucleus",
        "nuclear_success_fraction",
        "median_cell_area_px",
        "median_nuclear_area_px",
    ]

    total_cells = 0
    total_cells_with_nucleus = 0
    with open(object_csv, "w", newline="", encoding="utf-8") as object_handle, open(image_csv, "w", newline="", encoding="utf-8") as image_handle:
        object_writer = csv.DictWriter(object_handle, fieldnames=object_headers)
        image_writer = csv.DictWriter(image_handle, fieldnames=image_headers)
        object_writer.writeheader()
        image_writer.writeheader()

        for idx, row in enumerate(manifest_rows, start=1):
            key = row["image_key"]
            print(f"[{idx}/{len(manifest_rows)}] {key}", flush=True)
            try:
                channels = load_channels(raw_groups[key])
                cell_masks = segment_cell_masks(channels, model, args.cellpose_batch_size)
                if cell_masks.size == 0 or int(cell_masks.max()) == 0:
                    image_writer.writerow(summarize_image(row, [], "no_cells"))
                    continue
                object_rows, nuclear_masks = measure_objects(
                    cell_masks,
                    channels["nuclear01"],
                    channels["alive01"],
                    channels["dead01"],
                    args.min_cell_area_px,
                    args.min_nucleus_area_px,
                )
                for object_row in object_rows:
                    out = {name: row.get(name, "") for name in ("image_key", "cellLine", "ploidy", "G0", "glucose_bin", "time_bin", "hours")}
                    out.update(object_row)
                    object_writer.writerow(out)
                image_writer.writerow(summarize_image(row, object_rows, "ok"))
                total_cells += len(object_rows)
                total_cells_with_nucleus += sum(1 for object_row in object_rows if int(object_row["nuclear_area_px"]) > 0)

                if args.save_masks:
                    dtype = np.uint16 if int(cell_masks.max()) <= np.iinfo(np.uint16).max else np.uint32
                    tifffile.imwrite(os.path.join(run_dir, "masks", "cell", f"{safe_name(key)}_cell_masks.tif"), cell_masks.astype(dtype, copy=False))
                    tifffile.imwrite(os.path.join(run_dir, "masks", "nuclear", f"{safe_name(key)}_nuclear_masks.tif"), nuclear_masks.astype(dtype, copy=False))
                if idx <= args.qc_count:
                    write_qc_panel(os.path.join(run_dir, "qc", f"{idx:03d}_{safe_name(key)}.png"), row, channels, cell_masks, nuclear_masks)
                object_handle.flush()
                image_handle.flush()
            except Exception as exc:
                print(f"ERROR processing {key}: {exc}", flush=True)
                image_writer.writerow(summarize_image(row, [], "error", str(exc)))
                image_handle.flush()

    summary = {
        "n_images": len(manifest_rows),
        "n_cells": total_cells,
        "n_cells_with_nucleus": total_cells_with_nucleus,
        "nuclear_success_fraction": total_cells_with_nucleus / max(total_cells, 1),
        "object_features_csv": os.path.abspath(object_csv),
        "image_qc_csv": os.path.abspath(image_csv),
        "mask_dirs": {
            "cell": os.path.abspath(os.path.join(run_dir, "masks", "cell")),
            "nuclear": os.path.abspath(os.path.join(run_dir, "masks", "nuclear")),
        },
    }
    with open(os.path.join(run_dir, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Finished WP3 pilot. Wrote {object_csv} and {image_csv}", flush=True)


if __name__ == "__main__":
    main()
