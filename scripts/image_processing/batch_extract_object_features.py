#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import shutil
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import numpy as np
import torch
import tifffile

from cellpose.models import CellposeModel


MAX_WORKERS = 6
BG_EPS = 1e-8
VALID_EXTENSIONS = (".tif", ".tiff")
CHANNEL_RE = re.compile(r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$", re.IGNORECASE)
csv_writer_lock = threading.Lock()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Cellpose on raw images and write per-object area/live/dead features."
    )
    parser.add_argument("--raw_data_dir", required=True, help="Directory containing the raw image files to process.")
    parser.add_argument(
        "--output_root",
        required=True,
        help="Root directory under which a timestamped run folder will be created.",
    )
    parser.add_argument("--max_workers", type=int, default=MAX_WORKERS, help="Thread pool size.")
    parser.add_argument(
        "--sample_fraction",
        type=float,
        default=1.0,
        help="Optional fraction of image keys to process for smoke tests.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed used only when sample_fraction < 1.",
    )
    return parser.parse_args()


def make_run_dir(output_root: str) -> str:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(output_root, f"run_{stamp}")
    os.makedirs(run_dir, exist_ok=False)
    return run_dir


def shlex_quote(text: str) -> str:
    if text == "":
        return "''"
    if all(ch.isalnum() or ch in "-_./=:" for ch in text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def write_run_provenance(args: argparse.Namespace, run_dir: str) -> None:
    cmd = " ".join([shlex_quote(x) for x in sys.argv])
    with open(os.path.join(run_dir, "invocation.txt"), "w", encoding="utf-8") as handle:
        handle.write(cmd + "\n")

    with open(os.path.join(run_dir, "run_args.json"), "w", encoding="utf-8") as handle:
        json.dump(vars(args), handle, indent=2, sort_keys=True)
        handle.write("\n")

    meta = {
        "created_at": datetime.now().isoformat(),
        "hostname": os.uname().nodename,
        "python_executable": sys.executable,
        "torch_version": getattr(torch, "__version__", None),
        "cuda_available": bool(torch.cuda.is_available()),
    }
    with open(os.path.join(run_dir, "run_meta.json"), "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")

    shutil.copy2(os.path.abspath(__file__), os.path.join(run_dir, os.path.basename(__file__)))


def safe_channel_bg_stats(channel_img: np.ndarray, bg_mask: np.ndarray) -> tuple[float, float]:
    bg_values = channel_img[bg_mask]
    if bg_values.size == 0:
        bg_values = channel_img.ravel()
    bg_mean = float(np.mean(bg_values))
    bg_sd = float(np.std(bg_values))
    return bg_mean, max(bg_sd, BG_EPS)


def robust_normalize(channel_img: np.ndarray) -> np.ndarray:
    x = np.asarray(channel_img, dtype=np.float32)
    finite = np.isfinite(x)
    if not np.any(finite):
        return np.zeros_like(x, dtype=np.float32)
    vals = x[finite]
    lo = float(np.percentile(vals, 1))
    hi = float(np.percentile(vals, 99))
    if hi <= lo:
        hi = float(np.max(vals))
        lo = float(np.min(vals))
    if hi <= lo:
        return np.zeros_like(x, dtype=np.float32)
    out = (x - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0).astype(np.float32, copy=False)


def parse_channel_filename(filename: str) -> tuple[str, str] | None:
    match = CHANNEL_RE.match(filename)
    if match is None:
        return None
    channel = match.group("channel").lower()
    base_key = f"{match.group('prefix')}_{match.group('suffix')}"
    return base_key, channel


def build_raw_group_map(raw_data_dir: str) -> dict[str, dict[str, list[str]]]:
    raw_groups: dict[str, dict[str, list[str]]] = {}
    for fname in sorted(os.listdir(raw_data_dir)):
        if not fname.lower().endswith(VALID_EXTENSIONS):
            continue
        parsed = parse_channel_filename(fname)
        if parsed is None:
            continue
        base_key, channel = parsed
        channel_groups = raw_groups.setdefault(base_key, {"alive": [], "dead": [], "phase": []})
        full_path = os.path.join(raw_data_dir, fname)
        if channel.startswith("alive"):
            channel_groups["alive"].append(full_path)
        elif channel == "dead":
            channel_groups["dead"].append(full_path)
        elif channel == "phase":
            channel_groups["phase"].append(full_path)
    return raw_groups


def sum_channel_images(file_paths: list[str]) -> np.ndarray | None:
    if not file_paths:
        return None
    imgs = [tifffile.imread(path).astype(np.float32) for path in file_paths]
    ref_shape = imgs[0].shape
    for img in imgs[1:]:
        if img.shape != ref_shape:
            raise ValueError("Channel images do not share the same shape.")
    return np.sum(imgs, axis=0, dtype=np.float32)


def build_composite_and_cpose_input(channel_paths: dict[str, list[str]]) -> tuple[np.ndarray, np.ndarray]:
    alive_raw = sum_channel_images(channel_paths["alive"])
    dead_raw = sum_channel_images(channel_paths["dead"])
    phase_raw = sum_channel_images(channel_paths["phase"])

    ref = next((arr for arr in (dead_raw, alive_raw, phase_raw) if arr is not None), None)
    if ref is None:
        raise ValueError("No alive/dead/phase images were found for this key.")

    if dead_raw is None:
        dead_raw = np.zeros_like(ref, dtype=np.float32)
    if alive_raw is None:
        alive_raw = np.zeros_like(ref, dtype=np.float32)
    if phase_raw is None:
        phase_raw = np.zeros_like(ref, dtype=np.float32)

    dead_norm = robust_normalize(dead_raw)
    alive_norm = robust_normalize(alive_raw)
    phase_norm = robust_normalize(phase_raw)
    fluor_sum_norm = robust_normalize(alive_raw + dead_raw)

    composite = np.stack([dead_norm, alive_norm, phase_norm], axis=-1)
    cpose_input = np.stack([phase_norm, fluor_sum_norm], axis=0)
    return composite, cpose_input


def process_image_key(
    key: str,
    raw_groups: dict,
    cellpose_model: CellposeModel,
    csv_writer: csv.DictWriter,
) -> int:
    try:
        print(f"Processing: {key}", flush=True)

        composite_img, cpose_input_img = build_composite_and_cpose_input(raw_groups[key])
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None, channels=[2, 1])

        if masks is None or masks.size == 0 or masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.", flush=True)
            return 0

        dead_img = composite_img[:, :, 0].astype(np.float32, copy=False)
        live_img = composite_img[:, :, 1].astype(np.float32, copy=False)

        counts = np.bincount(masks.ravel())
        bg_mask = masks == 0
        live_bg_mean, live_bg_sd = safe_channel_bg_stats(live_img, bg_mask)
        dead_bg_mean, dead_bg_sd = safe_channel_bg_stats(dead_img, bg_mask)

        rows_to_write = []
        max_obj_id = int(masks.max())
        for obj_id in range(1, max_obj_id + 1):
            obj_mask = masks == obj_id
            area_px = int(counts[obj_id]) if obj_id < len(counts) else int(np.sum(obj_mask))
            if area_px == 0:
                continue

            live_mean = float(np.mean(live_img[obj_mask]))
            dead_mean = float(np.mean(dead_img[obj_mask]))
            rows_to_write.append(
                {
                    "image_key": key,
                    "object_id": obj_id,
                    "area_px": area_px,
                    "live_mean": f"{live_mean:.8f}",
                    "dead_mean": f"{dead_mean:.8f}",
                    "live_bg_z": f"{(live_mean - live_bg_mean) / live_bg_sd:.8f}",
                    "dead_bg_z": f"{(dead_mean - dead_bg_mean) / dead_bg_sd:.8f}",
                }
            )

        if rows_to_write:
            with csv_writer_lock:
                csv_writer.writerows(rows_to_write)

        return len(rows_to_write)

    except Exception as exc:
        print(f"ERROR processing key {key}: {exc}", flush=True)
        return 0


def main() -> None:
    args = parse_args()
    if not (0 < args.sample_fraction <= 1.0):
        raise ValueError("--sample_fraction must be in (0, 1].")

    run_dir = make_run_dir(args.output_root)
    write_run_provenance(args, run_dir)
    output_csv_path = os.path.join(run_dir, "object_features.csv")

    print(f"Run directory: {run_dir}", flush=True)
    print("--- Initializing Cellpose model ---", flush=True)
    cellpose_model = CellposeModel(gpu=torch.cuda.is_available())

    print("\n--- Scanning for raw image groups ---", flush=True)
    raw_groups = build_raw_group_map(args.raw_data_dir)
    image_keys = sorted(raw_groups.keys())
    if not image_keys:
        raise RuntimeError(f"No raw images found in '{args.raw_data_dir}'.")

    if args.sample_fraction < 1.0:
        rng = np.random.default_rng(args.seed)
        sample_n = max(1, int(len(image_keys) * args.sample_fraction))
        image_keys = sorted(rng.choice(image_keys, size=sample_n, replace=False).tolist())
        print(
            f"Sampling {len(image_keys)} / {len(raw_groups)} image keys "
            f"(fraction={args.sample_fraction:.4f}, seed={args.seed}).",
            flush=True,
        )
    else:
        print(f"Processing all {len(image_keys)} image keys.", flush=True)

    headers = [
        "image_key",
        "object_id",
        "area_px",
        "live_mean",
        "dead_mean",
        "live_bg_z",
        "dead_bg_z",
    ]

    total_objects_processed = 0
    print(f"\n--- Starting parallel feature extraction using {args.max_workers} workers ---", flush=True)
    with open(output_csv_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures = [
                executor.submit(process_image_key, key, raw_groups, cellpose_model, writer)
                for key in image_keys
            ]
            for future in futures:
                total_objects_processed += future.result()

    summary = {
        "n_image_keys": len(image_keys),
        "n_objects": total_objects_processed,
        "output_csv_path": output_csv_path,
    }
    with open(os.path.join(run_dir, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print("\n--- Feature extraction complete ---", flush=True)
    print(f"Successfully processed {total_objects_processed} objects.", flush=True)
    print(f"Results saved to '{output_csv_path}'", flush=True)


if __name__ == "__main__":
    main()
