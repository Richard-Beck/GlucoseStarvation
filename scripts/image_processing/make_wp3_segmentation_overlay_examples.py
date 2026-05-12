#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries


VALID_EXTENSIONS = (".tif", ".tiff")
CHANNEL_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make RGB overlay examples for WP3 CPSAM cell and nuclear segmentations."
    )
    parser.add_argument(
        "--run_dir",
        default="data/image_processing_runs/wp3_nuclear_size_pilot/run_20260511_143515",
        help="Completed WP3 nuclear-size run directory.",
    )
    parser.add_argument("--raw_data_dir", default="all_raw", help="Directory with raw channel TIFF symlinks.")
    parser.add_argument(
        "--output_dir",
        default="figures/wp3_nuclear_size/segmentation_examples",
        help="Directory for output PNG overlays.",
    )
    parser.add_argument("--max_examples", type=int, default=10, help="Maximum examples to write.")
    parser.add_argument(
        "--image_keys",
        default="",
        help="Optional comma-separated image_key list. If omitted, choose a balanced set from the run manifest.",
    )
    return parser.parse_args()


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


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


def parse_channel_filename(filename: str) -> Optional[tuple[str, str]]:
    match = CHANNEL_RE.match(filename)
    if match is None:
        return None
    channel = match.group("channel").lower()
    key = f"{match.group('prefix')}_{match.group('suffix')}"
    return key, channel


def build_raw_group_map(raw_data_dir: str, target_keys: set[str]) -> dict[str, dict[str, list[str]]]:
    groups = {key: {"alive": [], "dead": [], "phase": []} for key in target_keys}
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


def load_raw_rgb(channel_paths: dict[str, list[str]]) -> np.ndarray:
    alive = sum_channel(channel_paths["alive"])
    dead = sum_channel(channel_paths["dead"])
    phase = sum_channel(channel_paths["phase"])

    ref = next((x for x in (phase, alive, dead) if x is not None), None)
    if ref is None:
        raise ValueError("No raw channel files found.")
    if alive is None:
        alive = np.zeros_like(ref, dtype=np.float32)
    if dead is None:
        dead = np.zeros_like(ref, dtype=np.float32)
    if phase is None:
        phase = np.zeros_like(ref, dtype=np.float32)

    phase01 = robust01(phase)
    alive01 = robust01(alive)
    dead01 = robust01(dead)

    # Phase carries morphology; fluorescence channels keep the Incucyte live/dead signal visible.
    rgb = np.stack(
        [
            np.clip(0.72 * phase01 + 0.55 * dead01, 0, 1),
            np.clip(0.72 * phase01 + 0.55 * alive01, 0, 1),
            np.clip(0.72 * phase01 + 0.18 * alive01, 0, 1),
        ],
        axis=-1,
    )
    return rgb


def read_manifest(run_dir: str) -> list[dict]:
    manifest_path = os.path.join(run_dir, "pilot_manifest.csv")
    with open(manifest_path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def choose_examples(rows: list[dict], max_examples: int) -> list[dict]:
    selected = []
    seen = set()
    for row in rows:
        key = (row["cellLine"], row["ploidy"])
        if key in seen:
            continue
        selected.append(row)
        seen.add(key)
        if len(selected) >= max_examples:
            return selected

    for row in rows:
        if row in selected:
            continue
        selected.append(row)
        if len(selected) >= max_examples:
            break
    return selected


def load_mask(path: str) -> np.ndarray:
    arr = tifffile.imread(path)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D mask at {path}, got shape {arr.shape}")
    return arr


def overlay_boundaries(rgb: np.ndarray, cell_mask: np.ndarray, nuclear_mask: np.ndarray) -> np.ndarray:
    out = rgb.copy()
    cell_boundary = find_boundaries(cell_mask, mode="inner")
    nuclear_boundary = find_boundaries(nuclear_mask, mode="inner")

    # Thicken by one pixel for visibility in full-frame exports.
    for boundary, color in (
        (cell_boundary, np.array([0.0, 0.95, 1.0])),
        (nuclear_boundary, np.array([1.0, 0.0, 0.85])),
    ):
        yy, xx = np.where(boundary)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                y2 = np.clip(yy + dy, 0, out.shape[0] - 1)
                x2 = np.clip(xx + dx, 0, out.shape[1] - 1)
                out[y2, x2, :] = color
    return out


def write_one_panel(row: dict, raw_groups: dict[str, dict[str, list[str]]], run_dir: str, output_dir: str, index: int) -> str:
    key = row["image_key"]
    stem = safe_name(key)
    cell_path = os.path.join(run_dir, "masks", "cell", f"{stem}_cell_masks.tif")
    nuclear_path = os.path.join(run_dir, "masks", "nuclear", f"{stem}_nuclear_masks.tif")
    rgb = load_raw_rgb(raw_groups[key])
    cell_mask = load_mask(cell_path)
    nuclear_mask = load_mask(nuclear_path)
    overlay = overlay_boundaries(rgb, cell_mask, nuclear_mask)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    axes[0].imshow(rgb)
    axes[0].set_title("RGB input composite")
    axes[1].imshow(overlay)
    axes[1].set_title("CPSAM cells (cyan) + nuclei (magenta)")
    for ax in axes:
        ax.set_axis_off()
    fig.suptitle(
        f"{row['cellLine']} {row['ploidy']} G0={row['G0']} {row['time_bin']} h={row['hours']}",
        fontsize=11,
    )
    out_path = os.path.join(output_dir, f"{index:02d}_{row['cellLine']}_{row['ploidy']}_{row['glucose_bin']}_{row['time_bin']}.png")
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    rows = read_manifest(args.run_dir)
    if args.image_keys.strip():
        requested = {x.strip() for x in args.image_keys.split(",") if x.strip()}
        examples = [row for row in rows if row["image_key"] in requested]
    else:
        examples = choose_examples(rows, args.max_examples)
    examples = examples[: args.max_examples]
    if not examples:
        raise RuntimeError("No examples selected.")

    target_keys = {row["image_key"] for row in examples}
    raw_groups = build_raw_group_map(args.raw_data_dir, target_keys)
    missing = [key for key, groups in raw_groups.items() if not any(groups.values())]
    if missing:
        raise RuntimeError(f"Missing raw channels for {len(missing)} image_key values. First missing: {missing[0]}")

    written = []
    for idx, row in enumerate(examples, start=1):
        written.append(write_one_panel(row, raw_groups, args.run_dir, args.output_dir, idx))

    manifest_path = os.path.join(args.output_dir, "overlay_examples_manifest.csv")
    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        fieldnames = ["overlay_path", "image_key", "cellLine", "ploidy", "G0", "glucose_bin", "time_bin", "hours"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for out_path, row in zip(written, examples):
            writer.writerow({name: row.get(name, "") for name in fieldnames if name != "overlay_path"} | {"overlay_path": out_path})

    print(f"Wrote {len(written)} overlay examples to {args.output_dir}")
    print(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
