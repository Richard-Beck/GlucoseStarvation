#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
from collections import defaultdict
from statistics import median
from typing import Optional

import numpy as np
import tifffile


VALID_EXTENSIONS = (".tif", ".tiff")
CHANNEL_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$",
    re.IGNORECASE,
)
EXCLUDED_CELLLINES = {"MCF10A-ctrl", "MCF10A-hras"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a 90-image annotation stack stratified by "
            "cellLine x ploidy x glucose_bin x time_bin."
        )
    )
    parser.add_argument("--prelim_counts_csv", required=True, help="Path to prelim_counts_by_image.csv")
    parser.add_argument("--raw_data_dir", required=True, help="Directory with raw channel TIFFs.")
    parser.add_argument("--output_stack_tif", required=True, help="Output ImageJ stack TIFF path.")
    parser.add_argument("--output_manifest_csv", required=True, help="Output manifest CSV path.")
    parser.add_argument("--output_meta_json", required=True, help="Output metadata JSON path.")
    parser.add_argument("--context_size", type=int, default=400, help="Context crop size in pixels.")
    parser.add_argument("--target_size", type=int, default=200, help="Target region size in pixels.")
    parser.add_argument("--border_px", type=int, default=2, help="Target-box border thickness.")
    return parser.parse_args()


def parse_channel_filename(filename: str) -> Optional[tuple[str, str]]:
    match = CHANNEL_RE.match(filename)
    if match is None:
        return None
    channel = match.group("channel").lower()
    base_key = f"{match.group('prefix')}_{match.group('suffix')}"
    return base_key, channel


def build_image_key_groups(raw_data_dir: str) -> dict[str, dict[str, list[str]]]:
    groups: dict[str, dict[str, list[str]]] = {}
    for fname in sorted(os.listdir(raw_data_dir)):
        if not fname.lower().endswith(VALID_EXTENSIONS):
            continue
        parsed = parse_channel_filename(fname)
        if parsed is None:
            continue
        key, channel = parsed
        channel_groups = groups.setdefault(key, {"alive": [], "dead": [], "phase": []})
        full_path = os.path.join(raw_data_dir, fname)
        if channel.startswith("alive"):
            channel_groups["alive"].append(full_path)
        elif channel == "dead":
            channel_groups["dead"].append(full_path)
        elif channel == "phase":
            channel_groups["phase"].append(full_path)
    return groups


def sum_channel(paths: list[str]) -> Optional[np.ndarray]:
    if not paths:
        return None
    imgs = [tifffile.imread(p).astype(np.float32) for p in paths]
    ref_shape = imgs[0].shape
    for img in imgs[1:]:
        if img.shape != ref_shape:
            raise ValueError("Channel images do not share shape.")
    return np.sum(imgs, axis=0, dtype=np.float32)


def robust_normalize(img: np.ndarray) -> np.ndarray:
    x = np.asarray(img, dtype=np.float32)
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


def make_norm_channels(channels: dict[str, list[str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alive_raw = sum_channel(channels["alive"])
    dead_raw = sum_channel(channels["dead"])
    phase_raw = sum_channel(channels["phase"])

    ref = next((arr for arr in (dead_raw, alive_raw, phase_raw) if arr is not None), None)
    if ref is None:
        raise ValueError("No alive/dead/phase images found for image_key.")

    if alive_raw is None:
        alive_raw = np.zeros_like(ref, dtype=np.float32)
    if dead_raw is None:
        dead_raw = np.zeros_like(ref, dtype=np.float32)
    if phase_raw is None:
        phase_raw = np.zeros_like(ref, dtype=np.float32)

    dead_norm = robust_normalize(dead_raw)
    alive_norm = robust_normalize(alive_raw)
    phase_norm = robust_normalize(phase_raw)

    maxv = 65535.0
    phase_u16 = np.clip(phase_norm * maxv, 0, maxv).astype(np.uint16)
    dead_u16 = np.clip(dead_norm * maxv, 0, maxv).astype(np.uint16)
    alive_u16 = np.clip(alive_norm * maxv, 0, maxv).astype(np.uint16)
    return phase_u16, dead_u16, alive_u16


def bounded_crop_start(center: int, size: int, max_size: int) -> int:
    start = int(center - size // 2)
    if start < 0:
        return 0
    if start + size > max_size:
        return max(0, max_size - size)
    return start


def draw_target_box(ch_stack: np.ndarray, target_size: int, border_px: int) -> None:
    _, h, w = ch_stack.shape
    if target_size > h or target_size > w:
        raise ValueError("target_size cannot exceed crop dimensions.")
    y0 = (h - target_size) // 2
    y1 = y0 + target_size
    x0 = (w - target_size) // 2
    x1 = x0 + target_size
    b = max(1, int(border_px))

    ch_stack[0, y0:y0 + b, x0:x1] = 0
    ch_stack[0, y1 - b:y1, x0:x1] = 0
    ch_stack[0, y0:y1, x0:x0 + b] = 0
    ch_stack[0, y0:y1, x1 - b:x1] = 0
    ch_stack[1, y0:y0 + b, x0:x1] = 65535
    ch_stack[1, y1 - b:y1, x0:x1] = 65535
    ch_stack[1, y0:y1, x0:x0 + b] = 65535
    ch_stack[1, y0:y1, x1 - b:x1] = 65535
    ch_stack[2, y0:y0 + b, x0:x1] = 65535
    ch_stack[2, y1 - b:y1, x0:x1] = 65535
    ch_stack[2, y0:y1, x0:x0 + b] = 65535
    ch_stack[2, y0:y1, x1 - b:x1] = 65535


def get_imagej_luts_phase_dead_alive() -> list[np.ndarray]:
    lut_gray = np.zeros((3, 256), dtype=np.uint8)
    lut_gray[0, :] = np.arange(256, dtype=np.uint8)
    lut_gray[1, :] = np.arange(256, dtype=np.uint8)
    lut_gray[2, :] = np.arange(256, dtype=np.uint8)
    lut_red = np.zeros((3, 256), dtype=np.uint8)
    lut_red[0, :] = np.arange(256, dtype=np.uint8)
    lut_green = np.zeros((3, 256), dtype=np.uint8)
    lut_green[1, :] = np.arange(256, dtype=np.uint8)
    return [lut_gray, lut_red, lut_green]


def parse_glucose_bin(g0_value: str) -> Optional[str]:
    try:
        g0 = float(g0_value)
    except ValueError:
        return None
    if abs(g0 - 0.0) < 1e-9 or abs(g0 - 0.1) < 1e-9:
        return "low"
    if abs(g0 - 0.25) < 1e-9 or abs(g0 - 0.5) < 1e-9 or abs(g0 - 1.0) < 1e-9:
        return "mid"
    if abs(g0 - 5.0) < 1e-9 or abs(g0 - 25.0) < 1e-9:
        return "high"
    return None


def assign_time_bins(rows: list[dict]) -> None:
    rows_sorted = sorted(rows, key=lambda r: (float(r["hours"]), r["image_key"]))
    n = len(rows_sorted)
    for idx, row in enumerate(rows_sorted):
        row["time_bin_idx"] = min(2, (idx * 3) // n)
        row["time_bin"] = ("early", "mid", "late")[row["time_bin_idx"]]


def load_candidates(prelim_counts_csv: str) -> list[dict]:
    candidates = []
    with open(prelim_counts_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cell_line = row["cellLine"]
            if cell_line in EXCLUDED_CELLLINES:
                continue
            gbin = parse_glucose_bin(row["G0"])
            if gbin is None:
                continue
            total = int(float(row["total_count_kept"]))
            if total <= 0:
                continue
            candidates.append(
                {
                    "cellLine": cell_line,
                    "ploidy": row["ploidy"],
                    "G0": float(row["G0"]),
                    "glucose_bin": gbin,
                    "hours": float(row["hours"]),
                    "image_key": row["image_key"],
                    "base_key": row["base_key"],
                    "total_count_kept": total,
                    "dead_count": int(float(row["dead_count"])),
                    "alive_count": int(float(row["alive_count"])),
                }
            )
    return candidates


def select_90_images(candidates: list[dict]) -> list[dict]:
    by_core = defaultdict(list)
    for row in candidates:
        key = (row["cellLine"], row["ploidy"], row["glucose_bin"])
        by_core[key].append(row)

    for group_rows in by_core.values():
        assign_time_bins(group_rows)

    strata = defaultdict(list)
    for row in candidates:
        skey = (row["cellLine"], row["ploidy"], row["glucose_bin"], row["time_bin"])
        strata[skey].append(row)

    expected_line_ploidy = sorted({(r["cellLine"], r["ploidy"]) for r in candidates})
    expected_bins = ["low", "mid", "high"]
    expected_times = ["early", "mid", "late"]

    selected = []
    missing = []
    for cell_line, ploidy in expected_line_ploidy:
        for gbin in expected_bins:
            for tbin in expected_times:
                skey = (cell_line, ploidy, gbin, tbin)
                rows = strata.get(skey, [])
                if not rows:
                    missing.append(skey)
                    continue
                center_hour = median([r["hours"] for r in rows])
                pick = sorted(
                    rows,
                    key=lambda r: (
                        abs(r["hours"] - center_hour),
                        -r["total_count_kept"],
                        r["image_key"],
                    ),
                )[0]
                selected.append(pick)

    if missing:
        preview = "; ".join([f"{m[0]}|{m[1]}|{m[2]}|{m[3]}" for m in missing[:10]])
        raise RuntimeError(
            f"Missing strata ({len(missing)} total). First few: {preview}"
        )
    if len(selected) != 90:
        raise RuntimeError(f"Expected 90 selected images but got {len(selected)}")

    return sorted(selected, key=lambda r: (r["cellLine"], r["ploidy"], r["glucose_bin"], r["time_bin"], r["image_key"]))


def main() -> None:
    args = parse_args()
    if args.target_size > args.context_size:
        raise ValueError("--target_size must be <= --context_size")

    candidates = load_candidates(args.prelim_counts_csv)
    selected = select_90_images(candidates)

    groups = build_image_key_groups(args.raw_data_dir)

    stack_frames = []
    manifest_rows = []
    for idx, row in enumerate(selected, start=1):
        image_key = row["image_key"]
        if image_key not in groups:
            raise KeyError(f"image_key missing from raw_data_dir: {image_key}")

        phase_u16, dead_u16, alive_u16 = make_norm_channels(groups[image_key])
        src_h, src_w = phase_u16.shape
        cx = src_w // 2
        cy = src_h // 2
        x0 = bounded_crop_start(cx, args.context_size, src_w)
        y0 = bounded_crop_start(cy, args.context_size, src_h)
        x1 = x0 + args.context_size
        y1 = y0 + args.context_size

        crop = np.stack(
            [
                phase_u16[y0:y1, x0:x1],
                dead_u16[y0:y1, x0:x1],
                alive_u16[y0:y1, x0:x1],
            ],
            axis=0,
        ).copy()
        if crop.shape[1] != args.context_size or crop.shape[2] != args.context_size:
            raise RuntimeError(f"Bad crop shape for {image_key}: {crop.shape}")
        draw_target_box(crop, args.target_size, args.border_px)
        stack_frames.append(crop)

        target_x0 = x0 + (args.context_size - args.target_size) // 2
        target_y0 = y0 + (args.context_size - args.target_size) // 2
        manifest_rows.append(
            {
                "stack_index_1based": idx,
                "image_key": image_key,
                "cellLine": row["cellLine"],
                "ploidy": row["ploidy"],
                "G0": row["G0"],
                "glucose_bin": row["glucose_bin"],
                "time_bin": row["time_bin"],
                "hours": row["hours"],
                "base_key": row["base_key"],
                "crop_x0": x0,
                "crop_y0": y0,
                "crop_x1": x1,
                "crop_y1": y1,
                "target_x0": target_x0,
                "target_y0": target_y0,
                "target_x1": target_x0 + args.target_size,
                "target_y1": target_y0 + args.target_size,
            }
        )

    stack = np.stack(stack_frames, axis=0).astype(np.uint16, copy=False)  # T,C,Y,X

    os.makedirs(os.path.dirname(os.path.abspath(args.output_stack_tif)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_manifest_csv)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_meta_json)), exist_ok=True)

    tifffile.imwrite(
        args.output_stack_tif,
        stack,
        imagej=True,
        metadata={
            "axes": "TCYX",
            "mode": "composite",
            "LUTs": get_imagej_luts_phase_dead_alive(),
        },
    )

    manifest_fields = [
        "stack_index_1based",
        "image_key",
        "cellLine",
        "ploidy",
        "G0",
        "glucose_bin",
        "time_bin",
        "hours",
        "base_key",
        "crop_x0",
        "crop_y0",
        "crop_x1",
        "crop_y1",
        "target_x0",
        "target_y0",
        "target_x1",
        "target_y1",
    ]
    with open(args.output_manifest_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    meta = {
        "n_images": len(manifest_rows),
        "stratification": {
            "cellLine": sorted({m["cellLine"] for m in manifest_rows}),
            "ploidy": sorted({m["ploidy"] for m in manifest_rows}),
            "glucose_bins": ["low", "mid", "high"],
            "time_bins": ["early", "mid", "late"],
            "excluded_cell_lines": sorted(EXCLUDED_CELLLINES),
        },
        "crop_policy": {
            "context_size_px": args.context_size,
            "target_size_px": args.target_size,
            "annotation_rule": "Annotate objects whose center lies inside target box.",
            "crop_center_policy": "image center",
        },
        "imagej_channels": {
            "channel_order": ["phase", "dead", "alive"],
            "luts": {"phase": "grayscale", "dead": "red", "alive": "green"},
            "axes": "TCYX",
            "dtype": "uint16",
        },
        "outputs": {
            "stack_tif": os.path.abspath(args.output_stack_tif),
            "manifest_csv": os.path.abspath(args.output_manifest_csv),
            "meta_json": os.path.abspath(args.output_meta_json),
        },
    }
    with open(args.output_meta_json, "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote stack: {os.path.abspath(args.output_stack_tif)}")
    print(f"Wrote manifest: {os.path.abspath(args.output_manifest_csv)}")
    print(f"Wrote metadata: {os.path.abspath(args.output_meta_json)}")


if __name__ == "__main__":
    main()
