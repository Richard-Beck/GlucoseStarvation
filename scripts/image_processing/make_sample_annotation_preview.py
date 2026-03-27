#!/usr/bin/env python3
import argparse
import json
import os
import re
from typing import Optional

import numpy as np
import tifffile


VALID_EXTENSIONS = (".tif", ".tiff")
CHANNEL_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one annotation preview image with a 400x400 context crop and "
            "central 200x200 target region marked."
        )
    )
    parser.add_argument("--raw_data_dir", required=True, help="Directory with raw channel TIFFs.")
    parser.add_argument("--image_key", required=True, help="Image key (without channel token).")
    parser.add_argument(
        "--output_image",
        required=True,
        help="Output preview image path (TIFF).",
    )
    parser.add_argument(
        "--output_meta_json",
        required=True,
        help="Output JSON path with crop/target coordinates.",
    )
    parser.add_argument("--center_x", type=int, default=None, help="Optional crop center x in source image.")
    parser.add_argument("--center_y", type=int, default=None, help="Optional crop center y in source image.")
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
    # ch_stack expected shape: (C, H, W), channel order [phase, dead, alive]
    _, h, w = ch_stack.shape
    if target_size > h or target_size > w:
        raise ValueError("target_size cannot exceed crop dimensions.")
    y0 = (h - target_size) // 2
    y1 = y0 + target_size
    x0 = (w - target_size) // 2
    x1 = x0 + target_size

    b = max(1, int(border_px))

    # Draw yellow under LUT mapping by setting dead+alive high and phase low.
    # Channel order is [phase, dead, alive].
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
    # Channel 1 phase -> grayscale
    lut_gray = np.zeros((3, 256), dtype=np.uint8)
    lut_gray[0, :] = np.arange(256, dtype=np.uint8)
    lut_gray[1, :] = np.arange(256, dtype=np.uint8)
    lut_gray[2, :] = np.arange(256, dtype=np.uint8)
    # Channel 2 dead -> red
    lut_red = np.zeros((3, 256), dtype=np.uint8)
    lut_red[0, :] = np.arange(256, dtype=np.uint8)
    # Channel 3 alive -> green
    lut_green = np.zeros((3, 256), dtype=np.uint8)
    lut_green[1, :] = np.arange(256, dtype=np.uint8)
    return [lut_gray, lut_red, lut_green]


def main() -> None:
    args = parse_args()
    if args.target_size > args.context_size:
        raise ValueError("--target_size must be <= --context_size")

    groups = build_image_key_groups(args.raw_data_dir)
    if args.image_key not in groups:
        raise KeyError(f"image_key not found in raw_data_dir: {args.image_key}")

    phase_u8, dead_u8, alive_u8 = make_norm_channels(groups[args.image_key])
    src_h, src_w = phase_u8.shape

    cx = args.center_x if args.center_x is not None else src_w // 2
    cy = args.center_y if args.center_y is not None else src_h // 2

    x0 = bounded_crop_start(cx, args.context_size, src_w)
    y0 = bounded_crop_start(cy, args.context_size, src_h)
    x1 = x0 + args.context_size
    y1 = y0 + args.context_size

    crop = np.stack(
        [
            phase_u8[y0:y1, x0:x1],
            dead_u8[y0:y1, x0:x1],
            alive_u8[y0:y1, x0:x1],
        ],
        axis=0,
    ).copy()
    if crop.shape[1] != args.context_size or crop.shape[2] != args.context_size:
        raise RuntimeError("Computed crop shape is not the requested context size.")

    draw_target_box(crop, args.target_size, args.border_px)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_image)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_meta_json)), exist_ok=True)

    tifffile.imwrite(
        args.output_image,
        crop,
        imagej=True,
        metadata={
            "axes": "CYX",
            "mode": "composite",
            "LUTs": get_imagej_luts_phase_dead_alive(),
        },
    )

    target_x0 = x0 + (args.context_size - args.target_size) // 2
    target_y0 = y0 + (args.context_size - args.target_size) // 2
    meta = {
        "image_key": args.image_key,
        "source_shape_hw": [int(src_h), int(src_w)],
        "requested_center_xy": [int(cx), int(cy)],
        "context_size_px": int(args.context_size),
        "target_size_px": int(args.target_size),
        "crop_bounds_xyxy": [int(x0), int(y0), int(x1), int(y1)],
        "target_bounds_xyxy_in_source": [
            int(target_x0),
            int(target_y0),
            int(target_x0 + args.target_size),
            int(target_y0 + args.target_size),
        ],
        "annotation_rule": "Annotate objects whose center lies inside target_bounds_xyxy_in_source.",
        "channel_order": ["phase", "dead", "alive"],
        "luts": {"phase": "grayscale", "dead": "red", "alive": "green"},
        "output_image": os.path.abspath(args.output_image),
    }

    with open(args.output_meta_json, "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote preview image: {os.path.abspath(args.output_image)}")
    print(f"Wrote metadata: {os.path.abspath(args.output_meta_json)}")


if __name__ == "__main__":
    main()
