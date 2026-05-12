#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from datetime import datetime
from statistics import median
from typing import Iterable


EXCLUDED_CELLLINES_DEFAULT = ("MCF10A-ctrl", "MCF10A-hras")
GLUCOSE_BINS_DEFAULT = ("low", "mid", "high")
TIME_BINS = ("early", "mid", "late")
MAX_IMAGES_HARD_LIMIT = 1000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the bounded WP3 nuclear-size pilot manifest from prelim_counts_by_image.csv. "
            "The default selects one image per cellLine x ploidy x glucose_bin x time_bin stratum."
        )
    )
    parser.add_argument(
        "--prelim_counts_csv",
        default="data/image_processing_runs/run_20260324_233122/prelim_counts_by_image.csv",
        help="Input prelim_counts_by_image.csv with image metadata and per-image counts.",
    )
    parser.add_argument(
        "--output_csv",
        default="data/report_exports/wp3_nuclear_size/pilot_manifest.csv",
        help="Output manifest CSV.",
    )
    parser.add_argument(
        "--output_meta_json",
        default="data/report_exports/wp3_nuclear_size/pilot_manifest_meta.json",
        help="Output metadata JSON.",
    )
    parser.add_argument(
        "--images_per_stratum",
        type=int,
        default=1,
        help="Number of images to select per stratum before the max_images cap is applied.",
    )
    parser.add_argument(
        "--max_images",
        type=int,
        default=90,
        help="Hard cap on selected images. Must be <= 1000 unless --allow_over_1000 is used.",
    )
    parser.add_argument(
        "--glucose_bins",
        default=",".join(GLUCOSE_BINS_DEFAULT),
        help="Comma-separated glucose bins to include: low,mid,high.",
    )
    parser.add_argument(
        "--exclude_cell_lines",
        default=",".join(EXCLUDED_CELLLINES_DEFAULT),
        help="Comma-separated cellLine values to exclude.",
    )
    parser.add_argument(
        "--min_objects",
        type=int,
        default=1,
        help="Exclude images with fewer kept objects than this.",
    )
    parser.add_argument(
        "--allow_over_1000",
        action="store_true",
        help="Allow max_images above 1000. Off by default for the bounded WP3 pilot.",
    )
    return parser.parse_args()


def parse_csv_list(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def glucose_bin(g0: float) -> str | None:
    if abs(g0 - 0.0) < 1e-9 or abs(g0 - 0.1) < 1e-9:
        return "low"
    if abs(g0 - 0.25) < 1e-9 or abs(g0 - 0.5) < 1e-9 or abs(g0 - 1.0) < 1e-9:
        return "mid"
    if abs(g0 - 5.0) < 1e-9 or abs(g0 - 25.0) < 1e-9:
        return "high"
    return None


def assign_time_bins(rows: list[dict]) -> None:
    ordered = sorted(rows, key=lambda r: (float(r["hours"]), r["image_key"]))
    n = len(ordered)
    for idx, row in enumerate(ordered):
        row["time_bin"] = TIME_BINS[min(2, (idx * 3) // max(n, 1))]


def read_candidates(
    prelim_counts_csv: str,
    excluded_cell_lines: set[str],
    allowed_glucose_bins: set[str],
    min_objects: int,
) -> list[dict]:
    candidates: list[dict] = []
    with open(prelim_counts_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "cellLine",
            "experiment",
            "plateID",
            "pos",
            "ploidy",
            "G0",
            "hours",
            "image_key",
            "base_key",
            "alive_count",
            "dead_count",
            "total_count_kept",
            "total_area_px_all_objects",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")

        for row in reader:
            cell_line = row["cellLine"]
            if cell_line in excluded_cell_lines:
                continue
            try:
                g0 = float(row["G0"])
                total_count = int(float(row["total_count_kept"]))
                hours = float(row["hours"])
            except ValueError:
                continue
            if total_count < min_objects:
                continue
            gbin = glucose_bin(g0)
            if gbin is None or gbin not in allowed_glucose_bins:
                continue
            out = dict(row)
            out["G0"] = f"{g0:g}"
            out["hours"] = f"{hours:g}"
            out["glucose_bin"] = gbin
            out["total_count_kept"] = str(total_count)
            candidates.append(out)
    return candidates


def select_from_group(rows: list[dict], n: int) -> list[dict]:
    if n <= 0:
        return []
    if len(rows) <= n:
        return sorted(rows, key=lambda r: (float(r["hours"]), r["image_key"]))

    center = median(float(r["hours"]) for r in rows)
    ranked = sorted(
        rows,
        key=lambda r: (
            abs(float(r["hours"]) - center),
            -int(float(r["total_count_kept"])),
            r["image_key"],
        ),
    )
    return sorted(ranked[:n], key=lambda r: (float(r["hours"]), r["image_key"]))


def stride_cap(rows: list[dict], max_images: int) -> list[dict]:
    if len(rows) <= max_images:
        return rows
    ordered = sorted(rows, key=lambda r: (r["cellLine"], r["ploidy"], r["glucose_bin"], r["time_bin"], r["image_key"]))
    if max_images <= 0:
        return []
    step = len(ordered) / max_images
    picked = []
    used = set()
    for i in range(max_images):
        idx = min(len(ordered) - 1, int(i * step))
        while idx in used and idx + 1 < len(ordered):
            idx += 1
        if idx in used:
            break
        used.add(idx)
        picked.append(ordered[idx])
    return picked


def write_csv(path: str, rows: list[dict]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fieldnames = [
        "pilot_index",
        "cellLine",
        "experiment",
        "plateID",
        "pos",
        "ploidy",
        "G0",
        "glucose_bin",
        "time_bin",
        "hours",
        "image_key",
        "base_key",
        "alive_count",
        "dead_count",
        "total_count_kept",
        "total_area_px_all_objects",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for idx, row in enumerate(rows, start=1):
            out = {name: row.get(name, "") for name in fieldnames}
            out["pilot_index"] = idx
            writer.writerow(out)


def count_by(rows: Iterable[dict], keys: tuple[str, ...]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        label = " | ".join(row[k] for k in keys)
        counts[label] += 1
    return dict(sorted(counts.items()))


def main() -> None:
    args = parse_args()
    if args.images_per_stratum < 1:
        raise ValueError("--images_per_stratum must be >= 1.")
    if args.max_images < 1:
        raise ValueError("--max_images must be >= 1.")
    if args.max_images > MAX_IMAGES_HARD_LIMIT and not args.allow_over_1000:
        raise ValueError("--max_images must be <= 1000 for the bounded WP3 pilot.")

    excluded = set(parse_csv_list(args.exclude_cell_lines))
    allowed_bins = set(parse_csv_list(args.glucose_bins))
    unknown_bins = allowed_bins - set(GLUCOSE_BINS_DEFAULT)
    if unknown_bins:
        raise ValueError(f"Unknown glucose bins: {sorted(unknown_bins)}")

    candidates = read_candidates(args.prelim_counts_csv, excluded, allowed_bins, args.min_objects)
    by_core: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in candidates:
        by_core[(row["cellLine"], row["ploidy"], row["glucose_bin"])].append(row)
    for group_rows in by_core.values():
        assign_time_bins(group_rows)

    by_stratum: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in candidates:
        by_stratum[(row["cellLine"], row["ploidy"], row["glucose_bin"], row["time_bin"])].append(row)

    selected: list[dict] = []
    for stratum in sorted(by_stratum):
        selected.extend(select_from_group(by_stratum[stratum], args.images_per_stratum))
    selected = stride_cap(selected, args.max_images)
    selected = sorted(selected, key=lambda r: (r["cellLine"], r["ploidy"], r["glucose_bin"], r["time_bin"], float(r["hours"]), r["image_key"]))

    write_csv(args.output_csv, selected)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_meta_json)), exist_ok=True)
    meta = {
        "created_at": datetime.now().isoformat(),
        "input": os.path.abspath(args.prelim_counts_csv),
        "output_csv": os.path.abspath(args.output_csv),
        "rules": {
            "excluded_cell_lines": sorted(excluded),
            "glucose_bins": sorted(allowed_bins),
            "time_bins": list(TIME_BINS),
            "images_per_stratum": args.images_per_stratum,
            "max_images": args.max_images,
            "min_objects": args.min_objects,
            "selection": "nearest within-stratum median hour, then higher kept-object count, then image_key",
            "cap": "deterministic stride cap if selected images exceed max_images",
        },
        "summary": {
            "n_candidates": len(candidates),
            "n_selected": len(selected),
            "by_cellLine": count_by(selected, ("cellLine",)),
            "by_cellLine_ploidy": count_by(selected, ("cellLine", "ploidy")),
            "by_glucose_bin": count_by(selected, ("glucose_bin",)),
            "by_time_bin": count_by(selected, ("time_bin",)),
        },
    }
    with open(args.output_meta_json, "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote {len(selected)} pilot images to {args.output_csv}")
    print(f"Wrote metadata to {args.output_meta_json}")


if __name__ == "__main__":
    main()
