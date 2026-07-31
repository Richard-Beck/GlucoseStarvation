#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path


CHANNEL_SPECS = {
    "red": ("Red", re.compile(r"^Red_(?P<well>[FG]\d+)_(?P<position>[1-4])_(?P<timestamp>\d{2}d\d{2}h\d{2}m)\.tif$")),
    "green": ("Green", re.compile(r"^Green_(?P<well>[FG]\d+)_(?P<position>[1-4])_(?P<timestamp>\d{2}d\d{2}h\d{2}m)\.tif$")),
    "nir": ("NIR", re.compile(r"^NIR__(?P<well>[FG]\d+)_(?P<position>[1-4])_(?P<timestamp>\d{2}d\d{2}h\d{2}m)\.tif$")),
    "phase": ("phase", re.compile(r"^Phase_(?P<well>[FG]\d+)_(?P<position>[1-4])_(?P<timestamp>\d{2}d\d{2}h\d{2}m)\.tif$")),
}
GLUCOSE_BY_COLUMN = {
    2: 0.0,
    3: 0.0,
    4: 0.1,
    5: 0.1,
    6: 0.25,
    7: 0.25,
    8: 0.5,
    9: 0.5,
    10: 1.0,
    11: 1.0,
}
EXPECTED_WELLS = {f"{row}{column}" for row in ("F", "G") for column in range(2, 12)}
EXPECTED_POSITIONS = {1, 2, 3, 4}
EXPECTED_IMAGES_PER_SERIES = 15
EXPECTED_SERIES = 80
EXPECTED_IMAGES = EXPECTED_SERIES * EXPECTED_IMAGES_PER_SERIES
TIMESTAMP_RE = re.compile(r"^(?P<days>\d{2})d(?P<hours>\d{2})h(?P<minutes>\d{2})m$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the complete SUM-159 J00 competition image manifest.")
    parser.add_argument("--raw_data_dir", required=True, help="J00_IncucyteRawDataLiveDead_Competition directory.")
    parser.add_argument("--output_tsv", required=True, help="Destination manifest TSV.")
    parser.add_argument("--output_meta_json", required=True, help="Destination manifest metadata JSON.")
    return parser.parse_args()


def well_sort_key(well: str) -> tuple[str, int]:
    return well[0], int(well[1:])


def timestamp_hours(timestamp: str) -> float:
    match = TIMESTAMP_RE.match(timestamp)
    if match is None:
        raise ValueError(f"Invalid timestamp: {timestamp}")
    return (
        int(match.group("days")) * 24
        + int(match.group("hours"))
        + int(match.group("minutes")) / 60
    )


def scan_channel(raw_data_dir: Path, channel: str) -> dict[tuple[str, int, str], str]:
    dirname, pattern = CHANNEL_SPECS[channel]
    channel_dir = raw_data_dir / dirname
    if not channel_dir.is_dir():
        raise FileNotFoundError(f"Missing channel directory: {channel_dir}")

    found: dict[tuple[str, int, str], str] = {}
    for path in sorted(channel_dir.iterdir()):
        if not path.is_file():
            continue
        match = pattern.match(path.name)
        if match is None:
            continue
        key = (
            match.group("well"),
            int(match.group("position")),
            match.group("timestamp"),
        )
        if key in found:
            raise ValueError(f"Duplicate {channel} acquisition key: {key}")
        found[key] = str(path.resolve())
    return found


def main() -> None:
    args = parse_args()
    raw_data_dir = Path(args.raw_data_dir).resolve()
    channel_maps = {channel: scan_channel(raw_data_dir, channel) for channel in CHANNEL_SPECS}

    key_sets = {channel: set(paths) for channel, paths in channel_maps.items()}
    union_keys = set().union(*key_sets.values())
    incomplete = {
        key: sorted(channel for channel, keys in key_sets.items() if key not in keys)
        for key in union_keys
        if any(key not in keys for keys in key_sets.values())
    }
    if incomplete:
        first_key = sorted(incomplete)[0]
        raise RuntimeError(
            f"Found {len(incomplete)} incomplete image groups; first is {first_key}, "
            f"missing {incomplete[first_key]}"
        )

    series_to_keys: dict[tuple[str, int], list[tuple[str, int, str]]] = defaultdict(list)
    for key in union_keys:
        series_to_keys[(key[0], key[1])].append(key)

    observed_wells = {well for well, _ in series_to_keys}
    observed_positions = {position for _, position in series_to_keys}
    if observed_wells != EXPECTED_WELLS:
        raise RuntimeError(
            f"Well coverage mismatch: missing={sorted(EXPECTED_WELLS - observed_wells)}, "
            f"extra={sorted(observed_wells - EXPECTED_WELLS)}"
        )
    if observed_positions != EXPECTED_POSITIONS:
        raise RuntimeError(
            f"Position coverage mismatch: expected={sorted(EXPECTED_POSITIONS)}, "
            f"observed={sorted(observed_positions)}"
        )
    if len(series_to_keys) != EXPECTED_SERIES:
        raise RuntimeError(f"Expected {EXPECTED_SERIES} well-position series, found {len(series_to_keys)}")

    ordered_series = sorted(series_to_keys, key=lambda x: (well_sort_key(x[0]), x[1]))
    series_index = {series: idx for idx, series in enumerate(ordered_series)}
    rows: list[dict[str, object]] = []
    for series in ordered_series:
        keys = sorted(series_to_keys[series], key=lambda x: timestamp_hours(x[2]))
        if len(keys) != EXPECTED_IMAGES_PER_SERIES:
            raise RuntimeError(
                f"Expected {EXPECTED_IMAGES_PER_SERIES} images for {series}, found {len(keys)}"
            )
        well, position = series
        column = int(well[1:])
        for key in keys:
            timestamp = key[2]
            rows.append(
                {
                    "image_key": f"SUM159_J00_Competition_{well}_{position}_{timestamp}",
                    "series_id": f"{well}_{position}",
                    "series_index": series_index[series],
                    "well": well,
                    "position": position,
                    "timestamp": timestamp,
                    "hours": f"{timestamp_hours(timestamp):g}",
                    "G0_mM": f"{GLUCOSE_BY_COLUMN[column]:g}",
                    "red_path": channel_maps["red"][key],
                    "green_path": channel_maps["green"][key],
                    "nir_path": channel_maps["nir"][key],
                    "phase_path": channel_maps["phase"][key],
                }
            )

    if len(rows) != EXPECTED_IMAGES:
        raise RuntimeError(f"Expected {EXPECTED_IMAGES} complete image groups, found {len(rows)}")

    output_tsv = Path(args.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "created_at": datetime.now().isoformat(),
        "raw_data_dir": str(raw_data_dir),
        "n_images": len(rows),
        "n_series": len(ordered_series),
        "images_per_series": EXPECTED_IMAGES_PER_SERIES,
        "channel_counts": {channel: len(paths) for channel, paths in channel_maps.items()},
        "green_usage": "measured_only; excluded from segmentation and classifier",
        "segmentation_fluorescence": "red + NIR",
        "classifier_composite": ["NIR", "red", "phase"],
        "glucose_by_column_mM": {str(key): value for key, value in GLUCOSE_BY_COLUMN.items()},
    }
    output_meta = Path(args.output_meta_json)
    output_meta.parent.mkdir(parents=True, exist_ok=True)
    with output_meta.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote {len(rows)} complete image groups across {len(ordered_series)} series to {output_tsv}")


if __name__ == "__main__":
    main()
