#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp3_nuclear_size_pilot import build_raw_group_map, load_channels, safe_name, write_qc_panel  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest_csv", required=True)
    parser.add_argument("--mask_manifest_csv", required=True)
    parser.add_argument("--raw_data_dir", default="all_raw")
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def read_csv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    args = parse_args()
    manifest = read_csv(args.manifest_csv)
    mask_lookup = {row["image_key"]: row for row in read_csv(args.mask_manifest_csv)}
    selected = []
    for batch_id in ("C00", "I00", "M00b"):
        for ploidy in ("2N", "4N"):
            rows = [row for row in manifest if row["batch_id"] == batch_id and row["ploidy"] == ploidy]
            g0_values = sorted({float(row["G0"]) for row in rows})
            g0_selected = [g0_values[0], g0_values[len(g0_values) // 2], g0_values[-1]]
            for g0 in g0_selected:
                condition = [row for row in rows if float(row["G0"]) == g0]
                hours = sorted({float(row["hours"]) for row in condition})
                for hour in (hours[0], hours[-1]):
                    choices = [row for row in condition if float(row["hours"]) == hour]
                    choices.sort(key=lambda row: (row["plateID"], float(row["pos"]), row["image_key"]))
                    selected.append(choices[0])
    keys = {row["image_key"] for row in selected}
    raw_groups = build_raw_group_map(args.raw_data_dir, keys)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_rows = []
    for index, row in enumerate(selected, start=1):
        key = row["image_key"]
        masks = mask_lookup[key]
        cell_masks = tifffile.imread(masks["cell_mask"])
        nuclear_masks = tifffile.imread(masks["nuclear_mask"])
        channels = load_channels(raw_groups[key])
        output_path = output_dir / f"{index:03d}_{safe_name(key)}.png"
        write_qc_panel(str(output_path), row, channels, cell_masks, nuclear_masks)
        qc_rows.append({
            "qc_index": index,
            "output_png": str(output_path),
            **row,
            "cell_mask": masks["cell_mask"],
            "nuclear_mask": masks["nuclear_mask"],
        })
    fields = list(qc_rows[0])
    with open(output_dir / "qc_manifest.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(qc_rows)
    print(f"Wrote {len(qc_rows)} stratified QC panels to {output_dir}")


if __name__ == "__main__":
    main()
