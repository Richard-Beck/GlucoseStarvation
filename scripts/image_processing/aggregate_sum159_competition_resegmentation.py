#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path


EXPECTED_IMAGES = 1200
EXPECTED_SERIES = 80


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate and validate SUM-159 competition segmentation shards.")
    parser.add_argument("--manifest_tsv", required=True)
    parser.add_argument("--run_dir", required=True)
    return parser.parse_args()


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concatenate_csvs(paths: list[Path], output_path: Path) -> tuple[int, set[str]]:
    n_rows = 0
    image_keys: set[str] = set()
    expected_header: list[str] | None = None
    with output_path.open("w", newline="", encoding="utf-8") as output_handle:
        writer = None
        for path in paths:
            with path.open(newline="", encoding="utf-8") as input_handle:
                reader = csv.DictReader(input_handle)
                header = reader.fieldnames or []
                if expected_header is None:
                    expected_header = header
                    writer = csv.DictWriter(output_handle, fieldnames=expected_header)
                    writer.writeheader()
                elif header != expected_header:
                    raise RuntimeError(f"Header mismatch in {path}")
                assert writer is not None
                for row in reader:
                    writer.writerow(row)
                    n_rows += 1
                    image_keys.add(row["image_key"])
    return n_rows, image_keys


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    manifest_rows = read_tsv(args.manifest_tsv)
    if len(manifest_rows) != EXPECTED_IMAGES:
        raise RuntimeError(f"Manifest has {len(manifest_rows)} rows; expected {EXPECTED_IMAGES}")
    expected_keys = {row["image_key"] for row in manifest_rows}
    if len(expected_keys) != EXPECTED_IMAGES:
        raise RuntimeError("Manifest image_key values are not unique")

    object_parts: list[Path] = []
    image_parts: list[Path] = []
    total_shard_objects = 0
    for series_index in range(EXPECTED_SERIES):
        shard_dir = run_dir / "shards" / f"series_{series_index:03d}"
        status_path = shard_dir / "status.json"
        if not status_path.is_file():
            raise RuntimeError(f"Missing shard status: {status_path}")
        with status_path.open(encoding="utf-8") as handle:
            status = json.load(handle)
        if int(status["n_errors"]) != 0 or int(status["n_images"]) != 15:
            raise RuntimeError(f"Invalid shard status for series {series_index}: {status}")
        total_shard_objects += int(status["n_objects"])
        object_parts.append(shard_dir / "object_features_and_predictions.csv")
        image_parts.append(shard_dir / "image_summary.csv")

    n_objects, object_image_keys = concatenate_csvs(
        object_parts, run_dir / "object_features_and_predictions.csv"
    )
    n_image_rows, image_keys = concatenate_csvs(image_parts, run_dir / "image_summary.csv")
    if n_objects != total_shard_objects:
        raise RuntimeError(
            f"Object-row mismatch: concatenated={n_objects}, shard statuses={total_shard_objects}"
        )
    if n_image_rows != EXPECTED_IMAGES or image_keys != expected_keys:
        raise RuntimeError(
            f"Image-summary coverage mismatch: rows={n_image_rows}, "
            f"missing={len(expected_keys - image_keys)}, extra={len(image_keys - expected_keys)}"
        )
    with (run_dir / "image_summary.csv").open(newline="", encoding="utf-8") as handle:
        image_summary_rows = list(csv.DictReader(handle))
    no_cell_keys = {
        row["image_key"] for row in image_summary_rows if row["status"] == "no_cells"
    }
    expected_object_keys = expected_keys - no_cell_keys
    if object_image_keys != expected_object_keys:
        raise RuntimeError(
            f"Object-table image coverage mismatch after allowing no-cell images: "
            f"missing={len(expected_object_keys - object_image_keys)}, "
            f"extra={len(object_image_keys - expected_object_keys)}"
        )

    missing_cell_masks = []
    missing_nuclear_masks = []
    for key in sorted(expected_keys):
        safe_key = key
        cell_path = run_dir / "masks" / "cell" / f"{safe_key}_cell_masks.tif"
        nuclear_path = run_dir / "masks" / "nuclear" / f"{safe_key}_nuclear_masks.tif"
        if not cell_path.is_file():
            missing_cell_masks.append(str(cell_path))
        if not nuclear_path.is_file():
            missing_nuclear_masks.append(str(nuclear_path))
    if missing_cell_masks or missing_nuclear_masks:
        raise RuntimeError(
            f"Missing masks: cell={len(missing_cell_masks)}, nuclear={len(missing_nuclear_masks)}"
        )

    summary = {
        "status": "pass",
        "validated_at": datetime.now().isoformat(),
        "n_series": EXPECTED_SERIES,
        "n_images": EXPECTED_IMAGES,
        "n_objects": n_objects,
        "n_cell_masks": EXPECTED_IMAGES,
        "n_nuclear_masks": EXPECTED_IMAGES,
        "green_metrics": [
            "green_mean_raw",
            "green_median_raw",
            "green_p90_raw",
            "green_integrated_raw",
            "green_mean_robust01",
            "green_bg_z",
        ],
        "green_excluded_from_segmentation": True,
        "green_excluded_from_classifier": True,
        "qc_overlays_written": False,
    }
    with (run_dir / "validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
