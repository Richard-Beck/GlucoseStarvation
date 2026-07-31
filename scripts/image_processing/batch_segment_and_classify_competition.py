#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import tifffile
import torch
from cellpose.models import CellposeModel

sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier

from wp3_nuclear_size_pilot import measure_objects, robust01, segment_cell_masks


LABEL_MAP = {1: "alive", 2: "dead", 3: "junk"}
MIN_CELL_AREA_PX = 50
MIN_NUCLEUS_AREA_PX = 10
EPS = 1e-8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Segment and classify one SUM-159 competition well-position series, "
            "measure cytoplasmic GFP, and reproduce the manuscript nuclear-area method."
        )
    )
    parser.add_argument("--manifest_tsv", required=True)
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--series_index", required=True, type=int)
    parser.add_argument("--classifier_path", required=True)
    parser.add_argument("--cellpose_batch_size", type=int, default=8)
    return parser.parse_args()


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_manifest_series(path: str, series_index: int) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if int(row["series_index"]) == series_index
        ]
    rows.sort(key=lambda row: float(row["hours"]))
    if len(rows) != 15:
        raise RuntimeError(f"Series {series_index} has {len(rows)} rows; expected 15")
    return rows


def load_float32(path: str) -> np.ndarray:
    image = tifffile.imread(path).astype(np.float32)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D image at {path}, found shape {image.shape}")
    return image


def background_stats(image: np.ndarray, background_mask: np.ndarray) -> tuple[float, float]:
    values = image[background_mask]
    if values.size == 0:
        values = image.ravel()
    mean = float(np.mean(values))
    sd = max(float(np.std(values)), EPS)
    return mean, sd


def green_metrics(
    green_raw: np.ndarray,
    green01: np.ndarray,
    cell_masks: np.ndarray,
) -> dict[int, dict[str, object]]:
    background_mask = cell_masks == 0
    bg_mean, bg_sd = background_stats(green_raw, background_mask)
    metrics: dict[int, dict[str, object]] = {}
    for object_id in range(1, int(cell_masks.max()) + 1):
        object_mask = cell_masks == object_id
        values = green_raw[object_mask]
        if values.size == 0:
            continue
        mean_raw = float(np.mean(values))
        metrics[object_id] = {
            "green_mean_raw": f"{mean_raw:.8f}",
            "green_median_raw": f"{float(np.median(values)):.8f}",
            "green_p90_raw": f"{float(np.percentile(values, 90)):.8f}",
            "green_integrated_raw": f"{float(np.sum(values, dtype=np.float64)):.8f}",
            "green_mean_robust01": f"{float(np.mean(green01[object_mask])):.8f}",
            "green_bg_z": f"{((mean_raw - bg_mean) / bg_sd):.8f}",
        }
    return metrics


def summarize_image(
    row: dict[str, str],
    object_rows: list[dict[str, object]],
    status: str,
    elapsed_seconds: float,
    error: str = "",
) -> dict[str, object]:
    class_counts = {name: 0 for name in LABEL_MAP.values()}
    for object_row in object_rows:
        label = str(object_row.get("predicted_label_name", ""))
        if label in class_counts:
            class_counts[label] += 1
    n_nuclei = sum(int(object_row["nuclear_area_px"]) > 0 for object_row in object_rows)
    green_values = [
        float(object_row["green_mean_raw"])
        for object_row in object_rows
        if object_row.get("green_mean_raw", "") != ""
    ]
    return {
        "image_key": row["image_key"],
        "series_id": row["series_id"],
        "series_index": row["series_index"],
        "well": row["well"],
        "position": row["position"],
        "timestamp": row["timestamp"],
        "hours": row["hours"],
        "G0_mM": row["G0_mM"],
        "status": status,
        "error": error,
        "n_objects": len(object_rows),
        "n_alive": class_counts["alive"],
        "n_dead": class_counts["dead"],
        "n_junk": class_counts["junk"],
        "n_with_nucleus": n_nuclei,
        "nuclear_success_fraction": f"{n_nuclei / max(len(object_rows), 1):.8f}",
        "median_green_mean_raw": f"{float(np.median(green_values)):.8f}" if green_values else "",
        "elapsed_seconds": f"{elapsed_seconds:.3f}",
    }


def main() -> None:
    args = parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; submit through the competition GPU SLURM wrapper")

    manifest_rows = read_manifest_series(args.manifest_tsv, args.series_index)
    run_dir = Path(args.run_dir).resolve()
    shard_dir = run_dir / "shards" / f"series_{args.series_index:03d}"
    cell_mask_dir = run_dir / "masks" / "cell"
    nuclear_mask_dir = run_dir / "masks" / "nuclear"
    shard_dir.mkdir(parents=True, exist_ok=False)
    cell_mask_dir.mkdir(parents=True, exist_ok=True)
    nuclear_mask_dir.mkdir(parents=True, exist_ok=True)

    classifier = ObjectClassifier()
    classifier.load_state(args.classifier_path)
    if not classifier.is_trained:
        raise RuntimeError("The production object classifier did not load as trained")
    model = CellposeModel(gpu=True)

    object_headers = [
        "image_key", "series_id", "series_index", "well", "position", "timestamp", "hours", "G0_mM",
        "cell_mask_file", "nuclear_mask_file", "object_id",
        "cell_area_px", "cell_area_pass", "centroid_x", "centroid_y",
        "bbox_x0", "bbox_y0", "bbox_x1", "bbox_y1",
        "nuclear_area_px", "largest_nuclear_area_px", "nuclear_component_count",
        "nuclear_to_cell_area_ratio", "nuclear_mean_intensity", "nuclear_integrated_intensity",
        "cell_mean_nuclear_channel", "cell_mean_alive_channel", "cell_mean_dead_channel",
        "green_mean_raw", "green_median_raw", "green_p90_raw", "green_integrated_raw",
        "green_mean_robust01", "green_bg_z",
        "predicted_label_id", "predicted_label_name", "prob_alive", "prob_dead", "prob_junk",
    ]
    image_headers = [
        "image_key", "series_id", "series_index", "well", "position", "timestamp", "hours", "G0_mM",
        "status", "error", "n_objects", "n_alive", "n_dead", "n_junk", "n_with_nucleus",
        "nuclear_success_fraction", "median_green_mean_raw", "elapsed_seconds",
    ]
    object_csv = shard_dir / "object_features_and_predictions.csv"
    image_csv = shard_dir / "image_summary.csv"
    n_objects_total = 0
    n_errors = 0

    with object_csv.open("w", newline="", encoding="utf-8") as object_handle, image_csv.open(
        "w", newline="", encoding="utf-8"
    ) as image_handle:
        object_writer = csv.DictWriter(object_handle, fieldnames=object_headers)
        image_writer = csv.DictWriter(image_handle, fieldnames=image_headers)
        object_writer.writeheader()
        image_writer.writeheader()

        for image_number, row in enumerate(manifest_rows, start=1):
            start = time.monotonic()
            key = row["image_key"]
            print(f"[series {args.series_index:03d} {image_number:02d}/15] {key}", flush=True)
            object_rows: list[dict[str, object]] = []
            try:
                red_raw = load_float32(row["red_path"])
                green_raw = load_float32(row["green_path"])
                nir_raw = load_float32(row["nir_path"])
                phase_raw = load_float32(row["phase_path"])
                shapes = {red_raw.shape, green_raw.shape, nir_raw.shape, phase_raw.shape}
                if len(shapes) != 1:
                    raise ValueError(f"Channel shape mismatch: {sorted(shapes)}")

                channels = {
                    "alive_raw": red_raw,
                    "dead_raw": nir_raw,
                    "phase_raw": phase_raw,
                    "alive01": robust01(red_raw),
                    "dead01": robust01(nir_raw),
                    "phase01": robust01(phase_raw),
                    "nuclear01": robust01(red_raw + nir_raw),
                }
                cell_masks = segment_cell_masks(channels, model, args.cellpose_batch_size)
                if cell_masks.size == 0 or int(cell_masks.max()) == 0:
                    safe_key = safe_name(key)
                    empty_masks = np.zeros(red_raw.shape, dtype=np.uint16)
                    tifffile.imwrite(
                        cell_mask_dir / f"{safe_key}_cell_masks.tif",
                        empty_masks,
                    )
                    tifffile.imwrite(
                        nuclear_mask_dir / f"{safe_key}_nuclear_masks.tif",
                        empty_masks,
                    )
                    image_writer.writerow(
                        summarize_image(row, [], "no_cells", time.monotonic() - start)
                    )
                    image_handle.flush()
                    continue

                measured_rows, nuclear_masks = measure_objects(
                    cell_masks,
                    channels["nuclear01"],
                    channels["alive01"],
                    channels["dead01"],
                    MIN_CELL_AREA_PX,
                    MIN_NUCLEUS_AREA_PX,
                )
                green_by_object = green_metrics(green_raw, robust01(green_raw), cell_masks)

                classifier_composite = np.stack(
                    [channels["dead01"], channels["alive01"], channels["phase01"]],
                    axis=-1,
                ).astype(np.float32, copy=False)
                predictions = classifier.predict_with_probabilities(classifier_composite, cell_masks)

                safe_key = safe_name(key)
                cell_mask_file = f"masks/cell/{safe_key}_cell_masks.tif"
                nuclear_mask_file = f"masks/nuclear/{safe_key}_nuclear_masks.tif"
                mask_dtype = np.uint16 if int(cell_masks.max()) <= np.iinfo(np.uint16).max else np.uint32
                tifffile.imwrite(run_dir / cell_mask_file, cell_masks.astype(mask_dtype, copy=False))
                tifffile.imwrite(run_dir / nuclear_mask_file, nuclear_masks.astype(mask_dtype, copy=False))

                for measured in measured_rows:
                    object_id = int(measured["object_id"])
                    prediction = predictions.get(object_id)
                    if prediction is None:
                        raise RuntimeError(f"Classifier returned no result for object {object_id}")
                    output_row: dict[str, object] = {
                        name: row.get(name, "")
                        for name in (
                            "image_key", "series_id", "series_index", "well",
                            "position", "timestamp", "hours", "G0_mM",
                        )
                    }
                    output_row.update(
                        {
                            "cell_mask_file": cell_mask_file,
                            "nuclear_mask_file": nuclear_mask_file,
                        }
                    )
                    output_row.update(measured)
                    output_row.update(green_by_object.get(object_id, {}))
                    predicted_id = int(prediction["prediction"])
                    output_row.update(
                        {
                            "predicted_label_id": predicted_id,
                            "predicted_label_name": LABEL_MAP.get(predicted_id, "unknown"),
                            "prob_alive": f"{float(prediction['probabilities'][0]):.8f}",
                            "prob_dead": f"{float(prediction['probabilities'][1]):.8f}",
                            "prob_junk": f"{float(prediction['probabilities'][2]):.8f}",
                        }
                    )
                    object_writer.writerow(output_row)
                    object_rows.append(output_row)

                elapsed = time.monotonic() - start
                image_writer.writerow(summarize_image(row, object_rows, "ok", elapsed))
                n_objects_total += len(object_rows)
                object_handle.flush()
                image_handle.flush()
            except Exception as exc:
                n_errors += 1
                elapsed = time.monotonic() - start
                print(f"ERROR {key}: {exc}", flush=True)
                image_writer.writerow(summarize_image(row, [], "error", elapsed, str(exc)))
                image_handle.flush()

    status = {
        "completed_at": datetime.now().isoformat(),
        "series_index": args.series_index,
        "series_id": manifest_rows[0]["series_id"],
        "n_images": len(manifest_rows),
        "n_objects": n_objects_total,
        "n_errors": n_errors,
        "classifier_path": str(Path(args.classifier_path).resolve()),
        "classifier_sha256": sha256_file(args.classifier_path),
        "cellpose_backend": "CellposeModel / Cellpose v4 CPSAM",
        "segmentation_channels": ["phase", "red + NIR"],
        "classifier_channels": ["NIR", "red", "phase"],
        "green_usage": "per-object measurement only",
        "nuclear_method_source": "scripts/image_processing/wp3_nuclear_size_pilot.py",
        "min_cell_area_px": MIN_CELL_AREA_PX,
        "min_nucleus_area_px": MIN_NUCLEUS_AREA_PX,
    }
    with (shard_dir / "status.json").open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if n_errors:
        raise RuntimeError(f"Series {args.series_index} completed with {n_errors} image errors")


if __name__ == "__main__":
    main()
