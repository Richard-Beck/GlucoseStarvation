#!/usr/bin/env python3
"""Prepare and score the canonical 90-frame segmentation/classifier validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import tifffile


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RELEASE = PROJECT_ROOT / "data/image_processing_validation/segmentation_classifier_90frame/20260722"
DEFAULT_ANNOTATION_MANIFEST = PROJECT_ROOT / "data/image_processing_runs/run_20260324_233122/annotation_stack_90_manifest.csv"
DEFAULT_ANNOTATION_RESULTS = PROJECT_ROOT / "data/image_processing_runs/run_20260324_233122/annotation_stack_90_Results.csv"
DEFAULT_MASTER_MANIFEST = PROJECT_ROOT / "data/image_processing_runs/full_segmentation_classification_nuclear/run_20260721_163410/manifests/master_manifest.csv"
DEFAULT_FROZEN_RUNNER = PROJECT_ROOT / "data/image_processing_runs/full_segmentation_classification_nuclear/run_20260721_163410/code/run_shard.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "summarize"))
    parser.add_argument("--release-dir", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--annotation-manifest", type=Path, default=DEFAULT_ANNOTATION_MANIFEST)
    parser.add_argument("--annotation-results", type=Path, default=DEFAULT_ANNOTATION_RESULTS)
    parser.add_argument("--master-manifest", type=Path, default=DEFAULT_MASTER_MANIFEST)
    parser.add_argument("--pipeline-dir", type=Path)
    parser.add_argument("--frozen-runner", type=Path, default=DEFAULT_FROZEN_RUNNER)
    parser.add_argument("--point-snap-radius", type=int, default=3)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError(f"Cannot infer columns for empty table: {path}")
        fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepare(args: argparse.Namespace) -> None:
    release = args.release_dir.resolve()
    if release.exists():
        raise FileExistsError(f"Release directory already exists: {release}")
    annotation = read_csv(args.annotation_manifest)
    master = read_csv(args.master_manifest)
    master_by_key = {row["image_key"]: row for row in master}
    if len(master_by_key) != len(master):
        raise ValueError("Full-run master manifest contains duplicate image keys.")
    missing = [row["image_key"] for row in annotation if row["image_key"] not in master_by_key]
    if missing:
        raise ValueError(f"Annotation image keys absent from full-run manifest: {missing[:5]}")
    if len(annotation) != 90 or len({row["image_key"] for row in annotation}) != 90:
        raise ValueError("Expected exactly 90 unique annotation frames.")

    release.mkdir(parents=True)
    selected = []
    for ann in sorted(annotation, key=lambda row: int(row["stack_index_1based"])):
        row = dict(master_by_key[ann["image_key"]])
        row["manifest_index"] = ann["stack_index_1based"]
        row["shard_id"] = "validation_90"
        selected.append(row)
    write_csv(release / "pipeline_manifest.csv", selected)
    settings = {
        "created_at": utc_now(),
        "purpose": "90-frame manual-annotation validation using the exact 20260721 full-stack image-processing runner",
        "release_dir": str(release),
        "pipeline_output_dir": str(release / "pipeline"),
        "annotation_manifest": str(args.annotation_manifest.resolve()),
        "annotation_manifest_sha256": sha256_file(args.annotation_manifest),
        "annotation_results": str(args.annotation_results.resolve()),
        "annotation_results_sha256": sha256_file(args.annotation_results),
        "full_run_master_manifest": str(args.master_manifest.resolve()),
        "full_run_master_manifest_sha256": sha256_file(args.master_manifest),
        "frozen_run_shard": str(args.frozen_runner.resolve()),
        "frozen_run_shard_sha256": sha256_file(args.frozen_runner),
        "builder": str(Path(__file__).resolve()),
        "builder_sha256": sha256_file(Path(__file__)),
        "point_snap_radius_px": args.point_snap_radius,
        "scope_rule": "retain context-crop objects; metrics retain every manually positive object plus target-centroid negatives",
    }
    with (release / "settings.json").open("w", encoding="utf-8") as handle:
        json.dump(settings, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Prepared {len(selected)} validation fields at {release}")


def mask_path(pipeline: Path, image_key: str) -> Path:
    stem = hashlib.sha256(image_key.encode("utf-8")).hexdigest()[:20]
    return pipeline / "masks" / "cell" / f"{stem}.tif"


def point_to_object(labels: np.ndarray, x: float, y: float, radius: int) -> tuple[int, str]:
    height, width = labels.shape
    xi = min(max(int(math.floor(x)), 0), width - 1)
    yi = min(max(int(math.floor(y)), 0), height - 1)
    direct = int(labels[yi, xi])
    if direct or radius <= 0:
        return direct, "direct" if direct else "unmapped"
    x0, x1 = max(0, xi - radius), min(width, xi + radius + 1)
    y0, y1 = max(0, yi - radius), min(height, yi + radius + 1)
    patch = labels[y0:y1, x0:x1]
    ys, xs = np.where(patch > 0)
    if not len(ys):
        return 0, "unmapped"
    distances = (xs + x0 - x) ** 2 + (ys + y0 - y) ** 2
    nearest = int(np.argmin(distances))
    return int(patch[ys[nearest], xs[nearest]]), "snapped"


def auc_rank(y: np.ndarray, score: np.ndarray) -> float:
    positive = y == 1
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if not n_pos or not n_neg:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    ranks = np.empty(len(score), dtype=float)
    start = 0
    while start < len(score):
        end = start + 1
        while end < len(score) and sorted_score[end] == sorted_score[start]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2
        start = end
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def metric_value(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else float("nan")


def summarize(args: argparse.Namespace) -> None:
    release = args.release_dir.resolve()
    pipeline = (args.pipeline_dir or release / "pipeline").resolve()
    if not (pipeline / "COMPLETE").is_file():
        raise FileNotFoundError(f"Pipeline run is not complete: {pipeline}")
    images = read_csv(pipeline / "images.csv")
    if len(images) != 90 or any(row["processing_status"] != "ok" for row in images):
        raise ValueError("Expected 90 successfully processed validation fields.")
    annotation = read_csv(args.annotation_manifest)
    annotation_by_frame = {int(row["stack_index_1based"]): row for row in annotation}
    annotation_by_key = {row["image_key"]: row for row in annotation}
    objects = read_csv(pipeline / "objects.csv")

    points_by_frame: dict[int, list[dict]] = {frame: [] for frame in annotation_by_frame}
    for point_index, row in enumerate(read_csv(args.annotation_results), 1):
        if not row.get("Frame", "").strip():
            continue
        frame = int(float(row["Frame"]))
        ann = annotation_by_frame[frame]
        x_crop, y_crop = float(row["X"]), float(row["Y"])
        points_by_frame[frame].append({
            "manual_point_id": point_index, "frame": frame, "image_key": ann["image_key"],
            "x_crop": f"{x_crop:.6f}", "y_crop": f"{y_crop:.6f}",
            "x_full": f"{x_crop + float(ann['crop_x0']):.6f}",
            "y_full": f"{y_crop + float(ann['crop_y0']):.6f}",
        })

    manual_rows: list[dict] = []
    positive_ids: dict[str, set[int]] = {}
    for frame, ann in annotation_by_frame.items():
        labels = tifffile.imread(mask_path(pipeline, ann["image_key"]))
        positive_ids[ann["image_key"]] = set()
        for point in points_by_frame[frame]:
            object_id, mapping = point_to_object(labels, float(point["x_full"]), float(point["y_full"]), args.point_snap_radius)
            point["mapped_object_id"] = object_id
            point["mapping"] = mapping
            if object_id:
                positive_ids[ann["image_key"]].add(object_id)
            manual_rows.append(point)

    validation_objects: list[dict] = []
    for row in objects:
        ann = annotation_by_key[row["image_key"]]
        x, y = float(row["centroid_x"]), float(row["centroid_y"])
        in_context = float(ann["crop_x0"]) <= x < float(ann["crop_x1"]) and float(ann["crop_y0"]) <= y < float(ann["crop_y1"])
        if not in_context:
            continue
        in_target = float(ann["target_x0"]) <= x <= float(ann["target_x1"]) and float(ann["target_y0"]) <= y <= float(ann["target_y1"])
        object_id = int(row["object_id"])
        truth = int(object_id in positive_ids[row["image_key"]])
        keep_metric = bool(truth or in_target)
        validation_objects.append({
            "frame": int(ann["stack_index_1based"]), "image_key": row["image_key"],
            "object_id": object_id, "centroid_x_full": row["centroid_x"], "centroid_y_full": row["centroid_y"],
            "centroid_x_crop": f"{x - float(ann['crop_x0']):.3f}", "centroid_y_crop": f"{y - float(ann['crop_y0']):.3f}",
            "bbox_x0_full": row["bbox_x0"], "bbox_y0_full": row["bbox_y0"], "bbox_x1_full": row["bbox_x1"], "bbox_y1_full": row["bbox_y1"],
            "predicted_label_id": row["predicted_label_id"], "predicted_label_name": row["predicted_label_name"],
            "prob_alive": row["prob_alive"], "prob_dead": row["prob_dead"], "prob_junk": row["prob_junk"],
            "manual_alive": truth, "in_annotation_target": int(in_target), "keep_for_target_metrics": int(keep_metric),
            "segmented_area_px": row["segmented_area_px"], "cell_area_pass": row["cell_area_pass"],
            "nuclear_area_px": row["nuclear_area_px"],
        })

    metric_rows = [row for row in validation_objects if row["keep_for_target_metrics"]]
    y = np.asarray([row["manual_alive"] for row in metric_rows], dtype=int)
    predicted = np.asarray([int(row["predicted_label_id"]) == 1 for row in metric_rows], dtype=int)
    score = np.asarray([float(row["prob_alive"]) for row in metric_rows], dtype=float)
    tp = int(((y == 1) & (predicted == 1)).sum())
    fp = int(((y == 0) & (predicted == 1)).sum())
    tn = int(((y == 0) & (predicted == 0)).sum())
    fn = int(((y == 1) & (predicted == 0)).sum())
    precision, recall = metric_value(tp, tp + fp), metric_value(tp, tp + fn)
    specificity = metric_value(tn, tn + fp)
    metrics = [{
        "scope": "manual positives plus target-centroid negatives", "n_objects": len(metric_rows),
        "n_positive": int(y.sum()), "n_negative": int((y == 0).sum()),
        "roc_auc": f"{auc_rank(y, score):.8f}", "balanced_accuracy": f"{(recall + specificity) / 2:.8f}",
        "precision": f"{precision:.8f}", "recall": f"{recall:.8f}", "specificity": f"{specificity:.8f}",
        "f1": f"{metric_value(2 * tp, 2 * tp + fp + fn):.8f}", "tp": tp, "fp": fp, "tn": tn, "fn": fn,
    }]

    by_key: dict[str, list[dict]] = {}
    for row in validation_objects:
        by_key.setdefault(row["image_key"], []).append(row)
    frame_rows = []
    for frame, ann in sorted(annotation_by_frame.items()):
        rows = by_key.get(ann["image_key"], [])
        target = [row for row in rows if row["in_annotation_target"]]
        points = points_by_frame[frame]
        mapped = [row for row in points if row.get("mapped_object_id", 0)]
        frame_rows.append({
            "frame": frame, "image_key": ann["image_key"], "cellLine": ann["cellLine"], "ploidy": ann["ploidy"],
            "G0": ann["G0"], "glucose_bin": ann["glucose_bin"], "time_bin": ann["time_bin"], "hours": ann["hours"],
            "manual_point_count": len(points), "mapped_manual_point_count": len(mapped),
            "manual_alive_object_count": len(positive_ids[ann["image_key"]]),
            "target_segmented_object_count": len(target),
            "target_predicted_alive_count": sum(int(row["predicted_label_id"]) == 1 for row in target),
            "target_predicted_alive_probability_sum": f"{sum(float(row['prob_alive']) for row in target):.8f}",
        })

    write_csv(release / "manual_points.csv", manual_rows)
    write_csv(release / "object_predictions.csv", validation_objects)
    write_csv(release / "frame_summary.csv", frame_rows)
    write_csv(release / "object_metrics.csv", metrics)
    summary = {
        "completed_at": utc_now(), "n_frames": 90, "n_manual_points": len(manual_rows),
        "n_mapped_manual_points": sum(int(row["mapped_object_id"] != 0) for row in manual_rows),
        "n_unique_manual_alive_objects": sum(len(ids) for ids in positive_ids.values()),
        "n_context_objects": len(validation_objects), "metrics": metrics[0],
    }
    with (release / "validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    manifest_rows = []
    for path in sorted(item for item in release.rglob("*") if item.is_file() and item.name != "output_manifest.csv"):
        manifest_rows.append({"path": path.relative_to(release).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_csv(release / "output_manifest.csv", manifest_rows)
    (release / "COMPLETE").write_text("complete\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


def main() -> None:
    args = parse_args()
    if args.point_snap_radius < 0:
        raise ValueError("--point-snap-radius must be nonnegative.")
    if args.stage == "prepare":
        prepare(args)
    else:
        summarize(args)


if __name__ == "__main__":
    main()
