#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new_run_dir", required=True)
    parser.add_argument("--reference_run_dir", required=True)
    parser.add_argument("--output_json", required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def main() -> None:
    args = parse_args()
    new_dir = Path(args.new_run_dir)
    ref_dir = Path(args.reference_run_dir)
    new_images = read_csv(new_dir / "wp3_nuclear_image_qc.csv")
    keys = {row["image_key"] for row in new_images}
    if not keys:
        raise RuntimeError("New run has no image rows.")
    if any(row["status"] != "ok" for row in new_images):
        raise RuntimeError("New overlap run contains a non-ok image.")

    new_objects = [row for row in read_csv(new_dir / "wp3_nuclear_object_features.csv") if row["image_key"] in keys]
    ref_objects = [row for row in read_csv(ref_dir / "wp3_nuclear_object_features.csv") if row["image_key"] in keys]
    new_by_key = {(row["image_key"], int(row["object_id"])): row for row in new_objects}
    ref_by_key = {(row["image_key"], int(row["object_id"])): row for row in ref_objects}
    object_keys_match = set(new_by_key) == set(ref_by_key)
    common_columns = sorted(set(next(iter(new_by_key.values()))) & set(next(iter(ref_by_key.values()))))
    mismatched_rows = 0
    if object_keys_match:
        for key in new_by_key:
            if any(new_by_key[key][column] != ref_by_key[key][column] for column in common_columns):
                mismatched_rows += 1

    mask_checks = []
    mask_manifest = {row["image_key"]: row for row in read_csv(new_dir / "mask_manifest.csv")}
    for image_key in sorted(keys):
        stem = safe_name(image_key)
        for mask_type in ("cell", "nuclear"):
            new_path = Path(mask_manifest[image_key][f"{mask_type}_mask"])
            ref_path = ref_dir / "masks" / mask_type / f"{stem}_{mask_type}_masks.tif"
            match = new_path.exists() and ref_path.exists() and sha256(new_path) == sha256(ref_path)
            mask_checks.append({"image_key": image_key, "mask_type": mask_type, "exact_match": match})

    report = {
        "n_overlap_images": len(keys),
        "n_new_objects": len(new_objects),
        "n_reference_objects": len(ref_objects),
        "object_keys_match": object_keys_match,
        "n_object_rows_with_value_mismatch": mismatched_rows,
        "n_mask_files_checked": len(mask_checks),
        "n_exact_mask_matches": sum(row["exact_match"] for row in mask_checks),
        "passed": object_keys_match and mismatched_rows == 0 and all(row["exact_match"] for row in mask_checks),
        "mask_checks": mask_checks,
    }
    output = Path(args.output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, output)
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
