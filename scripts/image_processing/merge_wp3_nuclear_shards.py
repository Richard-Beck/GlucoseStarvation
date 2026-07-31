#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master_manifest", required=True)
    parser.add_argument("--shard_root", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--shard_prefix", default="shard")
    parser.add_argument("--n_shards", type=int, default=4)
    return parser.parse_args()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def atomic_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def main() -> None:
    args = parse_args()
    if args.n_shards != 4:
        raise ValueError("The approved run requires exactly four shards.")
    _, master = read_csv(Path(args.master_manifest))
    master_keys = [row["image_key"] for row in master]
    if len(master_keys) != len(set(master_keys)):
        raise RuntimeError("Master manifest contains duplicate image keys.")
    order = {key: index for index, key in enumerate(master_keys)}

    shard_manifest_keys: list[str] = []
    image_rows: list[dict[str, str]] = []
    object_rows: list[dict[str, str]] = []
    mask_rows: list[dict[str, str]] = []
    object_fields: list[str] | None = None
    image_fields: list[str] | None = None
    shard_summaries = []

    for shard_id in range(args.n_shards):
        run_dir = Path(args.shard_root) / f"{args.shard_prefix}_{shard_id:03d}"
        required = (
            run_dir / "pilot_manifest.csv",
            run_dir / "wp3_nuclear_object_features.csv",
            run_dir / "wp3_nuclear_image_qc.csv",
            run_dir / "summary.json",
        )
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise RuntimeError(f"Shard {shard_id} is incomplete: {missing}")
        _, shard_manifest = read_csv(required[0])
        shard_manifest_keys.extend(row["image_key"] for row in shard_manifest)
        fields, rows = read_csv(required[1])
        object_fields = object_fields or fields
        if fields != object_fields:
            raise RuntimeError(f"Object schema differs in shard {shard_id}.")
        object_rows.extend(rows)
        fields, rows = read_csv(required[2])
        image_fields = image_fields or fields
        if fields != image_fields:
            raise RuntimeError(f"Image schema differs in shard {shard_id}.")
        image_rows.extend(rows)
        with open(required[3], encoding="utf-8") as handle:
            shard_summaries.append(json.load(handle))
        for row in shard_manifest:
            key = row["image_key"]
            stem = "".join(character if character.isalnum() or character in "-_." else "_" for character in key)
            mask_rows.append({
                "image_key": key,
                "shard_id": str(shard_id),
                "cell_mask": str((run_dir / "masks" / "cell" / f"{stem}_cell_masks.tif").resolve()),
                "nuclear_mask": str((run_dir / "masks" / "nuclear" / f"{stem}_nuclear_masks.tif").resolve()),
            })

    if len(shard_manifest_keys) != len(set(shard_manifest_keys)):
        raise RuntimeError("Shard manifests overlap.")
    if set(shard_manifest_keys) != set(master_keys):
        raise RuntimeError("Shard-manifest union does not equal the master manifest.")
    image_keys = [row["image_key"] for row in image_rows]
    if len(image_keys) != len(set(image_keys)) or set(image_keys) != set(master_keys):
        raise RuntimeError("Image-QC rows do not map one-to-one to the master manifest.")
    errors = [row for row in image_rows if row["status"] == "error"]
    if errors:
        raise RuntimeError(f"Cannot merge {len(errors)} errored images; first is {errors[0]['image_key']}: {errors[0]['error']}")
    object_keys = [(row["image_key"], int(row["object_id"])) for row in object_rows]
    if len(object_keys) != len(set(object_keys)):
        raise RuntimeError("Duplicate (image_key, object_id) rows across shards.")
    missing_masks = [row for row in mask_rows if not Path(row["cell_mask"]).exists() or not Path(row["nuclear_mask"]).exists()]
    no_cell_keys = {row["image_key"] for row in image_rows if row["status"] == "no_cells"}
    missing_masks = [row for row in missing_masks if row["image_key"] not in no_cell_keys]
    if missing_masks:
        raise RuntimeError(f"Missing masks for {len(missing_masks)} successful images; first: {missing_masks[0]['image_key']}")

    image_rows.sort(key=lambda row: order[row["image_key"]])
    object_rows.sort(key=lambda row: (order[row["image_key"]], int(row["object_id"])))
    mask_rows.sort(key=lambda row: order[row["image_key"]])
    output_dir = Path(args.output_dir)
    atomic_csv(output_dir / "wp3_nuclear_image_qc.csv", image_fields or [], image_rows)
    atomic_csv(output_dir / "wp3_nuclear_object_features.csv", object_fields or [], object_rows)
    atomic_csv(output_dir / "mask_manifest.csv", ["image_key", "shard_id", "cell_mask", "nuclear_mask"], mask_rows)
    summary = {
        "n_images": len(image_rows),
        "n_images_ok": sum(row["status"] == "ok" for row in image_rows),
        "n_images_no_cells": len(no_cell_keys),
        "n_cells": len(object_rows),
        "n_cells_with_nucleus": sum(int(row["nuclear_area_px"]) > 0 for row in object_rows),
        "n_shards": args.n_shards,
        "shard_summaries": shard_summaries,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp = output_dir / "summary.json.tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, output_dir / "summary.json")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
