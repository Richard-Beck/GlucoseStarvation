#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from pathlib import Path


TARGETS = (
    ("SUM-159-fuse", "C00-IncucyteRawDataLiveDead-varyGlucose", "C00"),
    ("SUM-159-fuse", "I00-IncucyteRawDataLiveDead-varyGlucose", "I00"),
    ("SUM-159-chem", "M00b-IncucyteRawDataLiveDead-varyGlucose", "M00b"),
)
CHANNEL_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<channel>alive[^_]*|dead|phase)_(?P<suffix>[^.]+)\.tif{1,2}$",
    re.IGNORECASE,
)
OUTPUT_FIELDS = (
    "manifest_index", "shard_id", "cellLine", "experiment", "batch_id",
    "plateID", "pos", "ploidy", "G0", "glucose_bin", "time_bin",
    "hours", "image_key", "base_key", "alive_count", "dead_count",
    "total_count_kept", "total_area_px_all_objects",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prelim_counts_csv",
        default="data/image_processing_runs/run_20260324_233122/prelim_counts_by_image.csv",
    )
    parser.add_argument("--raw_data_dir", default="all_raw")
    parser.add_argument(
        "--prior_pilot_manifest",
        default="data/image_processing_runs/wp3_nuclear_size_pilot/run_20260511_143515/pilot_manifest.csv",
    )
    parser.add_argument(
        "--output_dir",
        default="data/report_exports/sum159_targeted_nuclear_area/manifest_20260721",
    )
    parser.add_argument("--timepoints_per_batch", type=int, default=5)
    parser.add_argument("--fields_per_condition", type=int, default=2)
    parser.add_argument("--n_shards", type=int, default=4)
    return parser.parse_args()


def glucose_bin(value: float) -> str:
    if value in (0.0, 0.1):
        return "low"
    if value in (0.25, 0.5, 1.0):
        return "mid"
    if value in (5.0, 25.0):
        return "high"
    raise ValueError(f"Unexpected G0 value: {value}")


def evenly_spaced(values: list[float], n: int) -> list[float]:
    values = sorted(set(values))
    if len(values) < n:
        raise ValueError(f"Need {n} shared timepoints but only found {len(values)}: {values}")
    if n == 1:
        return [values[len(values) // 2]]
    lo, hi = values[0], values[-1]
    targets = [lo + i * (hi - lo) / (n - 1) for i in range(n)]
    selected = [min(values, key=lambda value: (abs(value - target), value)) for target in targets]
    if len(set(selected)) != n:
        raise ValueError(f"Could not select {n} distinct evenly spaced times from {values}")
    return selected


def read_csv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def raw_channel_coverage(raw_data_dir: str, keys: set[str]) -> dict[str, set[str]]:
    coverage = {key: set() for key in keys}
    with os.scandir(raw_data_dir) as entries:
        for entry in entries:
            match = CHANNEL_RE.match(entry.name)
            if match is None:
                continue
            key = f"{match.group('prefix')}_{match.group('suffix')}"
            if key not in coverage:
                continue
            channel = match.group("channel").lower()
            coverage[key].add("alive" if channel.startswith("alive") else channel)
    return coverage


def main() -> None:
    args = parse_args()
    if args.timepoints_per_batch < 2 or args.fields_per_condition < 1 or args.n_shards != 4:
        raise ValueError("This approved design requires >=2 timepoints, >=1 field, and exactly four shards.")

    target_lookup = {(line, experiment): batch for line, experiment, batch in TARGETS}
    candidates = []
    for row in read_csv(args.prelim_counts_csv):
        batch = target_lookup.get((row["cellLine"], row["experiment"]))
        if batch is None or row["ploidy"] not in ("2N", "4N"):
            continue
        out = dict(row)
        out["batch_id"] = batch
        out["hours_num"] = float(row["hours"])
        out["G0_num"] = float(row["G0"])
        out["total_count_num"] = int(float(row["total_count_kept"]))
        candidates.append(out)
    if not candidates:
        raise RuntimeError("No target SUM-159 rows found.")

    selected: list[dict[str, object]] = []
    selected_hours: dict[str, list[float]] = {}
    shortages: list[dict[str, object]] = []
    for line, experiment, batch in TARGETS:
        batch_rows = [r for r in candidates if r["cellLine"] == line and r["experiment"] == experiment]
        conditions = sorted({(r["ploidy"], r["G0_num"]) for r in batch_rows}, key=lambda x: (x[0], x[1]))
        hour_sets = []
        for ploidy, g0 in conditions:
            hour_sets.append({r["hours_num"] for r in batch_rows if r["ploidy"] == ploidy and r["G0_num"] == g0})
        shared_hours = sorted(set.intersection(*hour_sets))
        chosen_hours = evenly_spaced(shared_hours, args.timepoints_per_batch)
        selected_hours[batch] = chosen_hours

        for ploidy, g0 in conditions:
            condition_rows = [
                r for r in batch_rows
                if r["ploidy"] == ploidy and r["G0_num"] == g0 and r["hours_num"] in chosen_hours
            ]
            by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in condition_rows:
                by_series[row["base_key"]].append(row)
            complete = []
            for base_key, rows in by_series.items():
                if {r["hours_num"] for r in rows} == set(chosen_hours) and len(rows) == len(chosen_hours):
                    complete.append((base_key, rows))
            by_well: dict[str, list[tuple[str, list[dict[str, str]]]]] = defaultdict(list)
            for item in complete:
                by_well[str(item[1][0]["plateID"])].append(item)
            well_ids = sorted(by_well)
            if len(well_ids) >= args.fields_per_condition:
                well_indices = [round(q * (len(well_ids) - 1)) for q in (0.25, 0.75)]
                if well_indices[0] == well_indices[1]:
                    well_indices[1] = min(len(well_ids) - 1, well_indices[0] + 1)
                preferred_positions = (2.0, 3.0)
                chosen_series = []
                for selection_index, well_index in enumerate(well_indices[: args.fields_per_condition]):
                    choices = by_well[well_ids[well_index]]
                    preferred = preferred_positions[min(selection_index, len(preferred_positions) - 1)]
                    choices.sort(
                        key=lambda item: (
                            abs(float(item[1][0]["pos"]) - preferred),
                            float(item[1][0]["pos"]),
                            item[0],
                        )
                    )
                    chosen_series.append(choices[0])
            else:
                complete.sort(key=lambda item: item[0])
                chosen_series = complete[: args.fields_per_condition]
            if len(chosen_series) < args.fields_per_condition:
                shortages.append({
                    "batch_id": batch, "ploidy": ploidy, "G0": g0,
                    "requested_fields": args.fields_per_condition,
                    "available_complete_fields": len(complete),
                })
            for _, rows in chosen_series:
                for row in rows:
                    hour_rank = chosen_hours.index(row["hours_num"])
                    out = dict(row)
                    out["G0"] = f"{g0:g}"
                    out["hours"] = f"{row['hours_num']:g}"
                    out["glucose_bin"] = glucose_bin(g0)
                    out["time_bin"] = "early" if hour_rank < 2 else ("mid" if hour_rank == 2 else "late")
                    selected.append(out)

    selected.sort(
        key=lambda r: (
            [x[2] for x in TARGETS].index(str(r["batch_id"])),
            str(r["ploidy"]), float(r["G0"]), float(r["hours"]), str(r["base_key"]),
        )
    )
    keys = [str(row["image_key"]) for row in selected]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Selected manifest contains duplicate image_key values.")

    coverage = raw_channel_coverage(args.raw_data_dir, set(keys))
    incomplete_raw = {key: sorted({"alive", "dead", "phase"} - channels) for key, channels in coverage.items() if channels != {"alive", "dead", "phase"}}
    if incomplete_raw:
        first = next(iter(incomplete_raw.items()))
        raise RuntimeError(f"Incomplete raw channels for {len(incomplete_raw)} selected images; first: {first}")

    for index, row in enumerate(selected, start=1):
        row["manifest_index"] = index
        row["shard_id"] = (index - 1) % args.n_shards

    output_dir = Path(args.output_dir)
    master_path = output_dir / "manifest.csv"
    write_csv(master_path, selected)
    shard_paths = []
    for shard_id in range(args.n_shards):
        path = output_dir / f"shard_{shard_id:03d}.csv"
        write_csv(path, [row for row in selected if row["shard_id"] == shard_id])
        shard_paths.append(path)

    prior_rows = read_csv(args.prior_pilot_manifest)
    smoke_rows = []
    for line, experiment, batch in TARGETS:
        for ploidy in ("2N", "4N"):
            rows = [
                r for r in prior_rows
                if r["cellLine"] == line and r["experiment"] == experiment and r["ploidy"] == ploidy
            ]
            if not rows:
                raise RuntimeError(f"No prior pilot rows for {batch} {ploidy}")
            rows.sort(key=lambda r: (float(r["hours"]), float(r["G0"]), r["image_key"]))
            for row in (rows[0], rows[-1]):
                out = dict(row)
                out["batch_id"] = batch
                smoke_rows.append(out)
    smoke_rows.sort(key=lambda r: str(r["image_key"]))
    for index, row in enumerate(smoke_rows, start=1):
        row["manifest_index"] = index
        row["shard_id"] = (index - 1) % args.n_shards
    smoke_path = output_dir / "smoke_manifest.csv"
    write_csv(smoke_path, smoke_rows)
    smoke_shard_paths = []
    for shard_id in range(args.n_shards):
        path = output_dir / f"smoke_shard_{shard_id:03d}.csv"
        write_csv(path, [row for row in smoke_rows if row["shard_id"] == shard_id])
        smoke_shard_paths.append(path)

    metadata = {
        "design": {
            "targets": [dict(cellLine=line, experiment=experiment, batch_id=batch) for line, experiment, batch in TARGETS],
            "timepoints_per_batch": args.timepoints_per_batch,
            "fields_per_condition": args.fields_per_condition,
            "n_shards": args.n_shards,
            "selected_hours": selected_hours,
        },
        "n_selected_images": len(selected),
        "n_smoke_images": len(smoke_rows),
        "n_shortages": len(shortages),
        "shortages": shortages,
        "n_incomplete_raw_channel_groups": len(incomplete_raw),
        "master_manifest": str(master_path),
        "master_manifest_sha256": sha256(master_path),
        "shards": [{"path": str(path), "n_images": sum(r["shard_id"] == i for r in selected), "sha256": sha256(path)} for i, path in enumerate(shard_paths)],
        "smoke_manifest": str(smoke_path),
        "smoke_shards": [{"path": str(path), "n_images": sum(r["shard_id"] == i for r in smoke_rows), "sha256": sha256(path)} for i, path in enumerate(smoke_shard_paths)],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    meta_path = output_dir / "manifest_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
