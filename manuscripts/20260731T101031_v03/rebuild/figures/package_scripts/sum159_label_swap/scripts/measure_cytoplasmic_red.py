#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from skimage import measure


def find_project_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "scripts" / "agentRrunner.sh").is_file():
            return candidate
    raise RuntimeError(f"Could not locate the project root from {start}")


PACKAGE_ROOT = Path(__file__).resolve().parent.parent
ROOT = find_project_root(PACKAGE_ROOT)
OUT = PACKAGE_ROOT / "derived_data" / "cytoplasmic_red_quantification"
SEED = 20260729
WORKERS = 1
MIN_CYTOPLASM_PX = 20
CELL_EROSION_PX = 1
NUCLEAR_DILATION_PX = 2
BACKGROUND_CELL_DILATION_PX = 3
BACKGROUND_SIGMA_PX = 75.0
BACKGROUND_DOWNSAMPLE = 4
COMMON_G0 = {0.1, 0.25, 0.5, 1.0}

MONO_RUN = (
    ROOT
    / "data/image_processing_runs/sum159_targeted_nuclear_area/run_20260721_114700"
)
MONO_MANIFEST = (
    ROOT
    / "data/report_exports/sum159_targeted_nuclear_area/manifest_20260721/manifest.csv"
)
MONO_MASK_MANIFEST = MONO_RUN / "merged/mask_manifest.csv"
MONO_OBJECTS = MONO_RUN / "merged/wp3_nuclear_object_features.csv"
FULL_RUN = (
    ROOT
    / "data/image_processing_runs/full_segmentation_classification_nuclear/run_20260721_163410"
)
FULL_MANIFEST = FULL_RUN / "manifests/master_manifest.csv"

MIX_RUN = (
    ROOT
    / "data/image_processing_runs/sum159_competition_resegmentation/run_20260720_143217"
)
MIX_MANIFEST = MIX_RUN / "input_manifest.tsv"
MIX_CALLS = (
    PACKAGE_ROOT
    / "derived_data"
    / "competition_ploidy_cytoplasmic_gfp_full"
    / "sum159_competition_object_ploidy_calls.csv.gz"
)


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mono_red_path(row: dict[str, str]) -> Path:
    suffix = f"{row['plateID']}_{row['pos']}_{row['image_key'].rsplit('_', 1)[-1]}"
    if not row["image_key"].endswith("_" + suffix):
        raise ValueError(f"Could not parse monoculture image key: {row['image_key']}")
    prefix = row["image_key"][: -(len(suffix) + 1)]
    matches = sorted((ROOT / "all_raw").glob(f"{prefix}_alive*_{suffix}.tif"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one red/alive file for {row['image_key']}, found {matches}")
    return matches[0]


def select_images() -> list[dict[str, str]]:
    mask_by_key = {
        row["image_key"]: row for row in read_csv(MONO_MASK_MANIFEST)
    }
    monoculture = []
    for row in read_csv(MONO_MANIFEST):
        if not (
            row["cellLine"] in {"SUM-159-fuse", "SUM-159-chem"}
            and row["batch_id"] in {"C00", "I00", "M00b"}
            and row["ploidy"] in {"2N", "4N"}
            and float(row["G0"]) in COMMON_G0
            and 24.0 <= float(row["hours"]) <= 48.0
        ):
            continue
        mask_row = mask_by_key.get(row["image_key"])
        if mask_row is None:
            raise RuntimeError(f"No maintained mask mapping for {row['image_key']}")
        out = dict(row)
        out.update(
            {
                "assay": f"{row['batch_id']} monoculture",
                "ploidy_state": "low" if row["ploidy"] == "2N" else "high",
                "red_path": str(mono_red_path(row)),
                "cell_mask_path": mask_row["cell_mask"],
                "nuclear_mask_path": mask_row["nuclear_mask"],
                "G0_mM": row["G0"],
                "well": row["plateID"],
                "position": row["pos"],
                "sampling_role": "all maintained targeted fields meeting criteria",
            }
        )
        monoculture.append(out)

    strata = sorted({(float(row["G0"]), float(row["hours"])) for row in monoculture})
    rng = random.Random(SEED)
    mix_rows = read_csv(MIX_MANIFEST, delimiter="\t")
    mixture = []
    for g0, hours in strata:
        pool = [
            row
            for row in mix_rows
            if float(row["G0_mM"]) == g0 and float(row["hours"]) == hours
        ]
        if len(pool) < 2:
            raise RuntimeError(f"Need two mixture fields at G0={g0}, h={hours}")
        for row in rng.sample(pool, 2):
            out = dict(row)
            out.update(
                {
                    "assay": "J00 mixture",
                    "ploidy": "mixed",
                    "ploidy_state": "mixed",
                    "cell_mask_path": str(
                        MIX_RUN / "masks/cell" / f"{row['image_key']}_cell_masks.tif"
                    ),
                    "nuclear_mask_path": str(
                        MIX_RUN / "masks/nuclear" / f"{row['image_key']}_nuclear_masks.tif"
                    ),
                    "sampling_role": "two random fields per matched condition-time stratum",
                }
            )
            mixture.append(out)

    selected = monoculture + mixture
    keys = [row["image_key"] for row in selected]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Selected image keys are not unique.")
    return selected


def load_monoculture_eligible(
    selected_keys: set[str],
) -> tuple[dict[str, dict[int, dict[str, str]]], dict[str, dict[int, dict[str, str]]]]:
    full_manifest = {
        row["image_key"]: row for row in read_csv(FULL_MANIFEST)
        if row["image_key"] in selected_keys
    }
    if set(full_manifest) != selected_keys:
        missing = sorted(selected_keys - set(full_manifest))
        raise RuntimeError(f"Selected monoculture images absent from full manifest: {missing[:3]}")

    keys_by_shard: dict[str, set[str]] = defaultdict(set)
    for key, row in full_manifest.items():
        keys_by_shard[row["shard_id"]].add(key)
    full_rows: dict[str, dict[int, dict[str, str]]] = {key: {} for key in selected_keys}
    for shard_id, keys in keys_by_shard.items():
        path = FULL_RUN / "shards" / f"shard_{int(shard_id):03d}" / "objects.csv"
        for row in read_csv(path):
            if row["image_key"] in keys:
                full_rows[row["image_key"]][int(row["object_id"])] = row

    targeted_rows: dict[str, dict[int, dict[str, str]]] = {key: {} for key in selected_keys}
    for row in read_csv(MONO_OBJECTS):
        if row["image_key"] in selected_keys:
            targeted_rows[row["image_key"]][int(row["object_id"])] = row
    return full_rows, targeted_rows


def load_mixture_calls(
    selected_keys: set[str],
) -> dict[str, dict[int, dict[str, str]]]:
    rows: dict[str, dict[int, dict[str, str]]] = {key: {} for key in selected_keys}
    with gzip.open(MIX_CALLS, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["image_key"] in selected_keys:
                rows[row["image_key"]][int(row["object_id"])] = row
    if any(not value for value in rows.values()):
        missing = [key for key, value in rows.items() if not value]
        raise RuntimeError(f"No mixture calls for selected images: {missing[:3]}")
    return rows


def match_monoculture_objects(
    targeted: dict[int, dict[str, str]],
    full: dict[int, dict[str, str]],
    max_centroid_distance_px: float = 5.0,
    max_relative_area_difference: float = 0.25,
) -> tuple[dict[int, dict[str, str]], int]:
    target_ids = sorted(targeted)
    full_ids = sorted(full)
    target_xy = np.array(
        [
            [
                float(targeted[object_id]["centroid_x"]),
                float(targeted[object_id]["centroid_y"]),
            ]
            for object_id in target_ids
        ],
        dtype=np.float64,
    )
    full_xy = np.array(
        [
            [
                float(full[object_id]["centroid_x"]),
                float(full[object_id]["centroid_y"]),
            ]
            for object_id in full_ids
        ],
        dtype=np.float64,
    )
    squared_distance = np.sum(
        (target_xy[:, None, :] - full_xy[None, :, :]) ** 2,
        axis=2,
    )
    target_index, full_index = linear_sum_assignment(squared_distance)
    matched: dict[int, dict[str, str]] = {}
    for target_i, full_i in zip(target_index, full_index):
        target_id = target_ids[int(target_i)]
        full_id = full_ids[int(full_i)]
        distance = float(np.sqrt(squared_distance[target_i, full_i]))
        target_area = float(targeted[target_id]["cell_area_px"])
        full_area = float(full[full_id]["segmented_area_px"])
        relative_area_difference = abs(target_area - full_area) / max(
            target_area, full_area, 1.0
        )
        if (
            distance <= max_centroid_distance_px
            and relative_area_difference <= max_relative_area_difference
        ):
            matched[target_id] = full[full_id]
    return matched, len(target_ids) - len(matched)


def local_background(
    raw: np.ndarray, cell_masks: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    excluded = ndi.binary_dilation(
        cell_masks > 0, iterations=BACKGROUND_CELL_DILATION_PX
    )
    background_mask = ~excluded
    if int(background_mask.sum()) < 1000:
        raise RuntimeError("Insufficient non-cell background pixels.")
    raw32 = raw.astype(np.float32, copy=False)
    factor = BACKGROUND_DOWNSAMPLE
    raw_small = raw32[::factor, ::factor]
    background_small = background_mask[::factor, ::factor]
    weights = background_small.astype(np.float32)
    weighted_sum = ndi.gaussian_filter(
        raw_small * weights,
        sigma=BACKGROUND_SIGMA_PX / factor,
        mode="reflect",
    )
    smoothed_weight = ndi.gaussian_filter(
        weights,
        sigma=BACKGROUND_SIGMA_PX / factor,
        mode="reflect",
    )
    fallback = float(np.median(raw32[background_mask]))
    estimate_small = np.divide(
        weighted_sum,
        smoothed_weight,
        out=np.full_like(weighted_sum, fallback, dtype=np.float32),
        where=smoothed_weight > 1e-6,
    )
    zoom = (
        raw32.shape[0] / estimate_small.shape[0],
        raw32.shape[1] / estimate_small.shape[1],
    )
    estimate = ndi.zoom(estimate_small, zoom=zoom, order=1, mode="nearest")
    estimate = estimate[: raw32.shape[0], : raw32.shape[1]]
    if estimate.shape != raw32.shape:
        estimate = np.pad(
            estimate,
            (
                (0, raw32.shape[0] - estimate.shape[0]),
                (0, raw32.shape[1] - estimate.shape[1]),
            ),
            mode="edge",
        )
    residual = raw32 - estimate
    bg_values = (raw_small - estimate_small)[background_small]
    center = float(np.median(bg_values))
    mad = float(np.median(np.abs(bg_values - center)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 1e-8:
        q25, q75 = np.percentile(bg_values, [25, 75])
        sigma = float((q75 - q25) / 1.349)
    if not np.isfinite(sigma) or sigma <= 1e-8:
        sigma = float(np.std(bg_values))
    if not np.isfinite(sigma) or sigma <= 1e-8:
        raise RuntimeError("Could not estimate a positive background-noise scale.")
    return estimate, residual, sigma


def measure_image(
    row: dict[str, str],
    mono_full: dict[str, dict[int, dict[str, str]]],
    mono_targeted: dict[str, dict[int, dict[str, str]]],
    mix_calls: dict[str, dict[int, dict[str, str]]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    key = row["image_key"]
    raw = np.asarray(tifffile.imread(row["red_path"]))
    cell_masks = np.asarray(tifffile.imread(row["cell_mask_path"]))
    nuclear_masks = np.asarray(tifffile.imread(row["nuclear_mask_path"]))
    if raw.ndim != 2 or raw.shape != cell_masks.shape or raw.shape != nuclear_masks.shape:
        raise ValueError(
            f"Image/mask shape mismatch for {key}: {raw.shape}, "
            f"{cell_masks.shape}, {nuclear_masks.shape}"
        )

    if row["assay"] == "J00 mixture":
        eligible = {
            object_id: call
            for object_id, call in mix_calls[key].items()
            if call["predicted_label_name"] == "alive"
            and call["cell_area_pass"] == "1"
            and int(float(call["nuclear_area_px"])) > 0
            and call["ploidy_call"] in {"low", "high"}
        }
    else:
        full = mono_full[key]
        targeted = mono_targeted[key]
        for object_id, target in targeted.items():
            mask_area = int(np.sum(cell_masks == object_id))
            if mask_area != int(target["cell_area_px"]):
                raise RuntimeError(f"Saved mask area mismatch for {key}, object {object_id}")
        matched, n_classifier_unmatched = match_monoculture_objects(targeted, full)
        eligible = {
            object_id: source
            for object_id, source in matched.items()
            if source["predicted_label_name"] == "alive"
            and source["cell_area_pass"] == "1"
            and int(float(targeted[object_id]["nuclear_area_px"])) > 0
        }
    if row["assay"] == "J00 mixture":
        n_classifier_unmatched = 0

    background, residual, background_sigma = local_background(raw, cell_masks)
    raw32 = raw.astype(np.float32, copy=False)
    p1, p99 = (float(value) for value in np.percentile(raw32, [1, 99]))
    image01 = np.clip((raw32 - p1) / max(p99 - p1, 1e-8), 0.0, 1.0)
    props = {int(prop.label): prop for prop in measure.regionprops(cell_masks)}

    object_rows = []
    n_missing_property = 0
    n_no_nucleus_mask = 0
    n_small_cytoplasm = 0
    for object_id, source in eligible.items():
        prop = props.get(object_id)
        if prop is None:
            n_missing_property += 1
            continue
        y0, x0, y1, x1 = (int(value) for value in prop.bbox)
        local_cell = prop.image.astype(bool, copy=False)
        eroded_cell = ndi.binary_erosion(
            local_cell, iterations=CELL_EROSION_PX, border_value=0
        )
        local_nucleus = nuclear_masks[y0:y1, x0:x1] == object_id
        if not np.any(local_nucleus):
            n_no_nucleus_mask += 1
            continue
        excluded_nucleus = ndi.binary_dilation(
            local_nucleus, iterations=NUCLEAR_DILATION_PX
        )
        cytoplasm = eroded_cell & ~excluded_nucleus
        n_cytoplasm = int(cytoplasm.sum())
        if n_cytoplasm < MIN_CYTOPLASM_PX:
            n_small_cytoplasm += 1
            continue
        local_raw = raw32[y0:y1, x0:x1][cytoplasm]
        local_background_values = background[y0:y1, x0:x1][cytoplasm]
        local_residual = residual[y0:y1, x0:x1][cytoplasm]
        local_image01 = image01[y0:y1, x0:x1][cytoplasm]
        ploidy_state = (
            source["ploidy_call"]
            if row["assay"] == "J00 mixture"
            else row["ploidy_state"]
        )
        object_rows.append(
            {
                "assay": row["assay"],
                "image_key": key,
                "object_id": object_id,
                "ploidy_state": ploidy_state,
                "G0_mM": float(row["G0_mM"]),
                "hours": float(row["hours"]),
                "well": row["well"],
                "position": row["position"],
                "cytoplasm_pixel_count": n_cytoplasm,
                "red_cytoplasm_raw_median": float(np.median(local_raw)),
                "red_local_background_median": float(np.median(local_background_values)),
                "red_cytoplasm_bg_corrected_median": float(np.median(local_residual)),
                "red_cytoplasm_bg_z": float(np.median(local_residual) / background_sigma),
                "red_cytoplasm_image_p1_p99_median": float(np.median(local_image01)),
                "image_background_noise_sigma": background_sigma,
                "image_red_p1": p1,
                "image_red_p99": p99,
            }
        )
    image_summary = {
        "assay": row["assay"],
        "image_key": key,
        "G0_mM": row["G0_mM"],
        "hours": row["hours"],
        "well": row["well"],
        "position": row["position"],
        "red_path": row["red_path"],
        "cell_mask_path": row["cell_mask_path"],
        "nuclear_mask_path": row["nuclear_mask_path"],
        "n_segmented": int(cell_masks.max()),
        "n_eligible_before_cytoplasm_gate": len(eligible),
        "n_measured": len(object_rows),
        "n_measured_low": sum(value["ploidy_state"] == "low" for value in object_rows),
        "n_measured_high": sum(value["ploidy_state"] == "high" for value in object_rows),
        "n_missing_property": n_missing_property,
        "n_classifier_unmatched": n_classifier_unmatched,
        "n_no_nucleus_mask": n_no_nucleus_mask,
        "n_small_cytoplasm": n_small_cytoplasm,
        "image_background_noise_sigma": background_sigma,
        "image_red_p1": p1,
        "image_red_p99": p99,
        "sampling_role": row["sampling_role"],
    }
    return object_rows, image_summary


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--n-chunks", type=int, default=6)
    parser.add_argument("--aggregate", action="store_true")
    return parser.parse_args()


def aggregate_chunks(n_chunks: int) -> None:
    all_cells: list[dict[str, str]] = []
    all_qc: list[dict[str, str]] = []
    all_selection: list[dict[str, str]] = []
    chunk_metadata = []
    for index in range(n_chunks):
        suffix = f".chunk_{index:02d}"
        cell_path = OUT / f"per_cell_cytoplasmic_red{suffix}.csv.gz"
        qc_path = OUT / f"image_measurement_qc{suffix}.csv"
        selection_path = OUT / f"sample_manifest{suffix}.csv"
        metadata_path = OUT / f"run_metadata{suffix}.json"
        for path in (cell_path, qc_path, selection_path, metadata_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        with gzip.open(cell_path, "rt", newline="", encoding="utf-8") as handle:
            all_cells.extend(csv.DictReader(handle))
        all_qc.extend(read_csv(qc_path))
        all_selection.extend(read_csv(selection_path))
        chunk_metadata.append(json.loads(metadata_path.read_text(encoding="utf-8")))
    all_cells.sort(
        key=lambda row: (row["assay"], row["image_key"], int(row["object_id"]))
    )
    all_qc.sort(
        key=lambda row: (
            row["assay"],
            float(row["G0_mM"]),
            float(row["hours"]),
            row["image_key"],
        )
    )
    all_selection.sort(
        key=lambda row: (
            row["assay"],
            float(row["G0_mM"]),
            float(row["hours"]),
            row["image_key"],
        )
    )
    with gzip.open(
        OUT / "per_cell_cytoplasmic_red.csv.gz",
        "wt",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_cells[0]))
        writer.writeheader()
        writer.writerows(all_cells)
    write_csv(OUT / "image_measurement_qc.csv", all_qc, list(all_qc[0]))
    write_csv(OUT / "sample_manifest.csv", all_selection, list(all_selection[0]))
    selected = select_images()
    metadata = dict(chunk_metadata[0])
    metadata.pop("chunk", None)
    metadata.update(
        {
            "workers": WORKERS,
            "n_images": len(all_selection),
            "n_monoculture_images": sum(
                row["assay"] != "J00 mixture" for row in all_selection
            ),
            "n_mixture_images": sum(
                row["assay"] == "J00 mixture" for row in all_selection
            ),
            "n_measured_cells": len(all_cells),
            "execution": {
                "n_chunks": n_chunks,
                "chunk_image_counts": [
                    value["n_images"] for value in chunk_metadata
                ],
            },
            "conditions": {
                "G0_mM": sorted(COMMON_G0),
                "hours_inclusive": [24, 48],
                "realized_strata": sorted(
                    {
                        f"G0={float(row['G0_mM']):g},h={float(row['hours']):g}"
                        for row in selected
                    }
                ),
            },
        }
    )
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Aggregated {len(all_cells)} cells from {len(all_selection)} images "
        f"across {n_chunks} chunks.",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    if args.aggregate:
        aggregate_chunks(args.n_chunks)
        return
    if args.chunk_index is not None and not 0 <= args.chunk_index < args.n_chunks:
        raise ValueError("--chunk-index must be in [0, --n-chunks).")
    selected = select_images()
    suffix = ""
    if args.chunk_index is not None:
        selected = [
            row for index, row in enumerate(selected)
            if index % args.n_chunks == args.chunk_index
        ]
        suffix = f".chunk_{args.chunk_index:02d}"
    if not selected:
        raise RuntimeError("This chunk contains no selected images.")
    mono_keys = {row["image_key"] for row in selected if row["assay"] != "J00 mixture"}
    mix_keys = {row["image_key"] for row in selected if row["assay"] == "J00 mixture"}
    mono_full, mono_targeted = (
        load_monoculture_eligible(mono_keys) if mono_keys else ({}, {})
    )
    mix_calls = load_mixture_calls(mix_keys) if mix_keys else {}

    per_cell: list[dict[str, object]] = []
    image_summaries: list[dict[str, object]] = []
    errors = []
    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = {
            executor.submit(
                measure_image, row, mono_full, mono_targeted, mix_calls
            ): row
            for row in selected
        }
        for future in as_completed(futures):
            row = futures[future]
            try:
                cells, summary = future.result()
            except Exception as exc:
                errors.append((row["image_key"], repr(exc)))
            else:
                per_cell.extend(cells)
                image_summaries.append(summary)
                print(
                    f"{row['assay']}: {row['image_key']} -> {len(cells)} cells",
                    flush=True,
                )
    if errors:
        raise RuntimeError(f"{len(errors)} image failures; first={errors[0]}")
    if not per_cell:
        raise RuntimeError("No cells were measured.")
    per_cell.sort(key=lambda row: (row["assay"], row["image_key"], int(row["object_id"])))
    image_summaries.sort(key=lambda row: (row["assay"], float(row["G0_mM"]), float(row["hours"]), row["image_key"]))

    cell_fields = list(per_cell[0])
    with gzip.open(
        OUT / f"per_cell_cytoplasmic_red{suffix}.csv.gz",
        "wt",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields)
        writer.writeheader()
        writer.writerows(per_cell)
    write_csv(
        OUT / f"image_measurement_qc{suffix}.csv",
        image_summaries,
        list(image_summaries[0]),
    )

    selection_rows = []
    summary_by_key = {row["image_key"]: row for row in image_summaries}
    for row in selected:
        summary = summary_by_key[row["image_key"]]
        selection_rows.append(
            {
                "assay": row["assay"],
                "image_key": row["image_key"],
                "ploidy": row["ploidy"],
                "G0_mM": row["G0_mM"],
                "hours": row["hours"],
                "well": row["well"],
                "position": row["position"],
                "red_path": row["red_path"],
                "cell_mask_path": row["cell_mask_path"],
                "nuclear_mask_path": row["nuclear_mask_path"],
                "sampling_role": row["sampling_role"],
                "n_measured": summary["n_measured"],
                "n_measured_low": summary["n_measured_low"],
                "n_measured_high": summary["n_measured_high"],
            }
        )
    write_csv(
        OUT / f"sample_manifest{suffix}.csv",
        selection_rows,
        list(selection_rows[0]),
    )

    metadata = {
        "random_seed": SEED,
        "workers": WORKERS,
        "n_images": len(selected),
        "n_monoculture_images": len(mono_keys),
        "n_mixture_images": len(mix_keys),
        "n_measured_cells": len(per_cell),
        "chunk": (
            {
                "index": args.chunk_index,
                "n_chunks": args.n_chunks,
            }
            if args.chunk_index is not None
            else None
        ),
        "conditions": {
            "G0_mM": sorted(COMMON_G0),
            "hours_inclusive": [24, 48],
            "realized_strata": sorted(
                {
                    f"G0={float(row['G0_mM']):g},h={float(row['hours']):g}"
                    for row in selected
                }
            ),
        },
        "compartment": {
            "cell_erosion_px": CELL_EROSION_PX,
            "nuclear_dilation_px": NUCLEAR_DILATION_PX,
            "minimum_cytoplasm_pixels": MIN_CYTOPLASM_PX,
        },
        "background_normalization": {
            "background_cell_dilation_px": BACKGROUND_CELL_DILATION_PX,
            "gaussian_sigma_px": BACKGROUND_SIGMA_PX,
            "background_grid_downsample": BACKGROUND_DOWNSAMPLE,
            "primary_metric": "median(red - local_background) / image background residual robust sigma",
            "robust_sigma": "1.4826 * MAD of non-cell residual pixels",
            "secondary_metrics": [
                "raw cytoplasm median",
                "background-corrected cytoplasm median",
                "per-image p1-p99 normalized cytoplasm median",
            ],
        },
        "inputs": {
            "monoculture_manifest": str(MONO_MANIFEST.relative_to(ROOT)),
            "monoculture_mask_manifest": str(MONO_MASK_MANIFEST.relative_to(ROOT)),
            "monoculture_full_classifier_run": str(FULL_RUN.relative_to(ROOT)),
            "mixture_manifest": str(MIX_MANIFEST.relative_to(ROOT)),
            "mixture_calls": str(MIX_CALLS.relative_to(ROOT)),
        },
        "input_hashes": {
            "monoculture_manifest": sha256(MONO_MANIFEST),
            "monoculture_mask_manifest": sha256(MONO_MASK_MANIFEST),
            "mixture_manifest": sha256(MIX_MANIFEST),
            "mixture_calls": sha256(MIX_CALLS),
        },
    }
    (OUT / f"run_metadata{suffix}.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
