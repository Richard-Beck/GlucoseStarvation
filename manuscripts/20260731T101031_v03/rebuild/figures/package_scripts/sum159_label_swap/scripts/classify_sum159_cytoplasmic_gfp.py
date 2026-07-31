#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import platform
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.special import logsumexp
from skimage.segmentation import find_boundaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify SUM-159 ploidy from locally background-corrected "
            "cytoplasmic GFP without assuming a 50:50 baseline mixture."
        )
    )
    parser.add_argument(
        "--run_dir",
        default=(
            "data/image_processing_runs/sum159_competition_resegmentation/"
            "run_20260720_143217"
        ),
    )
    parser.add_argument(
        "--output_dir",
        default=str(
            Path(__file__).resolve().parent.parent
            / "derived_data"
            / "competition_ploidy_cytoplasmic_gfp_full"
        ),
    )
    parser.add_argument(
        "--target_hours",
        default="24,60,120",
        help=(
            "Selected QC hours. Baseline images plus these representative "
            "images are processed unless --all_images is supplied."
        ),
    )
    parser.add_argument("--all_images", action="store_true")
    parser.add_argument("--background_sigma_px", type=float, default=75.0)
    parser.add_argument("--background_cell_dilation_px", type=int, default=3)
    parser.add_argument("--cell_erosion_px", type=int, default=1)
    parser.add_argument("--nuclear_dilation_px", type=int, default=2)
    parser.add_argument("--min_cytoplasm_px", type=int, default=20)
    parser.add_argument("--asinh_cofactor", type=float, default=1.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260721)
    return parser.parse_args()


def select_examples(manifest: pd.DataFrame, target_hours: list[float]) -> pd.DataFrame:
    selected: list[pd.Series] = []
    for g0 in sorted(manifest["G0_mM"].unique()):
        group = manifest.loc[manifest["G0_mM"] == g0]
        series_index = int(group["series_index"].min())
        series = group.loc[group["series_index"] == series_index]
        acquired = sorted(series["hours"].unique())
        for target in target_hours:
            actual = min(acquired, key=lambda hour: (abs(float(hour) - target), hour))
            match = series.loc[series["hours"] == actual]
            if len(match) != 1:
                raise RuntimeError(
                    f"Expected one image for G0={g0}, series={series_index}, "
                    f"hour={actual}; found {len(match)}"
                )
            selected.append(match.iloc[0])
    return pd.DataFrame(selected).reset_index(drop=True)


def read_image(path: Path) -> np.ndarray:
    image = tifffile.imread(path)
    if image.ndim != 2:
        raise ValueError(f"Expected a two-dimensional TIFF at {path}; got {image.shape}")
    return image


def normalized_gaussian_background(
    image: np.ndarray,
    background_mask: np.ndarray,
    sigma: float,
) -> tuple[np.ndarray, np.ndarray]:
    if sigma <= 0:
        raise ValueError("--background_sigma_px must be positive")
    image = np.asarray(image, dtype=np.float32)
    weights = np.asarray(background_mask, dtype=np.float32)
    if image.shape != weights.shape:
        raise ValueError("GFP image and background mask have different shapes")
    if not np.any(weights):
        raise ValueError("No non-cell pixels remain for GFP background estimation")
    weighted_sum = ndi.gaussian_filter(image * weights, sigma=sigma, mode="reflect")
    support = ndi.gaussian_filter(weights, sigma=sigma, mode="reflect")
    fallback = float(np.median(image[background_mask]))
    background = np.divide(
        weighted_sum,
        support,
        out=np.full_like(weighted_sum, fallback, dtype=np.float32),
        where=support > 1e-6,
    )
    return background, support


def label_median(values: np.ndarray, labels: np.ndarray, maximum: int) -> np.ndarray:
    indices = np.arange(1, maximum + 1, dtype=np.int32)
    result = np.asarray(ndi.median(values, labels=labels, index=indices), dtype=float)
    if result.ndim == 0:
        result = result.reshape(1)
    return result


def extract_image_features(task: dict[str, object]) -> pd.DataFrame:
    image_key = str(task["image_key"])
    green = read_image(Path(str(task["green_path"]))).astype(np.float32, copy=False)
    cell = read_image(Path(str(task["cell_path"]))).astype(np.int32, copy=False)
    nuclear = read_image(Path(str(task["nuclear_path"]))).astype(np.int32, copy=False)
    if green.shape != cell.shape or cell.shape != nuclear.shape:
        raise ValueError(f"Image/mask shape mismatch for {image_key}")

    maximum = int(cell.max())
    if maximum == 0:
        return pd.DataFrame()

    cell_pixels = cell > 0
    background_exclusion = cell_pixels
    dilation = int(task["background_cell_dilation_px"])
    if dilation > 0:
        background_exclusion = ndi.binary_dilation(
            background_exclusion,
            iterations=dilation,
        )
    background_mask = ~background_exclusion
    background, support = normalized_gaussian_background(
        green,
        background_mask,
        sigma=float(task["background_sigma_px"]),
    )
    corrected = green - background

    cell_boundary = find_boundaries(cell, mode="inner")
    erosion = int(task["cell_erosion_px"])
    if erosion <= 0:
        cell_interior = cell_pixels
    else:
        removed = cell_boundary
        if erosion > 1:
            removed = ndi.binary_dilation(removed, iterations=erosion - 1)
        cell_interior = cell_pixels & ~removed

    nuclear_exclusion = nuclear > 0
    nuclear_dilation = int(task["nuclear_dilation_px"])
    if nuclear_dilation > 0:
        nuclear_exclusion = ndi.binary_dilation(
            nuclear_exclusion,
            iterations=nuclear_dilation,
        )
    cytoplasm = cell_interior & ~nuclear_exclusion
    cytoplasm_labels = np.where(cytoplasm, cell, 0)

    counts = np.bincount(cytoplasm_labels.ravel(), minlength=maximum + 1)[1:]
    raw_median = label_median(green, cytoplasm_labels, maximum)
    corrected_median = label_median(corrected, cytoplasm_labels, maximum)
    background_median = label_median(background, cytoplasm_labels, maximum)
    support_median = label_median(support, cytoplasm_labels, maximum)
    missing = counts == 0
    raw_median[missing] = np.nan
    corrected_median[missing] = np.nan
    background_median[missing] = np.nan
    support_median[missing] = np.nan

    return pd.DataFrame(
        {
            "image_key": image_key,
            "object_id": np.arange(1, maximum + 1, dtype=np.int32),
            "cytoplasm_pixel_count": counts.astype(np.int32),
            "green_cytoplasmic_median_raw": raw_median,
            "green_local_background_median": background_median,
            "green_cytoplasmic_median_bg_corrected": corrected_median,
            "green_background_support_median": support_median,
        }
    )


def normal_log_density(x: np.ndarray, mean: float, variance: float) -> np.ndarray:
    return -0.5 * (
        math.log(2 * math.pi * variance) + ((x - mean) ** 2) / variance
    )


def fit_grouped_gaussian_mixture(
    values: np.ndarray,
    groups: np.ndarray,
    seed: int,
    max_iter: int = 500,
    tolerance: float = 1e-8,
) -> dict[str, object]:
    finite = np.isfinite(values)
    x = np.asarray(values[finite], dtype=float)
    group_names, group_index = np.unique(np.asarray(groups)[finite], return_inverse=True)
    if len(x) < 1000 or np.std(x) <= 0:
        raise ValueError("Insufficient baseline variation for a two-component mixture")

    overall_variance = float(np.var(x, ddof=1))
    variance_floor = max(overall_variance * 1e-6, np.finfo(float).eps)
    rng = np.random.default_rng(seed)
    starts = [
        np.quantile(x, [0.20, 0.80]),
        np.quantile(x, [0.30, 0.70]),
        np.array([np.mean(x[x <= np.median(x)]), np.mean(x[x > np.median(x)])]),
    ]
    for _ in range(3):
        starts.append(np.sort(rng.choice(x, size=2, replace=False)))

    group_counts = np.bincount(group_index).astype(float)
    fits: list[dict[str, object]] = []
    for initial_means in starts:
        means = np.sort(np.asarray(initial_means, dtype=float))
        variances = np.repeat(overall_variance / 2, 2)
        weights = np.repeat([[0.5, 0.5]], len(group_names), axis=0)
        previous = -np.inf
        converged = False
        collapsed = False

        for iteration in range(1, max_iter + 1):
            log_components = np.column_stack(
                [
                    normal_log_density(x, means[k], variances[k])
                    + np.log(weights[group_index, k])
                    for k in range(2)
                ]
            )
            denominator = logsumexp(log_components, axis=1)
            log_likelihood = float(denominator.sum())
            responsibilities = np.exp(log_components - denominator[:, None])
            effective_n = responsibilities.sum(axis=0)
            if np.any(effective_n < max(25, 0.005 * len(x))):
                collapsed = True
                break

            means = (responsibilities * x[:, None]).sum(axis=0) / effective_n
            variances = np.maximum(
                (
                    responsibilities
                    * ((x[:, None] - means[None, :]) ** 2)
                ).sum(axis=0)
                / effective_n,
                variance_floor,
            )
            for k in range(2):
                weights[:, k] = np.bincount(
                    group_index,
                    weights=responsibilities[:, k],
                    minlength=len(group_names),
                ) / group_counts
            weights = np.clip(weights, 1e-4, 1 - 1e-4)
            weights /= weights.sum(axis=1, keepdims=True)

            order = np.argsort(means)
            means = means[order]
            variances = variances[order]
            weights = weights[:, order]
            if (
                np.isfinite(previous)
                and abs(log_likelihood - previous)
                <= tolerance * (1 + abs(previous))
            ):
                converged = True
                break
            previous = log_likelihood
        else:
            log_likelihood = float(denominator.sum())

        if not collapsed and np.all(np.isfinite(means)) and means[0] < means[1]:
            fits.append(
                {
                    "means": means,
                    "variances": variances,
                    "weights": weights,
                    "group_names": group_names,
                    "log_likelihood": log_likelihood,
                    "iterations": iteration,
                    "converged": converged,
                }
            )
    if not fits:
        raise RuntimeError("All free-proportion mixture starts collapsed")
    return max(fits, key=lambda fit: float(fit["log_likelihood"]))


def central_likelihood_boundary(
    low_mean: float,
    low_variance: float,
    high_mean: float,
    high_variance: float,
) -> float:
    low_sd = math.sqrt(low_variance)
    high_sd = math.sqrt(high_variance)
    coefficients = np.array(
        [
            0.5 / low_variance - 0.5 / high_variance,
            -low_mean / low_variance + high_mean / high_variance,
            (
                low_mean**2 / (2 * low_variance)
                - high_mean**2 / (2 * high_variance)
                - math.log(high_sd / low_sd)
            ),
        ]
    )
    if abs(coefficients[0]) < 1e-12:
        roots = np.array([-coefficients[2] / coefficients[1]])
    else:
        roots = np.roots(coefficients)
    real_roots = np.real(roots[np.isreal(roots)])
    between = real_roots[
        (real_roots >= low_mean - 1e-10) & (real_roots <= high_mean + 1e-10)
    ]
    if len(between) != 1:
        raise RuntimeError(
            "Expected one equal-likelihood boundary between component means; "
            f"found {between.tolist()}"
        )
    return float(between[0])


def git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.workers < 1:
        raise ValueError("--workers must be at least one")
    if args.asinh_cofactor <= 0:
        raise ValueError("--asinh_cofactor must be positive")

    manifest_path = run_dir / "input_manifest.tsv"
    objects_path = run_dir / "object_features_and_predictions.csv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    objects = pd.read_csv(objects_path)
    calibration_hour = float(manifest["hours"].min())
    targets = [float(value) for value in args.target_hours.split(",") if value.strip()]
    examples = select_examples(manifest, targets)
    if args.all_images:
        selected_manifest = manifest.copy()
    else:
        selected_keys = set(examples["image_key"])
        selected_manifest = manifest.loc[
            (manifest["hours"] == calibration_hour)
            | manifest["image_key"].isin(selected_keys)
        ].copy()
    selected_manifest = selected_manifest.drop_duplicates("image_key")

    tasks: list[dict[str, object]] = []
    for row in selected_manifest.itertuples(index=False):
        tasks.append(
            {
                "image_key": row.image_key,
                "green_path": row.green_path,
                "cell_path": run_dir / "masks" / "cell" / f"{row.image_key}_cell_masks.tif",
                "nuclear_path": run_dir / "masks" / "nuclear" / f"{row.image_key}_nuclear_masks.tif",
                "background_sigma_px": args.background_sigma_px,
                "background_cell_dilation_px": args.background_cell_dilation_px,
                "cell_erosion_px": args.cell_erosion_px,
                "nuclear_dilation_px": args.nuclear_dilation_px,
            }
        )

    if args.workers == 1:
        feature_frames = [extract_image_features(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            feature_frames = list(executor.map(extract_image_features, tasks))
    features = pd.concat(feature_frames, ignore_index=True)
    if features.duplicated(["image_key", "object_id"]).any():
        raise RuntimeError("Extracted feature keys are not unique")

    selected_objects = objects.loc[
        objects["image_key"].isin(selected_manifest["image_key"])
    ].copy()
    calls = selected_objects.merge(
        features,
        on=["image_key", "object_id"],
        how="left",
        validate="one_to_one",
    )
    calls["eligible"] = (
        calls["predicted_label_name"].eq("alive")
        & calls["cell_area_pass"].eq(1)
        & calls["cytoplasm_pixel_count"].ge(args.min_cytoplasm_px)
        & np.isfinite(calls["green_cytoplasmic_median_bg_corrected"])
    )
    calls["ploidy_score"] = np.arcsinh(
        calls["green_cytoplasmic_median_bg_corrected"] / args.asinh_cofactor
    )

    calibration = calls.loc[calls["eligible"] & calls["hours"].eq(calibration_hour)]
    fit = fit_grouped_gaussian_mixture(
        calibration["ploidy_score"].to_numpy(),
        calibration["well"].astype(str).to_numpy(),
        seed=args.seed,
    )
    means = np.asarray(fit["means"], dtype=float)
    variances = np.asarray(fit["variances"], dtype=float)
    boundary = central_likelihood_boundary(
        means[0], variances[0], means[1], variances[1]
    )
    separation = float((means[1] - means[0]) / math.sqrt(np.mean(variances)))

    calls["log_likelihood_low"] = normal_log_density(
        calls["ploidy_score"].to_numpy(), means[0], variances[0]
    )
    calls["log_likelihood_high"] = normal_log_density(
        calls["ploidy_score"].to_numpy(), means[1], variances[1]
    )
    likelihood_difference = calls["log_likelihood_high"] - calls["log_likelihood_low"]
    calls["prob_4N"] = 1 / (1 + np.exp(-np.clip(likelihood_difference, -700, 700)))
    calls["ploidy_call"] = "not_eligible"
    calls.loc[calls["eligible"] & (calls["ploidy_score"] < boundary), "ploidy_call"] = "low"
    calls.loc[calls["eligible"] & (calls["ploidy_score"] >= boundary), "ploidy_call"] = "high"

    group_names = np.asarray(fit["group_names"]).astype(str)
    weights = np.asarray(fit["weights"], dtype=float)
    weight_table = pd.DataFrame(
        {
            "well": group_names,
            "baseline_weight_low": weights[:, 0],
            "baseline_weight_high": weights[:, 1],
        }
    ).sort_values("well")
    pooled_high = float(
        np.average(
            weight_table["baseline_weight_high"],
            weights=calibration.groupby("well").size().reindex(weight_table["well"]).to_numpy(),
        )
    )
    mixture_table = pd.DataFrame(
        {
            "component": ["low", "high"],
            "mean_ploidy_score": means,
            "variance_ploidy_score": variances,
            "pooled_baseline_weight": [1 - pooled_high, pooled_high],
            "standardized_mean_separation": separation,
            "equal_likelihood_boundary": boundary,
            "calibration_hour": calibration_hour,
            "n_calibration_objects": len(calibration),
            "em_log_likelihood": float(fit["log_likelihood"]),
            "em_iterations": int(fit["iterations"]),
            "em_converged": bool(fit["converged"]),
        }
    )

    export_columns = [
        "image_key", "well", "position", "hours", "G0_mM", "object_id",
        "predicted_label_name", "cell_area_px", "nuclear_area_px", "cell_area_pass",
        "cytoplasm_pixel_count", "green_cytoplasmic_median_raw",
        "green_local_background_median", "green_cytoplasmic_median_bg_corrected",
        "green_background_support_median", "ploidy_score", "prob_4N", "ploidy_call",
    ]
    calls_path = output_dir / "sum159_competition_object_ploidy_calls.csv.gz"
    calls.loc[calls["eligible"], export_columns].to_csv(
        calls_path,
        index=False,
        compression="gzip",
    )
    mixture_table.to_csv(output_dir / "sum159_competition_ploidy_mixture.csv", index=False)
    weight_table.to_csv(output_dir / "sum159_competition_baseline_well_weights.csv", index=False)
    call_summary = (
        calls.loc[calls["eligible"]]
        .groupby(["hours", "G0_mM", "ploidy_call"], dropna=False)
        .size()
        .rename("n_objects")
        .reset_index()
    )
    call_summary.to_csv(output_dir / "sum159_competition_call_summary.csv", index=False)

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command_parameters": vars(args),
        "run_dir": str(run_dir),
        "manifest_path": str(manifest_path),
        "objects_path": str(objects_path),
        "processed_image_count": int(len(selected_manifest)),
        "calibration_image_count": int((selected_manifest["hours"] == calibration_hour).sum()),
        "qc_image_count": int(len(examples)),
        "calibration_object_count": int(len(calibration)),
        "eligible_export_count": int(calls["eligible"].sum()),
        "git_commit": git_commit(),
        "python": platform.python_version(),
        "classification_rule": (
            "hard low/high call at the central equal-component-likelihood boundary; "
            "fitted baseline mixture weights do not enter the call boundary"
        ),
    }
    (output_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    readme = f"""# SUM-159 cytoplasmic GFP ploidy classification

- GFP background: normalized Gaussian estimate (sigma {args.background_sigma_px:g}px) from non-cell pixels after a {args.background_cell_dilation_px}px cell-mask dilation.
- Cytoplasm: cell mask with a {args.cell_erosion_px}px boundary erosion, excluding the nuclear mask after {args.nuclear_dilation_px}px dilation.
- Feature: asinh(background-corrected cytoplasmic GFP median / {args.asinh_cofactor:g}).
- Calibration: two Gaussian components with shared component means/variances across wells and freely estimated baseline mixture weights for every well.
- Hard call: central equal-component-likelihood boundary ({boundary:.8g}); fitted mixture weights do not move the boundary.
- Baseline pooled high-component weight: {pooled_high:.6f}; standardized component separation: {separation:.4f}.
- Scope: {'all manifest images' if args.all_images else 'all baseline images plus the selected QC images'} ({len(selected_manifest)} images).
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"Processed {len(selected_manifest)} images")
    print(f"Fitted {len(calibration)} eligible baseline cells")
    print(f"Boundary={boundary:.8g}; separation={separation:.4f}; pooled high weight={pooled_high:.6f}")
    print(f"Wrote calls to {calls_path}")


if __name__ == "__main__":
    main()
