#!/usr/bin/env python3
"""Normalize package-local outputs and require approved whole-figure hashes."""

from __future__ import annotations

import csv
import hashlib
import shutil
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ASSEMBLY_ROOT = SCRIPT_DIR.parents[1]
MANIFEST = SCRIPT_DIR / "figure_rebuild_manifest.tsv"

SPECS = [
    ("resegmentation_core", "figure_1.png", "F1"),
    ("resegmentation_core", "figure_2.png", "F2"),
    ("mechanistic_diagnostics", "figure_3.png", "F3"),
    ("posterior_size_context", "figure_4.png", "F4"),
    ("posterior_strategy_F5", "figure_5.png", "F5"),
    ("resegmentation_core", "figure_s1.png", "S1"),
    ("resegmentation_core", "figure_s2.png", "S2"),
    ("resegmentation_core", "figure_s3.png", "S3"),
    ("resegmentation_core", "figure_s4.png", "S4"),
    ("mechanistic_diagnostics", "figure_s5.png", "S5"),
    ("mechanistic_diagnostics", "figure_s6.png", "S6"),
    ("mechanistic_diagnostics", "figure_s7.png", "S7"),
    ("mechanistic_diagnostics", "figure_s8.png", "S8"),
    ("mechanistic_diagnostics", "figure_s9.png", "S9"),
    ("posterior_size_context", "figure_s10.png", "S10"),
    ("posterior_size_context", "figure_s13.png", "S11"),
    ("posterior_strategy_context", "posterior_strategy_context_support_median.png", "S12"),
    ("morphology_metrics", "candidate_supplement_a.png", "S13"),
    ("morphology_metrics", "candidate_supplement_b.png", "S14"),
    ("sum159_label_swap", "figure_1_all_timecourses.png", "S15"),
    ("sum159_label_swap", "figure_2_confluence_robustness.png", "S16"),
    ("sum159_label_swap", "figure_3_focused_distributions_and_same_2n.png", "S17"),
    ("sum159_label_swap", "figure_4_cytoplasmic_signal_and_multimodal_fields.png", "S18"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: normalize_rebuild_outputs.py PACKAGE_ROOT OUTPUT_DIR")
    package_root = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output_dir}")

    with MANIFEST.open(newline="") as handle:
        records = {row["figure_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    expected_ids = {figure_id for _, _, figure_id in SPECS}
    if set(records) != expected_ids:
        raise RuntimeError(
            f"rebuild manifest figure IDs differ: {sorted(set(records) ^ expected_ids)}"
        )

    for package, filename, figure_id in SPECS:
        source = package_root / package / "final_images" / filename
        if not source.is_file():
            raise FileNotFoundError(source)
        observed_hash = sha256(source)
        expected_hash = records[figure_id]["expected_sha256"]
        if observed_hash != expected_hash:
            raise RuntimeError(
                f"unverified rebuild for {figure_id}: {observed_hash} != {expected_hash}"
            )
        destination = output_dir / f"{figure_id}.png"
        shutil.copy2(source, destination)
        if sha256(destination) != expected_hash:
            raise RuntimeError(f"normalization changed bytes for {figure_id}")


if __name__ == "__main__":
    main()
