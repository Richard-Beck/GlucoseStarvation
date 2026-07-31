#!/usr/bin/env python3
"""Populate immutable v03 assembly inputs from approved upstream packages."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]

FIGURE_ROOT = PROJECT / "agent-dev/manuscript_integration/20260729_v03_figure_set_integration"
RESULTS_ROOT = PROJECT / "agent-dev/manuscript_results/20260731_v03_results_revision"
METHODS_ROOT = PROJECT / "agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction"
CLAIM_ROOT = PROJECT / "agent-dev/manuscript_claim_graph/20260730_v03_clean_figure_level_graph"
AID_ROOT = PROJECT / "agent-dev/manuscript_aid/20260731_v03_aid_serving"

FIGURE_IDS = [f"F{i}" for i in range(1, 6)] + [f"S{i}" for i in range(1, 19)]
RESULT_SIDECARS = [
    "results_measurement_foundation.md",
    "results_direct_features.md",
    "results_mechanistic_models.md",
    "results_posterior_size_context.md",
    "results_selection_simulations.md",
]

FIGURE_PACKAGES = {
    **{fid: "resegmentation_core" for fid in ["F1", "F2", "S1", "S2", "S3", "S4"]},
    **{fid: "mechanistic_diagnostics" for fid in ["F3", "S5", "S6", "S7", "S8", "S9"]},
    **{fid: "posterior_size_context" for fid in ["F4", "S10", "S11"]},
    "F5": "posterior_strategy_F5",
    "S12": "posterior_strategy_context",
    "S13": "morphology_metrics",
    "S14": "morphology_metrics",
    **{fid: "sum159_label_swap" for fid in ["S15", "S16", "S17", "S18"]},
}

PACKAGE_SCRIPTS = {
    "resegmentation_core": "scripts/polish_figures.R",
    "mechanistic_diagnostics": "polish_figures.R",
    "posterior_size_context": "polish_figures.R",
    "posterior_strategy_F5": "scripts/polish_figures.R",
    "posterior_strategy_context": "scripts/polish_figures.R",
    "morphology_metrics": "scripts/polish_figures.R",
    "sum159_label_swap": "scripts/polish_figures.R",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(
        p for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts
    ):
        rel = str(path.relative_to(root)).encode("utf-8")
        digest.update(len(rel).to_bytes(8, "big"))
        digest.update(rel)
        data_hash = bytes.fromhex(sha256(path))
        digest.update(data_hash)
    return digest.hexdigest()


def copy_file(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination)


def project_rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def main() -> None:
    guarded = [ROOT / name for name in ("source", "assets", "evidence", "traceability")]
    if any(path.exists() for path in guarded):
        raise SystemExit("refusing to overwrite an existing assembly population")

    # Manuscript-facing sources.
    copy_file(FIGURE_ROOT / "figure_set_manifest.csv", ROOT / "source/figure_set_manifest.csv")
    copy_file(
        FIGURE_ROOT / "integrated_figure_legends.md",
        ROOT / "source/integrated_figure_legends.md",
    )
    copy_file(
        RESULTS_ROOT / "combined_results_preview.md",
        ROOT / "source/manuscript_sections/results.md",
    )
    for filename in RESULT_SIDECARS:
        copy_file(
            RESULTS_ROOT / "manuscript_sections/results" / filename,
            ROOT / "source/manuscript_sections/results" / filename,
        )
    copy_file(METHODS_ROOT / "methods_text.md", ROOT / "source/manuscript_sections/methods.md")
    for section in ("abstract", "introduction", "discussion", "references"):
        copy_file(
            AID_ROOT / "served" / f"{section}.md",
            ROOT / "source/manuscript_sections" / f"{section}.md",
        )

    # Final assets and complete assembly-local polishing packages.
    for figure_id in FIGURE_IDS:
        copy_file(
            FIGURE_ROOT / "final_images" / f"{figure_id}.png",
            ROOT / "assets/figures" / f"{figure_id}.png",
        )
    copy_tree(
        FIGURE_ROOT / "final_figure_scripts/package_scripts",
        ROOT / "rebuild/figures/package_scripts",
    )

    # Concise evidence needed to audit manuscript strength and Methods provenance.
    for filename in (
        "claim_graph_integrated.json",
        "claim_reconciliation.md",
        "claim_graph_integration_report.md",
        "claim_graph_validation_report.txt",
    ):
        copy_file(CLAIM_ROOT / filename, ROOT / "evidence/claim_graph" / filename)
    for filename in (
        "target_figure_set.tsv",
        "locked_provenance_table.md",
        "provenance_lock_verification.md",
        "final_methods_draft_audit.md",
        "methods_handoff.md",
    ):
        copy_file(METHODS_ROOT / filename, ROOT / "evidence/methods" / filename)
    copy_file(
        RESULTS_ROOT / "validation_report.md",
        ROOT / "evidence/results/validation_report.md",
    )
    copy_file(
        RESULTS_ROOT / "recycling_audit.tsv",
        ROOT / "evidence/results/recycling_audit.tsv",
    )
    copy_file(
        FIGURE_ROOT / "legend_validation_report.md",
        ROOT / "evidence/legends/legend_validation_report.md",
    )
    copy_tree(AID_ROOT, ROOT / "evidence/abstract_introduction_discussion")

    # Figure-facing integration traceability.
    for filename in (
        "figure_rebuild_manifest.tsv",
        "figure_byte_identity_report.tsv",
        "figure_set_integration_report.md",
        "figure_set_validation_report.json",
        "whole_set_visual_qc.md",
        "whole_set_visual_qc.tsv",
        "omitted_packages.md",
    ):
        copy_file(FIGURE_ROOT / filename, ROOT / "traceability/figures" / filename)

    with (ROOT / "source/figure_set_manifest.csv").open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))
    if len(manifest_rows) != 79:
        raise RuntimeError(f"expected 79 panel rows, found {len(manifest_rows)}")

    # One assembly rebuild record per whole figure.
    rebuild_rows = []
    asset_rows = []
    for figure_id in FIGURE_IDS:
        rows = [row for row in manifest_rows if row["current_figure_name"] == figure_id]
        if not rows:
            raise RuntimeError(f"missing manifest rows for {figure_id}")
        expected_hashes = {row["final_image_sha256"] for row in rows}
        if len(expected_hashes) != 1:
            raise RuntimeError(f"non-unique expected hash for {figure_id}: {expected_hashes}")
        expected_hash = expected_hashes.pop()
        asset = ROOT / "assets/figures" / f"{figure_id}.png"
        if sha256(asset) != expected_hash:
            raise RuntimeError(f"asset hash mismatch for {figure_id}")
        package = FIGURE_PACKAGES[figure_id]
        script_rel = (
            f"rebuild/figures/package_scripts/{package}/{PACKAGE_SCRIPTS[package]}"
        )
        direct_inputs = sorted(
            {
                value
                for row in rows
                for value in row["data_input_paths"].split(";")
                if value and value != "NA"
            }
        )
        if figure_id == "F3":
            direct_inputs.append(
                "figures/user-approved-raster-figures/model_family_schematic_placeholder.png"
            )
            direct_inputs.sort()
        source_roots = sorted({row["source_root"] for row in rows})
        approved_raster = "yes" if figure_id == "F3" else "no"
        accepted_exception = (
            "Figure 3a uses the user-approved model-family schematic raster"
            if figure_id == "F3"
            else "none"
        )
        rebuild_rows.append(
            {
                "figure_id": figure_id,
                "asset_path": f"assets/figures/{figure_id}.png",
                "rebuild_output_path": f"OUTPUT_DIR/{figure_id}.png",
                "source_package": package,
                "polish_root": ";".join(source_roots),
                "polishing_script": script_rel,
                "rebuild_command": "rebuild/figures/run_all_figures.sh OUTPUT_DIR",
                "direct_inputs": ";".join(direct_inputs),
                "dependency_paths": (
                    f"rebuild/figures/package_scripts/{package};"
                    "rebuild/figures/normalize_rebuild_outputs.py;"
                    "rebuild/figures/run_all_figures.sh"
                ),
                "approved_raster": approved_raster,
                "expected_sha256": expected_hash,
                "accepted_exception": accepted_exception,
            }
        )
        asset_rows.append(
            {
                "asset_id": figure_id,
                "asset_path": f"assets/figures/{figure_id}.png",
                "manuscript_location": (
                    "main Results" if figure_id.startswith("F") else "supplemental figure index"
                ),
                "legend_label": (
                    f"Figure {figure_id[1:]}"
                    if figure_id.startswith("F")
                    else f"Figure {figure_id}"
                ),
                "rebuild_manifest_figure_id": figure_id,
                "sha256": expected_hash,
            }
        )

    rebuild_path = ROOT / "rebuild/figures/figure_rebuild_manifest.tsv"
    rebuild_path.parent.mkdir(parents=True, exist_ok=True)
    with rebuild_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rebuild_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rebuild_rows)
    inventory_path = ROOT / "assets/asset_inventory.tsv"
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    with inventory_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asset_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(asset_rows)

    # Component-level upstream register. Directory entries use a deterministic
    # tree hash over relative file names and file hashes.
    source_components = [
        ("approved figure integration", FIGURE_ROOT, "copied final images, legends, manifests, and rebuild packages"),
        ("approved Results", RESULTS_ROOT, "copied combined text, five sidecars, and audits"),
        ("approved Methods", METHODS_ROOT, "copied canonical prose and provenance evidence"),
        ("accepted claim graph", CLAIM_ROOT, "copied graph and reconciliation evidence"),
        ("served A/I/D", AID_ROOT, "copied served prose and ordered bundle evidence"),
    ]
    register_rows = []
    for role, path, rationale in source_components:
        register_rows.append(
            {
                "role": role,
                "upstream_path": project_rel(path),
                "object_type": "directory",
                "sha256": tree_sha256(path),
                "assembly_disposition": "selected manuscript-facing files copied; large upstream work referenced",
                "rationale": rationale,
            }
        )
    register_path = ROOT / "traceability/upstream_input_register.tsv"
    register_path.parent.mkdir(parents=True, exist_ok=True)
    with register_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(register_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(register_rows)

    preparation = {
        "schema_version": 1,
        "status": "pass",
        "figure_count": len(FIGURE_IDS),
        "panel_endpoint_count": len(manifest_rows),
        "result_sidecar_count": len(RESULT_SIDECARS),
        "assembly_root": str(ROOT.relative_to(PROJECT)),
        "source_roots": [project_rel(path) for _, path, _ in source_components],
    }
    (ROOT / "traceability/preparation_status.json").write_text(
        json.dumps(preparation, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
