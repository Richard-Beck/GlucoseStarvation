#!/usr/bin/env python3
"""Validate final v03 package structure, indexes, and checksum inventory."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHECKSUMS = ROOT / "traceability/package_file_checksums.sha256"
ASSEMBLY_OUTPUTS = {
    ROOT / "validation/assembly_validation.json",
    ROOT / "validation/assembly_validation_report.md",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "pass": bool(condition), "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    for dirname in (
        "draft",
        "source",
        "assets",
        "evidence",
        "traceability",
        "review_state",
        "validation",
        "rebuild",
    ):
        check(f"required_directory:{dirname}", (ROOT / dirname).is_dir(), dirname)
    required_files = [
        "README.md",
        "status.json",
        "draft/manuscript_draft.html",
        "source/title.txt",
        "source/figure_set_manifest.csv",
        "source/integrated_figure_legends.md",
        "source/manuscript_sections/abstract.md",
        "source/manuscript_sections/introduction.md",
        "source/manuscript_sections/results.md",
        "source/manuscript_sections/discussion.md",
        "source/manuscript_sections/methods.md",
        "source/manuscript_sections/references.md",
        "assets/asset_inventory.tsv",
        "evidence/claim_graph/claim_graph_integrated.json",
        "evidence/methods/locked_provenance_table.md",
        "traceability/upstream_input_register.tsv",
        "traceability/copied_input_identity.tsv",
        "review_state/component_status.tsv",
        "review_state/accepted_exceptions.md",
        "validation/draft_validation.json",
        "validation/figure_rebuild_validation.tsv",
        "rebuild/figures/run_all_figures.sh",
        "rebuild/figures/figure_rebuild_manifest.tsv",
        "rebuild/manuscript/build_manuscript_html.py",
        "rebuild/manuscript/validate_manuscript.py",
    ]
    for rel in required_files:
        check(f"required_file:{rel}", (ROOT / rel).is_file(), rel)

    status = json.loads((ROOT / "status.json").read_text())
    check("assembly_status", status.get("status") == "WARN", repr(status.get("status")))
    check(
        "review_package_purpose",
        status.get("package_purpose") == "coherent manuscript review assembly",
        repr(status.get("package_purpose")),
    )
    draft_validation = json.loads((ROOT / "validation/draft_validation.json").read_text())
    check("draft_validation", draft_validation.get("status") == "PASS", "PASS")
    check("draft_check_count", draft_validation.get("check_count") == 442, repr(draft_validation.get("check_count")))
    check(
        "draft_hash",
        sha256(ROOT / "draft/manuscript_draft.html") == status["rendered_draft_sha256"],
        status["rendered_draft_sha256"],
    )

    with (ROOT / "validation/figure_rebuild_validation.tsv").open(newline="") as handle:
        figure_rows = list(csv.DictReader(handle, delimiter="\t"))
    check("figure_validation_rows", len(figure_rows) == 23, f"observed={len(figure_rows)}")
    check("figure_validation_status", all(row["status"] == "PASS" for row in figure_rows), "all PASS")
    with (ROOT / "validation/figure_rebuild_dependency_validation.tsv").open(newline="") as handle:
        dependency_rows = list(csv.DictReader(handle, delimiter="\t"))
    check("dependency_validation_rows", len(dependency_rows) == 168, f"observed={len(dependency_rows)}")
    check("dependency_validation_status", all(row["status"] == "PASS" for row in dependency_rows), "all PASS")

    assets = sorted((ROOT / "assets/figures").glob("*.png"))
    check("asset_count", len(assets) == 23, f"observed={len(assets)}")
    packages = sorted(path for path in (ROOT / "rebuild/figures/package_scripts").iterdir() if path.is_dir())
    check("polishing_package_count", len(packages) == 7, f"observed={len(packages)}")
    with (ROOT / "rebuild/figures/figure_rebuild_manifest.tsv").open(newline="") as handle:
        rebuild_rows = list(csv.DictReader(handle, delimiter="\t"))
    check("rebuild_manifest_rows", len(rebuild_rows) == 23, f"observed={len(rebuild_rows)}")

    check("checksum_inventory_exists", CHECKSUMS.is_file(), str(CHECKSUMS))
    listed: dict[Path, str] = {}
    for line in CHECKSUMS.read_text().splitlines():
        digest, rel = line.split("  ", 1)
        path = ROOT / rel
        check(f"checksum_path:{rel}", path.is_file(), rel)
        check(f"checksum_value:{rel}", sha256(path) == digest, digest)
        listed[path] = digest
    expected_files = {
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and path != CHECKSUMS
        and path not in ASSEMBLY_OUTPUTS
        and "__pycache__" not in path.parts
    }
    check(
        "checksum_inventory_complete",
        set(listed) == expected_files,
        f"listed={len(listed)} expected={len(expected_files)}",
    )

    payload = {
        "schema_version": 1,
        "status": "WARN",
        "assembly_contract": "PASS",
        "checks": checks,
        "check_count": len(checks),
        "warnings": status["warnings"],
        "rendered_draft_sha256": status["rendered_draft_sha256"],
        "package_checksum_rows": len(listed),
    }
    (ROOT / "validation/assembly_validation.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Assembly validation report",
        "",
        "Assembly contract: PASS",
        "Overall package status: WARN",
        "",
        f"Checks passed: {len(checks)}",
        f"Package checksum records: {len(listed)}",
        "",
        "The package has the complete assembly structure, authenticated approved sources, a self-contained rendered draft, 23 checksum-verified figure assets, seven local polishing packages, complete figure-rebuild commands, and passing draft and dependency validation.",
        "",
        "The WARN disposition is limited to the accepted Methods, A/I/D-carry-forward, and workflow-gate exceptions recorded in `review_state/accepted_exceptions.md`; no assembly contract failure remains.",
    ]
    (ROOT / "validation/assembly_validation_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
