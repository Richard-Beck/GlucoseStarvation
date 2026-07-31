#!/usr/bin/env python3
"""Validate v03 source identity, rendering, and assembly-local figure replay."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
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
        digest.update(bytes.fromhex(sha256(path)))
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-dir", type=Path, required=True)
    args = parser.parse_args()
    rebuild_dir = args.rebuild_dir.resolve()

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "pass": bool(condition), "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    identity_specs: list[tuple[str, Path, Path]] = [
        (
            "figure manifest",
            ROOT / "source/figure_set_manifest.csv",
            FIGURE_ROOT / "figure_set_manifest.csv",
        ),
        (
            "figure legends",
            ROOT / "source/integrated_figure_legends.md",
            FIGURE_ROOT / "integrated_figure_legends.md",
        ),
        (
            "combined Results",
            ROOT / "source/manuscript_sections/results.md",
            RESULTS_ROOT / "combined_results_preview.md",
        ),
        (
            "Methods",
            ROOT / "source/manuscript_sections/methods.md",
            METHODS_ROOT / "methods_text.md",
        ),
        (
            "claim graph",
            ROOT / "evidence/claim_graph/claim_graph_integrated.json",
            CLAIM_ROOT / "claim_graph_integrated.json",
        ),
        (
            "Methods provenance",
            ROOT / "evidence/methods/locked_provenance_table.md",
            METHODS_ROOT / "locked_provenance_table.md",
        ),
    ]
    for filename in RESULT_SIDECARS:
        identity_specs.append(
            (
                f"Results sidecar {filename}",
                ROOT / "source/manuscript_sections/results" / filename,
                RESULTS_ROOT / "manuscript_sections/results" / filename,
            )
        )
    for section in ("abstract", "introduction", "discussion", "references"):
        identity_specs.append(
            (
                f"served {section}",
                ROOT / "source/manuscript_sections" / f"{section}.md",
                AID_ROOT / "served" / f"{section}.md",
            )
        )
    for figure_id in FIGURE_IDS:
        identity_specs.append(
            (
                f"figure {figure_id}",
                ROOT / "assets/figures" / f"{figure_id}.png",
                FIGURE_ROOT / "final_images" / f"{figure_id}.png",
            )
        )

    identity_rows = []
    for role, local, upstream in identity_specs:
        check(f"input_exists:{role}", local.is_file() and upstream.is_file(), str(upstream))
        local_hash = sha256(local)
        upstream_hash = sha256(upstream)
        check(f"input_identity:{role}", local_hash == upstream_hash, local_hash)
        identity_rows.append(
            {
                "role": role,
                "assembly_path": str(local.relative_to(PROJECT)),
                "upstream_path": str(upstream.relative_to(PROJECT)),
                "sha256": local_hash,
                "identity": "match",
            }
        )

    local_packages = ROOT / "rebuild/figures/package_scripts"
    upstream_packages = FIGURE_ROOT / "final_figure_scripts/package_scripts"
    local_package_hash = tree_sha256(local_packages)
    upstream_package_hash = tree_sha256(upstream_packages)
    check(
        "figure_package_tree_identity",
        local_package_hash == upstream_package_hash,
        local_package_hash,
    )

    with (ROOT / "traceability/copied_input_identity.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(identity_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(identity_rows)

    with (ROOT / "source/figure_set_manifest.csv").open(newline="") as handle:
        panel_rows = list(csv.DictReader(handle))
    observed_ids = []
    for row in panel_rows:
        if row["current_figure_name"] not in observed_ids:
            observed_ids.append(row["current_figure_name"])
    check("panel_endpoint_count", len(panel_rows) == 79, f"observed={len(panel_rows)}")
    check("figure_order", observed_ids == FIGURE_IDS, repr(observed_ids))
    check("figure_S19_absent", "S19" not in observed_ids, repr(observed_ids))

    with (ROOT / "rebuild/figures/figure_rebuild_manifest.tsv").open(newline="") as handle:
        rebuild_records = list(csv.DictReader(handle, delimiter="\t"))
    check("figure_rebuild_rows", len(rebuild_records) == 23, f"observed={len(rebuild_records)}")
    check(
        "figure_rebuild_ids",
        [row["figure_id"] for row in rebuild_records] == FIGURE_IDS,
        repr([row["figure_id"] for row in rebuild_records]),
    )

    dependency_rows: list[dict[str, str]] = []
    figure_validation_rows: list[dict[str, str]] = []
    for record in rebuild_records:
        figure_id = record["figure_id"]
        asset = ROOT / record["asset_path"]
        rebuilt = rebuild_dir / f"{figure_id}.png"
        expected_hash = record["expected_sha256"]
        check(f"asset_exists:{figure_id}", asset.is_file(), str(asset))
        check(f"rebuilt_exists:{figure_id}", rebuilt.is_file(), str(rebuilt))
        asset_hash = sha256(asset)
        rebuilt_hash = sha256(rebuilt)
        check(f"asset_hash:{figure_id}", asset_hash == expected_hash, asset_hash)
        check(f"rebuild_hash:{figure_id}", rebuilt_hash == expected_hash, rebuilt_hash)
        figure_validation_rows.append(
            {
                "figure_id": figure_id,
                "asset_path": record["asset_path"],
                "rebuild_path": str(rebuilt),
                "expected_sha256": expected_hash,
                "asset_sha256": asset_hash,
                "rebuild_sha256": rebuilt_hash,
                "status": "PASS",
            }
        )
        script = ROOT / record["polishing_script"]
        check(f"polishing_script:{figure_id}", script.is_file(), str(script))
        for kind, field, base in (
            ("direct_input", "direct_inputs", PROJECT),
            ("dependency", "dependency_paths", ROOT),
        ):
            for path_text in filter(None, record[field].split(";")):
                path = base / path_text
                exists = path.exists()
                dependency_rows.append(
                    {
                        "figure_id": figure_id,
                        "kind": kind,
                        "path": path_text,
                        "exists": "yes" if exists else "no",
                        "status": "PASS" if exists else "FAIL",
                    }
                )
                check(f"{kind}:{figure_id}:{path_text}", exists, str(path))

    raster_records = [row for row in rebuild_records if row["approved_raster"] == "yes"]
    check(
        "approved_raster_scope",
        len(raster_records) == 1 and raster_records[0]["figure_id"] == "F3",
        repr([(row["figure_id"], row["accepted_exception"]) for row in raster_records]),
    )

    with (ROOT / "validation/figure_rebuild_validation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(figure_validation_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(figure_validation_rows)
    with (ROOT / "validation/figure_rebuild_dependency_validation.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dependency_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(dependency_rows)

    legends = (ROOT / "source/integrated_figure_legends.md").read_text(encoding="utf-8")
    legend_labels = re.findall(r"^## Figure (S?\d+)\.", legends, flags=re.M)
    expected_labels = [str(i) for i in range(1, 6)] + [f"S{i}" for i in range(1, 19)]
    check("legend_top_heading", legends.startswith("# Figure Legends\n"), "heading")
    check("legend_blocks", legend_labels == expected_labels, repr(legend_labels))
    check("legend_S19_absent", "Figure S19" not in legends, "S19 absent")
    owner_legend_report = (ROOT / "evidence/legends/legend_validation_report.md").read_text()
    check("owner_legend_validation", "PASS" in owner_legend_report, "owner report passes")

    sidecars = sorted((ROOT / "source/manuscript_sections/results").glob("*.md"))
    check("result_sidecar_count", len(sidecars) == 5, f"observed={len(sidecars)}")
    check(
        "result_sidecar_contract",
        all("## Results Text" in path.read_text() and "## Revision Notes" in path.read_text() for path in sidecars),
        "five manuscript/audit sidecars",
    )
    check(
        "results_approved_validation",
        "Recommendation" in (ROOT / "evidence/results/validation_report.md").read_text(),
        "Results validation retained",
    )

    graph = json.loads((ROOT / "evidence/claim_graph/claim_graph_integrated.json").read_text())
    fixed_claims = [claim for claim in graph["claims"] if claim.get("user_fixed") is True]
    check("user_fixed_claim_count", len(fixed_claims) == 8, f"observed={len(fixed_claims)}")
    manuscript_source = "\n".join(
        (ROOT / "source/manuscript_sections" / filename).read_text().lower()
        for filename in (
            "abstract.md",
            "introduction.md",
            "results.md",
            "discussion.md",
            "methods.md",
        )
    )
    claim_patterns = {
        "C0": r"ploidy (?:modulates|changes).{0,80}(?:glucose|starvation)",
        "C1": r"glucose drawdown per (?:unit )?live-cell auc|greater glucose depletion after accounting",
        "C2": r"lower fitted per-cell yield|fewer additional cells under low glucose",
        "C3": r"more robust to starvation|maintained more viable cells",
        "C4": r"starvation selects for high ploidy|depletion favored high ploidy",
        "C7": r"cell size generally increased with ploidy|high-ploidy states had larger cells",
        "C11": r"abundance selects for low ploidy|replenishment favored low ploidy",
        "C12": r"hypothesis that high ploidy increases volumetric yield",
    }
    for claim_id, pattern in claim_patterns.items():
        check(
            f"user_fixed_claim_prose:{claim_id}",
            bool(re.search(pattern, manuscript_source, flags=re.S)),
            pattern,
        )

    references = (ROOT / "source/manuscript_sections/references.md").read_text()
    reference_numbers = [int(value) for value in re.findall(r"(?m)^(\d+)\.\s", references)]
    check(
        "reference_sequence",
        reference_numbers == list(range(1, max(reference_numbers) + 1)),
        f"count={len(reference_numbers)}",
    )
    aid_text = "\n".join(
        (ROOT / "source/manuscript_sections" / f"{name}.md").read_text()
        for name in ("abstract", "introduction", "discussion")
    )
    cited_numbers = [
        int(value)
        for group in re.findall(r"\[([0-9,–\- ]+)\]", aid_text)
        for value in re.findall(r"\d+", group)
    ]
    check(
        "citation_reference_bounds",
        cited_numbers and max(cited_numbers) <= max(reference_numbers),
        f"max_citation={max(cited_numbers)} references={max(reference_numbers)}",
    )

    html_path = ROOT / "draft/manuscript_draft.html"
    html_text = html_path.read_text(encoding="utf-8")
    image_srcs = re.findall(r'<img src="([^"]+)"', html_text)
    check("html_image_count", len(image_srcs) == 23, f"observed={len(image_srcs)}")
    check(
        "html_images_embedded",
        all(src.startswith("data:image/png;base64,") for src in image_srcs),
        "all images are PNG data URIs",
    )
    anchors = re.findall(r'id="([^"]+)"', html_text)
    check("unique_html_anchors", len(anchors) == len(set(anchors)), f"anchors={len(anchors)}")
    href_targets = re.findall(r'href="#([^"]+)"', html_text)
    unresolved_hrefs = sorted(set(href_targets) - set(anchors))
    check("internal_html_links", not unresolved_hrefs, repr(unresolved_hrefs))
    figure_blocks = re.findall(
        r'<figure[^>]*id="fig-([^"]+)"[^>]*>(.*?)</figure>', html_text, flags=re.S
    )
    block_ids = [figure_id.upper() for figure_id, _ in figure_blocks]
    check("html_figure_ids", block_ids == FIGURE_IDS, repr(block_ids))
    for (figure_id_lower, block), src in zip(figure_blocks, image_srcs):
        figure_id = figure_id_lower.upper()
        block_sources = re.findall(r'<img src="([^"]+)"', block)
        check(f"html_one_image:{figure_id}", len(block_sources) == 1, repr(block_sources))
        decoded = base64.b64decode(
            block_sources[0].removeprefix("data:image/png;base64,"), validate=True
        )
        decoded_hash = hashlib.sha256(decoded).hexdigest()
        check(
            f"html_image_identity:{figure_id}",
            decoded_hash == sha256(ROOT / "assets/figures" / f"{figure_id}.png"),
            decoded_hash,
        )
    check(
        "portable_html",
        not re.search(r'<(?:img|script)[^>]+src="(?:file:|https?:|/)', html_text, flags=re.I)
        and not re.search(r"<link[^>]+href=", html_text, flags=re.I),
        "no external image, script, or stylesheet dependency",
    )
    check("references_rendered", '<section id="references">' in html_text, "references section")
    check("audit_trail_rendered", '<section id="audit-trail">' in html_text, "audit section")
    journal_body = html_text.split('<section id="audit-trail">', 1)[0]
    check("results_notes_excluded", "Revision Notes" not in journal_body, "audit notes excluded")
    html_without_images = re.sub(
        r"data:image/png;base64,[A-Za-z0-9+/=]+", "[embedded PNG]", journal_body
    )
    leaked = [
        term
        for term in ("/share/", "agent-dev/", "manuscripts/2026", "polishing package")
        if term in html_without_images
    ]
    check("journal_provenance_hygiene", not leaked, repr(leaked))
    visible_callouts = re.findall(r"\bFigures?\s+(S?\d+)", html_without_images)
    invalid_callouts = sorted(set(visible_callouts) - set(expected_labels))
    check("figure_callout_targets", not invalid_callouts, repr(invalid_callouts))
    check("html_S19_absent", "Figure S19" not in html_without_images, "S19 absent")

    status_payload = {
        "schema_version": 1,
        "status": "PASS",
        "checks": checks,
        "check_count": len(checks),
        "figure_count": 23,
        "panel_endpoint_count": 79,
        "rendered_html_sha256": sha256(html_path),
        "figure_package_tree_sha256": local_package_hash,
        "figure_rebuild_directory": str(rebuild_dir),
    }
    (ROOT / "validation/draft_validation.json").write_text(
        json.dumps(status_payload, indent=2) + "\n", encoding="utf-8"
    )
    report = [
        "# Draft and component validation",
        "",
        "Status: PASS",
        "",
        f"Checks passed: {len(checks)}",
        "",
        "- 23 approved figures and all 79 panel endpoints authenticate against the integration root.",
        "- All 23 figures regenerated from assembly-local polishing packages and matched approved checksums.",
        "- The HTML contains exactly 23 embedded PNGs, unique anchors, resolved internal links, and no external asset dependencies.",
        "- The 23 legends map exactly to the rendered figures; Figure S19 is absent.",
        "- Five Results sidecars, canonical Methods, served A/I/D and References, and all eight user-fixed claims are represented.",
        "- Journal-facing content excludes Results revision notes and local provenance paths.",
        "",
        f"Rendered HTML SHA-256: `{sha256(html_path)}`",
    ]
    (ROOT / "validation/draft_validation_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
