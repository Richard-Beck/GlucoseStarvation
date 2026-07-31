#!/usr/bin/env python3
"""Render the manuscript claim graph as a Graphviz PNG.

The plot is claim-only by default. Evidence, method, and assumption nodes are
kept in the JSON for provenance but omitted from the visual graph so the
manuscript-level claim structure remains readable.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path


DEFAULT_GRAPH = Path(__file__).with_name("manuscript_claim_evidence_graph_v3_prototype.json")
DEFAULT_OUTPUT = Path(__file__).with_name("manuscript_claim_graph_v3_prototype.png")
DEFAULT_DOT = Path(__file__).with_name("manuscript_claim_graph_v3_prototype.dot")

TIER_ORDER = {
    "tier_0": 0,
    "tier_1": 1,
    "central": 0,
    "major_manuscript": 1,
    "supporting_result": 2,
    "boundary_or_undermining": 3,
}

TIER_PALETTE = {
    "tier_0": {"fill": "#5B2A86", "font": "#FFFFFF", "border": "#32104C"},
    "tier_1": {"fill": "#D7E9FF", "font": "#102A43", "border": "#2F6690"},
    "central": {"fill": "#5B2A86", "font": "#FFFFFF", "border": "#32104C"},
    "major_manuscript": {"fill": "#D7E9FF", "font": "#102A43", "border": "#2F6690"},
    "supporting_result": {"fill": "#E2F3E6", "font": "#173B23", "border": "#3A7D44"},
    "boundary_or_undermining": {"fill": "#FFF0D9", "font": "#402A0A", "border": "#C47F18"},
}

EDGE_STYLES = {
    "supports": {
        "color": "#2E7D32",
        "style": "solid",
        "penwidth": "2.1",
        "arrowhead": "normal",
        "label": "supports",
    },
    "qualifies": {
        "color": "#B36B00",
        "style": "dashed",
        "penwidth": "1.8",
        "arrowhead": "vee",
        "label": "qualifies",
    },
    "undermines": {
        "color": "#B00020",
        "style": "dotted",
        "penwidth": "2.2",
        "arrowhead": "tee",
        "label": "undermines",
    },
}


def load_graph(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def wrap_text(text: str, width: int) -> str:
    lines = textwrap.wrap(
        text,
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )
    return "\n".join(lines)


def dot_quote(value: str) -> str:
    return json.dumps(value)


def node_label(claim: dict, wrap_width: int) -> str:
    claim_id = claim.get("id", "")
    text = wrap_text(claim.get("text", ""), wrap_width)
    meta_parts = [claim.get("tier", ""), claim.get("status", "")]
    if claim.get("user_fixed") is True:
        meta_parts.append("user_fixed")
    meta = " | ".join(meta_parts)
    return f"{claim_id}\n{text}\n[{meta}]"


def collect_edges(claims: list[dict], complete_edges: bool = False) -> list[tuple[str, str, str]]:
    claims_by_id = {claim["id"]: claim for claim in claims}
    claim_ids = set(claims_by_id)
    edges: set[tuple[str, str, str]] = set()

    for claim in claims:
        claim_id = claim["id"]

        for target in claim.get("supports", []):
            if target in claim_ids:
                edges.add((claim_id, target, "supports"))
        for source in claim.get("supported_by_claims", []):
            if source in claim_ids:
                edges.add((source, claim_id, "supports"))

        for target in claim.get("qualifies", []):
            if target in claim_ids:
                edges.add((claim_id, target, "qualifies"))
        for source in claim.get("qualified_by_claims", []):
            if source in claim_ids:
                edges.add((source, claim_id, "qualifies"))

        for target in claim.get("undermines", []):
            if target in claim_ids:
                edges.add((claim_id, target, "undermines"))

    if not complete_edges:
        edges = {
            edge
            for edge in edges
            if include_in_presentation_view(edge, claims_by_id)
        }

    return sorted(edges, key=lambda edge: (edge[2], edge[1], edge[0]))


def include_in_presentation_view(edge: tuple[str, str, str], claims_by_id: dict[str, dict]) -> bool:
    source, target, relation = edge
    source_tier = claims_by_id[source].get("tier")
    target_tier = claims_by_id[target].get("tier")

    if relation == "undermines":
        return True

    if relation == "supports":
        if target_tier in {"central", "tier_0"}:
            return source_tier in {"major_manuscript", "tier_1"}
        return True

    if relation == "qualifies":
        return target_tier in {"central", "major_manuscript", "tier_0", "tier_1"}

    return True


def render_dot(graph: dict, wrap_width: int, edge_labels: bool, complete_edges: bool, cluster_tiers: bool) -> str:
    claims = sorted(
        graph.get("claims", []),
        key=lambda claim: (TIER_ORDER.get(claim.get("tier", ""), 99), claim.get("id", "")),
    )
    claims_by_tier: dict[str, list[dict]] = {}
    for claim in claims:
        claims_by_tier.setdefault(claim.get("tier", "other"), []).append(claim)

    lines = [
        "digraph manuscript_claim_graph {",
        "  graph [",
        '    rankdir="BT",',
        '    bgcolor="white",',
        '    margin="0.2",',
        '    pad="0.3",',
        '    nodesep="0.38",',
        '    ranksep="0.72",',
        '    splines="spline",',
        '    overlap="false",',
        '    concentrate="false",',
        '    labelloc="t",',
        '    label="Tiered manuscript claim graph prototype\\nsolid green = supports; dashed amber = qualifies; dotted red tee = undermines"',
        "  ];",
        '  node [shape="box", style="rounded,filled", fontname="Helvetica", fontsize="10", margin="0.09,0.07"];',
        '  edge [fontname="Helvetica", fontsize="8", arrowsize="0.8"];',
        "",
    ]

    if cluster_tiers:
        for tier, tier_claims in sorted(claims_by_tier.items(), key=lambda item: TIER_ORDER.get(item[0], 99)):
            lines.append(f"  subgraph cluster_{tier.replace('-', '_')} {{")
            lines.append(f'    label={dot_quote(tier.replace("_", " ").title())};')
            lines.append('    color="#D9DEE8";')
            lines.append('    fontname="Helvetica";')
            lines.append('    fontsize="12";')
            lines.append('    style="rounded";')
            for claim in tier_claims:
                lines.append(f"    {render_node(claim, wrap_width)}")
            lines.append("  }")
            lines.append("")
    else:
        for claim in claims:
            lines.append(f"  {render_node(claim, wrap_width)}")
        lines.append("")

    for source, target, relation in collect_edges(claims, complete_edges=complete_edges):
        style = EDGE_STYLES[relation]
        attrs = {
            "color": style["color"],
            "fontcolor": style["color"],
            "style": style["style"],
            "penwidth": style["penwidth"],
            "arrowhead": style["arrowhead"],
        }
        if edge_labels:
            attrs["label"] = style["label"]
        attr_text = ", ".join(f"{key}={dot_quote(value)}" for key, value in attrs.items())
        lines.append(f"  {dot_quote(source)} -> {dot_quote(target)} [{attr_text}];")

    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def render_node(claim: dict, wrap_width: int) -> str:
    tier = claim.get("tier", "other")
    palette = TIER_PALETTE.get(tier, {"fill": "#EEEEEE", "font": "#111111", "border": "#777777"})
    attrs = {
        "label": node_label(claim, wrap_width),
        "fillcolor": palette["fill"],
        "fontcolor": palette["font"],
        "color": palette["border"],
        "penwidth": "2.2" if tier in {"central", "major_manuscript", "tier_0", "tier_1"} else "1.4",
    }
    attr_text = ", ".join(f"{key}={dot_quote(value)}" for key, value in attrs.items())
    return f"{dot_quote(claim['id'])} [{attr_text}];"


def run_graphviz(dot_path: Path, output_path: Path, engine: str) -> None:
    executable = shutil.which(engine)
    if executable is None:
        raise RuntimeError(f"Graphviz executable not found: {engine}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [executable, "-Tpng", str(dot_path), "-o", str(output_path)],
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH, help="Input claim graph JSON.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PNG path.")
    parser.add_argument("--dot-output", type=Path, default=DEFAULT_DOT, help="Output DOT path.")
    parser.add_argument("--engine", default="dot", help="Graphviz layout executable to use.")
    parser.add_argument("--wrap-width", type=int, default=34, help="Approximate text wrap width for claim boxes.")
    parser.add_argument("--edge-labels", action="store_true", help="Print relation labels on every edge.")
    parser.add_argument("--complete-edges", action="store_true", help="Render every claim-claim edge from the JSON.")
    parser.add_argument("--cluster-tiers", action="store_true", help="Group node tiers into Graphviz clusters.")
    args = parser.parse_args()

    try:
        graph = load_graph(args.graph)
        dot_text = render_dot(graph, args.wrap_width, args.edge_labels, args.complete_edges, args.cluster_tiers)
        args.dot_output.write_text(dot_text, encoding="utf-8")
        run_graphviz(args.dot_output, args.output, args.engine)
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: Graphviz failed with exit code {exc.returncode}", file=sys.stderr)
        return exc.returncode or 1
    except Exception as exc:  # noqa: BLE001 - plotting script should report failures cleanly.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    edge_counts = {"supports": 0, "qualifies": 0, "undermines": 0}
    for _, _, relation in collect_edges(graph.get("claims", []), complete_edges=args.complete_edges):
        edge_counts[relation] += 1

    print(f"Wrote {args.output}")
    print(f"Wrote {args.dot_output}")
    print(
        "Edges: "
        f"{edge_counts['supports']} supports, "
        f"{edge_counts['qualifies']} qualifies, "
        f"{edge_counts['undermines']} undermines"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
