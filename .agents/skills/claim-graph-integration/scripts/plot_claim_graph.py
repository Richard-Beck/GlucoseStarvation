#!/usr/bin/env python3
"""Render a clean claim graph as Graphviz DOT and PNG."""

import argparse
import json
import shutil
import subprocess
import textwrap
from pathlib import Path


RELATION_STYLE = {
    "support": {"color": "#2E7D32", "arrowhead": "normal"},
    "undermine": {"color": "#C62828", "arrowhead": "tee"},
}
STRENGTH_STYLE = {
    "strong": {"penwidth": "2.5", "style": "solid"},
    "moderate": {"penwidth": "1.7", "style": "solid"},
    "weak": {"penwidth": "1.2", "style": "dashed"},
}


def dot_escape(value):
    return (
        str(value)
        .replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", "\\n")
    )


def wrapped(value, width=42):
    return "\n".join(textwrap.wrap(str(value), width=width))


def node_name(prefix, node_id):
    return f"{prefix}_{node_id}"


def graph_to_dot(graph):
    lines = [
        "digraph claim_graph {",
        '  graph [rankdir="LR", bgcolor="white", nodesep="0.25", ranksep="0.55"];',
        '  node [fontname="Helvetica", fontsize="9"];',
        '  edge [fontname="Helvetica", fontsize="8"];',
    ]

    claims = {claim["id"]: claim for claim in graph.get("claims", [])}
    evidence = {item["id"]: item for item in graph.get("evidence", [])}

    for claim_id, claim in claims.items():
        fixed = claim.get("user_fixed") is True
        label = f"{claim_id}\n{wrapped(claim.get('text', ''))}"
        attrs = {
            "shape": "ellipse",
            "style": "filled",
            "fillcolor": "#FFF2B2" if fixed else "#E8EEF7",
            "color": "#8A6D00" if fixed else "#506784",
            "penwidth": "2.2" if fixed else "1.2",
            "label": label,
        }
        rendered = ", ".join(
            f'{key}="{dot_escape(value)}"' for key, value in attrs.items()
        )
        lines.append(f'  "{dot_escape(node_name("claim", claim_id))}" [{rendered}];')

    for evidence_id, item in evidence.items():
        panels = ",".join(str(panel) for panel in item.get("panels", []))
        label = (
            f"{evidence_id}\n{item.get('figure_id', '')}/{panels}\n"
            f"{wrapped(item.get('observation', ''), width=36)}"
        )
        attrs = {
            "shape": "box",
            "style": "rounded,filled",
            "fillcolor": "#F5F5F5",
            "color": "#777777",
            "label": label,
        }
        rendered = ", ".join(
            f'{key}="{dot_escape(value)}"' for key, value in attrs.items()
        )
        lines.append(
            f'  "{dot_escape(node_name("evidence", evidence_id))}" [{rendered}];'
        )

    for relationship in graph.get("relationships", []):
        source_type = relationship["source_type"]
        source_id = relationship["source_id"]
        target_id = relationship["target_claim_id"]
        relation = relationship["relation"]
        strength = relationship["strength"]
        relation_style = RELATION_STYLE[relation]
        strength_style = STRENGTH_STYLE[strength]
        source = node_name("evidence" if source_type == "evidence" else "claim", source_id)
        target = node_name("claim", target_id)
        attrs = {
            "color": relation_style["color"],
            "arrowhead": relation_style["arrowhead"],
            "penwidth": strength_style["penwidth"],
            "style": strength_style["style"],
            "label": f"{relation} · {strength}",
        }
        rendered = ", ".join(
            f'{key}="{dot_escape(value)}"' for key, value in attrs.items()
        )
        lines.append(
            f'  "{dot_escape(source)}" -> "{dot_escape(target)}" [{rendered}];'
        )

    lines.append("}")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dot-output", required=True, type=Path)
    args = parser.parse_args()

    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    dot_text = graph_to_dot(graph)
    args.dot_output.parent.mkdir(parents=True, exist_ok=True)
    args.dot_output.write_text(dot_text, encoding="utf-8")

    dot = shutil.which("dot")
    if dot is None:
        raise SystemExit(
            f"Graphviz dot is unavailable; DOT written to {args.dot_output}"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [dot, "-Tpng", str(args.dot_output), "-o", str(args.output)],
        check=True,
    )
    print(f"Wrote {args.dot_output} and {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
