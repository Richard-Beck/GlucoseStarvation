#!/usr/bin/env python3
"""Render the plan-redraft workflow DAG JSON to a PNG with Graphviz."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path


def dot_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def wrapped_label(value: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False))


def load_dag(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    nodes = data.get("nodes")
    edges = data.get("edges")
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise ValueError("DAG JSON must contain list fields: nodes and edges")

    node_ids: set[str] = set()
    normalized_nodes: list[dict[str, str]] = []
    for node in nodes:
        if not isinstance(node, dict):
            raise ValueError("Each node must be an object")
        node_id = node.get("id")
        label = node.get("label")
        if not isinstance(node_id, str) or not isinstance(label, str):
            raise ValueError("Each node must have string id and label fields")
        if node_id in node_ids:
            raise ValueError(f"Duplicate node id: {node_id}")
        node_ids.add(node_id)
        normalized_nodes.append({"id": node_id, "label": label})

    normalized_edges: list[dict[str, str]] = []
    for edge in edges:
        if not isinstance(edge, dict):
            raise ValueError("Each edge must be an object")
        source = edge.get("from")
        target = edge.get("to")
        artifact = edge.get("artifact")
        if not all(isinstance(value, str) for value in (source, target, artifact)):
            raise ValueError("Each edge must have string from, to, and artifact fields")
        missing = [endpoint for endpoint in (source, target) if endpoint not in node_ids]
        if missing:
            raise ValueError(f"Edge references unknown node(s): {', '.join(missing)}")
        normalized_edges.append({"from": source, "to": target, "artifact": artifact})

    return normalized_nodes, normalized_edges


def build_dot(nodes: list[dict[str, str]], edges: list[dict[str, str]]) -> str:
    lines = [
        "digraph redraft_workflow {",
        '  graph [rankdir=TB, bgcolor="white", pad="0.25", nodesep="0.45", ranksep="0.75"];',
        '  node [shape=box, style="rounded,filled", fillcolor="#F8FAFC", color="#475569", penwidth=1.3, fontname="Helvetica", fontsize=12, margin="0.12,0.08"];',
        '  edge [color="#64748B", fontcolor="#334155", fontname="Helvetica", fontsize=10, arrowsize=0.75, penwidth=1.2];',
    ]
    for node in nodes:
        label = dot_escape(wrapped_label(node["label"], width=26))
        lines.append(f'  "{dot_escape(node["id"])}" [label="{label}"];')
    for edge in edges:
        source = dot_escape(edge["from"])
        target = dot_escape(edge["to"])
        label = dot_escape(wrapped_label(edge["artifact"], width=22))
        lines.append(f'  "{source}" -> "{target}" [label="{label}"];')
    lines.append("}")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    skill_dir = script_dir.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=skill_dir / "references" / "redraft_workflow_dag.json",
        help="Path to the workflow DAG JSON.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=skill_dir / "assets" / "redraft_workflow_dag.png",
        help="Path for the rendered PNG.",
    )
    parser.add_argument(
        "--dot-out",
        type=Path,
        default=None,
        help="Optional path to also write the intermediate DOT source.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if shutil.which("dot") is None:
        print("Graphviz 'dot' executable was not found on PATH.", file=sys.stderr)
        return 1

    nodes, edges = load_dag(args.input)
    dot_source = build_dot(nodes, edges)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.dot_out is not None:
        args.dot_out.parent.mkdir(parents=True, exist_ok=True)
        args.dot_out.write_text(dot_source, encoding="utf-8")

    subprocess.run(
        ["dot", "-Tpng", "-o", str(args.output)],
        input=dot_source,
        text=True,
        check=True,
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
