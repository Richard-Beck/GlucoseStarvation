#!/usr/bin/env python3
"""Build a computation-centered method-atom graph from classified provenance.

Inputs are:

* a Markdown provenance table with ``id | parent | what | why | comment``;
* a Markdown classification table with ``id | classification``.

The script treats each ``computation_step`` row as a method atom center. It then
walks upstream and downstream through non-computation nodes until another
computation step or graph boundary is reached. This produces an overlapping atom
graph: bridge artifacts can be both outputs of one atom and inputs to another.
"""

from __future__ import annotations

import argparse
import html
import re
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path


METHOD_HEADER = ["id", "parent", "what", "why", "comment"]
CLASSIFICATION_HEADER = ["id", "classification"]

VALID_CLASSIFICATIONS = {
    "display_node",
    "source_node",
    "data_artifact",
    "computation_step",
    "code_or_config_reference",
    "provenance_marker",
}

TERMINAL_VALUES = {
    "",
    "na",
    "n/a",
    "na (terminal)",
    "none",
    "null",
    "terminal",
}

DEFAULT_PANEL_LIKE_PATTERNS = [
    r"\bmanuscript-visible\b",
    r"\bmanuscript endpoint\b",
    r"\bmanifest row\b",
    r"\bmanifest-defined\b",
    r"\bmanifested panel\b",
    r"\bintegrated panel\b",
    r"\bregistered the manuscript\b",
    r"\bregistered .* manuscript endpoint\b",
    r"\bnormalized .* manuscript\b",
    r"\bfigure[_ ]s?\d+[a-z_]*_integrated_panel\b",
    r"\bfigure[_ ]s?\d+[a-z_]*_manuscript_panel\b",
]


@dataclass(frozen=True)
class Row:
    line_number: int
    row_index: int
    node_id: str
    parent_raw: str
    what: str
    why: str
    comment: str


@dataclass(frozen=True)
class BoundaryWalk:
    region_nodes: frozenset[str]
    boundary_computations: frozenset[str]


def clean_cell(text: str) -> str:
    text = html.unescape(text)
    text = text.replace("<br>", ";").replace("<br/>", ";").replace("<br />", ";")
    text = text.replace("&nbsp;", " ")
    text = text.replace("\\|", "|")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def strip_single_code_span(text: str) -> str:
    text = clean_cell(text)
    match = re.fullmatch(r"`([^`]+)`", text)
    if match:
        return match.group(1).strip()
    return text.strip()


def normalize_for_terminal(text: str) -> str:
    text = strip_single_code_span(text)
    text = text.strip().strip("\"'").strip()
    text = re.sub(r"\s+", " ", text)
    return text.lower()


def is_terminal_boundary(text: str) -> bool:
    norm = normalize_for_terminal(text).strip(" <>")
    return norm in TERMINAL_VALUES or norm.startswith("terminal:")


def split_markdown_row(line: str) -> list[str]:
    line = line.strip()
    if not (line.startswith("|") and line.endswith("|")):
        return []
    return [clean_cell(cell) for cell in line.strip("|").split("|")]


def parse_method_table(path: Path) -> list[Row]:
    rows: list[Row] = []
    saw_header = False
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.startswith("|"):
            continue
        cells = split_markdown_row(line)
        if len(cells) != 5:
            continue
        lower = [normalize_for_terminal(cell) for cell in cells]
        if lower == METHOD_HEADER:
            saw_header = True
            continue
        if not saw_header:
            continue
        if all(set(cell.replace(" ", "")) <= {"-"} for cell in cells):
            continue
        rows.append(
            Row(
                line_number=line_number,
                row_index=len(rows) + 1,
                node_id=strip_single_code_span(cells[0]),
                parent_raw=cells[1],
                what=cells[2],
                why=cells[3],
                comment=cells[4],
            )
        )
    return rows


def parse_classification_table(path: Path) -> dict[str, str]:
    classifications: dict[str, str] = {}
    saw_header = False
    duplicate_ids: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            continue
        cells = split_markdown_row(line)
        if len(cells) != 2:
            continue
        lower = [normalize_for_terminal(cell) for cell in cells]
        if lower == CLASSIFICATION_HEADER:
            saw_header = True
            continue
        if not saw_header:
            continue
        if all(set(cell.replace(" ", "")) <= {"-"} for cell in cells):
            continue
        node_id = strip_single_code_span(cells[0])
        classification = strip_single_code_span(cells[1])
        if node_id in classifications:
            duplicate_ids.append(node_id)
        classifications[node_id] = classification

    if duplicate_ids:
        shown = ", ".join(sorted(set(duplicate_ids))[:5])
        suffix = "" if len(set(duplicate_ids)) <= 5 else f" (+{len(set(duplicate_ids)) - 5} more)"
        raise ValueError(f"duplicate classification ids: {shown}{suffix}")
    invalid = sorted({kind for kind in classifications.values() if kind not in VALID_CLASSIFICATIONS})
    if invalid:
        raise ValueError(f"invalid classifications: {', '.join(invalid)}")
    return classifications


def split_parent_field(text: str) -> list[str]:
    text = clean_cell(text)
    if not text:
        return []
    parts = [part.strip() for part in text.split(";")]
    return [strip_single_code_span(part) for part in parts if strip_single_code_span(part)]


def build_graph(rows: list[Row]) -> tuple[dict[str, Row], dict[str, set[str]], dict[str, set[str]], list[str], dict[str, set[str]]]:
    by_id: dict[str, Row] = {}
    duplicate_ids: list[str] = []
    for row in rows:
        if row.node_id in by_id:
            duplicate_ids.append(row.node_id)
        by_id[row.node_id] = row

    ids = set(by_id)
    children: dict[str, set[str]] = defaultdict(set)
    parents: dict[str, set[str]] = defaultdict(set)
    external_parents: dict[str, set[str]] = defaultdict(set)

    for row in rows:
        for parent in split_parent_field(row.parent_raw):
            if parent in ids:
                children[parent].add(row.node_id)
                parents[row.node_id].add(parent)
            elif not is_terminal_boundary(parent):
                external_parents[row.node_id].add(parent)

    return by_id, children, parents, duplicate_ids, external_parents


def root_distances(by_id: dict[str, Row], children: dict[str, set[str]], parents: dict[str, set[str]]) -> dict[str, int]:
    roots = sorted(node_id for node_id in by_id if not parents.get(node_id))
    distance: dict[str, int] = {}
    queue: deque[tuple[str, int]] = deque((root, 0) for root in roots)
    while queue:
        node_id, dist = queue.popleft()
        if node_id in distance and distance[node_id] <= dist:
            continue
        distance[node_id] = dist
        for child in sorted(children.get(node_id, set())):
            queue.append((child, dist + 1))
    return distance


def descendants(start: str, children: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    queue: deque[str] = deque(sorted(children.get(start, set())))
    while queue:
        node_id = queue.popleft()
        if node_id in seen:
            continue
        seen.add(node_id)
        queue.extend(sorted(children.get(node_id, set()) - seen))
    return seen


def compile_patterns(patterns: list[str]) -> list[re.Pattern[str]]:
    return [re.compile(pattern, re.IGNORECASE) for pattern in patterns]


def is_panel_like_leaf(row: Row, children: dict[str, set[str]], patterns: list[re.Pattern[str]]) -> bool:
    if children.get(row.node_id):
        return False
    text = " ".join([row.node_id, row.what, row.why, row.comment]).lower()
    return any(pattern.search(text) for pattern in patterns)


def walk_boundary(
    center: str,
    adjacent: dict[str, set[str]],
    classifications: dict[str, str],
) -> BoundaryWalk:
    region_nodes: set[str] = set()
    boundary_computations: set[str] = set()
    queue: deque[str] = deque(sorted(adjacent.get(center, set())))
    while queue:
        node_id = queue.popleft()
        if classifications[node_id] == "computation_step":
            boundary_computations.add(node_id)
            continue
        if node_id in region_nodes:
            continue
        region_nodes.add(node_id)
        queue.extend(sorted(adjacent.get(node_id, set())))
    return BoundaryWalk(frozenset(region_nodes), frozenset(boundary_computations))


def bridge_nodes_between(
    upstream_center: str,
    downstream_center: str,
    downstream_region: set[str],
    parents: dict[str, set[str]],
    classifications: dict[str, str],
) -> set[str]:
    bridge_nodes: set[str] = set()
    queue: deque[str] = deque(sorted(parents.get(downstream_center, set())))
    while queue:
        node_id = queue.popleft()
        if node_id == upstream_center:
            continue
        if classifications[node_id] == "computation_step":
            continue
        if node_id not in downstream_region or node_id in bridge_nodes:
            continue
        bridge_nodes.add(node_id)
        queue.extend(sorted(parents.get(node_id, set())))
    return bridge_nodes


def bool_text(value: bool) -> str:
    return "yes" if value else "no"


def tsv_text(value: object) -> str:
    text = "" if value is None else str(value)
    return text.replace("\t", " ").replace("\n", " ").replace("\r", " ")


def join_ids(values: set[str] | list[str] | tuple[str, ...]) -> str:
    return "; ".join(sorted(values))


def nodes_of_kind(nodes: set[str], classifications: dict[str, str], kind: str) -> set[str]:
    return {node_id for node_id in nodes if classifications[node_id] == kind}


def write_tsv(path: Path, rows: list[dict[str, object]], headers: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        print("\t".join(headers), file=handle)
        for row in rows:
            print("\t".join(tsv_text(row.get(header, "")) for header in headers), file=handle)


def write_summary(
    path: Path,
    *,
    table_path: Path,
    classification_path: Path,
    node_count: int,
    edge_count: int,
    classification_counts: Counter[str],
    atom_rows: list[dict[str, object]],
    edge_rows: list[dict[str, object]],
    core_membership_rows: list[dict[str, object]],
    context_membership_rows: list[dict[str, object]],
    uncovered_rows: list[dict[str, object]],
    proxy_candidate_rows: list[dict[str, object]],
    panel_like_leaves: set[str],
) -> None:
    top_atoms = sorted(
        atom_rows,
        key=lambda row: (
            -int(row["manuscript_panel_like_leaf_descendants"]),
            -int(row["leaf_descendants"]),
            int(row["center_root_distance"]) if row["center_root_distance"] != "" else 10**9,
            str(row["center_id"]),
        ),
    )[:15]
    with path.open("w", encoding="utf-8", newline="") as handle:
        print("# Method Atom Graph Summary", file=handle)
        print("", file=handle)
        print(f"Input provenance table: `{table_path}`", file=handle)
        print(f"Input classification table: `{classification_path}`", file=handle)
        print("", file=handle)
        print(f"Parsed provenance nodes: {node_count}", file=handle)
        print(f"Parsed graph edges: {edge_count}", file=handle)
        print(f"Method atoms: {len(atom_rows)}", file=handle)
        print(f"Atom-to-atom edges: {len(edge_rows)}", file=handle)
        print(f"Atom core membership rows: {len(core_membership_rows)}", file=handle)
        print(f"Atom context membership rows: {len(context_membership_rows)}", file=handle)
        print(f"Uncovered non-computation nodes: {len(uncovered_rows)}", file=handle)
        print(f"Non-computation proxy candidate regions: {len(proxy_candidate_rows)}", file=handle)
        print(f"Manuscript/panel-like leaves: {len(panel_like_leaves)}", file=handle)
        print("", file=handle)
        print("Classification counts:", file=handle)
        for kind in sorted(classification_counts):
            print(f"- `{kind}`: {classification_counts[kind]}", file=handle)
        print("", file=handle)
        print(
            "Atom definition: each `computation_step` is an atom center. "
            "Core membership is limited to the center plus direct non-computation "
            "parents and children. Longer upstream/downstream non-computation "
            "walks are reported as context and bridge paths instead of being "
            "folded into the atom core. Boundary artifacts can therefore appear "
            "in more than one atom or edge context.",
            file=handle,
        )
        print("", file=handle)
        print("## Top Atoms By Manuscript/Panel-Like Leaf Reach", file=handle)
        print("", file=handle)
        headers = [
            "rank",
            "atom_id",
            "center_id",
            "center_root_distance",
            "upstream_atoms",
            "downstream_atoms",
            "leaf_descendants",
            "manuscript_panel_like_leaf_descendants",
        ]
        print("| " + " | ".join(headers) + " |", file=handle)
        print("|" + "|".join("---" for _ in headers) + "|", file=handle)
        for rank, row in enumerate(top_atoms, start=1):
            values = [
                str(rank),
                f"`{row['atom_id']}`",
                f"`{str(row['center_id']).replace('`', '')}`",
                str(row["center_root_distance"]),
                str(row["upstream_atom_count"]),
                str(row["downstream_atom_count"]),
                str(row["leaf_descendants"]),
                str(row["manuscript_panel_like_leaf_descendants"]),
            ]
            print("| " + " | ".join(values) + " |", file=handle)
        if proxy_candidate_rows:
            print("", file=handle)
            print("## Top Non-Computation Proxy Candidate Regions", file=handle)
            print("", file=handle)
            proxy_headers = [
                "rank",
                "proxy_atom_id",
                "candidate_kind",
                "node_count",
                "manuscript_panel_like_leaf_count",
                "anchor_node_ids",
            ]
            print("| " + " | ".join(proxy_headers) + " |", file=handle)
            print("|" + "|".join("---" for _ in proxy_headers) + "|", file=handle)
            for rank, row in enumerate(proxy_candidate_rows[:15], start=1):
                values = [
                    str(rank),
                    f"`{row['proxy_atom_id']}`",
                    f"`{row['candidate_kind']}`",
                    str(row["node_count"]),
                    str(row["manuscript_panel_like_leaf_count"]),
                    f"`{str(row['anchor_node_ids']).replace('`', '')}`",
                ]
                print("| " + " | ".join(values) + " |", file=handle)


def build_outputs(
    rows: list[Row],
    classifications: dict[str, str],
    by_id: dict[str, Row],
    children: dict[str, set[str]],
    parents: dict[str, set[str]],
    distances: dict[str, int],
    panel_like_leaves: set[str],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    computation_nodes = {node_id for node_id, kind in classifications.items() if kind == "computation_step"}
    atom_id_by_center = {
        row.node_id: f"atom_{row.row_index:04d}"
        for row in rows
        if row.node_id in computation_nodes
    }

    upstream_walks: dict[str, BoundaryWalk] = {}
    downstream_walks: dict[str, BoundaryWalk] = {}
    for center_id in sorted(computation_nodes):
        upstream_walks[center_id] = walk_boundary(center_id, parents, classifications)
        downstream_walks[center_id] = walk_boundary(center_id, children, classifications)

    atom_rows: list[dict[str, object]] = []
    core_membership_rows: list[dict[str, object]] = []
    context_membership_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    covered_nodes: set[str] = set()

    for center_id in sorted(computation_nodes, key=lambda node_id: by_id[node_id].row_index):
        atom_id = atom_id_by_center[center_id]
        row = by_id[center_id]
        upstream_walk = upstream_walks[center_id]
        downstream_walk = downstream_walks[center_id]
        upstream_context_nodes = set(upstream_walk.region_nodes)
        downstream_context_nodes = set(downstream_walk.region_nodes)
        direct_upstream_nodes = {
            node_id
            for node_id in parents.get(center_id, set())
            if classifications[node_id] != "computation_step"
        }
        direct_downstream_nodes = {
            node_id
            for node_id in children.get(center_id, set())
            if classifications[node_id] != "computation_step"
        }
        core_nodes = {center_id} | direct_upstream_nodes | direct_downstream_nodes
        covered_nodes.update({center_id} | upstream_context_nodes | downstream_context_nodes)

        downstream = descendants(center_id, children)
        leaf_descendants = {node_id for node_id in downstream if not children.get(node_id)}
        panel_descendants = leaf_descendants & panel_like_leaves

        upstream_atoms = {atom_id_by_center[node_id] for node_id in upstream_walk.boundary_computations}
        downstream_atoms = {atom_id_by_center[node_id] for node_id in downstream_walk.boundary_computations}

        atom_rows.append(
            {
                "atom_id": atom_id,
                "atom_kind": "computation_atom",
                "center_id": center_id,
                "center_line": row.line_number,
                "center_root_distance": distances.get(center_id, ""),
                "upstream_atom_count": len(upstream_atoms),
                "downstream_atom_count": len(downstream_atoms),
                "core_node_count": len(core_nodes),
                "direct_input_node_count": len(direct_upstream_nodes),
                "direct_output_node_count": len(direct_downstream_nodes),
                "upstream_context_node_count": len(upstream_context_nodes),
                "downstream_context_node_count": len(downstream_context_nodes),
                "descendant_rows": len(downstream),
                "leaf_descendants": len(leaf_descendants),
                "manuscript_panel_like_leaf_descendants": len(panel_descendants),
                "upstream_atom_ids": join_ids(upstream_atoms),
                "downstream_atom_ids": join_ids(downstream_atoms),
                "direct_input_source_nodes": join_ids(nodes_of_kind(direct_upstream_nodes, classifications, "source_node")),
                "direct_input_data_artifacts": join_ids(nodes_of_kind(direct_upstream_nodes, classifications, "data_artifact")),
                "direct_input_code_or_config_references": join_ids(nodes_of_kind(direct_upstream_nodes, classifications, "code_or_config_reference")),
                "direct_input_provenance_markers": join_ids(nodes_of_kind(direct_upstream_nodes, classifications, "provenance_marker")),
                "direct_input_display_nodes": join_ids(nodes_of_kind(direct_upstream_nodes, classifications, "display_node")),
                "direct_output_data_artifacts": join_ids(nodes_of_kind(direct_downstream_nodes, classifications, "data_artifact")),
                "direct_output_display_nodes": join_ids(nodes_of_kind(direct_downstream_nodes, classifications, "display_node")),
                "direct_output_code_or_config_references": join_ids(nodes_of_kind(direct_downstream_nodes, classifications, "code_or_config_reference")),
                "direct_output_provenance_markers": join_ids(nodes_of_kind(direct_downstream_nodes, classifications, "provenance_marker")),
                "direct_output_source_nodes": join_ids(nodes_of_kind(direct_downstream_nodes, classifications, "source_node")),
                "manuscript_panel_like_leaf_ids": join_ids(panel_descendants),
                "center_what": row.what,
                "center_why": row.why,
                "center_comment": row.comment,
            }
        )

        relation_specs = [
            ("center", {center_id}),
            ("direct_input", direct_upstream_nodes),
            ("direct_output", direct_downstream_nodes),
        ]
        for relation, node_ids in relation_specs:
            for node_id in sorted(node_ids, key=lambda value: by_id[value].row_index):
                core_membership_rows.append(
                    {
                        "atom_id": atom_id,
                        "center_id": center_id,
                        "relation": relation,
                        "node_id": node_id,
                        "classification": classifications[node_id],
                        "line": by_id[node_id].line_number,
                        "root_distance": distances.get(node_id, ""),
                    }
                )
        context_relation_specs = [
            ("upstream_context", upstream_context_nodes),
            ("downstream_context", downstream_context_nodes),
        ]
        for relation, node_ids in context_relation_specs:
            for node_id in sorted(node_ids, key=lambda value: by_id[value].row_index):
                context_membership_rows.append(
                    {
                        "atom_id": atom_id,
                        "center_id": center_id,
                        "relation": relation,
                        "node_id": node_id,
                        "classification": classifications[node_id],
                        "line": by_id[node_id].line_number,
                        "root_distance": distances.get(node_id, ""),
                    }
                )

    seen_edges: set[tuple[str, str]] = set()
    for center_id in sorted(computation_nodes, key=lambda node_id: by_id[node_id].row_index):
        source_atom = atom_id_by_center[center_id]
        downstream_region = set(downstream_walks[center_id].region_nodes)
        for target_center_id in sorted(downstream_walks[center_id].boundary_computations, key=lambda node_id: by_id[node_id].row_index):
            target_atom = atom_id_by_center[target_center_id]
            edge_key = (source_atom, target_atom)
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)
            bridge_nodes = bridge_nodes_between(
                center_id,
                target_center_id,
                downstream_region,
                parents,
                classifications,
            )
            edge_rows.append(
                {
                    "source_atom_id": source_atom,
                    "source_center_id": center_id,
                    "target_atom_id": target_atom,
                    "target_center_id": target_center_id,
                    "bridge_node_count": len(bridge_nodes),
                    "bridge_data_artifacts": join_ids(nodes_of_kind(bridge_nodes, classifications, "data_artifact")),
                    "bridge_display_nodes": join_ids(nodes_of_kind(bridge_nodes, classifications, "display_node")),
                    "bridge_code_or_config_references": join_ids(nodes_of_kind(bridge_nodes, classifications, "code_or_config_reference")),
                    "bridge_provenance_markers": join_ids(nodes_of_kind(bridge_nodes, classifications, "provenance_marker")),
                    "bridge_source_nodes": join_ids(nodes_of_kind(bridge_nodes, classifications, "source_node")),
                    "all_bridge_node_ids": join_ids(bridge_nodes),
                }
            )

    uncovered_rows: list[dict[str, object]] = []
    for row in rows:
        if row.node_id in covered_nodes:
            continue
        if classifications[row.node_id] == "computation_step":
            continue
        downstream = descendants(row.node_id, children)
        leaf_descendants = {node_id for node_id in downstream if not children.get(node_id)}
        panel_descendants = leaf_descendants & panel_like_leaves
        uncovered_rows.append(
            {
                "node_id": row.node_id,
                "classification": classifications[row.node_id],
                "line": row.line_number,
                "root_distance": distances.get(row.node_id, ""),
                "direct_parent_count": len(parents.get(row.node_id, set())),
                "direct_child_count": len(children.get(row.node_id, set())),
                "descendant_rows": len(downstream),
                "leaf_descendants": len(leaf_descendants),
                "manuscript_panel_like_leaf_descendants": len(panel_descendants),
                "manuscript_panel_like_leaf_ids": join_ids(panel_descendants),
                "what": row.what,
                "why": row.why,
            }
        )

    atom_rows.sort(
        key=lambda item: (
            -int(item["manuscript_panel_like_leaf_descendants"]),
            -int(item["leaf_descendants"]),
            int(item["center_root_distance"]) if item["center_root_distance"] != "" else 10**9,
            str(item["center_id"]),
        )
    )
    edge_rows.sort(key=lambda item: (str(item["source_atom_id"]), str(item["target_atom_id"])))
    core_membership_rows.sort(key=lambda item: (str(item["atom_id"]), str(item["relation"]), int(item["line"])))
    context_membership_rows.sort(key=lambda item: (str(item["atom_id"]), str(item["relation"]), int(item["line"])))
    uncovered_rows.sort(
        key=lambda item: (
            -int(item["manuscript_panel_like_leaf_descendants"]),
            -int(item["leaf_descendants"]),
            int(item["root_distance"]) if item["root_distance"] != "" else 10**9,
            str(item["node_id"]),
        )
    )
    return atom_rows, edge_rows, core_membership_rows, context_membership_rows, uncovered_rows


def candidate_kind_for_component(component: set[str], classifications: dict[str, str]) -> str:
    kinds = Counter(classifications[node_id] for node_id in component)
    if kinds["data_artifact"]:
        return "data_artifact_region"
    if kinds["code_or_config_reference"]:
        return "code_or_config_region"
    if kinds["source_node"]:
        return "source_region"
    if kinds["provenance_marker"]:
        return "provenance_marker_region"
    return "display_layout_region"


def proxy_anchor_sort_key(
    node_id: str,
    by_id: dict[str, Row],
    classifications: dict[str, str],
    children: dict[str, set[str]],
    panel_like_leaves: set[str],
    distances: dict[str, int],
) -> tuple[int, int, int, int, str]:
    priority = {
        "data_artifact": 0,
        "source_node": 1,
        "code_or_config_reference": 2,
        "provenance_marker": 3,
        "display_node": 4,
    }
    downstream = descendants(node_id, children)
    panel_count = len((downstream | ({node_id} if node_id in panel_like_leaves else set())) & panel_like_leaves)
    return (
        priority.get(classifications[node_id], 99),
        -panel_count,
        distances.get(node_id, 10**9),
        by_id[node_id].row_index,
        node_id,
    )


def build_proxy_candidates(
    uncovered_rows: list[dict[str, object]],
    by_id: dict[str, Row],
    classifications: dict[str, str],
    children: dict[str, set[str]],
    parents: dict[str, set[str]],
    distances: dict[str, int],
    panel_like_leaves: set[str],
) -> list[dict[str, object]]:
    uncovered_ids = {str(row["node_id"]) for row in uncovered_rows}
    visited: set[str] = set()
    components: list[set[str]] = []
    for node_id in sorted(uncovered_ids, key=lambda value: by_id[value].row_index):
        if node_id in visited:
            continue
        component: set[str] = set()
        queue: deque[str] = deque([node_id])
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
            visited.add(current)
            component.add(current)
            neighbors = (children.get(current, set()) | parents.get(current, set())) & uncovered_ids
            queue.extend(sorted(neighbors - visited))
        components.append(component)

    rows: list[dict[str, object]] = []
    for component_index, component in enumerate(components, start=1):
        support: set[str] = set()
        leaf_support: set[str] = set()
        for node_id in component:
            downstream = descendants(node_id, children)
            leaves = {value for value in downstream if not children.get(value)}
            if not children.get(node_id):
                leaves.add(node_id)
            leaf_support.update(leaves)
            support.update(leaves & panel_like_leaves)
        if not support:
            continue
        kinds = Counter(classifications[node_id] for node_id in component)
        anchors = sorted(
            component,
            key=lambda value: proxy_anchor_sort_key(
                value,
                by_id,
                classifications,
                children,
                panel_like_leaves,
                distances,
            ),
        )[:5]
        boundary_parents = sorted(
            {
                parent_id
                for node_id in component
                for parent_id in parents.get(node_id, set())
                if parent_id not in component
            },
            key=lambda value: by_id[value].row_index,
        )
        boundary_children = sorted(
            {
                child_id
                for node_id in component
                for child_id in children.get(node_id, set())
                if child_id not in component
            },
            key=lambda value: by_id[value].row_index,
        )
        root_values = [distances[node_id] for node_id in component if node_id in distances]
        rows.append(
            {
                "proxy_atom_id": f"proxy_{component_index:03d}",
                "candidate_kind": candidate_kind_for_component(component, classifications),
                "node_count": len(component),
                "source_node_count": kinds["source_node"],
                "data_artifact_count": kinds["data_artifact"],
                "code_or_config_reference_count": kinds["code_or_config_reference"],
                "display_node_count": kinds["display_node"],
                "provenance_marker_count": kinds["provenance_marker"],
                "min_root_distance": min(root_values) if root_values else "",
                "max_root_distance": max(root_values) if root_values else "",
                "leaf_count": len(leaf_support),
                "manuscript_panel_like_leaf_count": len(support),
                "manuscript_panel_like_leaf_ids": join_ids(support),
                "anchor_node_ids": join_ids(set(anchors)),
                "all_node_ids": join_ids(component),
                "boundary_parent_ids": join_ids(set(boundary_parents)),
                "boundary_child_ids": join_ids(set(boundary_children)),
            }
        )

    rows.sort(
        key=lambda item: (
            -int(item["manuscript_panel_like_leaf_count"]),
            str(item["candidate_kind"]) == "display_layout_region",
            -int(item["data_artifact_count"]),
            -int(item["node_count"]),
            str(item["anchor_node_ids"]),
        )
    )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Reduce a classified method-table provenance graph to a deterministic "
            "computation-centered method-atom graph."
        )
    )
    parser.add_argument("table", type=Path, help="Markdown table with id|parent|what|why|comment columns.")
    parser.add_argument("classification", type=Path, help="Markdown table with id|classification columns.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for method atom graph outputs. Defaults to METHOD_TABLE_DIR/method_atom_graph.",
    )
    parser.add_argument(
        "--panel-like-regex",
        action="append",
        default=[],
        help=(
            "Additional case-insensitive regex used to classify leaf rows as manuscript/panel-like. "
            "May be supplied multiple times."
        ),
    )
    args = parser.parse_args(argv)

    if not args.table.exists() or not args.table.is_file():
        parser.error(f"table is not a file: {args.table}")
    if not args.classification.exists() or not args.classification.is_file():
        parser.error(f"classification is not a file: {args.classification}")

    rows = parse_method_table(args.table)
    if not rows:
        parser.error(f"no method-table rows parsed from: {args.table}")

    by_id, children, parents, duplicate_ids, _external_parents = build_graph(rows)
    if duplicate_ids:
        shown = ", ".join(sorted(set(duplicate_ids))[:5])
        suffix = "" if len(set(duplicate_ids)) <= 5 else f" (+{len(set(duplicate_ids)) - 5} more)"
        parser.error(f"duplicate row ids are not supported: {shown}{suffix}")

    classifications = parse_classification_table(args.classification)
    table_ids = set(by_id)
    classification_ids = set(classifications)
    missing = sorted(table_ids - classification_ids)
    extra = sorted(classification_ids - table_ids)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"{len(missing)} table ids missing classifications, first: {missing[:5]}")
        if extra:
            parts.append(f"{len(extra)} classification ids absent from table, first: {extra[:5]}")
        parser.error("; ".join(parts))

    panel_patterns = compile_patterns(DEFAULT_PANEL_LIKE_PATTERNS + args.panel_like_regex)
    panel_like_leaves = {
        row.node_id for row in rows if is_panel_like_leaf(row, children, panel_patterns)
    }
    distances = root_distances(by_id, children, parents)
    atom_rows, edge_rows, core_membership_rows, context_membership_rows, uncovered_rows = build_outputs(
        rows,
        classifications,
        by_id,
        children,
        parents,
        distances,
        panel_like_leaves,
    )
    proxy_candidate_rows = build_proxy_candidates(
        uncovered_rows,
        by_id,
        classifications,
        children,
        parents,
        distances,
        panel_like_leaves,
    )

    output_dir = args.output_dir or args.table.parent / "method_atom_graph"
    output_dir.mkdir(parents=True, exist_ok=True)

    write_tsv(
        output_dir / "method_atoms.tsv",
        atom_rows,
        [
            "atom_id",
            "atom_kind",
            "center_id",
            "center_line",
            "center_root_distance",
            "upstream_atom_count",
            "downstream_atom_count",
            "core_node_count",
            "direct_input_node_count",
            "direct_output_node_count",
            "upstream_context_node_count",
            "downstream_context_node_count",
            "descendant_rows",
            "leaf_descendants",
            "manuscript_panel_like_leaf_descendants",
            "upstream_atom_ids",
            "downstream_atom_ids",
            "direct_input_source_nodes",
            "direct_input_data_artifacts",
            "direct_input_code_or_config_references",
            "direct_input_provenance_markers",
            "direct_input_display_nodes",
            "direct_output_data_artifacts",
            "direct_output_display_nodes",
            "direct_output_code_or_config_references",
            "direct_output_provenance_markers",
            "direct_output_source_nodes",
            "manuscript_panel_like_leaf_ids",
            "center_what",
            "center_why",
            "center_comment",
        ],
    )
    write_tsv(
        output_dir / "method_atom_edges.tsv",
        edge_rows,
        [
            "source_atom_id",
            "source_center_id",
            "target_atom_id",
            "target_center_id",
            "bridge_node_count",
            "bridge_data_artifacts",
            "bridge_display_nodes",
            "bridge_code_or_config_references",
            "bridge_provenance_markers",
            "bridge_source_nodes",
            "all_bridge_node_ids",
        ],
    )
    write_tsv(
        output_dir / "method_atom_node_membership.tsv",
        core_membership_rows,
        ["atom_id", "center_id", "relation", "node_id", "classification", "line", "root_distance"],
    )
    write_tsv(
        output_dir / "method_atom_context_membership.tsv",
        context_membership_rows,
        ["atom_id", "center_id", "relation", "node_id", "classification", "line", "root_distance"],
    )
    write_tsv(
        output_dir / "uncovered_noncomputation_nodes.tsv",
        uncovered_rows,
        [
            "node_id",
            "classification",
            "line",
            "root_distance",
            "direct_parent_count",
            "direct_child_count",
            "descendant_rows",
            "leaf_descendants",
            "manuscript_panel_like_leaf_descendants",
            "manuscript_panel_like_leaf_ids",
            "what",
            "why",
        ],
    )
    write_tsv(
        output_dir / "method_atom_proxy_candidates.tsv",
        proxy_candidate_rows,
        [
            "proxy_atom_id",
            "candidate_kind",
            "node_count",
            "source_node_count",
            "data_artifact_count",
            "code_or_config_reference_count",
            "display_node_count",
            "provenance_marker_count",
            "min_root_distance",
            "max_root_distance",
            "leaf_count",
            "manuscript_panel_like_leaf_count",
            "manuscript_panel_like_leaf_ids",
            "anchor_node_ids",
            "all_node_ids",
            "boundary_parent_ids",
            "boundary_child_ids",
        ],
    )
    write_summary(
        output_dir / "method_atom_graph_summary.md",
        table_path=args.table,
        classification_path=args.classification,
        node_count=len(rows),
        edge_count=sum(len(value) for value in children.values()),
        classification_counts=Counter(classifications.values()),
        atom_rows=atom_rows,
        edge_rows=edge_rows,
        core_membership_rows=core_membership_rows,
        context_membership_rows=context_membership_rows,
        uncovered_rows=uncovered_rows,
        proxy_candidate_rows=proxy_candidate_rows,
        panel_like_leaves=panel_like_leaves,
    )

    print(f"Wrote method atom graph outputs to {output_dir}", file=sys.stderr)
    print(f"Atoms: {len(atom_rows)}", file=sys.stderr)
    print(f"Atom edges: {len(edge_rows)}", file=sys.stderr)
    print(f"Uncovered non-computation nodes: {len(uncovered_rows)}", file=sys.stderr)
    print(f"Proxy candidate regions: {len(proxy_candidate_rows)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
