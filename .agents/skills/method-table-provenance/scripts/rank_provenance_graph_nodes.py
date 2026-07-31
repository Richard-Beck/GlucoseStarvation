#!/usr/bin/env python3
"""Rank method-table provenance nodes by downstream manuscript-panel reach.

The input is a Markdown table with columns:

    id | parent | what | why | comment

Edges are interpreted as ``parent -> id``. The script reports one row per
provenance node, ordered by decreasing count of downstream leaf nodes that look
like manuscript/panel endpoints.
"""

from __future__ import annotations

import argparse
import html
import re
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path


TABLE_HEADER = ["id", "parent", "what", "why", "comment"]

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
        if lower == TABLE_HEADER:
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


def split_parent_field(text: str) -> list[str]:
    text = clean_cell(text)
    if not text:
        return []
    parts = [part.strip() for part in text.split(";")]
    return [strip_single_code_span(part) for part in parts if strip_single_code_span(part)]


def build_graph(rows: list[Row]) -> tuple[dict[str, Row], dict[str, set[str]], dict[str, set[str]], list[str]]:
    by_id: dict[str, Row] = {}
    duplicate_ids: list[str] = []
    for row in rows:
        if row.node_id in by_id:
            duplicate_ids.append(row.node_id)
        by_id[row.node_id] = row

    ids = set(by_id)
    children: dict[str, set[str]] = defaultdict(set)
    parents: dict[str, set[str]] = defaultdict(set)

    for row in rows:
        for parent in split_parent_field(row.parent_raw):
            if parent in ids:
                children[parent].add(row.node_id)
                parents[row.node_id].add(parent)
            elif not is_terminal_boundary(parent):
                # Parent tokens not represented as rows are external boundaries for
                # root-distance purposes, but do not become nodes themselves.
                continue

    return by_id, children, parents, duplicate_ids


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


def is_panel_like_leaf(row: Row, children: dict[str, set[str]], patterns: list[re.Pattern[str]]) -> bool:
    if children.get(row.node_id):
        return False
    text = " ".join([row.node_id, row.what, row.why, row.comment]).lower()
    return any(pattern.search(text) for pattern in patterns)


def bool_text(value: bool) -> str:
    return "yes" if value else "no"


def markdown_escape(text: str) -> str:
    text = text.replace("\n", " ")
    text = text.replace("|", "\\|")
    return text


def make_metrics(
    rows: list[Row],
    by_id: dict[str, Row],
    children: dict[str, set[str]],
    parents: dict[str, set[str]],
    distances: dict[str, int],
    panel_patterns: list[re.Pattern[str]],
) -> list[dict[str, object]]:
    panel_like_leaves = {
        row.node_id for row in rows if is_panel_like_leaf(row, children, panel_patterns)
    }
    metrics: list[dict[str, object]] = []
    for row in rows:
        downstream = descendants(row.node_id, children)
        leaf_descendants = {node_id for node_id in downstream if not children.get(node_id)}
        panel_descendants = leaf_descendants & panel_like_leaves
        metrics.append(
            {
                "id": row.node_id,
                "line": row.line_number,
                "root_distance": distances.get(row.node_id),
                "is_root": not parents.get(row.node_id),
                "direct_children": len(children.get(row.node_id, set())),
                "descendant_rows": len(downstream),
                "leaf_descendants": len(leaf_descendants),
                "manuscript_panel_like_leaf_descendants": len(panel_descendants),
                "is_leaf": not children.get(row.node_id),
                "is_manuscript_panel_like_leaf": row.node_id in panel_like_leaves,
            }
        )
    metrics.sort(
        key=lambda item: (
            -int(item["manuscript_panel_like_leaf_descendants"]),
            -int(item["leaf_descendants"]),
            -int(item["descendant_rows"]),
            int(item["root_distance"]) if item["root_distance"] is not None else 10**9,
            str(item["id"]),
        )
    )
    return metrics


def write_markdown(
    metrics: list[dict[str, object]],
    output,
    input_path: Path,
    row_count: int,
    root_count: int,
    panel_like_leaf_count: int,
) -> None:
    print("# Provenance Graph Node Metrics", file=output)
    print("", file=output)
    print(f"Input table: `{input_path}`", file=output)
    print(f"Parsed rows: {row_count}", file=output)
    print(f"Roots: {root_count}", file=output)
    print(f"Manuscript/panel-like leaves: {panel_like_leaf_count}", file=output)
    print("", file=output)
    print(
        "`root_distance` is the shortest edge distance from any parsed root row. "
        "A root row has no nonterminal parent that is also represented as an `id` in the table.",
        file=output,
    )
    print(
        "`leaf_descendants` and `manuscript_panel_like_leaf_descendants` exclude the row itself.",
        file=output,
    )
    print("", file=output)

    headers = [
        "rank",
        "id",
        "line",
        "root_distance",
        "is_root",
        "direct_children",
        "descendant_rows",
        "leaf_descendants",
        "manuscript_panel_like_leaf_descendants",
        "is_leaf",
        "is_manuscript_panel_like_leaf",
    ]
    print("| " + " | ".join(headers) + " |", file=output)
    print("|" + "|".join("---" for _ in headers) + "|", file=output)
    for rank, metric in enumerate(metrics, start=1):
        row = [
            str(rank),
            "`" + markdown_escape(str(metric["id"]).replace("`", "")) + "`",
            str(metric["line"]),
            "" if metric["root_distance"] is None else str(metric["root_distance"]),
            bool_text(bool(metric["is_root"])),
            str(metric["direct_children"]),
            str(metric["descendant_rows"]),
            str(metric["leaf_descendants"]),
            str(metric["manuscript_panel_like_leaf_descendants"]),
            bool_text(bool(metric["is_leaf"])),
            bool_text(bool(metric["is_manuscript_panel_like_leaf"])),
        ]
        print("| " + " | ".join(row) + " |", file=output)


def write_tsv(metrics: list[dict[str, object]], output) -> None:
    headers = [
        "rank",
        "id",
        "line",
        "root_distance",
        "is_root",
        "direct_children",
        "descendant_rows",
        "leaf_descendants",
        "manuscript_panel_like_leaf_descendants",
        "is_leaf",
        "is_manuscript_panel_like_leaf",
    ]
    print("\t".join(headers), file=output)
    for rank, metric in enumerate(metrics, start=1):
        values = [
            str(rank),
            str(metric["id"]),
            str(metric["line"]),
            "" if metric["root_distance"] is None else str(metric["root_distance"]),
            bool_text(bool(metric["is_root"])),
            str(metric["direct_children"]),
            str(metric["descendant_rows"]),
            str(metric["leaf_descendants"]),
            str(metric["manuscript_panel_like_leaf_descendants"]),
            bool_text(bool(metric["is_leaf"])),
            bool_text(bool(metric["is_manuscript_panel_like_leaf"])),
        ]
        print("\t".join(values), file=output)


def compile_patterns(patterns: list[str]) -> list[re.Pattern[str]]:
    return [re.compile(pattern, re.IGNORECASE) for pattern in patterns]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank rows in a method-table provenance Markdown table by transitive "
            "downstream manuscript/panel-like leaf support."
        )
    )
    parser.add_argument("table", type=Path, help="Markdown table with id|parent|what|why|comment columns.")
    parser.add_argument("--output", type=Path, help="Optional output file. Defaults to stdout.")
    parser.add_argument(
        "--format",
        choices=["markdown", "tsv"],
        default="markdown",
        help="Output format. Default: markdown.",
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

    rows = parse_method_table(args.table)
    if not rows:
        parser.error(f"no method-table rows parsed from: {args.table}")

    by_id, children, parents, duplicate_ids = build_graph(rows)
    if duplicate_ids:
        shown = ", ".join(sorted(set(duplicate_ids))[:5])
        suffix = "" if len(set(duplicate_ids)) <= 5 else f" (+{len(set(duplicate_ids)) - 5} more)"
        parser.error(f"duplicate row ids are not supported: {shown}{suffix}")

    panel_patterns = compile_patterns(DEFAULT_PANEL_LIKE_PATTERNS + args.panel_like_regex)
    distances = root_distances(by_id, children, parents)
    metrics = make_metrics(rows, by_id, children, parents, distances, panel_patterns)
    root_count = sum(1 for node_id in by_id if not parents.get(node_id))
    panel_like_leaf_count = sum(
        1 for row in rows if is_panel_like_leaf(row, children, panel_patterns)
    )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            if args.format == "markdown":
                write_markdown(metrics, handle, args.table, len(rows), root_count, panel_like_leaf_count)
            else:
                write_tsv(metrics, handle)
    else:
        if args.format == "markdown":
            write_markdown(metrics, sys.stdout, args.table, len(rows), root_count, panel_like_leaf_count)
        else:
            write_tsv(metrics, sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main())
