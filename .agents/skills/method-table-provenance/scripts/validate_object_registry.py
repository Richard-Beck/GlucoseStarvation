#!/usr/bin/env python3
"""Validate object and dependency registries used during provenance closure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path


OBJECT_COLUMNS = [
    "object_id",
    "object_type",
    "canonical_path",
    "selector",
    "sha256",
    "resolution_status",
    "transformation",
    "methods_role",
    "evidence",
]
EDGE_COLUMNS = [
    "child_id",
    "parent_id",
    "relation",
    "dependency_status",
    "evidence",
    "comment",
]
OBJECT_TYPES = {
    "panel_endpoint",
    "script",
    "function",
    "action",
    "named_object",
    "data_file",
    "config",
    "raw_source",
    "external",
    "unresolved",
}
RESOLUTION_STATUSES = {"queued", "resolved", "terminal", "unresolved"}
DEPENDENCY_STATUSES = {
    "candidate",
    "confirmed",
    "unresolved",
}
LOCAL_FILE_TYPES = {"panel_endpoint", "script", "data_file", "config"}
MISSING_VALUES = {"", "na", "n/a", "none", "null"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def read_tsv(path: Path, expected_columns: list[str]) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as exc:
        return [], [f"cannot read {path}: {exc}"]

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_columns:
            return [], [
                f"{path} must have exactly these tab-separated columns: "
                + " | ".join(expected_columns)
            ]
        rows = [
            {column: (value or "").strip() for column, value in row.items()}
            for row in reader
        ]
    if not rows:
        errors.append(f"{path} contains no data rows")
    return rows, errors


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detect_cycles(parents: dict[str, set[str]]) -> list[str]:
    errors: list[str] = []
    state: dict[str, int] = {}
    stack: list[str] = []

    def visit(node_id: str) -> None:
        if state.get(node_id) == 2:
            return
        if state.get(node_id) == 1:
            start = stack.index(node_id)
            cycle = stack[start:] + [node_id]
            errors.append("confirmed dependency cycle: " + " -> ".join(cycle))
            return
        state[node_id] = 1
        stack.append(node_id)
        for parent_id in sorted(parents.get(node_id, set())):
            visit(parent_id)
        stack.pop()
        state[node_id] = 2

    for node_id in sorted(parents):
        visit(node_id)
    return errors


def validate(
    registry_path: Path,
    edge_path: Path,
    root: Path,
    require_closed: bool,
) -> tuple[list[str], dict[str, int]]:
    object_rows, object_errors = read_tsv(registry_path, OBJECT_COLUMNS)
    edge_rows, edge_errors = read_tsv(edge_path, EDGE_COLUMNS)
    errors = object_errors + edge_errors

    object_ids = [row["object_id"] for row in object_rows]
    object_counts = Counter(object_ids)
    for object_id, count in sorted(object_counts.items()):
        if not object_id:
            errors.append("object_registry.tsv contains an empty object_id")
        elif count > 1:
            errors.append(f"duplicate object_id: {object_id}")
    known_ids = {object_id for object_id in object_ids if object_id}

    for line_number, row in enumerate(object_rows, start=2):
        prefix = f"object row {line_number} ({row['object_id'] or 'empty id'})"
        if row["object_type"] not in OBJECT_TYPES:
            errors.append(f"{prefix}: invalid object_type {row['object_type']!r}")
        if row["resolution_status"] not in RESOLUTION_STATUSES:
            errors.append(
                f"{prefix}: invalid resolution_status {row['resolution_status']!r}"
            )
        if require_closed and row["resolution_status"] == "queued":
            errors.append(f"{prefix}: queued object remains in a closed registry")
        if not row["transformation"]:
            errors.append(f"{prefix}: transformation is empty")
        if not row["methods_role"]:
            errors.append(f"{prefix}: methods_role is empty")
        if not row["evidence"]:
            errors.append(f"{prefix}: evidence is empty")

        canonical_path = row["canonical_path"]
        path_missing = canonical_path.lower() in MISSING_VALUES
        if not path_missing and Path(canonical_path).is_absolute():
            errors.append(f"{prefix}: canonical_path must be repo-relative")
        if row["object_type"] in LOCAL_FILE_TYPES and path_missing:
            errors.append(f"{prefix}: local file object requires canonical_path")

        digest = row["sha256"].lower()
        digest_missing = digest in MISSING_VALUES
        if not digest_missing and not SHA256_RE.fullmatch(digest):
            errors.append(f"{prefix}: invalid sha256")
        if not digest_missing:
            if path_missing:
                errors.append(f"{prefix}: hashed object has no canonical_path")
            else:
                resolved = root / canonical_path
                if not resolved.is_file():
                    errors.append(f"{prefix}: hashed canonical_path is not a file")
                else:
                    current = file_sha256(resolved)
                    if current != digest:
                        errors.append(
                            f"{prefix}: sha256 mismatch; stored={digest}, current={current}"
                        )

    edge_keys: list[tuple[str, str, str, str]] = []
    confirmed_parents: dict[str, set[str]] = defaultdict(set)
    for line_number, row in enumerate(edge_rows, start=2):
        prefix = f"edge row {line_number} ({row['child_id']} -> {row['parent_id']})"
        for column in ("child_id", "parent_id", "relation", "dependency_status"):
            if not row[column]:
                errors.append(f"{prefix}: {column} is empty")
        if row["child_id"] and row["child_id"] not in known_ids:
            errors.append(f"{prefix}: unknown child_id")
        if row["parent_id"] and row["parent_id"] not in known_ids:
            errors.append(f"{prefix}: unknown parent_id")
        if row["child_id"] and row["child_id"] == row["parent_id"]:
            errors.append(f"{prefix}: self-dependency is not allowed")
        if row["dependency_status"] not in DEPENDENCY_STATUSES:
            errors.append(
                f"{prefix}: invalid dependency_status {row['dependency_status']!r}"
            )
        if require_closed and row["dependency_status"] == "candidate":
            errors.append(f"{prefix}: candidate edge remains in a closed registry")
        if row["dependency_status"] in {"confirmed", "unresolved"} and not row["evidence"]:
            errors.append(f"{prefix}: evidence is required")

        key = (
            row["child_id"],
            row["parent_id"],
            row["relation"],
            row["dependency_status"],
        )
        edge_keys.append(key)
        if (
            row["dependency_status"] == "confirmed"
            and row["child_id"] in known_ids
            and row["parent_id"] in known_ids
        ):
            confirmed_parents[row["child_id"]].add(row["parent_id"])

    for key, count in Counter(edge_keys).items():
        if count > 1:
            errors.append("duplicate dependency edge: " + " | ".join(key))
    errors.extend(detect_cycles(confirmed_parents))

    summary = {
        "objects": len(object_rows),
        "edges": len(edge_rows),
        "queued_objects": sum(
            row["resolution_status"] == "queued" for row in object_rows
        ),
        "candidate_edges": sum(
            row["dependency_status"] == "candidate" for row in edge_rows
        ),
        "errors": len(errors),
    }
    return errors, summary


def write_report(
    path: Path | None,
    registry_path: Path,
    edge_path: Path,
    require_closed: bool,
    errors: list[str],
    summary: dict[str, int],
) -> None:
    status = "PASS" if not errors else "FAIL"
    lines = [
        "# Object Registry Validation",
        "",
        f"- Status: `{status}`",
        f"- Object registry: `{registry_path}`",
        f"- Dependency edges: `{edge_path}`",
        f"- Require closed: `{str(require_closed).lower()}`",
        f"- Objects: {summary['objects']}",
        f"- Edges: {summary['edges']}",
        f"- Queued objects: {summary['queued_objects']}",
        f"- Candidate edges: {summary['candidate_edges']}",
        "",
        "## Errors",
        "",
    ]
    lines.extend(f"- {error}" for error in errors)
    if not errors:
        lines.append("- None.")
    text = "\n".join(lines) + "\n"
    if path is None:
        print(text, end="")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("object_registry", type=Path)
    parser.add_argument("dependency_edges", type=Path)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-closed", action="store_true")
    args = parser.parse_args(argv)

    errors, summary = validate(
        args.object_registry,
        args.dependency_edges,
        args.root.resolve(),
        args.require_closed,
    )
    write_report(
        args.output,
        args.object_registry,
        args.dependency_edges,
        args.require_closed,
        errors,
        summary,
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
