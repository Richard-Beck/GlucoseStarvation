#!/usr/bin/env python3
"""Validate an agent-authored consequence audit against a manuscript graph."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable


GRAPH_COLUMNS = ("upstream_class", "downstream_class", "relationship")
CLASSIFICATION_COLUMNS = ("artifact_path", "source", "artifact_class", "rationale")
CHANGE_COLUMNS = (
    "artifact_path",
    "kind",
    "change_type",
    "parent_sha256",
    "current_sha256",
    "recorded_at",
)
HASH_REFERENCE_COLUMNS = ("consumer_path", "parent_sha256", "changed_artifacts")
AUDIT_COLUMNS = (
    "upstream_class",
    "downstream_class",
    "relationship",
    "affected_artifacts",
    "decision",
    "rationale",
    "owner",
    "state",
    "resolution",
)
CLASSIFICATION_SOURCES = {"changed_artifact", "hash_backreference_consumer"}
DECISIONS = {"invalidated", "remains_valid", "unresolved"}
STATES = {"open", "resolved", "accepted_exception", "blocked"}
TERMINAL_CLASS = "terminal"


class AuditError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, columns: Iterable[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            missing = [column for column in columns if column not in fields]
            if missing:
                raise AuditError(f"{path} is missing column(s): {', '.join(missing)}")
            return list(reader)
    except OSError as exc:
        raise AuditError(f"Cannot read TSV {path}: {exc}") from exc


def safe_name(value: str, field: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value or ""):
        raise AuditError(f"Invalid {field}: {value!r}")
    return value


def safe_path(value: str, field: str) -> str:
    path = Path(value)
    if path.is_absolute() or not path.parts or path == Path(".") or ".." in path.parts:
        raise AuditError(f"{field} must be a safe assembly-relative path: {value!r}")
    return value


def split_paths(value: str, field: str) -> set[str]:
    return {
        safe_path(item, field)
        for item in (part.strip() for part in value.split(";"))
        if item
    }


def validate_graph(
    rows: list[dict[str, str]],
) -> tuple[dict[str, list[tuple[str, str, str]]], set[str]]:
    adjacency: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    nodes: set[str] = set()
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        upstream = safe_name(row["upstream_class"], "upstream_class")
        downstream = safe_name(row["downstream_class"], "downstream_class")
        relationship = safe_name(row["relationship"], "relationship")
        key = (upstream, downstream, relationship)
        if key in seen:
            raise AuditError(f"Duplicate dependency-graph edge: {key}")
        seen.add(key)
        adjacency[upstream].append(key)
        nodes.update((upstream, downstream))

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str) -> None:
        if node in visiting:
            raise AuditError(f"Dependency graph contains a cycle through {node!r}")
        if node in visited:
            return
        visiting.add(node)
        for _, downstream, _ in adjacency.get(node, []):
            visit(downstream)
        visiting.remove(node)
        visited.add(node)

    for node in sorted(nodes):
        visit(node)
    return adjacency, nodes


def validate_classifications(
    rows: list[dict[str, str]],
    graph_nodes: set[str],
    changed_paths: set[str],
    consumer_paths: set[str],
) -> dict[tuple[str, str], set[str]]:
    expected = {
        *( (path, "changed_artifact") for path in changed_paths ),
        *( (path, "hash_backreference_consumer") for path in consumer_paths ),
    }
    classified: dict[tuple[str, str], set[str]] = defaultdict(set)
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        path = safe_path(row["artifact_path"], "artifact_path")
        source = row["source"]
        if source not in CLASSIFICATION_SOURCES:
            raise AuditError(f"Unsupported classification source {source!r} for {path}")
        artifact_class = safe_name(row["artifact_class"], "artifact_class")
        if artifact_class not in graph_nodes and artifact_class != TERMINAL_CLASS:
            raise AuditError(
                f"Classification for {path!r} uses class {artifact_class!r}, which is absent from the graph"
            )
        if not row["rationale"].strip():
            raise AuditError(f"Classification rationale is required for {(path, source, artifact_class)}")
        key = (path, source, artifact_class)
        if key in seen:
            raise AuditError(f"Duplicate artifact classification: {key}")
        seen.add(key)
        classified[(path, source)].add(artifact_class)

    observed = set(classified)
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        details = []
        if missing:
            details.append("missing=" + ",".join(f"{path}:{source}" for path, source in missing))
        if extra:
            details.append("unexpected=" + ",".join(f"{path}:{source}" for path, source in extra))
        raise AuditError("Artifact classification does not match observed dependency inputs: " + "; ".join(details))
    mixed_terminal = sorted(
        key for key, classes in classified.items() if TERMINAL_CLASS in classes and len(classes) > 1
    )
    if mixed_terminal:
        raise AuditError(
            "The terminal class cannot be combined with downstream classes: "
            + ",".join(f"{path}:{source}" for path, source in mixed_terminal)
        )
    return classified


def validate_audit_files(
    graph_path: Path,
    classifications_path: Path,
    changes_path: Path,
    hash_references_path: Path,
    audit_path: Path,
) -> dict[str, object]:
    graph_rows = read_tsv(graph_path, GRAPH_COLUMNS)
    classification_rows = read_tsv(classifications_path, CLASSIFICATION_COLUMNS)
    change_rows = read_tsv(changes_path, CHANGE_COLUMNS)
    reference_rows = read_tsv(hash_references_path, HASH_REFERENCE_COLUMNS)
    audit_rows = read_tsv(audit_path, AUDIT_COLUMNS)

    adjacency, graph_nodes = validate_graph(graph_rows)
    changed_paths = {row["artifact_path"] for row in change_rows}
    consumer_paths = {row["consumer_path"] for row in reference_rows}
    classified = validate_classifications(
        classification_rows, graph_nodes, changed_paths, consumer_paths
    )

    audit_by_edge: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in audit_rows:
        upstream = safe_name(row["upstream_class"], "audit upstream_class")
        downstream = safe_name(row["downstream_class"], "audit downstream_class")
        relationship = safe_name(row["relationship"], "audit relationship")
        key = (upstream, downstream, relationship)
        if upstream not in graph_nodes or downstream not in graph_nodes:
            raise AuditError(f"Consequence-audit edge uses a class absent from the graph: {key}")
        if key in audit_by_edge:
            raise AuditError(f"Duplicate consequence-audit edge: {key}")
        decision = row["decision"]
        state = row["state"]
        if decision not in DECISIONS:
            raise AuditError(f"Unsupported audit decision {decision!r} for {key}")
        if state not in STATES:
            raise AuditError(f"Unsupported audit state {state!r} for {key}")
        if not row["rationale"].strip():
            raise AuditError(f"Audit rationale is required for {key}")
        if not row["owner"].strip():
            raise AuditError(f"Audit owner is required for {key}")
        safe_name(row["owner"], "audit owner")
        if decision == "remains_valid" and state != "resolved":
            raise AuditError(f"remains_valid must use state='resolved' for {key}")
        if decision == "unresolved" and state != "blocked":
            raise AuditError(f"unresolved must use state='blocked' for {key}")
        if decision == "invalidated" and state != "open" and not row["resolution"].strip():
            raise AuditError(f"Resolved or excepted invalidation requires a resolution for {key}")
        if state in {"accepted_exception", "blocked"} and not row["resolution"].strip():
            raise AuditError(f"State {state!r} requires a resolution for {key}")
        split_paths(row["affected_artifacts"], "affected_artifacts")
        audit_by_edge[key] = row

    classes_by_changed_path = {
        path: classified[(path, "changed_artifact")] - {TERMINAL_CLASS}
        for path in changed_paths
    }
    root_classes = set().union(*classes_by_changed_path.values()) if classes_by_changed_path else set()

    expected_artifacts: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for reference in reference_rows:
        consumer = reference["consumer_path"]
        consumer_classes = classified[(consumer, "hash_backreference_consumer")] - {TERMINAL_CLASS}
        for changed in split_paths(reference["changed_artifacts"], "changed_artifacts"):
            for upstream in classes_by_changed_path.get(changed, set()):
                for downstream in consumer_classes:
                    existing = [edge for edge in adjacency.get(upstream, []) if edge[1] == downstream]
                    keys = existing or [(upstream, downstream, "exact_hash_reference")]
                    for key in keys:
                        expected_artifacts[key].add(consumer)

    required: set[tuple[str, str, str]] = set()

    def traverse(artifact_class: str) -> None:
        for edge in adjacency.get(artifact_class, []):
            if edge in required:
                continue
            required.add(edge)
            row = audit_by_edge.get(edge)
            if row and row["decision"] == "invalidated":
                traverse(edge[1])

    for artifact_class in sorted(root_classes):
        traverse(artifact_class)
    for edge in sorted(expected_artifacts):
        if edge not in required:
            required.add(edge)
            row = audit_by_edge.get(edge)
            if row and row["decision"] == "invalidated":
                traverse(edge[1])

    missing_edges = sorted(required - set(audit_by_edge))
    extra_edges = sorted(set(audit_by_edge) - required)
    errors: list[str] = []
    if missing_edges:
        errors.append("missing=" + ",".join("/".join(edge) for edge in missing_edges))
    if extra_edges:
        errors.append("unreachable=" + ",".join("/".join(edge) for edge in extra_edges))
    for edge, artifacts in expected_artifacts.items():
        row = audit_by_edge.get(edge)
        if row:
            recorded = split_paths(row["affected_artifacts"], "affected_artifacts")
            missing_artifacts = sorted(artifacts - recorded)
            if missing_artifacts:
                errors.append(
                    "missing_exact_hash_consumers=" + "/".join(edge) + ":" + ",".join(missing_artifacts)
                )
    if errors:
        raise AuditError("Consequence audit is incomplete or inconsistent: " + "; ".join(errors))

    counts = {
        decision: sum(row["decision"] == decision for row in audit_rows)
        for decision in sorted(DECISIONS)
    }
    open_count = sum(row["state"] in {"open", "blocked"} for row in audit_rows)
    return {
        "status": "PASS",
        "graph_sha256": sha256_file(graph_path),
        "classifications_sha256": sha256_file(classifications_path),
        "changes_sha256": sha256_file(changes_path),
        "hash_references_sha256": sha256_file(hash_references_path),
        "audit_sha256": sha256_file(audit_path),
        "changed_artifacts": len(changed_paths),
        "hash_backreference_consumers": len(consumer_paths),
        "classified_rows": len(classification_rows),
        "root_classes": sorted(root_classes),
        "audited_edges": len(audit_rows),
        "decisions": counts,
        "open_or_blocked": open_count,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, required=True)
    parser.add_argument("--classifications", type=Path, required=True)
    parser.add_argument("--changes", type=Path, required=True)
    parser.add_argument("--hash-references", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = validate_audit_files(
        args.graph.resolve(),
        args.classifications.resolve(),
        args.changes.resolve(),
        args.hash_references.resolve(),
        args.audit.resolve(),
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AuditError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
