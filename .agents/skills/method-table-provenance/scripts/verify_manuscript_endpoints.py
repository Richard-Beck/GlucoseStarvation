#!/usr/bin/env python3
"""Validate a current panel-leaf provenance graph against a target manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from provenance_table import (
    MISSING_VALUES,
    ProvenanceRow,
    normalize_missing,
    normalize_sha256,
    parse_parent_ids,
    parse_provenance_table,
)


TARGET_COLUMNS = ["endpoint_id", "artifact_path", "lock_selector", "sha256"]
LOCK_KINDS = {
    "file",
    "code",
    "proxy_file",
    "representative_file",
    "terminal",
    "external",
    "unresolved",
}
HASHED_STATUS_KINDS = {
    "computed_self": {"file", "code"},
    "computed_downstream_proxy": {"proxy_file"},
    "computed_upstream_proxy": {"proxy_file"},
    "computed_code_proxy": {"code"},
    "computed_representative": {"representative_file"},
    "metadata_checksum": {"file", "code", "proxy_file", "representative_file"},
}
UNHASHED_STATUSES = {
    "not_applicable",
    "external",
    "missing",
    "ambiguous",
    "unresolved",
}
SELECTOR_RE = re.compile(r"^#panel_[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class TargetEndpoint:
    endpoint_id: str
    artifact_path: str
    lock_selector: str
    declared_sha256: str
    current_sha256: str
    artifact_status: str


@dataclass(frozen=True)
class EndpointResult:
    endpoint_id: str
    artifact_path: str
    graph_lock_target: str
    lock_selector: str
    graph_lock_selector: str
    current_sha256: str
    graph_sha256: str
    status: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_target_manifest(path: Path, root: Path) -> tuple[list[TargetEndpoint], list[str]]:
    errors: list[str] = []
    endpoints: list[TargetEndpoint] = []
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        return [], [f"cannot read target figure-set manifest: {exc}"]

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if columns != TARGET_COLUMNS:
            return [], [
                "target figure-set manifest must have exactly these tab-separated "
                f"columns: {' | '.join(TARGET_COLUMNS)}"
            ]

        for line_number, row in enumerate(reader, start=2):
            endpoint_id = (row.get("endpoint_id") or "").strip()
            artifact_path = (row.get("artifact_path") or "").strip()
            lock_selector = (row.get("lock_selector") or "").strip()
            declared_raw = (row.get("sha256") or "").strip()
            declared_sha256 = normalize_sha256(declared_raw)

            if not endpoint_id:
                errors.append(f"target manifest line {line_number}: empty endpoint_id")
            if not artifact_path:
                errors.append(f"target manifest line {line_number}: empty artifact_path")
            elif Path(artifact_path).is_absolute():
                errors.append(
                    f"target manifest line {line_number}: artifact_path must be repo-relative"
                )
            if not SELECTOR_RE.fullmatch(lock_selector):
                errors.append(
                    f"target manifest line {line_number}: lock_selector must match #panel_<id>"
                )
            if endpoint_id and lock_selector and not endpoint_id.endswith(lock_selector):
                errors.append(
                    f"target manifest line {line_number}: endpoint_id must end with lock_selector"
                )
            if declared_sha256 is None:
                errors.append(f"target manifest line {line_number}: invalid sha256")

            resolved = root / artifact_path if artifact_path else root
            if not artifact_path or not resolved.exists():
                current_sha256 = "NA"
                artifact_status = "current_artifact_missing"
            elif not resolved.is_file():
                current_sha256 = "NA"
                artifact_status = "current_artifact_not_file"
            else:
                try:
                    current_sha256 = sha256_file(resolved)
                    artifact_status = (
                        "ok"
                        if declared_sha256 == current_sha256
                        else "target_manifest_hash_mismatch"
                    )
                except OSError:
                    current_sha256 = "NA"
                    artifact_status = "current_artifact_unreadable"

            endpoints.append(
                TargetEndpoint(
                    endpoint_id=endpoint_id,
                    artifact_path=artifact_path,
                    lock_selector=lock_selector,
                    declared_sha256=declared_sha256 or declared_raw or "NA",
                    current_sha256=current_sha256,
                    artifact_status=artifact_status,
                )
            )

    if not endpoints:
        errors.append("target figure-set manifest contains no panel endpoints")

    endpoint_counts = Counter(endpoint.endpoint_id for endpoint in endpoints)
    for endpoint_id, count in endpoint_counts.items():
        if endpoint_id and count > 1:
            errors.append(f"duplicate target endpoint_id: {endpoint_id}")

    pair_counts = Counter(
        (endpoint.artifact_path, endpoint.lock_selector) for endpoint in endpoints
    )
    for (artifact_path, selector), count in pair_counts.items():
        if artifact_path and selector and count > 1:
            errors.append(
                f"duplicate target artifact/selector pair: {artifact_path} {selector}"
            )

    path_hashes: dict[str, set[str]] = defaultdict(set)
    hash_paths: dict[str, set[str]] = defaultdict(set)
    for endpoint in endpoints:
        if endpoint.artifact_path and normalize_sha256(endpoint.declared_sha256):
            path_hashes[endpoint.artifact_path].add(endpoint.declared_sha256)
            hash_paths[endpoint.declared_sha256].add(endpoint.artifact_path)
    for artifact_path, hashes in path_hashes.items():
        if len(hashes) != 1:
            errors.append(f"target artifact has conflicting hashes: {artifact_path}")
    for digest, paths in hash_paths.items():
        if len(paths) != 1:
            errors.append(
                f"target hash identifies multiple final-figure paths: {digest}"
            )
    if path_hashes and len(path_hashes) != len(hash_paths):
        errors.append(
            "target manifest must have one unique whole-figure hash per unique "
            "final-figure path"
        )
    return endpoints, errors


def validate_lock_fields(row: ProvenanceRow) -> list[str]:
    errors: list[str] = []
    prefix = f"line {row.line_number} ({row.node_id or 'empty id'})"
    target = row.value("lock_target")
    selector = row.value("lock_selector")
    lock_kind = row.value("lock_kind").lower()
    hash_status = row.value("hash_status").lower()
    digest = normalize_sha256(row.value("sha256"))
    target_missing = normalize_missing(target) in MISSING_VALUES
    selector_missing = normalize_missing(selector) in MISSING_VALUES

    if lock_kind not in LOCK_KINDS:
        errors.append(f"{prefix}: invalid lock_kind {lock_kind!r}")
    if hash_status not in set(HASHED_STATUS_KINDS) | UNHASHED_STATUSES:
        errors.append(f"{prefix}: invalid hash_status {hash_status!r}")

    if digest is not None:
        if target_missing:
            errors.append(f"{prefix}: hashed row has no lock_target")
        elif Path(target).is_absolute():
            errors.append(f"{prefix}: hashed lock_target must be repo-relative")
        elif any(character in target for character in "*?[]<>"):
            errors.append(f"{prefix}: hashed lock_target must identify one canonical file")
        allowed_kinds = HASHED_STATUS_KINDS.get(hash_status)
        if allowed_kinds is None:
            errors.append(f"{prefix}: stored sha256 is incompatible with {hash_status}")
        elif lock_kind not in allowed_kinds:
            errors.append(
                f"{prefix}: lock_kind {lock_kind!r} is incompatible with {hash_status}"
            )
    else:
        raw_digest = normalize_missing(row.value("sha256"))
        if raw_digest not in MISSING_VALUES:
            errors.append(f"{prefix}: malformed sha256")
        if hash_status in HASHED_STATUS_KINDS:
            errors.append(f"{prefix}: {hash_status} requires a sha256")
        if hash_status == "not_applicable" and lock_kind not in {
            "terminal",
            "external",
            "unresolved",
        }:
            errors.append(
                f"{prefix}: not_applicable requires terminal, external, or unresolved lock_kind"
            )
        if not target_missing and hash_status == "not_applicable":
            errors.append(f"{prefix}: not_applicable row should not name a lock_target")

    if not selector_missing and selector.startswith("#panel_") and not SELECTOR_RE.fullmatch(
        selector
    ):
        errors.append(f"{prefix}: malformed panel selector")
    return errors


def validate_graph(
    rows: list[ProvenanceRow], target_ids: set[str]
) -> tuple[dict[str, ProvenanceRow], dict[str, list[str]], list[str]]:
    errors: list[str] = []
    rows_by_id: dict[str, ProvenanceRow] = {}
    duplicate_ids: set[str] = set()
    parents_by_id: dict[str, list[str]] = {}

    for row in rows:
        errors.extend(validate_lock_fields(row))
        node_id = row.node_id
        if not node_id:
            errors.append(f"line {row.line_number}: empty id")
            continue
        if node_id in rows_by_id:
            duplicate_ids.add(node_id)
        else:
            rows_by_id[node_id] = row
        parent_ids = parse_parent_ids(row.values["parent"])
        repeated_parent_ids = {
            parent_id
            for parent_id, count in Counter(parent_ids).items()
            if count > 1
        }
        for parent_id in sorted(repeated_parent_ids):
            errors.append(
                f"duplicate parent reference for {node_id}: {parent_id}"
            )
        parents_by_id[node_id] = parent_ids

    for node_id in sorted(duplicate_ids):
        errors.append(f"duplicate provenance id: {node_id}")

    for node_id, parent_ids in parents_by_id.items():
        for parent_id in parent_ids:
            if parent_id == node_id:
                errors.append(f"self-parent edge: {node_id}")
            elif parent_id not in rows_by_id:
                errors.append(f"missing parent for {node_id}: {parent_id}")

    colors: dict[str, int] = {}
    stack: list[str] = []
    reported_cycles: set[tuple[str, ...]] = set()

    def visit(node_id: str) -> None:
        colors[node_id] = 1
        stack.append(node_id)
        for parent_id in parents_by_id.get(node_id, []):
            if parent_id not in rows_by_id:
                continue
            if colors.get(parent_id, 0) == 0:
                visit(parent_id)
            elif colors.get(parent_id) == 1:
                start = stack.index(parent_id)
                cycle = tuple(stack[start:] + [parent_id])
                if cycle not in reported_cycles:
                    reported_cycles.add(cycle)
                    errors.append(f"provenance cycle: {' -> '.join(cycle)}")
        stack.pop()
        colors[node_id] = 2

    for node_id in rows_by_id:
        if colors.get(node_id, 0) == 0:
            visit(node_id)

    referenced_parents = {
        parent_id for parent_ids in parents_by_id.values() for parent_id in parent_ids
    }
    leaves = {
        node_id: row for node_id, row in rows_by_id.items() if node_id not in referenced_parents
    }
    leaf_ids = set(leaves)
    for node_id in sorted(target_ids - leaf_ids):
        errors.append(f"current panel endpoint is not a graph leaf: {node_id}")
    for node_id in sorted(leaf_ids - target_ids):
        errors.append(f"graph leaf is absent from current panel scope: {node_id}")

    reachable: set[str] = set()
    pending = list(target_ids & leaf_ids)
    while pending:
        node_id = pending.pop()
        if node_id in reachable:
            continue
        reachable.add(node_id)
        pending.extend(
            parent_id
            for parent_id in parents_by_id.get(node_id, [])
            if parent_id in rows_by_id
        )
    for node_id in sorted(set(rows_by_id) - reachable):
        errors.append(f"row is not upstream-reachable from a current panel leaf: {node_id}")

    target_hashes: dict[str, set[str]] = defaultdict(set)
    direct_locks: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows_by_id.values():
        target = row.value("lock_target")
        selector = row.value("lock_selector")
        digest = normalize_sha256(row.value("sha256"))
        status = row.value("hash_status").lower()
        if digest is not None and normalize_missing(target) not in MISSING_VALUES:
            target_hashes[target].add(digest)
        if status == "computed_self" and digest is not None:
            normalized_selector = "" if normalize_missing(selector) in MISSING_VALUES else selector
            direct_locks[(target, normalized_selector)].append(row.node_id)

    for target, hashes in target_hashes.items():
        if len(hashes) > 1:
            errors.append(f"lock_target has conflicting sha256 values: {target}")
    for (target, selector), node_ids in direct_locks.items():
        if len(node_ids) > 1:
            selector_text = selector or "NA"
            errors.append(
                "duplicate direct-file object lock requires reconciliation: "
                f"{target} {selector_text} ({', '.join(sorted(node_ids))})"
            )
    return leaves, parents_by_id, errors


def compare_endpoints(
    targets: list[TargetEndpoint],
    leaves: dict[str, ProvenanceRow],
    parents_by_id: dict[str, list[str]],
) -> list[EndpointResult]:
    targets_by_id = {target.endpoint_id: target for target in targets if target.endpoint_id}
    results: list[EndpointResult] = []
    for endpoint_id in sorted(set(targets_by_id) | set(leaves)):
        target = targets_by_id.get(endpoint_id)
        leaf = leaves.get(endpoint_id)
        status = "endpoint_match"
        if target is None:
            status = "graph_leaf_not_in_current_scope"
        elif leaf is None:
            status = "current_endpoint_missing_from_graph_leaves"
        else:
            graph_target = leaf.value("lock_target")
            graph_selector = leaf.value("lock_selector")
            graph_sha = normalize_sha256(leaf.value("sha256"))
            if target.artifact_status != "ok":
                status = target.artifact_status
            elif graph_target != target.artifact_path:
                status = "endpoint_path_mismatch"
            elif graph_selector != target.lock_selector:
                status = "endpoint_selector_mismatch"
            elif graph_sha != target.current_sha256:
                status = "endpoint_hash_mismatch"
            elif leaf.value("lock_kind").lower() != "file":
                status = "endpoint_lock_kind_not_file"
            elif leaf.value("hash_status").lower() != "computed_self":
                status = "endpoint_hash_status_not_computed_self"
            elif not parents_by_id.get(endpoint_id):
                status = "endpoint_has_no_methods_parent"

        results.append(
            EndpointResult(
                endpoint_id=endpoint_id,
                artifact_path=target.artifact_path if target else "NA",
                graph_lock_target=leaf.value("lock_target") if leaf else "NA",
                lock_selector=target.lock_selector if target else "NA",
                graph_lock_selector=leaf.value("lock_selector") if leaf else "NA",
                current_sha256=target.current_sha256 if target else "NA",
                graph_sha256=(
                    normalize_sha256(leaf.value("sha256")) or leaf.value("sha256")
                    if leaf
                    else "NA"
                ),
                status=status,
            )
        )
    return results


def markdown_escape(value: object) -> str:
    return str(value).replace("\n", " ").replace("|", "\\|")


def format_report(
    manifest: Path,
    table: Path,
    root: Path,
    row_count: int,
    targets: list[TargetEndpoint],
    leaves: dict[str, ProvenanceRow],
    results: list[EndpointResult],
    errors: list[str],
) -> str:
    passed = bool(results) and not errors and all(
        result.status == "endpoint_match" for result in results
    )
    counts = Counter(result.status for result in results)
    target_hashes = {
        target.declared_sha256
        for target in targets
        if normalize_sha256(target.declared_sha256) is not None
    }
    leaf_hashes = {
        normalize_sha256(row.value("sha256"))
        for row in leaves.values()
        if normalize_sha256(row.value("sha256")) is not None
    }
    lines = [
        "# Manuscript Panel-Leaf Provenance Verification",
        "",
        f"- Status: `{'PASS' if passed else 'FAIL'}`",
        f"- Target figure-set manifest: `{manifest}`",
        f"- Canonical provenance table: `{table}`",
        f"- Root: `{root}`",
        f"- Provenance rows: {row_count}",
        f"- Target panel endpoints: {len(targets)}",
        f"- Graph leaves: {len(leaves)}",
        f"- Target final-figure hashes: {len(target_hashes)}",
        f"- Graph-leaf hashes: {len(leaf_hashes)}",
        "",
        "## Status Summary",
        "",
        "| status | n_endpoints |",
        "|---|---:|",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"| {markdown_escape(status)} | {count} |")

    if errors:
        lines.extend(["", "## Validation Errors", ""])
        lines.extend(f"- {markdown_escape(error)}" for error in errors)

    lines.extend(
        [
            "",
            "## Panel Endpoint Results",
            "",
            "| endpoint_id | artifact_path | graph_lock_target | lock_selector | "
            "graph_lock_selector | current_sha256 | graph_sha256 | status |",
            "|---|---|---|---|---|---|---|---|",
        ]
    )
    for result in results:
        lines.append(
            "| `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | {} |".format(
                markdown_escape(result.endpoint_id),
                markdown_escape(result.artifact_path),
                markdown_escape(result.graph_lock_target),
                markdown_escape(result.lock_selector),
                markdown_escape(result.graph_lock_selector),
                markdown_escape(result.current_sha256),
                markdown_escape(result.graph_sha256),
                markdown_escape(result.status),
            )
        )
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate one current canonical provenance graph and compare its "
            "structural leaves with a panel-level target figure-set manifest."
        ),
        epilog=(
            "Required target TSV columns: endpoint_id, artifact_path, "
            "lock_selector, sha256. Provenance endpoints are derived from graph leaves."
        ),
    )
    parser.add_argument("target_manifest", type=Path)
    parser.add_argument("locked_provenance_table", type=Path)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Root used to resolve repo-relative artifact paths. Default: current directory.",
    )
    parser.add_argument("--output", type=Path, help="Optional Markdown report path.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = args.root.resolve()
    targets, target_errors = parse_target_manifest(args.target_manifest, root)
    rows, table_errors = parse_provenance_table(args.locked_provenance_table)
    target_ids = {target.endpoint_id for target in targets if target.endpoint_id}
    leaves, parents_by_id, graph_errors = validate_graph(rows, target_ids)
    results = compare_endpoints(targets, leaves, parents_by_id)

    errors = target_errors + table_errors + graph_errors
    target_figure_count = len(
        {
            target.declared_sha256
            for target in targets
            if normalize_sha256(target.declared_sha256) is not None
        }
    )
    leaf_figure_count = len(
        {
            normalize_sha256(row.value("sha256"))
            for row in leaves.values()
            if normalize_sha256(row.value("sha256")) is not None
        }
    )
    if leaf_figure_count != target_figure_count:
        errors.append(
            "unique graph-leaf hash count does not match current final-figure count: "
            f"{leaf_figure_count} != {target_figure_count}"
        )

    report = format_report(
        args.target_manifest,
        args.locked_provenance_table,
        root,
        len(rows),
        targets,
        leaves,
        results,
        errors,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
    else:
        sys.stdout.write(report)

    return 0 if results and not errors and all(
        result.status == "endpoint_match" for result in results
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
