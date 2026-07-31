#!/usr/bin/env python3
"""Verify a declared manuscript figure set against locked provenance endpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


TARGET_COLUMNS = ["endpoint_id", "artifact_path", "sha256"]
LOCK_COLUMNS = ["endpoint_id", "lock_target", "sha256"]
SHA256_RE = re.compile(r"^(?:sha256:)?([0-9a-fA-F]{64})$")


@dataclass(frozen=True)
class TargetEndpoint:
    endpoint_id: str
    artifact_path: str
    declared_sha256: str
    current_sha256: str
    artifact_status: str


@dataclass(frozen=True)
class LockedEndpoint:
    endpoint_id: str
    lock_target: str
    sha256: str


@dataclass(frozen=True)
class EndpointResult:
    endpoint_id: str
    artifact_path: str
    lock_target: str
    declared_sha256: str
    current_sha256: str
    locked_sha256: str
    status: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_sha256(value: str) -> str | None:
    match = SHA256_RE.fullmatch(value.strip())
    return match.group(1).lower() if match else None


def clean_cell(value: str) -> str:
    value = html.unescape(value).strip()
    match = re.fullmatch(r"`([^`]+)`", value)
    return match.group(1).strip() if match else value


def split_markdown_row(line: str) -> list[str]:
    stripped = line.strip()
    if not (stripped.startswith("|") and stripped.endswith("|")):
        return []
    cells: list[str] = []
    current: list[str] = []
    escaped = False
    for character in stripped[1:-1]:
        if escaped:
            current.append(character)
            escaped = False
        elif character == "\\":
            escaped = True
        elif character == "|":
            cells.append(clean_cell("".join(current)))
            current = []
        else:
            current.append(character)
    if escaped:
        current.append("\\")
    cells.append(clean_cell("".join(current)))
    return cells


def is_separator(cells: list[str]) -> bool:
    return bool(cells) and all(
        cell and set(cell.replace(" ", "")) <= {"-", ":"} for cell in cells
    )


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
            declared_raw = (row.get("sha256") or "").strip()
            declared_sha256 = normalize_sha256(declared_raw)
            if not endpoint_id:
                errors.append(f"target manifest line {line_number}: empty endpoint_id")
            if not artifact_path:
                errors.append(f"target manifest line {line_number}: empty artifact_path")
            if declared_sha256 is None:
                errors.append(f"target manifest line {line_number}: invalid sha256")

            resolved = Path(artifact_path)
            if not resolved.is_absolute():
                resolved = root / resolved
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
                    declared_sha256=declared_sha256 or declared_raw or "NA",
                    current_sha256=current_sha256,
                    artifact_status=artifact_status,
                )
            )

    if not endpoints:
        errors.append("target figure-set manifest contains no endpoints")
    return endpoints, errors


def parse_locked_endpoints(path: Path) -> tuple[list[LockedEndpoint], list[str]]:
    errors: list[str] = []
    endpoints: list[LockedEndpoint] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        return [], [f"cannot read locked provenance table: {exc}"]

    header_count = 0
    in_table = False
    for line_number, line in enumerate(lines, start=1):
        cells = split_markdown_row(line)
        if cells == LOCK_COLUMNS:
            header_count += 1
            in_table = True
            continue
        if not in_table:
            continue
        if not cells:
            in_table = False
            continue
        if is_separator(cells):
            continue
        if len(cells) != len(LOCK_COLUMNS):
            errors.append(
                f"locked endpoint table line {line_number}: expected "
                f"{len(LOCK_COLUMNS)} cells, found {len(cells)}"
            )
            continue
        endpoint_id, lock_target, stored_raw = cells
        stored_sha256 = normalize_sha256(stored_raw)
        if not endpoint_id:
            errors.append(f"locked endpoint table line {line_number}: empty endpoint_id")
        if not lock_target:
            errors.append(f"locked endpoint table line {line_number}: empty lock_target")
        if stored_sha256 is None:
            errors.append(f"locked endpoint table line {line_number}: invalid sha256")
        endpoints.append(
            LockedEndpoint(
                endpoint_id=endpoint_id,
                lock_target=lock_target,
                sha256=stored_sha256 or stored_raw or "NA",
            )
        )

    if header_count != 1:
        errors.append(
            "locked provenance table must contain exactly one Markdown table with "
            f"this exact header: {' | '.join(LOCK_COLUMNS)}"
        )
    if header_count == 1 and not endpoints:
        errors.append("locked manuscript endpoint table contains no endpoints")
    return endpoints, errors


def compare_endpoints(
    targets: list[TargetEndpoint], locks: list[LockedEndpoint]
) -> list[EndpointResult]:
    target_groups: dict[str, list[TargetEndpoint]] = defaultdict(list)
    lock_groups: dict[str, list[LockedEndpoint]] = defaultdict(list)
    for target in targets:
        target_groups[target.endpoint_id].append(target)
    for lock in locks:
        lock_groups[lock.endpoint_id].append(lock)

    results: list[EndpointResult] = []
    for endpoint_id in sorted(set(target_groups) | set(lock_groups)):
        target_rows = target_groups[endpoint_id]
        lock_rows = lock_groups[endpoint_id]
        target = target_rows[0] if target_rows else None
        lock = lock_rows[0] if lock_rows else None

        if not target_rows:
            status = "table_endpoint_missing_from_current_scope"
        elif not lock_rows:
            status = "current_endpoint_missing_from_table"
        elif len(target_rows) != 1:
            status = "duplicate_target_endpoint"
        elif len(lock_rows) != 1:
            status = "duplicate_locked_endpoint"
        elif target.artifact_status != "ok":
            status = target.artifact_status
        elif normalize_sha256(lock.sha256) is None:
            status = "invalid_locked_sha256"
        elif target.current_sha256 != lock.sha256:
            status = "endpoint_hash_changed"
        else:
            status = "endpoint_hash_match"

        results.append(
            EndpointResult(
                endpoint_id=endpoint_id,
                artifact_path=target.artifact_path if target else "NA",
                lock_target=lock.lock_target if lock else "NA",
                declared_sha256=target.declared_sha256 if target else "NA",
                current_sha256=target.current_sha256 if target else "NA",
                locked_sha256=lock.sha256 if lock else "NA",
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
    results: list[EndpointResult],
    errors: list[str],
) -> str:
    passed = bool(results) and not errors and all(
        result.status == "endpoint_hash_match" for result in results
    )
    counts = Counter(result.status for result in results)
    lines = [
        "# Manuscript Endpoint Verification",
        "",
        f"- Status: `{'PASS' if passed else 'FAIL'}`",
        f"- Target figure-set manifest: `{manifest}`",
        f"- Locked provenance table: `{table}`",
        f"- Root: `{root}`",
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
            "## Endpoint Results",
            "",
            "| endpoint_id | artifact_path | lock_target | declared_sha256 | "
            "current_sha256 | locked_sha256 | status |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for result in results:
        lines.append(
            "| `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | {} |".format(
                markdown_escape(result.endpoint_id),
                markdown_escape(result.artifact_path),
                markdown_escape(result.lock_target),
                markdown_escape(result.declared_sha256),
                markdown_escape(result.current_sha256),
                markdown_escape(result.locked_sha256),
                markdown_escape(result.status),
            )
        )
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify a Methods-local target figure-set manifest against the "
            "manuscript-facing endpoint table in locked provenance."
        ),
        epilog=(
            "Required target TSV columns: endpoint_id, artifact_path, sha256. "
            "Required locked Markdown table header: endpoint_id | lock_target | sha256."
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
    locks, lock_errors = parse_locked_endpoints(args.locked_provenance_table)
    results = compare_endpoints(targets, locks)
    errors = target_errors + lock_errors
    report = format_report(
        args.target_manifest,
        args.locked_provenance_table,
        root,
        results,
        errors,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
    else:
        sys.stdout.write(report)

    return 0 if results and not errors and all(
        result.status == "endpoint_hash_match" for result in results
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
