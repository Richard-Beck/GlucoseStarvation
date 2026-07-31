#!/usr/bin/env python3
"""Recompute SHA256 locks in one canonical method provenance table."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from provenance_table import (
    MISSING_VALUES,
    ProvenanceRow,
    normalize_missing,
    normalize_sha256,
    parse_provenance_table,
)


HASH_REQUIRED_STATUSES = {
    "computed_self",
    "computed_downstream_proxy",
    "computed_upstream_proxy",
    "computed_code_proxy",
    "computed_representative",
    "metadata_checksum",
}
INTENTIONAL_UNHASHED_STATUSES = {
    "not_applicable",
    "external",
    "missing",
    "ambiguous",
    "unresolved",
}


@dataclass(frozen=True)
class CheckResult:
    row_index: int
    line_number: int
    node_id: str
    lock_target: str
    lock_kind: str
    source_hash_status: str
    status: str
    reason: str
    stored_sha256: str
    current_sha256: str
    resolved_path: str


def has_glob_meta(text: str) -> bool:
    return any(character in text for character in "*?[")


def resolve_target(raw_target: str, root: Path) -> tuple[Path | None, str]:
    if normalize_missing(raw_target) in MISSING_VALUES:
        return None, "no_lock_target"
    if "<" in raw_target or ">" in raw_target:
        return None, "template_lock_target"
    if has_glob_meta(raw_target):
        pattern = raw_target if Path(raw_target).is_absolute() else str(root / raw_target)
        matches = sorted(Path(match) for match in glob.glob(pattern))
        if not matches:
            return None, "glob_matched_no_files"
        if len(matches) > 1:
            return None, f"glob_matched_{len(matches)}_files"
        return matches[0], "resolved"
    path = Path(raw_target)
    return (path if path.is_absolute() else root / path), "resolved"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def result_for(
    row: ProvenanceRow,
    *,
    status: str,
    reason: str,
    stored_sha256: str = "NA",
    current_sha256: str = "NA",
    resolved_path: str = "NA",
) -> CheckResult:
    return CheckResult(
        row_index=row.row_index,
        line_number=row.line_number,
        node_id=row.node_id,
        lock_target=row.value("lock_target"),
        lock_kind=row.value("lock_kind"),
        source_hash_status=row.value("hash_status"),
        status=status,
        reason=reason,
        stored_sha256=stored_sha256,
        current_sha256=current_sha256,
        resolved_path=resolved_path,
    )


def check_row(row: ProvenanceRow, root: Path) -> CheckResult:
    raw_digest = row.value("sha256")
    digest = normalize_sha256(raw_digest)
    hash_status = row.value("hash_status").lower()
    if digest is None:
        if normalize_missing(raw_digest) not in MISSING_VALUES:
            return result_for(row, status="invalid", reason="invalid_stored_sha256")
        if hash_status in HASH_REQUIRED_STATUSES:
            return result_for(row, status="invalid", reason="required_sha256_missing")
        if hash_status in INTENTIONAL_UNHASHED_STATUSES:
            return result_for(row, status="not_applicable", reason=hash_status)
        return result_for(row, status="invalid", reason="invalid_hash_status_without_sha256")

    resolved, resolve_status = resolve_target(row.value("lock_target"), root)
    if resolved is None:
        status = "ambiguous" if resolve_status.startswith("glob_matched_") else "unresolved"
        if resolve_status == "glob_matched_no_files":
            status = "missing"
        return result_for(
            row,
            status=status,
            reason=resolve_status,
            stored_sha256=digest,
        )
    if not resolved.exists():
        return result_for(
            row,
            status="missing",
            reason="lock_target_missing",
            stored_sha256=digest,
            resolved_path=str(resolved),
        )
    if not resolved.is_file():
        return result_for(
            row,
            status="unresolved",
            reason="lock_target_not_file",
            stored_sha256=digest,
            resolved_path=str(resolved),
        )
    try:
        current_sha256 = sha256_file(resolved)
    except OSError as exc:
        return result_for(
            row,
            status="unresolved",
            reason=f"hash_error: {exc}",
            stored_sha256=digest,
            resolved_path=str(resolved),
        )
    status = "ok" if current_sha256 == digest else "changed"
    return result_for(
        row,
        status=status,
        reason="sha256_match" if status == "ok" else "sha256_mismatch",
        stored_sha256=digest,
        current_sha256=current_sha256,
        resolved_path=str(resolved),
    )


def markdown_escape(value: object) -> str:
    return str(value).replace("\n", " ").replace("|", "\\|")


def write_tsv(path: Path, results: list[CheckResult]) -> None:
    columns = [
        "row_index",
        "line_number",
        "node_id",
        "lock_target",
        "lock_kind",
        "source_hash_status",
        "status",
        "reason",
        "stored_sha256",
        "current_sha256",
        "resolved_path",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({column: getattr(result, column) for column in columns})


def format_markdown(
    table_path: Path,
    root: Path,
    results: list[CheckResult],
    errors: list[str],
    max_problems: int,
    show_ok: bool,
) -> str:
    problem_statuses = {"changed", "missing", "ambiguous", "unresolved", "invalid"}
    counts = Counter(result.status for result in results)
    passed = bool(results) and not errors and not any(
        result.status in problem_statuses for result in results
    )
    lines = [
        "# Provenance Lock Verification",
        "",
        f"- Status: `{'PASS' if passed else 'FAIL'}`",
        f"- Table: `{table_path}`",
        f"- Root: `{root}`",
        f"- Parsed rows: {len(results)}",
        "",
        "## Status Summary",
        "",
        "| status | n_rows |",
        "|---|---:|",
    ]
    for status in [
        "ok",
        "changed",
        "missing",
        "ambiguous",
        "unresolved",
        "invalid",
        "not_applicable",
    ]:
        lines.append(f"| {status} | {counts.get(status, 0)} |")

    if errors:
        lines.extend(["", "## Validation Errors", ""])
        lines.extend(f"- {markdown_escape(error)}" for error in errors)

    displayed = results if show_ok else [
        result for result in results if result.status not in {"ok", "not_applicable"}
    ]
    lines.extend(["", "## Row Results" if show_ok else "## Problem Rows", ""])
    if not displayed:
        lines.append("No rows to report.")
        return "\n".join(lines) + "\n"

    truncated_count = 0
    if max_problems >= 0:
        truncated_count = max(0, len(displayed) - max_problems)
        displayed = displayed[:max_problems]
    lines.extend(
        [
            "| row | line | status | reason | id | lock_target | stored_sha256 | current_sha256 |",
            "|---:|---:|---|---|---|---|---|---|",
        ]
    )
    for result in displayed:
        current = result.current_sha256
        stored = result.stored_sha256
        if len(current) > 16 and current != "NA":
            current = current[:12] + "..."
        if len(stored) > 16 and stored != "NA":
            stored = stored[:12] + "..."
        lines.append(
            "| {row} | {line} | {status} | {reason} | `{node}` | `{target}` | `{stored}` | `{current}` |".format(
                row=result.row_index,
                line=result.line_number,
                status=markdown_escape(result.status),
                reason=markdown_escape(result.reason),
                node=markdown_escape(result.node_id),
                target=markdown_escape(result.lock_target),
                stored=markdown_escape(stored),
                current=markdown_escape(current),
            )
        )
    if truncated_count:
        lines.extend(
            [
                "",
                f"Truncated {truncated_count} additional row(s); use `--max-problems -1` to show all.",
            ]
        )
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify SHA256 targets in one canonical method provenance table."
    )
    parser.add_argument("table", type=Path, help="Canonical Markdown provenance table.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Root used to resolve repo-relative lock_target paths. Default: current directory.",
    )
    parser.add_argument("--output", type=Path, help="Optional Markdown report path.")
    parser.add_argument("--details-tsv", type=Path, help="Optional TSV with one result per row.")
    parser.add_argument(
        "--max-problems",
        type=int,
        default=100,
        help="Maximum problem rows in Markdown. Use -1 for all. Default: 100.",
    )
    parser.add_argument("--show-ok", action="store_true", help="Include all rows in the report.")
    parser.add_argument(
        "--fail-on-problems",
        action="store_true",
        help="Exit 1 for parse errors or changed, missing, ambiguous, unresolved, or invalid rows.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = args.root.resolve()
    rows, errors = parse_provenance_table(args.table)
    results = [check_row(row, root) for row in rows]
    report = format_markdown(
        args.table,
        root,
        results,
        errors,
        args.max_problems,
        args.show_ok,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
    else:
        sys.stdout.write(report)
    if args.details_tsv:
        write_tsv(args.details_tsv, results)

    if args.fail_on_problems:
        problem_statuses = {"changed", "missing", "ambiguous", "unresolved", "invalid"}
        if errors or not results or any(
            result.status in problem_statuses for result in results
        ):
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
