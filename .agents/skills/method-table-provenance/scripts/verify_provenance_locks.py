#!/usr/bin/env python3
"""Verify lock-target SHA256 values in a locked method provenance table.

This script implements the first provenance-table invalidation check: locked
artifacts silently changed, moved, or disappeared on disk. It does not perform
manuscript endpoint scope checks.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import html
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


REQUIRED_COLUMNS = [
    "id",
    "parent",
    "what",
    "why",
    "comment",
    "lock_target",
    "lock_selector",
    "lock_kind",
    "sha256",
    "hash_status",
]

MISSING_VALUES = {"", "na", "n/a", "none", "null"}
HEX_SHA256_RE = re.compile(r"^(?:sha256:)?([0-9a-fA-F]{64})$")


@dataclass(frozen=True)
class LockedRow:
    line_number: int
    row_index: int
    values: dict[str, str]

    @property
    def node_id(self) -> str:
        return strip_single_code_span(self.values["id"])

    @property
    def lock_target(self) -> str:
        return strip_single_code_span(self.values["lock_target"])

    @property
    def lock_kind(self) -> str:
        return strip_single_code_span(self.values["lock_kind"])

    @property
    def stored_sha256_raw(self) -> str:
        return strip_single_code_span(self.values["sha256"])

    @property
    def source_hash_status(self) -> str:
        return strip_single_code_span(self.values["hash_status"])


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


def clean_cell(text: str) -> str:
    text = html.unescape(text)
    text = text.replace("<br>", ";").replace("<br/>", ";").replace("<br />", ";")
    text = text.replace("&nbsp;", " ")
    text = text.replace("\\|", "|")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def split_markdown_row(line: str) -> list[str]:
    line = line.rstrip("\n")
    stripped = line.strip()
    if not (stripped.startswith("|") and stripped.endswith("|")):
        return []
    inner = stripped[1:-1]
    cells: list[str] = []
    current: list[str] = []
    escaped = False
    for char in inner:
        if escaped:
            current.append(char)
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if char == "|":
            cells.append(clean_cell("".join(current)))
            current = []
            continue
        current.append(char)
    if escaped:
        current.append("\\")
    cells.append(clean_cell("".join(current)))
    return cells


def normalize_header(text: str) -> str:
    text = strip_single_code_span(text)
    return re.sub(r"\s+", " ", text).strip().lower()


def strip_single_code_span(text: str) -> str:
    text = clean_cell(text)
    match = re.fullmatch(r"`([^`]+)`", text)
    if match:
        return match.group(1).strip()
    return text.strip()


def is_separator_row(cells: list[str]) -> bool:
    return bool(cells) and all(set(cell.replace(" ", "")) <= {"-", ":"} for cell in cells)


def parse_locked_table(path: Path) -> tuple[list[LockedRow], list[str]]:
    rows: list[LockedRow] = []
    warnings: list[str] = []
    header: list[str] | None = None
    column_indices: dict[str, int] = {}

    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.lstrip().startswith("|"):
            continue
        cells = split_markdown_row(line)
        if not cells:
            continue
        normalized = [normalize_header(cell) for cell in cells]
        if all(column in normalized for column in REQUIRED_COLUMNS):
            header = normalized
            column_indices = {column: normalized.index(column) for column in REQUIRED_COLUMNS}
            continue
        if header is None:
            continue
        if is_separator_row(cells):
            continue
        if len(cells) != len(header):
            warnings.append(
                f"line {line_number}: expected {len(header)} cells from table header, found {len(cells)}"
            )
            continue
        values = {column: cells[index] for column, index in column_indices.items()}
        rows.append(LockedRow(line_number=line_number, row_index=len(rows) + 1, values=values))

    if header is None:
        warnings.append("no locked provenance table header found")
    return rows, warnings


def normalize_missing(text: str) -> str:
    return re.sub(r"\s+", " ", strip_single_code_span(text)).strip().lower()


def normalize_sha256(text: str) -> tuple[str | None, str | None]:
    cleaned = normalize_missing(text)
    if cleaned in MISSING_VALUES:
        return None, None
    match = HEX_SHA256_RE.fullmatch(cleaned)
    if not match:
        return None, "invalid_stored_sha256"
    return match.group(1).lower(), None


def has_glob_meta(text: str) -> bool:
    return any(char in text for char in "*?[")


def resolve_target(raw_target: str, root: Path) -> tuple[Path | None, str]:
    target = strip_single_code_span(raw_target)
    if normalize_missing(target) in MISSING_VALUES:
        return None, "no_lock_target"
    if "<" in target or ">" in target:
        return None, "template_lock_target"

    if has_glob_meta(target):
        pattern = target if Path(target).is_absolute() else str(root / target)
        matches = sorted(Path(match) for match in glob.glob(pattern))
        if not matches:
            return None, "glob_matched_no_files"
        if len(matches) > 1:
            return None, f"glob_matched_{len(matches)}_files"
        return matches[0], "resolved"

    path = Path(target)
    if not path.is_absolute():
        path = root / path
    return path, "resolved"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_row(row: LockedRow, root: Path) -> CheckResult:
    stored_sha256, stored_error = normalize_sha256(row.stored_sha256_raw)
    if stored_error:
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status="not_checked",
            reason=stored_error,
            stored_sha256=row.stored_sha256_raw,
            current_sha256="NA",
            resolved_path="NA",
        )
    if stored_sha256 is None:
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status="not_checked",
            reason="no_stored_sha256",
            stored_sha256="NA",
            current_sha256="NA",
            resolved_path="NA",
        )

    resolved, resolve_status = resolve_target(row.lock_target, root)
    if resolved is None:
        status = "ambiguous" if resolve_status.startswith("glob_matched_") else "unresolved"
        if resolve_status == "glob_matched_no_files":
            status = "missing"
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status=status,
            reason=resolve_status,
            stored_sha256=stored_sha256,
            current_sha256="NA",
            resolved_path="NA",
        )

    if not resolved.exists():
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status="missing",
            reason="lock_target_missing",
            stored_sha256=stored_sha256,
            current_sha256="NA",
            resolved_path=str(resolved),
        )
    if not resolved.is_file():
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status="unresolved",
            reason="lock_target_not_file",
            stored_sha256=stored_sha256,
            current_sha256="NA",
            resolved_path=str(resolved),
        )

    try:
        current_sha256 = sha256_file(resolved)
    except OSError as exc:
        return CheckResult(
            row_index=row.row_index,
            line_number=row.line_number,
            node_id=row.node_id,
            lock_target=row.lock_target,
            lock_kind=row.lock_kind,
            source_hash_status=row.source_hash_status,
            status="unresolved",
            reason=f"hash_error: {exc}",
            stored_sha256=stored_sha256,
            current_sha256="NA",
            resolved_path=str(resolved),
        )

    status = "ok" if current_sha256 == stored_sha256 else "changed"
    reason = "sha256_match" if status == "ok" else "sha256_mismatch"
    return CheckResult(
        row_index=row.row_index,
        line_number=row.line_number,
        node_id=row.node_id,
        lock_target=row.lock_target,
        lock_kind=row.lock_kind,
        source_hash_status=row.source_hash_status,
        status=status,
        reason=reason,
        stored_sha256=stored_sha256,
        current_sha256=current_sha256,
        resolved_path=str(resolved),
    )


def markdown_escape(text: object) -> str:
    return str(text).replace("\n", " ").replace("|", "\\|")


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
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({column: getattr(result, column) for column in columns})


def format_markdown(
    table_path: Path,
    root: Path,
    results: list[CheckResult],
    warnings: list[str],
    max_problems: int,
    show_ok: bool,
) -> str:
    counts = Counter(result.status for result in results)
    lines = [
        "# Provenance Lock Verification",
        "",
        f"- Table: `{table_path}`",
        f"- Root: `{root}`",
        f"- Parsed rows: {len(results)}",
        "",
        "## Status Summary",
        "",
        "| status | n_rows |",
        "|---|---:|",
    ]
    for status in ["ok", "changed", "missing", "ambiguous", "unresolved", "not_checked"]:
        lines.append(f"| {status} | {counts.get(status, 0)} |")

    if warnings:
        lines.extend(["", "## Parse Warnings", ""])
        for warning in warnings:
            lines.append(f"- {markdown_escape(warning)}")

    if show_ok:
        displayed = results
        heading = "## Row Results"
    else:
        displayed = [result for result in results if result.status != "ok"]
        heading = "## Non-OK Rows"

    lines.extend(["", heading, ""])
    if not displayed:
        lines.append("No rows to report.")
        return "\n".join(lines) + "\n"

    if max_problems >= 0:
        truncated_count = max(0, len(displayed) - max_problems)
        displayed = displayed[:max_problems]
    else:
        truncated_count = 0

    lines.extend(
        [
            "| row | line | status | reason | id | lock_target | stored_sha256 | current_sha256 |",
            "|---:|---:|---|---|---|---|---|---|",
        ]
    )
    for result in displayed:
        current = result.current_sha256
        if len(current) > 16 and current != "NA":
            current = current[:12] + "..."
        stored = result.stored_sha256
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
        lines.append("")
        lines.append(f"Truncated {truncated_count} additional row(s); use `--max-problems -1` to show all.")

    return "\n".join(lines) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify SHA256 lock targets in a locked method-table provenance Markdown table."
    )
    parser.add_argument("table", type=Path, help="Locked Markdown provenance table.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Root used to resolve repo-relative lock_target paths. Default: current directory.",
    )
    parser.add_argument("--output", type=Path, help="Optional Markdown report path. Defaults to stdout.")
    parser.add_argument("--details-tsv", type=Path, help="Optional TSV with one row per checked table row.")
    parser.add_argument(
        "--max-problems",
        type=int,
        default=100,
        help="Maximum non-OK rows to include in Markdown output. Use -1 for all. Default: 100.",
    )
    parser.add_argument("--show-ok", action="store_true", help="Include OK rows in the Markdown row table.")
    parser.add_argument(
        "--fail-on-problems",
        action="store_true",
        help="Exit with code 1 if any changed, missing, ambiguous, or unresolved rows are found.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    table_path = args.table
    root = args.root.resolve()

    rows, warnings = parse_locked_table(table_path)
    results = [check_row(row, root) for row in rows]

    report = format_markdown(
        table_path=table_path,
        root=root,
        results=results,
        warnings=warnings,
        max_problems=args.max_problems,
        show_ok=args.show_ok,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
    else:
        sys.stdout.write(report)

    if args.details_tsv:
        args.details_tsv.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(args.details_tsv, results)

    if args.fail_on_problems:
        problem_statuses = {"changed", "missing", "ambiguous", "unresolved"}
        if any(result.status in problem_statuses for result in results):
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
