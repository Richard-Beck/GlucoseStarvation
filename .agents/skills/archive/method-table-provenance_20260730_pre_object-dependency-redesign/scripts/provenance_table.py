#!/usr/bin/env python3
"""Shared parser and schema constants for canonical provenance tables."""

from __future__ import annotations

import html
import re
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

LEGACY_ENDPOINT_COLUMNS = ["endpoint_id", "lock_target", "sha256"]
MISSING_VALUES = {"", "na", "n/a", "none", "null"}
SHA256_RE = re.compile(r"^(?:sha256:)?([0-9a-fA-F]{64})$")


@dataclass(frozen=True)
class ProvenanceRow:
    line_number: int
    row_index: int
    values: dict[str, str]

    @property
    def node_id(self) -> str:
        return strip_single_code_span(self.values["id"])

    def value(self, column: str) -> str:
        return strip_single_code_span(self.values[column])


def clean_cell(text: str) -> str:
    text = html.unescape(text)
    text = text.replace("<br>", ";").replace("<br/>", ";").replace("<br />", ";")
    text = text.replace("&nbsp;", " ").replace("\\|", "|")
    return re.sub(r"\s+", " ", text).strip()


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


def strip_single_code_span(text: str) -> str:
    cleaned = clean_cell(text)
    match = re.fullmatch(r"`([^`]+)`", cleaned)
    return match.group(1).strip() if match else cleaned


def normalize_header(text: str) -> str:
    return re.sub(r"\s+", " ", strip_single_code_span(text)).strip().lower()


def normalize_missing(text: str) -> str:
    return re.sub(r"\s+", " ", strip_single_code_span(text)).strip().lower()


def normalize_sha256(text: str) -> str | None:
    cleaned = normalize_missing(text)
    if cleaned in MISSING_VALUES:
        return None
    match = SHA256_RE.fullmatch(cleaned)
    return match.group(1).lower() if match else None


def is_separator_row(cells: list[str]) -> bool:
    return bool(cells) and all(
        cell and set(cell.replace(" ", "")) <= {"-", ":"} for cell in cells
    )


def parse_parent_ids(text: str) -> list[str]:
    if normalize_missing(text) in MISSING_VALUES | {"terminal"}:
        return []
    return [
        strip_single_code_span(part)
        for part in re.split(r"\s*;\s*", clean_cell(text))
        if strip_single_code_span(part)
    ]


def parse_provenance_table(path: Path) -> tuple[list[ProvenanceRow], list[str]]:
    """Parse exactly one canonical 10-column Markdown provenance table."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        return [], [f"cannot read provenance table: {exc}"]

    rows: list[ProvenanceRow] = []
    errors: list[str] = []
    header_lines: list[int] = []
    legacy_header_lines: list[int] = []
    active = False
    parsed_first_table = False

    for line_number, line in enumerate(lines, start=1):
        cells = split_markdown_row(line)
        normalized = [normalize_header(cell) for cell in cells]

        if normalized == REQUIRED_COLUMNS:
            header_lines.append(line_number)
            active = not parsed_first_table
            parsed_first_table = True
            continue
        if normalized == LEGACY_ENDPOINT_COLUMNS:
            legacy_header_lines.append(line_number)
        if cells and set(REQUIRED_COLUMNS).issubset(normalized):
            errors.append(
                f"line {line_number}: provenance header must contain exactly these columns "
                f"in order: {' | '.join(REQUIRED_COLUMNS)}"
            )

        if not active:
            continue
        if not cells:
            active = False
            continue
        if is_separator_row(cells):
            continue
        if len(cells) != len(REQUIRED_COLUMNS):
            errors.append(
                f"line {line_number}: expected {len(REQUIRED_COLUMNS)} provenance cells, "
                f"found {len(cells)}"
            )
            continue
        rows.append(
            ProvenanceRow(
                line_number=line_number,
                row_index=len(rows) + 1,
                values=dict(zip(REQUIRED_COLUMNS, cells)),
            )
        )

    if len(header_lines) != 1:
        errors.append(
            "provenance document must contain exactly one canonical 10-column table; "
            f"found {len(header_lines)}"
        )
    if legacy_header_lines:
        joined = ", ".join(str(line) for line in legacy_header_lines)
        errors.append(
            "legacy three-column endpoint table is not allowed; derive endpoints from "
            f"canonical graph leaves (header line(s): {joined})"
        )
    if len(header_lines) == 1 and not rows:
        errors.append("canonical provenance table contains no rows")
    return rows, errors
