#!/usr/bin/env python3
"""Validate clean-room figure audit coverage and table contracts."""

import argparse
import json
import re
import sys
from pathlib import Path


HEADING = "### Observation–claim relationships"
ASCII_HEADING = "### Observation-claim relationships"
NO_RELATIONSHIP = "No material relationships identified."
RELATIONS = {"support", "undermine"}
STRENGTHS = {"strong", "moderate", "weak"}
EXPECTED_HEADER = [
    "Panel(s)",
    "Exact observation",
    "Claim",
    "Relation",
    "Strength",
    "Reason",
]


def read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"could not read JSON from {path}: {exc}") from exc


def clean_cell(value: str) -> str:
    return value.strip().strip("`").strip()


def parse_row(line: str):
    if not line.strip().startswith("|") or not line.strip().endswith("|"):
        return None
    return [cell.strip() for cell in line.strip()[1:-1].split("|")]


def claims_by_id(payload, errors):
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        errors.append("claims index schema_version must be 1")
        return {}
    claims = payload.get("claims")
    if not isinstance(claims, list):
        errors.append("claims index claims must be a list")
        return {}
    result = {}
    for claim in claims:
        if not isinstance(claim, dict):
            errors.append("claims index entries must be objects")
            continue
        claim_id = str(claim.get("id", "")).strip()
        text = str(claim.get("text", "")).strip()
        if not claim_id or not text:
            errors.append("claims index entries require id and text")
        elif claim_id in result:
            errors.append(f"duplicate authoritative claim id: {claim_id}")
        else:
            result[claim_id] = text
    return result


def figure_ids(payload, errors):
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        errors.append("figure index schema_version must be 1")
        return set()
    figures = payload.get("figures")
    if not isinstance(figures, list):
        errors.append("figure index figures must be a list")
        return set()
    result = set()
    for figure in figures:
        if not isinstance(figure, dict):
            errors.append("figure index entries must be objects")
            continue
        figure_id = str(figure.get("figure_id", "")).strip()
        if not figure_id:
            errors.append("figure index entries require figure_id")
        elif figure_id in result:
            errors.append(f"duplicate figure_id: {figure_id}")
        else:
            result.add(figure_id)
    return result


def audit_paths(payload, repo_root, errors):
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        errors.append("audit index schema_version must be 1")
        return {}
    audits = payload.get("audits")
    if not isinstance(audits, list):
        errors.append("audit index audits must be a list")
        return {}
    result = {}
    for audit in audits:
        if not isinstance(audit, dict):
            errors.append("audit index entries must be objects")
            continue
        figure_id = str(audit.get("figure_id", "")).strip()
        path_text = str(audit.get("path", "")).strip()
        if not figure_id or not path_text:
            errors.append("audit index entries require figure_id and path")
            continue
        if figure_id in result:
            errors.append(f"duplicate audit figure_id: {figure_id}")
            continue
        path = Path(path_text)
        if not path.is_absolute():
            path = repo_root / path
        result[figure_id] = path.resolve()
    return result


def validate_audit(figure_id, path, claims, errors):
    if not path.is_file():
        errors.append(f"{figure_id}: audit file is missing: {path}")
        return 0
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    try:
        heading_index = next(
            i
            for i, line in enumerate(lines)
            if line.strip() in {HEADING, ASCII_HEADING}
        )
    except StopIteration:
        errors.append(f"{figure_id}: missing observation-claim heading")
        return 0

    table_index = None
    for i in range(heading_index + 1, len(lines)):
        row = parse_row(lines[i])
        if row is not None:
            table_index = i
            break
        if lines[i].strip() and lines[i].strip() != NO_RELATIONSHIP:
            continue
    if table_index is None:
        errors.append(f"{figure_id}: missing observation-claim table")
        return 0

    header = parse_row(lines[table_index])
    if header != EXPECTED_HEADER:
        errors.append(
            f"{figure_id}: table header must be {' | '.join(EXPECTED_HEADER)}"
        )
        return 0
    if table_index + 1 >= len(lines):
        errors.append(f"{figure_id}: table separator is missing")
        return 0
    separator = parse_row(lines[table_index + 1])
    if (
        separator is None
        or len(separator) != len(EXPECTED_HEADER)
        or any(not re.fullmatch(r":?-{3,}:?", cell) for cell in separator)
    ):
        errors.append(f"{figure_id}: table separator is invalid")
        return 0

    row_count = 0
    for line in lines[table_index + 2 :]:
        row = parse_row(line)
        if row is None:
            if line.strip():
                break
            continue
        if len(row) != len(EXPECTED_HEADER):
            errors.append(f"{figure_id}: table row has {len(row)} cells, expected 6")
            continue
        cleaned = [
            clean_cell(cell) for cell in row
        ]
        if not any(cleaned):
            continue
        panel, observation, claim_cell, relation, strength, reason = cleaned
        row_count += 1
        if not panel:
            errors.append(f"{figure_id}: row {row_count} has no panel")
        if not observation:
            errors.append(f"{figure_id}: row {row_count} has no observation")
        if relation not in RELATIONS:
            errors.append(
                f"{figure_id}: row {row_count} has invalid relation {relation!r}"
            )
        if strength not in STRENGTHS:
            errors.append(
                f"{figure_id}: row {row_count} has invalid strength {strength!r}"
            )
        if not reason:
            errors.append(f"{figure_id}: row {row_count} has no reason")
        match = re.fullmatch(r"([^:]+):\s*(.+)", claim_cell)
        if not match:
            errors.append(
                f"{figure_id}: row {row_count} claim must be 'ID: exact text'"
            )
            continue
        claim_id = match.group(1).strip()
        claim_text = match.group(2).strip()
        if claim_id not in claims:
            errors.append(
                f"{figure_id}: row {row_count} references unknown claim {claim_id}"
            )
        elif claim_text != claims[claim_id]:
            errors.append(
                f"{figure_id}: row {row_count} changes authoritative text for "
                f"{claim_id}"
            )

    if row_count == 0 and NO_RELATIONSHIP not in text:
        errors.append(
            f"{figure_id}: empty table must include '{NO_RELATIONSHIP}'"
        )
    return row_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-index", required=True, type=Path)
    parser.add_argument("--claims-index", required=True, type=Path)
    parser.add_argument("--audit-index", required=True, type=Path)
    parser.add_argument("--repo-root", default=".", type=Path)
    args = parser.parse_args()

    errors = []
    repo_root = args.repo_root.resolve()
    try:
        figures = figure_ids(read_json(args.figure_index), errors)
        claims = claims_by_id(read_json(args.claims_index), errors)
        audits = audit_paths(read_json(args.audit_index), repo_root, errors)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    missing = sorted(figures.difference(audits))
    extra = sorted(set(audits).difference(figures))
    if missing:
        errors.append(f"missing figure audits: {', '.join(missing)}")
    if extra:
        errors.append(f"audits for unknown figures: {', '.join(extra)}")

    total_rows = 0
    for figure_id in sorted(figures.intersection(audits)):
        total_rows += validate_audit(
            figure_id,
            audits[figure_id],
            claims,
            errors,
        )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        f"OK: {len(figures)} figure audits, "
        f"{len(claims)} authoritative claims, "
        f"{total_rows} observation-claim relationships"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
