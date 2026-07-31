#!/usr/bin/env python3
"""Audit or apply the repository-local junk relocation policy."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import sys


ROOT = Path(__file__).resolve().parents[1]
KEEP_FILE = ROOT / "docs/current_manuscript_rebuild_keep_roots.txt"
PROVENANCE_LOCATIONS = (
    ROOT
    / "agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction"
    / "provenance_object_locations.tsv"
)
FIGURE_REBUILD_MANIFEST = (
    ROOT
    / "manuscripts/20260731T101031_v03/rebuild/figures"
    / "figure_rebuild_manifest.tsv"
)
EXCLUDED_ROOTS = {Path(".git"), Path("junk")}


def read_keep_roots() -> tuple[Path, ...]:
    roots: list[Path] = []
    for number, raw in enumerate(KEEP_FILE.read_text().splitlines(), start=1):
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or path in EXCLUDED_ROOTS:
            raise ValueError(f"unsafe keep path on line {number}: {value}")
        roots.append(path)
    if len(roots) != len(set(roots)):
        raise ValueError("duplicate paths in keep-root manifest")
    return tuple(sorted(roots, key=lambda path: path.parts))


def is_kept_or_ancestor(path: Path, keep_roots: tuple[Path, ...]) -> bool:
    return any(path == keep or path in keep.parents for keep in keep_roots)


def is_covered(path: Path, keep_roots: tuple[Path, ...]) -> bool:
    return any(path == keep or keep in path.parents for keep in keep_roots)


def verify_dependency_coverage(keep_roots: tuple[Path, ...]) -> tuple[int, int]:
    with PROVENANCE_LOCATIONS.open(newline="") as handle:
        provenance_rows = list(csv.DictReader(handle, delimiter="\t"))
    provenance_paths = {
        Path(row["canonical_path"])
        for row in provenance_rows
        if row["canonical_path"] != "NA"
    }
    uncovered_provenance = sorted(
        path for path in provenance_paths if not is_covered(path, keep_roots)
    )

    with FIGURE_REBUILD_MANIFEST.open(newline="") as handle:
        figure_rows = list(csv.DictReader(handle, delimiter="\t"))
    figure_inputs = {
        Path(value)
        for row in figure_rows
        for value in row["direct_inputs"].split(";")
        if value and not Path(value).is_absolute()
    }
    uncovered_figure_inputs = sorted(
        path for path in figure_inputs if not is_covered(path, keep_roots)
    )

    if uncovered_provenance or uncovered_figure_inputs:
        if uncovered_provenance:
            print("Uncovered provenance paths:", file=sys.stderr)
            for path in uncovered_provenance:
                print(f"  {path}", file=sys.stderr)
        if uncovered_figure_inputs:
            print("Uncovered declared figure inputs:", file=sys.stderr)
            for path in uncovered_figure_inputs:
                print(f"  {path}", file=sys.stderr)
        raise ValueError("keep-root manifest does not cover current dependencies")
    return len(provenance_paths), len(figure_inputs)


def collect_moves(directory: Path, keep_roots: tuple[Path, ...]) -> list[Path]:
    moves: list[Path] = []
    for entry in sorted(directory.iterdir(), key=lambda path: path.name):
        relative = entry.relative_to(ROOT)
        if relative in EXCLUDED_ROOTS:
            continue
        if relative in keep_roots:
            continue
        if entry.is_dir() and not entry.is_symlink() and is_kept_or_ancestor(relative, keep_roots):
            moves.extend(collect_moves(entry, keep_roots))
        else:
            moves.append(relative)
    return moves


def validate_keep_roots(keep_roots: tuple[Path, ...]) -> list[Path]:
    missing = [path for path in keep_roots if not (ROOT / path).exists()]
    nested = [
        path
        for path in keep_roots
        if any(parent in keep_roots for parent in path.parents if parent != Path("."))
    ]
    if nested:
        rendered = ", ".join(str(path) for path in nested)
        raise ValueError(f"redundant keep roots nested under kept directories: {rendered}")
    return missing


def move_to_junk(relative: Path) -> None:
    source = ROOT / relative
    destination = ROOT / "junk" / relative
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"junk destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--apply",
        action="store_true",
        help="move unexpected paths under junk/; default is a read-only audit",
    )
    args = parser.parse_args()

    keep_roots = read_keep_roots()
    missing = validate_keep_roots(keep_roots)
    if missing:
        print("Missing keep roots:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 2

    provenance_count, figure_input_count = verify_dependency_coverage(keep_roots)

    moves = collect_moves(ROOT, keep_roots)
    for relative in moves:
        print(relative)
    print(
        f"SUMMARY keep_roots={len(keep_roots)} move_roots={len(moves)} "
        f"provenance_paths={provenance_count} figure_inputs={figure_input_count} "
        f"mode={'apply' if args.apply else 'audit'}"
    )

    if args.apply:
        for relative in moves:
            move_to_junk(relative)
        remaining = collect_moves(ROOT, keep_roots)
        if remaining:
            print("Unexpected paths remain after relocation:", file=sys.stderr)
            for path in remaining:
                print(f"  {path}", file=sys.stderr)
            return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
