#!/usr/bin/env python3
"""Score overlap between per-subpanel method provenance chains.

The script reads Markdown tables with columns:

    id | parent | what | why | comment

and ranks chain pairs by weighted feature overlap. The score is intended to
prioritize pairwise LLM reconciliation tasks, not to decide node merges.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


TABLE_HEADER = ["id", "parent", "what", "why", "comment"]

PATH_RE = re.compile(
    r"(?<![A-Za-z0-9_./:-])"
    r"([A-Za-z0-9_./:*\-]+\."
    r"(?:csv|tsv|rds|rda|r|rmd|py|json|txt|md|png|tif|tiff|pkl|stan|sh|sbatch|xlsx))"
    r"(?![A-Za-z0-9_./:-])",
    re.IGNORECASE,
)
FUNC_RE = re.compile(r"\b([A-Za-z][A-Za-z0-9_.]*::[A-Za-z][A-Za-z0-9_.]*|[A-Za-z][A-Za-z0-9_.]*\(\))")
SHA_RE = re.compile(r"\bsha256[:= ]*([0-9a-f]{16,64})\b", re.IGNORECASE)
CODE_RE = re.compile(r"`([^`]+)`")

TERMINAL_VALUES = {
    "",
    "na",
    "n/a",
    "terminal",
    "na terminal",
    "none",
    "null",
    "fixed input",
}

LOW_VALUE_PARTS = (
    "final_images/",
    "/final_images/",
    "subpanels/figure_",
    "/subpanels/figure_",
    "layout/layout_plan.csv",
    "layout/subpanel_dimensions.csv",
    "figure_set_manifest.csv",
    "legend.md",
    "validation_report.json",
    "polishing/provenance.csv",
    "polishing/scripts/polish_figures.r",
)

LOW_VALUE_BASENAMES = {
    "figure_1.png",
    "figure_s1.png",
    "figure_s2.png",
    "layout_plan.csv",
    "subpanel_dimensions.csv",
    "validation_report.json",
    "provenance.csv",
    "manifest.csv",
    "legend.md",
}

HIGH_VALUE_EXTENSIONS = {
    ".csv",
    ".rds",
    ".rda",
    ".pkl",
    ".tif",
    ".tiff",
    ".stan",
    ".py",
    ".r",
    ".rmd",
}


@dataclass(frozen=True)
class Row:
    source_file: Path
    row_index: int
    id_raw: str
    parent_raw: str
    what: str
    why: str
    comment: str


@dataclass
class Chain:
    name: str
    path: Path
    rows: list[Row]
    features: Counter[str]
    node_ids: set[str]


def clean_cell(text: str) -> str:
    text = html.unescape(text)
    text = text.replace("<br>", ";").replace("<br/>", ";").replace("<br />", ";")
    text = text.replace("&nbsp;", " ")
    text = text.replace("\\|", "|")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def strip_markdown(text: str) -> str:
    text = clean_cell(text)
    text = text.replace("`", "")
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    return text.strip()


def normalize_token(text: str) -> str:
    text = strip_markdown(text).strip().strip("\"'")
    text = text.replace("\\", "/")
    text = text.replace("::panel_", "#panel_")
    text = re.sub(r"\s+", " ", text)
    text = re.sub(r"^\./", "", text)
    text = text.strip(" .;,")
    text = text.lower()
    text = text.replace("na (terminal)", "terminal")
    text = text.replace("`", "")
    return text


def is_terminal(text: str) -> bool:
    norm = normalize_token(text)
    return norm in TERMINAL_VALUES or norm.startswith("terminal:")


def is_low_value_path(token: str) -> bool:
    norm = normalize_token(token)
    base = Path(norm).name
    return base in LOW_VALUE_BASENAMES or any(part in norm for part in LOW_VALUE_PARTS)


def path_weight(token: str, base_weight: float) -> float:
    norm = normalize_token(token)
    suffix = Path(norm).suffix.lower()
    if is_low_value_path(norm):
        return min(base_weight, 0.35)
    if suffix in HIGH_VALUE_EXTENSIONS:
        return base_weight
    return max(0.5, base_weight * 0.7)


def add_feature(features: Counter[str], key: str, weight: float) -> None:
    key = normalize_token(key)
    if not key or is_terminal(key):
        return
    features[key] += weight


def split_markdown_row(line: str) -> list[str]:
    line = line.strip()
    if not (line.startswith("|") and line.endswith("|")):
        return []
    return [clean_cell(cell) for cell in line.strip("|").split("|")]


def parse_markdown_table(path: Path) -> list[Row]:
    rows: list[Row] = []
    saw_header = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            continue
        cells = split_markdown_row(line)
        if len(cells) != 5:
            continue
        lower = [normalize_token(cell) for cell in cells]
        if lower == TABLE_HEADER:
            saw_header = True
            continue
        if not saw_header:
            continue
        if all(set(cell.replace(" ", "")) <= {"-"} for cell in cells):
            continue
        rows.append(
            Row(
                source_file=path,
                row_index=len(rows) + 1,
                id_raw=cells[0],
                parent_raw=cells[1],
                what=cells[2],
                why=cells[3],
                comment=cells[4],
            )
        )
    return rows


def code_spans(text: str) -> list[str]:
    return [match.group(1) for match in CODE_RE.finditer(text)]


def path_tokens(text: str) -> list[str]:
    tokens = [match.group(1) for match in PATH_RE.finditer(strip_markdown(text))]
    return tokens


def function_tokens(text: str) -> list[str]:
    return [match.group(1) for match in FUNC_RE.finditer(strip_markdown(text))]


def sha_tokens(text: str) -> list[str]:
    return [f"sha256:{match.group(1).lower()}" for match in SHA_RE.finditer(strip_markdown(text))]


def split_parent_field(text: str) -> list[str]:
    raw = strip_markdown(text)
    parts = re.split(r"\s*;\s*|\s+<br\s*/?>\s*", raw)
    return [part.strip() for part in parts if part.strip()]


def row_id_candidates(row: Row) -> list[str]:
    candidates = [strip_markdown(row.id_raw)]
    candidates.extend(code_spans(row.id_raw))
    candidates.extend(path_tokens(row.id_raw))
    return candidates


def row_parent_candidates(row: Row) -> list[str]:
    candidates = split_parent_field(row.parent_raw)
    candidates.extend(code_spans(row.parent_raw))
    candidates.extend(path_tokens(row.parent_raw))
    return candidates


def basename_feature(path_token: str) -> str:
    return f"basename:{Path(normalize_token(path_token)).name}"


def build_chain(path: Path) -> Chain:
    rows = parse_markdown_table(path)
    features: Counter[str] = Counter()
    node_ids: set[str] = set()

    for row in rows:
        for candidate in row_id_candidates(row):
            norm = normalize_token(candidate)
            if not norm or is_terminal(norm):
                continue
            node_ids.add(norm)
            add_feature(features, f"id:{norm}", path_weight(norm, 5.0))
            if Path(norm).suffix:
                add_feature(features, f"path:{norm}", path_weight(norm, 4.0))
                add_feature(features, basename_feature(norm), path_weight(norm, 1.4))

        for candidate in row_parent_candidates(row):
            norm = normalize_token(candidate)
            if not norm or is_terminal(norm):
                continue
            add_feature(features, f"parent:{norm}", path_weight(norm, 3.0))
            if Path(norm).suffix:
                add_feature(features, f"path:{norm}", path_weight(norm, 2.5))
                add_feature(features, basename_feature(norm), path_weight(norm, 1.0))

        text = " ".join([row.what, row.why, row.comment])
        for path_token in path_tokens(text):
            norm = normalize_token(path_token)
            add_feature(features, f"path:{norm}", path_weight(norm, 1.2))
            add_feature(features, basename_feature(norm), path_weight(norm, 0.45))

        for func in function_tokens(text):
            add_feature(features, f"func:{func}", 0.9)

        for sha in sha_tokens(text):
            add_feature(features, sha, 2.0)

    return Chain(name=path.stem, path=path, rows=rows, features=features, node_ids=node_ids)


def weighted_jaccard(a: Counter[str], b: Counter[str]) -> tuple[float, float, float]:
    keys = set(a) | set(b)
    if not keys:
        return 0.0, 0.0, 0.0
    shared = sum(min(a[key], b[key]) for key in keys)
    union = sum(max(a[key], b[key]) for key in keys)
    if union <= 0:
        return 0.0, shared, union
    return shared / union, shared, union


def exact_node_matches(a: Chain, b: Chain) -> set[str]:
    matches = a.node_ids & b.node_ids
    return {match for match in matches if not is_low_value_path(match)}


def top_shared_features(a: Chain, b: Chain, limit: int) -> list[tuple[str, float]]:
    shared = []
    for key in set(a.features) & set(b.features):
        weight = min(a.features[key], b.features[key])
        if weight <= 0:
            continue
        shared.append((key, weight))
    shared.sort(key=lambda item: (-item[1], item[0]))
    return shared[:limit]


def pair_scores(chains: list[Chain]) -> list[dict[str, object]]:
    pairs: list[dict[str, object]] = []
    for i, left in enumerate(chains):
        for j in range(i + 1, len(chains)):
            right = chains[j]
            score, shared, union = weighted_jaccard(left.features, right.features)
            matches = exact_node_matches(left, right)
            pairs.append(
                {
                    "left": left.name,
                    "right": right.name,
                    "score": score,
                    "shared_weight": shared,
                    "union_weight": union,
                    "exact_node_matches": len(matches),
                    "matched_nodes": sorted(matches),
                    "top_features": top_shared_features(left, right, 8),
                }
            )
    pairs.sort(
        key=lambda item: (
            -float(item["score"]),
            -int(item["exact_node_matches"]),
            str(item["left"]),
            str(item["right"]),
        )
    )
    return pairs


def greedy_batch(pairs: list[dict[str, object]], batch_size: int) -> list[dict[str, object]]:
    selected = []
    used: set[str] = set()
    for pair in pairs:
        left = str(pair["left"])
        right = str(pair["right"])
        if left in used or right in used:
            continue
        selected.append(pair)
        used.add(left)
        used.add(right)
        if len(selected) >= batch_size:
            break
    return selected


def write_pairs_csv(path: Path, pairs: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left",
                "right",
                "score",
                "shared_weight",
                "union_weight",
                "exact_node_matches",
                "matched_nodes",
                "top_features",
            ],
        )
        writer.writeheader()
        for pair in pairs:
            row = dict(pair)
            row["score"] = f"{float(row['score']):.6f}"
            row["shared_weight"] = f"{float(row['shared_weight']):.3f}"
            row["union_weight"] = f"{float(row['union_weight']):.3f}"
            row["matched_nodes"] = json.dumps(row["matched_nodes"])
            row["top_features"] = json.dumps(row["top_features"])
            writer.writerow(row)


def write_matrix_csv(path: Path, chains: list[Chain], pairs: list[dict[str, object]]) -> None:
    lookup: dict[tuple[str, str], float] = {}
    for pair in pairs:
        left = str(pair["left"])
        right = str(pair["right"])
        score = float(pair["score"])
        lookup[(left, right)] = score
        lookup[(right, left)] = score
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chain", *[chain.name for chain in chains]])
        for left in chains:
            row = [left.name]
            for right in chains:
                if left.name == right.name:
                    row.append("1.000000")
                else:
                    row.append(f"{lookup.get((left.name, right.name), 0.0):.6f}")
            writer.writerow(row)


def format_features(features: list[tuple[str, float]]) -> str:
    if not features:
        return "-"
    return "; ".join(f"{name} ({weight:.1f})" for name, weight in features)


def print_pairs(pairs: list[dict[str, object]], limit: int) -> None:
    for rank, pair in enumerate(pairs[:limit], start=1):
        print(
            f"{rank:>2}. {pair['left']} <-> {pair['right']}  "
            f"score={float(pair['score']):.3f}  "
            f"shared={float(pair['shared_weight']):.1f}  "
            f"exact_nodes={pair['exact_node_matches']}"
        )
        matched_nodes = pair["matched_nodes"]
        if matched_nodes:
            shown = "; ".join(str(item) for item in list(matched_nodes)[:5])
            if len(matched_nodes) > 5:
                shown += f"; ... (+{len(matched_nodes) - 5})"
            print(f"    matched nodes: {shown}")
        print(f"    top overlap: {format_features(pair['top_features'])}")


def load_chains(input_dir: Path, pattern: str) -> list[Chain]:
    files = sorted(path for path in input_dir.glob(pattern) if path.is_file())
    chains = [build_chain(path) for path in files]
    return [chain for chain in chains if chain.rows]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Score overlap between method-provenance panel chains."
    )
    parser.add_argument("input_dir", type=Path, help="Directory containing per-panel Markdown files.")
    parser.add_argument("--glob", default="*.md", help="Input glob relative to input_dir. Default: *.md")
    parser.add_argument("--top", type=positive_int, default=20, help="Number of ranked pairs to print.")
    parser.add_argument(
        "--batch-size",
        type=positive_int,
        default=5,
        help="Number of disjoint top pairs to print as a parallel reconciliation batch.",
    )
    parser.add_argument("--pairs-csv", type=Path, help="Optional output CSV for all pair scores.")
    parser.add_argument("--matrix-csv", type=Path, help="Optional output CSV for the NxN score matrix.")
    args = parser.parse_args(argv)

    if not args.input_dir.exists() or not args.input_dir.is_dir():
        parser.error(f"input_dir is not a directory: {args.input_dir}")

    chains = load_chains(args.input_dir, args.glob)
    if len(chains) < 2:
        parser.error("need at least two parsed chains")

    pairs = pair_scores(chains)
    batch = greedy_batch(pairs, args.batch_size)

    total_rows = sum(len(chain.rows) for chain in chains)
    print(f"chains: {len(chains)}")
    print(f"rows: {total_rows}")
    print(f"pairs: {len(pairs)}")
    print(f"score_mean: {sum(float(pair['score']) for pair in pairs) / len(pairs):.3f}")
    print(f"score_max: {float(pairs[0]['score']):.3f}")
    print()
    print(f"Top {min(args.top, len(pairs))} pair scores")
    print_pairs(pairs, args.top)
    print()
    print(f"Greedy disjoint batch, m={args.batch_size}")
    print_pairs(batch, len(batch))

    if args.pairs_csv:
        write_pairs_csv(args.pairs_csv, pairs)
        print(f"\nwrote pairs CSV: {args.pairs_csv}")
    if args.matrix_csv:
        write_matrix_csv(args.matrix_csv, chains, pairs)
        print(f"wrote matrix CSV: {args.matrix_csv}")

    if any(math.isnan(float(pair["score"])) for pair in pairs):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
