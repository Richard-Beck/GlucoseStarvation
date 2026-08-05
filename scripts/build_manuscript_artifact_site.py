#!/usr/bin/env python3
"""Build the public manuscript artifact site from an assembly and report bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from build_manuscript_text_page import build_text_page


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ASSEMBLY = (
    PROJECT_ROOT
    / "manuscripts/20260731T101031_v03/revision_workspace/candidates"
    / "20260803_replace_s7_predictive_transfer_c02/assembly"
)
LANDING_MARKER = '<meta name="artifact-site" content="manuscript-landing">'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def staged_directory(site_dir: Path, name: str) -> Path:
    return Path(tempfile.mkdtemp(prefix=f".{name}.", dir=site_dir))


def install_directory(staged: Path, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    os.replace(staged, destination)


def resolve_source(assembly: Path, source_text: str) -> Path:
    source = Path(source_text)
    candidates = [source] if source.is_absolute() else [PROJECT_ROOT / source, assembly / source]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not resolve assembly source: {source_text}")


def page_css(*, wide: bool = False) -> str:
    width = "1440px" if wide else "940px"
    return f"""
    :root {{
      --ink: #20231f;
      --muted: #62685f;
      --line: #d7dbd2;
      --paper: #fbfaf6;
      --card: #ffffff;
      --accent: #176b73;
      --accent-wash: #e7f0ee;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: var(--paper);
      font-family: Georgia, "Times New Roman", serif;
      font-size: 17px;
      line-height: 1.62;
    }}
    .page {{ width: min({width}, calc(100% - 36px)); margin: 0 auto; padding: 58px 0 90px; }}
    header {{ padding-bottom: 30px; border-bottom: 1px solid var(--line); }}
    .eyebrow {{
      margin: 0 0 12px;
      color: var(--accent);
      font: 700 0.76rem/1.2 ui-monospace, "SFMono-Regular", Consolas, monospace;
      letter-spacing: 0.1em;
      text-transform: uppercase;
    }}
    h1 {{ margin: 0; max-width: 1100px; font-size: clamp(2.3rem, 7vw, 5rem); line-height: 1.03; }}
    h2 {{ margin: 42px 0 16px; font-size: 1.8rem; line-height: 1.2; }}
    h3 {{ margin: 26px 0 12px; font-size: 1.2rem; }}
    p {{ max-width: 880px; }}
    a {{ color: var(--accent); text-decoration-thickness: 0.08em; }}
    .lede {{ color: var(--muted); font-size: 1.12rem; max-width: 780px; }}
    .toolbar {{ display: flex; flex-wrap: wrap; gap: 10px 18px; margin: 22px 0; }}
    .button {{
      display: inline-block;
      padding: 8px 12px;
      border: 1px solid var(--line);
      background: var(--card);
      text-decoration: none;
      font: 650 0.85rem/1.3 ui-monospace, "SFMono-Regular", Consolas, monospace;
    }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(270px, 1fr)); gap: 16px; }}
    .card {{ padding: 20px; border: 1px solid var(--line); background: var(--card); }}
    .card h2, .card h3 {{ margin-top: 0; }}
    .meta {{ color: var(--muted); font: 0.82rem/1.5 ui-monospace, "SFMono-Regular", Consolas, monospace; }}
    details {{ margin: 12px 0; border: 1px solid var(--line); background: var(--card); }}
    summary {{ cursor: pointer; padding: 12px 14px; font-weight: 700; }}
    details > .content {{ padding: 0 14px 16px; }}
    pre {{ margin: 0; white-space: pre-wrap; overflow-wrap: anywhere; font-size: 0.84rem; line-height: 1.52; }}
    input[type="search"] {{
      width: min(560px, 100%);
      padding: 10px 12px;
      border: 1px solid var(--line);
      background: white;
      color: var(--ink);
      font: inherit;
    }}
    .table-wrap {{ width: 100%; overflow: auto; border: 1px solid var(--line); background: white; }}
    table {{ border-collapse: collapse; min-width: 1300px; font-size: 0.79rem; line-height: 1.4; }}
    th, td {{ padding: 8px 9px; border: 1px solid var(--line); text-align: left; vertical-align: top; }}
    th {{ position: sticky; top: 0; background: var(--accent-wash); }}
    ul.clean {{ padding-left: 20px; }}
    footer {{ margin-top: 48px; padding-top: 20px; border-top: 1px solid var(--line); color: var(--muted); }}
    @media (max-width: 640px) {{
      .page {{ padding-top: 32px; }}
      h1 {{ font-size: 2.4rem; }}
    }}
    """


def page(title: str, eyebrow: str, body: str, *, wide: bool = False) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>{page_css(wide=wide)}</style>
</head>
<body>
  <div class="page">
    <header>
      <p class="eyebrow">{html.escape(eyebrow)}</p>
      <h1>{html.escape(title)}</h1>
    </header>
    {body}
  </div>
</body>
</html>
"""


def prepare_report(site_dir: Path) -> tuple[Path, str]:
    root_index = site_dir / "index.html"
    existing_report = site_dir / "report" / "index.html"
    if root_index.is_file():
        root_text = root_index.read_text(encoding="utf-8")
    else:
        root_text = ""

    if root_text and LANDING_MARKER not in root_text:
        source_report = root_index
        report_text = root_text
        source_images = site_dir / "manuscript-images"
    elif existing_report.is_file():
        source_report = existing_report
        report_text = source_report.read_text(encoding="utf-8")
        source_images = existing_report.parent / "manuscript-images"
    else:
        raise FileNotFoundError("No illustrated manuscript report was found")

    if not source_images.is_dir():
        raise FileNotFoundError(f"Missing external manuscript images: {source_images}")

    staged = staged_directory(site_dir, "report")
    shutil.copy2(source_report, staged / "index.html")
    shutil.copytree(source_images, staged / "manuscript-images")
    install_directory(staged, site_dir / "report")

    if source_images.parent == site_dir and source_images.exists():
        shutil.rmtree(source_images)

    title_match = re.search(r"<h1>(.*?)</h1>", report_text, re.DOTALL)
    title = html.unescape(re.sub(r"<[^>]+>", "", title_match.group(1))).strip() if title_match else "Manuscript"
    return site_dir / "report" / "index.html", title


def prepare_text_page(site_dir: Path, report_path: Path) -> Path:
    staged = staged_directory(site_dir, "manuscript-text")
    output = build_text_page(report_path.read_text(encoding="utf-8"))
    write_text(staged / "index.html", output)
    install_directory(staged, site_dir / "manuscript-text")
    return site_dir / "manuscript-text" / "index.html"


def prepare_semantics(site_dir: Path, assembly: Path) -> tuple[Path, list[dict[str, str]]]:
    manifest_path = assembly / "source" / "figure_set_manifest.csv"
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        semantic_columns = [name for name in (reader.fieldnames or []) if "semantic" in name.lower()]
        if len(semantic_columns) != 1:
            raise ValueError(f"Expected one semantic column in {manifest_path}: {semantic_columns}")
        semantic_column = semantic_columns[0]
        rows = list(reader)

    staged = staged_directory(site_dir, "figure-semantics")
    panels_dir = staged / "panels"
    panels_dir.mkdir()
    entries: list[dict[str, str]] = []
    seen: set[str] = set()

    for row in rows:
        source_text = (row.get(semantic_column) or "").strip()
        if not source_text or source_text == "NA":
            continue
        figure_id = (row.get("current_figure_name") or "").strip()
        panel_id = (row.get("panel_id") or "").strip()
        if not figure_id or not panel_id:
            raise ValueError("Semantic manifest row lacks current_figure_name or panel_id")
        artifact_id = f"{figure_id}_{panel_id}"
        if artifact_id in seen:
            raise ValueError(f"Duplicate semantic interpretation: {artifact_id}")
        seen.add(artifact_id)
        source = resolve_source(assembly, source_text)
        destination = panels_dir / f"{artifact_id}.md"
        shutil.copy2(source, destination)
        entries.append(
            {
                "artifact_id": artifact_id,
                "figure_id": figure_id,
                "panel_id": panel_id,
                "source": source_text,
                "path": f"figure-semantics/panels/{artifact_id}.md",
            }
        )

    if not entries:
        raise ValueError("No semantic interpretations were selected by the manifest")

    combined = ["# Figure Semantic Interpretations", ""]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    detail_blocks: list[str] = []
    for entry in entries:
        grouped[entry["figure_id"]].append(entry)
        text = (panels_dir / f'{entry["artifact_id"]}.md').read_text(encoding="utf-8").strip()
        combined.extend([f'## {entry["artifact_id"]}', "", text, ""])
        detail_blocks.append(
            f'<details data-search="{html.escape((entry["artifact_id"] + " " + text).lower(), quote=True)}">'
            f'<summary>{html.escape(entry["artifact_id"])}</summary>'
            f'<div class="content"><p><a href="panels/{html.escape(entry["artifact_id"])}.md">Raw Markdown</a></p>'
            f'<pre>{html.escape(text)}</pre></div></details>'
        )
    write_text(staged / "all.md", "\n".join(combined))

    group_links = " ".join(
        f'<a class="button" href="#{html.escape(figure_id)}">{html.escape(figure_id)} ({len(items)})</a>'
        for figure_id, items in grouped.items()
    )
    grouped_details: list[str] = []
    offset = 0
    for figure_id, items in grouped.items():
        count = len(items)
        grouped_details.append(
            f'<section id="{html.escape(figure_id)}"><h2>{html.escape(figure_id)}</h2>'
            + "".join(detail_blocks[offset : offset + count])
            + "</section>"
        )
        offset += count
    body = f"""
    <p class="lede">Panel-level descriptions selected by the current figure-set manifest. Raw Markdown is preserved for direct retrieval.</p>
    <div class="toolbar"><a class="button" href="all.md">Combined Markdown</a>{group_links}</div>
    <p><input id="semantic-search" type="search" placeholder="Filter panels and interpretation text"></p>
    {''.join(grouped_details)}
    <script>
      const search = document.getElementById('semantic-search');
      const items = [...document.querySelectorAll('details[data-search]')];
      search.addEventListener('input', () => {{
        const query = search.value.trim().toLowerCase();
        items.forEach(item => item.hidden = query && !item.dataset.search.includes(query));
      }});
    </script>
    """
    write_text(staged / "index.html", page("Figure semantic interpretations", "Figure evidence", body))
    install_directory(staged, site_dir / "figure-semantics")
    return site_dir / "figure-semantics" / "index.html", entries


def parse_markdown_table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    start = next((index for index, line in enumerate(lines) if line.startswith("| id |")), None)
    if start is None:
        raise ValueError(f"Could not find provenance table in {path}")

    def cells(line: str) -> list[str]:
        return [cell.strip() for cell in line.strip().strip("|").split("|")]

    headers = cells(lines[start])
    rows = [cells(line) for line in lines[start + 2 :] if line.startswith("|")]
    if any(len(row) != len(headers) for row in rows):
        raise ValueError(f"Irregular provenance table rows in {path}")
    return headers, rows


def prepare_provenance(site_dir: Path, assembly: Path) -> tuple[Path, int]:
    source_dir = assembly / "evidence" / "methods"
    table_source = source_dir / "locked_provenance_table.md"
    verification_source = source_dir / "provenance_lock_verification.md"
    headers, rows = parse_markdown_table(table_source)
    staged = staged_directory(site_dir, "provenance")
    shutil.copy2(table_source, staged / table_source.name)
    shutil.copy2(verification_source, staged / verification_source.name)

    head_html = "".join(f"<th>{html.escape(value)}</th>" for value in headers)
    row_html = "\n".join(
        '<tr data-search="{}">{}</tr>'.format(
            html.escape(" ".join(row).lower(), quote=True),
            "".join(f"<td>{html.escape(value)}</td>" for value in row),
        )
        for row in rows
    )
    body = f"""
    <p class="lede">The locked recursive Methods provenance graph for manuscript-visible panel endpoints and their analytical dependencies.</p>
    <div class="toolbar">
      <a class="button" href="locked_provenance_table.md">Raw provenance Markdown</a>
      <a class="button" href="provenance_lock_verification.md">Lock verification</a>
    </div>
    <p class="meta">{len(rows)} provenance records</p>
    <p><input id="provenance-search" type="search" placeholder="Filter IDs, objects, roles, paths, or hashes"></p>
    <div class="table-wrap"><table><thead><tr>{head_html}</tr></thead><tbody>{row_html}</tbody></table></div>
    <script>
      const search = document.getElementById('provenance-search');
      const rows = [...document.querySelectorAll('tbody tr')];
      search.addEventListener('input', () => {{
        const query = search.value.trim().toLowerCase();
        rows.forEach(row => row.hidden = query && !row.dataset.search.includes(query));
      }});
    </script>
    """
    write_text(staged / "index.html", page("Methods provenance", "Traceability", body, wide=True))
    install_directory(staged, site_dir / "provenance")
    return site_dir / "provenance" / "index.html", len(rows)


def first_named_value(value: Any, names: tuple[str, ...]) -> str | None:
    if isinstance(value, dict):
        for name in names:
            candidate = value.get(name)
            if isinstance(candidate, (str, int)) and str(candidate).strip():
                return str(candidate).strip()
        for child in value.values():
            found = first_named_value(child, names)
            if found:
                return found
    return None


def render_value(value: Any) -> str:
    if isinstance(value, dict):
        return "<dl>" + "".join(
            f"<dt><strong>{html.escape(str(key))}</strong></dt><dd>{render_value(child)}</dd>"
            for key, child in value.items()
        ) + "</dl>"
    if isinstance(value, list):
        return "<ul class=" + '"clean">' + "".join(f"<li>{render_value(child)}</li>" for child in value) + "</ul>"
    if value is None:
        return '<span class="meta">null</span>'
    return html.escape(str(value))


def prepare_literature(site_dir: Path, assembly: Path) -> tuple[Path, dict[str, int]]:
    source = assembly / "evidence" / "abstract_introduction_discussion" / "old_bundle_1" / "literature_map.json"
    literature = json.loads(source.read_text(encoding="utf-8"))
    articles = literature.get("articles", [])
    if not isinstance(articles, list) or not articles:
        raise ValueError(f"No articles found in {source}")

    staged = staged_directory(site_dir, "literature-map")
    shutil.copy2(source, staged / "literature_map.json")
    article_blocks: list[str] = []
    for index, article in enumerate(articles, start=1):
        article_id = str(article.get("article_id") or f"article-{index}")
        title = first_named_value(article.get("article_info", {}), ("title", "article_title", "citation")) or article_id
        search_text = json.dumps(article, ensure_ascii=True).lower()
        article_blocks.append(
            f'<details data-search="{html.escape(search_text, quote=True)}">'
            f'<summary>{html.escape(title)}</summary>'
            f'<div class="content"><p class="meta">{html.escape(article_id)}</p>{render_value(article)}</div></details>'
        )

    counts = {
        "articles": len(articles),
        "anchors": len(literature.get("anchor_catalog", [])),
        "groups": len(literature.get("group_catalog", [])),
    }
    generated_on = str(literature.get("generated_on", "unknown"))
    body = f"""
    <p class="lede">A structured reference-discovery and claim-anchoring map inherited from the accepted prior A/I/D bundle.</p>
    <p class="meta">Generated {html.escape(generated_on)} | {counts['articles']} articles | {counts['anchors']} claim anchors | {counts['groups']} literature groups</p>
    <p><strong>Status:</strong> accepted inherited artifact; it was not regenerated for the current A/I/D context bundle.</p>
    <div class="toolbar"><a class="button" href="literature_map.json">Raw literature map JSON</a></div>
    <p><input id="literature-search" type="search" placeholder="Filter titles, claims, groups, or manuscript relevance"></p>
    <section>{''.join(article_blocks)}</section>
    <script>
      const search = document.getElementById('literature-search');
      const items = [...document.querySelectorAll('details[data-search]')];
      search.addEventListener('input', () => {{
        const query = search.value.trim().toLowerCase();
        items.forEach(item => item.hidden = query && !item.dataset.search.includes(query));
      }});
    </script>
    """
    write_text(staged / "index.html", page("Structured literature map", "A/I/D evidence", body))
    install_directory(staged, site_dir / "literature-map")
    return site_dir / "literature-map" / "index.html", counts


def prepare_claim_graph(site_dir: Path, assembly: Path) -> tuple[Path, dict[str, int]]:
    source = assembly / "evidence" / "claim_graph" / "claim_graph_integrated.json"
    graph = json.loads(source.read_text(encoding="utf-8"))
    collections = {
        "claims": graph.get("claims", []),
        "evidence": graph.get("evidence", []),
        "relationships": graph.get("relationships", []),
    }
    for name, records in collections.items():
        if not isinstance(records, list):
            raise ValueError(f"Claim graph collection is not a list: {name}")

    shutil.copy2(source, site_dir / "claim-graph.json")
    staged = staged_directory(site_dir, "claim-graph")
    collection_html: list[str] = []
    for name, records in collections.items():
        blocks: list[str] = []
        for index, record in enumerate(records, start=1):
            if not isinstance(record, dict):
                raise ValueError(f"Claim graph {name} record {index} is not an object")
            record_id = str(record.get("id") or f"{name}-{index}")
            label = str(record.get("text") or record.get("description") or record_id)
            search_text = json.dumps(record, ensure_ascii=True).lower()
            blocks.append(
                f'<details data-search="{html.escape(search_text, quote=True)}">'
                f'<summary>{html.escape(record_id)}: {html.escape(label)}</summary>'
                f'<div class="content">{render_value(record)}</div></details>'
            )
        collection_html.append(
            f'<section id="{html.escape(name)}"><h2>{html.escape(name.title())} ({len(records)})</h2>'
            + "".join(blocks)
            + "</section>"
        )

    metadata = graph.get("metadata", {})
    body = f"""
    <p class="lede">The existing integrated claim records, evidence records, and relationships, published without schema translation or scientific revalidation.</p>
    <div class="toolbar">
      <a class="button" href="../claim-graph.json">Complete claim graph JSON</a>
      <a class="button" href="#claims">Claims</a>
      <a class="button" href="#evidence">Evidence</a>
      <a class="button" href="#relationships">Relationships</a>
    </div>
    <details><summary>Graph metadata</summary><div class="content">{render_value(metadata)}</div></details>
    <p><input id="claim-search" type="search" placeholder="Filter claims, evidence, or relationships"></p>
    {''.join(collection_html)}
    <script>
      const search = document.getElementById('claim-search');
      const items = [...document.querySelectorAll('details[data-search]')];
      search.addEventListener('input', () => {{
        const query = search.value.trim().toLowerCase();
        items.forEach(item => item.hidden = query && !item.dataset.search.includes(query));
      }});
    </script>
    """
    write_text(staged / "index.html", page("Claim graph", "Claim and evidence structure", body))
    install_directory(staged, site_dir / "claim-graph")
    counts = {name: len(records) for name, records in collections.items()}
    return site_dir / "claim-graph" / "index.html", counts


def current_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=PROJECT_ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout.strip()


def prepare_release_identity(
    site_dir: Path,
    *,
    manuscript_version: str,
    repository: str,
    source_commit: str,
    release_tag: str,
) -> dict[str, str]:
    identity = {
        "site_schema_version": "1.1",
        "manuscript_version": manuscript_version,
        "repository": repository,
        "source_commit": source_commit,
        "release_tag": release_tag,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "unsealed_revision_preview",
    }
    write_text(site_dir / "release.json", json.dumps(identity, indent=2) + "\n")
    return identity


def prepare_landing(
    site_dir: Path,
    title: str,
    semantic_count: int,
    provenance_count: int,
    literature_counts: dict[str, int],
    claim_counts: dict[str, int],
    release_identity: dict[str, str],
) -> Path:
    repository_url = f'https://github.com/{release_identity["repository"]}'
    commit_url = f'{repository_url}/commit/{release_identity["source_commit"]}'
    release_url = f'{repository_url}/releases/tag/{release_identity["release_tag"]}'
    short_commit = release_identity["source_commit"][:12]
    body = f"""
    <meta name="artifact-site" content="manuscript-landing">
    <p class="lede">A public review and evidence package for the current manuscript, organized for both human inspection and machine retrieval.</p>
    <section class="grid" aria-label="Primary manuscript artifacts">
      <article class="card"><p class="eyebrow">Primary report</p><h2><a href="report/">Illustrated manuscript</a></h2><p>Complete review copy with all main and supplemental figures.</p></article>
      <article class="card"><p class="eyebrow">Lightweight text</p><h2><a href="manuscript-text/">Text-only manuscript</a></h2><p>Abstract, Introduction, Results, Discussion, Methods, and References without figures.</p></article>
      <article class="card"><p class="eyebrow">Figure evidence</p><h2><a href="figure-semantics/">Semantic interpretations</a></h2><p>{semantic_count} panel-level interpretations with individual and combined Markdown endpoints.</p></article>
      <article class="card"><p class="eyebrow">Claim structure</p><h2><a href="claim-graph/">Claim graph</a></h2><p>{claim_counts['claims']} claims, {claim_counts['evidence']} evidence records, and {claim_counts['relationships']} existing relationships.</p></article>
      <article class="card"><p class="eyebrow">Traceability</p><h2><a href="provenance/">Methods provenance</a></h2><p>{provenance_count} locked provenance records with raw Markdown and verification.</p></article>
      <article class="card"><p class="eyebrow">Literature</p><h2><a href="literature-map/">Structured literature map</a></h2><p>{literature_counts['articles']} articles, {literature_counts['anchors']} claim anchors, and {literature_counts['groups']} literature groups.</p></article>
    </section>
    <section>
      <h2>Machine-readable entry points</h2>
      <div class="toolbar">
        <a class="button" href="llms.txt">llms.txt</a>
        <a class="button" href="site-manifest.json">Site manifest JSON</a>
        <a class="button" href="release.json">Release identity JSON</a>
        <a class="button" href="claim-graph.json">Claim graph JSON</a>
        <a class="button" href="figure-semantics/all.md">All figure semantics</a>
        <a class="button" href="provenance/locked_provenance_table.md">Raw provenance</a>
        <a class="button" href="literature-map/literature_map.json">Raw literature map</a>
      </div>
    </section>
    <section>
      <h2>Package identity</h2>
      <p>This unsealed revision-preview package was assembled in the context of repository commit <a href="{html.escape(commit_url)}"><code>{html.escape(short_commit)}</code></a>. It is not represented as having been reproduced or scientifically validated from that commit.</p>
      <div class="toolbar">
        <a class="button" href="{html.escape(repository_url)}">Repository</a>
        <a class="button" href="{html.escape(commit_url)}">Associated commit</a>
        <a class="button" href="{html.escape(release_url)}">Release asset</a>
      </div>
      <p class="meta">Manuscript {html.escape(release_identity['manuscript_version'])} | generated {html.escape(release_identity['generated_at'])}</p>
    </section>
    <footer><p>This site publishes review artifacts from a versioned manuscript assembly. The structured literature map is explicitly retained from the accepted prior A/I/D bundle.</p></footer>
    """
    landing = page(title, "Manuscript artifact portal", body)
    landing = landing.replace("</head>", f"  {LANDING_MARKER}\n</head>")
    write_text(site_dir / "index.html", landing)
    return site_dir / "index.html"


def prepare_machine_indexes(
    site_dir: Path,
    assembly: Path,
    base_url: str,
    semantic_entries: list[dict[str, str]],
) -> tuple[Path, Path]:
    base_url = base_url.rstrip("/")
    llms = f"""# Glucose Starvation manuscript artifacts

> Public manuscript review, figure interpretation, provenance, and literature artifacts.

- Full illustrated manuscript: {base_url}/report/
- Text-only manuscript: {base_url}/manuscript-text/
- Figure semantic interpretation index: {base_url}/figure-semantics/
- Combined figure semantic interpretations: {base_url}/figure-semantics/all.md
- Claim graph browser: {base_url}/claim-graph/
- Complete claim graph JSON: {base_url}/claim-graph.json
- Methods provenance index: {base_url}/provenance/
- Raw locked provenance table: {base_url}/provenance/locked_provenance_table.md
- Structured literature map index: {base_url}/literature-map/
- Raw structured literature map: {base_url}/literature-map/literature_map.json
- Release identity: {base_url}/release.json
- Machine-readable site manifest: {base_url}/site-manifest.json
"""
    llms_path = site_dir / "llms.txt"
    write_text(llms_path, llms)

    artifact_specs = [
        ("index.html", "landing_page", None),
        ("report/index.html", "illustrated_manuscript", None),
        ("manuscript-text/index.html", "text_only_manuscript", None),
        ("figure-semantics/index.html", "semantic_interpretation_index", None),
        ("figure-semantics/all.md", "combined_semantic_interpretations", None),
        ("claim-graph/index.html", "claim_graph_index", None),
        ("claim-graph.json", "integrated_claim_graph", "evidence/claim_graph/claim_graph_integrated.json"),
        ("provenance/index.html", "provenance_index", None),
        ("provenance/locked_provenance_table.md", "locked_provenance_table", "evidence/methods/locked_provenance_table.md"),
        ("provenance/provenance_lock_verification.md", "provenance_verification", "evidence/methods/provenance_lock_verification.md"),
        ("literature-map/index.html", "literature_map_index", None),
        ("literature-map/literature_map.json", "structured_literature_map", "evidence/abstract_introduction_discussion/old_bundle_1/literature_map.json"),
        ("release.json", "release_identity", None),
        ("llms.txt", "llm_directory", None),
    ]
    artifacts: list[dict[str, Any]] = []
    for relative, role, source in artifact_specs:
        path = site_dir / relative
        entry: dict[str, Any] = {
            "path": relative,
            "url": f"{base_url}/{relative}",
            "role": role,
            "sha256": sha256(path),
            "byte_size": path.stat().st_size,
        }
        if source:
            entry["assembly_source"] = source
        artifacts.append(entry)
    for semantic in semantic_entries:
        path = site_dir / semantic["path"]
        artifacts.append(
            {
                "path": semantic["path"],
                "url": f'{base_url}/{semantic["path"]}',
                "role": "panel_semantic_interpretation",
                "artifact_id": semantic["artifact_id"],
                "assembly_source": semantic["source"],
                "sha256": sha256(path),
                "byte_size": path.stat().st_size,
            }
        )

    try:
        assembly_reference = str(assembly.relative_to(PROJECT_ROOT))
    except ValueError:
        assembly_reference = str(assembly)
    manifest = {
        "schema_version": "1.1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "base_url": base_url,
        "assembly_root": assembly_reference,
        "literature_map_status": "accepted_inherited_artifact_generated_2026-07-13",
        "artifacts": artifacts,
    }
    manifest_path = site_dir / "site-manifest.json"
    write_text(manifest_path, json.dumps(manifest, indent=2) + "\n")
    return llms_path, manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the manuscript artifact landing site.")
    parser.add_argument(
        "--site-dir",
        type=Path,
        default=PROJECT_ROOT / "latest_manuscript_bundle",
    )
    parser.add_argument("--assembly-root", type=Path, default=DEFAULT_ASSEMBLY)
    parser.add_argument(
        "--base-url",
        default="https://richard-beck.github.io/GlucoseStarvation",
    )
    parser.add_argument("--manuscript-version", default="20260731T101031_v03")
    parser.add_argument("--repository", default="Richard-Beck/GlucoseStarvation")
    parser.add_argument("--release-tag", default="pages-site")
    parser.add_argument("--source-commit")
    args = parser.parse_args()
    site_dir = args.site_dir.resolve()
    assembly = args.assembly_root.resolve()
    if not site_dir.is_dir():
        raise FileNotFoundError(f"Missing site directory: {site_dir}")
    if not assembly.is_dir():
        raise FileNotFoundError(f"Missing assembly root: {assembly}")

    report_path, manuscript_title = prepare_report(site_dir)
    prepare_text_page(site_dir, report_path)
    _, semantic_entries = prepare_semantics(site_dir, assembly)
    _, claim_counts = prepare_claim_graph(site_dir, assembly)
    _, provenance_count = prepare_provenance(site_dir, assembly)
    _, literature_counts = prepare_literature(site_dir, assembly)
    release_identity = prepare_release_identity(
        site_dir,
        manuscript_version=args.manuscript_version,
        repository=args.repository,
        source_commit=args.source_commit or current_commit(),
        release_tag=args.release_tag,
    )
    prepare_landing(
        site_dir,
        manuscript_title,
        len(semantic_entries),
        provenance_count,
        literature_counts,
        claim_counts,
        release_identity,
    )
    prepare_machine_indexes(site_dir, assembly, args.base_url, semantic_entries)
    print(
        f"Built artifact site with {len(semantic_entries)} semantic interpretations, "
        f"{provenance_count} provenance records, and "
        f"{literature_counts['articles']} literature records in {site_dir}"
    )


if __name__ == "__main__":
    main()
