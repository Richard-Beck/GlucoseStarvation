#!/usr/bin/env python3
"""Render the self-contained HTML draft from approved v03 assembly inputs."""

from __future__ import annotations

import base64
import csv
import html
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = ROOT.parents[1]
RUN_REL = ROOT.relative_to(PROJECT_ROOT)
MANIFEST = ROOT / "source" / "figure_set_manifest.csv"
LEGENDS = ROOT / "source" / "integrated_figure_legends.md"
LEGEND_VALIDATION = ROOT / "validation" / "draft_render_legend_validation_report.md"
OUT = ROOT / "draft" / "manuscript_draft.html"
RESULTS_DIR = ROOT / "source" / "manuscript_sections" / "results"
SECTIONS_DIR = ROOT / "source" / "manuscript_sections"
METHODS_SOURCE = ROOT / "source" / "manuscript_sections" / "methods.md"
TITLE_SOURCE = ROOT / "source" / "title.txt"

TITLE = TITLE_SOURCE.read_text(encoding="utf-8").strip()
if not TITLE:
    raise ValueError(f"Empty manuscript title: {TITLE_SOURCE}")

CONFIG = json.loads(Path(__file__).with_name("renderer_config.json").read_text())
EXPECTED_CONFIG = {
    "assembly_root": "manuscripts/20260731T101031_v03",
    "manifest": "source/figure_set_manifest.csv",
    "legends": "source/integrated_figure_legends.md",
    "results_dir": "source/manuscript_sections/results",
    "methods": "source/manuscript_sections/methods.md",
    "output": "draft/manuscript_draft.html",
    "expected_panel_rows": 79,
    "expected_figures": 23,
    "expected_results_sections": 5,
    "expected_user_fixed_claims": 8,
}
if CONFIG != EXPECTED_CONFIG:
    raise ValueError(f"renderer_config.json differs from the package contract: {CONFIG}")


def e(text: object) -> str:
    return html.escape(str(text), quote=True)


def display_label(figure_id: str) -> str:
    if re.match(r"^S\d+", figure_id):
        suffix = figure_id.removeprefix("S")
        if suffix.endswith("_continued"):
            return f"Figure S{suffix.removesuffix('_continued')} continued"
        return f"Figure S{suffix}"
    if re.match(r"^F\d+", figure_id):
        return f"Figure {figure_id.removeprefix('F')}"
    if figure_id.startswith("Figure_S"):
        suffix = figure_id.removeprefix("Figure_S")
        if suffix.endswith("_continued"):
            return f"Figure S{suffix.removesuffix('_continued')} continued"
        return f"Figure S{suffix}"
    if figure_id.startswith("Figure_"):
        return f"Figure {figure_id.removeprefix('Figure_')}"
    return figure_id.replace("_", " ")


def figure_anchor(figure_id: str) -> str:
    return f"fig-{figure_id.lower().replace('_', '-')}"


def label_to_figure_id(label: str) -> str:
    label = " ".join(label.strip().split())
    match = re.match(r"^Figure\s+(S?)(\d+)([A-Z]?)(?:\s+(continued))?$", label)
    if not match:
        raise ValueError(f"Could not parse figure label: {label}")
    is_supp, number, letter, continued = match.groups()
    if is_supp:
        figure_id = f"S{number}{letter}"
    else:
        figure_id = f"F{number}"
    if continued:
        figure_id += "_continued"
    return figure_id


def callout_to_figure_id(label: str, known_figures: set[str]) -> str:
    """Resolve a text callout to the rendered whole-figure anchor.

    Figure S8A/S8B are separate rendered figures, while callouts such as
    Figure S5A usually mean panel A within the rendered Figure S5.
    """

    figure_id = label_to_figure_id(label)
    if figure_id in known_figures:
        return figure_id
    if figure_id.startswith("Figure_S") and not figure_id.endswith("_continued"):
        suffix = figure_id.removeprefix("Figure_S")
        match = re.match(r"^(\d+)[A-Z]$", suffix)
        if match:
            parent = f"Figure_S{match.group(1)}"
            if parent in known_figures:
                return parent
    panel_parent = re.match(r"^(F|S)(\d+)[A-Z]$", figure_id)
    if panel_parent:
        parent = f"{panel_parent.group(1)}{panel_parent.group(2)}"
        if parent in known_figures:
            return parent
    return figure_id


def figure_sort_key(figure_id: str) -> tuple[int, int, str, int]:
    if re.match(r"^S\d+", figure_id):
        suffix = figure_id.removeprefix("S")
        continued = 1 if suffix.endswith("_continued") else 0
        suffix = suffix.removesuffix("_continued")
        match = re.match(r"^(\d+)([A-Z]?)$", suffix)
        if not match:
            return (1, 10_000, suffix, continued)
        return (1, int(match.group(1)), match.group(2), continued)
    if re.match(r"^F\d+", figure_id):
        return (0, int(figure_id.removeprefix("F")), "", 0)
    if figure_id.startswith("Figure_S"):
        suffix = figure_id.removeprefix("Figure_S")
        continued = 1 if suffix.endswith("_continued") else 0
        suffix = suffix.removesuffix("_continued")
        match = re.match(r"^(\d+)([A-Z]?)$", suffix)
        if not match:
            return (1, 10_000, suffix, continued)
        return (1, int(match.group(1)), match.group(2), continued)
    if figure_id.startswith("Figure_"):
        return (0, int(figure_id.removeprefix("Figure_")), "", 0)
    return (2, 10_000, figure_id, 0)


def repo_path_to_data_uri(path_text: str) -> str:
    repo_rel = Path(path_text)
    abs_path = PROJECT_ROOT / repo_rel
    if not abs_path.exists():
        raise FileNotFoundError(path_text)
    encoded = base64.b64encode(abs_path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def parse_manifest() -> dict[str, dict[str, object]]:
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    figures: dict[str, dict[str, object]] = {}
    for row in rows:
        figure_id = row["current_figure_name"]
        if figure_id == "NA":
            continue
        fig = figures.setdefault(
            figure_id,
            {
                "representative": row,
                "panels": [],
                "all_rows": [],
                "main_or_supplement": "main" if figure_id.startswith("F") else "supplement",
            },
        )
        fig["all_rows"].append(row)
        if row["panel_id"] not in {"whole_figure", "no_output", "NA"}:
            fig["panels"].append(row)
        fig["representative"] = row

    for figure_id, fig in figures.items():
        image_path = ROOT / "assets" / "figures" / f"{figure_id}.png"
        if not image_path.exists():
            raise FileNotFoundError(image_path)
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        fig["html_image_path"] = f"data:image/png;base64,{encoded}"
    return figures


def parse_legends() -> dict[str, dict[str, str]]:
    legends: dict[str, dict[str, str]] = {}
    current: str | None = None
    body_lines: list[str] = []

    def finish_current() -> None:
        if current is None:
            return
        body = "\n".join(body_lines).strip()
        if not body:
            raise ValueError(f"Empty legend body for {current}")
        legends[current]["body"] = body

    for line in LEGENDS.read_text().splitlines():
        header = re.match(r"^##\s+(Figure\s+S?\d+[A-Z]?(?:\s+continued)?)\.\s+(.+)$", line)
        if header:
            finish_current()
            current = label_to_figure_id(header.group(1))
            if current in legends:
                raise ValueError(f"Duplicate legend block for {current}")
            legends[current] = {"title": header.group(2).strip(), "body": ""}
            body_lines = []
            continue
        if current:
            body_lines.append(line)
    finish_current()
    return legends


def write_legend_validation_report(
    figures: dict[str, dict[str, object]], legends: dict[str, dict[str, str]]
) -> None:
    expected = set(figures)
    observed = set(legends)
    missing = sorted(expected - observed, key=figure_sort_key)
    extra = sorted(observed - expected, key=figure_sort_key)
    status = "PASS" if not missing and not extra else "FAIL"
    lines = [
        "# Legend Validation Report",
        "",
        f"Status: {status}",
        "",
        f"Expected rendered figure legends: {len(expected)}",
        f"Observed legend blocks: {len(observed)}",
        "",
    ]
    if missing:
        lines.extend(["Missing legend blocks:", *[f"- {display_label(fid)}" for fid in missing], ""])
    if extra:
        lines.extend(["Extra legend blocks:", *[f"- {display_label(fid)}" for fid in extra], ""])
    if status == "PASS":
        lines.append("Every rendered figure has exactly one legend block and no extra figure legend blocks were found.")
    LEGEND_VALIDATION.write_text("\n".join(lines) + "\n")
    if status != "PASS":
        raise ValueError("Legend validation failed; see legend_validation_report.md")


def link_figure_refs(html_text: str, known_figures: set[str]) -> str:
    """Link escaped figure callouts such as Figure 4a or Figure S10 continued c."""

    def replace(match: re.Match[str]) -> str:
        figure_label = match.group(1)
        panel = match.group(2) or ""
        figure_id = callout_to_figure_id(figure_label, known_figures)
        label = f"{figure_label}{(' ' + panel) if panel and 'continued' in figure_label else panel}"
        return f'<a href="#{figure_anchor(figure_id)}">{label}</a>'

    pattern = r"\b(Figure\s+S?\d+[A-Z]?(?:\s+continued)?)(?:\s*([a-z]))?\b"
    return re.sub(pattern, replace, html_text)


def format_math_expr(text: str) -> str:
    escaped = e(text)

    def subscript_repl(match: re.Match[str]) -> str:
        base, sub = match.groups()
        if base == "sum":
            return f'<span class="math-op">&sum;</span><sub>{e(sub)}</sub>'
        if base == "mu":
            base_html = "&mu;"
        elif base == "Phi":
            base_html = "&Phi;"
        elif base == "ae":
            return f"<i>a</i><sub>e,{e(sub)}</sub>"
        elif base == "ah":
            return f"<i>a</i><sub>h,{e(sub)}</sub>"
        else:
            base_html = e(base)
        return f"<i>{base_html}</i><sub>{e(sub)}</sub>"

    escaped = re.sub(r"\b([A-Za-z]+)_([A-Za-z0-9]+)\b", subscript_repl, escaped)
    escaped = re.sub(r"(?<!&)\bPhi\b(?!;)", r"<i>&Phi;</i>", escaped)
    escaped = re.sub(r"\b(log|exp)\b", r'<span class="math-fn">\1</span>', escaped)
    escaped = escaped.replace("*", "&middot;")
    escaped = re.sub(r"(?<=\s)-(?=\s|\w)", "&minus;", escaped)
    return escaped


def markdown_inline(text: str, known_figures: set[str], *, enable_math: bool = False) -> str:
    code_spans: list[str] = []

    def stash_code(match: re.Match[str]) -> str:
        code_spans.append(f"<code>{e(match.group(1))}</code>")
        return f"@@CODE{len(code_spans) - 1}@@"

    text = re.sub(r"`([^`]+)`", stash_code, text)
    escaped = e(text)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = link_figure_refs(escaped, known_figures)
    if enable_math:
        math_patterns = [
            r"Lum = aG \+ b",
            r"log\[\(N_high \+ 1e-8\)/\(N_low \+ 1e-8\)\]",
            r"dN/dt = gN; dR_i/dt = -u_i N; and dW/dt = K u N",
        ]
        for pattern in math_patterns:
            escaped = re.sub(
                pattern,
                lambda match: f'<span class="math-inline">{format_math_expr(match.group(0))}</span>',
                escaped,
            )
    for idx, code_html in enumerate(code_spans):
        escaped = escaped.replace(f"@@CODE{idx}@@", code_html)
    return escaped


def is_math_block(block: str) -> bool:
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    if len(lines) != 1:
        return False
    line = lines[0].rstrip(".,")
    return bool(re.match(r"^(?:[A-Za-z]+(?:_[A-Za-z0-9]+)?|d[A-Za-z]/dt)\s*=", line))


def markdown_blocks(
    text: str,
    known_figures: set[str],
    *,
    enable_math: bool = False,
    extended_markdown: bool = False,
) -> str:
    chunks: list[str] = []
    for block in [p.strip() for p in re.split(r"\n\s*\n", text.strip()) if p.strip()]:
        lines = [line.strip() for line in block.splitlines() if line.strip()]
        is_table = (
            extended_markdown
            and len(lines) >= 2
            and "|" in lines[0]
            and bool(re.match(r"^\|?\s*:?-+", lines[1]))
        )
        if is_table:
            def cells(line: str) -> list[str]:
                return [cell.strip() for cell in line.strip().strip("|").split("|")]

            header = cells(lines[0])
            body_rows = [cells(line) for line in lines[2:]]
            head_html = "".join(
                f"<th>{markdown_inline(cell, known_figures, enable_math=enable_math)}</th>"
                for cell in header
            )
            body_html = "\n".join(
                "<tr>"
                + "".join(
                    f"<td>{markdown_inline(cell, known_figures, enable_math=enable_math)}</td>"
                    for cell in row
                )
                + "</tr>"
                for row in body_rows
            )
            chunks.append(
                '<div class="table-wrap"><table><thead><tr>'
                + head_html
                + "</tr></thead><tbody>"
                + body_html
                + "</tbody></table></div>"
            )
            continue
        if extended_markdown:
            heading = re.match(r"^(#{1,4})\s+(.+)$", block)
        else:
            heading = re.match(r"^(#{3,4})\s+(.+)$", block)
        if heading:
            level = len(heading.group(1))
            if extended_markdown:
                level = max(2, min(4, level))
            chunks.append(
                f"<h{level}>{markdown_inline(heading.group(2), known_figures, enable_math=enable_math)}</h{level}>"
            )
        elif extended_markdown and all(
            line.startswith("- ") for line in block.splitlines() if line.strip()
        ):
            items = "\n".join(
                f"<li>{markdown_inline(line.strip()[2:], known_figures, enable_math=enable_math)}</li>"
                for line in block.splitlines()
                if line.strip()
            )
            chunks.append(f'<ul class="manuscript-list">\n{items}\n</ul>')
        elif extended_markdown and all(
            re.match(r"^\d+\.\s+", line) for line in lines
        ):
            items = "\n".join(
                f"<li>{markdown_inline(line.split('. ', 1)[1], known_figures, enable_math=enable_math)}</li>"
                for line in lines
            )
            chunks.append(f'<ol class="manuscript-list">\n{items}\n</ol>')
        elif enable_math and is_math_block(block):
            formula = block.strip()
            trailing = ""
            if formula[-1:] in {".", ","}:
                trailing = e(formula[-1])
                formula = formula[:-1]
            chunks.append(f'<div class="math-block">{format_math_expr(formula)}{trailing}</div>')
        else:
            chunks.append(f"<p>{markdown_inline(block, known_figures, enable_math=enable_math)}</p>")
    return "\n".join(chunks)


def figure_block(
    figure_id: str,
    figures: dict[str, dict[str, object]],
    legends: dict[str, dict[str, str]],
    *,
    compact: bool = False,
) -> str:
    fig = figures[figure_id]
    representative = fig["representative"]
    assert isinstance(representative, dict)
    legend = legends[figure_id]
    title = legend["title"]
    body = legend["body"]
    src = fig["html_image_path"]
    panel_count = len(fig["panels"])
    panel_word = "panel" if panel_count == 1 else "panels"
    css_class = "figure-card compact" if compact else "figure-card"
    legend_html = markdown_blocks(body, set(figures))
    return f"""
<figure class="{css_class}" id="{figure_anchor(figure_id)}">
  <img src="{e(src)}" alt="{e(display_label(figure_id) + '. ' + title)}" loading="lazy">
  <figcaption>
    <strong>{e(display_label(figure_id))}. {e(title)}</strong>
    <div class="legend-body">{legend_html}</div>
  </figcaption>
</figure>
""".strip()


RESULT_SECTIONS = [
    {
        "id": "results_measurement_foundation",
        "filename": "results_measurement_foundation.md",
        "fallback_title": "Measurement foundations establish ploidy-linked starvation responses",
        "figure": "F1",
    },
    {
        "id": "results_direct_features",
        "filename": "results_direct_features.md",
        "fallback_title": "Direct trajectory features separate starvation robustness from cell-number yield costs",
        "figure": "F2",
    },
    {
        "id": "results_mechanistic_models",
        "filename": "results_mechanistic_models.md",
        "fallback_title": "Hidden-state model families explain trajectory structure but limit transferability",
        "figure": "F3",
    },
    {
        "id": "results_posterior_size_context",
        "filename": "results_posterior_size_context.md",
        "fallback_title": "High ploidy links higher apparent uptake, lower per-cell yield, and starvation robustness in a size-aware context",
        "figure": "F4",
    },
    {
        "id": "results_selection_simulations",
        "filename": "results_selection_simulations.md",
        "fallback_title": "Fitted simulations predict resource-regime-dependent selection for high or low ploidy",
        "figure": "F5",
    },
]


def parse_front_matter(text: str, path: Path) -> tuple[dict[str, str], str]:
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError(f"Missing front matter in {path}")
    end = None
    for idx, line in enumerate(lines[1:], start=1):
        if line.strip() == "---":
            end = idx
            break
    if end is None:
        raise ValueError(f"Unclosed front matter in {path}")

    metadata: dict[str, str] = {}
    for line in lines[1:end]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"Malformed front matter line in {path}: {line}")
        key, value = line.split(":", 1)
        metadata[key.strip()] = value.strip().strip('"')
    return metadata, "\n".join(lines[end + 1 :])


def extract_markdown_section(text: str, heading: str, path: Path) -> str:
    lines = text.splitlines()
    start = None
    pattern = re.compile(rf"^##\s+{re.escape(heading)}\s*$")
    for idx, line in enumerate(lines):
        if pattern.match(line):
            start = idx + 1
            break
    if start is None:
        raise ValueError(f"Missing '## {heading}' in {path}")

    end = len(lines)
    for idx in range(start, len(lines)):
        if lines[idx].startswith("## "):
            end = idx
            break
    body = "\n".join(lines[start:end]).strip()
    if not body:
        raise ValueError(f"Empty '## {heading}' in {path}")
    return body


def load_manuscript_section(section_id: str) -> dict[str, str]:
    path = SECTIONS_DIR / f"{section_id}.md"
    if not path.exists():
        raise FileNotFoundError(f"Missing manuscript section: {path}")
    raw_text = path.read_text().strip()
    if raw_text.startswith("---"):
        metadata, body = parse_front_matter(raw_text, path)
        if metadata.get("section_id") != section_id:
            raise ValueError(f"Section {path} has section_id={metadata.get('section_id')!r}")
        title = metadata.get("title", section_id.title())
        text = extract_markdown_section(body, "Manuscript Text", path)
    else:
        expected_heading = "References cited" if section_id == "references" else section_id.title()
        text = extract_markdown_section(raw_text, expected_heading, path)
        title = expected_heading
    return {
        "id": section_id,
        "title": title,
        "text": text,
        "path": str(path.relative_to(PROJECT_ROOT)),
    }


def load_methods_section() -> dict[str, str]:
    if not METHODS_SOURCE.exists():
        raise FileNotFoundError(f"Missing Methods source: {METHODS_SOURCE}")
    raw_text = METHODS_SOURCE.read_text().strip()
    if not raw_text:
        raise ValueError(f"Empty Methods source: {METHODS_SOURCE}")
    if raw_text.startswith("---"):
        metadata, body = parse_front_matter(raw_text, METHODS_SOURCE)
        text = extract_markdown_section(body, "Manuscript Text", METHODS_SOURCE)
        title = metadata.get("title", "Methods")
    else:
        text = raw_text
        title = "Methods"
    return {
        "id": "methods",
        "title": title,
        "text": text,
        "path": str(METHODS_SOURCE.relative_to(ROOT)),
    }


def load_result_sections() -> list[dict[str, object]]:
    loaded: list[dict[str, object]] = []
    for section in RESULT_SECTIONS:
        section_id = str(section["id"])
        path = RESULTS_DIR / str(section["filename"])
        if not path.exists():
            raise FileNotFoundError(f"Missing Results sidecar for {section_id}: {path}")
        metadata, body = parse_front_matter(path.read_text(), path)
        sidecar_section_id = metadata.get("section_id")
        if sidecar_section_id != section_id:
            raise ValueError(
                f"Results sidecar {path} has section_id={sidecar_section_id!r}, expected {section_id!r}"
            )
        if "## Revision Notes" not in body:
            raise ValueError(f"Missing Revision Notes review notes in {path}")
        loaded.append(
            {
                **section,
                "title": metadata.get("title") or section["fallback_title"],
                "results_text": extract_markdown_section(body, "Results Text", path),
                "sidecar_path": str(path.relative_to(PROJECT_ROOT)),
            }
        )
    return loaded


def render_toc(result_sections: list[dict[str, object]]) -> str:
    result_links = "\n".join(
        f'  <a href="#{e(section["id"])}">{e(section["title"])}</a>'
        for section in result_sections
    )
    return f"""
<nav class="toc" aria-label="Table of contents">
  <strong>Contents</strong>
  <a href="#abstract">Abstract</a>
  <a href="#introduction">Introduction</a>
{result_links}
  <a href="#discussion">Discussion</a>
  <a href="#methods">Methods</a>
  <a href="#references">References</a>
  <a href="#figures">Figure Index</a>
  <a href="#audit-trail">Audit Trail</a>
</nav>
""".strip()


def render_results(
    figures: dict[str, dict[str, object]],
    legends: dict[str, dict[str, str]],
    result_sections: list[dict[str, object]],
) -> str:
    chunks: list[str] = []
    for section in result_sections:
        results_text = markdown_blocks(str(section["results_text"]), set(figures))
        chunks.append(
            f"""
<section id="{e(section["id"])}">
  <h3>{e(section["title"])}</h3>
  {results_text}
  {figure_block(str(section["figure"]), figures, legends)}
</section>
""".strip()
        )
    return "\n\n".join(chunks)


def render_figure_index(
    figures: dict[str, dict[str, object]],
    legends: dict[str, dict[str, str]],
) -> str:
    main = []
    supplements = []
    for fid, fig in figures.items():
        if fig["main_or_supplement"] == "main":
            main.append(fid)
        else:
            supplements.append(fid)
    main_links = "\n".join(
        f'<li><a href="#{figure_anchor(fid)}">{e(display_label(fid))}</a></li>'
        for fid in sorted(main, key=figure_sort_key)
    )
    supplement_blocks = "\n".join(
        figure_block(fid, figures, legends, compact=True)
        for fid in sorted(supplements, key=figure_sort_key)
    )
    return f"""
<section id="figures">
  <h2>Figure Index</h2>
  <p>Main figures are cited in the Results. Supplemental figures provide additional raw-data, quality-control, model-diagnostic, posterior, and sensitivity analyses.</p>
  <ul class="main-figure-links">{main_links}</ul>
  <h3>Supplemental Figures</h3>
  <div class="supplement-list">{supplement_blocks}</div>
</section>
""".strip()


def render_audit_trail(result_sections: list[dict[str, object]], methods_path: str) -> str:
    sources = [
        Path("rebuild/manuscript/build_manuscript_html.py"),
        Path("source/figure_set_manifest.csv"),
        Path("source/integrated_figure_legends.md"),
        Path("traceability/upstream_input_register.tsv"),
        Path("review_state/accepted_exceptions.md"),
        Path("evidence/claim_graph/claim_graph_integrated.json"),
        Path("source/manuscript_sections/abstract.md"),
        Path("source/manuscript_sections/introduction.md"),
        Path("source/manuscript_sections/discussion.md"),
        Path(methods_path),
    ]
    sources.extend(Path(str(section["sidecar_path"])) for section in result_sections)
    source_items = "\n".join(f"<li><code>{e(path)}</code></li>" for path in sources)
    return f"""
<section id="audit-trail">
  <h2>Audit Trail</h2>
  <p>This review copy was generated mechanically from the accepted current figure, claim, Results, Methods, and legend inputs. Provenance, feedback, and validation notes are kept out of the journal-facing body and retained in the files below.</p>
  <ul class="provenance-list">{source_items}</ul>
</section>
""".strip()


def render_html() -> str:
    figures = parse_manifest()
    legends = parse_legends()
    write_legend_validation_report(figures, legends)
    result_sections = load_result_sections()
    abstract = load_manuscript_section("abstract")
    introduction = load_manuscript_section("introduction")
    discussion = load_manuscript_section("discussion")
    references = load_manuscript_section("references")
    methods = load_methods_section()
    methods_has_heading = bool(re.match(r"^\s*#{1,2}\s+Methods\s*$", methods["text"].splitlines()[0]))
    methods_heading = "" if methods_has_heading else f'        <h2>{e(methods["title"])}</h2>'

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{e(TITLE)}</title>
  <style>
    :root {{
      --ink: #1b1f23;
      --muted: #5d6874;
      --line: #d7dde3;
      --paper: #ffffff;
      --wash: #f5f7f9;
      --accent: #1f6f78;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: var(--paper);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 17px;
      line-height: 1.58;
    }}
    a {{ color: var(--accent); text-decoration-thickness: 0.08em; }}
    .page {{
      display: grid;
      grid-template-columns: minmax(190px, 250px) minmax(0, 1fr);
      gap: 32px;
      max-width: 1320px;
      margin: 0 auto;
      padding: 32px 28px 72px;
    }}
    header {{
      grid-column: 1 / -1;
      border-bottom: 1px solid var(--line);
      padding: 36px 0 28px;
    }}
    .kicker {{
      color: var(--accent);
      font-size: 0.82rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      margin: 0 0 12px;
    }}
    h1 {{
      margin: 0;
      max-width: 980px;
      font-size: clamp(2.1rem, 4.6vw, 4.2rem);
      line-height: 1.04;
      font-weight: 760;
      letter-spacing: 0;
    }}
    h2 {{
      margin: 0 0 18px;
      font-size: 1.75rem;
      line-height: 1.18;
      letter-spacing: 0;
    }}
    h3 {{
      margin: 0 0 14px;
      font-size: 1.28rem;
      line-height: 1.25;
      letter-spacing: 0;
    }}
    p {{ max-width: 900px; margin: 0 0 16px; }}
    .toc {{
      position: sticky;
      top: 20px;
      align-self: start;
      padding: 16px 0;
      font-size: 0.92rem;
    }}
    .toc strong {{
      display: block;
      color: var(--muted);
      margin: 0 0 8px;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .toc a {{
      display: block;
      padding: 6px 0;
      color: var(--ink);
      text-decoration: none;
      border-bottom: 1px solid transparent;
    }}
    .toc a:hover {{ color: var(--accent); border-bottom-color: var(--line); }}
    main {{ min-width: 0; }}
    section {{
      padding: 34px 0;
      border-bottom: 1px solid var(--line);
    }}
    .figure-card {{
      margin: 28px 0 8px;
      border: 1px solid var(--line);
      border-radius: 6px;
      overflow: hidden;
      background: #fff;
    }}
    .figure-card img {{
      display: block;
      width: 100%;
      max-height: 86vh;
      object-fit: contain;
      background: white;
      border-bottom: 1px solid var(--line);
    }}
    .figure-card figcaption {{
      padding: 14px 16px 16px;
      color: var(--muted);
      font-size: 0.9rem;
      line-height: 1.45;
    }}
    .figure-card figcaption strong {{
      display: block;
      color: var(--ink);
      margin-bottom: 8px;
      font-size: 0.98rem;
    }}
    .legend-body p {{ max-width: none; margin: 0 0 10px; }}
    .figure-meta {{
      max-width: none;
      margin: 10px 0 0;
      color: var(--muted);
      font-size: 0.78rem;
    }}
    .figure-card.compact {{ margin: 0; }}
    .figure-card.compact img {{ max-height: 620px; }}
    .supplement-list {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(310px, 1fr));
      gap: 18px;
      margin-top: 18px;
    }}
    .main-figure-links {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px 14px;
      margin: 0 0 24px;
      padding: 0;
      list-style: none;
    }}
    .main-figure-links a {{
      display: inline-block;
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 5px 9px;
      text-decoration: none;
      background: var(--wash);
    }}
    code {{
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace;
      font-size: 0.88em;
      background: var(--wash);
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 0.08rem 0.26rem;
    }}
    .math-block {{
      max-width: 900px;
      margin: 8px 0 18px;
      padding-left: 18px;
      border-left: 3px solid var(--line);
      font-family: Georgia, "Times New Roman", serif;
      font-size: 1.04rem;
      font-style: italic;
      line-height: 1.45;
    }}
    .math-inline {{
      font-family: Georgia, "Times New Roman", serif;
      font-style: italic;
      white-space: nowrap;
    }}
    .math-fn,
    .math-op {{
      font-style: normal;
    }}
    .manuscript-list {{
      max-width: 900px;
      margin: 0 0 16px 22px;
      padding: 0;
    }}
    .manuscript-list li {{ margin: 0 0 8px; }}
    .table-wrap {{ max-width: 100%; overflow-x: auto; margin: 0 0 20px; }}
    table {{ border-collapse: collapse; min-width: 680px; font-size: 0.88rem; }}
    th, td {{ border: 1px solid var(--line); padding: 7px 9px; text-align: left; vertical-align: top; }}
    th {{ background: var(--wash); }}
    .provenance-list {{
      margin: 0 0 0 20px;
      padding: 0;
      color: var(--muted);
      font-size: 0.95rem;
    }}
    @media (max-width: 900px) {{
      .page {{
        display: block;
        padding: 22px 18px 54px;
      }}
      .toc {{
        position: static;
        border-bottom: 1px solid var(--line);
        margin-bottom: 16px;
      }}
      h1 {{ font-size: 2.2rem; }}
      .supplement-list {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <header>
      <p class="kicker">Manuscript draft</p>
      <h1>{e(TITLE)}</h1>
    </header>
    {render_toc(result_sections)}
    <main>
      <section id="abstract">
        <h2>{e(abstract["title"])}</h2>
        {markdown_blocks(abstract["text"], set(figures))}
      </section>
      <section id="introduction">
        <h2>{e(introduction["title"])}</h2>
        {markdown_blocks(introduction["text"], set(figures))}
      </section>
      <section id="results">
        <h2>Results</h2>
      </section>
      {render_results(figures, legends, result_sections)}
      <section id="discussion">
        <h2>{e(discussion["title"])}</h2>
        {markdown_blocks(discussion["text"], set(figures))}
      </section>
      <section id="methods">
{methods_heading}
        {markdown_blocks(methods["text"], set(figures), enable_math=True, extended_markdown=True)}
      </section>
      <section id="references">
        <h2>{e(references["title"])}</h2>
        {markdown_blocks(references["text"], set(figures), extended_markdown=True)}
      </section>
      {render_figure_index(figures, legends)}
      {render_audit_trail(result_sections, methods["path"])}
    </main>
  </div>
</body>
</html>
"""


def validate_output(html_text: str, figures: dict[str, dict[str, object]]) -> None:
    srcs = re.findall(r'<img src="([^"]+)"', html_text)
    if len(srcs) != len(figures):
        raise ValueError(f"Rendered {len(srcs)} images, expected {len(figures)}")
    non_embedded = [src for src in srcs if not src.startswith("data:image/png;base64,")]
    if non_embedded:
        raise ValueError("Non-embedded image sources: " + ", ".join(non_embedded))
    for src in srcs:
        base64.b64decode(src.removeprefix("data:image/png;base64,"), validate=True)
    anchors = re.findall(r'id="([^"]+)"', html_text)
    duplicates = sorted({anchor for anchor in anchors if anchors.count(anchor) > 1})
    if duplicates:
        raise ValueError("Duplicate HTML anchors: " + ", ".join(duplicates))


def main() -> None:
    figures = parse_manifest()
    html_text = render_html()
    validate_output(html_text, figures)
    OUT.write_text(html_text)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
