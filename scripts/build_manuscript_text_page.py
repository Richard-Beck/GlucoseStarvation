#!/usr/bin/env python3
"""Build a lightweight, figure-free manuscript webpage from the report HTML."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path


TITLE_PATTERN = re.compile(r"<h1>(?P<title>.*?)</h1>", re.DOTALL)
SECTION_PATTERN = re.compile(
    r'<section id="(?P<section_id>[^"]+)">(?P<body>.*?)</section>',
    re.DOTALL,
)
FIGURE_PATTERN = re.compile(r"<figure\b.*?</figure>", re.DOTALL)
FIGURE_LINK_PATTERN = re.compile(
    r'<a href="#fig-[^"]+">(?P<label>.*?)</a>',
    re.DOTALL,
)
REQUIRED_SECTIONS = {
    "abstract",
    "introduction",
    "results",
    "discussion",
    "methods",
    "references",
}


def build_text_page(report_html: str) -> str:
    title_match = TITLE_PATTERN.search(report_html)
    if title_match is None:
        raise ValueError("Could not find the manuscript title")
    title = title_match.group("title").strip()

    sections: list[tuple[str, str]] = []
    observed: set[str] = set()
    for match in SECTION_PATTERN.finditer(report_html):
        section_id = match.group("section_id")
        if section_id not in REQUIRED_SECTIONS and not section_id.startswith("results_"):
            continue
        body = FIGURE_PATTERN.sub("", match.group("body"))
        body = FIGURE_LINK_PATTERN.sub(lambda item: item.group("label"), body)
        sections.append((section_id, body.strip()))
        observed.add(section_id)

    missing = sorted(REQUIRED_SECTIONS - observed)
    if missing:
        raise ValueError("Missing manuscript sections: " + ", ".join(missing))
    if not any(section_id.startswith("results_") for section_id, _ in sections):
        raise ValueError("No Results subsections were found")

    body_html = "\n\n".join(
        f'<section id="{section_id}">\n{body}\n</section>'
        for section_id, body in sections
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title} | Text-only manuscript</title>
  <style>
    :root {{
      --ink: #20231f;
      --muted: #62685f;
      --line: #d8dbd3;
      --paper: #fbfaf6;
      --accent: #1f6f78;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: var(--paper);
      font-family: Georgia, "Times New Roman", serif;
      font-size: 18px;
      line-height: 1.68;
    }}
    .page {{
      width: min(860px, calc(100% - 36px));
      margin: 0 auto;
      padding: 64px 0 96px;
    }}
    header {{ padding-bottom: 34px; border-bottom: 1px solid var(--line); }}
    .format-label {{
      margin: 0 0 14px;
      color: var(--accent);
      font: 700 0.76rem/1.2 ui-monospace, "SFMono-Regular", Consolas, monospace;
      letter-spacing: 0.1em;
      text-transform: uppercase;
    }}
    h1 {{ margin: 0; font-size: clamp(2.2rem, 7vw, 4.8rem); line-height: 1.03; }}
    h2 {{ margin: 0 0 18px; font-size: 2rem; line-height: 1.18; }}
    h3 {{ margin: 0 0 16px; font-size: 1.38rem; line-height: 1.3; }}
    h4 {{ margin: 26px 0 12px; font-size: 1.08rem; }}
    p {{ margin: 0 0 18px; }}
    nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px 18px;
      padding: 22px 0;
      border-bottom: 1px solid var(--line);
      font: 600 0.85rem/1.4 ui-monospace, "SFMono-Regular", Consolas, monospace;
    }}
    a {{ color: var(--accent); }}
    nav a {{ text-decoration: none; }}
    section {{ padding: 42px 0 18px; }}
    section[id^="results_"] {{ padding-top: 30px; }}
    code {{ font-size: 0.88em; }}
    .math-block {{ margin: 14px 0 22px; padding-left: 18px; border-left: 3px solid var(--line); }}
    .table-wrap {{ max-width: 100%; overflow-x: auto; margin: 0 0 24px; }}
    table {{ border-collapse: collapse; min-width: 680px; font-size: 0.84rem; }}
    th, td {{ border: 1px solid var(--line); padding: 7px 9px; text-align: left; vertical-align: top; }}
    th {{ background: #f0f1eb; }}
    @media (max-width: 620px) {{
      .page {{ padding-top: 34px; }}
      body {{ font-size: 17px; }}
      h1 {{ font-size: 2.35rem; }}
      h2 {{ font-size: 1.65rem; }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <header>
      <p class="format-label">Text-only manuscript</p>
      <h1>{title}</h1>
    </header>
    <nav aria-label="Manuscript sections">
      <a href="#abstract">Abstract</a>
      <a href="#introduction">Introduction</a>
      <a href="#results">Results</a>
      <a href="#discussion">Discussion</a>
      <a href="#methods">Methods</a>
      <a href="#references">References</a>
    </nav>
    <main>
{body_html}
    </main>
  </div>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create manuscript-text/index.html from a report bundle."
    )
    parser.add_argument(
        "site_dir",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "latest_manuscript_bundle",
    )
    args = parser.parse_args()
    site_dir = args.site_dir.resolve()
    source = site_dir / "index.html"
    destination_dir = site_dir / "manuscript-text"
    destination = destination_dir / "index.html"

    if not source.is_file():
        raise FileNotFoundError(f"Missing report entry point: {source}")
    output = build_text_page(source.read_text(encoding="utf-8"))
    destination_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=destination_dir,
        prefix=".index.",
        suffix=".html",
        delete=False,
    ) as handle:
        handle.write(output)
        temporary = Path(handle.name)
    os.replace(temporary, destination)
    print(f"Wrote {destination}")


if __name__ == "__main__":
    main()
