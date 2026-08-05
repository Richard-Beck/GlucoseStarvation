#!/usr/bin/env python3
"""Extract embedded manuscript PNGs and replace their HTML data URIs."""

from __future__ import annotations

import argparse
import base64
import os
import re
import shutil
import tempfile
from pathlib import Path


IMAGE_PATTERN = re.compile(
    r'(?P<prefix><figure\b[^>]*\bid="(?P<figure_id>[^"]+)"[^>]*>'
    r'.*?<img\b[^>]*\bsrc=")'
    r'data:image/png;base64,(?P<payload>[^"]+)'
    r'(?P<suffix>")',
    re.DOTALL,
)
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def safe_filename(figure_id: str) -> str:
    stem = re.sub(r"[^a-zA-Z0-9._-]+", "-", figure_id).strip("-.")
    if not stem:
        raise ValueError(f"Figure ID cannot form a filename: {figure_id!r}")
    return f"{stem}.png"


def externalize(site_dir: Path) -> int:
    index_path = site_dir / "index.html"
    if not index_path.is_file():
        raise FileNotFoundError(f"Missing report entry point: {index_path}")

    html_text = index_path.read_text(encoding="utf-8")
    image_dir = site_dir / "manuscript-images"
    temporary_image_dir = Path(
        tempfile.mkdtemp(prefix=".manuscript-images.", dir=site_dir)
    )
    temporary_html: Path | None = None
    filenames: set[str] = set()

    try:
        def replace(match: re.Match[str]) -> str:
            filename = safe_filename(match.group("figure_id"))
            if filename in filenames:
                raise ValueError(f"Duplicate manuscript figure filename: {filename}")
            filenames.add(filename)

            payload = re.sub(r"\s+", "", match.group("payload"))
            image = base64.b64decode(payload, validate=True)
            if not image.startswith(PNG_SIGNATURE):
                raise ValueError(f"Embedded figure is not a PNG: {filename}")
            (temporary_image_dir / filename).write_bytes(image)

            relative_path = f"manuscript-images/{filename}"
            return match.group("prefix") + relative_path + match.group("suffix")

        externalized_html, count = IMAGE_PATTERN.subn(replace, html_text)
        if count == 0:
            raise ValueError(f"No embedded manuscript PNGs found in {index_path}")

        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=site_dir,
            prefix=".index.",
            suffix=".html",
            delete=False,
        ) as handle:
            handle.write(externalized_html)
            temporary_html = Path(handle.name)

        if image_dir.exists():
            shutil.rmtree(image_dir)
        os.replace(temporary_image_dir, image_dir)
        os.replace(temporary_html, index_path)
        temporary_html = None
        return count
    finally:
        if temporary_image_dir.exists():
            shutil.rmtree(temporary_image_dir)
        if temporary_html is not None and temporary_html.exists():
            temporary_html.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract embedded PNGs from a manuscript site's index.html."
    )
    parser.add_argument(
        "site_dir",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "latest_manuscript_bundle",
    )
    args = parser.parse_args()
    count = externalize(args.site_dir.resolve())
    print(f"Externalized {count} PNGs in {args.site_dir.resolve()}")


if __name__ == "__main__":
    main()
