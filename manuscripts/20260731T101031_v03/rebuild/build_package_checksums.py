#!/usr/bin/env python3
"""Write the non-circular v03 package file checksum inventory."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "traceability/package_file_checksums.sha256"
EXCLUDED = {
    OUT,
    ROOT / "validation/assembly_validation.json",
    ROOT / "validation/assembly_validation_report.md",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file() and path not in EXCLUDED and "__pycache__" not in path.parts
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{sha256(path)}  {path.relative_to(ROOT)}" for path in files]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
