#!/usr/bin/env python3
"""Build and validate the v03 ordered A/I/D serving bundle."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


OUTPUT_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = OUTPUT_ROOT.parents[2]
PREVIOUS_OUTPUT_ROOT = (
    PROJECT_ROOT
    / "manuscripts/20260713T163237_v02/evidence/abstract_introduction_discussion"
)
CONTEXT_SEPARATOR = b"\n\n"
CONTEXT_SOURCES = [
    {
        "path": PROJECT_ROOT
        / "agent-dev/manuscript_results/20260731_v03_results_revision/combined_results_preview.md",
        "role": "approved_results",
    },
    {
        "path": PROJECT_ROOT
        / "agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction/methods_text.md",
        "role": "canonical_methods",
    },
    {
        "path": PROJECT_ROOT
        / "agent-dev/manuscript_integration/20260729_v03_figure_set_integration/integrated_figure_legends.md",
        "role": "approved_figure_legends",
    },
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def payload_record(path: Path, role: str) -> dict[str, object]:
    data = path.read_bytes()
    return {
        "path": path.name,
        "role": role,
        "byte_size": len(data),
        "sha256": sha256_bytes(data),
    }


def bundle_id(payload_files: list[dict[str, object]]) -> str:
    identity = json.dumps(
        payload_files, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return f"aid-bundle-sha256:{sha256_bytes(identity)}"


def ordered_bundle_dirs(root: Path) -> list[Path]:
    current = root / "current_bundle"
    if not current.is_dir():
        return []
    old_dirs: dict[int, Path] = {}
    for candidate in root.glob("old_bundle_*"):
        match = re.fullmatch(r"old_bundle_(\d+)", candidate.name)
        if match and candidate.is_dir():
            old_dirs[int(match.group(1))] = candidate
    if old_dirs and sorted(old_dirs) != list(range(1, max(old_dirs) + 1)):
        return []
    return [current] + [old_dirs[i] for i in sorted(old_dirs)]


def validate_bundle(bundle: Path) -> tuple[bool, list[str], dict | None]:
    failures: list[str] = []
    manifest_path = bundle / "manifest.json"
    context_path = bundle / "context.txt"
    if not context_path.is_file():
        failures.append(f"{bundle.name}: missing context.txt")
    if not manifest_path.is_file():
        failures.append(f"{bundle.name}: missing manifest.json")
        return False, failures, None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        failures.append(f"{bundle.name}: invalid manifest.json: {exc}")
        return False, failures, None
    records = manifest.get("payload_files")
    if not isinstance(records, list):
        failures.append(f"{bundle.name}: payload_files is not a list")
        return False, failures, manifest
    seen: set[str] = set()
    roles: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            failures.append(f"{bundle.name}: malformed payload record")
            continue
        rel = record.get("path")
        role = record.get("role")
        if not isinstance(rel, str) or Path(rel).name != rel:
            failures.append(f"{bundle.name}: invalid payload path {rel!r}")
            continue
        if rel in seen:
            failures.append(f"{bundle.name}: duplicate payload path {rel}")
        seen.add(rel)
        if isinstance(role, str):
            roles.add(role)
        payload = bundle / rel
        if not payload.is_file():
            failures.append(f"{bundle.name}: missing payload {rel}")
            continue
        data = payload.read_bytes()
        if record.get("byte_size") != len(data):
            failures.append(f"{bundle.name}: byte-size mismatch for {rel}")
        if record.get("sha256") != sha256_bytes(data):
            failures.append(f"{bundle.name}: hash mismatch for {rel}")
    if "context.txt" not in seen or "context" not in roles:
        failures.append(f"{bundle.name}: context is not registered")
    for optional_name, optional_role in (
        ("literature_map.json", "literature_map"),
        ("aid_source.md", "aid_source"),
    ):
        if (bundle / optional_name).exists() and (
            optional_name not in seen or optional_role not in roles
        ):
            failures.append(f"{bundle.name}: orphan {optional_name}")
    bundle_value = manifest.get("bundle_id")
    if not isinstance(bundle_value, str) or not bundle_value.startswith(
        "aid-bundle-sha256:"
    ):
        failures.append(f"{bundle.name}: invalid bundle_id")
    return not failures, failures, manifest


def extract_sections(source: bytes) -> dict[str, bytes]:
    specifications = {
        "abstract.md": "Abstract",
        "introduction.md": "Introduction",
        "discussion.md": "Discussion",
        "references.md": "References cited",
    }
    output: dict[str, bytes] = {}
    for filename, heading in specifications.items():
        start_match = re.search(
            rb"(?m)^## " + re.escape(heading.encode("utf-8")) + rb"\r?$", source
        )
        if not start_match:
            if filename == "references.md":
                continue
            raise ValueError(f"missing required A/I/D heading: {heading}")
        next_match = re.search(rb"(?m)^## [^\r\n]+\r?$", source[start_match.end() :])
        end = len(source)
        if next_match:
            end = start_match.end() + next_match.start()
        # Preserve the complete heading-delimited byte range, including the
        # blank line immediately preceding the next level-two heading.
        output[filename] = source[start_match.start() : end]
    return output


def validate_output(root: Path) -> list[str]:
    failures: list[str] = []
    bundles = ordered_bundle_dirs(root)
    if not bundles:
        failures.append("bundle ordering is missing, incomplete, or non-unique")
        return failures
    manifests: list[dict] = []
    for bundle in bundles:
        valid, bundle_failures, manifest = validate_bundle(bundle)
        failures.extend(bundle_failures)
        if valid and manifest is not None:
            manifests.append(manifest)
    ids = [manifest.get("bundle_id") for manifest in manifests]
    if len(ids) != len(set(ids)):
        failures.append("bundle IDs are not unique")

    served_root = root / "served"
    served_manifest_path = served_root / "manifest.json"
    if not served_manifest_path.is_file():
        failures.append("missing served/manifest.json")
        return failures
    try:
        served_manifest = json.loads(served_manifest_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        failures.append(f"invalid served manifest: {exc}")
        return failures
    served_records = served_manifest.get("served_files", [])
    for record in served_records:
        path = served_root / record["path"]
        if not path.is_file():
            failures.append(f"missing served file {record['path']}")
            continue
        data = path.read_bytes()
        if record.get("byte_size") != len(data):
            failures.append(f"served byte-size mismatch for {record['path']}")
        if record.get("sha256") != sha256_bytes(data):
            failures.append(f"served hash mismatch for {record['path']}")

    source_rel = served_manifest.get("source_path")
    if source_rel:
        source_path = (served_root / source_rel).resolve()
        if not source_path.is_file():
            failures.append("served source_path does not resolve")
        else:
            expected = extract_sections(source_path.read_bytes())
            for filename, data in expected.items():
                served_path = served_root / filename
                if not served_path.is_file() or served_path.read_bytes() != data:
                    failures.append(f"served block is not byte-preserving: {filename}")
    else:
        for filename in ("abstract.md", "introduction.md", "discussion.md"):
            path = served_root / filename
            if not path.is_file() or path.stat().st_size != 0:
                failures.append(f"invalid explicit empty state for {filename}")
        if (served_root / "references.md").exists():
            failures.append("references.md present in explicit empty state")
    return failures


def validate_previous(root: Path) -> list[str]:
    failures = validate_output(root)
    status_path = root / "status.json"
    if not status_path.is_file():
        failures.append("missing status.json")
    else:
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
            if status.get("status") != "pass":
                failures.append("prior status is not pass")
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            failures.append(f"invalid prior status.json: {exc}")
    return failures


def main() -> None:
    generated_targets = [
        OUTPUT_ROOT / "current_bundle",
        OUTPUT_ROOT / "served",
        OUTPUT_ROOT / "status.json",
    ]
    if any(path.exists() for path in generated_targets) or list(
        OUTPUT_ROOT.glob("old_bundle_*")
    ):
        raise SystemExit("refusing to overwrite an existing generated bundle")

    source_records: list[dict[str, object]] = []
    source_data: list[bytes] = []
    for item in CONTEXT_SOURCES:
        path = item["path"]
        if not path.is_file():
            raise SystemExit(f"missing context source: {path}")
        data = path.read_bytes()
        source_data.append(data.rstrip(b"\r\n"))
        source_records.append(
            {
                "path": str(path.relative_to(PROJECT_ROOT)),
                "role": item["role"],
                "byte_size": len(data),
                "sha256": sha256_bytes(data),
            }
        )

    current = OUTPUT_ROOT / "current_bundle"
    current.mkdir(parents=True)
    context = CONTEXT_SEPARATOR.join(source_data) + b"\n"
    (current / "context.txt").write_bytes(context)
    current_payloads = [payload_record(current / "context.txt", "context")]
    current_manifest = {
        "schema_version": 1,
        "bundle_id": bundle_id(current_payloads),
        "hash_algorithm": "sha256",
        "payload_files": current_payloads,
        "context_separator": "one blank line (LF LF)",
        "context_sources": source_records,
    }
    (current / "manifest.json").write_text(
        json.dumps(current_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    previous_failures = validate_previous(PREVIOUS_OUTPUT_ROOT)
    previous_disposition: dict[str, object]
    if previous_failures:
        previous_disposition = {
            "supplied": True,
            "path": str(PREVIOUS_OUTPUT_ROOT.relative_to(PROJECT_ROOT)),
            "disposition": "ignored_invalid",
            "validation_failures": previous_failures,
        }
    else:
        prior_bundles = ordered_bundle_dirs(PREVIOUS_OUTPUT_ROOT)
        for index, source_bundle in enumerate(prior_bundles, start=1):
            shutil.copytree(source_bundle, OUTPUT_ROOT / f"old_bundle_{index}")
        previous_disposition = {
            "supplied": True,
            "path": str(PREVIOUS_OUTPUT_ROOT.relative_to(PROJECT_ROOT)),
            "disposition": "accepted",
            "validation_failures": [],
        }

    bundles = ordered_bundle_dirs(OUTPUT_ROOT)
    served_root = OUTPUT_ROOT / "served"
    served_root.mkdir()
    source_bundle: Path | None = None
    source_manifest: dict | None = None
    for bundle in bundles:
        if (bundle / "aid_source.md").is_file():
            source_bundle = bundle
            source_manifest = json.loads(
                (bundle / "manifest.json").read_text(encoding="utf-8")
            )
            break

    served_records: list[dict[str, object]] = []
    if source_bundle is None:
        for filename, section in (
            ("abstract.md", "Abstract"),
            ("introduction.md", "Introduction"),
            ("discussion.md", "Discussion"),
        ):
            path = served_root / filename
            path.write_bytes(b"")
            served_records.append(
                {
                    "path": filename,
                    "section": section,
                    "byte_size": 0,
                    "sha256": sha256_bytes(b""),
                }
            )
        source_path = None
        source_bundle_id = None
        extraction = "explicit zero-byte empty state"
    else:
        extracted = extract_sections((source_bundle / "aid_source.md").read_bytes())
        section_names = {
            "abstract.md": "Abstract",
            "introduction.md": "Introduction",
            "discussion.md": "Discussion",
            "references.md": "References cited",
        }
        for filename, data in extracted.items():
            path = served_root / filename
            path.write_bytes(data)
            served_records.append(
                {
                    "path": filename,
                    "section": section_names[filename],
                    "byte_size": len(data),
                    "sha256": sha256_bytes(data),
                }
            )
        source_path = f"../{source_bundle.name}/aid_source.md"
        source_bundle_id = source_manifest["bundle_id"]
        extraction = "heading-delimited byte-preserving blocks"

    served_manifest = {
        "schema_version": 1,
        "source_bundle_id": source_bundle_id,
        "source_path": source_path,
        "extraction": extraction,
        "served_files": served_records,
    }
    (served_root / "manifest.json").write_text(
        json.dumps(served_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    validation_failures = validate_output(OUTPUT_ROOT)
    status = {
        "schema_version": 1,
        "status": "pass" if not validation_failures else "fail",
        "generated_at": datetime.now(ZoneInfo("America/New_York")).isoformat(
            timespec="seconds"
        ),
        "current_bundle_id": current_manifest["bundle_id"],
        "bundle_order": [bundle.name for bundle in bundles],
        "previous_output_root": previous_disposition,
        "context_sources": source_records,
        "served_manifest": "served/manifest.json",
        "validation": {
            "bundle_payload_hashes_verified": not validation_failures,
            "served_file_hashes_verified": not validation_failures,
            "served_blocks_byte_preserving": not validation_failures,
            "failures": validation_failures,
        },
    }
    (OUTPUT_ROOT / "status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if validation_failures:
        raise SystemExit("bundle validation failed: " + "; ".join(validation_failures))


if __name__ == "__main__":
    main()
