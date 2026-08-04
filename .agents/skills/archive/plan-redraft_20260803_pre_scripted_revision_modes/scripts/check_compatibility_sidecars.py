#!/usr/bin/env python3
"""Validate and compare consumer-authored manuscript compatibility sidecars."""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List


FACT_STATUSES = {"validated", "unresolved"}
REQUIRED_CONSUMERS = {
    "analysis",
    "manuscript-figure-workflow",
    "claim-graph-integration",
    "results-text",
    "method-table-provenance",
    "manuscript-legend-writing",
    "serve-manuscript-abstract-introduction-discussion",
}


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise ValueError("missing file: %s" % path)
    except json.JSONDecodeError as exc:
        raise ValueError("invalid JSON %s: %s" % (path, exc))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def canonical_assertion(fact: Dict[str, Any]) -> str:
    payload = {
        "value": fact.get("value"),
        "unit": fact.get("unit", ""),
        "scope": fact.get("scope", ""),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def validate_source_artifacts(
    repo_root: Path,
    consumer: str,
    fact_id: str,
    entries: Any,
    errors: List[str],
) -> None:
    if not isinstance(entries, list) or not entries:
        errors.append("%s:%s requires nonempty source_artifacts" % (consumer, fact_id))
        return
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or not str(entry.get("path", "")).strip():
            errors.append("%s:%s source_artifacts[%d] lacks path" % (consumer, fact_id, index))
            continue
        path = resolve(repo_root, str(entry["path"]))
        if not path.exists():
            errors.append("%s:%s source artifact is missing: %s" % (consumer, fact_id, path))
            continue
        expected = str(entry.get("sha256", "")).strip()
        if expected:
            if not path.is_file():
                errors.append("%s:%s checksummed source is not a file: %s" % (consumer, fact_id, path))
            elif sha256(path) != expected:
                errors.append("%s:%s source checksum mismatch: %s" % (consumer, fact_id, path))


def markdown_report(report: Dict[str, Any]) -> str:
    lines = [
        "# Pre-Assembly Compatibility Gate",
        "",
        "Status: **%s**" % report["status"],
        "",
        "- Sidecars checked: %d" % report["sidecars_checked"],
        "- Facts checked: %d" % report["facts_checked"],
        "- Conflicts: %d" % len(report["conflicts"]),
        "- Unresolved facts: %d" % len(report["unresolved_facts"]),
        "- Accepted exceptions: %d" % len(report["accepted_exceptions_used"]),
        "- Errors: %d" % len(report["errors"]),
    ]
    for title, key in (
        ("Errors", "errors"),
        ("Conflicts", "conflicts"),
        ("Unresolved Facts", "unresolved_facts"),
        ("Accepted Exceptions", "accepted_exceptions_used"),
    ):
        entries = report[key]
        if not entries:
            continue
        lines.extend(["", "## %s" % title, ""])
        for entry in entries:
            if isinstance(entry, str):
                lines.append("- %s" % entry)
            else:
                lines.append("- `%s`: %s" % (entry.get("fact_id", "unknown"), entry.get("summary", "")))
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--report-json", required=True)
    parser.add_argument("--report-md", required=True)
    parser.add_argument("--repo-root", default=".")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    registry_path = resolve(repo_root, args.registry).resolve()
    errors = []  # type: List[str]
    conflicts = []  # type: List[Dict[str, Any]]
    unresolved = []  # type: List[Dict[str, Any]]
    accepted_used = []  # type: List[Dict[str, Any]]
    assertions = {}  # type: Dict[str, List[Dict[str, Any]]]
    sidecars_checked = 0
    facts_checked = 0

    try:
        registry = read_json(registry_path)
    except ValueError as exc:
        registry = {}
        errors.append(str(exc))
    if not isinstance(registry, dict) or registry.get("schema_version") != 1:
        errors.append("registry schema_version must be 1")

    accepted = {}  # type: Dict[str, Dict[str, Any]]
    accepted_entries = registry.get("accepted_exceptions", []) if isinstance(registry, dict) else []
    if not isinstance(accepted_entries, list):
        errors.append("accepted_exceptions must be a list")
        accepted_entries = []
    for entry in accepted_entries:
        if not isinstance(entry, dict):
            errors.append("accepted exception entries must be objects")
            continue
        fact_id = str(entry.get("fact_id", "")).strip()
        reason = str(entry.get("reason", "")).strip()
        authorized_by = str(entry.get("authorized_by", "")).strip()
        if not fact_id or not reason or not authorized_by:
            errors.append("accepted exceptions require fact_id, reason, and authorized_by")
            continue
        accepted[fact_id] = entry

    sidecar_entries = registry.get("sidecars", []) if isinstance(registry, dict) else []
    if not isinstance(sidecar_entries, list) or not sidecar_entries:
        errors.append("registry requires a nonempty sidecars list")
        sidecar_entries = []

    seen_consumers = set()
    for entry in sidecar_entries:
        if not isinstance(entry, dict):
            errors.append("sidecar registry entries must be objects")
            continue
        expected_consumer = str(entry.get("consumer", "")).strip()
        path_text = str(entry.get("path", "")).strip()
        if not expected_consumer or not path_text:
            errors.append("sidecar registry entries require consumer and path")
            continue
        if expected_consumer in seen_consumers:
            errors.append("duplicate sidecar registry consumer: %s" % expected_consumer)
            continue
        seen_consumers.add(expected_consumer)
        sidecar_path = resolve(repo_root, path_text).resolve()
        try:
            sidecar = read_json(sidecar_path)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        sidecars_checked += 1
        if not isinstance(sidecar, dict) or sidecar.get("schema_version") != 1:
            errors.append("%s sidecar schema_version must be 1" % expected_consumer)
            continue
        consumer = str(sidecar.get("consumer", "")).strip()
        if consumer != expected_consumer:
            errors.append(
                "sidecar consumer mismatch for %s: found %s" % (expected_consumer, consumer or "missing")
            )
        if sidecar.get("completion_attestation") is not True:
            errors.append("%s sidecar lacks completion_attestation=true" % expected_consumer)
        if not str(sidecar.get("bundle_id", "")).strip():
            errors.append("%s sidecar lacks bundle_id" % expected_consumer)
        facts = sidecar.get("facts")
        if not isinstance(facts, list):
            errors.append("%s sidecar facts must be a list" % expected_consumer)
            continue
        local_ids = set()
        for fact in facts:
            facts_checked += 1
            if not isinstance(fact, dict):
                errors.append("%s sidecar fact entries must be objects" % expected_consumer)
                continue
            fact_id = str(fact.get("fact_id", "")).strip()
            if not fact_id:
                errors.append("%s sidecar fact lacks fact_id" % expected_consumer)
                continue
            if fact_id in local_ids:
                errors.append("%s sidecar repeats fact_id %s" % (expected_consumer, fact_id))
                continue
            local_ids.add(fact_id)
            status = str(fact.get("status", "")).strip()
            if status not in FACT_STATUSES:
                errors.append("%s:%s has invalid status %s" % (expected_consumer, fact_id, status))
            if "value" not in fact:
                errors.append("%s:%s lacks value" % (expected_consumer, fact_id))
            validate_source_artifacts(
                repo_root,
                expected_consumer,
                fact_id,
                fact.get("source_artifacts"),
                errors,
            )
            record = {
                "consumer": expected_consumer,
                "fact_id": fact_id,
                "status": status,
                "assertion": canonical_assertion(fact),
                "value": fact.get("value"),
                "unit": fact.get("unit", ""),
                "scope": fact.get("scope", ""),
            }
            assertions.setdefault(fact_id, []).append(record)
            if status == "unresolved":
                unresolved.append(
                    {
                        "fact_id": fact_id,
                        "summary": "%s reported this fact as unresolved" % expected_consumer,
                        "consumer": expected_consumer,
                    }
                )

    missing_consumers = sorted(REQUIRED_CONSUMERS.difference(seen_consumers))
    if missing_consumers:
        errors.append("registry is missing required consumers: %s" % missing_consumers)
    unexpected_consumers = sorted(seen_consumers.difference(REQUIRED_CONSUMERS))
    if unexpected_consumers:
        errors.append("registry contains unexpected consumers: %s" % unexpected_consumers)

    for fact_id, records in sorted(assertions.items()):
        variants = {record["assertion"] for record in records}
        if len(variants) > 1:
            conflicts.append(
                {
                    "fact_id": fact_id,
                    "summary": "conflicting assertions from %s"
                    % ", ".join(record["consumer"] for record in records),
                    "assertions": records,
                }
            )

    blocking_conflicts = []
    for entry in conflicts:
        if entry["fact_id"] in accepted:
            used = dict(accepted[entry["fact_id"]])
            used["summary"] = "accepted cross-consumer conflict: %s" % used["reason"]
            accepted_used.append(used)
        else:
            blocking_conflicts.append(entry)
    blocking_unresolved = []
    for entry in unresolved:
        if entry["fact_id"] in accepted:
            if not any(used.get("fact_id") == entry["fact_id"] for used in accepted_used):
                used = dict(accepted[entry["fact_id"]])
                used["summary"] = "accepted unresolved fact: %s" % used["reason"]
                accepted_used.append(used)
        else:
            blocking_unresolved.append(entry)

    status = "BLOCKED" if errors or blocking_conflicts or blocking_unresolved else (
        "WARN" if accepted_used else "PASS"
    )
    report = {
        "schema_version": 1,
        "status": status,
        "registry": str(registry_path),
        "sidecars_checked": sidecars_checked,
        "facts_checked": facts_checked,
        "required_consumers": sorted(REQUIRED_CONSUMERS),
        "errors": errors,
        "conflicts": blocking_conflicts,
        "unresolved_facts": blocking_unresolved,
        "accepted_exceptions_used": accepted_used,
    }
    report_json = resolve(repo_root, args.report_json)
    report_md = resolve(repo_root, args.report_md)
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_md.write_text(markdown_report(report), encoding="utf-8")
    return 0 if status in {"PASS", "WARN"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
