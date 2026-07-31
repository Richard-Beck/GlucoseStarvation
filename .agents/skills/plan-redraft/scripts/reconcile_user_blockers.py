#!/usr/bin/env python3
"""Fingerprint, deduplicate, defer, and resolve redraft user-decision blockers."""

import argparse
import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple


SCHEMA_VERSION = 1
CONTRACT_VERSION = "redraft-user-blockers/v1"
MODES = {"interactive", "defer_to_assembly"}
DISPOSITIONS = {"awaiting_user", "deferred_to_assembly", "resolved"}
BLOCKER_RE = re.compile(r"^BLOCKER-[0-9a-f]{64}$")
KEY_RE = re.compile(r"^[A-Za-z0-9_.:-]+$")


def now_text() -> str:
    return dt.datetime.now().astimezone().replace(microsecond=0).isoformat()


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise SystemExit("Missing JSON file: %s" % path)
    except json.JSONDecodeError as exc:
        raise SystemExit("%s: invalid JSON: %s" % (path, exc))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.%s" % os.getpid())
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(str(temporary), str(path))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def resolve_path(repo_root: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def display_path(repo_root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def require_schema(payload: Any, label: str) -> Dict[str, Any]:
    if not isinstance(payload, dict) or payload.get("schema_version") != SCHEMA_VERSION:
        raise SystemExit("%s schema_version must be 1" % label)
    return payload


def load_policy(path: Path) -> Tuple[str, str]:
    if not path.exists():
        return "interactive", "default: missing run policy"
    payload = require_schema(read_json(path), "run policy")
    mode = str(payload.get("user_blocker_mode", "")).strip()
    if mode not in MODES:
        raise SystemExit("run policy user_blocker_mode must be interactive or defer_to_assembly")
    return mode, str(path)


def validate_artifacts(repo_root: Path, entries: Any, context: str) -> List[Dict[str, str]]:
    if not isinstance(entries, list):
        raise SystemExit("%s relevant_artifacts must be a list" % context)
    output = []
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise SystemExit("%s artifact entries must be objects" % context)
        path_text = str(entry.get("path", "")).strip()
        declared = str(entry.get("sha256", "")).strip().lower()
        if not path_text or not re.fullmatch(r"[0-9a-f]{64}", declared):
            raise SystemExit("%s artifacts require path and full sha256" % context)
        path = resolve_path(repo_root, path_text)
        if not path.is_file():
            raise SystemExit("%s artifact is not a regular file: %s" % (context, path))
        actual = sha256_file(path)
        if actual != declared:
            raise SystemExit("%s artifact checksum mismatch: %s" % (context, path))
        display = display_path(repo_root, path)
        if display in seen:
            raise SystemExit("%s repeats artifact path: %s" % (context, display))
        seen.add(display)
        output.append({"path": display, "sha256": actual})
    return sorted(output, key=lambda item: (item["sha256"], item["path"]))


def normalize_sidecar(repo_root: Path, path: Path) -> Tuple[str, str, List[Dict[str, Any]]]:
    payload = require_schema(read_json(path), "blocker sidecar")
    consumer = str(payload.get("consumer", "")).strip()
    bundle_id = str(payload.get("bundle_id", "")).strip()
    if not consumer or not bundle_id:
        raise SystemExit("blocker sidecar requires consumer and bundle_id")
    if payload.get("completion_attestation") is not True:
        raise SystemExit("blocker sidecar requires completion_attestation=true")
    blockers = payload.get("blockers")
    if not isinstance(blockers, list):
        raise SystemExit("blocker sidecar blockers must be a list")
    output = []
    seen_keys = set()
    for entry in blockers:
        if not isinstance(entry, dict):
            raise SystemExit("blocker entries must be objects")
        key = str(entry.get("blocker_key", "")).strip()
        question = str(entry.get("question", "")).strip()
        scope = str(entry.get("blocking_scope", "")).strip()
        provisional = str(entry.get("provisional_handling", "")).strip()
        if not key or not KEY_RE.fullmatch(key):
            raise SystemExit("blocker_key must use letters, digits, dot, colon, underscore, or hyphen")
        if key in seen_keys:
            raise SystemExit("blocker sidecar repeats blocker_key %s" % key)
        seen_keys.add(key)
        if not question or not provisional:
            raise SystemExit("blocker %s requires question and provisional_handling" % key)
        if scope != "user-decision":
            raise SystemExit("blocker %s blocking_scope must be user-decision" % key)
        options = entry.get("options", [])
        if not isinstance(options, list) or not all(
            isinstance(option, str) and option.strip() for option in options
        ):
            raise SystemExit("blocker %s options must be nonempty strings" % key)
        artifacts = validate_artifacts(
            repo_root, entry.get("relevant_artifacts"), "blocker %s" % key
        )
        fingerprint_payload = {
            "contract_version": CONTRACT_VERSION,
            "consumer": consumer,
            "blocker_key": key,
            "artifact_sha256": sorted(item["sha256"] for item in artifacts),
        }
        fingerprint = "BLOCKER-" + hashlib.sha256(
            canonical_json(fingerprint_payload)
        ).hexdigest()
        output.append(
            {
                "fingerprint": fingerprint,
                "consumer": consumer,
                "blocker_key": key,
                "question": question,
                "blocking_scope": scope,
                "provisional_handling": provisional,
                "options": options,
                "relevant_artifacts": artifacts,
                "bundle_id": bundle_id,
                "source_sidecar": display_path(repo_root, path),
            }
        )
    return consumer, bundle_id, output


def load_ledger(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {
            "schema_version": SCHEMA_VERSION,
            "contract_version": CONTRACT_VERSION,
            "blockers": [],
        }
    payload = require_schema(read_json(path), "blocker ledger")
    if payload.get("contract_version") != CONTRACT_VERSION:
        raise SystemExit("blocker ledger contract_version mismatch")
    blockers = payload.get("blockers")
    if not isinstance(blockers, list):
        raise SystemExit("blocker ledger blockers must be a list")
    seen = set()
    for record in blockers:
        if not isinstance(record, dict):
            raise SystemExit("blocker ledger entries must be objects")
        fingerprint = str(record.get("fingerprint", ""))
        if not BLOCKER_RE.fullmatch(fingerprint) or fingerprint in seen:
            raise SystemExit("blocker ledger has invalid or duplicate fingerprint")
        seen.add(fingerprint)
        if record.get("disposition") not in DISPOSITIONS:
            raise SystemExit("blocker ledger has invalid disposition for %s" % fingerprint)
    return payload


def action_payload(ledger: Dict[str, Any], ledger_path: Path) -> Dict[str, Any]:
    awaiting = [
        record
        for record in ledger["blockers"]
        if record.get("disposition") == "awaiting_user"
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "awaiting_user" if awaiting else "clear",
        "generated_at": now_text(),
        "ledger": str(ledger_path.resolve()),
        "blocker_fingerprints": sorted(record["fingerprint"] for record in awaiting),
        "blockers": [
            {
                "fingerprint": record["fingerprint"],
                "consumer": record["consumer"],
                "question": record["question"],
                "options": record.get("options", []),
                "provisional_handling": record["provisional_handling"],
                "relevant_artifacts": record["relevant_artifacts"],
            }
            for record in sorted(awaiting, key=lambda item: item["fingerprint"])
        ],
    }


def markdown_report(report: Dict[str, Any]) -> str:
    lines = [
        "# User-Decision Blocker Reconciliation",
        "",
        "Status: **%s**" % report["status"],
        "",
        "- Mode: `%s`" % report["user_blocker_mode"],
        "- Current blockers: %d" % report["current_blocker_count"],
        "- New fingerprints: %d" % len(report["new_fingerprints"]),
        "- Repeated fingerprints: %d" % len(report["repeated_fingerprints"]),
        "- Awaiting user: %d" % len(report["awaiting_user_fingerprints"]),
        "- Deferred to Assembly: %d" % len(report["deferred_fingerprints"]),
    ]
    for record in report["current_blockers"]:
        lines.extend(
            [
                "",
                "## %s" % record["fingerprint"],
                "",
                "- Consumer: `%s`" % record["consumer"],
                "- Disposition: `%s`" % record["disposition"],
                "- Question: %s" % record["question"],
                "- Provisional handling: %s" % record["provisional_handling"],
            ]
        )
    lines.append("")
    return "\n".join(lines)


def reconcile(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    sidecar_path = resolve_path(repo_root, args.sidecar)
    ledger_path = resolve_path(repo_root, args.ledger)
    action_path = resolve_path(repo_root, args.user_action_file)
    policy_path = resolve_path(repo_root, args.run_policy)
    mode, policy_source = load_policy(policy_path)
    consumer, bundle_id, current = normalize_sidecar(repo_root, sidecar_path)
    ledger = load_ledger(ledger_path)
    by_fingerprint = {record["fingerprint"]: record for record in ledger["blockers"]}
    stamp = now_text()
    if mode == "defer_to_assembly":
        for record in by_fingerprint.values():
            if record["disposition"] == "awaiting_user":
                record["disposition"] = "deferred_to_assembly"
                record["deferred_at"] = stamp
                record["deferred_by"] = "run_policy"
    new = []
    repeated = []
    current_records = []
    for blocker in current:
        fingerprint = blocker["fingerprint"]
        record = by_fingerprint.get(fingerprint)
        if record is None:
            disposition = (
                "deferred_to_assembly" if mode == "defer_to_assembly" else "awaiting_user"
            )
            record = dict(blocker)
            record.update(
                {
                    "first_seen_at": stamp,
                    "last_seen_at": stamp,
                    "occurrence_count": 1,
                    "disposition": disposition,
                }
            )
            if disposition == "deferred_to_assembly":
                record["deferred_at"] = stamp
                record["deferred_by"] = "run_policy"
            by_fingerprint[fingerprint] = record
            new.append(fingerprint)
        else:
            repeated.append(fingerprint)
            record.update(blocker)
            record["last_seen_at"] = stamp
            record["occurrence_count"] = int(record.get("occurrence_count", 0)) + 1
            if mode == "defer_to_assembly" and record["disposition"] == "awaiting_user":
                record["disposition"] = "deferred_to_assembly"
                record["deferred_at"] = stamp
                record["deferred_by"] = "run_policy"
        current_records.append(record)
    ledger["blockers"] = sorted(by_fingerprint.values(), key=lambda item: item["fingerprint"])
    ledger["updated_at"] = stamp
    write_json(ledger_path, ledger)
    action = action_payload(ledger, ledger_path)
    write_json(action_path, action)
    awaiting = [
        record["fingerprint"]
        for record in ledger["blockers"]
        if record["disposition"] == "awaiting_user"
    ]
    deferred = [
        record["fingerprint"]
        for record in ledger["blockers"]
        if record["disposition"] == "deferred_to_assembly"
    ]
    report = {
        "schema_version": SCHEMA_VERSION,
        "status": "USER_ACTION_REQUIRED" if awaiting else "CONTINUE",
        "user_blocker_mode": mode,
        "policy_source": policy_source,
        "consumer": consumer,
        "bundle_id": bundle_id,
        "sidecar": display_path(repo_root, sidecar_path),
        "ledger": display_path(repo_root, ledger_path),
        "user_action_file": display_path(repo_root, action_path),
        "current_blocker_count": len(current),
        "new_fingerprints": sorted(new),
        "repeated_fingerprints": sorted(repeated),
        "awaiting_user_fingerprints": sorted(awaiting),
        "deferred_fingerprints": sorted(deferred),
        "current_blockers": [
            {
                "fingerprint": record["fingerprint"],
                "consumer": record["consumer"],
                "question": record["question"],
                "provisional_handling": record["provisional_handling"],
                "disposition": record["disposition"],
            }
            for record in current_records
        ],
    }
    report_json = resolve_path(repo_root, args.report_json)
    report_md = resolve_path(repo_root, args.report_md)
    write_json(report_json, report)
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_md.write_text(markdown_report(report), encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def fingerprint_command(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    sidecar_path = resolve_path(repo_root, args.sidecar)
    consumer, bundle_id, blockers = normalize_sidecar(repo_root, sidecar_path)
    output = {
        "schema_version": SCHEMA_VERSION,
        "contract_version": CONTRACT_VERSION,
        "consumer": consumer,
        "bundle_id": bundle_id,
        "sidecar": display_path(repo_root, sidecar_path),
        "blockers": [
            {
                "blocker_key": blocker["blocker_key"],
                "fingerprint": blocker["fingerprint"],
            }
            for blocker in blockers
        ],
    }
    output_path = resolve_path(repo_root, args.output)
    write_json(output_path, output)
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


def resolve_command(args: argparse.Namespace) -> int:
    ledger_path = Path(args.ledger).resolve()
    action_path = Path(args.user_action_file).resolve()
    ledger = load_ledger(ledger_path)
    payload = require_schema(read_json(Path(args.resolution_file).resolve()), "resolution file")
    resolutions = payload.get("resolutions")
    if not isinstance(resolutions, list) or not resolutions:
        raise SystemExit("resolution file requires a nonempty resolutions list")
    records = {record["fingerprint"]: record for record in ledger["blockers"]}
    stamp = now_text()
    changed = []
    for resolution in resolutions:
        if not isinstance(resolution, dict):
            raise SystemExit("resolution entries must be objects")
        fingerprint = str(resolution.get("fingerprint", ""))
        disposition = str(resolution.get("disposition", ""))
        text = str(resolution.get("resolution", "")).strip()
        if not BLOCKER_RE.fullmatch(fingerprint) or fingerprint not in records:
            raise SystemExit("resolution references an unknown fingerprint: %s" % fingerprint)
        if resolution.get("authorized_by") != "user":
            raise SystemExit("resolution entries require authorized_by: user")
        if disposition not in {"resolved", "deferred_to_assembly"} or not text:
            raise SystemExit("resolutions require disposition resolved/deferred_to_assembly and text")
        record = records[fingerprint]
        record["disposition"] = disposition
        record["resolution"] = text
        record["resolved_at"] = stamp
        record["resolved_by"] = "user"
        changed.append(fingerprint)
    ledger["updated_at"] = stamp
    write_json(ledger_path, ledger)
    write_json(action_path, action_payload(ledger, ledger_path))
    summary = {
        "schema_version": SCHEMA_VERSION,
        "updated_fingerprints": sorted(changed),
        "remaining_awaiting_user": sorted(
            record["fingerprint"]
            for record in ledger["blockers"]
            if record["disposition"] == "awaiting_user"
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    fingerprint_parser = subparsers.add_parser("fingerprint")
    fingerprint_parser.add_argument("--sidecar", required=True)
    fingerprint_parser.add_argument("--output", required=True)
    fingerprint_parser.add_argument("--repo-root", default=".")
    reconcile_parser = subparsers.add_parser("reconcile")
    reconcile_parser.add_argument("--sidecar", required=True)
    reconcile_parser.add_argument("--ledger", required=True)
    reconcile_parser.add_argument("--run-policy", required=True)
    reconcile_parser.add_argument("--user-action-file", required=True)
    reconcile_parser.add_argument("--report-json", required=True)
    reconcile_parser.add_argument("--report-md", required=True)
    reconcile_parser.add_argument("--repo-root", default=".")
    resolve_parser = subparsers.add_parser("resolve")
    resolve_parser.add_argument("--ledger", required=True)
    resolve_parser.add_argument("--resolution-file", required=True)
    resolve_parser.add_argument("--user-action-file", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "fingerprint":
        return fingerprint_command(args)
    return reconcile(args) if args.command == "reconcile" else resolve_command(args)


if __name__ == "__main__":
    raise SystemExit(main())
