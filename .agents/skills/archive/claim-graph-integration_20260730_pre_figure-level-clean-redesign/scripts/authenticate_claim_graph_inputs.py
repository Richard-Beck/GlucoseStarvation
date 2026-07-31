#!/usr/bin/env python3
"""Snapshot and compare exact declared Claim Graph inputs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple


CONTRACT_VERSION = "claim-graph-input-authentication/v1"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit("%s: invalid JSON: %s" % (path, exc))


def display_path(repo_root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def hash_path(path: Path) -> Tuple[str, str, int]:
    if not path.exists() and not path.is_symlink():
        raise SystemExit("Declared input does not exist: %s" % path)
    if path.is_symlink():
        target = os.readlink(path)
        return "symlink", sha256_bytes(target.encode("utf-8")), len(target.encode("utf-8"))
    if path.is_file():
        return "file", sha256_file(path), path.stat().st_size
    if not path.is_dir():
        raise SystemExit("Unsupported declared input type: %s" % path)
    entries = []  # type: List[Dict[str, Any]]
    total_size = 0
    for child in sorted(path.rglob("*"), key=lambda item: str(item.relative_to(path))):
        if child.is_dir() and not child.is_symlink():
            continue
        rel = str(child.relative_to(path))
        kind, digest, size = hash_path(child)
        entries.append({"path": rel, "kind": kind, "sha256": digest, "size": size})
        total_size += size
    return "directory", sha256_bytes(canonical_json(entries)), total_size


def parse_evidence_input(value: str) -> Tuple[str, Path]:
    if "=" not in value:
        raise SystemExit("--evidence-input must use label=path: %s" % value)
    label, path_text = value.split("=", 1)
    label = label.strip()
    path_text = path_text.strip()
    if not label or not path_text:
        raise SystemExit("--evidence-input requires nonempty label and path")
    return label, Path(path_text)


def item_hashes(index_path: Path) -> Tuple[str, Dict[str, str], int]:
    payload = read_json(index_path)
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise SystemExit("%s: item index schema_version must be 1" % index_path)
    items = payload.get("items")
    if not isinstance(items, list):
        raise SystemExit("%s: items must be a list" % index_path)
    hashes = {}  # type: Dict[str, str]
    canonical_items = []
    for item in items:
        if not isinstance(item, dict):
            raise SystemExit("%s: item entries must be objects" % index_path)
        item_id = str(item.get("item_id", "")).strip()
        if not item_id or "canonical_input" not in item:
            raise SystemExit("%s: items require item_id and canonical_input" % index_path)
        if item_id in hashes:
            raise SystemExit("%s: duplicate item_id %s" % (index_path, item_id))
        digest = sha256_bytes(canonical_json(item["canonical_input"]))
        hashes[item_id] = digest
        canonical_items.append({"item_id": item_id, "sha256": digest})
    canonical_items.sort(key=lambda item: item["item_id"])
    return sha256_bytes(canonical_json(canonical_items)), hashes, len(items)


def audit_records(repo_root: Path, index_path: Path, item_ids: set) -> Dict[str, Dict[str, str]]:
    payload = read_json(index_path)
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise SystemExit("%s: audit index schema_version must be 1" % index_path)
    audits = payload.get("audits")
    if not isinstance(audits, list):
        raise SystemExit("%s: audits must be a list" % index_path)
    records = {}  # type: Dict[str, Dict[str, str]]
    for audit in audits:
        if not isinstance(audit, dict):
            raise SystemExit("%s: audit entries must be objects" % index_path)
        item_id = str(audit.get("item_id", "")).strip()
        path_text = str(audit.get("path", "")).strip()
        if not item_id or not path_text:
            raise SystemExit("%s: audits require item_id and path" % index_path)
        if item_id in records:
            raise SystemExit("%s: duplicate audit item_id %s" % (index_path, item_id))
        if item_ids and item_id not in item_ids:
            raise SystemExit("%s: audit item_id is absent from item index: %s" % (index_path, item_id))
        path = Path(path_text)
        if not path.is_absolute():
            path = repo_root / path
        path = path.resolve()
        if not path.is_file():
            raise SystemExit("Audit is not a file: %s" % path)
        digest = sha256_file(path)
        declared_digest = str(audit.get("sha256", "")).strip()
        if declared_digest and declared_digest != digest:
            raise SystemExit("Audit checksum mismatch: %s" % path)
        records[item_id] = {
            "path": display_path(repo_root, path),
            "sha256": digest,
        }
    return dict(sorted(records.items()))


def authenticated_prior_audits(
    repo_root: Path,
    prior: Dict[str, Any],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    valid = {}  # type: Dict[str, Dict[str, str]]
    invalid = {}  # type: Dict[str, str]
    records = prior.get("audit_records") or {}
    if not isinstance(records, dict):
        return valid, {"*": "prior audit_records is not an object"}
    for item_id, record in records.items():
        item_id = str(item_id)
        if not isinstance(record, dict):
            invalid[item_id] = "prior audit record is not an object"
            continue
        path_text = str(record.get("path", "")).strip()
        expected = str(record.get("sha256", "")).strip()
        if not path_text or not expected:
            invalid[item_id] = "prior audit record lacks path or sha256"
            continue
        path = Path(path_text)
        if not path.is_absolute():
            path = repo_root / path
        if not path.is_file():
            invalid[item_id] = "prior audit file is missing"
            continue
        if sha256_file(path) != expected:
            invalid[item_id] = "prior audit checksum mismatch"
            continue
        valid[item_id] = {
            "path": display_path(repo_root, path),
            "sha256": expected,
        }
    return valid, invalid


def snapshot(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    starting_graph = Path(args.starting_graph).resolve()
    if not starting_graph.is_file():
        raise SystemExit("Starting graph is not a file: %s" % starting_graph)
    declared = []  # type: List[Dict[str, Any]]
    seen_labels = set()
    for raw in args.evidence_input:
        label, path = parse_evidence_input(raw)
        if label in seen_labels:
            raise SystemExit("Duplicate evidence-input label: %s" % label)
        seen_labels.add(label)
        path = path.resolve()
        kind, digest, size = hash_path(path)
        declared.append(
            {
                "label": label,
                "path": display_path(repo_root, path),
                "kind": kind,
                "sha256": digest,
                "size": size,
            }
        )
    declared.sort(key=lambda entry: entry["label"])
    evidence_package_hash = sha256_bytes(canonical_json(declared))

    index_record = None
    items = {}  # type: Dict[str, str]
    if args.item_index:
        index_path = Path(args.item_index).resolve()
        index_hash, items, item_count = item_hashes(index_path)
        index_record = {
            "path": display_path(repo_root, index_path),
            "sha256": index_hash,
            "item_count": item_count,
        }

    audits = {}  # type: Dict[str, Dict[str, str]]
    if args.audit_index:
        audits = audit_records(repo_root, Path(args.audit_index).resolve(), set(items))

    result_graph = None
    if args.result_graph:
        result_path = Path(args.result_graph).resolve()
        if not result_path.is_file():
            raise SystemExit("Result graph is not a file: %s" % result_path)
        result_graph = {
            "path": display_path(repo_root, result_path),
            "sha256": sha256_file(result_path),
        }

    whole_payload = {
        "contract_version": CONTRACT_VERSION,
        "evidence_package_sha256": evidence_package_hash,
        "item_index_sha256": index_record["sha256"] if index_record else None,
    }
    output = {
        "schema_version": 1,
        "contract_version": CONTRACT_VERSION,
        "starting_graph": {
            "path": display_path(repo_root, starting_graph),
            "sha256": sha256_file(starting_graph),
        },
        "evidence_inputs": declared,
        "evidence_package_sha256": evidence_package_hash,
        "item_index": index_record,
        "item_hashes": dict(sorted(items.items())),
        "audit_records": audits,
        "whole_evidence_input_sha256": sha256_bytes(canonical_json(whole_payload)),
        "result_graph": result_graph,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def compare(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    current = read_json(Path(args.current))
    prior = read_json(Path(args.prior))
    for label, payload in (("current", current), ("prior", prior)):
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise SystemExit("%s authentication schema_version must be 1" % label)
    version_match = current.get("contract_version") == prior.get("contract_version") == CONTRACT_VERSION
    current_start = current.get("starting_graph", {}).get("sha256")
    prior_start = prior.get("starting_graph", {}).get("sha256")
    prior_result = (prior.get("result_graph") or {}).get("sha256")
    graph_match = bool(version_match and current_start and current_start in {prior_start, prior_result})
    whole_match = bool(
        graph_match
        and current.get("whole_evidence_input_sha256")
        == prior.get("whole_evidence_input_sha256")
    )
    current_items = current.get("item_hashes") or {}
    prior_items = prior.get("item_hashes") or {}
    if not isinstance(current_items, dict) or not isinstance(prior_items, dict):
        raise SystemExit("item_hashes must be objects")
    valid_audits, invalid_audits = authenticated_prior_audits(repo_root, prior)
    input_identical = sorted(
        item_id
        for item_id, digest in current_items.items()
        if graph_match and prior_items.get(item_id) == digest
    )
    for item_id in input_identical:
        if item_id not in valid_audits:
            invalid_audits.setdefault(item_id, "prior audit record is missing")
    reusable = sorted(item_id for item_id in input_identical if item_id in valid_audits)
    changed = sorted(
        item_id
        for item_id, digest in current_items.items()
        if item_id in prior_items and prior_items.get(item_id) != digest
    )
    new = sorted(set(current_items).difference(prior_items))
    removed = sorted(set(prior_items).difference(current_items))
    misses = sorted(set(current_items).difference(reusable))
    audit_complete = set(current_items).issubset(reusable)
    if whole_match and (not current_items or audit_complete):
        mode = "reuse_complete"
    elif graph_match and reusable:
        mode = "reuse_partial"
    else:
        mode = "fresh_full"
        reusable = []
        misses = sorted(current_items)
    report = {
        "schema_version": 1,
        "contract_version": CONTRACT_VERSION,
        "mode": mode,
        "contract_version_match": version_match,
        "starting_graph_match": graph_match,
        "whole_evidence_input_match": whole_match,
        "audit_receipts_complete": audit_complete,
        "reusable_item_ids": reusable,
        "reused_audits": {item_id: valid_audits[item_id] for item_id in reusable},
        "cache_miss_item_ids": misses,
        "invalid_prior_audits": {
            item_id: reason
            for item_id, reason in sorted(invalid_audits.items())
            if item_id in current_items or item_id == "*"
        },
        "changed_item_ids": changed,
        "new_item_ids": new,
        "removed_item_ids": removed,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    snapshot_parser = subparsers.add_parser("snapshot")
    snapshot_parser.add_argument("--starting-graph", required=True)
    snapshot_parser.add_argument("--evidence-input", action="append", required=True)
    snapshot_parser.add_argument("--item-index")
    snapshot_parser.add_argument("--audit-index")
    snapshot_parser.add_argument("--result-graph")
    snapshot_parser.add_argument("--output", required=True)
    snapshot_parser.add_argument("--repo-root", default=".")
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--current", required=True)
    compare_parser.add_argument("--prior", required=True)
    compare_parser.add_argument("--output", required=True)
    compare_parser.add_argument("--repo-root", default=".")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return snapshot(args) if args.command == "snapshot" else compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
