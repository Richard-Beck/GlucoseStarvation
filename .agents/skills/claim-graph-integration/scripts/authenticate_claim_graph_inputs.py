#!/usr/bin/env python3
"""Authenticate figure-level claim-audit inputs and reusable audit receipts."""

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple


CONTRACT_VERSION = "claim-graph-input-authentication/v2"


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
        raise SystemExit(f"{path}: invalid JSON: {exc}")


def display_path(repo_root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def hash_path(path: Path) -> Tuple[str, str, int]:
    if not path.exists() and not path.is_symlink():
        raise SystemExit(f"Declared input does not exist: {path}")
    if path.is_symlink():
        target = os.readlink(path)
        encoded = target.encode("utf-8")
        return "symlink", sha256_bytes(encoded), len(encoded)
    if path.is_file():
        return "file", sha256_file(path), path.stat().st_size
    if not path.is_dir():
        raise SystemExit(f"Unsupported declared input type: {path}")
    entries: List[Dict[str, Any]] = []
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
        raise SystemExit(f"--evidence-input must use label=path: {value}")
    label, path_text = value.split("=", 1)
    label = label.strip()
    path_text = path_text.strip()
    if not label or not path_text:
        raise SystemExit("--evidence-input requires nonempty label and path")
    return label, Path(path_text)


def figure_hashes(index_path: Path) -> Tuple[str, Dict[str, str], int]:
    payload = read_json(index_path)
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise SystemExit(f"{index_path}: figure index schema_version must be 1")
    figures = payload.get("figures")
    if not isinstance(figures, list):
        raise SystemExit(f"{index_path}: figures must be a list")
    hashes: Dict[str, str] = {}
    canonical_figures = []
    for figure in figures:
        if not isinstance(figure, dict):
            raise SystemExit(f"{index_path}: figure entries must be objects")
        figure_id = str(figure.get("figure_id", "")).strip()
        if not figure_id or "canonical_input" not in figure:
            raise SystemExit(
                f"{index_path}: figures require figure_id and canonical_input"
            )
        if figure_id in hashes:
            raise SystemExit(f"{index_path}: duplicate figure_id {figure_id}")
        digest = sha256_bytes(canonical_json(figure["canonical_input"]))
        hashes[figure_id] = digest
        canonical_figures.append({"figure_id": figure_id, "sha256": digest})
    canonical_figures.sort(key=lambda item: item["figure_id"])
    return (
        sha256_bytes(canonical_json(canonical_figures)),
        hashes,
        len(figures),
    )


def claims_hash(index_path: Path) -> Tuple[str, int]:
    payload = read_json(index_path)
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise SystemExit(f"{index_path}: claims index schema_version must be 1")
    claims = payload.get("claims")
    if not isinstance(claims, list):
        raise SystemExit(f"{index_path}: claims must be a list")
    seen = set()
    canonical_claims = []
    for claim in claims:
        if not isinstance(claim, dict):
            raise SystemExit(f"{index_path}: claim entries must be objects")
        claim_id = str(claim.get("id", "")).strip()
        text = str(claim.get("text", "")).strip()
        if not claim_id or not text:
            raise SystemExit(f"{index_path}: claims require id and text")
        if claim_id in seen:
            raise SystemExit(f"{index_path}: duplicate claim id {claim_id}")
        seen.add(claim_id)
        canonical_claims.append(claim)
    canonical_claims.sort(key=lambda item: str(item["id"]))
    return sha256_bytes(canonical_json(canonical_claims)), len(claims)


def audit_records(
    repo_root: Path,
    index_path: Path,
    figure_ids: set,
) -> Dict[str, Dict[str, str]]:
    payload = read_json(index_path)
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise SystemExit(f"{index_path}: audit index schema_version must be 1")
    audits = payload.get("audits")
    if not isinstance(audits, list):
        raise SystemExit(f"{index_path}: audits must be a list")
    records: Dict[str, Dict[str, str]] = {}
    for audit in audits:
        if not isinstance(audit, dict):
            raise SystemExit(f"{index_path}: audit entries must be objects")
        figure_id = str(audit.get("figure_id", "")).strip()
        path_text = str(audit.get("path", "")).strip()
        if not figure_id or not path_text:
            raise SystemExit(f"{index_path}: audits require figure_id and path")
        if figure_id in records:
            raise SystemExit(f"{index_path}: duplicate audit figure_id {figure_id}")
        if figure_id not in figure_ids:
            raise SystemExit(
                f"{index_path}: audit figure_id is absent from figure index: "
                f"{figure_id}"
            )
        path = Path(path_text)
        if not path.is_absolute():
            path = repo_root / path
        path = path.resolve()
        if not path.is_file():
            raise SystemExit(f"Audit is not a file: {path}")
        digest = sha256_file(path)
        declared_digest = str(audit.get("sha256", "")).strip()
        if declared_digest and declared_digest != digest:
            raise SystemExit(f"Audit checksum mismatch: {path}")
        records[figure_id] = {
            "path": display_path(repo_root, path),
            "sha256": digest,
        }
    return dict(sorted(records.items()))


def authenticated_prior_audits(
    repo_root: Path,
    prior: Dict[str, Any],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str]]:
    valid: Dict[str, Dict[str, str]] = {}
    invalid: Dict[str, str] = {}
    records = prior.get("audit_records") or {}
    if not isinstance(records, dict):
        return valid, {"*": "prior audit_records is not an object"}
    for figure_id, record in records.items():
        figure_id = str(figure_id)
        if not isinstance(record, dict):
            invalid[figure_id] = "prior audit record is not an object"
            continue
        path_text = str(record.get("path", "")).strip()
        expected = str(record.get("sha256", "")).strip()
        if not path_text or not expected:
            invalid[figure_id] = "prior audit record lacks path or sha256"
            continue
        path = Path(path_text)
        if not path.is_absolute():
            path = repo_root / path
        if not path.is_file():
            invalid[figure_id] = "prior audit file is missing"
            continue
        if sha256_file(path) != expected:
            invalid[figure_id] = "prior audit checksum mismatch"
            continue
        valid[figure_id] = {
            "path": display_path(repo_root, path),
            "sha256": expected,
        }
    return valid, invalid


def authenticated_prior_result(
    repo_root: Path,
    prior: Dict[str, Any],
) -> Tuple[bool, str]:
    record = prior.get("result_graph")
    if not isinstance(record, dict):
        return False, "prior result graph record is missing"
    path_text = str(record.get("path", "")).strip()
    expected = str(record.get("sha256", "")).strip()
    if not path_text or not expected:
        return False, "prior result graph record lacks path or sha256"
    path = Path(path_text)
    if not path.is_absolute():
        path = repo_root / path
    if not path.is_file():
        return False, "prior result graph file is missing"
    if sha256_file(path) != expected:
        return False, "prior result graph checksum mismatch"
    return True, ""


def snapshot(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    figure_index = Path(args.figure_index).resolve()
    claims_index = Path(args.claims_index).resolve()
    audit_contract = Path(args.audit_contract).resolve()
    for label, path in (
        ("figure index", figure_index),
        ("claims index", claims_index),
        ("audit contract", audit_contract),
    ):
        if not path.is_file():
            raise SystemExit(f"{label} is not a file: {path}")

    figure_index_hash, figures, figure_count = figure_hashes(figure_index)
    authoritative_claims_hash, claim_count = claims_hash(claims_index)
    audit_contract_hash = sha256_file(audit_contract)
    audit_context_hash = sha256_bytes(
        canonical_json(
            {
                "contract_version": CONTRACT_VERSION,
                "authoritative_claims_sha256": authoritative_claims_hash,
                "audit_contract_sha256": audit_contract_hash,
            }
        )
    )

    declared = []
    seen_labels = set()
    for raw in args.evidence_input:
        label, path = parse_evidence_input(raw)
        if label in seen_labels:
            raise SystemExit(f"Duplicate evidence-input label: {label}")
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

    audits: Dict[str, Dict[str, str]] = {}
    if args.audit_index:
        audits = audit_records(
            repo_root,
            Path(args.audit_index).resolve(),
            set(figures),
        )

    result_graph = None
    if args.result_graph:
        result_path = Path(args.result_graph).resolve()
        if not result_path.is_file():
            raise SystemExit(f"Result graph is not a file: {result_path}")
        result_graph = {
            "path": display_path(repo_root, result_path),
            "sha256": sha256_file(result_path),
        }

    whole_payload = {
        "contract_version": CONTRACT_VERSION,
        "figure_index_sha256": figure_index_hash,
        "authoritative_claims_sha256": authoritative_claims_hash,
        "audit_contract_sha256": audit_contract_hash,
        "evidence_package_sha256": evidence_package_hash,
    }
    output = {
        "schema_version": 2,
        "contract_version": CONTRACT_VERSION,
        "figure_index": {
            "path": display_path(repo_root, figure_index),
            "sha256": figure_index_hash,
            "figure_count": figure_count,
        },
        "figure_hashes": dict(sorted(figures.items())),
        "authoritative_claims": {
            "path": display_path(repo_root, claims_index),
            "sha256": authoritative_claims_hash,
            "claim_count": claim_count,
        },
        "audit_contract": {
            "path": display_path(repo_root, audit_contract),
            "sha256": audit_contract_hash,
        },
        "audit_context_sha256": audit_context_hash,
        "evidence_inputs": declared,
        "evidence_package_sha256": evidence_package_hash,
        "whole_input_sha256": sha256_bytes(canonical_json(whole_payload)),
        "audit_records": audits,
        "result_graph": result_graph,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


def compare(args: argparse.Namespace) -> int:
    repo_root = Path(args.repo_root).resolve()
    current = read_json(Path(args.current))
    prior = read_json(Path(args.prior))
    for label, payload in (("current", current), ("prior", prior)):
        if not isinstance(payload, dict) or payload.get("schema_version") != 2:
            raise SystemExit(f"{label} authentication schema_version must be 2")

    version_match = (
        current.get("contract_version")
        == prior.get("contract_version")
        == CONTRACT_VERSION
    )
    context_match = bool(
        version_match
        and current.get("audit_context_sha256")
        == prior.get("audit_context_sha256")
    )
    whole_match = bool(
        context_match
        and current.get("whole_input_sha256") == prior.get("whole_input_sha256")
    )
    current_figures = current.get("figure_hashes") or {}
    prior_figures = prior.get("figure_hashes") or {}
    if not isinstance(current_figures, dict) or not isinstance(prior_figures, dict):
        raise SystemExit("figure_hashes must be objects")

    valid_audits, invalid_audits = authenticated_prior_audits(repo_root, prior)
    identical = sorted(
        figure_id
        for figure_id, digest in current_figures.items()
        if context_match and prior_figures.get(figure_id) == digest
    )
    for figure_id in identical:
        if figure_id not in valid_audits:
            invalid_audits.setdefault(
                figure_id,
                "prior audit record is missing",
            )
    reusable = sorted(
        figure_id for figure_id in identical if figure_id in valid_audits
    )
    changed = sorted(
        figure_id
        for figure_id, digest in current_figures.items()
        if figure_id in prior_figures and prior_figures.get(figure_id) != digest
    )
    new = sorted(set(current_figures).difference(prior_figures))
    removed = sorted(set(prior_figures).difference(current_figures))
    misses = sorted(set(current_figures).difference(reusable))
    audits_complete = set(current_figures).issubset(reusable)
    result_valid, result_error = authenticated_prior_result(repo_root, prior)

    if whole_match and audits_complete and result_valid:
        mode = "reuse_complete"
    elif reusable:
        mode = "reuse_partial"
    else:
        mode = "fresh_full"
        reusable = []
        misses = sorted(current_figures)

    report = {
        "schema_version": 2,
        "contract_version": CONTRACT_VERSION,
        "mode": mode,
        "contract_version_match": version_match,
        "audit_context_match": context_match,
        "whole_input_match": whole_match,
        "audit_receipts_complete": audits_complete,
        "prior_result_graph_authenticated": result_valid,
        "prior_result_graph_error": result_error,
        "reusable_figure_ids": reusable,
        "reused_audits": {
            figure_id: valid_audits[figure_id] for figure_id in reusable
        },
        "cache_miss_figure_ids": misses,
        "invalid_prior_audits": {
            figure_id: reason
            for figure_id, reason in sorted(invalid_audits.items())
            if figure_id in current_figures or figure_id == "*"
        },
        "changed_figure_ids": changed,
        "new_figure_ids": new,
        "removed_figure_ids": removed,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    snapshot_parser = subparsers.add_parser("snapshot")
    snapshot_parser.add_argument("--figure-index", required=True)
    snapshot_parser.add_argument("--claims-index", required=True)
    snapshot_parser.add_argument("--audit-contract", required=True)
    snapshot_parser.add_argument("--evidence-input", action="append", default=[])
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
