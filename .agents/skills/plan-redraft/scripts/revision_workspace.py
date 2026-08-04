#!/usr/bin/env python3
"""Manage staged revisions of an immutable manuscript assembly.

The controller is deliberately mechanical.  It authenticates an existing
checksum inventory, manages ordered batch overlays, records observed changes
and exact-hash references, validates agent-authored consequence audits,
materializes review candidates, and invokes a configured preview renderer.  It
never edits the baseline and never adjudicates semantic validity.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from validate_consequence_audit import (
    AUDIT_COLUMNS as AUDIT_INPUT_COLUMNS,
    AuditError,
    CLASSIFICATION_COLUMNS,
    GRAPH_COLUMNS,
    validate_audit_files,
    validate_graph,
)

SCHEMA_VERSION = 2
CHANGE_COLUMNS = (
    "artifact_path",
    "kind",
    "change_type",
    "parent_sha256",
    "current_sha256",
    "recorded_at",
)
CONSEQUENCE_COLUMNS = (
    "consequence_id",
    "triggering_batch",
    *AUDIT_INPUT_COLUMNS,
)
SNAPSHOT_COLUMNS = ("path", "sha256", "size", "mtime_ns")
BATCH_INDEX_COLUMNS = (
    "batch_id",
    "state",
    "created_at",
    "predecessor_batches",
    "prompt_path",
    "assembly_contract_disposition",
    "accepted_overlaps",
)
CONTRACT_CHANGE_COLUMNS = ("artifact_path", "matched_pattern")
HASH_REFERENCE_COLUMNS = ("consumer_path", "parent_sha256", "changed_artifacts")


class ControllerError(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(encoded)


def safe_id(value: str, field: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value or ""):
        raise ControllerError(f"{field} must contain only letters, digits, '.', '_', and '-': {value!r}")
    return value


def clean_relative(value: str, field: str) -> Path:
    path = Path(value)
    if path.is_absolute() or not path.parts or path == Path(".") or ".." in path.parts:
        raise ControllerError(f"{field} must be a safe relative path: {value!r}")
    return path


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ControllerError(f"Cannot read JSON {path}: {exc}") from exc


def contract_approval_source(value: str | None, disposition: str) -> Path | None:
    if disposition == "preserved":
        return None
    if not value:
        raise ControllerError(
            f"Changing contract disposition to {disposition!r} requires --approval-file"
        )
    source = Path(value).resolve()
    if not source.is_file():
        raise ControllerError(f"Contract approval file does not exist: {source}")
    try:
        content = source.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise ControllerError(f"Contract approval must be a readable UTF-8 text file: {source}") from exc
    if not content.strip():
        raise ControllerError(f"Contract approval file is empty: {source}")
    return source


def store_contract_approval(root: Path, source: Path, disposition: str) -> str:
    decisions = root / "decisions"
    decisions.mkdir(parents=True, exist_ok=True)
    sequence = len(list(decisions.glob("*_contract_disposition_*.md"))) + 1
    destination = decisions / f"{sequence:04d}_contract_disposition_{disposition}.md"
    shutil.copy2(source, destination)
    return str(destination.relative_to(root))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], columns: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in writer.fieldnames})
    temporary.replace(path)


def read_tsv(path: Path, columns: Iterable[str], *, allow_empty: bool = True) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            missing = [column for column in columns if column not in fields]
            if missing:
                raise ControllerError(f"{path} is missing column(s): {', '.join(missing)}")
            rows = list(reader)
    except OSError as exc:
        raise ControllerError(f"Cannot read TSV {path}: {exc}") from exc
    if not allow_empty and not rows:
        raise ControllerError(f"{path} has no rows")
    return rows


def config_path_value(value: str, config_path: Path, field: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def load_initial_config(path: Path) -> dict[str, Any]:
    path = path.resolve()
    config = read_json(path)
    if config.get("schema_version") != SCHEMA_VERSION:
        raise ControllerError(f"Unsupported schema_version: {config.get('schema_version')!r}")
    if config.get("mode") != "in_place":
        raise ControllerError("The deterministic workspace manager supports mode='in_place' only")

    assembly_value = config.get("assembly_root")
    if not isinstance(assembly_value, str) or not assembly_value:
        raise ControllerError("assembly_root must be a non-empty path")
    assembly_root = config_path_value(assembly_value, path, "assembly_root")
    if not assembly_root.is_dir():
        raise ControllerError(f"assembly_root is not a directory: {assembly_root}")

    workspace_rel = clean_relative(config.get("workspace", ""), "workspace")
    if len(workspace_rel.parts) != 1:
        raise ControllerError("workspace must be one top-level directory within the assembly root")

    baseline = config.get("baseline")
    if not isinstance(baseline, dict):
        raise ControllerError("baseline must be an object")
    checksum_rel = clean_relative(baseline.get("checksum_manifest", ""), "baseline.checksum_manifest")
    status_value = baseline.get("status_file")
    status_rel = clean_relative(status_value, "baseline.status_file") if status_value else None
    allowed = baseline.get("allowed_unsealed_paths", [])
    if not isinstance(allowed, list) or not all(isinstance(value, str) and value for value in allowed):
        raise ControllerError("baseline.allowed_unsealed_paths must be a list of path patterns")

    renderer = config.get("preview_renderer")
    if not isinstance(renderer, dict):
        raise ControllerError("preview_renderer must be an object")
    argv = renderer.get("argv")
    if not isinstance(argv, list) or not argv or not all(isinstance(value, str) and value for value in argv):
        raise ControllerError("preview_renderer.argv must be a non-empty string array")
    output_rel = clean_relative(renderer.get("output", ""), "preview_renderer.output")
    cwd_value = renderer.get("cwd", "{candidate_root}")
    if not isinstance(cwd_value, str) or not cwd_value:
        raise ControllerError("preview_renderer.cwd must be a non-empty string")
    environment = renderer.get("environment", {})
    if not isinstance(environment, dict) or not all(
        isinstance(key, str) and isinstance(value, str) for key, value in environment.items()
    ):
        raise ControllerError("preview_renderer.environment must be a string-to-string object")

    audit_config = config.get("consequence_audit")
    if not isinstance(audit_config, dict):
        raise ControllerError("consequence_audit must be an object")
    graph_value = audit_config.get("dependency_graph")
    if not isinstance(graph_value, str) or not graph_value:
        raise ControllerError("consequence_audit.dependency_graph must name a TSV graph")
    graph_path = config_path_value(graph_value, path, "consequence_audit.dependency_graph")
    graph_rows = read_tsv(graph_path, GRAPH_COLUMNS)
    try:
        validate_graph(graph_rows)
    except AuditError as exc:
        raise ControllerError(str(exc)) from exc

    sensitive_patterns = config.get("contract_sensitive_patterns")
    if not isinstance(sensitive_patterns, list) or not sensitive_patterns or not all(
        isinstance(value, str) and value for value in sensitive_patterns
    ):
        raise ControllerError("contract_sensitive_patterns must be a non-empty list of path patterns")
    for pattern in sensitive_patterns:
        if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
            raise ControllerError(f"Invalid contract-sensitive pattern: {pattern!r}")
        if any(marker in pattern for marker in ("*", "?", "[")):
            raise ControllerError(
                "contract_sensitive_patterns must name exact control-plane paths, not globs: "
                f"{pattern!r}"
            )

    normalized = {
        "schema_version": SCHEMA_VERSION,
        "mode": "in_place",
        "assembly_root": str(assembly_root),
        "workspace": str(workspace_rel),
        "preview_policy": config.get("preview_policy", "batch_complete"),
        "baseline": {
            "checksum_manifest": str(checksum_rel),
            "status_file": str(status_rel) if status_rel else None,
            "require_complete_inventory": bool(baseline.get("require_complete_inventory", True)),
            "allowed_unsealed_paths": allowed,
        },
        "preview_renderer": {
            "argv": argv,
            "cwd": cwd_value,
            "output": str(output_rel),
            "environment": environment,
            "timeout_seconds": int(renderer.get("timeout_seconds", 300)),
        },
        "consequence_audit": {
            "dependency_graph": "controller/manuscript_dependency_graph.tsv",
        },
        "contract_sensitive_patterns": sensitive_patterns,
    }
    if normalized["preview_renderer"]["timeout_seconds"] <= 0:
        raise ControllerError("preview_renderer.timeout_seconds must be positive")
    normalized["_source_config"] = str(path)
    normalized["_source_dependency_graph"] = str(graph_path)
    return normalized


def parse_checksum_manifest(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ControllerError(f"Cannot read checksum manifest {path}: {exc}") from exc
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        parts = line.split("  ", 1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-fA-F]{64}", parts[0]):
            raise ControllerError(f"Malformed checksum row {path}:{line_number}")
        rel = str(clean_relative(parts[1], f"checksum row {line_number}"))
        if rel in records:
            raise ControllerError(f"Duplicate checksum path: {rel}")
        records[rel] = parts[0].lower()
    if not records:
        raise ControllerError(f"Checksum manifest has no records: {path}")
    return records


def matches_any(path: str, patterns: Iterable[str]) -> bool:
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)


def snapshot_baseline(assembly_root: Path, workspace: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(value for value in assembly_root.rglob("*") if value.is_file()):
        if path == workspace or workspace in path.parents:
            continue
        stat = path.stat()
        rows.append(
            {
                "path": str(path.relative_to(assembly_root)),
                "sha256": sha256_file(path),
                "size": str(stat.st_size),
                "mtime_ns": str(stat.st_mtime_ns),
            }
        )
    return rows


def describe_snapshot_drift(
    stored_rows: list[dict[str, str]], current_rows: list[dict[str, str]]
) -> str:
    stored = {row["path"]: row for row in stored_rows}
    current = {row["path"]: row for row in current_rows}
    messages: list[str] = []
    for path in sorted(set(stored) | set(current)):
        if path not in stored:
            messages.append(f"added {path}")
        elif path not in current:
            messages.append(f"deleted {path}")
        elif stored[path] != current[path]:
            fields = [field for field in SNAPSHOT_COLUMNS[1:] if stored[path][field] != current[path][field]]
            messages.append(f"changed {path} ({', '.join(fields)})")
    return "; ".join(messages[:12])


def authenticate(config: dict[str, Any], *, require_existing_workspace_receipt: bool = False) -> dict[str, Any]:
    assembly_root = Path(config["assembly_root"])
    workspace_rel = clean_relative(config["workspace"], "workspace")
    workspace = assembly_root / workspace_rel
    baseline = config["baseline"]
    manifest = assembly_root / clean_relative(baseline["checksum_manifest"], "checksum_manifest")
    if not manifest.is_file():
        raise ControllerError(f"Missing baseline checksum manifest: {manifest}")
    records = parse_checksum_manifest(manifest)

    mismatches: list[dict[str, str]] = []
    for rel, expected in records.items():
        path = assembly_root / rel
        if not path.is_file():
            mismatches.append({"path": rel, "expected": expected, "observed": "MISSING"})
            continue
        observed = sha256_file(path)
        if observed != expected:
            mismatches.append({"path": rel, "expected": expected, "observed": observed})
    if mismatches:
        preview = "; ".join(f"{row['path']} ({row['observed']})" for row in mismatches[:8])
        raise ControllerError(f"Baseline checksum authentication failed for {len(mismatches)} file(s): {preview}")

    manifest_rel = str(manifest.relative_to(assembly_root))
    allowed_patterns = list(baseline.get("allowed_unsealed_paths", []))
    observed_unsealed: list[str] = []
    unexpected: list[str] = []
    if baseline.get("require_complete_inventory", True):
        for path in sorted(value for value in assembly_root.rglob("*") if value.is_file()):
            if path == manifest or path == workspace or workspace in path.parents:
                continue
            rel = str(path.relative_to(assembly_root))
            if rel in records:
                continue
            observed_unsealed.append(rel)
            if not matches_any(rel, allowed_patterns):
                unexpected.append(rel)
        if unexpected:
            raise ControllerError(
                "Baseline contains file(s) outside the checksum inventory and allowed exclusions: "
                + ", ".join(unexpected[:12])
            )

    status_rel = baseline.get("status_file")
    status_sha = None
    if status_rel:
        status_path = assembly_root / clean_relative(status_rel, "status_file")
        if not status_path.is_file():
            raise ControllerError(f"Missing baseline status file: {status_path}")
        status_sha = sha256_file(status_path)

    current_snapshot = snapshot_baseline(assembly_root, workspace)
    baseline_tree_id = canonical_hash(current_snapshot)

    identity_payload = {
        "checksum_manifest_sha256": sha256_file(manifest),
        "status_sha256": status_sha,
        "verified_file_count": len(records),
        "baseline_tree_id": baseline_tree_id,
    }
    receipt = {
        "schema_version": SCHEMA_VERSION,
        "authenticated_at": now(),
        "assembly_root": str(assembly_root),
        "workspace": str(workspace_rel),
        "checksum_manifest": manifest_rel,
        "checksum_manifest_sha256": identity_payload["checksum_manifest_sha256"],
        "status_file": status_rel,
        "status_sha256": status_sha,
        "verified_file_count": len(records),
        "require_complete_inventory": bool(baseline.get("require_complete_inventory", True)),
        "allowed_unsealed_paths": allowed_patterns,
        "observed_unsealed_paths": observed_unsealed,
        "baseline_file_count": len(current_snapshot),
        "baseline_tree_id": baseline_tree_id,
        "baseline_id": canonical_hash(identity_payload),
        "result": "PASS",
    }

    if require_existing_workspace_receipt:
        stored_path = workspace / "baseline_receipt.json"
        stored = read_json(stored_path)
        snapshot_path = workspace / "controller/baseline_snapshot.tsv"
        stored_snapshot = read_tsv(snapshot_path, SNAPSHOT_COLUMNS)
        if stored_snapshot != current_snapshot:
            detail = describe_snapshot_drift(stored_snapshot, current_snapshot)
            raise ControllerError(f"Immutable baseline snapshot changed: {detail}")
        if stored.get("baseline_id") != receipt["baseline_id"]:
            raise ControllerError("Current baseline identity differs from the workspace receipt")
    else:
        receipt["_baseline_snapshot"] = current_snapshot
    return receipt


def workspace_from_arg(value: str) -> Path:
    path = Path(value).resolve()
    if not path.is_dir():
        raise ControllerError(f"Workspace is not a directory: {path}")
    return path


def load_workspace(value: str, *, verify_baseline: bool = True) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    workspace = workspace_from_arg(value)
    config = read_json(workspace / "controller/controller_config.json")
    if config.get("schema_version") != SCHEMA_VERSION or config.get("mode") != "in_place":
        raise ControllerError("Unsupported workspace configuration")
    state = read_json(workspace / "working_state.json")
    if verify_baseline:
        authenticate(config, require_existing_workspace_receipt=True)
    return workspace, config, state


def batch_dir(workspace: Path, batch_id: str) -> Path:
    return workspace / "batches" / safe_id(batch_id, "batch_id")


def read_batch_index(workspace: Path) -> list[dict[str, str]]:
    return read_tsv(workspace / "batch_index.tsv", BATCH_INDEX_COLUMNS)


def write_batch_index(workspace: Path, rows: list[dict[str, str]]) -> None:
    write_tsv(workspace / "batch_index.tsv", rows, BATCH_INDEX_COLUMNS)


def update_state(workspace: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = now()
    write_json(workspace / "working_state.json", state)


def mark_batch_active(workspace: Path, batch_id: str, state: dict[str, Any]) -> None:
    metadata_path = batch_dir(workspace, batch_id) / "batch.json"
    metadata = read_json(metadata_path)
    if metadata.get("state") == "abandoned":
        raise ControllerError("An abandoned batch is not editable")
    metadata["state"] = "active"
    metadata["updated_at"] = now()
    write_json(metadata_path, metadata)
    rows = read_batch_index(workspace)
    for row in rows:
        if row["batch_id"] == batch_id:
            row["state"] = "active"
    write_batch_index(workspace, rows)
    state["latest_preview_current"] = False
    update_state(workspace, state)


def init_workspace(args: argparse.Namespace) -> int:
    config = load_initial_config(args.config)
    assembly_root = Path(config["assembly_root"])
    workspace = assembly_root / config["workspace"]
    receipt = authenticate(config)
    baseline_snapshot = receipt.pop("_baseline_snapshot")
    if workspace.exists():
        if not (workspace / "controller/controller_config.json").is_file():
            raise ControllerError(f"Workspace exists but is not initialized: {workspace}")
        _, existing, _ = load_workspace(str(workspace))
        comparable = {key: value for key, value in config.items() if not key.startswith("_source_")}
        if existing != comparable:
            raise ControllerError(f"Workspace already exists with a different configuration: {workspace}")
        print(json.dumps({"status": "EXISTS", "workspace": str(workspace), "baseline_id": receipt["baseline_id"]}, indent=2))
        return 0

    for rel in ("controller", "renderer", "batches", "candidates", "previews"):
        (workspace / rel).mkdir(parents=True, exist_ok=False if rel == "controller" else True)

    persistent = {key: value for key, value in config.items() if not key.startswith("_source_")}
    write_json(workspace / "controller/controller_config.json", persistent)
    shutil.copy2(
        Path(config["_source_dependency_graph"]),
        workspace / "controller/manuscript_dependency_graph.tsv",
    )
    shutil.copy2(Path(__file__).resolve(), workspace / "controller/revision_workspace.py")
    shutil.copy2(
        Path(__file__).resolve().with_name("validate_consequence_audit.py"),
        workspace / "controller/validate_consequence_audit.py",
    )
    write_tsv(workspace / "controller/baseline_snapshot.tsv", baseline_snapshot, SNAPSHOT_COLUMNS)
    write_json(workspace / "renderer/renderer_contract.json", persistent["preview_renderer"])
    write_json(workspace / "baseline_receipt.json", receipt)
    write_tsv(workspace / "batch_index.tsv", [], BATCH_INDEX_COLUMNS)
    write_tsv(workspace / "cumulative_consequences.tsv", [], CONSEQUENCE_COLUMNS)
    state = {
        "schema_version": SCHEMA_VERSION,
        "mode": "in_place",
        "baseline_id": receipt["baseline_id"],
        "preview_policy": persistent["preview_policy"],
        "active_batches": [],
        "current_candidate": None,
        "latest_preview": None,
        "latest_preview_current": False,
        "created_at": now(),
    }
    update_state(workspace, state)
    print(json.dumps({"status": "INITIALIZED", "workspace": str(workspace), "baseline_id": receipt["baseline_id"]}, indent=2))
    return 0


def command_authenticate(args: argparse.Namespace) -> int:
    config = load_initial_config(args.config)
    receipt = authenticate(config)
    receipt.pop("_baseline_snapshot", None)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


def begin_batch(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    rows = read_batch_index(workspace)
    batch_id = safe_id(args.batch_id or datetime.now().strftime("%Y%m%dT%H%M%S_revision"), "batch_id")
    if any(row["batch_id"] == batch_id for row in rows) or batch_dir(workspace, batch_id).exists():
        raise ControllerError(f"Batch already exists: {batch_id}")
    prompt_path = Path(args.prompt_file).resolve()
    if not prompt_path.is_file():
        raise ControllerError(f"Prompt file does not exist: {prompt_path}")
    content = prompt_path.read_text(encoding="utf-8")
    if not content.strip():
        raise ControllerError("Prompt file is empty")
    accepted_overlaps = [safe_id(value, "allow_overlap_with") for value in args.allow_overlap_with]
    known_ids = {row["batch_id"] for row in rows}
    unknown = sorted(set(accepted_overlaps) - known_ids)
    if unknown:
        raise ControllerError(f"Accepted-overlap batch does not exist: {', '.join(unknown)}")
    approval_source = contract_approval_source(args.contract_approval_file, args.contract_disposition)

    root = batch_dir(workspace, batch_id)
    for rel in ("prompts", "decisions", "checks", "files", "generated"):
        (root / rel).mkdir(parents=True, exist_ok=True)
    (root / "prompts/0001.md").write_text(content, encoding="utf-8")
    approval_path = (
        store_contract_approval(root, approval_source, args.contract_disposition)
        if approval_source
        else None
    )
    write_tsv(root / "scope.tsv", [], ("artifact", "intent"))
    write_tsv(root / "deletions.tsv", [], ("artifact_path",))
    write_tsv(root / "artifact_changes.tsv", [], CHANGE_COLUMNS)
    write_tsv(root / "artifact_classification.tsv", [], CLASSIFICATION_COLUMNS)
    write_tsv(root / "consequence_delta.tsv", [], CONSEQUENCE_COLUMNS)
    write_tsv(root / "consequence_audit.tsv", [], AUDIT_INPUT_COLUMNS)
    write_json(root / "consequence_audit_validation.json", {"status": "NOT_RECORDED"})
    write_tsv(root / "contract_sensitive_changes.tsv", [], CONTRACT_CHANGE_COLUMNS)
    write_tsv(root / "hash_backreferences.tsv", [], HASH_REFERENCE_COLUMNS)
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "batch_id": batch_id,
        "state": "active",
        "created_at": now(),
        "predecessor_batches": list(state["active_batches"]),
        "accepted_overlaps": accepted_overlaps,
        "assembly_contract_disposition": args.contract_disposition,
        "contract_disposition_approval": approval_path,
        "consequence_audit_status": "not_recorded",
        "prompt_count": 1,
    }
    write_json(root / "batch.json", metadata)
    rows.append(
        {
            "batch_id": batch_id,
            "state": "active",
            "created_at": metadata["created_at"],
            "predecessor_batches": ";".join(metadata["predecessor_batches"]),
            "prompt_path": f"batches/{batch_id}/prompts/0001.md",
            "assembly_contract_disposition": args.contract_disposition,
            "accepted_overlaps": ";".join(accepted_overlaps),
        }
    )
    write_batch_index(workspace, rows)
    state["active_batches"].append(batch_id)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(json.dumps({"status": "BATCH_STARTED", "batch_id": batch_id, "overlay": str(root)}, indent=2))
    return 0


def append_prompt(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if metadata.get("state") == "abandoned":
        raise ControllerError("Cannot append a prompt to an abandoned batch")
    source = Path(args.prompt_file).resolve()
    if not source.is_file() or not source.read_text(encoding="utf-8").strip():
        raise ControllerError(f"Prompt file is missing or empty: {source}")
    count = int(metadata.get("prompt_count", 1)) + 1
    destination = root / "prompts" / f"{count:04d}_clarification.md"
    shutil.copy2(source, destination)
    metadata["prompt_count"] = count
    metadata["updated_at"] = now()
    write_json(root / "batch.json", metadata)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(destination)
    return 0


def active_before(state: dict[str, Any], batch_id: str) -> list[str]:
    try:
        index = state["active_batches"].index(batch_id)
    except ValueError as exc:
        raise ControllerError(f"Batch is not active: {batch_id}") from exc
    return list(state["active_batches"][:index])


def overlay_candidates(workspace: Path, batch_id: str, rel: Path) -> list[Path]:
    root = batch_dir(workspace, batch_id)
    return [root / "files" / rel, root / "generated" / rel]


def deleted_paths(workspace: Path, batch_id: str) -> list[Path]:
    rows = read_tsv(batch_dir(workspace, batch_id) / "deletions.tsv", ("artifact_path",))
    return [clean_relative(row["artifact_path"], "deleted artifact") for row in rows]


def is_deleted(workspace: Path, batch_id: str, rel: Path) -> bool:
    for deletion in deleted_paths(workspace, batch_id):
        if rel == deletion or deletion in rel.parents:
            return True
    return False


def effective_file(workspace: Path, config: dict[str, Any], batches: list[str], rel: Path) -> Path | None:
    current = Path(config["assembly_root"]) / rel
    result: Path | None = current if current.is_file() else None
    for batch_id in batches:
        if is_deleted(workspace, batch_id, rel):
            result = None
        for candidate in overlay_candidates(workspace, batch_id, rel):
            if candidate.is_file():
                result = candidate
    return result


def merge_tree(source: Path, destination: Path) -> None:
    if source.is_file():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        return
    if not source.is_dir():
        raise ControllerError(f"Cannot merge missing path: {source}")
    destination.mkdir(parents=True, exist_ok=True)
    for path in sorted(source.rglob("*")):
        rel = path.relative_to(source)
        target = destination / rel
        if path.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        elif path.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)


def stage_path(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    rel = clean_relative(args.path, "path")
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if metadata.get("state") not in {"active", "recorded", "awaiting_review"}:
        raise ControllerError(f"Batch is not editable: {metadata.get('state')}")
    destination = root / "files" / rel
    if destination.exists():
        raise ControllerError(f"Path is already staged: {destination}")
    predecessors = active_before(state, args.batch)
    baseline_source = Path(config["assembly_root"]) / rel
    if baseline_source.is_file():
        source = effective_file(workspace, config, predecessors, rel)
        if source is None:
            raise ControllerError(f"Effective file does not exist: {rel}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    elif baseline_source.is_dir():
        shutil.copytree(baseline_source, destination)
        for predecessor in predecessors:
            for overlay_root in (batch_dir(workspace, predecessor) / "files", batch_dir(workspace, predecessor) / "generated"):
                source = overlay_root / rel
                if source.exists():
                    merge_tree(source, destination)
            for deletion in deleted_paths(workspace, predecessor):
                if rel == deletion:
                    if destination.exists():
                        shutil.rmtree(destination)
                elif rel in deletion.parents:
                    target = destination / deletion.relative_to(rel)
                    if target.is_dir():
                        shutil.rmtree(target)
                    elif target.exists():
                        target.unlink()
    else:
        source = effective_file(workspace, config, predecessors, rel)
        if source is None:
            raise ControllerError(f"Path does not exist in the effective parent: {rel}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    mark_batch_active(workspace, args.batch, state)
    print(destination)
    return 0


def register_generated(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    source = Path(args.source).resolve()
    if not source.is_file():
        raise ControllerError(f"Generated source is not a file: {source}")
    rel = clean_relative(args.path, "path")
    destination = batch_dir(workspace, args.batch) / "generated" / rel
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not args.replace:
        raise ControllerError(f"Generated artifact is already registered: {destination}")
    shutil.copy2(source, destination)
    mark_batch_active(workspace, args.batch, state)
    print(destination)
    return 0


def register_check(args: argparse.Namespace) -> int:
    workspace, _, _ = load_workspace(args.workspace)
    source = Path(args.source).resolve()
    if not source.is_file():
        raise ControllerError(f"Check or execution receipt is not a file: {source}")
    rel = clean_relative(args.name, "name")
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if metadata.get("state") == "abandoned":
        raise ControllerError("Cannot register evidence for an abandoned batch")
    destination = root / "checks" / rel
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not args.replace:
        raise ControllerError(f"Check or execution receipt is already registered: {destination}")
    shutil.copy2(source, destination)
    print(destination)
    return 0


def register_consequence_audit(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if metadata.get("state") not in {"recorded", "awaiting_review", "complete"}:
        raise ControllerError("Record the batch before registering its consequence audit")
    if metadata.get("recorded_overlay_fingerprint") != batch_overlay_fingerprint(workspace, args.batch):
        raise ControllerError("Batch overlay changed after recording; re-run record before auditing")
    classification_source = Path(args.classification_source).resolve()
    audit_source = Path(args.audit_source).resolve()
    if not classification_source.is_file():
        raise ControllerError(f"Artifact classification is not a file: {classification_source}")
    if not audit_source.is_file():
        raise ControllerError(f"Consequence audit is not a file: {audit_source}")
    classification_destination = root / "artifact_classification.tsv"
    audit_destination = root / "consequence_audit.tsv"
    if classification_source != classification_destination.resolve():
        shutil.copy2(classification_source, classification_destination)
    if audit_source != audit_destination.resolve():
        shutil.copy2(audit_source, audit_destination)

    audit_config = config["consequence_audit"]
    try:
        receipt = validate_audit_files(
            workspace / audit_config["dependency_graph"],
            classification_destination,
            root / "artifact_changes.tsv",
            root / "hash_backreferences.tsv",
            audit_destination,
        )
    except AuditError as exc:
        write_json(
            root / "consequence_audit_validation.json",
            {
                "status": "BLOCK",
                "validated_at": now(),
                "error": str(exc),
                "recorded_overlay_fingerprint": metadata.get("recorded_overlay_fingerprint"),
            },
        )
        metadata["consequence_audit_status"] = "block"
        metadata["updated_at"] = now()
        write_json(root / "batch.json", metadata)
        state["latest_preview_current"] = False
        update_state(workspace, state)
        raise ControllerError(str(exc)) from exc

    audit_rows = read_tsv(audit_destination, AUDIT_INPUT_COLUMNS)
    delta: list[dict[str, str]] = []
    for row in audit_rows:
        identity = {
            "triggering_batch": args.batch,
            "upstream_class": row["upstream_class"],
            "downstream_class": row["downstream_class"],
            "relationship": row["relationship"],
        }
        delta.append(
            {
                "consequence_id": "con_" + canonical_hash(identity)[:16],
                "triggering_batch": args.batch,
                **row,
            }
        )
    delta.sort(key=lambda row: (row["upstream_class"], row["downstream_class"], row["relationship"]))
    cumulative_path = workspace / "cumulative_consequences.tsv"
    cumulative = [
        row
        for row in read_tsv(cumulative_path, CONSEQUENCE_COLUMNS)
        if row["triggering_batch"] != args.batch
    ]
    cumulative.extend(delta)
    cumulative.sort(key=lambda row: (row["triggering_batch"], row["consequence_id"]))
    write_tsv(root / "consequence_delta.tsv", delta, CONSEQUENCE_COLUMNS)
    write_tsv(cumulative_path, cumulative, CONSEQUENCE_COLUMNS)

    receipt.update(
        {
            "validated_at": now(),
            "batch_id": args.batch,
            "recorded_overlay_fingerprint": metadata.get("recorded_overlay_fingerprint"),
        }
    )
    write_json(root / "consequence_audit_validation.json", receipt)
    metadata["consequence_audit_status"] = "pass"
    metadata["consequence_count"] = len(delta)
    metadata["updated_at"] = now()
    write_json(root / "batch.json", metadata)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


def register_deletion(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    rel = clean_relative(args.path, "path")
    predecessors = active_before(state, args.batch)
    baseline = Path(config["assembly_root"]) / rel
    if not baseline.exists() and effective_file(workspace, config, predecessors, rel) is None:
        raise ControllerError(f"Cannot delete an absent artifact: {rel}")
    path = batch_dir(workspace, args.batch) / "deletions.tsv"
    rows = read_tsv(path, ("artifact_path",))
    if any(row["artifact_path"] == str(rel) for row in rows):
        raise ControllerError(f"Deletion is already registered: {rel}")
    rows.append({"artifact_path": str(rel)})
    write_tsv(path, rows, ("artifact_path",))
    mark_batch_active(workspace, args.batch, state)
    print(rel)
    return 0


def overlay_files(root: Path) -> dict[str, tuple[str, Path]]:
    found: dict[str, tuple[str, Path]] = {}
    for kind in ("files", "generated"):
        base = root / kind
        for path in sorted(value for value in base.rglob("*") if value.is_file()):
            rel = str(path.relative_to(base))
            if rel in found:
                raise ControllerError(f"Artifact exists in both files/ and generated/: {rel}")
            found[rel] = (kind, path)
    return found


def file_hash_references(path: Path, hashes: Iterable[str]) -> set[str]:
    needles = {value.encode("ascii"): value for value in hashes}
    if not needles:
        return set()
    found: set[str] = set()
    overlap = max(len(value) for value in needles) - 1
    tail = b""
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                content = tail + block
                for needle, value in needles.items():
                    if value not in found and needle in content:
                        found.add(value)
                if len(found) == len(needles):
                    break
                tail = content[-overlap:] if overlap else b""
    except OSError as exc:
        raise ControllerError(f"Cannot inspect exact-hash consumer {path}: {exc}") from exc
    return found


def find_hash_backreferences(
    workspace: Path,
    config: dict[str, Any],
    predecessors: list[str],
    changes: list[dict[str, str]],
) -> list[dict[str, str]]:
    artifacts_by_hash: dict[str, set[str]] = {}
    changed_paths = {row["artifact_path"] for row in changes}
    for row in changes:
        parent_sha = row["parent_sha256"]
        if parent_sha:
            artifacts_by_hash.setdefault(parent_sha, set()).add(row["artifact_path"])
    if not artifacts_by_hash:
        return []

    assembly_root = Path(config["assembly_root"])
    candidate_paths = {
        str(path.relative_to(assembly_root))
        for path in assembly_root.rglob("*")
        if path.is_file() and workspace not in path.parents
    }
    for batch_id in predecessors:
        candidate_paths.update(overlay_files(batch_dir(workspace, batch_id)))

    rows: list[dict[str, str]] = []
    for rel_text in sorted(candidate_paths - changed_paths):
        rel = clean_relative(rel_text, "effective parent artifact")
        source = effective_file(workspace, config, predecessors, rel)
        if source is None:
            continue
        for parent_sha in sorted(file_hash_references(source, artifacts_by_hash)):
            rows.append(
                {
                    "consumer_path": rel_text,
                    "parent_sha256": parent_sha,
                    "changed_artifacts": ";".join(sorted(artifacts_by_hash[parent_sha])),
                }
            )
    return rows


def record_batch(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    root = batch_dir(workspace, args.batch)
    predecessors = active_before(state, args.batch)
    changes: list[dict[str, str]] = []
    timestamp = now()
    for rel_text, (kind, path) in overlay_files(root).items():
        rel = clean_relative(rel_text, "overlay path")
        parent = effective_file(workspace, config, predecessors, rel)
        parent_sha = sha256_file(parent) if parent else ""
        current_sha = sha256_file(path)
        change_type = "added" if parent is None else ("unchanged" if parent_sha == current_sha else "modified")
        if change_type != "unchanged":
            changes.append(
                {
                    "artifact_path": rel_text,
                    "kind": kind,
                    "change_type": change_type,
                    "parent_sha256": parent_sha,
                    "current_sha256": current_sha,
                    "recorded_at": timestamp,
                }
            )
    for rel in deleted_paths(workspace, args.batch):
        parent = effective_file(workspace, config, predecessors, rel)
        if parent:
            changes.append(
                {
                    "artifact_path": str(rel),
                    "kind": "deletion",
                    "change_type": "deleted",
                    "parent_sha256": sha256_file(parent),
                    "current_sha256": "",
                    "recorded_at": timestamp,
                }
            )
    changes.sort(key=lambda row: (row["artifact_path"], row["kind"]))
    write_tsv(root / "artifact_changes.tsv", changes, CHANGE_COLUMNS)
    changed = sorted({row["artifact_path"] for row in changes})
    sensitive_rows = [
        {"artifact_path": artifact, "matched_pattern": pattern}
        for artifact in changed
        for pattern in config["contract_sensitive_patterns"]
        if fnmatch.fnmatchcase(artifact, pattern)
    ]
    write_tsv(root / "contract_sensitive_changes.tsv", sensitive_rows, CONTRACT_CHANGE_COLUMNS)
    metadata = read_json(root / "batch.json")
    metadata["contract_sensitive_changed_artifact_count"] = len(
        {row["artifact_path"] for row in sensitive_rows}
    )
    write_json(root / "batch.json", metadata)
    if sensitive_rows and metadata.get("assembly_contract_disposition") == "preserved":
        paths = sorted({row["artifact_path"] for row in sensitive_rows})
        raise ControllerError(
            "Contract-sensitive scope expansion requires an explicit batch disposition: "
            + ", ".join(paths[:12])
        )
    hash_references = find_hash_backreferences(workspace, config, predecessors, changes)
    write_tsv(root / "hash_backreferences.tsv", hash_references, HASH_REFERENCE_COLUMNS)

    cumulative_path = workspace / "cumulative_consequences.tsv"
    cumulative = read_tsv(cumulative_path, CONSEQUENCE_COLUMNS)
    cumulative = [row for row in cumulative if row["triggering_batch"] != args.batch]
    cumulative.sort(key=lambda row: (row["triggering_batch"], row["consequence_id"]))
    write_tsv(root / "consequence_delta.tsv", [], CONSEQUENCE_COLUMNS)
    write_tsv(cumulative_path, cumulative, CONSEQUENCE_COLUMNS)
    write_tsv(root / "artifact_classification.tsv", [], CLASSIFICATION_COLUMNS)
    write_tsv(root / "consequence_audit.tsv", [], AUDIT_INPUT_COLUMNS)
    audit_required = bool(changes)
    write_json(
        root / "consequence_audit_validation.json",
        {
            "status": "PENDING" if audit_required else "NOT_REQUIRED",
            "recorded_overlay_fingerprint": batch_overlay_fingerprint(workspace, args.batch),
        },
    )

    metadata["state"] = "recorded"
    metadata["recorded_at"] = timestamp
    metadata["changed_artifact_count"] = len(changed)
    metadata["consequence_count"] = 0
    metadata["consequence_audit_status"] = "pending" if audit_required else "not_required"
    metadata["hash_backreference_consumer_count"] = len(
        {row["consumer_path"] for row in hash_references}
    )
    metadata["recorded_overlay_fingerprint"] = batch_overlay_fingerprint(workspace, args.batch)
    write_json(root / "batch.json", metadata)
    index = read_batch_index(workspace)
    for row in index:
        if row["batch_id"] == args.batch:
            row["state"] = "recorded"
    write_batch_index(workspace, index)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(
        json.dumps(
            {
                "status": "RECORDED",
                "batch_id": args.batch,
                "changed_artifacts": len(changed),
                "consequence_audit": "PENDING" if audit_required else "NOT_REQUIRED",
                "hash_backreference_consumers": metadata["hash_backreference_consumer_count"],
            },
            indent=2,
        )
    )
    return 0


def set_contract_disposition(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if metadata.get("state") == "abandoned":
        raise ControllerError("Cannot change the disposition of an abandoned batch")
    approval_source = contract_approval_source(args.approval_file, args.disposition)
    approval_path = (
        store_contract_approval(root, approval_source, args.disposition)
        if approval_source
        else None
    )
    metadata["assembly_contract_disposition"] = args.disposition
    metadata["contract_disposition_approval"] = approval_path
    metadata["updated_at"] = now()
    write_json(root / "batch.json", metadata)
    rows = read_batch_index(workspace)
    for row in rows:
        if row["batch_id"] == args.batch:
            row["assembly_contract_disposition"] = args.disposition
    write_batch_index(workspace, rows)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(json.dumps({"batch_id": args.batch, "assembly_contract_disposition": args.disposition}, indent=2))
    return 0


def batch_overlay_fingerprint(workspace: Path, batch_id: str) -> str:
    root = batch_dir(workspace, batch_id)
    files = [
        {"path": rel, "kind": kind, "sha256": sha256_file(path)}
        for rel, (kind, path) in sorted(overlay_files(root).items())
    ]
    prompts = [
        {"path": str(path.relative_to(root / "prompts")), "sha256": sha256_file(path)}
        for path in sorted(value for value in (root / "prompts").rglob("*") if value.is_file())
    ]
    decisions = [
        {"path": str(path.relative_to(root / "decisions")), "sha256": sha256_file(path)}
        for path in sorted(value for value in (root / "decisions").rglob("*") if value.is_file())
    ]
    deletions = [str(value) for value in deleted_paths(workspace, batch_id)]
    metadata = read_json(root / "batch.json")
    return canonical_hash(
        {
            "batch_id": batch_id,
            "files": files,
            "prompts": prompts,
            "decisions": decisions,
            "deletions": deletions,
            "assembly_contract_disposition": metadata.get("assembly_contract_disposition"),
        }
    )


def working_fingerprint(workspace: Path, state: dict[str, Any], batches: list[str] | None = None) -> str:
    selected = batches if batches is not None else list(state["active_batches"])
    return canonical_hash(
        [{"batch_id": batch_id, "overlay_fingerprint": batch_overlay_fingerprint(workspace, batch_id)} for batch_id in selected]
    )


def consequence_fingerprint(rows: list[dict[str, str]], active_batches: Iterable[str]) -> str:
    active = set(active_batches)
    selected = [row for row in rows if row["triggering_batch"] in active]
    selected.sort(key=lambda row: row["consequence_id"])
    return canonical_hash(selected)


def dependency_state_fingerprint(
    workspace: Path, rows: list[dict[str, str]], active_batches: Iterable[str]
) -> str:
    active = list(active_batches)
    classifications = {
        batch_id: read_tsv(
            batch_dir(workspace, batch_id) / "artifact_classification.tsv",
            CLASSIFICATION_COLUMNS,
        )
        for batch_id in active
    }
    audit_rows = {
        batch_id: read_tsv(
            batch_dir(workspace, batch_id) / "consequence_audit.tsv",
            AUDIT_INPUT_COLUMNS,
        )
        for batch_id in active
    }
    audits = {
        batch_id: read_json(batch_dir(workspace, batch_id) / "consequence_audit_validation.json")
        for batch_id in active
    }
    return canonical_hash(
        {
            "consequences": consequence_fingerprint(rows, active),
            "artifact_classifications": classifications,
            "consequence_audits": audit_rows,
            "audit_validation": audits,
        }
    )


def detect_overlaps(workspace: Path, batches: list[str]) -> None:
    seen: dict[str, str] = {}
    for batch_id in batches:
        metadata = read_json(batch_dir(workspace, batch_id) / "batch.json")
        accepted = set(metadata.get("accepted_overlaps", []))
        changes = read_tsv(batch_dir(workspace, batch_id) / "artifact_changes.tsv", CHANGE_COLUMNS)
        for row in changes:
            if row["change_type"] == "unchanged":
                continue
            prior = seen.get(row["artifact_path"])
            if prior and prior not in accepted:
                raise ControllerError(
                    f"Unresolved overlapping write to {row['artifact_path']}: {prior} and {batch_id}"
                )
            seen[row["artifact_path"]] = batch_id


def copy_baseline(config: dict[str, Any], destination: Path) -> None:
    source = Path(config["assembly_root"])
    excluded = config["workspace"]

    def ignore(directory: str, names: list[str]) -> set[str]:
        if Path(directory).resolve() == source.resolve() and excluded in names:
            return {excluded}
        return set()

    shutil.copytree(source, destination, ignore=ignore)


def apply_batch(workspace: Path, batch_id: str, assembly: Path) -> None:
    root = batch_dir(workspace, batch_id)
    for kind in ("files", "generated"):
        merge_tree(root / kind, assembly)
    for rel in deleted_paths(workspace, batch_id):
        target = assembly / rel
        if target.is_dir():
            shutil.rmtree(target)
        elif target.exists() or target.is_symlink():
            target.unlink()


def materialize(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    batches = list(args.batches or state["active_batches"])
    known = {row["batch_id"]: row for row in read_batch_index(workspace)}
    unknown = [batch for batch in batches if batch not in known]
    if unknown:
        raise ControllerError(f"Unknown batch(es): {', '.join(unknown)}")
    not_recorded = [
        batch
        for batch in batches
        if known[batch]["state"] not in {"recorded", "awaiting_review", "complete"}
    ]
    if not_recorded:
        raise ControllerError(f"Record batch changes before materialization: {', '.join(not_recorded)}")
    stale_records = [
        batch
        for batch in batches
        if read_json(batch_dir(workspace, batch) / "batch.json").get("recorded_overlay_fingerprint")
        != batch_overlay_fingerprint(workspace, batch)
    ]
    if stale_records:
        raise ControllerError(f"Batch overlay changed after recording; re-run record: {', '.join(stale_records)}")
    detect_overlaps(workspace, batches)
    prior_candidates = []
    for composition_path in sorted((workspace / "candidates").glob("*/composition.json")):
        composition = read_json(composition_path)
        if composition.get("active_batches") == batches:
            prior_candidates.append(composition.get("candidate_id", composition_path.parent.name))
    if prior_candidates and not args.additional_reason:
        raise ControllerError(
            "A candidate already exists for this batch set; provide --additional-reason after a user change "
            "or genuine render failure: " + ", ".join(prior_candidates)
        )
    candidate_id = safe_id(args.candidate_id or datetime.now().strftime("%Y%m%dT%H%M%S_candidate"), "candidate_id")
    root = workspace / "candidates" / candidate_id
    if root.exists():
        raise ControllerError(f"Candidate already exists: {candidate_id}")
    assembly = root / "assembly"
    root.mkdir(parents=True)
    copy_baseline(config, assembly)
    for batch in batches:
        apply_batch(workspace, batch, assembly)
    fingerprint = working_fingerprint(workspace, state, batches)
    composition = {
        "schema_version": SCHEMA_VERSION,
        "candidate_id": candidate_id,
        "created_at": now(),
        "baseline_id": state["baseline_id"],
        "active_batches": batches,
        "working_fingerprint": fingerprint,
        "additional_candidate_reason": args.additional_reason,
        "assembly_root": str(assembly),
        "sealed": False,
    }
    write_json(root / "composition.json", composition)
    state["current_candidate"] = candidate_id
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(json.dumps({"status": "MATERIALIZED", "candidate_id": candidate_id, "assembly": str(assembly)}, indent=2))
    return 0


def render_value(value: str, replacements: dict[str, str]) -> str:
    for marker, replacement in replacements.items():
        value = value.replace(marker, replacement)
    return value


def write_preview_consequences(path: Path, rows: list[dict[str, str]]) -> None:
    write_tsv(path, rows, CONSEQUENCE_COLUMNS)


def add_preview_banner(
    html_text: str, candidate_id: str, open_count: int, pending_audit_count: int
) -> str:
    banner = (
        '<div style="position:sticky;top:0;z-index:9999;padding:10px 16px;'
        'background:#7f1d1d;color:white;font:600 14px sans-serif;text-align:center">'
        f'UNSEALED REVISION PREVIEW — candidate {candidate_id} — {open_count} open consequence(s)'
        f' — {pending_audit_count} pending consequence audit(s)'
        "</div>"
    )
    match = re.search(r"<body(?:\s[^>]*)?>", html_text, flags=re.I)
    if match:
        return html_text[: match.end()] + "\n" + banner + html_text[match.end() :]
    return banner + "\n" + html_text


def preview(args: argparse.Namespace) -> int:
    workspace, config, state = load_workspace(args.workspace)
    candidate_id = safe_id(args.candidate or state.get("current_candidate") or "", "candidate_id")
    candidate_root = workspace / "candidates" / candidate_id
    composition = read_json(candidate_root / "composition.json")
    assembly = candidate_root / "assembly"
    current_fingerprint = working_fingerprint(workspace, state, composition["active_batches"])
    if current_fingerprint != composition["working_fingerprint"]:
        raise ControllerError("Candidate is stale relative to its batch overlays; materialize a new candidate")

    renderer = config["preview_renderer"]
    output = assembly / clean_relative(renderer["output"], "preview_renderer.output")
    replacements = {
        "{python}": sys.executable,
        "{candidate_root}": str(assembly),
        "{workspace}": str(workspace),
        "{assembly_root}": config["assembly_root"],
        "{output_path}": str(output),
    }
    argv = [render_value(value, replacements) for value in renderer["argv"]]
    cwd = Path(render_value(renderer["cwd"], replacements))
    environment = os.environ.copy()
    environment.update({key: render_value(value, replacements) for key, value in renderer["environment"].items()})

    preview_id = safe_id(args.preview_id or datetime.now().strftime("%Y%m%dT%H%M%S_preview"), "preview_id")
    preview_root = workspace / "previews" / preview_id
    if preview_root.exists():
        raise ControllerError(f"Preview already exists: {preview_id}")
    preview_root.mkdir(parents=True)
    started = now()
    result = subprocess.run(
        argv,
        cwd=str(cwd),
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=renderer["timeout_seconds"],
    )
    cumulative = read_tsv(workspace / "cumulative_consequences.tsv", CONSEQUENCE_COLUMNS)
    active = set(composition["active_batches"])
    open_rows = [
        row
        for row in cumulative
        if row["triggering_batch"] in active and row["state"] in {"open", "blocked"}
    ]
    pending_audits = [
        batch_id
        for batch_id in composition["active_batches"]
        if read_json(batch_dir(workspace, batch_id) / "batch.json").get("consequence_audit_status")
        not in {"pass", "not_required"}
    ]
    receipt = {
        "schema_version": SCHEMA_VERSION,
        "preview_id": preview_id,
        "candidate_id": candidate_id,
        "baseline_id": state["baseline_id"],
        "active_batches": composition["active_batches"],
        "started_at": started,
        "completed_at": now(),
        "renderer_argv": argv,
        "renderer_cwd": str(cwd),
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "open_consequence_count": len(open_rows),
        "pending_consequence_audits": pending_audits,
        "working_fingerprint": composition["working_fingerprint"],
        "consequence_fingerprint": consequence_fingerprint(cumulative, composition["active_batches"]),
        "dependency_state_fingerprint": dependency_state_fingerprint(
            workspace, cumulative, composition["active_batches"]
        ),
        "sealed": False,
    }
    if result.returncode != 0:
        receipt["status"] = "BLOCK"
        write_json(preview_root / "render_receipt.json", receipt)
        raise ControllerError(f"Preview renderer returned {result.returncode}; see {preview_root / 'render_receipt.json'}")
    if not output.is_file():
        receipt["status"] = "BLOCK"
        write_json(preview_root / "render_receipt.json", receipt)
        raise ControllerError(f"Preview renderer did not create configured output: {output}")

    source_hash = sha256_file(output)
    preview_html = add_preview_banner(
        output.read_text(encoding="utf-8"), candidate_id, len(open_rows), len(pending_audits)
    )
    preview_path = preview_root / "manuscript_preview.html"
    preview_path.write_text(preview_html, encoding="utf-8")
    receipt.update(
        {
            "status": "PASS",
            "renderer_output": str(output),
            "renderer_output_sha256": source_hash,
            "preview_output": str(preview_path),
            "preview_output_sha256": sha256_file(preview_path),
        }
    )
    write_json(preview_root / "render_receipt.json", receipt)
    write_tsv(
        preview_root / "active_batches.tsv",
        [{"batch_id": value} for value in composition["active_batches"]],
        ("batch_id",),
    )
    write_preview_consequences(preview_root / "open_consequences.tsv", open_rows)
    index = read_batch_index(workspace)
    for batch_id in composition["active_batches"]:
        metadata_path = batch_dir(workspace, batch_id) / "batch.json"
        metadata = read_json(metadata_path)
        if metadata.get("state") == "recorded":
            metadata["state"] = "awaiting_review"
            metadata["updated_at"] = now()
            write_json(metadata_path, metadata)
            for row in index:
                if row["batch_id"] == batch_id:
                    row["state"] = "awaiting_review"
    write_batch_index(workspace, index)
    state["latest_preview"] = preview_id
    state["latest_preview_current"] = True
    update_state(workspace, state)
    print(
        json.dumps(
            {
                "status": "PREVIEW_RENDERED",
                "preview_id": preview_id,
                "html": str(preview_path),
                "open_consequences": len(open_rows),
                "pending_consequence_audits": pending_audits,
            },
            indent=2,
        )
    )
    return 0


def resolve_consequence(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    path = workspace / "cumulative_consequences.tsv"
    rows = read_tsv(path, CONSEQUENCE_COLUMNS)
    found: dict[str, str] | None = None
    for row in rows:
        if row["consequence_id"] == args.consequence:
            if row["decision"] == "unresolved" and args.state == "resolved":
                raise ControllerError(
                    "Register a revised consequence audit before resolving an unresolved decision"
                )
            row["state"] = args.state
            row["resolution"] = args.resolution
            found = row
    if not found:
        raise ControllerError(f"Unknown consequence: {args.consequence}")
    write_tsv(path, rows, CONSEQUENCE_COLUMNS)
    delta_path = batch_dir(workspace, found["triggering_batch"]) / "consequence_delta.tsv"
    delta = read_tsv(delta_path, CONSEQUENCE_COLUMNS)
    for row in delta:
        if row["consequence_id"] == args.consequence:
            row["state"] = args.state
            row["resolution"] = args.resolution
    write_tsv(delta_path, delta, CONSEQUENCE_COLUMNS)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(args.consequence)
    return 0


def set_batch_state(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    root = batch_dir(workspace, args.batch)
    metadata = read_json(root / "batch.json")
    if args.state == "complete" and metadata.get("state") not in {
        "recorded",
        "awaiting_review",
        "complete",
    }:
        raise ControllerError("Record the batch before marking it complete")
    if args.state == "complete" and metadata.get("consequence_audit_status") not in {
        "pass",
        "not_required",
    }:
        raise ControllerError("Register a passing consequence audit before marking the batch complete")
    metadata["state"] = args.state
    metadata["updated_at"] = now()
    write_json(root / "batch.json", metadata)
    rows = read_batch_index(workspace)
    for row in rows:
        if row["batch_id"] == args.batch:
            row["state"] = args.state
    write_batch_index(workspace, rows)
    if args.state == "abandoned" and args.batch in state["active_batches"]:
        state["active_batches"].remove(args.batch)
    state["latest_preview_current"] = False
    update_state(workspace, state)
    print(json.dumps({"batch_id": args.batch, "state": args.state}, indent=2))
    return 0


def status(args: argparse.Namespace) -> int:
    workspace, _, state = load_workspace(args.workspace)
    rows = read_batch_index(workspace)
    consequences = read_tsv(workspace / "cumulative_consequences.tsv", CONSEQUENCE_COLUMNS)
    active = set(state["active_batches"])
    open_rows = [
        row
        for row in consequences
        if row["triggering_batch"] in active and row["state"] in {"open", "blocked"}
    ]
    audit_states = {
        batch_id: read_json(batch_dir(workspace, batch_id) / "batch.json").get(
            "consequence_audit_status", "unknown"
        )
        for batch_id in state["active_batches"]
    }
    stale_recorded_batches = []
    for batch_id in state["active_batches"]:
        metadata = read_json(batch_dir(workspace, batch_id) / "batch.json")
        if metadata.get("state") in {"recorded", "awaiting_review", "complete"} and metadata.get(
            "recorded_overlay_fingerprint"
        ) != batch_overlay_fingerprint(workspace, batch_id):
            stale_recorded_batches.append(batch_id)
    preview_current = False
    latest_preview = state.get("latest_preview")
    if latest_preview:
        receipt = read_json(workspace / "previews" / latest_preview / "render_receipt.json")
        preview_current = (
            receipt.get("status") == "PASS"
            and receipt.get("active_batches") == state["active_batches"]
            and receipt.get("working_fingerprint") == working_fingerprint(workspace, state)
            and receipt.get("consequence_fingerprint")
            == consequence_fingerprint(consequences, state["active_batches"])
            and receipt.get("dependency_state_fingerprint")
            == dependency_state_fingerprint(workspace, consequences, state["active_batches"])
        )
    if bool(state.get("latest_preview_current")) != preview_current:
        state["latest_preview_current"] = preview_current
        update_state(workspace, state)
    by_decision = {
        decision: sum(row["decision"] == decision for row in open_rows)
        for decision in ("invalidated", "unresolved")
    }
    payload = {
        "status": "READY",
        "workspace": str(workspace),
        "baseline_id": state["baseline_id"],
        "active_batches": state["active_batches"],
        "batch_states": {row["batch_id"]: row["state"] for row in rows},
        "current_candidate": state.get("current_candidate"),
        "latest_preview": latest_preview,
        "latest_preview_current": preview_current,
        "stale_recorded_batches": stale_recorded_batches,
        "open_consequences": len(open_rows),
        "open_consequences_by_decision": by_decision,
        "consequence_audit_states": audit_states,
        "pending_consequence_audits": sorted(
            batch_id
            for batch_id, audit_state in audit_states.items()
            if audit_state not in {"pass", "not_required"}
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    command = subparsers.add_parser("authenticate", help="verify a baseline without creating a workspace")
    command.add_argument("--config", type=Path, required=True)
    command.set_defaults(func=command_authenticate)

    command = subparsers.add_parser("init", help="authenticate a baseline and initialize its revision workspace")
    command.add_argument("--config", type=Path, required=True)
    command.set_defaults(func=init_workspace)

    command = subparsers.add_parser("begin-batch", help="register a coherent revision request")
    command.add_argument("--workspace", required=True)
    command.add_argument("--prompt-file", required=True)
    command.add_argument("--batch-id")
    command.add_argument(
        "--contract-disposition",
        choices=("preserved", "possibly_changed", "replacement_required"),
        default="preserved",
    )
    command.add_argument(
        "--contract-approval-file",
        help="verbatim user approval required for a non-preserved starting disposition",
    )
    command.add_argument("--allow-overlap-with", action="append", default=[])
    command.set_defaults(func=begin_batch)

    command = subparsers.add_parser("append-prompt", help="append a verbatim clarification to a batch")
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--prompt-file", required=True)
    command.set_defaults(func=append_prompt)

    command = subparsers.add_parser("stage-path", help="copy an effective baseline file or directory into a batch overlay")
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--path", required=True)
    command.set_defaults(func=stage_path)

    command = subparsers.add_parser("register-generated", help="copy a generated file into a batch overlay")
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--source", required=True)
    command.add_argument("--path", required=True)
    command.add_argument("--replace", action="store_true")
    command.set_defaults(func=register_generated)

    command = subparsers.add_parser(
        "register-check",
        help="preserve a batch-local replay, validation, or execution receipt",
    )
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--source", required=True)
    command.add_argument("--name", required=True)
    command.add_argument("--replace", action="store_true")
    command.set_defaults(func=register_check)

    command = subparsers.add_parser(
        "register-consequence-audit",
        help="validate and register an agent-authored graph traversal audit",
    )
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--classification-source", required=True)
    command.add_argument("--audit-source", required=True)
    command.set_defaults(func=register_consequence_audit)

    command = subparsers.add_parser("delete-path", help="register deletion of an effective artifact")
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument("--path", required=True)
    command.set_defaults(func=register_deletion)

    command = subparsers.add_parser(
        "record",
        help="record actual changes and exact-hash backreferences without semantic adjudication",
    )
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.set_defaults(func=record_batch)

    command = subparsers.add_parser(
        "set-contract-disposition",
        help="approve or revise a batch's assembly-contract disposition",
    )
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument(
        "--disposition",
        choices=("preserved", "possibly_changed", "replacement_required"),
        required=True,
    )
    command.add_argument(
        "--approval-file",
        help="verbatim user approval required when leaving the preserved disposition",
    )
    command.set_defaults(func=set_contract_disposition)

    command = subparsers.add_parser("materialize", help="compose a candidate from the baseline and selected batches")
    command.add_argument("--workspace", required=True)
    command.add_argument("--candidate-id")
    command.add_argument("--batches", nargs="*")
    command.add_argument("--additional-reason")
    command.set_defaults(func=materialize)

    command = subparsers.add_parser("preview", help="invoke the configured renderer on a materialized candidate")
    command.add_argument("--workspace", required=True)
    command.add_argument("--candidate")
    command.add_argument("--preview-id")
    command.set_defaults(func=preview)

    command = subparsers.add_parser("resolve", help="record a consequence disposition")
    command.add_argument("--workspace", required=True)
    command.add_argument("--consequence", required=True)
    command.add_argument("--state", choices=("resolved", "accepted_exception", "blocked"), required=True)
    command.add_argument("--resolution", required=True)
    command.set_defaults(func=resolve_consequence)

    command = subparsers.add_parser("set-batch-state", help="mark a recorded batch complete or abandoned")
    command.add_argument("--workspace", required=True)
    command.add_argument("--batch", required=True)
    command.add_argument(
        "--state",
        choices=("recorded", "awaiting_review", "complete", "abandoned"),
        required=True,
    )
    command.set_defaults(func=set_batch_state)

    command = subparsers.add_parser("status", help="authenticate the baseline and report workspace state")
    command.add_argument("--workspace", required=True)
    command.set_defaults(func=status)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ControllerError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    except subprocess.TimeoutExpired as exc:
        print(f"ERROR: preview renderer timed out after {exc.timeout} seconds", file=sys.stderr)
        sys.exit(2)
