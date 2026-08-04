"""Shared dependency-contract implementation for manuscript assembly tools.

The project-independent entrypoint consumes an assembly-local JSON contract.
``materialize`` builds a reviewable live-dependency candidate manifest from
declared sources. ``validate`` copies the assembly into a sparse temporary
project tree and exposes only locked live and runtime dependencies. Lineage
sources are recorded for audit and deliberately remain unavailable.

The script is compatible with the system Python 3.6 available on the HPC login
nodes.  Commands are represented as JSON argument arrays; no command is passed
through a shell.
"""

from __future__ import print_function

import csv
import grp
import hashlib
import json
import os
import pwd
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


DEPENDENCY_COLUMNS = (
    "dependency_id",
    "source_path",
    "mount_path",
    "path_type",
    "expected_sha256",
    "rationale",
)
ENTRYPOINT_COLUMNS = ("step_id", "argv_json", "expected_exit")
CODE_SUFFIXES = {".py", ".r", ".sh"}
CONFIG_SCHEMA_VERSION = 1


class GateError(RuntimeError):
    pass


def normalized_path(path):
    return Path(os.path.abspath(os.path.normpath(str(path))))


def resolve_existing_path(value, base, field_name):
    raw = Path(value)
    if not raw.is_absolute():
        raw = base / clean_relative_path(value, field_name)
    raw = normalized_path(raw)
    if not raw.exists():
        raise GateError("{} does not exist: {}".format(field_name, raw))
    return raw


def load_config(config_path, project_root):
    config_path = normalized_path(config_path)
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise GateError("Cannot read assembly config {}: {}".format(config_path, exc))
    if config.get("schema_version") != CONFIG_SCHEMA_VERSION:
        raise GateError(
            "Unsupported assembly config schema_version: {!r}".format(
                config.get("schema_version")
            )
        )
    assembly_value = config.get("assembly_root")
    if not isinstance(assembly_value, str) or not assembly_value:
        raise GateError("assembly_root must be a non-empty project-relative path")
    assembly_root = resolve_existing_path(
        assembly_value, project_root, "assembly_root"
    )
    if not assembly_root.is_dir() or not is_relative_to(assembly_root, project_root):
        raise GateError("assembly_root must be a directory inside project_root")
    return config, assembly_root


def assembly_path(config, key, assembly_root, required=True):
    value = config.get(key)
    if value is None and not required:
        return None
    if not isinstance(value, str) or not value:
        raise GateError("{} must be a non-empty assembly-relative path".format(key))
    path = assembly_root / clean_relative_path(value, key)
    if required and not path.exists():
        raise GateError("{} does not exist: {}".format(key, path))
    return path


def read_delimited(path, delimiter, required_columns):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = tuple(reader.fieldnames or ())
        missing = [column for column in required_columns if column not in fields]
        if missing:
            raise GateError(
                "{} is missing column(s): {}".format(path, ", ".join(missing))
            )
        return list(reader)


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def is_relative_to(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def clean_relative_path(value, field_name):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise GateError("{} must be a safe relative path: {}".format(field_name, value))
    if not path.parts:
        raise GateError("{} may not be empty".format(field_name))
    return path


def read_tsv(path, required_columns):
    rows = read_delimited(path, "\t", required_columns)
    if not rows:
        raise GateError("Manifest has no records: {}".format(path))
    return rows


def resolve_source(value, project_root):
    raw = Path(value)
    if not raw.is_absolute():
        raw = project_root / clean_relative_path(value, "source_path")
    if not raw.exists():
        raise GateError("Declared source does not exist: {}".format(raw))
    return raw.resolve()


def resolve_mount(value, project_root):
    raw = Path(value)
    if raw.is_absolute():
        normalized = Path(os.path.normpath(str(raw)))
        if ".." in raw.parts:
            raise GateError("mount_path may not contain '..': {}".format(value))
        return normalized
    return project_root / clean_relative_path(value, "mount_path")


def dependency_records(
    path, dependency_class, project_root, assembly_root, lineage_roots
):
    rows = read_tsv(path, DEPENDENCY_COLUMNS)
    records = []
    ids = set()
    targets = set()
    for row in rows:
        dependency_id = row["dependency_id"].strip()
        if not dependency_id or dependency_id in ids:
            raise GateError("Invalid or duplicate dependency_id: {!r}".format(dependency_id))
        ids.add(dependency_id)
        source = resolve_source(row["source_path"].strip(), project_root)
        target = resolve_mount(row["mount_path"].strip(), project_root)
        path_type = row["path_type"].strip().lower()
        if path_type not in ("file", "directory"):
            raise GateError(
                "{} has unsupported path_type: {}".format(dependency_id, path_type)
            )
        if path_type == "file" and not source.is_file():
            raise GateError("{} is not a file: {}".format(dependency_id, source))
        if path_type == "directory" and not source.is_dir():
            raise GateError("{} is not a directory: {}".format(dependency_id, source))
        target_key = str(target)
        if target_key in targets:
            raise GateError("Duplicate mount target: {}".format(target))
        targets.add(target_key)
        if is_relative_to(target, assembly_root):
            raise GateError(
                "Dependencies may not replace assembly-local content: {}".format(target)
            )
        for lineage_root in lineage_roots:
            if source == lineage_root or is_relative_to(source, lineage_root):
                raise GateError(
                    "{} dependency {} overlaps a declared lineage source: {}".format(
                        dependency_class.capitalize(), dependency_id, source
                    )
                )
        expected = row["expected_sha256"].strip().lower()
        observed = sha256_file(source) if source.is_file() else None
        if expected and expected not in ("-", "na"):
            if path_type != "file":
                raise GateError(
                    "expected_sha256 is supported only for files: {}".format(dependency_id)
                )
            if expected != observed:
                raise GateError(
                    "Checksum mismatch for {}: expected {}, observed {}".format(
                        dependency_id, expected, observed
                    )
                )
        records.append(
            {
                "dependency_id": dependency_id,
                "dependency_class": dependency_class,
                "source_path": str(source),
                "mount_path": str(target),
                "path_type": path_type,
                "expected_sha256": expected or None,
                "observed_sha256": observed,
                "rationale": row["rationale"].strip(),
                "environment": row.get("environment", "").strip(),
            }
        )
    return records


def entrypoint_records(path, replacements):
    rows = read_tsv(path, ENTRYPOINT_COLUMNS)
    records = []
    ids = set()
    for row in rows:
        step_id = row["step_id"].strip()
        if not step_id or step_id in ids:
            raise GateError("Invalid or duplicate step_id: {!r}".format(step_id))
        ids.add(step_id)
        try:
            argv = json.loads(row["argv_json"])
        except ValueError as exc:
            raise GateError("Invalid argv_json for {}: {}".format(step_id, exc))
        if not isinstance(argv, list) or not argv or not all(
            isinstance(value, str) and value for value in argv
        ):
            raise GateError("argv_json must be a non-empty string array: {}".format(step_id))
        expanded = []
        for value in argv:
            for marker, replacement in replacements.items():
                value = value.replace(marker, replacement)
            expanded.append(value)
        try:
            expected_exit = int(row["expected_exit"])
        except ValueError:
            raise GateError("expected_exit must be an integer for {}".format(step_id))
        records.append(
            {"step_id": step_id, "argv": expanded, "expected_exit": expected_exit}
        )
    return records


def prepare_target(sparse_project, project_root, record):
    target = Path(record["mount_path"])
    if is_relative_to(target, project_root):
        relative = target.relative_to(project_root)
        sparse_target = sparse_project / relative
        sparse_target.parent.mkdir(parents=True, exist_ok=True)
        if record["path_type"] == "directory":
            sparse_target.mkdir(parents=True, exist_ok=True)
        elif not sparse_target.exists():
            sparse_target.touch()


def parent_dirs_for_external_target(target, project_root):
    target = Path(target)
    if is_relative_to(target, project_root):
        return []
    dirs = []
    current = target.parent
    while str(current) != "/":
        dirs.append(str(current))
        current = current.parent
    return list(reversed(dirs))


def lineage_records(config, project_root, assembly_root):
    specs = config.get("lineage_sources", [])
    if not isinstance(specs, list):
        raise GateError("lineage_sources must be a list")
    records = []
    seen = set()
    for spec in specs:
        if not isinstance(spec, dict):
            raise GateError("Each lineage_sources record must be an object")
        manifest = assembly_path(spec, "manifest", assembly_root)
        path_column = spec.get("path_column")
        if not isinstance(path_column, str) or not path_column:
            raise GateError("lineage_sources.path_column must be non-empty")
        delimiter = spec.get("delimiter", "\t")
        if delimiter == "tsv":
            delimiter = "\t"
        elif delimiter == "csv":
            delimiter = ","
        rows = read_delimited(manifest, delimiter, (path_column,))
        for row in rows:
            recorded = row[path_column].strip()
            if not recorded or recorded.lower() in ("na", "none"):
                continue
            raw = Path(recorded)
            resolved = normalized_path(raw if raw.is_absolute() else project_root / raw)
            key = str(resolved)
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "recorded_path": recorded,
                    "resolved_path": str(resolved),
                    "manifest": str(manifest),
                }
            )
    return records


def scan_lineage_references(assembly_root, lineage):
    findings = []
    markers = []
    for record in lineage:
        markers.append((record["recorded_path"], record))
        markers.append((record["resolved_path"], record))
    for path in sorted(assembly_root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in CODE_SUFFIXES:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for line_number, line in enumerate(lines, 1):
            for marker, record in markers:
                if marker in line:
                    findings.append(
                        {
                            "path": str(path.relative_to(assembly_root)),
                            "line": line_number,
                            "lineage_path": record["recorded_path"],
                            "text": line.strip()[:500],
                        }
                    )
    return findings


def display_path(path, project_root):
    path = normalized_path(path)
    if is_relative_to(path, project_root):
        return str(path.relative_to(project_root))
    return str(path)


def candidate_path(value, base, project_root, assembly_root):
    raw = Path(value)
    if raw.is_absolute():
        logical = normalized_path(raw)
    else:
        roots = [base, project_root, assembly_root]
        logical = None
        for root in roots:
            test = normalized_path(root / raw)
            if test.exists() or test.is_symlink():
                logical = test
                break
        if logical is None:
            return None
    if not logical.exists() and not logical.is_symlink():
        return None
    if is_relative_to(logical, assembly_root):
        return None
    return logical


def quoted_existing_paths(code_path, project_root, assembly_root):
    try:
        text = code_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    candidates = []
    for match in re.finditer(r"file\.path\(\s*project_root\s*,([^\)]*)\)", text):
        parts = re.findall(r"[\"']([^\"']+)[\"']", match.group(1))
        if parts:
            candidates.append("/".join(parts))
    for value in re.findall(r"[\"']([^\"'\r\n]+)[\"']", text):
        if "/" in value or value.lower().endswith(tuple(CODE_SUFFIXES)):
            candidates.append(value)
    resolved = []
    for value in candidates:
        path = candidate_path(value, code_path.parent, project_root, assembly_root)
        if path is not None and path.resolve().is_file():
            resolved.append(path)
    return resolved


def materialize_candidates(config, project_root, assembly_root, output_override=None):
    materialization = config.get("materialization")
    if not isinstance(materialization, dict):
        raise GateError("Config has no materialization object")
    discovered = {}
    excluded_roots = [
        Path(record["resolved_path"])
        for record in lineage_records(config, project_root, assembly_root)
    ]
    dependency_config = config.get("dependency_manifests", {})
    runtime_manifest = assembly_path(
        dependency_config, "runtime", assembly_root
    )
    for row in read_tsv(runtime_manifest, DEPENDENCY_COLUMNS):
        for field in ("source_path", "mount_path"):
            raw = Path(row[field].strip())
            excluded_roots.append(
                normalized_path(raw if raw.is_absolute() else project_root / raw)
            )

    def add(path, reason):
        path = normalized_path(path)
        if is_relative_to(path, assembly_root):
            return
        resolved = path.resolve()
        for root in excluded_roots:
            if path == root or is_relative_to(path, root):
                return
            if resolved == root or is_relative_to(resolved, root):
                return
        key = str(path)
        discovered.setdefault(key, {"logical": path, "reasons": set()})[
            "reasons"
        ].add(reason)

    table_sources = materialization.get("table_sources", [])
    if not isinstance(table_sources, list):
        raise GateError("materialization.table_sources must be a list")
    for spec in table_sources:
        if not isinstance(spec, dict):
            raise GateError("Each materialization.table_sources record must be an object")
        table = assembly_path(spec, "path", assembly_root)
        columns = spec.get("columns")
        if not isinstance(columns, list) or not columns:
            raise GateError("materialization table source requires columns")
        delimiter = spec.get("delimiter", "\t")
        if delimiter == "tsv":
            delimiter = "\t"
        elif delimiter == "csv":
            delimiter = ","
        separator = spec.get("value_separator", ";")
        base_name = spec.get("base", "project")
        if base_name == "project":
            base = project_root
        elif base_name == "assembly":
            base = assembly_root
        else:
            raise GateError("materialization table source base must be project or assembly")
        for row in read_delimited(table, delimiter, columns):
            for column in columns:
                for value in [item.strip() for item in row[column].split(separator)]:
                    if not value or value.lower() in ("na", "none"):
                        continue
                    path = candidate_path(value, base, project_root, assembly_root)
                    if path is not None:
                        add(path, "{}:{}".format(table.relative_to(assembly_root), column))

    extra = materialization.get("extra_candidates", [])
    if not isinstance(extra, list):
        raise GateError("materialization.extra_candidates must be a list")
    for record in extra:
        if isinstance(record, str):
            value = record
            reason = "config:extra_candidates"
        elif isinstance(record, dict):
            value = record.get("path")
            reason = record.get("rationale", "config:extra_candidates")
        else:
            raise GateError("extra candidate records must be strings or objects")
        if not isinstance(value, str) or not value:
            raise GateError("extra candidate path must be non-empty")
        path = candidate_path(value, project_root, project_root, assembly_root)
        if path is None:
            raise GateError("Extra candidate does not exist: {}".format(value))
        add(path, reason)

    code_queue = []
    for value in materialization.get("code_roots", []):
        root = assembly_root / clean_relative_path(value, "materialization.code_roots")
        if not root.exists():
            raise GateError("Code root does not exist: {}".format(root))
        if root.is_file():
            code_queue.append(root)
        else:
            code_queue.extend(
                path
                for path in sorted(root.rglob("*"))
                if path.is_file() and path.suffix.lower() in CODE_SUFFIXES
            )

    scanned = set()
    while code_queue:
        code_path = normalized_path(code_queue.pop(0))
        if str(code_path) in scanned:
            continue
        scanned.add(str(code_path))
        for path in quoted_existing_paths(code_path, project_root, assembly_root):
            add(path, "code:{}".format(display_path(code_path, project_root)))
            if path.is_file() and path.suffix.lower() in CODE_SUFFIXES:
                code_queue.append(path)

    rows = []
    for index, item in enumerate(
        sorted(discovered.values(), key=lambda value: str(value["logical"])), 1
    ):
        logical = item["logical"]
        source = logical.resolve()
        rows.append(
            {
                "dependency_id": "candidate_{:04d}".format(index),
                "source_path": display_path(source, project_root),
                "mount_path": display_path(logical, project_root),
                "path_type": "directory" if source.is_dir() else "file",
                "expected_sha256": sha256_file(source) if source.is_file() else "",
                "rationale": "; ".join(sorted(item["reasons"])),
            }
        )

    output_value = output_override or materialization.get("candidate_manifest")
    if not output_value:
        raise GateError("Candidate output is not configured; pass --output")
    output = Path(output_value)
    if not output.is_absolute():
        output = assembly_root / clean_relative_path(
            output_value, "materialization.candidate_manifest"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DEPENDENCY_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print("WROTE {} candidate dependencies: {}".format(len(rows), output))
    return 0


def materialize_identity_files(work_root, account):
    identity_root = work_root / "runtime_identity"
    identity_root.mkdir()
    passwd_text = Path("/etc/passwd").read_text(encoding="utf-8")
    if not passwd_text.endswith("\n"):
        passwd_text += "\n"
    uid_marker = ":{}:".format(account.pw_uid)
    if not any(uid_marker in line for line in passwd_text.splitlines()):
        passwd_text += "{}:x:{}:{}:{}:{}:{}\n".format(
            account.pw_name,
            account.pw_uid,
            account.pw_gid,
            account.pw_gecos.replace(":", ""),
            account.pw_dir,
            account.pw_shell or "/bin/bash",
        )
    group = grp.getgrgid(account.pw_gid)
    group_text = Path("/etc/group").read_text(encoding="utf-8")
    if not group_text.endswith("\n"):
        group_text += "\n"
    gid_marker = ":{}:".format(account.pw_gid)
    if not any(gid_marker in line for line in group_text.splitlines()):
        group_text += "{}:x:{}:{}\n".format(
            group.gr_name, account.pw_gid, account.pw_name
        )
    passwd_path = identity_root / "passwd"
    group_path = identity_root / "group"
    passwd_path.write_text(passwd_text, encoding="utf-8")
    group_path.write_text(group_text, encoding="utf-8")
    return passwd_path, group_path


def bwrap_base(
    project_root,
    sparse_project,
    records,
    user_name,
    user_home,
    passwd_path,
    group_path,
):
    command = [
        "bwrap",
        "--unshare-all",
        "--die-with-parent",
        "--new-session",
        "--proc",
        "/proc",
        "--dev",
        "/dev",
        "--dev-bind",
        "/dev/fuse",
        "/dev/fuse",
        "--tmpfs",
        "/tmp",
        "--tmpfs",
        "/run",
    ]

    external_dirs = set()
    project_parent = project_root.parent
    current = project_parent
    while str(current) != "/":
        external_dirs.add(str(current))
        current = current.parent
    for record in records:
        external_dirs.update(
            parent_dirs_for_external_target(record["mount_path"], project_root)
        )
    external_dirs.add("/home")
    external_dirs.add(user_home)
    external_dirs.add("/mnt")
    external_dirs.add("/var/tmp")
    for directory in sorted(external_dirs, key=lambda value: (value.count("/"), value)):
        command.extend(["--dir", directory])

    command.extend(["--bind", str(sparse_project), str(project_root)])
    for record in records:
        command.extend(
            ["--ro-bind", record["source_path"], record["mount_path"]]
        )
    command.extend(["--ro-bind", str(passwd_path), "/etc/passwd"])
    command.extend(["--ro-bind", str(group_path), "/etc/group"])

    environment = {
        "HOME": user_home,
        "USER": user_name,
        "USER_NAME": user_name,
        "LOGNAME": user_name,
        "PATH": "/usr/local/bin:/usr/bin:/bin",
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "TMPDIR": "/tmp",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    for record in records:
        declarations = record.get("environment", "")
        for declaration in [value for value in declarations.split(";") if value]:
            if "=" not in declaration:
                raise GateError(
                    "Runtime environment declaration must use NAME=value: {}".format(
                        declaration
                    )
                )
            key, value = declaration.split("=", 1)
            if not key or not key.replace("_", "A").isalnum() or key[0].isdigit():
                raise GateError("Unsafe runtime environment name: {}".format(key))
            if key in environment and environment[key] != value:
                raise GateError("Conflicting runtime environment value for {}".format(key))
            environment[key] = value
    for key, value in environment.items():
        command.extend(["--setenv", key, value])
    return command


def shortened(text, limit=24000):
    if len(text) <= limit:
        return text
    half = limit // 2
    return text[:half] + "\n... OUTPUT TRUNCATED ...\n" + text[-half:]


def write_reports(json_path, markdown_path, payload):
    json_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Scaffolding-independence validation",
        "",
        "Status: **{}**".format(payload["status"]),
        "",
        "Assembly: `{}`".format(payload["assembly_root"]),
        "",
        "The assembly was copied into a sparse project tree and replayed with only the listed live and runtime dependencies mounted read-only. Declared lineage sources were not mounted.",
        "",
        "## Summary",
        "",
        "- Live dependency mounts: {}".format(payload["live_dependency_count"]),
        "- Runtime mounts: {}".format(payload["runtime_dependency_count"]),
        "- Entrypoints attempted: {} of {}".format(
            len(payload["entrypoint_results"]), payload["entrypoint_count"]
        ),
        "- References to declared lineage sources in executable code: {}".format(
            len(payload["lineage_reference_findings"])
        ),
        "",
    ]
    if payload["lineage_reference_findings"]:
        lines.extend(["## References to lineage sources in code", ""])
        lines.append(
            "These references are audit findings, not failures by themselves. The isolation gate blocks only if a terminal entrypoint requires an unavailable lineage source."
        )
        lines.append("")
        for finding in payload["lineage_reference_findings"]:
            lines.append(
                "- `{path}:{line}` references lineage source `{lineage_path}`: `{text}`".format(**finding)
            )
        lines.append("")
    lines.extend(["## Entrypoints", ""])
    for result in payload["entrypoint_results"]:
        lines.append(
            "- **{status}** `{step_id}` (exit {returncode}; expected {expected_exit}; {elapsed_seconds:.1f} s)".format(
                **result
            )
        )
        if result.get("failure_tail"):
            lines.extend(["", "  Failure tail:", "", "```text"])
            lines.extend(result["failure_tail"].splitlines())
            lines.extend(["```", ""])
    if payload.get("blocked_reason"):
        lines.extend(["", "Blocked reason: {}".format(payload["blocked_reason"]), ""])
    markdown_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def load_context(project_root_value, config_value):
    project_root = Path(project_root_value).resolve()
    if not project_root.is_dir():
        raise GateError("Project root is not a directory: {}".format(project_root))
    config_path = resolve_existing_path(config_value, project_root, "config")
    config, assembly_root = load_config(config_path, project_root)
    return project_root, config_path, config, assembly_root


def validate_independence(
    project_root, config_path, config, assembly_root, timeout_override=None,
    keep_work_root=False
):
    if shutil.which("bwrap") is None:
        raise GateError("bubblewrap (bwrap) is required")
    dependency_config = config.get("dependency_manifests")
    if not isinstance(dependency_config, dict):
        raise GateError("dependency_manifests must be an object")
    isolation = config.get("isolation")
    if not isinstance(isolation, dict):
        raise GateError("isolation must be an object")
    reports = config.get("reports")
    if not isinstance(reports, dict):
        raise GateError("reports must be an object")

    lineage = lineage_records(config, project_root, assembly_root)
    lineage_roots = [Path(record["resolved_path"]) for record in lineage]
    live_manifest = assembly_path(
        dependency_config, "live", assembly_root
    ).resolve()
    runtime_manifest = assembly_path(
        dependency_config, "runtime", assembly_root
    ).resolve()
    entrypoints_manifest = assembly_path(
        isolation, "entrypoints_manifest", assembly_root
    ).resolve()
    report_json = assembly_path(
        reports, "json", assembly_root, required=False
    ).resolve()
    report_md = assembly_path(
        reports, "markdown", assembly_root, required=False
    ).resolve()
    timeout_seconds = timeout_override
    if timeout_seconds is None:
        timeout_seconds = isolation.get("timeout_seconds", 1800)
    if not isinstance(timeout_seconds, int) or timeout_seconds <= 0:
        raise GateError("--timeout-seconds must be positive")

    live = dependency_records(
        live_manifest, "live", project_root, assembly_root, lineage_roots
    )
    runtime = dependency_records(
        runtime_manifest, "runtime", project_root, assembly_root, lineage_roots
    )
    all_records = live + runtime
    targets = [record["mount_path"] for record in all_records]
    if len(targets) != len(set(targets)):
        raise GateError("Live and runtime manifests contain duplicate mount targets")

    assembly_relative = assembly_root.relative_to(project_root)
    lineage_findings = scan_lineage_references(assembly_root, lineage)
    work_root = Path(tempfile.mkdtemp(prefix="manuscript-scaffolding-gate."))
    sparse_project = work_root / "project"
    sparse_project.mkdir()
    sparse_assembly = sparse_project / assembly_relative
    sparse_assembly.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(str(assembly_root), str(sparse_assembly), symlinks=True)
    scratch_root = sparse_project / ".scaffolding_gate"
    scratch_root.mkdir()
    for record in all_records:
        prepare_target(sparse_project, project_root, record)

    replacements = {
        "{project_root}": str(project_root),
        "{assembly_root}": str(project_root / assembly_relative),
        "{scratch_root}": str(project_root / ".scaffolding_gate"),
    }
    entrypoints = entrypoint_records(entrypoints_manifest, replacements)
    account = pwd.getpwuid(os.getuid())
    user_name = account.pw_name
    user_home = account.pw_dir
    passwd_path, group_path = materialize_identity_files(work_root, account)
    base = bwrap_base(
        project_root,
        sparse_project,
        all_records,
        user_name,
        user_home,
        passwd_path,
        group_path,
    )

    results = []
    blocked_reason = None
    try:
        for entrypoint in entrypoints:
            command = base + ["--chdir", str(project_root), "--"] + entrypoint["argv"]
            outer_environment = {
                "PATH": "/usr/local/bin:/usr/bin:/bin",
                "LANG": "C.UTF-8",
                "LC_ALL": "C.UTF-8",
            }
            started = time.time()
            try:
                completed = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    timeout=timeout_seconds,
                    env=outer_environment,
                )
                returncode = completed.returncode
                stdout = completed.stdout
                stderr = completed.stderr
            except subprocess.TimeoutExpired as exc:
                returncode = 124
                stdout = exc.stdout or ""
                stderr = (exc.stderr or "") + "\nEntrypoint timed out."
            elapsed = time.time() - started
            passed = returncode == entrypoint["expected_exit"]
            result = {
                "step_id": entrypoint["step_id"],
                "argv": entrypoint["argv"],
                "expected_exit": entrypoint["expected_exit"],
                "returncode": returncode,
                "elapsed_seconds": round(elapsed, 3),
                "status": "PASS" if passed else "BLOCK",
                "stdout": shortened(stdout),
                "stderr": shortened(stderr),
            }
            if not passed:
                combined = (stderr or stdout).strip().splitlines()
                result["failure_tail"] = "\n".join(combined[-20:])
                blocked_reason = "Entrypoint {} returned {} (expected {})".format(
                    entrypoint["step_id"], returncode, entrypoint["expected_exit"]
                )
            results.append(result)
            if not passed:
                break
    finally:
        payload = {
            "schema_version": 1,
            "status": "PASS" if blocked_reason is None else "BLOCK",
            "assembly_root": str(assembly_root),
            "project_root": str(project_root),
            "isolation_engine": "bubblewrap",
            "network_unshared": True,
            "assembly_copy_writable": True,
            "dependency_mounts_read_only": True,
            "synthetic_runtime_identity": True,
            "assembly_config": str(config_path),
            "live_dependency_count": len(live),
            "runtime_dependency_count": len(runtime),
            "lineage_source_count": len(lineage),
            "entrypoint_count": len(entrypoints),
            "entrypoint_results": results,
            "lineage_sources": lineage,
            "lineage_reference_findings": lineage_findings,
            "blocked_reason": blocked_reason,
            "live_dependencies": live,
            "runtime_dependencies": runtime,
            "temporary_work_root_retained": bool(keep_work_root),
            "temporary_work_root": str(work_root) if keep_work_root else None,
        }
        write_reports(report_json, report_md, payload)
        if not keep_work_root:
            shutil.rmtree(str(work_root))

    print("{}: {}".format(payload["status"], report_md))
    return 0 if payload["status"] == "PASS" else 1
