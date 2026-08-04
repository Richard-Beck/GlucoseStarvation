#!/usr/bin/env python3
"""Disposable end-to-end tests for the in-place revision workspace manager."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


CONTROLLER = Path(__file__).resolve().parents[1] / "scripts/revision_workspace.py"
AUDIT_VALIDATOR = Path(__file__).resolve().parents[1] / "scripts/validate_consequence_audit.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class RevisionWorkspaceTest(unittest.TestCase):
    def run_controller(self, *args: str, controller: Path = CONTROLLER) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(controller), *args],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            self.fail(f"controller failed ({result.returncode}): {result.stderr}\n{result.stdout}")
        return result

    def make_fixture(self, root: Path) -> tuple[Path, Path]:
        assembly = root / "assembly"
        (assembly / "source").mkdir(parents=True)
        (assembly / "draft").mkdir()
        (assembly / "rebuild").mkdir()
        (assembly / "trace").mkdir()
        (assembly / "validation").mkdir()
        (assembly / "source/manuscript.md").write_text("baseline\n", encoding="utf-8")
        (assembly / "draft/manuscript.html").write_text("<html><body>baseline</body></html>\n", encoding="utf-8")
        (assembly / "status.json").write_text('{"status":"PASS"}\n', encoding="utf-8")
        renderer = assembly / "rebuild/render.py"
        renderer.write_text(
            "from pathlib import Path\n"
            "root=Path(__file__).resolve().parents[1]\n"
            "text=(root/'source/manuscript.md').read_text()\n"
            "(root/'draft/manuscript.html').write_text('<html><body>'+text+'</body></html>')\n",
            encoding="utf-8",
        )
        unsealed = assembly / "validation/working_receipt.json"
        unsealed.write_text('{"note":"excluded from canonical seal"}\n', encoding="utf-8")
        (assembly / "trace/source_identity.tsv").write_text(
            "artifact\tsha256\nsource/manuscript.md\t"
            + sha256(assembly / "source/manuscript.md")
            + "\n",
            encoding="utf-8",
        )
        manifest = assembly / "trace/checksums.sha256"
        files = sorted(
            path
            for path in assembly.rglob("*")
            if path.is_file() and path not in {manifest, unsealed}
        )
        manifest.write_text(
            "".join(f"{sha256(path)}  {path.relative_to(assembly)}\n" for path in files),
            encoding="utf-8",
        )
        graph = root / "graph.tsv"
        graph.write_text(
            "upstream_class\tdownstream_class\trelationship\n"
            "manuscript_source\trendered_manuscript\trendered_into\n"
            "rendered_manuscript\tseal_record\tsealed_by\n",
            encoding="utf-8",
        )
        config = root / "config.json"
        config.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "mode": "in_place",
                    "assembly_root": str(assembly),
                    "workspace": "revision_work",
                    "preview_policy": "batch_complete",
                    "baseline": {
                        "checksum_manifest": "trace/checksums.sha256",
                        "status_file": "status.json",
                        "require_complete_inventory": True,
                        "allowed_unsealed_paths": ["validation/working_receipt.json"],
                    },
                    "consequence_audit": {
                        "dependency_graph": str(graph),
                    },
                    "contract_sensitive_patterns": [
                        "rebuild/render.py",
                        "validation/assembly_validation.json",
                        "status.json",
                        "trace/checksums.sha256",
                    ],
                    "preview_renderer": {
                        "argv": ["{python}", "{candidate_root}/rebuild/render.py"],
                        "cwd": "{candidate_root}",
                        "output": "draft/manuscript.html",
                        "environment": {},
                        "timeout_seconds": 30,
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return assembly, config

    def test_end_to_end_overlay_and_preview_preserve_baseline(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_test_") as value:
            root = Path(value)
            assembly, config = self.make_fixture(root)
            baseline_text_hash = sha256(assembly / "source/manuscript.md")
            baseline_draft_hash = sha256(assembly / "draft/manuscript.html")

            auth = json.loads(self.run_controller("authenticate", "--config", str(config)).stdout)
            self.assertEqual(auth["result"], "PASS")
            self.assertEqual(auth["verified_file_count"], 5)

            initialized = json.loads(self.run_controller("init", "--config", str(config)).stdout)
            workspace = Path(initialized["workspace"])
            copied_controller = workspace / "controller/revision_workspace.py"
            prompt = root / "prompt.md"
            prompt.write_text("Change the manuscript text.\n", encoding="utf-8")
            self.run_controller(
                "begin-batch",
                "--workspace",
                str(workspace),
                "--batch-id",
                "edit_text",
                "--prompt-file",
                str(prompt),
                controller=copied_controller,
            )
            staged = Path(
                self.run_controller(
                    "stage-path",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "edit_text",
                    "--path",
                    "source/manuscript.md",
                    controller=copied_controller,
                ).stdout.strip()
            )
            staged.write_text("revised\n", encoding="utf-8")
            recorded = json.loads(
                self.run_controller(
                    "record",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "edit_text",
                    controller=copied_controller,
                ).stdout
            )
            self.assertEqual(recorded["changed_artifacts"], 1)
            self.assertEqual(recorded["consequence_audit"], "PENDING")
            self.assertEqual(recorded["hash_backreference_consumers"], 2)
            references = (workspace / "batches/edit_text/hash_backreferences.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("trace/checksums.sha256", references)
            self.assertIn("trace/source_identity.tsv", references)
            classification = root / "artifact_classification.tsv"
            classification.write_text(
                "artifact_path\tsource\tartifact_class\trationale\n"
                "source/manuscript.md\tchanged_artifact\tmanuscript_source\tThe edited file is the manuscript source.\n"
                "trace/checksums.sha256\thash_backreference_consumer\tseal_record\tThe manifest records sealed artifact hashes.\n"
                "trace/source_identity.tsv\thash_backreference_consumer\tseal_record\tThe identity table records source hashes.\n",
                encoding="utf-8",
            )
            audit = root / "consequence_audit.tsv"
            audit.write_text(
                "upstream_class\tdownstream_class\trelationship\taffected_artifacts\tdecision\trationale\towner\tstate\tresolution\n"
                "manuscript_source\trendered_manuscript\trendered_into\tdraft/manuscript.html\tinvalidated\tThe rendered draft consumes the changed source.\tassembly\topen\t\n"
                "rendered_manuscript\tseal_record\tsealed_by\ttrace/checksums.sha256\tinvalidated\tThe seal records the rendered manuscript state.\tassembly\topen\t\n"
                "manuscript_source\tseal_record\texact_hash_reference\ttrace/checksums.sha256;trace/source_identity.tsv\tinvalidated\tBoth records contain the previous source hash.\tassembly\topen\t\n",
                encoding="utf-8",
            )
            registered_audit = json.loads(
                self.run_controller(
                    "register-consequence-audit",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "edit_text",
                    "--classification-source",
                    str(classification),
                    "--audit-source",
                    str(audit),
                    controller=copied_controller,
                ).stdout
            )
            self.assertEqual(registered_audit["status"], "PASS")
            self.assertEqual(registered_audit["audited_edges"], 3)
            materialized = json.loads(
                self.run_controller(
                    "materialize",
                    "--workspace",
                    str(workspace),
                    "--candidate-id",
                    "candidate_1",
                    controller=copied_controller,
                ).stdout
            )
            candidate = Path(materialized["assembly"])
            self.assertEqual((candidate / "source/manuscript.md").read_text(), "revised\n")
            rendered = json.loads(
                self.run_controller(
                    "preview",
                    "--workspace",
                    str(workspace),
                    "--candidate",
                    "candidate_1",
                    "--preview-id",
                    "preview_1",
                    controller=copied_controller,
                ).stdout
            )
            preview = Path(rendered["html"])
            self.assertIn("UNSEALED REVISION PREVIEW", preview.read_text())
            self.assertIn("revised", preview.read_text())
            state = json.loads(
                self.run_controller("status", "--workspace", str(workspace), controller=copied_controller).stdout
            )
            self.assertTrue(state["latest_preview_current"])
            self.assertEqual(state["open_consequences"], 3)
            self.assertEqual(state["consequence_audit_states"]["edit_text"], "pass")
            self.assertEqual(state["batch_states"]["edit_text"], "awaiting_review")

            receipt_source = root / "figure_replay.tsv"
            receipt_source.write_text("figure\tsha256\nF1\tabc\n", encoding="utf-8")
            receipt_path = Path(
                self.run_controller(
                    "register-check",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "edit_text",
                    "--source",
                    str(receipt_source),
                    "--name",
                    "figure_replay.tsv",
                    controller=copied_controller,
                ).stdout.strip()
            )
            self.assertEqual(receipt_path.read_text(encoding="utf-8"), receipt_source.read_text(encoding="utf-8"))
            state = json.loads(
                self.run_controller("status", "--workspace", str(workspace), controller=copied_controller).stdout
            )
            self.assertTrue(state["latest_preview_current"])

            duplicate_candidate = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "materialize",
                    "--workspace",
                    str(workspace),
                    "--candidate-id",
                    "candidate_duplicate",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(duplicate_candidate.returncode, 2)
            self.assertIn("--additional-reason", duplicate_candidate.stderr)

            staged.write_text("revised again\n", encoding="utf-8")
            state = json.loads(
                self.run_controller("status", "--workspace", str(workspace), controller=copied_controller).stdout
            )
            self.assertFalse(state["latest_preview_current"])
            self.assertEqual(state["stale_recorded_batches"], ["edit_text"])
            stale = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "materialize",
                    "--workspace",
                    str(workspace),
                    "--candidate-id",
                    "candidate_stale",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(stale.returncode, 2)
            self.assertIn("re-run record", stale.stderr)

            self.assertEqual(sha256(assembly / "source/manuscript.md"), baseline_text_hash)
            self.assertEqual(sha256(assembly / "draft/manuscript.html"), baseline_draft_hash)

    def test_workspace_rejects_metadata_drift_in_checksum_excluded_baseline_file(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_test_") as value:
            root = Path(value)
            assembly, config = self.make_fixture(root)
            initialized = json.loads(self.run_controller("init", "--config", str(config)).stdout)
            workspace = Path(initialized["workspace"])
            copied_controller = workspace / "controller/revision_workspace.py"
            unsealed = assembly / "validation/working_receipt.json"
            before = unsealed.stat()
            os.utime(unsealed, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
            result = subprocess.run(
                [sys.executable, str(copied_controller), "status", "--workspace", str(workspace)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("Immutable baseline snapshot changed", result.stderr)
            self.assertIn("mtime_ns", result.stderr)

    def test_contract_sensitive_change_requires_explicit_disposition(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_test_") as value:
            root = Path(value)
            _, config = self.make_fixture(root)
            initialized = json.loads(self.run_controller("init", "--config", str(config)).stdout)
            workspace = Path(initialized["workspace"])
            copied_controller = workspace / "controller/revision_workspace.py"
            prompt = root / "prompt.md"
            prompt.write_text("Make a narrow source edit.\n", encoding="utf-8")
            unapproved_start = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "begin-batch",
                    "--workspace",
                    str(workspace),
                    "--batch-id",
                    "unapproved_start",
                    "--prompt-file",
                    str(prompt),
                    "--contract-disposition",
                    "possibly_changed",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(unapproved_start.returncode, 2)
            self.assertIn("requires --approval-file", unapproved_start.stderr)
            self.run_controller(
                "begin-batch",
                "--workspace",
                str(workspace),
                "--batch-id",
                "scope_gate",
                "--prompt-file",
                str(prompt),
                controller=copied_controller,
            )
            staged = Path(
                self.run_controller(
                    "stage-path",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "scope_gate",
                    "--path",
                    "rebuild/render.py",
                    controller=copied_controller,
                ).stdout.strip()
            )
            staged.write_text(staged.read_text(encoding="utf-8") + "# changed\n", encoding="utf-8")
            blocked = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "record",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "scope_gate",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(blocked.returncode, 2)
            self.assertIn("Contract-sensitive scope expansion", blocked.stderr)
            missing_approval = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "set-contract-disposition",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "scope_gate",
                    "--disposition",
                    "possibly_changed",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(missing_approval.returncode, 2)
            self.assertIn("requires --approval-file", missing_approval.stderr)
            approval = root / "approval.md"
            approval.write_text("User approves exploring a possible contract change.\n", encoding="utf-8")
            self.run_controller(
                "set-contract-disposition",
                "--workspace",
                str(workspace),
                "--batch",
                "scope_gate",
                "--disposition",
                "possibly_changed",
                "--approval-file",
                str(approval),
                controller=copied_controller,
            )
            copied_approvals = list((workspace / "batches/scope_gate/decisions").glob("*.md"))
            self.assertEqual(len(copied_approvals), 1)
            self.assertEqual(copied_approvals[0].read_text(encoding="utf-8"), approval.read_text(encoding="utf-8"))
            recorded = json.loads(
                self.run_controller(
                    "record",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "scope_gate",
                    controller=copied_controller,
                ).stdout
            )
            self.assertEqual(recorded["changed_artifacts"], 1)

    def test_configuration_requires_exact_contract_sensitive_paths(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_test_") as value:
            root = Path(value)
            _, config = self.make_fixture(root)
            payload = json.loads(config.read_text(encoding="utf-8"))
            payload["contract_sensitive_patterns"] = ["rebuild/*"]
            config.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(CONTROLLER), "authenticate", "--config", str(config)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("exact control-plane paths", result.stderr)

    def test_batch_completion_requires_a_passing_agent_audit(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_complete_test_") as value:
            root = Path(value)
            _, config = self.make_fixture(root)
            initialized = json.loads(self.run_controller("init", "--config", str(config)).stdout)
            workspace = Path(initialized["workspace"])
            copied_controller = workspace / "controller/revision_workspace.py"
            prompt = root / "prompt.md"
            prompt.write_text("Revise the source.\n", encoding="utf-8")
            self.run_controller(
                "begin-batch",
                "--workspace",
                str(workspace),
                "--batch-id",
                "audit_gate",
                "--prompt-file",
                str(prompt),
                controller=copied_controller,
            )
            staged = Path(
                self.run_controller(
                    "stage-path",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "audit_gate",
                    "--path",
                    "source/manuscript.md",
                    controller=copied_controller,
                ).stdout.strip()
            )
            staged.write_text("revised\n", encoding="utf-8")
            self.run_controller(
                "record",
                "--workspace",
                str(workspace),
                "--batch",
                "audit_gate",
                controller=copied_controller,
            )
            blocked = subprocess.run(
                [
                    sys.executable,
                    str(copied_controller),
                    "set-batch-state",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "audit_gate",
                    "--state",
                    "complete",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(blocked.returncode, 2)
            self.assertIn("passing consequence audit", blocked.stderr)

            classification = root / "classification.tsv"
            classification.write_text(
                "artifact_path\tsource\tartifact_class\trationale\n"
                "source/manuscript.md\tchanged_artifact\tmanuscript_source\tThe edited file is manuscript source.\n"
                "trace/checksums.sha256\thash_backreference_consumer\tseal_record\tThe checksum manifest is a seal record.\n"
                "trace/source_identity.tsv\thash_backreference_consumer\tseal_record\tThe identity table is a seal record.\n",
                encoding="utf-8",
            )
            audit = root / "audit.tsv"
            audit.write_text(
                "upstream_class\tdownstream_class\trelationship\taffected_artifacts\tdecision\trationale\towner\tstate\tresolution\n"
                "manuscript_source\trendered_manuscript\trendered_into\tdraft/manuscript.html\tinvalidated\tThe draft consumes the changed source.\tassembly\topen\t\n"
                "rendered_manuscript\tseal_record\tsealed_by\ttrace/checksums.sha256\tinvalidated\tThe seal records the rendered manuscript state.\tassembly\topen\t\n"
                "manuscript_source\tseal_record\texact_hash_reference\ttrace/checksums.sha256;trace/source_identity.tsv\tinvalidated\tThe records contain the previous hash.\tassembly\topen\t\n",
                encoding="utf-8",
            )
            self.run_controller(
                "register-consequence-audit",
                "--workspace",
                str(workspace),
                "--batch",
                "audit_gate",
                "--classification-source",
                str(classification),
                "--audit-source",
                str(audit),
                controller=copied_controller,
            )
            completed = json.loads(
                self.run_controller(
                    "set-batch-state",
                    "--workspace",
                    str(workspace),
                    "--batch",
                    "audit_gate",
                    "--state",
                    "complete",
                    controller=copied_controller,
                ).stdout
            )
            self.assertEqual(completed["state"], "complete")

    def test_audit_terminates_at_remains_valid_and_requires_invalidated_descendants(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_audit_test_") as value:
            root = Path(value)
            graph = root / "graph.tsv"
            graph.write_text(
                "upstream_class\tdownstream_class\trelationship\n"
                "source\tinterpretation\tdescribed_by\n"
                "interpretation\tresults\tdescribed_in\n",
                encoding="utf-8",
            )
            classifications = root / "classifications.tsv"
            classifications.write_text(
                "artifact_path\tsource\tartifact_class\trationale\n"
                "source/F2.png\tchanged_artifact\tsource\tThis changed file is a figure source artifact.\n",
                encoding="utf-8",
            )
            changes = root / "changes.tsv"
            changes.write_text(
                "artifact_path\tkind\tchange_type\tparent_sha256\tcurrent_sha256\trecorded_at\n"
                "source/F2.png\tgenerated\tmodified\told\tnew\tnow\n",
                encoding="utf-8",
            )
            references = root / "references.tsv"
            references.write_text(
                "consumer_path\tparent_sha256\tchanged_artifacts\n",
                encoding="utf-8",
            )
            audit = root / "audit.tsv"
            audit.write_text(
                "upstream_class\tdownstream_class\trelationship\taffected_artifacts\tdecision\trationale\towner\tstate\tresolution\n"
                "source\tinterpretation\tdescribed_by\t\tremains_valid\tThe visual meaning is unchanged.\tfigures\tresolved\t\n",
                encoding="utf-8",
            )
            command = [
                sys.executable,
                str(AUDIT_VALIDATOR),
                "--graph",
                str(graph),
                "--classifications",
                str(classifications),
                "--changes",
                str(changes),
                "--hash-references",
                str(references),
                "--audit",
                str(audit),
            ]
            passed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(passed.returncode, 0, passed.stderr)
            self.assertEqual(json.loads(passed.stdout)["audited_edges"], 1)

            audit.write_text(
                "upstream_class\tdownstream_class\trelationship\taffected_artifacts\tdecision\trationale\towner\tstate\tresolution\n"
                "source\tinterpretation\tdescribed_by\t\tinvalidated\tThe interpretation must change.\tfigures\topen\t\n",
                encoding="utf-8",
            )
            blocked = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(blocked.returncode, 2)
            self.assertIn("interpretation/results/described_in", blocked.stderr)

    def test_authentication_rejects_a_mutated_baseline(self) -> None:
        with tempfile.TemporaryDirectory(prefix="plan_redraft_test_") as value:
            root = Path(value)
            assembly, config = self.make_fixture(root)
            (assembly / "source/manuscript.md").write_text("mutated\n", encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(CONTROLLER), "authenticate", "--config", str(config)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("Baseline checksum authentication failed", result.stderr)


if __name__ == "__main__":
    unittest.main()
