#!/usr/bin/env python3
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_ROOT / "scripts" / "reconcile_user_blockers.py"


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


class UserBlockerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.repo = Path(self.temp.name)
        self.source = self.repo / "evidence.json"
        self.source.write_text('{"value": 1}\n', encoding="utf-8")
        self.sidecar = self.repo / "user_blockers.json"
        self.policy = self.repo / "run_policy.json"
        self.ledger = self.repo / "ledger.json"
        self.action = self.repo / "user_action_required.json"

    def tearDown(self):
        self.temp.cleanup()

    def source_hash(self):
        return hashlib.sha256(self.source.read_bytes()).hexdigest()

    def write_policy(self, mode):
        write_json(
            self.policy,
            {"schema_version": 1, "user_blocker_mode": mode, "source": "test"},
        )

    def write_sidecar(self, question="Choose a status.", sha=None, blockers=True):
        entries = []
        if blockers:
            entries.append(
                {
                    "blocker_key": "S13a.choose-support-status",
                    "question": question,
                    "blocking_scope": "user-decision",
                    "relevant_artifacts": [
                        {"path": "evidence.json", "sha256": sha or self.source_hash()}
                    ],
                    "provisional_handling": "Retain the current status pending review.",
                    "options": ["Retain", "Change"],
                }
            )
        write_json(
            self.sidecar,
            {
                "schema_version": 1,
                "consumer": "claim-graph-integration",
                "bundle_id": "bundle-1",
                "completion_attestation": True,
                "blockers": entries,
            },
        )

    def reconcile(self, check=True):
        report_json = self.repo / "report.json"
        result = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "reconcile",
                "--sidecar",
                str(self.sidecar),
                "--ledger",
                str(self.ledger),
                "--run-policy",
                str(self.policy),
                "--user-action-file",
                str(self.action),
                "--report-json",
                str(report_json),
                "--report-md",
                str(self.repo / "report.md"),
                "--repo-root",
                str(self.repo),
            ],
            capture_output=True,
            text=True,
        )
        if check:
            self.assertEqual(result.returncode, 0, result.stderr)
            return result, json.loads(report_json.read_text(encoding="utf-8"))
        return result, None

    def fingerprint(self):
        output = self.repo / "fingerprints.json"
        result = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "fingerprint",
                "--sidecar",
                str(self.sidecar),
                "--output",
                str(output),
                "--repo-root",
                str(self.repo),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(output.read_text(encoding="utf-8"))

    def test_defer_mode_deduplicates_reworded_question(self):
        self.write_policy("defer_to_assembly")
        self.write_sidecar()
        _, first = self.reconcile()
        self.assertEqual(first["status"], "CONTINUE")
        fingerprint = first["new_fingerprints"][0]
        self.assertEqual(json.loads(self.action.read_text())["status"], "clear")
        self.write_sidecar(question="Please choose the support status.")
        _, second = self.reconcile()
        self.assertEqual(second["repeated_fingerprints"], [fingerprint])
        ledger = json.loads(self.ledger.read_text())
        self.assertEqual(len(ledger["blockers"]), 1)
        self.assertEqual(ledger["blockers"][0]["occurrence_count"], 2)
        self.assertEqual(ledger["blockers"][0]["disposition"], "deferred_to_assembly")

    def test_consumer_fingerprint_matches_planner_reconciliation(self):
        self.write_policy("defer_to_assembly")
        self.write_sidecar()
        fingerprint = self.fingerprint()["blockers"][0]["fingerprint"]
        _, report = self.reconcile()
        self.assertEqual(report["new_fingerprints"], [fingerprint])

    def test_changed_artifact_creates_new_fingerprint(self):
        self.write_policy("defer_to_assembly")
        self.write_sidecar()
        _, first = self.reconcile()
        self.source.write_text('{"value": 2}\n', encoding="utf-8")
        self.write_sidecar()
        _, second = self.reconcile()
        self.assertNotEqual(first["new_fingerprints"], second["new_fingerprints"])
        self.assertEqual(len(json.loads(self.ledger.read_text())["blockers"]), 2)

    def test_interactive_mode_waits_once_and_resolution_clears(self):
        self.write_policy("interactive")
        self.write_sidecar()
        _, first = self.reconcile()
        fingerprint = first["new_fingerprints"][0]
        self.assertEqual(first["status"], "USER_ACTION_REQUIRED")
        self.assertEqual(json.loads(self.action.read_text())["status"], "awaiting_user")
        _, second = self.reconcile()
        self.assertEqual(second["new_fingerprints"], [])
        self.assertEqual(second["repeated_fingerprints"], [fingerprint])
        resolution = self.repo / "resolution.json"
        write_json(
            resolution,
            {
                "schema_version": 1,
                "resolutions": [
                    {
                        "fingerprint": fingerprint,
                        "authorized_by": "user",
                        "disposition": "resolved",
                        "resolution": "Retain the current support status.",
                    }
                ],
            },
        )
        result = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "resolve",
                "--ledger",
                str(self.ledger),
                "--resolution-file",
                str(resolution),
                "--user-action-file",
                str(self.action),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(self.action.read_text())["status"], "clear")
        self.assertEqual(json.loads(self.ledger.read_text())["blockers"][0]["disposition"], "resolved")

    def test_bad_artifact_checksum_rejects_without_ledger(self):
        self.write_policy("interactive")
        self.write_sidecar(sha="0" * 64)
        result, _ = self.reconcile(check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum mismatch", result.stderr)
        self.assertFalse(self.ledger.exists())

    def test_defer_mode_converts_prior_awaiting_blocker(self):
        self.write_policy("interactive")
        self.write_sidecar()
        self.reconcile()
        self.write_policy("defer_to_assembly")
        self.write_sidecar(blockers=False)
        _, report = self.reconcile()
        self.assertEqual(report["status"], "CONTINUE")
        self.assertEqual(json.loads(self.action.read_text())["status"], "clear")
        self.assertEqual(
            json.loads(self.ledger.read_text())["blockers"][0]["disposition"],
            "deferred_to_assembly",
        )


if __name__ == "__main__":
    unittest.main()
