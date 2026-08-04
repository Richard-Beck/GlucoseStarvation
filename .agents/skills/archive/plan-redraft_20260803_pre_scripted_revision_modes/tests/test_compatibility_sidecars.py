#!/usr/bin/env python3
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
CHECKER = SKILL_ROOT / "scripts" / "check_compatibility_sidecars.py"


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


class CompatibilitySidecarTests(unittest.TestCase):
    REQUIRED_CONSUMERS = [
        "analysis",
        "manuscript-figure-workflow",
        "claim-graph-integration",
        "results-text",
        "method-table-provenance",
        "manuscript-legend-writing",
        "serve-manuscript-abstract-introduction-discussion",
    ]

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.repo = Path(self.temp.name)
        self.source = self.repo / "source.txt"
        self.source.write_text("source\n", encoding="utf-8")
        self.source_hash = hashlib.sha256(self.source.read_bytes()).hexdigest()

    def tearDown(self):
        self.temp.cleanup()

    def sidecar(self, consumer, value, status="validated"):
        return {
            "schema_version": 1,
            "consumer": consumer,
            "bundle_id": consumer + "-bundle",
            "completion_attestation": True,
            "facts": [
                {
                    "fact_id": "S13a.qc_pass_count.MDA-MB-231",
                    "value": value,
                    "unit": "selected model fits",
                    "scope": "posterior sampling QC; denominator 5",
                    "status": status,
                    "source_artifacts": [
                        {"path": "source.txt", "sha256": self.source_hash}
                    ],
                    "manuscript_artifacts_using_fact": [],
                    "notes": "",
                }
            ],
        }

    def empty_sidecar(self, consumer):
        return {
            "schema_version": 1,
            "consumer": consumer,
            "bundle_id": consumer + "-bundle",
            "completion_attestation": True,
            "facts": [],
        }

    def run_gate(self, left_value=5, right_value=5, right_status="validated", accepted=None):
        left = self.repo / "figures.json"
        right = self.repo / "legends.json"
        write_json(left, self.sidecar("manuscript-figure-workflow", left_value))
        write_json(right, self.sidecar("manuscript-legend-writing", right_value, right_status))
        sidecars = [
            {"consumer": "manuscript-figure-workflow", "path": "figures.json"},
            {"consumer": "manuscript-legend-writing", "path": "legends.json"},
        ]
        for consumer in self.REQUIRED_CONSUMERS:
            if consumer in {"manuscript-figure-workflow", "manuscript-legend-writing"}:
                continue
            path = consumer + ".json"
            write_json(self.repo / path, self.empty_sidecar(consumer))
            sidecars.append({"consumer": consumer, "path": path})
        registry = self.repo / "registry.json"
        write_json(
            registry,
            {
                "schema_version": 1,
                "sidecars": sidecars,
                "accepted_exceptions": accepted or [],
            },
        )
        report_json = self.repo / "report.json"
        result = subprocess.run(
            [
                "python3",
                str(CHECKER),
                "--registry",
                str(registry),
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
        return result, json.loads(report_json.read_text())

    def test_matching_facts_pass(self):
        result, report = self.run_gate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(report["status"], "PASS")

    def test_conflicting_facts_block(self):
        result, report = self.run_gate(left_value=5, right_value=3)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertEqual(len(report["conflicts"]), 1)

    def test_authorized_exception_warns(self):
        accepted = [
            {
                "fact_id": "S13a.qc_pass_count.MDA-MB-231",
                "reason": "Project owner deferred reconciliation.",
                "authorized_by": "user",
            }
        ]
        result, report = self.run_gate(
            left_value=5,
            right_value=3,
            right_status="unresolved",
            accepted=accepted,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(report["status"], "WARN")
        self.assertEqual(len(report["accepted_exceptions_used"]), 1)

    def test_missing_required_consumer_blocks(self):
        result, report = self.run_gate()
        registry_path = self.repo / "registry.json"
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
        registry["sidecars"] = [
            entry for entry in registry["sidecars"] if entry["consumer"] != "analysis"
        ]
        write_json(registry_path, registry)
        result = subprocess.run(
            [
                "python3",
                str(CHECKER),
                "--registry",
                str(registry_path),
                "--report-json",
                str(self.repo / "missing-report.json"),
                "--report-md",
                str(self.repo / "missing-report.md"),
                "--repo-root",
                str(self.repo),
            ],
            capture_output=True,
            text=True,
        )
        report = json.loads((self.repo / "missing-report.json").read_text(encoding="utf-8"))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("analysis" in error for error in report["errors"]))


if __name__ == "__main__":
    unittest.main()
