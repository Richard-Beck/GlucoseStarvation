#!/usr/bin/env python3

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
AUTH = SKILL_ROOT / "scripts" / "authenticate_claim_graph_inputs.py"


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


class InputAuthenticationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.figure_index = self.root / "figures.json"
        self.claims_index = self.root / "claims.json"
        self.audit_contract = self.root / "prompt.md"
        self.evidence = self.root / "evidence.json"
        self.audit_index = self.root / "audit_index.json"
        self.result_graph = self.root / "claim_graph.json"
        self.audit_contract.write_text("prompt version 1\n", encoding="utf-8")
        write_json(self.evidence, {"package": "v1"})
        write_json(self.result_graph, {"metadata": {"schema_version": "claim-graph/v4"}})
        self.write_figures({"F1": {"value": 1}, "F5": {"value": 5}})
        self.write_claims(
            [
                {"id": "C0", "text": "Claim zero."},
                {"id": "C1", "text": "Claim one."},
            ]
        )

    def tearDown(self):
        self.temp.cleanup()

    def write_figures(self, values):
        write_json(
            self.figure_index,
            {
                "schema_version": 1,
                "figures": [
                    {"figure_id": figure_id, "canonical_input": value}
                    for figure_id, value in sorted(values.items())
                ],
            },
        )

    def write_claims(self, claims):
        write_json(
            self.claims_index,
            {"schema_version": 1, "claims": claims},
        )

    def write_audits(self, figure_ids):
        audit_dir = self.root / "audits"
        audit_dir.mkdir(exist_ok=True)
        audits = []
        for figure_id in figure_ids:
            path = audit_dir / f"{figure_id}.md"
            path.write_text(f"audit for {figure_id}\n", encoding="utf-8")
            audits.append({"figure_id": figure_id, "path": str(path)})
        write_json(self.audit_index, {"schema_version": 1, "audits": audits})

    def snapshot(self, output, with_audits=False, result_graph=None):
        command = [
            "python3",
            str(AUTH),
            "snapshot",
            "--figure-index",
            str(self.figure_index),
            "--claims-index",
            str(self.claims_index),
            "--audit-contract",
            str(self.audit_contract),
            "--evidence-input",
            f"canonical={self.evidence}",
            "--output",
            str(output),
            "--repo-root",
            str(self.root),
        ]
        if with_audits:
            command.extend(["--audit-index", str(self.audit_index)])
        if result_graph:
            command.extend(["--result-graph", str(result_graph)])
        subprocess.run(command, check=True, capture_output=True, text=True)

    def compare(self, current, prior):
        output = self.root / "comparison.json"
        subprocess.run(
            [
                "python3",
                str(AUTH),
                "compare",
                "--current",
                str(current),
                "--prior",
                str(prior),
                "--output",
                str(output),
                "--repo-root",
                str(self.root),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(output.read_text(encoding="utf-8"))

    def test_identical_inputs_reuse_complete(self):
        prior = self.root / "prior.json"
        current = self.root / "current.json"
        self.write_audits(["F1", "F5"])
        self.snapshot(
            prior,
            with_audits=True,
            result_graph=self.result_graph,
        )
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_complete")
        self.assertEqual(report["cache_miss_figure_ids"], [])

    def test_one_changed_figure_reuses_other_figure(self):
        prior = self.root / "prior.json"
        self.write_audits(["F1", "F5"])
        self.snapshot(
            prior,
            with_audits=True,
            result_graph=self.result_graph,
        )
        self.write_figures({"F1": {"value": 1}, "F5": {"value": 6}})
        current = self.root / "current.json"
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_partial")
        self.assertEqual(report["reusable_figure_ids"], ["F1"])
        self.assertEqual(report["cache_miss_figure_ids"], ["F5"])

    def test_changed_claims_invalidate_every_figure_audit(self):
        prior = self.root / "prior.json"
        self.write_audits(["F1", "F5"])
        self.snapshot(
            prior,
            with_audits=True,
            result_graph=self.result_graph,
        )
        self.write_claims(
            [
                {"id": "C0", "text": "Changed claim zero."},
                {"id": "C1", "text": "Claim one."},
            ]
        )
        current = self.root / "current.json"
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "fresh_full")
        self.assertEqual(report["cache_miss_figure_ids"], ["F1", "F5"])

    def test_changed_prompt_contract_invalidates_every_audit(self):
        prior = self.root / "prior.json"
        self.write_audits(["F1", "F5"])
        self.snapshot(
            prior,
            with_audits=True,
            result_graph=self.result_graph,
        )
        self.audit_contract.write_text("prompt version 2\n", encoding="utf-8")
        current = self.root / "current.json"
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "fresh_full")
        self.assertFalse(report["audit_context_match"])

    def test_missing_cached_audit_is_a_partial_miss(self):
        prior = self.root / "prior.json"
        current = self.root / "current.json"
        self.write_audits(["F1", "F5"])
        self.snapshot(
            prior,
            with_audits=True,
            result_graph=self.result_graph,
        )
        (self.root / "audits" / "F1.md").unlink()
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_partial")
        self.assertEqual(report["reusable_figure_ids"], ["F5"])
        self.assertEqual(report["cache_miss_figure_ids"], ["F1"])


if __name__ == "__main__":
    unittest.main()
