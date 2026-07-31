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
        self.graph = self.root / "graph.json"
        self.result_graph = self.root / "result_graph.json"
        self.evidence = self.root / "evidence.json"
        self.index = self.root / "items.json"
        self.audit_index = self.root / "audit_index.json"
        write_json(self.graph, {"claims": [], "evidence": []})
        write_json(self.result_graph, {"claims": [], "evidence": [], "validated": True})
        write_json(self.evidence, {"package": "v1"})
        self.write_index({"F1_a": {"value": 1}, "S13_a": {"value": 5}})

    def tearDown(self):
        self.temp.cleanup()

    def write_index(self, values):
        write_json(
            self.index,
            {
                "schema_version": 1,
                "items": [
                    {"item_id": item_id, "canonical_input": value}
                    for item_id, value in sorted(values.items())
                ],
            },
        )

    def write_audits(self, item_ids):
        audit_dir = self.root / "audits"
        audit_dir.mkdir(exist_ok=True)
        audits = []
        for item_id in item_ids:
            audit_path = audit_dir / (item_id + ".md")
            audit_path.write_text("audit for %s\n" % item_id, encoding="utf-8")
            audits.append({"item_id": item_id, "path": str(audit_path)})
        write_json(self.audit_index, {"schema_version": 1, "audits": audits})

    def snapshot(self, output, graph=None, result_graph=None, with_audits=False):
        command = [
            "python3",
            str(AUTH),
            "snapshot",
            "--starting-graph",
            str(graph or self.graph),
            "--evidence-input",
            "canonical=%s" % self.evidence,
            "--item-index",
            str(self.index),
            "--output",
            str(output),
            "--repo-root",
            str(self.root),
        ]
        if result_graph:
            command.extend(["--result-graph", str(result_graph)])
        if with_audits:
            command.extend(["--audit-index", str(self.audit_index)])
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
        return json.loads(output.read_text())

    def test_identical_inputs_reuse_complete(self):
        prior = self.root / "prior.json"
        current = self.root / "current.json"
        self.write_audits(["F1_a", "S13_a"])
        self.snapshot(prior, result_graph=self.result_graph, with_audits=True)
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_complete")
        self.assertEqual(report["cache_miss_item_ids"], [])

    def test_one_changed_item_reuses_other_item(self):
        prior = self.root / "prior.json"
        self.write_audits(["F1_a", "S13_a"])
        self.snapshot(prior, result_graph=self.result_graph, with_audits=True)
        self.write_index({"F1_a": {"value": 1}, "S13_a": {"value": 3}})
        current = self.root / "current.json"
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_partial")
        self.assertEqual(report["reusable_item_ids"], ["F1_a"])
        self.assertEqual(report["cache_miss_item_ids"], ["S13_a"])

    def test_changed_starting_graph_forces_full_audit(self):
        prior = self.root / "prior.json"
        self.write_audits(["F1_a", "S13_a"])
        self.snapshot(prior, result_graph=self.result_graph, with_audits=True)
        other_graph = self.root / "other_graph.json"
        write_json(other_graph, {"claims": [{"id": "changed"}], "evidence": []})
        current = self.root / "current.json"
        self.snapshot(current, graph=other_graph)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "fresh_full")
        self.assertEqual(report["reusable_item_ids"], [])

    def test_missing_cached_audit_is_a_miss(self):
        prior = self.root / "prior.json"
        current = self.root / "current.json"
        self.write_audits(["F1_a", "S13_a"])
        self.snapshot(prior, result_graph=self.result_graph, with_audits=True)
        (self.root / "audits" / "F1_a.md").unlink()
        self.snapshot(current)
        report = self.compare(current, prior)
        self.assertEqual(report["mode"], "reuse_partial")
        self.assertEqual(report["reusable_item_ids"], ["S13_a"])
        self.assertEqual(report["cache_miss_item_ids"], ["F1_a"])
        self.assertEqual(
            report["invalid_prior_audits"]["F1_a"], "prior audit file is missing"
        )


if __name__ == "__main__":
    unittest.main()
