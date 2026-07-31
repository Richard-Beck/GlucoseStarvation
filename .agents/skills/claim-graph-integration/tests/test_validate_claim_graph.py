#!/usr/bin/env python3

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = SKILL_ROOT / "scripts" / "validate_claim_graph.py"
PLOTTER = SKILL_ROOT / "scripts" / "plot_claim_graph.py"


def valid_graph():
    return {
        "metadata": {
            "schema_version": "claim-graph/v4",
            "relation_values": ["support", "undermine"],
            "strength_values": ["strong", "moderate", "weak"],
            "authoritative_claim_contract": {
                "claims": [
                    {
                        "id": "C11",
                        "text": "Abundance selects for low ploidy.",
                    }
                ]
            },
        },
        "claims": [
            {
                "id": "C11",
                "text": "Abundance selects for low ploidy.",
                "user_fixed": True,
            }
        ],
        "evidence": [
            {
                "id": "E_F5_e_1",
                "figure_id": "F5",
                "panels": ["e"],
                "observation": "Higher-resource regions favor low ploidy.",
                "source": "Figure F5 semantic audit.",
            }
        ],
        "relationships": [
            {
                "id": "R1",
                "source_type": "evidence",
                "source_id": "E_F5_e_1",
                "target_claim_id": "C11",
                "relation": "support",
                "strength": "moderate",
                "reason": "The modeled growth contrast is negative.",
            }
        ],
    }


class ClaimGraphValidatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.graph = self.root / "graph.json"

    def tearDown(self):
        self.temp.cleanup()

    def run_validator(self, payload):
        self.graph.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return subprocess.run(
            ["python3", str(VALIDATOR), str(self.graph)],
            capture_output=True,
            text=True,
        )

    def test_valid_clean_graph_passes(self):
        result = self.run_validator(valid_graph())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("1 relationships", result.stdout)

    def test_qualification_relation_fails(self):
        payload = valid_graph()
        payload["relationships"][0]["relation"] = "qualify"
        result = self.run_validator(payload)
        self.assertEqual(result.returncode, 1)
        self.assertIn("must be support or undermine", result.stderr)

    def test_embedded_qualification_fields_fail(self):
        payload = valid_graph()
        payload["claims"][0]["qualifies"] = []
        result = self.run_validator(payload)
        self.assertEqual(result.returncode, 1)
        self.assertIn("obsolete embedded relationship fields", result.stderr)

    def test_unlinked_evidence_fails(self):
        payload = valid_graph()
        payload["relationships"] = []
        result = self.run_validator(payload)
        self.assertEqual(result.returncode, 1)
        self.assertIn("evidence nodes without relationships", result.stderr)

    @unittest.skipUnless(shutil.which("dot"), "Graphviz dot is unavailable")
    def test_valid_clean_graph_renders(self):
        payload = valid_graph()
        self.graph.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        dot_path = self.root / "graph.dot"
        png_path = self.root / "graph.png"
        result = subprocess.run(
            [
                "python3",
                str(PLOTTER),
                "--graph",
                str(self.graph),
                "--dot-output",
                str(dot_path),
                "--output",
                str(png_path),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(dot_path.is_file())
        self.assertGreater(png_path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
