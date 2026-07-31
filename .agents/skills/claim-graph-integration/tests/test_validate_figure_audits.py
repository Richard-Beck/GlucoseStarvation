#!/usr/bin/env python3

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = SKILL_ROOT / "scripts" / "validate_figure_audits.py"


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


VALID_AUDIT = """### Observation–claim relationships

| Panel(s) | Exact observation | Claim | Relation | Strength | Reason |
|---|---|---|---|---|---|
| e | Higher resources favor low ploidy. | `C11: Abundance selects for low ploidy.` | `support` | `moderate` | The growth contrast is negative. |
"""

EMPTY_AUDIT = """### Observation–claim relationships

| Panel(s) | Exact observation | Claim | Relation | Strength | Reason |
|---|---|---|---|---|---|
| | | | | | |

No material relationships identified.
"""


class FigureAuditValidatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.figures = self.root / "figures.json"
        self.claims = self.root / "claims.json"
        self.audits = self.root / "audits.json"
        self.f5 = self.root / "F5.md"
        self.s1 = self.root / "S1.md"
        write_json(
            self.figures,
            {
                "schema_version": 1,
                "figures": [
                    {"figure_id": "F5", "canonical_input": {}},
                    {"figure_id": "S1", "canonical_input": {}},
                ],
            },
        )
        write_json(
            self.claims,
            {
                "schema_version": 1,
                "claims": [
                    {
                        "id": "C11",
                        "text": "Abundance selects for low ploidy.",
                    }
                ],
            },
        )
        self.f5.write_text(VALID_AUDIT, encoding="utf-8")
        self.s1.write_text(EMPTY_AUDIT, encoding="utf-8")
        write_json(
            self.audits,
            {
                "schema_version": 1,
                "audits": [
                    {"figure_id": "F5", "path": str(self.f5)},
                    {"figure_id": "S1", "path": str(self.s1)},
                ],
            },
        )

    def tearDown(self):
        self.temp.cleanup()

    def run_validator(self):
        return subprocess.run(
            [
                "python3",
                str(VALIDATOR),
                "--figure-index",
                str(self.figures),
                "--claims-index",
                str(self.claims),
                "--audit-index",
                str(self.audits),
                "--repo-root",
                str(self.root),
            ],
            capture_output=True,
            text=True,
        )

    def test_valid_relationship_and_empty_audit_pass(self):
        result = self.run_validator()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("2 figure audits", result.stdout)
        self.assertIn("1 observation-claim relationships", result.stdout)

    def test_qualification_relation_fails(self):
        self.f5.write_text(
            VALID_AUDIT.replace("`support`", "`qualify`"),
            encoding="utf-8",
        )
        result = self.run_validator()
        self.assertEqual(result.returncode, 1)
        self.assertIn("invalid relation", result.stderr)

    def test_changed_authoritative_text_fails(self):
        self.f5.write_text(
            VALID_AUDIT.replace(
                "Abundance selects for low ploidy.",
                "Abundance sometimes selects for low ploidy.",
            ),
            encoding="utf-8",
        )
        result = self.run_validator()
        self.assertEqual(result.returncode, 1)
        self.assertIn("changes authoritative text", result.stderr)


if __name__ == "__main__":
    unittest.main()
