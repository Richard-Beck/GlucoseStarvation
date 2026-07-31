#!/usr/bin/env python3
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SERVE = SKILL_ROOT / "scripts" / "serve_feedback"
SIGNOFF = SKILL_ROOT / "scripts" / "feedback_signoff"


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_jsonl(path, records):
    path.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in records),
        encoding="utf-8",
    )


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FeedbackSignoffTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.repo = Path(self.temp.name)
        self.raw = self.repo / "raw_feedback.txt"
        self.raw.write_text("Revise the manuscript.\n", encoding="utf-8")
        self.tracking = self.repo / "feedback_tracking" / "review"
        self.tracking.mkdir(parents=True)
        write_jsonl(
            self.tracking / "feedback_spans.jsonl",
            [
                {
                    "span_id": "SPAN-1",
                    "raw_source": "raw_feedback.txt",
                    "start_line": 1,
                    "end_line": 1,
                    "created_by": "feedback-manager",
                    "created_at": "2026-07-13T00:00:00-04:00",
                    "consumer_responses": [],
                }
            ],
        )
        write_jsonl(
            self.tracking / "feedback_items.jsonl",
            [
                {
                    "item_id": "ITEM-1",
                    "source_spans": ["SPAN-1"],
                    "created_by_consumer": "analysis",
                    "created_at": "2026-07-13T00:00:00-04:00",
                    "active": True,
                    "item_update_note": "Canonical item note.",
                    "consumer_responses": [
                        {
                            "response_id": "RESPONSE-1-analysis",
                            "consumer_name": "analysis",
                            "responded_at": "2026-07-13T00:01:00-04:00",
                            "status": "signed_off",
                            "response_type": "scope_satisfied",
                            "response": "Analysis complete.",
                        },
                        {
                            "response_id": "RESPONSE-2-figures",
                            "consumer_name": "manuscript-figure-workflow",
                            "responded_at": "2026-07-13T00:02:00-04:00",
                            "status": "signed_off",
                            "response_type": "scope_satisfied",
                            "response": "Figures complete.",
                        },
                    ],
                }
            ],
        )
        self.served = self.repo / "served"
        subprocess.run(
            [
                str(SERVE),
                "--consumer_name",
                "manuscript-draft-rendering",
                "--tracking_root",
                str(self.tracking),
                "--output_dir",
                str(self.served),
                "--repo_root",
                str(self.repo),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    def tearDown(self):
        self.temp.cleanup()

    def submit(self, payload):
        response_path = self.repo / "response.json"
        write_json(response_path, payload)
        return subprocess.run(
            [str(SIGNOFF), str(self.served), str(response_path)],
            capture_output=True,
            text=True,
        )

    def response_fields(self):
        return {
            "response_type": "scope_satisfied",
            "status": "signed_off",
            "response": "Current consumer scope is complete.",
            "checked_artifacts": [],
            "evidence_artifacts": [],
            "changed_artifacts": [],
            "residual_risk": "",
        }

    def test_generated_form_separates_signoff_from_modification(self):
        template = json.loads(
            (self.served / "feedback_response_template.json").read_text(encoding="utf-8")
        )
        self.assertNotIn("item_updates", template)
        self.assertEqual(template["item_signoffs"], [])
        self.assertEqual(template["item_modifications"], [])
        self.assertFalse(
            {"source_spans", "active", "item_update_note"}.intersection(
                template["item_signoff_template"]
            )
        )
        self.assertEqual(
            set(template["item_modification_template"]["set"]),
            {"source_spans", "active", "item_update_note"},
        )

    def test_legacy_false_modifies_item_note_is_rejected_without_write(self):
        ledger = self.tracking / "feedback_items.jsonl"
        before = sha256(ledger)
        for note in (
            "Renderer replay is complete.",
            "Claim Graph replay is complete.",
        ):
            update = {"item_id": "ITEM-1", "modifies_item": False, "item_update_note": note}
            update.update(self.response_fields())
            result = self.submit(
                {
                    "consumer_name": "manuscript-draft-rendering",
                    "checked_artifacts": [],
                    "span_signoffs": [],
                    "new_items": [],
                    "item_updates": [update],
                }
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("item_updates is no longer supported", result.stderr)
            self.assertEqual(sha256(ledger), before)

    def test_item_signoff_cannot_contain_durable_item_fields(self):
        ledger = self.tracking / "feedback_items.jsonl"
        before = sha256(ledger)
        signoff = {"item_id": "ITEM-1", "item_update_note": "Stage complete."}
        signoff.update(self.response_fields())
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [signoff],
                "item_modifications": [],
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cannot modify item fields", result.stderr)
        self.assertEqual(sha256(ledger), before)

    def test_plain_signoff_preserves_item_and_prior_signoffs(self):
        signoff = {"item_id": "ITEM-1"}
        signoff.update(self.response_fields())
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [signoff],
                "item_modifications": [],
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        item = json.loads((self.tracking / "feedback_items.jsonl").read_text())
        self.assertEqual(item["item_update_note"], "Canonical item note.")
        self.assertEqual(
            [response["status"] for response in item["consumer_responses"]],
            ["signed_off", "signed_off", "signed_off"],
        )

    def test_real_modification_reopens_prior_signoffs(self):
        modification = {
            "item_id": "ITEM-1",
            "modification_reason": "The durable work-unit definition changed.",
            "set": {"item_update_note": "Revised canonical note."},
        }
        modification.update(self.response_fields())
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [],
                "item_modifications": [modification],
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        item = json.loads((self.tracking / "feedback_items.jsonl").read_text())
        self.assertEqual(item["item_update_note"], "Revised canonical note.")
        self.assertEqual(
            [response["status"] for response in item["consumer_responses"]],
            ["reopened", "reopened", "signed_off"],
        )

    def test_noop_modification_is_rejected(self):
        ledger = self.tracking / "feedback_items.jsonl"
        before = sha256(ledger)
        modification = {
            "item_id": "ITEM-1",
            "modification_reason": "Attempted repair.",
            "set": {"item_update_note": "Canonical item note."},
        }
        modification.update(self.response_fields())
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [],
                "item_modifications": [modification],
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not change the durable item record", result.stderr)
        self.assertEqual(sha256(ledger), before)

    def test_item_cannot_be_signed_off_and_modified_together(self):
        ledger = self.tracking / "feedback_items.jsonl"
        before = sha256(ledger)
        signoff = {"item_id": "ITEM-1"}
        signoff.update(self.response_fields())
        modification = {
            "item_id": "ITEM-1",
            "modification_reason": "The durable note changed.",
            "set": {"item_update_note": "Changed note."},
        }
        modification.update(self.response_fields())
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [signoff],
                "item_modifications": [modification],
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cannot be both signed off and modified", result.stderr)
        self.assertEqual(sha256(ledger), before)

    def test_deferred_by_run_policy_is_a_nonmutating_signoff(self):
        fingerprint = "BLOCKER-" + ("a" * 64)
        signoff = {
            "item_id": "ITEM-1",
            "response_type": "deferred_by_run_policy",
            "status": "signed_off",
            "response": "Current consumer work is complete under the routed run policy.",
            "checked_artifacts": [],
            "evidence_artifacts": [],
            "changed_artifacts": [],
            "residual_risk": fingerprint,
        }
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [signoff],
                "item_modifications": [],
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        item = json.loads((self.tracking / "feedback_items.jsonl").read_text())
        self.assertEqual(item["item_update_note"], "Canonical item note.")
        self.assertEqual(item["consumer_responses"][-1]["response_type"], "deferred_by_run_policy")
        self.assertEqual(item["consumer_responses"][-1]["residual_risk"], fingerprint)

    def test_deferred_by_run_policy_requires_full_fingerprint(self):
        ledger = self.tracking / "feedback_items.jsonl"
        before = sha256(ledger)
        signoff = {
            "item_id": "ITEM-1",
            "response_type": "deferred_by_run_policy",
            "status": "signed_off",
            "response": "Deferred.",
            "checked_artifacts": [],
            "evidence_artifacts": [],
            "changed_artifacts": [],
            "residual_risk": "BLOCKER-short",
        }
        result = self.submit(
            {
                "consumer_name": "manuscript-draft-rendering",
                "checked_artifacts": [],
                "span_signoffs": [],
                "new_items": [],
                "item_signoffs": [signoff],
                "item_modifications": [],
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("full BLOCKER-<sha256>", result.stderr)
        self.assertEqual(sha256(ledger), before)


if __name__ == "__main__":
    unittest.main()
