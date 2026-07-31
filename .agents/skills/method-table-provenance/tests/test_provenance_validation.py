#!/usr/bin/env python3
"""Regression tests for the panel-leaf provenance validation contract."""

from __future__ import annotations

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from verify_manuscript_endpoints import main as endpoint_main  # noqa: E402
from verify_provenance_locks import main as lock_main  # noqa: E402


COLUMNS = [
    "id",
    "parent",
    "what",
    "why",
    "comment",
    "lock_target",
    "lock_selector",
    "lock_kind",
    "sha256",
    "hash_status",
]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class ProvenanceValidationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        (self.root / "figures").mkdir()
        (self.root / "derived").mkdir()
        self.f1_bytes = b"whole figure one"
        self.f2_bytes = b"whole figure two"
        self.shared_bytes = b"shared upstream table"
        (self.root / "figures/F1.png").write_bytes(self.f1_bytes)
        (self.root / "figures/F2.png").write_bytes(self.f2_bytes)
        (self.root / "derived/shared.csv").write_bytes(self.shared_bytes)
        self.f1_hash = digest(self.f1_bytes)
        self.f2_hash = digest(self.f2_bytes)
        self.shared_hash = digest(self.shared_bytes)
        self.manifest = self.root / "target_figure_set.tsv"
        self.table = self.root / "locked_provenance_table.md"
        self.report = self.root / "report.md"
        self.write_manifest()
        self.rows = self.valid_rows()
        self.write_table(self.rows)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def write_manifest(self) -> None:
        rows = [
            ["F1#panel_a", "figures/F1.png", "#panel_a", self.f1_hash],
            ["F1#panel_b", "figures/F1.png", "#panel_b", self.f1_hash],
            ["F2#panel_a", "figures/F2.png", "#panel_a", self.f2_hash],
        ]
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["endpoint_id", "artifact_path", "lock_selector", "sha256"])
            writer.writerows(rows)

    def valid_rows(self) -> list[list[str]]:
        return [
            [
                "F1#panel_a",
                "action_a",
                "Expose the current manuscript panel.",
                "Anchor the panel to its last Methods-relevant transformation.",
                "Panel endpoint.",
                "figures/F1.png",
                "#panel_a",
                "file",
                self.f1_hash,
                "computed_self",
            ],
            [
                "F1#panel_b",
                "action_b",
                "Expose the current manuscript panel.",
                "Anchor the panel to its last Methods-relevant transformation.",
                "Panel endpoint.",
                "figures/F1.png",
                "#panel_b",
                "file",
                self.f1_hash,
                "computed_self",
            ],
            [
                "F2#panel_a",
                "action_c",
                "Expose the current manuscript panel.",
                "Anchor the panel to its last Methods-relevant transformation.",
                "Panel endpoint.",
                "figures/F2.png",
                "#panel_a",
                "file",
                self.f2_hash,
                "computed_self",
            ],
            [
                "action_a",
                "raw_data",
                "Compute the first displayed quantity.",
                "Define the panel value.",
                "The final figure is a downstream proxy lock.",
                "figures/F1.png",
                "NA",
                "proxy_file",
                self.f1_hash,
                "computed_downstream_proxy",
            ],
            [
                "action_b",
                "raw_data",
                "Compute the second displayed quantity.",
                "Define the panel value.",
                "The same final figure is a valid shared proxy lock.",
                "figures/F1.png",
                "NA",
                "proxy_file",
                self.f1_hash,
                "computed_downstream_proxy",
            ],
            [
                "action_c",
                "raw_data",
                "Compute the third displayed quantity.",
                "Define the panel value.",
                "The final figure is a downstream proxy lock.",
                "figures/F2.png",
                "NA",
                "proxy_file",
                self.f2_hash,
                "computed_downstream_proxy",
            ],
            [
                "raw_data",
                "NA",
                "Load the durable analysis input.",
                "Provide observations for all three transformations.",
                "Current durable input.",
                "derived/shared.csv",
                "NA",
                "file",
                self.shared_hash,
                "computed_self",
            ],
        ]

    def write_table(self, rows: list[list[str]], suffix: str = "") -> None:
        lines = [
            "# Locked provenance",
            "",
            "| " + " | ".join(COLUMNS) + " |",
            "|" + "|".join(["---"] * len(COLUMNS)) + "|",
        ]
        lines.extend("| " + " | ".join(row) + " |" for row in rows)
        if suffix:
            lines.extend(["", suffix])
        self.table.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def endpoint_status(self) -> int:
        return endpoint_main(
            [
                str(self.manifest),
                str(self.table),
                "--root",
                str(self.root),
                "--output",
                str(self.report),
            ]
        )

    def test_valid_panel_graph_and_shared_proxy_locks_pass(self) -> None:
        self.assertEqual(self.endpoint_status(), 0)
        self.assertIn("- Status: `PASS`", self.report.read_text(encoding="utf-8"))
        self.assertEqual(
            lock_main(
                [
                    str(self.table),
                    "--root",
                    str(self.root),
                    "--fail-on-problems",
                    "--output",
                    str(self.report),
                ]
            ),
            0,
        )

    def test_legacy_endpoint_table_is_rejected(self) -> None:
        legacy = "\n".join(
            [
                "| endpoint_id | lock_target | sha256 |",
                "|---|---|---|",
                f"| F1 | figures/F1.png | {self.f1_hash} |",
            ]
        )
        self.write_table(self.rows, legacy)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("legacy three-column endpoint table", self.report.read_text())

    def test_multiple_canonical_tables_are_rejected(self) -> None:
        second_table = "\n".join(
            [
                "| " + " | ".join(COLUMNS) + " |",
                "|" + "|".join(["---"] * len(COLUMNS)) + "|",
                "| old | NA | Old object. | Historical. | NA | NA | NA | terminal | NA | not_applicable |",
            ]
        )
        self.write_table(self.rows, second_table)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("found 2", self.report.read_text())

    def test_one_hash_cannot_identify_two_final_figure_paths(self) -> None:
        (self.root / "figures/F2.png").write_bytes(self.f1_bytes)
        self.f2_hash = self.f1_hash
        self.write_manifest()
        for row in self.rows:
            if row[5] == "figures/F2.png":
                row[8] = self.f1_hash
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("hash identifies multiple final-figure paths", self.report.read_text())

    def test_disconnected_historical_component_is_rejected(self) -> None:
        self.rows.append(
            [
                "legacy_overlay",
                "NA",
                "Retain an old overlay.",
                "Historical archive.",
                "Not current.",
                "NA",
                "NA",
                "terminal",
                "NA",
                "not_applicable",
            ]
        )
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        report = self.report.read_text()
        self.assertIn("graph leaf is absent from current panel scope", report)
        self.assertIn("not upstream-reachable", report)

    def test_duplicate_direct_file_object_is_rejected(self) -> None:
        for row in self.rows:
            if row[0] == "action_c":
                row[1] = "raw_data; raw_data_duplicate"
        self.rows.append(
            [
                "raw_data_duplicate",
                "NA",
                "Duplicate the same durable analysis input.",
                "Invalid duplicate object node.",
                "Must be reconciled.",
                "derived/shared.csv",
                "NA",
                "file",
                self.shared_hash,
                "computed_self",
            ]
        )
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("duplicate direct-file object lock", self.report.read_text())

    def test_duplicate_parent_reference_is_rejected(self) -> None:
        for row in self.rows:
            if row[0] == "action_c":
                row[1] = "raw_data; raw_data"
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn(
            "duplicate parent reference for action_c: raw_data",
            self.report.read_text(),
        )

    def test_cycle_is_rejected(self) -> None:
        for row in self.rows:
            if row[0] == "raw_data":
                row[1] = "action_a"
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("provenance cycle", self.report.read_text())

    def test_panel_selector_is_mandatory(self) -> None:
        for row in self.rows:
            if row[0] == "F1#panel_a":
                row[6] = "NA"
        self.write_table(self.rows)
        self.assertEqual(self.endpoint_status(), 1)
        self.assertIn("endpoint_selector_mismatch", self.report.read_text())


if __name__ == "__main__":
    unittest.main()
