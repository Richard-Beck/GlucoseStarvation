#!/usr/bin/env python3
"""Regression tests for recursive object-registry validation."""

from __future__ import annotations

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from validate_object_registry import EDGE_COLUMNS, OBJECT_COLUMNS, main  # noqa: E402


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class ObjectRegistryValidationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        (self.root / "figures").mkdir()
        (self.root / "scripts").mkdir()
        (self.root / "data").mkdir()
        (self.root / "figures/F1.png").write_bytes(b"figure")
        (self.root / "scripts/plot.R").write_bytes(b"plot code")
        (self.root / "data/raw.csv").write_bytes(b"raw")
        self.registry = self.root / "object_registry.tsv"
        self.edges = self.root / "dependency_edges.tsv"
        self.report = self.root / "report.md"
        self.objects = [
            [
                "F1#panel_a",
                "panel_endpoint",
                "figures/F1.png",
                "#panel_a",
                digest(b"figure"),
                "resolved",
                "Expose the current panel.",
                "Anchor the manuscript endpoint.",
                "current figure manifest",
            ],
            [
                "plot_action",
                "action",
                "scripts/plot.R",
                "::build_panel",
                digest(b"plot code"),
                "resolved",
                "Summarize raw values for plotting.",
                "Define the displayed endpoint.",
                "scripts/plot.R:1",
            ],
            [
                "raw_table",
                "data_file",
                "data/raw.csv",
                "NA",
                digest(b"raw"),
                "terminal",
                "Store raw observations.",
                "Provide the measurement source.",
                "data/raw.csv",
            ],
        ]
        self.dependencies = [
            [
                "F1#panel_a",
                "plot_action",
                "generated_by",
                "confirmed",
                "scripts/plot.R:1",
                "",
            ],
            [
                "plot_action",
                "raw_table",
                "consumes",
                "confirmed",
                "scripts/plot.R:1",
                "",
            ],
        ]
        self.write_inputs()

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def write_inputs(self) -> None:
        with self.registry.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(OBJECT_COLUMNS)
            writer.writerows(self.objects)
        with self.edges.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(EDGE_COLUMNS)
            writer.writerows(self.dependencies)

    def status(self, require_closed: bool = True) -> int:
        args = [
            str(self.registry),
            str(self.edges),
            "--root",
            str(self.root),
            "--output",
            str(self.report),
        ]
        if require_closed:
            args.append("--require-closed")
        return main(args)

    def test_closed_registry_passes(self) -> None:
        self.assertEqual(self.status(), 0)
        self.assertIn("- Status: `PASS`", self.report.read_text())

    def test_duplicate_object_fails(self) -> None:
        self.objects.append(self.objects[-1])
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("duplicate object_id", self.report.read_text())

    def test_unknown_parent_fails(self) -> None:
        self.dependencies[1][1] = "missing_object"
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("unknown parent_id", self.report.read_text())

    def test_invented_object_type_fails(self) -> None:
        self.objects[2][1] = "data_table"
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("invalid object_type", self.report.read_text())

    def test_terminal_dependency_status_fails(self) -> None:
        self.dependencies[1][3] = "terminal"
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("invalid dependency_status", self.report.read_text())

    def test_rejected_dependency_status_fails(self) -> None:
        self.dependencies[1][3] = "display_only"
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("invalid dependency_status", self.report.read_text())

    def test_confirmed_external_scientific_boundary_passes(self) -> None:
        self.objects.append(
            [
                "external_reference_measurements",
                "external",
                "NA",
                "NA",
                "NA",
                "terminal",
                "Provide externally supplied reference measurements.",
                "Anchor an external scientific input boundary.",
                "study accession and retrieval record",
            ]
        )
        self.dependencies.append(
            [
                "plot_action",
                "external_reference_measurements",
                "consumes",
                "confirmed",
                "scripts/plot.R:1",
                "External scientific boundary.",
            ]
        )
        self.write_inputs()
        self.assertEqual(self.status(), 0)

    def test_open_queue_fails_only_when_closed_is_required(self) -> None:
        self.objects[1][5] = "queued"
        self.dependencies[1][3] = "candidate"
        self.write_inputs()
        self.assertEqual(self.status(require_closed=False), 0)
        self.assertEqual(self.status(require_closed=True), 1)
        report = self.report.read_text()
        self.assertIn("queued object remains", report)
        self.assertIn("candidate edge remains", report)

    def test_confirmed_cycle_fails(self) -> None:
        self.dependencies.append(
            [
                "raw_table",
                "plot_action",
                "generated_by",
                "confirmed",
                "scripts/plot.R:1",
                "",
            ]
        )
        self.write_inputs()
        self.assertEqual(self.status(), 1)
        self.assertIn("confirmed dependency cycle", self.report.read_text())

    def test_hash_drift_fails(self) -> None:
        (self.root / "data/raw.csv").write_bytes(b"changed")
        self.assertEqual(self.status(), 1)
        self.assertIn("sha256 mismatch", self.report.read_text())


if __name__ == "__main__":
    unittest.main()
