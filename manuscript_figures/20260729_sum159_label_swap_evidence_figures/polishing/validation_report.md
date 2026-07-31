# Polishing validation report

Status: **PASS** on 2026-07-29.

- The package contract requires panel-level provenance, regenerated subpanels,
  layout optimization/QC, and a single invoked R entrypoint.
- All 14 panel-map entries have locally generated audit subpanels, dimension
  records, optimizer placements, and provenance rows.
- All four expected final PNGs are readable and were assembled from live R
  objects.
- Figure 4d was regenerated from raw TIFF channels and masks; no approved,
  audit, intermediate, or final PNG was used as a composite input.
- Two final-phase replays produced identical SHA-256 hashes for every output.
- Independent visual QC passed after a clipped Figure 4d column header was
  repaired in the source constructor and the full chain was rerun.
- The machine-readable postflight report is `validation_report.json`; its
  overall status is `OK`. The only validator warning is the expected use of
  raw `grid` viewports. `notes.md` records that optimizer coordinates are
  consumed directly as lower-left origins without y inversion.
