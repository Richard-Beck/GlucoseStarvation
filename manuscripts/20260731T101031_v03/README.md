# Glucose-starvation manuscript v03

Status: **WARN — coherent review assembly; not submission-complete**

This package assembles the project-owner-approved v03 figures, legends,
Results, and Methods with the ordered, verbatim-carried v02 Abstract,
Introduction, Discussion, and References. It contains 23 figures and 79
manuscript-visible panel endpoints; Figure S19 is intentionally omitted.

## Start here

- `draft/manuscript_draft.html` — self-contained manuscript review copy with
  all 23 PNGs embedded.
- `source/` — editable manuscript text, legends, title, and figure manifest.
- `review_state/review_summary.md` — approvals and current review status.
- `review_state/accepted_exceptions.md` — accepted and deferred limitations.
- `validation/assembly_validation_report.md` — final package verdict.
- `validation/draft_validation_report.md` — detailed render and source checks.
- `rebuild/README.md` — commands for rebuilding the draft and figures.

## Validation outcome

The assembly-local replay regenerated all 23 figures from the seven accepted
polishing families and matched every approved PNG checksum. The HTML contains
exactly 23 embedded figures, 23 matched legends, five injected Results
sidecars, canonical Methods, served A/I/D and References, unique anchors,
resolved internal links, and all eight user-fixed claims.

The `WARN` status reflects approved, manuscript-visible limitations rather than
an assembly failure: wet-lab protocol and biological-replication records,
several software/numerical identities, and archival release metadata remain
incomplete. See `evidence/methods/final_methods_draft_audit.md`.
