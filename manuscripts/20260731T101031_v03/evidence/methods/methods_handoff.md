# Methods Handoff

## Status

Reviewer-cleared and ready for project-owner review.

The canonical manuscript prose is:

`agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction/methods_text.md`

SHA-256:

`b766a8019ad718462614fcb07f4c45633a7ee9a40334018fa68f275c998a0ad7`

The 7,493-word canonical file is an exact extraction of the final
writer-authored draft. The Overseer did not revise its prose.

## Canonical Methods package

- `target_figure_set.tsv`
- `manuscript_endpoint_verification.md`
- `methods_text.md`
- `locked_provenance_table.md`
- `provenance_lock_verification.md`
- `node_classification.csv`
- `methods_spine.md`
- `final_methods_draft_audit.md`
- `methods_handoff.md`

The verified provenance graph contains 271 objects and 417 edges and exactly
matches all 79 panel endpoints across the 23-figure set.

## Drafting and review record

- `drafting/writer_round_01.md`: first spine-based draft and first request set.
- `drafting/inspector_round_01.md`: ten current-source fact cards.
- `drafting/writer_round_02.md`: fact-card revision and review-ready draft.
- `drafting/reviewer_round_01.md`: 14 reviewer points, including v02
  omission-check points.
- `drafting/writer_round_03.md`: reviewer-response revision and second request
  set.
- `drafting/inspector_round_02.md`: eleven current-source fact cards.
- `drafting/writer_round_04.md`: final writer-authored draft and complete
  dispositions.
- `drafting/reviewer_round_02.md`: no remaining reviewer points.

Writers and reviewers used `gpt-5.6-sol` at xhigh reasoning. Each
source-inspection round used one `gpt-5.6-sol` medium inspector. The v02 Methods
text was visible only to reviewers and was used only to identify potentially
missing details.

## Material project-owner decisions

The final draft deliberately exposes, rather than repairs or conceals:

- incomplete wet-lab culture, imaging, assay, and replicate documentation;
- two likelihood contributions for each ploidy-unspecified glucose entry;
- retained nuclear-mask overshoots in the 35 rerun fields;
- unpinned classifier, CPSAM, preliminary Cellpose, and software identities;
- operational rather than biologically validated SUM-159 green-component
  labels;
- MAP selection without a convergence/return-code rejection rule;
- posterior diagnostics without encoded acceptance thresholds;
- posterior-prediction and intervention integration using LSODA defaults;
- absence of a formal intervention-schedule ranking analysis; and
- pending code/data availability identifiers and release metadata.

See `final_methods_draft_audit.md` and the final writer's `## Open Issues` for
the complete decision list.

## Consumer guidance

Use `methods_text.md` as the sole manuscript-prose input. Do not consume
writer packets, fact cards, provenance tables, or review reports as manuscript
text.

If project-owner review changes scientific procedures, analysis scope, or any
material limitation, revise through a new writer turn and return the modified
text to a Methods reviewer. Pure journal-formatting changes do not require
provenance reconstruction.
