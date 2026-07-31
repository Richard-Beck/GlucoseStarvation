# Clean claim-graph integration report

## Status

Reconciliation complete. Deterministic audit/graph/package validation and
rendering are performed after this build by `validate_clean_return.py`.
The completed production return passes that validator.

## Inputs

- Figure integration: `agent-dev/manuscript_integration/20260729_v03_figure_set_integration`
- Current figures: 23 (F1–F5 and S1–S18)
- Semantic panel descriptions: 79
- Authoritative claims: 8
- Audit prompt: `figure_audit_prompt_contract.md`
- Reuse mode: `fresh_full`

## Audit execution

One fresh-context GPT-5.6-sol medium auditor assessed each figure. Auditors
received only one figure's semantic descriptions, the eight authoritative
claims, and the audit prompt contract. Execution used at most eight simultaneous
auditors. The 23 normalized audit tables contain 103 candidate
relationships and pass the bundled figure-audit validator.

Formatting normalization changed only Markdown/table contract details and exact
claim-cell expansion; semantic observations, relations, and strengths were not
altered during audit materialization.

## Reconciliation

The main agent accepted 81 worker rows into
40 grouped evidence nodes and rejected
22 rows with explicit reasons. Five lean umbrella claims
replace panel-level claim proliferation. The resulting graph has
13 claims and 56 relationships. No prior
non-authoritative structure and no qualification edge were imported.

## Boundary

This return does not update manuscript prose, Results, Methods, rendering, or
any downstream consumer of the old graph contract. Those consumer updates
remain explicitly deferred.
