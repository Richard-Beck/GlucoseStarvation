# Legend validation report

**Status: PASS — USER APPROVED (2026-07-31)**

## Scope

- Integration root: `agent-dev/manuscript_integration/20260729_v03_figure_set_integration`
- Expected rendered figures: 23
- Present legend blocks: 23
- Expected and explicitly described panels: 79
- Main-then-supplement order: Figure 1-Figure 5, then Figure S1-Figure S18
- Figure S19: absent by user decision
- Renderer-facing output: `integrated_figure_legends.md`

## Sources inspected

- `figure_set_manifest.csv`, `semantic_routing.tsv`, and the 79 routed panel semantic descriptions
- The 23 accepted final PNGs and package-local panel maps/notes for the replaced figure packages
- The prior v02 integrated legend set as a recycling baseline
- The current reviewer-cleared Methods text for definitions and statistical essentials
- The current Results preview as a non-authoritative consistency guardrail
- The SUM-159 polishing-package legend and source descriptions for Figures S15-S18

## Deterministic checks

| Check | Result |
|---|---|
| One `## Figure ...` block per final PNG | PASS (23/23) |
| Manifest order matches legend order | PASS |
| Every manifest panel letter is explicitly described in its figure block | PASS (79/79) |
| Nonempty title and body for every block | PASS |
| No extra figure labels or Figure S19 | PASS |
| No source paths, commands, hashes, provenance notes, or workflow prose in legend bodies | PASS |
| No stale 448-schedule/legacy Figure 5 text | PASS |
| No stale learned-environment or leave-one-line Figure S11-S14 text | PASS |
| Current S2 validation values present | PASS |
| Current Figure S7 cases and scores (+494, -802) present | PASS |
| Current Figure S10 diagnostic exception present | PASS |

## Reconciliation notes

- Figures 1-4 and S1-S10 retain prior wording where the current panels still support it, with targeted factual corrections.
- Figure 5 and Figures S11-S14 were rebuilt to match the replacement posterior-strategy, single-line posterior, fit-context, and morphology panels.
- Figures S15-S18 were revised to distinguish operational class labels from biological interpretation and to preserve the unresolved cytoplasmic-red boundary.
- Figure S6 and Figure S9 now correctly describe posterior-median trajectories rather than selected optimization draws.
- Figure S13 and Figure S14 are identified as exploratory maximum-a-posteriori comparisons without posterior or model-selection uncertainty.

## Accepted caveats and unresolved decisions

- Figure 1a and Figure S18d have no calibrated physical scale; Figure S18d reports field dimensions in pixels.
- The SUM-159 lower/2N and higher/4N labels are operational and retain the project-owner-accepted provenance limitation.
- The cytoplasmic-red signal is not assigned a biological reporter interpretation.
- No blocking legend decisions remain. Final journal style, title preference, and copyediting remain available for user review or manuscript rendering.

## User disposition

On 2026-07-31, the project owner declined further legend review and approved
the legend package as presented. It is ready for downstream manuscript use.
