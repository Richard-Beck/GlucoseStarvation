# Results revision validation report

## Status

PASS — ready for project-owner review.

The package contains five complete renderer-ready Results sidecars and a
manuscript-facing combined preview. Legend writing completed separately and is
user-approved; manuscript rendering and assembly have not started.

## Inputs

- Prior language baseline: `manuscripts/20260713T163237_v02/source/manuscript_sections/results/`
- Current figure evidence: `agent-dev/manuscript_integration/20260729_v03_figure_set_integration/`
- Current accepted claim graph: `agent-dev/manuscript_claim_graph/20260730_v03_clean_figure_level_graph/`
- Style authority: `.agents/skills/results-text/references/results_style_exemplars.md`

The v02 Results were used only as a language and organizational baseline.
Current semantic figure interpretations were treated as the primary scientific
evidence.

## Delegation and reconciliation

Five disjoint `gpt-5.6-sol` xhigh writer threads revised one existing Results
section each. A later dedicated Results-writer thread revised only the first
section in response to three inline project-owner comments. The main agent then
checked the complete five-section return for contract conformance, current
figure grounding, claim coverage, obsolete content, cross-section continuity,
and manuscript-facing hygiene, and made one narrow wording repair to avoid an
unsupported comparison of glucose-decline magnitude.

## Contract checks

- Expected sidecars: 5; present: 5.
- Every sidecar contains YAML routing metadata, a nonempty conclusion-first
  title, `## Results Text`, and `## Revision Notes`.
- Every assigned main and supplemental figure resolves in the current manifest.
- The combined Results body contains 2,893 words in 23 paragraphs.
- The combined preview contains only the manuscript-facing Results bodies;
  revision notes remain in the individual sidecars.

## Figure coverage

All 23 current rendered figures are cited at least once across the five-section
Results package: Figures 1-5 and Figures S1-S18. No retired Figure S19 or
non-current figure identifier is cited.

Panel citations were retained only when they support a specific statement.
Exhaustive subpanel-by-subpanel narration was not required.

## Authoritative-claim coverage

| Claim | Results handling | Assessment |
|---|---|---|
| C0: Ploidy modulates response to glucose starvation. | Directly synthesized in the measurement, direct-feature, posterior, and selection sections. | Pass |
| C1: High ploidy increases glucose consumption rate. | Stated directly for four-line culture-level glucose drawdown per live-cell AUC and separately as higher model-inferred apparent uptake. | Pass with current scope boundary |
| C2: High ploidy decreases yield (per cell). | Direct features retain the distinction from unnormalized total-cell output; posterior fits state lower apparent per-cell yield across current F4 comparisons. | Pass as model-supported |
| C3: High ploidy cells are more robust to starvation. | Stated directly from zero-glucose live-cell persistence; raw-trajectory and modeled countercontexts remain local. | Pass with current context boundary |
| C4: Starvation selects for high ploidy. | Stated directly as a fitted-simulation prediction, not empirical competition evidence. | Pass as model-supported |
| C7: Cell size is proportional to ploidy. | The general positive association is retained while repeated SUM-159-fuse and other size reversals are stated explicitly. | Pass at qualified strength |
| C11: Abundance selects for low ploidy. | Stated directly as a fitted-simulation prediction with schedule, line, model, and fit-context exceptions. | Pass as model-supported |
| C12: High ploidy increases yield (volumetric). | Retained explicitly as a partially supported proxy-based working hypothesis rather than a measured result. | Pass as working hypothesis |

## Recycling outcome

The prior Results contained 25 paragraphs and 3,363 words. The current Results
contain 23 paragraphs and 2,893 words.

- 2 prior paragraphs were retained.
- 11 were lightly revised.
- 7 were substantially revised.
- 5 were removed or replaced.
- 2 focused new paragraphs were added.

Thus, 13 of 25 prior paragraphs were retained or lightly revised. The largest
necessary changes were removal of the obsolete Figure S14 environment-selection
transfer branch, replacement of old S13/S14 interpretations, correction of the
current F4 parameter-direction account, reconstruction of the selection section
around current F5/S11/S12, and the project-owner-directed tightening of the
SUM-159 identity concern and raw-trajectory synthesis in the first section.

## Scientific and hygiene checks

- Targeted numerical statements were checked against the current semantic
  interpretations, including measurement performance, F3/S7 model examples,
  F4 distances and parameter directions, S10 diagnostics, S13 morphology-AIC
  comparisons, and current F5/S12 support displays.
- Old numerical claims not authenticated by the current semantic package were
  removed rather than inherited.
- The obsolete 9-of-10/10-of-10, model-native/common-R1, target-inclusive, and
  35-of-50 transfer discussion is absent from the manuscript-facing text.
- Current Results bodies contain no local paths, workflow labels, claim-graph
  labels, provenance language, or generation bookkeeping.
- Model-implied R2/W1 states are not presented as measured biological entities.
- Culture-level glucose drawdown is not presented as direct single-cell uptake
  or metabolic flux.
- Fitted selection predictions are not presented as empirical competition.
- The biological meaning of the Figure S18 cytoplasmic-red measurement remains
  explicitly unresolved.
- All three inline `USER_COMMENT` annotations were addressed in the first
  sidecar and are absent from the manuscript-facing combined preview.

## Remaining review surface

No independent whole-Results reviewer was required or used. Project-owner
review should therefore focus on scientific emphasis and desired narrative
density, particularly the new SUM-159 paragraph, the revised current F4
parameter interpretation, and the reconstructed Figure 5 selection section.

## Recommendation

Use `combined_results_preview.md` for narrative review and the individual
sidecars for section-specific revision notes. After project-owner acceptance,
the current Results can be included in the next A/I/D context bundle; manuscript
rendering and assembly should remain downstream of that handoff.
