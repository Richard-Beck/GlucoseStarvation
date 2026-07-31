# v03 figure-set integration report

## Outcome

PASS. This root is the canonical integration of Figures 1–5 and Figures
S1–S18. It contains 23 normalized whole-figure PNGs and 79 manuscript-visible
panel records. Figure S19 is absent by explicit user decision.

All 23 accepted polishing PNGs were regenerated from promoted
integration-local code and fixed source inputs. Every polishing-final,
integration-rebuild, wrapper-rebuild, and normalized-final SHA256 comparison
matches.

## Integrated generation families

| Family | Canonical figures |
|---|---|
| resegmentation core | F1, F2, S1–S4 |
| mechanistic diagnostics | F3, S5–S9 |
| posterior size context | F4, S10, S11 |
| posterior strategy main | F5 |
| posterior strategy context | S12 |
| morphology metrics | S13, S14 |
| SUM-159 label-swap evidence | S15–S18 |

The review-manifest generator entry for S10 was corrected during integration:
its actual generator is the posterior-size-context entrypoint, not the
mechanistic-diagnostics entrypoint.

Canonical panel routing also expands under-declared source-package provenance
for the SUM figures: S16–S17 now name the field-metadata and mixture-call inputs
actually consumed by their plotting code, and S18d names the seven selected raw
TIFF/mask input sets behind its multimodal grid. The later permanent-package
repair promoted the active SUM inputs and generators into that polishing
package; its accepted PNG content and hashes remain unchanged.

## Reproducibility and raster policy

`final_figure_scripts/run_all_figures.sh` reruns all seven promoted generation
families through `scripts/agentRrunner.sh`, writes their outputs beneath
`package_rebuilds/`, verifies them against the accepted polishing outputs, and
normalizes only matching integration-local rebuilds into `final_images/`.

No accepted final or subpanel raster was used as a shortcut for local
regeneration. Figure 3a retains the explicit user-approved immutable schematic
raster exception. The SUM-159 microscopy panels are reconstructed locally from
their declared source image channels, masks, and live plotting/assembly code.

The manuscript-level wrapper completed two clean integration replays during
this run: the initial promoted-code verification and the final release
validation replay.

## Semantic and visual contracts

Fresh semantic interpretations were produced for all 79 panels. In accordance
with the user's supervision rule, these were the only delegated tasks and used
`gpt-5.6-terra` at medium reasoning. The main agent checked required section
completeness and linkage but did not reinterpret panel semantics.

Whole-set visual QC is inherited from the user's explicit acceptance of the
review set because all 23 integrated finals are byte-identical to those
reviewed PNGs. The structured record is `whole_set_visual_qc.tsv`.

## Validation

`figure_set_validation_report.json` records:

- status: `PASS`
- errors: 0
- final images: 23
- canonical panel rows: 79
- semantic interpretation artifacts: 79
- rebuild records: 46
- whole-image visual-QC records: 23
- S19: absent

The two recorded warnings are intentional scope notes, not unresolved
integration defects: package subagents were prohibited by the user, and a full
prior-subpanel lineage audit was not requested. Minimal prior correspondence
fields are retained in the canonical figure index where recoverable.

## Canonical outputs

- `final_images/`
- `figure_set_manifest.csv`
- `figure_rebuild_manifest.tsv`
- `figure_byte_identity_report.tsv`
- `semantic_interpretations/`
- `whole_set_visual_qc.tsv`
- `figure_set_validation_report.json`
- `package_inventories/`
- `omitted_packages.md`
