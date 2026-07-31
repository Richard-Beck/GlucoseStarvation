# Polishing and raster-compliance notes

This permanent package supersedes the temporary development package at
`tmp/20260729_sum159_label_swap_evidence_figures/`. The temporary root remains
untouched as historical development provenance.

The package also owns the previously external SUM-159 analysis inputs required
by its active figure constructor. The GFP-defined competition calls,
cytoplasmic-red measurements and sampling manifest, and full-scope platemap
summary now live under `derived_data/`; their generators are retained under
`scripts/`. The package does not read either `tmp/` or the v02 manuscript
development tree during figure regeneration. The promoted numerical files are
byte-identical to the accepted inputs. Upstream classification and measurement
analyses were not rerun during this path-only provenance repair.

The original promotion was invalid because it allowed split R/Python assembly,
used intermediate rendered PNGs as figure inputs, and restored an accepted
whole-figure PNG after rebuilding. This repair removes those paths. The active
package has one invoked R entrypoint, `scripts/polish_figures.R`, with
`subpanels` and `final` phases. Its package-local sourced helper constructs the
quantitative plots and derived tables; the final phase rebuilds all live panel
objects and draws them directly into optimizer-defined viewports.
The optimizer's `x_npc` and `y_npc` coordinates are consumed as lower-left
origins without y-axis inversion.

All 14 displayed panels have fresh audit subpanels. Figures 1–3 and Figure
4a–c are regenerated from numeric measurements. Figure 4d is regenerated in R
from the original TIFF channels and segmentation masks identified by
`derived_data/multimodal_field_manifest.csv`: crop extraction, channel
normalization, composites, and class-colored object boundaries are recomputed
locally. Raw microscopy channels are data inputs; no whole microscopy
subpanel, final PNG, approved-reference PNG, or audit PNG is copied or reread
into a final composite. No raster exception is used by this package.

The jittered quantitative panels use fixed seeds. Two successive final-phase
runs produced byte-identical outputs for all four figures. The historical v3
PNGs are retained only in `approved_reference_images/` for visual comparison;
they are not build inputs and byte identity to them is neither expected nor
claimed.

The independent visual audit identified a clipped rightmost header in Figure
4d. The header was wrapped in the source constructor and the complete
subpanel/layout/final chain was rerun before postflight validation.

No manuscript numbering is fixed inside this package. The current review set
maps these four source-local figures to S15–S18.

Project-map decision: `docs/project_map.txt` is updated because the package's
active regeneration contract materially changed from copied accepted rasters
to complete local regeneration.
