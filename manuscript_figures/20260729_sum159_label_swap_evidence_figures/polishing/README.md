# SUM-159 label-swap evidence — permanent polishing package

The four figures in `final_images/` are the authoritative package outputs:

1. `figure_1_all_timecourses.png` — complete segmented-cell and nuclear-area
   time courses.
2. `figure_2_confluence_robustness.png` — robustness to field cell count and
   total segmented area.
3. `figure_3_focused_distributions_and_same_2n.png` — focused distributions and
   direct low/2N comparisons.
4. `figure_4_cytoplasmic_signal_and_multimodal_fields.png` — cytoplasmic-red
   contrasts and raw-channel multimodal fields.

`figure_story.md` contains the panel-referenced interpretation, `legends.md`
defines the panels, and `provenance.csv` maps every displayed panel to its
local generator, inputs, audit subpanel, layout plan, and final image.

## Reproduction

Run the full package-local chain from the project root:

```bash
scripts/agentRrunner.sh manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/scripts/polish_figures.R --phase subpanels
scripts/agentRrunner.sh .agents/skills/manuscript-figure-workflow/scripts/optimize_panel_layout.R --input manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/layout/subpanel_dimensions.csv --output-dir manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/layout --target-width 7 --max-height 9.25 --gap 0.08 --sx-lo 0.9 --sx-hi 1.05 --sy-lo 0.9 --sy-hi 1.05
scripts/agentRrunner.sh manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/scripts/polish_figures.R --phase final
```

The final phase reconstructs every panel as a live R object and assembles the
figures in memory. It does not import the audit subpanel PNGs or the visual-only
files in `approved_reference_images/`.
