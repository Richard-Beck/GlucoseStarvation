# Panel layout optimization report

- Input dimensions: `manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/layout/subpanel_dimensions.csv`
- Layout plan: `manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/layout/layout_plan.csv`
- Target width: 7.00 in
- Maximum height: 9.25 in
- Gap: 0.08 in
- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.

## Figures

- `Candidate supplement A`: 6.99 x 5.53 in, wasted 2.0%, tree `[[a | b] / c]`
- `Candidate supplement B`: 7.00 x 6.13 in, wasted 1.3%, tree `[a / b]`

## Scale recommendations

Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly.
