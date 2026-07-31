# Panel layout optimization report

- Input dimensions: `manuscript_figures/20260728_posterior_strategy_figure_set/F5/polishing/layout/subpanel_dimensions.csv`
- Layout plan: `manuscript_figures/20260728_posterior_strategy_figure_set/F5/polishing/layout/layout_plan.csv`
- Target width: 7.00 in
- Maximum height: 9.25 in
- Gap: 0.00 in
- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.

## Figures

- `Figure 5`: 7.00 x 9.25 in, wasted 0.0%, tree `[[a / [b / [c / d]]] | [e / [f / g]]]`

## Scale recommendations

Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly.
