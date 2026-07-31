# Panel layout optimization report

- Input dimensions: `manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/layout/subpanel_dimensions.csv`
- Layout plan: `manuscript_figures/20260729_sum159_label_swap_evidence_figures/polishing/layout/layout_plan.csv`
- Target width: 7.00 in
- Maximum height: 9.25 in
- Gap: 0.08 in
- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.

## Figures

- `SUM159 evidence 1`: 7.00 x 9.18 in, wasted 0.9%, tree `[a / b]`
- `SUM159 evidence 2`: 7.00 x 9.22 in, wasted 2.0%, tree `[[a | b] / [c | d]]`
- `SUM159 evidence 3`: 6.92 x 7.98 in, wasted 2.1%, tree `[[a | b] / [c | d]]`
- `SUM159 evidence 4`: 6.97 x 9.15 in, wasted 1.4%, tree `[[[a | b] | c] / d]`

## Scale recommendations

Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly.
