# Panel layout optimization report

- Input dimensions: `manuscript_figures/20260722T141236_resegmentation_fix/layout/subpanel_dimensions.csv`
- Layout plan: `manuscript_figures/20260722T141236_resegmentation_fix/layout/layout_plan.csv`
- Target width: 7.00 in
- Maximum height: 9.25 in
- Gap: 0.08 in
- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.

## Figures

- `Figure 1`: 6.97 x 9.09 in, wasted 2.9%, tree `[a / [b | [[c | d] / [e | f]]]]`
- `Figure S1`: 7.00 x 6.85 in, wasted 0.0%, tree `a`
- `Figure S2`: 6.71 x 8.47 in, wasted 3.8%, tree `[[[a / [b / c]] | [[d / e] / f]] / g]`
- `Figure 2`: 7.00 x 7.81 in, wasted 2.0%, tree `[a / [b / c]]`
- `Figure S3`: 7.00 x 4.25 in, wasted 0.0%, tree `a`
- `Figure S4`: 7.00 x 8.20 in, wasted 0.0%, tree `a`

## Scale recommendations

Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly.
