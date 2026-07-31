# Selection and drafting record

## Decision

Promote two supplementary candidates, without manuscript numbering:

1. **Candidate supplement A** — covariate definition, likelihood-AIC comparison, and localization of fit changes.
2. **Candidate supplement B** — the fitted per-unit covariate effects and their realized elevated-versus-baseline consequences for v1, K1, y1, and m.

This is preferred over a main-plus-supplement split. The nuclear-area result is biologically suggestive, but it is an exploratory MAP comparison without parameter uncertainty or validation analysis, and the improvement does not motivate changing the primary ploidy encoding used throughout the manuscript.

## Disposition of rough figures

- Retain and revise the covariate-value display, Pareto-model AIC comparison, condition-level likelihood decomposition, and realized-parameter effects.
- Replace the 12-model raw-id AIC display with the five manuscript-retained Pareto aliases.
- Replace the separate mean-condition and SUM-159-fuse heatmaps with one integrated glucose-by-cell-line display faceted into five model columns.
- Promote a clarified per-unit coefficient panel with explicit model rows, because the fitted coefficient is global rather than cell-line-specific.
- Retain the realized elevated-versus-baseline parameter panel, showing the five model-level MAP points plus their median and range for each cell line.
- Do not promote the absolute AIC-versus-complexity plot or the within-family AIC heatmap.
- Do not promote the 15 full trajectory atlases. They remain useful fit audits but are too repetitive for the incremental claim and do not localize the likelihood gain as directly as the condition decomposition.

## Interpretation boundary

Likelihood AIC is computed from observation log likelihood at the best MAP solution with the same nominal parameter count within each model family. Negative delta AIC favors the area metric. These fixed-case MAP comparisons do not include posterior uncertainty or model-selection uncertainty. Condition-level log-likelihood changes are descriptive decompositions of the same fits and should not be interpreted as independent tests.

Per-unit coefficients and realized pair effects answer different questions. The former is the global fitted log2 parameter change for a +1 shift in the chosen covariate and is therefore displayed by model, not cell line. The latter multiplies that fitted relationship by each line's observed elevated-versus-baseline covariate contrast and is therefore displayed by cell line.
