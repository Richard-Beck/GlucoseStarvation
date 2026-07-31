# Polishing notes

## Recommended placement

Both outputs are recommended as new supplementary figures. A main-plus-supplement split was considered and rejected: the nuclear-area encoding improves the optimized all-line observation likelihood consistently, but this is an exploratory MAP comparison without posterior uncertainty, holdout assessment, or a preregistered model-selection decision. The result is therefore better used to refine biological interpretation than to replace the manuscript-wide ploidy analysis.

## Figure content

- **Candidate supplement A** establishes the covariates and the fit comparison. Panel a shows the elevated member of each experimental pair; all within-line baseline values are zero. Panel b compares likelihood AIC against the existing ploidy fit for the five manuscript-retained Pareto models. Panel c integrates the nuclear-area-minus-ploidy observation log-likelihood difference across starting glucose, cell line, and five model columns.
- **Candidate supplement B** separates the fitted coefficient from its realized consequence. Panel a shows the global fitted log2 parameter change per +1 unit of ploidy or nuclear-area covariate for each Pareto model; these coefficients do not vary by cell line. Panel b applies the relevant line-specific elevated-versus-baseline covariate contrast: faint points are the five model-level MAP effects and the emphasized point and interval are their median and range.

## Scientific interpretation

The nuclear-area encoding is favored over ploidy in each of the five retained model structures, whereas the cell-area encoding is disfavored. Nuclear area therefore provides intriguing support for a size-linked explanation of the shared ploidy effect. The effect is not uniform across conditions: the integrated likelihood decomposition includes both gains and losses, including within SUM-159-fuse. The parameter-effect directions under nuclear area remain consistent with those under ploidy across the four displayed parameters and five retained models. Per-unit coefficient magnitudes should not be compared as though the two covariates had identical scales; the realized panel provides that line-specific comparison. Larger nuclear-area contrasts in several lines yield larger realized effects, whereas the near-zero SUM-159-fuse contrast attenuates its inferred effect. This combination supports retaining ploidy as the primary manuscript covariate while presenting nuclear size as a plausible mechanistic refinement.

Likelihood AIC is calculated from the observation log likelihood at the best MAP solution with the same nominal parameter count within a model family. Negative delta AIC favors an area encoding. This is conventional likelihood-AIC arithmetic applied to fixed MAP solutions, not a posterior model probability or an uncertainty-aware model comparison. The condition decompositions are descriptive parts of the same fitted likelihood and are not independent statistical tests.

## Data and provenance

The area covariates were derived exclusively from alive objects in 24--48 hour images from the corrected classification run:

`data/image_processing_runs/full_segmentation_classification_nuclear/run_20260721_163410/`

The four package-local source tables are immutable copies of the rough-plot tables derived from the resulting all-five-line cell-area and nuclear-area optimization release and the existing ploidy optimization release. Their generating scripts and upstream paths are recorded in `source_tables/README.md`. No new fits, segmentation, or image classification were run during polishing.

The maintained generators now obtain the five accepted model identifiers, aliases, and order from the current `red_a30_counts_20260722` all-lines optimization assessment. Existing source tables and figures were not regenerated for this selector-only correction because those consumed fields are exactly identical to the previously used selector; the equivalence and regeneration waiver are recorded in the current Methods provenance package.

The source-table generators are now retained inside this permanent package as
`scripts/make_condition_source_tables.R` and
`scripts/make_parameter_effect_source_tables.R`; neither the figure constructor
nor its provenance contract depends on `tmp/`. The historical rough PNGs are
not inputs to the polished figures.

## Rebuild

Run, in order:

1. `scripts/agentRrunner.sh manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/scripts/polish_figures.R --phase subpanels`
2. `scripts/agentRrunner.sh .agents/skills/manuscript-figure-workflow/scripts/optimize_panel_layout.R --input manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/layout/subpanel_dimensions.csv --output-dir manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/layout --target-width 7 --max-height 9.25`
3. `scripts/agentRrunner.sh manuscript_figures/20260722T141236_resegmentation_fix/morphology_metric_exploration/polishing/scripts/polish_figures.R --phase final`

Project-map decision: `docs/project_map.txt` was updated because this is a new promotion-ready figure package adjacent to the current canonical manuscript rebuild. The candidates are not yet integrated into or numbered within the canonical figure set.
