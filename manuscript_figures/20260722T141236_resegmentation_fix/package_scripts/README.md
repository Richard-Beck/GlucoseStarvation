# Figure-package scripts

`scripts/polish_figures.R` remains the completed FG1/FG2 rebuild entrypoint.
The four remaining packages have one self-contained entrypoint each:

- `FG3_mechanistic_diagnostics/polish_figures.R`
- `FG4_posterior_size_context/polish_figures.R`
- `FG6_selection_simulation_revision/polish_figures.R`
- `FG7_prediction_transfer/polish_figures.R`

Each entrypoint reads canonical project data, performs derived transformations
in memory, and writes only its whole-figure PNGs to `../../final_images/` and
its operation timings to `../../timings/`. Persisted transformed tables,
subpanel rasters, layout previews, caches, and other intermediate artifacts are
not permitted. Session-temporary files must use `tempdir()` and are not package
outputs.
