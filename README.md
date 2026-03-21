# GlucoseStarvation

This project analyzes how changes in ploidy alter cancer cell behavior under glucose limitation.

The scientific aim is not predictive transfer as an end in itself. The central question is whether ploidy-associated effects on growth, death, glucose use, and yield are generalizable enough across cell lines to support strategies that would selectively favor or disfavor high-ploidy states, including in cell lines not directly observed during model fitting.

The repo approaches that question in three connected ways:

- a hierarchical mechanistic `gpath` modeling workflow built from live/dead count trajectories plus glucose measurements,
- a reduced model-free summary workflow that extracts interpretable ploidy-response features from the same processed dataset,
- an exploratory morphology side branch based on image-derived area quantiles.

Directional transfer workflows are included because they stress-test whether inferred ploidy effects are stable across cell lines. In this project, those transfer benchmarks are used as evidence about cross-line generalizability for intervention or selection ideas, not as the primary scientific objective.

For project organization and the current active analysis path, start with [docs/project_map.txt](docs/project_map.txt).
