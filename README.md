# GlucoseStarvation

This project analyzes how changes in ploidy alter cancer cell behavior under glucose limitation.

The scientific aim is not predictive transfer as an end in itself. The central question is whether ploidy-associated effects on growth, death, glucose use, and yield are generalizable enough across cell lines to support strategies that would selectively favor or disfavor high-ploidy states, including in cell lines not directly observed during model fitting.

The repo approaches that question in three connected ways:

- a hierarchical mechanistic `gpath` modeling workflow built from live/dead count trajectories plus glucose measurements,
- a reduced model-free summary workflow that extracts interpretable ploidy-response features from the same processed dataset,
- an exploratory morphology side branch based on image-derived area quantiles.

Directional transfer workflows are included because they stress-test whether inferred ploidy effects are stable across cell lines. In this project, those transfer benchmarks are used as evidence about cross-line generalizability for intervention or selection ideas, not as the primary scientific objective.

The current reproducibility target is
[`manuscripts/20260731T101031_v03`](manuscripts/20260731T101031_v03). For its
rebuild boundary and validation commands, see
[`docs/current_manuscript_rebuild.md`](docs/current_manuscript_rebuild.md).
For compact project navigation, see
[`docs/project_map.txt`](docs/project_map.txt).

Superseded, exploratory, and transient material is retained locally under the
Git-ignored `junk/` tree with its original repository-relative layout. Active
code must not depend on `junk/`.
