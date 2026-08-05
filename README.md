# GlucoseStarvation

This project analyzes how changes in ploidy alter cancer cell behavior under glucose limitation.

The scientific aim is not predictive transfer as an end in itself. The central question is whether ploidy-associated effects on growth, death, glucose use, and yield are generalizable enough across cell lines to support strategies that would selectively favor or disfavor high-ploidy states, including in cell lines not directly observed during model fitting.

The repo approaches that question in three connected ways:

- a hierarchical mechanistic `gpath` modeling workflow built from live/dead count trajectories plus glucose measurements,
- a reduced model-free summary workflow that extracts interpretable ploidy-response features from the same processed dataset,
- an exploratory morphology side branch based on image-derived area quantiles.

Directional transfer workflows are included because they stress-test whether inferred ploidy effects are stable across cell lines. In this project, those transfer benchmarks are used as evidence about cross-line generalizability for intervention or selection ideas, not as the primary scientific objective.

## Current manuscript audit package

The current public manuscript and audit artifacts are available at
[richard-beck.github.io/GlucoseStarvation](https://richard-beck.github.io/GlucoseStarvation/).

The package includes:

- illustrated and text-only manuscript views;
- panel-level figure semantics;
- the integrated claim graph;
- locked Methods provenance;
- the structured literature map;
- an artifact manifest and release identity.

The published package is a rapidly iterated, unsealed review artifact distributed
separately from the Git repository because the rendered manuscript bundle is
large. Publication does not imply that a complete isolated rebuild or manuscript
validation gate has passed.

Associated repository commit:
[`3643d579f771b64bb10e00cbed6e790f5e23f380`](https://github.com/Richard-Beck/GlucoseStarvation/commit/3643d579f771b64bb10e00cbed6e790f5e23f380).
The exact published payload is distributed through the
[`pages-site` release](https://github.com/Richard-Beck/GlucoseStarvation/releases/tag/pages-site).

The proposed long-term target-workspace policy is recorded in
[`docs/target_policy_DRAFT.md`](docs/target_policy_DRAFT.md); it is not yet an
active policy.

Superseded, exploratory, and transient material is retained locally under the
Git-ignored `junk/` tree with its original repository-relative layout. Active
code must not depend on `junk/`.
