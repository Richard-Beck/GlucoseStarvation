# Resegmentation figure rebuild notes

Figures 1, 2 and S1-S4 were regenerated from the 20260722 canonical all-line Stan data, resegmentation area summaries, model-free feature release, and rerun classifier-validation release.
The accepted v02 plotting constructors, themes, facet orientations, panel dimensions, and Figure 1 manual layout were retained. Only their data-loading and derivation layer was replaced; no old report_exports tables are runtime inputs.
The count adapter uses N_obs directly as live cells and defines total cells as N_obs + D_obs. Decimal G0 labels retain their leading zero (0.1, 0.25, and 0.5 mM).
Cell and nuclear area intervals are empirical replicate-level interquartile ranges. Glucose calibration panels are derived directly from the fixed calibration fields in the canonical Stan data.
Visual inspection was skipped by explicit user instruction; deterministic file, schema, layout, and provenance validation was still performed.
Project-map decision: no update to docs/project_map.txt; this timestamped manuscript figure package does not change maintained project organization.
