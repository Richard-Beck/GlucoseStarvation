# Optimization release record

- Approval: user requested dispatch of all-line, five leave-one-line-out, and five single-line datasets on 2026-07-22.
- Input: canonical `run_20260721_163410/summaries/counts.Rds`, SHA256 `9d17ded4eb7787d02e2917a111eecf770fd8d05ff176b7a09c58ead2ffd35b4c`; arithmetic mean per imaged field including zero-object fields, referenced directly rather than copied.
- Wave 1 scope: 11 datasets x 12 assessment models x 500 starts = 66,000 starts.
- Resources: 8 CPUs/8 GB per start; 30 minutes for all-line/LOO, 15 minutes for single-line; `xxlarge` QOS. TIMEOUT shards retry once at 3x walltime.
- Outputs: `data/modelling/gpath_v1/red_a30_counts_20260722/`.
- Validation: rebuilt release passed 110 prepare-phase checks and contains all 47 declared Stan datasets.
- Wave 1 status: submitted on 2026-07-22 under `jobs/optim_20260722_154917/`.
- Optimization wave 2: 36 directional and matched-null datasets x 12 models x 500 starts = 216,000 starts; submitted under `jobs/optim_20260722_163144/`.
- NUTS scope: five all-line Pareto models across all-line, four-line no-SUM-159-fuse, and five single-line datasets (35 fits; four chains each).
- NUTS settings: 1,000 warmup, 1,500 sampling, `adapt_delta=0.995`, maximum treedepth 14, dense metric, 16 CPUs/32 GB/12 hours per chain, `xxlarge` QOS.
- NUTS status: all 35 fits/140 chains completed successfully under `jobs/nuts_20260722_163134/`; all expected draw and diagnostic files are present.
- Derived status: optimization assessment written and checked at `derived/optimization/assessment.Rds` (564 fits, 282,000 starts, 230 Pareto rows, 216 transfer/null rows). Posterior QC, seven parameter shards, and 35 prediction shards were submitted under `jobs/derived_20260723_141554/`.
