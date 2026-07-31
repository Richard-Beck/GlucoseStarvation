# Wave 1 launch validation

- Release preparation: PASS (107 checks).
- Fit declaration: PASS (47 datasets, 12 models, 564 fits, 282,000 starts; zero NUTS tasks).
- Wave 1 selection: PASS (11 datasets, 132 fits, 66,000 starts).
- Smoke: PASS. Jobs `19194594` and `19194596` produced 10 complete start triplets; the reducer combined them successfully.
- Dispatch: PASS. Job root `data/modelling/gpath_v1/red_a30_counts_20260722/jobs/optim_20260722_101722/` records 132 optimization arrays and 132 timeout-recovery coordinators.
- Initial execution check: PASS. Per-start outputs and successful task logs appeared immediately; no shell-level failure markers were found.
- Completion: PASS. All 66,000 start triplets and 132 combined-fit triplets are present; seven timeout retries completed.
- Retroactive repair: PASS. All 47 Stan-data RDS files and 396 combined optimization RDS files remained byte-identical after relocation.
- Canonical release validation: PASS (108 prepare checks). Active manifests use the canonical release root and direct canonical count source; only immutable historical SLURM logs retain the original launch path.
- Optimization wave 2 configuration: PASS. The 36-dataset filter selects exactly 432 fit families and 216,000 starts, all with 8 CPUs/8 GB/30 minutes and `xxlarge` QOS; no pre-existing wave-2 shards were found.
- Optimization wave 2 dispatch: PASS. The job root `jobs/optim_20260722_145619/` records 432 arrays and 432 timeout-recovery coordinators; initial per-start outputs appeared successfully.
- NUTS configuration: PASS. The selected five-model front is crossed with seven dataset contexts, producing exactly 35 fits and 140 chain tasks; all optimizer initializers are present.
- NUTS release validation: PASS (110 prepare checks). All chain configs parse and carry the approved sampling and `xxlarge` resource settings.
- NUTS dispatch: PASS. Job arrays `19265005` through `19265039` were submitted under `jobs/nuts_20260722_144933/`; initial inspection found 99 running and 41 pending tasks.
