# Parameter estimation

One release config builds immutable, dataset-first model inputs and fit plans under:

`data/modelling/<model>_<version>/<release_id>/datasets/<dataset>/{optim,nuts}/<model_id>/<fit_id>/`

The preparation chain is:

1. count input -> validated direct reference or release-local `uncorrected.Rds`;
2. corrected `R/prepare_data.R` -> `all_lines/stan_data.Rds`;
3. declared subsets/holdouts -> one `stan_data.Rds` per dataset;
4. fit specification -> optimization and NUTS manifests.

The optional `morphology_metrics` release block enables
`morphology_cell_area` and `morphology_nuclear_area` dataset transforms. These
replace the model's baseline-relative `ploidy_metric` covariate with a
baseline-relative log2 morphology metric while preserving `ploidy_abs` as the
experimental pair identity. Metric values and source provenance are written
under the release-local `inputs/morphology/` directory.

Optimization starts use CmdStan's default random initialization radius of 2
unless `GPATH_OPTIM_INIT_RADIUS` is set to a smaller positive value. The
effective radius is recorded in each optimization run manifest.

`GPATH_OPTIM_INIT_SOURCE_ROOT` may point to an authenticated all-line
optimization tree. In that mode each task initializes from the correspondingly
ranked completed start for its model ID, and records the resolved source
directory in its run manifest. This is useful for matched covariate-only refits.
`GPATH_OPTIM_INIT_OFFSET` selects a later ranked source start while retaining
the output task ID; `repair_missing_optim.sh` uses this for targeted replacement
of missing non-timeout shards before recombination.

Run from the project root:

```bash
scripts/parameter_estimation/run.sh prepare CONFIG.json
scripts/parameter_estimation/run.sh submit-optim CONFIG.json --dry-run
scripts/parameter_estimation/run.sh submit-optim CONFIG.json
scripts/parameter_estimation/run.sh submit-optim CONFIG.json --datasets DATASET_IDS.tsv
scripts/parameter_estimation/run.sh submit-nuts CONFIG.json --dry-run
scripts/parameter_estimation/run.sh submit-nuts CONFIG.json
scripts/parameter_estimation/run.sh submit-nuts CONFIG.json --defer-after-optim JOB_IDS.tsv
scripts/parameter_estimation/run.sh status CONFIG.json
scripts/parameter_estimation/run.sh validate CONFIG.json complete
scripts/parameter_estimation/run.sh derive-plan CONFIG.json
scripts/parameter_estimation/run.sh derive-optim CONFIG.json
scripts/parameter_estimation/run.sh submit-derived CONFIG.json --dry-run
scripts/parameter_estimation/run.sh submit-derived CONFIG.json
scripts/parameter_estimation/run.sh derived-status CONFIG.json
scripts/parameter_estimation/run.sh validate-derived CONFIG.json
scripts/parameter_estimation/run.sh strategy-plan CONFIG.json
scripts/parameter_estimation/run.sh submit-strategies CONFIG.json
scripts/parameter_estimation/run.sh strategy-status CONFIG.json
scripts/parameter_estimation/run.sh validate-strategies CONFIG.json
scripts/parameter_estimation/run.sh reduce-strategy-endpoints CONFIG.json
```

Optimization arrays retain image-independent start shards. After every array, a small coordinator checks `sacct`; only missing tasks whose SLURM state is `TIMEOUT` are automatically retried once with three times the original wall time. Non-timeout failures stop recovery and are reported for inspection.

`--datasets` restricts an optimization submission wave to the listed dataset IDs while retaining a single immutable release and output hierarchy. The submitted job directory snapshots both the selected optimization plan and the filter.

`--defer-after-optim` queues NUTS before its optimization wave has fully drained.
Each NUTS fit receives a scheduler gate tied to its matching recovery coordinator;
the gate releases its chains only after all three combined optimization objects
exist.

Preparation refuses to overwrite an existing release. Validation records the exact count input, Stan-data SHA256, builder SHA256/mtime, and enforces the corrected SUM-159-fuse glucose mapping (912 glucose rows; affected batch 60 low/60 high; 44 censored).

The derived-data layer is release-local. `derive-optim` consolidates optimization
assessment metrics. After the selected NUTS fits complete, `submit-derived` runs
seven context-level reconstructed-parameter tasks and 35 context/model prediction
tasks, then validates them and writes `derived/manifest.tsv`. Prediction shards
contain chain-balanced trajectories and growth surfaces; figure-specific summaries
are intentionally left to in-memory transforms.

Posterior strategy simulations are an optional release-local stage. They write
compressed context/model/line qs files under
`derived/posterior/strategy_simulations/` without changing the completed generic
posterior-prediction layer. After validation, `reduce-strategy-endpoints` writes
one compact endpoint-only RDS shard per context under
`derived/posterior/strategy_simulations/endpoints/`; figure workflows should use
these reductions for schedule ranking and read full qs trajectories only for
the selected displayed schedules.

For a canonical external count object, set `counts.materialize` to `false` and pin `counts.expected_sha256`. Preparation then validates and records the source directly without copying or reserializing it.

For covariate-only releases that must remain exactly comparable to an existing
Stan input, configure `base_stan_data.input_path` and
`base_stan_data.expected_sha256`. Preparation copies that authenticated object
byte-for-byte instead of rebuilding it with the current `R/prepare_data.R`.
