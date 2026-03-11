# Slurm Layout

This directory is the non-`ecology` home for cluster submission code.

Pattern:

- `jobs/` contains stable job scripts that do one thing on the cluster.
- `submit_*.sh` scripts are thin local wrappers that create manifests, submit jobs, wire dependencies, and record job IDs.
- `runs/<pipeline>/<run_id>/<timestamp>/` stores generated manifests and submission metadata for a specific launch.
- `logs/` is the shared default target for Slurm stdout/stderr.

Recommended conventions for future pipelines:

- Use manifest-driven array jobs when task grids come from data rather than simple integer formulas.
- Keep model logic in `R/` or `scripts/`; Slurm scripts should only decode task metadata and invoke those entry points.
- Write one `submit_*.sh` wrapper per pipeline so a full model run can be launched with a single `sbatch`-free command.
- Save `job_ids.txt` or similar metadata under `slurm/runs/...` so jobs can be traced after submission.
