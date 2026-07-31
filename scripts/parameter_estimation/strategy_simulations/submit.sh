#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT=${1:?release root required}
CONFIG=${2:?release config required}
DRY_RUN=${3:-0}
STRATEGY_CONFIG=${4:-"$(dirname "$CONFIG")/strategy_simulations.json"}
export GPATH_STRATEGY_CONFIG="$STRATEGY_CONFIG"
[[ -f "$STRATEGY_CONFIG" ]] || { echo "Missing strategy config: $STRATEGY_CONFIG" >&2; exit 1; }

read -r cpus mem_gb walltime qos validator_cpus validator_mem validator_time max_concurrent output_root < <(
  scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/read_resources.R "$CONFIG"
)
PLAN="${output_root}/plans/tasks.tsv"
[[ -f "$PLAN" ]] || { echo "Missing strategy plan: $PLAN" >&2; exit 1; }
n_tasks=$(( $(wc -l < "$PLAN") - 1 ))
if (( max_concurrent > n_tasks )); then max_concurrent=$n_tasks; fi
if [[ "$DRY_RUN" == 1 ]]; then
  echo "DRY_RUN strategy_tasks=${n_tasks} max_concurrent=${max_concurrent} cpus=${cpus} mem_gb=${mem_gb} time=${walltime} qos=${qos} output_root=${output_root}"
  exit 0
fi

stamp=$(date +%Y%m%d_%H%M%S)
job_root="${RELEASE_ROOT}/jobs/strategy_simulations_${stamp}"
log_root="${job_root}/logs"
mkdir -p "$log_root"
cp "$PLAN" "$STRATEGY_CONFIG" "$job_root/"
mkdir -p "$job_root/source_snapshot"
cp \
  scripts/parameter_estimation/strategy_simulations/make_plan.R \
  scripts/parameter_estimation/strategy_simulations/simulate_task.R \
  scripts/parameter_estimation/strategy_simulations/validate.R \
  scripts/parameter_estimation/strategy_simulations/submit.sh \
  scripts/parameter_estimation/strategy_simulations/read_resources.R \
  R/gpath_posterior_strategy_utils.R \
  R/selection_strategy_utils.R \
  "$job_root/source_snapshot/"

array_id=$(sbatch --parsable --job-name=pe_strategy \
  --array="1-${n_tasks}%${max_concurrent}" --cpus-per-task="$cpus" --mem="${mem_gb}G" \
  --time="$walltime" --qos="$qos" \
  --export="ALL,GPATH_STRATEGY_CONFIG=${STRATEGY_CONFIG}" \
  --output="${log_root}/%x_%A_%a.out" --error="${log_root}/%x_%A_%a.err" \
  slurm/jobs/gpath_derived_array.sh \
  scripts/parameter_estimation/strategy_simulations/simulate_task.R "$CONFIG")

validator_id=$(sbatch --parsable --dependency="afterok:${array_id}" \
  --job-name=pe_strategy_validate --cpus-per-task="$validator_cpus" \
  --mem="${validator_mem}G" --time="$validator_time" --qos="$qos" \
  --export="ALL,GPATH_STRATEGY_CONFIG=${STRATEGY_CONFIG}" \
  --output="${log_root}/%x_%j.out" --error="${log_root}/%x_%j.err" \
  slurm/jobs/gpath_derived_array.sh \
  scripts/parameter_estimation/strategy_simulations/validate.R "$CONFIG")

printf 'component\tjob_id\narray\t%s\nvalidator\t%s\n' \
  "$array_id" "$validator_id" > "$job_root/job_ids.tsv"
echo "strategy_job_root=${job_root}"
