#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT=${1:?release root required}
CONFIG=${2:?release config required}
DRY_RUN=${3:-0}
PLAN_ROOT="${RELEASE_ROOT}/derived/plans"
PARAM_PLAN="${PLAN_ROOT}/parameter_tasks.tsv"
PRED_PLAN="${PLAN_ROOT}/prediction_tasks.tsv"
RESOURCES="${PLAN_ROOT}/resources.tsv"
[[ -f "$PARAM_PLAN" && -f "$PRED_PLAN" && -f "$RESOURCES" ]] || { echo "Missing derived plans" >&2; exit 1; }

read_resource() {
  awk -F '\t' -v component="$1" 'NR>1 && $1==component {print $2"\t"$3"\t"$4"\t"$5}' "$RESOURCES"
}

n_param=$(( $(wc -l < "$PARAM_PLAN") - 1 ))
n_pred=$(( $(wc -l < "$PRED_PLAN") - 1 ))
if [[ "$DRY_RUN" == 1 ]]; then
  echo "DRY_RUN parameter_tasks=${n_param} prediction_tasks=${n_pred}"
  exit 0
fi

FIT_PLAN="${RELEASE_ROOT}/manifests/fit_plan.tsv"
while IFS=$'\t' read -r fit_id dataset_id model_id stan_path stan_sha run_optim optim_starts optim_dir run_nuts nuts_chains nuts_dir; do
  [[ "$fit_id" == "fit_id" || "$run_nuts" != 1 ]] && continue
  for ((chain=1; chain<=nuts_chains; chain++)); do
    path=$(printf '%s/nuts_draws_chain%02d.Rds' "$nuts_dir" "$chain")
    [[ -f "$path" ]] || { echo "Missing NUTS chain: $path" >&2; exit 1; }
  done
done < "$FIT_PLAN"

stamp=$(date +%Y%m%d_%H%M%S)
job_root="${RELEASE_ROOT}/jobs/derived_${stamp}"
log_root="${job_root}/logs"
mkdir -p "$log_root"
cp "$PARAM_PLAN" "$PRED_PLAN" "$RESOURCES" "$job_root/"

IFS=$'\t' read -r qc_cpus qc_mem qc_time qc_qos <<< "$(read_resource qc)"
IFS=$'\t' read -r par_cpus par_mem par_time par_qos <<< "$(read_resource parameters)"
IFS=$'\t' read -r pred_cpus pred_mem pred_time pred_qos <<< "$(read_resource predictions)"
IFS=$'\t' read -r fin_cpus fin_mem fin_time fin_qos <<< "$(read_resource finalize)"

qc_id=$(sbatch --parsable --job-name=pe_derived_qc \
  --cpus-per-task="$qc_cpus" --mem="${qc_mem}G" --time="$qc_time" --qos="$qc_qos" \
  --output="${log_root}/%x_%j.out" --error="${log_root}/%x_%j.err" \
  slurm/jobs/gpath_derived_array.sh scripts/parameter_estimation/build_nuts_qc.R "$CONFIG")
param_id=$(sbatch --parsable --job-name=pe_derived_parameters --array="1-${n_param}" \
  --cpus-per-task="$par_cpus" --mem="${par_mem}G" --time="$par_time" --qos="$par_qos" \
  --output="${log_root}/%x_%A_%a.out" --error="${log_root}/%x_%A_%a.err" \
  slurm/jobs/gpath_derived_array.sh scripts/parameter_estimation/build_posterior_parameters.R "$CONFIG")
pred_id=$(sbatch --parsable --job-name=pe_derived_predictions --array="1-${n_pred}" \
  --cpus-per-task="$pred_cpus" --mem="${pred_mem}G" --time="$pred_time" --qos="$pred_qos" \
  --output="${log_root}/%x_%A_%a.out" --error="${log_root}/%x_%A_%a.err" \
  slurm/jobs/gpath_derived_array.sh scripts/parameter_estimation/build_posterior_predictions.R "$CONFIG")
final_id=$(sbatch --parsable --dependency="afterok:${qc_id}:${param_id}:${pred_id}" --job-name=pe_derived_finalize \
  --cpus-per-task="$fin_cpus" --mem="${fin_mem}G" --time="$fin_time" --qos="$fin_qos" \
  --output="${log_root}/%x_%j.out" --error="${log_root}/%x_%j.err" \
  slurm/jobs/gpath_derived_array.sh scripts/parameter_estimation/validate_derived.R "$CONFIG")

printf 'component\tjob_id\nqc\t%s\nparameters\t%s\npredictions\t%s\nfinalize\t%s\n' \
  "$qc_id" "$param_id" "$pred_id" "$final_id" > "$job_root/job_ids.tsv"
echo "derived_job_root=${job_root}"
