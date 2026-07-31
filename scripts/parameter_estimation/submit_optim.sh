#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT=${1:?release root required}
DRY_RUN=${2:-0}
DATASET_FILTER=${3:-}
PLAN="${RELEASE_ROOT}/manifests/optim_plan.tsv"
[[ -f "$PLAN" ]] || { echo "Missing optimization plan: $PLAN" >&2; exit 1; }

expected=$'fit_id\tdataset_id\tmodel_id\tn_starts\tthreads\tstan_data_path\toutput_dir\tcpus\tmem_gb\ttime\tqos'
[[ "$(head -n 1 "$PLAN")" == "$expected" ]] || { echo "Unexpected optimization plan schema" >&2; exit 1; }

stamp=$(date +%Y%m%d_%H%M%S)
job_root="${RELEASE_ROOT}/jobs/optim_${stamp}"
log_root="${job_root}/logs"
mkdir -p "$log_root"
selected_plan="${job_root}/optim_plan.tsv"
if [[ -n "$DATASET_FILTER" ]]; then
  [[ -f "$DATASET_FILTER" ]] || { echo "Dataset filter does not exist: $DATASET_FILTER" >&2; exit 1; }
  declare -A requested=()
  while IFS=$'\t' read -r dataset_id _; do
    [[ -z "$dataset_id" || "$dataset_id" == "dataset_id" || "$dataset_id" == \#* ]] && continue
    requested["$dataset_id"]=1
  done < "$DATASET_FILTER"
  ((${#requested[@]})) || { echo "Dataset filter is empty: $DATASET_FILTER" >&2; exit 1; }
  head -n 1 "$PLAN" > "$selected_plan"
  declare -A matched=()
  while IFS=$'\t' read -r fit_id dataset_id rest; do
    [[ "$fit_id" == "fit_id" ]] && continue
    if [[ -n "${requested[$dataset_id]:-}" ]]; then
      printf '%s\t%s\t%s\n' "$fit_id" "$dataset_id" "$rest" >> "$selected_plan"
      matched["$dataset_id"]=1
    fi
  done < "$PLAN"
  missing=()
  for dataset_id in "${!requested[@]}"; do
    [[ -n "${matched[$dataset_id]:-}" ]] || missing+=("$dataset_id")
  done
  ((${#missing[@]} == 0)) || { echo "Dataset filter contains unknown/unplanned ids: ${missing[*]}" >&2; exit 1; }
  cp "$DATASET_FILTER" "${job_root}/dataset_filter.tsv"
else
  cp "$PLAN" "$selected_plan"
fi
printf 'fit_id\tarray_job_id\trecovery_job_id\toutput_dir\n' > "${job_root}/job_ids.tsv"

n_fits=$(( $(wc -l < "$selected_plan") - 1 ))
n_starts=$(awk -F '\t' 'NR > 1 {n += $4} END {print n + 0}' "$selected_plan")
echo "selected_optimization_families=${n_fits} selected_starts=${n_starts}"

while IFS=$'\t' read -r fit_id dataset_id model_id n_starts threads stan_data_path output_dir cpus mem_gb walltime qos; do
  [[ "$fit_id" == "fit_id" ]] && continue
  [[ -f "$stan_data_path" ]] || { echo "Missing Stan data: $stan_data_path" >&2; exit 1; }
  if compgen -G "${output_dir}/optim_*_[0-9]*.Rds" >/dev/null; then
    echo "Refusing to overwrite existing optimization shards: $output_dir" >&2
    exit 1
  fi
  mkdir -p "$output_dir"
  if [[ "$DRY_RUN" == 1 ]]; then
    printf 'DRY_RUN fit=%s array=1-%s time=%s timeout_retry_time=3x\n' "$fit_id" "$n_starts" "$walltime"
    continue
  fi
  array_id=$(sbatch --parsable \
    --job-name="pe_opt" --array="1-${n_starts}" \
    --cpus-per-task="$cpus" --mem="${mem_gb}G" --time="$walltime" --qos="$qos" \
    --output="${log_root}/%x_%A_%a.out" --error="${log_root}/%x_%A_%a.err" \
    slurm/jobs/gpath_optim_array.sh \
    gpath v1 "$model_id" "$threads" "$stan_data_path" "$dataset_id" "$output_dir")
  record_path="${job_root}/recovery_${fit_id}.txt"
  recovery_id=$(sbatch --parsable \
    --dependency="afterany:${array_id}" --job-name="pe_opt_recover" \
    --cpus-per-task=1 --mem=1G --time=00:10:00 --qos="$qos" \
    --output="${log_root}/%x_%j.out" --error="${log_root}/%x_%j.err" \
    scripts/parameter_estimation/recover_optim_timeouts.sh \
    "$array_id" "$n_starts" "$walltime" "$output_dir" "$log_root" \
    gpath v1 "$model_id" "$threads" "$stan_data_path" "$dataset_id" \
    "$cpus" "$mem_gb" "$qos" "$record_path")
  printf '%s\t%s\t%s\t%s\n' "$fit_id" "$array_id" "$recovery_id" "$output_dir" >> "${job_root}/job_ids.tsv"
done < "$selected_plan"

echo "optimization_job_root=${job_root}"
