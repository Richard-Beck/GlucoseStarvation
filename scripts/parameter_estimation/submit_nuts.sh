#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT=${1:?release root required}
DRY_RUN=${2:-0}
DEFERRED_OPTIM_JOBS=${3:-}
PLAN="${RELEASE_ROOT}/manifests/nuts_plan.tsv"
FIT_PLAN="${RELEASE_ROOT}/manifests/fit_plan.tsv"
[[ -f "$PLAN" && -f "$FIT_PLAN" ]] || { echo "Missing NUTS/fit plan" >&2; exit 1; }
[[ -z "$DEFERRED_OPTIM_JOBS" || -f "$DEFERRED_OPTIM_JOBS" ]] || {
  echo "Missing optimization job table: $DEFERRED_OPTIM_JOBS" >&2
  exit 1
}

stamp=$(date +%Y%m%d_%H%M%S)
job_root="${RELEASE_ROOT}/jobs/nuts_${stamp}"
log_root="${job_root}/logs"
mkdir -p "$log_root" "${job_root}/manifests"
printf 'fit_id\tarray_job_id\tgate_job_id\tn_tasks\n' > "${job_root}/job_ids.tsv"

cut -f2 "$PLAN" | tail -n +2 | sort -u | while read -r fit_id; do
  [[ -n "$fit_id" ]] || continue
  optim_dir=$(awk -F '\t' -v id="$fit_id" 'NR>1 && $1==id {print $8}' "$FIT_PLAN")
  optim_ready=0
  [[ -f "${optim_dir}/optim_draws_all.Rds" && -f "${optim_dir}/optim_lp_all.Rds" && -f "${optim_dir}/optim_rc_all.Rds" ]] && optim_ready=1
  recovery_id=""
  if [[ "$optim_ready" == 0 ]]; then
    [[ -n "$DEFERRED_OPTIM_JOBS" ]] || {
      echo "Optimization is not combined for $fit_id: $optim_dir" >&2
      exit 1
    }
    recovery_id=$(awk -F '\t' -v id="$fit_id" 'NR>1 && $1==id {print $3; exit}' "$DEFERRED_OPTIM_JOBS")
    [[ -n "$recovery_id" ]] || {
      echo "No recovery job found for deferred NUTS fit $fit_id" >&2
      exit 1
    }
  fi
  subset="${job_root}/manifests/${fit_id}.tsv"
  awk -F '\t' -v id="$fit_id" 'BEGIN{OFS="\t"} NR==1 {print "task_id","run_id","chain_id","config_path","output_dir","run_tag"; next} $2==id {n++; print n,$3,$4,$5,$6,$7}' "$PLAN" > "$subset"
  n_tasks=$(awk 'END{print NR-1}' "$subset")
  first=$(awk -F '\t' -v id="$fit_id" 'NR>1 && $2==id {print $8"\t"$9"\t"$10"\t"$11; exit}' "$PLAN")
  IFS=$'\t' read -r cpus mem_gb walltime qos <<< "$first"
  if [[ "$DRY_RUN" == 1 ]]; then
    printf 'DRY_RUN fit=%s chains=%s time=%s optim_ready=%s recovery_job=%s\n' \
      "$fit_id" "$n_tasks" "$walltime" "$optim_ready" "${recovery_id:-none}"
    continue
  fi
  dependency_args=()
  gate_id=""
  if [[ "$optim_ready" == 0 ]]; then
    gate_id=$(sbatch --parsable \
      --dependency="afterany:${recovery_id}" --job-name="pe_nuts_gate" \
      --cpus-per-task=1 --mem=1G --time=12:00:00 --qos="$qos" \
      --output="${log_root}/%x_%j.out" --error="${log_root}/%x_%j.err" \
      scripts/parameter_estimation/wait_for_optim_combine.sh "$optim_dir")
    dependency_args+=(--dependency="afterok:${gate_id}")
  fi
  array_id=$(sbatch --parsable "${dependency_args[@]}" \
    --job-name="pe_nuts" --array="1-${n_tasks}" \
    --cpus-per-task="$cpus" --mem="${mem_gb}G" --time="$walltime" --qos="$qos" \
    --output="${log_root}/%x_%A_%a.out" --error="${log_root}/%x_%A_%a.err" \
    slurm/jobs/gpath_full_oracle_array.sh "$subset" scripts/run_job.R)
  printf '%s\t%s\t%s\t%s\n' "$fit_id" "$array_id" "$gate_id" "$n_tasks" >> "${job_root}/job_ids.tsv"
done

echo "nuts_job_root=${job_root}"
