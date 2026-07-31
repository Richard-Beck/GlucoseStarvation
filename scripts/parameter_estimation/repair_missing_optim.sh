#!/usr/bin/env bash
set -euo pipefail

release_root=${1:?release root required}
init_source_root=${2:?matched optimization init source root required}
init_offset=${3:-100}
plan="${release_root}/manifests/optim_plan.tsv"

if [[ ! -f "${plan}" ]]; then
  echo "Optimization plan not found: ${plan}" >&2
  exit 1
fi
if ! [[ "${init_offset}" =~ ^[0-9]+$ ]]; then
  echo "init_offset must be a non-negative integer" >&2
  exit 1
fi

stamp=$(date +%Y%m%d_%H%M%S)
job_root="${release_root}/jobs/optim_repair_${stamp}"
log_root="${job_root}/logs"
mkdir -p "${log_root}"
printf 'fit_id\tmissing_tasks\trepair_job_id\tcombine_job_id\toutput_dir\n' > "${job_root}/job_ids.tsv"

while IFS=$'\t' read -r fit_id dataset_id model_id n_starts threads stan_data_path output_dir cpus mem_gb walltime qos; do
  if [[ "${fit_id}" == "fit_id" ]]; then
    continue
  fi

  missing=()
  for ((task_id=1; task_id<=n_starts; task_id++)); do
    if [[ ! -f "${output_dir}/optim_draws_${task_id}.Rds" ||
          ! -f "${output_dir}/optim_lp_${task_id}.Rds" ||
          ! -f "${output_dir}/optim_rc_${task_id}.Rds" ]]; then
      missing+=("${task_id}")
    fi
  done
  if [[ ${#missing[@]} -eq 0 ]]; then
    continue
  fi

  missing_csv=$(IFS=,; echo "${missing[*]}")
  repair_job_id=$(
    sbatch --parsable \
      --job-name="repair_${fit_id}" \
      --array="${missing_csv}" \
      --cpus-per-task="${cpus}" \
      --mem="${mem_gb}G" \
      --time="${walltime}" \
      --qos="${qos}" \
      --output="${log_root}/${fit_id}_%A_%a.out" \
      --error="${log_root}/${fit_id}_%A_%a.err" \
      --export="ALL,GPATH_OPTIM_INIT_SOURCE_ROOT=${init_source_root},GPATH_OPTIM_INIT_OFFSET=${init_offset},GPATH_OPTIM_INIT_RADIUS=0.2" \
      slurm/jobs/gpath_optim_array.sh \
      gpath v1 "${model_id}" "${threads}" "${stan_data_path}" "${dataset_id}" "${output_dir}"
  )
  combine_job_id=$(
    sbatch --parsable \
      --job-name="combine_${fit_id}" \
      --dependency="afterok:${repair_job_id}" \
      --cpus-per-task=2 \
      --mem=8G \
      --time=00:30:00 \
      --qos=small \
      --output="${log_root}/${fit_id}_combine_%j.out" \
      --error="${log_root}/${fit_id}_combine_%j.err" \
      --wrap="scripts/agentRrunner.sh pipelines/gpath/optim/combine.R '${output_dir}' 1"
  )

  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${fit_id}" "${missing_csv}" "${repair_job_id}" "${combine_job_id}" "${output_dir}" \
    >> "${job_root}/job_ids.tsv"
done < "${plan}"

echo "${job_root}"
