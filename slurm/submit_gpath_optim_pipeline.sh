#!/bin/bash

set -euo pipefail

SPEC_PATH=${1:?Usage: bash slurm/submit_gpath_optim_pipeline.sh <spec_file.tsv>}

if [[ ! -f "$SPEC_PATH" ]]; then
  echo "Spec file not found: $SPEC_PATH" >&2
  exit 1
fi

mkdir -p slurm/logs
find slurm/logs -maxdepth 1 -type f -delete

BATCH_LABEL=$(date +%Y%m%d_%H%M%S)
RUN_DIR="slurm/runs/gpath_optim_batch/${BATCH_LABEL}"
mkdir -p "$RUN_DIR"
tr -d '\r' < "$SPEC_PATH" > "${RUN_DIR}/spec.tsv"
SPEC_COPY="${RUN_DIR}/spec.tsv"

JOB_TABLE="${RUN_DIR}/job_ids.tsv"
printf "row_index\trun_id\trun_label\toutput_dir\tarray_job_id\tcombine_job_id\n" > "$JOB_TABLE"

echo "PWD: $(pwd)"
echo "Submitting optimisation batch from ${SPEC_PATH}"
echo "Batch metadata directory: ${RUN_DIR}"

trim() {
  local x="$1"
  x="${x#"${x%%[![:space:]]*}"}"
  x="${x%"${x##*[![:space:]]}"}"
  printf '%s' "$x"
}

nz_or() {
  local value
  value=$(trim "${1-}")
  local default="$2"
  if [[ -z "$value" ]]; then
    printf '%s' "$default"
  else
    printf '%s' "$value"
  fi
}

header_expected=$'enabled\tmodel_name\tmodel_version\trun_id\tn_starts\tnum_threads\tstan_data_path\tdataset_label\trun_label\tarray_cpus\tarray_mem_gb\tarray_time\tcombine_mem_gb\tcombine_time\tqos\tdelete_after'
header_actual=$(head -n 1 "$SPEC_COPY")
if [[ "$header_actual" != "$header_expected" ]]; then
  echo "Spec file header does not match expected columns." >&2
  echo "Expected: $header_expected" >&2
  echo "Actual:   $header_actual" >&2
  exit 1
fi

row_index=0
tail -n +2 "$SPEC_COPY" | while IFS=$'\t' read -r enabled model_name model_version run_id n_starts num_threads stan_data_path dataset_label run_label array_cpus array_mem_gb array_time combine_mem_gb combine_time qos delete_after; do
  row_index=$((row_index + 1))

  enabled=$(nz_or "$enabled" "1")
  run_id=$(nz_or "$run_id" "")
  if [[ "$enabled" == "0" || -z "$run_id" ]]; then
    continue
  fi

  model_name=$(nz_or "$model_name" "gpath")
  model_version=$(nz_or "$model_version" "v1")
  n_starts=$(nz_or "$n_starts" "250")
  num_threads=$(nz_or "$num_threads" "8")
  stan_data_path=$(nz_or "$stan_data_path" "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds")
  dataset_label=$(nz_or "$dataset_label" "gstarvation_v1")
  run_label=$(nz_or "$run_label" "$BATCH_LABEL")
  array_cpus=$(nz_or "$array_cpus" "$num_threads")
  array_mem_gb=$(nz_or "$array_mem_gb" "8")
  array_time=$(nz_or "$array_time" "00:15:00")
  combine_mem_gb=$(nz_or "$combine_mem_gb" "4")
  combine_time=$(nz_or "$combine_time" "00:15:00")
  qos=$(nz_or "$qos" "normal")
  delete_after=$(nz_or "$delete_after" "1")

  output_dir="data/runs/${model_name}/${model_version}/optim/${dataset_label}/${run_id}/${run_label}"

  mkdir -p "$output_dir"

  echo "Submitting row ${row_index}: ${run_id} (${n_starts} starts) -> ${output_dir}"

  array_job_id=$(
    sbatch --parsable \
      --array="1-${n_starts}" \
      --cpus-per-task="${array_cpus}" \
      --mem="${array_mem_gb}G" \
      --time="${array_time}" \
      --qos="${qos}" \
      slurm/jobs/gpath_optim_array.sh \
      "$model_name" \
      "$model_version" \
      "$run_id" \
      "$num_threads" \
      "$stan_data_path" \
      "$dataset_label" \
      "$output_dir"
  )

  combine_job_id=$(
    sbatch --parsable \
      --dependency="afterany:${array_job_id}" \
      --mem="${combine_mem_gb}G" \
      --time="${combine_time}" \
      --qos="${qos}" \
      slurm/jobs/gpath_optim_combine.sh \
      "$output_dir" \
      "$delete_after"
  )

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$row_index" "$run_id" "$run_label" "$output_dir" "$array_job_id" "$combine_job_id" >> "$JOB_TABLE"
done

echo "Submitted optimisation batch. Job table: ${JOB_TABLE}"
