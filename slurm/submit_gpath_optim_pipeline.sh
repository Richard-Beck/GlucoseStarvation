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

Rscript - "$SPEC_COPY" "$BATCH_LABEL" <<'EOF' | while IFS='|' read -r row_index enabled model_name model_version run_id n_starts num_threads stan_data_path dataset_label run_label array_cpus array_mem_gb array_time combine_mem_gb combine_time qos delete_after; do
args <- commandArgs(trailingOnly = TRUE)
spec_path <- args[1]
batch_label <- args[2]

spec <- read.delim(
  spec_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_cols <- c(
  "enabled", "model_name", "model_version", "run_id", "n_starts", "num_threads",
  "stan_data_path", "dataset_label", "run_label", "array_cpus", "array_mem_gb",
  "array_time", "combine_mem_gb", "combine_time", "qos", "delete_after"
)

missing_cols <- setdiff(required_cols, names(spec))
if (length(missing_cols)) {
  stop(sprintf("Spec file is missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

nz_or <- function(x, default) {
  x <- trimws(as.character(x))
  ifelse(is.na(x) | x == "", default, x)
}

for (i in seq_len(nrow(spec))) {
  row <- spec[i, , drop = FALSE]

  enabled <- nz_or(row$enabled, "1")
  run_id <- nz_or(row$run_id, "")
  if (enabled == "0" || run_id == "") {
    next
  }

  vals <- c(
    i,
    enabled,
    nz_or(row$model_name, "gpath"),
    nz_or(row$model_version, "v1"),
    run_id,
    nz_or(row$n_starts, "250"),
    nz_or(row$num_threads, "8"),
    nz_or(row$stan_data_path, "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"),
    nz_or(row$dataset_label, "gstarvation_v1"),
    nz_or(row$run_label, batch_label),
    nz_or(row$array_cpus, nz_or(row$num_threads, "8")),
    nz_or(row$array_mem_gb, "8"),
    nz_or(row$array_time, "00:15:00"),
    nz_or(row$combine_mem_gb, "4"),
    nz_or(row$combine_time, "00:15:00"),
    nz_or(row$qos, "normal"),
    nz_or(row$delete_after, "1")
  )

  cat(paste(vals, collapse = "|"), "\n", sep = "")
}
EOF
  output_dir=$(
    MODEL_NAME="$model_name" \
    MODEL_VERSION="$model_version" \
    DATASET_LABEL="$dataset_label" \
    RUN_ID="$run_id" \
    RUN_LABEL="$run_label" \
    Rscript -e 'source("R/project_paths.R"); cat(get_run_output_dir(model_name = Sys.getenv("MODEL_NAME"), model_version = Sys.getenv("MODEL_VERSION"), pipeline_name = "optim", dataset_label = Sys.getenv("DATASET_LABEL"), run_id = Sys.getenv("RUN_ID"), run_label = Sys.getenv("RUN_LABEL")))' \
  )

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
