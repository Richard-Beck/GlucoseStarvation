#!/bin/bash

set -euo pipefail

SPEC_PATH=${1:?Usage: bash slurm/submit_gpath_transfer_screen_pipeline.sh <spec_file.tsv>}

if [[ ! -f "$SPEC_PATH" ]]; then
  echo "Spec file not found: $SPEC_PATH" >&2
  exit 1
fi

mkdir -p slurm/logs
find slurm/logs -maxdepth 1 -type f -delete

BATCH_LABEL=$(date +%Y%m%d_%H%M%S)
RUN_DIR="slurm/runs/gpath_transfer_screen/${BATCH_LABEL}"
mkdir -p "$RUN_DIR"
tr -d '\r' < "$SPEC_PATH" > "${RUN_DIR}/spec.tsv"
SPEC_COPY="${RUN_DIR}/spec.tsv"

JOB_TABLE="${RUN_DIR}/job_ids.tsv"
printf "row_index\trun_id\ttransfer_run_dir\tarray_job_id\tsummary_job_id\tbatch_run_dir\n" > "$JOB_TABLE"

echo "PWD: $(pwd)"
echo "Submitting transfer screening batch from ${SPEC_PATH}"
echo "Batch metadata directory: ${RUN_DIR}"

Rscript - "$SPEC_COPY" <<'EOF' | while IFS='|' read -r row_index model_name run_id line_spec direction_spec fit_spec start_spec iter_warmup iter_sampling adapt_delta max_treedepth num_threads stan_data_path output_root ploidy_effect_mask_spec run_prefit prefit_chains init_mode qos; do
args <- commandArgs(trailingOnly = TRUE)
spec_path <- args[1]

spec <- read.delim(
  spec_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_cols <- c(
  "enabled", "model_name", "run_id", "line_spec", "direction_spec", "fit_spec", "start_spec",
  "iter_warmup", "iter_sampling", "adapt_delta", "max_treedepth", "num_threads",
  "stan_data_path", "output_root", "ploidy_effect_mask_spec", "run_prefit", "prefit_chains", "init_mode", "qos"
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
  if (enabled == "0") {
    next
  }

  vals <- c(
    i,
    nz_or(row$model_name, "gpath"),
    nz_or(row$run_id, ""),
    nz_or(row$line_spec, "all"),
    nz_or(row$direction_spec, "low_to_high,high_to_low"),
    nz_or(row$fit_spec, "null,transfer,oracle"),
    nz_or(row$start_spec, "1:50"),
    nz_or(row$iter_warmup, "500"),
    nz_or(row$iter_sampling, "1000"),
    nz_or(row$adapt_delta, "0.99"),
    nz_or(row$max_treedepth, "12"),
    nz_or(row$num_threads, "16"),
    nz_or(row$stan_data_path, "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"),
    nz_or(row$output_root, "data/gpath_transfer_cv"),
    nz_or(row$ploidy_effect_mask_spec, "all"),
    nz_or(row$run_prefit, "0"),
    nz_or(row$prefit_chains, "4"),
    nz_or(row$init_mode, "random"),
    nz_or(row$qos, "medium")
  )

  if (vals[3] == "") {
    next
  }

  cat(paste(vals, collapse = "|"), "\n", sep = "")
}
EOF
  run_parts=$(Rscript -e "source('models/gpath/v1/model.R'); x <- parse_run_id('$run_id'); cat(x\$R, x\$P, x\$W, x\$C, x\$M)")
  read -r r_val p_val w_val c_val m_val <<< "$run_parts"

  echo "Submitting transfer screening row ${row_index}: ${run_id} (${start_spec}, init_mode=${init_mode}, qos=${qos}, mask=${ploidy_effect_mask_spec})"

  submit_output=$(bash slurm/submit_gpath_transfer_pipeline.sh \
    "$model_name" \
    "$r_val" \
    "$p_val" \
    "$w_val" \
    "$c_val" \
    "$m_val" \
    "$line_spec" \
    "$direction_spec" \
    "$fit_spec" \
    "$start_spec" \
    "$iter_warmup" \
    "$iter_sampling" \
    "$adapt_delta" \
    "$max_treedepth" \
    "$num_threads" \
    "$stan_data_path" \
    "$output_root" \
    "$ploidy_effect_mask_spec" \
    "$run_prefit" \
    "$prefit_chains" \
    "$init_mode" \
    "$qos" \
    0)

  transfer_run_dir=$(printf "%s\n" "$submit_output" | awk -F= '/^run_dir=/{print $2; exit}')
  if [[ -z "$transfer_run_dir" ]]; then
    transfer_run_dir=$(printf "%s\n" "$submit_output" | awk '/Run metadata written to /{print $NF; exit}')
    transfer_run_dir=$(dirname "$transfer_run_dir")
  fi

  array_job_id=$(printf "%s\n" "$submit_output" | awk -F': ' '/^Array job submitted:/{print $2; exit}')
  summary_job_id=$(printf "%s\n" "$submit_output" | awk -F': ' '/^Summary job submitted:/{print $2; exit}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$row_index" "$run_id" "$transfer_run_dir" "$array_job_id" "$summary_job_id" "$RUN_DIR" >> "$JOB_TABLE"
done

echo "Submitted transfer screening batch. Job table: ${JOB_TABLE}"
