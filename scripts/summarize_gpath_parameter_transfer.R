#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

MODEL_ID <- if (length(args) >= 1) args[1] else "1R_1P_0W_C0_M1"
OUTPUT_ROOT <- if (length(args) >= 2) args[2] else "data/gpath_transfer_cv"
STAN_DATA_PATH <- if (length(args) >= 3) args[3] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"
MODEL_NAME <- if (length(args) >= 4) args[4] else "gpath"

source("R/parameter_transfer_utils.R")

cat(sprintf(">>> summarize_gpath_parameter_transfer.R cwd: %s\n", getwd()))

tables <- build_parameter_transfer_tables(
  model_id = MODEL_ID,
  output_root = OUTPUT_ROOT,
  stan_data_path = STAN_DATA_PATH,
  model_name = MODEL_NAME
)

base_dir <- file.path(OUTPUT_ROOT, MODEL_NAME, MODEL_ID)
saveRDS(tables$parameter_states, file.path(base_dir, "parameter_transfer_states.Rds"))
saveRDS(tables$parameter_shifts, file.path(base_dir, "parameter_transfer_shifts.Rds"))
saveRDS(tables$parameter_comparison, file.path(base_dir, "parameter_transfer_comparison.Rds"))
saveRDS(tables$comparison_summary, file.path(base_dir, "parameter_transfer_summary.Rds"))

print(tables$comparison_summary[order(tables$comparison_summary$direction, tables$comparison_summary$line_id), ], row.names = FALSE)
