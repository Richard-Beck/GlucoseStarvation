#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

MODEL_ID <- if (length(args) >= 1) args[1] else "1R_1P_0W_C0_M1"
LINE_ID <- if (length(args) >= 2) as.integer(args[2]) else 1L
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else "data/gpath_transfer_cv"
STAN_DATA_PATH <- if (length(args) >= 4) args[4] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"

source("R/transfer_cv_plot_utils.R")

cat(sprintf(">>> plot_gpath_transfer_cv.R cwd: %s\n", getwd()))

transfer_data <- generate_transfer_overlay_data(
  model_id = MODEL_ID,
  line_id = LINE_ID,
  output_root = OUTPUT_ROOT,
  stan_data_path = STAN_DATA_PATH
)

plot_transfer_line_trajectories(
  transfer_data = transfer_data,
  line_id = LINE_ID,
  model_id = MODEL_ID
)
