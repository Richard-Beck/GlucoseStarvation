#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

MODEL_ID <- if (length(args) >= 1) args[1] else "1R_1P_0W_C0_M1"
LINE_ID <- if (length(args) >= 2) as.integer(args[2]) else 1L
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else "data/gpath_transfer_cv"
STAN_DATA_PATH <- if (length(args) >= 4) args[4] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"

source("R/transfer_cv_plot_utils.R")

cat(sprintf(">>> plot_gpath_transfer_cv.R cwd: %s\n", getwd()))

plot_data <- build_transfer_overlay_plot_data(
  model_id = MODEL_ID,
  line_id = LINE_ID,
  output_root = OUTPUT_ROOT,
  stan_data_path = STAN_DATA_PATH
)

plot_fit_overlays(
  sim_df = plot_data$sim_df,
  obs_df = plot_data$obs_df,
  color_by = "group_1",
  linetype_by = "group_2",
  color_values = c(null = "#6f6f6f", transfer = "#1b9e77", oracle = "#d95f02"),
  linetype_values = c(low_to_high = "solid", high_to_low = "22"),
  color_label = "Fit",
  linetype_label = "Direction",
  title = sprintf(
    "Transfer CV Fits | Model: %s | Cell Line: %s",
    MODEL_ID,
    plot_data$context$line_name
  ),
  subtitle = plot_data$score_text
)
