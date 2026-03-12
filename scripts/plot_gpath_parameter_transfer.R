#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

MODEL_ID <- if (length(args) >= 1) args[1] else "1R_1P_0W_C0_M1"
OUTPUT_ROOT <- if (length(args) >= 2) args[2] else "data/gpath_transfer_cv"
MODEL_NAME <- if (length(args) >= 3) args[3] else "gpath"

source("R/parameter_transfer_plot_utils.R")

cat(sprintf(">>> plot_gpath_parameter_transfer.R cwd: %s\n", getwd()))

p1 <- plot_parameter_transfer_improvement(
  model_id = MODEL_ID,
  output_root = OUTPUT_ROOT,
  model_name = MODEL_NAME
)
print(p1)

p2 <- plot_parameter_transfer_summary(
  model_id = MODEL_ID,
  output_root = OUTPUT_ROOT,
  model_name = MODEL_NAME
)
print(p2)
