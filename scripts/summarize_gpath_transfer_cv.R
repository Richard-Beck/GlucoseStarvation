#!/usr/bin/env Rscript

library(posterior)

args <- commandArgs(trailingOnly = TRUE)

MODEL_NAME <- if (length(args) >= 1) args[1] else "gpath"
RUN_ID <- if (length(args) >= 2) args[2] else "1R_1P_0W_C0_M1"
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else "data/gpath_transfer_cv"

source("R/gpath_run_utils.R")
source("R/posterior_io_utils.R")
source("R/elpd_transfer_utils.R")
source("R/parameter_transfer_utils.R")

cat(sprintf(">>> summarize_gpath_transfer_cv.R cwd: %s\n", getwd()))

base_dir <- file.path(OUTPUT_ROOT, MODEL_NAME, RUN_ID)
if (!dir.exists(base_dir)) {
  stop(sprintf("Transfer CV output directory not found: %s", base_dir))
}

dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
dirs <- dirs[dirs != base_dir]

rows <- lapply(dirs, function(dir_path) {
  start_files <- list.files(dir_path, pattern = "^start_summary_.*\\.Rds$", full.names = TRUE)
  if (!length(start_files)) {
    return(NULL)
  }

  do.call(rbind, lapply(start_files, readRDS))
})

rows <- Filter(Negate(is.null), rows)
if (!length(rows)) {
  stop(sprintf("No chain summaries found under %s", base_dir))
}

summary_df <- do.call(rbind, rows)

best_idx <- unlist(lapply(
  split(seq_len(nrow(summary_df)), interaction(summary_df$line_id, summary_df$direction, summary_df$fit_type, drop = TRUE)),
  function(idx) idx[which.max(summary_df$lp__[idx])]
))
best_df <- summary_df[sort(best_idx), ]

wide <- reshape(
  best_df[, c("line_id", "direction", "fit_type", "elpd")],
  idvar = c("line_id", "direction"),
  timevar = "fit_type",
  direction = "wide"
)

names(wide) <- sub("^elpd\\.", "", names(wide))
wide$transfer_gain <- wide$transfer - wide$null
wide$transfer_regret <- wide$oracle - wide$transfer

saveRDS(summary_df, file.path(base_dir, "transfer_start_summaries.Rds"))
saveRDS(best_df, file.path(base_dir, "transfer_best_start_summary.Rds"))
saveRDS(wide, file.path(base_dir, "transfer_comparison_summary.Rds"))

  param_tables <- build_parameter_transfer_tables(
  model_id = RUN_ID,
  output_root = OUTPUT_ROOT,
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  model_name = MODEL_NAME
)
saveRDS(param_tables$parameter_states, file.path(base_dir, "parameter_transfer_states.Rds"))
saveRDS(param_tables$parameter_shifts, file.path(base_dir, "parameter_transfer_shifts.Rds"))
saveRDS(param_tables$parameter_comparison, file.path(base_dir, "parameter_transfer_comparison.Rds"))
saveRDS(param_tables$comparison_summary, file.path(base_dir, "parameter_transfer_summary.Rds"))

print(wide[order(wide$direction, wide$line_id), ], row.names = FALSE)
