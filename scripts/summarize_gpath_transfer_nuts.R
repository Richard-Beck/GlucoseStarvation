#!/usr/bin/env Rscript

library(posterior)

args <- commandArgs(trailingOnly = TRUE)

MODEL_NAME <- if (length(args) >= 1) args[1] else "gpath"
RUN_ID <- if (length(args) >= 2) args[2] else "1R_1P_0W_C0_M1"
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else "data/gpath_transfer_cv_nuts"

source("R/posterior_io_utils.R")
source("R/elpd_transfer_utils.R")

cat(sprintf(">>> summarize_gpath_transfer_nuts.R cwd: %s\n", getwd()))

base_dir <- file.path(OUTPUT_ROOT, MODEL_NAME, RUN_ID)
if (!dir.exists(base_dir)) {
  stop(sprintf("Transfer CV NUTS output directory not found: %s", base_dir))
}

dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
dirs <- dirs[dirs != base_dir]

condition_rows <- list()

for (dir_path in dirs) {
  draw_files <- list.files(dir_path, pattern = "^nuts_draws_.*\\.Rds$", full.names = TRUE)
  meta_files <- list.files(dir_path, pattern = "^split_meta_.*\\.Rds$", full.names = TRUE)
  if (!length(draw_files) || !length(meta_files)) next

  draw_tags <- sub("^nuts_draws_(.*)\\.Rds$", "\\1", basename(draw_files))
  meta_tags <- sub("^split_meta_(.*)\\.Rds$", "\\1", basename(meta_files))

  common_tags <- intersect(draw_tags, meta_tags)
  if (!length(common_tags)) {
    next
  }

  group_keys <- sub("_chain[0-9]+$", "", common_tags)

  for (group_key in unique(group_keys)) {
    idx <- which(group_keys == group_key)
    group_draw_files <- draw_files[match(common_tags[idx], draw_tags)]
    group_meta_files <- meta_files[match(common_tags[idx], meta_tags)]

    dr <- bind_rds_draws(group_draw_files, along = "chain")
    meta <- readRDS(group_meta_files[[1]])
    ss <- summarize_transfer_draws(dr, meta)
    ss$n_chains <- length(group_draw_files)

    diags <- suppressWarnings(summarise_draws(dr, rhat = rhat, ess_bulk = ess_bulk))
    ss$bad_rhat <- sum(diags$rhat > 1.01, na.rm = TRUE)
    ss$bad_ess <- sum(diags$ess_bulk < 400, na.rm = TRUE)

    condition_rows[[length(condition_rows) + 1L]] <- ss
  }
}

if (!length(condition_rows)) {
  stop(sprintf("No NUTS transfer outputs found under %s", base_dir))
}

condition_df <- do.call(rbind, condition_rows)

wide <- reshape(
  condition_df[, c("line_id", "direction", "fit_type", "elpd")],
  idvar = c("line_id", "direction"),
  timevar = "fit_type",
  direction = "wide"
)
names(wide) <- sub("^elpd\\.", "", names(wide))
wide$transfer_gain <- wide$transfer - wide$null
wide$transfer_regret <- wide$oracle - wide$transfer

saveRDS(condition_df, file.path(base_dir, "transfer_condition_summaries.Rds"))
saveRDS(wide, file.path(base_dir, "transfer_comparison_summary.Rds"))

print(wide[order(wide$direction, wide$line_id), ], row.names = FALSE)
