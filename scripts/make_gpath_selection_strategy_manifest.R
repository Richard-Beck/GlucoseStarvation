#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
CONFIG_PATH <- if (length(args) >= 1) args[1] else file.path("scripts", "selection_strategy_config.R")
OUTPUT_PATH <- if (length(args) >= 2) args[2] else stop("output manifest path required")

source("R/selection_strategy_utils.R")

cfg <- load_selection_config(CONFIG_PATH)
stan_data <- resolve_selection_stan(cfg)
tasks <- build_selection_tasks(cfg, stan_data)

task_df <- data.frame(task_id = seq_len(nrow(tasks)), tasks, stringsAsFactors = FALSE)

dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  task_df,
  file = OUTPUT_PATH,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Wrote %d tasks to %s\n", nrow(task_df), OUTPUT_PATH))
