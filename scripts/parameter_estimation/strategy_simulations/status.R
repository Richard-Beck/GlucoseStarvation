#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: status.R RELEASE_CONFIG", call. = FALSE)
source("scripts/parameter_estimation/common.R")
source("R/gpath_posterior_strategy_utils.R")
cfg <- pe_read_config(args[[1]])
strategy_cfg <- gps_read_strategy_config(cfg)
root <- gps_output_root(cfg, strategy_cfg)
tasks <- pe_read_tsv(file.path(root, "plans", "tasks.tsv"), c("task_id", "output_path"))
exists <- file.exists(tasks$output_path)
sizes <- rep(0, nrow(tasks))
sizes[exists] <- file.info(tasks$output_path[exists])$size
cat(sprintf("strategy_outputs=%d/%d\n", sum(exists), length(exists)))
cat(sprintf("bytes_written=%.0f\n", sum(sizes, na.rm = TRUE)))
if (any(!exists)) {
  cat("missing_task_ids=", paste(tasks$task_id[!exists], collapse = ","), "\n", sep = "")
}
if (!identical(as.character(strategy_cfg$storage$mode), "compact_observation_quantiles")) {
  endpoint_paths <- file.path(
    root, "endpoints",
    paste0(unique(as.character(tasks$dataset_id)), ".Rds")
  )
  cat(sprintf("endpoint_shards=%d/%d\n",
              sum(file.exists(endpoint_paths)), length(endpoint_paths)))
}
