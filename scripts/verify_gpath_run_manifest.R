#!/usr/bin/env Rscript

library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- if (length(args) >= 1) args[1] else stop("run_dir required")

source("R/run_manifest_utils.R")

manifest_path <- file.path(RUN_DIR, "run_manifest.json")
if (!file.exists(manifest_path)) {
  stop(sprintf("Manifest not found: %s", manifest_path))
}

manifest <- read_json(manifest_path, simplifyVector = TRUE)

current <- data.frame(
  field = c("stan_file_md5", "stan_data_md5", "script_md5"),
  manifest_value = c(manifest$stan_file_md5, manifest$stan_data_md5, manifest$script_md5),
  current_value = c(
    file_md5_or_na(manifest$stan_file_path),
    file_md5_or_na(manifest$stan_data_path),
    file_md5_or_na(manifest$script_path)
  ),
  stringsAsFactors = FALSE
)
current$match <- current$manifest_value == current$current_value

print(current, row.names = FALSE)
