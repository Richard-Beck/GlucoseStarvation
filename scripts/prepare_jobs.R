#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

JOB_CONFIG_PATH <- if (length(args) >= 1) args[1] else stop("job config path required")
OUTPUT_DIR <- if (length(args) >= 2) args[2] else NULL

source("R/job_config_utils.R")

result <- prepare_jobs(JOB_CONFIG_PATH, output_dir = OUTPUT_DIR)

cat(sprintf("Prepared job configs in %s\n", result$config_dir))
cat(sprintf("Manifest: %s\n", result$manifest_path))
