#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

CONFIG_PATH <- if (length(args) >= 1) args[1] else stop("config path required")

source("R/job_config_utils.R")
run_job(CONFIG_PATH)
