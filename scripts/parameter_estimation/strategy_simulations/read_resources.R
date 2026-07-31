#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: read_resources.R RELEASE_CONFIG", call. = FALSE)
source("scripts/parameter_estimation/common.R")
source("R/gpath_posterior_strategy_utils.R")
cfg <- pe_read_config(args[[1]])
strategy_cfg <- gps_read_strategy_config(cfg)
x <- strategy_cfg$resources
cat(
  as.integer(x$cpus), as.integer(x$mem_gb), as.character(x$time),
  as.character(x$qos), as.integer(x$validator_cpus),
  as.integer(x$validator_mem_gb), as.character(x$validator_time),
  if (is.null(x$max_concurrent)) 1000001L else as.integer(x$max_concurrent),
  gps_output_root(cfg, strategy_cfg), "\n"
)
