#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")
cat(pe_release_root(pe_read_config(args[[1]])))
