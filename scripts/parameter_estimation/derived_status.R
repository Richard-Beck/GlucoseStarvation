#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
derived_root <- pe_derived_root(cfg)
param_path <- file.path(derived_root, "plans", "parameter_tasks.tsv")
pred_path <- file.path(derived_root, "plans", "prediction_tasks.tsv")
if (!file.exists(param_path) || !file.exists(pred_path)) pe_fail("Derived plans are absent; run derive-plan")
parameters <- pe_read_tsv(param_path, c("output_path"))
predictions <- pe_read_tsv(pred_path, c("output_path"))
cat(sprintf("optimization_assessment=%s\n", file.exists(file.path(derived_root, "optimization", "assessment.Rds"))))
cat(sprintf("nuts_qc=%s\n", file.exists(file.path(derived_root, "posterior", "qc.Rds"))))
cat(sprintf("parameter_shards=%d/%d\n", sum(file.exists(parameters$output_path)), nrow(parameters)))
cat(sprintf("prediction_shards=%d/%d\n", sum(file.exists(predictions$output_path)), nrow(predictions)))
cat(sprintf("validated=%s\n", file.exists(file.path(derived_root, "validation.json")) && isTRUE(jsonlite::fromJSON(file.path(derived_root, "validation.json"))$validation_passed)))
