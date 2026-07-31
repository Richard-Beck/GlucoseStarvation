#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
derived_cfg <- pe_read_derived_config(cfg)
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
fit_plan <- pe_read_tsv(file.path(root, "manifests", "fit_plan.tsv"), c(
  "fit_id", "dataset_id", "model_id", "stan_data_path", "run_nuts", "nuts_chains", "nuts_dir"
))
nuts_fits <- fit_plan[as.integer(fit_plan$run_nuts) == 1L, , drop = FALSE]
if (!nrow(nuts_fits)) pe_fail("Fit plan contains no NUTS fits")
if (anyDuplicated(nuts_fits[, c("dataset_id", "model_id")])) pe_fail("NUTS plan repeats a dataset/model pair")

dataset_ids <- unique(nuts_fits$dataset_id)
parameter_tasks <- do.call(rbind, lapply(seq_along(dataset_ids), function(i) {
  dataset_id <- dataset_ids[[i]]
  rows <- nuts_fits[nuts_fits$dataset_id == dataset_id, , drop = FALSE]
  data.frame(
    task_id = i,
    dataset_id = dataset_id,
    n_models = nrow(rows),
    stan_data_path = unique(rows$stan_data_path)[[1]],
    output_path = file.path(derived_root, "posterior", "parameters", sprintf("%s.Rds", dataset_id)),
    stringsAsFactors = FALSE
  )
}))

prediction_tasks <- nuts_fits[, c(
  "fit_id", "dataset_id", "model_id", "stan_data_path", "nuts_chains", "nuts_dir"
), drop = FALSE]
prediction_tasks$task_id <- seq_len(nrow(prediction_tasks))
prediction_tasks$output_path <- file.path(
  derived_root, "posterior", "predictions", prediction_tasks$dataset_id,
  sprintf("%s.Rds", prediction_tasks$model_id)
)
prediction_tasks <- prediction_tasks[, c(
  "task_id", "fit_id", "dataset_id", "model_id", "stan_data_path", "nuts_chains", "nuts_dir", "output_path"
)]

resources <- do.call(rbind, lapply(names(derived_cfg$resources), function(component) {
  x <- derived_cfg$resources[[component]]
  data.frame(
    component = component,
    cpus = as.integer(x$cpus),
    mem_gb = as.integer(x$mem_gb),
    time = as.character(x$time),
    qos = as.character(x$qos),
    stringsAsFactors = FALSE
  )
}))

plan_root <- file.path(derived_root, "plans")
dir.create(plan_root, recursive = TRUE, showWarnings = FALSE)
utils::write.table(parameter_tasks, file.path(plan_root, "parameter_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
utils::write.table(prediction_tasks, file.path(plan_root, "prediction_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
utils::write.table(resources, file.path(plan_root, "resources.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
invisible(file.copy(derived_cfg$config_path, file.path(plan_root, "derived_config.json"), overwrite = TRUE))

cat(sprintf("Planned %d parameter shards and %d prediction shards under %s\n",
            nrow(parameter_tasks), nrow(prediction_tasks), derived_root))
