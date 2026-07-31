#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: make_plan.R RELEASE_CONFIG", call. = FALSE)
suppressPackageStartupMessages(library(jsonlite))
source("scripts/parameter_estimation/common.R")
source("R/project_paths.R")
source("R/gpath_posterior_strategy_utils.R")

cfg <- pe_read_config(args[[1]])
strategy_cfg <- gps_read_strategy_config(cfg)
release_root <- pe_release_root(cfg)
output_root <- gps_output_root(cfg, strategy_cfg)
fit_plan <- pe_read_tsv(
  file.path(release_root, "manifests", "fit_plan.tsv"),
  c("fit_id", "dataset_id", "model_id", "stan_data_path", "run_nuts", "nuts_chains", "nuts_dir")
)

fit_plan <- fit_plan[
  as.integer(fit_plan$run_nuts) == 1L &
    fit_plan$dataset_id %in% strategy_cfg$dataset_ids &
    fit_plan$model_id %in% strategy_cfg$model_ids,
  , drop = FALSE
]
expected_fits <- length(strategy_cfg$dataset_ids) * length(strategy_cfg$model_ids)
if (nrow(fit_plan) != expected_fits) {
  pe_fail("Expected %d strategy NUTS fits; found %d", expected_fits, nrow(fit_plan))
}

source(get_model_r_path(cfg$model_name, cfg$model_version))
source("R/gpath_run_utils.R")
source("R/gpath_derived_utils.R")
source("R/selection_strategy_utils.R")

rows <- list()
ptr <- 1L
for (fit_i in seq_len(nrow(fit_plan))) {
  fit <- fit_plan[fit_i, , drop = FALSE]
  stan_data <- add_group_structure(readRDS(fit$stan_data_path[[1]]))
  lines <- gpd_line_table(stan_data)
  if (anyNA(lines$line_name)) pe_fail("Unresolved line names in %s", fit$stan_data_path[[1]])
  for (line_i in seq_len(nrow(lines))) {
    line_name <- lines$line_name[[line_i]]
    rows[[ptr]] <- data.frame(
      task_id = ptr,
      fit_id = fit$fit_id[[1]],
      dataset_id = fit$dataset_id[[1]],
      model_id = fit$model_id[[1]],
      line_id = lines$line_id[[line_i]],
      line_name = line_name,
      stan_data_path = fit$stan_data_path[[1]],
      nuts_chains = as.integer(fit$nuts_chains[[1]]),
      nuts_dir = fit$nuts_dir[[1]],
      output_path = file.path(
        output_root, fit$dataset_id[[1]], fit$model_id[[1]],
        paste0(gps_slug(line_name), ".qs")
      ),
      stringsAsFactors = FALSE
    )
    ptr <- ptr + 1L
  }
}
tasks <- do.call(rbind, rows)
if (nrow(tasks) != 70L || anyDuplicated(tasks[, c("dataset_id", "model_id", "line_id")])) {
  pe_fail("Expected 70 unique context/model/line tasks; found %d", nrow(tasks))
}

schedule <- gps_build_strategy_panel(strategy_cfg)
expected_schedules <- if (is.null(strategy_cfg$expected_schedules)) {
  448L
} else {
  as.integer(strategy_cfg$expected_schedules)
}
if (nrow(schedule) != expected_schedules) {
  pe_fail("Expected %d schedules; found %d", expected_schedules, nrow(schedule))
}

plan_root <- file.path(output_root, "plans")
dir.create(plan_root, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  tasks, file.path(plan_root, "tasks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
utils::write.table(
  schedule, file.path(output_root, "schedule_index.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
file.copy(strategy_cfg$config_path, file.path(plan_root, "strategy_config.json"), overwrite = TRUE)

cat(sprintf("Planned %d strategy tasks (%d schedules each) under %s\n",
            nrow(tasks), nrow(schedule), output_root))
