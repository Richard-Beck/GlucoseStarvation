#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

TARGET_LINE_ID <- if (length(args) >= 1) as.integer(args[1]) else stop("line_id required")
OUTPUT_PATH <- if (length(args) >= 2) args[2] else file.path(
  "data", "inputs", "stan", "gstarvation_v1_single_line",
  sprintf("line_%02d", TARGET_LINE_ID),
  "stan_ready_data.Rds"
)
INPUT_PATH <- if (length(args) >= 3) args[3] else file.path("data", "stan_ready_data.Rds")

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")

cat(sprintf(">>> make_gpath_line_stan_data.R cwd: %s\n", getwd()))
cat(sprintf(">>> Input: %s\n", INPUT_PATH))
cat(sprintf(">>> Output: %s\n", OUTPUT_PATH))
cat(sprintf(">>> Target line_id: %d\n", TARGET_LINE_ID))

stan_data_path <- resolve_stan_data_path(INPUT_PATH)
stan_data <- readRDS(stan_data_path)
line_data <- subset_stan_data_to_line(stan_data, TARGET_LINE_ID)

dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
saveRDS(line_data, OUTPUT_PATH)

cat(sprintf(">>> Saved single-line stan_data with %d wells, %d count obs, %d glucose obs.\n",
            line_data$N_wells, line_data$N_obs_count, line_data$N_obs_gluc))
if ("subset_source_line_name" %in% names(line_data)) {
  cat(sprintf(">>> Source line name: %s\n", line_data$subset_source_line_name))
}
