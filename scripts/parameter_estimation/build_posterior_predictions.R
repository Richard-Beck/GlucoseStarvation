#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: build_posterior_predictions.R CONFIG TASK_ID", call. = FALSE)
suppressPackageStartupMessages({library(posterior); library(deSolve)})
source("scripts/parameter_estimation/common.R")
source("R/project_paths.R")

cfg <- pe_read_config(args[[1]])
derived_cfg <- pe_read_derived_config(cfg)
task_id <- as.integer(args[[2]])
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
tasks <- pe_read_tsv(file.path(derived_root, "plans", "prediction_tasks.tsv"), c(
  "task_id", "fit_id", "dataset_id", "model_id", "stan_data_path", "nuts_chains", "nuts_dir", "output_path"
))
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) pe_fail("Unknown prediction task id: %s", args[[2]])
if (file.exists(task$output_path[[1]])) pe_fail("Refusing to overwrite prediction shard: %s", task$output_path[[1]])

source(get_model_r_path(cfg$model_name, cfg$model_version))
source("R/gpath_run_utils.R")
source("R/gpath_derived_utils.R")

stan_data <- add_group_structure(readRDS(task$stan_data_path[[1]]))
max_draws <- as.integer(derived_cfg$prediction_draws)
surface_draws <- as.integer(derived_cfg$surface_draws)
draw_info <- gpd_read_fit_draws(task$nuts_dir[[1]], task$nuts_chains[[1]], max_draws = max_draws)
times <- if (identical(as.character(derived_cfg$prediction_time_grid), "stan")) as.numeric(stan_data$t_grid) else {
  stop("Only prediction_time_grid='stan' is currently supported", call. = FALSE)
}
cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

simulate_one <- function(i) {
  draw_vec <- draw_info$draws[i, , drop = TRUE]
  gpd_simulate_draw(draw_vec, stan_data, task$model_id[[1]], times)
}
state_parts <- if (cores > 1L) parallel::mclapply(seq_len(nrow(draw_info$draws)), simulate_one, mc.cores = cores) else lapply(seq_len(nrow(draw_info$draws)), simulate_one)
states <- gpd_bind_first_dimension(state_parts)
dimnames(states)[[1]] <- as.character(draw_info$index$draw_id)
names(dimnames(states))[1] <- "draw"

surface_cfg <- derived_cfg$growth_surface
glucose_grid <- exp(seq(log(as.numeric(surface_cfg$glucose_min)), log(as.numeric(surface_cfg$glucose_max)), length.out = as.integer(surface_cfg$glucose_points)))
resource2_grid <- seq(as.numeric(surface_cfg$resource2_min), as.numeric(surface_cfg$resource2_max), length.out = as.integer(surface_cfg$resource2_points))
surface_idx <- unique(as.integer(round(seq(1, nrow(draw_info$draws), length.out = min(surface_draws, nrow(draw_info$draws))))))
surface_one <- function(i) {
  draw_vec <- draw_info$draws[i, , drop = TRUE]
  gpd_growth_delta_surface(
    draw_vec, stan_data, task$model_id[[1]], glucose_grid, resource2_grid,
    extra_resource_value = as.numeric(surface_cfg$extra_resource_value)
  )
}
surface_parts <- if (cores > 1L) parallel::mclapply(surface_idx, surface_one, mc.cores = cores) else lapply(surface_idx, surface_one)
growth_delta <- gpd_bind_first_dimension(surface_parts)
dimnames(growth_delta)[[1]] <- as.character(draw_info$index$draw_id[surface_idx])
names(dimnames(growth_delta))[1] <- "draw"

lum_mean <- matrix(NA_real_, nrow = dim(states)[1], ncol = stan_data$N_obs_gluc)
r1_idx <- match("R1", dimnames(states)$state)
for (obs in seq_len(stan_data$N_obs_gluc)) {
  well <- stan_data$well_idx_gluc[[obs]]
  time <- stan_data$grid_idx_gluc[[obs]]
  exp_id <- stan_data$exp_id[[well]]
  lum_mean[, obs] <- stan_data$calib_a_fixed[[exp_id]] * states[, well, time, r1_idx] * stan_data$dilution[[obs]] + stan_data$calib_b_fixed[[exp_id]]
}

artifact <- list(
  schema_version = 1L,
  metadata = list(
    release_id = cfg$release_id,
    fit_id = task$fit_id[[1]], dataset_id = task$dataset_id[[1]], model_id = task$model_id[[1]],
    stan_data_path = task$stan_data_path[[1]], stan_data_sha256 = pe_sha256(task$stan_data_path[[1]]),
    source_chain_sha256 = setNames(lapply(draw_info$chain_paths, pe_sha256), draw_info$chain_paths),
    draw_policy = sprintf("%d chain-balanced evenly spaced sampling draws", nrow(draw_info$draws)),
    generated_at = format(Sys.time(), usetz = TRUE)
  ),
  draw_index = draw_info$index,
  well_metadata = gpd_well_metadata(stan_data),
  time = times,
  state = states,
  glucose_observation = list(
    well_idx = stan_data$well_idx_gluc,
    grid_idx = stan_data$grid_idx_gluc,
    dilution = stan_data$dilution,
    expected_luminescence = lum_mean
  ),
  growth_surface = list(
    draw_index = draw_info$index[surface_idx, , drop = FALSE],
    glucose = glucose_grid,
    resource2 = resource2_grid,
    high_minus_low_growth = growth_delta
  )
)
gpd_atomic_save_rds(artifact, task$output_path[[1]], compress = "gzip")
cat(sprintf("Wrote %s (%d trajectory draws; %d surface draws)\n",
            task$output_path[[1]], nrow(draw_info$draws), length(surface_idx)))
