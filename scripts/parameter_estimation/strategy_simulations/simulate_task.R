#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: simulate_task.R RELEASE_CONFIG TASK_ID [MAX_DRAWS] [MAX_SCHEDULES] [OUTPUT_OVERRIDE]",
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(deSolve)
  library(jsonlite)
  library(posterior)
  library(qs)
})
source("scripts/parameter_estimation/common.R")
source("R/project_paths.R")
source("R/gpath_posterior_strategy_utils.R")

cfg <- pe_read_config(args[[1]])
strategy_cfg <- gps_read_strategy_config(cfg)
output_root <- gps_output_root(cfg, strategy_cfg)
tasks <- pe_read_tsv(
  file.path(output_root, "plans", "tasks.tsv"),
  c("task_id", "dataset_id", "model_id", "line_id", "line_name",
    "stan_data_path", "nuts_chains", "nuts_dir", "output_path")
)
task_id <- as.integer(args[[2]])
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) pe_fail("Unknown strategy task: %s", args[[2]])

max_draws <- if (length(args) >= 3L) as.integer(args[[3]]) else as.integer(strategy_cfg$posterior_draws)
max_schedules <- if (length(args) >= 4L) as.integer(args[[4]]) else Inf
output_path <- if (length(args) >= 5L) args[[5]] else task$output_path[[1]]
if (file.exists(output_path)) pe_fail("Refusing to overwrite strategy output: %s", output_path)

source(get_model_r_path(cfg$model_name, cfg$model_version))
source("R/gpath_run_utils.R")
source("R/gpath_derived_utils.R")
source("R/selection_strategy_utils.R")

started <- proc.time()[["elapsed"]]
stan_data <- add_group_structure(readRDS(task$stan_data_path[[1]]))
draw_info <- gpd_read_fit_draws(
  task$nuts_dir[[1]], as.integer(task$nuts_chains[[1]]), max_draws = max_draws
)
strategy_grid <- gps_build_strategy_panel(strategy_cfg)
if (is.finite(max_schedules)) {
  strategy_grid <- strategy_grid[seq_len(min(nrow(strategy_grid), max_schedules)), , drop = FALSE]
}

line_id <- as.integer(task$line_id[[1]])
model_id <- task$model_id[[1]]
state_names <- gps_state_names(model_id)
endpoint_names <- gps_endpoint_names()
time_hours <- seq(
  0,
  as.numeric(strategy_cfg$interval_hours) * as.integer(strategy_cfg$segments),
  by = as.numeric(strategy_cfg$time_step_hours)
)
ndraw <- nrow(draw_info$draws)
nschedule <- nrow(strategy_grid)
ntime <- length(time_hours)
cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
batch_size <- min(cores, ndraw)
compact_mode <- identical(
  as.character(strategy_cfg$storage$mode),
  "compact_observation_quantiles"
)
if (compact_mode && !requireNamespace("matrixStats", quietly = TRUE)) {
  stop("Compact strategy simulation requires matrixStats", call. = FALSE)
}
if (compact_mode) {
  density_multipliers <- as.numeric(strategy_cfg$seed_total_multipliers)
  occasions <- c(
    "day0_post_seed", "day2_pre_action", "day2_post_action",
    "day4_pre_action", "day4_post_action", "day6_end"
  )
  dense_state <- array(
    NA_real_,
    dim = c(ndraw, length(density_multipliers), nschedule, ntime, length(state_names))
  )
  observation_state <- array(
    NA_real_,
    dim = c(ndraw, length(density_multipliers), nschedule, length(occasions), length(state_names))
  )
  fitted_seed_total_n <- numeric(ndraw)
  actual_seed_total_n <- matrix(NA_real_, ndraw, length(density_multipliers))
  draw_elapsed_seconds <- numeric(ndraw)
} else {
  draw_results <- vector("list", ndraw)
}

for (batch_start in seq(1L, ndraw, by = batch_size)) {
  batch <- batch_start:min(ndraw, batch_start + batch_size - 1L)
  results <- if (length(batch) > 1L && cores > 1L) {
    parallel::mclapply(
      batch,
      if (compact_mode) gps_simulate_compact_posterior_draw else gps_simulate_posterior_draw,
      draw_matrix = draw_info$draws,
      stan_data = stan_data,
      model_id = model_id,
      line_id = line_id,
      strategy_grid = strategy_grid,
      strategy_cfg = strategy_cfg,
      mc.cores = min(cores, length(batch)),
      mc.preschedule = TRUE
    )
  } else {
    lapply(
      batch,
      if (compact_mode) gps_simulate_compact_posterior_draw else gps_simulate_posterior_draw,
      draw_matrix = draw_info$draws,
      stan_data = stan_data,
      model_id = model_id,
      line_id = line_id,
      strategy_grid = strategy_grid,
      strategy_cfg = strategy_cfg
    )
  }
  if (any(vapply(results, inherits, logical(1), "try-error"))) {
    stop("One or more posterior-draw simulations failed", call. = FALSE)
  }
  for (result in results) {
    i <- result$draw_index
    if (compact_mode) {
      dense_state[i, , , , ] <- result$trajectory
      observation_state[i, , , , ] <- result$snapshots
      fitted_seed_total_n[[i]] <- result$fitted_seed_total_n
      actual_seed_total_n[i, ] <- result$actual_seed_total_n
      draw_elapsed_seconds[[i]] <- result$elapsed_seconds
    } else {
      draw_results[[i]] <- result
    }
  }
  cat(sprintf("task=%d completed_draws=%d/%d\n", task_id, max(batch), ndraw))
}

common_metadata <- list(
  release_id = cfg$release_id,
  dataset_id = task$dataset_id[[1]],
  model_id = model_id,
  line_id = line_id,
  line_name = task$line_name[[1]],
  stan_data_path = task$stan_data_path[[1]],
  stan_data_sha256 = pe_sha256(task$stan_data_path[[1]]),
  time_step_hours = as.numeric(strategy_cfg$time_step_hours),
  initial_high_fraction = as.numeric(strategy_cfg$initial_high_fraction),
  seed_total_n = "arithmetic mean of draw-specific fitted low/high N0",
  refresh_resets_total_n = isTRUE(strategy_cfg$refresh_resets_total_n),
  refresh_resets_latent_resources = isTRUE(strategy_cfg$refresh_resets_latent_resources),
  refresh_clears_waste = isTRUE(strategy_cfg$refresh_clears_waste),
  elapsed_seconds = proc.time()[["elapsed"]] - started,
  generated_at = format(Sys.time(), usetz = TRUE)
)

if (compact_mode) {
  probabilities <- as.numeric(strategy_cfg$storage$quantiles)
  dense_matrix <- matrix(dense_state, nrow = ndraw)
  trajectory_quantiles <- matrixStats::colQuantiles(
    dense_matrix, probs = probabilities, na.rm = FALSE, drop = FALSE
  )
  if (nrow(trajectory_quantiles) != length(probabilities)) {
    trajectory_quantiles <- t(trajectory_quantiles)
  }
  trajectory_quantiles <- array(
    trajectory_quantiles,
    dim = c(
      length(probabilities), length(density_multipliers), nschedule,
      ntime, length(state_names)
    ),
    dimnames = list(
      quantile = paste0("q", probabilities * 100),
      density_multiplier = as.character(density_multipliers),
      schedule = strategy_grid$strategy_code,
      time_hours = as.character(time_hours),
      state = state_names
    )
  )
  dimnames(observation_state) <- list(
    draw = draw_info$index$draw_id,
    density_multiplier = as.character(density_multipliers),
    schedule = strategy_grid$strategy_code,
    occasion = occasions,
    state = state_names
  )
  dimnames(actual_seed_total_n) <- list(
    draw = draw_info$index$draw_id,
    density_multiplier = as.character(density_multipliers)
  )
  artifact <- list(
    schema_version = 2L,
    metadata = c(common_metadata, list(
      storage_mode = "all_draws_at_observation_times_plus_dense_pointwise_quantiles",
      topup_definition = "add X mM to existing G1; no medium replacement or other state reset"
    )),
    draw_index = draw_info$index,
    schedule_index = strategy_grid,
    density_multipliers = density_multipliers,
    observation_occasions = occasions,
    time_hours = time_hours,
    state_names = state_names,
    observation_state = observation_state,
    trajectory_quantile_probabilities = probabilities,
    trajectory_quantiles = trajectory_quantiles,
    fitted_seed_total_n = fitted_seed_total_n,
    actual_seed_total_n = actual_seed_total_n,
    draw_elapsed_seconds = draw_elapsed_seconds
  )
  rm(dense_state, dense_matrix)
} else {
  endpoint <- array(
    NA_real_,
    dim = c(ndraw, nschedule, length(endpoint_names)),
    dimnames = list(draw = draw_info$index$draw_id,
                    schedule = strategy_grid$strategy_code,
                    metric = endpoint_names)
  )
  for (i in seq_len(ndraw)) endpoint[i, , ] <- draw_results[[i]]$endpoints
  artifact <- list(
    schema_version = 1L,
    metadata = common_metadata,
    draw_index = draw_info$index,
    schedule_index = strategy_grid,
    time_hours = time_hours,
    state_names = state_names,
    state = lapply(draw_results, `[[`, "trajectory"),
    endpoint = endpoint,
    seed_total_n = vapply(draw_results, `[[`, numeric(1), "seed_total_n"),
    draw_elapsed_seconds = vapply(draw_results, `[[`, numeric(1), "elapsed_seconds")
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
tmp_path <- tempfile(pattern = paste0(".", basename(output_path), "."), tmpdir = dirname(output_path))
on.exit(unlink(tmp_path), add = TRUE)
qs::qsave(
  artifact, tmp_path,
  preset = as.character(strategy_cfg$storage$preset),
  algorithm = as.character(strategy_cfg$storage$algorithm),
  nthreads = min(cores, 8L)
)
check <- qs::qread(tmp_path, nthreads = min(cores, 8L))
if (compact_mode) {
  expected_obs <- c(
    ndraw, length(density_multipliers), nschedule,
    length(occasions), length(state_names)
  )
  if (!identical(dim(check$observation_state), expected_obs)) {
    stop("Serialized compact artifact failed immediate readback validation", call. = FALSE)
  }
} else if (length(check$state) != ndraw ||
           !identical(dim(check$endpoint), c(ndraw, nschedule, length(endpoint_names)))) {
  stop("Serialized strategy artifact failed immediate readback validation", call. = FALSE)
}
rm(check)
if (!file.rename(tmp_path, output_path)) pe_fail("Could not install %s", output_path)
cat(sprintf("Wrote %s (%d draws x %d densities x %d schedules x %d times)\n",
            output_path, ndraw,
            if (compact_mode) length(density_multipliers) else 1L,
            nschedule, ntime))
