gps_strategy_config_path <- function(release_cfg, explicit_path = NULL) {
  candidate <- explicit_path
  if (is.null(candidate) || !nzchar(candidate)) {
    candidate <- Sys.getenv("GPATH_STRATEGY_CONFIG", unset = "")
  }
  if (!nzchar(candidate)) {
    candidate <- file.path(dirname(release_cfg$config_path), "strategy_simulations.json")
  }
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

gps_read_strategy_config <- function(release_cfg, explicit_path = NULL) {
  path <- gps_strategy_config_path(release_cfg, explicit_path)
  if (!file.exists(path)) stop("Missing strategy-simulation config: ", path, call. = FALSE)
  out <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  out$config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out
}

gps_output_root <- function(release_cfg, strategy_cfg = NULL) {
  root <- file.path(pe_derived_root(release_cfg), "posterior", "strategy_simulations")
  subdir <- if (is.null(strategy_cfg$output_subdir)) "" else as.character(strategy_cfg$output_subdir)
  if (!nzchar(subdir)) return(root)
  if (grepl("^/", subdir) || any(strsplit(subdir, "/", fixed = TRUE)[[1]] == "..")) {
    stop("strategy output_subdir must be a safe relative path", call. = FALSE)
  }
  file.path(root, subdir)
}

gps_slug <- function(x) {
  out <- tolower(gsub("[^A-Za-z0-9]+", "_", as.character(x)))
  gsub("^_+|_+$", "", out)
}

gps_prior_n0_center <- function(stan_data) {
  center <- as.numeric(stan_data$prior_N0_center)
  if (length(center) != 1L || !is.finite(center) || center <= 0) {
    stop("stan_data$prior_N0_center must be positive and finite", call. = FALSE)
  }
  center
}

gps_estimate_state_n0 <- function(draw_vec, stan_data, line_id, ploidy_metric, tol = 1e-12) {
  group_id <- get_group_id_for_state(
    stan_data, line_id = line_id, ploidy_metric = ploidy_metric, tol = tol
  )
  raw_name <- sprintf("raw_N0[%d]", group_id)
  raw_val <- as.numeric(draw_vec[[raw_name]])
  if (length(raw_val) != 1L || !is.finite(raw_val)) {
    stop("Posterior draw lacks finite ", raw_name, call. = FALSE)
  }
  exp(log(gps_prior_n0_center(stan_data)) + raw_val)
}

gps_estimate_mixture_seed_total_n <- function(
  draw_vec, stan_data, line_id, low_ploidy, high_ploidy
) {
  mean(c(
    gps_estimate_state_n0(draw_vec, stan_data, line_id, low_ploidy),
    gps_estimate_state_n0(draw_vec, stan_data, line_id, high_ploidy)
  ))
}

gps_state_names <- function(model_id) {
  dims <- parse_run_id(model_id)
  c(
    "N_low", "N_high", "G1",
    if (dims$R > 1L) sprintf("R%d", 2:dims$R),
    if (dims$W > 0L) sprintf("W%d", seq_len(dims$W))
  )
}

gps_endpoint_names <- function() {
  c(
    "N_low", "N_high", "total_live", "high_fraction",
    "log_ratio_high_low", "low_fold_change", "high_fold_change"
  )
}

gps_build_strategy_panel <- function(strategy_cfg) {
  design <- if (is.null(strategy_cfg$schedule_design)) {
    "legacy_full_grid"
  } else {
    as.character(strategy_cfg$schedule_design)
  }
  glucose <- as.numeric(strategy_cfg$glucose_mM)
  if (identical(design, "legacy_full_grid")) {
    return(build_strategy_grid(glucose))
  }
  if (!identical(design, "simple_xxx_xcc_addx")) {
    stop("Unknown strategy schedule_design: ", design, call. = FALSE)
  }

  glucose_label <- function(x) format(x, scientific = FALSE, trim = TRUE)
  rows <- lapply(seq_along(glucose), function(i) {
    x <- glucose[[i]]
    label <- glucose_label(x)
    data.frame(
      day0_idx = i,
      day2_idx = i,
      day4_idx = i,
      strategy_code = c(
        sprintf("%s,%s,%s", label, label, label),
        sprintf("%s,C,C", label),
        sprintf("%s,A%s,A%s", label, label, label)
      ),
      g0_day0 = x,
      g0_day2 = x,
      g0_day4 = x,
      day2_action = c("refresh", "carry", "add_glucose"),
      day4_action = c("refresh", "carry", "add_glucose"),
      strategy_family = c("XXX", "XCC", "XA_xA_x"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gps_simulate_compact_posterior_draw <- function(
  draw_index,
  draw_matrix,
  stan_data,
  model_id,
  line_id,
  strategy_grid,
  strategy_cfg
) {
  started <- proc.time()[["elapsed"]]
  draw_vec <- draw_matrix[draw_index, , drop = TRUE]
  states <- get_line_ploidy_states(stan_data, line_id)
  comp <- reconstruct_mixture_components(
    draw_vec = draw_vec,
    model_id = model_id,
    line_id = line_id,
    low_ploidy = states$low_value,
    high_ploidy = states$high_value
  )
  fitted_seed_total <- gps_estimate_mixture_seed_total_n(
    draw_vec, stan_data, line_id, states$low_value, states$high_value
  )
  multipliers <- as.numeric(strategy_cfg$seed_total_multipliers)
  state_names <- gps_state_names(model_id)
  time_hours <- seq(
    0,
    as.numeric(strategy_cfg$interval_hours) * as.integer(strategy_cfg$segments),
    by = as.numeric(strategy_cfg$time_step_hours)
  )
  occasions <- c(
    "day0_post_seed", "day2_pre_action", "day2_post_action",
    "day4_pre_action", "day4_post_action", "day6_end"
  )
  trajectory <- array(
    NA_real_,
    dim = c(length(multipliers), nrow(strategy_grid), length(time_hours), length(state_names))
  )
  snapshots <- array(
    NA_real_,
    dim = c(length(multipliers), nrow(strategy_grid), length(occasions), length(state_names))
  )

  for (density_i in seq_along(multipliers)) {
    actual_seed <- fitted_seed_total * multipliers[[density_i]]
    for (strategy_i in seq_len(nrow(strategy_grid))) {
      sim <- simulate_strategy(
        strategy_row = strategy_grid[strategy_i, , drop = FALSE],
        parms_mix = comp$parms_mix,
        seed_total_n = actual_seed,
        initial_high_fraction = as.numeric(strategy_cfg$initial_high_fraction),
        interval_hours = as.numeric(strategy_cfg$interval_hours),
        time_step_hours = as.numeric(strategy_cfg$time_step_hours),
        reset_total_n_on_refresh = isTRUE(strategy_cfg$refresh_resets_total_n),
        detailed = TRUE
      )
      detail <- sim$detail
      if (nrow(detail) != length(time_hours) ||
          max(abs(as.numeric(detail$time_global) - time_hours)) > 1e-8) {
        stop("Unexpected time grid for draw ", draw_index, call. = FALSE)
      }
      values <- as.matrix(detail[, state_names, drop = FALSE])
      transition <- as.matrix(sim$transition_states[, state_names, drop = FALSE])
      if (any(!is.finite(values)) || any(!is.finite(transition))) {
        stop("Non-finite compact trajectory for draw ", draw_index, call. = FALSE)
      }
      trajectory[density_i, strategy_i, , ] <- values
      snapshots[density_i, strategy_i, , ] <- transition
    }
  }

  list(
    draw_index = draw_index,
    trajectory = trajectory,
    snapshots = snapshots,
    fitted_seed_total_n = fitted_seed_total,
    actual_seed_total_n = fitted_seed_total * multipliers,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}

gps_simulate_posterior_draw <- function(
  draw_index,
  draw_matrix,
  stan_data,
  model_id,
  line_id,
  strategy_grid,
  strategy_cfg
) {
  started <- proc.time()[["elapsed"]]
  draw_vec <- draw_matrix[draw_index, , drop = TRUE]
  states <- get_line_ploidy_states(stan_data, line_id)
  comp <- reconstruct_mixture_components(
    draw_vec = draw_vec,
    model_id = model_id,
    line_id = line_id,
    low_ploidy = states$low_value,
    high_ploidy = states$high_value
  )
  seed_total_n <- gps_estimate_mixture_seed_total_n(
    draw_vec, stan_data, line_id, states$low_value, states$high_value
  )

  state_names <- gps_state_names(model_id)
  expected_time <- seq(
    0,
    as.numeric(strategy_cfg$interval_hours) * as.integer(strategy_cfg$segments),
    by = as.numeric(strategy_cfg$time_step_hours)
  )
  trajectory <- array(
    NA_real_,
    dim = c(nrow(strategy_grid), length(expected_time), length(state_names))
  )
  endpoints <- matrix(
    NA_real_,
    nrow = nrow(strategy_grid),
    ncol = length(gps_endpoint_names()),
    dimnames = list(NULL, gps_endpoint_names())
  )

  for (strategy_i in seq_len(nrow(strategy_grid))) {
    sim <- simulate_strategy(
      strategy_row = strategy_grid[strategy_i, , drop = FALSE],
      parms_mix = comp$parms_mix,
      seed_total_n = seed_total_n,
      initial_high_fraction = as.numeric(strategy_cfg$initial_high_fraction),
      interval_hours = as.numeric(strategy_cfg$interval_hours),
      time_step_hours = as.numeric(strategy_cfg$time_step_hours),
      reset_total_n_on_refresh = isTRUE(strategy_cfg$refresh_resets_total_n),
      detailed = TRUE
    )
    detail <- sim$detail
    if (nrow(detail) != length(expected_time) ||
        max(abs(as.numeric(detail$time_global) - expected_time)) > 1e-8) {
      stop("Unexpected detailed time grid for draw ", draw_index,
           ", schedule ", strategy_grid$strategy_code[[strategy_i]], call. = FALSE)
    }
    values <- as.matrix(detail[, state_names, drop = FALSE])
    if (any(!is.finite(values))) {
      stop("Non-finite trajectory for draw ", draw_index,
           ", schedule ", strategy_grid$strategy_code[[strategy_i]], call. = FALSE)
    }
    trajectory[strategy_i, , ] <- values
    final <- sim$final_summary[1, , drop = FALSE]
    endpoints[strategy_i, ] <- c(
      final$N_low, final$N_high, final$total_live, final$high_fraction_final,
      final$log_ratio_final, final$low_fold_change, final$high_fold_change
    )
  }

  list(
    draw_index = draw_index,
    trajectory = trajectory,
    endpoints = endpoints,
    seed_total_n = seed_total_n,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}
