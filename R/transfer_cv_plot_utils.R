library(deSolve)
library(dplyr)
library(tidyr)

source("R/project_paths.R")
source("R/plot_utils.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/elpd_transfer_utils.R")

get_stan_line_name <- function(stan_data, line_id) {
  if (!("line_map" %in% names(stan_data)) || is.null(stan_data$line_map) || !length(stan_data$line_map)) {
    return(sprintf("line %d", as.integer(line_id)))
  }

  hit <- names(stan_data$line_map)[as.integer(unname(stan_data$line_map)) == as.integer(line_id)]
  if (!length(hit)) {
    return(sprintf("line %d", as.integer(line_id)))
  }

  hit[[1]]
}

load_transfer_best_fit <- function(
  model_id,
  line_id,
  direction,
  fit_type,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  best_path <- file.path(output_root, model_name, model_id, "transfer_best_start_summary.Rds")
  if (!file.exists(best_path)) {
    stop(sprintf("Best-start summary not found: %s", best_path))
  }

  best_df <- readRDS(best_path)
  direction <- normalize_transfer_direction(direction)
  fit_type <- normalize_fit_type(fit_type)

  hit <- best_df[
    best_df$line_id == line_id &
      best_df$direction == direction &
      best_df$fit_type == fit_type,
  ]

  if (nrow(hit) != 1L) {
    stop(sprintf(
      "Expected exactly one best-start row for model=%s line=%d direction=%s fit=%s",
      model_id,
      line_id,
      direction,
      fit_type
    ))
  }

  run_tag <- hit$run_tag[[1]]
  draws_path <- file.path(
    output_root,
    model_name,
    model_id,
    direction,
    fit_type,
    sprintf("optim_draws_%s.Rds", run_tag)
  )
  meta_path <- file.path(
    output_root,
    model_name,
    model_id,
    direction,
    fit_type,
    sprintf("split_meta_%s.Rds", run_tag)
  )

  if (!file.exists(draws_path)) {
    stop(sprintf("Optim draws not found: %s", draws_path))
  }

  list(
    summary = hit,
    draws = readRDS(draws_path),
    split_meta = readRDS(meta_path),
    run_tag = run_tag
  )
}

make_transfer_obs_df <- function(stan_data, line_id) {
  obs_n <- data.frame(
    well_idx = stan_data$well_idx_count,
    time = stan_data$t_grid[stan_data$grid_idx_count],
    variable = "N",
    value = as.numeric(stan_data$N_obs),
    stringsAsFactors = FALSE
  )

  exps_g <- stan_data$exp_id[stan_data$well_idx_gluc]
  lum <- as.numeric(stan_data$lum_obs)
  dil <- as.numeric(stan_data$dilution)
  obs_g <- data.frame(
    well_idx = stan_data$well_idx_gluc,
    time = stan_data$t_grid[stan_data$grid_idx_gluc],
    variable = "R1",
    value = pmax(0, (lum - stan_data$calib_b_fixed[exps_g]) / (stan_data$calib_a_fixed[exps_g] * dil + 1e-12)),
    stringsAsFactors = FALSE
  )

  obs_df <- bind_rows(obs_n, obs_g)
  obs_df$line_id <- as.integer(stan_data$line_id[obs_df$well_idx])
  obs_df <- obs_df[obs_df$line_id == as.integer(line_id), , drop = FALSE]
  obs_df$line_name <- get_stan_line_name(stan_data, line_id)
  obs_df$ploidy <- as.numeric(stan_data$ploidy_abs[obs_df$well_idx])
  obs_df$G0 <- as.numeric(stan_data$G0_per_well[obs_df$well_idx])

  split_low_high <- get_line_ploidy_states(stan_data, line_id)
  obs_df$group_1 <- ifelse(
    abs(stan_data$ploidy_metric[obs_df$well_idx] - split_low_high$low_value) < split_low_high$tol,
    "low",
    "high"
  )
  obs_df$group_2 <- NA_character_

  obs_df
}

simulate_transfer_fit <- function(
  model_id,
  line_id,
  direction,
  fit_type,
  output_root = "data/gpath_transfer_cv",
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  times = NULL
) {
  dims <- parse_run_id(model_id)
  R <- dims$R
  P <- dims$P
  W <- dims$W
  strict_spec <- (dims$C == 1)
  M <- dims$M
  L <- 3 * R + (P - 1) * R + W * R + 1

  stan_data <- readRDS(resolve_stan_data_path(stan_data_path))
  stan_data <- add_group_structure(stan_data)
  line_name <- get_stan_line_name(stan_data, line_id)

  if (is.null(times)) {
    times <- seq(0, max(stan_data$t_grid), by = 0.5)
  }

  best_fit <- load_transfer_best_fit(
    model_id = model_id,
    line_id = line_id,
    direction = direction,
    fit_type = fit_type,
    output_root = output_root
  )

  draw_vals <- best_fit$draws
  if (is.matrix(draw_vals)) {
    draw_vals <- draw_vals[1, , drop = TRUE]
  }

  line_wells <- which(stan_data$line_id == line_id)

  sim_one_well <- function(well_idx) {
    group_id <- stan_data$group_id[well_idx]
    target_line_id <- stan_data$line_id[well_idx]
    p_metric <- stan_data$ploidy_metric[well_idx]
    if (fit_type == "null") {
      p_metric <- 0
    }

    N0_raw <- draw_vals[sprintf("raw_N0[%d]", group_id)]
    N0 <- exp(log(500) + N0_raw * 1.0)
    G1_0 <- draw_vals[sprintf("G1_0[%d]", stan_data$g1_id[well_idx])]

    raw_theta_line <- draw_vals[sprintf("raw_theta_line[%d,%d]", 1:L, target_line_id)]
    raw_theta_ploidy <- draw_vals[sprintf("raw_theta_ploidy[%d]", 1:L)]

    parms <- reconstruct_parms(
      R = R,
      P = P,
      W = W,
      strict_spec = strict_spec,
      M = M,
      base_priors = base_priors,
      raw_theta_line = raw_theta_line,
      raw_theta_ploidy = raw_theta_ploidy,
      ploidy_metric = p_metric,
      ploidy_effect_mask = best_fit$split_meta$ploidy_effect_mask
    )

    y0 <- c(N = max(as.numeric(N0), 1e-6), R1 = as.numeric(G1_0))
    if (R > 1) {
      y0 <- c(y0, setNames(rep(1.0, R - 1), paste0("R", 2:R)))
    }
    if (W > 0) {
      y0 <- c(y0, setNames(rep(0.0, W), paste0("W", 1:W)))
    }

    out <- as.data.frame(ode(y = y0, times = times, func = rhs, parms = parms, method = "lsoda"))
    out_long <- pivot_longer(out, cols = -time, names_to = "variable", values_to = "value")
    out_long$well_idx <- as.integer(well_idx)
    out_long$line_id <- as.integer(target_line_id)
    out_long$line_name <- line_name
    out_long$ploidy <- as.numeric(stan_data$ploidy_abs[well_idx])
    out_long$G0 <- as.numeric(stan_data$G0_per_well[well_idx])
    out_long$group_1 <- fit_type
    out_long$group_2 <- direction
    out_long$run_tag <- best_fit$run_tag
    out_long
  }

  sim_df <- bind_rows(lapply(line_wells, sim_one_well))
  split_meta <- best_fit$split_meta
  sim_df$heldout_state <- ifelse(
    sim_df$well_idx %in% split_meta$holdout_wells,
    "heldout",
    "observed"
  )

  list(
    sim_df = sim_df,
    obs_df = make_transfer_obs_df(stan_data, line_id),
    summary = best_fit$summary,
    split_meta = split_meta
  )
}

build_transfer_overlay_plot_data <- function(
  model_id,
  line_id,
  directions = c("low_to_high", "high_to_low"),
  fit_types = c("null", "transfer", "oracle"),
  output_root = "data/gpath_transfer_cv",
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  times = NULL,
  include_score_text = TRUE,
  model_name = "gpath"
) {
  runs <- list()
  idx <- 1L

  for (direction in directions) {
    for (fit_type in fit_types) {
      runs[[idx]] <- simulate_transfer_fit(
        model_id = model_id,
        line_id = line_id,
        direction = direction,
        fit_type = fit_type,
        output_root = output_root,
        stan_data_path = stan_data_path,
        times = times
      )
      idx <- idx + 1L
    }
  }

  sim_df <- bind_rows(lapply(runs, `[[`, "sim_df"))
  obs_df <- runs[[1]]$obs_df
  line_name <- unique(obs_df$line_name)
  if (!length(line_name)) {
    line_name <- sprintf("line %d", as.integer(line_id))
  } else {
    line_name <- line_name[[1]]
  }

  score_text <- NULL
  if (isTRUE(include_score_text)) {
    score_text <- tryCatch(
      format_transfer_overlay_scores(
        model_id = model_id,
        line_id = line_id,
        output_root = output_root,
        model_name = model_name
      ),
      error = function(e) NULL
    )
  }

  list(
    sim_df = sim_df,
    obs_df = obs_df,
    context = list(
      model_id = model_id,
      line_id = as.integer(line_id),
      line_name = line_name
    ),
    score_text = score_text
  )
}

generate_transfer_overlay_data <- function(
  model_id,
  line_id,
  directions = c("low_to_high", "high_to_low"),
  fit_types = c("null", "transfer", "oracle"),
  output_root = "data/gpath_transfer_cv",
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  times = NULL
) {
  build_transfer_overlay_plot_data(
    model_id = model_id,
    line_id = line_id,
    directions = directions,
    fit_types = fit_types,
    output_root = output_root,
    stan_data_path = stan_data_path,
    times = times,
    include_score_text = FALSE
  )
}

load_transfer_comparison_summary <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  summary_path <- file.path(output_root, model_name, model_id, "transfer_comparison_summary.Rds")
  if (!file.exists(summary_path)) {
    stop(sprintf("Transfer comparison summary not found: %s", summary_path))
  }

  readRDS(summary_path)
}

format_transfer_overlay_scores <- function(
  model_id,
  line_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath",
  digits = 2
) {
  df <- load_transfer_comparison_summary(
    model_id = model_id,
    output_root = output_root,
    model_name = model_name
  )

  df <- df[df$line_id == line_id, , drop = FALSE]
  if (!nrow(df)) {
    return(NULL)
  }

  df$direction <- factor(df$direction, levels = c("low_to_high", "high_to_low"))
  df <- df[order(df$direction), , drop = FALSE]

  apply(df, 1, function(row) {
    sprintf(
      "%s: null=%.*f | transfer=%.*f | oracle=%.*f | gain=%.*f | regret=%.*f",
      row[["direction"]],
      digits, as.numeric(row[["null"]]),
      digits, as.numeric(row[["transfer"]]),
      digits, as.numeric(row[["oracle"]]),
      digits, as.numeric(row[["transfer_gain"]]),
      digits, as.numeric(row[["transfer_regret"]])
    )
  }) %>%
    paste(collapse = "\n")
}

plot_transfer_line_trajectories <- function(
  transfer_data,
  line_id,
  model_id,
  line_name = NULL,
  score_text = NULL
) {
  if (is.null(line_name) && !is.null(transfer_data$context$line_name)) {
    line_name <- transfer_data$context$line_name
  }
  if (is.null(line_name)) {
    line_name <- sprintf("line %d", as.integer(line_id))
  }
  if (is.null(score_text) && !is.null(transfer_data$score_text)) {
    score_text <- transfer_data$score_text
  }

  plot_fit_overlays(
    sim_df = transfer_data$sim_df,
    obs_df = transfer_data$obs_df,
    color_by = "group_1",
    linetype_by = "group_2",
    color_values = c(null = "#6f6f6f", transfer = "#1b9e77", oracle = "#d95f02"),
    linetype_values = c(low_to_high = "solid", high_to_low = "22"),
    color_label = "Fit",
    linetype_label = "Direction",
    title = sprintf("Transfer CV Fits | Model: %s | Cell Line: %s", model_id, line_name),
    subtitle = score_text
  )
}
