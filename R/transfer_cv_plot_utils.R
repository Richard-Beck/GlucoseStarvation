library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/elpd_transfer_utils.R")

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
  obsN <- data.frame(
    well_idx = stan_data$well_idx_count,
    time = stan_data$t_grid[stan_data$grid_idx_count],
    Value = as.numeric(stan_data$N_obs),
    Variable = "N"
  )

  exps_G <- stan_data$exp_id[stan_data$well_idx_gluc]
  lum <- as.numeric(stan_data$lum_obs)
  dil <- as.numeric(stan_data$dilution)
  obsG <- data.frame(
    well_idx = stan_data$well_idx_gluc,
    time = stan_data$t_grid[stan_data$grid_idx_gluc],
    Value = pmax(0, (lum - stan_data$calib_b_fixed[exps_G]) / (stan_data$calib_a_fixed[exps_G] * dil + 1e-12)),
    Variable = "R1"
  )

  obs_df <- bind_rows(obsN, obsG)
  obs_df$ploidy_abs <- stan_data$ploidy_abs[obs_df$well_idx]
  obs_df$G0 <- stan_data$G0_per_well[obs_df$well_idx]
  obs_df$line_id <- stan_data$line_id[obs_df$well_idx]
  obs_df <- obs_df[obs_df$line_id == line_id, ]

  split_low_high <- get_line_ploidy_states(stan_data, line_id)
  obs_df$ploidy_role <- ifelse(
    abs(stan_data$ploidy_metric[obs_df$well_idx] - split_low_high$low_value) < split_low_high$tol,
    "low",
    "high"
  )

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
      ploidy_metric = p_metric
    )

    y0 <- c(N = max(as.numeric(N0), 1e-6), R1 = as.numeric(G1_0))
    if (R > 1) {
      y0 <- c(y0, setNames(rep(1.0, R - 1), paste0("R", 2:R)))
    }
    if (W > 0) {
      y0 <- c(y0, setNames(rep(0.0, W), paste0("W", 1:W)))
    }

    out <- as.data.frame(ode(y = y0, times = times, func = rhs, parms = parms, method = "lsoda"))
    out_long <- pivot_longer(out, cols = -time, names_to = "Variable", values_to = "Value")
    out_long$well_idx <- well_idx
    out_long$ploidy_abs <- stan_data$ploidy_abs[well_idx]
    out_long$G0 <- stan_data$G0_per_well[well_idx]
    out_long$line_id <- target_line_id
    out_long$fit_type <- fit_type
    out_long$direction <- direction
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

generate_transfer_overlay_data <- function(
  model_id,
  line_id,
  directions = c("low_to_high", "high_to_low"),
  fit_types = c("null", "transfer", "oracle"),
  output_root = "data/gpath_transfer_cv",
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  times = NULL
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

  list(
    sim_df = sim_df,
    obs_df = obs_df
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
  if (is.null(line_name)) {
    stan_data <- readRDS(resolve_stan_data_path("data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"))
    line_name <- if (!is.null(stan_data$line_map)) names(stan_data$line_map)[line_id] else sprintf("line %d", line_id)
  }

  sim_df <- transfer_data$sim_df
  obs_df <- transfer_data$obs_df

  g0_levels <- sort(unique(sim_df$G0))
  g0_labels <- paste0("G[0]=", g0_levels)
  sim_df$G0_label <- factor(paste0("G[0]=", sim_df$G0), levels = g0_labels)
  obs_df$G0_label <- factor(paste0("G[0]=", obs_df$G0), levels = g0_labels)

  sim_df$fit_type <- factor(sim_df$fit_type, levels = c("null", "transfer", "oracle"))
  sim_df$direction <- factor(sim_df$direction, levels = c("low_to_high", "high_to_low"))

  vars_to_plot <- unique(sim_df$Variable)

  plot_list <- lapply(vars_to_plot, function(v) {
    sub_mean <- sim_df[sim_df$Variable == v, ]
    sub_obs <- obs_df[obs_df$Variable == v, ]

    title_str <- v
    if (v == "N") {
      title_str <- "N (alive cells)"
    }
    if (v == "R1") {
      title_str <- "R1: Glucose"
    }

    ggplot() +
      geom_line(
        data = sub_mean,
        aes(x = time, y = Value, color = fit_type, linetype = direction, group = interaction(well_idx, fit_type, direction)),
        linewidth = 0.9,
        alpha = 0.85
      ) +
      geom_point(
        data = sub_obs,
        aes(x = time, y = Value),
        color = "black",
        size = 1.1,
        alpha = 0.75
      ) +
      facet_grid(rows = vars(G0_label), cols = vars(ploidy_abs), scales = "free") +
      scale_color_manual(values = c(null = "#6f6f6f", transfer = "#1b9e77", oracle = "#d95f02")) +
      scale_linetype_manual(values = c(low_to_high = "solid", high_to_low = "22")) +
      theme_minimal() +
      labs(title = title_str, x = "Time", y = "", color = "Fit", linetype = "Direction") +
      theme(
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
  })

  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = sprintf("Transfer CV Fits | Model: %s | Cell Line: %s", model_id, line_name),
      subtitle = score_text,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  print(combined_plot)
  invisible(combined_plot)
}
