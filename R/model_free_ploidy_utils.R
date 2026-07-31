source("R/project_paths.R")
source("R/gpath_run_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

g0_display_levels <- function() {
  c(0, 0.1, 0.25, 0.5, 1, 5, 25)
}

as_g0_factor <- function(x) {
  factor(as.character(x), levels = as.character(g0_display_levels()))
}

nearest_value_at_or_after <- function(hours, values, target_time) {
  ok <- is.finite(hours) & is.finite(values)
  hours <- hours[ok]
  values <- values[ok]

  if (!length(hours)) {
    return(NA_real_)
  }

  ord <- order(hours)
  hours <- hours[ord]
  values <- values[ord]

  idx <- which(hours >= target_time)
  if (length(idx)) {
    return(values[min(idx)])
  }

  values[length(values)]
}

safe_trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

calc_rate_extrema <- function(hours, value, log1p_scale = FALSE) {
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (length(hours) < 2L) {
    return(list(
      max_rate = NA_real_,
      min_rate = NA_real_,
      time_at_max_rate = NA_real_,
      time_at_min_rate = NA_real_
    ))
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  x <- if (log1p_scale) log1p(pmax(value, 0)) else value
  dt <- diff(hours)
  dy <- diff(x)

  keep <- is.finite(dt) & dt > 0 & is.finite(dy)
  dt <- dt[keep]
  dy <- dy[keep]

  if (!length(dt)) {
    return(list(
      max_rate = NA_real_,
      min_rate = NA_real_,
      time_at_max_rate = NA_real_,
      time_at_min_rate = NA_real_
    ))
  }

  rate <- dy / dt
  mids <- (hours[-1L] + hours[-length(hours)]) / 2
  mids <- mids[keep]

  idx_max <- which.max(rate)
  idx_min <- which.min(rate)

  list(
    max_rate = rate[idx_max],
    min_rate = rate[idx_min],
    time_at_max_rate = mids[idx_max],
    time_at_min_rate = mids[idx_min]
  )
}

calc_smoothed_rate_extrema <- function(hours, value, log1p_scale = FALSE, spar = 0.55, n_eval = 200L) {
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (length(hours) < 4L || length(unique(hours)) < 4L) {
    return(calc_rate_extrema(hours, value, log1p_scale = log1p_scale))
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]
  y <- if (log1p_scale) log1p(pmax(value, 0)) else value

  fit <- try(stats::smooth.spline(x = hours, y = y, spar = spar), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(calc_rate_extrema(hours, value, log1p_scale = log1p_scale))
  }

  eval_hours <- seq(min(hours), max(hours), length.out = n_eval)
  deriv <- predict(fit, x = eval_hours, deriv = 1)$y

  if (!length(deriv) || all(!is.finite(deriv))) {
    return(calc_rate_extrema(hours, value, log1p_scale = log1p_scale))
  }

  idx_max <- which.max(deriv)
  idx_min <- which.min(deriv)

  list(
    max_rate = deriv[idx_max],
    min_rate = deriv[idx_min],
    time_at_max_rate = eval_hours[idx_max],
    time_at_min_rate = eval_hours[idx_min]
  )
}

fit_exponential_phase <- function(hours, live_cells, frac_of_max = 0.7) {
  ok <- is.finite(hours) & is.finite(live_cells)
  hours <- hours[ok]
  live_cells <- live_cells[ok]

  if (length(hours) < 3L) {
    return(list(
      rate = NA_real_,
      intercept = NA_real_,
      end_time = NA_real_,
      n_fit_points = 0L,
      fit_r_squared = NA_real_
    ))
  }

  ord <- order(hours)
  hours <- hours[ord]
  live_cells <- live_cells[ord]

  max_live <- max(live_cells, na.rm = TRUE)
  threshold <- frac_of_max * max_live
  hit_idx <- which(live_cells >= threshold)
  if (!length(hit_idx)) {
    return(list(
      rate = NA_real_,
      intercept = NA_real_,
      end_time = NA_real_,
      n_fit_points = 0L,
      fit_r_squared = NA_real_
    ))
  }

  end_idx <- min(hit_idx)
  fit_hours <- hours[seq_len(end_idx)]
  fit_live <- live_cells[seq_len(end_idx)]

  if (length(fit_hours) < 3L || length(unique(fit_hours)) < 2L) {
    return(list(
      rate = NA_real_,
      intercept = NA_real_,
      end_time = hours[end_idx],
      n_fit_points = length(fit_hours),
      fit_r_squared = NA_real_
    ))
  }

  fit <- lm(log1p(fit_live) ~ fit_hours)
  coef_fit <- coef(fit)

  list(
    rate = unname(coef_fit[["fit_hours"]]),
    intercept = unname(coef_fit[["(Intercept)"]]),
    end_time = hours[end_idx],
    n_fit_points = length(fit_hours),
    fit_r_squared = summary(fit)$r.squared
  )
}

time_to_threshold <- function(hours, value, threshold, direction = c("below", "above")) {
  direction <- match.arg(direction)
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (!length(hours)) {
    return(NA_real_)
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  hit_idx <- if (direction == "below") {
    which(value <= threshold)
  } else {
    which(value >= threshold)
  }

  if (!length(hit_idx)) {
    return(NA_real_)
  }

  hours[min(hit_idx)]
}

build_model_free_tables <- function(stan_data_path = NULL) {
  stan_data_path <- resolve_stan_data_path(stan_data_path)
  stan_data <- readRDS(stan_data_path)

  line_names <- names(stan_data$line_map)[match(stan_data$line_id, unname(stan_data$line_map))]

  well_meta <- tibble(
    well_idx = seq_len(stan_data$N_wells),
    line_id = as.integer(stan_data$line_id),
    cellLine = line_names,
    ploidy_metric = as.numeric(stan_data$ploidy_metric),
    ploidy_abs = as.numeric(stan_data$ploidy_abs),
    G0 = as.numeric(stan_data$G0_per_well),
    exp_id = as.integer(stan_data$exp_id),
    has_starvation = as.integer(stan_data$has_starvation)
  )

  count_obs <- tibble(
    well_idx = as.integer(stan_data$well_idx_count),
    hours = as.numeric(stan_data$t_grid[stan_data$grid_idx_count]),
    rep_id = as.character(stan_data$rep_id_count),
    live_cells = as.numeric(stan_data$N_obs),
    dead_cells = as.numeric(stan_data$D_obs)
  ) %>%
    mutate(
      total_cells = live_cells + dead_cells,
      dead_fraction = if_else(total_cells > 0, dead_cells / total_cells, NA_real_)
    )

  count_summary <- count_obs %>%
    group_by(well_idx, hours) %>%
    summarise(
      n_count_reps = n(),
      live_cells = median(live_cells, na.rm = TRUE),
      total_cells = median(total_cells, na.rm = TRUE),
      dead_cells = median(dead_cells, na.rm = TRUE),
      dead_fraction = median(dead_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(well_meta, by = "well_idx")

  gluc_obs <- tibble(
    well_idx = as.integer(stan_data$well_idx_gluc),
    hours = as.numeric(stan_data$t_grid[stan_data$grid_idx_gluc]),
    lum_obs = as.numeric(stan_data$lum_obs),
    dilution = as.numeric(stan_data$dilution),
    is_censored = as.integer(stan_data$is_censored)
  ) %>%
    mutate(
      exp_id = well_meta$exp_id[well_idx],
      calib_a = as.numeric(stan_data$calib_a_fixed[exp_id]),
      calib_b = as.numeric(stan_data$calib_b_fixed[exp_id]),
      glucose_hat = pmax(0, (lum_obs - calib_b) / (calib_a * pmax(dilution, 1e-12)))
    )

  gluc_summary <- gluc_obs %>%
    group_by(well_idx, hours) %>%
    summarise(
      n_glucose_reps = n(),
      glucose_hat = {
        uncens <- glucose_hat[is_censored == 0L]
        if (length(uncens)) {
          median(uncens, na.rm = TRUE)
        } else {
          median(glucose_hat, na.rm = TRUE)
        }
      },
      glucose_hat_lo = quantile(glucose_hat, probs = 0.25, na.rm = TRUE, names = FALSE),
      glucose_hat_hi = quantile(glucose_hat, probs = 0.75, na.rm = TRUE, names = FALSE),
      any_glucose_censored = as.integer(any(is_censored == 1L)),
      .groups = "drop"
    ) %>%
    left_join(well_meta, by = "well_idx")

  list(
    stan_data_path = stan_data_path,
    stan_data = stan_data,
    well_meta = well_meta,
    count_obs = count_obs,
    count_summary = count_summary,
    glucose_obs = gluc_obs,
    glucose_summary = gluc_summary
  )
}

summarize_one_well_features <- function(well_meta_row, count_df, glucose_df, glucose_floor = 0.1, min_glucose_drawdown_for_yield = 0.05) {
  count_df <- count_df %>% arrange(hours)
  glucose_df <- glucose_df %>% arrange(hours)

  live_rate_raw <- calc_rate_extrema(count_df$hours, count_df$live_cells, log1p_scale = TRUE)
  dead_rate_raw <- calc_rate_extrema(count_df$hours, count_df$dead_cells, log1p_scale = TRUE)
  live_rate <- calc_smoothed_rate_extrema(count_df$hours, count_df$live_cells, log1p_scale = TRUE)
  dead_rate <- calc_smoothed_rate_extrema(count_df$hours, count_df$dead_cells, log1p_scale = TRUE)
  exp_phase_fit <- fit_exponential_phase(count_df$hours, count_df$live_cells, frac_of_max = 0.7)
  dead_frac_rate <- calc_rate_extrema(count_df$hours, count_df$dead_fraction, log1p_scale = FALSE)
  glucose_rate <- calc_rate_extrema(glucose_df$hours, glucose_df$glucose_hat, log1p_scale = FALSE)

  peak_live_idx <- if (nrow(count_df)) which.max(count_df$live_cells) else integer(0)
  nadir_live_idx <- if (nrow(count_df)) which.min(count_df$live_cells) else integer(0)
  peak_dead_idx <- if (nrow(count_df)) which.max(count_df$dead_cells) else integer(0)
  min_glucose_idx <- if (nrow(glucose_df)) which.min(glucose_df$glucose_hat) else integer(0)
  glucose_start_time <- if (nrow(glucose_df)) min(glucose_df$hours, na.rm = TRUE) else NA_real_
  glucose_end_time <- if (nrow(glucose_df)) max(glucose_df$hours, na.rm = TRUE) else NA_real_
  count_glucose_window <- if (is.finite(glucose_start_time) && is.finite(glucose_end_time)) {
    count_df %>% filter(hours >= glucose_start_time, hours <= glucose_end_time)
  } else {
    count_df[0, ]
  }
  total_peak_to_glucose_end <- if (nrow(count_glucose_window)) {
    max(count_glucose_window$total_cells, na.rm = TRUE)
  } else {
    NA_real_
  }
  total_initial_at_glucose_start <- if (nrow(count_glucose_window)) {
    nearest_value_at_or_after(count_glucose_window$hours, count_glucose_window$total_cells, glucose_start_time)
  } else {
    NA_real_
  }
  total_net_gain_to_glucose_end <- if (is.finite(total_peak_to_glucose_end) && is.finite(total_initial_at_glucose_start)) {
    total_peak_to_glucose_end - total_initial_at_glucose_start
  } else {
    NA_real_
  }
  live_auc_glucose_window <- if (nrow(count_glucose_window)) {
    safe_trapz(count_glucose_window$hours, count_glucose_window$live_cells)
  } else {
    NA_real_
  }

  out <- tibble(
    well_idx = well_meta_row$well_idx,
    cellLine = well_meta_row$cellLine,
    line_id = well_meta_row$line_id,
    ploidy_metric = well_meta_row$ploidy_metric,
    ploidy_abs = well_meta_row$ploidy_abs,
    G0 = well_meta_row$G0,
    exp_id = well_meta_row$exp_id,
    has_starvation = well_meta_row$has_starvation,
    n_count_times = nrow(count_df),
    n_glucose_times = nrow(glucose_df),
    live_initial = if (nrow(count_df)) count_df$live_cells[1] else NA_real_,
    live_final = if (nrow(count_df)) dplyr::last(count_df$live_cells) else NA_real_,
    live_peak = if (length(peak_live_idx)) count_df$live_cells[peak_live_idx] else NA_real_,
    live_nadir = if (length(nadir_live_idx)) count_df$live_cells[nadir_live_idx] else NA_real_,
    live_peak_time = if (length(peak_live_idx)) count_df$hours[peak_live_idx] else NA_real_,
    live_nadir_time = if (length(nadir_live_idx)) count_df$hours[nadir_live_idx] else NA_real_,
    live_net_change = if (nrow(count_df)) dplyr::last(count_df$live_cells) - count_df$live_cells[1] else NA_real_,
    live_fold_change = if (nrow(count_df)) (dplyr::last(count_df$live_cells) + 1) / (count_df$live_cells[1] + 1) else NA_real_,
    live_auc = safe_trapz(count_df$hours, count_df$live_cells),
    live_auc_glucose_window = live_auc_glucose_window,
    exp_growth_rate_70 = exp_phase_fit$rate,
    exp_growth_intercept_70 = exp_phase_fit$intercept,
    exp_growth_fit_end_time_70 = exp_phase_fit$end_time,
    exp_growth_fit_n_points_70 = exp_phase_fit$n_fit_points,
    exp_growth_fit_r2_70 = exp_phase_fit$fit_r_squared,
    max_growth_rate_raw = live_rate_raw$max_rate,
    max_decline_rate_raw = live_rate_raw$min_rate,
    max_growth_rate = live_rate$max_rate,
    max_decline_rate = live_rate$min_rate,
    time_max_growth_rate = live_rate$time_at_max_rate,
    time_max_decline_rate = live_rate$time_at_min_rate,
    dead_initial = if (nrow(count_df)) count_df$dead_cells[1] else NA_real_,
    dead_final = if (nrow(count_df)) dplyr::last(count_df$dead_cells) else NA_real_,
    dead_peak = if (length(peak_dead_idx)) count_df$dead_cells[peak_dead_idx] else NA_real_,
    dead_peak_time = if (length(peak_dead_idx)) count_df$hours[peak_dead_idx] else NA_real_,
    dead_net_change = if (nrow(count_df)) dplyr::last(count_df$dead_cells) - count_df$dead_cells[1] else NA_real_,
    dead_auc = safe_trapz(count_df$hours, count_df$dead_cells),
    max_death_rate_raw = dead_rate_raw$max_rate,
    max_death_rate = dead_rate$max_rate,
    dead_fraction_initial = if (nrow(count_df)) count_df$dead_fraction[1] else NA_real_,
    dead_fraction_final = if (nrow(count_df)) dplyr::last(count_df$dead_fraction) else NA_real_,
    dead_fraction_peak = if (nrow(count_df)) max(count_df$dead_fraction, na.rm = TRUE) else NA_real_,
    max_dead_fraction_rate = dead_frac_rate$max_rate,
    glucose_initial = if (nrow(glucose_df)) glucose_df$glucose_hat[1] else NA_real_,
    glucose_final = if (nrow(glucose_df)) dplyr::last(glucose_df$glucose_hat) else NA_real_,
    glucose_min = if (length(min_glucose_idx)) glucose_df$glucose_hat[min_glucose_idx] else NA_real_,
    glucose_min_time = if (length(min_glucose_idx)) glucose_df$hours[min_glucose_idx] else NA_real_,
    glucose_drawdown = if (nrow(glucose_df)) glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat) else NA_real_,
    glucose_drawdown_fraction = if (nrow(glucose_df) && glucose_df$glucose_hat[1] > 0) {
      (glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat)) / glucose_df$glucose_hat[1]
    } else {
      NA_real_
    },
    glucose_auc = safe_trapz(glucose_df$hours, glucose_df$glucose_hat),
    max_glucose_drawdown_rate = if (is.finite(glucose_rate$min_rate)) -glucose_rate$min_rate else NA_real_,
    time_max_glucose_drawdown = glucose_rate$time_at_min_rate,
    time_to_glucose_floor = time_to_threshold(glucose_df$hours, glucose_df$glucose_hat, glucose_floor, direction = "below"),
    live_gain_per_glucose = if (nrow(count_df) && nrow(glucose_df)) {
      (dplyr::last(count_df$live_cells) - count_df$live_cells[1]) / max(glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat), 1e-8)
    } else {
      NA_real_
    },
    dead_gain_per_glucose = if (nrow(count_df) && nrow(glucose_df)) {
      (dplyr::last(count_df$dead_cells) - count_df$dead_cells[1]) / max(glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat), 1e-8)
    } else {
      NA_real_
    },
    total_peak = if (nrow(count_df)) max(count_df$total_cells, na.rm = TRUE) else NA_real_,
    total_peak_time = if (nrow(count_df)) count_df$hours[which.max(count_df$total_cells)] else NA_real_,
    total_auc = safe_trapz(count_df$hours, count_df$total_cells),
    glucose_start_time = glucose_start_time,
    glucose_end_time = glucose_end_time,
    total_peak_to_glucose_end = total_peak_to_glucose_end,
    total_initial_at_glucose_start = total_initial_at_glucose_start,
    total_net_gain_to_glucose_end = total_net_gain_to_glucose_end,
    peak_total_yield_per_glucose = if (nrow(glucose_df)) {
      drawdown <- glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat)
      if (is.finite(drawdown) && drawdown >= min_glucose_drawdown_for_yield) {
        pmax(total_net_gain_to_glucose_end, 0) / drawdown
      } else {
        NA_real_
      }
    } else {
      NA_real_
    }
  )

  out
}

build_feature_panel <- function(stan_data_path = NULL, glucose_floor = 0.1, min_glucose_drawdown_for_yield = 0.05) {
  tables <- build_model_free_tables(stan_data_path = stan_data_path)

  feature_panel <- tables$well_meta %>%
    split(.$well_idx) %>%
    map_dfr(function(meta_row) {
      well_idx <- meta_row$well_idx[[1]]
      count_df <- tables$count_summary %>% filter(well_idx == !!well_idx)
      glucose_df <- tables$glucose_summary %>% filter(well_idx == !!well_idx)
      summarize_one_well_features(
        meta_row[1, ],
        count_df,
        glucose_df,
        glucose_floor = glucose_floor,
        min_glucose_drawdown_for_yield = min_glucose_drawdown_for_yield
      )
    }) %>%
    arrange(cellLine, ploidy_metric, G0)

  list(
    stan_data_path = tables$stan_data_path,
    well_meta = tables$well_meta,
    count_summary = tables$count_summary,
    glucose_summary = tables$glucose_summary,
    feature_panel = feature_panel
  )
}

fit_log_glucose_response <- function(g0, y, pred_grid = c(0, 0.1, 0.25, 1, 5, 25), max_g0 = Inf) {
  pred_names <- c("pred_0", "pred_0p1", "pred_0p25", "pred_1", "pred_5", "pred_25")
  ok <- is.finite(g0) & is.finite(y) & g0 <= max_g0
  g0 <- g0[ok]
  y <- y[ok]

  if (!length(y)) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      r_squared = NA_real_,
      n = 0L,
      pred = setNames(rep(NA_real_, length(pred_grid)), pred_names)
    ))
  }

  if (length(unique(g0)) < 2L || length(y) < 3L) {
    med <- median(y, na.rm = TRUE)
    return(list(
      intercept = med,
      slope = 0,
      r_squared = NA_real_,
      n = length(y),
      pred = setNames(rep(med, length(pred_grid)), pred_names)
    ))
  }

  fit <- lm(y ~ log1p(g0))
  pred <- stats::predict(fit, newdata = data.frame(g0 = pred_grid))

  list(
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2]),
    r_squared = summary(fit)$r.squared,
    n = length(y),
    pred = setNames(as.numeric(pred), pred_names)
  )
}

fit_glucose_per_live_auc <- function(
  g0,
  live_auc,
  glucose_drawdown,
  max_g0 = 1,
  include_intercept = TRUE,
  conf_level = 0.95
) {
  ok <- is.finite(g0) & is.finite(live_auc) & is.finite(glucose_drawdown) & g0 <= max_g0
  fit_data <- data.frame(
    g0 = as.numeric(g0[ok]),
    live_auc = as.numeric(live_auc[ok]),
    glucose_drawdown = as.numeric(glucose_drawdown[ok])
  )

  if (nrow(fit_data) < 3L || length(unique(fit_data$live_auc)) < 2L) {
    return(list(
      intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
      slope_conf_low = NA_real_, slope_conf_high = NA_real_,
      r_squared = NA_real_, n = nrow(fit_data),
      auc_min = if (nrow(fit_data)) min(fit_data$live_auc) else NA_real_,
      auc_max = if (nrow(fit_data)) max(fit_data$live_auc) else NA_real_,
      include_intercept = include_intercept, max_g0 = max_g0,
      fit = NULL
    ))
  }

  fit <- if (include_intercept) {
    stats::lm(glucose_drawdown ~ live_auc, data = fit_data)
  } else {
    stats::lm(glucose_drawdown ~ 0 + live_auc, data = fit_data)
  }
  coef_table <- summary(fit)$coefficients
  slope_row <- coef_table["live_auc", , drop = FALSE]
  slope_ci <- stats::confint(fit, "live_auc", level = conf_level)

  list(
    intercept = if (include_intercept) unname(stats::coef(fit)["(Intercept)"]) else 0,
    slope = unname(stats::coef(fit)["live_auc"]),
    slope_se = unname(slope_row[1, "Std. Error"]),
    slope_conf_low = unname(slope_ci[1]),
    slope_conf_high = unname(slope_ci[2]),
    r_squared = summary(fit)$r.squared,
    n = nrow(fit_data),
    auc_min = min(fit_data$live_auc),
    auc_max = max(fit_data$live_auc),
    include_intercept = include_intercept,
    max_g0 = max_g0,
    fit = fit
  )
}

get_model_free_feature_catalog <- function() {
  tibble(
    feature = c(
      "growth_lowG_median",
      "growth_highG_median",
      "death_lowG_median",
      "death_highG_median",
      "yield_alive_auc_intercept",
      "yield_alive_auc_slope",
      "peak_total_yield_intercept",
      "peak_total_yield_slope"
    ),
    short_label = c(
      "Growth Low-G Median",
      "Growth High-G Median",
      "Death Low-G Median",
      "Death High-G Median",
      "Alive AUC Regression Intercept",
      "Alive AUC Regression Slope",
      "Peak Total Yield Regression Intercept",
      "Peak Total Yield Regression Slope"
    ),
    category = c(
      "Growth",
      "Growth",
      "Death",
      "Death",
      "Yield",
      "Yield",
      "Yield",
      "Yield"
    ),
    definition = c(
      "Per condition, compute the smoothed maximum upward derivative of log1p(live_cells) over time, then take the median of those per-condition values across G0 <= 0.25 within each cell line and ploidy.",
      "Restrict to wells with G0 equal to 5 or 25. For each such well, fit log1p(live_cells) ~ hours from time 0 through the first timepoint where live_cells reaches 70% of that well's maximum observed live-cell count. The feature is the median fitted slope within each cell line and ploidy.",
      "Per condition, compute the smoothed maximum upward derivative of log1p(dead_cells) over time, then take the median of those per-condition values across G0 <= 0.25 within each cell line and ploidy.",
      "Per condition, compute the smoothed maximum upward derivative of log1p(dead_cells) over time, then take the median of those per-condition values across G0 >= 1 within each cell line and ploidy.",
      "For each condition, live_auc_glucose_window is the trapezoidal area under the live-cell trajectory between the first and last glucose measurement times. The feature is the intercept from regressing that quantity on log1p(G0), using only G0 <= 1.",
      "For each condition, live_auc_glucose_window is the trapezoidal area under the live-cell trajectory between the first and last glucose measurement times. The feature is the slope from regressing that quantity on log1p(G0), using only G0 <= 1.",
      "For each condition, peak_total_yield_per_glucose is max(live_cells + dead_cells) up to the last glucose measurement, minus the initial total cell count at the first glucose measurement, divided by glucose_initial - glucose_final. Negative net gain is floored at zero before division, and the feature is the intercept from regressing that quantity on log1p(G0).",
      "For each condition, peak_total_yield_per_glucose is max(live_cells + dead_cells) up to the last glucose measurement, minus the initial total cell count at the first glucose measurement, divided by glucose_initial - glucose_final. Negative net gain is floored at zero before division, and the feature is the slope from regressing that quantity on log1p(G0)."
    ),
    rationale = c(
      "Median maximum live-cell growth rate under severe glucose limitation (G0 <= 0.25). Tests whether ploidy changes the ability to grow when glucose is scarce.",
      "Median exponential-phase live-cell growth rate from high-glucose wells (G0 = 5 or 25). Targets the pre-depletion, pre-confluence proliferation constant rather than the steepest local derivative.",
      "Median maximum dead-cell accumulation rate under severe glucose limitation (G0 <= 0.25). Targets starvation-associated death pressure.",
      "Median maximum dead-cell accumulation rate when glucose is more available (G0 >= 1). Checks whether any ploidy-associated death effect persists outside the strongest starvation regime.",
      "Intercept from regression of alive-cell AUC on log1p(G0). Captures baseline cumulative alive biomass at near-zero glucose.",
      "Slope from regression of alive-cell AUC on log1p(G0). Captures how strongly cumulative alive biomass increases as glucose supply rises.",
      "Intercept from regression of peak total-cell yield per glucose consumed on log1p(G0). Estimates baseline conversion efficiency of glucose into total cell output near zero glucose.",
      "Slope from regression of peak total-cell yield per glucose consumed on log1p(G0). Tests whether ploidy changes how glucose-to-cell conversion efficiency scales with available glucose."
    ),
    computation = c(
      "Median of per-condition smoothed max_growth_rate across G0 <= 0.25.",
      "Median of per-well exp_growth_rate_70 across G0 %in% c(5, 25). Each exp_growth_rate_70 comes from lm(log1p(live_cells) ~ hours) fit from t = 0 to the first time live_cells reaches 70% of that well's maximum.",
      "Median of per-condition smoothed max_death_rate across G0 <= 0.25.",
      "Median of per-condition smoothed max_death_rate across G0 >= 1.",
      "Intercept from lm(live_auc_glucose_window ~ log1p(G0)) fit only on G0 <= 1.",
      "Slope from lm(live_auc_glucose_window ~ log1p(G0)) fit only on G0 <= 1.",
      "Intercept from lm(peak_total_yield_per_glucose ~ log1p(G0)).",
      "Slope from lm(peak_total_yield_per_glucose ~ log1p(G0))."
    ),
    interpretation = c(
      "Higher values mean faster best-case live-cell expansion under very low glucose.",
      "Higher values mean a faster estimated exponential-phase live-cell growth constant in high-glucose wells.",
      "Higher values mean faster dead-cell accumulation under very low glucose.",
      "Higher values mean faster dead-cell accumulation when glucose is not strongly limiting.",
      "Higher values mean higher baseline cumulative live-cell burden near zero glucose.",
      "Higher values mean alive-cell burden rises more steeply as initial glucose increases.",
      "Higher values mean better baseline conversion of consumed glucose into total cell output near zero glucose.",
      "Higher values mean glucose-to-total-cell conversion efficiency improves more strongly as initial glucose increases."
    )
  )
}

median_in_g0_band <- function(g0, y, lower = -Inf, upper = Inf) {
  keep <- is.finite(g0) & is.finite(y) & g0 >= lower & g0 <= upper
  if (!any(keep)) {
    return(NA_real_)
  }
  median(y[keep], na.rm = TRUE)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

summarize_glucose_signatures <- function(feature_panel) {
  groups <- split(feature_panel, paste(feature_panel$cellLine, feature_panel$ploidy_metric, sep = "||"))

  sig_rows <- lapply(groups, function(df) {
    df <- df %>% arrange(G0)
    g0 <- df$G0

    fit_alive_auc <- fit_log_glucose_response(g0, df$live_auc_glucose_window, max_g0 = 1)
    fit_peak_total_yield <- fit_log_glucose_response(g0, df$peak_total_yield_per_glucose)

    tibble(
      cellLine = df$cellLine[1],
      line_id = df$line_id[1],
      ploidy_metric = df$ploidy_metric[1],
      ploidy_abs = df$ploidy_abs[1],
      has_starvation = df$has_starvation[1],
      n_glucose_conditions = nrow(df),
      growth_lowG_median = median_in_g0_band(g0, df$max_growth_rate, upper = 0.25),
      growth_highG_median = safe_median(df$exp_growth_rate_70[g0 %in% c(5, 25)]),
      growth_logG_slope = fit_log_glucose_response(g0, df$max_growth_rate)$slope,
      growth_peak_time_highG = median_in_g0_band(g0, df$time_max_growth_rate, lower = 1),
      death_lowG_median = median_in_g0_band(g0, df$max_death_rate, upper = 0.25),
      death_highG_median = median_in_g0_band(g0, df$max_death_rate, lower = 1),
      death_logG_slope = fit_log_glucose_response(g0, df$max_death_rate)$slope,
      yield_alive_auc_intercept = fit_alive_auc$intercept,
      yield_alive_auc_slope = fit_alive_auc$slope,
      yield_alive_auc_pred_0 = fit_alive_auc$pred[["pred_0"]],
      yield_alive_auc_pred_1 = fit_alive_auc$pred[["pred_1"]],
      yield_alive_auc_pred_25 = fit_alive_auc$pred[["pred_25"]],
      peak_total_yield_intercept = fit_peak_total_yield$intercept,
      peak_total_yield_slope = fit_peak_total_yield$slope,
      peak_total_yield_pred_0 = fit_peak_total_yield$pred[["pred_0"]],
      peak_total_yield_pred_1 = fit_peak_total_yield$pred[["pred_1"]],
      peak_total_yield_pred_25 = fit_peak_total_yield$pred[["pred_25"]]
    )
  })

  bind_rows(sig_rows) %>% arrange(cellLine, ploidy_metric)
}

compute_empirical_effects <- function(signature_panel) {
  feature_cols <- get_model_free_feature_catalog()$feature

  paired <- signature_panel %>%
    group_by(cellLine) %>%
    filter(n() == 2L) %>%
    arrange(ploidy_metric, .by_group = TRUE) %>%
    summarise(
      observed_ploidy = first(ploidy_metric),
      holdout_ploidy = last(ploidy_metric),
      delta_ploidy = last(ploidy_metric) - first(ploidy_metric),
      across(all_of(feature_cols), ~ (dplyr::last(.) - dplyr::first(.)), .names = "{.col}_delta"),
      .groups = "drop"
    )

  long <- paired %>%
    pivot_longer(
      cols = ends_with("_delta"),
      names_to = "feature",
      values_to = "delta_value"
    ) %>%
    mutate(
      feature = sub("_delta$", "", feature),
      effect_per_ploidy = delta_value / delta_ploidy
    )

  effect_matrix <- long %>%
    select(cellLine, feature, effect_per_ploidy) %>%
    tidyr::pivot_wider(names_from = feature, values_from = effect_per_ploidy) %>%
    arrange(cellLine)

  pca <- NULL
  score_df <- tibble()
  loading_df <- tibble()

  mat <- effect_matrix %>%
    ungroup() %>%
    select(-cellLine) %>%
    as.matrix()

  keep_cols <- apply(mat, 2, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  if (nrow(mat) >= 3L && sum(keep_cols) >= 2L) {
    mat_sub <- mat[, keep_cols, drop = FALSE]
    pca <- prcomp(mat_sub, center = TRUE, scale. = TRUE)
    score_df <- bind_cols(
      effect_matrix %>% select(cellLine),
      as_tibble(pca$x, .name_repair = "minimal")
    )
    loading_df <- tibble(
      feature = rownames(pca$rotation),
      PC1 = pca$rotation[, 1],
      PC2 = if (ncol(pca$rotation) >= 2L) pca$rotation[, 2] else NA_real_
    )
  }

  list(
    feature_cols = feature_cols,
    paired_effects_wide = paired,
    paired_effects_long = long,
    effect_matrix = effect_matrix,
    pca = pca,
    pca_scores = score_df,
    pca_loadings = loading_df
  )
}

safe_mad <- function(x, center = median(x, na.rm = TRUE), constant = 1.4826) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  mad_val <- stats::mad(x, center = center, constant = constant, na.rm = TRUE)
  if (is.finite(mad_val) && mad_val > 0) {
    return(mad_val)
  }

  iqr_val <- stats::IQR(x, na.rm = TRUE) / 1.349
  if (is.finite(iqr_val) && iqr_val > 0) {
    return(iqr_val)
  }

  sd_val <- stats::sd(x, na.rm = TRUE)
  if (is.finite(sd_val) && sd_val > 0) {
    return(sd_val)
  }

  1.0
}

compute_signature_outlier_scores <- function(signature_panel, sign_tol = 1e-8) {
  effects <- compute_empirical_effects(signature_panel)
  effect_matrix <- effects$effect_matrix

  if (!nrow(effect_matrix)) {
    return(list(
      feature_stats = tibble(),
      line_scores = tibble()
    ))
  }

  feature_tbl <- effect_matrix %>%
    pivot_longer(cols = -cellLine, names_to = "feature", values_to = "effect_per_ploidy") %>%
    group_by(feature) %>%
    summarise(
      center = median(effect_per_ploidy, na.rm = TRUE),
      scale = safe_mad(effect_per_ploidy),
      .groups = "drop"
    )

  line_scores <- effect_matrix %>%
    pivot_longer(cols = -cellLine, names_to = "feature", values_to = "effect_per_ploidy") %>%
    left_join(feature_tbl, by = "feature") %>%
    mutate(
      robust_z = if_else(
        is.finite(effect_per_ploidy) & is.finite(center) & is.finite(scale) & scale > 0,
        (effect_per_ploidy - center) / scale,
        NA_real_
      ),
      sign_discordant = if_else(
        is.finite(effect_per_ploidy) & is.finite(center) &
          abs(effect_per_ploidy) > sign_tol & abs(center) > sign_tol,
        as.integer(sign(effect_per_ploidy) != sign(center)),
        NA_integer_
      )
    ) %>%
    group_by(cellLine) %>%
    summarise(
      n_features = sum(is.finite(robust_z)),
      robust_distance_sq = sum(robust_z[is.finite(robust_z)] ^ 2),
      mean_abs_robust_z = mean(abs(robust_z), na.rm = TRUE),
      max_abs_robust_z = max(abs(robust_z), na.rm = TRUE),
      sign_discordance_rate = mean(sign_discordant, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(robust_distance_sq), desc(mean_abs_robust_z), desc(sign_discordance_rate), cellLine)

  list(
    feature_stats = feature_tbl,
    line_scores = line_scores
  )
}

evaluate_transfer_predictions <- function(signature_panel, excluded_training_lines = NULL) {
  feature_cols <- get_model_free_feature_catalog()$feature
  excluded_training_lines <- unique(excluded_training_lines[!is.na(excluded_training_lines) & nzchar(excluded_training_lines)])

  rows <- list()
  idx <- 1L

  for (line_name in sort(unique(signature_panel$cellLine))) {
    line_df <- signature_panel %>% filter(cellLine == line_name)
    states <- sort(unique(line_df$ploidy_metric))

    if (length(states) != 2L) {
      next
    }

    for (direction in c("low_to_high", "high_to_low")) {
      observed_ploidy <- if (direction == "low_to_high") states[1] else states[2]
      holdout_ploidy <- if (direction == "low_to_high") states[2] else states[1]
      delta_ploidy <- holdout_ploidy - observed_ploidy

      obs_row <- signature_panel %>%
        filter(cellLine == line_name, abs(ploidy_metric - observed_ploidy) < 1e-12)

      hold_row <- signature_panel %>%
        filter(cellLine == line_name, abs(ploidy_metric - holdout_ploidy) < 1e-12)

      train_pairs <- signature_panel %>%
        filter(cellLine != line_name, !(cellLine %in% excluded_training_lines)) %>%
        group_by(cellLine) %>%
        filter(n() == 2L) %>%
        summarise(
          low_ploidy = min(ploidy_metric),
          high_ploidy = max(ploidy_metric),
          delta_ploidy_train = high_ploidy - low_ploidy,
          across(all_of(feature_cols), ~ .[which.min(ploidy_metric)], .names = "low__{.col}"),
          across(all_of(feature_cols), ~ .[which.max(ploidy_metric)], .names = "high__{.col}"),
          .groups = "drop"
        )

      if (!nrow(obs_row) || !nrow(hold_row) || !nrow(train_pairs)) {
        next
      }

      for (feature in feature_cols) {
        obs_val <- obs_row[[feature]][1]
        true_val <- hold_row[[feature]][1]

        low_col <- paste0("low__", feature)
        high_col <- paste0("high__", feature)
        low_vals <- train_pairs[[low_col]]
        high_vals <- train_pairs[[high_col]]
        slope_vals <- (high_vals - low_vals) / train_pairs$delta_ploidy_train

        effect_hat <- median(slope_vals[is.finite(slope_vals)], na.rm = TRUE)
        if (!is.finite(effect_hat) || !is.finite(obs_val) || !is.finite(true_val)) {
          next
        }

        pred_null <- obs_val
        pred_transfer <- obs_val + effect_hat * delta_ploidy
        true_effect <- (true_val - obs_val) / delta_ploidy
        scale_ref <- median(abs(c(low_vals, high_vals, obs_val, true_val)), na.rm = TRUE)
        scale_ref <- max(scale_ref, 1e-8)
        abs_err_null <- abs(pred_null - true_val)
        abs_err_transfer <- abs(pred_transfer - true_val)

        rows[[idx]] <- tibble(
          cellLine = line_name,
          direction = direction,
          feature = feature,
          observed_ploidy = observed_ploidy,
          holdout_ploidy = holdout_ploidy,
          delta_ploidy = delta_ploidy,
          observed_value = obs_val,
          true_value = true_val,
          true_effect = true_effect,
          transfer_effect = effect_hat,
          pred_null = pred_null,
          pred_transfer = pred_transfer,
          abs_err_null = abs_err_null,
          abs_err_transfer = abs_err_transfer,
          scaled_abs_err_null = abs_err_null / scale_ref,
          scaled_abs_err_transfer = abs_err_transfer / scale_ref,
          sq_err_null = (pred_null - true_val) ^ 2,
          sq_err_transfer = (pred_transfer - true_val) ^ 2,
          abs_err_improvement = abs_err_null - abs_err_transfer,
          scaled_abs_err_improvement = (abs_err_null - abs_err_transfer) / scale_ref,
          sign_match = as.integer(sign(true_effect) == sign(effect_hat)),
          n_training_lines = sum(is.finite(slope_vals)),
          excluded_training_lines = if (length(excluded_training_lines)) {
            paste(sort(excluded_training_lines), collapse = "|")
          } else {
            "none"
          }
        )
        idx <- idx + 1L
      }
    }
  }

  predictions <- bind_rows(rows)

  feature_summary <- predictions %>%
    group_by(feature) %>%
    summarise(
      n_cases = n(),
      mae_null = mean(abs_err_null, na.rm = TRUE),
      mae_transfer = mean(abs_err_transfer, na.rm = TRUE),
      scaled_mae_null = mean(scaled_abs_err_null, na.rm = TRUE),
      scaled_mae_transfer = mean(scaled_abs_err_transfer, na.rm = TRUE),
      rmse_null = sqrt(mean(sq_err_null, na.rm = TRUE)),
      rmse_transfer = sqrt(mean(sq_err_transfer, na.rm = TRUE)),
      mean_abs_err_improvement = mean(abs_err_improvement, na.rm = TRUE),
      median_abs_err_improvement = median(abs_err_improvement, na.rm = TRUE),
      mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
      median_scaled_err_improvement = median(scaled_abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_scaled_err_improvement), scaled_mae_transfer)

  case_summary <- predictions %>%
    group_by(cellLine, direction) %>%
    summarise(
      n_features = n(),
      mean_abs_err_improvement = mean(abs_err_improvement, na.rm = TRUE),
      mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cellLine, direction)

  list(
    feature_cols = feature_cols,
    predictions = predictions,
    feature_summary = feature_summary,
    case_summary = case_summary
  )
}

summarize_transfer_performance <- function(predictions) {
  if (!nrow(predictions)) {
    return(tibble(
      n_predictions = integer(),
      n_cell_lines = integer(),
      mean_scaled_err_improvement = numeric(),
      median_scaled_err_improvement = numeric(),
      transfer_win_rate = numeric(),
      sign_accuracy = numeric()
    ))
  }

  tibble(
    n_predictions = nrow(predictions),
    n_cell_lines = dplyr::n_distinct(predictions$cellLine),
    mean_scaled_err_improvement = mean(predictions$scaled_abs_err_improvement, na.rm = TRUE),
    median_scaled_err_improvement = median(predictions$scaled_abs_err_improvement, na.rm = TRUE),
    transfer_win_rate = mean(predictions$abs_err_transfer < predictions$abs_err_null, na.rm = TRUE),
    sign_accuracy = mean(predictions$sign_match, na.rm = TRUE)
  )
}

analyze_model_free_outlier_robustness <- function(signature_panel) {
  outlier_scores <- compute_signature_outlier_scores(signature_panel)
  baseline <- evaluate_transfer_predictions(signature_panel)
  baseline_predictions <- baseline$predictions
  line_ids <- sort(unique(signature_panel$cellLine))

  baseline_overall <- summarize_transfer_performance(baseline_predictions) %>%
    mutate(summary_scope = "all_lines")

  exclusion_summary <- bind_rows(lapply(line_ids, function(excluded_line) {
    excl <- evaluate_transfer_predictions(
      signature_panel,
      excluded_training_lines = excluded_line
    )

    remaining_baseline <- baseline_predictions %>%
      filter(cellLine != excluded_line)
    remaining_excl <- excl$predictions %>%
      filter(cellLine != excluded_line)

    base_sum <- summarize_transfer_performance(remaining_baseline)
    excl_sum <- summarize_transfer_performance(remaining_excl)

    tibble(
      excluded_line = excluded_line,
      n_holdout_lines_evaluated = dplyr::n_distinct(remaining_excl$cellLine),
      n_predictions_evaluated = nrow(remaining_excl),
      baseline_mean_scaled_err_improvement = base_sum$mean_scaled_err_improvement,
      excluded_mean_scaled_err_improvement = excl_sum$mean_scaled_err_improvement,
      delta_mean_scaled_err_improvement = excl_sum$mean_scaled_err_improvement - base_sum$mean_scaled_err_improvement,
      baseline_median_scaled_err_improvement = base_sum$median_scaled_err_improvement,
      excluded_median_scaled_err_improvement = excl_sum$median_scaled_err_improvement,
      delta_median_scaled_err_improvement = excl_sum$median_scaled_err_improvement - base_sum$median_scaled_err_improvement,
      baseline_transfer_win_rate = base_sum$transfer_win_rate,
      excluded_transfer_win_rate = excl_sum$transfer_win_rate,
      delta_transfer_win_rate = excl_sum$transfer_win_rate - base_sum$transfer_win_rate,
      baseline_sign_accuracy = base_sum$sign_accuracy,
      excluded_sign_accuracy = excl_sum$sign_accuracy,
      delta_sign_accuracy = excl_sum$sign_accuracy - base_sum$sign_accuracy
    )
  })) %>%
    left_join(outlier_scores$line_scores, by = c("excluded_line" = "cellLine")) %>%
    mutate(
      robust_distance_rank = min_rank(desc(robust_distance_sq)),
      influence_rank = min_rank(desc(delta_mean_scaled_err_improvement)),
      composite_rank = robust_distance_rank + influence_rank
    ) %>%
    arrange(desc(delta_mean_scaled_err_improvement), desc(delta_transfer_win_rate), desc(robust_distance_sq), excluded_line)

  feature_exclusion_summary <- bind_rows(lapply(line_ids, function(excluded_line) {
    excl <- evaluate_transfer_predictions(
      signature_panel,
      excluded_training_lines = excluded_line
    )

    base_feat <- baseline_predictions %>%
      filter(cellLine != excluded_line) %>%
      group_by(feature) %>%
      summarise(
        baseline_mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
        baseline_transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
        baseline_sign_accuracy = mean(sign_match, na.rm = TRUE),
        .groups = "drop"
      )

    excl_feat <- excl$predictions %>%
      filter(cellLine != excluded_line) %>%
      group_by(feature) %>%
      summarise(
        excluded_mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
        excluded_transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
        excluded_sign_accuracy = mean(sign_match, na.rm = TRUE),
        .groups = "drop"
      )

    base_feat %>%
      inner_join(excl_feat, by = "feature") %>%
      mutate(
        excluded_line = excluded_line,
        delta_mean_scaled_err_improvement = excluded_mean_scaled_err_improvement - baseline_mean_scaled_err_improvement,
        delta_transfer_win_rate = excluded_transfer_win_rate - baseline_transfer_win_rate,
        delta_sign_accuracy = excluded_sign_accuracy - baseline_sign_accuracy
      )
  })) %>%
    arrange(desc(delta_mean_scaled_err_improvement), desc(delta_transfer_win_rate), excluded_line, feature)

  recommended_exclusion <- exclusion_summary %>%
    arrange(composite_rank, influence_rank, robust_distance_rank, excluded_line) %>%
    slice_head(n = 1) %>%
    mutate(
      recommendation_flag = is.finite(delta_mean_scaled_err_improvement) &
        delta_mean_scaled_err_improvement > 0
    )

  list(
    baseline = baseline,
    baseline_overall = baseline_overall,
    outlier_scores = outlier_scores$line_scores,
    outlier_feature_stats = outlier_scores$feature_stats,
    exclusion_summary = exclusion_summary,
    feature_exclusion_summary = feature_exclusion_summary,
    recommended_exclusion = recommended_exclusion
  )
}
