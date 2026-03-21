knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.height = 6
)
knitr::opts_knit$set(root.dir = normalizePath("../"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

resolve_area_stan_data_path <- function() {
  candidates <- c(
    file.path("data", "stan_ready_data_with_area.Rds"),
    file.path("dev_data", "stan_ready_data_with_area.Rds"),
    file.path("data", "stan_ready_data.Rds"),
    file.path("ecology", "stan_ready_data.Rds")
  )

  for (path in candidates) {
    if (file.exists(path)) {
      return(path)
    }
  }

  stop(
    "Could not find an area-aware stan data object. Tried: ",
    paste(candidates, collapse = ", ")
  )
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

smooth_derivative <- function(hours, value, spar = 0.55) {
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (length(hours) < 4L || length(unique(hours)) < 4L) {
    return(tibble(
      hours = hours,
      value_smooth = value,
      deriv = c(NA_real_, diff(value) / diff(hours))
    ))
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  fit <- try(stats::smooth.spline(hours, value, spar = spar), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(
      hours = hours,
      value_smooth = value,
      deriv = c(NA_real_, diff(value) / diff(hours))
    ))
  }

  tibble(
    hours = hours,
    value_smooth = as.numeric(stats::predict(fit, x = hours, deriv = 0)$y),
    deriv = as.numeric(stats::predict(fit, x = hours, deriv = 1)$y)
  )
}

lagged_correlations <- function(hours, x, y, lag_hours = seq(-24, 24, by = 8)) {
  out <- lapply(lag_hours, function(lag_h) {
    x_shift <- approx(
      x = hours + lag_h,
      y = x,
      xout = hours,
      rule = 1
    )$y

    tibble(
      lag_hours = lag_h,
      corr = suppressWarnings(cor(x_shift, y, use = "complete.obs")),
      n_overlap = sum(is.finite(x_shift) & is.finite(y))
    )
  })

  bind_rows(out)
}

fit_cost_model <- function(y, x_build, x_maint) {
  ok <- is.finite(y) & is.finite(x_build) & is.finite(x_maint)
  y <- y[ok]
  x_build <- x_build[ok]
  x_maint <- x_maint[ok]

  if (length(y) < 4L) {
    return(c(a = NA_real_, m = NA_real_, rmse = NA_real_, n = length(y)))
  }

  obj <- function(par) {
    pred <- par[1] * x_build + par[2] * x_maint
    sum((y - pred)^2)
  }

  fit0 <- try(lm(y ~ 0 + x_build + x_maint), silent = TRUE)
  par0 <- if (inherits(fit0, "try-error")) {
    c(0.1, 0.001)
  } else {
    pmax(as.numeric(coef(fit0)), 0)
  }
  if (length(par0) != 2L || any(!is.finite(par0))) {
    par0 <- c(0.1, 0.001)
  }

  opt <- optim(
    par = par0,
    fn = obj,
    method = "L-BFGS-B",
    lower = c(0, 0)
  )

  pred <- opt$par[1] * x_build + opt$par[2] * x_maint
  rmse <- sqrt(mean((y - pred)^2))

  c(a = unname(opt$par[1]), m = unname(opt$par[2]), rmse = rmse, n = length(y))
}

stan_path <- resolve_area_stan_data_path()
stan <- readRDS(stan_path)

if (is.null(stan$area_alive_quantiles)) {
  stop("Loaded stan data does not contain area_alive_quantiles.")
}

line_names <- names(stan$line_map)[match(stan$line_id, unname(stan$line_map))]
well_meta <- tibble(
  well_idx = seq_len(stan$N_wells),
  cellLine = line_names,
  line_id = as.integer(stan$line_id),
  ploidy_metric = as.numeric(stan$ploidy_metric),
  ploidy_abs = as.numeric(stan$ploidy_abs),
  G0 = as.numeric(stan$G0_per_well),
  exp_id = as.integer(stan$exp_id),
  has_starvation = as.integer(stan$has_starvation)
) %>%
  mutate(
    ploidy_label = if_else(ploidy_metric == ave(ploidy_metric, cellLine, FUN = min), "baseline", "elevated"),
    G0_label = factor(as.character(G0), levels = as.character(sort(unique(G0))))
  )

area_q <- as.data.frame(stan$area_alive_quantiles, check.names = FALSE)

count_obs <- tibble(
  well_idx = as.integer(stan$well_idx_count),
  rep_id = as.character(stan$rep_id_count),
  hours = as.numeric(stan$t_grid[stan$grid_idx_count]),
  live_cells = as.numeric(stan$N_obs),
  dead_cells = as.numeric(stan$D_obs),
  n_area_alive = as.numeric(stan$n_area_alive %||% rep(NA_real_, length(stan$N_obs))),
  q0 = area_q[["q0"]],
  q01 = area_q[["q0.01"]],
  q05 = area_q[["q0.05"]],
  q10 = area_q[["q0.1"]],
  q25 = area_q[["q0.25"]],
  q50 = area_q[["q0.5"]],
  q75 = area_q[["q0.75"]],
  q90 = area_q[["q0.9"]],
  q95 = area_q[["q0.95"]],
  q99 = area_q[["q0.99"]],
  q100 = area_q[["q1"]]
) %>%
  mutate(
    total_cells = live_cells + dead_cells,
    dead_fraction = if_else(total_cells > 0, dead_cells / total_cells, NA_real_),
    tail_width = q99 - q50,
    tail_ratio = q99 / pmax(q50, 1e-8)
  ) %>%
  left_join(well_meta, by = "well_idx")

glucose_obs <- tibble(
  well_idx = as.integer(stan$well_idx_gluc),
  hours = as.numeric(stan$t_grid[stan$grid_idx_gluc]),
  lum = as.numeric(stan$lum_obs),
  dilution = as.numeric(stan$dilution),
  is_censored = as.integer(stan$is_censored)
) %>%
  left_join(well_meta %>% select(well_idx, exp_id), by = "well_idx") %>%
  mutate(
    calib_a = as.numeric(stan$calib_a_fixed[exp_id]),
    calib_b = as.numeric(stan$calib_b_fixed[exp_id]),
    glucose_hat = pmax(0, (lum - calib_b) / (calib_a * pmax(dilution, 1e-12)))
) %>%
  group_by(well_idx, hours) %>%
  summarise(
    glucose_hat = {
      uncens <- glucose_hat[is_censored == 0L]
      if (length(uncens)) median(uncens, na.rm = TRUE) else NA_real_
    },
    any_glucose_censored = any(is_censored == 1L, na.rm = TRUE),
    n_glucose_uncensored = sum(is_censored == 0L, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(well_meta, by = "well_idx")

condition_counts <- count_obs %>%
  group_by(cellLine, ploidy_metric, ploidy_label, G0, G0_label, hours) %>%
  summarise(
    n_reps = n(),
    live_cells = median(live_cells, na.rm = TRUE),
    dead_cells = median(dead_cells, na.rm = TRUE),
    total_cells = median(total_cells, na.rm = TRUE),
    dead_fraction = median(dead_fraction, na.rm = TRUE),
    q50 = median(q50, na.rm = TRUE),
    q90 = median(q90, na.rm = TRUE),
    q99 = median(q99, na.rm = TRUE),
    tail_width = median(tail_width, na.rm = TRUE),
    tail_ratio = median(tail_ratio, na.rm = TRUE),
    .groups = "drop"
  )

condition_glucose <- glucose_obs %>%
  group_by(cellLine, ploidy_metric, ploidy_label, G0, G0_label, hours) %>%
  summarise(
    glucose_hat = median(glucose_hat, na.rm = TRUE),
    .groups = "drop"
  )

count_glucose_joint <- condition_counts %>%
  left_join(
    condition_glucose,
    by = c("cellLine", "ploidy_metric", "ploidy_label", "G0", "G0_label", "hours")
  )

paired_ploidy <- condition_counts %>%
  group_by(cellLine, G0, G0_label, hours) %>%
  filter(n_distinct(ploidy_metric) == 2L) %>%
  arrange(ploidy_metric, .by_group = TRUE) %>%
  summarise(
    delta_live = live_cells[2] - live_cells[1],
    delta_dead_fraction = dead_fraction[2] - dead_fraction[1],
    delta_q50 = q50[2] - q50[1],
    delta_q90 = q90[2] - q90[1],
    delta_q99 = q99[2] - q99[1],
    delta_tail_ratio = tail_ratio[2] - tail_ratio[1],
    .groups = "drop"
  )

rate_obs <- count_obs %>%
  group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_label, G0, G0_label) %>%
  group_modify(function(df, keys) {
    death_fit <- smooth_derivative(df$hours, df$dead_cells)
    area_fit <- smooth_derivative(df$hours, df$q50)
    total_ref <- df %>%
      select(hours, total_cells) %>%
      distinct()

    full_join(
      death_fit %>% rename(death_cells_smooth = value_smooth, death_deriv = deriv),
      area_fit %>% rename(area_q50_smooth = value_smooth, area_rate = deriv),
      by = "hours"
    ) %>%
      left_join(total_ref, by = "hours") %>%
      transmute(
        hours = hours,
        death_cells_smooth = death_cells_smooth,
        death_rate_per_cell = pmax(death_deriv, 0) / pmax(total_cells, 1e-8),
        area_q50_smooth = area_q50_smooth,
        area_rate = area_rate
      )
  }) %>%
  ungroup()

lag_profiles <- rate_obs %>%
  group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_label, G0, G0_label) %>%
  group_modify(function(df, keys) {
    lagged_correlations(df$hours, df$area_rate, df$death_rate_per_cell)
  }) %>%
  ungroup()

lag_summary <- lag_profiles %>%
  filter(n_overlap >= 4) %>%
  group_by(cellLine, ploidy_label, G0, G0_label, lag_hours) %>%
  summarise(
    median_corr = median(corr, na.rm = TRUE),
    mean_corr = mean(corr, na.rm = TRUE),
    n_profiles = sum(is.finite(corr)),
    .groups = "drop"
  )

best_lag_tbl <- lag_profiles %>%
  filter(n_overlap >= 4, is.finite(corr)) %>%
  group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_label, G0, G0_label) %>%
  summarise(
    best_abs_corr = corr[which.max(abs(corr))],
    best_lag_hours = lag_hours[which.max(abs(corr))],
    .groups = "drop"
  )

glucose_last <- glucose_obs %>%
  group_by(well_idx) %>%
  summarise(
    t_first_glucose = if (any(is.finite(glucose_hat))) min(hours[is.finite(glucose_hat)], na.rm = TRUE) else NA_real_,
    t_last_glucose = if (any(is.finite(glucose_hat))) max(hours[is.finite(glucose_hat)], na.rm = TRUE) else NA_real_,
    glucose_initial = if (is.finite(t_first_glucose)) glucose_hat[match(t_first_glucose, hours)] else NA_real_,
    glucose_last = if (is.finite(t_last_glucose)) glucose_hat[match(t_last_glucose, hours)] else NA_real_,
    glucose_consumed = pmax(0, glucose_initial - glucose_last),
    n_uncensored_timepoints = sum(is.finite(glucose_hat)),
    .groups = "drop"
  )

replicate_cost_df <- count_obs %>%
  left_join(glucose_last, by = "well_idx") %>%
  filter(is.finite(t_last_glucose), hours <= t_last_glucose) %>%
  group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_abs, ploidy_label, G0, G0_label, glucose_consumed) %>%
  summarise(
    B0 = total_cells[which.min(hours)],
    B_end = total_cells[which.max(hours)],
    dB = pmax(0, B_end - B0),
    AUC_live = safe_trapz(hours, live_cells),
    area_24h = safe_median(q50[hours >= 20 & hours <= 28]),
    .groups = "drop"
  ) %>%
  filter(is.finite(glucose_consumed), is.finite(dB), is.finite(AUC_live), is.finite(area_24h)) %>%
  mutate(
    vol_flat = area_24h,
    vol_sphere = area_24h^1.5
  )

cost_fits <- replicate_cost_df %>%
  filter(G0 != 25) %>%
  group_by(cellLine, ploidy_metric, ploidy_abs, ploidy_label) %>%
  summarise(
    fit = list(fit_cost_model(glucose_consumed, dB, AUC_live)),
    .groups = "drop"
  ) %>%
  mutate(
    a = map_dbl(fit, ~ .x[["a"]]),
    m = map_dbl(fit, ~ .x[["m"]]),
    rmse = map_dbl(fit, ~ .x[["rmse"]]),
    n = map_dbl(fit, ~ .x[["n"]])
  ) %>%
  select(-fit) %>%
  left_join(
    replicate_cost_df %>%
      group_by(cellLine, ploidy_metric, ploidy_abs, ploidy_label) %>%
      summarise(
        area_24h = median(area_24h, na.rm = TRUE),
        vol_flat = median(vol_flat, na.rm = TRUE),
        vol_sphere = median(vol_sphere, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("cellLine", "ploidy_metric", "ploidy_abs", "ploidy_label")
  ) %>%
  group_by(cellLine) %>%
  mutate(
    a_base = a[which.min(ploidy_metric)],
    vol_flat_base = vol_flat[which.min(ploidy_metric)],
    vol_sphere_base = vol_sphere[which.min(ploidy_metric)],
    ratio_cost = a / a_base,
    ratio_vol_flat = vol_flat / vol_flat_base,
    ratio_vol_sphere = vol_sphere / vol_sphere_base,
    cost_per_flat_vol = a / vol_flat,
    cost_per_sphere_vol = a / vol_sphere
  ) %>%
  ungroup()

cost_fit_diagnostics <- replicate_cost_df %>%
  left_join(
    cost_fits %>%
      select(cellLine, ploidy_metric, a, m),
    by = c("cellLine", "ploidy_metric")
  ) %>%
  mutate(
    glucose_pred = a * dB + m * AUC_live,
    resid = glucose_consumed - glucose_pred
  )

glucose_init_check <- glucose_last %>%
  left_join(well_meta, by = "well_idx")

rate_summary <- rate_obs %>%
  group_by(cellLine, ploidy_label, G0, G0_label, hours) %>%
  summarise(
    death_rate_per_cell = median(death_rate_per_cell, na.rm = TRUE),
    area_rate = median(area_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(death_rate_per_cell, area_rate),
    names_to = "signal",
    values_to = "value"
  ) %>%
  mutate(
    signal = recode(
      signal,
      death_rate_per_cell = "1/N dD/dt",
      area_rate = "dA/dt"
    )
  )

cellline_ploidy_capacity <- count_obs %>%
  group_by(cellLine, ploidy_metric, ploidy_label) %>%
  summarise(max_total_cells = max(total_cells, na.rm = TRUE), .groups = "drop")

count_glucose_aligned <- count_obs %>%
  group_by(well_idx, rep_id) %>%
  group_modify(function(df, keys) {
    gluc_df <- glucose_obs %>%
      filter(well_idx == keys$well_idx[[1]]) %>%
      arrange(hours)

    if (!nrow(gluc_df)) {
      return(df %>%
        mutate(
          glucose_hat = NA_real_,
          glucose_lag24 = NA_real_,
          glucose_lag48 = NA_real_
        ))
    }

    approx_glucose <- function(target_hours) {
      ok <- is.finite(gluc_df$hours) & is.finite(gluc_df$glucose_hat)
      if (sum(ok) < 2L) {
        return(rep(NA_real_, length(target_hours)))
      }
      approx(
        x = gluc_df$hours[ok],
        y = gluc_df$glucose_hat[ok],
        xout = target_hours,
        rule = 1
      )$y
    }

    df %>%
      arrange(hours) %>%
      mutate(
        glucose_hat = approx_glucose(hours),
        glucose_lag24 = approx_glucose(hours - 24),
        glucose_lag48 = approx_glucose(hours - 48)
      )
  }) %>%
  ungroup() %>%
  left_join(
    cellline_ploidy_capacity,
    by = c("cellLine", "ploidy_metric", "ploidy_label")
  ) %>%
  filter(total_cells <= 0.7 * max_total_cells)

morph_glucose_summary <- count_glucose_aligned %>%
  mutate(
    glucose_hat_plot = pmax(glucose_hat, 0.02),
    glucose_lag24_plot = pmax(glucose_lag24, 0.02),
    glucose_lag48_plot = pmax(glucose_lag48, 0.02)
  ) %>%
  group_by(cellLine, ploidy_metric, ploidy_label, G0, G0_label, hours) %>%
  summarise(
    q50 = median(q50, na.rm = TRUE),
    glucose_hat_plot = median(glucose_hat_plot, na.rm = TRUE),
    glucose_lag24_plot = median(glucose_lag24_plot, na.rm = TRUE),
    glucose_lag48_plot = median(glucose_lag48_plot, na.rm = TRUE),
    .groups = "drop"
  )

cost_fit_timecourse <- count_glucose_aligned %>%
  left_join(
    cost_fits %>%
      select(cellLine, ploidy_metric, a, m),
    by = c("cellLine", "ploidy_metric")
  ) %>%
  filter(is.finite(glucose_hat), is.finite(a), is.finite(m)) %>%
  group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_label, G0, G0_label, a, m) %>%
  arrange(hours, .by_group = TRUE) %>%
  mutate(
    glucose_initial_timecourse = glucose_hat[which.min(hours)],
    glucose_consumed_to_t = pmax(0, glucose_initial_timecourse - glucose_hat),
    AUC_live_to_t = cumsum(c(
      0,
      diff(hours) * (head(live_cells, -1) + tail(live_cells, -1)) / 2
    )),
    B0 = total_cells[which.min(hours)],
    dB_pred = pmax((glucose_consumed_to_t - m * AUC_live_to_t) / pmax(a, 1e-8), 0),
    total_cells_pred = B0 + dB_pred
  ) %>%
  ungroup()

cost_fit_timecourse_summary <- cost_fit_timecourse %>%
  group_by(cellLine, ploidy_metric, ploidy_label, G0, G0_label, hours) %>%
  summarise(
    total_cells_obs = median(total_cells, na.rm = TRUE),
    total_cells_pred = median(total_cells_pred, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(total_cells_obs, total_cells_pred),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      total_cells_obs = "Observed total cells",
      total_cells_pred = "Predicted total cells"
    )
  )

coverage_tbl <- tibble(
  metric = c("count rows", "rows with any area summary", "rows missing all area quantiles"),
  value = c(
    nrow(count_obs),
    sum(rowSums(is.na(stan$area_alive_quantiles)) < ncol(stan$area_alive_quantiles)),
    sum(rowSums(is.na(stan$area_alive_quantiles)) == ncol(stan$area_alive_quantiles))
  )
)

knitr::kable(coverage_tbl)

plot_traj_counts <- condition_counts %>%
  pivot_longer(
    cols = c(live_cells, dead_cells, total_cells),
    names_to = "state",
    values_to = "value"
  ) %>%
  mutate(
    state = recode(
      state,
      live_cells = "Live",
      dead_cells = "Dead",
      total_cells = "Total"
    )
  ) %>%
  ggplot(
    aes(hours, value, color = ploidy_label, group = interaction(ploidy_metric, G0))
  ) +
  geom_line(linewidth = 0.7) +
  facet_grid(state + G0_label ~ cellLine, scales = "free_y") +
  scale_color_manual(values = c("baseline" = "#0072B2", "elevated" = "#D55E00")) +
  theme_bw() +
  labs(
    title = "Count Trajectories",
    x = "Hours",
    y = "Cells",
    color = "Ploidy"
  )

plot_traj_q50 <- ggplot(
  condition_counts,
  aes(hours, q50, color = ploidy_label, group = interaction(ploidy_metric, G0))
) +
  geom_line(linewidth = 0.7) +
  facet_grid(G0_label ~ cellLine, scales = "free_y") +
  scale_color_manual(values = c("baseline" = "#0072B2", "elevated" = "#D55E00")) +
  theme_bw() +
  labs(
    title = "Median Alive-Cell Area Trajectories",
    x = "Hours",
    y = "Alive-cell area q50",
    color = "Ploidy"
  )

plot_traj_tail <- ggplot(
  condition_counts,
  aes(hours, tail_ratio, color = ploidy_label, group = interaction(ploidy_metric, G0))
) +
  geom_line(linewidth = 0.7) +
  facet_grid(G0_label ~ cellLine, scales = "free_y") +
  scale_color_manual(values = c("baseline" = "#0072B2", "elevated" = "#D55E00")) +
  theme_bw() +
  labs(
    title = "Upper-Tail Inflation Over Time",
    subtitle = "tail_ratio = q99 / q50",
    x = "Hours",
    y = "q99 / q50",
    color = "Ploidy"
  )

(plot_traj_counts / plot_traj_q50 / plot_traj_tail) + plot_layout(heights = c(1.6, 1.2, 1.2))

size_summary_tbl <- condition_counts %>%
  group_by(cellLine, ploidy_label, G0) %>%
  summarise(
    peak_live = max(live_cells, na.rm = TRUE),
    peak_q50 = max(q50, na.rm = TRUE),
    peak_tail_ratio = max(tail_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cellLine, G0, ploidy_label)

knitr::kable(head(size_summary_tbl, 18), digits = 2)

contrast_long <- paired_ploidy %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      delta_live = "Live Cells",
      delta_dead_fraction = "Dead Fraction",
      delta_q50 = "Area q50",
      delta_q90 = "Area q90",
      delta_q99 = "Area q99",
      delta_tail_ratio = "Tail Ratio"
    )
  )

ggplot(contrast_long, aes(hours, value, color = metric)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_grid(G0_label ~ cellLine, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Within-Line Ploidy Contrasts Over Time",
    subtitle = "elevated ploidy minus baseline ploidy",
    x = "Hours",
    y = "Contrast"
  ) +
  theme(legend.position = "bottom")

contrast_summary_tbl <- paired_ploidy %>%
  group_by(cellLine, G0) %>%
  summarise(
    max_abs_delta_live = max(abs(delta_live), na.rm = TRUE),
    max_abs_delta_q50 = max(abs(delta_q50), na.rm = TRUE),
    max_abs_delta_tail_ratio = max(abs(delta_tail_ratio), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_delta_q50))

knitr::kable(head(contrast_summary_tbl, 15), digits = 2)

p_lag_profile <- ggplot(
  lag_summary,
  aes(lag_hours, median_corr, color = ploidy_label)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_grid(G0_label ~ cellLine, scales = "free_y") +
  scale_color_manual(values = c("baseline" = "#0072B2", "elevated" = "#D55E00")) +
  theme_bw() +
  labs(
    title = "Lagged Coupling Between Area Change and Death Rate",
    subtitle = "Positive lag means area-rate signal leads death-rate signal",
    x = "Lag (hours)",
    y = "Median replicate-level correlation",
    color = "Ploidy"
  )

rate_summary_scaled <- rate_summary %>%
  group_by(cellLine, G0, G0_label, signal) %>%
  mutate(value_scaled = value / max(abs(value), na.rm = TRUE)) %>%
  ungroup()

p_rate_signals <- ggplot(
  rate_summary_scaled,
  aes(hours, value_scaled, color = ploidy_label, linetype = signal, group = interaction(ploidy_label, signal, G0))
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.8) +
  facet_grid(G0_label ~ cellLine, scales = "free_y") +
  scale_color_manual(values = c("baseline" = "#0072B2", "elevated" = "#D55E00")) +
  theme_bw() +
  labs(
    title = "Smoothed Signals Used for Lag Analysis",
    subtitle = "Signals are scaled within each cellLine x G0 x signal block to aid visual comparison",
    x = "Hours",
    y = "Scaled signal",
    color = "Ploidy",
    linetype = "Signal"
  )

(p_lag_profile / p_rate_signals) + plot_layout(heights = c(1.1, 1.4))

lag_corr_tbl <- best_lag_tbl %>%
  group_by(cellLine) %>%
  summarise(
    median_best_corr = median(best_abs_corr, na.rm = TRUE),
    median_best_lag_hours = median(best_lag_hours, na.rm = TRUE),
    .groups = "drop"
  )

knitr::kable(lag_corr_tbl, digits = 3)

p_gluc_q50 <- ggplot(
  morph_glucose_summary,
  aes(glucose_hat_plot, q50, color = hours, shape = ploidy_label, group = interaction(ploidy_metric, G0))
) +
  geom_point(alpha = 0.8) +
  geom_line(alpha = 0.7) +
  facet_wrap(~cellLine, scales = "free") +
  scale_x_log10() +
  scale_color_viridis_c(option = "C") +
  theme_bw() +
  labs(
    title = "Median Area Versus Remaining Glucose",
    x = "Estimated glucose",
    y = "Area q50",
    color = "Hours",
    shape = "Ploidy"
  )

p_gluc_lag24 <- ggplot(
  morph_glucose_summary,
  aes(glucose_lag24_plot, q50, color = hours, shape = ploidy_label, group = interaction(ploidy_metric, G0))
) +
  geom_point(alpha = 0.8) +
  geom_line(alpha = 0.7) +
  facet_wrap(~cellLine, scales = "free") +
  scale_x_log10() +
  scale_color_viridis_c(option = "C") +
  theme_bw() +
  labs(
    title = "Median Area Versus Glucose Lagged 24h",
    x = "Estimated glucose 24h earlier",
    y = "Area q50",
    color = "Hours",
    shape = "Ploidy"
  )

p_gluc_lag48 <- ggplot(
  morph_glucose_summary,
  aes(glucose_lag48_plot, q50, color = hours, shape = ploidy_label, group = interaction(ploidy_metric, G0))
) +
  geom_point(alpha = 0.8) +
  geom_line(alpha = 0.7) +
  facet_wrap(~cellLine, scales = "free") +
  scale_x_log10() +
  scale_color_viridis_c(option = "C") +
  theme_bw() +
  labs(
    title = "Median Area Versus Glucose Lagged 48h",
    x = "Estimated glucose 48h earlier",
    y = "Area q50",
    color = "Hours",
    shape = "Ploidy"
  )

(p_gluc_q50 / p_gluc_lag24 / p_gluc_lag48) + plot_layout(heights = c(1.2, 1.2, 1.2))

morph_glucose_assoc <- morph_glucose_summary %>%
  group_by(cellLine) %>%
  summarise(
    cor_q50_vs_glucose = suppressWarnings(cor(log10(glucose_hat_plot), q50, use = "complete.obs")),
    cor_q50_vs_glucose_lag24 = suppressWarnings(cor(log10(glucose_lag24_plot), q50, use = "complete.obs")),
    cor_q50_vs_glucose_lag48 = suppressWarnings(cor(log10(glucose_lag48_plot), q50, use = "complete.obs")),
    .groups = "drop"
  )

knitr::kable(morph_glucose_assoc, digits = 3)

knitr::kable(cost_fits, digits = 3)

p_gluc_init <- ggplot(
  glucose_init_check,
  aes(G0, glucose_initial, color = ploidy_label)
) +
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Estimated Initial Glucose Versus Nominal G0",
    x = "Nominal initial glucose (G0)",
    y = "Estimated glucose at first measurement",
    color = "Ploidy"
  )

p_build <- ggplot(
  cost_fit_diagnostics,
  aes(dB, glucose_consumed, color = ploidy_label)
) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Glucose Consumed Versus Biomass Built",
    x = "dB = net total-cell gain up to last glucose measurement",
    y = "Estimated glucose consumed",
    color = "Ploidy"
  )

p_auc <- ggplot(
  cost_fit_diagnostics,
  aes(AUC_live, glucose_consumed, color = ploidy_label)
) +
  geom_point(alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Glucose Consumed Versus Live-Cell AUC",
    x = "Live-cell AUC up to last glucose measurement",
    y = "Estimated glucose consumed",
    color = "Ploidy"
  )

p_gluc_init / p_build / p_auc

p_cost <- ggplot(cost_fits, aes(factor(ploidy_metric), a, fill = cellLine)) +
  geom_col() +
  facet_wrap(~cellLine, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Estimated Construction Cost by Ploidy",
    subtitle = "Parameter a from glucose_consumed ~ a*dB + m*AUC_live",
    x = "Ploidy metric",
    y = "Construction cost"
  )

scaling_data <- bind_rows(
  cost_fits %>%
    transmute(cellLine, ploidy_metric, assumption = "Flat", volume_ratio = ratio_vol_flat, cost_ratio = ratio_cost),
  cost_fits %>%
    transmute(cellLine, ploidy_metric, assumption = "Sphere", volume_ratio = ratio_vol_sphere, cost_ratio = ratio_cost)
) %>%
  filter(is.finite(volume_ratio), is.finite(cost_ratio), volume_ratio > 0, cost_ratio > 0)

p_scaling <- ggplot(scaling_data, aes(volume_ratio, cost_ratio, color = assumption)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 2.5) +
  geom_line(aes(group = interaction(cellLine, assumption)), alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Does Glucose Cost Scale With Size?",
    subtitle = "Reference line is isometric scaling",
    x = "Volume ratio relative to baseline ploidy",
    y = "Construction-cost ratio relative to baseline ploidy",
    color = "Volume proxy"
  )

p_cost / p_scaling

p_fit_quality <- ggplot(
  cost_fit_diagnostics,
  aes(glucose_consumed, glucose_pred, color = ploidy_label)
) +
  geom_point(alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Observed Versus Predicted Glucose Consumption",
    x = "Observed glucose consumed",
    y = "Predicted glucose consumed",
    color = "Ploidy"
  )

p_resid <- ggplot(
  cost_fit_diagnostics,
  aes(dB, resid, color = ploidy_label)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_point(alpha = 0.85) +
  facet_wrap(~cellLine, scales = "free") +
  theme_bw() +
  labs(
    title = "Cost-Model Residuals Versus Biomass Built",
    x = "dB",
    y = "Observed - predicted glucose consumed",
    color = "Ploidy"
  )

p_fit_quality / p_resid

cellline_timecourse_plots <- cost_fit_timecourse_summary %>%
  split(.$cellLine) %>%
  lapply(function(df) {
    ggplot(
      df,
      aes(hours, value, color = series, linetype = series)
    ) +
      geom_line(linewidth = 0.8) +
      facet_grid(G0_label ~ ploidy_label, scales = "free_y") +
      theme_bw() +
      labs(
        title = unique(df$cellLine),
        x = "Hours",
        y = "Total cells",
        color = NULL,
        linetype = NULL
      )
  })

for (p in cellline_timecourse_plots) {
  print(p)
}

cost_fit_summary_tbl <- cost_fit_timecourse_summary %>%
  pivot_wider(names_from = series, values_from = value) %>%
  group_by(cellLine, ploidy_label, G0) %>%
  summarise(
    rmse_timecourse = sqrt(mean((`Observed total cells` - `Predicted total cells`)^2, na.rm = TRUE)),
    cor_timecourse = suppressWarnings(cor(`Observed total cells`, `Predicted total cells`, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  arrange(cellLine, G0, ploidy_label)

knitr::kable(head(cost_fit_summary_tbl, 20), digits = 3)

efficiency_long <- cost_fits %>%
  select(cellLine, ploidy_metric, cost_per_flat_vol, cost_per_sphere_vol) %>%
  pivot_longer(
    cols = starts_with("cost_per_"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      cost_per_flat_vol = "Cost per flat-volume proxy",
      cost_per_sphere_vol = "Cost per sphere-volume proxy"
    )
  )

ggplot(efficiency_long, aes(factor(ploidy_metric), value, fill = factor(ploidy_metric))) +
  geom_col(show.legend = FALSE) +
  facet_grid(metric ~ cellLine, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Glucose Cost Per Unit Size Proxy",
    x = "Ploidy metric",
    y = "Cost / size"
  )
