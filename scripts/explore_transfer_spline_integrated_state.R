suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(mgcv)
})

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/model_free_ploidy_utils.R")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

cumtrapz_vec <- function(x, y) {
  n <- length(x)
  out <- rep(0, n)
  if (n < 2L) {
    return(out)
  }
  for (i in 2:n) {
    dx <- x[i] - x[i - 1L]
    out[i] <- out[i - 1L] + dx * (y[i] + y[i - 1L]) / 2
  }
  out
}

lag1_resid_cor <- function(df, resid_col = "resid") {
  vals <- df[[resid_col]]
  if (length(vals) < 3L || sd(vals, na.rm = TRUE) == 0 || any(!is.finite(vals))) {
    return(NA_real_)
  }
  stats::cor(vals[-length(vals)], vals[-1L])
}

fit_positive_spline <- function(df, spar = 0.72, min_unique_x = 4L, n_eval = 200L) {
  df <- df %>%
    filter(is.finite(hours), is.finite(value)) %>%
    arrange(hours)

  unique_x <- sort(unique(df$hours))
  if (nrow(df) < min_unique_x || length(unique_x) < min_unique_x) {
    return(NULL)
  }

  fit <- try(
    stats::smooth.spline(
      x = df$hours,
      y = log1p(pmax(df$value, 0)),
      spar = spar,
      cv = FALSE
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  eval_hours <- seq(min(df$hours), max(df$hours), length.out = n_eval)
  pred <- predict(fit, x = eval_hours)
  deriv <- predict(fit, x = eval_hours, deriv = 1)

  fit_value <- pmax(expm1(pred$y), 0)
  fit_deriv <- pmax(exp(pred$y), 1) * deriv$y

  tibble(
    hours = eval_hours,
    fit_value = fit_value,
    fit_deriv = fit_deriv
  )
}

fit_panel_splines <- function(df, spar_map = c(N = 0.45, R1 = 0.70)) {
  nested <- df %>%
    group_by(cellLine, line_id, ploidy_abs, ploidy_metric, G0, variable, well_idx) %>%
    nest()

  fits <- nested %>%
    mutate(
      fit_spar = unname(spar_map[as.character(variable)]),
      fit_spar = ifelse(is.na(fit_spar), 0.65, fit_spar),
      fit = map2(data, fit_spar, fit_positive_spline),
      fit_ok = map_lgl(fit, ~ !is.null(.x))
    )

  fits %>%
    filter(fit_ok) %>%
    select(-data, -fit_ok, -fit_spar) %>%
    unnest(fit)
}

build_rate_relationships <- function(fit_df, min_n = 1) {
  n_fit <- fit_df %>%
    filter(variable == "N") %>%
    transmute(
      cellLine,
      line_id,
      ploidy_abs,
      ploidy_metric,
      G0,
      well_idx,
      hours,
      N_fit = fit_value,
      dN_dt = fit_deriv
    )

  r1_fit <- fit_df %>%
    filter(variable == "R1") %>%
    transmute(
      cellLine,
      line_id,
      ploidy_abs,
      ploidy_metric,
      G0,
      well_idx,
      hours,
      R1_fit = fit_value,
      dR1_dt = fit_deriv
    )

  inner_join(
    n_fit,
    r1_fit,
    by = c("cellLine", "line_id", "ploidy_abs", "ploidy_metric", "G0", "well_idx", "hours")
  ) %>%
    filter(is.finite(N_fit), is.finite(R1_fit), N_fit >= min_n) %>%
    mutate(
      uptake_per_cell = dR1_dt / pmax(N_fit, min_n),
      growth_per_cell = dN_dt / pmax(N_fit, min_n)
    )
}

fit_rate_model <- function(df, response, model_formula, model_name) {
  dat <- df %>%
    ungroup() %>%
    transmute(
      y = .data[[response]],
      log_R1,
      log_N,
      hours,
      log_G0,
      gluc_drawdown,
      frac_gluc_remaining,
      lag_log_R1,
      lag_log_N,
      cum_drawdown_auc,
      cum_N_auc,
      cellLine,
      ploidy_abs,
      well_idx
    ) %>%
    filter(
      is.finite(y),
      is.finite(log_R1),
      is.finite(log_N),
      is.finite(hours),
      is.finite(log_G0),
      is.finite(gluc_drawdown),
      is.finite(frac_gluc_remaining),
      is.finite(lag_log_R1),
      is.finite(lag_log_N),
      is.finite(cum_drawdown_auc),
      is.finite(cum_N_auc)
    )

  fit <- mgcv::gam(model_formula, data = dat, method = "REML")
  resid <- residuals(fit, type = "response")
  pred <- fitted(fit)

  acf_by_well <- dat %>%
    mutate(resid = resid) %>%
    group_by(well_idx) %>%
    summarise(acf1 = lag1_resid_cor(pick(everything()), "resid"), .groups = "drop")

  tibble(
    model_name = model_name,
    rmse = sqrt(mean((dat$y - pred)^2)),
    mae = mean(abs(dat$y - pred)),
    dev_expl = summary(fit)$dev.expl %||% NA_real_,
    r_sq = summary(fit)$r.sq %||% NA_real_,
    aic = AIC(fit),
    median_abs_acf1 = median(abs(acf_by_well$acf1), na.rm = TRUE),
    mean_abs_acf1 = mean(abs(acf_by_well$acf1), na.rm = TRUE),
    n = nrow(dat),
    n_wells = n_distinct(dat$well_idx)
  )
}

model_formulas <- list(
  baseline = y ~ s(log_R1, k = 6),
  plus_logN = y ~ s(log_R1, k = 6) + s(log_N, k = 6),
  plus_time = y ~ s(log_R1, k = 6) + s(hours, k = 6),
  plus_logG0 = y ~ s(log_R1, k = 6) + s(log_G0, k = 6),
  plus_drawdown = y ~ s(log_R1, k = 6) + s(gluc_drawdown, k = 6),
  plus_frac_remaining = y ~ s(log_R1, k = 6) + s(frac_gluc_remaining, k = 6),
  plus_lag_log_R1 = y ~ s(log_R1, k = 6) + s(lag_log_R1, k = 6),
  plus_lag_log_N = y ~ s(log_R1, k = 6) + s(lag_log_N, k = 6),
  plus_cum_drawdown_auc = y ~ s(log_R1, k = 6) + s(cum_drawdown_auc, k = 6),
  plus_cum_N_auc = y ~ s(log_R1, k = 6) + s(cum_N_auc, k = 6),
  plus_integrated_state = y ~ s(log_R1, k = 6) + s(log_N, k = 6) + s(cum_drawdown_auc, k = 6) + s(cum_N_auc, k = 6),
  plus_logN_time = y ~ s(log_R1, k = 6) + s(log_N, k = 6) + s(hours, k = 6)
)

message("Building model-free tables and uncensored well set...")
tables <- build_model_free_tables()

excluded_censored_wells <- tables$glucose_obs %>%
  group_by(well_idx) %>%
  summarise(any_censored = any(is_censored == 1L), .groups = "drop") %>%
  filter(any_censored) %>%
  pull(well_idx)

included_wells <- tables$well_meta %>%
  filter(!well_idx %in% excluded_censored_wells) %>%
  pull(well_idx)

count_df <- tables$count_obs %>%
  filter(well_idx %in% included_wells) %>%
  left_join(
    tables$well_meta %>% select(well_idx, cellLine, line_id, ploidy_abs, ploidy_metric, G0),
    by = "well_idx"
  ) %>%
  transmute(cellLine, line_id, ploidy_abs, ploidy_metric, G0, well_idx, hours, variable = "N", value = total_cells)

glucose_df <- tables$glucose_summary %>%
  filter(well_idx %in% included_wells, any_glucose_censored == 0L) %>%
  transmute(cellLine, line_id, ploidy_abs, ploidy_metric, G0, well_idx, hours, variable = "R1", value = glucose_hat)

obs_df <- bind_rows(count_df, glucose_df) %>%
  mutate(value = pmax(value, 0))

message("Fitting splines and building rate table...")
fit_df <- fit_panel_splines(obs_df)
rate_df <- build_rate_relationships(fit_df, min_n = 1) %>%
  mutate(
    log_R1 = log10(pmax(R1_fit, 1e-6)),
    log_N = log10(pmax(N_fit, 1)),
    log_G0 = log10(pmax(G0, 1e-4)),
    gluc_drawdown = pmax(G0 - R1_fit, 0),
    frac_gluc_remaining = if_else(G0 > 0, pmin(pmax(R1_fit / G0, 0), 5), NA_real_)
  ) %>%
  group_by(cellLine, ploidy_abs, G0, well_idx) %>%
  arrange(hours, .by_group = TRUE) %>%
  mutate(
    lag_log_R1 = lag(log_R1),
    lag_log_N = lag(log_N),
    cum_drawdown_auc = cumtrapz_vec(hours, gluc_drawdown),
    cum_N_auc = cumtrapz_vec(hours, N_fit)
  ) %>%
  ungroup()

group_keys <- rate_df %>%
  distinct(cellLine, ploidy_abs) %>%
  arrange(cellLine, ploidy_abs)

responses <- c("uptake_per_cell", "growth_per_cell")

message("Running grouped GAM screen...")
model_fit_tbl <- bind_rows(lapply(seq_len(nrow(group_keys)), function(i) {
  cell_line <- group_keys$cellLine[[i]]
  ploidy_value <- group_keys$ploidy_abs[[i]]
  df_group <- rate_df %>%
    filter(cellLine == cell_line, ploidy_abs == ploidy_value)

  bind_rows(lapply(responses, function(response_name) {
    bind_rows(lapply(names(model_formulas), function(model_name) {
      fit_rate_model(df_group, response_name, model_formulas[[model_name]], model_name) %>%
        mutate(
          cellLine = cell_line,
          ploidy_abs = ploidy_value,
          response = response_name,
          .before = 1
        )
    }))
  }))
}))

summary_by_response <- model_fit_tbl %>%
  group_by(response, model_name) %>%
  summarise(
    n_groups = n(),
    median_rmse = median(rmse, na.rm = TRUE),
    median_mae = median(mae, na.rm = TRUE),
    median_dev_expl = median(dev_expl, na.rm = TRUE),
    median_r_sq = median(r_sq, na.rm = TRUE),
    median_abs_acf1 = median(median_abs_acf1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(response, median_rmse)

single_model_names <- c(
  "plus_logN", "plus_time", "plus_logG0", "plus_drawdown", "plus_frac_remaining",
  "plus_lag_log_R1", "plus_lag_log_N", "plus_cum_drawdown_auc", "plus_cum_N_auc"
)

single_var_summary <- summary_by_response %>%
  filter(model_name %in% single_model_names) %>%
  group_by(response) %>%
  arrange(median_rmse, .by_group = TRUE) %>%
  mutate(rank_rmse = row_number()) %>%
  ungroup()

single_var_wins <- model_fit_tbl %>%
  filter(model_name %in% single_model_names) %>%
  group_by(response, cellLine, ploidy_abs) %>%
  arrange(rmse, desc(dev_expl), median_abs_acf1, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  count(response, model_name, sort = TRUE)

pairwise_vs_time <- model_fit_tbl %>%
  filter(model_name %in% c("plus_time", "plus_cum_drawdown_auc", "plus_cum_N_auc", "plus_integrated_state", "plus_logN_time")) %>%
  select(cellLine, ploidy_abs, response, model_name, rmse, dev_expl, median_abs_acf1) %>%
  pivot_wider(names_from = model_name, values_from = c(rmse, dev_expl, median_abs_acf1), values_fn = dplyr::first) %>%
  transmute(
    cellLine,
    ploidy_abs,
    response,
    rmse_gain_cum_drawdown_vs_time = rmse_plus_time - rmse_plus_cum_drawdown_auc,
    rmse_gain_cum_N_vs_time = rmse_plus_time - rmse_plus_cum_N_auc,
    rmse_gain_integrated_vs_time = rmse_plus_time - rmse_plus_integrated_state,
    rmse_gain_logN_time_vs_integrated = rmse_plus_integrated_state - rmse_plus_logN_time,
    dev_gain_integrated_vs_time = dev_expl_plus_integrated_state - dev_expl_plus_time,
    acf_gain_integrated_vs_time = median_abs_acf1_plus_time - median_abs_acf1_plus_integrated_state
  )

pairwise_summary <- pairwise_vs_time %>%
  group_by(response) %>%
  summarise(
    median_rmse_gain_cum_drawdown_vs_time = median(rmse_gain_cum_drawdown_vs_time, na.rm = TRUE),
    median_rmse_gain_cum_N_vs_time = median(rmse_gain_cum_N_vs_time, na.rm = TRUE),
    median_rmse_gain_integrated_vs_time = median(rmse_gain_integrated_vs_time, na.rm = TRUE),
    median_rmse_gain_logN_time_vs_integrated = median(rmse_gain_logN_time_vs_integrated, na.rm = TRUE),
    n_groups_integrated_beats_time = sum(rmse_gain_integrated_vs_time > 0, na.rm = TRUE),
    n_groups_integrated_beats_logN_time = sum(rmse_gain_logN_time_vs_integrated > 0, na.rm = TRUE),
    median_dev_gain_integrated_vs_time = median(dev_gain_integrated_vs_time, na.rm = TRUE),
    median_acf_gain_integrated_vs_time = median(acf_gain_integrated_vs_time, na.rm = TRUE),
    .groups = "drop"
  )

out_dir <- file.path("data", "report_exports", "transfer_spline_integrated_state")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(model_fit_tbl, file.path(out_dir, "model_fit_by_group.csv"), row.names = FALSE)
write.csv(summary_by_response, file.path(out_dir, "model_summary_by_response.csv"), row.names = FALSE)
write.csv(single_var_summary, file.path(out_dir, "single_variable_summary.csv"), row.names = FALSE)
write.csv(single_var_wins, file.path(out_dir, "single_variable_wins_by_group.csv"), row.names = FALSE)
write.csv(pairwise_vs_time, file.path(out_dir, "pairwise_vs_time_by_group.csv"), row.names = FALSE)
write.csv(pairwise_summary, file.path(out_dir, "pairwise_vs_time_summary.csv"), row.names = FALSE)

message("Wrote outputs to: ", out_dir)
print(summary_by_response)
print(single_var_summary)
print(single_var_wins)
print(pairwise_summary)
