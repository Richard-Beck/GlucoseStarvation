#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: build_model_free_feature_dataset.R STAN_DATA_RDS OUTPUT_RDS",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

stan_data_recorded_path <- args[[1]]
stan_data_path <- normalizePath(stan_data_recorded_path, mustWork = TRUE)
output_path <- args[[2]]
if (file.exists(output_path)) {
  stop("Refusing to overwrite existing output: ", output_path, call. = FALSE)
}

utility_paths <- normalizePath(
  c(
    "R/model_free_ploidy_utils.R",
    "R/project_paths.R",
    "R/gpath_run_utils.R"
  ),
  mustWork = TRUE
)
source(utility_paths[[1]])

glucose_floor <- 0.1
min_glucose_drawdown_for_yield <- 0.05
growth_low_g0_max <- 0.25
growth_high_g0 <- c(5, 25)
live_auc_fit_g0_max <- 1
glucose_per_live_auc_g0_max <- 1
fixed_yield_g0 <- 1
confidence_level <- 0.95
robust_derivative_spar <- 0.62
robust_derivative_n_eval <- 260L
robust_derivative_edge_fraction <- 0.06
robust_derivative_quantile <- 0.90

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, call. = FALSE)
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_args)) {
  normalizePath(sub("^--file=", "", script_args[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

finite_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

count_finite <- function(x) {
  sum(is.finite(as.numeric(x)))
}

compute_robust_derivative <- function(
  df,
  spar = 0.62,
  n_eval = 260L,
  edge_fraction = 0.06,
  derivative_quantile = 0.90
) {
  df <- df %>% arrange(hours)
  keep <- is.finite(df$hours) & is.finite(df$live_cells)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 2L) {
    stop("Robust derivative requires at least two finite trajectory points.", call. = FALSE)
  }

  eval_hours <- seq(min(df$hours), max(df$hours), length.out = n_eval)
  log_live <- log1p(pmax(df$live_cells, 0))
  fit <- if (length(unique(df$hours)) >= 4L) {
    try(stats::smooth.spline(df$hours, log_live, spar = spar), silent = TRUE)
  } else {
    structure(list(), class = "try-error")
  }

  derivative <- if (inherits(fit, "try-error")) {
    smooth_log_live <- stats::approx(
      df$hours,
      log_live,
      xout = eval_hours,
      ties = "ordered",
      rule = 2
    )$y
    c(diff(smooth_log_live) / diff(eval_hours), NA_real_)
  } else {
    stats::predict(fit, x = eval_hours, deriv = 1)$y
  }

  span <- diff(range(eval_hours))
  inner <- eval_hours >= min(eval_hours) + edge_fraction * span &
    eval_hours <= max(eval_hours) - edge_fraction * span
  finite_inner <- is.finite(derivative) & inner
  if (!any(finite_inner)) {
    stop("Robust derivative has no finite interior evaluation points.", call. = FALSE)
  }

  robust_value <- unname(stats::quantile(
    derivative[finite_inner],
    probs = derivative_quantile,
    na.rm = TRUE,
    names = FALSE
  ))
  top_quantile <- finite_inner & derivative >= robust_value
  robust_time <- if (any(top_quantile)) {
    stats::median(eval_hours[top_quantile])
  } else {
    NA_real_
  }

  tibble(
    well_idx = df$well_idx[[1]],
    robust_max_derivative = robust_value,
    robust_derivative_time = robust_time,
    robust_derivative_n_points = sum(finite_inner),
    robust_derivative_n_top_quantile_points = sum(top_quantile),
    robust_derivative_rule = paste0(
      derivative_quantile * 100,
      "th percentile of interior smooth-spline derivative"
    )
  )
}

fit_diagnostics <- function(
  data,
  fit_id,
  response,
  predictor,
  formula_text,
  feature_term,
  g0_max,
  include_intercept = TRUE,
  conf_level = 0.95
) {
  x <- as.numeric(data[[predictor]])
  y <- as.numeric(data[[response]])
  keep <- is.finite(x) & is.finite(y)
  fit_data <- data.frame(x = x[keep], y = y[keep])

  empty <- tibble(
    fit_id = fit_id,
    formula = formula_text,
    response = response,
    predictor = predictor,
    feature_term = feature_term,
    include_intercept = include_intercept,
    g0_max = g0_max,
    confidence_level = conf_level,
    intercept = NA_real_,
    intercept_se = NA_real_,
    intercept_conf_low = NA_real_,
    intercept_conf_high = NA_real_,
    slope = NA_real_,
    slope_se = NA_real_,
    slope_conf_low = NA_real_,
    slope_conf_high = NA_real_,
    r_squared = NA_real_,
    adjusted_r_squared = NA_real_,
    residual_sigma = NA_real_,
    residual_df = NA_integer_,
    n = nrow(fit_data),
    predictor_min = if (nrow(fit_data)) min(fit_data$x) else NA_real_,
    predictor_max = if (nrow(fit_data)) max(fit_data$x) else NA_real_,
    response_min = if (nrow(fit_data)) min(fit_data$y) else NA_real_,
    response_max = if (nrow(fit_data)) max(fit_data$y) else NA_real_,
    n_negative_response = if (nrow(fit_data)) sum(fit_data$y < 0) else 0L
  )

  if (nrow(fit_data) < 3L || length(unique(fit_data$x)) < 2L) {
    return(empty)
  }

  fit <- if (include_intercept) {
    stats::lm(y ~ x, data = fit_data)
  } else {
    stats::lm(y ~ 0 + x, data = fit_data)
  }
  fit_summary <- summary(fit)
  coefficients <- fit_summary$coefficients
  slope_row <- coefficients["x", , drop = FALSE]
  slope_ci <- stats::confint(fit, "x", level = conf_level)

  intercept <- 0
  intercept_se <- NA_real_
  intercept_ci <- c(NA_real_, NA_real_)
  if (include_intercept) {
    intercept <- unname(stats::coef(fit)["(Intercept)"])
    intercept_se <- unname(coefficients["(Intercept)", "Std. Error"])
    intercept_ci <- as.numeric(stats::confint(fit, "(Intercept)", level = conf_level))
  }

  empty %>%
    mutate(
      intercept = .env$intercept,
      intercept_se = .env$intercept_se,
      intercept_conf_low = .env$intercept_ci[[1]],
      intercept_conf_high = .env$intercept_ci[[2]],
      slope = unname(stats::coef(fit)["x"]),
      slope_se = unname(slope_row[1, "Std. Error"]),
      slope_conf_low = as.numeric(slope_ci[[1]]),
      slope_conf_high = as.numeric(slope_ci[[2]]),
      r_squared = fit_summary$r.squared,
      adjusted_r_squared = fit_summary$adj.r.squared,
      residual_sigma = fit_summary$sigma,
      residual_df = as.integer(fit$df.residual)
    )
}

message("Building condition-level empirical features from ", stan_data_path)
built <- build_feature_panel(
  stan_data_path = stan_data_path,
  glucose_floor = glucose_floor,
  min_glucose_drawdown_for_yield = min_glucose_drawdown_for_yield
)

robust_derivatives <- built$count_summary %>%
  split(.$well_idx) %>%
  map_dfr(
    compute_robust_derivative,
    spar = robust_derivative_spar,
    n_eval = robust_derivative_n_eval,
    edge_fraction = robust_derivative_edge_fraction,
    derivative_quantile = robust_derivative_quantile
  )

condition_features <- built$feature_panel %>%
  left_join(robust_derivatives, by = "well_idx") %>%
  mutate(
    peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)
  ) %>%
  arrange(line_id, ploidy_metric, exp_id, G0, well_idx)

required_condition_columns <- c(
  "well_idx", "cellLine", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id",
  "max_growth_rate", "robust_max_derivative", "exp_growth_rate_70", "live_auc_glucose_window",
  "glucose_initial", "glucose_final", "glucose_drawdown",
  "total_net_gain_to_glucose_end", "peak_total_yield_net"
)
missing_condition_columns <- setdiff(required_condition_columns, names(condition_features))
if (length(missing_condition_columns)) {
  stop(
    "Condition feature table lacks required columns: ",
    paste(missing_condition_columns, collapse = ", "),
    call. = FALSE
  )
}
if (!nrow(condition_features) || anyDuplicated(condition_features$well_idx)) {
  stop("Condition feature table must contain one unique row per Stan well.", call. = FALSE)
}
drawdown_expected <- condition_features$glucose_initial - condition_features$glucose_final
drawdown_matches <- dplyr::near(condition_features$glucose_drawdown, drawdown_expected)
if (any(is.na(drawdown_matches)) || any(!drawdown_matches)) {
  stop("Glucose drawdown identity check failed.", call. = FALSE)
}

group_columns <- c("cellLine", "line_id", "ploidy_metric", "ploidy_abs")

feature_fits <- condition_features %>%
  group_by(across(all_of(group_columns))) %>%
  group_modify(function(df, key) {
    low_g <- df %>%
      filter(G0 <= live_auc_fit_g0_max) %>%
      mutate(log1p_G0 = log1p(G0))

    bind_rows(
      fit_diagnostics(
        data = low_g,
        fit_id = "yield_alive_auc_logG",
        response = "live_auc_glucose_window",
        predictor = "log1p_G0",
        formula_text = "live_auc_glucose_window ~ log1p(G0)",
        feature_term = "intercept",
        g0_max = live_auc_fit_g0_max,
        include_intercept = TRUE,
        conf_level = confidence_level
      ),
      fit_diagnostics(
        data = low_g,
        fit_id = "glucose_per_live_auc",
        response = "glucose_drawdown",
        predictor = "live_auc_glucose_window",
        formula_text = "glucose_drawdown ~ live_auc_glucose_window",
        feature_term = "slope",
        g0_max = glucose_per_live_auc_g0_max,
        include_intercept = TRUE,
        conf_level = confidence_level
      )
    )
  }) %>%
  ungroup() %>%
  arrange(line_id, ploidy_metric, fit_id)

legacy_signatures <- summarize_glucose_signatures(condition_features)

direct_features <- condition_features %>%
  group_by(across(all_of(group_columns))) %>%
  summarise(
    n_conditions = n(),
    n_experiments = n_distinct(exp_id),
    growth_lowG_median = finite_median(robust_max_derivative[G0 <= growth_low_g0_max]),
    live_auc_0mM = finite_median(live_auc_glucose_window[abs(G0) < 1e-8]),
    peak_total_yield_1mM = finite_median(peak_total_yield_net[abs(G0 - fixed_yield_g0) < 1e-8]),
    n_growth_lowG = count_finite(robust_max_derivative[G0 <= growth_low_g0_max]),
    n_growth_highG = count_finite(exp_growth_rate_70[G0 %in% growth_high_g0]),
    n_live_auc_0mM = count_finite(live_auc_glucose_window[abs(G0) < 1e-8]),
    n_peak_total_yield_1mM = count_finite(peak_total_yield_net[abs(G0 - fixed_yield_g0) < 1e-8]),
    .groups = "drop"
  )

fit_features <- feature_fits %>%
  transmute(
    across(all_of(group_columns)),
    fit_id,
    feature_value = if_else(feature_term == "intercept", intercept, slope),
    n
  ) %>%
  pivot_wider(
    names_from = fit_id,
    values_from = c(feature_value, n),
    names_glue = "{.value}_{fit_id}"
  ) %>%
  rename(
    yield_alive_auc_intercept = feature_value_yield_alive_auc_logG,
    glucose_per_live_auc_slope = feature_value_glucose_per_live_auc,
    n_yield_alive_auc_fit = n_yield_alive_auc_logG,
    n_glucose_per_live_auc_fit = n_glucose_per_live_auc
  )

line_ploidy_features <- legacy_signatures %>%
  select(
    all_of(group_columns), has_starvation,
    growth_highG_median,
    legacy_yield_alive_auc_intercept = yield_alive_auc_intercept
  ) %>%
  left_join(direct_features, by = group_columns) %>%
  left_join(fit_features, by = group_columns) %>%
  mutate(
    auc_intercept_check = abs(yield_alive_auc_intercept - legacy_yield_alive_auc_intercept)
  ) %>%
  arrange(line_id, ploidy_metric)

if (any(!is.finite(line_ploidy_features$auc_intercept_check)) ||
    any(line_ploidy_features$auc_intercept_check > 1e-8)) {
  stop("Live-AUC intercept does not match the maintained utility calculation.", call. = FALSE)
}
line_ploidy_features <- line_ploidy_features %>%
  select(-legacy_yield_alive_auc_intercept, -auc_intercept_check)

feature_ids <- c(
  "growth_lowG_median",
  "growth_highG_median",
  "live_auc_0mM",
  "yield_alive_auc_intercept",
  "glucose_per_live_auc_slope",
  "peak_total_yield_1mM"
)

feature_catalog <- tibble(
  feature_id = feature_ids,
  feature_class = c("Growth", "Growth", "Alive AUC", "Alive AUC", "Glucose use", "Total-cell yield"),
  display_label = c(
    "Low-glucose growth",
    "High-glucose exponential growth",
    "Live-cell AUC at 0 mM",
    "Live-cell AUC regression intercept",
    "Glucose drawdown per live-cell AUC",
    "Net total-cell yield at 1 mM"
  ),
  units = c(
    "log1p(cells) per hour",
    "log1p(cells) per hour",
    "live-cell hours",
    "live-cell hours",
    "glucose concentration per live-cell hour",
    "cells"
  ),
  condition_source = c(
    "robust_max_derivative",
    "exp_growth_rate_70",
    "live_auc_glucose_window",
    "live_auc_glucose_window",
    "glucose_drawdown and live_auc_glucose_window",
    "peak_total_yield_net"
  ),
  definition = c(
    "Median condition-level robust derivative over G0 <= 0.25 mM. Each robust derivative is the 90th percentile of the finite first derivative after smoothing log1p(live cells) with spar = 0.62 on 260 points and excluding the outer 6% of time.",
    "Median condition-level exponential-phase slope over G0 in {5, 25} mM; each slope uses observations through the first time live cells reach 70% of their observed maximum.",
    "Median trapezoidal live-cell AUC within the glucose-measurement window among G0 = 0 mM conditions.",
    "Intercept from an intercept-bearing regression of live-cell AUC on log1p(G0), using finite conditions with G0 <= 1 mM.",
    "Slope from an intercept-bearing regression of glucose drawdown on live-cell AUC, using finite conditions with G0 <= 1 mM.",
    "Median nonnegative peak-minus-initial total-cell count within the glucose-measurement window among G0 = 1 mM conditions."
  ),
  larger_value_means = c(
    "faster maximum growth under low glucose",
    "faster exponential-phase growth under high glucose",
    "greater cumulative live-cell persistence without supplied glucose",
    "greater fitted baseline cumulative live-cell persistence",
    "greater culture-level glucose-concentration drawdown per live-cell hour",
    "greater net total-cell production at 1 mM starting glucose"
  )
)

ploidy_effects <- line_ploidy_features %>%
  select(cellLine, line_id, ploidy_metric, all_of(feature_ids)) %>%
  pivot_longer(
    cols = all_of(feature_ids),
    names_to = "feature_id",
    values_to = "value"
  ) %>%
  group_by(cellLine, line_id, feature_id) %>%
  arrange(ploidy_metric, .by_group = TRUE) %>%
  summarise(
    n_ploidy_states = n(),
    low_ploidy = first(ploidy_metric),
    high_ploidy = last(ploidy_metric),
    delta_ploidy = high_ploidy - low_ploidy,
    low_value = first(value),
    high_value = last(value),
    delta_value = high_value - low_value,
    effect_per_ploidy = delta_value / delta_ploidy,
    .groups = "drop"
  ) %>%
  left_join(feature_catalog %>% select(feature_id, feature_class, display_label, units), by = "feature_id") %>%
  arrange(line_id, match(feature_id, feature_ids))

if (any(ploidy_effects$n_ploidy_states != 2L) ||
    any(!is.finite(ploidy_effects$delta_ploidy)) ||
    any(ploidy_effects$delta_ploidy <= 0)) {
  stop("Each cell line must have exactly two ordered, distinct ploidy states.", call. = FALSE)
}

feature_matrix <- as.matrix(line_ploidy_features[, feature_ids, drop = FALSE])
if (any(!is.finite(feature_matrix))) {
  stop("Manuscript-facing line/ploidy feature values must all be finite.", call. = FALSE)
}
if (nrow(feature_fits) != 2L * nrow(line_ploidy_features) ||
    any(feature_fits$n < 3L) ||
    any(!is.finite(feature_fits$intercept)) ||
    any(!is.finite(feature_fits$slope))) {
  stop("Feature-fit coverage or coefficient validation failed.", call. = FALSE)
}
if (nrow(ploidy_effects) != length(feature_ids) * n_distinct(line_ploidy_features$cellLine)) {
  stop("Ploidy-effect table has an unexpected number of rows.", call. = FALSE)
}

stan_data <- readRDS(stan_data_path)
git_commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]],
  error = function(e) NA_character_
)

code_paths <- c(script_path, utility_paths)
code_paths <- code_paths[!is.na(code_paths)]
settings <- list(
  schema_version = 1L,
  generated_at = format(Sys.time(), usetz = TRUE),
  description = "Model-free empirical trajectory features used by manuscript Figures 2, S3, and S4.",
  source = list(
    stan_data_path = stan_data_recorded_path,
    stan_data_sha256 = sha256_file(stan_data_path),
    n_wells = as.integer(stan_data$N_wells),
    n_count_observations = length(stan_data$N_obs),
    n_glucose_observations = length(stan_data$lum_obs)
  ),
  code = tibble(
    path = code_paths,
    sha256 = vapply(code_paths, sha256_file, character(1))
  ),
  parameters = list(
    maintained_condition_rate_spline_spar = 0.55,
    maintained_condition_rate_derivative_evaluation_points = 200L,
    robust_derivative_spar = robust_derivative_spar,
    robust_derivative_evaluation_points = robust_derivative_n_eval,
    robust_derivative_excluded_edge_fraction = robust_derivative_edge_fraction,
    robust_derivative_quantile = robust_derivative_quantile,
    exponential_phase_fraction_of_observed_maximum = 0.70,
    glucose_floor = glucose_floor,
    min_glucose_drawdown_for_legacy_yield_ratio = min_glucose_drawdown_for_yield,
    growth_low_g0_max = growth_low_g0_max,
    growth_high_g0 = growth_high_g0,
    live_auc_fit_g0_max = live_auc_fit_g0_max,
    glucose_per_live_auc_g0_max = glucose_per_live_auc_g0_max,
    glucose_per_live_auc_include_intercept = TRUE,
    fixed_yield_g0 = fixed_yield_g0,
    confidence_level = confidence_level
  ),
  git_commit = git_commit,
  r_version = R.version.string,
  package_versions = c(
    dplyr = as.character(utils::packageVersion("dplyr")),
    purrr = as.character(utils::packageVersion("purrr")),
    tibble = as.character(utils::packageVersion("tibble")),
    tidyr = as.character(utils::packageVersion("tidyr"))
  ),
  row_counts = list(
    condition_features = nrow(condition_features),
    line_ploidy_features = nrow(line_ploidy_features),
    feature_fits = nrow(feature_fits),
    ploidy_effects = nrow(ploidy_effects),
    feature_catalog = nrow(feature_catalog)
  )
)

result <- list(
  condition_features = as.data.frame(condition_features),
  line_ploidy_features = as.data.frame(line_ploidy_features),
  feature_fits = as.data.frame(feature_fits),
  ploidy_effects = as.data.frame(ploidy_effects),
  feature_catalog = as.data.frame(feature_catalog),
  settings = settings
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_path, compress = "xz")

round_trip <- readRDS(output_path)
expected_names <- c(
  "condition_features", "line_ploidy_features", "feature_fits",
  "ploidy_effects", "feature_catalog", "settings"
)
if (!identical(names(round_trip), expected_names)) {
  stop("Output RDS failed list-schema validation.", call. = FALSE)
}
if (!identical(round_trip$settings$row_counts, settings$row_counts) ||
    !identical(round_trip$settings$source$stan_data_sha256, settings$source$stan_data_sha256)) {
  stop("Output RDS failed round-trip provenance validation.", call. = FALSE)
}

message("Wrote ", output_path)
message("  condition features: ", nrow(condition_features))
message("  line/ploidy features: ", nrow(line_ploidy_features))
message("  feature fits: ", nrow(feature_fits))
message("  ploidy effects: ", nrow(ploidy_effects))
