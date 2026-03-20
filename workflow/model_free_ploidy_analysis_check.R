## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE
)
knitr::opts_knit$set(root.dir = normalizePath("../"))


## ----libraries, include=FALSE-------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)
library(forcats)
library(purrr)

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/model_free_ploidy_utils.R")

analysis_root <- file.path("data", "model_free_ploidy")

paths <- list(
  stan_data = resolve_stan_data_path(),
  feature_panel = file.path(analysis_root, "feature_panel.Rds"),
  count_summary = file.path(analysis_root, "count_summary.Rds"),
  glucose_summary = file.path(analysis_root, "glucose_summary.Rds"),
  agent_text = file.path("workflow", "model_free_ploidy_analysis.txt")
)

safe_read_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

feature_panel <- safe_read_rds(paths$feature_panel)
count_summary <- safe_read_rds(paths$count_summary)
glucose_summary <- safe_read_rds(paths$glucose_summary)
feature_catalog <- get_model_free_feature_catalog()

if (is.null(feature_panel) && !file.exists(paths$stan_data)) {
  stop("Neither feature_panel nor stan_data is available. Run the extraction script or provide stan_ready_data before rendering this notebook.")
}

feature_panel_required_cols <- c(
  "exp_growth_rate_70",
  "exp_growth_fit_end_time_70",
  "exp_growth_fit_n_points_70",
  "exp_growth_fit_r2_70",
  "max_growth_rate_raw",
  "max_growth_rate",
  "max_death_rate_raw",
  "max_death_rate"
)

need_rebuild_feature_panel <- is.null(feature_panel) ||
  !all(feature_panel_required_cols %in% names(feature_panel))

if (need_rebuild_feature_panel || is.null(count_summary) || is.null(glucose_summary)) {
  rebuilt <- build_feature_panel(paths$stan_data)
  if (need_rebuild_feature_panel) {
    feature_panel <- rebuilt$feature_panel
  }
  if (is.null(count_summary)) {
    count_summary <- rebuilt$count_summary
  }
  if (is.null(glucose_summary)) {
    glucose_summary <- rebuilt$glucose_summary
  }
}

if (!all(feature_panel_required_cols %in% names(feature_panel))) {
  stop(
    "feature_panel is missing required rate columns even after rebuild: ",
    paste(setdiff(feature_panel_required_cols, names(feature_panel)), collapse = ", ")
  )
}

signature_panel <- summarize_glucose_signatures(feature_panel)
signature_effects <- compute_empirical_effects(signature_panel)$paired_effects_long
transfer <- evaluate_transfer_predictions(signature_panel)
transfer_predictions <- transfer$predictions
transfer_feature_summary <- transfer$feature_summary
transfer_case_summary <- transfer$case_summary

exp_fit_qc <- feature_panel %>%
  filter(G0 %in% c(5, 25)) %>%
  group_by(cellLine, ploidy_metric) %>%
  summarise(
    n_highG_wells = sum(is.finite(exp_growth_rate_70)),
    median_exp_growth_rate_70 = safe_median(exp_growth_rate_70),
    median_fit_end_time_70 = safe_median(exp_growth_fit_end_time_70),
    median_fit_n_points_70 = safe_median(exp_growth_fit_n_points_70),
    median_fit_r2_70 = safe_median(exp_growth_fit_r2_70),
    median_initial_live_highG = safe_median(live_initial),
    qc_flag = dplyr::case_when(
      n_highG_wells < 2 ~ "too_few_highG_wells",
      !is.finite(median_fit_r2_70) ~ "missing_fit_r2",
      median_fit_r2_70 < 0.75 ~ "low_fit_r2",
      median_fit_n_points_70 < 5 ~ "short_fit_window",
      TRUE ~ "ok"
    ),
    .groups = "drop"
  )

transfer_feature_ranked <- transfer_feature_summary %>%
  left_join(feature_catalog, by = "feature") %>%
  select(feature, short_label, category, n_cases, scaled_mae_null, scaled_mae_transfer, mean_scaled_err_improvement, transfer_win_rate, sign_accuracy) %>%
  arrange(desc(mean_scaled_err_improvement), desc(transfer_win_rate), short_label)

signature_effect_ranked <- signature_effects %>%
  left_join(feature_catalog, by = "feature") %>%
  arrange(desc(abs(effect_per_ploidy)), category, short_label, cellLine)

print_feature_details <- function(features) {
  rows <- feature_catalog %>% filter(feature %in% features)
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    cat("**", row$short_label, "** (`", row$feature, "`)\n\n", sep = "")
    cat("Definition: ", row$definition, "\n\n", sep = "")
    cat("Rationale: ", row$rationale, "\n\n", sep = "")
    cat("Computation: ", row$computation, "\n\n", sep = "")
    cat("How to read it: ", row$interpretation, "\n\n", sep = "")
  }
}

write_model_free_agent_text_report <- function(
  path = paths$agent_text,
  feature_catalog,
  feature_panel,
  signature_panel,
  signature_effects,
  transfer_feature_summary,
  transfer_case_summary,
  exp_fit_qc
) {
  top_transfer <- transfer_feature_summary %>%
    left_join(feature_catalog %>% select(feature, short_label, category), by = "feature") %>%
    select(feature, short_label, category, n_cases, mean_scaled_err_improvement, transfer_win_rate, sign_accuracy, scaled_mae_null, scaled_mae_transfer) %>%
    arrange(desc(mean_scaled_err_improvement), desc(transfer_win_rate), short_label)

  top_effects <- signature_effects %>%
    left_join(feature_catalog %>% select(feature, short_label, category), by = "feature") %>%
    mutate(direction = if_else(effect_per_ploidy >= 0, "positive", "negative")) %>%
    arrange(desc(abs(effect_per_ploidy)), short_label, cellLine)

  qc_flags <- exp_fit_qc %>%
    arrange(qc_flag != "ok", median_fit_r2_70)

  signature_long <- signature_panel %>%
    select(cellLine, ploidy_metric, all_of(feature_catalog$feature)) %>%
    pivot_longer(cols = all_of(feature_catalog$feature), names_to = "feature", values_to = "value") %>%
    left_join(feature_catalog %>% select(feature, short_label, category), by = "feature") %>%
    arrange(category, short_label, cellLine, ploidy_metric)

  feature_definitions <- feature_catalog %>%
    select(feature, short_label, category, definition, interpretation)

  lines <- c(
    "MODEL_FREE_PLOIDY_ANALYSIS",
    sprintf("generated\t%s", Sys.time()),
    sprintf("stan_data\t%s", paths$stan_data),
    sprintf("feature_panel_source\t%s", if (file.exists(paths$feature_panel)) paths$feature_panel else "rebuilt_from_stan_data"),
    "",
    "SECTION\tKEY_RESULTS",
    capture.output(write.table(as.data.frame(top_transfer), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tHIGH_G_EXP_GROWTH_QC",
    capture.output(write.table(as.data.frame(qc_flags), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tPAIRED_PLOIDY_EFFECTS",
    capture.output(write.table(as.data.frame(top_effects), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tTRANSFER_CASE_SUMMARY",
    capture.output(write.table(as.data.frame(transfer_case_summary %>% arrange(desc(mean_scaled_err_improvement), desc(transfer_win_rate))), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tFEATURE_DEFINITIONS",
    capture.output(write.table(as.data.frame(feature_definitions), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tPER_LINE_SIGNATURE_VALUES",
    capture.output(write.table(as.data.frame(signature_long), sep = "\t", row.names = FALSE, quote = FALSE))
  )

  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

plot_feature_scatter <- function(raw_col, title, ylab) {
  ggplot(feature_panel, aes(as_g0_factor(G0), .data[[raw_col]], color = factor(ploidy_metric))) +
    geom_point(size = 2) +
    geom_line(aes(group = interaction(cellLine, ploidy_metric)), alpha = 0.5) +
    facet_wrap(~cellLine, scales = "free_y") +
    scale_color_brewer(palette = "Set1", name = "ploidy") +
    labs(title = title, x = "initial glucose (G0)", y = ylab) +
    theme_bw()
}

build_fit_df <- function(raw_col, max_g0 = Inf) {
  feature_panel %>%
    filter(is.finite(.data[[raw_col]]), G0 <= max_g0) %>%
    group_by(cellLine, ploidy_metric) %>%
    group_modify(~ {
      if (nrow(.x) < 2) {
        return(tibble())
      }
      fit <- lm(reformulate("log1p(G0)", response = raw_col), data = .x)
      tibble(
        G0 = sort(unique(.x$G0)),
        fitted = predict(fit, newdata = data.frame(G0 = sort(unique(.x$G0))))
      )
    }) %>%
    ungroup()
}

plot_feature_regression <- function(raw_col, title, ylab, max_g0 = Inf) {
  fit_df <- build_fit_df(raw_col, max_g0 = max_g0)

  ggplot(feature_panel %>% filter(is.finite(.data[[raw_col]]), G0 <= max_g0), aes(as_g0_factor(G0), .data[[raw_col]], color = factor(ploidy_metric))) +
    geom_point(size = 2) +
    geom_line(data = fit_df, aes(x = as_g0_factor(G0), y = fitted, color = factor(ploidy_metric), group = interaction(cellLine, ploidy_metric)), linewidth = 0.8) +
    facet_wrap(~cellLine, scales = "free_y") +
    scale_color_brewer(palette = "Set1", name = "ploidy") +
    labs(title = title, x = "initial glucose (G0)", y = ylab) +
    theme_bw()
}

plot_signature_values <- function(sig_col, title, ylab) {
  ggplot(signature_panel, aes(factor(ploidy_metric), .data[[sig_col]], fill = factor(ploidy_metric))) +
    geom_col() +
    facet_wrap(~cellLine, scales = "free_y") +
    scale_fill_brewer(palette = "Set1", name = "ploidy") +
    labs(title = title, x = "ploidy metric", y = ylab) +
    theme_bw()
}

build_live_reference_df <- function() {
  feature_panel %>%
    filter(G0 %in% c(5, 25)) %>%
    group_by(cellLine, ploidy_metric) %>%
    summarise(
      initial_live_median = median(live_initial, na.rm = TRUE),
      max_live_obs = max(live_peak, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      signature_panel %>%
        select(cellLine, ploidy_metric, growth_highG_median),
      by = c("cellLine", "ploidy_metric")
    ) %>%
    inner_join(
      count_summary %>%
        group_by(cellLine, ploidy_metric) %>%
        summarise(hours = list(sort(unique(hours))), .groups = "drop") %>%
        tidyr::unnest(hours),
      by = c("cellLine", "ploidy_metric")
    ) %>%
    group_by(cellLine, ploidy_metric) %>%
    mutate(
      start_hour = min(hours, na.rm = TRUE),
      ref_live_cells = pmax(
        expm1(log1p(initial_live_median) + growth_highG_median * (hours - start_hour)),
        0
      ),
      ref_live_cells = if_else(ref_live_cells <= max_live_obs, ref_live_cells, NA_real_)
    ) %>%
    ungroup()
}

plot_trajectory_view <- function(value_col, title, ylab, add_live_reference = FALSE, swap_facets = FALSE) {
  count_long <- feature_panel %>%
    select(cellLine, ploidy_metric, G0, well_idx)

  p <- count_summary %>%
    inner_join(count_long, by = c("well_idx", "cellLine", "ploidy_metric", "G0")) %>%
    ggplot(aes(hours, .data[[value_col]], color = as_g0_factor(G0), group = interaction(well_idx, G0))) +
    geom_line(alpha = 0.7) +
    scale_color_brewer(palette = "Spectral", name = "G0") +
    labs(title = title, x = "hours", y = ylab) +
    theme_bw()

  if (swap_facets) {
    p <- p + facet_grid(rows = vars(cellLine), cols = vars(factor(ploidy_metric)), scales = "free_y")
  } else {
    p <- p + facet_grid(rows = vars(factor(ploidy_metric)), cols = vars(cellLine), scales = "free_y")
  }

  if (add_live_reference && identical(value_col, "live_cells")) {
    p <- p +
      geom_line(
        data = build_live_reference_df(),
        aes(hours, ref_live_cells, group = interaction(cellLine, ploidy_metric)),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.9,
        linetype = "22"
      )
  }

  p
}

plot_rate_processing <- function(raw_col, smooth_col, title, ylab) {
  feature_panel %>%
    select(cellLine, ploidy_metric, G0, all_of(raw_col), all_of(smooth_col)) %>%
    pivot_longer(cols = c(all_of(raw_col), all_of(smooth_col)), names_to = "rate_type", values_to = "value") %>%
    mutate(rate_type = recode(rate_type, !!raw_col := "raw", !!smooth_col := "smoothed")) %>%
    ggplot(aes(as_g0_factor(G0), value, color = factor(ploidy_metric), shape = rate_type, group = interaction(cellLine, ploidy_metric, rate_type))) +
    geom_point(size = 2) +
    geom_line(alpha = 0.6) +
    facet_wrap(~cellLine, scales = "free_y") +
    scale_color_brewer(palette = "Set1", name = "ploidy") +
    labs(title = title, x = "initial glucose (G0)", y = ylab, shape = "rate source") +
    theme_bw()
}

plot_yield_inputs <- function() {
  feature_panel %>%
    select(cellLine, ploidy_metric, G0, glucose_drawdown, total_net_gain_to_glucose_end, peak_total_yield_per_glucose) %>%
    pivot_longer(
      cols = c(glucose_drawdown, total_net_gain_to_glucose_end, peak_total_yield_per_glucose),
      names_to = "component",
      values_to = "value"
    ) %>%
    ggplot(aes(as_g0_factor(G0), value, color = factor(ploidy_metric), group = interaction(cellLine, ploidy_metric))) +
    geom_point(size = 2) +
    geom_line(alpha = 0.6) +
    facet_grid(component ~ cellLine, scales = "free_y") +
    scale_color_brewer(palette = "Set1", name = "ploidy") +
    labs(title = "Peak Yield Inputs By Condition", x = "initial glucose (G0)", y = NULL) +
    theme_bw()
}


## ----executive-transfer-------------------------------------------------------
transfer_feature_ranked


## ----executive-qc-------------------------------------------------------------
exp_fit_qc %>%
  arrange(qc_flag != "ok", median_fit_r2_70)


## ----transfer-bar, fig.height=5, fig.width=9----------------------------------
transfer_feature_ranked %>%
  mutate(short_label = forcats::fct_reorder(short_label, mean_scaled_err_improvement)) %>%
  ggplot(aes(x = short_label, y = mean_scaled_err_improvement, fill = transfer_win_rate)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Transfer Gain By Curated Feature",
    x = NULL,
    y = "mean scaled error improvement",
    fill = "transfer\nwin rate"
  ) +
  theme_bw()


## ----transfer-summary---------------------------------------------------------
transfer_feature_ranked


## ----transfer-case-summary----------------------------------------------------
transfer_case_summary %>%
  arrange(desc(mean_scaled_err_improvement))


## ----effect-heatmap, fig.height=5, fig.width=10-------------------------------
signature_effects %>%
  left_join(feature_catalog, by = "feature") %>%
  group_by(feature) %>%
  mutate(effect_z = as.numeric(scale(effect_per_ploidy))) %>%
  ungroup() %>%
  ggplot(aes(x = short_label, y = cellLine, fill = effect_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(
    title = "Standardized Ploidy Effects Across Curated Features",
    x = NULL,
    y = NULL,
    fill = "within-feature\nz-score"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----feature-catalog----------------------------------------------------------
feature_catalog %>%
  select(feature, short_label, category, definition, rationale, computation, interpretation)


## ----growth-lowg-text, results='asis'-----------------------------------------
print_feature_details(c("growth_lowG_median", "growth_highG_median"))


## ----growth-plot, fig.height=6, fig.width=10----------------------------------
plot_trajectory_view("live_cells", "Underlying Live-Cell Trajectories With High-G Growth Reference", "live cells", add_live_reference = TRUE, swap_facets = TRUE)


## ----growth-processing, fig.height=6, fig.width=10----------------------------
plot_rate_processing(
  raw_col = "max_growth_rate_raw",
  smooth_col = "max_growth_rate",
  title = "Raw vs Smoothed Max Growth Rate By Condition",
  ylab = "max growth rate"
)


## ----growth-signatures, fig.height=7, fig.width=10----------------------------
plot_signature_values("growth_lowG_median", "Aggregated Feature: Growth Low-G Median", "median smoothed max growth rate at G0 <= 0.25") /
  plot_signature_values("growth_highG_median", "Aggregated Feature: Growth High-G Median", "median exponential-phase growth rate from G0 = 5, 25")


## ----death-lowg-text, results='asis'------------------------------------------
print_feature_details(c("death_lowG_median", "death_highG_median"))


## ----death-plot, fig.height=6, fig.width=10-----------------------------------
plot_trajectory_view("dead_cells", "Underlying Dead-Cell Trajectories", "dead cells")


## ----death-processing, fig.height=6, fig.width=10-----------------------------
plot_rate_processing(
  raw_col = "max_death_rate_raw",
  smooth_col = "max_death_rate",
  title = "Raw vs Smoothed Max Death Rate By Condition",
  ylab = "max death rate"
)


## ----death-signatures, fig.height=7, fig.width=10-----------------------------
plot_signature_values("death_lowG_median", "Aggregated Feature: Death Low-G Median", "median smoothed max death rate at G0 <= 0.25") /
  plot_signature_values("death_highG_median", "Aggregated Feature: Death High-G Median", "median smoothed max death rate at G0 >= 1")


## ----aliveauc-int-text, results='asis'----------------------------------------
print_feature_details(c("yield_alive_auc_intercept", "yield_alive_auc_slope"))


## ----aliveauc-reg-plot, fig.height=6, fig.width=10----------------------------
plot_feature_regression(
  raw_col = "live_auc_glucose_window",
  title = "Alive AUC Over The Glucose-Measurement Window vs Initial Glucose",
  ylab = "alive AUC over glucose window",
  max_g0 = 1
)


## ----aliveauc-raw-view, fig.height=6, fig.width=10----------------------------
plot_feature_scatter(
  raw_col = "live_auc_glucose_window",
  title = "Alive AUC Over The Glucose-Measurement Window By Condition",
  ylab = "alive AUC over glucose window"
)


## ----aliveauc-signatures, fig.height=7, fig.width=10--------------------------
plot_signature_values("yield_alive_auc_intercept", "Aggregated Feature: Alive AUC Regression Intercept", "intercept of alive AUC ~ log1p(G0), fit on G0 <= 1") /
  plot_signature_values("yield_alive_auc_slope", "Aggregated Feature: Alive AUC Regression Slope", "slope of alive AUC ~ log1p(G0), fit on G0 <= 1")


## ----peakyield-int-text, results='asis'---------------------------------------
print_feature_details(c("peak_total_yield_intercept", "peak_total_yield_slope"))


## ----peakyield-reg-plot, fig.height=6, fig.width=10---------------------------
plot_feature_regression(
  raw_col = "peak_total_yield_per_glucose",
  title = "Peak Total Yield Per Glucose Regressed Against Initial Glucose",
  ylab = "net gain in total cells / glucose drawdown"
)


## ----peakyield-inputs, fig.height=8, fig.width=10-----------------------------
plot_yield_inputs()


## ----peakyield-signatures, fig.height=7, fig.width=10-------------------------
plot_signature_values("peak_total_yield_intercept", "Aggregated Feature: Peak Yield Regression Intercept", "intercept of peak yield ~ log1p(G0)") /
  plot_signature_values("peak_total_yield_slope", "Aggregated Feature: Peak Yield Regression Slope", "slope of peak yield ~ log1p(G0)")


## ----per-condition-preview----------------------------------------------------
feature_panel %>%
  select(
    cellLine, ploidy_metric, G0,
    max_growth_rate, max_death_rate,
    exp_growth_rate_70, exp_growth_fit_r2_70,
    live_auc_glucose_window,
    total_net_gain_to_glucose_end,
    glucose_initial, glucose_final,
    peak_total_yield_per_glucose
  ) %>%
  arrange(cellLine, ploidy_metric, G0)


## ----output-check-------------------------------------------------------------
tibble(
  artifact = c(names(paths), "agent_text"),
  path = unlist(c(paths, agent_text = paths$agent_text)),
  exists = file.exists(unlist(c(paths, agent_text = paths$agent_text)))
)


## ----agent-text-export, include=FALSE-----------------------------------------
write_model_free_agent_text_report(
  path = paths$agent_text,
  feature_catalog = feature_catalog,
  feature_panel = feature_panel,
  signature_panel = signature_panel,
  signature_effects = signature_effects,
  transfer_feature_summary = transfer_feature_summary,
  transfer_case_summary = transfer_case_summary,
  exp_fit_qc = exp_fit_qc
)

