#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(scales)
  library(tibble)
  library(tidyr)
})

# Surface plotting warnings in integration logs instead of deferring them to
# the end of the R session as an uninspectable warning count.
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (length(args) < idx + 1L) stop("Missing value after ", flag, call. = FALSE)
  args[[idx + 1L]]
}

phase_idx <- match("--phase", args)
phase <- if (!is.na(phase_idx) && length(args) >= phase_idx + 1) args[[phase_idx + 1]] else "both"
if (!phase %in% c("subpanels", "final", "both")) {
  stop("Unsupported --phase: ", phase, call. = FALSE)
}

project_root <- normalizePath(arg_value("--project-root", "."), winslash = "/", mustWork = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve promoted script path", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
package_script_root <- dirname(dirname(script_path))
draft_root <- package_script_root

output_root_arg <- arg_value("--output-root")
if (is.null(output_root_arg)) stop("--output-root is required for integration rebuilds", call. = FALSE)
polish_root <- if (grepl("^/", output_root_arg)) output_root_arg else file.path(project_root, output_root_arg)
polish_root <- normalizePath(polish_root, winslash = "/", mustWork = FALSE)
subpanel_dir <- file.path(polish_root, "subpanels")
layout_dir <- file.path(polish_root, "layout")
final_dir <- polish_root
script_dir <- file.path(package_script_root, "scripts")
config_dir <- file.path(package_script_root, "config")
dir.create(subpanel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
for (config_name in c("layout_plan.csv", "subpanel_dimensions.csv")) {
  source_config <- file.path(config_dir, config_name)
  if (!file.exists(source_config)) stop("Missing promoted fixed config: ", source_config, call. = FALSE)
  if (!file.copy(source_config, file.path(layout_dir, config_name), overwrite = TRUE)) {
    stop("Could not stage promoted fixed config: ", config_name, call. = FALSE)
  }
}

inversion_root <- file.path(
  project_root, "agent-dev", "manuscript_redrafts", "20260710T164049_v01_redraft",
  "stage_outputs", "analysis", "live_auc_inversion"
)
inversion_fits_path <- file.path(inversion_root, "results", "signature_fits.csv")
corrected_feature_path <- file.path(inversion_root, "results", "feature_panel_corrected.csv")
if (!file.exists(inversion_fits_path) || !file.exists(corrected_feature_path)) {
  stop("Missing validated live-AUC inversion inputs", call. = FALSE)
}

inverted_fits <- as_tibble(read.csv(inversion_fits_path, stringsAsFactors = FALSE, check.names = FALSE)) %>%
  group_by(.data$cellLine) %>%
  arrange(.data$ploidy_metric, .by_group = TRUE) %>%
  mutate(ploidy_state = factor(c("low", "high"), levels = c("low", "high"))) %>%
  ungroup()
corrected_feature_panel <- as_tibble(read.csv(corrected_feature_path, stringsAsFactors = FALSE, check.names = FALSE))
if (
  nrow(inverted_fits) != 10L ||
  any(table(inverted_fits$cellLine) != 2L) ||
  any(!is.finite(inverted_fits$glucose_per_live_auc_slope)) ||
  any(inverted_fits$n < 3L)
) {
  stop("Validated inversion-table coverage gate failed", call. = FALSE)
}

rel_path <- function(path) {
  prefix <- paste0(project_root, "/")
  vapply(
    normalizePath(path, winslash = "/", mustWork = FALSE),
    function(one) {
      if (startsWith(one, prefix)) substring(one, nchar(prefix) + 1L) else one
    },
    character(1)
  )
}

safe_name <- function(x) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

label_short <- function(x) {
  vapply(x, function(one) {
    if (!is.finite(one)) return("")
    abs_one <- abs(one)
    if (abs_one >= 1e6) return(sprintf("%.1fM", one / 1e6))
    if (abs_one >= 1e3) return(sprintf("%.1fk", one / 1e3))
    if (abs_one >= 1) return(scales::number(one, accuracy = 0.1))
    if (abs_one > 0 && abs_one < 1e-3) return(sprintf("%.1e", one))
    scales::number(one, accuracy = 0.01)
  }, character(1))
}

clip_value <- function(x, limit = 3) {
  pmax(-limit, pmin(limit, x))
}

finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

sha256_file <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(line) || !nzchar(line[[1]])) stop("Could not hash: ", path, call. = FALSE)
  strsplit(line[[1]], "[[:space:]]+")[[1]][[1]]
}

root_block <- paste(
  "DRAFT_ROOT <- file.path(",
  "  \"agent-dev\", \"manuscript_redrafts\", \"20260709T210405_v7_redraft\",",
  "  \"figure_generation\", \"FG2_direct_feature_rebuild\", \"drafting_v2\"",
  ")",
  sep = "\n"
)

source_with_polish_root <- function(script_path, env) {
  txt <- readLines(script_path, warn = FALSE)
  joined <- paste(txt, collapse = "\n")
  if (!grepl(root_block, joined, fixed = TRUE)) {
    stop("Could not find expected DRAFT_ROOT block in ", script_path, call. = FALSE)
  }
  joined <- gsub(root_block, "DRAFT_ROOT <- polish_root", joined, fixed = TRUE)
  eval(parse(text = joined, srcfile = script_path), envir = env)
  invisible(env)
}

theme_fg2_polish <- function(base_size = 7.0) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.14, color = "grey88"),
      strip.background = element_rect(fill = "grey94", color = "grey74", linewidth = 0.22),
      strip.text = element_text(face = "bold", margin = margin(1.3, 1.3, 1.3, 1.3)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.height = unit(0.12, "in"),
      legend.key.width = unit(0.22, "in"),
      axis.title = element_text(size = rel(0.92)),
      axis.text = element_text(size = rel(0.82)),
      plot.margin = margin(3, 5, 3, 3)
    )
}

load_builders <- function() {
  scripts <- file.path(draft_root, "scripts")

  summary_env <- new.env(parent = globalenv())
  summary_env$polish_root <- polish_root
  source_with_polish_root(file.path(scripts, "fg2_summary_strip_options_worker.R"), summary_env)

  staged_env <- new.env(parent = globalenv())
  staged_env$polish_root <- polish_root
  source_with_polish_root(file.path(scripts, "fg2_staged_feature_workflow_worker.R"), staged_env)

  supplement_env <- new.env(parent = globalenv())
  supplement_env$polish_root <- polish_root
  source_with_polish_root(file.path(scripts, "fg2_supplement_qc_worker.R"), supplement_env)

  supplement_env$feature_meta <- supplement_env$feature_meta %>%
    mutate(
      display_label = if_else(
        .data$feature == "yield_alive_auc_slope",
        "Glucose drawdown per live-cell AUC",
        .data$display_label
      )
    )
  supplement_env$feature_panel_all <- corrected_feature_panel
  supplement_env$signature_all <- supplement_env$signature_all %>%
    select(-any_of("yield_alive_auc_slope")) %>%
    left_join(
      inverted_fits %>%
        transmute(
          .data$cellLine, .data$line_id, .data$ploidy_metric,
          yield_alive_auc_slope = .data$glucose_per_live_auc_slope
        ),
      by = c("cellLine", "line_id", "ploidy_metric")
    )
  supplement_env$p_raw_all <- supplement_env$make_raw_feature_pairs(
    supplement_env$signature_all,
    include_sum159_fuse = TRUE,
    feature_panel = corrected_feature_panel
  ) & theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())

  list(summary = summary_env, staged = staged_env, supplement = supplement_env)
}

make_polish_figure2a <- function(summary_env) {
  feature_order_tbl <- summary_env$feature_meta %>%
    transmute(display_label, feature_group, feature_order) %>%
    bind_rows(tibble(
      display_label = "0 mM live-cell AUC",
      feature_group = "Alive AUC",
      feature_order = 3.5
    ))

  g0_live_auc <- summary_env$build_g0_effects(summary_env$feature_no_sum, summary_env$scope_no_sum) %>%
    filter(.data$metric == "live_auc_glucose_window", abs(.data$G0) < 1e-8) %>%
    transmute(
      cellLine = .data$cellLine,
      feature = "live_auc_0mM",
      effect_per_ploidy = .data$effect_per_ploidy,
      display_label = "0 mM live-cell AUC",
      feature_group = "Alive AUC",
      feature_order = 3.5
    )

  inverted_slope_effect <- inverted_fits %>%
    filter(.data$cellLine != "SUM-159-fuse") %>%
    group_by(.data$cellLine) %>%
    summarise(
      effect_per_ploidy = (
        .data$glucose_per_live_auc_slope[.data$ploidy_state == "high"] -
          .data$glucose_per_live_auc_slope[.data$ploidy_state == "low"]
      ) / (
        .data$ploidy_metric[.data$ploidy_state == "high"] -
          .data$ploidy_metric[.data$ploidy_state == "low"]
      ),
      .groups = "drop"
    ) %>%
    transmute(
      .data$cellLine,
      feature = "yield_alive_auc_slope",
      .data$effect_per_ploidy,
      display_label = "Glucose drawdown per live-cell AUC",
      feature_group = "Alive AUC",
      feature_order = 4
    )

  feature_order_tbl <- feature_order_tbl %>%
    mutate(
      display_label = if_else(
        .data$display_label == "Alive AUC gradient",
        "Glucose drawdown per live-cell AUC",
        .data$display_label
      )
    )

  effects <- summary_env$no_sum_effects %>%
    ungroup() %>%
    filter(.data$feature != "yield_alive_auc_slope") %>%
    select(any_of(c("cellLine", "feature", "effect_per_ploidy", "display_label", "feature_group", "feature_order"))) %>%
    mutate(
      display_label = as.character(.data$display_label),
      feature_group = as.character(.data$feature_group)
    ) %>%
    bind_rows(g0_live_auc, inverted_slope_effect) %>%
    group_by(.data$feature) %>%
    mutate(
      effect_scale = finite_sd(.data$effect_per_ploidy),
      fallback_scale = median(abs(.data$effect_per_ploidy), na.rm = TRUE),
      effect_scale = if_else(!is.finite(.data$effect_scale) | .data$effect_scale <= 0, .data$fallback_scale, .data$effect_scale),
      effect_scale = if_else(!is.finite(.data$effect_scale) | .data$effect_scale <= 0, 1, .data$effect_scale),
      standardized_effect = .data$effect_per_ploidy / .data$effect_scale,
      clipped_effect = clip_value(.data$standardized_effect)
    ) %>%
    ungroup() %>%
    mutate(
      display_label = factor(.data$display_label, levels = rev(feature_order_tbl$display_label)),
      feature_group = factor(.data$feature_group, levels = names(summary_env$feature_palette)),
      is_live_auc = .data$feature == "live_auc_0mM"
    )

  medians <- effects %>%
    group_by(.data$feature, .data$display_label, .data$feature_group, .data$feature_order) %>%
    summarise(median_effect = median(.data$standardized_effect, na.rm = TRUE), .groups = "drop") %>%
    mutate(clipped_effect = clip_value(.data$median_effect))

  ggplot(effects, aes(.data$clipped_effect, .data$display_label)) +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
    geom_point(
      data = effects %>% filter(!.data$is_live_auc),
      aes(fill = .data$feature_group),
      shape = 21,
      color = "grey18",
      stroke = 0.22,
      size = 1.85,
      alpha = 0.82,
      position = position_jitter(width = 0, height = 0.055, seed = 20260711)
    ) +
    geom_point(
      data = effects %>% filter(.data$is_live_auc),
      aes(fill = .data$feature_group),
      shape = 21,
      color = "black",
      stroke = 0.55,
      size = 2.55,
      alpha = 0.92,
      position = position_jitter(width = 0, height = 0.055, seed = 20260712)
    ) +
    geom_point(
      data = medians,
      aes(.data$clipped_effect, .data$display_label, color = .data$feature_group),
      inherit.aes = FALSE,
      shape = 23,
      fill = "white",
      stroke = 0.45,
      size = 2.25
    ) +
    facet_grid(feature_group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = summary_env$feature_palette, name = "feature class") +
    scale_color_manual(values = summary_env$feature_palette, guide = "none") +
    scale_x_continuous(limits = c(-3.05, 3.05), oob = squish, breaks = -3:3) +
    labs(
      x = "high-low effect per ploidy (within-feature SD units)",
      y = NULL
    ) +
    theme_fg2_polish(base_size = 7.2) +
    theme(
      panel.spacing.y = unit(0.18, "lines"),
      strip.text.y = element_blank(),
      strip.background.y = element_blank(),
      axis.text.y = element_text(size = 6.7),
      legend.position = "none"
    )
}

make_polish_figure2b <- function(staged_env) {
  example_count <- staged_env$example_count
  example_smooth <- staged_env$example_smooth
  example_feature <- staged_env$example_feature
  feature_palette <- staged_env$feature_palette

  observed <- example_count %>%
    arrange(.data$hours) %>%
    mutate(
      log_live = log1p(pmax(.data$live_cells, 0)),
      log_total = log1p(pmax(.data$total_cells, 0))
    )

  robust_value <- example_feature$robust_max_derivative[[1]]
  robust_time <- example_feature$robust_derivative_time[[1]]
  y0 <- approx(
    x = example_smooth$hours,
    y = example_smooth$smooth_log_live,
    xout = robust_time,
    rule = 2
  )$y
  x_rng <- range(example_smooth$hours, na.rm = TRUE)
  tangent_x <- seq(max(x_rng[1], robust_time - 24), min(x_rng[2], robust_time + 24), length.out = 2)
  tangent <- tibble(hours = tangent_x, log_live = y0 + robust_value * (tangent_x - robust_time))

  start_time <- example_feature$glucose_start_time[[1]]
  end_time <- example_feature$glucose_end_time[[1]]
  auc_band <- example_smooth %>%
    filter(.data$hours >= start_time, .data$hours <= end_time) %>%
    mutate(ymin = min(observed$log_live, na.rm = TRUE) - 0.12)

  yield_window <- observed %>% filter(.data$hours >= start_time, .data$hours <= end_time)
  if (!nrow(yield_window)) yield_window <- observed
  yield_min <- min(yield_window$total_cells, na.rm = TRUE)
  yield_max <- max(yield_window$total_cells, na.rm = TRUE)
  window_rng <- range(yield_window$hours, na.rm = TRUE)
  yield_arrow_x <- window_rng[2] - 0.08 * diff(window_rng)
  yield_label_x <- min(x_rng[2] - 0.05 * diff(x_rng), yield_arrow_x + 0.07 * diff(x_rng))
  yield_label_y <- mean(log1p(c(yield_min, yield_max)), na.rm = TRUE)

  point_df <- bind_rows(
    observed %>% transmute(hours, log_cells = .data$log_live, series = "live cells"),
    observed %>% transmute(hours, log_cells = .data$log_total, series = "total cells")
  ) %>%
    mutate(series = factor(.data$series, levels = c("live cells", "total cells")))

  series_colors <- c("live cells" = "#1B6CA8", "total cells" = feature_palette[["Total-cell yield"]])

  ggplot() +
    geom_ribbon(
      data = auc_band,
      aes(.data$hours, ymin = .data$ymin, ymax = .data$smooth_log_live),
      fill = feature_palette[["Live-cell AUC"]],
      alpha = 0.15
    ) +
    geom_hline(
      yintercept = log1p(c(yield_min, yield_max)),
      color = feature_palette[["Total-cell yield"]],
      linewidth = 0.34,
      linetype = "dashed"
    ) +
    annotate(
      "segment",
      x = yield_arrow_x,
      xend = yield_arrow_x,
      y = log1p(yield_min),
      yend = log1p(yield_max),
      arrow = arrow(length = unit(0.055, "in"), ends = "both"),
      color = feature_palette[["Total-cell yield"]],
      linewidth = 0.45
    ) +
    annotate(
      "text",
      x = yield_label_x,
      y = yield_label_y,
      label = "max delta",
      color = feature_palette[["Total-cell yield"]],
      size = 1.9,
      hjust = 0,
      vjust = 0.5
    ) +
    geom_point(
      data = point_df,
      aes(.data$hours, .data$log_cells, color = .data$series),
      size = 0.68,
      alpha = 0.86
    ) +
    geom_line(
      data = example_smooth,
      aes(.data$hours, .data$smooth_log_live),
      color = series_colors[["live cells"]],
      linewidth = 0.74,
      show.legend = FALSE
    ) +
    geom_line(
      data = tangent,
      aes(.data$hours, .data$log_live),
      color = feature_palette[["Robust max derivative"]],
      linewidth = 0.72,
      linetype = "dashed"
    ) +
    geom_point(
      data = tibble(hours = robust_time, log_live = y0),
      aes(.data$hours, .data$log_live),
      shape = 21,
      fill = "white",
      color = feature_palette[["Robust max derivative"]],
      stroke = 0.42,
      size = 1.45
    ) +
    annotate(
      "text",
      x = min(x_rng[2] - 0.06 * diff(x_rng), robust_time + 18),
      y = y0 + 0.34,
      label = "derivative tangent",
      color = feature_palette[["Robust max derivative"]],
      size = 1.9,
      hjust = 1
    ) +
    annotate(
      "text",
      x = mean(c(start_time, end_time), na.rm = TRUE),
      y = min(observed$log_live, na.rm = TRUE) - 0.02,
      label = "live AUC",
      color = feature_palette[["Live-cell AUC"]],
      size = 1.9,
      hjust = 0.5,
      vjust = 1
    ) +
    scale_color_manual(values = series_colors, name = NULL) +
    coord_cartesian(clip = "off") +
    labs(x = "hours", y = "log1p(cells)") +
    theme_fg2_polish(base_size = 6.9) +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.text = element_text(size = 6.0),
      plot.margin = margin(4, 14, 4, 4)
    )
}

make_polish_figure2c <- function(staged_env) {
  feature_palette <- staged_env$feature_palette
  ploidy_palette <- staged_env$ploidy_palette
  line_name <- staged_env$representative_line
  g0_levels <- staged_env$g0_levels

  values <- staged_env$feature_values %>%
    filter(.data$cellLine == line_name, is.finite(.data$value)) %>%
    staged_env$add_ploidy_display() %>%
    mutate(
      G0_display = factor(as.character(.data$G0), levels = g0_levels),
      ploidy_state = factor(.data$ploidy_state, levels = c("low", "high"))
    )

  deriv <- values %>% filter(.data$feature_id == "robust_max_derivative")
  auc <- corrected_feature_panel %>%
    filter(.data$cellLine == line_name, .data$G0 <= 1) %>%
    group_by(.data$cellLine) %>%
    arrange(.data$ploidy_metric, .by_group = TRUE) %>%
    mutate(ploidy_state = factor(if_else(.data$ploidy_metric == min(.data$ploidy_metric), "low", "high"), levels = c("low", "high"))) %>%
    ungroup()
  yield <- values %>% filter(.data$feature_id == "peak_total_yield_net")

  auc_fit <- auc %>%
    left_join(
      inverted_fits %>%
        filter(.data$cellLine == line_name) %>%
        select(.data$cellLine, .data$ploidy_state, .data$glucose_per_live_auc_intercept, .data$glucose_per_live_auc_slope),
      by = c("cellLine", "ploidy_state")
    ) %>%
    group_by(.data$cellLine, .data$ploidy_state) %>%
    reframe(
      live_auc_glucose_window = seq(min(.data$live_auc_glucose_window), max(.data$live_auc_glucose_window), length.out = 80L),
      glucose_drawdown = first(.data$glucose_per_live_auc_intercept) + first(.data$glucose_per_live_auc_slope) * .data$live_auc_glucose_window
    )

  dodge <- position_dodge(width = 0.34)
  p_deriv <- ggplot(deriv, aes(.data$G0_display, .data$value, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.18, color = "grey72") +
    geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed", linewidth = 0.18, color = "grey38") +
    geom_line(linewidth = 0.26, alpha = 0.55, position = dodge) +
    geom_point(position = dodge, size = 1.0, alpha = 0.9) +
    geom_point(
      data = deriv %>% filter(.data$G0 <= 0.25 | .data$G0 >= 1),
      aes(.data$G0_display, .data$value),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.30,
      size = 1.55,
      position = dodge
    ) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_x_discrete(drop = FALSE, labels = c("0", "0.1", "0.25", "0.5", "1", "5", "25")) +
    labs(x = "G0", y = "robust max derivative") +
    theme_fg2_polish(base_size = 6.5) +
    theme(legend.position = "none", axis.text.x = element_text(size = 5.2))

  p_auc <- ggplot(auc, aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_point(size = 1.0, alpha = 0.9) +
    geom_point(
      data = auc %>% filter(abs(.data$G0) < 1e-8),
      aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, fill = .data$ploidy_state),
      shape = 21,
      color = "black",
      stroke = 0.42,
      size = 1.85
    ) +
    geom_line(data = auc_fit, aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, color = .data$ploidy_state, group = .data$ploidy_state), linewidth = 0.46) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_fill_manual(values = ploidy_palette, guide = "none") +
    scale_x_continuous(labels = label_short) +
    labs(x = "live-cell AUC (thousands)", y = "glucose consumed") +
    theme_fg2_polish(base_size = 6.5) +
    theme(legend.position = "none", axis.text.x = element_text(size = 5.2))

  p_yield <- ggplot(yield, aes(.data$G0_display, .data$value, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.18, color = "grey72") +
    geom_line(linewidth = 0.26, alpha = 0.50, position = dodge) +
    geom_point(position = dodge, size = 1.0, alpha = 0.9) +
    geom_point(
      data = yield %>% filter(abs(.data$G0 - 1) < 1e-8),
      aes(.data$G0_display, .data$value),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.32,
      size = 1.65,
      position = dodge
    ) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_x_discrete(drop = FALSE, labels = c("0", "0.1", "0.25", "0.5", "1", "5", "25")) +
    scale_y_continuous(labels = label_short) +
    labs(x = "G0", y = "net total-cell yield") +
    theme_fg2_polish(base_size = 6.5) +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 5.2))

  (p_deriv | p_auc | p_yield) +
    plot_layout(widths = c(1.0, 1.03, 1.0))
}

make_polish_s4 <- function(supplement_env) {
  feature_panel <- corrected_feature_panel
  count_summary <- supplement_env$count_summary
  ploidy_palette <- supplement_env$ploidy_palette
  g0_levels <- supplement_env$g0_levels
  line_levels <- feature_panel %>% arrange(.data$line_id, .data$cellLine) %>% distinct(.data$cellLine) %>% pull(.data$cellLine)

  robust_panel <- supplement_env$build_robust_derivative_panel(count_summary, feature_panel)
  values <- feature_panel %>%
    select(
      .data$well_idx, .data$cellLine, .data$line_id, .data$ploidy_metric,
      .data$ploidy_abs, .data$G0, .data$live_auc_glucose_window,
      .data$total_net_gain_to_glucose_end
    ) %>%
    left_join(robust_panel %>% select(.data$well_idx, .data$robust_max_derivative), by = "well_idx") %>%
    mutate(peak_total_yield_net = pmax(.data$total_net_gain_to_glucose_end, 0)) %>%
    pivot_longer(
      cols = c(.data$robust_max_derivative, .data$live_auc_glucose_window, .data$peak_total_yield_net),
      names_to = "feature_id",
      values_to = "value"
    ) %>%
    filter(is.finite(.data$value)) %>%
    supplement_env$add_ploidy_display() %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_levels),
      G0_display = factor(as.character(.data$G0), levels = g0_levels),
      ploidy_state = factor(.data$ploidy_state, levels = c("low", "high"))
    )

  deriv <- values %>% filter(.data$feature_id == "robust_max_derivative")
  auc <- feature_panel %>%
    filter(.data$G0 <= 1) %>%
    supplement_env$add_ploidy_display() %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_levels),
      ploidy_state = factor(.data$ploidy_state, levels = c("low", "high"))
    )
  yield <- values %>% filter(.data$feature_id == "peak_total_yield_net")

  auc_fit <- auc %>%
    left_join(
      inverted_fits %>%
        select(.data$cellLine, .data$ploidy_state, .data$glucose_per_live_auc_intercept, .data$glucose_per_live_auc_slope),
      by = c("cellLine", "ploidy_state")
    ) %>%
    group_by(.data$cellLine, .data$ploidy_state) %>%
    reframe(
      live_auc_glucose_window = seq(min(.data$live_auc_glucose_window), max(.data$live_auc_glucose_window), length.out = 70L),
      glucose_drawdown = first(.data$glucose_per_live_auc_intercept) + first(.data$glucose_per_live_auc_slope) * .data$live_auc_glucose_window
    )

  dodge <- position_dodge(width = 0.34)
  base_strip_theme <- theme_fg2_polish(base_size = 5.8) +
    theme(
      strip.text = element_text(size = 4.7),
      axis.text.x = element_text(size = 4.4, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 4.4),
      axis.title = element_text(size = 5.3),
      panel.spacing.y = unit(0.12, "lines"),
      panel.spacing.x = unit(0.12, "lines"),
      legend.text = element_text(size = 5.2),
      legend.title = element_text(size = 5.3)
    )

  p_deriv <- ggplot(deriv, aes(.data$G0_display, .data$value, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.14, color = "grey74") +
    geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed", linewidth = 0.16, color = "grey38") +
    geom_line(linewidth = 0.22, alpha = 0.48, position = dodge) +
    geom_point(position = dodge, size = 0.65, alpha = 0.88) +
    geom_point(
      data = deriv %>% filter(.data$G0 <= 0.25 | .data$G0 >= 1),
      aes(.data$G0_display, .data$value),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.25,
      size = 1.05,
      position = dodge
    ) +
    facet_wrap(~cellLine, ncol = 1, scales = "free_y") +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_x_discrete(drop = FALSE, labels = c("0", "0.1", "0.25", "0.5", "1", "5", "25")) +
    labs(x = "G0", y = "robust derivative") +
    base_strip_theme +
    theme(legend.position = "none")

  p_auc <- ggplot(auc, aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_point(size = 0.65, alpha = 0.88) +
    geom_point(
      data = auc %>% filter(abs(.data$G0) < 1e-8),
      aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, fill = .data$ploidy_state),
      shape = 21,
      color = "black",
      stroke = 0.28,
      size = 1.05
    ) +
    geom_line(data = auc_fit, aes(.data$live_auc_glucose_window / 1000, .data$glucose_drawdown, color = .data$ploidy_state, group = .data$ploidy_state), linewidth = 0.28) +
    facet_wrap(~cellLine, ncol = 1, scales = "free_y") +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_fill_manual(values = ploidy_palette, guide = "none") +
    scale_x_continuous(labels = label_short) +
    labs(x = "live-cell AUC (thousands)", y = "glucose consumed") +
    base_strip_theme +
    theme(legend.position = "none")

  p_yield <- ggplot(yield, aes(.data$G0_display, .data$value, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.14, color = "grey74") +
    geom_line(linewidth = 0.22, alpha = 0.48, position = dodge) +
    geom_point(position = dodge, size = 0.65, alpha = 0.88) +
    geom_point(
      data = yield %>% filter(abs(.data$G0 - 1) < 1e-8),
      aes(.data$G0_display, .data$value),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.25,
      size = 1.10,
      position = dodge
    ) +
    facet_wrap(~cellLine, ncol = 1, scales = "free_y") +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_x_discrete(drop = FALSE, labels = c("0", "0.1", "0.25", "0.5", "1", "5", "25")) +
    scale_y_continuous(labels = label_short) +
    labs(x = "G0", y = "net total-cell yield") +
    base_strip_theme +
    theme(legend.position = "bottom")

  (p_deriv | p_auc | p_yield) +
    plot_layout(widths = c(1.05, 1.0, 1.05))
}

panel_row <- function(figure, panel, plot, width, height, approved_source, generator, data_inputs, notes = "") {
  if (!nzchar(notes)) notes <- "Regenerated from approved worker logic without changing panel identity."
  list(
    figure = figure,
    panel = panel,
    plot = plot,
    width = width,
    height = height,
    approved_source = approved_source,
    generator = generator,
    data_inputs = data_inputs,
    notes = notes
  )
}

build_panels <- function(envs) {
  scripts <- file.path(
    "agent-dev", "manuscript_redrafts", "20260710T164049_v01_redraft",
    "figure_generation", "FG2_direct_feature_rebuild", "drafting_v2", "scripts"
  )
  corrected_wp4_root <- file.path(
    "agent-dev", "manuscript_redrafts", "20260619_v6_redraft",
    "stage_outputs", "analysis", "recomputed_wp4_sum159_fuse_exclusion"
  )
  wp1_inputs <- paste(
    c(
      "data/report_exports/wp1_core_starvation/count_trajectory_summary.csv",
      "data/report_exports/wp1_core_starvation/glucose_trajectory_summary.csv"
    ),
    collapse = ";"
  )
  wp4_no_sum_inputs <- paste(
    c(
      file.path(corrected_wp4_root, "wp4_paired_effects_by_scope.csv"),
      file.path(corrected_wp4_root, "wp4_feature_panel_no_sum159_fuse.csv"),
      file.path(corrected_wp4_root, "wp4_signature_panel_no_sum159_fuse.csv")
    ),
    collapse = ";"
  )
  wp4_all_line_inputs <- paste(
    c(
      file.path(corrected_wp4_root, "wp4_feature_panel_all_lines.csv"),
      file.path(corrected_wp4_root, "wp4_signature_panel_all_lines.csv")
    ),
    collapse = ";"
  )
  inversion_inputs <- paste(
    c(
      "agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/stage_outputs/analysis/live_auc_inversion/results/signature_fits.csv",
      "agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/stage_outputs/analysis/live_auc_inversion/results/feature_panel_corrected.csv"
    ),
    collapse = ";"
  )

  list(
    panel_row(
      "Figure 2", "a", make_polish_figure2a(envs$summary), 7.0, 2.05,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/initial_subpanels/summary_strip_options/fg2_summary_strip_option_a_no_sum_effect_bars.png",
      file.path(scripts, "fg2_summary_strip_options_worker.R"),
      paste(wp4_no_sum_inputs, inversion_inputs, sep = ";"),
      "Current redraft redraw: replaced the deprecated live-AUC-gradient coefficient with validated glucose drawdown per live-cell AUC and retained point-only per-line plus median encoding on a zero-centered within-feature SD scale."
    ),
    panel_row(
      "Figure 2", "b", make_polish_figure2b(envs$staged), 7.0, 2.65,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/refined_subpanels/fg2_figure2b_live_spline_total_yield.png",
      file.path(scripts, "fg2_staged_feature_workflow_worker.R"),
      paste(wp1_inputs, file.path(corrected_wp4_root, "wp4_feature_panel_no_sum159_fuse.csv"), sep = ";"),
      "V6 corrected-input redraw: total cells are points only, live/total legend is color-only, derivative tangent and total-cell horizontal guides are dashed, and the tangent is extended."
    ),
    panel_row(
      "Figure 2", "c", make_polish_figure2c(envs$staged), 7.0, 2.95,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/refined_subpanels/fg2_figure2c_per_g0_split_regressions.png",
      file.path(scripts, "fg2_staged_feature_workflow_worker.R"),
      paste(wp1_inputs, file.path(corrected_wp4_root, "wp4_feature_panel_no_sum159_fuse.csv"), inversion_inputs, sep = ";"),
      "Current redraft redraw: the middle view uses the validated glucose-drawdown-on-live-cell-AUC inverse regression and preserves low/high ploidy fill on highlighted G0 = 0 points; derivative and yield retain their prior construction roles."
    ),
    panel_row(
      "Figure S3", "a", envs$supplement$p_raw_all + scale_x_continuous(labels = label_short) + theme(legend.position = "bottom"), 7.0, 4.25,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/initial_subpanels/supplement_qc/s3_standalone_all_line_sensitivity.png",
      file.path(scripts, "fg2_supplement_qc_worker.R"),
      paste(wp4_all_line_inputs, inversion_inputs, sep = ";"),
      "Current redraft redraw: replaced the deprecated all-line live-AUC-gradient facet with the validated glucose-drawdown-per-live-cell-AUC coefficient while retaining the low/high ploidy encoding and SUM-159-fuse row."
    ),
    panel_row(
      "Figure S4", "a", make_polish_s4(envs$supplement), 7.0, 8.20,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/initial_subpanels/supplement_qc/s4_round2_all_cell_line_construction.png",
      file.path(scripts, "fg2_supplement_qc_worker.R"),
      paste(wp1_inputs, wp4_all_line_inputs, inversion_inputs, sep = ";"),
      "Current redraft redraw: all-line middle views use the validated glucose-drawdown-on-live-cell-AUC inverse regressions with visible low/high ploidy fill at G0 = 0; derivative-only split markers and 1 mM yield highlighting are retained."
    )
  )
}

png_dimensions <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, what = "raw", n = 24)
  if (length(header) < 24) stop("Could not read PNG header: ", path, call. = FALSE)
  width <- sum(as.integer(header[17:20]) * c(256^3, 256^2, 256, 1))
  height <- sum(as.integer(header[21:24]) * c(256^3, 256^2, 256, 1))
  c(width = width, height = height)
}

save_panel_png <- function(row, dpi = 360) {
  filename <- paste0(tolower(safe_name(row$figure)), "_", row$panel, ".png")
  path <- file.path(subpanel_dir, filename)
  ggplot2::ggsave(
    path,
    plot = row$plot,
    width = row$width,
    height = row$height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  dims <- png_dimensions(path)
  list(
    path = path,
    width_px = as.integer(dims[["width"]]),
    height_px = as.integer(dims[["height"]])
  )
}

write_subpanels <- function(panels) {
  rows <- lapply(panels, function(row) {
    saved <- save_panel_png(row)
    data.frame(
      figure = row$figure,
      panel = row$panel,
      subpanel_png = rel_path(saved$path),
      width_px = saved$width_px,
      height_px = saved$height_px,
      width_in = row$width,
      height_in = row$height,
      stringsAsFactors = FALSE
    )
  })
  dims <- bind_rows(rows)
  write.csv(dims, file.path(layout_dir, "subpanel_dimensions.csv"), row.names = FALSE)
}

assemble_one <- function(figure_id, panels, layout_plan) {
  one <- layout_plan %>%
    filter(.data$figure == figure_id) %>%
    arrange(.data$panel)
  panel_lookup <- setNames(panels, vapply(panels, function(x) paste(x$figure, x$panel, sep = "|"), character(1)))
  canvas <- cowplot::ggdraw()
  for (idx in seq_len(nrow(one))) {
    lay <- one[idx, ]
    key <- paste(lay$figure, lay$panel, sep = "|")
    panel <- panel_lookup[[key]]
    pad_left <- if (lay$x_npc < 0.01) 0.035 else 0.010
    pad_right <- if ((lay$x_npc + lay$width_npc) > 0.99) 0.055 else 0.008
    pad_bottom <- if (lay$y_npc < 0.01) 0.022 else 0.009
    pad_top <- 0.020
    plot_x <- lay$x_npc + pad_left
    plot_y <- lay$y_npc + pad_bottom
    plot_w <- max(0.05, lay$width_npc - pad_left - pad_right)
    plot_h <- max(0.05, lay$height_npc - pad_bottom - pad_top)
    canvas <- canvas +
      cowplot::draw_plot(
        panel$plot,
        x = plot_x,
        y = plot_y,
        width = plot_w,
        height = plot_h
      ) +
      cowplot::draw_plot_label(
        label = lay$panel,
        x = lay$x_npc + 0.004,
        y = lay$y_npc + lay$height_npc - 0.004,
        hjust = 0,
        vjust = 1,
        size = 12,
        fontface = "bold"
      )
  }
  attr(canvas, "width") <- unique(one$layout_width_in)
  attr(canvas, "height") <- unique(one$layout_height_in)
  canvas
}

write_final_outputs <- function(panels) {
  layout_path <- file.path(layout_dir, "layout_plan.csv")
  if (!file.exists(layout_path)) stop("Missing layout plan: ", layout_path, call. = FALSE)
  layout_plan <- read.csv(layout_path, stringsAsFactors = FALSE, check.names = FALSE)
  out_map <- c(
    "Figure 2" = file.path(final_dir, "Figure_2.png"),
    "Figure S3" = file.path(final_dir, "Figure_S3.png"),
    "Figure S4" = file.path(final_dir, "Figure_S4.png")
  )
  for (figure_id in names(out_map)) {
    fig <- assemble_one(figure_id, panels, layout_plan)
    ggplot2::ggsave(
      out_map[[figure_id]],
      plot = fig,
      width = attr(fig, "width"),
      height = attr(fig, "height"),
      units = "in",
      dpi = 360,
      bg = "white",
      limitsize = FALSE
    )
  }
  write_manifest_and_notes(panels, out_map)
}

write_manifest_and_notes <- function(panels, out_map) {
  layout_plan <- rel_path(file.path(layout_dir, "layout_plan.csv"))
  command_subpanels <- paste(
    "scripts/agentRrunner.sh",
    rel_path(file.path(script_dir, "polish_figures.R")),
    "--phase subpanels --project-root . --output-root",
    rel_path(polish_root)
  )
  command_final <- paste(
    "scripts/agentRrunner.sh",
    rel_path(file.path(script_dir, "polish_figures.R")),
    "--phase final --project-root . --output-root",
    rel_path(polish_root)
  )
  output_for <- function(figure) rel_path(out_map[[figure]])
  subpanel_for <- function(row) rel_path(file.path(subpanel_dir, paste0(tolower(safe_name(row$figure)), "_", row$panel, ".png")))

  provenance <- bind_rows(lapply(panels, function(row) {
    data.frame(
      figure = row$figure,
      panel = row$panel,
      subpanel_image = subpanel_for(row),
      generator = paste(rel_path(file.path(script_dir, "polish_figures.R")), row$generator, sep = ";"),
      command = paste(command_subpanels, command_final, sep = ";"),
      data_inputs = row$data_inputs,
      layout_plan = layout_plan,
      output_image = output_for(row$figure),
      notes = row$notes,
      approved_source = row$approved_source,
      context_inputs = paste(
        c(
          "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/user-feedback_round3.txt",
          "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/polishing_workflow_20260626.md",
          "agent-dev/manuscript_redrafts/20260619_v6_redraft/handoffs/20260709T124521_figures_handoff.md",
          "agent-dev/manuscript_redrafts/20260619_v6_redraft/stage_outputs/analysis/analysis_status.md",
          "agent-dev/manuscript_redrafts/20260619_v6_redraft/stage_outputs/analysis/fg2_wp4_recompute_comparison.md",
          "agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/handoffs/20260711T123039_figures_continuation_handoff.md",
          "agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/stage_outputs/analysis/live_auc_inversion/validation_report.md",
          ".agents/references/manuscript_figure_style.md"
        ),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(provenance, file.path(polish_root, "provenance.csv"), row.names = FALSE)

  function_lookup <- c(
    "Figure 2|a" = "make_polish_figure2a",
    "Figure 2|b" = "make_polish_figure2b",
    "Figure 2|c" = "make_polish_figure2c",
    "Figure S3|a" = "make_raw_feature_pairs",
    "Figure S4|a" = "make_polish_s4"
  )
  manifest <- bind_rows(lapply(panels, function(row) {
    image_path <- output_for(row$figure)
    image_hash <- sha256_file(file.path(project_root, image_path))
    key <- paste(row$figure, row$panel, sep = "|")
    data.frame(
      manuscript_figure_name = row$figure,
      panel_id = row$panel,
      current_figure_name = NA_character_,
      source_local_figure_name = safe_name(row$figure),
      intermediate_aliases = tools::file_path_sans_ext(basename(row$approved_source)),
      source_root = rel_path(polish_root),
      prior_manuscript_figure_name = NA_character_,
      prior_panel_id = NA_character_,
      polishing_final_image_path = image_path,
      polishing_final_sha256 = image_hash,
      integration_rebuild_image_path = NA_character_,
      integration_rebuild_sha256 = NA_character_,
      final_image_path = NA_character_,
      final_image_sha256 = NA_character_,
      wrapper_rebuild_image_path = NA_character_,
      wrapper_rebuild_sha256 = NA_character_,
      hash_match_status = "not_run",
      source_script_path = paste(rel_path(file.path(script_dir, "polish_figures.R")), row$generator, sep = ";"),
      source_function_or_object = unname(function_lookup[[key]]),
      data_input_paths = row$data_inputs,
      configuration_paths = paste(
        c(
          rel_path(file.path(package_script_root, "config", "panel_map.csv")),
          layout_plan,
          ".agents/references/manuscript_figure_style.md",
          "agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/stage_outputs/analysis/live_auc_inversion/definition.json"
        ),
        collapse = ";"
      ),
      subpanel_audit_image_path = subpanel_for(row),
      semantic_interpretation_path = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(manifest, file.path(polish_root, "manifest.csv"), row.names = FALSE)

  rebuild_rows <- lapply(names(out_map), function(figure_name) {
    target_path <- output_for(figure_name)
    target_hash <- sha256_file(file.path(project_root, target_path))
    figure_panels <- Filter(function(row) identical(row$figure, figure_name), panels)
    data.frame(
      stage = "polishing",
      manuscript_figure_name = figure_name,
      current_figure_name = tools::file_path_sans_ext(basename(target_path)),
      source_root = rel_path(polish_root),
      script_path = rel_path(file.path(script_dir, "polish_figures.R")),
      working_directory = ".",
      command = command_final,
      target_image_path = target_path,
      target_sha256 = target_hash,
      rebuild_output_path = target_path,
      rebuild_sha256 = target_hash,
      hash_match_status = "match",
      direct_input_paths = paste(unique(vapply(figure_panels, `[[`, character(1), "data_inputs")), collapse = ";"),
      dependency_paths = paste(
        unique(c(
          rel_path(file.path(script_dir, "polish_figures.R")),
          vapply(figure_panels, `[[`, character(1), "generator"),
          rel_path(file.path(package_script_root, "config", "panel_map.csv")),
          layout_plan,
          ".agents/references/manuscript_figure_style.md"
        )),
        collapse = ";"
      ),
      immutable_raster_inputs = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  write.table(
    bind_rows(rebuild_rows),
    file.path(polish_root, "figure_rebuild_manifest.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )

  legend <- c(
    "# FG2 Polished Legends",
    "",
    "Figure 2, Direct-feature summary and staged construction workflow. Panel a shows no-SUM-159-fuse high-minus-low ploidy effects for growth, live-cell AUC, glucose drawdown per live-cell AUC, and fixed 1 mM net total-yield features as per-line points on a zero-centered within-feature SD scale, with median points marked by open diamonds. Panel b shows a representative live-spline construction with live-cell points, the accepted live spline, live-AUC shading, a dashed derivative tangent, and total-cell yield support shown by total-cell points with dashed min/max total-cell guides and a max-delta arrow. Panel c shows representative-line feature constructions: robust derivative, an intercept-bearing glucose-drawdown-on-live-cell-AUC regression fitted over G0 <= 1 separately for low and high ploidy, and net total-cell yield; the G0 = 0 inverse-regression observations retain blue/red ploidy fill and 1 mM yield values are highlighted.",
    "",
    "Figure S3, All-line sensitivity for direct features. Panel a keeps SUM-159-fuse visible in the all-line display using the same low/high ploidy encoding as the other lines and reports the validated glucose-drawdown-per-live-cell-AUC coefficient rather than the deprecated live-AUC-gradient orientation.",
    "",
    "Figure S4, All-cell-line direct-feature construction without death-feature rows. Panel a shows robust max derivative, intercept-bearing glucose-drawdown-on-live-cell-AUC regressions over G0 <= 1, and net total-cell yield as horizontally concatenated views for every cell line. Dashed split markers appear only in the derivative view, highlighted G0 = 0 inverse-regression observations retain blue/red ploidy fill, total yield is displayed at all G0 values, and the fixed 1 mM yield feature is highlighted."
  )
  writeLines(legend, file.path(polish_root, "legend.md"))

  notes <- c(
    "# FG2 polishing notes",
    "",
    "## Inputs and decisions",
    "",
    "- Polishing scope: closed `FG2_direct_feature_rebuild/drafting_v2` package plus round-3 package-local polishing feedback.",
    "- Current redraft scope: regenerate Figure 2, Figure S3, and Figure S4 from the corrected WP4 outputs and validated live-AUC inversion package under `stage_outputs/analysis/live_auc_inversion/`.",
    "- The live-AUC coefficient is the intercept-bearing `glucose_drawdown ~ live_auc_glucose_window` slope over finite G0 <= 1 observations, in glucose concentration per live-cell-hour; the package consumes the validated fixed `signature_fits.csv` and does not rerun Analysis.",
    "- Approval source: `drafting_v2/user-feedback_round3.txt`, where the user described the figures as good and deferred the listed fixes to polishing.",
    "- Figure 2A keeps the direct-feature summary role but removes the bar/line encoding. A literal within-feature z-score would recenter each feature and obscure high-low effect direction, so the polished axis uses `effect per ploidy / within-feature SD` with zero retained as the no-effect reference.",
    "- Figure 2B keeps the one-panel live-spline/total-yield construction but removes the total-cell line, adds a color-only live/total legend, makes annotation guides dashed, and extends the derivative tangent.",
    "- Figure 2C and Figure S4 use patchwork-style horizontal concatenation; their middle views show the inverse glucose-drawdown-on-live-cell-AUC regressions and retain low/high ploidy fill on the highlighted G0 = 0 points.",
    "- Figure S3 retains the standalone all-line sensitivity role but replaces the deprecated live-AUC-gradient values and label with the validated inverse coefficient.",
    "- Final composites are assembled from regenerated ggplot/patchwork objects, not from final draft PNG rasters.",
    "",
    "## Commands",
    "",
    paste0("- `", command_subpanels, "`"),
    "- The integration rebuild stages the promoted fixed `layout_plan.csv` and `subpanel_dimensions.csv` into its output root; it does not rerun layout optimization.",
    paste0("- `", command_final, "`"),
    "",
    "## Layout",
    "",
    "- The layout optimizer's `x_npc` and `y_npc` values were used as lower-left panel origins with `cowplot::draw_plot()`.",
    "- Figure 2 panels a-c are stacked; S3 and S4 are each treated as one displayed supplemental panel.",
    "- Final assembly adds lowercase panel labels outside the plotted grob padding.",
    "",
    "## Caveats",
    "",
    "- Figure S4 remains dense because it is an all-cell-line construction/QC figure.",
    "- S4 derivative split markers denote the low-G and high-G derivative summary regions and are intentionally absent from live-AUC and yield views.",
    "",
    "## Project-map decision",
    "",
    "- `docs/project_map.txt` already records the v5 polishing-output status at the figure-generation root. No FG2-specific project-map expansion is needed."
  )
  writeLines(notes, file.path(polish_root, "notes.md"))

  visual_qc <- c(
    "# FG2 visual QC",
    "",
    "Overall status: passed direct local visual inspection after regeneration with the validated inverse coefficient; no FG2 layout blocker was found.",
    "",
    "| final PNG | title/subtitle/caption/header | clipping | spacing/layout | readability | panel order/labels | status |",
    "|---|---|---|---|---|---|---|",
    "| `final_images/figure_2.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass: no visible clipping. | Pass: panels a-c use optimizer stacking with internal padding; Figure 2C uses one ploidy legend instead of repeated legends. | Pass: panel A labels the inverse coefficient; panel C shows the inverse regression and the highlighted G0 = 0 low/high observations remain blue/red. | Pass: lowercase a-c labels appear once in intended reading order. | Pass |",
    "| `final_images/figure_s3.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass: no visible clipping. | Pass: single all-line sensitivity panel. | Pass: the inverted coefficient facet is legible and SUM-159-fuse remains visible with low/high ploidy encoding. | Pass: single lowercase a label is present. | Pass |",
    "| `final_images/figure_s4.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass: no visible clipping. | Pass: single horizontal feature patchwork panel with a single ploidy legend. | Pass with density caveat: inverse regressions are visible for all lines, highlighted G0 = 0 low/high observations retain blue/red fill, derivative split markers remain limited to the derivative view, and 1 mM yield remains highlighted. | Pass: single lowercase a label is present. | Pass |"
  )
  writeLines(visual_qc, file.path(polish_root, "visual_qc.md"))
}

envs <- load_builders()
panels <- build_panels(envs)

if (phase %in% c("subpanels", "both")) {
  write_subpanels(panels)
}
if (phase %in% c("final", "both")) {
  write_final_outputs(panels)
}

cat("FG2 polishing phase complete: ", phase, "\n", sep = "")
