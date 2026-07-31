#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

DRAFT_ROOT <- file.path(
  "agent-dev", "manuscript_redrafts", "20260709T210405_v7_redraft",
  "figure_generation", "FG2_direct_feature_rebuild", "drafting_v2"
)
INITIAL_DIR <- file.path(DRAFT_ROOT, "initial_subpanels", "staged_feature_workflow")
REFINED_PATH <- file.path(DRAFT_ROOT, "refined_subpanels", "fg2_staged_feature_workflow_contact_sheet.png")
REFINED_DIR <- dirname(REFINED_PATH)
NOTES_PATH <- file.path(DRAFT_ROOT, "worker_notes", "fg2_staged_feature_workflow_coverage.md")
ROUND3_NOTES_PATH <- file.path(DRAFT_ROOT, "worker_notes", "fg2_round3_figure2bc_coverage.md")

dir.create(INITIAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REFINED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(NOTES_PATH), recursive = TRUE, showWarnings = FALSE)

png_dpi <- 420
representative_line <- "MDA-MB-231"
g0_levels <- c("0", "0.1", "0.25", "0.5", "1", "5", "25")

ploidy_palette <- c(low = "#1B6CA8", high = "#C43C39")
feature_palette <- c(
  `Robust max derivative` = "#009E73",
  `Live-cell AUC` = "#7A5195",
  `Total-cell yield` = "#D55E00"
)

count_summary_path <- file.path("data", "report_exports", "wp1_core_starvation", "count_trajectory_summary.csv")
glucose_summary_path <- file.path("data", "report_exports", "wp1_core_starvation", "glucose_trajectory_summary.csv")
corrected_wp4_root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v6_redraft",
  "stage_outputs", "analysis", "recomputed_wp4_sum159_fuse_exclusion"
)
feature_panel_path <- file.path(corrected_wp4_root, "wp4_feature_panel_no_sum159_fuse.csv")
signature_panel_path <- file.path(corrected_wp4_root, "wp4_signature_panel_no_sum159_fuse.csv")
feedback_path <- file.path(
  "agent-dev", "manuscript_redrafts", "20260709T210405_v7_redraft",
  "handoffs", "20260709T124521_figures_handoff.md"
)

read_csv_tbl <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...))
}

write_csv_tbl <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(as.data.frame(x), path, row.names = FALSE)
  invisible(path)
}

save_png <- function(plot, path, width, height, dpi = png_dpi) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

pretty_axis <- function(x, accuracy = NULL, blank_missing = TRUE) {
  vapply(x, function(one) {
    if (!is.finite(one)) {
      return(if (blank_missing) "" else "NA")
    }
    abs_one <- abs(one)
    if (abs_one >= 1e6) {
      return(sprintf("%.1fM", one / 1e6))
    }
    if (abs_one >= 1e3) {
      return(sprintf("%.1fk", one / 1e3))
    }
    if (is.null(accuracy)) {
      accuracy <- if (abs(one) < 0.1) 0.001 else 0.01
    }
    number(one, accuracy = accuracy)
  }, character(1))
}

short_number <- function() {
  function(x) pretty_axis(x)
}

pretty_num <- function(x, accuracy = NULL) {
  pretty_axis(x, accuracy = accuracy, blank_missing = FALSE)
}

format_ploidy_abs <- function(x) {
  out <- ifelse(abs(x - round(x)) < 1e-8, sprintf("%.0fN", x), sprintf("%.1fN", x))
  out[!is.finite(x)] <- NA_character_
  out
}

g0_display_factor <- function(g0) {
  factor(as.character(g0), levels = g0_levels)
}

add_ploidy_display <- function(df) {
  out <- df %>%
    group_by(cellLine) %>%
    mutate(
      ploidy_state = if_else(ploidy_metric == min(ploidy_metric, na.rm = TRUE), "low", "high"),
      ploidy_abs_label = format_ploidy_abs(ploidy_abs),
      ploidy_display = paste0(ploidy_state, " ", ploidy_abs_label)
    ) %>%
    ungroup() %>%
    mutate(
      ploidy_state = factor(ploidy_state, levels = c("low", "high")),
      ploidy_display = factor(ploidy_display, levels = unique(ploidy_display))
    )
  if ("G0" %in% names(out)) {
    out <- out %>%
      mutate(
        G0_display = g0_display_factor(G0),
        G0_index = match(as.character(G0), g0_levels)
      )
  }
  out
}

wp2_theme <- function(base_size = 7.0) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey93", color = "grey72", linewidth = 0.25),
      strip.text = element_text(face = "bold", margin = margin(2, 2, 2, 2)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.height = unit(0.12, "in"),
      legend.key.width = unit(0.22, "in"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(face = "bold", size = rel(1.0), margin = margin(0, 0, 2, 0)),
      plot.subtitle = element_text(size = rel(0.82), color = "grey24", margin = margin(0, 0, 2, 0)),
      plot.caption = element_text(size = rel(0.72), color = "grey35", hjust = 0)
    )
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_lm_grid_log_g0 <- function(df, y_col = "value", n = 120L) {
  df <- df[is.finite(df$G0) & is.finite(df[[y_col]]), , drop = FALSE]
  if (nrow(df) < 2L || length(unique(df$G0)) < 2L) {
    return(tibble())
  }
  fit <- lm(stats::as.formula(paste(y_col, "~ log1p(G0)")), data = df)
  grid <- tibble(G0 = seq(min(df$G0, na.rm = TRUE), max(df$G0, na.rm = TRUE), length.out = n))
  grid$pred <- as.numeric(predict(fit, newdata = grid))
  grid
}

safe_lm_summary_log_g0 <- function(df) {
  df <- df[is.finite(df$G0) & is.finite(df$value), , drop = FALSE]
  if (nrow(df) < 2L || length(unique(df$G0)) < 2L) {
    med <- safe_median(df$value)
    return(tibble(intercept = med, slope = NA_real_, n = nrow(df), r_squared = NA_real_))
  }
  fit <- lm(value ~ log1p(G0), data = df)
  coef_fit <- coef(fit)
  tibble(
    intercept = unname(coef_fit[[1]]),
    slope = unname(coef_fit[[2]]),
    n = nrow(df),
    r_squared = summary(fit)$r.squared
  )
}

compute_spline_one <- function(df, value_col = "live_cells", spar = 0.62, n_eval = 260L, edge_frac = 0.06) {
  df <- df %>% arrange(hours)
  ok <- is.finite(df$hours) & is.finite(df[[value_col]])
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 2L) {
    return(list(smooth = tibble(), feature = tibble()))
  }

  x <- df$hours
  y <- log1p(pmax(df[[value_col]], 0))
  eval_hours <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n_eval)

  if (length(unique(x)) >= 4L) {
    fit <- try(stats::smooth.spline(x = x, y = y, spar = spar), silent = TRUE)
  } else {
    fit <- structure(list(), class = "try-error")
  }

  if (inherits(fit, "try-error")) {
    smooth_y <- approx(x = x, y = y, xout = eval_hours, ties = "ordered", rule = 2)$y
    derivative <- c(diff(smooth_y) / diff(eval_hours), NA_real_)
  } else {
    smooth_y <- predict(fit, x = eval_hours)$y
    derivative <- predict(fit, x = eval_hours, deriv = 1)$y
  }

  span <- diff(range(eval_hours, na.rm = TRUE))
  inner <- eval_hours >= min(eval_hours, na.rm = TRUE) + edge_frac * span &
    eval_hours <= max(eval_hours, na.rm = TRUE) - edge_frac * span
  finite_inner <- is.finite(derivative) & inner
  robust_max <- if (any(finite_inner)) {
    unname(quantile(derivative[finite_inner], probs = 0.9, na.rm = TRUE, names = FALSE))
  } else {
    NA_real_
  }
  top_decile <- is.finite(derivative) & inner & is.finite(robust_max) & derivative >= robust_max
  mark_hour <- if (any(top_decile)) {
    median(eval_hours[top_decile], na.rm = TRUE)
  } else if (any(finite_inner)) {
    eval_hours[which.max(derivative * finite_inner)]
  } else {
    NA_real_
  }

  meta_cols <- intersect(
    c("well_idx", "cellLine", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id", "has_starvation"),
    names(df)
  )
  meta <- df[1, meta_cols, drop = FALSE]

  list(
    smooth = bind_cols(
      meta[rep(1L, length(eval_hours)), , drop = FALSE],
      tibble(
        hours = eval_hours,
        smooth_log_live = smooth_y,
        derivative = derivative,
        is_inner = inner,
        is_top_decile = top_decile
      )
    ),
    feature = bind_cols(
      meta,
      tibble(
        robust_max_derivative = robust_max,
        robust_derivative_time = mark_hour,
        robust_rule = "90th percentile of inner spline derivative",
        n_derivative_points = sum(finite_inner),
        n_top_decile_points = sum(top_decile)
      )
    )
  )
}

build_spline_bundle <- function(count_summary, feature_panel) {
  count_keep <- count_summary %>%
    semi_join(feature_panel %>% distinct(well_idx), by = "well_idx")
  pieces <- count_keep %>%
    split(.$well_idx) %>%
    lapply(compute_spline_one)

  smooth <- bind_rows(lapply(pieces, `[[`, "smooth"))
  robust <- bind_rows(lapply(pieces, `[[`, "feature"))

  list(smooth = smooth, robust = robust)
}

make_feature_values <- function(feature_panel, robust_panel) {
  feature_panel %>%
    select(
      well_idx, cellLine, line_id, ploidy_metric, ploidy_abs, G0, exp_id, has_starvation,
      live_auc_glucose_window,
      glucose_start_time, glucose_end_time, total_initial_at_glucose_start,
      total_peak_to_glucose_end, total_net_gain_to_glucose_end, glucose_initial,
      glucose_final, glucose_drawdown
    ) %>%
    left_join(
      robust_panel %>% select(well_idx, robust_max_derivative, robust_derivative_time, robust_rule),
      by = "well_idx"
    ) %>%
    mutate(peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)) %>%
    pivot_longer(
      cols = c(robust_max_derivative, live_auc_glucose_window, peak_total_yield_net),
      names_to = "feature_id",
      values_to = "value"
    ) %>%
    mutate(
      feature_label = recode(
        feature_id,
        robust_max_derivative = "Robust max derivative",
        live_auc_glucose_window = "Live-cell AUC",
        peak_total_yield_net = "Total-cell yield"
      ),
      units = recode(
        feature_id,
        robust_max_derivative = "d log1p(live) / hour",
        live_auc_glucose_window = "live cell-hours",
        peak_total_yield_net = "total cells"
      ),
      feature_label = factor(feature_label, levels = names(feature_palette))
    )
}

summarize_line_features <- function(feature_values, line_name = representative_line) {
  x <- feature_values %>%
    filter(cellLine == line_name, is.finite(value)) %>%
    add_ploidy_display() %>%
    mutate(log1p_G0 = log1p(G0))

  deriv_summary <- x %>%
    filter(feature_id == "robust_max_derivative") %>%
    group_by(cellLine, ploidy_state, ploidy_display) %>%
    summarise(
      `Low-G median` = safe_median(value[G0 <= 0.25]),
      `High-G median` = safe_median(value[G0 >= 1]),
      .groups = "drop"
    ) %>%
    pivot_longer(c(`Low-G median`, `High-G median`), names_to = "summary_label", values_to = "summary_value") %>%
    mutate(
      feature_label = "Robust max derivative",
      summary_rule = if_else(
        summary_label == "Low-G median",
        "Median per-G0 robust derivatives across G0 <= 0.25",
        "Median per-G0 robust derivatives across G0 >= 1"
      )
    )

  auc_summary <- x %>%
    filter(feature_id == "live_auc_glucose_window", G0 <= 1) %>%
    group_by(cellLine, feature_id, feature_label, ploidy_state, ploidy_display) %>%
    group_modify(~safe_lm_summary_log_g0(.x)) %>%
    ungroup() %>%
    transmute(
      cellLine,
      ploidy_state,
      ploidy_display,
      feature_label = as.character(feature_label),
      `Regression intercept` = intercept,
      `Regression slope` = slope,
      summary_rule = paste0("G0 <= 1 lm(live-cell AUC ~ log1p(G0)), displayed on a linear G0 axis; n=", n)
    ) %>%
    pivot_longer(c(`Regression intercept`, `Regression slope`), names_to = "summary_label", values_to = "summary_value")

  yield_summary <- x %>%
    filter(feature_id == "peak_total_yield_net", abs(G0 - 1) < 1e-8) %>%
    group_by(cellLine, feature_label, ploidy_state, ploidy_display) %>%
    summarise(
      `1 mM net yield` = safe_median(value),
      .groups = "drop"
    ) %>%
    mutate(summary_rule = "Net peak total-cell yield at G0 = 1: max(total_peak_to_glucose_end - total_initial_at_glucose_start, 0).") %>%
    pivot_longer(`1 mM net yield`, names_to = "summary_label", values_to = "summary_value")

  bind_rows(deriv_summary, auc_summary, yield_summary) %>%
    mutate(
      feature_label = factor(feature_label, levels = names(feature_palette)),
      summary_display = case_when(
        feature_label == "Robust max derivative" & summary_label == "Low-G median" ~ "Derivative\nlow-G median",
        feature_label == "Robust max derivative" & summary_label == "High-G median" ~ "Derivative\nhigh-G median",
        feature_label == "Live-cell AUC" & summary_label == "Regression intercept" ~ "Live AUC\nintercept",
        feature_label == "Live-cell AUC" & summary_label == "Regression slope" ~ "Live AUC\ngradient",
        feature_label == "Total-cell yield" & summary_label == "1 mM net yield" ~ "1 mM\nnet yield",
        TRUE ~ paste(feature_label, summary_label, sep = "\n")
      ),
      summary_display = factor(
        summary_display,
        levels = c(
          "Derivative\nlow-G median",
          "Derivative\nhigh-G median",
          "Live AUC\nintercept",
          "Live AUC\ngradient",
          "1 mM\nnet yield"
        )
      )
    )
}

pick_example <- function(feature_panel) {
  x <- feature_panel %>%
    add_ploidy_display() %>%
    filter(
      cellLine == representative_line,
      is.finite(live_auc_glucose_window),
      is.finite(total_net_gain_to_glucose_end)
    )
  preferred <- x %>%
    filter(ploidy_state == "low", G0 == 1) %>%
    slice(1)
  if (nrow(preferred)) {
    return(preferred)
  }
  x %>%
    arrange(abs(G0 - 1), ploidy_metric) %>%
    slice(1)
}

make_spline_fit_panel <- function(
    example_count,
    example_smooth,
    example_feature,
    base_size = 6.1) {
  observed <- example_count %>%
    arrange(hours) %>%
    mutate(
      log_live = log1p(pmax(live_cells, 0)),
      log_total = log1p(pmax(total_cells, 0))
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
  tangent_x <- seq(max(x_rng[1], robust_time - 14), min(x_rng[2], robust_time + 14), length.out = 2)
  tangent <- tibble(
    hours = tangent_x,
    log_live = y0 + robust_value * (tangent_x - robust_time)
  )

  start_time <- example_feature$glucose_start_time[[1]]
  end_time <- example_feature$glucose_end_time[[1]]
  auc_band <- example_smooth %>%
    filter(hours >= start_time, hours <= end_time) %>%
    mutate(ymin = min(observed$log_live, na.rm = TRUE) - 0.15)

  yield_window <- observed %>%
    filter(hours >= start_time, hours <= end_time)
  if (!nrow(yield_window)) {
    yield_window <- observed
  }
  yield_min <- min(yield_window$total_cells, na.rm = TRUE)
  yield_max <- max(yield_window$total_cells, na.rm = TRUE)
  yield_arrow_x <- max(yield_window$hours, na.rm = TRUE) - 0.09 * diff(range(yield_window$hours, na.rm = TRUE))
  yield_label_x <- yield_arrow_x + 0.03 * diff(range(observed$hours, na.rm = TRUE))
  yield_label_y <- mean(log1p(c(yield_min, yield_max)), na.rm = TRUE)

  ggplot() +
    geom_ribbon(
      data = auc_band,
      aes(hours, ymin = ymin, ymax = smooth_log_live),
      fill = feature_palette[["Live-cell AUC"]],
      alpha = 0.15
    ) +
    geom_line(data = observed, aes(hours, log_total), color = feature_palette[["Total-cell yield"]], linewidth = 0.48, alpha = 0.78) +
    geom_point(data = observed, aes(hours, log_total), color = feature_palette[["Total-cell yield"]], size = 0.48, alpha = 0.62) +
    geom_hline(
      yintercept = log1p(c(yield_min, yield_max)),
      color = feature_palette[["Total-cell yield"]],
      linewidth = 0.3,
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
    annotate(
      "text",
      x = max(yield_window$hours, na.rm = TRUE) - 0.04 * diff(range(observed$hours, na.rm = TRUE)),
      y = log1p(yield_max) + 0.06,
      label = "total cells",
      color = feature_palette[["Total-cell yield"]],
      size = 1.9,
      hjust = 1,
      vjust = 0
    ) +
    geom_point(data = observed, aes(hours, log_live), color = "grey25", size = 0.58, alpha = 0.82) +
    geom_line(data = example_smooth, aes(hours, smooth_log_live), color = "#1B6CA8", linewidth = 0.76) +
    geom_line(data = tangent, aes(hours, log_live), color = feature_palette[["Robust max derivative"]], linewidth = 0.72) +
    geom_point(
      data = tibble(hours = robust_time, log_live = y0),
      aes(hours, log_live),
      shape = 21,
      fill = "white",
      color = feature_palette[["Robust max derivative"]],
      stroke = 0.42,
      size = 1.45
    ) +
    annotate(
      "text",
      x = robust_time,
      y = y0 + 0.24,
      label = "derivative tangent",
      color = feature_palette[["Robust max derivative"]],
      size = 1.9,
      hjust = 0.5
    ) +
    annotate(
      "text",
      x = mean(c(start_time, end_time), na.rm = TRUE),
      y = min(observed$log_live, na.rm = TRUE),
      label = "live AUC",
      color = feature_palette[["Live-cell AUC"]],
      size = 1.9,
      hjust = 0.5
    ) +
    labs(x = "hours", y = "log1p(cells)", title = NULL) +
    wp2_theme(base_size = base_size) +
    theme(legend.position = "none")
}

make_total_yield_panel <- function(example_count, example_feature, base_size = 6.1) {
  start_time <- example_feature$glucose_start_time[[1]]
  end_time <- example_feature$glucose_end_time[[1]]
  observed <- example_count %>%
    filter(hours >= start_time, hours <= end_time) %>%
    arrange(hours)
  if (!nrow(observed)) {
    observed <- example_count %>% arrange(hours)
  }

  eval_hours <- seq(min(observed$hours, na.rm = TRUE), max(observed$hours, na.rm = TRUE), length.out = 180L)
  fit <- try(stats::smooth.spline(observed$hours, log1p(pmax(observed$total_cells, 0)), spar = 0.62), silent = TRUE)
  smooth_total <- if (inherits(fit, "try-error")) {
    tibble(
      hours = eval_hours,
      smooth_total = expm1(approx(observed$hours, log1p(pmax(observed$total_cells, 0)), xout = eval_hours, rule = 2)$y)
    )
  } else {
    tibble(hours = eval_hours, smooth_total = expm1(predict(fit, x = eval_hours)$y))
  }

  peak_idx <- which.max(observed$total_cells)
  start_mark <- tibble(
    point = "start",
    hours = start_time,
    total_cells = example_feature$total_initial_at_glucose_start[[1]]
  )
  peak_mark <- tibble(
    point = "peak",
    hours = observed$hours[[peak_idx]],
    total_cells = example_feature$total_peak_to_glucose_end[[1]]
  )
  marks <- bind_rows(start_mark, peak_mark)
  label_text <- paste0(
    "yield = max(delta total, 0)\n/ glucose drawdown"
  )

  ggplot() +
    geom_point(data = observed, aes(hours, total_cells), color = "grey25", size = 0.56, alpha = 0.82) +
    geom_line(data = smooth_total, aes(hours, smooth_total), color = feature_palette[["Total-cell yield"]], linewidth = 0.72) +
    geom_segment(
      data = marks %>% summarise(
        x = start_mark$hours[[1]],
        xend = peak_mark$hours[[1]],
        y = start_mark$total_cells[[1]],
        yend = peak_mark$total_cells[[1]]
      ),
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.045, "in")),
      color = feature_palette[["Total-cell yield"]],
      linewidth = 0.55
    ) +
    geom_point(data = marks, aes(hours, total_cells), shape = 21, fill = "white", color = "grey18", size = 1.15, stroke = 0.28) +
    annotate(
      "label",
      x = min(observed$hours, na.rm = TRUE) + 0.04 * diff(range(observed$hours, na.rm = TRUE)),
      y = max(observed$total_cells, na.rm = TRUE),
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 1.7,
      label.size = 0.14,
      fill = "white"
    ) +
    scale_y_continuous(labels = short_number()) +
    labs(x = "hours", y = "total cells", title = "example: total-cell yield") +
    wp2_theme(base_size = base_size) +
    theme(legend.position = "none")
}

make_derivative_panel <- function(example_smooth, example_spline_feature, base_size = 6.2) {
  robust_value <- example_spline_feature$robust_max_derivative[[1]]
  robust_time <- example_spline_feature$robust_derivative_time[[1]]
  y_rng <- range(example_smooth$derivative, robust_value, na.rm = TRUE)
  x_rng <- range(example_smooth$hours, na.rm = TRUE)

  ggplot(example_smooth, aes(hours, derivative)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.2) +
    geom_line(color = "grey38", linewidth = 0.38) +
    geom_point(
      data = example_smooth %>% filter(is_top_decile),
      aes(hours, derivative),
      color = feature_palette[["Robust max derivative"]],
      size = 0.72,
      alpha = 0.85
    ) +
    geom_hline(yintercept = robust_value, color = feature_palette[["Robust max derivative"]], linewidth = 0.32, linetype = "dashed") +
    geom_point(
      data = tibble(hours = robust_time, derivative = robust_value),
      aes(hours, derivative),
      shape = 24,
      fill = "white",
      color = feature_palette[["Robust max derivative"]],
      stroke = 0.35,
      size = 1.4
    ) +
    annotate(
      "label",
      x = x_rng[2] - 0.02 * diff(x_rng),
      y = robust_value,
      label = paste0("robust max\n", pretty_num(robust_value)),
      hjust = 1,
      vjust = -0.12,
      size = 1.9,
      label.size = 0.14,
      fill = "white"
    ) +
    coord_cartesian(ylim = y_rng + c(-0.08, 0.14) * diff(y_rng)) +
    labs(x = "hours", y = "spline derivative", title = "2a robust max derivative") +
    wp2_theme(base_size = base_size) +
    theme(legend.position = "none")
}

make_auc_panel <- function(example_count, example_feature, base_size = 6.2) {
  start_time <- example_feature$glucose_start_time[[1]]
  end_time <- example_feature$glucose_end_time[[1]]
  window <- example_count %>%
    filter(hours >= start_time, hours <= end_time) %>%
    arrange(hours)
  auc_poly <- bind_rows(
    window %>% slice(1) %>% mutate(live_cells = 0),
    window,
    window %>% slice(n()) %>% mutate(live_cells = 0)
  )

  ggplot(window, aes(hours, live_cells)) +
    geom_polygon(data = auc_poly, aes(hours, live_cells), fill = feature_palette[["Live-cell AUC"]], alpha = 0.26, color = NA) +
    geom_line(color = feature_palette[["Live-cell AUC"]], linewidth = 0.46) +
    geom_point(color = feature_palette[["Live-cell AUC"]], size = 0.62) +
    geom_vline(xintercept = c(start_time, end_time), color = "grey62", linewidth = 0.22, linetype = "dotted") +
    annotate(
      "label",
      x = min(window$hours, na.rm = TRUE) + 0.04 * diff(range(window$hours, na.rm = TRUE)),
      y = max(window$live_cells, na.rm = TRUE),
      label = paste0("AUC ", pretty_num(example_feature$live_auc_glucose_window[[1]])),
      hjust = 0,
      vjust = 1,
      size = 1.9,
      label.size = 0.14,
      fill = "white"
    ) +
    labs(x = "hours", y = "live cells", title = "2b live-cell AUC") +
    scale_y_continuous(labels = short_number()) +
    wp2_theme(base_size = base_size) +
    theme(legend.position = "none")
}

make_yield_panel <- function(example_count, example_glucose, example_feature, base_size = 6.2) {
  start_time <- example_feature$glucose_start_time[[1]]
  end_time <- example_feature$glucose_end_time[[1]]
  count_window <- example_count %>%
    filter(hours >= start_time, hours <= end_time) %>%
    arrange(hours)
  glucose_window <- example_glucose %>%
    filter(hours >= start_time, hours <= end_time) %>%
    arrange(hours)

  peak_idx <- if (nrow(count_window)) which.max(count_window$total_cells) else integer()
  total_marks <- tibble(
    measurement = "total cells",
    hours = c(start_time, if (length(peak_idx)) count_window$hours[peak_idx] else NA_real_),
    value = c(example_feature$total_initial_at_glucose_start[[1]], example_feature$total_peak_to_glucose_end[[1]]),
    point = c("start", "peak")
  )
  glucose_marks <- tibble(
    measurement = "glucose",
    hours = c(start_time, end_time),
    value = c(example_feature$glucose_initial[[1]], example_feature$glucose_final[[1]]),
    point = c("start", "end")
  )
  marks <- bind_rows(total_marks, glucose_marks) %>%
    filter(is.finite(hours), is.finite(value))

  long <- bind_rows(
    count_window %>% transmute(measurement = "total cells", hours, value = total_cells),
    glucose_window %>% transmute(measurement = "glucose", hours, value = glucose_hat)
  ) %>%
    mutate(measurement = factor(measurement, levels = c("total cells", "glucose")))

  yield_value <- pmax(example_feature$total_net_gain_to_glucose_end[[1]], 0)
  yield_label <- paste0("yield ", pretty_num(yield_value), "\nmax(delta total,0)")
  yield_label_df <- tibble(
    measurement = factor("total cells", levels = c("total cells", "glucose")),
    hours = min(long$hours, na.rm = TRUE) + 0.04 * diff(range(long$hours, na.rm = TRUE)),
    value = max(count_window$total_cells, na.rm = TRUE),
    label = yield_label
  )

  ggplot(long, aes(hours, value)) +
    geom_line(color = feature_palette[["Total-cell yield"]], linewidth = 0.42) +
    geom_point(color = feature_palette[["Total-cell yield"]], size = 0.52, alpha = 0.9) +
    geom_point(data = marks, aes(hours, value, shape = point), fill = "white", color = "grey20", size = 1.25, stroke = 0.3) +
    facet_wrap(~measurement, ncol = 1, scales = "free_y") +
    scale_y_continuous(labels = short_number()) +
    scale_shape_manual(values = c(start = 21, peak = 24, end = 25), guide = "none") +
    geom_label(
      data = yield_label_df,
      aes(hours, value, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 1.75,
      label.size = 0.14,
      fill = "white"
    ) +
    labs(x = "hours", y = NULL, title = "2c total-cell yield") +
    wp2_theme(base_size = base_size) +
    theme(legend.position = "none", strip.text = element_text(size = rel(0.82)))
}

make_all_g0_grid <- function(feature_values, line_name = representative_line, compact = FALSE, example_well = NULL) {
  raw_values <- feature_values %>%
    filter(cellLine == line_name, is.finite(value)) %>%
    add_ploidy_display() %>%
    mutate(
      feature_label = factor(as.character(feature_label), levels = names(feature_palette)),
      plot_feature_label = if (compact) {
        recode(
          as.character(feature_label),
          `Robust max derivative` = "Derivative",
          `Live-cell AUC` = "Live AUC",
          `Total-cell yield` = "1 mM yield"
        )
      } else {
        as.character(feature_label)
      },
      plot_feature_label = factor(
        plot_feature_label,
        levels = if (compact) c("Derivative", "Live AUC", "1 mM yield") else names(feature_palette)
      ),
      highlighted_derivative_region = feature_id == "robust_max_derivative" & (G0 <= 0.25 | G0 >= 1),
      highlighted_yield_point = feature_id == "peak_total_yield_net" & abs(G0 - 1) < 1e-8,
      example_point = !is.null(example_well) & well_idx == example_well
    )

  line_values <- raw_values %>%
    filter(case_when(
      feature_id == "live_auc_glucose_window" ~ G0 <= 1,
      feature_id == "peak_total_yield_net" ~ abs(G0 - 1) < 1e-8,
      TRUE ~ TRUE
    ))

  fit_lines <- line_values %>%
    filter(feature_id == "live_auc_glucose_window") %>%
    group_by(plot_feature_label, ploidy_state) %>%
    group_modify(~safe_lm_grid_log_g0(.x, y_col = "value")) %>%
    ungroup()

  split_df <- line_values %>%
    filter(feature_id == "robust_max_derivative") %>%
    distinct(plot_feature_label) %>%
    mutate(xintercept = mean(c(0.25, 1)))

  p <- ggplot(line_values, aes(G0, value, color = ploidy_state, group = ploidy_state)) +
    geom_line(
      data = line_values %>% filter(feature_id == "robust_max_derivative"),
      linewidth = 0.28,
      alpha = 0.42
    ) +
    geom_line(
      data = fit_lines,
      aes(G0, pred, group = ploidy_state),
      linewidth = 0.44,
      alpha = 0.82
    ) +
    geom_vline(
      data = split_df,
      aes(xintercept = xintercept),
      linetype = "dashed",
      linewidth = 0.22,
      color = "grey35"
    ) +
    geom_point(size = if (compact) 0.88 else 1.15, alpha = 0.92) +
    geom_point(
      data = line_values %>% filter(highlighted_derivative_region | highlighted_yield_point),
      aes(G0, value),
      shape = 21,
      fill = NA,
      color = "black",
      stroke = 0.35,
      size = if (compact) 1.45 else 1.8
    ) +
    geom_point(
      data = line_values %>% filter(example_point),
      aes(G0, value),
      shape = 21,
      fill = NA,
      color = "black",
      stroke = 0.65,
      size = if (compact) 1.8 else 2.1
    ) +
    facet_wrap(~plot_feature_label, scales = "free", ncol = if (compact) 3 else 1) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_x_continuous(
      breaks = as.numeric(g0_levels),
      labels = if (compact) c("0", "", ".25", "", "1", "5", "25") else g0_levels,
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_y_continuous(labels = short_number()) +
    labs(
      x = "starting glucose G0",
      y = "per-G0 feature value",
      title = NULL
    ) +
    wp2_theme(base_size = if (compact) 5.8 else 6.5) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = if (compact) 5.1 else 5.9),
      panel.spacing = unit(if (compact) 0.25 else 0.34, "lines"),
      axis.text.x = element_text(
        size = if (compact) 4.3 else 5.2,
        angle = if (compact) 45 else 0,
        hjust = if (compact) 1 else 0.5
      )
    )

  p
}

make_summary_pair_panel <- function(line_summary, ncol = 3L, base_size = 6.2, show_title = TRUE) {
  x <- line_summary %>%
    filter(is.finite(summary_value)) %>%
    mutate(summary_display = factor(summary_display, levels = levels(line_summary$summary_display)))

  pair_df <- x %>%
    group_by(feature_label, summary_display, summary_label) %>%
    summarise(
      low_value = summary_value[ploidy_state == "low"][1],
      high_value = summary_value[ploidy_state == "high"][1],
      .groups = "drop"
    )

  ggplot() +
    geom_segment(
      data = pair_df,
      aes(x = low_value, xend = high_value, y = 1, yend = 1),
      color = "grey54",
      linewidth = 0.36
    ) +
    geom_point(
      data = x,
      aes(summary_value, 1, fill = ploidy_state, shape = ploidy_state),
      color = "grey16",
      size = 1.65,
      stroke = 0.26
    ) +
    facet_wrap(~summary_display, scales = "free_x", ncol = ncol) +
    scale_fill_manual(values = ploidy_palette, name = "ploidy") +
    scale_shape_manual(values = c(low = 21, high = 24), name = "ploidy") +
    scale_x_continuous(labels = short_number()) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(
      x = "final feature summary value",
      title = if (show_title) "3 aggregate per-G0 values to the cell-line summary" else NULL
    ) +
    wp2_theme(base_size = base_size) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = rel(0.82)),
      panel.grid.major.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = rel(0.78))
    )
}

make_feature_rule_table <- function(base_size = 6.2) {
  rules <- tibble(
    y = factor(
      c("Robust max derivative", "Live-cell AUC", "Total-cell yield"),
      levels = rev(c("Robust max derivative", "Live-cell AUC", "Total-cell yield"))
    ),
    x = 1,
    label = c(
      "spline derivative -> 90th percentile inside time range",
      "live-cell AUC across glucose window, G0 <= 1",
      "max(delta total cells,0) at G0 = 1"
    ),
    fill = names(feature_palette)
  )

  ggplot(rules, aes(x, y)) +
    geom_tile(aes(fill = fill), width = 0.92, height = 0.72, alpha = 0.16, color = "grey74", linewidth = 0.22) +
    geom_text(aes(label = label), hjust = 0.5, size = 1.9, lineheight = 0.9) +
    scale_fill_manual(values = feature_palette, guide = "none") +
    scale_x_continuous(NULL, breaks = NULL, limits = c(0.5, 1.5)) +
    labs(y = NULL, title = "feature rules") +
    wp2_theme(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = rel(0.82)),
      axis.ticks = element_blank()
    )
}

count_summary <- read_csv_tbl(count_summary_path)
glucose_summary <- read_csv_tbl(glucose_summary_path)
feature_panel <- read_csv_tbl(feature_panel_path)
signature_panel <- read_csv_tbl(signature_panel_path)

spline_bundle <- build_spline_bundle(count_summary, feature_panel)
feature_values <- make_feature_values(feature_panel, spline_bundle$robust)
line_summary <- summarize_line_features(feature_values, representative_line)

example_feature <- pick_example(feature_panel) %>%
  left_join(
    spline_bundle$robust %>% select(well_idx, robust_max_derivative, robust_derivative_time, robust_rule),
    by = "well_idx"
  ) %>%
  add_ploidy_display()
example_well <- example_feature$well_idx[[1]]
example_count <- count_summary %>%
  filter(well_idx == example_well)
example_glucose <- glucose_summary %>%
  filter(well_idx == example_well)
example_smooth <- spline_bundle$smooth %>%
  filter(well_idx == example_well)
example_spline_feature <- spline_bundle$robust %>%
  filter(well_idx == example_well)
stage_spline <- make_spline_fit_panel(
  example_count,
  example_smooth,
  example_feature
)
stage_example_construction <- stage_spline
stage_derivative <- make_derivative_panel(example_smooth, example_spline_feature)
stage_auc <- make_auc_panel(example_count, example_feature)
stage_yield <- make_yield_panel(example_count, example_glucose, example_feature)
stage_rules <- make_feature_rule_table()
stage_grid <- make_all_g0_grid(feature_values, representative_line, compact = FALSE, example_well = example_well)
stage_grid_compact <- make_all_g0_grid(feature_values, representative_line, compact = TRUE, example_well = example_well)
stage_summary <- make_summary_pair_panel(line_summary, ncol = 3L)
stage_summary_compact <- make_summary_pair_panel(line_summary, ncol = 3L, base_size = 5.7, show_title = FALSE)

single_trajectory_features <- (stage_derivative | stage_auc | stage_yield) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

option_a <- wrap_elements(full = stage_example_construction) /
  wrap_elements(full = single_trajectory_features) /
  wrap_elements(full = stage_grid) +
  plot_layout(heights = c(0.92, 1.12, 1.02), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

option_b <- (stage_grid | stage_summary) +
  plot_layout(widths = c(1.15, 1.0), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

option_c <- (wrap_elements(full = stage_example_construction) | wrap_elements(full = stage_grid_compact)) +
  plot_layout(widths = c(1.02, 1.28), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

promoted_figure2_construction <- (wrap_elements(full = stage_example_construction) | wrap_elements(full = stage_grid_compact)) +
  plot_layout(widths = c(1.02, 1.28), guides = "collect") &
  theme(legend.position = "bottom")

panel_bc_pair <- (wrap_elements(full = stage_example_construction) | wrap_elements(full = stage_grid_compact)) +
  plot_layout(widths = c(1.02, 1.28), guides = "collect") +
  plot_annotation(
    tag_levels = list(c("b", "c")),
    theme = theme(plot.tag = element_text(size = 9, face = "bold"))
  ) &
  theme(legend.position = "bottom")

contact_label <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = label, hjust = 0, vjust = 0.5, size = 3.0, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, 4, 0, 4))
}

option_paths <- tibble(
  option_id = c(
    "three_stage_single_trajectory_to_summary",
    "all_g0_feature_grid_plus_summary",
    "compact_figure2_construction_strip",
    "promoted_figure2_construction_strip_no_tags",
    "figure2b_live_spline_total_yield",
    "figure2c_per_g0_split_regressions",
    "figure2bc_pair_with_tags"
  ),
  png_path = file.path(
    INITIAL_DIR,
    c(
      "fg2_staged_workflow_option_a_three_stage.png",
      "fg2_staged_workflow_option_b_all_g0_summary_grid.png",
      "fg2_staged_workflow_option_c_compact_figure2_strip.png",
      "fg2_staged_workflow_promoted_figure2_construction_strip.png",
      "fg2_figure2b_live_spline_total_yield.png",
      "fg2_figure2c_per_g0_split_regressions.png",
      "fg2_figure2bc_construction_pair.png"
    )
  ),
  width_in = c(7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0),
  height_in = c(8.1, 6.7, 3.25, 3.25, 3.15, 3.15, 3.25),
  description = c(
    "Full workflow: annotated example trajectory construction, single-trajectory feature extraction checks, and all-G0 feature values.",
    "All starting glucose feature grid beside the final summaries for the representative cell line.",
    "Compact construction strip candidate for pairing with a results-first Figure 2A-like display; separates the live-spline/total-yield example from the per-G0 regression grid.",
    "Promoted compact Figure 2 construction strip without internal patchwork tags; live spline plus total-yield example beside the per-G0 scatter/regression panel.",
    "Round-3 Figure 2B candidate: one example trajectory with live AUC shading, derivative tangent, total-cell-yield overlay, dashed min/max total-cell lines, and a vertical max-delta arrow.",
    "Round-3 Figure 2C candidate: per-G0 split/regression panel moved out of Figure 2B; live AUC shows only G0 <= 1 with the canonical log-G0 regression on a linear G0 axis, and yield is the net total-cell gain at G0 = 1 with no regression.",
    "Round-3 Figure 2B/C paired review strip with lowercase b/c patchwork tags."
  )
)

plots <- list(
  option_a,
  option_b,
  option_c,
  promoted_figure2_construction,
  stage_example_construction,
  stage_grid_compact,
  panel_bc_pair
)
for (i in seq_along(plots)) {
  save_png(plots[[i]], option_paths$png_path[[i]], option_paths$width_in[[i]], option_paths$height_in[[i]])
}

refined_paths <- tibble(
  option_id = c(
    "refined_figure2b_live_spline_total_yield",
    "refined_figure2c_per_g0_split_regressions",
    "refined_figure2bc_pair_with_tags"
  ),
  png_path = file.path(
    REFINED_DIR,
    c(
      "fg2_figure2b_live_spline_total_yield.png",
      "fg2_figure2c_per_g0_split_regressions.png",
      "fg2_figure2bc_construction_pair.png"
    )
  ),
  width_in = c(7.0, 7.0, 7.0),
  height_in = c(3.15, 3.15, 3.25),
  description = c(
    "Refined round-3 Figure 2B candidate generated directly from the staged workflow worker.",
    "Refined round-3 Figure 2C candidate generated directly from the staged workflow worker.",
    "Refined round-3 Figure 2B/C paired review strip generated directly from the staged workflow worker."
  )
)
refined_plots <- list(stage_example_construction, stage_grid_compact, panel_bc_pair)
for (i in seq_along(refined_plots)) {
  save_png(refined_plots[[i]], refined_paths$png_path[[i]], refined_paths$width_in[[i]], refined_paths$height_in[[i]])
}

contact_sheet <-
  contact_label("Option A: full three-stage workflow") /
  wrap_elements(full = option_a) /
  contact_label("Option B: all-G0 feature grid plus final summaries") /
  wrap_elements(full = option_b) /
  contact_label("Option C: compact construction strip for Figure 2") /
  wrap_elements(full = option_c) +
  plot_layout(heights = c(0.04, 1.0, 0.04, 0.70, 0.04, 0.48), guides = "collect") &
  theme(legend.position = "bottom")
save_png(contact_sheet, REFINED_PATH, 7.0, 16.5)

feature_values_out <- feature_values %>%
  mutate(
    source_count_summary = count_summary_path,
    source_feature_panel = feature_panel_path,
    robust_derivative_note = if_else(
      feature_id == "robust_max_derivative",
      "Computed locally for this drafting option as 90th percentile of inner spline derivative from exported count trajectories.",
      ""
    )
  )

line_summary_out <- line_summary %>%
  mutate(
    representative_line = representative_line,
    source_feature_values = "initial_subpanels/staged_feature_workflow/fg2_staged_feature_values.csv"
  )

write_csv_tbl(feature_values_out, file.path(INITIAL_DIR, "fg2_staged_feature_values.csv"))
write_csv_tbl(line_summary_out, file.path(INITIAL_DIR, "fg2_staged_feature_summary.csv"))

manifest <- option_paths %>%
  mutate(
    status = if_else(file.exists(png_path), "written", "missing"),
    source_script = file.path(DRAFT_ROOT, "scripts", "fg2_staged_feature_workflow_worker.R"),
    data_inputs = paste(
      count_summary_path,
      glucose_summary_path,
      feature_panel_path,
      signature_panel_path,
      sep = ";"
    ),
    feedback_source = feedback_path
  ) %>%
  bind_rows(refined_paths %>%
    mutate(
      status = if_else(file.exists(png_path), "written", "missing"),
      source_script = file.path(DRAFT_ROOT, "scripts", "fg2_staged_feature_workflow_worker.R"),
      data_inputs = paste(
        count_summary_path,
        glucose_summary_path,
        feature_panel_path,
        signature_panel_path,
        sep = ";"
      ),
      feedback_source = feedback_path
    )) %>%
  bind_rows(tibble(
    option_id = "refined_contact_sheet",
    png_path = REFINED_PATH,
    width_in = 7.0,
    height_in = 16.5,
    description = "Contact sheet combining the three staged workflow options for reviewer choice.",
    status = if_else(file.exists(REFINED_PATH), "written", "missing"),
    source_script = file.path(DRAFT_ROOT, "scripts", "fg2_staged_feature_workflow_worker.R"),
    data_inputs = paste(
      count_summary_path,
      glucose_summary_path,
      feature_panel_path,
      signature_panel_path,
      sep = ";"
    ),
    feedback_source = feedback_path
  ))
write_csv_tbl(manifest, file.path(INITIAL_DIR, "fg2_staged_feature_workflow_manifest.csv"))

source_status <- tibble(
  source_path = c(count_summary_path, glucose_summary_path, feature_panel_path, signature_panel_path, feedback_path),
  role = c(
    "exported live/total count trajectories for spline overlays",
    "exported glucose trajectories for yield drawdown display",
    "existing no-SUM per-well direct features for AUC/yield and metadata",
    "existing no-SUM line-level signatures, read for provenance continuity",
    "current direct user feedback"
  ),
  status = if_else(file.exists(source_path), "present", "missing")
)
write_csv_tbl(source_status, file.path(INITIAL_DIR, "fg2_staged_feature_workflow_sources.csv"))

notes <- c(
  "# FG2 staged feature-workflow coverage",
  "",
  "Worker: `scripts/fg2_staged_feature_workflow_worker.R`.",
  "",
  "Command: `scripts/agentRrunner.sh agent-dev/manuscript_redrafts/20260709T210405_v7_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/scripts/fg2_staged_feature_workflow_worker.R`.",
  "",
  "Scope: staged workflow construction outputs for Figure 2B/C. This worker regenerates local initial/refined PNGs only; final figure assembly, `report_manifest.csv`, `review_report.html`, and `scripts/make_fg2_v2_review_report.R` remain out of scope for this round.",
  "",
  "## Data sources",
  "",
  paste0("- `", count_summary_path, "`: exported live/dead/total count trajectory summaries used for local spline overlays."),
  paste0("- `", glucose_summary_path, "`: exported glucose trajectory summaries used to show glucose drawdown in total-cell yield."),
  paste0("- `", feature_panel_path, "`: existing no-SUM per-well direct-feature export used for live-cell AUC, total-cell yield, glucose-window metadata, and per-G0 values."),
  paste0("- `", signature_panel_path, "`: existing no-SUM signature export read for provenance continuity with recommended Figure 2A-like feature summaries."),
  "",
  "No major analysis, upstream refit, segmentation rerun, classifier retraining, or SLURM job was run. The new robust-derivative values are local drafting-option overlays computed from exported count trajectories.",
  "",
  "## Directive coverage",
  "",
  "- FG2-ROUND3-D01, remove the in-plot example-condition text label: addressed. The live spline example no longer includes the `MDA-MB-231 lower 3.5 mM G0` style label.",
  "- FG2-ROUND3-D02, keep live AUC shading and allow the derivative tangent: addressed. The live AUC shaded region and derivative tangent remain in Figure 2B.",
  "- FG2-ROUND3-D03, remove spline fits from other initial glucose conditions: addressed. The Figure 2B example shows only the selected live spline and no grey all-G0 spline overlays.",
  "- FG2-ROUND3-D04, put total cell yield in the same plot as the live spline: addressed. Figure 2B overlays log1p(total cells) on the same time axis and y scale as the live spline.",
  "- FG2-ROUND3-D05, represent max delta total with two horizontal dashed lines and a vertical double-ended arrow: addressed. Figure 2B uses dashed min/max total-cell lines plus a vertical arrow labeled `max delta`.",
  "- FG2-ROUND3-D06, break the per-G0 value split/regressions out of panel B so it can become panel C: addressed. The worker writes separate Figure 2B, Figure 2C, and paired B/C review PNGs under initial and refined outputs.",
  "- Current redraft, replace the misleading yield-per-glucose regression: addressed for Figure 2C. The yield facet now shows only net peak total-cell gain at G0 = 1 (`max(total_peak_to_glucose_end - total_initial_at_glucose_start, 0)`) and draws no regression line.",
  "- Current redraft, make the live-AUC construction match the canonical feature: addressed for Figure 2C. The live-AUC facet now omits G0 > 1 points and overlays the true `lm(live_auc_glucose_window ~ log1p(G0))` fit on a linear G0 axis.",
  "- FG2-STAGED-D04, aggregate values to a final feature summary for a given cell line: addressed without the rejected final-summary mini-plot. The per-G0 panels show the derivative low/high-G split, halo the points contributing to derivative medians, and overlay all-G0 regression lines for AUC and yield gradients.",
  "- FG2-ROUND3-D08, Figure S4 curved-regression concern: not edited in this worker because the explicit write scope for this task is Figure 2B/C staged workflow outputs only.",
  "",
  "## Outputs",
  "",
  paste0("- `", option_paths$png_path, "`"),
  paste0("- `", refined_paths$png_path, "`"),
  paste0("- `", REFINED_PATH, "`"),
  paste0("- `", file.path(INITIAL_DIR, "fg2_staged_feature_values.csv"), "`"),
  paste0("- `", file.path(INITIAL_DIR, "fg2_staged_feature_summary.csv"), "`"),
  paste0("- `", file.path(INITIAL_DIR, "fg2_staged_feature_workflow_manifest.csv"), "`"),
  paste0("- `", file.path(INITIAL_DIR, "fg2_staged_feature_workflow_sources.csv"), "`"),
  "",
  "## Visual QC",
  "",
  "- PNGs are 7 inches wide at 420 dpi and contain no figure-number headers or manuscript captions.",
  "- Option PNGs use lowercase patchwork tags where they are multi-panel review options; the promoted Figure 2 construction strip intentionally removes those tags for downstream assembly compatibility.",
  "- The promoted Figure 2 construction strip does not include the small final-summary mini-plot rejected in round-2 feedback.",
  "- Figure 2B/C are now available as separate refined PNGs, plus a paired review strip with lowercase b/c tags.",
  "- The contact sheet is review-facing and denser than a final figure by design; individual option PNGs remain separately inspectable.",
  "- Death-derived features are absent from this staged workflow.",
  "",
  "## Caveats and blockers",
  "",
  "- No blocker for drafting or promoting the requested current-feedback candidates.",
  "- The robust max derivative is a candidate construction rule for this feedback round; it is promoted here as a review candidate, not as a new upstream analysis/refit.",
  "- Yield values at G0 levels with insufficient glucose drawdown remain absent rather than imputed.",
  "- The existing HTML report is regenerated separately by `scripts/make_fg2_v2_review_report.R`."
)
writeLines(notes, NOTES_PATH, useBytes = TRUE)
writeLines(c(
  "# FG2 round-3 Figure 2B/C coverage",
  "",
  "Worker: `scripts/fg2_staged_feature_workflow_worker.R`.",
  "",
  "Command: `scripts/agentRrunner.sh agent-dev/manuscript_redrafts/20260709T210405_v7_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/scripts/fg2_staged_feature_workflow_worker.R`.",
  "",
  "Feedback source: `agent-dev/manuscript_redrafts/20260619_v6_redraft/handoffs/20260709T124521_figures_handoff.md`.",
  "",
  "## Changed construction outputs",
  "",
  "- Figure 2B is split into `initial_subpanels/staged_feature_workflow/fg2_figure2b_live_spline_total_yield.png` and `refined_subpanels/fg2_figure2b_live_spline_total_yield.png`.",
  "- Figure 2C is split into `initial_subpanels/staged_feature_workflow/fg2_figure2c_per_g0_split_regressions.png` and `refined_subpanels/fg2_figure2c_per_g0_split_regressions.png`.",
  "- A paired review strip is written to `initial_subpanels/staged_feature_workflow/fg2_figure2bc_construction_pair.png` and `refined_subpanels/fg2_figure2bc_construction_pair.png`.",
  "- The existing downstream-compatible promoted strip path is also regenerated: `initial_subpanels/staged_feature_workflow/fg2_staged_workflow_promoted_figure2_construction_strip.png`.",
  "",
  "## Directive coverage",
  "",
  "- Removed the internal example-condition label from the live spline panel.",
  "- Preserved live AUC shading and derivative tangent annotation.",
  "- Removed spline overlays from other initial glucose conditions.",
  "- Moved total-cell yield into the live spline plot as log1p(total cells) on the same time axis.",
  "- Replaced the slanted total-yield arrow with dashed min/max total-cell lines and a vertical double-ended `max delta` arrow.",
  "- Broke the per-G0 split/regression display out as Figure 2C.",
  "- Current feature-definition redraft: Figure 2C now omits G0 > 1 live-AUC points and draws the canonical `lm(live_auc_glucose_window ~ log1p(G0))` fit on a linear G0 axis; yield is shown only as fixed G0 = 1 net total-cell gain with no regression line.",
  "- Did not revise Figure S4 in this worker because this task's write scope is Figure 2B/C staged workflow outputs only.",
  "",
  "## Project-map decision",
  "",
  "- No `docs/project_map.txt` update needed: this is a local drafting-output revision, not a maintained workflow or canonical-path change."
), ROUND3_NOTES_PATH, useBytes = TRUE)

cat(sprintf("Wrote %d staged workflow option PNGs and contact sheet under %s\n", nrow(option_paths), DRAFT_ROOT))
