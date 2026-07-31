#!/usr/bin/env Rscript

DRAFT_ROOT <- file.path(
  "agent-dev", "manuscript_redrafts", "20260709T210405_v7_redraft",
  "figure_generation", "FG2_direct_feature_rebuild", "drafting_v2"
)
INITIAL_DIR <- file.path(DRAFT_ROOT, "initial_subpanels", "supplement_qc")
REFINED_DIR <- file.path(DRAFT_ROOT, "refined_subpanels", "supplement_qc")
WORKER_NOTES_DIR <- file.path(DRAFT_ROOT, "worker_notes")

dir.create(INITIAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REFINED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORKER_NOTES_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(scales)
  library(grid)
})

sum_fuse_line <- "SUM-159-fuse"
representative_line <- "MDA-MB-231"
png_dpi <- 420

ploidy_palette <- c(low = "#1B6CA8", high = "#C43C39")
feature_palette <- c(Growth = "#009E73", `Alive AUC` = "#7A5195", `Total-cell yield` = "#0072B2")
g0_levels <- c("0", "0.1", "0.25", "0.5", "1", "5", "25")

source_manifest_rows <- list()

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

record_source <- function(source_path, role) {
  source_manifest_rows[[length(source_manifest_rows) + 1L]] <<- tibble(
    source_path = source_path,
    role = role,
    status = if_else(file.exists(source_path), "present", "missing")
  )
  invisible(NULL)
}

save_plot_png <- function(plot, path, width, height, dpi = png_dpi) {
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

save_review_png <- function(plot, dir, filename, width, height, dpi = png_dpi) {
  png_path <- file.path(dir, filename)
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(png_path)
}

save_subpanel <- function(plot, figure, panel, filename, width, height, dpi = png_dpi) {
  png_path <- file.path(INITIAL_DIR, filename)
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  tibble(
    figure = figure,
    panel = panel,
    subpanel_png = png_path,
    width_px = as.integer(round(width * dpi)),
    height_px = as.integer(round(height * dpi)),
    width_in = width,
    height_in = height
  )
}

wp2_theme <- function(base_size = 7.2) {
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
      plot.title = element_text(face = "bold", size = rel(1.02)),
      plot.subtitle = element_text(size = rel(0.88), color = "grey22"),
      plot.caption = element_text(size = rel(0.76), color = "grey35")
    )
}

wrap_label <- function(x, width = 18) {
  vapply(as.character(x), function(one) paste(strwrap(one, width = width), collapse = "\n"), character(1))
}

format_ploidy_abs <- function(x) {
  out <- ifelse(abs(x - round(x)) < 1e-8, sprintf("%.0fN", x), sprintf("%.1fN", x))
  out[!is.finite(x)] <- NA_character_
  out
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
    out <- out %>% mutate(G0_display = factor(as.character(G0), levels = g0_levels))
  }
  out
}

short_number <- function() {
  label_number(scale_cut = cut_short_scale())
}

feature_meta <- tibble(
  feature = c(
    "growth_lowG_median",
    "growth_highG_median",
    "yield_alive_auc_intercept",
    "yield_alive_auc_slope",
    "peak_total_yield_1mM"
  ),
  display_label = c(
    "Low-G growth",
    "High-G exp. growth",
    "Alive AUC intercept",
    "Alive AUC gradient",
    "1 mM net total yield"
  ),
  feature_class = c("Growth", "Growth", "Alive AUC", "Alive AUC", "Total-cell yield")
)

add_fixed_yield_signature <- function(signature_panel, feature_panel) {
  yield_1mm <- feature_panel %>%
    filter(abs(G0 - 1) < 1e-8) %>%
    mutate(peak_total_yield_1mM = pmax(total_net_gain_to_glucose_end, 0)) %>%
    group_by(cellLine, line_id, ploidy_metric) %>%
    summarise(peak_total_yield_1mM = median(peak_total_yield_1mM, na.rm = TRUE), .groups = "drop") %>%
    mutate(peak_total_yield_1mM = if_else(is.nan(peak_total_yield_1mM), NA_real_, peak_total_yield_1mM))

  signature_panel %>%
    select(-any_of(c("peak_total_yield_intercept", "peak_total_yield_slope", "peak_total_yield_1mM"))) %>%
    left_join(yield_1mm, by = c("cellLine", "line_id", "ploidy_metric"))
}

smooth_one_well <- function(df, value_col, spar = 0.55, n_eval = 180L) {
  df <- df %>% arrange(hours)
  ok <- is.finite(df$hours) & is.finite(df[[value_col]])
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 2L) {
    return(list(smooth = tibble(), point = tibble()))
  }

  x <- df$hours
  y_plot <- log1p(pmax(df[[value_col]], 0))
  eval_hours <- seq(min(x), max(x), length.out = n_eval)

  if (length(unique(x)) >= 4L) {
    fit <- try(stats::smooth.spline(x = x, y = y_plot, spar = spar), silent = TRUE)
  } else {
    fit <- structure(list(), class = "try-error")
  }

  if (inherits(fit, "try-error")) {
    smooth_y <- approx(x = x, y = y_plot, xout = eval_hours, ties = "ordered", rule = 2)$y
    deriv <- c(diff(smooth_y) / diff(eval_hours), NA_real_)
  } else {
    smooth_y <- predict(fit, x = eval_hours)$y
    deriv <- predict(fit, x = eval_hours, deriv = 1)$y
  }

  idx <- if (all(!is.finite(deriv))) which.max(smooth_y) else which.max(deriv)
  meta <- df[1, c("well_idx", "cellLine", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id", "has_starvation"), drop = FALSE]

  list(
    smooth = bind_cols(
      meta[rep(1L, length(eval_hours)), ],
      tibble(hours = eval_hours, smooth_value = smooth_y, derivative = deriv)
    ),
    point = bind_cols(
      meta,
      tibble(hours = eval_hours[idx], smooth_value = smooth_y[idx], max_rate = deriv[idx])
    )
  )
}

build_smooth_bundle <- function(count_df, value_col) {
  pieces <- count_df %>%
    split(.$well_idx) %>%
    lapply(smooth_one_well, value_col = value_col)
  list(
    smooth = bind_rows(lapply(pieces, `[[`, "smooth")),
    extrema = bind_rows(lapply(pieces, `[[`, "point"))
  )
}

make_lm_grid <- function(df, x_col, y_col, n = 80L) {
  ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 2L || length(unique(df[[x_col]])) < 2L) {
    return(tibble())
  }
  fit <- lm(stats::as.formula(paste(y_col, "~", x_col)), data = df)
  grid <- tibble(x_value = seq(min(df[[x_col]]), max(df[[x_col]]), length.out = n))
  names(grid) <- x_col
  bind_cols(grid, tibble(pred = predict(fit, newdata = grid)))
}

make_log_glucose_lm_curve <- function(df, y_col = "value", n = 80L) {
  ok <- is.finite(df$G0) & is.finite(df[[y_col]])
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 2L || length(unique(df$G0)) < 2L) {
    return(tibble())
  }
  fit <- lm(stats::as.formula(paste(y_col, "~ log1p(G0)")), data = df)
  grid <- tibble(G0 = seq(min(df$G0), max(df$G0), length.out = n))
  grid$log1p_G0 <- log1p(grid$G0)
  grid$pred <- predict(fit, newdata = grid)
  grid
}

make_low_g_growth_panel <- function(count_summary, feature_panel) {
  traj <- count_summary %>%
    filter(cellLine == representative_line, G0 <= 0.25) %>%
    add_ploidy_display() %>%
    mutate(log_live = log1p(pmax(live_cells, 0)))

  bundle <- build_smooth_bundle(traj, "live_cells")
  smooth <- bundle$smooth %>% add_ploidy_display()
  extrema <- bundle$extrema %>% add_ploidy_display()

  rate_summary <- feature_panel %>%
    filter(cellLine == representative_line, G0 <= 0.25, is.finite(max_growth_rate)) %>%
    add_ploidy_display()

  log_plot <- ggplot(traj, aes(hours, log_live, color = ploidy_state)) +
    geom_line(aes(group = well_idx), linewidth = 0.24, alpha = 0.38) +
    geom_line(data = smooth, aes(y = smooth_value, group = well_idx), linewidth = 0.48, alpha = 0.9) +
    geom_point(size = 0.35, alpha = 0.45) +
    facet_wrap(~paste0("G0=", G0_display), nrow = 1, ncol = 3) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(x = NULL, y = "log1p(live cells)") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", strip.text = element_text(size = 5.5), panel.spacing = unit(0.2, "lines"))

  deriv_plot <- ggplot(smooth, aes(hours, derivative, color = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey58") +
    geom_line(aes(group = well_idx), linewidth = 0.35, alpha = 0.58) +
    geom_point(data = extrema, aes(y = max_rate), shape = 24, fill = "white", stroke = 0.28, size = 1.1) +
    facet_wrap(~paste0("G0=", G0_display), nrow = 1, ncol = 3) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(x = "hours", y = "d log1p(live) / dt") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", strip.text = element_blank(), panel.spacing = unit(0.2, "lines"))

  rate_plot <- ggplot(rate_summary, aes(G0_display, max_growth_rate, color = ploidy_state, group = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey58") +
    geom_point(position = position_dodge(width = 0.38), size = 1.25) +
    stat_summary(fun = median, geom = "crossbar", width = 0.42, linewidth = 0.22, position = position_dodge(width = 0.38)) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(x = "G0", y = "max derivative") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))

  wrap_plots(log_plot, deriv_plot, rate_plot, ncol = 1, heights = c(0.88, 0.82, 0.62)) &
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
}

make_high_g_growth_panel <- function(count_summary, feature_panel) {
  traj <- count_summary %>%
    filter(cellLine == representative_line, G0 %in% c(5, 25)) %>%
    add_ploidy_display() %>%
    mutate(log_live = log1p(pmax(live_cells, 0)))

  fit_rows <- traj %>%
    group_by(well_idx) %>%
    group_modify(function(df, key) {
      feat <- feature_panel %>% filter(well_idx == key$well_idx[[1]])
      end_time <- feat$exp_growth_fit_end_time_70[[1]]
      if (!is.finite(end_time)) return(tibble())
      fit_df <- df %>% filter(hours <= end_time)
      if (nrow(fit_df) < 3L || length(unique(fit_df$hours)) < 2L) return(tibble())
      fit <- lm(log1p(pmax(live_cells, 0)) ~ hours, data = fit_df)
      pred_hours <- seq(min(fit_df$hours), max(fit_df$hours), length.out = 60L)
      meta <- df %>% slice(1) %>% select(cellLine, line_id, ploidy_metric, ploidy_abs, G0, ploidy_state, ploidy_display, G0_display)
      bind_cols(meta[rep(1L, length(pred_hours)), ], tibble(hours = pred_hours, fit_log_live = predict(fit, newdata = tibble(hours = pred_hours))))
    }) %>%
    ungroup()

  slopes <- feature_panel %>%
    filter(cellLine == representative_line, G0 %in% c(5, 25), is.finite(exp_growth_rate_70)) %>%
    add_ploidy_display()

  traj_plot <- ggplot(traj, aes(hours, log_live, color = ploidy_state)) +
    geom_line(aes(group = well_idx), linewidth = 0.28, alpha = 0.45) +
    geom_point(size = 0.36, alpha = 0.5) +
    geom_line(data = fit_rows, aes(y = fit_log_live, group = well_idx), color = "#D55E00", linewidth = 0.58, alpha = 0.9) +
    facet_wrap(~paste0("G0=", G0_display), nrow = 1) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(x = "hours", y = "log1p(live cells)") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", strip.text = element_text(size = 5.5), panel.spacing = unit(0.25, "lines"))

  slope_plot <- ggplot(slopes, aes(G0_display, exp_growth_rate_70, color = ploidy_state, group = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey58") +
    geom_point(position = position_dodge(width = 0.38), size = 1.35) +
    stat_summary(fun = median, geom = "crossbar", width = 0.42, linewidth = 0.22, position = position_dodge(width = 0.38)) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(x = "G0", y = "fitted slope") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none")

  wrap_plots(traj_plot, slope_plot, ncol = 1, heights = c(0.95, 0.58)) &
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
}

make_auc_panel <- function(count_summary, feature_panel) {
  reg_data <- feature_panel %>%
    filter(cellLine == representative_line, G0 <= 1, is.finite(live_auc_glucose_window)) %>%
    add_ploidy_display() %>%
    mutate(log1p_G0 = log1p(G0))

  example_row <- reg_data %>%
    filter(G0 == 1) %>%
    arrange(ploidy_metric, well_idx) %>%
    slice(1)
  example_traj <- count_summary %>%
    filter(well_idx == example_row$well_idx[[1]], hours >= example_row$glucose_start_time[[1]], hours <= example_row$glucose_end_time[[1]]) %>%
    arrange(hours)
  auc_poly <- bind_rows(
    example_traj %>% slice(1) %>% mutate(live_cells = 0),
    example_traj,
    example_traj %>% slice(n()) %>% mutate(live_cells = 0)
  )

  fit_df <- reg_data %>%
    group_by(ploidy_state, ploidy_display) %>%
    group_modify(~make_lm_grid(.x, "log1p_G0", "live_auc_glucose_window")) %>%
    ungroup()

  coef_df <- reg_data %>%
    group_by(ploidy_state, ploidy_display) %>%
    summarise(
      intercept = coef(lm(live_auc_glucose_window ~ log1p_G0))[["(Intercept)"]],
      slope = coef(lm(live_auc_glucose_window ~ log1p_G0))[["log1p_G0"]],
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(ploidy_display, ": int ", comma(round(intercept)), "; slope ", comma(round(slope))),
      x = min(reg_data$log1p_G0, na.rm = TRUE) + 0.02,
      y = max(reg_data$live_auc_glucose_window, na.rm = TRUE) - row_number() * 0.13 * diff(range(reg_data$live_auc_glucose_window, na.rm = TRUE))
    )

  auc_plot <- ggplot(example_traj, aes(hours, live_cells)) +
    geom_polygon(data = auc_poly, aes(hours, live_cells), fill = "#80B1D3", alpha = 0.34, color = NA) +
    geom_line(color = "#1B6CA8", linewidth = 0.42) +
    geom_point(size = 0.55, color = "#1B6CA8") +
    labs(x = "hours", y = "live cells") +
    wp2_theme(base_size = 6.2) +
    theme(plot.title = element_text(size = 6.4))

  reg_plot <- ggplot(reg_data, aes(log1p_G0, live_auc_glucose_window, color = ploidy_state)) +
    geom_point(size = 1.05) +
    geom_line(data = fit_df, aes(y = pred, group = ploidy_state), linewidth = 0.5) +
    geom_text(data = coef_df, aes(x = x, y = y, label = label), hjust = 0, size = 1.8, show.legend = FALSE) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_y_continuous(labels = short_number()) +
    labs(x = "log1p(G0)", y = "live-cell AUC") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", plot.title = element_text(size = 6.4))

  wrap_plots(auc_plot, reg_plot, nrow = 1, widths = c(0.72, 1.25)) &
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
}

make_yield_panel <- function(feature_panel) {
  yield_data <- feature_panel %>%
    mutate(peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)) %>%
    filter(cellLine == representative_line, abs(G0 - 1) < 1e-8, is.finite(peak_total_yield_net)) %>%
    add_ploidy_display()

  example_row <- yield_data %>%
    arrange(ploidy_metric, well_idx) %>%
    slice(1)

  component_df <- tibble(
    component = factor("total-cell gain", levels = "total-cell gain"),
    start = example_row$total_initial_at_glucose_start,
    end = example_row$total_peak_to_glucose_end,
    value = example_row$peak_total_yield_net,
    label = paste0("gain ", comma(round(example_row$peak_total_yield_net)))
  )

  component_plot <- ggplot(component_df) +
    geom_segment(
      aes(x = start, xend = end, y = component, yend = component),
      linewidth = 0.7,
      color = "grey42",
      arrow = grid::arrow(length = unit(0.055, "in"))
    ) +
    geom_point(aes(start, component), size = 1.15, color = "#1B6CA8") +
    geom_point(aes(end, component), size = 1.15, color = "#C43C39") +
    facet_wrap(~component, scales = "free_x", ncol = 1) +
    scale_x_continuous(labels = short_number(), expand = expansion(mult = c(0.08, 0.4))) +
    labs(x = NULL, y = NULL) +
    wp2_theme(base_size = 6.2) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 5.4),
      plot.title = element_text(size = 6.4)
    )

  yield_plot <- ggplot(yield_data, aes(ploidy_state, peak_total_yield_net, color = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey58") +
    geom_point(size = 1.35) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_y_continuous(labels = short_number()) +
    labs(x = NULL, y = "1 mM net total yield") +
    wp2_theme(base_size = 6.2) +
    theme(legend.position = "none", plot.title = element_text(size = 6.4))

  wrap_plots(component_plot, yield_plot, nrow = 1, widths = c(0.72, 1.25)) &
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
}

make_raw_feature_pairs <- function(signature_panel, include_sum159_fuse, feature_panel = NULL) {
  line_levels <- signature_panel %>%
    arrange(line_id, cellLine) %>%
    distinct(cellLine) %>%
    pull(cellLine)

  x <- signature_panel %>%
    add_ploidy_display() %>%
    pivot_longer(all_of(feature_meta$feature), names_to = "feature", values_to = "raw_value") %>%
    left_join(feature_meta, by = "feature") %>%
    mutate(
      cellLine = factor(cellLine, levels = rev(line_levels)),
      display_label = factor(display_label, levels = feature_meta$display_label),
      feature_class = factor(feature_class, levels = names(feature_palette))
    )

  if (!is.null(feature_panel)) {
    auc0 <- feature_panel %>%
      filter(abs(G0) < 1e-8, is.finite(live_auc_glucose_window)) %>%
      add_ploidy_display() %>%
      transmute(
        cellLine = factor(cellLine, levels = rev(line_levels)),
        ploidy_state,
        feature = "live_auc_0mM",
        raw_value = live_auc_glucose_window,
        display_label = factor("0 mM live-cell AUC", levels = c(as.character(feature_meta$display_label), "0 mM live-cell AUC")),
        feature_class = factor("Alive AUC", levels = names(feature_palette))
      )
    x <- bind_rows(
      x %>% mutate(display_label = factor(as.character(display_label), levels = c(as.character(feature_meta$display_label), "0 mM live-cell AUC"))),
      auc0
    )
  }

  pair_df <- x %>%
    group_by(cellLine, feature, display_label, feature_class) %>%
    summarise(
      low_value = raw_value[ploidy_state == "low"][1],
      high_value = raw_value[ploidy_state == "high"][1],
      .groups = "drop"
    )

  p <- ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey68") +
    geom_segment(
      data = pair_df,
      aes(x = low_value, xend = high_value, y = cellLine, yend = cellLine),
      linewidth = 0.42,
      color = "grey55"
    ) +
    geom_point(
      data = x,
      aes(raw_value, cellLine, fill = ploidy_state, shape = ploidy_state),
      color = "grey18",
      size = 1.62,
      stroke = 0.22
    ) +
    facet_wrap(~display_label, scales = "free_x", ncol = 3) +
    scale_fill_manual(values = ploidy_palette, name = "ploidy state") +
    scale_shape_manual(values = c(low = 21, high = 24), name = "ploidy state") +
    scale_x_continuous(labels = short_number()) +
    labs(x = "raw feature value", y = NULL) +
    wp2_theme(base_size = 6.8) +
    theme(
      strip.text = element_text(size = 6),
      axis.text.x = element_text(size = 5.5),
      panel.spacing = unit(0.38, "lines")
    )

  if (include_sum159_fuse) {
    p <- p +
      geom_point(
        data = x %>% filter(cellLine == sum_fuse_line),
        aes(raw_value, cellLine, fill = ploidy_state, shape = ploidy_state),
        inherit.aes = FALSE,
        color = "black",
        stroke = 0.55,
        size = 2.35
      )
  }

  p
}

make_diagnostic_derivation_figure <- function(feature_panel) {
  x <- feature_panel %>% add_ploidy_display() %>% mutate(log1p_G0 = log1p(G0))
  line_levels <- x %>% arrange(line_id, cellLine) %>% distinct(cellLine) %>% pull(cellLine)

  diagnostics <- bind_rows(
    x %>%
      filter(G0 <= 0.25, is.finite(max_growth_rate)) %>%
      transmute(cellLine, line_id, ploidy_state, ploidy_display, G0, log1p_G0, value = max_growth_rate, feature_group = "Low-G growth\nmax derivative", fit_type = "points"),
    x %>%
      filter(G0 %in% c(5, 25), is.finite(exp_growth_rate_70)) %>%
      transmute(cellLine, line_id, ploidy_state, ploidy_display, G0, log1p_G0, value = exp_growth_rate_70, feature_group = "High-G growth\nexp. slope", fit_type = "points"),
    x %>%
      filter(G0 <= 1, is.finite(live_auc_glucose_window)) %>%
      transmute(cellLine, line_id, ploidy_state, ploidy_display, G0, log1p_G0, value = live_auc_glucose_window / 1000, feature_group = "Alive AUC / 1k\nintercept + slope", fit_type = "lm"),
    x %>%
      mutate(peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)) %>%
      filter(abs(G0 - 1) < 1e-8, is.finite(peak_total_yield_net)) %>%
      transmute(cellLine, line_id, ploidy_state, ploidy_display, G0, log1p_G0, value = peak_total_yield_net / 1000, feature_group = "1 mM net yield / 1k", fit_type = "points"),
    tibble()
  ) %>%
    mutate(
      cellLine = factor(cellLine, levels = line_levels),
      feature_group = factor(feature_group, levels = c(
        "Low-G growth\nmax derivative",
        "High-G growth\nexp. slope",
        "Alive AUC / 1k\nintercept + slope",
        "1 mM net yield / 1k"
      )),
      panel_label = factor(
        paste(feature_group, cellLine, sep = "\n"),
        levels = as.vector(t(outer(levels(feature_group), line_levels, paste, sep = "\n")))
      )
    )

  fit_lines <- diagnostics %>%
    filter(fit_type == "lm") %>%
    group_by(feature_group, cellLine, ploidy_state) %>%
    group_modify(function(df, key) {
      make_log_glucose_lm_curve(df, n = 60L)
    }) %>%
    ungroup() %>%
    mutate(
      panel_label = factor(
        paste(feature_group, cellLine, sep = "\n"),
        levels = levels(diagnostics$panel_label)
      )
    )

  ggplot(diagnostics, aes(G0, value, color = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.18, color = "grey70") +
    geom_line(
      data = diagnostics %>% filter(fit_type == "points", feature_group != "1 mM net yield / 1k"),
      aes(group = ploidy_state),
      linewidth = 0.22,
      alpha = 0.42
    ) +
    geom_line(
      data = fit_lines,
      aes(G0, pred, color = ploidy_state, group = ploidy_state),
      linewidth = 0.35,
      alpha = 0.8,
      inherit.aes = FALSE
    ) +
    geom_point(size = 0.72, alpha = 0.86) +
    facet_wrap(~panel_label, ncol = length(line_levels), scales = "free") +
    scale_color_manual(values = ploidy_palette, name = "ploidy state") +
    scale_x_continuous(breaks = c(0, 0.25, 1, 5, 25), labels = c("0", ".25", "1", "5", "25")) +
    labs(x = "starting glucose G0", y = "feature-construction value (row-specific units)") +
    wp2_theme(base_size = 6.0) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 4.7),
      axis.text.x = element_text(size = 4.8, angle = 0),
      axis.text.y = element_text(size = 4.8),
      panel.spacing.x = unit(0.18, "lines"),
      panel.spacing.y = unit(0.22, "lines")
    )
}

compute_robust_derivative_one <- function(df, spar = 0.62, n_eval = 260L, edge_frac = 0.06) {
  df <- df %>% arrange(hours)
  ok <- is.finite(df$hours) & is.finite(df$live_cells)
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 2L) {
    return(tibble())
  }

  eval_hours <- seq(min(df$hours, na.rm = TRUE), max(df$hours, na.rm = TRUE), length.out = n_eval)
  y <- log1p(pmax(df$live_cells, 0))
  fit <- if (length(unique(df$hours)) >= 4L) {
    try(stats::smooth.spline(df$hours, y, spar = spar), silent = TRUE)
  } else {
    structure(list(), class = "try-error")
  }
  derivative <- if (inherits(fit, "try-error")) {
    smooth_y <- approx(df$hours, y, xout = eval_hours, rule = 2)$y
    c(diff(smooth_y) / diff(eval_hours), NA_real_)
  } else {
    predict(fit, x = eval_hours, deriv = 1)$y
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

  meta_cols <- intersect(
    c("well_idx", "cellLine", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id", "has_starvation"),
    names(df)
  )
  bind_cols(
    df[1, meta_cols, drop = FALSE],
    tibble(
      robust_max_derivative = robust_max,
      robust_rule = "90th percentile of inner smooth-spline derivative"
    )
  )
}

build_robust_derivative_panel <- function(count_summary, feature_panel) {
  count_summary %>%
    semi_join(feature_panel %>% distinct(well_idx), by = "well_idx") %>%
    split(.$well_idx) %>%
    lapply(compute_robust_derivative_one) %>%
    bind_rows()
}

make_s4_round2_all_cell_line_construction <- function(count_summary, feature_panel) {
  robust_panel <- build_robust_derivative_panel(count_summary, feature_panel)
  line_levels <- feature_panel %>% arrange(line_id, cellLine) %>% distinct(cellLine) %>% pull(cellLine)

	  x <- feature_panel %>%
	    select(
	      well_idx, cellLine, line_id, ploidy_metric, ploidy_abs, G0, exp_id, has_starvation,
	      live_auc_glucose_window, total_net_gain_to_glucose_end
	    ) %>%
	    left_join(robust_panel %>% select(well_idx, robust_max_derivative), by = "well_idx") %>%
	    mutate(peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)) %>%
	    pivot_longer(
	      cols = c(robust_max_derivative, live_auc_glucose_window, peak_total_yield_net),
	      names_to = "feature_id",
	      values_to = "value"
	    ) %>%
	    filter(case_when(
	      feature_id == "live_auc_glucose_window" ~ G0 <= 1,
	      feature_id == "peak_total_yield_net" ~ abs(G0 - 1) < 1e-8,
	      TRUE ~ TRUE
	    )) %>%
	    filter(is.finite(value)) %>%
	    add_ploidy_display() %>%
	    mutate(
	      cellLine = factor(cellLine, levels = line_levels),
	      feature_label = recode(
	        feature_id,
	        robust_max_derivative = "Robust max derivative",
	        live_auc_glucose_window = "Live AUC\nintercept/gradient",
	        peak_total_yield_net = "1 mM net\ntotal yield"
	      ),
	      feature_label = factor(
	        feature_label,
	        levels = c(
	          "Robust max derivative",
	          "Live AUC\nintercept/gradient",
	          "1 mM net\ntotal yield"
	        )
	      ),
	      panel_label = factor(
	        paste(feature_label, cellLine, sep = "\n"),
	        levels = as.vector(t(outer(levels(feature_label), line_levels, paste, sep = "\n")))
	      ),
	      derivative_summary_point = feature_id == "robust_max_derivative" & (G0 <= 0.25 | G0 >= 1),
	      fixed_yield_point = feature_id == "peak_total_yield_net"
	    )

	  fit_lines <- x %>%
	    filter(feature_id == "live_auc_glucose_window") %>%
	    group_by(feature_label, cellLine, ploidy_state) %>%
	    group_modify(function(df, key) {
	      make_log_glucose_lm_curve(df, n = 70L)
	    }) %>%
	    ungroup() %>%
	    mutate(
	      panel_label = factor(
	        paste(feature_label, cellLine, sep = "\n"),
	        levels = levels(x$panel_label)
	      )
	    )

  split_df <- x %>%
	    filter(feature_id == "robust_max_derivative") %>%
	    distinct(feature_label, cellLine) %>%
	    mutate(xintercept = mean(c(0.25, 1)))

	  ggplot(x, aes(G0, value, color = ploidy_state, group = ploidy_state)) +
	    geom_hline(yintercept = 0, linewidth = 0.16, color = "grey74") +
	    geom_line(
	      data = x %>% filter(feature_id == "robust_max_derivative"),
	      linewidth = 0.22,
	      alpha = 0.42
	    ) +
	    geom_line(
	      data = fit_lines,
	      aes(G0, pred, color = ploidy_state, group = ploidy_state),
	      linewidth = 0.34,
	      alpha = 0.82,
	      inherit.aes = FALSE
	    ) +
    geom_vline(
      data = split_df,
      aes(xintercept = xintercept),
      linetype = "dashed",
      linewidth = 0.18,
      color = "grey35"
    ) +
	    geom_point(size = 0.74, alpha = 0.88) +
	    geom_point(
	      data = x %>% filter(feature_id == "live_auc_glucose_window", abs(G0) < 1e-8),
	      aes(G0, value),
	      shape = 21,
	      fill = "white",
	      color = "black",
	      stroke = 0.30,
	      size = 1.36
	    ) +
	    geom_point(
	      data = x %>% filter(derivative_summary_point | fixed_yield_point),
	      aes(G0, value),
	      shape = 21,
	      fill = NA,
	      color = "black",
      stroke = 0.28,
	      size = 1.32
	    ) +
	    facet_wrap(~panel_label, ncol = length(line_levels), scales = "free") +
	    scale_color_manual(values = ploidy_palette, name = "ploidy state") +
	    scale_x_continuous(
	      breaks = c(0, 0.25, 1, 5, 25),
	      labels = c("0", ".25", "1", "5", "25"),
	      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_y_continuous(labels = short_number()) +
    labs(
      x = "starting glucose G0",
      y = "feature value"
    ) +
    wp2_theme(base_size = 5.7) +
    theme(
	      legend.position = "bottom",
	      strip.text = element_text(size = 4.7),
      axis.text.x = element_text(size = 4.7),
      axis.text.y = element_text(size = 4.5),
      panel.spacing.x = unit(0.16, "lines"),
      panel.spacing.y = unit(0.20, "lines")
    )
}

make_range_inclusion_panel <- function(feature_panel) {
  g0_values <- sort(unique(feature_panel$G0[is.finite(feature_panel$G0)]))
  x <- tibble(G0 = g0_values) %>%
    tidyr::crossing(
      feature_group = factor(
        c(
          "Low-G growth\nmax derivative",
          "High-G growth\nexp. slope",
          "Alive AUC\nregression",
          "1 mM net\ntotal yield"
        ),
        levels = c(
          "Low-G growth\nmax derivative",
          "High-G growth\nexp. slope",
          "Alive AUC\nregression",
          "1 mM net\ntotal yield"
        )
      )
    ) %>%
    mutate(
      included = case_when(
        feature_group == "Low-G growth\nmax derivative" ~ G0 <= 0.25,
        feature_group == "High-G growth\nexp. slope" ~ G0 %in% c(5, 25),
        feature_group == "Alive AUC\nregression" ~ G0 <= 1,
        feature_group == "1 mM net\ntotal yield" ~ abs(G0 - 1) < 1e-8,
        TRUE ~ FALSE
      ),
      G0_display = factor(as.character(G0), levels = g0_levels),
      status = if_else(included, "included", "not used")
    )

  ggplot(x, aes(G0_display, feature_group, fill = status)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = if_else(included, "in", "")), size = 2.1, color = "grey18") +
    scale_fill_manual(values = c(included = "#D9EAD3", `not used` = "#F2F2F2"), name = "range") +
    labs(x = "starting glucose G0", y = NULL) +
    wp2_theme(base_size = 6.4) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 5.4),
      axis.text.y = element_text(size = 5.8),
      panel.grid = element_blank()
    )
}

make_spline_qc_panel <- function(count_summary, feature_panel) {
  traj <- count_summary %>%
    filter(cellLine == representative_line, G0 <= 0.25) %>%
    add_ploidy_display() %>%
    mutate(log_live = log1p(pmax(live_cells, 0)))

  bundle <- build_smooth_bundle(traj, "live_cells")
  smooth <- bundle$smooth %>% add_ploidy_display()
  extrema <- bundle$extrema %>% add_ploidy_display()

  ggplot(smooth, aes(hours, derivative, color = ploidy_state)) +
    geom_hline(yintercept = 0, linewidth = 0.22, color = "grey58") +
    geom_line(aes(group = well_idx), linewidth = 0.32, alpha = 0.56) +
    geom_point(data = extrema, aes(y = max_rate), shape = 24, fill = "white", stroke = 0.28, size = 1.1) +
    facet_wrap(~paste0("G0=", G0_display), nrow = 1, ncol = 3) +
    scale_color_manual(values = ploidy_palette, name = "ploidy state") +
    scale_x_continuous(breaks = c(0, 60, 120)) +
    labs(x = "hours", y = "spline derivative") +
    wp2_theme(base_size = 6.4) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 5.7),
      panel.spacing = unit(0.2, "lines")
    )
}

make_s4_construction_qc_figure <- function(count_summary, feature_panel) {
  diagnostic <- make_diagnostic_derivation_figure(feature_panel)
  spline_qc <- make_spline_qc_panel(count_summary, feature_panel)
  range_qc <- make_range_inclusion_panel(feature_panel)

  wrap_elements(full = diagnostic) /
    (wrap_elements(full = spline_qc) | wrap_elements(full = range_qc)) +
    plot_layout(heights = c(1.65, 1.0), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")
}

wp1_dir <- file.path("data", "report_exports", "wp1_core_starvation")
wp4_dir <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v6_redraft",
  "stage_outputs", "analysis", "recomputed_wp4_sum159_fuse_exclusion"
)

count_summary_path <- file.path(wp1_dir, "count_trajectory_summary.csv")
feature_panel_all_path <- file.path(wp4_dir, "wp4_feature_panel_all_lines.csv")
signature_all_path <- file.path(wp4_dir, "wp4_signature_panel_all_lines.csv")
user_feedback_path <- file.path("figures", "manuscript_draft_v3", "wp2_drafting", "user-feedback.txt")
fg2_feedback_path <- file.path(
  "agent-dev", "manuscript_redrafts", "20260709T210405_v7_redraft",
  "handoffs", "20260709T124521_figures_handoff.md"
)

for (path in c(count_summary_path, feature_panel_all_path, signature_all_path, user_feedback_path, fg2_feedback_path)) {
  record_source(path, "FG2 v5 supplement-QC input")
}

count_summary <- read_csv_tbl(count_summary_path)
feature_panel_all <- read_csv_tbl(feature_panel_all_path)
signature_all <- add_fixed_yield_signature(read_csv_tbl(signature_all_path), feature_panel_all)

p_raw_all <- make_raw_feature_pairs(signature_all, include_sum159_fuse = TRUE, feature_panel = feature_panel_all) &
  theme(legend.position = "bottom")

p_s4_diagnostic <- make_diagnostic_derivation_figure(feature_panel_all)
p_s4_round2_all <- make_s4_round2_all_cell_line_construction(count_summary, feature_panel_all)
p_s4_spline <- make_spline_qc_panel(count_summary, feature_panel_all)
p_s4_range <- make_range_inclusion_panel(feature_panel_all)
p_s4_assembled <- make_s4_construction_qc_figure(count_summary, feature_panel_all)

outputs <- tibble(
  figure = c("Figure S3", "Figure S4", "Figure S4", "Figure S4", "Figure S4", "Figure S4"),
  panel = c("standalone", "round2", "a", "b", "c", "assembled"),
  filename = c(
    "s3_standalone_all_line_sensitivity.png",
    "s4_round2_all_cell_line_construction.png",
    "s4_regression_feature_qc_no_death.png",
    "s4_spline_derivative_example.png",
    "s4_glucose_range_inclusion.png",
    "s4_construction_qc_support.png"
  ),
  width = c(7.0, 7.0, 7.0, 4.3, 3.0, 7.0),
  height = c(4.25, 6.2, 5.75, 2.65, 2.65, 8.5),
  directive_ids = c(
    "FG2-D04;FG2-R2-D10",
    "FG2-D05;FG2-D06;FG2-D18;FG2-R2-D11",
    "FG2-D05;FG2-D06",
    "FG2-D06",
    "FG2-D06",
    "FG2-D05;FG2-D06"
  )
)

plots <- list(p_raw_all, p_s4_round2_all, p_s4_diagnostic, p_s4_spline, p_s4_range, p_s4_assembled)

purrr::walk2(plots, seq_along(plots), function(plot, i) {
  spec <- outputs[i, ]
  save_review_png(plot, INITIAL_DIR, spec$filename, spec$width, spec$height)
  save_review_png(plot, REFINED_DIR, spec$filename, spec$width, spec$height)
})

write_csv_tbl(
  outputs %>%
    mutate(
      initial_png = file.path(INITIAL_DIR, filename),
      refined_png = file.path(REFINED_DIR, filename),
      active_script = file.path(DRAFT_ROOT, "scripts", "fg2_supplement_qc_worker.R"),
      data_inputs = case_when(
        figure == "Figure S3" ~ signature_all_path,
        TRUE ~ paste(count_summary_path, feature_panel_all_path, sep = ";")
      ),
      status = if_else(file.exists(initial_png) & file.exists(refined_png), "written", "missing")
    ),
  file.path(INITIAL_DIR, "supplement_qc_outputs.csv")
)

write_csv_tbl(bind_rows(source_manifest_rows), file.path(INITIAL_DIR, "supplement_qc_source_manifest.csv"))

round3_notes <- c(
  "# FG2 round-3 Figure S4 regression-display coverage",
  "",
  "Worker: `scripts/fg2_supplement_qc_worker.R`.",
  "",
  "Command: `scripts/agentRrunner.sh agent-dev/manuscript_redrafts/20260709T210405_v7_redraft/figure_generation/FG2_direct_feature_rebuild/drafting_v2/scripts/fg2_supplement_qc_worker.R`.",
  "",
  "Raw user-feedback source: `agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/user-feedback-20260625145000.txt`.",
  "",
  "## Directive coverage",
  "",
  "- Current redraft, replace total-yield regression with fixed 1 mM net yield: addressed. The active S4 all-cell-line construction view now shows `max(total_peak_to_glucose_end - total_initial_at_glucose_start, 0)` only at `G0 == 1` and draws no yield regression line.",
  "- Current redraft, make live-AUC construction match the canonical feature: addressed. The S4 live-AUC row now omits G0 > 1 points and overlays the true `lm(live_auc_glucose_window ~ log1p(G0))` fit on a linear G0 axis.",
  "- S4 role/content are otherwise preserved: all cell lines remain shown, death-feature rows remain absent, robust-derivative summary points remain highlighted, and the glucose-range inclusion support display now records the G0 <= 1 AUC and G0 = 1 yield restrictions.",
  "",
  "## Output paths",
  "",
  paste0("- `", file.path(INITIAL_DIR, outputs$filename[outputs$figure == "Figure S4"]), "`"),
  paste0("- `", file.path(REFINED_DIR, outputs$filename[outputs$figure == "Figure S4"]), "`"),
  "",
  "## Visual QC",
  "",
  "- Targeted code path uses literal linear starting-glucose axes with free x-scales per feature/line panel so the G0 <= 1 AUC and G0 = 1 yield views remain inspectable.",
  "- Regression overlays are present only for live AUC and are generated from the canonical log-G0 model; yield is a point feature at 1 mM.",
  "- The S4 panels remain dense supplemental QC displays; no new figure-level title, caption, or manuscript prose was added to the PNGs."
)
writeLines(round3_notes, file.path(WORKER_NOTES_DIR, "fg2_round3_s4_regression_coverage.md"), useBytes = TRUE)

cat(sprintf("Wrote FG2 supplement-QC PNGs to %s and %s\n", INITIAL_DIR, REFINED_DIR))
