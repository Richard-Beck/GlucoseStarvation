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

root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v5_redraft",
  "figure_generation", "FG1_measurement_foundation", "drafting_v2"
)
out_dir <- file.path(root, "initial_subpanels", "trajectory_size")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- list(
  raw_live = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "raw_live_trajectories", "wp1_raw_live_trajectories_raw_points.csv"
  ),
  live_summary = file.path("data", "report_exports", "wp1_core_starvation", "count_trajectory_summary.csv"),
  glucose_summary = file.path("data", "report_exports", "wp1_core_starvation", "glucose_trajectory_summary.csv"),
  stan_data = file.path("data", "stan_ready_data.Rds"),
  area_observations = file.path("data", "report_exports", "wp2_morphology_volume", "wp2_area_count_observations.csv"),
  adjusted_area = file.path("data", "report_exports", "wp2_morphology_volume", "wp2_adjusted_ploidy_effects.csv"),
  nuclear_high_low = file.path("data", "report_exports", "wp3_nuclear_size", "wp3_nuclear_high_vs_low_effects.csv")
)

for (path in unlist(input_paths)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}

read_tbl <- function(path) {
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), .name_repair = "unique")
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 8))
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

format_g0 <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  out <- ifelse(
    abs(y) < 1e-8,
    "0",
    ifelse(abs(y - round(y)) < 1e-8, sprintf("%.0f", y), sprintf("%.2f", y))
  )
  out <- ifelse(grepl("\\.", out), sub("0$", "", out), out)
  out <- sub("\\.$", "", out)
  out[is.na(y)] <- as.character(x)[is.na(y)]
  out
}

line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
no_sum_lines <- setdiff(line_order, "SUM-159-fuse")
g0_desc <- c("25", "5", "1", "0.5", "0.25", "0.1", "0")
g0_main <- c("25", "1", "0.25", "0")
g0_label_levels <- paste0("G0 ", g0_desc, " mM")
s1_g0_label_levels <- paste0("G0 ", rev(g0_desc), " mM")
ploidy_colors <- c(low = "#2166AC", high = "#B2182B")

add_factors <- function(df) {
  df %>%
    mutate(
      cellLine = factor(as.character(.data$cellLine), levels = line_order),
      G0_display = factor(format_g0(.data$G0), levels = g0_desc),
      G0_label = factor(paste0("G0 ", as.character(.data$G0_display), " mM"), levels = g0_label_levels),
      ploidy_state = factor(as.character(.data$ploidy_state), levels = c("low", "high"))
    )
}

add_ploidy_state <- function(df) {
  df %>%
    group_by(.data$cellLine) %>%
    mutate(ploidy_state = if_else(.data$ploidy_metric == min(.data$ploidy_metric, na.rm = TRUE), "low", "high")) %>%
    ungroup()
}

theme_fg1 <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.14, color = "grey88"),
      strip.background = element_rect(fill = "grey94", color = "grey75", linewidth = 0.22),
      strip.text = element_text(face = "bold", margin = margin(1.2, 1.2, 1.2, 1.2)),
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.9)),
      legend.key.width = unit(0.22, "in"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

save_png <- function(plot, filename, width, height, dpi = 300) {
  path <- file.path(out_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(path)
}

raw_live <- read_tbl(input_paths$raw_live) %>%
  add_factors() %>%
  mutate(hours_plot = .data$hours + .data$jitter_offset_hours * 0.42)

live_summary <- read_tbl(input_paths$live_summary) %>%
  add_ploidy_state() %>%
  add_factors()

glucose_summary <- read_tbl(input_paths$glucose_summary) %>%
  add_ploidy_state() %>%
  add_factors()

load_raw_glucose_all <- function(glucose_summary) {
  stan_data <- readRDS(input_paths$stan_data)
  well_idx <- as.integer(stan_data$well_idx_gluc)
  exp_id <- as.integer(stan_data$exp_id[well_idx])
  raw <- tibble(
    raw_obs_id = seq_along(well_idx),
    well_idx = well_idx,
    hours = as.numeric(stan_data$t_grid[as.integer(stan_data$grid_idx_gluc)]),
    lum_obs = as.numeric(stan_data$lum_obs),
    dilution = as.numeric(stan_data$dilution),
    is_censored = as.integer(stan_data$is_censored),
    calib_exp_id = exp_id,
    calib_a = as.numeric(stan_data$calib_a_fixed[exp_id]),
    calib_b = as.numeric(stan_data$calib_b_fixed[exp_id])
  ) %>%
    mutate(glucose_hat_raw = pmax(0, (.data$lum_obs - .data$calib_b) / (.data$calib_a * pmax(.data$dilution, 1e-12)))) %>%
    group_by(.data$well_idx, .data$hours) %>%
    mutate(
      raw_rep_idx = row_number(),
      n_raw_reps = n(),
      hours_plot = if_else(
        .data$n_raw_reps > 1,
        .data$hours + ((.data$raw_rep_idx - 1) / pmax(.data$n_raw_reps - 1, 1) - 0.5) * 0.18,
        .data$hours
      )
    ) %>%
    ungroup()

  meta <- glucose_summary %>%
    select(.data$well_idx, .data$hours, .data$line_id, .data$cellLine, .data$ploidy_metric,
           .data$ploidy_abs, .data$G0, .data$exp_id, .data$has_starvation, .data$ploidy_state)
  raw %>%
    left_join(meta, by = c("well_idx", "hours")) %>%
    filter(!is.na(.data$cellLine)) %>%
    add_factors()
}

raw_glucose <- load_raw_glucose_all(glucose_summary)

live_summary_all <- live_summary %>%
  select(.data$well_idx, .data$hours, .data$live_cells, .data$dead_cells, .data$cellLine,
         .data$line_id, .data$ploidy_metric, .data$ploidy_abs, .data$G0, .data$ploidy_state,
         .data$G0_display, .data$G0_label)

glucose_summary_all <- glucose_summary %>%
  select(.data$well_idx, .data$hours, .data$glucose_hat, .data$glucose_hat_lo, .data$glucose_hat_hi,
         .data$any_glucose_censored, .data$cellLine, .data$line_id, .data$ploidy_metric,
         .data$ploidy_abs, .data$G0, .data$ploidy_state, .data$G0_display, .data$G0_label)

make_live_line_block <- function(line) {
  pts <- raw_live %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$G0_display) %in% g0_main)
  lines <- live_summary_all %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$G0_display) %in% g0_main)

  ggplot() +
    geom_point(
      data = pts,
      aes(.data$hours_plot, .data$live_cells, color = .data$ploidy_state),
      size = 0.42, alpha = 1, stroke = 0, na.rm = TRUE
    ) +
    geom_line(
      data = lines,
      aes(.data$hours, .data$live_cells, color = .data$ploidy_state, group = interaction(.data$well_idx, .data$ploidy_state)),
      linewidth = 0.26, alpha = 0.90, na.rm = TRUE
    ) +
    facet_wrap(~G0_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = ploidy_colors, drop = FALSE, name = "ploidy", labels = c(low = "low", high = "high")) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = "hours", y = "live cells", title = line) +
    theme_fg1(base_size = 6.0) +
    theme(
      plot.title = element_text(face = "bold", size = 7.2, hjust = 0.5),
      strip.text = element_text(size = 5.1),
      axis.text = element_text(size = 4.9),
      axis.title = element_text(size = 5.8),
      panel.spacing = unit(0.28, "lines")
    )
}

make_glucose_line_block <- function(line) {
  pts <- raw_glucose %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$G0_display) %in% g0_main)
  lines <- glucose_summary_all %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$G0_display) %in% g0_main)

  ggplot() +
    geom_point(
      data = pts,
      aes(.data$hours_plot, .data$glucose_hat_raw, color = .data$ploidy_state, shape = factor(.data$is_censored)),
      size = 0.44, alpha = 1, stroke = 0.16, na.rm = TRUE
    ) +
    geom_line(
      data = lines,
      aes(.data$hours, .data$glucose_hat, color = .data$ploidy_state, group = interaction(.data$well_idx, .data$ploidy_state)),
      linewidth = 0.26, alpha = 0.92, na.rm = TRUE
    ) +
    facet_wrap(~G0_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = ploidy_colors, drop = FALSE, name = "ploidy", labels = c(low = "low", high = "high")) +
    scale_shape_manual(values = c("0" = 16, "1" = 4), drop = FALSE, name = "censored") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = "hours", y = "glucose (mM)", title = line) +
    theme_fg1(base_size = 6.0) +
    theme(
      plot.title = element_text(face = "bold", size = 7.2, hjust = 0.5),
      strip.text = element_text(size = 5.1),
      axis.text = element_text(size = 4.9),
      axis.title = element_text(size = 5.8),
      panel.spacing = unit(0.28, "lines")
    )
}

make_line_patchwork <- function(block_fun, lines_keep = no_sum_lines) {
  wrap_plots(lapply(lines_keep, block_fun), ncol = 2, guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
}

add_s1_display_factors <- function(df) {
  df %>%
    mutate(
      G0_s1_label = factor(paste0("G0 ", as.character(.data$G0_display), " mM"), levels = s1_g0_label_levels),
      readout = factor(.data$readout, levels = c("Live cells", "Dead cells", "Glucose"))
    )
}

s1_raw_by_readout <- bind_rows(
  raw_live %>% transmute(cellLine, G0_display, ploidy_state, hours_plot, value = live_cells, readout = "Live cells"),
  raw_live %>% transmute(cellLine, G0_display, ploidy_state, hours_plot, value = dead_cells, readout = "Dead cells"),
  raw_glucose %>% transmute(cellLine, G0_display, ploidy_state, hours_plot, value = glucose_hat_raw, readout = "Glucose")
) %>%
  add_s1_display_factors()

s1_summary_by_readout <- bind_rows(
  live_summary_all %>% transmute(cellLine, G0_display, ploidy_state, well_idx, hours, value = live_cells, readout = "Live cells"),
  live_summary_all %>% transmute(cellLine, G0_display, ploidy_state, well_idx, hours, value = dead_cells, readout = "Dead cells"),
  glucose_summary_all %>% transmute(cellLine, G0_display, ploidy_state, well_idx, hours, value = glucose_hat, readout = "Glucose")
) %>%
  add_s1_display_factors()

make_s1_readout_panel <- function(line, readout_name, show_y = TRUE) {
  raw <- s1_raw_by_readout %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$readout) == readout_name)
  summaries <- s1_summary_by_readout %>%
    filter(as.character(.data$cellLine) == line, as.character(.data$readout) == readout_name)

  ggplot() +
    geom_line(
      data = summaries,
      aes(.data$hours, .data$value, color = .data$ploidy_state, group = interaction(.data$well_idx, .data$ploidy_state)),
      linewidth = 0.15, alpha = 0.80, na.rm = TRUE
    ) +
    geom_point(data = raw, aes(.data$hours_plot, .data$value), color = "black", size = 0.20, alpha = 1, stroke = 0, na.rm = TRUE) +
    facet_grid(G0_s1_label ~ readout, scales = "free_y") +
    scale_color_manual(values = ploidy_colors, drop = FALSE, name = "ploidy", labels = c(low = "low", high = "high")) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = "hours", y = if (show_y) NULL else NULL) +
    theme_fg1(base_size = 4.45) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 3.55, margin = margin(t = 1)),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 2.85),
      axis.ticks = element_line(linewidth = 0.14),
      strip.text.y = element_text(size = 3.05, angle = 0),
      strip.text.x = element_text(size = 3.45),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2),
      panel.spacing = unit(0.22, "lines")
    )
}

make_s1_line_block <- function(line) {
  readout_grid <- make_s1_readout_panel(line, "Live cells") +
    make_s1_readout_panel(line, "Dead cells") +
    make_s1_readout_panel(line, "Glucose") +
    plot_layout(ncol = 3, widths = c(1, 1, 1), guides = "collect")

  cowplot::ggdraw() +
    cowplot::draw_plot(readout_grid, x = 0, y = 0, width = 1, height = 0.925) +
    cowplot::draw_label(line, x = 0.01, y = 0.995, hjust = 0, vjust = 1, size = 6.4, fontface = "bold")
}

make_ploidy_legend <- function() {
  legend_lines <- tibble(
    ploidy_state = factor(c("low", "high"), levels = c("low", "high")),
    x = c(0.42, 0.53),
    xend = c(0.47, 0.58),
    label_x = c(0.485, 0.595),
    label = c("low", "high")
  )
  ggplot(legend_lines) +
    geom_segment(
      aes(x = .data$x, xend = .data$xend, y = 0.5, yend = 0.5, color = .data$ploidy_state),
      linewidth = 0.35
    ) +
    geom_text(aes(.data$label_x, 0.5, label = .data$label), hjust = 0, vjust = 0.5, size = 2.3) +
    annotate("text", x = 0.38, y = 0.5, label = "ploidy", hjust = 1, vjust = 0.5, size = 2.3) +
    scale_color_manual(values = ploidy_colors, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void()
}

make_s1_atlas <- function() {
  blocks <- lapply(line_order, make_s1_line_block)
  atlas <- wrap_plots(blocks[[1]], blocks[[2]], blocks[[3]], blocks[[4]], blocks[[5]], plot_spacer(), ncol = 2) +
    plot_layout(widths = c(1, 1), heights = c(1, 1, 1))
  cowplot::plot_grid(atlas, make_ploidy_legend(), ncol = 1, rel_heights = c(1, 0.085))
}

write_s1_line_pages <- function() {
  for (line in line_order) {
    file_line <- gsub("[^A-Za-z0-9]+", "_", tolower(line))
    page <- cowplot::plot_grid(make_s1_line_block(line), make_ploidy_legend(), ncol = 1, rel_heights = c(1, 0.12))
    save_png(page, paste0("fg1_v2_figure_s1_review_page_", file_line, ".png"), 7.0, 2.45, dpi = 300)
  }
}

area_obs <- read_tbl(input_paths$area_observations) %>%
  mutate(
    cellLine = factor(as.character(.data$cellLine), levels = line_order),
    G0_display = factor(format_g0(.data$G0), levels = rev(g0_desc)),
    ploidy_state = factor(as.character(.data$ploidy_state), levels = c("baseline", "elevated"))
  )

adjusted_area <- read_tbl(input_paths$adjusted_area)
nuclear_high_low <- read_tbl(input_paths$nuclear_high_low)

make_q50_by_glucose <- function() {
  q50_summary <- area_obs %>%
    filter(is.finite(.data$q50), .data$n_area_alive >= 30) %>%
    group_by(.data$cellLine, .data$G0_display, .data$ploidy_state) %>%
    summarise(
      median_q50 = finite_median(.data$q50),
      q25_q50 = safe_quantile(.data$q50, 0.25),
      q75_q50 = safe_quantile(.data$q50, 0.75),
      n = n(),
      .groups = "drop"
    )

  ggplot(q50_summary, aes(.data$G0_display, .data$median_q50, color = .data$ploidy_state, group = .data$ploidy_state)) +
    geom_linerange(aes(ymin = .data$q25_q50, ymax = .data$q75_q50), linewidth = 0.28, alpha = 0.82) +
    geom_line(linewidth = 0.30) +
    geom_point(size = 1.05) +
    facet_wrap(~cellLine, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c(baseline = ploidy_colors[["low"]], elevated = ploidy_colors[["high"]]), name = "ploidy") +
    labs(x = "starting glucose (mM)", y = "cell area q50 (px)") +
    theme_fg1(base_size = 6.4) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.2), strip.text = element_text(size = 5.7))
}

make_ratio_summary <- function() {
  nuclear_summary <- nuclear_high_low %>%
    group_by(.data$cellLine) %>%
    summarise(
      estimate_log2 = finite_median(.data$log2_nuclear_area_ratio_elevated_vs_baseline),
      lower_log2 = safe_quantile(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.25),
      upper_log2 = safe_quantile(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.75),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(measure = "Nuclear area")

  body_summary <- adjusted_area %>%
    filter(.data$metric == "Area q50") %>%
    transmute(
      cellLine,
      measure = "Cell area q50",
      estimate_log2 = .data$log2_ratio,
      lower_log2 = .data$lower_log2_ratio,
      upper_log2 = .data$upper_log2_ratio,
      n = .data$n_model_rows
    )

  ratio_summary <- bind_rows(body_summary, nuclear_summary) %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = rev(line_order)),
      measure = factor(.data$measure, levels = c("Cell area q50", "Nuclear area")),
      sum_context = if_else(as.character(.data$cellLine) == "SUM-159-fuse", "SUM-159-fuse", "other lines")
    )

  ggplot(ratio_summary, aes(.data$estimate_log2, .data$cellLine, color = .data$measure, shape = .data$sum_context)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey45") +
    geom_errorbarh(aes(xmin = .data$lower_log2, xmax = .data$upper_log2), height = 0.13, linewidth = 0.38, na.rm = TRUE) +
    geom_point(size = 1.75, stroke = 0.65, fill = "white") +
    scale_color_manual(values = c("Cell area q50" = "#2166AC", "Nuclear area" = "#B2182B"), name = "measure") +
    scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    labs(x = "log2(elevated / baseline ploidy)", y = NULL) +
    theme_fg1(base_size = 7.0)
}

make_cell_vs_nuclear <- function() {
  condition_points <- nuclear_high_low %>%
    filter(
      is.finite(.data$cell_area_ratio_elevated_vs_baseline),
      is.finite(.data$log2_nuclear_area_ratio_elevated_vs_baseline)
    ) %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_order),
      log2_cell_area_ratio = log2(.data$cell_area_ratio_elevated_vs_baseline),
      sum_context = if_else(as.character(.data$cellLine) == "SUM-159-fuse", "SUM-159-fuse", "other lines")
    )

  summary_points <- condition_points %>%
    group_by(.data$cellLine, .data$sum_context) %>%
    summarise(
      cell_area = finite_median(.data$log2_cell_area_ratio),
      nuclear_area = finite_median(.data$log2_nuclear_area_ratio_elevated_vs_baseline),
      .groups = "drop"
    )

  ggplot(condition_points, aes(.data$log2_cell_area_ratio, .data$log2_nuclear_area_ratio_elevated_vs_baseline)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.20, color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.20, color = "grey55") +
    geom_point(color = "grey72", size = 0.75) +
    geom_point(
      data = summary_points,
      aes(.data$cell_area, .data$nuclear_area, shape = .data$sum_context),
      color = "grey10",
      fill = "white",
      size = 1.85,
      stroke = 0.70
    ) +
    geom_text(
      data = summary_points,
      aes(.data$cell_area, .data$nuclear_area, label = .data$cellLine),
      hjust = -0.08,
      size = 2.1,
      color = "grey15"
    ) +
    scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    coord_cartesian(clip = "off") +
    labs(x = "cell-area ratio, log2(elevated / baseline)", y = "nuclear-area ratio, log2(elevated / baseline)") +
    theme_fg1(base_size = 7.0)
}

write_manifest <- function(paths) {
  manifest <- tibble(
    path = paths,
    scope = "trajectory_size_worker",
    directive_ids = c(
      "FG1-D04",
      "FG1-D04",
      "FG1-D05",
      rep("FG1-D05", length(line_order)),
      "FG1-D06",
      "FG1-D06",
      "FG1-D06"
    ),
    notes = c(
      "Figure 1D live trajectory draft with per-line patchwork blocks, free y scales, opaque raw points, and x-only jitter.",
      "Figure 1E glucose trajectory draft with matched per-line patchwork blocks, free y scales, opaque raw points, and x-only jitter.",
      "Standalone Figure S1 raw atlas draft with no visible S1A label and visible low/high ploidy legend.",
      rep("Per-line S1 review page for density/readability inspection.", length(line_order)),
      "Early size-rationale option: q50 cell-area trajectories by starting glucose.",
      "Early size-rationale option: simple cell-area q50 and nuclear-area high/low ratio summary.",
      "Early size-rationale option: condition-level cell-area versus nuclear-area ratio scatter."
    )
  )
  write.csv(as.data.frame(manifest), file.path(out_dir, "trajectory_size_manifest.csv"), row.names = FALSE)
}

paths <- c(
  save_png(make_line_patchwork(make_live_line_block), "fg1_v2_figure1d_live_free_line_scales.png", 7.0, 5.8),
  save_png(make_line_patchwork(make_glucose_line_block), "fg1_v2_figure1e_glucose_free_line_scales.png", 7.0, 5.8),
  save_png(make_s1_atlas(), "fg1_v2_figure_s1_raw_atlas_standalone.png", 7.0, 7.25),
  {
    write_s1_line_pages()
    file.path(out_dir, paste0("fg1_v2_figure_s1_review_page_", gsub("[^A-Za-z0-9]+", "_", tolower(line_order)), ".png"))
  },
  save_png(make_q50_by_glucose(), "fg1_v2_early_size_q50_by_glucose.png", 7.0, 2.35),
  save_png(make_ratio_summary(), "fg1_v2_early_size_area_nuclear_ratio_summary.png", 4.9, 2.45),
  save_png(make_cell_vs_nuclear(), "fg1_v2_early_size_cell_vs_nuclear_ratio.png", 4.9, 3.0)
)

write_manifest(paths)
cat("Wrote trajectory/size drafts to ", out_dir, "\n", sep = "")
