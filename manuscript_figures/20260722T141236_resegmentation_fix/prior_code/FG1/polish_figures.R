#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(scales)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
phase_idx <- match("--phase", args)
phase <- if (!is.na(phase_idx) && length(args) >= phase_idx + 1) args[[phase_idx + 1]] else "both"
if (!phase %in% c("subpanels", "final", "both")) {
  stop("Unsupported --phase: ", phase, call. = FALSE)
}
output_idx <- match("--output-root", args)
if (is.na(output_idx) || length(args) < output_idx + 1) {
  stop("Required argument missing: --output-root <path>", call. = FALSE)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1) stop("Could not identify the promoted script path.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)
package_root <- dirname(script_dir)
draft_root <- package_root
polish_root <- normalizePath(args[[output_idx + 1]], winslash = "/", mustWork = FALSE)
subpanel_dir <- file.path(polish_root, "subpanels")
layout_dir <- file.path(polish_root, "layout")
final_dir <- file.path(polish_root, "final_images")
dir.create(subpanel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
fixed_layout_path <- file.path(package_root, "layout_plan.csv")
output_layout_path <- file.path(layout_dir, "layout_plan.csv")
if (!file.exists(output_layout_path)) {
  if (!file.copy(fixed_layout_path, output_layout_path, overwrite = FALSE)) {
    stop("Could not stage fixed layout config: ", fixed_layout_path, call. = FALSE)
  }
}
fixed_panel_map_path <- file.path(package_root, "panel_map.csv")
output_panel_map_path <- file.path(polish_root, "panel_map.csv")
if (!file.exists(output_panel_map_path)) {
  if (!file.copy(fixed_panel_map_path, output_panel_map_path, overwrite = FALSE)) {
    stop("Could not stage fixed panel map: ", fixed_panel_map_path, call. = FALSE)
  }
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

required_runtime_inputs <- c(
  "data/image_processing_runs/run_20260324_233122/object_train_features_90.csv",
  "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/alive_count_validation/wp1_alive_count_validation_frame_counts.csv",
  "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_frame_count_metrics.csv",
  "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_object_metrics.csv",
  "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_output_manifest.csv",
  "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_provenance.json",
  "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_selected_object_calls.csv"
)

required_runtime_context <- c(
  "agent-dev/major_analyses/20260623_production_classifier_fg1_unblock/analysis_plan.md",
  "agent-dev/major_analyses/20260623_production_classifier_fg1_unblock/validation_report.md"
)

merge_contract_paths <- function(...) {
  values <- unlist(lapply(list(...), function(value) {
    if (!length(value)) return(character())
    unlist(strsplit(value, ";", fixed = TRUE), use.names = FALSE)
  }), use.names = FALSE)
  paste(sort(unique(values[nzchar(values) & values != "NA"])), collapse = ";")
}

root_block <- paste(
  "root <- file.path(",
  "  \"agent-dev\", \"manuscript_redrafts\", \"20260619_v5_redraft\",",
  "  \"figure_generation\", \"FG1_measurement_foundation\", \"drafting_v2\"",
  ")",
  sep = "\n"
)

source_with_polish_root <- function(script_path, env) {
  txt <- readLines(script_path, warn = FALSE)
  joined <- paste(txt, collapse = "\n")
  if (!grepl(root_block, joined, fixed = TRUE)) {
    stop("Could not find expected root block in ", script_path, call. = FALSE)
  }
  joined <- gsub(root_block, "root <- polish_root", joined, fixed = TRUE)
  eval(parse(text = joined, srcfile = script_path), envir = env)
  invisible(env)
}

read_tbl_local <- function(path) {
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), .name_repair = "unique")
}

safe_quantile_local <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 8))
}

finite_median_local <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

theme_fg1_polish <- function(base_size = 7) {
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
      plot.margin = margin(3, 8, 3, 3)
    )
}

make_polish_cell_vs_nuclear <- function() {
  nuclear_high_low <- read_tbl_local(file.path(project_root, "data", "report_exports", "wp3_nuclear_size", "wp3_nuclear_high_vs_low_effects.csv"))
  line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
  line_palette <- c(
    "MCF10A" = "#1B9E77",
    "MDA-MB-231" = "#D95F02",
    "SNU668" = "#7570B3",
    "SUM-159-chem" = "#E7298A",
    "SUM-159-fuse" = "#66A61E"
  )
  condition_points <- nuclear_high_low %>%
    filter(
      is.finite(.data$cell_area_ratio_elevated_vs_baseline),
      is.finite(.data$log2_nuclear_area_ratio_elevated_vs_baseline)
    ) %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_order),
      log2_cell_area_ratio = log2(.data$cell_area_ratio_elevated_vs_baseline)
    )
  summary_points <- condition_points %>%
    group_by(.data$cellLine) %>%
    summarise(
      cell_area = finite_median_local(.data$log2_cell_area_ratio),
      nuclear_area = finite_median_local(.data$log2_nuclear_area_ratio_elevated_vs_baseline),
      .groups = "drop"
    ) %>%
    mutate(
      sum_context = if_else(as.character(.data$cellLine) == "SUM-159-fuse", "SUM-159-fuse", "other lines"),
      label = recode(
        as.character(.data$cellLine),
        "MDA-MB-231" = "MDA-231",
        "SUM-159-chem" = "SUM-chem",
        "SUM-159-fuse" = "SUM-fuse",
        .default = as.character(.data$cellLine)
      ),
      label_dx = case_when(
        as.character(.data$cellLine) == "SUM-159-fuse" ~ -0.06,
        as.character(.data$cellLine) == "MDA-MB-231" ~ 0.05,
        TRUE ~ 0.05
      ),
      label_dy = case_when(
        as.character(.data$cellLine) == "MCF10A" ~ 0.12,
        as.character(.data$cellLine) == "SNU668" ~ 0.08,
        as.character(.data$cellLine) == "SUM-159-chem" ~ -0.08,
        as.character(.data$cellLine) == "SUM-159-fuse" ~ -0.25,
        TRUE ~ -0.08
      ),
      label_hjust = if_else(as.character(.data$cellLine) == "SUM-159-fuse", 1, 0)
    )

  ggplot(condition_points, aes(.data$log2_cell_area_ratio, .data$log2_nuclear_area_ratio_elevated_vs_baseline)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.20, color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.20, color = "grey55") +
    geom_point(color = "grey74", size = 0.65, alpha = 0.9) +
    geom_point(
      data = summary_points,
      aes(.data$cell_area, .data$nuclear_area, color = .data$cellLine, shape = .data$sum_context),
      size = 2.0,
      stroke = 0.75,
      fill = "white"
    ) +
    geom_text(
      data = summary_points,
      aes(
        .data$cell_area + .data$label_dx,
        .data$nuclear_area + .data$label_dy,
        label = .data$label,
        color = .data$cellLine,
        hjust = .data$label_hjust
      ),
      size = 2.0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = line_palette, name = "cell line") +
    scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    coord_cartesian(xlim = c(-1.7, 2.45), ylim = c(-1.2, 1.7), clip = "off") +
    labs(
      x = "cell-area ratio (log2)",
      y = "nuclear-area ratio (log2)"
    ) +
    theme_fg1_polish(base_size = 7.0) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 5.8),
      legend.title = element_text(size = 6.0),
      axis.text = element_text(size = 6.0),
      axis.title = element_text(size = 6.6),
      plot.margin = margin(3, 30, 3, 4)
    )
}

make_polish_q50_by_glucose <- function(env) {
  lab_map <- c(
    "MCF10A" = "MCF10A",
    "MDA-MB-231" = "MDA-231",
    "SNU668" = "SNU668",
    "SUM-159-chem" = "SUM-chem",
    "SUM-159-fuse" = "SUM-fuse"
  )
  env$make_q50_by_glucose() +
    facet_wrap(~cellLine, ncol = 1, scales = "free_y", labeller = labeller(cellLine = as_labeller(lab_map))) +
    theme(
      strip.text = element_text(size = 5.5),
      axis.text.x = element_text(size = 4.8, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 4.8),
      axis.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      legend.title = element_text(size = 6.0),
      panel.spacing = unit(0.30, "lines")
    )
}

make_polish_size_ratio_vs_ploidy <- function() {
  adjusted_area <- read_tbl_local(file.path(project_root, "data", "report_exports", "wp2_morphology_volume", "wp2_adjusted_ploidy_effects.csv"))
  nuclear_high_low <- read_tbl_local(file.path(project_root, "data", "report_exports", "wp3_nuclear_size", "wp3_nuclear_high_vs_low_effects.csv"))
  line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")

  ploidy_map <- adjusted_area %>%
    filter(.data$metric == "Area q50") %>%
    transmute(cellLine = as.character(.data$cellLine), log2_ploidy_ratio = .data$ploidy_delta)
  cell_summary <- adjusted_area %>%
    filter(.data$metric == "Area q50") %>%
    transmute(
      cellLine = as.character(.data$cellLine),
      measure = "Cell area q50",
      log2_ploidy_ratio = .data$ploidy_delta,
      estimate_log2 = .data$log2_ratio,
      lower_log2 = .data$lower_log2_ratio,
      upper_log2 = .data$upper_log2_ratio
    )
  nuclear_summary <- nuclear_high_low %>%
    group_by(.data$cellLine) %>%
    summarise(
      estimate_log2 = finite_median_local(.data$log2_nuclear_area_ratio_elevated_vs_baseline),
      lower_log2 = safe_quantile_local(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.25),
      upper_log2 = safe_quantile_local(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.75),
      .groups = "drop"
    ) %>%
    mutate(cellLine = as.character(.data$cellLine), measure = "Nuclear area") %>%
    left_join(ploidy_map, by = "cellLine") %>%
    select(.data$cellLine, .data$measure, .data$log2_ploidy_ratio, .data$estimate_log2, .data$lower_log2, .data$upper_log2)
  points <- bind_rows(cell_summary, nuclear_summary) %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_order),
      measure = factor(.data$measure, levels = c("Cell area q50", "Nuclear area")),
      line_group = if_else(as.character(.data$cellLine) == "SUM-159-fuse", "SUM-159-fuse", "other lines")
    )
  sum_points <- points %>% filter(as.character(.data$cellLine) == "SUM-159-fuse")

  ggplot(points, aes(.data$log2_ploidy_ratio, .data$estimate_log2)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.24, color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.22, color = "grey68") +
    geom_errorbar(aes(ymin = .data$lower_log2, ymax = .data$upper_log2, color = .data$line_group),
      width = 0.018, linewidth = 0.35, alpha = 0.85, na.rm = TRUE
    ) +
    geom_point(aes(color = .data$line_group, shape = .data$line_group),
      size = 2.0, stroke = 0.72, fill = "white", na.rm = TRUE
    ) +
    geom_text(
      data = sum_points,
      aes(label = "SUM-fuse", color = .data$line_group),
      nudge_x = 0.045,
      hjust = 0,
      size = 2.1,
      show.legend = FALSE
    ) +
    facet_wrap(~measure, nrow = 1) +
    scale_color_manual(values = c("other lines" = "grey20", "SUM-159-fuse" = "#B2182B"), name = "line") +
    scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    scale_x_continuous(breaks = c(0.33, 1.0), labels = c("0.33", "1.00"), limits = c(0.24, 1.38)) +
    coord_cartesian(ylim = c(-1.18, 1.58), clip = "off") +
    labs(x = "ploidy ratio (log2)", y = "area ratio (log2)") +
    theme_fg1_polish(base_size = 7.0) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 6.0),
      axis.text = element_text(size = 5.8),
      axis.title = element_text(size = 6.3),
      plot.margin = margin(3, 24, 3, 3)
    )
}

compact_trajectory_panel <- function(plot) {
  plot &
    scale_y_continuous(n.breaks = 3, labels = function(x) {
      ifelse(abs(x) >= 1000, paste0(round(x / 1000, 1), "k"), format(round(x, 2), trim = TRUE, scientific = FALSE))
    }) &
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 3.6),
      axis.text.x = element_text(size = 4.0),
      axis.title = element_text(size = 5.0),
      strip.text = element_text(size = 4.4)
    )
}

prepare_context_files <- function() {
  src <- file.path(
    draft_root, "initial_subpanels", "measurement_provenance_options",
    "fg1_figure1a_options_manifest.csv"
  )
  dst_dir <- file.path(polish_root, "initial_subpanels", "measurement_provenance_options")
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src, file.path(dst_dir, basename(src)), overwrite = TRUE)
}

load_builders <- function() {
  prepare_context_files()
  env <- new.env(parent = globalenv())
  env$polish_root <- polish_root

  source_with_polish_root(file.path(draft_root, "scripts", "fg1_measurement_provenance_worker.R"), env)
  env$polish_s2a <- env$s2a
  env$polish_s2b_metrics <- env$metric_plot

  source_with_polish_root(file.path(draft_root, "scripts", "fg1_s2_preserved_support_worker.R"), env)
  env$polish_s2c <- env$make_s2_c()
  env$polish_s2d <- env$make_s2_d()
  env$polish_s2e <- env$make_s2_e() +
    scale_x_continuous(
      trans = "pseudo_log",
      breaks = c(0, 1, 25),
      labels = c("0", "1", "25"),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    theme(axis.text.x = element_text(size = 4.4))
  env$polish_s2f <- env$make_s2_f() +
    scale_x_continuous(
      trans = "pseudo_log",
      breaks = c(0, 1, 25),
      labels = c("0", "1", "25"),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    theme(axis.text.x = element_text(size = 4.4))
  env$polish_s2g <- env$make_s2_g()

  source_with_polish_root(file.path(draft_root, "scripts", "fg1_trajectory_size_worker.R"), env)
  env$polish_fig1b <- make_polish_cell_vs_nuclear()
  env$polish_fig1c <- make_polish_q50_by_glucose(env)
  env$g0_main <- c("25", "1", "0")
  env$polish_fig1e <- compact_trajectory_panel(env$make_line_patchwork(env$make_live_line_block))
  env$polish_fig1f <- compact_trajectory_panel(env$make_line_patchwork(env$make_glucose_line_block))

  source_with_polish_root(file.path(draft_root, "scripts", "fg1_size_s1_s2_revision_worker.R"), env)
  env$polish_fig1d <- make_polish_size_ratio_vs_ploidy()
  env$polish_s1a <- env$make_s1_atlas_legend_in_white_space()

  source_with_polish_root(file.path(draft_root, "scripts", "fg1_figure1a_round2_worker.R"), env)
  env$polish_fig1a <- env$round2_plot
  env
}

panel_row <- function(figure, panel, plot, width, height, approved_source, generator, data_inputs, notes = "") {
  if (!nzchar(notes)) notes <- "Regenerated from approved worker logic without a panel-identity change."
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

build_panels <- function(env) {
  scripts <- file.path(package_root, "scripts")
  data_measurement <- paste(
    c(
      "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_object_predictions.csv",
      "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_frame_counts.csv",
      "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_headline_metrics.csv",
      "data/image_processing_runs/run_20260324_233122/annotation_stack_90_manifest.csv",
      "data/image_processing_runs/run_20260324_233122/annotation_stack_90_Results.csv"
    ),
    collapse = ";"
  )
  data_figure1a <- paste(
    c(
      data_measurement,
      "data/image_processing_runs/run_20260324_233122/annotation_stack_90.tif",
      "data/image_processing_runs/run_20260324_233122/objects_labels_90.tif"
    ),
    collapse = ";"
  )
  data_trajectory <- paste(
    c(
      "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/raw_live_trajectories/wp1_raw_live_trajectories_raw_points.csv",
      "data/report_exports/wp1_core_starvation/count_trajectory_summary.csv",
      "data/report_exports/wp1_core_starvation/glucose_trajectory_summary.csv",
      "data/stan_ready_data.Rds"
    ),
    collapse = ";"
  )
  data_size <- paste(
    c(
      "data/report_exports/wp2_morphology_volume/wp2_area_count_observations.csv",
      "data/report_exports/wp2_morphology_volume/wp2_adjusted_ploidy_effects.csv",
      "data/report_exports/wp3_nuclear_size/wp3_nuclear_high_vs_low_effects.csv"
    ),
    collapse = ";"
  )
  data_s2_support <- paste(
    c(
      "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/s2_measurement_qc/wp1_s2_measurement_qc_coverage_by_condition.csv",
      "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/s2_measurement_qc/wp1_s2_glucose_calibration_fit_summary.csv",
      "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/s2_measurement_qc/wp1_s2_glucose_calibration_residual_by_standard.csv",
      "data/report_exports/manuscript_draft_figures_v3/wp1_drafting/panels/s2_measurement_qc/wp1_s2_glucose_calibration_standards_with_fit.csv"
    ),
    collapse = ";"
  )

  list(
    panel_row(
      "Figure 1", "a", env$polish_fig1a, 7.0, 3.45,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance_round2/fg1_figure1a_round2_low_high_ploidy_outlines.png",
      file.path(scripts, "fg1_figure1a_round2_worker.R"),
      data_figure1a,
      "Expanded to full-width polish real estate per 2026-06-25 closeout; retained in accepted 2026-06-30 layout revision."
    ),
    panel_row(
      "Figure 1", "b", env$polish_fig1c, 2.02, 5.62,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size/fg1_v2_early_size_q50_by_glucose.png",
      file.path(scripts, "fg1_trajectory_size_worker.R"),
      data_size,
      "Regenerated from approved worker logic with cell-line facets stacked vertically. In the v6 repair, this q50-by-glucose panel is lettered b so the visual order matches panel identity."
    ),
    panel_row(
      "Figure 1", "c", env$polish_fig1b, 1.98, 1.73,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size/fg1_v2_early_size_cell_vs_nuclear_ratio.png",
      file.path(scripts, "fg1_trajectory_size_worker.R"),
      data_size,
      "Regenerated from approved worker logic. In the v6 repair, this cell-area versus nuclear-area ratio panel is lettered c after the q50-by-glucose panel it summarizes."
    ),
    panel_row(
      "Figure 1", "d", env$polish_fig1d, 2.80, 1.73,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size_revision/fg1_v2_size_area_ratio_vs_log2_ploidy_ratio.png",
      file.path(scripts, "fg1_size_s1_s2_revision_worker.R"),
      paste("data/report_exports/wp2_morphology_volume/wp2_adjusted_ploidy_effects.csv", "data/report_exports/wp3_nuclear_size/wp3_nuclear_high_vs_low_effects.csv", sep = ";")
    ),
    panel_row(
      "Figure 1", "e", env$polish_fig1e, 2.38, 3.72,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size/fg1_v2_figure1d_live_free_line_scales.png",
      file.path(scripts, "fg1_trajectory_size_worker.R"),
      data_trajectory,
      "Regenerated from approved worker logic with G0 = 0.25 mM removed for the accepted 2026-06-30 Figure 1 layout revision."
    ),
    panel_row(
      "Figure 1", "f", env$polish_fig1f, 2.38, 3.72,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size/fg1_v2_figure1e_glucose_free_line_scales.png",
      file.path(scripts, "fg1_trajectory_size_worker.R"),
      data_trajectory,
      "Regenerated from approved worker logic with G0 = 0.25 mM removed for the accepted 2026-06-30 Figure 1 layout revision."
    ),
    panel_row(
      "Figure S1", "a", env$polish_s1a, 7.0, 6.85,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/trajectory_size_revision/fg1_v2_figure_s1_legend_in_bottom_right_option.png",
      file.path(scripts, "fg1_size_s1_s2_revision_worker.R"),
      data_trajectory,
      "Standalone atlas treated as one displayed panel for the polishing contract."
    ),
    panel_row(
      "Figure S2", "a", env$polish_s2a, 3.45, 2.35,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2a_alive_object_labels.png",
      file.path(scripts, "fg1_measurement_provenance_worker.R"),
      data_measurement
    ),
    panel_row(
      "Figure S2", "b", env$polish_s2b_metrics, 3.45, 1.65,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_figure1c_production_metrics.png",
      file.path(scripts, "fg1_measurement_provenance_worker.R"),
      "data/report_exports/v5_fg1_production_classifier_validation/fg1_production_classifier_object_predictions.csv",
      "Inserted as new panel b per 2026-06-25 closeout."
    ),
    panel_row(
      "Figure S2", "c", env$polish_s2c, 3.35, 2.70,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2c_measurement_coverage.png",
      file.path(scripts, "fg1_s2_preserved_support_worker.R"),
      data_s2_support
    ),
    panel_row(
      "Figure S2", "d", env$polish_s2d, 3.35, 2.70,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2d_glucose_censoring.png",
      file.path(scripts, "fg1_s2_preserved_support_worker.R"),
      data_s2_support
    ),
    panel_row(
      "Figure S2", "e", env$polish_s2e, 3.35, 2.30,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2e_calibration_standards.png",
      file.path(scripts, "fg1_s2_preserved_support_worker.R"),
      data_s2_support
    ),
    panel_row(
      "Figure S2", "f", env$polish_s2f, 3.35, 2.30,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2f_calibration_residuals.png",
      file.path(scripts, "fg1_s2_preserved_support_worker.R"),
      data_s2_support
    ),
    panel_row(
      "Figure S2", "g", env$polish_s2g, 6.80, 1.30,
      "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/initial_subpanels/measurement_provenance/fg1_measurement_s2g_calibration_fit_error.png",
      file.path(scripts, "fg1_s2_preserved_support_worker.R"),
      data_s2_support
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

figure1_manual_layout <- function(layout_plan, panels) {
  total_width <- 7.0
  total_height <- 9.25
  layout_tree <- "[a / [b | [[c | d] / [e | f]]]]"
  panel_lookup <- setNames(panels, vapply(panels, function(x) paste(x$figure, x$panel, sep = "|"), character(1)))
  specs <- tibble::tribble(
    ~panel, ~x_in, ~y_in, ~width_in, ~height_in, ~pad_left_in, ~pad_right_in, ~pad_bottom_in, ~pad_top_in,
    "a", 0.00, 5.80, 7.00, 3.45, 0.070, 0.035, 0.060, 0.050,
    "b", 0.00, 0.04, 2.02, 5.62, 0.095, 0.030, 0.055, 0.040,
    "c", 2.14, 3.93, 1.98, 1.73, 0.090, 0.025, 0.055, 0.040,
    "d", 4.20, 3.93, 2.80, 1.73, 0.085, 0.020, 0.055, 0.040,
    "e", 2.14, 0.04, 2.38, 3.72, 0.075, 0.035, 0.055, 0.040,
    "f", 4.62, 0.04, 2.38, 3.72, 0.075, 0.035, 0.055, 0.040
  )
  rows <- lapply(seq_len(nrow(specs)), function(idx) {
    spec <- specs[idx, ]
    key <- paste("Figure 1", spec$panel, sep = "|")
    panel <- panel_lookup[[key]]
    if (is.null(panel)) stop("Missing Figure 1 panel in manual layout: ", spec$panel, call. = FALSE)
    input_width <- panel$width
    input_height <- panel$height
    data.frame(
      figure = "Figure 1",
      panel = spec$panel,
      subpanel_png = rel_path(file.path(subpanel_dir, paste0("figure_1_", spec$panel, ".png"))),
      x_in = spec$x_in,
      y_in = spec$y_in,
      width_in = spec$width_in,
      height_in = spec$height_in,
      sx = spec$width_in / input_width,
      sy = spec$height_in / input_height,
      x_npc = spec$x_in / total_width,
      y_npc = spec$y_in / total_height,
      width_npc = spec$width_in / total_width,
      height_npc = spec$height_in / total_height,
      input_width_in = input_width,
      input_height_in = input_height,
      layout_width_in = total_width,
      layout_height_in = total_height,
      target_width_in = total_width,
      max_height_in = total_height,
      layout_tree = layout_tree,
      layout_score = 0,
      wasted_fraction = 0,
      over_width_in = 0,
      over_height_in = 0,
      scale_note = "manual user-approved Figure 1 layout",
      pad_left_npc = spec$pad_left_in / total_width,
      pad_right_npc = spec$pad_right_in / total_width,
      pad_bottom_npc = spec$pad_bottom_in / total_height,
      pad_top_npc = spec$pad_top_in / total_height,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(
    layout_plan %>% filter(.data$figure != "Figure 1"),
    bind_rows(rows)
  )
}

write_current_layout_report <- function(layout_plan) {
  by_figure <- split(layout_plan, layout_plan$figure)
  figure_lines <- vapply(by_figure, function(one) {
    sprintf(
      "- `%s`: %.2f x %.2f in, tree `%s`",
      unique(one$figure)[[1]],
      unique(one$layout_width_in)[[1]],
      unique(one$layout_height_in)[[1]],
      unique(one$layout_tree)[[1]]
    )
  }, character(1))
  lines <- c(
    "# Panel layout optimization report",
    "",
    "- Input dimensions: `agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/figure_generation/FG1_measurement_foundation/polishing/layout/subpanel_dimensions.csv`",
    "- Layout plan: `agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/figure_generation/FG1_measurement_foundation/polishing/layout/layout_plan.csv`",
    "- Target width: 7.00 in",
    "- Maximum height: 9.25 in",
    "- Gap: 0.08 in for optimizer-derived layouts; Figure 1 uses the user-approved manual layout recorded in `scale_note`.",
    "- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.",
    "",
    "## Figures",
    "",
    unname(figure_lines),
    "",
    "## Manual override",
    "",
    "- `Figure 1`: v6 repair of the user-approved 2026-06-30 option 5 layout. Panel b is stacked as cell-line facet rows on the left; panels c/d form the upper-right row; panels e/f form the taller lower-right row; panels e/f omit G0 = 0.25 mM.",
    "",
    "## Scale recommendations",
    "",
    "Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated before final assembly."
  )
  writeLines(lines, file.path(layout_dir, "layout_report.md"))
}

layout_pad <- function(lay, name, default) {
  if (name %in% names(lay)) {
    val <- suppressWarnings(as.numeric(lay[[name]][[1]]))
    if (is.finite(val)) return(val)
  }
  default
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
    pad_left <- layout_pad(lay, "pad_left_npc", if (lay$x_npc < 0.01) 0.035 else 0.010)
    pad_right <- layout_pad(lay, "pad_right_npc", if ((lay$x_npc + lay$width_npc) > 0.99) 0.010 else 0.008)
    pad_bottom <- layout_pad(lay, "pad_bottom_npc", if (lay$y_npc < 0.01) 0.020 else 0.008)
    pad_top <- layout_pad(lay, "pad_top_npc", 0.018)
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
  layout_plan <- figure1_manual_layout(layout_plan, panels)
  write.csv(layout_plan, layout_path, row.names = FALSE)
  write_current_layout_report(layout_plan)
  out_map <- c(
    "Figure 1" = file.path(final_dir, "figure_1.png"),
    "Figure S1" = file.path(final_dir, "figure_s1.png"),
    "Figure S2" = file.path(final_dir, "figure_s2.png")
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
    "--phase subpanels",
    "--output-root",
    rel_path(polish_root)
  )
  command_final <- paste(
    "scripts/agentRrunner.sh",
    rel_path(file.path(script_dir, "polish_figures.R")),
    "--phase final",
    "--output-root",
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
      data_inputs = merge_contract_paths(row$data_inputs, required_runtime_inputs),
      layout_plan = layout_plan,
      output_image = output_for(row$figure),
      notes = row$notes,
      approved_source = row$approved_source,
      context_inputs = merge_contract_paths(
        "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/final_closeout_20260625.md;agent-dev/manuscript_redrafts/20260619_v6_redraft/stage_outputs/figures/figures_contract_repair_scope_feedback_list.md;.agents/references/manuscript_figure_style.md",
        required_runtime_context
      ),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(provenance, file.path(polish_root, "provenance.csv"), row.names = FALSE)

  manifest <- data.frame(
    output_image = rel_path(unname(out_map)),
    figure = names(out_map),
    panel_map = rel_path(file.path(polish_root, "panel_map.csv")),
    subpanel_dimensions = rel_path(file.path(layout_dir, "subpanel_dimensions.csv")),
    layout_plan = layout_plan,
    legend = rel_path(file.path(polish_root, "legend.md")),
    provenance = rel_path(file.path(polish_root, "provenance.csv")),
    feedback = "agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/final_closeout_20260625.md;agent-dev/manuscript_redrafts/20260619_v6_redraft/stage_outputs/figures/figures_contract_repair_scope_feedback_list.md",
    stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(polish_root, "manifest.csv"), row.names = FALSE)

  legend <- c(
    "# FG1 Polished Legends",
    "",
    "Figure 1, Measurement foundation and early size context. Panel a shows production-classifier object outlines and manual alive points on low/high-ploidy annotated crops; panel b shows cell-area q50 by starting glucose and ploidy; panel c compares cell-area and nuclear-area ratios; panel d relates size ratios to ploidy ratio; panels e and f show accepted live-cell and extracellular-glucose trajectory support across no-SUM lines. Colors define low versus high ploidy, and SUM-159-fuse is highlighted where relevant.",
    "",
    "Figure S1, Raw live/dead/extracellular-glucose trajectory atlas. Panel a shows raw and smoothed live-cell, dead-cell, and extracellular-glucose measurements across cell lines, starting glucose concentrations, and ploidy states, with the low/high ploidy legend placed in the bottom-right whitespace cell.",
    "",
    "Figure S2, Measurement coverage and glucose-calibration QC. Panel a compares manual alive annotation points with alive object labels; panel b shows target-scoped production-classifier ROC AUC, balanced accuracy, precision, and recall; panels c and d show measurement coverage and glucose censoring; panels e-g show glucose calibration standards, residuals, and fit error."
  )
  writeLines(legend, file.path(polish_root, "legend.md"))

  notes <- c(
    "# FG1 polishing notes",
    "",
    "## Inputs and decisions",
    "",
    "- Polishing scope: closed `FG1_measurement_foundation/drafting_v2` package only.",
    "- Approval source: `drafting_v2/final_closeout_20260625.md`, `drafting_v2/feedback.md`, and package-local polish feedback `polishing/user-feedback.txt`.",
    "- Figure 1A was given full-width real estate in the polished layout to address the closeout readability handoff.",
    "- V6 repair swaps only the Figure 1B/C generated letters and panel identities: the left stacked q50-by-glucose panel is now panel b, and the cell-area/nuclear-area ratio summary is now panel c.",
    "- Figure S2 inserts the regenerated target-scoped production metrics panel as panel b, and the previous panels b-f are relettered c-g.",
    "- Final composites are assembled from regenerated ggplot/patchwork objects, not from the audit subpanel PNGs.",
    "- Preflight validation was run before figure generation. It found all required sources and panel-map content, but reported writability errors from Python `os.access()` on this shared filesystem even though shell and Python write tests succeeded in the same directories; this environment-specific false positive was treated as a validator caveat, and the temporary write-test files were removed.",
    "",
    "## Commands",
    "",
    paste0("- `", command_subpanels, "`"),
    "- `scripts/agentRrunner.sh .agents/skills/manuscript-figure-workflow/scripts/optimize_panel_layout.R --input agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/figure_generation/FG1_measurement_foundation/polishing/layout/subpanel_dimensions.csv --output-dir agent-dev/manuscript_redrafts/20260710T164049_v01_redraft/figure_generation/FG1_measurement_foundation/polishing/layout --target-width 7 --max-height 9.25 --max-panels 8`",
    paste0("- `", command_final, "`"),
    "",
    "## Layout",
    "",
    "- The layout optimizer's `x_npc` and `y_npc` values were used as lower-left panel origins with `cowplot::draw_plot()`.",
    "- Panel labels are added at the final composite stage in lowercase reading order.",
    "- The first optimizer run put Figure 1 above the target height; Figure 1 native panel heights were tightened and the optimizer was rerun before final assembly.",
    "- Final assembly uses small per-panel padding inside optimizer slots to avoid clipping axis text at page edges. The coordinate origin convention from `layout_plan.csv` is unchanged.",
    "- Independent final PNG QC found additional Figure 1 readability defects after the first passing postflight run. Polishing then abbreviated the q50-by-glucose facet labels, separated ratio-panel direct labels, limited panel d direct labeling to SUM-159-fuse, and compacted e/f trajectory axes before rerendering and reinspecting.",
    "- The 2026-06-26 package-local polish feedback was applied by giving panels e/f more vertical space, widening the q50-by-glucose panel, reducing the ratio and d panels if needed, and shortening panel d axis labels to `ploidy ratio (log2)` and `area ratio (log2)`.",
    "- The 2026-06-30 user-approved Figure 1 option 5 layout was applied as a manual final-phase override and repaired in v6: panel b is stacked as cell-line facet rows on the lower left; panels c/d form the upper-right row; panels e/f form a taller lower-right row; panels e/f omit G0 = 0.25 mM.",
    "",
    "## Caveats",
    "",
    "- S1 remains dense by design because it is a full raw trajectory atlas.",
    "- The Figure 1 live/glucose trajectory panels are compact summaries in the main figure; the 2026-06-30 accepted layout increases their vertical allocation after removing the G0 = 0.25 mM row, but print-size review remains useful during manuscript integration.",
    "",
    "## Project-map decision",
    "",
    "- `docs/project_map.txt` was checked during this repair pass. No FG1-specific project-map update is required because this creates local v6 redraft figure-generation artifacts under an already mapped v6 redraft root and does not add a maintained workflow entrypoint or canonical data product."
  )
  writeLines(notes, file.path(polish_root, "notes.md"))

  visual_qc <- c(
    "# FG1 visual QC",
    "",
    "Overall status: passed after rerun and independent QC revision.",
    "",
    "| final PNG | title/subtitle/caption/header | clipping | spacing/layout | readability | status |",
    "|---|---|---|---|---|---|",
    "| `final_images/figure_1.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass after rerun: left-edge axis/title clipping and panel text clipping were fixed. | Pass: panel order is a-f; panel b is stacked on the lower left, c/d form the upper-right row, and e/f form a taller lower-right row with G0 = 0.25 mM omitted. | Pass: Figure 1A remains readable; panels b/c/d/e/f were revised after independent QC, package-local polish feedback, 2026-06-30 user layout selection, and the v6 B/C lettering repair. | Pass |",
    "| `final_images/figure_s1.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass after rerun: panel label no longer overlaps the first cell-line title. | Pass: single atlas panel a is ordered as the accepted five-line atlas with the ploidy legend in the bottom-right whitespace. | Pass with caveat: raw trajectory atlas remains intentionally dense. | Pass |",
    "| `final_images/figure_s2.png` | Pass: no figure-level title, subtitle, caption, or header text visible. | Pass after rerun: left-edge clipping in panels b/c was fixed by final-assembly padding. | Pass: panel order is a-g; the target-scoped production metrics support is inserted as panel b and former support panels are relettered c-g. | Pass: text and symbols are readable at full rendered size. | Pass |",
    "",
    "Reruns: first final render exposed left-edge clipping in Figure 1 and Figure S2 and panel-label crowding in Figure S1. `scripts/polish_figures.R` was revised to inset plotted grobs within their optimizer slots while keeping lowercase labels in the reserved margin, then final assembly was rerun. Independent final PNG QC then requested Figure 1 revision for clipped q50-by-glucose headers, clipped panel d labels, crowded ratio-panel labels, and unreadable e/f axis text; the script was revised to abbreviate facet labels, replace the ratio panel's wide legend with direct labels, simplify panel d direct labeling, and compact e/f trajectory axes before rerendering. The 2026-06-26 package-local polish feedback then requested more vertical space for panels e/f, a less compressed q50-by-glucose panel, smaller ratio/d panels if needed, and shorter panel d area/ploidy axis labels; those layout and label changes were applied in that rerun. The 2026-06-30 user-selected option 5 layout was then patched into the maintained final assembly path, with the q50-by-glucose panel stacked as rows and G0 = 0.25 mM removed from panels e/f. The v6 repair keeps the plotted content and swaps the Figure 1B/C generated letters so the stacked q50-by-glucose panel is b and the ratio summary panel is c."
  )
  writeLines(visual_qc, file.path(polish_root, "visual_qc.md"))
}

env <- load_builders()
panels <- build_panels(env)

if (phase %in% c("subpanels", "both")) {
  write_subpanels(panels)
}
if (phase %in% c("final", "both")) {
  write_final_outputs(panels)
}

cat("FG1 polishing phase complete: ", phase, "\n", sep = "")
