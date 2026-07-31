#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(jsonlite)
  library(patchwork)
  library(scales)
  library(tibble)
  library(tidyr)
})

root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v5_redraft",
  "figure_generation", "FG1_measurement_foundation", "drafting_v2"
)
out_dir <- file.path(root, "initial_subpanels", "trajectory_size_revision")
refined_dir <- file.path(root, "refined_subpanels")
worker_dir <- file.path(root, "worker_notes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(refined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(worker_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- list(
  raw_live = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "raw_live_trajectories", "wp1_raw_live_trajectories_raw_points.csv"
  ),
  live_summary = file.path("data", "report_exports", "wp1_core_starvation", "count_trajectory_summary.csv"),
  glucose_summary = file.path("data", "report_exports", "wp1_core_starvation", "glucose_trajectory_summary.csv"),
  stan_data = file.path("data", "stan_ready_data.Rds"),
  adjusted_area = file.path("data", "report_exports", "wp2_morphology_volume", "wp2_adjusted_ploidy_effects.csv"),
  nuclear_high_low = file.path("data", "report_exports", "wp3_nuclear_size", "wp3_nuclear_high_vs_low_effects.csv"),
  annotation_manifest = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_manifest.csv"),
  annotation_points = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_Results.csv"),
  object_features = file.path("data", "image_processing_runs", "run_20260324_233122", "object_train_features_90.csv"),
  objects_labels = file.path("data", "image_processing_runs", "run_20260324_233122", "objects_labels_90.tif"),
  production_frames = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_frame_counts.csv"),
  headline = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_headline_metrics.csv"),
  object_metrics = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_object_metrics.csv"),
  frame_count_metrics = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_frame_count_metrics.csv"),
  object_predictions = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_object_predictions.csv"),
  output_manifest = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_output_manifest.csv"),
  provenance = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_provenance.json"),
  production_plan = file.path("agent-dev", "major_analyses", "20260623_production_classifier_fg1_unblock", "analysis_plan.md"),
  validation_report = file.path("agent-dev", "major_analyses", "20260623_production_classifier_fg1_unblock", "validation_report.md")
)

missing <- names(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing)) {
  stop("Missing required input(s): ", paste(missing, collapse = ", "), call. = FALSE)
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
g0_desc <- c("25", "5", "1", "0.5", "0.25", "0.1", "0")
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

save_initial_png <- function(plot, filename, width, height, dpi = 300) {
  path <- file.path(out_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(path)
}

save_refined_png <- function(plot, filename, width, height, dpi = 300) {
  path <- file.path(refined_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(path)
}

make_size_ratio_vs_ploidy_ratio <- function() {
  adjusted_area <- read_tbl(input_paths$adjusted_area)
  nuclear_high_low <- read_tbl(input_paths$nuclear_high_low)

  ploidy_map <- adjusted_area %>%
    filter(.data$metric == "Area q50") %>%
    transmute(
      cellLine = as.character(.data$cellLine),
      log2_ploidy_ratio = .data$ploidy_delta
    )

  cell_summary <- adjusted_area %>%
    filter(.data$metric == "Area q50") %>%
    transmute(
      cellLine = as.character(.data$cellLine),
      measure = "Cell area q50",
      log2_ploidy_ratio = .data$ploidy_delta,
      estimate_log2 = .data$log2_ratio,
      lower_log2 = .data$lower_log2_ratio,
      upper_log2 = .data$upper_log2_ratio,
      n = .data$n_model_rows
    )

  nuclear_summary <- nuclear_high_low %>%
    group_by(.data$cellLine) %>%
    summarise(
      estimate_log2 = finite_median(.data$log2_nuclear_area_ratio_elevated_vs_baseline),
      lower_log2 = safe_quantile(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.25),
      upper_log2 = safe_quantile(.data$log2_nuclear_area_ratio_elevated_vs_baseline, 0.75),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(cellLine = as.character(.data$cellLine), measure = "Nuclear area") %>%
    left_join(ploidy_map, by = "cellLine") %>%
    select(cellLine, measure, log2_ploidy_ratio, estimate_log2, lower_log2, upper_log2, n)

  points <- bind_rows(cell_summary, nuclear_summary) %>%
    mutate(
      cellLine = factor(.data$cellLine, levels = line_order),
      measure = factor(.data$measure, levels = c("Cell area q50", "Nuclear area")),
      line_group = if_else(as.character(.data$cellLine) == "SUM-159-fuse", "SUM-159-fuse", "other lines"),
      label = recode(
        as.character(.data$cellLine),
        "SUM-159-chem" = "SUM-chem",
        "SUM-159-fuse" = "SUM-fuse",
        .default = as.character(.data$cellLine)
      ),
      label_x = .data$log2_ploidy_ratio + 0.018,
      label_y = .data$estimate_log2 + case_when(
        as.character(.data$cellLine) == "MCF10A" ~ 0.070,
        as.character(.data$cellLine) == "SNU668" ~ -0.070,
        TRUE ~ 0
      )
    )

  if (any(!is.finite(points$log2_ploidy_ratio))) {
    stop("Could not determine log2 ploidy ratio for all size-context rows.", call. = FALSE)
  }

  ggplot(points, aes(.data$log2_ploidy_ratio, .data$estimate_log2)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.24, color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.22, color = "grey68") +
    geom_errorbar(aes(ymin = .data$lower_log2, ymax = .data$upper_log2, color = .data$line_group),
                  width = 0.018, linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
    geom_point(aes(color = .data$line_group, shape = .data$line_group),
               size = 2.05, stroke = 0.72, fill = "white", na.rm = TRUE) +
    geom_text(aes(.data$label_x, .data$label_y, label = .data$label, color = .data$line_group),
              hjust = 0, size = 2.05, show.legend = FALSE, na.rm = TRUE) +
    facet_wrap(~measure, nrow = 1) +
    scale_color_manual(values = c("other lines" = "grey20", "SUM-159-fuse" = "#B2182B"), name = "line") +
    scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    scale_x_continuous(
      breaks = c(0.33, 1.0),
      labels = c("0.33", "1.00"),
      limits = c(0.24, 1.21),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = c(-1.18, 1.58), clip = "off") +
    labs(
      x = "log2 ploidy ratio, elevated / baseline",
      y = "log2 area ratio, elevated / baseline"
    ) +
    theme_fg1(base_size = 7.0) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 6.3),
      axis.text = element_text(size = 6.0),
      axis.title = element_text(size = 6.7),
      plot.margin = margin(3, 24, 3, 3)
    )
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
    select(well_idx, hours, line_id, cellLine, ploidy_metric, ploidy_abs, G0, exp_id, has_starvation, ploidy_state)
  raw %>%
    left_join(meta, by = c("well_idx", "hours")) %>%
    filter(!is.na(.data$cellLine)) %>%
    add_factors()
}

raw_glucose <- load_raw_glucose_all(glucose_summary)

live_summary_all <- live_summary %>%
  select(well_idx, hours, live_cells, dead_cells, cellLine, line_id, ploidy_metric, ploidy_abs,
         G0, ploidy_state, G0_display, G0_label)

glucose_summary_all <- glucose_summary %>%
  select(well_idx, hours, glucose_hat, glucose_hat_lo, glucose_hat_hi, any_glucose_censored,
         cellLine, line_id, ploidy_metric, ploidy_abs, G0, ploidy_state, G0_display, G0_label)

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

make_s1_readout_panel <- function(line, readout_name) {
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
    labs(x = "hours", y = NULL) +
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

make_ploidy_legend_compact <- function() {
  legend_lines <- tibble(
    ploidy_state = factor(c("low", "high"), levels = c("low", "high")),
    x = c(0.30, 0.55),
    xend = c(0.43, 0.68),
    label_x = c(0.455, 0.705),
    label = c("low", "high")
  )
  ggplot(legend_lines) +
    geom_segment(
      aes(x = .data$x, xend = .data$xend, y = 0.5, yend = 0.5, color = .data$ploidy_state),
      linewidth = 0.42
    ) +
    geom_text(aes(.data$label_x, 0.5, label = .data$label), hjust = 0, vjust = 0.5, size = 2.5) +
    annotate("text", x = 0.24, y = 0.5, label = "ploidy", hjust = 1, vjust = 0.5, size = 2.5) +
    scale_color_manual(values = ploidy_colors, guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(2, 2, 2, 2))
}

make_s1_atlas_legend_in_white_space <- function() {
  blocks <- lapply(line_order, make_s1_line_block)
  legend_cell <- cowplot::ggdraw() +
    cowplot::draw_plot(make_ploidy_legend_compact(), x = 0.08, y = 0.38, width = 0.84, height = 0.24)

  wrap_plots(blocks[[1]], blocks[[2]], blocks[[3]], blocks[[4]], blocks[[5]], legend_cell, ncol = 2) +
    plot_layout(widths = c(1, 1), heights = c(1, 1, 1))
}

write_manifest <- function(paths) {
  manifest <- tibble(
    path = paths,
    scope = "fg1_size_s1_s2_revision_worker",
    directive_ids = c("FG1-current-size-log2-ploidy", "FG1-current-s1-legend", "FG1-current-contact-sheet"),
    notes = c(
      "Exploratory size-context subpanel plotting cell-area and nuclear-area log2 ratios against log2 ploidy ratio.",
      "Local Figure S1 option with the ploidy legend placed in the bottom-right white-space cell; final recommended PNG is regenerated from this promoted option.",
      "Contact sheet collecting the two local revision images for reviewer triage."
    )
  )
  write.csv(as.data.frame(manifest), file.path(out_dir, "trajectory_size_revision_manifest.csv"), row.names = FALSE)
}

write_s2ab_explainer <- function() {
  frames <- read_tbl(input_paths$production_frames)
  headline <- read_tbl(input_paths$headline)
  obj_metrics <- read_tbl(input_paths$object_metrics)
  frame_metrics <- read_tbl(input_paths$frame_count_metrics)
  output_manifest <- read_tbl(input_paths$output_manifest)
  provenance <- jsonlite::fromJSON(input_paths$provenance, simplifyVector = TRUE)
  manifest <- read_tbl(input_paths$annotation_manifest)

  all_obj <- obj_metrics %>% filter(.data$metric_scope == "all") %>% slice(1)
  all_frame <- frame_metrics %>% filter(.data$metric_scope == "all") %>% slice(1)
  test_obj <- obj_metrics %>% filter(.data$metric_scope == "legacy_split_test") %>% slice(1)
  train_obj <- obj_metrics %>% filter(.data$metric_scope == "legacy_split_train") %>% slice(1)
  selected_rows <- output_manifest %>%
    filter(.data$role %in% c("object_predictions", "frame_counts", "headline_metrics", "selected_object_calls", "provenance_json"))

  fmt <- function(x, digits = 3) {
    formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
  }
  fmt_int <- function(x) comma(as.integer(round(as.numeric(x))))

  lines <- c(
    "# Figure S2A/B provenance explainer",
    "",
    "## Short answer",
    "",
    "- Figure S2A and Figure S2B use the same 90-frame annotation crop set, not two unrelated biological datasets.",
    "- In Figure S2A, `alive object labels` means manual alive labels after manual points were mapped onto segmented object records. These are manual labels on objects, not production-classifier calls.",
    "- In Figure S2B, `production-predicted alive objects` means the production classifier's thresholded class-1 alive calls on those same scored object records.",
    "- The user's intuition that S2A is cleaner than S2B is consistent with the exports: S2A is a point-to-object mapping QC, whereas S2B exposes that the production classifier overcalls thresholded alive objects on this 90-frame annotation set.",
    "",
    "## Direct datasets used",
    "",
    paste0("- Annotation image stack: `", input_paths$annotation_manifest, "` describes `", nrow(manifest), "` frames from `data/image_processing_runs/run_20260324_233122/annotation_stack_90.tif`."),
    paste0("- Manual points: `", input_paths$annotation_points, "` supplies point annotations by frame."),
    paste0("- Object masks and labels: `", input_paths$objects_labels, "` plus `", input_paths$object_features, "` supply kept segmented objects, per-object features, and manual alive/not-alive labels."),
    paste0("- Production classifier: `data/train/classifier_training_outputs/object_classifier.pkl`, observed sha256 `", provenance$classifier$observed_sha256, "`."),
    paste0("- Production scoring exports: `", dirname(input_paths$headline), "/` contains object predictions, frame counts, metrics, selected overlay objects, output manifest, and provenance JSON."),
    paste0("- First-pass classifier outputs enter only through `", provenance$inputs$first_pass_predictions_csv, "` to carry the legacy train/test split colors. They are not the predictions plotted in S2B."),
    "",
    "## Processing chain",
    "",
    "1. A 90-frame annotation crop stack was selected across cell-line, ploidy, starting-glucose, and time-bin strata.",
    "2. Object masks were generated for the same stack. Manual point annotations were mapped to object ids, and the kept-object feature table records `label_int = 1` for manual alive objects and `label_int = 0` for manual not-alive objects.",
    "3. The production `ObjectClassifier` rescored the kept objects. The annotation stack is stored as phase/dead/alive channels; the scoring entrypoint converts it to dead/alive/phase before classifier prediction.",
    "4. Production class 1 is mapped to alive; classes 2 and 3 are mapped to not alive for binary comparison against manual labels.",
    "5. Frame-level counts are then computed from the scored object table: manual point counts from the annotation CSV, manual alive object counts from `sum(label_int)`, and production-predicted alive counts from `sum(production_y_pred_alive)`.",
    "",
    "## What each point means",
    "",
    "Each point in S2A or S2B is one annotation-stack frame/crop.",
    "",
    paste0("- S2A x-axis, `manual alive annotation points`: number of manual point annotations in that frame. Across all frames: ", fmt_int(headline$manual_annotation_points_total[[1]]), "."),
    paste0("- S2A y-axis, `alive object labels`: number of kept segmented objects in that frame with manual alive object label. Across all frames: ", fmt_int(headline$manual_alive_objects_total[[1]]), "."),
    paste0("- S2A therefore tests point-to-object recovery. The aggregate recovery is ", fmt(headline$point_to_object_recovery[[1]], 3), "; differences reflect points that did not become one kept alive object label after object mapping/filtering."),
    "- S2B x-axis, `manual alive objects`: the same manual object-label count used on the S2A y-axis.",
    paste0("- S2B y-axis, `production-predicted alive objects`: thresholded production class-1 alive calls in the same frame. Across all frames: ", fmt_int(headline$predicted_alive_threshold_total[[1]]), " thresholded alive calls; probability-summed alive count: ", fmt(headline$predicted_alive_probability_total[[1]], 1), "."),
    "- S2B color, `legacy frame split`: train/test labels copied from the old first-pass classifier table only for context. This color is not a different production pipeline.",
    "",
    "## Evaluated pipelines and metrics",
    "",
    "- Manual object-label pipeline: annotation points plus object masks/features. This creates the ground-truth binary object labels used in both S2A and S2B.",
    paste0("- Production classifier pipeline: RandomForest classifier with sha256 `", headline$classifier_sha256[[1]], "` rescored ", fmt_int(headline$n_objects_scored[[1]]), " kept objects across ", fmt_int(headline$n_frames[[1]]), " frames."),
    paste0("- Object-level all-frame metrics: precision ", fmt(all_obj$precision[[1]], 3), ", recall ", fmt(all_obj$recall[[1]], 3), ", balanced accuracy ", fmt(all_obj$balanced_accuracy[[1]], 3), ", ROC AUC ", fmt(all_obj$roc_auc[[1]], 3), "."),
    paste0("- Object-level split context: train precision/recall ", fmt(train_obj$precision[[1]], 3), "/", fmt(train_obj$recall[[1]], 3), "; test precision/recall ", fmt(test_obj$precision[[1]], 3), "/", fmt(test_obj$recall[[1]], 3), "."),
    paste0("- Frame-level count comparison: thresholded count MAE ", fmt(all_frame$threshold_mae[[1]], 2), ", RMSE ", fmt(all_frame$threshold_rmse[[1]], 2), ", Pearson r ", fmt(all_frame$threshold_pearson_r[[1]], 3), ", Spearman r ", fmt(all_frame$threshold_spearman_r[[1]], 3), "."),
    "- First-pass classifier pipeline: not evaluated as an S2B predictor in the current v2 figure. It supplies only the legacy frame split labels.",
    "",
    "## Interpretation for the next figure decision",
    "",
    "- S2A is appropriate as a manual point-to-object mapping QC panel.",
    "- S2B is provenance-correct as a production-classifier validation panel, but it should be labeled as a caveat-bearing validation result. It should not be used as evidence that thresholded production alive counts are accurate.",
    "- If S2B remains in the figure, the legend or report should explicitly state that the production classifier has high recall and strong frame-level rank correlation but low precision and substantial thresholded alive-count overprediction on this annotation set.",
    "- If a cleaner reviewer-facing S2B is desired later, a probability-summed alive count or rank/correlation view may communicate the production validation result with less apparent threshold-count inflation.",
    "",
    "## Output files backing this explanation",
    "",
    paste0("- `", paste(selected_rows$path, collapse = "`\n- `"), "`")
  )

  writeLines(lines, file.path(worker_dir, "fg1_s2ab_provenance_explainer.md"))
}

write_coverage_note <- function(size_path, s1_path, contact_path) {
  lines <- c(
    "# FG1 size/S1/S2 revision coverage",
    "",
    "## Scope",
    "",
    "This sidecar pass addresses the non-Figure-1A items in the current `drafting_v2/user-feedback.txt`. The final redrafter promotes the S1 legend-placement option through `scripts/assemble_fg1_v2_drafts.R` while preserving the accepted atlas content.",
    "",
    "Owned paths used:",
    "",
    "- `scripts/fg1_size_s1_s2_revision_worker.R`",
    "- `initial_subpanels/trajectory_size_revision/`",
    "- `worker_notes/fg1_s2ab_provenance_explainer.md`",
    "- `worker_notes/fg1_size_s1_revision_coverage.md`",
    "- `refined_subpanels/fg1_size_revision_contact_sheet.png`",
    "",
    "## Directive coverage",
    "",
    "| feedback item | status | output / disposition |",
    "|---|---|---|",
    paste0("| Size-context plot versus log2 ploidy ratio | addressed | `", size_path, "` plots cell-area and nuclear-area log2 ratios against log2 ploidy ratio. SUM-159-fuse is highlighted and falls below the other 4N/2N lines in the cell-area facet. |"),
    paste0("| Figure S1 bottom-right legend move | addressed and promoted | `", s1_path, "` places the ploidy legend in the existing bottom-right blank atlas cell; `final_figures/recommended/Figure_S1_FG1_raw_trajectory_atlas_v2_recommended.png` is regenerated from this promoted option. |"),
    "| Figure S2A/B provenance explainer | addressed | `worker_notes/fg1_s2ab_provenance_explainer.md` details datasets, processing, evaluated pipelines, and the point meanings. |",
    paste0("| Local review triage | addressed | `", contact_path, "` combines the new size-context plot and S1 legend option for quick visual inspection. |"),
    "",
    "## S1 legend-placement inspection",
    "",
    "The legend can be moved into the bottom-right white-space cell because the current five-line atlas is arranged as a 2-column by 3-row grid with the sixth cell unused. The local option removes the separate bottom legend band and uses that sixth cell for the low/high ploidy swatches. This is more space-efficient, but still dense; final adoption should be checked at the intended print size.",
    "",
    "## Visual QC",
    "",
    "- Size-ratio plot: no figure-level title/caption; axis labels and facet labels are readable; SUM-159-fuse is directly labeled/highlighted; the panel is exploratory and not inserted into a final figure.",
    "- S1 legend option: no `S1A` label or figure-level title; legend is visible in the blank sixth cell; raw atlas remains dense but not obviously clipped.",
    "- Contact sheet: for review triage only, not a final figure candidate.",
    "",
    "## Command",
    "",
    "```bash",
    "scripts/agentRrunner.sh agent-dev/manuscript_redrafts/20260619_v5_redraft/figure_generation/FG1_measurement_foundation/drafting_v2/scripts/fg1_size_s1_s2_revision_worker.R",
    "```",
    "",
    "## Project-map decision",
    "",
    "`docs/project_map.txt` was read. No update is needed because this pass creates local review sidecar outputs inside the already documented v5 redraft FG1 drafting package and does not alter a maintained workflow, shared entrypoint, or stable data product."
  )

  writeLines(lines, file.path(worker_dir, "fg1_size_s1_revision_coverage.md"))
}

size_plot <- make_size_ratio_vs_ploidy_ratio()
s1_option <- make_s1_atlas_legend_in_white_space()

size_path <- save_initial_png(size_plot, "fg1_v2_size_area_ratio_vs_log2_ploidy_ratio.png", 5.9, 2.85, dpi = 360)
s1_path <- save_initial_png(s1_option, "fg1_v2_figure_s1_legend_in_bottom_right_option.png", 7.0, 6.85, dpi = 300)
contact <- cowplot::plot_grid(
  size_plot,
  s1_option,
  ncol = 1,
  rel_heights = c(0.33, 0.67)
)
contact_path <- save_refined_png(contact, "fg1_size_revision_contact_sheet.png", 7.0, 9.2, dpi = 300)

write_manifest(c(size_path, s1_path, contact_path))
write_s2ab_explainer()
write_coverage_note(size_path, s1_path, contact_path)

cat("Wrote FG1 size/S1/S2 revision sidecar outputs under ", root, "\n", sep = "")
