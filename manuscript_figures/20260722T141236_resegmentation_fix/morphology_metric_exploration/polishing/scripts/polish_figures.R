#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(scales)
  library(tibble)
})
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (length(args) < i + 1L) stop("Missing value after ", flag, call. = FALSE)
  args[[i + 1L]]
}
phase <- arg_value("--phase", "both")
if (!phase %in% c("subpanels", "final", "both")) {
  stop("Unsupported --phase: ", phase, call. = FALSE)
}
project_root <- normalizePath(arg_value("--project-root", "."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R/model_analysis_utils.R"))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script path", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
polish_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
subpanel_dir <- file.path(polish_root, "subpanels")
layout_dir <- file.path(polish_root, "layout")
final_dir <- file.path(polish_root, "final_images")
dir.create(subpanel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  comparison = file.path(project_root, "data/report_exports/morphology_vs_ploidy_five_line_20260729/morphology_vs_ploidy_by_model.csv"),
  pareto_assessment = file.path(project_root, "data/modelling/gpath_v1/red_a30_counts_20260722/derived/optimization/assessment.Rds"),
  metric_values = file.path(polish_root, "source_tables/ploidy_effect_metric_values.csv"),
  conditions = file.path(polish_root, "source_tables/condition_delta_loglik_summary.csv"),
  effects = file.path(polish_root, "source_tables/realized_parameter_effects_by_metric.csv"),
  per_unit = file.path(polish_root, "source_tables/per_unit_parameter_effect_coefficients.csv")
)
if (!all(file.exists(unlist(paths)))) {
  stop("One or more declared source tables are missing", call. = FALSE)
}

rel_path <- function(path) {
  sub(paste0("^", project_root, "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
}
sha256_file <- function(path) {
  strsplit(system2("sha256sum", path, stdout = TRUE)[[1]], "[[:space:]]+")[[1]][[1]]
}
line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
line_short <- c(
  "MCF10A" = "MCF10A", "MDA-MB-231" = "MDA-231", "SNU668" = "SNU668",
  "SUM-159-chem" = "SUM-chem", "SUM-159-fuse" = "SUM-fuse"
)
metric_colors <- c("Ploidy" = "#3B3B3B", "Cell area" = "#D55E00", "Nuclear area" = "#0072B2")

theme_paper <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, colour = "grey90"),
      strip.background = element_rect(fill = "grey94", colour = "grey72", linewidth = 0.25),
      strip.text = element_text(face = "bold", margin = margin(1.5, 1.5, 1.5, 1.5)),
      axis.title = element_text(size = rel(0.95)),
      axis.text = element_text(size = rel(0.86)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.height = unit(0.13, "in"),
      legend.key.width = unit(0.23, "in"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

pareto_assessment <- readRDS(paths$pareto_assessment)$fit_summary
pareto <- pareto_assessment |>
  filter(dataset_id == "all_lines", pareto_member %in% TRUE, is.finite(deviance)) |>
  transmute(
    model_id,
    pareto_rank = as.integer(pareto_rank),
    alias = vapply(model_id, build_model_alias, character(1), format = "text")
  ) |>
  arrange(pareto_rank, model_id)
if (nrow(pareto) != 5L || anyDuplicated(pareto$model_id)) {
  stop("Expected five unique all-lines Pareto models in the current count assessment")
}
pareto_levels <- pareto$alias

metric_values <- read.csv(paths$metric_values, stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(pair_state == "elevated") |>
  mutate(
    cellLine = factor(cellLine, levels = rev(line_order)),
    metric_label = recode(metric, ploidy = "Ploidy", cell_area = "Cell area", nuclear_area = "Nuclear area"),
    metric_label = factor(metric_label, levels = c("Ploidy", "Cell area", "Nuclear area"))
  )

comparison <- read.csv(paths$comparison, stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(model_id %in% pareto$model_id, metric %in% c("cell_area", "nuclear_area")) |>
  left_join(pareto[, c("model_id", "alias", "pareto_rank")], by = "model_id") |>
  mutate(
    model_alias = factor(alias, levels = rev(pareto_levels)),
    metric_label = recode(metric, cell_area = "Cell area", nuclear_area = "Nuclear area")
  )

conditions <- read.csv(paths$conditions, stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(metric == "nuclear_area") |>
  mutate(
    cell_line = factor(cell_line, levels = rev(line_order), labels = rev(unname(line_short[line_order]))),
    model_alias = factor(model_alias, levels = pareto_levels),
    glucose_label = factor(
      format(glucose_mM, trim = TRUE, scientific = FALSE),
      levels = format(sort(unique(glucose_mM)), trim = TRUE, scientific = FALSE)
    )
  )
effects <- read.csv(paths$effects, stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(metric_label %in% c("Ploidy", "Nuclear area")) |>
  mutate(
    metric_label = factor(metric_label, levels = c("Ploidy", "Nuclear area")),
    cell_line = factor(cell_line, levels = rev(line_order), labels = rev(unname(line_short[line_order]))),
    parameter_label = factor(parameter_label, levels = c("v1", "K1", "y1", "m")),
    model_alias = factor(model_alias, levels = pareto_levels),
    effect = log2_elevated_vs_baseline
  )
effect_summary <- effects |>
  group_by(metric_label, cell_line, parameter_label) |>
  summarise(median = median(effect), low = min(effect), high = max(effect), .groups = "drop")
per_unit <- read.csv(paths$per_unit, stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(metric_label %in% c("Ploidy", "Nuclear area")) |>
  mutate(
    metric_label = factor(metric_label, levels = c("Ploidy", "Nuclear area")),
    model_alias = factor(model_alias, levels = rev(pareto_levels)),
    parameter_label = factor(parameter_label, levels = c("v1", "K1", "y1", "m")),
    effect = log2_effect_per_metric_unit
  )

make_a_a <- function() {
  lim <- max(abs(metric_values$ploidy_effect_metric), na.rm = TRUE)
  ggplot(metric_values, aes(metric_label, cellLine, fill = ploidy_effect_metric)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(aes(label = sprintf("%.2f", ploidy_effect_metric)), size = 2.3) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-lim, lim), name = "Covariate\nvalue"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) unname(line_short[x])) +
    labs(x = NULL, y = NULL) +
    theme_paper(7.4) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

make_a_b <- function() {
  ggplot(comparison, aes(delta_likelihood_AIC_vs_ploidy, model_alias, colour = metric_label)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
    geom_segment(
      aes(x = 0, xend = delta_likelihood_AIC_vs_ploidy, yend = model_alias),
      colour = "grey78", linewidth = 0.45
    ) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = metric_colors[c("Cell area", "Nuclear area")], name = NULL) +
    labs(x = expression(Delta * " likelihood AIC vs ploidy"), y = NULL) +
    theme_paper(7.4) +
    theme(legend.position = "bottom")
}

make_a_c <- function() {
  lim <- max(abs(conditions$delta_log_lik), na.rm = TRUE)
  ggplot(conditions, aes(glucose_label, cell_line, fill = delta_log_lik)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    facet_grid(cols = vars(model_alias)) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-lim, lim), name = expression(Delta*log~L)
    ) +
    labs(x = "Starting glucose (mM)", y = NULL) +
    theme_paper(7.4) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

make_b_a <- function() {
  ggplot(per_unit, aes(effect, model_alias, colour = metric_label)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_line(aes(group = model_alias), colour = "grey78", linewidth = 0.45) +
    geom_point(size = 1.65) +
    facet_wrap(~parameter_label, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = metric_colors[c("Ploidy", "Nuclear area")], name = NULL) +
    labs(x = expression(log[2]("parameter change per +1 covariate unit")), y = NULL) +
    theme_paper(7.5) +
    theme(legend.position = "bottom")
}

make_b_b <- function() {
  dodge <- position_dodge(width = 0.55)
  ggplot(effects, aes(effect, cell_line, colour = metric_label)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_point(position = dodge, alpha = 0.34, size = 0.75) +
    geom_errorbar(
      data = effect_summary,
      aes(x = median, xmin = low, xmax = high, y = cell_line),
      orientation = "y", position = dodge, linewidth = 0.52, width = 0.18
    ) +
    geom_point(
      data = effect_summary, aes(x = median, y = cell_line),
      position = dodge, size = 1.65
    ) +
    facet_wrap(~parameter_label, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = metric_colors[c("Ploidy", "Nuclear area")], name = NULL) +
    labs(x = expression(log[2]("elevated / baseline parameter")), y = NULL) +
    theme_paper(7.5) +
    theme(legend.position = "bottom")
}

panel_meta <- tribble(
  ~figure, ~panel, ~width_in, ~height_in, ~filename, ~order,
  "Candidate supplement A", "a", 3.15, 2.60, "candidate_supplement_a_a.png", 1,
  "Candidate supplement A", "b", 3.75, 2.60, "candidate_supplement_a_b.png", 2,
  "Candidate supplement A", "c", 7.00, 2.85, "candidate_supplement_a_c.png", 3,
  "Candidate supplement B", "a", 7.00, 2.50, "candidate_supplement_b_a.png", 1,
  "Candidate supplement B", "b", 7.00, 3.55, "candidate_supplement_b_b.png", 2
)

build_panels <- function() {
  list(
    "Candidate supplement A::a" = make_a_a(),
    "Candidate supplement A::b" = make_a_b(),
    "Candidate supplement A::c" = make_a_c(),
    "Candidate supplement B::a" = make_b_a(),
    "Candidate supplement B::b" = make_b_b()
  )
}

save_subpanels <- function(panels) {
  dpi <- 600
  dims <- panel_meta |>
    rowwise() |>
    mutate(
      subpanel_png = rel_path(file.path(subpanel_dir, filename)),
      width_px = as.integer(width_in * dpi),
      height_px = as.integer(height_in * dpi)
    ) |>
    ungroup() |>
    select(figure, panel, subpanel_png, width_px, height_px, width_in, height_in, order)
  for (i in seq_len(nrow(panel_meta))) {
    key <- paste(panel_meta$figure[[i]], panel_meta$panel[[i]], sep = "::")
    ggsave(
      file.path(subpanel_dir, panel_meta$filename[[i]]), panels[[key]],
      width = panel_meta$width_in[[i]], height = panel_meta$height_in[[i]],
      units = "in", dpi = dpi, bg = "white", limitsize = FALSE
    )
  }
  write.csv(dims, file.path(layout_dir, "subpanel_dimensions.csv"), row.names = FALSE)
}

final_filename <- c(
  "Candidate supplement A" = "candidate_supplement_a.png",
  "Candidate supplement B" = "candidate_supplement_b.png"
)

assemble_final <- function(panels) {
  layout_path <- file.path(layout_dir, "layout_plan.csv")
  if (!file.exists(layout_path)) {
    stop("Missing layout plan: run the optimizer between subpanel and final phases")
  }
  plan <- read.csv(layout_path, stringsAsFactors = FALSE, check.names = FALSE)
  for (fig in unique(panel_meta$figure)) {
    intended_order <- panel_meta$panel[panel_meta$figure == fig]
    one <- plan[plan$figure == fig, , drop = FALSE]
    one <- one[match(intended_order, one$panel), , drop = FALSE]
    if (anyNA(one$panel)) stop("Layout plan missing one or more panels for ", fig)
    canvas <- ggdraw()
    for (i in seq_len(nrow(one))) {
      key <- paste(fig, one$panel[[i]], sep = "::")
      canvas <- canvas +
        draw_plot(
          panels[[key]],
          x = one$x_npc[[i]], y = one$y_npc[[i]],
          width = one$width_npc[[i]], height = one$height_npc[[i]]
        ) +
        draw_label(
          one$panel[[i]],
          x = one$x_npc[[i]] + 0.002,
          y = one$y_npc[[i]] + one$height_npc[[i]] - 0.002,
          hjust = 0, vjust = 1, fontface = "bold", size = 9
        )
    }
    ggsave(
      file.path(final_dir, final_filename[[fig]]), canvas,
      width = unique(one$layout_width_in)[[1]],
      height = unique(one$layout_height_in)[[1]],
      units = "in", dpi = 600, bg = "white", limitsize = FALSE
    )
  }
}

write_package_records <- function() {
  data_by_panel <- c(
    "Candidate supplement A::a" = rel_path(paths$metric_values),
    "Candidate supplement A::b" = paste(rel_path(paths$comparison), rel_path(paths$pareto_assessment), sep = ";"),
    "Candidate supplement A::c" = rel_path(paths$conditions),
    "Candidate supplement B::a" = rel_path(paths$per_unit),
    "Candidate supplement B::b" = rel_path(paths$effects)
  )
  provenance <- panel_meta |>
    mutate(
      generator = rel_path(script_path),
      command = paste("scripts/agentRrunner.sh", rel_path(script_path), "--phase final"),
      data_inputs = unname(data_by_panel[paste(figure, panel, sep = "::")]),
      subpanel_image = rel_path(file.path(subpanel_dir, filename)),
      layout_plan = rel_path(file.path(layout_dir, "layout_plan.csv")),
      output_image = rel_path(file.path(final_dir, unname(final_filename[figure]))),
      notes = "Exploratory all-five-line MAP comparison; no posterior or model-selection uncertainty."
    ) |>
    select(figure, panel, generator, command, data_inputs, subpanel_image, layout_plan, output_image, notes)
  write.csv(provenance, file.path(polish_root, "provenance.csv"), row.names = FALSE)

  output_paths <- c(
    file.path(subpanel_dir, panel_meta$filename),
    file.path(final_dir, unname(final_filename)),
    file.path(layout_dir, "subpanel_dimensions.csv"),
    file.path(layout_dir, "layout_plan.csv"),
    file.path(layout_dir, "layout_report.md")
  )
  output_paths <- output_paths[file.exists(output_paths)]
  manifest <- tibble(
    path = vapply(output_paths, rel_path, character(1)),
    size_bytes = file.info(output_paths)$size,
    sha256 = vapply(output_paths, sha256_file, character(1))
  )
  write.csv(manifest, file.path(polish_root, "manifest.csv"), row.names = FALSE)

  finals <- file.path(final_dir, unname(final_filename))
  rebuild <- tibble(
    figure = names(final_filename),
    stage = "polishing",
    output = vapply(finals, rel_path, character(1)),
    output_sha256 = vapply(finals, sha256_file, character(1)),
    generator = rel_path(script_path),
    command = paste("scripts/agentRrunner.sh", rel_path(script_path), "--phase final"),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write.table(
    rebuild, file.path(polish_root, "figure_rebuild_manifest.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

panels <- build_panels()
if (phase %in% c("subpanels", "both")) save_subpanels(panels)
if (phase %in% c("final", "both")) {
  assemble_final(panels)
  write_package_records()
}
cat("Completed phase ", phase, " in ", polish_root, "\n", sep = "")
