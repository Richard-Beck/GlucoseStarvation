#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(scales)
  library(tibble)
})

root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v5_redraft",
  "figure_generation", "FG1_measurement_foundation", "drafting_v2"
)
initial_dir <- file.path(root, "initial_subpanels", "measurement_provenance")
refined_dir <- file.path(root, "refined_subpanels")
worker_dir <- file.path(root, "worker_notes")
dir.create(initial_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(refined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(worker_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  coverage = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "s2_measurement_qc", "wp1_s2_measurement_qc_coverage_by_condition.csv"
  ),
  calibration_fit = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "s2_measurement_qc", "wp1_s2_glucose_calibration_fit_summary.csv"
  ),
  calibration_resid = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "s2_measurement_qc", "wp1_s2_glucose_calibration_residual_by_standard.csv"
  ),
  calibration_standards = file.path(
    "data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting",
    "panels", "s2_measurement_qc", "wp1_s2_glucose_calibration_standards_with_fit.csv"
  )
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) {
  stop("Missing required preserved S2 support input(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

read_tbl <- function(path) {
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), .name_repair = "unique")
}

save_png <- function(plot, filename, width, height, dpi = 360) {
  out <- file.path(initial_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  file.copy(out, file.path(refined_dir, filename), overwrite = TRUE)
  invisible(out)
}

line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
g0_desc <- c("25", "5", "1", "0.5", "0.25", "0.1", "0")
g0_label_levels <- paste0("G0 ", g0_desc, " mM")

format_g0 <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  out <- ifelse(
    abs(y) < 1e-8,
    "0",
    ifelse(abs(y - round(y)) < 1e-8, sprintf("%.0f", y), sprintf("%.2f", y))
  )
  out <- ifelse(out == "0", out, sub("0$", "", out))
  out <- sub("\\.$", "", out)
  out[is.na(y)] <- as.character(x)[is.na(y)]
  out
}

add_factors <- function(df) {
  df %>%
    mutate(
      cellLine = factor(as.character(cellLine), levels = line_order),
      G0_display = factor(format_g0(G0), levels = g0_desc),
      G0_label = factor(paste0("G0 ", as.character(G0_display), " mM"), levels = g0_label_levels),
      ploidy_state = factor(as.character(ploidy_state), levels = c("low", "high"))
    )
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
      legend.title = element_text(size = ggplot2::rel(0.9)),
      plot.margin = margin(3, 3, 3, 3)
    )
}

coverage <- read_tbl(paths$coverage) %>% add_factors()
cal_fit <- read_tbl(paths$calibration_fit)
cal_resid <- read_tbl(paths$calibration_resid)
standards <- read_tbl(paths$calibration_standards)

make_s2_c <- function() {
  coverage %>%
    mutate(
      condition_label = factor(condition_label, levels = rev(unique(condition_label))),
      G0_axis = factor(as.character(G0_display), levels = rev(g0_desc))
    ) %>%
    ggplot(aes(G0_axis, condition_label, fill = count_time_coverage)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = paste0(n_count_times, "c/", n_glucose_times, "g")), size = 1.55) +
    scale_fill_gradient(low = "white", high = "#386CB0", limits = c(0, 1), name = "timepoint\ncoverage") +
    labs(x = "starting glucose (mM)", y = NULL) +
    theme_fg1(base_size = 5.8) +
    theme(axis.text.y = element_text(size = 4.2), legend.position = "right")
}

make_s2_d <- function() {
  coverage %>%
    mutate(
      condition_label = factor(condition_label, levels = rev(unique(condition_label))),
      G0_axis = factor(as.character(G0_display), levels = rev(g0_desc))
    ) %>%
    ggplot(aes(G0_axis, condition_label, fill = glucose_censored_fraction)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = if_else(glucose_censored_observations > 0, as.character(glucose_censored_observations), "")), size = 1.65) +
    scale_fill_gradient(low = "white", high = "#D95F02", name = "censored\nfraction") +
    labs(x = "starting glucose (mM)", y = NULL) +
    theme_fg1(base_size = 5.8) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")
}

make_s2_e <- function() {
  ggplot(standards, aes(G, Lum, color = calibration_set)) +
    geom_point(size = 0.65, alpha = 0.85) +
    geom_line(aes(y = fitted_lum), color = "black", linewidth = 0.28) +
    facet_wrap(~ calibration_set, scales = "free_y", ncol = 2) +
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 0.1, 1, 5, 25), labels = c("0", "0.1", "1", "5", "25")) +
    labs(x = "known glucose", y = "luminescence") +
    theme_fg1(base_size = 5.9) +
    theme(legend.position = "none", strip.text = element_text(size = 4.8), axis.text = element_text(size = 4.6))
}

make_s2_f <- function() {
  ggplot(cal_resid, aes(G, residual_median, color = calibration_set)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey45") +
    geom_linerange(aes(ymin = residual_q25, ymax = residual_q75), alpha = 0.45, linewidth = 0.30) +
    geom_point(size = 0.65) +
    facet_wrap(~ calibration_set, ncol = 2) +
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 0.1, 1, 5, 25), labels = c("0", "0.1", "1", "5", "25")) +
    labs(x = "known glucose", y = "log residual") +
    theme_fg1(base_size = 5.9) +
    theme(legend.position = "none", strip.text = element_text(size = 4.8), axis.text = element_text(size = 4.6))
}

make_s2_g <- function() {
  ggplot(cal_fit, aes(reorder(calibration_set, log_rmse), log_rmse)) +
    geom_col(width = 0.62, fill = "#6B7FAF") +
    coord_flip() +
    labs(x = NULL, y = "log residual RMSE") +
    theme_fg1(base_size = 6.7)
}

outputs <- tibble(
  filename = c(
    "fg1_measurement_s2c_measurement_coverage.png",
    "fg1_measurement_s2d_glucose_censoring.png",
    "fg1_measurement_s2e_calibration_standards.png",
    "fg1_measurement_s2f_calibration_residuals.png",
    "fg1_measurement_s2g_calibration_fit_error.png"
  ),
  panel_id = c(
    "figure_s2c_measurement_coverage",
    "figure_s2d_glucose_censoring",
    "figure_s2e_calibration_standards",
    "figure_s2f_calibration_residuals",
    "figure_s2g_calibration_fit_error"
  ),
  width = c(3.35, 3.35, 3.35, 3.35, 6.80),
  height = c(2.70, 2.70, 2.30, 2.30, 1.30),
  data_source = c(paths$coverage, paths$coverage, paths$calibration_standards, paths$calibration_resid, paths$calibration_fit)
)

plots <- list(make_s2_c(), make_s2_d(), make_s2_e(), make_s2_f(), make_s2_g())
for (idx in seq_len(nrow(outputs))) {
  save_png(plots[[idx]], outputs$filename[[idx]], outputs$width[[idx]], outputs$height[[idx]])
}

notes <- c(
  "# Preserved S2C-G support worker notes",
  "",
  "Generated local drafting_v2 copies of the v4 polished Figure S2C-G support panels from the same measurement-QC and glucose-calibration exports used by the prior reviewed v4 polishing script.",
  "",
  "Outputs:",
  "",
  "- `fg1_measurement_s2c_measurement_coverage.png`: count/glucose timepoint coverage by condition.",
  "- `fg1_measurement_s2d_glucose_censoring.png`: localized glucose censoring by condition.",
  "- `fg1_measurement_s2e_calibration_standards.png`: glucose standard curves with fitted luminescence.",
  "- `fg1_measurement_s2f_calibration_residuals.png`: glucose-calibration residual summaries.",
  "- `fg1_measurement_s2g_calibration_fit_error.png`: calibration fit-error summary.",
  "",
  "Fidelity:",
  "",
  "- These panels are classified as `inherited_preserve` in `prior_code_fidelity.csv`.",
  "- The active script is a local extraction of the reviewed v4 `make_s2_c()` through `make_s2_g()` plotting logic, with output paths relocated into `drafting_v2/` and panel labels left to the final v2 assembler.",
  "- No data source, filtering rule, or visual encoding change was intended for S2C-G."
)
writeLines(notes, file.path(worker_dir, "s2_preserved_support_coverage.md"))

manifest_out <- outputs %>%
  transmute(
    output_png = file.path(initial_dir, filename),
    refined_copy = file.path(refined_dir, filename),
    directive_ids = "FG1-D08",
    data_source = data_source,
    panel_id = panel_id
  )
write.csv(manifest_out, file.path(initial_dir, "s2_preserved_support_manifest.csv"), row.names = FALSE)

cat("Wrote preserved FG1 drafting_v2 S2C-G support panels under ", root, "\n", sep = "")
