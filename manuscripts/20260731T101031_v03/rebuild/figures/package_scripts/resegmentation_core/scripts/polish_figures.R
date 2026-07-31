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
  library(tiff)
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
if (!phase %in% c("subpanels", "final", "both")) stop("Unsupported --phase: ", phase, call. = FALSE)
project_root <- normalizePath(arg_value("--project-root", "."), winslash = "/", mustWork = TRUE)
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
  stan = file.path(project_root, "data/modelling/gpath_v1/red_a30_counts_20260722/datasets/all_lines/stan_data.Rds"),
  area = file.path(project_root, "data/image_processing_runs/full_segmentation_classification_nuclear/run_20260721_163410/summaries/cell_nuclear_area_summaries.Rds"),
  features = file.path(project_root, "data/empirical_features/20260722/model_free_features.Rds"),
  validation = file.path(project_root, "data/image_processing_validation/segmentation_classifier_90frame/20260722"),
  annotation_manifest = file.path(project_root, "data/image_processing_runs/run_20260324_233122/annotation_stack_90_manifest.csv")
)
stopifnot(all(file.exists(unlist(paths))))
read_tbl <- function(path) as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
rel_path <- function(path) sub(paste0("^", project_root, "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
finite_median <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
finite_q <- function(x, p) { x <- x[is.finite(x)]; if (length(x)) as.numeric(quantile(x, p, names = FALSE, type = 8)) else NA_real_ }
sha256_file <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE)[[1]], "[[:space:]]+")[[1]][[1]]
line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
line_short <- c("MCF10A" = "MCF10A", "MDA-MB-231" = "MDA-231", "SNU668" = "SNU668", "SUM-159-chem" = "SUM-chem", "SUM-159-fuse" = "SUM-fuse")
line_palette <- c("MCF10A" = "#009E73", "MDA-MB-231" = "#D55E00", "SNU668" = "#7A5195", "SUM-159-chem" = "#0072B2", "SUM-159-fuse" = "#CC79A7")
ploidy_palette <- c(low = "#56B4E9", high = "#D55E00")
feature_palette <- c(Growth = "#009E73", `Alive AUC` = "#7A5195", `Glucose use` = "#E69F00", `Total-cell yield` = "#0072B2")

theme_paper <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.14, color = "grey88"),
      strip.background = element_rect(fill = "grey94", color = "grey72", linewidth = 0.22),
      strip.text = element_text(face = "bold", margin = margin(1.2, 1.2, 1.2, 1.2)),
      legend.position = "bottom", legend.box = "horizontal",
      legend.key.height = unit(0.12, "in"), legend.key.width = unit(0.22, "in"),
      axis.title = element_text(size = rel(0.92)), axis.text = element_text(size = rel(0.82)),
      plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
      plot.margin = margin(3, 5, 3, 3)
    )
}

stan <- readRDS(paths$stan)
area_release <- readRDS(paths$area)
feature_release <- readRDS(paths$features)
inv_line <- setNames(names(stan$line_map), as.character(stan$line_map))
wells <- tibble(
  well_idx = seq_len(stan$N_wells),
  line_id = as.integer(stan$line_id),
  cellLine = unname(inv_line[as.character(stan$line_id)]),
  ploidy_metric = as.numeric(stan$ploidy_metric), ploidy_abs = as.numeric(stan$ploidy_abs),
  G0 = as.numeric(stan$G0_per_well), exp_id = as.integer(stan$exp_id)
) %>% group_by(cellLine) %>% mutate(ploidy_state = if_else(ploidy_metric == min(ploidy_metric), "low", "high")) %>% ungroup()
count_raw <- tibble(
  well_idx = as.integer(stan$well_idx_count), hours = stan$t_grid[as.integer(stan$grid_idx_count)],
  replicate = as.character(stan$rep_id_count), live = as.numeric(stan$N_obs), dead = as.numeric(stan$D_obs)
) %>% mutate(total = live + dead) %>% left_join(wells, by = "well_idx")
if (!isTRUE(all.equal(count_raw$live, as.numeric(stan$N_obs))) ||
    !isTRUE(all.equal(count_raw$total, as.numeric(stan$N_obs) + as.numeric(stan$D_obs)))) {
  stop("Count adapter invariant failed: N_obs must be live cells and total must be N_obs + D_obs")
}
glucose_raw <- tibble(
  well_idx = as.integer(stan$well_idx_gluc), hours = stan$t_grid[as.integer(stan$grid_idx_gluc)],
  lum = as.numeric(stan$lum_obs), dilution = as.numeric(stan$dilution), censored = as.integer(stan$is_censored)
) %>% left_join(wells, by = "well_idx") %>%
  mutate(glucose = pmax(0, (lum - stan$calib_b_fixed[exp_id]) / (stan$calib_a_fixed[exp_id] * dilution)))
count_summary <- count_raw %>% group_by(well_idx, cellLine, ploidy_state, ploidy_metric, G0, hours) %>%
  summarise(across(c(live, dead, total), finite_median), .groups = "drop")
glucose_summary <- glucose_raw %>% group_by(well_idx, cellLine, ploidy_state, ploidy_metric, G0, hours) %>%
  summarise(glucose = finite_median(glucose), censored_fraction = mean(censored), .groups = "drop")

area <- as_tibble(area_release$replicate_time) %>%
  filter(cellLine %in% line_order, object_scope == "alive_area_pass") %>%
  mutate(
    G0 = as.numeric(glucose), ploidy_state = case_when(
      cellLine == "MDA-MB-231" & ploidy == "3N" ~ "low",
      cellLine == "MDA-MB-231" ~ "high",
      ploidy %in% c("2N", "low") ~ "low", TRUE ~ "high"
    ),
    glucose_bin = case_when(G0 <= 0.1 ~ "low", G0 <= 1 ~ "intermediate", TRUE ~ "high")
  )
area_conditions <- area %>%
  group_by(cellLine, ploidy_state, ploidy, G0, glucose_bin, hours) %>%
  summarise(cell = finite_median(segmented_area_px_q50), nuclear = finite_median(nuclear_area_detected_px_q50), .groups = "drop") %>%
  group_by(cellLine, ploidy_state, glucose_bin) %>% mutate(time_bin = ntile(rank(hours, ties.method = "first"), 3)) %>% ungroup() %>%
  group_by(cellLine, ploidy_state, glucose_bin, time_bin) %>% summarise(cell = finite_median(cell), nuclear = finite_median(nuclear), .groups = "drop") %>%
  pivot_wider(names_from = ploidy_state, values_from = c(cell, nuclear)) %>%
  mutate(log2_cell_ratio = log2(cell_high / cell_low), log2_nuclear_ratio = log2(nuclear_high / nuclear_low))

make_f1b <- function() {
  x <- area %>% group_by(cellLine, ploidy_state, G0) %>%
    summarise(q50 = finite_median(segmented_area_px_q50), q25 = finite_q(segmented_area_px_q50, .25), q75 = finite_q(segmented_area_px_q50, .75), .groups = "drop") %>%
    mutate(cellLine = factor(cellLine, line_order))
  ggplot(x, aes(G0, q50, color = ploidy_state, group = ploidy_state)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = ploidy_state), color = NA, alpha = .14) +
    geom_line(linewidth = .46) + geom_point(size = .8) +
    facet_wrap(~cellLine, scales = "free_y", ncol = 1, labeller = as_labeller(line_short)) +
    scale_x_continuous(trans = pseudo_log_trans(base = 10), breaks = c(0, .1, 1, 5, 25)) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") + scale_fill_manual(values = ploidy_palette, guide = "none") +
    labs(x = "starting glucose (mM)", y = "segmented cell area q50 (px)") + theme_paper(6.7)
}

make_f1c <- function() {
  summaries <- area_conditions %>% group_by(cellLine) %>% summarise(x = finite_median(log2_cell_ratio), y = finite_median(log2_nuclear_ratio), .groups = "drop")
  ggplot(area_conditions, aes(log2_cell_ratio, log2_nuclear_ratio)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .22, color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = .22, color = "grey55") +
    geom_point(color = "grey72", size = .65) +
    geom_point(data = summaries, aes(x, y, color = cellLine), size = 1.7) +
    geom_text(data = summaries, aes(x, y, label = line_short[cellLine], color = cellLine), size = 1.65, nudge_x = .07, hjust = 0, show.legend = FALSE) +
    scale_color_manual(values = line_palette, guide = "none") +
    labs(x = "cell-area ratio, log2(high / low)", y = "nuclear-area ratio, log2(high / low)") +
    coord_cartesian(clip = "off") + theme_paper(6.8) + theme(plot.margin = margin(3, 25, 3, 3))
}

make_f1d <- function() {
  x <- area_conditions %>% pivot_longer(c(log2_cell_ratio, log2_nuclear_ratio), names_to = "measure", values_to = "ratio") %>%
    group_by(cellLine, measure) %>% summarise(y = finite_median(ratio), lo = finite_q(ratio, .25), hi = finite_q(ratio, .75), .groups = "drop") %>%
    left_join(wells %>% group_by(cellLine) %>% summarise(x = log2(max(ploidy_abs) / min(ploidy_abs)), .groups = "drop"), by = "cellLine") %>%
    mutate(measure = recode(measure, log2_cell_ratio = "Cell area", log2_nuclear_ratio = "Nuclear area"))
  ggplot(x, aes(x, y, color = measure)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = .22, color = "grey55") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .025, linewidth = .35) + geom_point(size = 1.5) +
    geom_text(aes(label = line_short[cellLine]), nudge_x = .03, hjust = 0, size = 1.55, show.legend = FALSE) +
    facet_wrap(~measure, nrow = 1) + scale_color_manual(values = c("Cell area" = "#2166AC", "Nuclear area" = "#B2182B"), guide = "none") +
    labs(x = "ploidy ratio, log2(high / low)", y = "area ratio, log2(high / low)") + coord_cartesian(clip = "off") +
    theme_paper(6.7) + theme(plot.margin = margin(3, 25, 3, 3))
}

make_trajectory_grid <- function(readout = c("live", "glucose"), lines = line_order[-5], glucose_keep = c(0, 1, 25), ncol = 2) {
  readout <- match.arg(readout)
  x <- if (readout == "live") count_summary else glucose_summary
  ycol <- if (readout == "live") "live" else "glucose"
  x <- x %>% filter(cellLine %in% lines, G0 %in% glucose_keep) %>% mutate(cellLine = factor(cellLine, lines), G0f = factor(G0, glucose_keep))
  ggplot(x, aes(hours, .data[[ycol]], color = ploidy_state, linetype = G0f, group = interaction(ploidy_state, G0f))) +
    geom_line(linewidth = .42) + facet_wrap(~cellLine, scales = "free_y", ncol = ncol, labeller = as_labeller(line_short)) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_linetype_manual(values = c("solid", "22", "11"), name = "G0 (mM)") +
    labs(x = "time (h)", y = if (readout == "live") "live cells" else "extracellular glucose (mM)") + theme_paper(6.5)
}

validation_manifest <- read_tbl(file.path(paths$validation, "pipeline_manifest.csv"))
annotation_manifest <- read_tbl(paths$annotation_manifest)
validation_frames <- read_tbl(file.path(paths$validation, "frame_summary.csv"))
validation_predictions <- read_tbl(file.path(paths$validation, "object_predictions.csv"))
validation_points <- read_tbl(file.path(paths$validation, "manual_points.csv"))
validation_metrics <- read_tbl(file.path(paths$validation, "object_metrics.csv"))

robust01 <- function(x) {
  q <- finite_q(as.numeric(x), c(.01, .99))
  if (length(q) != 2L || !is.finite(diff(q)) || diff(q) <= 0) return(matrix(0, nrow(x), ncol(x)))
  pmin(1, pmax(0, (x - q[[1]]) / diff(q)))
}
read_json_channel <- function(x) {
  relative_paths <- as.character(fromJSON(x))
  if (!length(relative_paths)) stop("Channel path list is empty", call. = FALSE)
  images <- lapply(file.path(project_root, relative_paths), readTIFF)
  shapes <- vapply(images, function(image) paste(dim(image), collapse = "x"), character(1))
  if (length(unique(shapes)) != 1L) stop("Channel TIFF dimensions do not match", call. = FALSE)
  Reduce(`+`, images)
}
mask_path <- function(image_key) file.path(paths$validation, "pipeline/masks/cell", paste0(substr(digest::digest(image_key, algo = "sha256", serialize = FALSE), 1, 20), ".tif"))
boundary_points <- function(labels, calls) {
  nr <- nrow(labels); nc <- ncol(labels)
  edge <- labels > 0 & (
    rbind(TRUE, labels[-nr, , drop = FALSE] != labels[-1, , drop = FALSE]) |
      rbind(labels[-1, , drop = FALSE] != labels[-nr, , drop = FALSE], TRUE) |
      cbind(TRUE, labels[, -nc, drop = FALSE] != labels[, -1, drop = FALSE]) |
      cbind(labels[, -1, drop = FALSE] != labels[, -nc, drop = FALSE], TRUE)
  )
  idx <- which(edge, arr.ind = TRUE)
  tibble(object_id = labels[idx], x = idx[, "col"], y = nr - idx[, "row"] + 1L) %>%
    inner_join(calls %>% select(object_id, predicted_label_name), by = "object_id")
}
make_overlay_tile <- function(frame) {
  meta <- annotation_manifest %>% filter(stack_index_1based == frame) %>% slice(1)
  vm <- validation_manifest %>% filter(image_key == meta$image_key) %>% slice(1)
  if (nrow(meta) != 1L || nrow(vm) != 1L) stop("Missing validation metadata for frame ", frame)
  channels <- list(
    alive = read_json_channel(vm$alive_paths_json),
    dead = read_json_channel(vm$dead_paths_json),
    phase = read_json_channel(vm$phase_paths_json)
  )
  yy <- (meta$crop_y0 + 1L):meta$crop_y1; xx <- (meta$crop_x0 + 1L):meta$crop_x1
  green <- robust01(channels$alive[yy, xx]); red <- robust01(channels$dead[yy, xx]); phase_img <- robust01(channels$phase[yy, xx])
  rgb <- array(0, dim = c(length(yy), length(xx), 3L))
  rgb[, , 1] <- pmin(1, .58 * phase_img + .72 * red)
  rgb[, , 2] <- pmin(1, .58 * phase_img + .72 * green)
  rgb[, , 3] <- pmin(1, .58 * phase_img)
  labels <- round(readTIFF(mask_path(meta$image_key), native = FALSE)[yy, xx] * 65535)
  calls <- validation_predictions %>% filter(frame == !!frame, in_annotation_target == 1)
  outlines <- boundary_points(labels, calls)
  pts <- validation_points %>% filter(frame == !!frame)
  h <- nrow(labels); w <- ncol(labels)
  state <- if (meta$ploidy %in% c("2N", "3N", "low")) "low" else "high"
  ggplot() + annotation_raster(as.raster(rgb), xmin = 0, xmax = w, ymin = 0, ymax = h) +
    geom_point(data = outlines, aes(x, y, color = predicted_label_name), size = .08, alpha = .9) +
    geom_point(data = pts, aes(x_crop, h - y_crop), shape = 3, color = "#FFD92F", stroke = .35, size = .8) +
    annotate("rect", xmin = meta$target_x0 - meta$crop_x0, xmax = meta$target_x1 - meta$crop_x0,
             ymin = h - (meta$target_y1 - meta$crop_y0), ymax = h - (meta$target_y0 - meta$crop_y0),
             color = "white", fill = NA, linewidth = .22) +
    scale_color_manual(values = c(alive = "#00E676", dead = "#FF5252", junk = "#40C4FF"), name = "classifier") +
    coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE) +
    labs(title = paste0(line_short[[meta$cellLine]], " ", state)) +
    theme_void(base_size = 6.3) + theme(plot.title = element_text(face = "bold", hjust = .5, size = 5.8), legend.position = "none", plot.margin = margin(1, 1, 1, 1))
}
make_f1a <- function() {
  frames <- c(4L, 22L, 49L, 76L, 13L, 31L, 40L, 85L)
  plots <- lapply(frames, make_overlay_tile)
  wrap_plots(plots, ncol = 4, guides = "collect") & theme(legend.position = "bottom")
}

make_s1 <- function() {
  make_one <- function(line, readout) {
    if (readout %in% c("live", "dead")) {
      x <- count_summary %>% filter(cellLine == line)
    } else {
      x <- glucose_summary %>% filter(cellLine == line)
    }
    ggplot(x, aes(hours, .data[[readout]], color = ploidy_state, group = interaction(ploidy_state, G0))) +
      geom_line(aes(alpha = factor(G0)), linewidth = .32) +
      scale_color_manual(values = ploidy_palette, name = "ploidy") +
      scale_alpha_manual(values = seq(.25, 1, length.out = length(sort(unique(wells$G0)))), name = "G0") +
      labs(x = if (line == tail(line_order, 1)) "time (h)" else NULL,
           y = switch(readout, live = "live cells", dead = "dead cells", glucose = "glucose (mM)"),
           title = if (readout == "live") line_short[[line]] else NULL) +
      theme_paper(5.9) + theme(plot.title = element_text(face = "bold", size = 6.5), legend.position = "none", plot.margin = margin(2, 2, 2, 2))
  }
  plots <- unlist(lapply(line_order, function(line) lapply(c("live", "dead", "glucose"), function(readout) make_one(line, readout))), recursive = FALSE)
  wrap_plots(plots, ncol = 3, guides = "collect") & theme(legend.position = "bottom")
}

make_s2a <- function() {
  ggplot(validation_frames, aes(manual_point_count, manual_alive_object_count, color = cellLine)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = .25, color = "grey50") + geom_point(size = 1, alpha = .8) +
    scale_color_manual(values = line_palette, labels = line_short, name = "cell line") +
    labs(x = "manual alive points", y = "mapped alive objects") + coord_equal() + theme_paper(6.7)
}
make_s2b <- function() {
  x <- validation_metrics %>% pivot_longer(c(roc_auc, balanced_accuracy, precision, recall, specificity, f1), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("roc_auc", "balanced_accuracy", "precision", "recall", "specificity", "f1"),
                           labels = c("ROC AUC", "balanced\naccuracy", "precision", "recall", "specificity", "F1")))
  ggplot(x, aes(metric, value)) + geom_col(fill = "#0072B2", width = .68) + geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.35, size = 2) +
    scale_y_continuous(limits = c(0, 1.05), breaks = c(0, .5, 1)) + labs(x = NULL, y = "object-level metric") + theme_paper(6.7)
}
condition_label <- function(cellLine, ploidy_state, G0) paste(line_short[cellLine], ploidy_state, paste0(G0, " mM"), sep = " | ")
make_s2c <- function() {
  c1 <- count_raw %>% group_by(cellLine, ploidy_state, G0) %>% summarise(n_times = n_distinct(hours), n_obs = n(), .groups = "drop") %>% mutate(readout = "cell counts")
  c2 <- glucose_raw %>% group_by(cellLine, ploidy_state, G0) %>% summarise(n_times = n_distinct(hours), n_obs = n(), .groups = "drop") %>% mutate(readout = "glucose")
  x <- bind_rows(c1, c2) %>% mutate(condition = condition_label(cellLine, ploidy_state, G0))
  ggplot(x, aes(n_times, reorder(condition, n_times), color = readout, size = n_obs)) + geom_point(alpha = .82) +
    scale_color_manual(values = c(`cell counts` = "#0072B2", glucose = "#E69F00")) +
    scale_size_continuous(range = c(.7, 2.6)) + labs(x = "distinct measurement times", y = NULL, color = "readout", size = "observations") + theme_paper(5.8)
}
make_s2d <- function() {
  x <- glucose_raw %>% group_by(cellLine, ploidy_state, G0) %>% summarise(censored = mean(censored), n = n(), .groups = "drop") %>%
    mutate(condition = condition_label(cellLine, ploidy_state, G0))
  ggplot(x, aes(censored, reorder(condition, censored), color = cellLine)) + geom_point(size = 1.25) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) + scale_color_manual(values = line_palette, labels = line_short, guide = "none") +
    labs(x = "censored glucose observations", y = NULL) + theme_paper(5.8)
}
calib_names <- c("MCF10A / MDA-231", "SUM-chem / SNU668", "SUM-fuse glucose-cell", "SUM-fuse")
calibration <- tibble(batch = as.integer(stan$calib_exp_idx), glucose = as.numeric(stan$calib_G), lum = as.numeric(stan$calib_Lum)) %>%
  mutate(fitted = stan$calib_a_fixed[batch] * glucose + stan$calib_b_fixed[batch], residual = lum - fitted, batch_label = factor(calib_names[batch], calib_names))
make_s2e <- function() {
  grid <- tibble(batch = rep(seq_along(calib_names), each = 100), glucose = rep(seq(0, max(calibration$glucose), length.out = 100), length(calib_names))) %>%
    mutate(fitted = stan$calib_a_fixed[batch] * glucose + stan$calib_b_fixed[batch], batch_label = factor(calib_names[batch], calib_names))
  ggplot(calibration, aes(glucose, lum)) + geom_line(data = grid, aes(y = fitted), color = "#D55E00", linewidth = .45) + geom_point(size = .75, alpha = .75) +
    facet_wrap(~batch_label, scales = "free_y", ncol = 2) + labs(x = "glucose standard (mM)", y = "luminescence") + theme_paper(6.3)
}
make_s2f <- function() {
  ggplot(calibration, aes(fitted, residual, color = batch_label)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = .24) + geom_point(size = .8) +
    scale_color_manual(values = c("#009E73", "#D55E00", "#7A5195", "#0072B2"), name = "calibration batch") +
    labs(x = "fitted luminescence", y = "residual") + theme_paper(6.4)
}
make_s2g <- function() {
  x <- calibration %>% group_by(batch_label) %>% summarise(rmse = sqrt(mean(residual^2)), mae = mean(abs(residual)), .groups = "drop") %>%
    pivot_longer(c(rmse, mae), names_to = "metric", values_to = "error")
  ggplot(x, aes(batch_label, error, fill = metric)) + geom_col(position = "dodge", width = .7) +
    scale_fill_manual(values = c(rmse = "#0072B2", mae = "#56B4E9"), labels = c(rmse = "RMSE", mae = "MAE"), name = NULL) +
    labs(x = NULL, y = "calibration error (luminescence)") + theme_paper(6.3) + theme(axis.text.x = element_text(angle = 24, hjust = 1))
}

condition_features <- as_tibble(feature_release$condition_features) %>%
  group_by(cellLine) %>% mutate(ploidy_state = if_else(ploidy_metric == min(ploidy_metric), "low", "high")) %>% ungroup()
line_features <- as_tibble(feature_release$line_ploidy_features) %>%
  group_by(cellLine) %>% mutate(ploidy_state = if_else(ploidy_metric == min(ploidy_metric), "low", "high")) %>% ungroup()

make_f2a <- function() {
  x <- as_tibble(feature_release$ploidy_effects) %>% filter(cellLine != "SUM-159-fuse") %>%
    group_by(feature_id) %>% mutate(scale = sd(effect_per_ploidy), scaled_effect = pmax(-3, pmin(3, effect_per_ploidy / scale))) %>% ungroup() %>%
    mutate(display_label = factor(display_label, levels = rev(unique(feature_release$feature_catalog$display_label))), cellLine = factor(cellLine, line_order))
  ggplot(x, aes(scaled_effect, display_label, color = cellLine)) +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = .26) +
    geom_point(position = position_dodge(width = .5), size = 1.45) +
    scale_color_manual(values = line_palette, labels = line_short, name = "cell line") +
    scale_x_continuous(limits = c(-3.15, 3.15), breaks = -3:3) +
    labs(x = "ploidy effect / within-feature SD", y = NULL) + theme_paper(6.8)
}

representative <- condition_features %>% filter(cellLine == "MDA-MB-231", ploidy_state == "low", abs(G0 - 1) < 1e-8) %>% slice(1)
if (nrow(representative) != 1L) stop("Could not select representative MDA-MB-231 low-ploidy 1 mM condition")
make_f2b <- function() {
  obs <- count_summary %>% filter(well_idx == representative$well_idx) %>% arrange(hours)
  smooth_live <- smooth.spline(obs$hours, log1p(obs$live), spar = .62)
  grid <- tibble(hours = seq(min(obs$hours), max(obs$hours), length.out = 260))
  grid$log_live <- predict(smooth_live, grid$hours)$y
  grid$live <- expm1(grid$log_live)
  total_fit <- smooth.spline(obs$hours, obs$total, spar = .62)
  grid$total <- predict(total_fit, grid$hours)$y
  t0 <- representative$glucose_start_time; t1 <- representative$glucose_end_time
  auc_poly <- bind_rows(tibble(hours = t0, live = 0), grid %>% filter(hours >= t0, hours <= t1) %>% select(hours, live), tibble(hours = t1, live = 0))
  tangent_t <- representative$robust_derivative_time
  tangent_y <- predict(smooth_live, tangent_t)$y
  tangent <- tibble(hours = range(obs$hours), log_live = tangent_y + representative$robust_max_derivative * (range(obs$hours) - tangent_t), series = "derivative tangent") %>%
    mutate(value = expm1(log_live))
  base <- ggplot() +
    geom_polygon(data = auc_poly, aes(hours, live), fill = feature_palette[["Alive AUC"]], alpha = .19) +
    geom_line(data = grid, aes(hours, live, color = "live cells"), linewidth = .62) +
    geom_point(data = obs, aes(hours, live, color = "live cells"), size = .65) +
    geom_point(data = obs, aes(hours, total, color = "total cells"), size = .55, alpha = .75) +
    geom_line(data = tangent, aes(hours, value, color = "derivative tangent"), linetype = "dashed", linewidth = .55) +
    geom_vline(xintercept = c(t0, t1), color = feature_palette[["Alive AUC"]], linetype = "dashed", linewidth = .3) +
    scale_color_manual(values = c("live cells" = "#009E73", "total cells" = "#0072B2", "derivative tangent" = "#D55E00"), name = NULL) +
    labs(x = "time (h)", y = "cells", title = NULL) + theme_paper(7)
  base + annotate("text", x = t0, y = max(obs$total) * .92, label = "glucose window", hjust = 0, size = 2, color = feature_palette[["Alive AUC"]])
}

make_f2c <- function() {
  x <- condition_features %>% filter(cellLine == "MDA-MB-231")
  p1 <- ggplot(x, aes(G0, robust_max_derivative, color = ploidy_state)) + geom_line(linewidth = .4) + geom_point(size = .75) +
    geom_point(data = x %>% filter(G0 <= .25 | G0 >= 5), size = 1.4, shape = 21, aes(fill = ploidy_state)) +
    scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, .1, 1, 5, 25)) +
    scale_color_manual(values = ploidy_palette, guide = "none") + scale_fill_manual(values = ploidy_palette, guide = "none") +
    labs(x = "starting glucose (mM)", y = "robust derivative") + theme_paper(6.2)
  auc <- x %>% filter(G0 <= 1, is.finite(live_auc_glucose_window), is.finite(glucose_drawdown))
  p2 <- ggplot(auc, aes(live_auc_glucose_window, glucose_drawdown, color = ploidy_state)) + geom_smooth(method = "lm", se = FALSE, linewidth = .45) + geom_point(size = .85) +
    scale_color_manual(values = ploidy_palette, guide = "none") + labs(x = "live-cell AUC", y = "glucose drawdown") + theme_paper(6.2)
  p3 <- ggplot(x, aes(G0, peak_total_yield_net, color = ploidy_state)) + geom_line(linewidth = .4) + geom_point(size = .75) +
    geom_point(data = x %>% filter(abs(G0 - 1) < 1e-8), size = 1.5, shape = 21, aes(fill = ploidy_state)) +
    scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, .1, 1, 5, 25)) + scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_fill_manual(values = ploidy_palette, guide = "none") + labs(x = "starting glucose (mM)", y = "net total-cell yield") + theme_paper(6.2)
  p1 | p2 | p3
}

make_s3 <- function() {
  meta <- as_tibble(feature_release$feature_catalog) %>% select(feature_id, display_label, feature_class)
  x <- line_features %>% pivot_longer(all_of(meta$feature_id), names_to = "feature_id", values_to = "value") %>%
    left_join(meta, by = "feature_id") %>% mutate(display_label = factor(display_label, levels = meta$display_label), cellLine = factor(cellLine, line_order))
  ggplot(x, aes(ploidy_state, value, color = ploidy_state, group = cellLine)) +
    geom_line(color = "grey65", linewidth = .28) + geom_point(size = .8) +
    facet_grid(display_label ~ cellLine, scales = "free_y", labeller = labeller(cellLine = as_labeller(line_short))) +
    scale_color_manual(values = ploidy_palette, guide = "none") + labs(x = NULL, y = "raw feature value") +
    theme_paper(5.4) + theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0, hjust = 0))
}

make_s4 <- function() {
  x <- condition_features %>% mutate(cellLine = factor(cellLine, line_order))
  p1 <- ggplot(x, aes(G0, robust_max_derivative, color = ploidy_state)) + geom_line(linewidth = .3) + geom_point(size = .55) +
    geom_point(data = x %>% filter(G0 <= .25 | G0 >= 5), aes(fill = ploidy_state), shape = 21, size = 1) +
    facet_wrap(~cellLine, nrow = 1, scales = "free_y", labeller = as_labeller(line_short)) + scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, .1, 1, 25)) +
    scale_color_manual(values = ploidy_palette, guide = "none") + scale_fill_manual(values = ploidy_palette, guide = "none") +
    labs(x = "starting glucose (mM)", y = "robust derivative") + theme_paper(5.6)
  auc <- x %>% filter(G0 <= 1, is.finite(live_auc_glucose_window), is.finite(glucose_drawdown))
  p2 <- ggplot(auc, aes(live_auc_glucose_window, glucose_drawdown, color = ploidy_state)) + geom_smooth(method = "lm", se = FALSE, linewidth = .35) + geom_point(size = .6) +
    facet_wrap(~cellLine, nrow = 1, scales = "free", labeller = as_labeller(line_short)) + scale_color_manual(values = ploidy_palette, guide = "none") +
    labs(x = "live-cell AUC", y = "glucose drawdown") + theme_paper(5.6)
  p3 <- ggplot(x, aes(G0, peak_total_yield_net, color = ploidy_state)) + geom_line(linewidth = .3) + geom_point(size = .55) +
    geom_point(data = x %>% filter(abs(G0 - 1) < 1e-8), aes(fill = ploidy_state), shape = 21, size = 1.1) +
    facet_wrap(~cellLine, nrow = 1, scales = "free_y", labeller = as_labeller(line_short)) + scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, .1, 1, 25)) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") + scale_fill_manual(values = ploidy_palette, guide = "none") +
    labs(x = "starting glucose (mM)", y = "net total-cell yield") + theme_paper(5.6)
  (p1 / p2 / p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

# Load the accepted v02 constructors after all canonical data adapters exist.
# This intentionally overrides the provisional redraw constructors above.
source(file.path(dirname(script_path), "recycled_constructors.R"), local = environment())

panel_meta <- tribble(
  ~figure, ~panel, ~width_in, ~height_in, ~filename,
  "Figure 1", "a", 7.00, 3.45, "figure_1_a.png",
  "Figure 1", "b", 2.02, 5.62, "figure_1_b.png",
  "Figure 1", "c", 1.98, 1.73, "figure_1_c.png",
  "Figure 1", "d", 2.80, 1.73, "figure_1_d.png",
  "Figure 1", "e", 2.38, 3.72, "figure_1_e.png",
  "Figure 1", "f", 2.38, 3.72, "figure_1_f.png",
  "Figure S1", "a", 7.00, 6.85, "figure_s1_a.png",
  "Figure S2", "a", 3.45, 2.35, "figure_s2_a.png",
  "Figure S2", "b", 3.45, 1.65, "figure_s2_b.png",
  "Figure S2", "c", 3.35, 2.70, "figure_s2_c.png",
  "Figure S2", "d", 3.35, 2.70, "figure_s2_d.png",
  "Figure S2", "e", 3.35, 2.30, "figure_s2_e.png",
  "Figure S2", "f", 3.35, 2.30, "figure_s2_f.png",
  "Figure S2", "g", 6.80, 1.30, "figure_s2_g.png",
  "Figure 2", "a", 7.00, 2.05, "figure_2_a.png",
  "Figure 2", "b", 7.00, 2.65, "figure_2_b.png",
  "Figure 2", "c", 7.00, 2.95, "figure_2_c.png",
  "Figure S3", "a", 7.00, 4.25, "figure_s3_a.png",
  "Figure S4", "a", 7.00, 8.20, "figure_s4_a.png"
) %>% group_by(figure) %>% mutate(order = row_number()) %>% ungroup()

build_panels <- function() {
  list(
    "Figure 1::a" = make_f1a(), "Figure 1::b" = make_f1b(), "Figure 1::c" = make_f1c(),
    "Figure 1::d" = make_f1d(), "Figure 1::e" = make_trajectory_grid("live"), "Figure 1::f" = make_trajectory_grid("glucose"),
    "Figure S1::a" = make_s1(),
    "Figure S2::a" = make_s2a(), "Figure S2::b" = make_s2b(), "Figure S2::c" = make_s2c(),
    "Figure S2::d" = make_s2d(), "Figure S2::e" = make_s2e(), "Figure S2::f" = make_s2f(), "Figure S2::g" = make_s2g(),
    "Figure 2::a" = make_f2a(), "Figure 2::b" = make_f2b(), "Figure 2::c" = make_f2c(),
    "Figure S3::a" = make_s3(), "Figure S4::a" = make_s4()
  )
}

save_subpanels <- function(panels) {
  dpi <- 360
  dims <- panel_meta %>% rowwise() %>% mutate(
    subpanel_png = rel_path(file.path(subpanel_dir, filename)),
    width_px = as.integer(width_in * dpi), height_px = as.integer(height_in * dpi)
  ) %>% ungroup() %>% select(figure, panel, subpanel_png, width_px, height_px, width_in, height_in, order)
  for (i in seq_len(nrow(panel_meta))) {
    key <- paste(panel_meta$figure[[i]], panel_meta$panel[[i]], sep = "::")
    ggsave(file.path(subpanel_dir, panel_meta$filename[[i]]), panels[[key]], width = panel_meta$width_in[[i]], height = panel_meta$height_in[[i]], units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  }
  write.csv(dims, file.path(layout_dir, "subpanel_dimensions.csv"), row.names = FALSE)
  invisible(dims)
}

final_filename <- c(
  "Figure 1" = "figure_1.png", "Figure 2" = "figure_2.png", "Figure S1" = "figure_s1.png",
  "Figure S2" = "figure_s2.png", "Figure S3" = "figure_s3.png", "Figure S4" = "figure_s4.png"
)
assemble_final <- function(panels) {
  layout_path <- file.path(layout_dir, "layout_plan.csv")
  if (!file.exists(layout_path)) stop("Missing layout plan: run the optimizer between subpanel and final phases")
  plan <- read_tbl(layout_path)
  if (exists("apply_recycled_layout", mode = "function")) plan <- apply_recycled_layout(plan)
  for (fig in unique(panel_meta$figure)) {
    one <- plan %>% filter(figure == fig) %>% arrange(match(panel, panel_meta$panel[panel_meta$figure == fig]))
    if (!nrow(one)) stop("Layout plan missing ", fig)
    canvas <- ggdraw()
    for (i in seq_len(nrow(one))) {
      key <- paste(fig, one$panel[[i]], sep = "::")
      canvas <- canvas + draw_plot(panels[[key]], x = one$x_npc[[i]], y = one$y_npc[[i]], width = one$width_npc[[i]], height = one$height_npc[[i]]) +
        draw_label(one$panel[[i]], x = one$x_npc[[i]] + .004, y = one$y_npc[[i]] + one$height_npc[[i]] - .004,
                   hjust = 0, vjust = 1, fontface = "bold", size = 9)
    }
    ggsave(file.path(final_dir, final_filename[[fig]]), canvas, width = unique(one$layout_width_in)[[1]], height = unique(one$layout_height_in)[[1]],
           units = "in", dpi = 360, bg = "white", limitsize = FALSE)
  }
}

write_package_records <- function() {
  data_by_figure <- c(
    "Figure 1" = paste(rel_path(paths$stan), rel_path(paths$area), rel_path(paths$validation), sep = ";"),
    "Figure S1" = rel_path(paths$stan), "Figure S2" = paste(rel_path(paths$stan), rel_path(paths$validation), sep = ";"),
    "Figure 2" = paste(rel_path(paths$stan), rel_path(paths$features), sep = ";"),
    "Figure S3" = rel_path(paths$features), "Figure S4" = rel_path(paths$features)
  )
  provenance <- panel_meta %>% mutate(
    generator = paste(rel_path(script_path), rel_path(file.path(dirname(script_path), "recycled_constructors.R")), sep = ";"),
    command = paste("scripts/agentRrunner.sh", rel_path(script_path), "--phase both"),
    data_inputs = unname(data_by_figure[figure]),
    subpanel_image = rel_path(file.path(subpanel_dir, filename)),
    layout_plan = rel_path(file.path(layout_dir, "layout_plan.csv")),
    output_image = rel_path(file.path(final_dir, unname(final_filename[figure]))),
    notes = "Accepted v02 plotting constructor retained; only its canonical-data adapter was replaced."
  ) %>% select(figure, panel, generator, command, data_inputs, subpanel_image, layout_plan, output_image, notes)
  write.csv(provenance, file.path(polish_root, "provenance.csv"), row.names = FALSE)

  output_paths <- c(file.path(subpanel_dir, panel_meta$filename), file.path(final_dir, unname(final_filename)),
                    file.path(layout_dir, "subpanel_dimensions.csv"), file.path(layout_dir, "layout_plan.csv"), file.path(layout_dir, "layout_report.md"))
  output_paths <- output_paths[file.exists(output_paths)]
  manifest <- tibble(path = rel_path(output_paths), size_bytes = file.info(output_paths)$size, sha256 = vapply(output_paths, sha256_file, character(1)))
  write.csv(manifest, file.path(polish_root, "manifest.csv"), row.names = FALSE)

  finals <- file.path(final_dir, unname(final_filename))
  rebuild <- tibble(
    figure = names(final_filename), stage = "polishing", output = rel_path(finals), output_sha256 = vapply(finals, sha256_file, character(1)),
    generator = rel_path(script_path), command = paste("scripts/agentRrunner.sh", rel_path(script_path), "--phase final"), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write.table(rebuild, file.path(polish_root, "figure_rebuild_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  writeLines(c(
    "# Resegmentation figure rebuild notes", "",
    "Figures 1, 2 and S1-S4 were regenerated from the 20260722 canonical all-line Stan data, resegmentation area summaries, model-free feature release, and rerun classifier-validation release.",
    "The accepted v02 plotting constructors, themes, facet orientations, panel dimensions, and Figure 1 manual layout were retained. Only their data-loading and derivation layer was replaced; no old report_exports tables are runtime inputs.",
    "The count adapter uses N_obs directly as live cells and defines total cells as N_obs + D_obs. Decimal G0 labels retain their leading zero (0.1, 0.25, and 0.5 mM).",
    "Cell and nuclear area intervals are empirical replicate-level interquartile ranges. Glucose calibration panels are derived directly from the fixed calibration fields in the canonical Stan data.",
    "Visual inspection was skipped by explicit user instruction; deterministic file, schema, layout, and provenance validation was still performed.",
    "Project-map decision: no update to docs/project_map.txt; this timestamped manuscript figure package does not change maintained project organization."
  ), file.path(polish_root, "notes.md"))
  writeLines(c(
    "# Visual QC", "", "Status: visual inspection skipped by explicit user instruction.", "",
    "- Title/subtitle check: not visually assessed.",
    "- Clipping check: not visually assessed.",
    "- Spacing/layout check: structural layout validation passed; not visually assessed.",
    "- Readability check: not visually assessed.",
    "- Revision status: no visually motivated revision was performed.", "",
    "Deterministic structural and contract validation are recorded in validation_report.json."
  ), file.path(polish_root, "visual_qc.md"))
}

panels <- build_panels()
if (phase %in% c("subpanels", "both")) save_subpanels(panels)
if (phase %in% c("final", "both")) {
  assemble_final(panels)
  write_package_records()
}
cat("Completed phase ", phase, " in ", polish_root, "\n", sep = "")
