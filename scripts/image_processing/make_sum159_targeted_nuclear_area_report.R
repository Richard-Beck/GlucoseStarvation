#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: make_sum159_targeted_nuclear_area_report.R MERGED_DIR MANIFEST_CSV OUTPUT_DIR FIGURE_DIR [ALIVE_FEATURE_CSV]", call. = FALSE)
}
merged_dir <- args[[1]]
manifest_csv <- args[[2]]
output_dir <- args[[3]]
figure_dir <- args[[4]]
alive_feature_csv <- if (length(args) >= 5L) args[[5]] else file.path(
  "data", "image_processing_runs", "run_20260324_233122", "object_features.csv"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

objects <- fread(file.path(merged_dir, "wp3_nuclear_object_features.csv"), showProgress = FALSE)
images <- fread(file.path(merged_dir, "wp3_nuclear_image_qc.csv"), showProgress = FALSE)
manifest <- fread(manifest_csv, showProgress = FALSE)
stopifnot(
  uniqueN(manifest$image_key) == nrow(manifest),
  uniqueN(images$image_key) == nrow(manifest),
  setequal(images$image_key, manifest$image_key),
  all(images$status == "ok")
)

manifest_columns <- c("image_key", "manifest_index", "shard_id", "cellLine", "experiment", "batch_id", "plateID", "pos", "ploidy", "G0", "glucose_bin", "time_bin", "hours")
manifest_small <- manifest[, ..manifest_columns]
objects[, object_id := as.integer(object_id)]
objects[, cell_area_pass := as.integer(cell_area_pass)]
objects[, `:=`(
  cell_area_px = as.numeric(cell_area_px),
  nuclear_area_px = as.numeric(nuclear_area_px),
  largest_nuclear_area_px = as.numeric(largest_nuclear_area_px),
  nuclear_component_count = as.integer(nuclear_component_count),
  nuclear_to_cell_area_ratio = as.numeric(nuclear_to_cell_area_ratio)
)]
objects[, c("cellLine", "ploidy", "G0", "glucose_bin", "time_bin", "hours") := NULL]
objects <- merge(objects, manifest_small, by = "image_key", all.x = TRUE, sort = FALSE)

message("Reading production alive/dead scores from ", alive_feature_csv)
alive_scores <- fread(
  alive_feature_csv,
  select = c("image_key", "object_id", "live_bg_z", "dead_bg_z"),
  showProgress = FALSE
)[image_key %chin% manifest$image_key]
alive_scores[, object_id := as.integer(object_id)]
setkey(alive_scores, image_key, object_id)
setkey(objects, image_key, object_id)
objects <- alive_scores[objects]
match_fraction <- mean(is.finite(objects$live_bg_z) & is.finite(objects$dead_bg_z))
if (!is.finite(match_fraction) || match_fraction < 0.98) {
  stop(sprintf("Only %.2f%% of re-segmented objects matched production alive scores.", 100 * match_fraction), call. = FALSE)
}
objects[, `:=`(
  alive_score = live_bg_z - dead_bg_z,
  is_alive = is.finite(live_bg_z) & is.finite(dead_bg_z) & live_bg_z - dead_bg_z >= 0,
  has_nucleus = is.finite(nuclear_area_px) & nuclear_area_px > 0,
  ploidy_state = factor(ifelse(ploidy == "2N", "baseline", "elevated"), levels = c("baseline", "elevated"))
)]

measured <- objects[
  cell_area_pass == 1L & is_alive & has_nucleus &
    is.finite(nuclear_area_px) & nuclear_area_px > 0
]
if (!nrow(measured)) stop("No alive, area-passing cells with nuclei remain.", call. = FALSE)

trajectory <- measured[
  ,
  .(
    n_images = uniqueN(image_key),
    n_alive_cells_with_nucleus = .N,
    nuclear_area_q10_px = as.numeric(quantile(nuclear_area_px, 0.10, type = 7)),
    nuclear_area_median_px = as.numeric(median(nuclear_area_px)),
    nuclear_area_q90_px = as.numeric(quantile(nuclear_area_px, 0.90, type = 7)),
    largest_nuclear_area_median_px = as.numeric(median(largest_nuclear_area_px)),
    nuclear_component_count_median = as.numeric(median(nuclear_component_count)),
    nuclear_to_cell_area_ratio_median = as.numeric(median(nuclear_to_cell_area_ratio))
  ),
  by = .(cellLine, experiment, batch_id, ploidy, ploidy_state, G0, glucose_bin, time_bin, hours)
][order(batch_id, G0, hours, ploidy)]
trajectory[, estimable := n_images == 2L & n_alive_cells_with_nucleus >= 10L]
trajectory[, `:=`(
  plot_nuclear_area_q10_px = fifelse(estimable, nuclear_area_q10_px, NA_real_),
  plot_nuclear_area_median_px = fifelse(estimable, nuclear_area_median_px, NA_real_),
  plot_nuclear_area_q90_px = fifelse(estimable, nuclear_area_q90_px, NA_real_)
)]

image_summary <- measured[
  ,
  .(
    n_alive_cells_with_nucleus = .N,
    nuclear_area_median_px = as.numeric(median(nuclear_area_px)),
    largest_nuclear_area_median_px = as.numeric(median(largest_nuclear_area_px)),
    nuclear_component_count_median = as.numeric(median(nuclear_component_count)),
    nuclear_to_cell_area_ratio_median = as.numeric(median(nuclear_to_cell_area_ratio))
  ),
  by = .(image_key, cellLine, experiment, batch_id, plateID, pos, ploidy, ploidy_state, G0, glucose_bin, time_bin, hours)
][order(batch_id, G0, hours, ploidy, plateID, pos)]

field_reproducibility <- image_summary[
  ,
  .(
    n_fields = .N,
    min_field_cells = min(n_alive_cells_with_nucleus),
    max_field_cells = max(n_alive_cells_with_nucleus),
    nuclear_area_field_range_px = if (.N == 2L) diff(range(nuclear_area_median_px)) else NA_real_,
    nuclear_area_field_log2_ratio = if (.N == 2L) log2(max(nuclear_area_median_px) / min(nuclear_area_median_px)) else NA_real_
  ),
  by = .(cellLine, experiment, batch_id, ploidy, ploidy_state, G0, hours)
][order(batch_id, G0, hours, ploidy)]

fwrite(measured, file.path(output_dir, "sum159_nuclear_objects_alive.csv.gz"))
fwrite(trajectory, file.path(output_dir, "sum159_nuclear_area_time_quantiles.csv"))
fwrite(image_summary, file.path(output_dir, "sum159_nuclear_area_image_summary.csv"))
fwrite(field_reproducibility, file.path(output_dir, "sum159_nuclear_area_field_reproducibility.csv"))

batch_levels <- c("C00", "I00", "M00b")
batch_labels <- c(
  C00 = "SUM-159-fuse\nGlucCell batch (C00)",
  I00 = "SUM-159-fuse\nluminescence batch (I00)",
  M00b = "SUM-159-chem\nluminescence batch (M00b)"
)
g0_levels <- c(0, 0.1, 0.25, 0.5, 1, 5, 25)
format_g0 <- function(x) {
  ifelse(abs(x - round(x)) < 1e-10, sprintf("%.0f mM", x), paste0(format(x, trim = TRUE), " mM"))
}
trajectory[, `:=`(
  batch_display = factor(batch_id, levels = batch_levels, labels = unname(batch_labels[batch_levels])),
  G0_label = factor(G0, levels = g0_levels, labels = format_g0(g0_levels))
)]
insufficient <- trajectory[estimable == FALSE]
insufficient[, hours_marker := hours + fifelse(ploidy == "2N", -1.5, 1.5)]

missing_facets <- CJ(batch_id = batch_levels, G0 = g0_levels)[
  !unique(trajectory[, .(batch_id, G0)]), on = .(batch_id, G0)
]
missing_facets[, `:=`(
  batch_display = factor(batch_id, levels = batch_levels, labels = unname(batch_labels[batch_levels])),
  G0_label = factor(G0, levels = g0_levels, labels = format_g0(g0_levels)),
  hours = median(trajectory$hours),
  nuclear_area_median_px = 0.52 * max(trajectory$nuclear_area_q90_px)
)]

colors <- c(baseline = "#2166AC", elevated = "#B2182B")
p <- ggplot(
  trajectory,
  aes(
    x = hours,
    y = plot_nuclear_area_median_px,
    color = ploidy_state,
    fill = ploidy_state,
    group = ploidy_state
  )
) +
  geom_ribbon(aes(ymin = plot_nuclear_area_q10_px, ymax = plot_nuclear_area_q90_px), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.1, na.rm = TRUE) +
  geom_point(
    data = insufficient,
    aes(x = hours_marker, y = 0, color = ploidy_state),
    inherit.aes = FALSE,
    shape = 4,
    size = 1.6,
    stroke = 0.6
  ) +
  geom_text(
    data = missing_facets,
    aes(x = hours, y = nuclear_area_median_px, label = "not measured"),
    inherit.aes = FALSE,
    color = "grey55",
    size = 2.5
  ) +
  facet_grid(rows = vars(G0_label), cols = vars(batch_display), drop = FALSE) +
  scale_color_manual(values = colors, labels = c("Low ploidy (2N)", "High ploidy (4N)"), name = NULL) +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(breaks = breaks_width(24), expand = expansion(mult = c(0.01, 0.025))) +
  scale_y_continuous(labels = label_number(big.mark = ","), limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "SUM-159 nuclear area trajectories by experimental batch",
    subtitle = "Two longitudinal well-position tracks per condition; medians and pooled 10th-90th percentiles; x marks insufficient endpoints",
    x = "Time (hours)",
    y = "Nuclear area (mask pixels)",
    caption = "Alive, area-passing cells with a detected nucleus; alive gate uses production live-minus-dead z-scores. Estimates require both fields and at least 10 cells."
  ) +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.22),
    strip.background = element_rect(fill = "grey94", color = "grey72", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "top",
    legend.justification = "left",
    panel.spacing = grid::unit(0.8, "lines"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 7, hjust = 0),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.margin = margin(7, 8, 7, 7)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.1)))

ggsave(file.path(figure_dir, "sum159_nuclear_area_trajectories.png"), p, width = 12, height = 14, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "sum159_nuclear_area_trajectories.pdf"), p, width = 12, height = 14, units = "in", bg = "white")

summary_lines <- c(
  "# Targeted SUM-159 nuclear-area report",
  "",
  sprintf("- Images: %s", comma(nrow(images))),
  sprintf("- Re-segmented objects: %s", comma(nrow(objects))),
  sprintf("- Production alive-score join: %.2f%%", 100 * match_fraction),
  sprintf("- Alive, area-passing cells with nuclei: %s", comma(nrow(measured))),
  sprintf("- Trajectory rows: %s", comma(nrow(trajectory))),
  sprintf("- Estimable trajectory rows (both fields and at least 10 cells): %s/%s", comma(sum(trajectory$estimable)), comma(nrow(trajectory))),
  sprintf("- Insufficient late-starvation endpoints: %s", comma(sum(!trajectory$estimable))),
  sprintf("- Median absolute field log2 ratio: %.3f", median(abs(field_reproducibility$nuclear_area_field_log2_ratio), na.rm = TRUE))
)
writeLines(summary_lines, file.path(output_dir, "run_summary.md"))
cat(paste(summary_lines, collapse = "\n"), "\n")
