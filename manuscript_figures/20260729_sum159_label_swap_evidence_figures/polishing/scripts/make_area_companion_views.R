#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
script_dir <- dirname(script_path)
package_root <- dirname(script_dir)
find_project_root <- function(start) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "scripts", "agentRrunner.sh"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the project root from ", start, call. = FALSE)
    }
    current <- parent
  }
}
project_root <- find_project_root(package_root)
output_dir <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  file.path(package_root, "derived_data", "area_companion_views")
}
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

processing_root <- file.path(
  project_root, "data", "image_processing_runs",
  "full_segmentation_classification_nuclear",
  "run_20260721_163410"
)
object_path <- file.path(processing_root, "merged", "objects.csv")
metadata_path <- file.path(
  processing_root, "summaries", "cell_nuclear_area_summaries.Rds"
)
if (!file.exists(object_path) || !file.exists(metadata_path)) {
  stop("Required maintained segmentation/classification inputs are missing.", call. = FALSE)
}

batch_key <- tibble::tribble(
  ~cellLine, ~experiment, ~batch_id, ~batch_display,
  "SUM-159-fuse", "C00-IncucyteRawDataLiveDead-varyGlucose",
  "fuse_gluccell_C00", "SUM-159-fuse / GlucCell batch (C00)",
  "SUM-159-fuse", "I00-IncucyteRawDataLiveDead-varyGlucose",
  "fuse_luminescence_I00", "SUM-159-fuse / luminescence batch (I00)",
  "SUM-159-chem", "M00b-IncucyteRawDataLiveDead-varyGlucose",
  "chem_luminescence_M00b", "SUM-159-chem / luminescence batch (M00b)"
)
batch_levels <- batch_key$batch_id
batch_labels <- c(
  fuse_gluccell_C00 = "SUM-159-fuse\nGlucCell (C00)",
  fuse_luminescence_I00 = "SUM-159-fuse\nluminescence (I00)",
  chem_luminescence_M00b = "SUM-159-chem\nluminescence (M00b)"
)
matched_g0 <- c(0.1, 0.5, 1)
selected_hours <- c(24, 32, 40, 48)
all_g0_levels <- c(0, 0.1, 0.25, 0.5, 1, 5, 25)
ploidy_colors <- c(low = "#2166AC", high = "#B2182B")

format_g0 <- function(x) {
  ifelse(abs(x - round(x)) < 1e-10, sprintf("%.0f mM", x),
         paste0(format(x, trim = TRUE, scientific = FALSE), " mM"))
}

metadata_bundle <- readRDS(metadata_path)
if (!is.list(metadata_bundle) || is.null(metadata_bundle$field_qc) ||
    is.null(metadata_bundle$settings)) {
  stop("Metadata RDS lacks the expected field_qc/settings schema.", call. = FALSE)
}
metadata_object_path <- metadata_bundle$settings$provenance$objects_path
if (is.null(metadata_object_path) ||
    normalizePath(metadata_object_path, mustWork = TRUE) !=
      normalizePath(object_path, mustWork = TRUE)) {
  stop("Metadata provenance does not identify the selected object table.", call. = FALSE)
}

metadata_all <- as_tibble(metadata_bundle$field_qc) %>%
  distinct(
    .data$image_key, .data$cellLine, .data$experiment, .data$plateID,
    .data$ploidy, .data$glucose, .data$hours
  ) %>%
  inner_join(batch_key, by = c("cellLine", "experiment")) %>%
  transmute(
    image_key = as.character(.data$image_key),
    cell_line = as.character(.data$cellLine),
    experiment = as.character(.data$experiment),
    batch_id = as.character(.data$batch_id),
    batch_display = as.character(.data$batch_display),
    replicate_well = as.character(.data$plateID),
    ploidy_raw = as.character(.data$ploidy),
    G0_mM = suppressWarnings(as.numeric(as.character(.data$glucose))),
    hours = as.numeric(.data$hours)
  )
metadata_hist <- metadata_all %>%
  filter(
    .data$G0_mM %in% matched_g0,
    .data$hours >= 24,
    .data$hours <= 48
  )

if (!setequal(sort(unique(metadata_hist$hours)), selected_hours)) {
  stop(
    "The selected metadata does not have the expected shared hours: ",
    paste(sort(unique(metadata_hist$hours)), collapse = ", "),
    call. = FALSE
  )
}
coverage_check <- metadata_hist %>%
  distinct(.data$batch_id, .data$ploidy_raw, .data$G0_mM, .data$hours) %>%
  count(.data$batch_id, .data$ploidy_raw, .data$G0_mM, name = "n_hours")
if (nrow(coverage_check) != length(batch_levels) * 2L * length(matched_g0) ||
    any(coverage_check$n_hours != length(selected_hours))) {
  stop("Matched batch/ploidy/G0/hour coverage is incomplete.", call. = FALSE)
}
if (anyDuplicated(metadata_all$image_key)) {
  stop("Target metadata contains duplicate image keys.", call. = FALSE)
}
well_token_matches_filename <- mapply(
  function(image_key, well) {
    grepl(paste0("_", well, "_"), image_key, fixed = TRUE)
  },
  metadata_all$image_key,
  metadata_all$replicate_well,
  USE.NAMES = FALSE
)
if (!all(well_token_matches_filename)) {
  stop("At least one metadata well token does not match its image filename.", call. = FALSE)
}

# Prefix filtering avoids materializing unrelated rows from the 4.1-GB object table.
read_command <- paste(
  "grep -E", shQuote("^(image_key|SUM-159_)"), shQuote(object_path)
)
area_raw <- fread(
  cmd = read_command,
  select = c(
    "image_key", "predicted_label_name", "segmented_area_px", "nuclear_area_px"
  ),
  showProgress = FALSE
)
area_raw[, image_key := as.character(image_key)]
metadata_dt <- as.data.table(metadata_all)
setkey(metadata_dt, image_key)
setkey(area_raw, image_key)
area_all <- metadata_dt[area_raw, nomatch = 0L]
area_all <- area_all[
  predicted_label_name == "alive" &
    is.finite(segmented_area_px) &
    segmented_area_px > 0
]
area_all[, ploidy_state := fcase(
  ploidy_raw == "2N", "low",
  ploidy_raw == "4N", "high",
  default = NA_character_
)]
if (!nrow(area_all) || anyNA(area_all$ploidy_state)) {
  stop("No eligible area data or unexpected ploidy labels.", call. = FALSE)
}
area_hist <- area_all[
  G0_mM %in% matched_g0 &
    hours >= 24 &
    hours <= 48
]

group_columns <- c(
  "cell_line", "experiment", "batch_id", "batch_display",
  "ploidy_state", "G0_mM", "hours"
)
summarize_strata <- function(data) {
  result <- data[
    ,
    .(
      n_alive_cells = .N,
      n_images = uniqueN(image_key),
      n_replicate_wells = uniqueN(replicate_well),
      total_segmented_area_px = sum(segmented_area_px),
      cell_area_q10_px = as.numeric(
        quantile(segmented_area_px, 0.10, type = 7, names = FALSE)
      ),
      cell_area_median_px = as.numeric(
        quantile(segmented_area_px, 0.50, type = 7, names = FALSE)
      ),
      cell_area_q90_px = as.numeric(
        quantile(segmented_area_px, 0.90, type = 7, names = FALSE)
      ),
      n_alive_cells_with_nucleus = sum(
        is.finite(nuclear_area_px) & nuclear_area_px > 0
      ),
      nuclear_area_q10_px = as.numeric(
        quantile(
          nuclear_area_px[is.finite(nuclear_area_px) & nuclear_area_px > 0],
          0.10, type = 7, names = FALSE
        )
      ),
      nuclear_area_median_px = as.numeric(
        quantile(
          nuclear_area_px[is.finite(nuclear_area_px) & nuclear_area_px > 0],
          0.50, type = 7, names = FALSE
        )
      ),
      nuclear_area_q90_px = as.numeric(
        quantile(
          nuclear_area_px[is.finite(nuclear_area_px) & nuclear_area_px > 0],
          0.90, type = 7, names = FALSE
        )
      )
    ),
    by = group_columns
  ]
  result[, batch_order__ := match(batch_id, batch_levels)]
  setorderv(result, c("batch_order__", "G0_mM", "ploidy_state", "hours"))
  result[, batch_order__ := NULL]
  result
}

stratum_summary_matched <- summarize_strata(area_hist)
stratum_summary_full <- summarize_strata(area_all)
fwrite(
  stratum_summary_matched,
  file.path(output_dir, "sum159_matched_24_48h_stratum_summary.csv")
)
fwrite(
  stratum_summary_full,
  file.path(output_dir, "sum159_full_scope_stratum_summary.csv")
)

make_histogram_bins <- function(data, value_column, metric) {
  values <- data[[value_column]]
  limits <- range(values[is.finite(values) & values > 0])
  breaks <- 10^seq(
    floor(log10(limits[[1]])),
    ceiling(log10(limits[[2]])),
    length.out = 81L
  )
  result <- copy(data[is.finite(get(value_column)) & get(value_column) > 0])
  result[, bin_id__ := findInterval(
    get(value_column), breaks, rightmost.closed = TRUE, all.inside = TRUE
  )]
  result <- result[
    ,
    .(n_cells = .N),
    by = .(batch_id, batch_display, G0_mM, ploidy_state, bin_id__)
  ]
  result[, `:=`(
    bin_left_px = breaks[bin_id__],
    bin_right_px = breaks[pmin(bin_id__ + 1L, length(breaks))],
    metric = metric
  )]
  result[
    ,
    fraction_cells := n_cells / sum(n_cells),
    by = .(batch_id, G0_mM, ploidy_state)
  ]
  medians <- data[
    is.finite(get(value_column)) & get(value_column) > 0,
    .(median_px = as.numeric(
      quantile(get(value_column), 0.50, type = 7, names = FALSE)
    )),
    by = .(batch_id, G0_mM, ploidy_state)
  ]
  result <- merge(
    result,
    medians,
    by = c("batch_id", "G0_mM", "ploidy_state"),
    all.x = TRUE,
    sort = FALSE
  )
  result[, bin_id__ := NULL]
  result
}

cell_hist <- make_histogram_bins(
  area_hist, "segmented_area_px", "cell_area"
)
nuclear_hist <- make_histogram_bins(
  area_hist[is.finite(nuclear_area_px) & nuclear_area_px > 0],
  "nuclear_area_px",
  "nuclear_area"
)
histogram_bins <- rbindlist(list(cell_hist, nuclear_hist), use.names = TRUE)
histogram_bins[, batch_order__ := match(batch_id, batch_levels)]
setorderv(
  histogram_bins,
  c("metric", "batch_order__", "G0_mM", "ploidy_state", "bin_left_px")
)
histogram_bins[, batch_order__ := NULL]
fwrite(
  histogram_bins,
  file.path(output_dir, "sum159_matched_24_48h_histogram_bins.csv")
)

hist_plot <- function(hist_data, title, x_label) {
  plot_data <- as_tibble(hist_data) %>%
    mutate(
      batch_display = factor(
        .data$batch_id,
        levels = batch_levels,
        labels = unname(batch_labels[batch_levels])
      ),
      G0_label = factor(
        .data$G0_mM,
        levels = all_g0_levels,
        labels = format_g0(all_g0_levels)
      ),
      ploidy_state = factor(.data$ploidy_state, levels = c("low", "high")),
      ploidy_display = factor(
        .data$ploidy_state,
        levels = c("low", "high"),
        labels = c("Low ploidy (2N)", "High ploidy (4N)")
      )
    )
  median_data <- plot_data %>%
    distinct(
      .data$batch_display, .data$G0_label, .data$ploidy_state,
      .data$ploidy_display, .data$median_px
    )
  ggplot(plot_data) +
    geom_rect(
      aes(
        xmin = .data$bin_left_px,
        xmax = .data$bin_right_px,
        ymin = 0,
        ymax = .data$fraction_cells,
        fill = .data$ploidy_state
      ),
      alpha = 0.42,
      color = NA
    ) +
    geom_vline(
      data = median_data,
      aes(
        xintercept = .data$median_px,
        color = .data$ploidy_state
      ),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    facet_grid(
      rows = vars(batch_display),
      cols = vars(G0_label),
      scales = "free_y"
    ) +
    scale_x_log10(
      labels = label_number(big.mark = ","),
      breaks = trans_breaks("log10", function(x) 10^x)
    ) +
    scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
    scale_fill_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low ploidy (2N)", "High ploidy (4N)"),
      name = NULL
    ) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low ploidy (2N)", "High ploidy (4N)"),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = paste(
        "Hard-classified alive cells pooled over matched G0 = 0.1, 0.5, and 1 mM;",
        "24, 32, 40, and 48 h"
      ),
      x = x_label,
      y = "Fraction of cells per log-spaced bin",
      caption = paste(
        "Low and high ploidy are fully overlaid with transparency. Dashed vertical lines mark medians.",
        "Each histogram is normalized within batch, G0, and ploidy; x-axis and bins are logarithmic."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.22),
      strip.background = element_rect(fill = "grey94", color = "grey72"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "top",
      legend.justification = "left",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 7, hjust = 0),
      axis.text = element_text(size = 8),
      panel.spacing = grid::unit(0.7, "lines"),
      plot.margin = margin(8, 9, 8, 8)
    ) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 0.55)),
      color = "none"
    )
}

prepare_trajectory_data <- function(summary_data, metric_prefix) {
  as_tibble(summary_data) %>%
    transmute(
      batch_id = .data$batch_id,
      batch_display = factor(
        .data$batch_id,
        levels = batch_levels,
        labels = unname(batch_labels[batch_levels])
      ),
      G0_label = factor(
        .data$G0_mM,
        levels = all_g0_levels,
        labels = format_g0(all_g0_levels)
      ),
      ploidy_state = factor(.data$ploidy_state, levels = c("low", "high")),
      hours = .data$hours,
      n_alive_cells = .data$n_alive_cells,
      total_segmented_area_px = .data$total_segmented_area_px,
      q10_px = .data[[paste0(metric_prefix, "_q10_px")]],
      median_px = .data[[paste0(metric_prefix, "_median_px")]],
      q90_px = .data[[paste0(metric_prefix, "_q90_px")]]
    )
}

alternative_x_plot <- function(
    summary_data, metric_prefix, x_column, title, x_label, y_label,
    caption_detail) {
  plot_data <- prepare_trajectory_data(summary_data, metric_prefix)
  ggplot(
    plot_data,
    aes(
      x = .data[[x_column]],
      y = .data$median_px,
      color = .data$ploidy_state,
      group = .data$ploidy_state
    )
  ) +
    geom_path(
      linewidth = 0.65,
      alpha = 0.8,
      arrow = grid::arrow(
        length = grid::unit(0.065, "inches"),
        type = "closed"
      )
    ) +
    geom_errorbar(
      aes(ymin = .data$q10_px, ymax = .data$q90_px),
      width = 0,
      linewidth = 0.38,
      alpha = 0.48
    ) +
    geom_point(size = 1.55) +
    facet_grid(
      rows = vars(G0_label),
      cols = vars(batch_display),
      scales = "free_x"
    ) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low ploidy (2N)", "High ploidy (4N)"),
      name = NULL,
      drop = FALSE
    ) +
    scale_x_continuous(
      labels = label_number(big.mark = ","),
      expand = expansion(mult = c(0.08, 0.10))
    ) +
    scale_y_continuous(
      labels = label_number(big.mark = ","),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      title = title,
      subtitle = paste(
        "All available G0 values and acquisition times;",
        "paths follow time toward the arrow; bars span pooled 10th-90th percentiles"
      ),
      x = x_label,
      y = y_label,
      caption = paste(
        caption_detail,
        "Counts and areas are pooled across the displayed stratum's replicate wells and imaging fields."
      )
    ) +
    theme_bw(base_size = 9) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.22),
      strip.background = element_rect(fill = "grey94", color = "grey72"),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = grid::unit(0.48, "in"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 7, hjust = 0),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 9),
      plot.margin = margin(7, 8, 7, 7)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.1)))
}

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    bg = "white", limitsize = FALSE
  )
}

save_plot(
  hist_plot(
    cell_hist,
    "SUM-159 cell-area distributions by batch and ploidy",
    "Segmented cell area (mask pixels; log scale)"
  ),
  "sum159_cell_area_histograms_matched_24_48h", 10.5, 9
)
save_plot(
  hist_plot(
    nuclear_hist,
    "SUM-159 nuclear-area distributions by batch and ploidy",
    "Nuclear area (mask pixels; log scale)"
  ),
  "sum159_nuclear_area_histograms_matched_24_48h", 10.5, 9
)

cell_count_plot <- alternative_x_plot(
  stratum_summary_full, "cell_area", "n_alive_cells",
  "SUM-159 alive-cell area versus cell count",
  "Hard-classified alive-cell count",
  "Cell area (mask pixels)",
  "The cell count includes the same alive objects used for the cell-area quantiles."
)
nuclear_count_plot <- alternative_x_plot(
  stratum_summary_full, "nuclear_area", "n_alive_cells",
  "SUM-159 nuclear area versus cell count",
  "Hard-classified alive-cell count",
  "Nuclear area (mask pixels)",
  paste(
    "Nuclear quantiles use alive cells with nuclear area > 0;",
    "the x-axis counts all hard-classified alive cells."
  )
)
cell_total_area_plot <- alternative_x_plot(
  stratum_summary_full, "cell_area", "total_segmented_area_px",
  "SUM-159 alive-cell area versus total segmented area",
  "Total segmented alive-cell area (mask pixels)",
  "Cell area (mask pixels)",
  "Total segmented area sums the cell-mask pixels of all hard-classified alive cells."
)
nuclear_total_area_plot <- alternative_x_plot(
  stratum_summary_full, "nuclear_area", "total_segmented_area_px",
  "SUM-159 nuclear area versus total segmented area",
  "Total segmented alive-cell area (mask pixels)",
  "Nuclear area (mask pixels)",
  paste(
    "Nuclear quantiles use alive cells with nuclear area > 0;",
    "total segmented area sums all hard-classified alive cell masks."
  )
)

save_plot(
  cell_count_plot,
  "sum159_cell_area_vs_cell_count_full_scope", 12, 14
)
save_plot(
  nuclear_count_plot,
  "sum159_nuclear_area_vs_cell_count_full_scope", 12, 14
)
save_plot(
  cell_total_area_plot,
  "sum159_cell_area_vs_total_segmented_area_full_scope", 12, 14
)
save_plot(
  nuclear_total_area_plot,
  "sum159_nuclear_area_vs_total_segmented_area_full_scope", 12, 14
)

make_batch_summary <- function(data, strata) {
  result <- data[
  ,
  .(
    n_alive_cells = .N,
    n_alive_cells_with_nucleus = sum(
      is.finite(nuclear_area_px) & nuclear_area_px > 0
    ),
    n_unique_images = uniqueN(image_key),
    n_unique_replicate_wells = uniqueN(replicate_well)
  ),
  by = .(cell_line, experiment, batch_id, batch_display, ploidy_state)
  ]
  stratum_range <- strata[
  ,
  .(
    minimum_cells_per_stratum = min(n_alive_cells),
    maximum_cells_per_stratum = max(n_alive_cells)
  ),
  by = .(cell_line, experiment, batch_id, batch_display, ploidy_state)
  ]
  result <- merge(
    result,
    stratum_range,
    by = c(
      "cell_line", "experiment", "batch_id", "batch_display", "ploidy_state"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  result[, batch_order__ := match(batch_id, batch_levels)]
  setorderv(result, c("batch_order__", "ploidy_state"))
  result[, batch_order__ := NULL]
  result
}
batch_summary_matched <- make_batch_summary(area_hist, stratum_summary_matched)
batch_summary_full <- make_batch_summary(area_all, stratum_summary_full)
fwrite(
  batch_summary_matched,
  file.path(output_dir, "sum159_matched_24_48h_batch_ploidy_summary.csv")
)
fwrite(
  batch_summary_full,
  file.path(output_dir, "sum159_full_scope_batch_ploidy_summary.csv")
)

well_design <- as.data.table(metadata_all)[
  ,
  .(
    cell_line = unique(cell_line),
    experiment = unique(experiment),
    batch_display = unique(batch_display),
    n_ploidy_values = uniqueN(ploidy_raw),
    n_G0_values = uniqueN(G0_mM),
    ploidy_raw = unique(ploidy_raw),
    G0_mM = unique(G0_mM),
    n_images = uniqueN(image_key),
    first_hour = min(hours),
    last_hour = max(hours)
  ),
  by = .(batch_id, replicate_well)
]
if (any(well_design$n_ploidy_values != 1L) ||
    any(well_design$n_G0_values != 1L)) {
  stop("At least one batch/well has inconsistent ploidy or G0 metadata.", call. = FALSE)
}
if (any(!grepl("^[A-Za-z]+[0-9]+$", well_design$replicate_well))) {
  stop("At least one well name cannot be parsed into plate row and column.", call. = FALSE)
}
well_design[, `:=`(
  plate_row = toupper(sub("[0-9]+$", "", replicate_well)),
  plate_col = as.integer(sub("^[A-Za-z]+", "", replicate_well)),
  ploidy_state = fcase(
    ploidy_raw == "2N", "low",
    ploidy_raw == "4N", "high",
    default = NA_character_
  )
)]
if (anyNA(well_design$ploidy_state) ||
    any(!well_design$plate_row %in% LETTERS[1:8]) ||
    any(!well_design$plate_col %in% 1:12)) {
  stop("Parsed wells fall outside an A-H by 1-12 plate or have unknown ploidy.", call. = FALSE)
}
if (well_design[, anyDuplicated(paste(batch_id, plate_row, plate_col))] > 0L) {
  stop("Parsed plate locations are duplicated within a batch.", call. = FALSE)
}

area_platemap <- area_all[hours >= 24 & hours <= 48]
well_area <- area_platemap[
  ,
  .(
    n_alive_cells = .N,
    n_area_images_24_48h = uniqueN(image_key),
    area_first_hour = min(hours),
    area_last_hour = max(hours),
    n_alive_cells_with_nucleus = sum(
      is.finite(nuclear_area_px) & nuclear_area_px > 0
    ),
    median_segmented_area_px = as.numeric(
      quantile(segmented_area_px, 0.50, type = 7, names = FALSE)
    ),
    median_nuclear_area_px = as.numeric(
      quantile(
        nuclear_area_px[is.finite(nuclear_area_px) & nuclear_area_px > 0],
        0.50, type = 7, names = FALSE
      )
    )
  ),
  by = .(batch_id, replicate_well)
]
well_summary <- merge(
  well_design,
  well_area,
  by = c("batch_id", "replicate_well"),
  all.x = TRUE,
  sort = FALSE
)
well_summary[, batch_order__ := match(batch_id, batch_levels)]
setorderv(well_summary, c("batch_order__", "plate_row", "plate_col"))
well_summary[, batch_order__ := NULL]
fwrite(
  well_summary,
  file.path(output_dir, "sum159_full_scope_well_platemap_summary.csv")
)

base_platemap <- function(plot_data, title, subtitle) {
  ggplot(
    plot_data,
    aes(x = .data$plate_col, y = .data$plate_row)
  ) +
    geom_tile(
      aes(fill = .data$fill_value),
      color = "white",
      linewidth = 0.75,
      width = 0.96,
      height = 0.96
    ) +
    scale_x_continuous(
      breaks = 1:12,
      limits = c(0.5, 12.5),
      expand = c(0, 0),
      position = "top"
    ) +
    scale_y_discrete(
      limits = rev(LETTERS[1:8]),
      drop = FALSE,
      expand = c(0, 0)
    ) +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      text = element_text(family = "sans"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 8, color = "grey20"),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5, color = "grey35"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_batch_platemap <- function(batch_id_value) {
  plot_data <- as_tibble(well_summary[batch_id == batch_id_value]) %>%
    mutate(
      plate_row = factor(.data$plate_row, levels = LETTERS[1:8]),
      G0_factor = factor(
        .data$G0_mM,
        levels = all_g0_levels,
        labels = format_g0(all_g0_levels)
      ),
      ploidy_factor = factor(
        .data$ploidy_state,
        levels = c("low", "high"),
        labels = c("Low (2N)", "High (4N)")
      )
    )
  batch_title <- unique(plot_data$batch_display)
  if (length(batch_title) != 1L) {
    stop("Batch display label is not unique for plate map.", call. = FALSE)
  }

  g0_data <- plot_data %>%
    mutate(
      fill_value = .data$G0_factor,
      tile_label = format(.data$G0_mM, trim = TRUE, scientific = FALSE)
    )
  g0_plot <- base_platemap(
    g0_data,
    "Starting glucose (G0)",
    "Labels are mM"
  ) +
    geom_text(aes(label = .data$tile_label), size = 2.55, fontface = "bold") +
    scale_fill_viridis_d(
      option = "D",
      na.value = "grey92",
      name = "G0"
    )

  ploidy_data <- plot_data %>%
    mutate(
      fill_value = .data$ploidy_factor,
      tile_label = ifelse(.data$ploidy_state == "low", "Lo", "Hi")
    )
  ploidy_plot <- base_platemap(
    ploidy_data,
    "Ploidy identity",
    "Low = 2N; High = 4N"
  ) +
    geom_text(
      aes(label = .data$tile_label),
      size = 2.55,
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c("Low (2N)" = ploidy_colors[["low"]],
                 "High (4N)" = ploidy_colors[["high"]]),
      na.value = "grey92",
      name = "Ploidy"
    )

  cell_data <- plot_data %>%
    mutate(
      fill_value = .data$median_segmented_area_px,
      tile_label = comma(round(.data$median_segmented_area_px))
    )
  cell_plot <- base_platemap(
    cell_data,
    "Median segmented-cell area",
    "Hard-classified alive cells acquired from 24-48 h"
  ) +
    geom_text(aes(label = .data$tile_label), size = 2.25, fontface = "bold") +
    scale_fill_viridis_c(
      option = "C",
      na.value = "grey92",
      labels = label_number(big.mark = ","),
      name = "Pixels"
    )

  nuclear_data <- plot_data %>%
    mutate(
      fill_value = .data$median_nuclear_area_px,
      tile_label = comma(round(.data$median_nuclear_area_px))
    )
  nuclear_plot <- base_platemap(
    nuclear_data,
    "Median nuclear area",
    "Alive cells with detected nuclei acquired from 24-48 h"
  ) +
    geom_text(aes(label = .data$tile_label), size = 2.25, fontface = "bold") +
    scale_fill_viridis_c(
      option = "B",
      na.value = "grey92",
      labels = label_number(big.mark = ","),
      name = "Pixels"
    )

  gridExtra::arrangeGrob(
    g0_plot,
    ploidy_plot,
    cell_plot,
    nuclear_plot,
    ncol = 2,
    top = grid::textGrob(
      paste0(batch_title, " plate map"),
      gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 15)
    ),
    bottom = grid::textGrob(
      "Grey wells were not assigned in this batch. Area medians pool hard-classified alive objects acquired from 24 through 48 hours.",
      gp = grid::gpar(fontfamily = "sans", fontsize = 8),
      x = 0.01,
      hjust = 0
    )
  )
}

for (batch_id_value in batch_levels) {
  plate_grob <- make_batch_platemap(batch_id_value)
  stem <- paste0("sum159_", batch_id_value, "_platemap")
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plate_grob,
    width = 14,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plate_grob,
    width = 14,
    height = 9,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

run_metadata <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  object_source = object_path,
  metadata_source = metadata_path,
  common_filters = list(
    batches = unname(batch_key$experiment),
    object_class = "alive",
    ploidy_mapping = list(`2N` = "low", `4N` = "high")
  ),
  histogram_filters = list(
    G0_mM = matched_g0,
    hours_inclusive = c(24, 48),
    realized_hours = selected_hours
  ),
  trajectory_scope = list(
    G0_mM_by_batch = lapply(
      split(metadata_all$G0_mM, metadata_all$batch_id),
      function(x) sort(unique(x))
    ),
    hour_range_by_batch = lapply(
      split(metadata_all$hours, metadata_all$batch_id),
      range
    )
  ),
  counts = list(
    histogram_images = nrow(metadata_hist),
    histogram_alive_objects = nrow(area_hist),
    histogram_alive_objects_with_nucleus = area_hist[
      is.finite(nuclear_area_px) & nuclear_area_px > 0, .N
    ],
    histogram_strata = nrow(stratum_summary_matched),
    full_scope_images = nrow(metadata_all),
    full_scope_alive_objects = nrow(area_all),
    full_scope_alive_objects_with_nucleus = area_all[
      is.finite(nuclear_area_px) & nuclear_area_px > 0, .N
    ],
    full_scope_strata = nrow(stratum_summary_full),
    plate_map_wells = nrow(well_summary),
    plate_map_alive_objects_24_48h = nrow(area_platemap),
    plate_map_alive_objects_with_nucleus_24_48h = area_platemap[
      is.finite(nuclear_area_px) & nuclear_area_px > 0, .N
    ]
  ),
  histogram = list(
    bins = 80,
    spacing = "equal width on log10 area scale",
    normalization = "fraction of cells within batch, G0, and ploidy",
    overlay = "low and high ploidy overlaid with alpha 0.42",
    median_marker = "dashed vertical line"
  ),
  alternative_x_axes = list(
    cell_count = "pooled hard-classified alive-object count per batch/ploidy/G0/hour stratum",
    total_segmented_area = "sum of segmented_area_px over the same alive objects",
    path_order = "chronological by hours"
  ),
  plate_maps = list(
    layout = "one A-H by 1-12 plate per batch; four panels in a 2x2 grid",
    panels = c(
      "G0_mM", "ploidy_state", "median_segmented_area_px",
      "median_nuclear_area_px"
    ),
    area_scope = "acquisitions from 24 through 48 hours inclusive within each well"
  )
)
jsonlite::write_json(
  run_metadata,
  file.path(output_dir, "run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "Wrote matched histograms for ", comma(nrow(area_hist)),
  " alive objects and full-scope trajectories for ", comma(nrow(area_all)),
  " alive objects to: ",
  normalizePath(output_dir)
)
