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
package_root <- dirname(dirname(script_path))
figure_dir <- file.path(package_root, "intermediate_images")
data_dir <- file.path(package_root, "derived_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

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
processing_root <- file.path(
  project_root, "data", "image_processing_runs",
  "full_segmentation_classification_nuclear", "run_20260721_163410"
)
object_path <- file.path(processing_root, "merged", "objects.csv")
metadata_path <- file.path(
  processing_root, "summaries", "cell_nuclear_area_summaries.Rds"
)
mixture_calls_path <- file.path(
  data_dir, "competition_ploidy_cytoplasmic_gfp_full",
  "sum159_competition_object_ploidy_calls.csv.gz"
)
cytoplasm_path <- file.path(
  data_dir, "cytoplasmic_red_quantification",
  "per_cell_cytoplasmic_red.csv.gz"
)
platemap_path <- file.path(
  data_dir, "area_companion_views",
  "sum159_full_scope_well_platemap_summary.csv"
)
required <- c(
  object_path, metadata_path, mixture_calls_path, cytoplasm_path, platemap_path
)
if (any(!file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

low_color <- "#2166AC"
high_color <- "#B2182B"
ploidy_colors <- c(low = low_color, high = high_color)
context_levels <- c("C00 fuse", "I00 fuse", "M00b chem", "J00 mixture")
context_labels <- c(
  "C00 fuse" = "SUM-159-fuse\nC00",
  "I00 fuse" = "SUM-159-fuse\nI00",
  "M00b chem" = "SUM-159-chem\nM00b",
  "J00 mixture" = "Mixture\nGFP-defined"
)
focus_g0 <- c(0.1, 0.5, 1)
dynamic_g0 <- c(0.5, 1)

format_g0 <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-10,
    sprintf("%.0f mM", x),
    paste0(format(x, trim = TRUE, scientific = FALSE), " mM")
  )
}

figure_theme <- theme_bw(base_size = 8.5) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.22),
    strip.background = element_rect(fill = "grey94", color = "grey70", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 7.4),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = grid::unit(0.42, "in"),
    legend.margin = margin(0, 0, 1, 0),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.8),
    panel.spacing = grid::unit(0.52, "lines"),
    plot.margin = margin(4, 5, 4, 5),
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.002, 0.995)
  )

metadata_bundle <- readRDS(metadata_path)
metadata <- as_tibble(metadata_bundle$field_qc) %>%
  distinct(
    .data$image_key, .data$cellLine, .data$experiment, .data$plateID,
    .data$ploidy, .data$glucose, .data$hours
  ) %>%
  filter(
    (.data$cellLine == "SUM-159-fuse" &
       .data$experiment %in% c(
         "C00-IncucyteRawDataLiveDead-varyGlucose",
         "I00-IncucyteRawDataLiveDead-varyGlucose"
       )) |
      (.data$cellLine == "SUM-159-chem" &
         .data$experiment == "M00b-IncucyteRawDataLiveDead-varyGlucose")
  ) %>%
  transmute(
    image_key = as.character(.data$image_key),
    well = as.character(.data$plateID),
    context = case_when(
      grepl("^C00-", .data$experiment) ~ "C00 fuse",
      grepl("^I00-", .data$experiment) ~ "I00 fuse",
      grepl("^M00b-", .data$experiment) ~ "M00b chem",
      TRUE ~ NA_character_
    ),
    ploidy_state = case_when(
      as.character(.data$ploidy) == "2N" ~ "low",
      as.character(.data$ploidy) == "4N" ~ "high",
      TRUE ~ NA_character_
    ),
    G0_mM = as.numeric(as.character(.data$glucose)),
    hours = as.numeric(.data$hours)
  )
if (anyNA(metadata$context) || anyNA(metadata$ploidy_state) ||
    anyDuplicated(metadata$image_key)) {
  stop("Monoculture metadata mapping failed.")
}

read_command <- paste(
  "grep -E", shQuote("^(image_key|SUM-159_)"), shQuote(object_path)
)
mono_objects <- fread(
  cmd = read_command,
  select = c(
    "image_key", "predicted_label_name", "segmented_area_px",
    "nuclear_area_px"
  ),
  showProgress = FALSE
)
mono_objects[, image_key := as.character(image_key)]
metadata_dt <- as.data.table(metadata)
setkey(metadata_dt, image_key)
setkey(mono_objects, image_key)
mono <- metadata_dt[mono_objects, nomatch = 0L][
  predicted_label_name == "alive" &
    is.finite(segmented_area_px) &
    segmented_area_px > 0
]
mono <- mono[, .(
  context, ploidy_state, G0_mM, hours, well, image_key,
  cell_area_px = as.numeric(segmented_area_px),
  nuclear_area_px = as.numeric(nuclear_area_px)
)]

mixture <- fread(mixture_calls_path, showProgress = FALSE)[
  predicted_label_name == "alive" &
    cell_area_pass == 1 &
    ploidy_call %chin% c("low", "high") &
    is.finite(cell_area_px) &
    cell_area_px > 0
]
mixture <- mixture[, .(
  context = "J00 mixture",
  ploidy_state = as.character(ploidy_call),
  G0_mM = as.numeric(G0_mM),
  hours = as.numeric(hours),
  well = as.character(well),
  image_key = as.character(image_key),
  cell_area_px = as.numeric(cell_area_px),
  nuclear_area_px = as.numeric(nuclear_area_px)
)]

area <- rbindlist(list(mono, mixture), use.names = TRUE)
area[, context := factor(context, levels = context_levels)]
area[, ploidy_state := factor(ploidy_state, levels = c("low", "high"))]
if (anyNA(area$context) || anyNA(area$ploidy_state)) {
  stop("Combined area table contains unmapped context or ploidy.")
}

make_log_hist_bins <- function(data, value_column, metric) {
  values <- data[[value_column]]
  values <- values[is.finite(values) & values > 0]
  breaks <- 10^seq(
    floor(log10(min(values))),
    ceiling(log10(max(values))),
    length.out = 71L
  )
  work <- copy(data[is.finite(get(value_column)) & get(value_column) > 0])
  work[, bin_id__ := findInterval(
    get(value_column), breaks, rightmost.closed = TRUE, all.inside = TRUE
  )]
  bins <- work[, .(n_cells = .N), by = .(
    context, G0_mM, ploidy_state, bin_id__
  )]
  bins[, `:=`(
    bin_left_px = breaks[bin_id__],
    bin_right_px = breaks[pmin(bin_id__ + 1L, length(breaks))],
    metric = metric
  )]
  bins[, fraction_cells := n_cells / sum(n_cells),
       by = .(context, G0_mM, ploidy_state)]
  medians <- work[, .(
    median_px = as.numeric(quantile(
      get(value_column), 0.5, type = 7, names = FALSE
    ))
  ), by = .(context, G0_mM, ploidy_state)]
  bins <- merge(
    bins, medians,
    by = c("context", "G0_mM", "ploidy_state"),
    all.x = TRUE, sort = FALSE
  )
  bins[, bin_id__ := NULL]
  bins
}

hist_source <- area[
  G0_mM %in% focus_g0 &
    hours >= 24 &
    hours <= 48
]
cell_hist <- make_log_hist_bins(hist_source, "cell_area_px", "cell")
nuclear_hist <- make_log_hist_bins(
  hist_source[is.finite(nuclear_area_px) & nuclear_area_px > 0],
  "nuclear_area_px", "nuclear"
)
hist_bins <- rbindlist(list(cell_hist, nuclear_hist), use.names = TRUE)
fwrite(hist_bins, file.path(data_dir, "area_histogram_bins.csv"))

make_hist_plot <- function(data, x_label, tag) {
  plot_data <- as_tibble(data) %>%
    mutate(
      context_label = factor(
        as.character(.data$context),
        levels = context_levels,
        labels = unname(context_labels[context_levels])
      ),
      G0_label = factor(
        .data$G0_mM,
        levels = focus_g0,
        labels = format_g0(focus_g0)
      )
    )
  median_data <- plot_data %>%
    distinct(
      .data$context_label, .data$G0_label, .data$ploidy_state,
      .data$median_px
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
      aes(xintercept = .data$median_px, color = .data$ploidy_state),
      linewidth = 0.58,
      linetype = "dashed"
    ) +
    facet_grid(
      rows = vars(context_label),
      cols = vars(G0_label),
      scales = "free_y"
    ) +
    scale_x_log10(labels = function(x) {
      ifelse(
        x >= 1e6,
        paste0(format(x / 1e6, trim = TRUE), "M"),
        ifelse(
          x >= 1e3,
          paste0(format(x / 1e3, trim = TRUE), "k"),
          comma(x)
        )
      )
    }) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    scale_fill_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = NULL
    ) +
    scale_color_manual(values = ploidy_colors, guide = "none") +
    labs(
      tag = tag,
      x = x_label,
      y = "Cells per log-area bin"
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.55))) +
    figure_theme
}

cell_hist_plot <- make_hist_plot(
  cell_hist, "Segmented-cell area (pixels; log scale)", "a"
)
nuclear_hist_plot <- make_hist_plot(
  nuclear_hist, "Nuclear area (pixels; log scale)", "b"
) + theme(legend.position = "none")
hist_figure <- arrangeGrob(
  cell_hist_plot, nuclear_hist_plot,
  ncol = 1, heights = c(1, 1)
)
ggsave(
  file.path(figure_dir, "area_distributions_24_48h.png"),
  hist_figure, width = 7, height = 9.2, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

summarize_area <- function(data, value_column, prefix) {
  result <- data[
    is.finite(get(value_column)) & get(value_column) > 0,
    .(
      n_objects = .N,
      n_images = uniqueN(image_key),
      n_wells = uniqueN(well),
      q10 = as.numeric(quantile(
        get(value_column), 0.10, type = 7, names = FALSE
      )),
      median = as.numeric(quantile(
        get(value_column), 0.50, type = 7, names = FALSE
      )),
      q90 = as.numeric(quantile(
        get(value_column), 0.90, type = 7, names = FALSE
      ))
    ),
    by = .(context, G0_mM, hours, ploidy_state)
  ]
  result[, metric := prefix]
  result
}

dynamic_source <- area[G0_mM %in% dynamic_g0]
cell_summary <- summarize_area(dynamic_source, "cell_area_px", "cell")
nuclear_summary <- summarize_area(
  dynamic_source, "nuclear_area_px", "nuclear"
)
area_summary <- rbindlist(list(cell_summary, nuclear_summary), use.names = TRUE)
fwrite(area_summary, file.path(data_dir, "area_time_quantiles.csv"))

prepare_summary_plot <- function(data) {
  as_tibble(data) %>%
    mutate(
      context_label = factor(
        as.character(.data$context),
        levels = context_levels,
        labels = unname(context_labels[context_levels])
      ),
      G0_label = factor(
        .data$G0_mM,
        levels = dynamic_g0,
        labels = format_g0(dynamic_g0)
      )
    )
}

make_time_plot <- function(data, y_label, tag) {
  plot_data <- prepare_summary_plot(data)
  ggplot(
    plot_data,
    aes(
      x = .data$hours,
      y = .data$median,
      color = .data$ploidy_state,
      fill = .data$ploidy_state,
      group = .data$ploidy_state
    )
  ) +
    geom_ribbon(
      aes(ymin = .data$q10, ymax = .data$q90),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(linewidth = 0.62) +
    facet_grid(rows = vars(context_label), cols = vars(G0_label)) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = NULL
    ) +
    scale_fill_manual(values = ploidy_colors, guide = "none") +
    scale_x_continuous(breaks = breaks_width(24)) +
    scale_y_continuous(
      labels = label_number(big.mark = ","),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(tag = tag, x = "Time (h)", y = y_label) +
    figure_theme
}

time_cell_plot <- make_time_plot(
  cell_summary, "Segmented-cell area (pixels)", "a"
)
time_nuclear_plot <- make_time_plot(
  nuclear_summary, "Nuclear area (pixels)", "b"
) + theme(legend.position = "none")
time_figure <- arrangeGrob(
  time_cell_plot, time_nuclear_plot,
  ncol = 1, heights = c(1, 1)
)
ggsave(
  file.path(figure_dir, "area_over_time_with_mixture.png"),
  time_figure, width = 7, height = 9.2, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

density_x <- dynamic_source[, .(
  total_segmented_area_px = sum(cell_area_px),
  n_cells = .N
), by = .(context, G0_mM, hours, ploidy_state)]
density_summary <- merge(
  area_summary,
  density_x,
  by = c("context", "G0_mM", "hours", "ploidy_state"),
  all.x = TRUE,
  sort = FALSE
)
fwrite(density_summary, file.path(data_dir, "area_vs_total_segmented_area.csv"))

make_density_plot <- function(data, y_label, tag) {
  plot_data <- prepare_summary_plot(data)
  ggplot(
    plot_data,
    aes(
      x = .data$total_segmented_area_px,
      y = .data$median,
      color = .data$ploidy_state,
      group = .data$ploidy_state
    )
  ) +
    geom_errorbar(
      aes(ymin = .data$q10, ymax = .data$q90),
      width = 0,
      linewidth = 0.26,
      alpha = 0.42
    ) +
    geom_path(
      linewidth = 0.56,
      arrow = grid::arrow(
        length = grid::unit(0.045, "inches"),
        type = "closed"
      )
    ) +
    geom_point(size = 1.15) +
    facet_grid(rows = vars(context_label), cols = vars(G0_label)) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = NULL
    ) +
    scale_x_log10(labels = label_number(big.mark = ",")) +
    scale_y_continuous(
      labels = label_number(big.mark = ","),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      tag = tag,
      x = "Total segmented area (pixels; log scale)",
      y = y_label
    ) +
    figure_theme
}

density_cell_plot <- make_density_plot(
  density_summary[metric == "cell"],
  "Segmented-cell area (pixels)", "a"
)
density_nuclear_plot <- make_density_plot(
  density_summary[metric == "nuclear"],
  "Nuclear area (pixels)", "b"
) + theme(legend.position = "none")
density_figure <- arrangeGrob(
  density_cell_plot, density_nuclear_plot,
  ncol = 1, heights = c(1, 1)
)
ggsave(
  file.path(figure_dir, "area_vs_total_segmented_area_with_mixture.png"),
  density_figure, width = 7, height = 9.2, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

platemap <- fread(platemap_path)
platemap[, context := fcase(
  batch_id == "fuse_gluccell_C00", "C00 fuse",
  batch_id == "fuse_luminescence_I00", "I00 fuse",
  batch_id == "chem_luminescence_M00b", "M00b chem",
  default = NA_character_
)]
platemap[, context := factor(context, levels = context_levels[1:3])]
platemap[, ploidy_state := factor(ploidy_state, levels = c("low", "high"))]
if (anyNA(platemap$context) || anyNA(platemap$ploidy_state) ||
    any(platemap$area_first_hour != 24) || any(platemap$area_last_hour != 48)) {
  stop("Plate-map summary does not satisfy the expected 24-48 h contract.")
}

make_platemap_panel <- function(value_column, legend_title, tag) {
  plot_data <- as_tibble(platemap) %>%
    mutate(
      context_label = factor(
        as.character(.data$context),
        levels = context_levels[1:3],
        labels = unname(context_labels[context_levels[1:3]])
      ),
      plate_row = factor(.data$plate_row, levels = LETTERS[1:8]),
      value = .data[[value_column]]
    )
  ggplot(
    plot_data,
    aes(x = .data$plate_col, y = .data$plate_row)
  ) +
    geom_tile(
      aes(fill = .data$value, color = .data$ploidy_state),
      width = 0.94,
      height = 0.94,
      linewidth = 0.72
    ) +
    geom_text(
      aes(label = comma(round(.data$value))),
      size = 1.8,
      fontface = "bold"
    ) +
    facet_grid(rows = vars(context_label)) +
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
    scale_fill_viridis_c(
      option = "C",
      name = legend_title,
      labels = label_number(big.mark = ",")
    ) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = "Tile border"
    ) +
    coord_equal(clip = "off") +
    labs(tag = tag, x = NULL, y = NULL) +
    figure_theme +
    theme(
      axis.text = element_text(size = 6),
      panel.spacing = grid::unit(0.35, "lines"),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}

platemap_cell <- make_platemap_panel(
  "median_segmented_area_px", "Median cell area", "a"
)
platemap_nuclear <- make_platemap_panel(
  "median_nuclear_area_px", "Median nuclear area", "b"
)
platemap_figure <- arrangeGrob(
  platemap_cell, platemap_nuclear,
  ncol = 2, widths = c(1, 1)
)
ggsave(
  file.path(figure_dir, "platemap_area_qc_24_48h.png"),
  platemap_figure, width = 7, height = 6.6, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

cyto <- fread(cytoplasm_path)
cyto[, context := fcase(
  assay == "C00 monoculture", "C00 fuse",
  assay == "I00 monoculture", "I00 fuse",
  assay == "M00b monoculture", "M00b chem",
  assay == "J00 mixture", "J00 mixture",
  default = NA_character_
)]
cyto[, context := factor(context, levels = context_levels)]
cyto[, ploidy_state := factor(ploidy_state, levels = c("low", "high"))]
if (anyNA(cyto$context) || anyNA(cyto$ploidy_state)) {
  stop("Cytoplasmic-red context mapping failed.")
}
cyto_image <- cyto[, .(
  n_cells = .N,
  image_median_bg_z = median(red_cytoplasm_bg_z)
), by = .(
  context, assay, image_key, ploidy_state, G0_mM, hours, position
)]
cyto_image[, pair_id := fifelse(
  context == "J00 mixture",
  image_key,
  paste(context, G0_mM, hours, position, sep = "|")
)]
cyto_contrast <- dcast(
  cyto_image,
  context + G0_mM + hours + pair_id ~ ploidy_state,
  value.var = "image_median_bg_z"
)
cyto_contrast[, high_minus_low_bg_z := high - low]
fwrite(cyto_image, file.path(data_dir, "cytoplasmic_red_image_medians.csv"))
fwrite(cyto_contrast, file.path(data_dir, "cytoplasmic_red_matched_contrasts.csv"))

cyto_plot_data <- as_tibble(cyto_image) %>%
  mutate(
    context_label = factor(
      as.character(.data$context),
      levels = context_levels,
      labels = unname(context_labels[context_levels])
    )
  )
cyto_panel_a <- ggplot(
  cyto_plot_data,
  aes(
    x = .data$ploidy_state,
    y = .data$image_median_bg_z,
    color = .data$ploidy_state
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    color = "grey35",
    fill = "white",
    linewidth = 0.38
  ) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 157),
    size = 1.55, alpha = 0.82
  ) +
  facet_wrap(vars(context_label), nrow = 1) +
  scale_color_manual(values = ploidy_colors, guide = "none") +
  scale_x_discrete(labels = c(low = "Low", high = "High")) +
  labs(
    tag = "a",
    x = NULL,
    y = "Image-median cytoplasmic red\n(robust background units)"
  ) +
  figure_theme +
  theme(strip.text = element_text(size = 6.8))

contrast_plot_data <- as_tibble(cyto_contrast) %>%
  mutate(
    context_label = factor(
      as.character(.data$context),
      levels = context_levels,
      labels = unname(context_labels[context_levels])
    ),
    G0_label = factor(
      .data$G0_mM,
      levels = sort(unique(.data$G0_mM)),
      labels = format_g0(sort(unique(.data$G0_mM)))
    )
  )
cyto_panel_b <- ggplot(
  contrast_plot_data,
  aes(
    x = .data$context_label,
    y = .data$high_minus_low_bg_z,
    color = .data$G0_label
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(
    aes(group = .data$context_label),
    width = 0.48,
    outlier.shape = NA,
    color = "grey35",
    fill = "grey96",
    linewidth = 0.38
  ) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 158),
    size = 1.55, alpha = 0.85
  ) +
  scale_color_viridis_d(name = "G0") +
  labs(
    tag = "b",
    x = NULL,
    y = "High - low matched image median"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 24, hjust = 1),
    legend.key.width = grid::unit(0.22, "in")
  )

cyto_figure <- arrangeGrob(
  cyto_panel_a, cyto_panel_b,
  ncol = 2, widths = c(1.35, 1)
)
ggsave(
  file.path(figure_dir, "cytoplasmic_red_quantification.png"),
  cyto_figure, width = 7, height = 3.5, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

context_summary <- area[, .(
  n_objects = .N,
  n_images = uniqueN(image_key),
  n_wells = uniqueN(well),
  min_hour = min(hours),
  max_hour = max(hours),
  G0_values = paste(sort(unique(G0_mM)), collapse = ",")
), by = context]
fwrite(context_summary, file.path(data_dir, "area_context_summary.csv"))

cat(
  "Wrote five quantitative figure composites from ",
  format(nrow(area), big.mark = ","), " area observations and ",
  format(nrow(cyto), big.mark = ","), " cytoplasmic-red observations.\n",
  sep = ""
)
