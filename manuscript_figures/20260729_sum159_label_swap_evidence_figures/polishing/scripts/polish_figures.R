#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(scales)
  library(tiff)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
package_root <- dirname(dirname(script_path))
phase_index <- match("--phase", args)
phase <- if (!is.na(phase_index) && length(args) >= phase_index + 1L) {
  args[[phase_index + 1L]]
} else {
  "final"
}
if (!phase %in% c("subpanels", "final")) {
  stop("Unknown --phase value: ", phase)
}
source(file.path(package_root, "scripts", "make_evidence_figures.R"), local = FALSE)

v3_dir <- file.path(package_root, "final_images")
dir.create(v3_dir, recursive = TRUE, showWarnings = FALSE)
subpanel_dir <- file.path(package_root, "subpanels")
layout_dir <- file.path(package_root, "layout")
dir.create(subpanel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
shared_g0 <- c(0.1, 0.5, 1)
all_g0 <- sort(unique(area$G0_mM))
all_g0_labels <- vapply(all_g0, function(x) {
  base <- if (abs(x - round(x)) < 1e-10) sprintf("%.0f mM", x) else
    paste0(format(x, trim = TRUE, scientific = FALSE), " mM")
  if (x %in% shared_g0) paste0("★ ", base) else base
}, character(1))

context_display <- c(
  "C00 fuse" = "SUM-159-fuse C00",
  "I00 fuse" = "SUM-159-fuse I00",
  "M00b chem" = "SUM-159-chem M00b",
  "J00 mixture" = "Mixture GFP-defined"
)

# ---- Figure 1: complete time-course overview --------------------------------
time_all <- area[
  is.finite(cell_area_px) & cell_area_px > 0 &
    is.finite(nuclear_area_px) & nuclear_area_px > 0,
  .(
    n_cells = .N,
    median_cell_area_px = median(cell_area_px),
    cell_q10 = as.numeric(quantile(cell_area_px, 0.10)),
    cell_q90 = as.numeric(quantile(cell_area_px, 0.90)),
    median_nuclear_area_px = median(nuclear_area_px),
    nuclear_q10 = as.numeric(quantile(nuclear_area_px, 0.10)),
    nuclear_q90 = as.numeric(quantile(nuclear_area_px, 0.90))
  ),
  by = .(context, G0_mM, hours, ploidy_state)
]
time_all[, context_label := factor(
  as.character(context),
  levels = context_levels,
  labels = unname(context_display[context_levels])
)]
time_all[, G0_label := factor(G0_mM, levels = all_g0, labels = all_g0_labels)]

make_all_time_plot <- function(median_col, q10_col, q90_col, y_label, tag) {
  ggplot(time_all, aes(
    x = hours, y = .data[[median_col]],
    color = ploidy_state, fill = ploidy_state
  )) +
    annotate(
      "rect", xmin = 24, xmax = 48, ymin = -Inf, ymax = Inf,
      fill = "#F2C94C", alpha = 0.16
    ) +
    geom_ribbon(
      aes(ymin = .data[[q10_col]], ymax = .data[[q90_col]]),
      alpha = 0.12, color = NA
    ) +
    geom_line(linewidth = 0.48) +
    facet_grid(rows = vars(context_label), cols = vars(G0_label)) +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = NULL
    ) +
    scale_fill_manual(values = ploidy_colors, guide = "none") +
    scale_x_continuous(breaks = c(0, 48, 96, 144)) +
    scale_y_continuous(
      labels = comma,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(tag = tag, x = "Time (h)", y = y_label) +
    figure_theme +
    theme(
      strip.text.x = element_text(size = 6.5),
      strip.text.y = element_text(size = 6.5),
      panel.spacing = grid::unit(0.35, "lines")
    )
}

f1a <- make_all_time_plot(
  "median_cell_area_px", "cell_q10", "cell_q90",
  "Segmented-cell area (pixels)", "a"
)
f1b <- make_all_time_plot(
  "median_nuclear_area_px", "nuclear_q10", "nuclear_q90",
  "Nuclear area (pixels)", "b"
) + theme(legend.position = "none")

# ---- Figure 2: field-level confluence robustness ----------------------------
highlight_area <- area[
  hours >= 24 & hours <= 48 & G0_mM %in% shared_g0
]
field_totals <- highlight_area[, .(
  field_cell_count = .N,
  field_total_segmented_area_px = sum(cell_area_px)
), by = .(context, G0_mM, hours, image_key)]
field_ploidy <- highlight_area[, .(
  n_group_cells = .N,
  median_cell_area_px = median(cell_area_px),
  median_nuclear_area_px = median(nuclear_area_px)
), by = .(context, G0_mM, hours, image_key, ploidy_state)]
field_summary <- merge(
  field_ploidy, field_totals,
  by = c("context", "G0_mM", "hours", "image_key"),
  all.x = TRUE
)
field_summary[, context_label := factor(
  as.character(context),
  levels = context_levels,
  labels = unname(context_display[context_levels])
)]
field_summary[, G0_label := factor(
  G0_mM, levels = shared_g0,
  labels = c("0.1 mM", "0.5 mM", "1 mM")
)]

make_confluence_plot <- function(x_col, y_col, x_label, y_label, tag) {
  x_axis_labels <- if (grepl("total_segmented", x_col)) {
    function(x) {
      ifelse(
        x >= 1e6,
        paste0(format(x / 1e6, trim = TRUE), "M"),
        ifelse(
          x >= 1e3,
          paste0(format(x / 1e3, trim = TRUE), "k"),
          comma(x)
        )
      )
    }
  } else {
    comma
  }
  ggplot(field_summary, aes(
    x = .data[[x_col]], y = .data[[y_col]], color = ploidy_state
  )) +
    geom_point(alpha = 0.22, size = 0.65) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55) +
    facet_grid(rows = vars(context_label), cols = vars(G0_label), scales = "free_x") +
    scale_color_manual(
      values = ploidy_colors,
      breaks = c("low", "high"),
      labels = c("Low / 2N", "High / 4N"),
      name = NULL
    ) +
    scale_x_log10(labels = x_axis_labels) +
    scale_y_continuous(labels = comma, limits = c(0, NA)) +
    labs(tag = tag, x = x_label, y = y_label) +
    figure_theme +
    theme(
      strip.text = element_text(size = 6.5),
      panel.spacing = grid::unit(0.36, "lines")
    )
}

f2a <- make_confluence_plot(
  "field_cell_count", "median_cell_area_px",
  "Field cell count (log scale)", "Median cell area (pixels)", "a"
)
f2b <- make_confluence_plot(
  "field_cell_count", "median_nuclear_area_px",
  "Field cell count (log scale)", "Median nuclear area (pixels)", "b"
) + theme(legend.position = "none")
f2c <- make_confluence_plot(
  "field_total_segmented_area_px", "median_cell_area_px",
  "Field total segmented area (pixels; log scale)",
  "Median cell area (pixels)", "c"
) + theme(legend.position = "none")
f2d <- make_confluence_plot(
  "field_total_segmented_area_px", "median_nuclear_area_px",
  "Field total segmented area (pixels; log scale)",
  "Median nuclear area (pixels)", "d"
) + theme(legend.position = "none")

# ---- Figure 3: highlighted distributions and nominally same 2N comparison ---
broad_focus <- area[
  hours >= 24 & hours <= 48 & G0_mM %in% shared_g0,
  .(
    median_cell_area_px = median(cell_area_px),
    median_nuclear_area_px = median(nuclear_area_px)
  ),
  by = .(context, ploidy_state, G0_mM)
]
chem_reference <- broad_focus[
  as.character(context) == "M00b chem" & as.character(ploidy_state) == "low",
  .(
    G0_mM,
    chem_cell_area_px = median_cell_area_px,
    chem_nuclear_area_px = median_nuclear_area_px
  )
]
same_2n <- merge(
  broad_focus[as.character(ploidy_state) == "low"],
  chem_reference, by = "G0_mM"
)
same_2n[, cell_ratio_to_chem_2n := median_cell_area_px / chem_cell_area_px]
same_2n[, nuclear_ratio_to_chem_2n :=
  median_nuclear_area_px / chem_nuclear_area_px]
same_2n[, context_label := factor(
  as.character(context),
  levels = context_levels,
  labels = unname(context_display[context_levels])
)]
same_2n[, G0_label := factor(
  G0_mM, levels = shared_g0,
  labels = c("0.1 mM", "0.5 mM", "1 mM")
)]
context_colors <- c(
  "SUM-159-fuse C00" = "#0072B2",
  "SUM-159-fuse I00" = "#56B4E9",
  "SUM-159-chem M00b" = "#555555",
  "Mixture GFP-defined" = "#009E73"
)
make_ratio_plot <- function(value_col, y_label, tag, legend = TRUE) {
  ggplot(same_2n, aes(
    x = G0_label, y = .data[[value_col]],
    color = context_label, group = context_label
  )) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2.2) +
    scale_color_manual(values = context_colors, name = NULL) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.1, suffix = "×"),
      limits = c(0.75, 2.0), breaks = seq(0.75, 2, 0.25)
    ) +
    labs(tag = tag, x = "Initial glucose, G0", y = y_label) +
    figure_theme +
    theme(legend.position = if (legend) "top" else "none")
}

f3a <- cell_hist_plot + labs(tag = "a")
f3b <- nuclear_hist_plot + labs(tag = "b") + theme(legend.position = "none")
f3c <- make_ratio_plot(
  "cell_ratio_to_chem_2n",
  "Low / 2N cell area\n(relative to chem 2N)", "c", TRUE
) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.text = element_text(size = 6),
    legend.key.width = grid::unit(0.23, "in")
  )
f3d <- make_ratio_plot(
  "nuclear_ratio_to_chem_2n",
  "Low / 2N nuclear area\n(relative to chem 2N)", "d", FALSE
)
# ---- Figure 4a-c: cytoplasmic red with cell-size standardization -------------
exact_mask_area_path <- file.path(
  data_dir, "cytoplasmic_cell_area_from_exact_masks.csv.gz"
)
if (!file.exists(exact_mask_area_path)) {
  stop("Missing exact-mask area map: ", exact_mask_area_path)
}
area_map <- fread(exact_mask_area_path)
area_map[, `:=`(
  image_key = as.character(image_key),
  object_id = as.integer(object_id),
  cell_area_px = as.numeric(cell_area_px)
)]
cyto_size <- merge(
  copy(cyto), area_map,
  by = c("image_key", "object_id"),
  all.x = TRUE
)
if (anyNA(cyto_size$cell_area_px)) {
  stop("Failed to attach cell area to all cytoplasmic-red observations.")
}
cyto_size[, log_cell_area := log10(cell_area_px)]
cyto_size[, size_decile := {
  cutpoints <- quantile(
    log_cell_area, probs = seq(0.1, 0.9, by = 0.1),
    na.rm = TRUE, names = FALSE
  )
  findInterval(log_cell_area, cutpoints) + 1L
}, by = assay]
cyto_size[, pair_id := fifelse(
  as.character(context) == "J00 mixture",
  image_key,
  paste(as.character(context), G0_mM, hours, position, sep = "|")
)]

cyto_stratum_image <- cyto_size[, .(
  n_cells = .N,
  stratum_median_bg_z = median(red_cytoplasm_bg_z)
), by = .(
  context, assay, G0_mM, hours, pair_id, image_key,
  ploidy_state, size_decile
)]
cyto_stratum_image <- cyto_stratum_image[n_cells >= 10]
cyto_stratum_pair <- dcast(
  cyto_stratum_image,
  context + assay + G0_mM + hours + pair_id + size_decile ~ ploidy_state,
  value.var = c("stratum_median_bg_z", "n_cells")
)
cyto_stratum_pair <- cyto_stratum_pair[
  is.finite(stratum_median_bg_z_low) &
    is.finite(stratum_median_bg_z_high) &
    n_cells_low >= 10 & n_cells_high >= 10
]
cyto_stratum_pair[, high_minus_low_bg_z :=
  stratum_median_bg_z_high - stratum_median_bg_z_low]
cyto_size_standardized <- cyto_stratum_pair[, .(
  n_shared_size_deciles = .N,
  size_standardized_high_minus_low = mean(high_minus_low_bg_z)
), by = .(context, assay, G0_mM, hours, pair_id)]

size_plot_data <- as_tibble(cyto_size_standardized) %>%
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
cyto_panel_c <- ggplot(
  size_plot_data,
  aes(
    x = context_label,
    y = size_standardized_high_minus_low,
    color = G0_label
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(
    aes(group = context_label),
    width = 0.48, outlier.shape = NA,
    color = "grey35", fill = "grey96", linewidth = 0.38
  ) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 159),
    size = 1.55, alpha = 0.85
  ) +
  scale_color_viridis_d(name = "G0") +
  labs(
    tag = "c", x = NULL,
    y = "High - low matched median\n(size-decile standardized)"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 24, hjust = 1),
    legend.position = "none"
  )

short_context_labels <- setNames(
  c("C00", "I00", "M00b", "Mixture"),
  unname(context_labels[context_levels])
)
cyto_panel_a <- cyto_panel_a +
  facet_wrap(
    vars(context_label),
    nrow = 1,
    labeller = labeller(context_label = short_context_labels)
  ) +
  theme(
    strip.text = element_text(size = 6.2),
    axis.title.y = element_text(size = 7)
  )
cyto_panel_b <- cyto_panel_b +
  scale_x_discrete(labels = short_context_labels) +
  labs(y = "High - low\nmatched median") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
    axis.title.y = element_text(size = 7),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    legend.key.width = grid::unit(0.16, "in")
  )
cyto_panel_c <- cyto_panel_c +
  scale_x_discrete(labels = short_context_labels) +
  labs(y = "High - low\n(size-decile standardized)") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
    axis.title.y = element_text(size = 7)
  )

fwrite(time_all, file.path(data_dir, "v3_all_g0_time_summary.csv"))
fwrite(field_summary, file.path(data_dir, "v3_field_confluence_summary_24_48h.csv"))
fwrite(same_2n, file.path(data_dir, "v3_same_2n_ratios_shared_g0.csv"))
fwrite(cyto_stratum_pair, file.path(data_dir, "v3_cytoplasmic_size_stratum_contrasts.csv"))
fwrite(
  cyto_size_standardized,
  file.path(data_dir, "v3_cytoplasmic_size_standardized_pair_contrasts.csv")
)

# ---- Figure 4d: rebuild microscopy fields directly from raw TIFF data -------
read_tiff_matrix <- function(path, labels = FALSE) {
  image <- suppressWarnings(
    tiff::readTIFF(path, native = FALSE, convert = FALSE)
  )
  if (length(dim(image)) > 2L) image <- image[, , 1L]
  image <- as.matrix(image)
  if (labels) image <- round(image * 65535)
  image
}

crop_matrix <- function(image, y0, x0, size) {
  image[
    seq.int(as.integer(y0) + 1L, as.integer(y0) + as.integer(size)),
    seq.int(as.integer(x0) + 1L, as.integer(x0) + as.integer(size)),
    drop = FALSE
  ]
}

normalize_channel <- function(image, limits) {
  if (any(!is.finite(limits)) || limits[[2]] <= limits[[1]]) {
    return(matrix(0, nrow(image), ncol(image)))
  }
  normalized <- (image - limits[[1]]) / (limits[[2]] - limits[[1]])
  normalized[normalized < 0] <- 0
  normalized[normalized > 1] <- 1
  dim(normalized) <- dim(image)
  normalized
}

mask_boundary <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  boundary <- matrix(FALSE, nr, nc)
  boundary[-nr, ] <- boundary[-nr, ] |
    (mask[-nr, ] > 0 & mask[-nr, ] != mask[-1L, ])
  boundary[-1L, ] <- boundary[-1L, ] |
    (mask[-1L, ] > 0 & mask[-1L, ] != mask[-nr, ])
  boundary[, -nc] <- boundary[, -nc] |
    (mask[, -nc] > 0 & mask[, -nc] != mask[, -1L])
  boundary[, -1L] <- boundary[, -1L] |
    (mask[, -1L] > 0 & mask[, -1L] != mask[, -nc])
  boundary
}

overlay_boundaries <- function(rgb, mask, object_states, alpha = 0.95) {
  boundary <- mask_boundary(mask)
  state <- matrix("unknown", nrow(mask), ncol(mask))
  if (length(object_states) == 1L && is.null(names(object_states))) {
    state[mask > 0] <- object_states[[1]]
  } else {
    object_ids <- as.integer(names(object_states))
    matched <- match(as.vector(mask), object_ids)
    mapped <- unname(object_states)[matched]
    mapped[is.na(mapped)] <- "unknown"
    state[] <- mapped
  }
  colors <- list(
    low = c(0x00, 0xD9, 0xFF) / 255,
    high = c(0xFF, 0xB0, 0x00) / 255,
    unknown = rep(0.75, 3)
  )
  for (state_name in names(colors)) {
    pixels <- boundary & state == state_name
    if (!any(pixels)) next
    for (channel in seq_len(3L)) {
      layer <- rgb[, , channel]
      layer[pixels] <- (1 - alpha) * layer[pixels] +
        alpha * colors[[state_name]][[channel]]
      rgb[, , channel] <- layer
    }
  }
  rgb
}

make_field_grob <- function(
  rgb, row_label = NULL, missing_label = FALSE
) {
  array_grob <- grid::rasterGrob
  children <- list(
    array_grob(
      rgb, x = 0, y = 0, width = 1, height = 1,
      just = c("left", "bottom"), interpolate = FALSE
    )
  )
  if (missing_label) {
    children <- c(
      children,
      list(textGrob(
        "not acquired", x = 0.5, y = 0.5,
        gp = gpar(col = "white", fontsize = 7.5)
      ))
    )
  }
  if (!is.null(row_label)) {
    children <- c(
      children,
      list(
        rectGrob(
          x = 0.015, y = 0.985, width = 0.67, height = 0.18,
          just = c("left", "top"),
          gp = gpar(fill = "white", col = NA, alpha = 0.82)
        ),
        textGrob(
          row_label, x = 0.03, y = 0.965,
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 6.5, fontface = "bold")
        )
      )
    )
  }
  do.call(grobTree, children)
}

build_microscopy_panel <- function() {
  field_manifest <- fread(
    file.path(data_dir, "multimodal_field_manifest.csv"),
    encoding = "UTF-8"
  )
  if (nrow(field_manifest) != 7L) {
    stop("Expected seven selected microscopy fields; found ", nrow(field_manifest))
  }

  fields <- lapply(seq_len(nrow(field_manifest)), function(index) {
    row <- field_manifest[index]
    size <- as.integer(row$crop_size_pixels)
    mask <- crop_matrix(
      read_tiff_matrix(row$cell_mask_path, labels = TRUE),
      row$crop_y0, row$crop_x0, size
    )
    gfp_available <- !is.na(row$gfp_path) &&
      row$gfp_path != "not acquired" &&
      file.exists(row$gfp_path)
    list(
      row = row,
      mask = mask,
      phase = crop_matrix(
        read_tiff_matrix(row$phase_path),
        row$crop_y0, row$crop_x0, size
      ),
      red = crop_matrix(
        read_tiff_matrix(row$red_path),
        row$crop_y0, row$crop_x0, size
      ),
      nir = crop_matrix(
        read_tiff_matrix(row$nir_path),
        row$crop_y0, row$crop_x0, size
      ),
      gfp = if (gfp_available) crop_matrix(
        read_tiff_matrix(row$gfp_path),
        row$crop_y0, row$crop_x0, size
      ) else NULL
    )
  })

  assay_names <- unique(field_manifest$assay)
  red_limits <- setNames(vector("list", length(assay_names)), assay_names)
  for (assay_name in assay_names) {
    values <- unlist(lapply(
      fields,
      function(field) {
        if (field$row$assay == assay_name) as.vector(field$red) else numeric()
      }
    ), use.names = FALSE)
    red_limits[[assay_name]] <- as.numeric(
      quantile(values, c(0.01, 0.96), na.rm = TRUE, names = FALSE)
    )
  }

  titles <- c(
    "Phase", "mCherry / red", "NIR", "GFP",
    "Composite +\ncell outlines"
  )
  title_grobs <- lapply(
    titles,
    function(label) textGrob(label, gp = gpar(fontsize = 8.5))
  )
  field_grobs <- list()
  for (field_index in seq_along(fields)) {
    field <- fields[[field_index]]
    row <- field$row
    if (row$ploidy_state == "mixed") {
      state_rows <- cyto[
        image_key == row$image_key,
        .(object_id = as.integer(object_id), ploidy_state = as.character(ploidy_state))
      ]
      object_states <- setNames(
        state_rows$ploidy_state,
        as.character(state_rows$object_id)
      )
    } else {
      object_states <- as.character(row$ploidy_state)
    }

    phase_norm <- normalize_channel(
      field$phase,
      as.numeric(quantile(field$phase, c(0.01, 0.99), na.rm = TRUE))
    )
    red_norm <- normalize_channel(field$red, red_limits[[row$assay]])
    nir_norm <- normalize_channel(
      field$nir,
      as.numeric(quantile(field$nir, c(0.01, 0.995), na.rm = TRUE))
    )
    gfp_norm <- if (is.null(field$gfp)) {
      matrix(0, nrow(field$red), ncol(field$red))
    } else {
      normalize_channel(
        field$gfp,
        as.numeric(quantile(field$gfp, c(0.01, 0.995), na.rm = TRUE))
      )
    }

    grayscale <- list(phase_norm, red_norm, nir_norm, gfp_norm)
    for (column_index in seq_len(4L)) {
      gray_rgb <- array(
        rep(grayscale[[column_index]], 3L),
        dim = c(nrow(field$mask), ncol(field$mask), 3L)
      )
      field_grobs[[length(field_grobs) + 1L]] <- make_field_grob(
        overlay_boundaries(gray_rgb, field$mask, object_states),
        row_label = if (column_index == 1L) row$row_label else NULL,
        missing_label = column_index == 4L && is.null(field$gfp)
      )
    }

    base <- 0.22 * phase_norm
    composite <- array(0, dim = c(nrow(base), ncol(base), 3L))
    composite[, , 1L] <- base + red_norm
    composite[, , 2L] <- base + gfp_norm
    composite[, , 3L] <- base + 0.72 * red_norm + nir_norm
    composite[composite < 0] <- 0
    composite[composite > 1] <- 1
    field_grobs[[length(field_grobs) + 1L]] <- make_field_grob(
      overlay_boundaries(composite, field$mask, object_states)
    )
  }

  layout_matrix <- matrix(seq_len(40L), nrow = 8L, ncol = 5L, byrow = TRUE)
  body <- arrangeGrob(
    grobs = c(title_grobs, field_grobs),
    layout_matrix = layout_matrix,
    heights = c(0.25, rep(1, 7)),
    widths = rep(1, 5)
  )
  arrangeGrob(
    body,
    top = textGrob(
      "d", x = 0.005, just = "left",
      gp = gpar(fontsize = 12, fontface = "bold")
    ),
    bottom = textGrob(
      paste(
        "All fields: 1 mM G0, 512 × 512 pixels.",
        "Cyan outlines = low / 2N; amber outlines = high / 4N."
      ),
      gp = gpar(fontsize = 7.5)
    )
  )
}

f4d <- build_microscopy_panel()

# ---- Shared panel registry and audit exports --------------------------------
panel_specs <- data.table(
  figure = c(
    rep("SUM159 evidence 1", 2),
    rep("SUM159 evidence 2", 4),
    rep("SUM159 evidence 3", 4),
    rep("SUM159 evidence 4", 4)
  ),
  panel = c(
    "a", "b",
    "a", "b", "c", "d",
    "a", "b", "c", "d",
    "a", "b", "c", "d"
  ),
  width_in = c(
    7, 7,
    rep(1037 / 300, 4),
    3.42, 3.42, 3.42, 3.42,
    2.75, 2.02, 2.02, 7
  ),
  height_in = c(
    4.55, 4.55,
    rep(4.568, 4),
    5.15, 5.15, 2.75, 2.75,
    2.02, 2.02, 2.02, 7.05
  ),
  order = seq_len(14L)
)
panel_specs[, subpanel_png := file.path(
  subpanel_dir,
  sprintf(
    "sum159_evidence_%s%s.png",
    sub(".* ", "", figure),
    panel
  )
)]
panel_grobs <- list(
  "SUM159 evidence 1::a" = ggplotGrob(f1a),
  "SUM159 evidence 1::b" = ggplotGrob(f1b),
  "SUM159 evidence 2::a" = ggplotGrob(f2a),
  "SUM159 evidence 2::b" = ggplotGrob(f2b),
  "SUM159 evidence 2::c" = ggplotGrob(f2c),
  "SUM159 evidence 2::d" = ggplotGrob(f2d),
  "SUM159 evidence 3::a" = ggplotGrob(f3a),
  "SUM159 evidence 3::b" = ggplotGrob(f3b),
  "SUM159 evidence 3::c" = ggplotGrob(f3c),
  "SUM159 evidence 3::d" = ggplotGrob(f3d),
  "SUM159 evidence 4::a" = ggplotGrob(cyto_panel_a),
  "SUM159 evidence 4::b" = ggplotGrob(cyto_panel_b),
  "SUM159 evidence 4::c" = ggplotGrob(cyto_panel_c),
  "SUM159 evidence 4::d" = f4d
)

save_grob_png <- function(grob, path, width, height, dpi) {
  png(
    path, width = width, height = height, units = "in", res = dpi,
    type = "cairo", bg = "white"
  )
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  grid.draw(grob)
}

if (phase == "subpanels") {
  audit_dpi <- 300
  for (row_index in seq_len(nrow(panel_specs))) {
    row <- panel_specs[row_index]
    key <- paste(row$figure, row$panel, sep = "::")
    save_grob_png(
      panel_grobs[[key]], row$subpanel_png,
      row$width_in, row$height_in, audit_dpi
    )
  }
  dimensions <- panel_specs[, .(
    figure,
    panel,
    subpanel_png,
    width_px = as.integer(round(width_in * audit_dpi)),
    height_px = as.integer(round(height_in * audit_dpi)),
    width_in,
    height_in,
    order
  )]
  fwrite(dimensions, file.path(layout_dir, "subpanel_dimensions.csv"))
  message("Regenerated 14 audit subpanels and wrote subpanel dimensions.")
  quit(save = "no", status = 0)
}

# ---- Final composites from live grobs and optimizer coordinates -------------
layout_path <- file.path(layout_dir, "layout_plan.csv")
if (!file.exists(layout_path)) {
  stop("Missing layout plan; run the bundled layout optimizer first: ", layout_path)
}
layout_plan <- fread(layout_path)

compose_from_plan <- function(figure_id, output_name, dpi = 450) {
  plan <- layout_plan[figure == figure_id]
  if (nrow(plan) == 0L) stop("No layout rows for ", figure_id)
  width <- unique(plan$layout_width_in)
  height <- unique(plan$layout_height_in)
  if (length(width) != 1L || length(height) != 1L) {
    stop("Inconsistent layout dimensions for ", figure_id)
  }
  output_path <- file.path(v3_dir, output_name)
  png(
    output_path, width = width, height = height,
    units = "in", res = dpi, type = "cairo", bg = "white"
  )
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  for (row_index in seq_len(nrow(plan))) {
    row <- plan[row_index]
    key <- paste(figure_id, row$panel, sep = "::")
    pushViewport(viewport(
      x = row$x_npc,
      y = row$y_npc,
      width = row$width_npc,
      height = row$height_npc,
      just = c("left", "bottom"),
      name = sprintf("%s_%s", gsub(" ", "_", figure_id), row$panel)
    ))
    grid.draw(panel_grobs[[key]])
    upViewport()
  }
}

compose_from_plan("SUM159 evidence 1", "figure_1_all_timecourses.png")
compose_from_plan("SUM159 evidence 2", "figure_2_confluence_robustness.png")
compose_from_plan(
  "SUM159 evidence 3",
  "figure_3_focused_distributions_and_same_2n.png"
)
compose_from_plan(
  "SUM159 evidence 4",
  "figure_4_cytoplasmic_signal_and_multimodal_fields.png"
)

message("Regenerated all four final composites locally from live panel objects.")
