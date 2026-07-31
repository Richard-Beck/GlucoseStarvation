#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(scales)
  library(tibble)
})

root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v5_redraft",
  "figure_generation", "FG1_measurement_foundation", "drafting_v2"
)
round2_dir <- file.path(root, "initial_subpanels", "measurement_provenance_round2")
refined_dir <- file.path(root, "refined_subpanels")
worker_dir <- file.path(root, "worker_notes")
dir.create(round2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(refined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(worker_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  annotation_stack = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90.tif"),
  object_labels = file.path("data", "image_processing_runs", "run_20260324_233122", "objects_labels_90.tif"),
  annotation_manifest = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_manifest.csv"),
  annotation_points = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_Results.csv"),
  object_predictions = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_object_predictions.csv"),
  frame_counts = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_frame_counts.csv"),
  provenance = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_provenance.json"),
  prior_options_manifest = file.path(root, "initial_subpanels", "measurement_provenance_options", "fg1_figure1a_options_manifest.csv")
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) {
  stop("Missing required Figure 1A round-2 inputs: ", paste(missing, collapse = ", "), call. = FALSE)
}
if (!requireNamespace("tiff", quietly = TRUE)) {
  stop("The R package 'tiff' is required for annotation and label-mask TIFF rendering.", call. = FALSE)
}

read_tbl <- function(path) {
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), .name_repair = "unique_quiet")
}

read_tiff_all <- function(path, raw_labels = FALSE) {
  if (raw_labels) {
    out <- tryCatch(
      suppressWarnings(tiff::readTIFF(path, all = TRUE, as.is = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  suppressWarnings(tiff::readTIFF(path, all = TRUE))
}

annotation_pages <- read_tiff_all(paths$annotation_stack)
label_pages <- read_tiff_all(paths$object_labels, raw_labels = TRUE)
if (!is.list(annotation_pages)) annotation_pages <- list(annotation_pages)
if (!is.list(label_pages)) label_pages <- list(label_pages)

manifest <- read_tbl(paths$annotation_manifest)
points <- read_tbl(paths$annotation_points)
object_predictions <- read_tbl(paths$object_predictions)
frame_counts <- read_tbl(paths$frame_counts)
prior_options_manifest <- read_tbl(paths$prior_options_manifest)

line_palette <- c(
  alive = "#0f9e9a",
  dead = "#c73e5a",
  junk = "#7b61ff"
)
manual_color <- "#ffd33d"
target_color <- "#f2c94c"

robust01 <- function(x) {
  arr <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
  finite <- is.finite(arr)
  if (!any(finite)) return(matrix(0, nrow = nrow(arr), ncol = ncol(arr)))
  vals <- arr[finite]
  lo <- as.numeric(stats::quantile(vals, 0.01, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(vals, 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    lo <- min(vals, na.rm = TRUE)
    hi <- max(vals, na.rm = TRUE)
  }
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(matrix(0, nrow = nrow(arr), ncol = ncol(arr)))
  pmin(pmax((arr - lo) / (hi - lo), 0), 1)
}

frame_channels <- function(frame) {
  first <- annotation_pages[[1]]
  if (length(dim(first)) == 3 && length(annotation_pages) >= frame) {
    arr <- annotation_pages[[frame]]
    if (dim(arr)[[3]] >= 3) {
      return(list(phase = arr[, , 1], dead = arr[, , 2], alive = arr[, , 3]))
    }
  }
  idx <- ((frame - 1L) * 3L + 1L):(frame * 3L)
  if (max(idx) > length(annotation_pages)) {
    stop("Could not read three annotation-stack channels for frame ", frame, call. = FALSE)
  }
  list(
    phase = annotation_pages[[idx[[1]]]],
    dead = annotation_pages[[idx[[2]]]],
    alive = annotation_pages[[idx[[3]]]]
  )
}

annotation_frame_rgb <- function(frame) {
  channels <- frame_channels(frame)
  phase <- robust01(channels$phase)
  dead <- robust01(channels$dead)
  alive <- robust01(channels$alive)
  base <- 0.58 * phase
  rgb <- array(0, dim = c(nrow(phase), ncol(phase), 3))
  rgb[, , 1] <- pmin(base + 0.75 * dead, 1)
  rgb[, , 2] <- pmin(base + 0.82 * alive, 1)
  rgb[, , 3] <- pmin(base + 0.12 * dead, 1)
  as.raster(rgb)
}

label_frame <- function(frame, max_object_id) {
  first <- label_pages[[1]]
  labels <- NULL
  if (length(label_pages) >= frame && length(dim(label_pages[[frame]])) == 2) {
    labels <- label_pages[[frame]]
  } else if (length(label_pages) == 1 && length(dim(first)) == 3) {
    if (dim(first)[[3]] >= frame) {
      labels <- first[, , frame]
    } else if (dim(first)[[1]] >= frame) {
      labels <- first[frame, , ]
    }
  }
  if (is.null(labels)) {
    stop("Could not extract label-mask frame ", frame, " from ", paths$object_labels, call. = FALSE)
  }
  labels <- matrix(as.numeric(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (max(labels, na.rm = TRUE) <= 1 && max_object_id > 1) {
    labels <- round(labels * 65535)
  } else {
    labels <- round(labels)
  }
  labels
}

label_boundary <- function(labels) {
  nr <- nrow(labels)
  nc <- ncol(labels)
  boundary <- matrix(FALSE, nr, nc)
  if (nr > 1) {
    boundary[-nr, ] <- boundary[-nr, ] | labels[-nr, ] != labels[-1, ]
    boundary[-1, ] <- boundary[-1, ] | labels[-1, ] != labels[-nr, ]
  }
  if (nc > 1) {
    boundary[, -nc] <- boundary[, -nc] | labels[, -nc] != labels[, -1]
    boundary[, -1] <- boundary[, -1] | labels[, -1] != labels[, -nc]
  }
  boundary & labels > 0
}

production_state <- function(predicted_label_name) {
  case_when(
    predicted_label_name == "alive" ~ "alive",
    predicted_label_name == "dead" ~ "dead",
    predicted_label_name == "junk" ~ "junk",
    TRUE ~ "dead"
  )
}

outline_df_for_frame <- function(frame, calls, crop_h) {
  labels <- label_frame(frame, max(calls$object_id, na.rm = TRUE))
  calls <- calls %>%
    mutate(
      object_id = as.integer(.data$object_id),
      production_state = production_state(.data$predicted_label_name)
    ) %>%
    distinct(.data$object_id, .data$production_state)
  boundary <- label_boundary(labels) & labels %in% calls$object_id
  idx <- which(boundary, arr.ind = TRUE)
  if (!nrow(idx)) {
    return(tibble(x = numeric(), y = numeric(), production_state = character()))
  }
  object_id <- as.integer(labels[idx])
  tibble(
    x = idx[, "col"] - 0.5,
    y = crop_h - (idx[, "row"] - 0.5),
    object_id = object_id
  ) %>%
    left_join(calls, by = "object_id") %>%
    filter(!is.na(.data$production_state))
}

fmt_num <- function(x) {
  out <- formatC(as.numeric(x), format = "fg", digits = 4)
  sub("\\.0$", "", out)
}

tile_label <- function(row) {
  paste0(row$cellLine, " ", row$ploidy, "\n", fmt_num(row$G0), " mM, ", fmt_num(row$hours), " h")
}

make_tile <- function(frame) {
  m <- manifest %>% filter(.data$stack_index_1based == frame) %>% slice(1)
  calls <- object_predictions %>% filter(.data$frame == !!frame)
  if (nrow(m) != 1 || nrow(calls) == 0) {
    stop("Missing tile input for frame ", frame, call. = FALSE)
  }
  crop_h <- as.numeric(m$crop_y1 - m$crop_y0)
  crop_w <- as.numeric(m$crop_x1 - m$crop_x0)
  target <- tibble(
    xmin = as.numeric(m$target_x0 - m$crop_x0),
    xmax = as.numeric(m$target_x1 - m$crop_x0),
    ymin = crop_h - as.numeric(m$target_y1 - m$crop_y0),
    ymax = crop_h - as.numeric(m$target_y0 - m$crop_y0)
  )
  point_frame <- points %>%
    filter(.data$Frame == frame) %>%
    transmute(x = .data$X, y = crop_h - .data$Y)
  outlines <- outline_df_for_frame(frame, calls, crop_h)
  label_text <- tile_label(m)

  ggplot() +
    annotation_raster(annotation_frame_rgb(frame), xmin = 0, xmax = crop_w, ymin = 0, ymax = crop_h) +
    geom_tile(
      data = outlines,
      aes(x = .data$x, y = .data$y, fill = .data$production_state),
      width = 1.15,
      height = 1.15,
      alpha = 0.92,
      show.legend = FALSE
    ) +
    geom_point(
      data = point_frame,
      aes(.data$x, .data$y),
      color = "#1a1a1a",
      fill = manual_color,
      shape = 21,
      stroke = 0.32,
      size = 1.08,
      alpha = 1
    ) +
    geom_rect(data = target, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = NA, color = target_color, linewidth = 0.34
    ) +
    annotate("rect", xmin = 0, xmax = crop_w, ymin = 0, ymax = crop_h, fill = NA, color = "#3b424d", linewidth = 0.16) +
    annotate(
      "label",
      x = 5,
      y = crop_h - 5,
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 1.55,
      lineheight = 0.86,
      fontface = "bold",
      color = "white",
      fill = alpha("#20242b", 0.70),
      label.size = 0,
      label.r = unit(0.025, "lines")
    ) +
    scale_fill_manual(values = line_palette, guide = "none") +
    coord_fixed(xlim = c(0, crop_w), ylim = c(0, crop_h), expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(0.4, 0.4, 0.4, 0.4),
      panel.background = element_rect(fill = "#20242b", color = NA)
    )
}

legend_panel <- function() {
  ggplot() +
    annotate("segment", x = 0.025, xend = 0.080, y = 0.52, yend = 0.52, color = line_palette[["alive"]], linewidth = 1.2) +
    annotate("text", x = 0.090, y = 0.52, label = "alive", hjust = 0, vjust = 0.5, size = 2.2) +
    annotate("segment", x = 0.205, xend = 0.260, y = 0.52, yend = 0.52, color = line_palette[["dead"]], linewidth = 1.2) +
    annotate("text", x = 0.270, y = 0.52, label = "dead", hjust = 0, vjust = 0.5, size = 2.2) +
    annotate("segment", x = 0.385, xend = 0.440, y = 0.52, yend = 0.52, color = line_palette[["junk"]], linewidth = 1.2) +
    annotate("text", x = 0.450, y = 0.52, label = "junk", hjust = 0, vjust = 0.5, size = 2.2) +
    annotate("point", x = 0.565, y = 0.52, shape = 21, size = 2.5, stroke = 0.35, color = "#1a1a1a", fill = manual_color) +
    annotate("text", x = 0.590, y = 0.52, label = "manual alive", hjust = 0, vjust = 0.5, size = 2.2) +
    annotate("rect", xmin = 0.795, xmax = 0.850, ymin = 0.34, ymax = 0.70, fill = NA, color = target_color, linewidth = 0.55) +
    annotate("text", x = 0.860, y = 0.52, label = "target", hjust = 0, vjust = 0.5, size = 2.2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 1, 0, 1))
}

ploidy_map <- tribble(
  ~cellLine, ~ploidy_group, ~ploidy,
  "MCF10A", "low", "2N",
  "MCF10A", "high", "4N",
  "MDA-MB-231", "low", "parental",
  "MDA-MB-231", "high", "3N",
  "SNU668", "low", "low",
  "SNU668", "high", "high",
  "SUM-159-fuse", "low", "2N",
  "SUM-159-fuse", "high", "4N"
) %>%
  mutate(
    line_order = match(.data$cellLine, c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-fuse")),
    ploidy_order = match(.data$ploidy_group, c("low", "high"))
  )

selected_frames <- bind_rows(lapply(seq_len(nrow(ploidy_map)), function(i) {
  line <- ploidy_map$cellLine[[i]]
  ploidy_value <- ploidy_map$ploidy[[i]]
  frame <- frame_counts %>%
    filter(
      .data$cellLine == line,
      .data$ploidy == ploidy_value,
      .data$glucose_bin == "low",
      .data$time_bin == "early"
    ) %>%
    arrange(.data$hours, .data$G0, .data$frame) %>%
    slice(1) %>%
    pull(.data$frame)
  ploidy_map[i, ] %>%
    mutate(frame = if (length(frame)) as.integer(frame[[1]]) else NA_integer_)
})) %>%
  arrange(.data$ploidy_order, .data$line_order)

missing_pairs <- selected_frames %>% filter(is.na(.data$frame))
expansion_feasible <- nrow(missing_pairs) == 0

if (!expansion_feasible) {
  fallback_frames <- c(13L, 22L, 40L, 85L)
  selected_frames <- frame_counts %>%
    filter(.data$frame %in% fallback_frames) %>%
    transmute(
      cellLine = .data$cellLine,
      ploidy_group = "high",
      ploidy = .data$ploidy,
      line_order = match(.data$cellLine, c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-fuse")),
      ploidy_order = 1L,
      frame = .data$frame
    ) %>%
    arrange(.data$line_order)
}

frames <- selected_frames$frame
tile_plots <- lapply(frames, make_tile)

if (length(tile_plots) == 8) {
  body <- wrap_plots(tile_plots, ncol = 4)
  plot_height <- 3.78
  layout_note <- "Eight-image low/high-ploidy option: four cell lines by two ploidy rows."
} else {
  body <- wrap_plots(tile_plots, ncol = 4)
  plot_height <- 1.98
  layout_note <- "Fallback compact four-line option: high/low expansion was not feasible from available exports."
}

round2_plot <- body / legend_panel() +
  plot_layout(heights = c(1, 0.055)) &
  theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4))

standalone_path <- file.path(round2_dir, "fg1_figure1a_round2_low_high_ploidy_outlines.png")
contact_path <- file.path(refined_dir, "fg1_figure1a_round2_contact_sheet.png")
ggsave(standalone_path, plot = round2_plot, width = 7.10, height = plot_height, units = "in", dpi = 360, bg = "white", limitsize = FALSE)
ggsave(contact_path, plot = round2_plot, width = 7.10, height = plot_height, units = "in", dpi = 360, bg = "white", limitsize = FALSE)

frame_summary <- frame_counts %>%
  filter(.data$frame %in% frames) %>%
  select(all_of(c(
    "frame", "image_key", "cellLine", "ploidy", "G0", "glucose_bin",
    "time_bin", "hours", "scored_objects", "manual_alive_objects",
    "predicted_alive_count", "predicted_dead_count", "predicted_junk_count"
  ))) %>%
  left_join(
    selected_frames %>% select(all_of(c("frame", "ploidy_group", "line_order", "ploidy_order"))),
    by = "frame"
  ) %>%
  arrange(.data$ploidy_order, .data$line_order)

manifest_out <- tibble(
  option_id = "round2_low_high_ploidy_outlines",
  filename = basename(standalone_path),
  output_png = standalone_path,
  refined_contact_sheet = contact_path,
  frames = paste(frames, collapse = ";"),
  overlay = "outlines",
  layout = if (length(tile_plots) == 8) "4 columns x 2 rows plus compact legend" else "4 columns x 1 row plus compact legend",
  rationale = paste(
    "Refines preferred round-1 option 3 using the same production-classifier validation stack,",
    "keeps the liked outline style, reduces gutters/margins, shortens legend labels,",
    "and expands to low/high ploidy when those frames are present."
  ),
  expansion_feasible = expansion_feasible,
  directive_ids = "FG1A-R2-D01;FG1A-R2-D02;FG1A-R2-D03;FG1A-R2-D04;FG1A-R2-D05",
  data_source = paste(
    paths$object_predictions,
    paths$annotation_points,
    paths$annotation_stack,
    paths$object_labels,
    paths$frame_counts,
    sep = ";"
  )
)
write.csv(manifest_out, file.path(round2_dir, "fg1_figure1a_round2_manifest.csv"), row.names = FALSE)
write.csv(frame_summary, file.path(round2_dir, "fg1_figure1a_round2_frame_summary.csv"), row.names = FALSE)

prior_option3 <- prior_options_manifest %>%
  filter(.data$option_id == "option03") %>%
  slice(1)
prior_frames <- if (nrow(prior_option3)) prior_option3$frames[[1]] else "13;22;40;85"

feasibility_text <- if (expansion_feasible) {
  "Low/high-ploidy expansion was feasible from the existing frame-count and object-prediction exports; no new image processing was needed."
} else {
  paste0(
    "Low/high-ploidy expansion was blocked from existing exports. Missing pairs: ",
    paste(paste(missing_pairs$cellLine, missing_pairs$ploidy_group, missing_pairs$ploidy, sep = "/"), collapse = "; "),
    ". The worker generated the compact four-line option-3-style fallback instead."
  )
}

notes <- c(
  "# Figure 1A round-2 compact option notes",
  "",
  "Scope: targeted round-2 Figure 1A option refinement only. Existing final figures, report manifest, review report, feedback intake, and prior scripts were not edited.",
  "",
  "## Source continuity",
  "",
  paste0("- Starting point: prior option 3 frames `", prior_frames, "` from `initial_subpanels/measurement_provenance_options/fg1_figure1a_options_manifest.csv`."),
  "- Data/provenance reused from the current option worker: 90-frame annotation stack, object label masks, manual annotation points, production classifier object predictions, and frame-count metadata.",
  "- Object outlines are still reconstructed from `objects_labels_90.tif` joined to production predictions by `frame` and `object_id`.",
  "",
  "## Round-2 change",
  "",
  paste0("- ", layout_note),
  "- Uses the liked outline style from options 2-4 and does not regenerate the disliked box-style option 1.",
  "- Reduces white space through 0.4-point tile margins, a 4-column layout, no option header, and a short legend strip.",
  "- Shortens legend labels to `alive`, `dead`, `junk`, `manual alive`, and `target`.",
  paste0("- ", feasibility_text),
  "",
  "## Generated files",
  "",
  paste0("- Standalone option: `", standalone_path, "`."),
  paste0("- Refined contact sheet: `", contact_path, "`."),
  paste0("- Manifest: `", file.path(round2_dir, "fg1_figure1a_round2_manifest.csv"), "`."),
  paste0("- Frame summary: `", file.path(round2_dir, "fg1_figure1a_round2_frame_summary.csv"), "`."),
  "",
  "## Visual QC",
  "",
  "- The rendered PNG was visually inspected after generation.",
  "- No figure-level title, subtitle, caption, or `Figure 1` header is present.",
  "- The eight image panels are in reading order: low-ploidy row, then high-ploidy row; columns are MCF10A, MDA-MB-231, SNU668, and SUM-159-fuse.",
  "- Labels and legend are short enough to fit at the rendered size; no obvious clipping or overlap was observed.",
  "- The option remains dense by design, but the gutters and legend are substantially tighter than the prior option contact sheet."
)
writeLines(notes, file.path(round2_dir, "notes.md"))

coverage <- c(
  "# FG1 Figure 1A round-2 coverage",
  "",
  "## Directive coverage",
  "",
  "| directive_id | request | status | output/disposition |",
  "|---|---|---|---|",
  "| FG1A-R2-D01 | Drop/disfavor option 1 and keep the liked option 2-4 outline plotting style. | addressed | The worker generates only outline-based panels; no box-style option 1 is regenerated. |",
  "| FG1A-R2-D02 | Prefer option 3 because it shows four cell lines. | addressed | The refined option preserves the same four-cell-line scope: MCF10A, MDA-MB-231, SNU668, and SUM-159-fuse. |",
  "| FG1A-R2-D03 | Reduce white space and move image panels closer together. | addressed | `fg1_figure1a_round2_low_high_ploidy_outlines.png` and `fg1_figure1a_round2_contact_sheet.png` use 0.4-point tile margins, a 4-column layout, and a compact legend strip. |",
  "| FG1A-R2-D04 | Shorten in-image legend labels to brief terms such as alive/dead/junk and manual alive. | addressed | Legend labels are `alive`, `dead`, `junk`, `manual alive`, and `target`. |",
  paste0("| FG1A-R2-D05 | If feasible, show low and high ploidy so the four-cell-line option becomes eight images. | ", if (expansion_feasible) "addressed" else "blocked", " | ", if (expansion_feasible) "Eight low-glucose early-time images were generated from existing exports without new processing." else feasibility_text, " |"),
  "",
  "## Frame selection",
  "",
  paste0("- Frames: `", paste(frames, collapse = ";"), "`."),
  "- Selection rule: one low-glucose, early-time validation frame for each requested cell line and low/high ploidy pairing, using the existing `fg1_production_classifier_frame_counts.csv` export.",
  "",
  "## Files",
  "",
  paste0("- `", standalone_path, "`"),
  paste0("- `", contact_path, "`"),
  paste0("- `", file.path(round2_dir, "fg1_figure1a_round2_manifest.csv"), "`"),
  paste0("- `", file.path(round2_dir, "fg1_figure1a_round2_frame_summary.csv"), "`"),
  paste0("- `", file.path(round2_dir, "notes.md"), "`"),
  "",
  "## Scope and caveats",
  "",
  "- No classifier retraining, segmentation rerun, annotation rebuild, Stan/model-free/posterior rerun, report regeneration, or recommended figure regeneration was performed.",
  "- Existing final figures, `report_manifest.csv`, `review_report.html`, `feedback_intake.md`, and existing scripts were not edited.",
  "- `docs/project_map.txt` was read as instructed; this narrowly scoped review-option addition does not need a project-map update."
)
writeLines(coverage, file.path(worker_dir, "fg1_figure1a_round2_coverage.md"))

cat("Wrote FG1 Figure 1A round-2 option to ", standalone_path, "\n", sep = "")
cat("Wrote FG1 Figure 1A round-2 contact sheet to ", contact_path, "\n", sep = "")
