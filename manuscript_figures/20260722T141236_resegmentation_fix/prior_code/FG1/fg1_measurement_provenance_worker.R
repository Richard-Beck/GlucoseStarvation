#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(jsonlite)
  library(patchwork)
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
  annotation_stack = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90.tif"),
  annotation_manifest = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_manifest.csv"),
  annotation_points = file.path("data", "image_processing_runs", "run_20260324_233122", "annotation_stack_90_Results.csv"),
  production_calls = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_selected_object_calls.csv"),
  production_frames = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_frame_counts.csv"),
  object_predictions = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_object_predictions.csv"),
  frame_count_metrics = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_frame_count_metrics.csv"),
  headline = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_headline_metrics.csv"),
  object_metrics = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_object_metrics.csv"),
  output_manifest = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_output_manifest.csv"),
  provenance = file.path("data", "report_exports", "v5_fg1_production_classifier_validation", "fg1_production_classifier_provenance.json"),
  prior_frame_counts = file.path("data", "report_exports", "manuscript_draft_figures_v3", "wp1_drafting", "panels", "alive_count_validation", "wp1_alive_count_validation_frame_counts.csv")
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) {
  stop("Missing required production/provenance inputs: ", paste(missing, collapse = ", "), call. = FALSE)
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

theme_fg1 <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, color = "grey88"),
      strip.background = element_rect(fill = "grey94", color = "grey75", linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.margin = margin(4, 4, 4, 4)
    )
}

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

annotation_frame_rgb <- function(frame) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("The R package 'tiff' is required to render annotation-stack overlays.", call. = FALSE)
  }
  imgs <- suppressWarnings(tiff::readTIFF(paths$annotation_stack, all = TRUE))
  if (!is.list(imgs)) imgs <- list(imgs)
  first <- imgs[[1]]
  if (length(dim(first)) == 3 && dim(first)[[3]] >= 3) {
    phase <- first[, , 1]
    dead <- first[, , 2]
    alive <- first[, , 3]
  } else {
    idx <- ((frame - 1L) * 3L + 1L):(frame * 3L)
    if (max(idx) > length(imgs)) stop("Could not read three channels for frame ", frame, call. = FALSE)
    phase <- imgs[[idx[[1]]]]
    dead <- imgs[[idx[[2]]]]
    alive <- imgs[[idx[[3]]]]
  }
  phase <- robust01(phase)
  dead <- robust01(dead)
  alive <- robust01(alive)
  base <- 0.58 * phase
  rgb <- array(0, dim = c(nrow(phase), ncol(phase), 3))
  rgb[, , 1] <- pmin(base + 0.75 * dead, 1)
  rgb[, , 2] <- pmin(base + 0.82 * alive, 1)
  rgb[, , 3] <- pmin(base + 0.12 * dead, 1)
  as.raster(rgb)
}

manifest <- read_tbl(paths$annotation_manifest)
points <- read_tbl(paths$annotation_points)
production_calls <- read_tbl(paths$production_calls)
production_frames <- read_tbl(paths$production_frames)
object_predictions <- read_tbl(paths$object_predictions)
frame_count_metrics <- read_tbl(paths$frame_count_metrics)
headline <- read_tbl(paths$headline)
object_metrics <- read_tbl(paths$object_metrics)
output_manifest <- read_tbl(paths$output_manifest)
provenance <- jsonlite::fromJSON(paths$provenance, simplifyVector = TRUE)
prior_frame_counts <- read_tbl(paths$prior_frame_counts)

expected_sha <- "e20436c7a5ef141f06996b47982bcfaf5f17bb1a799f601d297008da6d4b3aed"
if (!identical(as.character(headline$classifier_sha256[[1]]), expected_sha)) {
  stop("Production headline SHA does not match expected classifier SHA.", call. = FALSE)
}
if (!identical(as.character(provenance$classifier$observed_sha256), expected_sha)) {
  stop("Production provenance SHA does not match expected classifier SHA.", call. = FALSE)
}

make_overlay_tile <- function(frame, mode) {
  m <- manifest %>% filter(.data$stack_index_1based == frame) %>% slice(1)
  calls <- production_calls %>% filter(.data$frame == !!frame)
  if (nrow(m) != 1 || nrow(calls) == 0) stop("Missing overlay input for frame ", frame, call. = FALSE)
  crop_h <- as.numeric(m$crop_y1 - m$crop_y0)
  crop_w <- as.numeric(m$crop_x1 - m$crop_x0)
  target <- tibble(
    xmin = as.numeric(m$target_x0 - m$crop_x0),
    xmax = as.numeric(m$target_x1 - m$crop_x0),
    ymin = crop_h - as.numeric(m$target_y1 - m$crop_y0),
    ymax = crop_h - as.numeric(m$target_y0 - m$crop_y0)
  )
  rects <- calls %>%
    mutate(
      xmin = .data$bbox_x0,
      xmax = .data$bbox_x1,
      ymin = crop_h - .data$bbox_y1,
      ymax = crop_h - .data$bbox_y0,
      manual_state = if_else(.data$label == "alive", "manual alive", "manual not alive"),
      production_state = case_when(
        .data$predicted_label_name == "alive" ~ "production alive",
        .data$predicted_label_name == "dead" ~ "production dead",
        .data$predicted_label_name == "junk" ~ "production junk",
        TRUE ~ "production not alive"
      )
    )
  point_frame <- points %>%
    filter(.data$Frame == frame) %>%
    transmute(x = .data$X, y = crop_h - .data$Y)
  label <- paste0(m$cellLine, " ", m$ploidy, ", G0 ", m$G0, " mM")
  mode_label <- if (mode == "manual") "manual labels" else "production calls"
  base <- ggplot() +
    annotation_raster(annotation_frame_rgb(frame), xmin = 0, xmax = crop_w, ymin = 0, ymax = crop_h) +
    geom_rect(aes(xmin = 0, xmax = crop_w, ymin = 0, ymax = crop_h), fill = NA, color = "#3b424d", linewidth = 0.22) +
    geom_rect(data = target, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = "#f2c94c", linewidth = 0.48) +
    annotate("text", x = 9, y = crop_h - 10, label = label, hjust = 0, vjust = 1, size = 2.15, fontface = "bold", color = "white") +
    annotate("label", x = 9, y = 10, label = mode_label, hjust = 0, vjust = 0, size = 1.8, fontface = "bold", color = "white", fill = "#20242b", alpha = 0.72, label.size = 0) +
    coord_fixed(xlim = c(0, crop_w), ylim = c(0, crop_h), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1), panel.background = element_rect(fill = "#20242b", color = NA))
  if (mode == "manual") {
    return(
      base +
        geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = manual_state), fill = NA, linewidth = 0.19, alpha = 0.82) +
        geom_point(data = point_frame, aes(x, y), color = "#ffd33d", fill = "#ffd33d", size = 0.55, alpha = 0.88) +
        scale_color_manual(values = c("manual alive" = "#20a568", "manual not alive" = "#d9dee7"), guide = "none")
    )
  }
  base +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = production_state), fill = NA, linewidth = 0.19, alpha = 0.86) +
    scale_color_manual(values = c("production alive" = "#0f9e9a", "production dead" = "#c73e5a", "production junk" = "#7b61ff", "production not alive" = "#c73e5a"), guide = "none")
}

frame_ids <- sort(unique(production_calls$frame))
overlay_rows <- lapply(frame_ids, function(frame) {
  make_overlay_tile(frame, "manual") | make_overlay_tile(frame, "production")
})
overlay_plot <- wrap_plots(overlay_rows, ncol = 1)
save_png(overlay_plot, "fg1_measurement_figure1a_production_overlay.png", width = 6.8, height = 3.15)

auc_rank <- function(y_true, y_prob) {
  ok <- is.finite(y_true) & is.finite(y_prob)
  y_true <- as.integer(y_true[ok])
  y_prob <- as.numeric(y_prob[ok])
  n_pos <- sum(y_true == 1L)
  n_neg <- sum(y_true == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(y_prob, ties.method = "average")
  (sum(ranks[y_true == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

target_bounds <- manifest %>%
  transmute(
    frame = as.integer(.data$stack_index_1based),
    target_x0_local = as.numeric(.data$target_x0 - .data$crop_x0),
    target_x1_local = as.numeric(.data$target_x1 - .data$crop_x0),
    target_y0_local = as.numeric(.data$target_y0 - .data$crop_y0),
    target_y1_local = as.numeric(.data$target_y1 - .data$crop_y0)
  )
target_scoped_predictions <- object_predictions %>%
  left_join(target_bounds, by = "frame") %>%
  mutate(
    in_annotation_target = .data$centroid_x >= .data$target_x0_local &
      .data$centroid_x <= .data$target_x1_local &
      .data$centroid_y >= .data$target_y0_local &
      .data$centroid_y <= .data$target_y1_local,
    keep_for_s2a_matched_metrics = .data$label_int == 1L | .data$in_annotation_target,
    y_true = as.integer(.data$label_int == 1L),
    y_pred = as.integer(.data$production_y_pred_alive == 1L),
    y_prob = as.numeric(.data$production_y_prob_alive)
  ) %>%
  filter(.data$keep_for_s2a_matched_metrics)

tp <- sum(target_scoped_predictions$y_true == 1L & target_scoped_predictions$y_pred == 1L, na.rm = TRUE)
fp <- sum(target_scoped_predictions$y_true == 0L & target_scoped_predictions$y_pred == 1L, na.rm = TRUE)
tn <- sum(target_scoped_predictions$y_true == 0L & target_scoped_predictions$y_pred == 0L, na.rm = TRUE)
fn <- sum(target_scoped_predictions$y_true == 1L & target_scoped_predictions$y_pred == 0L, na.rm = TRUE)
target_metrics <- tibble(
  n_objects = nrow(target_scoped_predictions),
  precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
  recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
  specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
  balanced_accuracy = (recall + specificity) / 2,
  roc_auc = auc_rank(target_scoped_predictions$y_true, target_scoped_predictions$y_prob),
  tp = tp,
  fp = fp,
  tn = tn,
  fn = fn
)

metric_plot_data <- tibble(
  metric = c("ROC AUC", "Balanced accuracy", "Precision", "Recall"),
  value = c(
    target_metrics$roc_auc[[1]],
    target_metrics$balanced_accuracy[[1]],
    target_metrics$precision[[1]],
    target_metrics$recall[[1]]
  )
) %>%
  mutate(metric = factor(metric, levels = rev(metric)))

metric_plot <- ggplot(metric_plot_data, aes(.data$value, .data$metric, fill = .data$metric)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = number(.data$value, accuracy = 0.01)), hjust = -0.12, size = 2.8) +
  coord_cartesian(xlim = c(0, 1.03), clip = "off") +
  scale_fill_manual(values = c("ROC AUC" = "#5A6FAE", "Balanced accuracy" = "#26938E", "Precision" = "#D65F5F", "Recall" = "#5BA867")) +
  labs(
    x = "target-scoped classifier value",
    y = NULL
  ) +
  theme_fg1(base_size = 8) +
  theme(axis.text.y = element_text(size = 7.2), plot.subtitle = element_text(size = 7.2), plot.margin = margin(4, 12, 4, 4))
save_png(metric_plot, "fg1_measurement_figure1c_production_metrics.png", width = 3.8, height = 1.55)

count_plot_data <- bind_rows(
  production_frames %>%
    transmute(
      manual_alive_objects = .data$manual_alive_objects,
      classifier_alive_objects = .data$predicted_alive_prob_sum,
      count_type = "probability-summed alive count"
    ),
  production_frames %>%
    transmute(
      manual_alive_objects = .data$manual_alive_objects,
      classifier_alive_objects = .data$predicted_alive_count,
      count_type = "thresholded alive calls"
    )
) %>%
  mutate(
    count_type = factor(
      .data$count_type,
      levels = c("probability-summed alive count", "thresholded alive calls")
    )
  )

max_count <- max(
  count_plot_data$manual_alive_objects,
  count_plot_data$classifier_alive_objects,
  na.rm = TRUE
)
prob_r <- stats::cor(
  production_frames$manual_alive_objects,
  production_frames$predicted_alive_prob_sum,
  use = "complete.obs",
  method = "pearson"
)
count_label <- paste0(
  "prob-sum r=", number(prob_r, accuracy = 0.01),
  "\n",
  "threshold overpredicts\n",
  comma(headline$predicted_alive_threshold_total[[1]]),
  " vs manual ", comma(headline$manual_alive_objects_total[[1]])
)
count_plot <- ggplot(
  count_plot_data,
  aes(.data$manual_alive_objects, .data$classifier_alive_objects, color = .data$count_type)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_point(size = 1.15, alpha = 0.74) +
  annotate(
    "text",
    x = max_count * 0.98,
    y = max_count * 0.07,
    label = count_label,
    hjust = 1,
    vjust = 0,
    size = 2.1,
    lineheight = 0.93,
    color = "#202124"
  ) +
  coord_equal(xlim = c(0, max_count * 1.02), ylim = c(0, max_count * 1.02), expand = FALSE, clip = "off") +
  scale_color_manual(
    values = c("probability-summed alive count" = "#1B9E77", "thresholded alive calls" = "#D95F02"),
    name = NULL
  ) +
  labs(
    x = "manual alive objects",
    y = "classifier alive count"
  ) +
  theme_fg1(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.16, "in"),
    legend.text = element_text(size = 6.0),
    axis.text = element_text(size = 6.2),
    axis.title = element_text(size = 6.9),
    plot.margin = margin(4, 4, 4, 4)
  )
save_png(count_plot, "fg1_measurement_figure1b_count_agreement.png", width = 3.65, height = 2.0)

line_palette <- c(
  "MCF10A" = "#1B9E77",
  "MDA-MB-231" = "#D95F02",
  "SNU668" = "#7570B3",
  "SUM-159-chem" = "#E7298A",
  "SUM-159-fuse" = "#66A61E"
)

s2a <- ggplot(production_frames, aes(.data$manual_annotation_points, .data$manual_alive_objects, color = .data$cellLine, shape = .data$glucose_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_point(size = 1.55, alpha = 0.9) +
  coord_equal(expand = TRUE) +
  scale_color_manual(values = line_palette, name = "cell line") +
  labs(
    x = "manual alive annotation points",
    y = "alive object labels",
    shape = "G0 bin"
  ) +
  theme_fg1(base_size = 8) +
  theme(legend.position = "right", legend.key.size = unit(0.13, "in"), legend.text = element_text(size = 5.5), legend.title = element_text(size = 6.2))
save_png(s2a, "fg1_measurement_s2a_alive_object_labels.png", width = 3.6, height = 2.45)

s2b <- ggplot(production_frames, aes(.data$manual_alive_objects, .data$predicted_alive_count, color = .data$legacy_split)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_point(size = 1.55, alpha = 0.9) +
  coord_cartesian(expand = TRUE) +
  scale_color_manual(values = c(train = "grey55", test = "#D55E00"), name = "legacy frame split") +
  labs(
    x = "manual alive objects",
    y = "production-predicted alive objects"
  ) +
  theme_fg1(base_size = 8) +
  theme(legend.position = "right", legend.key.size = unit(0.13, "in"), legend.text = element_text(size = 6.0), legend.title = element_text(size = 6.2))
save_png(s2b, "fg1_measurement_s2b_production_predictions.png", width = 3.6, height = 2.45)

source_counts <- tibble(
  source = c("annotation frames", "manual alive points", "crop objects scored", "production overlay objects", "full images in count lineage"),
  n = c(headline$n_frames, headline$manual_annotation_points_total, headline$n_objects_scored, nrow(production_calls), 33600),
  status = c("available", "available", "production scored", "production scored", "downstream production lineage")
) %>%
  mutate(source = factor(source, levels = rev(source)))

source_plot <- ggplot(source_counts, aes(log10(.data$n), .data$source, fill = .data$status)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = comma(.data$n)), hjust = -0.06, size = 2.6) +
  coord_cartesian(xlim = c(0, 5.15), clip = "off") +
  scale_fill_manual(values = c("available" = "#5A6FAE", "production scored" = "#26938E", "downstream production lineage" = "#B9812C"), name = NULL) +
  labs(x = "log10 count", y = NULL) +
  theme_fg1(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.7), legend.position = "bottom", plot.margin = margin(4, 14, 4, 4))
save_png(source_plot, "fg1_measurement_provenance_source_counts.png", width = 4.1, height = 1.9)

lineage_text <- paste(
  "Production classifier: data/train/classifier_training_outputs/object_classifier.pkl",
  paste0("sha256: ", expected_sha),
  "Class mapping: 1 alive; 2 dead; 3 junk. Binary comparison maps 2 and 3 to not alive.",
  paste0("Scored export: ", headline$n_frames, " frames, ", comma(headline$n_objects_scored), " objects."),
  "Downstream analyses already trace through batch_results_v2.csv and stan_ready_data.Rds; no downstream rerun is required for this display correction.",
  sep = "\n"
)
lineage_plot <- ggplot() +
  annotate("label", x = 0, y = 1, label = lineage_text, hjust = 0, vjust = 1, size = 2.65, label.size = 0.25, fill = "#f7f8fa", color = "#202124", label.r = unit(0.03, "lines")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))
save_png(lineage_plot, "fg1_measurement_production_lineage_tile.png", width = 6.5, height = 1.45)

notes <- c(
  "# Measurement/provenance production-classifier worker notes",
  "",
  "Generated production-classifier replacements from `data/report_exports/v5_fg1_production_classifier_validation/`.",
  "",
  "Outputs:",
  "",
  "- `fg1_measurement_figure1a_production_overlay.png`: production classifier overlay for the two selected prior overlay frames (`MCF10A` and `SUM-159-fuse`).",
  "- `fg1_measurement_figure1b_count_agreement.png`: recommended replacement for the confusing Figure 1B metric bars; it shows frame-level manual alive objects against probability-summed and thresholded production alive counts.",
  "- `fg1_measurement_figure1c_production_metrics.png`: round-3 regenerated production metrics matched to the current Figure S2A target scope. The metric set keeps all S2A alive object labels and restricts comparison negatives to the annotated 200x200 target region, so it no longer uses the stale all-context 400x400 object metric scope.",
  "- `fg1_measurement_s2a_alive_object_labels.png`: relabeled manual alive annotation points versus alive object labels; avoids the unexplained word `kept`.",
  "- `fg1_measurement_s2b_production_predictions.png`: production-predicted alive objects versus manual alive objects, colored by legacy train/test frame split.",
  "- `fg1_measurement_provenance_source_counts.png` and `fg1_measurement_production_lineage_tile.png`: provenance/context panels for review.",
  "",
  "Production artifact validation:",
  "",
  paste0("- Observed classifier sha256 matches expected `", expected_sha, "`."),
  "- Production class `1` is mapped to alive; classes `2` and `3` are mapped to not alive for binary manual-label comparison.",
  paste0("- Performance caveat on the 90-frame annotation set: precision ", number(headline$object_precision[[1]], accuracy = 0.001),
         ", recall ", number(headline$object_recall[[1]], accuracy = 0.001),
         ", target-scoped precision ", number(target_metrics$precision[[1]], accuracy = 0.001),
         ", recall ", number(target_metrics$recall[[1]], accuracy = 0.001),
         ", balanced accuracy ", number(target_metrics$balanced_accuracy[[1]], accuracy = 0.001),
         ", and ROC AUC ", number(target_metrics$roc_auc[[1]], accuracy = 0.001),
         " across ", comma(target_metrics$n_objects[[1]]), " target-scoped objects."),
  "- Round-3 feedback requested production metrics that match the current recommended Figure S2A; the regenerated metric panel therefore recomputes classifier metrics from `fg1_production_classifier_object_predictions.csv` after keeping all alive labels used by S2A and limiting not-alive comparison objects to the 200x200 annotated target region.",
  "- The production replacement supports artifact provenance and high recall/frame-level rank agreement, but it should not be read as strong thresholded alive-count accuracy.",
  "- The first-pass classifier outputs are not reused as fixed replacement panels. They are retained only as prior-code/provenance context.",
  "- No Stan/model-free/posterior/selection rerun is required for this display correction because the v4 production-classifier issue note validated the downstream count lineage.",
  "",
  "Directive coverage:",
  "",
  "- `FG1-D01`: addressed with a production overlay; the optional multi-line suggestion is partially addressed by using the two selected overlay frames available in the production selected-call export.",
  "- `FG1-D02`: addressed with production classifier metrics from `fg1_production_classifier_headline_metrics.csv`.",
  "- `FG1-R3-D01`: addressed by regenerating `fg1_measurement_figure1c_production_metrics.png` from target-scoped `fg1_production_classifier_object_predictions.csv` so the panel matches the current recommended Figure S2A annotation scope.",
  "- `FG1-D03`: addressed for S2A labels and S2B production predictions/split notation.",
  "- `FG1-R2-D01`: addressed with a recommended count-agreement display from `fg1_production_classifier_frame_counts.csv`; the metric bars are retained only as comparison/context.",
  "",
  "Visual QC:",
  "",
  "- PNGs have no figure-level titles or captions.",
  "- Axis labels use explicit quantities for S2A/S2B.",
  "- The production overlay is dense but readable at review size; selected-object rectangles are intentionally thin to preserve image context."
)
writeLines(notes, file.path(worker_dir, "measurement_provenance_coverage.md"))

manifest_out <- tibble(
  output_png = file.path(initial_dir, c(
    "fg1_measurement_figure1a_production_overlay.png",
    "fg1_measurement_figure1b_count_agreement.png",
    "fg1_measurement_figure1c_production_metrics.png",
    "fg1_measurement_s2a_alive_object_labels.png",
    "fg1_measurement_s2b_production_predictions.png",
    "fg1_measurement_provenance_source_counts.png",
    "fg1_measurement_production_lineage_tile.png"
  )),
  refined_copy = file.path(refined_dir, basename(output_png)),
  directive_ids = c("FG1-D01", "FG1-R2-D01", "FG1-D02;FG1-R3-D01", "FG1-D03", "FG1-D03", "FG1-D01;FG1-D02;FG1-D03", "FG1-D01;FG1-D02;FG1-D03"),
  data_source = c(
    paths$production_calls,
    paths$production_frames,
    paths$object_predictions,
    paths$production_frames,
    paths$production_frames,
    paths$headline,
    paths$provenance
  )
)
write.csv(manifest_out, file.path(initial_dir, "measurement_provenance_manifest.csv"), row.names = FALSE)

cat("Wrote FG1 drafting_v2 measurement/provenance production-classifier panels under ", root, "\n", sep = "")
