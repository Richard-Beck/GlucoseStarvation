# These constructors preserve the accepted v02 plotting code. Only the data
# frames supplied to them are rebuilt from the 20260722 canonical objects.

old_ploidy_colors <- c(low = "#2166AC", high = "#B2182B")
fg2_ploidy_palette <- c(low = "#1B6CA8", high = "#C43C39")
g0_levels <- c("0", "0.1", "0.25", "0.5", "1", "5", "25")
format_g0 <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  out <- ifelse(abs(y) < 1e-8, "0", ifelse(abs(y - round(y)) < 1e-8, sprintf("%.0f", y), sprintf("%.2f", y)))
  out <- ifelse(grepl("\\.", out), sub("0$", "", out), out)
  sub("\\.$", "", out)
}
short_number <- function() label_number(scale_cut = cut_short_scale())
label_short <- function(x) {
  vapply(x, function(one) {
    if (!is.finite(one)) return("")
    if (abs(one) >= 1e6) return(sprintf("%.1fM", one / 1e6))
    if (abs(one) >= 1e3) return(sprintf("%.1fk", one / 1e3))
    if (abs(one) >= 1) return(number(one, accuracy = .1))
    number(one, accuracy = .01)
  }, character(1))
}

theme_fg1_recycled <- function(base_size = 7) {
  theme_bw(base_size = base_size) + theme(
    text = element_text(family = "sans"), panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = .14, color = "grey88"),
    strip.background = element_rect(fill = "grey94", color = "grey75", linewidth = .22),
    strip.text = element_text(face = "bold", margin = margin(1.2, 1.2, 1.2, 1.2)),
    legend.position = "bottom", legend.title = element_text(size = rel(.9)),
    legend.key.width = unit(.22, "in"), plot.margin = margin(3, 3, 3, 3),
    plot.subtitle = element_blank(), plot.caption = element_blank()
  )
}
theme_fg2_recycled <- function(base_size = 7) {
  theme_bw(base_size = base_size) + theme(
    text = element_text(family = "sans"), panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = .14, color = "grey88"),
    strip.background = element_rect(fill = "grey94", color = "grey74", linewidth = .22),
    strip.text = element_text(face = "bold", margin = margin(1.3, 1.3, 1.3, 1.3)),
    legend.position = "bottom", legend.box = "horizontal", legend.key.height = unit(.12, "in"),
    legend.key.width = unit(.22, "in"), axis.title = element_text(size = rel(.92)),
    axis.text = element_text(size = rel(.82)), plot.margin = margin(3, 5, 3, 3),
    plot.subtitle = element_blank(), plot.caption = element_blank()
  )
}

make_overlay_tile_recycled <- function(frame) {
  meta <- annotation_manifest %>% filter(stack_index_1based == frame) %>% slice(1)
  vm <- validation_manifest %>% filter(image_key == meta$image_key) %>% slice(1)
  channels <- list(
    alive = read_json_channel(vm$alive_paths_json),
    dead = read_json_channel(vm$dead_paths_json),
    phase = read_json_channel(vm$phase_paths_json)
  )
  yy <- (meta$crop_y0 + 1L):meta$crop_y1; xx <- (meta$crop_x0 + 1L):meta$crop_x1
  green <- robust01(channels$alive[yy, xx]); red <- robust01(channels$dead[yy, xx]); phase_img <- robust01(channels$phase[yy, xx])
  rgb <- array(0, c(length(yy), length(xx), 3L)); rgb[, , 1] <- pmin(1, .58 * phase_img + .72 * red); rgb[, , 2] <- pmin(1, .58 * phase_img + .72 * green); rgb[, , 3] <- pmin(1, .58 * phase_img)
  labels <- round(readTIFF(mask_path(meta$image_key), native = FALSE)[yy, xx] * 65535)
  calls <- validation_predictions %>% filter(frame == !!frame)
  outlines <- boundary_points(labels, calls) %>% rename(production_state = predicted_label_name)
  pts <- validation_points %>% filter(frame == !!frame)
  h <- nrow(labels); w <- ncol(labels)
  target <- tibble(xmin = meta$target_x0 - meta$crop_x0, xmax = meta$target_x1 - meta$crop_x0,
                   ymin = h - (meta$target_y1 - meta$crop_y0), ymax = h - (meta$target_y0 - meta$crop_y0))
  label_text <- paste0(meta$cellLine, " ", meta$ploidy, "\n", format_g0(meta$G0), " mM, ", format_g0(meta$hours), " h")
  ggplot() + annotation_raster(as.raster(rgb), xmin = 0, xmax = w, ymin = 0, ymax = h) +
    geom_tile(data = outlines, aes(x, y, fill = production_state), width = 1.15, height = 1.15, alpha = .92, show.legend = FALSE) +
    geom_point(data = pts, aes(x_crop, h - y_crop), color = "#1a1a1a", fill = "#FFD92F", shape = 21, stroke = .32, size = 1.08) +
    geom_rect(data = target, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = "#f2c94c", linewidth = .34) +
    annotate("rect", xmin = 0, xmax = w, ymin = 0, ymax = h, fill = NA, color = "#3b424d", linewidth = .16) +
    annotate("label", x = 5, y = h - 5, label = label_text, hjust = 0, vjust = 1, size = 1.55, lineheight = .86, fontface = "bold", color = "white", fill = alpha("#20242b", .70), label.size = 0, label.r = unit(.025, "lines")) +
    scale_fill_manual(values = c(alive = "#0f9e9a", dead = "#c73e5a", junk = "#7b61ff"), guide = "none") +
    coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE) + theme_void() +
    theme(plot.margin = margin(.4, .4, .4, .4), panel.background = element_rect(fill = "#20242b", color = NA))
}
overlay_legend_recycled <- function() {
  ggplot() +
    annotate("segment", x = .025, xend = .080, y = .52, yend = .52, color = "#0f9e9a", linewidth = 1.2) + annotate("text", x = .090, y = .52, label = "alive", hjust = 0, size = 2.2) +
    annotate("segment", x = .205, xend = .260, y = .52, yend = .52, color = "#c73e5a", linewidth = 1.2) + annotate("text", x = .270, y = .52, label = "dead", hjust = 0, size = 2.2) +
    annotate("segment", x = .385, xend = .440, y = .52, yend = .52, color = "#7b61ff", linewidth = 1.2) + annotate("text", x = .450, y = .52, label = "junk", hjust = 0, size = 2.2) +
    annotate("point", x = .565, y = .52, shape = 21, size = 2.5, stroke = .35, color = "#1a1a1a", fill = "#FFD92F") + annotate("text", x = .590, y = .52, label = "manual alive", hjust = 0, size = 2.2) +
    annotate("rect", xmin = .795, xmax = .850, ymin = .34, ymax = .70, fill = NA, color = "#f2c94c", linewidth = .55) + annotate("text", x = .860, y = .52, label = "target", hjust = 0, size = 2.2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void()
}
make_f1a <- function() {
  frames <- c(4L, 22L, 49L, 76L, 13L, 31L, 40L, 85L)
  wrap_plots(lapply(frames, make_overlay_tile_recycled), ncol = 4) / overlay_legend_recycled() + plot_layout(heights = c(1, .055))
}

make_f1b <- function() {
  x <- area %>% mutate(G0_display = factor(format_g0(G0), levels = rev(c("25", "5", "1", "0.5", "0.25", "0.1", "0"))),
                       ploidy_state = factor(if_else(ploidy_state == "low", "baseline", "elevated"), c("baseline", "elevated"))) %>%
    group_by(cellLine, G0_display, ploidy_state) %>% summarise(median_q50 = finite_median(segmented_area_px_q50), q25_q50 = finite_q(segmented_area_px_q50, .25), q75_q50 = finite_q(segmented_area_px_q50, .75), .groups = "drop") %>%
    mutate(cellLine = factor(cellLine, line_order))
  ggplot(x, aes(G0_display, median_q50, color = ploidy_state, group = ploidy_state)) +
    geom_linerange(aes(ymin = q25_q50, ymax = q75_q50), linewidth = .28, alpha = .82) + geom_line(linewidth = .30) + geom_point(size = 1.05) +
    facet_wrap(~cellLine, ncol = 1, scales = "free_y", labeller = as_labeller(line_short)) +
    scale_color_manual(values = c(baseline = old_ploidy_colors[["low"]], elevated = old_ploidy_colors[["high"]]), name = "ploidy") +
    labs(x = "starting glucose (mM)", y = "cell area q50 (px)") + theme_fg1_recycled(6.4) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 4.8), strip.text = element_text(size = 5.5), axis.text.y = element_text(size = 4.8), axis.title = element_text(size = 6.2), legend.text = element_text(size = 5.8), panel.spacing = unit(.30, "lines"))
}

make_f1c <- function() {
  summaries <- area_conditions %>% group_by(cellLine) %>% summarise(cell_area = finite_median(log2_cell_ratio), nuclear_area = finite_median(log2_nuclear_ratio), .groups = "drop") %>%
    mutate(label = recode(cellLine, `MDA-MB-231` = "MDA-231", `SUM-159-chem` = "SUM-chem", `SUM-159-fuse` = "SUM-fuse", .default = cellLine),
           label_dx = if_else(cellLine == "SUM-159-fuse", -.06, .05), label_hjust = if_else(cellLine == "SUM-159-fuse", 1, 0))
  ggplot(area_conditions, aes(log2_cell_ratio, log2_nuclear_ratio)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = .20, color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = .20, color = "grey55") + geom_point(color = "grey74", size = .65, alpha = .9) +
    geom_point(data = summaries, aes(cell_area, nuclear_area, color = cellLine), size = 2, stroke = .75) +
    geom_text(data = summaries, aes(cell_area + label_dx, nuclear_area, label = label, color = cellLine, hjust = label_hjust), size = 2, show.legend = FALSE) +
    scale_color_manual(values = line_palette, guide = "none") + coord_cartesian(xlim = c(-1.7, 2.45), ylim = c(-1.2, 1.7), clip = "off") +
    labs(x = "cell-area ratio (log2)", y = "nuclear-area ratio (log2)") + theme_fg1_recycled(7) +
    theme(legend.position = "none", axis.text = element_text(size = 6), axis.title = element_text(size = 6.6), plot.margin = margin(3, 30, 3, 4))
}

make_f1d <- function() {
  x <- area_conditions %>% pivot_longer(c(log2_cell_ratio, log2_nuclear_ratio), names_to = "measure", values_to = "ratio") %>%
    group_by(cellLine, measure) %>% summarise(estimate = finite_median(ratio), lo = finite_q(ratio, .25), hi = finite_q(ratio, .75), .groups = "drop") %>%
    left_join(wells %>% group_by(cellLine) %>% summarise(ploidy_ratio = log2(max(ploidy_abs) / min(ploidy_abs)), .groups = "drop"), by = "cellLine") %>%
    mutate(measure = recode(measure, log2_cell_ratio = "Cell area q50", log2_nuclear_ratio = "Nuclear area"), line_group = if_else(cellLine == "SUM-159-fuse", "SUM-159-fuse", "other lines"))
  ggplot(x, aes(ploidy_ratio, estimate)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = .24, color = "grey50") + geom_abline(slope = 1, linetype = "dotted", linewidth = .22, color = "grey68") +
    geom_errorbar(aes(ymin = lo, ymax = hi, color = line_group), width = .018, linewidth = .35, alpha = .85) + geom_point(aes(color = line_group, shape = line_group), size = 2, stroke = .72, fill = "white") +
    geom_text(data = x %>% filter(cellLine == "SUM-159-fuse"), aes(label = "SUM-fuse", color = line_group), nudge_x = .045, hjust = 0, size = 2.1, show.legend = FALSE) +
    facet_wrap(~measure, nrow = 1) + scale_color_manual(values = c("other lines" = "grey20", "SUM-159-fuse" = "#B2182B"), name = "line") + scale_shape_manual(values = c("other lines" = 16, "SUM-159-fuse" = 21), name = "line") +
    scale_x_continuous(breaks = c(.33, 1), labels = c("0.33", "1.00"), limits = c(.24, 1.38)) + coord_cartesian(ylim = c(-1.18, 1.58), clip = "off") +
    labs(x = "ploidy ratio (log2)", y = "area ratio (log2)") + theme_fg1_recycled(7) + theme(strip.text = element_text(size = 6), axis.text = element_text(size = 5.8), axis.title = element_text(size = 6.3), plot.margin = margin(3, 24, 3, 3))
}

raw_live_recycled <- count_raw %>% group_by(well_idx, hours) %>% mutate(n_reps = n(), jitter_rank = row_number(), hours_plot = hours + if_else(n_reps > 1, ((jitter_rank - 1) / (n_reps - 1) - .5) * .42, 0)) %>% ungroup()
raw_glucose_recycled <- glucose_raw %>% group_by(well_idx, hours) %>% mutate(n_reps = n(), jitter_rank = row_number(), hours_plot = hours + if_else(n_reps > 1, ((jitter_rank - 1) / (n_reps - 1) - .5) * .18, 0)) %>% ungroup()
add_traj_factors <- function(x) x %>% mutate(G0_display = factor(format_g0(G0), levels = c("25", "5", "1", "0.5", "0.25", "0.1", "0")), G0_label = factor(paste0("G0 ", G0_display, " mM"), levels = paste0("G0 ", c("25", "5", "1", "0.5", "0.25", "0.1", "0"), " mM")), ploidy_state = factor(ploidy_state, c("low", "high")))
if (anyNA(add_traj_factors(raw_live_recycled)$G0_label) || anyNA(add_traj_factors(raw_glucose_recycled)$G0_label)) {
  stop("G0 adapter invariant failed: a trajectory condition mapped to an NA display label")
}

make_trajectory_block_recycled <- function(line, readout) {
  if (readout == "live") {
    pts <- raw_live_recycled %>% filter(cellLine == line, G0 %in% c(25, 1, 0)) %>% add_traj_factors()
    lines <- count_summary %>% filter(cellLine == line, G0 %in% c(25, 1, 0)) %>% add_traj_factors()
    p <- ggplot() + geom_point(data = pts, aes(hours_plot, live, color = ploidy_state), size = .42, stroke = 0) +
      geom_line(data = lines, aes(hours, live, color = ploidy_state, group = interaction(well_idx, ploidy_state)), linewidth = .26, alpha = .9) + labs(y = "live cells")
  } else {
    pts <- raw_glucose_recycled %>% filter(cellLine == line, G0 %in% c(25, 1, 0)) %>% add_traj_factors()
    lines <- glucose_summary %>% filter(cellLine == line, G0 %in% c(25, 1, 0)) %>% add_traj_factors()
    p <- ggplot() + geom_point(data = pts, aes(hours_plot, glucose, color = ploidy_state, shape = factor(censored)), size = .44, stroke = .16) +
      geom_line(data = lines, aes(hours, glucose, color = ploidy_state, group = interaction(well_idx, ploidy_state)), linewidth = .26, alpha = .92) + scale_shape_manual(values = c("0" = 16, "1" = 4), name = "censored") + labs(y = "glucose (mM)")
  }
  p + facet_wrap(~G0_label, ncol = 1, scales = "free_y") + scale_color_manual(values = old_ploidy_colors, name = "ploidy") +
    scale_x_continuous(expand = expansion(mult = c(.01, .03))) + labs(x = "hours", title = line) + theme_fg1_recycled(6) +
    theme(plot.title = element_text(face = "bold", size = 7.2, hjust = .5), strip.text = element_text(size = 4.4), axis.text.y = element_text(size = 3.6), axis.text.x = element_text(size = 4), axis.title = element_text(size = 5), panel.spacing = unit(.28, "lines"), legend.position = "none")
}
make_trajectory_grid <- function(readout = c("live", "glucose"), lines = line_order[-5], glucose_keep = c(0, 1, 25), ncol = 2) {
  readout <- match.arg(readout)
  wrap_plots(lapply(lines, function(line) make_trajectory_block_recycled(line, readout)), ncol = 2, guides = "collect") & theme(legend.position = "bottom")
}

make_s1_cell <- function(line, readout) {
  g0_rev <- paste0("G0 ", c("0", "0.1", "0.25", "0.5", "1", "5", "25"), " mM")
  if (readout %in% c("Live cells", "Dead cells")) {
    value_name <- if (readout == "Live cells") "live" else "dead"
    pts <- raw_live_recycled %>% transmute(cellLine, G0, ploidy_state, hours_plot, value = .data[[value_name]])
    lines <- count_summary %>% transmute(cellLine, G0, ploidy_state, well_idx, hours, value = .data[[value_name]])
  } else {
    pts <- raw_glucose_recycled %>% transmute(cellLine, G0, ploidy_state, hours_plot, value = glucose)
    lines <- glucose_summary %>% transmute(cellLine, G0, ploidy_state, well_idx, hours, value = glucose)
  }
  pts <- pts %>% filter(cellLine == line) %>% mutate(G0_label = factor(paste0("G0 ", format_g0(G0), " mM"), g0_rev))
  lines <- lines %>% filter(cellLine == line) %>% mutate(G0_label = factor(paste0("G0 ", format_g0(G0), " mM"), g0_rev))
  ggplot() + geom_line(data = lines, aes(hours, value, color = ploidy_state, group = interaction(well_idx, ploidy_state)), linewidth = .15, alpha = .8) +
    geom_point(data = pts, aes(hours_plot, value), color = "black", size = .2, stroke = 0) + facet_grid(G0_label ~ ., scales = "free_y") +
    scale_color_manual(values = old_ploidy_colors, guide = "none") + labs(x = "hours", y = NULL) + theme_fg1_recycled(4.45) +
    theme(legend.position = "none", axis.title.x = element_text(size = 3.55), axis.title.y = element_blank(), axis.text = element_text(size = 2.85), strip.text.y = element_text(size = 3.05, angle = 0), plot.margin = margin(1.2, 1.2, 1.2, 1.2), panel.spacing = unit(.22, "lines"))
}
make_s1_line_block_recycled <- function(line) {
  grid <- make_s1_cell(line, "Live cells") + make_s1_cell(line, "Dead cells") + make_s1_cell(line, "Glucose") + plot_layout(ncol = 3)
  ggdraw() + draw_plot(grid, x = 0, y = 0, width = 1, height = .925) + draw_label(line, x = .01, y = .995, hjust = 0, vjust = 1, size = 6.4, fontface = "bold")
}
make_s1 <- function() {
  blocks <- lapply(line_order, make_s1_line_block_recycled)
  atlas <- wrap_plots(blocks[[1]], blocks[[2]], blocks[[3]], blocks[[4]], blocks[[5]], plot_spacer(), ncol = 2)
  legend <- ggplot(tibble(x = c(.42, .53), xend = c(.47, .58), state = factor(c("low", "high"), c("low", "high")), label = c("low", "high"))) +
    geom_segment(aes(x = x, xend = xend, y = .5, yend = .5, color = state), linewidth = .35) + geom_text(aes(xend + .015, .5, label = label), hjust = 0, size = 2.3) +
    annotate("text", x = .38, y = .5, label = "ploidy", hjust = 1, size = 2.3) + scale_color_manual(values = old_ploidy_colors, guide = "none") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
  plot_grid(atlas, legend, ncol = 1, rel_heights = c(1, .085))
}

make_s2a <- function() {
  ggplot(validation_frames, aes(manual_point_count, manual_alive_object_count, color = cellLine, shape = glucose_bin)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = .35) + geom_point(size = 1.55, alpha = .9) + coord_equal(expand = TRUE) +
    scale_color_manual(values = c("MCF10A" = "#1B9E77", "MDA-MB-231" = "#D95F02", "SNU668" = "#7570B3", "SUM-159-chem" = "#E7298A", "SUM-159-fuse" = "#66A61E"), name = "cell line") +
    labs(x = "manual alive annotation points", y = "alive object labels", shape = "G0 bin") + theme_fg1_recycled(8) +
    theme(legend.position = "right", legend.key.size = unit(.13, "in"), legend.text = element_text(size = 5.5), legend.title = element_text(size = 6.2))
}
make_s2b <- function() {
  x <- validation_metrics %>% transmute(`ROC AUC` = roc_auc, `Balanced accuracy` = balanced_accuracy, Precision = precision, Recall = recall) %>% pivot_longer(everything(), names_to = "metric", values_to = "value") %>% mutate(metric = factor(metric, levels = rev(c("ROC AUC", "Balanced accuracy", "Precision", "Recall"))))
  ggplot(x, aes(value, metric, fill = metric)) + geom_col(width = .62, show.legend = FALSE) + geom_text(aes(label = number(value, accuracy = .01)), hjust = -.12, size = 2.8) +
    coord_cartesian(xlim = c(0, 1.03), clip = "off") + scale_fill_manual(values = c("ROC AUC" = "#5A6FAE", "Balanced accuracy" = "#26938E", Precision = "#D65F5F", Recall = "#5BA867")) +
    labs(x = "target-scoped classifier value", y = NULL) + theme_fg1_recycled(8) + theme(axis.text.y = element_text(size = 7.2), plot.margin = margin(4, 12, 4, 4))
}

coverage_recycled <- {
  counts <- count_raw %>% group_by(cellLine, ploidy_state, G0) %>% summarise(n_count_times = n_distinct(hours), .groups = "drop")
  gluc <- glucose_raw %>% group_by(cellLine, ploidy_state, G0) %>% summarise(n_glucose_times = n_distinct(hours), glucose_censored_observations = sum(censored), glucose_censored_fraction = mean(censored), .groups = "drop")
  full_join(counts, gluc, by = c("cellLine", "ploidy_state", "G0")) %>% mutate(
    count_time_coverage = n_count_times / max(n_count_times, na.rm = TRUE),
    condition_label = paste(line_short[cellLine], ploidy_state, sep = " | "),
    G0_display = factor(format_g0(G0), levels = c("25", "5", "1", "0.5", "0.25", "0.1", "0"))
  )
}
make_s2c <- function() {
  coverage_recycled %>% mutate(condition_label = factor(condition_label, levels = rev(unique(condition_label))), G0_axis = factor(as.character(G0_display), levels = rev(levels(G0_display)))) %>%
    ggplot(aes(G0_axis, condition_label, fill = count_time_coverage)) + geom_tile(color = "white", linewidth = .2) + geom_text(aes(label = paste0(n_count_times, "c/", n_glucose_times, "g")), size = 1.55) +
    scale_fill_gradient(low = "white", high = "#386CB0", limits = c(0, 1), name = "timepoint\ncoverage") + labs(x = "starting glucose (mM)", y = NULL) + theme_fg1_recycled(5.8) + theme(axis.text.y = element_text(size = 4.2), legend.position = "right")
}
make_s2d <- function() {
  coverage_recycled %>% mutate(condition_label = factor(condition_label, levels = rev(unique(condition_label))), G0_axis = factor(as.character(G0_display), levels = rev(levels(G0_display)))) %>%
    ggplot(aes(G0_axis, condition_label, fill = glucose_censored_fraction)) + geom_tile(color = "white", linewidth = .2) + geom_text(aes(label = if_else(glucose_censored_observations > 0, as.character(glucose_censored_observations), "")), size = 1.65) +
    scale_fill_gradient(low = "white", high = "#D95F02", name = "censored\nfraction") + labs(x = "starting glucose (mM)", y = NULL) + theme_fg1_recycled(5.8) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")
}
calib_recycled <- calibration %>% mutate(calibration_set = batch_label, G = glucose, Lum = lum, fitted_lum = fitted, log_residual = log(pmax(Lum, 1e-12)) - log(pmax(fitted_lum, 1e-12)))
make_s2e <- function() ggplot(calib_recycled, aes(G, Lum, color = calibration_set)) + geom_point(size = .65, alpha = .85) + geom_line(aes(y = fitted_lum), color = "black", linewidth = .28) + facet_wrap(~calibration_set, scales = "free_y", ncol = 2) + scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 25)) + labs(x = "known glucose", y = "luminescence") + theme_fg1_recycled(5.9) + theme(legend.position = "none", strip.text = element_text(size = 4.8), axis.text = element_text(size = 4.4))
make_s2f <- function() {
  x <- calib_recycled %>% group_by(calibration_set, G) %>% summarise(residual_median = finite_median(log_residual), residual_q25 = finite_q(log_residual, .25), residual_q75 = finite_q(log_residual, .75), .groups = "drop")
  ggplot(x, aes(G, residual_median, color = calibration_set)) + geom_hline(yintercept = 0, linewidth = .25, color = "grey45") + geom_linerange(aes(ymin = residual_q25, ymax = residual_q75), alpha = .45, linewidth = .3) + geom_point(size = .65) +
    facet_wrap(~calibration_set, ncol = 2) + scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 25)) + labs(x = "known glucose", y = "log residual") + theme_fg1_recycled(5.9) + theme(legend.position = "none", strip.text = element_text(size = 4.8), axis.text = element_text(size = 4.4))
}
make_s2g <- function() {
  x <- calib_recycled %>% group_by(calibration_set) %>% summarise(log_rmse = sqrt(mean(log_residual^2)), .groups = "drop")
  ggplot(x, aes(reorder(calibration_set, log_rmse), log_rmse)) + geom_col(width = .62, fill = "#6B7FAF") + coord_flip() + labs(x = NULL, y = "log residual RMSE") + theme_fg1_recycled(6.7)
}

make_f2a <- function() {
  effects <- as_tibble(feature_release$ploidy_effects) %>% filter(cellLine != "SUM-159-fuse") %>% group_by(feature_id) %>%
    mutate(effect_scale = sd(effect_per_ploidy), standardized_effect = effect_per_ploidy / if_else(is.finite(effect_scale) & effect_scale > 0, effect_scale, 1), clipped_effect = pmax(-3, pmin(3, standardized_effect))) %>% ungroup() %>%
    mutate(display_label = factor(display_label, levels = rev(feature_release$feature_catalog$display_label)), feature_class = factor(feature_class, levels = names(feature_palette)))
  medians <- effects %>% group_by(feature_id, display_label, feature_class) %>% summarise(clipped_effect = pmax(-3, pmin(3, median(standardized_effect))), .groups = "drop")
  ggplot(effects, aes(clipped_effect, display_label)) + geom_vline(xintercept = 0, linewidth = .25, color = "grey55") +
    geom_point(aes(fill = feature_class), shape = 21, color = "grey18", stroke = .22, size = 1.85, alpha = .82, position = position_jitter(width = 0, height = .055, seed = 20260711)) +
    geom_point(data = medians, aes(clipped_effect, display_label, color = feature_class), inherit.aes = FALSE, shape = 23, fill = "white", stroke = .45, size = 2.25) +
    facet_grid(feature_class ~ ., scales = "free_y", space = "free_y") + scale_fill_manual(values = feature_palette, guide = "none") + scale_color_manual(values = feature_palette, guide = "none") +
    scale_x_continuous(limits = c(-3.05, 3.05), oob = squish, breaks = -3:3) + labs(x = "high-low effect per ploidy (within-feature SD units)", y = NULL) +
    theme_fg2_recycled(7.2) + theme(panel.spacing.y = unit(.18, "lines"), strip.text.y = element_blank(), strip.background.y = element_blank(), axis.text.y = element_text(size = 6.7), legend.position = "none")
}

representative_recycled <- condition_features %>% filter(cellLine == "MDA-MB-231", ploidy_state == "low", abs(G0 - 1) < 1e-8) %>% slice(1)
make_f2b <- function() {
  observed <- count_summary %>% filter(well_idx == representative_recycled$well_idx) %>% arrange(hours) %>% mutate(live_cells = live, total_cells = total, log_live = log1p(pmax(live, 0)), log_total = log1p(pmax(total, 0)))
  fit <- smooth.spline(observed$hours, observed$log_live, spar = .62)
  example_smooth <- tibble(hours = seq(min(observed$hours), max(observed$hours), length.out = 260)) %>% mutate(smooth_log_live = predict(fit, hours)$y)
  robust_value <- representative_recycled$robust_max_derivative; robust_time <- representative_recycled$robust_derivative_time
  y0 <- approx(example_smooth$hours, example_smooth$smooth_log_live, xout = robust_time, rule = 2)$y
  x_rng <- range(example_smooth$hours); tangent_x <- seq(max(x_rng[1], robust_time - 24), min(x_rng[2], robust_time + 24), length.out = 2)
  tangent <- tibble(hours = tangent_x, log_live = y0 + robust_value * (tangent_x - robust_time))
  start_time <- representative_recycled$glucose_start_time; end_time <- representative_recycled$glucose_end_time
  auc_band <- example_smooth %>% filter(hours >= start_time, hours <= end_time) %>% mutate(ymin = min(observed$log_live) - .12)
  yield_window <- observed %>% filter(hours >= start_time, hours <= end_time); yield_min <- min(yield_window$total_cells); yield_max <- max(yield_window$total_cells); wr <- range(yield_window$hours)
  arrow_x <- wr[2] - .08 * diff(wr); label_x <- min(x_rng[2] - .05 * diff(x_rng), arrow_x + .07 * diff(x_rng)); label_y <- mean(log1p(c(yield_min, yield_max)))
  point_df <- bind_rows(observed %>% transmute(hours, log_cells = log_live, series = "live cells"), observed %>% transmute(hours, log_cells = log_total, series = "total cells"))
  series_colors <- c("live cells" = "#1B6CA8", "total cells" = feature_palette[["Total-cell yield"]])
  ggplot() + geom_ribbon(data = auc_band, aes(hours, ymin = ymin, ymax = smooth_log_live), fill = feature_palette[["Alive AUC"]], alpha = .15) +
    geom_hline(yintercept = log1p(c(yield_min, yield_max)), color = feature_palette[["Total-cell yield"]], linewidth = .34, linetype = "dashed") +
    annotate("segment", x = arrow_x, xend = arrow_x, y = log1p(yield_min), yend = log1p(yield_max), arrow = arrow(length = unit(.055, "in"), ends = "both"), color = feature_palette[["Total-cell yield"]], linewidth = .45) +
    annotate("text", x = label_x, y = label_y, label = "max delta", color = feature_palette[["Total-cell yield"]], size = 1.9, hjust = 0) +
    geom_point(data = point_df, aes(hours, log_cells, color = series), size = .68, alpha = .86) + geom_line(data = example_smooth, aes(hours, smooth_log_live), color = series_colors[["live cells"]], linewidth = .74) +
    geom_line(data = tangent, aes(hours, log_live), color = feature_palette[["Growth"]], linewidth = .72, linetype = "dashed") + geom_point(data = tibble(hours = robust_time, log_live = y0), aes(hours, log_live), shape = 21, fill = "white", color = feature_palette[["Growth"]], stroke = .42, size = 1.45) +
    annotate("text", x = min(x_rng[2] - .06 * diff(x_rng), robust_time + 18), y = y0 + .34, label = "derivative tangent", color = feature_palette[["Growth"]], size = 1.9, hjust = 1) +
    annotate("text", x = mean(c(start_time, end_time)), y = min(observed$log_live) - .02, label = "live AUC", color = feature_palette[["Alive AUC"]], size = 1.9, vjust = 1) +
    scale_color_manual(values = series_colors, name = NULL) + coord_cartesian(clip = "off") + labs(x = "hours", y = "log1p(cells)") + theme_fg2_recycled(6.9) + theme(legend.position = "bottom", legend.justification = "left", legend.text = element_text(size = 6), plot.margin = margin(4, 14, 4, 4))
}

make_f2c <- function() {
  x <- condition_features %>% filter(cellLine == "MDA-MB-231") %>% mutate(G0_display = factor(format_g0(G0), levels = g0_levels), ploidy_state = factor(ploidy_state, c("low", "high")))
  dodge <- position_dodge(width = .34)
  deriv <- x
  p1 <- ggplot(deriv, aes(G0_display, robust_max_derivative, color = ploidy_state, group = ploidy_state)) + geom_hline(yintercept = 0, linewidth = .18, color = "grey72") + geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed", linewidth = .18, color = "grey38") + geom_line(linewidth = .26, alpha = .55, position = dodge) + geom_point(position = dodge, size = 1) +
    geom_point(data = deriv %>% filter(G0 <= .25 | G0 >= 1), shape = 21, fill = "white", color = "black", stroke = .3, size = 1.55, position = dodge) + scale_color_manual(values = fg2_ploidy_palette, guide = "none") + labs(x = "G0", y = "robust max derivative") + theme_fg2_recycled(6.5)
  auc <- x %>% filter(G0 <= 1); fits <- auc %>% group_by(ploidy_state) %>% group_modify(~{m <- lm(glucose_drawdown ~ live_auc_glucose_window, data = .x); xx <- seq(min(.x$live_auc_glucose_window), max(.x$live_auc_glucose_window), length.out = 80); tibble(live_auc_glucose_window = xx, glucose_drawdown = predict(m, newdata = data.frame(live_auc_glucose_window = xx)))})
  p2 <- ggplot(auc, aes(live_auc_glucose_window / 1000, glucose_drawdown, color = ploidy_state)) + geom_point(size = 1) + geom_point(data = auc %>% filter(abs(G0) < 1e-8), aes(fill = ploidy_state), shape = 21, color = "black", stroke = .42, size = 1.85) +
    geom_line(data = fits, aes(live_auc_glucose_window / 1000, glucose_drawdown, color = ploidy_state), linewidth = .46) + scale_color_manual(values = fg2_ploidy_palette, guide = "none") + scale_fill_manual(values = fg2_ploidy_palette, guide = "none") + labs(x = "live-cell AUC (thousands)", y = "glucose consumed") + theme_fg2_recycled(6.5)
  p3 <- ggplot(x, aes(G0_display, peak_total_yield_net, color = ploidy_state, group = ploidy_state)) + geom_hline(yintercept = 0, linewidth = .18, color = "grey72") + geom_line(linewidth = .26, alpha = .5, position = dodge) + geom_point(position = dodge, size = 1) + geom_point(data = x %>% filter(abs(G0 - 1) < 1e-8), shape = 21, fill = "white", color = "black", stroke = .32, size = 1.65, position = dodge) + scale_color_manual(values = fg2_ploidy_palette, name = "ploidy") + scale_y_continuous(labels = label_short) + labs(x = "G0", y = "net total-cell yield") + theme_fg2_recycled(6.5) + theme(legend.position = "bottom")
  (p1 | p2 | p3) + plot_layout(widths = c(1, 1.03, 1))
}

make_s3 <- function() {
  meta <- as_tibble(feature_release$feature_catalog) %>% select(feature_id, display_label, feature_class)
  x <- line_features %>% pivot_longer(all_of(meta$feature_id), names_to = "feature_id", values_to = "raw_value") %>% left_join(meta, by = "feature_id") %>%
    mutate(cellLine = factor(cellLine, levels = rev(line_order)), display_label = factor(display_label, levels = meta$display_label), ploidy_state = factor(ploidy_state, c("low", "high")))
  pairs <- x %>% group_by(cellLine, feature_id, display_label, feature_class) %>% summarise(low_value = raw_value[ploidy_state == "low"][1], high_value = raw_value[ploidy_state == "high"][1], .groups = "drop")
  ggplot() + geom_vline(xintercept = 0, linewidth = .2, color = "grey68") + geom_segment(data = pairs, aes(x = low_value, xend = high_value, y = cellLine, yend = cellLine), linewidth = .42, color = "grey55") +
    geom_point(data = x, aes(raw_value, cellLine, fill = ploidy_state, shape = ploidy_state), color = "grey18", size = 1.62, stroke = .22) +
    geom_point(data = x %>% filter(cellLine == "SUM-159-fuse"), aes(raw_value, cellLine, fill = ploidy_state, shape = ploidy_state), color = "black", stroke = .55, size = 2.35) +
    facet_wrap(~display_label, scales = "free_x", ncol = 3) + scale_fill_manual(values = fg2_ploidy_palette, name = "ploidy state") + scale_shape_manual(values = c(low = 21, high = 24), name = "ploidy state") +
    scale_x_continuous(labels = label_short) + labs(x = "raw feature value", y = NULL) + theme_fg2_recycled(6.8) + theme(strip.text = element_text(size = 6), axis.text.x = element_text(size = 5.5), panel.spacing = unit(.38, "lines"), legend.position = "bottom")
}

make_s4 <- function() {
  x <- condition_features %>% mutate(cellLine = factor(cellLine, line_order), G0_display = factor(format_g0(G0), levels = g0_levels), ploidy_state = factor(ploidy_state, c("low", "high")))
  dodge <- position_dodge(width = .34)
  base_theme <- theme_fg2_recycled(5.8) + theme(strip.text = element_text(size = 4.7), axis.text.x = element_text(size = 4.4, angle = 35, hjust = 1), axis.text.y = element_text(size = 4.4), axis.title = element_text(size = 5.3), panel.spacing = unit(.12, "lines"), legend.text = element_text(size = 5.2), legend.title = element_text(size = 5.3))
  p1 <- ggplot(x, aes(G0_display, robust_max_derivative, color = ploidy_state, group = ploidy_state)) + geom_hline(yintercept = 0, linewidth = .14, color = "grey74") + geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed", linewidth = .16, color = "grey38") + geom_line(linewidth = .22, alpha = .48, position = dodge) + geom_point(position = dodge, size = .65) +
    geom_point(data = x %>% filter(G0 <= .25 | G0 >= 1), shape = 21, fill = "white", color = "black", stroke = .25, size = 1.05, position = dodge) + facet_wrap(~cellLine, ncol = 1, scales = "free_y", labeller = as_labeller(line_short)) + scale_color_manual(values = fg2_ploidy_palette, guide = "none") + labs(x = "G0", y = "robust derivative") + base_theme
  auc <- x %>% filter(G0 <= 1); fits <- auc %>% group_by(cellLine, ploidy_state) %>% group_modify(~{m <- lm(glucose_drawdown ~ live_auc_glucose_window, data = .x); xx <- seq(min(.x$live_auc_glucose_window), max(.x$live_auc_glucose_window), length.out = 70); tibble(live_auc_glucose_window = xx, glucose_drawdown = predict(m, newdata = data.frame(live_auc_glucose_window = xx)))})
  p2 <- ggplot(auc, aes(live_auc_glucose_window / 1000, glucose_drawdown, color = ploidy_state)) + geom_point(size = .65, alpha = .88) + geom_point(data = auc %>% filter(abs(G0) < 1e-8), aes(fill = ploidy_state), shape = 21, color = "black", stroke = .28, size = 1.05) + geom_line(data = fits, aes(live_auc_glucose_window / 1000, glucose_drawdown, color = ploidy_state), linewidth = .28) + facet_wrap(~cellLine, ncol = 1, scales = "free_y", labeller = as_labeller(line_short)) + scale_color_manual(values = fg2_ploidy_palette, guide = "none") + scale_fill_manual(values = fg2_ploidy_palette, guide = "none") + scale_x_continuous(labels = label_short) + labs(x = "live-cell AUC (thousands)", y = "glucose consumed") + base_theme
  p3 <- ggplot(x, aes(G0_display, peak_total_yield_net, color = ploidy_state, group = ploidy_state)) + geom_hline(yintercept = 0, linewidth = .14, color = "grey74") + geom_line(linewidth = .22, alpha = .48, position = dodge) + geom_point(position = dodge, size = .65) + geom_point(data = x %>% filter(abs(G0 - 1) < 1e-8), shape = 21, fill = "white", color = "black", stroke = .25, size = 1.1, position = dodge) + facet_wrap(~cellLine, ncol = 1, scales = "free_y", labeller = as_labeller(line_short)) + scale_color_manual(values = fg2_ploidy_palette, name = "ploidy") + scale_y_continuous(labels = label_short) + labs(x = "G0", y = "net total-cell yield") + base_theme + theme(legend.position = "bottom")
  (p1 | p2 | p3) + plot_layout(widths = c(1.05, 1, 1.05))
}

apply_recycled_layout <- function(plan) {
  specs <- tribble(
    ~panel, ~x_in, ~y_in, ~width_in, ~height_in,
    "a", 0.00, 5.80, 7.00, 3.45,
    "b", 0.00, 0.04, 2.02, 5.62,
    "c", 2.14, 3.93, 1.98, 1.73,
    "d", 4.20, 3.93, 2.80, 1.73,
    "e", 2.14, 0.04, 2.38, 3.72,
    "f", 4.62, 0.04, 2.38, 3.72
  ) %>% mutate(
    figure = "Figure 1", subpanel_png = rel_path(file.path(subpanel_dir, paste0("figure_1_", panel, ".png"))),
    sx = 1, sy = 1, x_npc = x_in / 7, y_npc = y_in / 9.25, width_npc = width_in / 7, height_npc = height_in / 9.25,
    layout_width_in = 7, layout_height_in = 9.25
  )
  common <- intersect(names(plan), names(specs))
  missing <- setdiff(names(plan), names(specs))
  for (nm in missing) specs[[nm]] <- if (is.numeric(plan[[nm]])) 0 else "manual accepted v02 layout"
  bind_rows(plan %>% filter(figure != "Figure 1"), specs %>% select(all_of(names(plan))))
}
