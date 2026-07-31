#!/usr/bin/env Rscript

# FG4 current-data rebuild.
#
# This script intentionally persists no intermediate objects. Every posterior,
# size-context, ordination, distance, diagnostic, and layout transform is done
# in memory. The only outputs are the three requested final PNGs and the package
# timing TSV.
#
# Figure S13 semantic adaptation:
# The historical panel summarized five leave-one-line-out NUTS fits. Those fits
# were not rerun for red_a30_counts_20260722. The closest current posterior
# evidence is therefore used: five independently fitted single-line posteriors.
# The heatmap keeps the accepted parameter-by-line identity and reports, for
# each line, direction counts across the same five model structures plus the
# median posterior log2(high/low) effect. It must be interpreted as single-line
# posterior robustness, not leave-one-line-out robustness.

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
default_output_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
output_root <- normalizePath(
  arg_value("--output-root", default_output_root),
  winslash = "/",
  mustWork = TRUE
)
final_dir <- file.path(output_root, "final_images")
timing_dir <- file.path(output_root, "timings")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)

package_id <- "FG4_posterior_size_context"
run_started <- proc.time()[["elapsed"]]
timing_rows <- list()
add_timing <- function(operation, started, details) {
  timing_rows[[length(timing_rows) + 1L]] <<- data.frame(
    package_id = package_id,
    operation = operation,
    elapsed_seconds = round(proc.time()[["elapsed"]] - started, 3),
    details = details,
    stringsAsFactors = FALSE
  )
}

model_ids <- c(
  "1R_1P_0W_C0_M1",
  "1R_1P_1W_C0_M0",
  "2R_1P_0W_C0_M1",
  "2R_2P_0W_C0_M1",
  "2R_2P_1W_C0_M0"
)
model_alias <- c(
  "1R_1P_0W_C0_M1" = "1R",
  "1R_1P_1W_C0_M0" = "1R,W(m)",
  "2R_1P_0W_C0_M1" = "2R(a)",
  "2R_2P_0W_C0_M1" = "2R(f)",
  "2R_2P_1W_C0_M0" = "2R(f),W(m)"
)
model_order <- unname(model_alias[model_ids])
line_order <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
four_lines <- line_order[1:4]
param_map <- c("Y_R[1]" = "y1", "ae[1]" = "v1", "ah[1]" = "K1", "m" = "m")
param_order <- c("y1", "v1", "K1", "m")
param_labels <- c(
  y1 = "y1 yield",
  v1 = "v1 uptake",
  K1 = "K1 response",
  m = "m maintenance/\ndeath"
)
line_colors <- c(
  "MCF10A" = "#0072B2",
  "MDA-MB-231" = "#D55E00",
  "SNU668" = "#009E73",
  "SUM-159-chem" = "#CC79A7",
  "SUM-159-fuse" = "#E69F00"
)
parameter_colors <- c(y1 = "#0072B2", v1 = "#D55E00", K1 = "#009E73", m = "#CC79A7")
direction_colors <- c(positive = "#B2182B", overlaps_zero = "#8C8C8C", negative = "#2166AC")

model_root <- file.path(project_root, "data/modelling/gpath/v1")
if (!dir.exists(model_root)) {
  model_root <- file.path(project_root, "data/modelling/gpath_v1")
}
release_root <- file.path(model_root, "red_a30_counts_20260722")
posterior_root <- file.path(release_root, "derived/posterior")
parameter_root <- file.path(posterior_root, "parameters")
qc_path <- file.path(posterior_root, "qc.Rds")
optimization_path <- file.path(release_root, "derived/optimization/assessment.Rds")
area_path <- file.path(
  project_root,
  "data/image_processing_runs/full_segmentation_classification_nuclear",
  "run_20260721_163410/summaries/cell_nuclear_area_summaries.Rds"
)
stan_paths <- file.path(
  release_root,
  "datasets",
  c(
    "all_lines", "loo_exclude_sum_159_fuse",
    "single_mcf10a", "single_mda_mb_231", "single_snu668",
    "single_sum_159_chem", "single_sum_159_fuse"
  ),
  "stan_data.Rds"
)
required_inputs <- c(
  file.path(parameter_root, c(
    "all_lines.Rds", "loo_exclude_sum_159_fuse.Rds",
    "single_mcf10a.Rds", "single_mda_mb_231.Rds", "single_snu668.Rds",
    "single_sum_159_chem.Rds", "single_sum_159_fuse.Rds"
  )),
  qc_path, optimization_path, area_path, stan_paths
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing canonical FG4 input(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

quant <- function(x, p) {
  as.numeric(stats::quantile(x[is.finite(x)], p, names = FALSE, type = 8))
}
classify_interval <- function(lo, hi) {
  fifelse(lo > 0, "positive", fifelse(hi < 0, "negative", "overlaps_zero"))
}
theme_panel <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, color = "grey90"),
      strip.background = element_rect(fill = "grey94", color = "grey72", linewidth = 0.25),
      strip.text = element_text(face = "bold", color = "grey12"),
      axis.text = element_text(color = "grey15"),
      axis.title = element_text(color = "grey15"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(2, 3, 2, 3)
    )
}

read_posterior_reduced <- function(path, keep_draws = TRUE) {
  obj <- readRDS(path)
  if (!identical(obj$schema_version, 1L) || is.null(obj$draws)) {
    stop("Unexpected posterior parameter schema: ", path, call. = FALSE)
  }
  dt <- as.data.table(obj$draws)
  dt <- dt[
    model_id %chin% model_ids &
      parameter %chin% names(param_map) &
      is.finite(value) & value > 0,
    .(model_id, draw_id, line_name, ploidy, parameter, value)
  ]
  dt[, param := unname(param_map[parameter])]
  dt[, log2_value := log2(value)]
  state <- dt[, .(
    median_log2 = median(log2_value),
    q025_log2 = quant(log2_value, 0.025),
    q975_log2 = quant(log2_value, 0.975)
  ), by = .(model_id, line_name, ploidy, param)]
  effect <- dcast(
    dt,
    model_id + draw_id + line_name + param ~ ploidy,
    value.var = "log2_value"
  )
  if (!all(c("low", "high") %in% names(effect))) {
    stop("Posterior states low/high were not both present in ", path, call. = FALSE)
  }
  effect[, effect := high - low]
  effect <- effect[, .(model_id, draw_id, line_name, param, effect)]
  effect_summary <- effect[, .(
    median = median(effect),
    q025 = quant(effect, 0.025),
    q25 = quant(effect, 0.25),
    q75 = quant(effect, 0.75),
    q975 = quant(effect, 0.975)
  ), by = .(model_id, line_name, param)]
  effect_summary[, direction := classify_interval(q025, q975)]
  rm(obj, dt)
  invisible(gc())
  list(
    state = state,
    effect = if (keep_draws) effect else NULL,
    summary = effect_summary
  )
}

load_started <- proc.time()[["elapsed"]]
no_sum <- read_posterior_reduced(
  file.path(parameter_root, "loo_exclude_sum_159_fuse.Rds"),
  keep_draws = TRUE
)
all_line <- read_posterior_reduced(
  file.path(parameter_root, "all_lines.Rds"),
  keep_draws = FALSE
)
qc <- readRDS(qc_path)
optimization <- readRDS(optimization_path)
area_obj <- readRDS(area_path)
stan_data <- lapply(stan_paths, readRDS)
names(stan_data) <- basename(dirname(stan_paths))
add_timing(
  "load_inputs",
  load_started,
  "Loaded canonical posterior parameter RDS files, posterior QC, optimization assessment, seven stan_data RDS files, and the 20260721 area summary."
)

if (!all(model_ids %in% optimization$pareto$model_id)) {
  stop("The five accepted model structures are not all present in the current optimization Pareto assessment.", call. = FALSE)
}
expected_maps <- list(
  all_lines = line_order,
  loo_exclude_sum_159_fuse = four_lines,
  single_mcf10a = "MCF10A",
  single_mda_mb_231 = "MDA-MB-231",
  single_snu668 = "SNU668",
  single_sum_159_chem = "SUM-159-chem",
  single_sum_159_fuse = "SUM-159-fuse"
)
for (nm in names(expected_maps)) {
  observed <- names(stan_data[[nm]]$line_map)
  if (!setequal(observed, expected_maps[[nm]])) {
    stop("Unexpected line map in datasets/", nm, "/stan_data.Rds", call. = FALSE)
  }
}
rm(stan_data)
invisible(gc())

posterior_transform_started <- proc.time()[["elapsed"]]
no_sum$summary[, alias := factor(unname(model_alias[model_id]), levels = model_order)]
all_line$summary[, alias := factor(unname(model_alias[model_id]), levels = model_order)]
no_sum$state[, alias := factor(unname(model_alias[model_id]), levels = model_order)]
no_sum$effect[, alias := factor(unname(model_alias[model_id]), levels = model_order)]
add_timing(
  "transform_posterior",
  posterior_transform_started,
  "Reduced all-line and SUM-159-fuse-excluded posterior draws to state summaries and paired high-minus-low effects for five accepted model structures."
)

area_transform_started <- proc.time()[["elapsed"]]
area_dt <- as.data.table(area_obj$replicate_time)
area_dt <- area_dt[
  object_scope == "alive_area_pass" &
    hours == 0 &
    cellLine %chin% four_lines &
    is.finite(segmented_area_px_q50) &
    segmented_area_px_q50 > 0,
  .(
    cellLine, experiment, plateID, glucose, hours, ploidy,
    area_q50 = segmented_area_px_q50
  )
]
area_dt[, state := fcase(
  cellLine == "MCF10A" & ploidy == "4N", "high",
  cellLine == "MCF10A" & ploidy == "2N", "low",
  cellLine == "MDA-MB-231" & ploidy == "parental", "high",
  cellLine == "MDA-MB-231" & ploidy == "3N", "low",
  cellLine == "SNU668" & ploidy == "high", "high",
  cellLine == "SNU668" & ploidy == "low", "low",
  cellLine == "SUM-159-chem" & ploidy == "4N", "high",
  cellLine == "SUM-159-chem" & ploidy == "2N", "low",
  default = NA_character_
)]
area_dt <- area_dt[!is.na(state)]
# Ploidy series occupy distinct plates in the canonical summary, so a literal
# plate-paired ratio is unavailable. Use the central high-minus-low log2 area
# contrast and a conservative empirical envelope (high q025 - low q975 through
# high q975 - low q025). This is a direct descriptive size context, not the
# historical crowding-adjusted morphology model.
area_state <- area_dt[, .(
  center = median(log2(area_q50)),
  q025 = quant(log2(area_q50), 0.025),
  q975 = quant(log2(area_q50), 0.975),
  n_contexts = .N
), by = .(cellLine, state)]
area_state_wide <- dcast(
  area_state,
  cellLine ~ state,
  value.var = c("center", "q025", "q975", "n_contexts")
)
area_summary <- area_state_wide[, .(
  line_name = cellLine,
  area_shift = center_high - center_low,
  area_lo = q025_high - q975_low,
  area_hi = q975_high - q025_low,
  n_contexts_per_state = pmin(n_contexts_high, n_contexts_low)
)]
if (!setequal(area_summary$line_name, four_lines) || any(area_summary$n_contexts_per_state < 2L)) {
  stop("Canonical area summary did not provide baseline contexts for all four FG4 lines.", call. = FALSE)
}
rm(area_obj, area_dt, area_state, area_state_wide)
invisible(gc())
add_timing(
  "transform_area_context",
  area_transform_started,
  "Computed unpaired baseline high-minus-low log2 segmented-area-q50 contrasts and conservative empirical envelopes in memory; volume proxy uses area_ratio^(3/2)."
)

figure4_transform_started <- proc.time()[["elapsed"]]
coords <- dcast(
  no_sum$state,
  model_id + alias + line_name + ploidy ~ param,
  value.var = "median_log2"
)
setorderv(coords, c("model_id", "line_name", "ploidy"))
param_sd <- vapply(param_order, function(p) stats::sd(coords[[p]]), numeric(1))
if (any(!is.finite(param_sd)) || any(param_sd <= 0)) {
  stop("Non-finite global parameter scale in Figure 4 ordination.", call. = FALSE)
}
for (p in param_order) {
  coords[, (paste0(p, "_z")) := (get(p) - mean(get(p))) / param_sd[[p]]]
}
zmat <- as.matrix(coords[, paste0(param_order, "_z"), with = FALSE])
ordination <- cmdscale(stats::dist(zmat), k = 2, eig = TRUE, add = TRUE)
coords[, `:=`(axis_1 = ordination$points[, 1], axis_2 = ordination$points[, 2])]
if (stats::cor(coords$K1_z, coords$axis_1) < 0) coords[, axis_1 := -axis_1]
if (stats::cor(coords$m_z, coords$axis_2) < 0) coords[, axis_2 := -axis_2]
positive_eigen <- ordination$eig[ordination$eig > 0]
axis_var <- positive_eigen[1:2] / sum(positive_eigen)

low_coords <- coords[ploidy == "low"]
high_coords <- coords[ploidy == "high"]
segments <- merge(
  low_coords[, .(model_id, alias, line_name, axis_1_low = axis_1, axis_2_low = axis_2)],
  high_coords[, .(model_id, alias, line_name, axis_1_high = axis_1, axis_2_high = axis_2)],
  by = c("model_id", "alias", "line_name")
)
centroids <- coords[, .(axis_1 = mean(axis_1), axis_2 = mean(axis_2)), by = alias]

axis_fit_1 <- lm(coords$axis_1 ~ zmat)
axis_fit_2 <- lm(coords$axis_2 ~ zmat)
projection <- cbind(
  axis_1 = coef(axis_fit_1)[-1],
  axis_2 = coef(axis_fit_2)[-1]
)
effect_wide <- dcast(
  no_sum$effect,
  model_id + alias + draw_id + line_name ~ param,
  value.var = "effect"
)
delta_z <- sweep(as.matrix(effect_wide[, ..param_order]), 2, param_sd, "/")
delta_axis <- delta_z %*% projection
effect_wide[, `:=`(delta_axis_1 = delta_axis[, 1], delta_axis_2 = delta_axis[, 2])]
ellipse_points <- function(df, level = 0.80, n = 81L) {
  center <- high_coords[
    model_id == df$model_id[[1]] & line_name == df$line_name[[1]],
    c("axis_1", "axis_2")
  ]
  covariance <- stats::cov(cbind(df$delta_axis_1, df$delta_axis_2))
  eig <- eigen(covariance, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  theta <- seq(0, 2 * pi, length.out = n)
  offsets <- eig$vectors %*% diag(sqrt(eig$values), 2) %*%
    rbind(cos(theta), sin(theta)) * sqrt(qchisq(level, 2))
  data.table(
    model_id = df$model_id[[1]],
    alias = df$alias[[1]],
    line_name = df$line_name[[1]],
    axis_1 = center$axis_1[[1]] + offsets[1, ],
    axis_2 = center$axis_2[[1]] + offsets[2, ]
  )
}
ellipses <- rbindlist(lapply(
  split(effect_wide, interaction(effect_wide$model_id, effect_wide$line_name, drop = TRUE)),
  ellipse_points
))
ellipses[, alias := factor(as.character(alias), levels = model_order)]

z_cols <- paste0(param_order, "_z")
line_profiles <- dcast(
  melt(
    coords,
    id.vars = c("model_id", "alias", "line_name", "ploidy"),
    measure.vars = z_cols,
    variable.name = "param",
    value.name = "z"
  )[, param := sub("_z$", "", param)],
  model_id + alias + line_name ~ ploidy + param,
  value.var = "z"
)
profile_cols <- setdiff(names(line_profiles), c("model_id", "alias", "line_name"))
line_distance_rows <- list()
k <- 1L
for (mid in model_ids) {
  one <- line_profiles[model_id == mid]
  for (li in four_lines) for (lj in four_lines) {
    vi <- as.numeric(one[line_name == li, ..profile_cols])
    vj <- as.numeric(one[line_name == lj, ..profile_cols])
    line_distance_rows[[k]] <- data.table(
      model_id = mid, line_i = li, line_j = lj,
      distance = sqrt(mean((vi - vj)^2))
    )
    k <- k + 1L
  }
}
line_distances <- rbindlist(line_distance_rows)[, .(
  median_distance = median(distance)
), by = .(line_i, line_j)]

ploidy_distances <- merge(
  low_coords[, c("model_id", "alias", "line_name", z_cols), with = FALSE],
  high_coords[, c("model_id", "alias", "line_name", z_cols), with = FALSE],
  by = c("model_id", "alias", "line_name"),
  suffixes = c("_low", "_high")
)
ploidy_distances[, distance := sqrt(rowMeans(vapply(
  param_order,
  function(p) (get(paste0(p, "_z_high")) - get(paste0(p, "_z_low")))^2,
  numeric(.N)
)))]

model_centers <- coords[, lapply(.SD, mean), by = alias, .SDcols = z_cols]
centroid_rows <- list()
k <- 1L
for (mi in model_order) for (mj in model_order) {
  vi <- as.numeric(model_centers[alias == mi, ..z_cols])
  vj <- as.numeric(model_centers[alias == mj, ..z_cols])
  centroid_rows[[k]] <- data.table(
    model_i = mi, model_j = mj,
    distance = sqrt(mean((vi - vj)^2))
  )
  k <- k + 1L
}
centroid_distances <- rbindlist(centroid_rows)

components <- merge(
  low_coords[, c("model_id", "alias", "line_name", z_cols), with = FALSE],
  high_coords[, c("model_id", "alias", "line_name", z_cols), with = FALSE],
  by = c("model_id", "alias", "line_name"),
  suffixes = c("_low", "_high")
)
components <- rbindlist(lapply(param_order, function(p) {
  components[, .(
    model_id, alias, line_name, param = p,
    contribution = get(paste0(p, "_z_high")) - get(paste0(p, "_z_low"))
  )]
}))
add_timing(
  "transform_figure_4",
  figure4_transform_started,
  "Recomputed global-z PCoA, posterior displacement ellipses, exact line/ploidy/model distances, and signed parameter contributions from current no-SUM-equivalent posterior data."
)

axis_limits <- function(x) {
  r <- range(x, finite = TRUE)
  p <- max(diff(r) * 0.10, 0.15)
  r + c(-p, p)
}
x_lim <- axis_limits(c(coords$axis_1, ellipses$axis_1))
y_lim <- axis_limits(c(coords$axis_2, ellipses$axis_2))
model_panel <- function(alias_value) {
  seg <- segments[alias == alias_value]
  el <- ellipses[alias == alias_value]
  cen <- centroids[alias == alias_value]
  ggplot() +
    geom_polygon(
      data = el,
      aes(axis_1, axis_2, group = line_name, fill = line_name, color = line_name),
      alpha = 0.15, linewidth = 0.45
    ) +
    geom_segment(
      data = seg,
      aes(axis_1_low, axis_2_low, xend = axis_1_high, yend = axis_2_high, color = line_name),
      linewidth = 0.65,
      arrow = arrow(length = unit(0.055, "in"), type = "closed")
    ) +
    geom_point(data = cen, aes(axis_1, axis_2), size = 1.8) +
    annotate("text", x = x_lim[1] + 0.04 * diff(x_lim), y = y_lim[2] - 0.04 * diff(y_lim),
             label = alias_value, hjust = 0, vjust = 1, fontface = "bold", size = 2.7) +
    scale_color_manual(values = line_colors, guide = "none") +
    scale_fill_manual(values = line_colors, guide = "none") +
    scale_x_continuous(sprintf("PCoA 1 (%.1f%%)", 100 * axis_var[1]), limits = x_lim) +
    scale_y_continuous(sprintf("PCoA 2 (%.1f%%)", 100 * axis_var[2]), limits = y_lim) +
    coord_equal() +
    theme_panel(6.7) +
    theme(axis.text = element_text(size = 5.8), axis.title = element_text(size = 6.3))
}
legend_df <- data.table(
  line_name = factor(four_lines, levels = four_lines),
  y = rev(seq_along(four_lines))
)
ellipse_key <- data.table(
  x = 0.23 + 0.15 * cos(seq(0, 2 * pi, length.out = 81)),
  y = -1 + 0.20 * sin(seq(0, 2 * pi, length.out = 81))
)
p_legend <- ggplot(legend_df) +
  annotate("text", x = 0, y = 5.0, label = "Line context", hjust = 0, fontface = "bold", size = 2.8) +
  geom_segment(
    aes(x = 0.05, xend = 0.42, y = y, yend = y, color = line_name),
    linewidth = 0.7,
    arrow = arrow(length = unit(0.05, "in"), type = "closed")
  ) +
  geom_text(aes(x = 0.55, y = y, label = line_name), hjust = 0, size = 2.35) +
  annotate("point", x = 0.23, y = 0, size = 1.8) +
  annotate("text", x = 0.55, y = 0, label = "Model centroid", hjust = 0, size = 2.35) +
  geom_polygon(
    data = ellipse_key,
    aes(x, y),
    inherit.aes = FALSE,
    fill = "grey70",
    color = "grey45",
    alpha = 0.22,
    linewidth = 0.35
  ) +
  annotate("text", x = 0.55, y = -1, label = "80% posterior displacement ellipse",
           hjust = 0, size = 2.15) +
  scale_color_manual(values = line_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 2.8), ylim = c(-1.45, 5.3), clip = "off") +
  theme_void()
p4a <- plot_grid(
  plotlist = c(lapply(model_order, model_panel), list(p_legend)),
  ncol = 2,
  align = "hv"
)

vector_rows <- rbindlist(lapply(param_order, function(p) {
  data.table(
    param = p,
    axis_1 = cor(coords[[paste0(p, "_z")]], coords$axis_1),
    axis_2 = cor(coords[[paste0(p, "_z")]], coords$axis_2)
  )
}))
p4b <- ggplot(vector_rows) +
  annotate("path", x = cos(seq(0, 2 * pi, length.out = 201)),
           y = sin(seq(0, 2 * pi, length.out = 201)), color = "grey80", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_segment(
    aes(x = 0, y = 0, xend = axis_1, yend = axis_2, color = param),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.06, "in"), type = "closed")
  ) +
  geom_text(aes(axis_1 * 1.08, axis_2 * 1.08, label = param, color = param),
            fontface = "bold", size = 2.7) +
  scale_color_manual(values = parameter_colors, guide = "none") +
  coord_equal(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", 100 * axis_var[1]),
       y = sprintf("PCoA 2 (%.1f%%)", 100 * axis_var[2])) +
  theme_panel(7)

tile_text_color <- function(x, threshold) ifelse(x >= threshold, "white", "grey10")
line_distances[, `:=`(
  line_i = factor(line_i, levels = four_lines),
  line_j = factor(line_j, levels = rev(four_lines)),
  text_color = tile_text_color(median_distance, 0.62 * max(median_distance))
)]
p4c <- ggplot(line_distances, aes(line_i, line_j, fill = median_distance)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", median_distance), color = text_color),
            fontface = "bold", size = 2.25) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F4F7FA", high = "#2F6B67", name = "Exact 8D\nRMS distance") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 32, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.8)
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(1.15, "in"), barheight = unit(0.07, "in")))

ploidy_distances[, `:=`(
  alias = factor(as.character(alias), levels = rev(model_order)),
  line_name = factor(line_name, levels = four_lines),
  text_color = tile_text_color(distance, 0.62 * max(distance))
)]
p4d <- ggplot(ploidy_distances, aes(line_name, alias, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", distance), color = text_color),
            fontface = "bold", size = 2.15) +
  scale_color_identity() +
  scale_fill_gradient(low = "#FAF7FA", high = "#8C2D60", name = "Exact 4D\nRMS distance") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 32, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 5.8),
    legend.text = element_text(size = 5.5)
  )

centroid_distances[, `:=`(
  model_i = factor(model_i, levels = model_order),
  model_j = factor(model_j, levels = rev(model_order)),
  text_color = tile_text_color(distance, 0.62 * max(distance))
)]
p4e <- ggplot(centroid_distances, aes(model_i, model_j, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", distance), color = text_color),
            fontface = "bold", size = 2.05) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F7FBFF", high = "#174A7E", name = "Exact 4D\nRMS distance") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 32, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 5.8),
    legend.text = element_text(size = 5.5)
  )

components[, `:=`(
  alias = factor(as.character(alias), levels = rev(model_order)),
  line_name = factor(line_name, levels = four_lines),
  param = factor(param, levels = param_order)
)]
component_limit <- max(abs(components$contribution))
p4f <- ggplot(components, aes(line_name, alias, fill = contribution)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", contribution)), size = 1.85) +
  facet_grid(. ~ param) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#B5433F",
    midpoint = 0, limits = c(-component_limit, component_limit),
    name = "Standardized\nhigh - low"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.7) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 32, hjust = 1, size = 5.8),
    axis.text.y = element_text(size = 5.8),
    legend.title = element_text(size = 5.8),
    legend.text = element_text(size = 5.5)
  )

figure_4 <- ggdraw() +
  draw_plot(p4a, x = 0.015, y = 0.43, width = 0.59, height = 0.555) +
  draw_plot(p4b, x = 0.625, y = 0.68, width = 0.36, height = 0.295) +
  draw_plot(p4c, x = 0.625, y = 0.43, width = 0.36, height = 0.235) +
  draw_plot(p4d, x = 0.015, y = 0.245, width = 0.47, height = 0.18) +
  draw_plot(p4e, x = 0.515, y = 0.245, width = 0.47, height = 0.18) +
  draw_plot(p4f, x = 0.015, y = 0.015, width = 0.97, height = 0.225) +
  draw_plot_label(c("a", "b", "c", "d", "e", "f"),
                  x = c(0.003, 0.61, 0.61, 0.003, 0.50, 0.003),
                  y = c(0.995, 0.995, 0.675, 0.435, 0.435, 0.245),
                  size = 13, fontface = "bold")

s10_transform_started <- proc.time()[["elapsed"]]
size_base <- copy(no_sum$summary[line_name %chin% four_lines & param %chin% c("y1", "v1")])
size_base <- merge(size_base, area_summary, by = "line_name", all.x = TRUE)
size_rows <- rbindlist(list(
  size_base[, `:=`(quantity = "Per-cell\nposterior", shift = 0, shift_lo = 0, shift_hi = 0)],
  copy(size_base)[, `:=`(
    quantity = "Area\nfilled",
    shift = area_shift, shift_lo = area_lo, shift_hi = area_hi
  )],
  copy(size_base)[, `:=`(
    quantity = "Volume\nproxy",
    shift = 1.5 * area_shift, shift_lo = 1.5 * area_lo, shift_hi = 1.5 * area_hi
  )]
), use.names = TRUE)
size_rows[, `:=`(
  median_filled = median + shift,
  q025_filled = q025 + shift_lo,
  q975_filled = q975 + shift_hi,
  interval_class = classify_interval(q025 + shift_lo, q975 + shift_hi)
)]
size_summary <- size_rows[, .(
  median = median(median_filled),
  q25 = quant(median_filled, 0.25),
  q75 = quant(median_filled, 0.75),
  lo = min(q025_filled),
  hi = max(q975_filled)
), by = .(line_name, param, quantity)]
size_summary[, direction := classify_interval(lo, hi)]
size_rows[, `:=`(
  line_name = factor(line_name, levels = rev(four_lines)),
  param_label = factor(unname(param_labels[param]), levels = unname(param_labels[c("y1", "v1")])),
  quantity = factor(quantity, levels = c("Per-cell\nposterior", "Area\nfilled", "Volume\nproxy")),
  interval_display = factor(
    interval_class,
    levels = c("negative", "overlaps_zero", "positive"),
    labels = c("< 0", "overlaps 0", "> 0")
  )
)]
size_summary[, `:=`(
  line_name = factor(line_name, levels = rev(four_lines)),
  param_label = factor(unname(param_labels[param]), levels = unname(param_labels[c("y1", "v1")])),
  quantity = factor(quantity, levels = c("Per-cell\nposterior", "Area\nfilled", "Volume\nproxy"))
)]
p10a <- ggplot(size_summary, aes(median, line_name)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_segment(aes(x = lo, xend = hi, yend = line_name), linewidth = 0.5, color = "grey72") +
  geom_segment(aes(x = q25, xend = q75, yend = line_name), linewidth = 1.25, color = "grey22") +
  geom_point(
    data = size_rows,
    aes(x = median_filled, y = line_name, color = interval_display),
    inherit.aes = FALSE, size = 1.35, alpha = 0.75
  ) +
  geom_point(aes(fill = direction), shape = 23, size = 2.4, color = "white", stroke = 0.35) +
  facet_grid(quantity ~ param_label) +
  scale_color_manual(
    values = c("< 0" = "#2166AC", "overlaps 0" = "#8C8C8C", "> 0" = "#B2182B"),
    name = NULL
  ) +
  scale_fill_manual(values = direction_colors, guide = "none") +
  labs(x = "Displayed effect, log2 high / low\n(size-filled rows add area or volume ratio)", y = NULL) +
  theme_panel(6.4) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 6.6),
    strip.text.y = element_text(size = 6.0, angle = 0),
    axis.text = element_text(size = 5.8),
    axis.title.x = element_text(size = 6.2)
  ) +
  guides(color = guide_legend(nrow = 1))

compare_rows <- rbindlist(list(
  copy(no_sum$summary)[, scope := "no-SUM primary\n4 lines"],
  copy(all_line$summary)[, scope := "all-line sensitivity\n5 lines"]
), use.names = TRUE)
compare_rows <- compare_rows[param %chin% param_order]
compare_consensus <- compare_rows[, .(
  median = median(median),
  q25 = quant(median, 0.25),
  q75 = quant(median, 0.75)
), by = .(scope, line_name, param)]
compare_consensus[, `:=`(
  scope = factor(scope, levels = c("no-SUM primary\n4 lines", "all-line sensitivity\n5 lines")),
  line_name = factor(line_name, levels = rev(line_order)),
  param_label = factor(
    unname(c(y1 = "y1 yield", v1 = "v1 uptake", K1 = "K1 response", m = "m")[param]),
    levels = c("y1 yield", "v1 uptake", "K1 response", "m")
  )
)]
p10b <- ggplot(compare_consensus, aes(median, line_name, color = scope)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0, linewidth = 0.8,
                 position = position_dodge(width = 0.38)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.38)) +
  facet_grid(. ~ param_label) +
  scale_color_manual(values = c(
    "no-SUM primary\n4 lines" = "#176D6A",
    "all-line sensitivity\n5 lines" = "#7A4FA3"
  ), name = NULL) +
  labs(x = "Posterior effect, log2 high / low", y = NULL) +
  theme_panel(6.3) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 6.0),
    axis.text = element_text(size = 5.5),
    axis.title.x = element_text(size = 6.0),
    legend.text = element_text(size = 5.5)
  )

fit_qc <- as.data.table(qc$fit_summary)[
  model_id %chin% model_ids &
    dataset_id %chin% c(
      "loo_exclude_sum_159_fuse", "all_lines",
      "single_mcf10a", "single_mda_mb_231", "single_snu668",
      "single_sum_159_chem", "single_sum_159_fuse"
    )
]
fit_qc[, `:=`(
  alias = factor(unname(model_alias[model_id]), levels = rev(model_order)),
  inventory_group = factor(
    dataset_id,
    levels = c(
      "loo_exclude_sum_159_fuse", "all_lines",
      "single_mcf10a", "single_mda_mb_231", "single_snu668",
      "single_sum_159_chem", "single_sum_159_fuse"
    ),
    labels = c(
      "no-SUM\nposterior", "all-line\nposterior",
      "MCF10A\nsingle-line", "MDA-MB-231\nsingle-line", "SNU668\nsingle-line",
      "SUM-159-chem\nsingle-line", "SUM-159-fuse\nsingle-line"
    )
  ),
  rhat_flag = as.integer(max_rhat > 1.01),
  bulk_flag = as.integer(min_ess_bulk < 400),
  tail_flag = as.integer(min_ess_tail < 400)
)]
inventory <- melt(
  fit_qc,
  id.vars = c("dataset_id", "inventory_group", "model_id", "alias"),
  measure.vars = c("rhat_flag", "bulk_flag", "tail_flag"),
  variable.name = "metric",
  value.name = "value"
)
candidate_inventory <- unique(fit_qc[, .(
  dataset_id, inventory_group, model_id, alias,
  metric = factor("candidate", levels = c("candidate", "rhat_flag", "bulk_flag", "tail_flag")),
  value = 0L
)])
inventory[, metric := factor(metric, levels = c("candidate", "rhat_flag", "bulk_flag", "tail_flag"))]
inventory <- rbindlist(list(candidate_inventory, inventory), use.names = TRUE)
inventory[, `:=`(
  display_metric = factor(
    metric,
    levels = c("candidate", "rhat_flag", "bulk_flag", "tail_flag"),
    labels = c("Cand", "Rhat", "B-ESS", "T-ESS")
  ),
  label = fifelse(metric == "candidate", "yes", as.character(value)),
  status = fcase(
    metric == "candidate", "candidate",
    value == 0, "0 flags",
    value == 1, "minor",
    default = "flag"
  )
)]
p10c <- ggplot(inventory, aes(display_metric, alias, fill = status)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 2.75) +
  facet_grid(inventory_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(candidate = "#78A870", "0 flags" = "#F4F4F4", minor = "#F2CC5C", flag = "#B2182B"),
    breaks = c("candidate", "0 flags", "minor", "flag"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8.0) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 0, size = 8.0, lineheight = 0.86),
    axis.text.x = element_text(size = 8.0, angle = 40, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.0),
    legend.position = "bottom",
    legend.text = element_text(size = 8.0),
    panel.spacing.y = unit(0.045, "in"),
    plot.margin = margin(2, 3, 2, 8)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

figure_s10 <- ggdraw() +
  draw_plot(p10a, x = 0.015, y = 0.405, width = 0.55, height = 0.58) +
  draw_plot(p10b, x = 0.015, y = 0.02, width = 0.55, height = 0.36) +
  draw_plot(p10c, x = 0.585, y = 0.02, width = 0.40, height = 0.965) +
  draw_plot_label("a", x = 0.006, y = 0.985, hjust = 0, vjust = 1,
                  size = 13, fontface = "bold") +
  draw_plot_label("b", x = 0.006, y = 0.395, hjust = 0, vjust = 1,
                  size = 13, fontface = "bold") +
  draw_plot_label("c", x = 0.57, y = 0.985, hjust = 0, vjust = 1,
                  size = 13, fontface = "bold")
add_timing(
  "transform_figure_s10",
  s10_transform_started,
  "Built size-context effects from current no-SUM posterior plus canonical area ratios, all-line sensitivity, and current posterior-QC inventory."
)

s13_transform_started <- proc.time()[["elapsed"]]
single_ids <- c(
  "single_mcf10a", "single_mda_mb_231", "single_snu668",
  "single_sum_159_chem", "single_sum_159_fuse"
)
single_paths <- file.path(parameter_root, paste0(single_ids, ".Rds"))
single_effect_rows <- list()
for (i in seq_along(single_ids)) {
  reduced <- read_posterior_reduced(single_paths[[i]], keep_draws = FALSE)
  one <- reduced$summary
  one[, dataset_id := single_ids[[i]]]
  single_effect_rows[[i]] <- one
  rm(reduced)
  invisible(gc())
}
single_effects <- rbindlist(single_effect_rows)
single_cells <- single_effects[, .(
  positive = sum(direction == "positive"),
  overlap = sum(direction == "overlaps_zero"),
  negative = sum(direction == "negative"),
  median_effect = median(median)
), by = .(dataset_id, line_name, param)]
single_qc <- fit_qc[dataset_id %chin% single_ids, .(
  n_pass = sum(
    divergences == 0 &
      max_rhat <= 1.01 &
      min_ess_bulk >= 400 &
      min_ess_tail >= 400
  )
), by = dataset_id]
single_cells <- merge(single_cells, single_qc, by = "dataset_id")
single_cells[, `:=`(
  line_label = factor(
    paste0(line_name, "\n", n_pass, "/5 pass"),
    levels = rev(paste0(line_order, "\n", single_qc$n_pass[match(
      c("single_mcf10a", "single_mda_mb_231", "single_snu668",
        "single_sum_159_chem", "single_sum_159_fuse"),
      single_qc$dataset_id
    )], "/5 pass"))
  ),
  param_label = factor(unname(param_labels[param]), levels = unname(param_labels)),
  cell_label = sprintf("%d/%d/%d\n%+.2f", positive, overlap, negative, median_effect)
)]
s13_limit <- max(abs(single_cells$median_effect))
p13 <- ggplot(single_cells, aes(param_label, line_label, fill = median_effect)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = cell_label), size = 2.65, lineheight = 0.92) +
  scale_fill_gradient2(
    low = "#3B95C3", mid = "#BDBDBD", high = "#D85C4A",
    midpoint = 0, limits = c(-s13_limit, s13_limit),
    name = "Median posterior\nlog2 high / low"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 32, hjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.5),
    legend.position = "bottom",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2),
    plot.margin = margin(4, 4, 3, 5)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    barwidth = unit(1.5, "in"),
    barheight = unit(0.08, "in")
  ))
figure_s13 <- ggdraw() +
  draw_plot(p13, x = 0.055, y = 0.04, width = 0.93, height = 0.94) +
  draw_plot_label("a", x = 0.01, y = 0.985, size = 13, fontface = "bold")
add_timing(
  "transform_figure_s13",
  s13_transform_started,
  "Replaced unavailable five-LOO NUTS evidence with five current single-line posteriors; cells show +/0/- model counts and median posterior effect."
)

render_started <- proc.time()[["elapsed"]]
f4_started <- proc.time()[["elapsed"]]
ggsave(
  file.path(final_dir, "figure_4.png"),
  figure_4,
  width = 7,
  height = 9.25,
  units = "in",
  dpi = 360,
  bg = "white",
  limitsize = FALSE
)
add_timing("render_figure_4", f4_started, "Rendered 7.00 x 9.25 inch PNG at 360 dpi.")

s10_started <- proc.time()[["elapsed"]]
ggsave(
  file.path(final_dir, "figure_s10.png"),
  figure_s10,
  width = 7,
  height = 9.25,
  units = "in",
  dpi = 360,
  bg = "white",
  limitsize = FALSE
)
add_timing("render_figure_s10", s10_started, "Rendered 7.00 x 9.25 inch PNG at 360 dpi.")

s13_started <- proc.time()[["elapsed"]]
ggsave(
  file.path(final_dir, "figure_s13.png"),
  figure_s13,
  width = 4.2,
  height = 3.3,
  units = "in",
  dpi = 360,
  bg = "white",
  limitsize = FALSE
)
add_timing("render_figure_s13", s13_started, "Rendered 4.20 x 3.30 inch PNG at 360 dpi.")
add_timing("render_all", render_started, "Rendered the three final composites directly from in-memory plot objects.")

add_timing(
  "total",
  run_started,
  "Complete FG4 current-data rebuild; no intermediate artifacts persisted."
)
timing_path <- file.path(timing_dir, paste0(package_id, ".tsv"))
fwrite(rbindlist(timing_rows), timing_path, sep = "\t", quote = TRUE)

cat(
  "FG4 rebuild complete: ",
  round(proc.time()[["elapsed"]] - run_started, 3),
  " seconds\n",
  sep = ""
)
