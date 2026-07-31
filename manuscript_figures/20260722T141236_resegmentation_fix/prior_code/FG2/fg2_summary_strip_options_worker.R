#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

DRAFT_ROOT <- file.path(
  "agent-dev", "manuscript_redrafts", "20260709T210405_v7_redraft",
  "figure_generation", "FG2_direct_feature_rebuild", "drafting_v2"
)
INITIAL_DIR <- file.path(DRAFT_ROOT, "initial_subpanels", "summary_strip_options")
REFINED_CONTACT <- file.path(DRAFT_ROOT, "refined_subpanels", "fg2_summary_strip_options_contact_sheet.png")
NOTES_PATH <- file.path(DRAFT_ROOT, "worker_notes", "fg2_summary_strip_options_coverage.md")

dir.create(INITIAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REFINED_CONTACT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(NOTES_PATH), recursive = TRUE, showWarnings = FALSE)

png_dpi <- 420
sum_fuse_line <- "SUM-159-fuse"
scope_no_sum <- "SUM-159-fuse excluded"
scope_all <- "All lines"

feature_palette <- c(
  Growth = "#009E73",
  `Alive AUC` = "#7A5195",
  `Total-cell yield` = "#0072B2"
)

effect_fill <- c(
  negative = "#3B6EA8",
  neutral = "white",
  positive = "#B3423F"
)

feature_meta <- tibble(
  feature = c(
    "growth_lowG_median",
    "growth_highG_median",
    "yield_alive_auc_intercept",
    "yield_alive_auc_slope",
    "peak_total_yield_1mM"
  ),
  display_label = c(
    "Low-G growth",
    "High-G exp. growth",
    "Alive AUC intercept",
    "Alive AUC gradient",
    "1 mM net total yield"
  ),
  feature_group = c(
    "Growth", "Growth",
    "Alive AUC", "Alive AUC",
    "Total-cell yield"
  ),
  feature_order = seq_len(5)
)

g0_feature_meta <- tibble(
  metric = c("max_growth_rate", "live_auc_glucose_window", "peak_total_yield_net"),
  display_label = c("Max growth derivative", "Live-cell AUC", "Net total yield"),
  feature_group = c("Growth", "Alive AUC", "Total-cell yield"),
  metric_order = seq_len(3)
)

read_csv_tbl <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }
  as_tibble(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...))
}

write_csv_tbl <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(as.data.frame(x), path, row.names = FALSE)
  invisible(path)
}

save_png <- function(plot, path, width, height, dpi = png_dpi) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

strip_theme <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey94", color = "grey74", linewidth = 0.25),
      strip.text = element_text(face = "bold", margin = margin(2, 2, 2, 2)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.height = unit(0.12, "in"),
      legend.key.width = unit(0.24, "in"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank()
    )
}

clip_effect <- function(x, limit = 3) {
  pmax(-limit, pmin(limit, x))
}

scale_effects <- function(df, group_cols, value_col = "effect_per_ploidy") {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(
      scale_denom = median(abs(.data[[value_col]]), na.rm = TRUE),
      scale_denom = if_else(!is.finite(scale_denom) | scale_denom <= 0, 1, scale_denom),
      scaled_effect = .data[[value_col]] / scale_denom,
      clipped_scaled_effect = clip_effect(scaled_effect)
    ) %>%
    ungroup()
}

line_order_from <- function(signature_panel) {
  signature_panel %>%
    arrange(line_id, cellLine) %>%
    distinct(cellLine) %>%
    pull(cellLine)
}

add_ploidy_state <- function(df) {
  df %>%
    group_by(cellLine) %>%
    mutate(
      ploidy_state = if_else(ploidy_metric == min(ploidy_metric, na.rm = TRUE), "low", "high")
    ) %>%
    ungroup()
}

build_g0_effects <- function(feature_panel, scope_name) {
  feature_panel %>%
    mutate(peak_total_yield_net = pmax(total_net_gain_to_glucose_end, 0)) %>%
    add_ploidy_state() %>%
    select(cellLine, line_id, G0, ploidy_state, ploidy_metric, all_of(g0_feature_meta$metric)) %>%
    pivot_longer(
      cols = all_of(g0_feature_meta$metric),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(cellLine, line_id, G0, metric, ploidy_state) %>%
    summarise(
      value = median(value, na.rm = TRUE),
      ploidy_metric = median(ploidy_metric, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(value = if_else(is.nan(value), NA_real_, value)) %>%
    pivot_wider(
      names_from = ploidy_state,
      values_from = c(value, ploidy_metric)
    ) %>%
    mutate(
      effect_per_ploidy = (value_high - value_low) / (ploidy_metric_high - ploidy_metric_low),
      scope = scope_name
    ) %>%
    left_join(g0_feature_meta, by = "metric") %>%
    filter(case_when(
      metric == "live_auc_glucose_window" ~ G0 <= 1,
      metric == "peak_total_yield_net" ~ abs(G0 - 1) < 1e-8,
      TRUE ~ TRUE
    )) %>%
    filter(is.finite(effect_per_ploidy))
}

build_fixed_yield_effects <- function(feature_panel, scope_name) {
  feature_panel %>%
    mutate(peak_total_yield_1mM = pmax(total_net_gain_to_glucose_end, 0)) %>%
    filter(abs(G0 - 1) < 1e-8) %>%
    add_ploidy_state() %>%
    group_by(cellLine, line_id, ploidy_state) %>%
    summarise(
      value = median(peak_total_yield_1mM, na.rm = TRUE),
      ploidy_metric = median(ploidy_metric, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(value = if_else(is.nan(value), NA_real_, value)) %>%
    pivot_wider(
      names_from = ploidy_state,
      values_from = c(value, ploidy_metric)
    ) %>%
    mutate(
      feature = "peak_total_yield_1mM",
      effect_per_ploidy = (value_high - value_low) / (ploidy_metric_high - ploidy_metric_low),
      scope = scope_name
    ) %>%
    select(scope, cellLine, feature, effect_per_ploidy) %>%
    filter(is.finite(effect_per_ploidy))
}

corrected_wp4_root <- file.path(
  "agent-dev", "manuscript_redrafts", "20260619_v6_redraft",
  "stage_outputs", "analysis", "recomputed_wp4_sum159_fuse_exclusion"
)
paired_effects_path <- file.path(corrected_wp4_root, "wp4_paired_effects_by_scope.csv")
feature_no_sum_path <- file.path(corrected_wp4_root, "wp4_feature_panel_no_sum159_fuse.csv")
feature_all_path <- file.path(corrected_wp4_root, "wp4_feature_panel_all_lines.csv")
signature_no_sum_path <- file.path(corrected_wp4_root, "wp4_signature_panel_no_sum159_fuse.csv")
signature_all_path <- file.path(corrected_wp4_root, "wp4_signature_panel_all_lines.csv")

feature_no_sum <- read_csv_tbl(feature_no_sum_path)
feature_all <- read_csv_tbl(feature_all_path)
signature_no_sum <- read_csv_tbl(signature_no_sum_path)
signature_all <- read_csv_tbl(signature_all_path)
line_levels <- line_order_from(signature_all)

paired_effects <- bind_rows(
  read_csv_tbl(paired_effects_path) %>%
    filter(feature %in% setdiff(feature_meta$feature, "peak_total_yield_1mM")) %>%
    select(scope, cellLine, feature, effect_per_ploidy),
  build_fixed_yield_effects(feature_all, scope_all),
  build_fixed_yield_effects(feature_no_sum, scope_no_sum)
) %>%
  filter(scope %in% c(scope_no_sum, scope_all)) %>%
  select(-any_of(c("display_label", "category"))) %>%
  left_join(feature_meta, by = "feature") %>%
  mutate(
    scope = factor(scope, levels = c(scope_no_sum, scope_all)),
    display_label = factor(display_label, levels = rev(feature_meta$display_label)),
    feature_group = factor(feature_group, levels = names(feature_palette))
  ) %>%
  scale_effects(group_cols = c("feature"))

effect_consistency <- paired_effects %>%
  group_by(scope, feature, display_label, feature_order) %>%
  summarise(
    n_lines = n(),
    n_positive = sum(effect_per_ploidy > 0, na.rm = TRUE),
    n_negative = sum(effect_per_ploidy < 0, na.rm = TRUE),
    consensus_n = pmax(n_positive, n_negative),
    consensus_direction = case_when(
      n_positive > n_negative ~ "positive",
      n_negative > n_positive ~ "negative",
      TRUE ~ "tie"
    ),
    median_effect_per_ploidy = median(effect_per_ploidy, na.rm = TRUE),
    .groups = "drop"
  )

effect_matrix <- paired_effects %>%
  select(scope, cellLine, feature, effect_per_ploidy) %>%
  left_join(feature_meta, by = "feature") %>%
  mutate(
    scope = factor(scope, levels = c(scope_no_sum, scope_all)),
    display_label = factor(display_label, levels = rev(feature_meta$display_label)),
    feature_group = factor(feature_group, levels = names(feature_palette))
  ) %>%
  scale_effects(group_cols = c("feature"))

no_sum_effects <- paired_effects %>%
  filter(scope == scope_no_sum) %>%
  mutate(cellLine = factor(cellLine, levels = line_levels))

no_sum_summary <- no_sum_effects %>%
  group_by(feature, display_label, feature_group, feature_order) %>%
  summarise(
    median_scaled_effect = median(scaled_effect, na.rm = TRUE),
    mean_scaled_effect = mean(scaled_effect, na.rm = TRUE),
    n_positive = sum(effect_per_ploidy > 0, na.rm = TRUE),
    n_negative = sum(effect_per_ploidy < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    clipped_median_scaled_effect = clip_effect(median_scaled_effect),
    display_label = factor(display_label, levels = rev(feature_meta$display_label))
  )

option_a <- ggplot(no_sum_effects, aes(x = clipped_scaled_effect, y = display_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_segment(
    data = no_sum_summary,
    aes(
      x = 0,
      xend = clipped_median_scaled_effect,
      y = display_label,
      yend = display_label,
      color = feature_group
    ),
    inherit.aes = FALSE,
    linewidth = 2.2,
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = feature_group),
    shape = 21,
    color = "grey20",
    stroke = 0.22,
    size = 1.75,
    position = position_jitter(width = 0, height = 0.055)
  ) +
  facet_grid(feature_group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = feature_palette, guide = "none") +
  scale_fill_manual(values = feature_palette, name = "feature class") +
  scale_x_continuous(
    limits = c(-3.05, 3.05),
    oob = squish,
    breaks = -3:3
  ) +
  labs(
    x = "high-low effect per ploidy, scaled within feature (SUM-159-fuse excluded)",
    y = NULL
  ) +
  strip_theme(base_size = 7.2) +
  theme(
    panel.spacing.y = unit(0.18, "lines"),
    strip.text.y = element_text(angle = 0, size = 6.2),
    axis.text.y = element_text(size = 6.6),
    legend.position = "none"
  )

option_b_df <- effect_matrix %>%
  filter(scope == scope_no_sum) %>%
  mutate(
    cellLine = factor(cellLine, levels = rev(line_levels)),
    display_label = factor(display_label, levels = feature_meta$display_label)
  )

option_b <- ggplot(option_b_df, aes(display_label, cellLine)) +
  geom_tile(aes(fill = clipped_scaled_effect), color = "white", linewidth = 0.36) +
  geom_point(aes(size = pmin(abs(scaled_effect), 3)), shape = 21, color = "grey25", fill = "white", stroke = 0.24, alpha = 0.86) +
  facet_grid(. ~ feature_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = effect_fill[["negative"]],
    mid = effect_fill[["neutral"]],
    high = effect_fill[["positive"]],
    midpoint = 0,
    limits = c(-3, 3),
    oob = squish,
    name = "scaled\nhigh-low"
  ) +
  scale_size_continuous(range = c(0.7, 2.6), guide = "none") +
  labs(
    x = NULL,
    y = "cell line"
  ) +
  strip_theme(base_size = 7.2) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6.6),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

scope_medians <- paired_effects %>%
  group_by(scope, feature, display_label, feature_group, feature_order) %>%
  summarise(
    median_scaled_effect = median(scaled_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    marker = if_else(scope == scope_no_sum, "no-SUM median", "all-line median"),
    marker_value = clip_effect(median_scaled_effect)
  )

all_line_points <- paired_effects %>%
  filter(scope == scope_all) %>%
  mutate(
    marker = if_else(cellLine == sum_fuse_line, "SUM-159-fuse line", "other all-line line"),
    marker_value = clipped_scaled_effect
  )

option_c_df <- bind_rows(
  all_line_points %>%
    select(display_label, feature_group, marker, marker_value, feature_order),
  scope_medians %>%
    select(display_label, feature_group, marker, marker_value, feature_order)
) %>%
  mutate(
    marker = factor(
      marker,
      levels = c("other all-line line", "SUM-159-fuse line", "all-line median", "no-SUM median")
    ),
    display_label = factor(display_label, levels = rev(feature_meta$display_label))
  )

option_c <- ggplot(option_c_df, aes(marker_value, display_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_point(
    aes(shape = marker, fill = marker, color = marker),
    size = 2.05,
    stroke = 0.34,
    position = position_jitter(width = 0, height = 0.045)
  ) +
  facet_grid(feature_group ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(
    values = c(
      "other all-line line" = 21,
      "SUM-159-fuse line" = 23,
      "all-line median" = 22,
      "no-SUM median" = 24
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "other all-line line" = "grey82",
      "SUM-159-fuse line" = "#D55E00",
      "all-line median" = "black",
      "no-SUM median" = "white"
    ),
    name = NULL
  ) +
  scale_color_manual(
    values = c(
      "other all-line line" = "grey35",
      "SUM-159-fuse line" = "black",
      "all-line median" = "black",
      "no-SUM median" = "black"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-3.05, 3.05),
    oob = squish,
    breaks = -3:3
  ) +
  labs(
    x = "high-low effect per ploidy, scaled within feature",
    y = NULL
  ) +
  strip_theme(base_size = 7.2) +
  theme(
    panel.spacing.y = unit(0.18, "lines"),
    strip.text.y = element_text(angle = 0, size = 6.2),
    axis.text.y = element_text(size = 6.6),
    legend.position = "bottom"
  )

g0_effects <- bind_rows(
  build_g0_effects(feature_no_sum, scope_no_sum),
  build_g0_effects(feature_all, scope_all)
) %>%
  mutate(
    scope = factor(scope, levels = c(scope_no_sum, scope_all)),
    G0_display = factor(as.character(G0), levels = c("0", "0.1", "0.25", "0.5", "1", "5", "25")),
    cellLine = factor(cellLine, levels = rev(line_levels)),
    display_label = factor(display_label, levels = g0_feature_meta$display_label)
  ) %>%
  scale_effects(group_cols = c("metric"))

option_d <- ggplot(g0_effects, aes(G0_display, cellLine)) +
  geom_tile(aes(fill = clipped_scaled_effect), color = "white", linewidth = 0.25) +
  facet_grid(scope ~ display_label, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = effect_fill[["negative"]],
    mid = effect_fill[["neutral"]],
    high = effect_fill[["positive"]],
    midpoint = 0,
    limits = c(-3, 3),
    oob = squish,
    name = "scaled\nhigh-low"
  ) +
  labs(
    x = "starting glucose (mM)",
    y = "cell line"
  ) +
  strip_theme(base_size = 6.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8),
    axis.text.y = element_text(size = 6.2),
    strip.text.x = element_text(size = 6.2),
    strip.text.y = element_text(angle = 90, size = 6.0),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

option_paths <- tibble(
  option_id = c("A", "B", "C", "D"),
  target = c("Figure 2A", "Figure 2A", "Figure S3", "Figure 2A and Figure S3"),
  option_summary = c(
    "No-SUM feature-effect bar strip: per-line points with the median high-low effect highlighted for each displayed feature.",
    "No-SUM line-by-feature consistency tile: compact matrix showing which lines/features have positive or negative high-low effects.",
    "All-line/SUM sensitivity strip: all-line points with SUM-159-fuse highlighted plus no-SUM and all-line medians.",
    "Starting-glucose delta strip: high-low effects for growth across G0, live-cell AUC across G0 <= 1, and net total yield at G0 = 1."
  ),
  png_path = file.path(
    INITIAL_DIR,
    c(
      "fg2_summary_strip_option_a_no_sum_effect_bars.png",
      "fg2_summary_strip_option_b_no_sum_line_consistency_tiles.png",
      "fg2_summary_strip_option_c_s3_all_line_scope_comparison.png",
      "fg2_summary_strip_option_d_starting_glucose_delta_tiles.png"
    )
  )
)

save_png(option_a, option_paths$png_path[[1]], width = 7.0, height = 2.25)
save_png(option_b, option_paths$png_path[[2]], width = 7.0, height = 2.45)
save_png(option_c, option_paths$png_path[[3]], width = 7.0, height = 2.55)
save_png(option_d, option_paths$png_path[[4]], width = 7.0, height = 3.35)

contact_sheet <- ((option_a | option_b) / (option_c | option_d)) +
  plot_layout(guides = "collect", heights = c(1, 1.18)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.tag.position = c(0.012, 0.985),
    legend.position = "bottom"
  )

save_png(contact_sheet, REFINED_CONTACT, width = 12.4, height = 8.6, dpi = 360)

effect_value_path <- file.path(INITIAL_DIR, "fg2_summary_strip_effect_values.csv")
g0_value_path <- file.path(INITIAL_DIR, "fg2_summary_strip_starting_glucose_effect_values.csv")
manifest_path <- file.path(INITIAL_DIR, "fg2_summary_strip_option_manifest.csv")

write_csv_tbl(
  paired_effects %>%
    select(
      scope, cellLine, feature, display_label, feature_group,
      effect_per_ploidy, scaled_effect, clipped_scaled_effect
    ),
  effect_value_path
)

write_csv_tbl(
  g0_effects %>%
    select(
      scope, cellLine, line_id, G0, metric, display_label, feature_group,
      effect_per_ploidy, scaled_effect, clipped_scaled_effect
    ),
  g0_value_path
)

write_csv_tbl(
  option_paths %>%
    mutate(
      contact_sheet = REFINED_CONTACT,
      effect_values_csv = effect_value_path,
      starting_glucose_values_csv = g0_value_path
    ),
  manifest_path
)

no_sum_consensus <- effect_consistency %>%
  filter(scope == scope_no_sum) %>%
  arrange(feature_order) %>%
  transmute(
    line = sprintf(
      "- %s: %s consensus (%s/%s lines), median effect per ploidy %.4g.",
      display_label, consensus_direction, consensus_n, n_lines, median_effect_per_ploidy
    )
  ) %>%
  pull(line)

coverage_lines <- c(
  "# FG2 summary-strip options coverage",
  "",
  "## Scope",
  "",
  "This worker addresses the current request to redraft the Figure 2A-compatible summary strip so yield is net peak total-cell gain at 1 mM rather than a yield-per-glucose regression feature.",
  "",
  "## Direct feedback coverage",
  "",
  "| directive | status | output |",
  "|---|---|---|",
  "| Replace yield-per-glucose intercept/gradient with net peak total-cell yield at 1 mM. | addressed | Option A/B/C and effect audit CSV now use `peak_total_yield_1mM`, computed as `max(total_peak_to_glucose_end - total_initial_at_glucose_start, 0)` at `G0 == 1`. |",
  "| Keep live-AUC regression features restricted to G0 <= 1. | addressed | Option A/B/C retain the exported canonical live-AUC intercept/gradient features; option D only shows live-AUC starting-glucose effects through G0 <= 1. |",
  "| Draft a few alternatives for reviewer choice rather than picking one final solution. | addressed | Four alternatives were generated: effect bars, line-feature tiles, all-line scope comparison, and starting-glucose delta tiles. |",
  "| Use existing data/exports only; no major analysis or refit. | addressed | Inputs were corrected redraft-local WP4 CSV exports under `stage_outputs/analysis/recomputed_wp4_sum159_fuse_exclusion`. No refit, segmentation, classifier, or SLURM job was run. |",
  "| Regenerate local review surfaces consistently. | addressed | This worker writes summary-strip subpanels/CSVs for downstream final assembly and report regeneration. |",
  "",
  "## Outputs",
  "",
  paste0("- Option ", option_paths$option_id, " (", option_paths$target, "): `", option_paths$png_path, "` - ", option_paths$option_summary),
  paste0("- Contact sheet: `", REFINED_CONTACT, "`."),
  paste0("- Option manifest: `", manifest_path, "`."),
  paste0("- Effect value audit CSV: `", effect_value_path, "`."),
  paste0("- Starting-glucose value audit CSV: `", g0_value_path, "`."),
  "",
  "## Data sources",
  "",
  paste0("- `", paired_effects_path, "`: existing high-low paired effects by scope for growth and live-AUC features."),
  paste0("- `", feature_no_sum_path, "`: existing no-SUM per-condition feature panel."),
  paste0("- `", feature_all_path, "`: existing all-line per-condition feature panel."),
  paste0("- `", signature_no_sum_path, "` and `", signature_all_path, "`: existing line order/signature context."),
  "- Fixed-yield effects are computed locally from the feature panels at G0 = 1 as `pmax(total_net_gain_to_glucose_end, 0)`.",
  "",
  "## No-SUM effect context",
  "",
  no_sum_consensus,
  "",
  "## Visual QC",
  "",
  "- PNGs are option strips/contact material, not replacement final figures.",
  "- Images avoid figure-level manuscript titles/captions and use only axis, facet, legend, and panel-tag labeling.",
  "- All option strips omit death features.",
  "- Option A gives the fastest Figure 2A readout of which retained features have nonzero no-SUM effects.",
  "- Option B emphasizes cross-line consistency for Figure 2A.",
  "- Option C is targeted to Figure S3 and keeps SUM-159-fuse visually explicit rather than hidden in the all-line context.",
  "- Option D is the only option that explicitly summarizes effects across starting glucose; it is useful for explaining aggregation but is denser than A-C.",
  "",
  "## Caveats and blockers",
  "",
  "- Values are existing-export summaries, with scaled effects used only for display across heterogeneous feature units. The underlying unscaled values are written to the audit CSVs.",
  "- Starting-glucose tiles use high-low paired condition summaries from existing feature-panel exports; they are not a new spline-fitting or feature-engineering pass.",
  "- No blockers encountered."
)

writeLines(coverage_lines, NOTES_PATH)

message("Wrote summary-strip options to: ", INITIAL_DIR)
message("Wrote contact sheet to: ", REFINED_CONTACT)
message("Wrote coverage notes to: ", NOTES_PATH)
