args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp1_core_starvation")
}

figure_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp1_core_starvation")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/model_free_ploidy_utils.R")

sum_fuse_line <- "SUM-159-fuse"
representative_g0 <- c(0, 0.25, 1, 25)
sign_tol <- 1e-8

g0_palette <- c(
  "0" = "#000000",
  "0.1" = "#E69F00",
  "0.25" = "#56B4E9",
  "0.5" = "#009E73",
  "1" = "#F0E442",
  "5" = "#0072B2",
  "25" = "#D55E00"
)

wp1_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey72"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

save_plot_pair <- function(plot, basename, width, height, dpi = 450) {
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))

  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)

  invisible(c(pdf = pdf_path, png = png_path))
}

format_ploidy_abs <- function(x) {
  out <- ifelse(abs(x - round(x)) < 1e-8, sprintf("%.0fN", x), sprintf("%.1fN", x))
  out[!is.finite(x)] <- NA_character_
  out
}

wrap_label <- function(x, width = 22) {
  vapply(
    x,
    function(one) paste(strwrap(one, width = width), collapse = "\n"),
    character(1)
  )
}

ordered_g0_factor <- function(x) {
  factor(as.character(x), levels = names(g0_palette))
}

add_display_columns <- function(df) {
  df %>%
    mutate(
      G0_display = ordered_g0_factor(G0),
      ploidy_abs_label = format_ploidy_abs(ploidy_abs),
      ploidy_state = if_else(ploidy_metric == min(ploidy_metric, na.rm = TRUE), "low", "high"),
      ploidy_display = paste0(
        ploidy_state,
        ": ",
        ploidy_abs_label
      )
    )
}

make_trajectory_long <- function(count_summary, glucose_summary) {
  count_live <- count_summary %>%
    transmute(
      well_idx, cellLine, line_id, ploidy_metric, ploidy_abs, G0, hours,
      measurement = "Live cells",
      value = live_cells
    )

  count_dead <- count_summary %>%
    transmute(
      well_idx, cellLine, line_id, ploidy_metric, ploidy_abs, G0, hours,
      measurement = "Dead cells",
      value = dead_cells
    )

  glucose <- glucose_summary %>%
    transmute(
      well_idx, cellLine, line_id, ploidy_metric, ploidy_abs, G0, hours,
      measurement = "Glucose",
      value = glucose_hat
    )

  out <- bind_rows(count_live, count_dead, glucose) %>%
    mutate(
      measurement = factor(measurement, levels = c("Live cells", "Dead cells", "Glucose"))
    ) %>%
    group_by(cellLine) %>%
    add_display_columns() %>%
    ungroup()

  ploidy_levels <- out %>%
    distinct(ploidy_state, ploidy_abs, ploidy_display) %>%
    arrange(factor(ploidy_state, levels = c("low", "high")), ploidy_abs) %>%
    pull(ploidy_display)

  out %>%
    mutate(
      ploidy_state = factor(ploidy_state, levels = c("low", "high")),
      ploidy_display = factor(ploidy_display, levels = unique(ploidy_levels))
    )
}

make_experiment_schematic <- function(feature_panel, count_summary, glucose_summary) {
  design <- tibble(
    step = factor(
      c("Cell lines", "Ploidy pairs", "Glucose starts", "Readouts"),
      levels = c("Cell lines", "Ploidy pairs", "Glucose starts", "Readouts")
    ),
    xpos = seq_len(4),
    detail = c(
      sprintf("%d engineered\ncomparisons", n_distinct(feature_panel$cellLine)),
      "paired low/high\nploidy states",
      paste(g0_display_levels(), collapse = ", "),
      sprintf(
        "live/dead: %s h\n glucose: %s h",
        paste(sort(unique(diff(sort(unique(count_summary$hours))))), collapse = "/"),
        paste(sort(unique(diff(sort(unique(glucose_summary$hours))))), collapse = "/")
      )
    )
  )

  arrows <- tibble(x = c(1.43, 2.43, 3.43), xend = c(1.57, 2.57, 3.57), y = 1, yend = 1)

  ggplot(design, aes(xpos, 1)) +
    geom_segment(
      data = arrows,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.45,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.08, "in"))
    ) +
    geom_label(
      aes(label = paste(step, detail, sep = "\n")),
      label.size = 0.25,
      label.padding = grid::unit(0.12, "lines"),
      fill = "white",
      color = "grey15",
      size = 2.5,
      lineheight = 0.92
    ) +
    coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.65, 1.35), clip = "off") +
    labs(title = "A. Experimental design") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_trajectory_plot <- function(trajectory_long, g0_keep = g0_display_levels(), title = NULL) {
  trajectory_long %>%
    filter(G0 %in% g0_keep) %>%
    mutate(G0_display = ordered_g0_factor(G0)) %>%
    ggplot(aes(hours, value, color = G0_display, group = interaction(well_idx, measurement))) +
    geom_line(linewidth = 0.35, alpha = 0.85) +
    facet_grid(measurement ~ cellLine + ploidy_display, scales = "free_y") +
    scale_color_manual(values = g0_palette, drop = TRUE, name = "starting glucose") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    labs(title = title, x = "hours", y = NULL) +
    wp1_theme(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 0),
      strip.text.x = element_text(size = 5.8),
      strip.text.y = element_text(size = 6.5),
      panel.spacing.x = grid::unit(0.45, "lines"),
      panel.spacing.y = grid::unit(0.35, "lines")
    )
}

make_line_trajectory_plot <- function(trajectory_long, line_name) {
  trajectory_long %>%
    filter(cellLine == line_name) %>%
    mutate(G0_display = ordered_g0_factor(G0)) %>%
    ggplot(aes(hours, value, color = G0_display, group = interaction(well_idx, measurement))) +
    geom_line(linewidth = 0.55, alpha = 0.9) +
    facet_grid(measurement ~ ploidy_display, scales = "free_y") +
    scale_color_manual(values = g0_palette, drop = TRUE, name = "starting glucose") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    labs(
      title = paste("Supplemental trajectory grid:", line_name),
      x = "hours",
      y = NULL
    ) +
    wp1_theme(base_size = 9) +
    theme(panel.spacing = grid::unit(0.6, "lines"))
}

feature_catalog_for_wp1 <- function() {
  get_model_free_feature_catalog() %>%
    mutate(
      display_label = case_when(
        feature == "growth_lowG_median" ~ "Low-G growth",
        feature == "growth_highG_median" ~ "High-G exponential growth",
        feature == "death_lowG_median" ~ "Low-G death",
        feature == "death_highG_median" ~ "High-G death",
        feature == "yield_alive_auc_intercept" ~ "Alive AUC baseline",
        feature == "yield_alive_auc_slope" ~ "Alive AUC glucose response",
        feature == "peak_total_yield_intercept" ~ "Glucose-yield baseline",
        feature == "peak_total_yield_slope" ~ "Glucose-yield response",
        TRUE ~ short_label
      ),
      display_label = factor(display_label, levels = display_label),
      category = factor(category, levels = c("Growth", "Death", "Yield"))
    )
}

build_effect_scope_tables <- function(signature_panel) {
  catalog <- feature_catalog_for_wp1()

  paired_effects <- compute_empirical_effects(signature_panel)$paired_effects_long %>%
    left_join(catalog, by = "feature") %>%
    mutate(
      effect_direction = case_when(
        effect_per_ploidy > sign_tol ~ "higher in high ploidy",
        effect_per_ploidy < -sign_tol ~ "lower in high ploidy",
        TRUE ~ "near zero"
      )
    )

  effect_scale <- paired_effects %>%
    group_by(feature) %>%
    summarise(
      center = mean(effect_per_ploidy, na.rm = TRUE),
      spread = sd(effect_per_ploidy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(spread = if_else(is.finite(spread) & spread > 0, spread, 1))

  scoped_effects <- bind_rows(
    paired_effects %>% mutate(scope = "All lines"),
    paired_effects %>% filter(cellLine != sum_fuse_line) %>% mutate(scope = "SUM-159-fuse excluded")
  ) %>%
    left_join(effect_scale, by = "feature") %>%
    mutate(
      scope = factor(scope, levels = c("All lines", "SUM-159-fuse excluded")),
      effect_z = (effect_per_ploidy - center) / spread
    )

  consistency <- scoped_effects %>%
    group_by(scope, feature, display_label, category) %>%
    summarise(
      n_lines = sum(is.finite(effect_per_ploidy)),
      n_positive = sum(effect_per_ploidy > sign_tol, na.rm = TRUE),
      n_negative = sum(effect_per_ploidy < -sign_tol, na.rm = TRUE),
      n_near_zero = sum(abs(effect_per_ploidy) <= sign_tol | !is.finite(effect_per_ploidy), na.rm = TRUE),
      consensus_n = max(n_positive, n_negative),
      consensus_fraction = if_else(n_lines > 0, consensus_n / n_lines, NA_real_),
      consensus_direction = case_when(
        n_positive > n_negative ~ "positive",
        n_negative > n_positive ~ "negative",
        TRUE ~ "tie"
      ),
      .groups = "drop"
    ) %>%
    arrange(scope, category, display_label)

  pass_table <- consistency %>%
    mutate(
      scope_key = recode(
        as.character(scope),
        "All lines" = "all_lines",
        "SUM-159-fuse excluded" = "no_sum159_fuse"
      )
    ) %>%
    select(scope_key, feature, display_label, n_lines, consensus_n, consensus_fraction, consensus_direction) %>%
    tidyr::pivot_wider(
      names_from = scope_key,
      values_from = c(n_lines, consensus_n, consensus_fraction, consensus_direction),
      names_glue = "{.value}__{scope_key}"
    ) %>%
    mutate(
      wp1_sign_pass = consensus_n__all_lines >= 4 &
        consensus_n__no_sum159_fuse >= 4 &
        consensus_fraction__no_sum159_fuse >= consensus_fraction__all_lines
    )

  list(
    paired_effects = paired_effects,
    scoped_effects = scoped_effects,
    consistency = consistency,
    pass_table = pass_table
  )
}

make_effect_heatmap <- function(scoped_effects) {
  scoped_effects %>%
    mutate(
      display_label_wrapped = wrap_label(as.character(display_label), width = 16),
      display_label_wrapped = factor(display_label_wrapped, levels = unique(display_label_wrapped))
    ) %>%
    ggplot(aes(display_label_wrapped, cellLine, fill = effect_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_point(
      aes(size = abs(effect_per_ploidy), shape = effect_direction),
      color = "grey12",
      fill = "white",
      stroke = 0.35,
      alpha = 0.85
    ) +
    facet_wrap(~scope, nrow = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "within-feature\nz-score"
    ) +
    scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
    scale_shape_manual(
      values = c("higher in high ploidy" = 24, "lower in high ploidy" = 25, "near zero" = 21),
      name = "raw effect sign"
    ) +
    labs(
      title = "A. Paired high-minus-low ploidy contrasts",
      x = NULL,
      y = NULL
    ) +
    wp1_theme(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
}

make_consistency_plot <- function(consistency) {
  consistency %>%
    mutate(
      display_label_wrapped = wrap_label(as.character(display_label), width = 16),
      display_label_wrapped = factor(display_label_wrapped, levels = unique(display_label_wrapped)),
      consensus_label = paste0(consensus_n, "/", n_lines)
    ) %>%
    ggplot(aes(display_label_wrapped, consensus_fraction, fill = consensus_direction)) +
    geom_col(width = 0.72, color = "grey20", linewidth = 0.18) +
    geom_text(aes(label = consensus_label), vjust = -0.25, size = 2.4) +
    facet_wrap(~scope, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.03))) +
    scale_fill_manual(
      values = c(positive = "#B2182B", negative = "#2166AC", tie = "grey70"),
      name = "consensus\nsign"
    ) +
    labs(
      title = "B. Sign consistency across paired cell-line contrasts",
      x = NULL,
      y = "same-sign line fraction"
    ) +
    wp1_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
}

write_wp1_text_report <- function(
  path,
  stan_data_path,
  figure_dir,
  output_dir,
  feature_panel,
  signature_panel,
  consistency,
  pass_table
) {
  passing_features <- pass_table %>%
    filter(wp1_sign_pass) %>%
    arrange(display_label)

  line_summary <- feature_panel %>%
    group_by(cellLine) %>%
    summarise(
      n_ploidy_states = n_distinct(ploidy_metric),
      n_glucose_conditions = n_distinct(G0),
      min_hours = min(c(0, live_peak_time), na.rm = TRUE),
      max_hours = max(c(glucose_end_time, live_peak_time, total_peak_time), na.rm = TRUE),
      .groups = "drop"
    )

  lines <- c(
    "WP1_CORE_STARVATION_PHENOTYPE",
    sprintf("generated\t%s", Sys.time()),
    sprintf("stan_data\t%s", stan_data_path),
    sprintf("figure_dir\t%s", figure_dir),
    sprintf("table_dir\t%s", output_dir),
    "",
    "SECTION\tPASS_CRITERION_SUMMARY",
    sprintf("n_features_passing_sign_gate\t%d", nrow(passing_features)),
    if (nrow(passing_features)) {
      capture.output(write.table(
        as.data.frame(passing_features %>% select(feature, display_label, starts_with("consensus_n"), starts_with("consensus_fraction"), starts_with("consensus_direction"))),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      ))
    } else {
      "No curated phenotype met the >=4/5 all-line and 4/4 no-SUM-159-fuse sign-consistency gate."
    },
    "",
    "SECTION\tDESIGN_SUMMARY_BY_LINE",
    capture.output(write.table(as.data.frame(line_summary), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tCONSISTENCY_BY_SCOPE",
    capture.output(write.table(
      as.data.frame(consistency %>% arrange(scope, category, display_label)),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "SECTION\tSIGNATURE_PANEL",
    capture.output(write.table(
      as.data.frame(signature_panel %>% arrange(cellLine, ploidy_metric)),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    ))
  )

  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

stan_data_path <- resolve_stan_data_path()
message("Building WP1 model-free tables from ", stan_data_path)
tables <- build_feature_panel(stan_data_path = stan_data_path)

count_summary <- tables$count_summary
glucose_summary <- tables$glucose_summary
feature_panel <- tables$feature_panel
signature_panel <- summarize_glucose_signatures(feature_panel)
trajectory_long <- make_trajectory_long(count_summary, glucose_summary)
effect_tables <- build_effect_scope_tables(signature_panel)

write.csv(count_summary, file.path(output_dir, "count_trajectory_summary.csv"), row.names = FALSE)
write.csv(glucose_summary, file.path(output_dir, "glucose_trajectory_summary.csv"), row.names = FALSE)
write.csv(feature_panel, file.path(output_dir, "wp1_feature_panel.csv"), row.names = FALSE)
write.csv(signature_panel, file.path(output_dir, "wp1_signature_panel.csv"), row.names = FALSE)
write.csv(effect_tables$paired_effects, file.path(output_dir, "wp1_paired_effects.csv"), row.names = FALSE)
write.csv(effect_tables$consistency, file.path(output_dir, "wp1_effect_consistency.csv"), row.names = FALSE)
write.csv(effect_tables$pass_table, file.path(output_dir, "wp1_pass_gate.csv"), row.names = FALSE)

saveRDS(count_summary, file.path(output_dir, "count_trajectory_summary.Rds"))
saveRDS(glucose_summary, file.path(output_dir, "glucose_trajectory_summary.Rds"))
saveRDS(feature_panel, file.path(output_dir, "wp1_feature_panel.Rds"))
saveRDS(signature_panel, file.path(output_dir, "wp1_signature_panel.Rds"))
saveRDS(effect_tables$paired_effects, file.path(output_dir, "wp1_paired_effects.Rds"))
saveRDS(effect_tables$consistency, file.path(output_dir, "wp1_effect_consistency.Rds"))
saveRDS(effect_tables$pass_table, file.path(output_dir, "wp1_pass_gate.Rds"))

schematic_plot <- make_experiment_schematic(feature_panel, count_summary, glucose_summary)
rep_trajectory_plot <- make_trajectory_plot(
  trajectory_long,
  g0_keep = representative_g0,
  title = "B. Representative live, dead, and glucose trajectories"
)

figure_1 <- schematic_plot / rep_trajectory_plot +
  plot_layout(heights = c(0.42, 2.7)) +
  plot_annotation(
    title = "Figure 1. Core Glucose-Starvation Phenotype",
    caption = "Trajectory panels show median object-count summaries by well and calibrated glucose estimates. Representative G0 values are 0, 0.25, 1, and 25."
  )

save_plot_pair(figure_1, "figure_1_core_starvation_trajectories", width = 15.5, height = 8.6)

figure_2 <- make_effect_heatmap(effect_tables$scoped_effects) /
  make_consistency_plot(effect_tables$consistency) +
  plot_layout(heights = c(1.7, 1)) +
  plot_annotation(
    title = "Figure 2. Paired Ploidy Effects With And Without SUM-159-fuse",
    caption = "Effects are high-ploidy minus low-ploidy contrasts per ploidy-metric unit. Tile color is feature-wise standardized for visualization; raw values are exported in wp1_paired_effects.csv."
  )

save_plot_pair(figure_2, "figure_2_paired_ploidy_effects", width = 11.8, height = 8.6)

supplement_pdf <- file.path(figure_dir, "supplement_per_line_trajectory_grids.pdf")
line_plots <- lapply(sort(unique(trajectory_long$cellLine)), function(line_name) {
  p <- make_line_trajectory_plot(trajectory_long, line_name)
  save_plot_pair(
    p,
    paste0("supplement_trajectory_grid_", gsub("[^A-Za-z0-9]+", "_", line_name)),
    width = 8.4,
    height = 6.8
  )
  p
})

grDevices::pdf(supplement_pdf, width = 8.4, height = 6.8, onefile = TRUE)
for (p in line_plots) {
  print(p)
}
grDevices::dev.off()

write_wp1_text_report(
  path = file.path(output_dir, "wp1_summary.txt"),
  stan_data_path = stan_data_path,
  figure_dir = figure_dir,
  output_dir = output_dir,
  feature_panel = feature_panel,
  signature_panel = signature_panel,
  consistency = effect_tables$consistency,
  pass_table = effect_tables$pass_table
)

cat(sprintf("Wrote WP1 figures to %s\n", figure_dir))
cat(sprintf("Wrote WP1 audit tables to %s\n", output_dir))
