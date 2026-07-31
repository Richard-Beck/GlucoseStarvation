#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(ggplot2)
  library(qs)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
phase_i <- match("--phase", args)
phase <- if (!is.na(phase_i) && phase_i < length(args)) {
  args[[phase_i + 1L]]
} else {
  sub("^--phase=", "", args[grepl("^--phase=", args)][1L])
}
if (is.na(phase) || !nzchar(phase)) phase <- "subpanels"
if (!phase %in% c("subpanels", "final")) {
  stop("--phase must be subpanels or final", call. = FALSE)
}

project_root <- normalizePath(".", mustWork = TRUE)
polish_root <- file.path(
  project_root, "manuscript_figures",
  "20260728_posterior_strategy_figure_set",
  "S_strategy_context_support", "polishing"
)
subpanel_dir <- file.path(polish_root, "subpanels")
layout_dir <- file.path(polish_root, "layout")
final_dir <- file.path(polish_root, "final_images")
derived_dir <- file.path(polish_root, "derived")
for (path in c(subpanel_dir, layout_dir, final_dir, derived_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

simulation_root <- file.path(
  project_root, "data", "modelling", "gpath_v1",
  "red_a30_counts_20260722", "derived", "posterior",
  "strategy_simulations", "runs",
  "seed3density_7g_xxx_xcc_axax_20260727"
)
task_path <- file.path(simulation_root, "plans", "tasks.tsv")
schedule_path <- file.path(simulation_root, "schedule_index.tsv")
validation_path <- file.path(simulation_root, "validation.json")
cache_path <- file.path(derived_dir, "support_median_summary.qs")

model_labels <- c(
  "1R_1P_0W_C0_M1" = "1R1P",
  "1R_1P_1W_C0_M0" = "1R1PW",
  "2R_1P_0W_C0_M1" = "2R1P",
  "2R_2P_0W_C0_M1" = "2R2P",
  "2R_2P_1W_C0_M0" = "2R2PW"
)
line_labels <- c(
  "MCF10A" = "MCF10A",
  "MDA-MB-231" = "MDA231",
  "SNU668" = "SNU668",
  "SUM-159-chem" = "SUMchem",
  "SUM-159-fuse" = "SUMfuse"
)
line_order <- unname(line_labels)
context_order <- c("Four-line", "All-line", "Single-line")
day_map <- c(day4_pre_action = "Day 4", day6_end = "Day 6")
omitted_strategy <- "0,A0,A0"

context_label <- function(dataset_id) {
  fifelse(
    dataset_id == "loo_exclude_sum_159_fuse", "Four-line",
    fifelse(dataset_id == "all_lines", "All-line", "Single-line")
  )
}

compact_strategy_label <- function(x) {
  first <- sub(",.*$", "", x)
  family <- fifelse(
    grepl(",C,C$", x), "CC",
    fifelse(grepl(",A[^,]+,A[^,]+$", x), "AA", "XX")
  )
  paste0(first, family)
}

build_summary <- function() {
  validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
  if (!isTRUE(validation$validation_passed)) {
    stop("The compact posterior strategy run is not validated.", call. = FALSE)
  }
  tasks <- fread(task_path)
  schedules <- fread(schedule_path)
  if (
    nrow(tasks) != 70L ||
    nrow(schedules) != 21L ||
    !all(file.exists(file.path(project_root, tasks$output_path)))
  ) {
    stop("Unexpected compact simulation task coverage.", call. = FALSE)
  }
  strategy_order <- schedules[
    strategy_code != omitted_strategy,
    strategy_code
  ]
  parts <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i]
    sim <- qread(
      file.path(project_root, task$output_path),
      nthreads = 2L
    )
    density_i <- match(1, as.numeric(sim$density_multipliers))
    occasion_i <- match(names(day_map), sim$observation_occasions)
    strategy_i <- match(strategy_order, sim$schedule_index$strategy_code)
    low_i <- match("N_low", sim$state_names)
    high_i <- match("N_high", sim$state_names)
    if (anyNA(c(density_i, occasion_i, strategy_i, low_i, high_i))) {
      stop("Required state or index missing from ", task$output_path,
           call. = FALSE)
    }
    low <- sim$observation_state[
      , density_i, strategy_i, occasion_i, low_i, drop = FALSE
    ]
    high <- sim$observation_state[
      , density_i, strategy_i, occasion_i, high_i, drop = FALSE
    ]
    dim(low) <- dim(high) <- c(
      dim(sim$observation_state)[[1L]],
      length(strategy_i),
      length(occasion_i)
    )
    grid <- as.data.table(expand.grid(
      draw = seq_len(dim(low)[[1L]]),
      strategy = strategy_order,
      day = unname(day_map),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
    grid[, `:=`(
      log_high_low = log(
        pmax(as.vector(high), 1e-12) /
          pmax(as.vector(low), 1e-12)
      ),
      context = context_label(task$dataset_id),
      dataset_id = task$dataset_id,
      model = unname(model_labels[[task$model_id]]),
      line = unname(line_labels[[task$line_name]])
    )]
    parts[[i]] <- grid
  }
  draws <- rbindlist(parts)
  summary <- draws[, .(
    support_high = mean(log_high_low > 0),
    median_log_ratio = median(log_high_low)
  ), by = .(
    context, dataset_id, model, line, strategy, day
  )]
  effect_limit <- as.numeric(quantile(
    abs(summary$median_log_ratio), 0.99, na.rm = TRUE
  ))
  bundle <- list(
    summary = summary,
    effect_limit = effect_limit,
    strategy_order = strategy_order,
    metadata = list(
      density_multiplier = 1,
      days = unname(day_map),
      omitted_strategy = omitted_strategy,
      posterior_draws = 100L,
      contexts = context_order
    )
  )
  qsave(bundle, cache_path, preset = "high", nthreads = 2L)
  bundle
}

theme_heatmap <- function(show_y, show_x, show_legend) {
  theme_minimal(base_size = 6, base_family = "sans") +
    theme(
      text = element_text(colour = "#202020"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.28, "mm"),
      strip.background = element_rect(
        fill = "#ECEDE7", colour = "#B4B4AE", linewidth = 0.20
      ),
      strip.text.x = element_text(face = "bold", size = 5.1),
      strip.text.y = element_text(face = "bold", size = 5.0),
      axis.title = element_blank(),
      axis.text.x = if (show_x) {
        element_text(
          angle = 90, hjust = 1, vjust = 0.5,
          size = 4.0, colour = "#252525"
        )
      } else {
        element_blank()
      },
      axis.text.y = if (show_y) {
        element_text(size = 4.4, colour = "#252525")
      } else {
        element_blank()
      },
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.direction = "vertical",
      legend.title = element_text(size = 5.4, face = "bold"),
      legend.text = element_text(size = 4.9),
      legend.key.height = unit(17, "mm"),
      legend.key.width = unit(2.3, "mm"),
      legend.spacing.x = unit(0.7, "mm"),
      plot.margin = margin(0.8, 1.2, 0.8, 1.2)
    )
}

wrap_header <- function(plot, header) {
  ggdraw() +
    draw_grob(
      grid::rectGrob(
        gp = grid::gpar(fill = "#ECEDE7", col = NA)
      ),
      x = 0, y = 0.965, width = 1, height = 0.035
    ) +
    draw_plot(plot, x = 0, y = 0, width = 1, height = 0.965) +
    draw_label(
      header, x = 0.5, y = 0.982,
      hjust = 0.5, vjust = 0.5,
      size = 5.9, fontface = "bold", fontfamily = "sans"
    )
}

make_panel <- function(
  bundle, context, metric,
  show_y, show_x, show_legend
) {
  current_context <- context
  data <- copy(bundle$summary[context == current_context])
  context_lines <- if (context == "Four-line") {
    setdiff(line_order, "SUMfuse")
  } else {
    line_order
  }
  data[, `:=`(
    line = factor(line, levels = context_lines),
    model = factor(model, levels = unname(model_labels)),
    strategy = factor(
      strategy, levels = rev(bundle$strategy_order)
    ),
    day = factor(day, levels = unname(day_map))
  )]
  fill_column <- if (metric == "support") {
    "support_high"
  } else {
    "median_log_ratio"
  }
  plot <- ggplot(
    data,
    aes(x = model, y = strategy, fill = .data[[fill_column]])
  ) +
    geom_tile(colour = "#F4F4F1", linewidth = 0.05) +
    facet_grid(day ~ line, drop = FALSE) +
    scale_y_discrete(labels = compact_strategy_label) +
    theme_heatmap(show_y, show_x, show_legend)
  if (metric == "support") {
    plot <- plot + scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0.5, limits = c(0, 1),
      breaks = c(0, 1),
      name = "Pr(H/L > 1)"
    ) +
      guides(fill = guide_colourbar(
        direction = "vertical", title.position = "top",
        barheight = unit(17, "mm"), barwidth = unit(2.3, "mm")
      ))
    header <- paste(context, "support", sep = " · ")
  } else {
    plot <- plot + scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-bundle$effect_limit, bundle$effect_limit),
      oob = squish,
      breaks = c(-bundle$effect_limit, bundle$effect_limit),
      labels = label_number(accuracy = 0.1),
      name = "Median\nlog(H/L)"
    ) +
      guides(fill = guide_colourbar(
        direction = "vertical", title.position = "top",
        barheight = unit(17, "mm"), barwidth = unit(2.3, "mm")
      ))
    header <- paste(context, "median", sep = " · ")
  }
  wrap_header(plot, header)
}

make_panels <- function(bundle) {
  list(
    a = make_panel(
      bundle, "Four-line", "support",
      show_y = TRUE, show_x = FALSE, show_legend = TRUE
    ),
    b = make_panel(
      bundle, "Four-line", "median",
      show_y = FALSE, show_x = FALSE, show_legend = TRUE
    ),
    c = make_panel(
      bundle, "All-line", "support",
      show_y = TRUE, show_x = FALSE, show_legend = FALSE
    ),
    d = make_panel(
      bundle, "All-line", "median",
      show_y = FALSE, show_x = FALSE, show_legend = FALSE
    ),
    e = make_panel(
      bundle, "Single-line", "support",
      show_y = TRUE, show_x = TRUE, show_legend = FALSE
    ),
    f = make_panel(
      bundle, "Single-line", "median",
      show_y = FALSE, show_x = TRUE, show_legend = FALSE
    )
  )
}

write_subpanels <- function(panels) {
  relative_dir <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "S_strategy_context_support", "polishing", "subpanels"
  )
  dimensions <- data.table(
    figure = "Posterior strategy context supplement",
    panel = letters[1:6],
    subpanel_png = file.path(
      relative_dir, paste0("panel_", letters[1:6], ".png")
    ),
    width_in = rep(3.47, 6),
    height_in = rep(3.043, 6),
    order = 1:6
  )
  dpi <- 600
  for (i in seq_len(nrow(dimensions))) {
    ggsave(
      file.path(project_root, dimensions$subpanel_png[[i]]),
      panels[[dimensions$panel[[i]]]],
      width = dimensions$width_in[[i]],
      height = dimensions$height_in[[i]],
      dpi = dpi,
      bg = "white"
    )
  }
  dimensions[, `:=`(
    width_px = as.integer(round(width_in * dpi)),
    height_px = as.integer(floor(height_in * dpi + 1e-7))
  )]
  setcolorder(dimensions, c(
    "figure", "panel", "subpanel_png",
    "width_px", "height_px", "width_in", "height_in", "order"
  ))
  fwrite(dimensions, file.path(layout_dir, "subpanel_dimensions.csv"))
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

write_metadata <- function(final_path) {
  relative_final <- sub(
    paste0("^", project_root, "/"), "", normalizePath(final_path)
  )
  script_path <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "S_strategy_context_support", "polishing", "scripts",
    "polish_figures.R"
  )
  layout_path <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "S_strategy_context_support", "polishing", "layout",
    "layout_plan.csv"
  )
  subpanels <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "S_strategy_context_support", "polishing", "subpanels",
    paste0("panel_", letters[1:6], ".png")
  )
  data_inputs <- c(
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/posterior/strategy_simulations/runs/seed3density_7g_xxx_xcc_axax_20260727/validation.json",
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/posterior/strategy_simulations/runs/seed3density_7g_xxx_xcc_axax_20260727/plans/tasks.tsv",
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/posterior/strategy_simulations/runs/seed3density_7g_xxx_xcc_axax_20260727/schedule_index.tsv"
  )
  fwrite(
    data.table(
      figure = "Posterior strategy context supplement",
      source_local_name = relative_final,
      source_root = dirname(relative_final),
      panel_labels = paste(letters[1:6], collapse = ";"),
      polished_final_image = relative_final,
      polished_image_sha256 = sha256_file(final_path),
      data_inputs = paste(data_inputs, collapse = ";"),
      subpanel_audit_paths = paste(subpanels, collapse = ";"),
      polishing_script = script_path,
      review_status = "visual_qc_complete"
    ),
    file.path(polish_root, "manifest.csv")
  )
  fwrite(
    data.table(
      figure = "Posterior strategy context supplement",
      panel = letters[1:6],
      subpanel_image = subpanels,
      generator = script_path,
      command = paste(
        "scripts/agentRrunner.sh", script_path, "--phase final"
      ),
      data_inputs = paste(data_inputs, collapse = ";"),
      layout_plan = layout_path,
      output_image = relative_final,
      notes = c(
        "Four-line support at days 4 and 6.",
        "Four-line median log(H/L) at days 4 and 6.",
        "All-line support at days 4 and 6.",
        "All-line median log(H/L) at days 4 and 6.",
        "Unified single-line support at days 4 and 6.",
        "Unified single-line median log(H/L) at days 4 and 6."
      )
    ),
    file.path(polish_root, "provenance.csv")
  )
  fwrite(
    data.table(
      figure = "Posterior strategy context supplement",
      stage = "polishing",
      output_image = relative_final,
      output_sha256 = sha256_file(final_path),
      generator = script_path,
      command = paste(
        "scripts/agentRrunner.sh", script_path, "--phase final"
      ),
      data_inputs = paste(data_inputs, collapse = ";"),
      layout_plan = layout_path,
      status = "visual_qc_complete"
    ),
    file.path(polish_root, "figure_rebuild_manifest.tsv"),
    sep = "\t"
  )
  writeLines(
    c(
      "# Supplement polishing notes",
      "",
      "- Regenerated directly from the validated 70-task compact simulation run; no tmp raster or tmp table is embedded.",
      "- Data slice: 1x density, days 4 and 6, 100 posterior draws, five Pareto models, and 20 strategies after omitting 0AA.",
      "- Support is Pr(H/L > 1); median is posterior median log(H/L).",
      "- The median color limit is one shared symmetric 99th-percentile limit across all contexts.",
      "- Four-line excludes SUM-159-fuse; all-line and unified single-line panels include it.",
      "- Visual QC reduced structural typography, compacted strategy labels, moved both color guides to the right, and retained only endpoint guide labels.",
      "- Project-map decision: no project_map update is needed while this remains a review-stage supplementary candidate."
    ),
    file.path(polish_root, "notes.md")
  )
}

if (phase == "subpanels") {
  bundle <- if (file.exists(cache_path)) {
    qread(cache_path, nthreads = 2L)
  } else {
    build_summary()
  }
  panels <- make_panels(bundle)
  write_subpanels(panels)
  cat("Wrote six regenerated audit subpanels.\n")
}

if (phase == "final") {
  bundle <- if (file.exists(cache_path)) {
    qread(cache_path, nthreads = 2L)
  } else {
    build_summary()
  }
  panels <- make_panels(bundle)
  layout <- fread(file.path(layout_dir, "layout_plan.csv"))
  layout <- layout[
    figure == "Posterior strategy context supplement"
  ][match(letters[1:6], panel)]
  if (nrow(layout) != 6L || anyNA(layout$panel)) {
    stop("Layout plan does not contain panels a-f.", call. = FALSE)
  }
  canvas <- ggdraw()
  for (i in seq_len(nrow(layout))) {
    row <- layout[i]
    canvas <- canvas +
      draw_plot(
        panels[[row$panel]],
        x = row$x_npc,
        y = row$y_npc,
        width = row$width_npc,
        height = row$height_npc
      ) +
      draw_label(
        row$panel,
        x = row$x_npc + 0.003,
        y = row$y_npc + row$height_npc - 0.003,
        hjust = 0,
        vjust = 1,
        size = 8,
        fontface = "bold",
        fontfamily = "sans"
      )
  }
  final_path <- file.path(
    final_dir, "posterior_strategy_context_support_median.png"
  )
  ggsave(
    final_path,
    canvas,
    width = unique(layout$layout_width_in),
    height = unique(layout$layout_height_in),
    dpi = 600,
    bg = "white"
  )
  write_metadata(final_path)
  cat("Wrote supplementary composite: ", final_path, "\n", sep = "")
}
