#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/optim_utils.R")
source("R/parameter_transfer_utils.R")
source("R/model_analysis_utils.R")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
package_root <- dirname(dirname(script_path))
output_dir <- file.path(package_root, "source_tables")
current_root <- "data/modelling/gpath_v1/morphology_metrics_a30_20260729"
reference_root <- "data/modelling/gpath_v1/red_a30_counts_20260722"

reference_assessment <- readRDS(
  file.path(reference_root, "derived/optimization/assessment.Rds")
)$fit_summary
pareto <- reference_assessment[
  reference_assessment$dataset_id == "all_lines" &
    reference_assessment$pareto_member %in% TRUE &
    is.finite(reference_assessment$deviance),
  c("model_id", "pareto_rank"),
  drop = FALSE
]
pareto$alias <- vapply(pareto$model_id, build_model_alias, character(1), format = "text")
pareto <- pareto[order(pareto$pareto_rank, pareto$model_id), , drop = FALSE]
if (nrow(pareto) != 5L || anyDuplicated(pareto$model_id)) {
  stop("Expected five unique all-lines Pareto models in the current count assessment")
}
pair_map <- read.delim(
  "config/parameter_estimation/morphology_metrics_a30_20260729/pair_map.tsv",
  stringsAsFactors = FALSE
)

cases <- data.frame(
  metric = c("ploidy", "cell_area", "nuclear_area"),
  metric_label = c("Ploidy", "Cell area", "Nuclear area"),
  root = c(reference_root, current_root, current_root),
  dataset_id = c("all_lines", "all_lines_cell_area", "all_lines_nuclear_area"),
  stringsAsFactors = FALSE
)
focus_parameters <- c("ae[1]", "ah[1]", "Y_R[1]", "m")
parameter_labels <- c(
  "ae[1]" = "v1",
  "ah[1]" = "K1",
  "Y_R[1]" = "y1",
  "m" = "m"
)
line_levels <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
line_colors <- c(
  "MCF10A" = "#1B9E77",
  "MDA-MB-231" = "#D95F02",
  "SNU668" = "#7570B3",
  "SUM-159-chem" = "#E7298A",
  "SUM-159-fuse" = "#66A61E"
)
metric_colors <- c(
  "Ploidy" = "#333333",
  "Cell area" = "#D55E00",
  "Nuclear area" = "#0072B2"
)

read_best_fit <- function(root, dataset_id, model_id) {
  assessment <- readRDS(
    file.path(root, "derived/optimization/assessment.Rds")
  )$fit_summary
  fit <- assessment[
    assessment$dataset_id == dataset_id & assessment$model_id == model_id,
    , drop = FALSE
  ]
  if (nrow(fit) != 1L) stop(sprintf("Missing assessment row: %s / %s", dataset_id, model_id))
  plan <- read.delim(
    file.path(root, "manifests/fit_plan.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  plan_row <- plan[plan$fit_id == fit$fit_id, , drop = FALSE]
  if (nrow(plan_row) != 1L) stop(sprintf("Missing plan row: %s", fit$fit_id))
  draws <- readRDS(file.path(plan_row$optim_dir, "optim_draws_all.Rds"))
  draw <- extract_draw_vector(draws[[fit$best_start]])
  rm(draws)
  gc(verbose = FALSE)
  list(
    draw = draw,
    stan_data = readRDS(plan_row$stan_data_path),
    fit_id = fit$fit_id,
    best_start = fit$best_start
  )
}

get_line_states <- function(stan_data, line_id, line_name) {
  idx <- which(as.integer(stan_data$line_id) == as.integer(line_id))
  states <- unique(data.frame(
    ploidy_abs = as.numeric(stan_data$ploidy_abs[idx]),
    metric_value = as.numeric(stan_data$ploidy_metric[idx])
  ))
  states <- merge(
    states,
    pair_map[pair_map$cellLine == line_name, c("ploidy_abs", "pair_state")],
    by = "ploidy_abs",
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(states) != 2L || anyNA(states$pair_state)) {
    stop(sprintf("Could not resolve pair states for %s", line_name))
  }
  states
}

effect_rows <- list()
unit_rows <- list()
effect_ptr <- 1L
unit_ptr <- 1L

for (case_i in seq_len(nrow(cases))) {
  metric <- cases$metric[[case_i]]
  metric_label <- cases$metric_label[[case_i]]
  for (model_i in seq_len(nrow(pareto))) {
    model_id <- pareto$model_id[[model_i]]
    model_alias <- pareto$alias[[model_i]]
    fit <- read_best_fit(
      cases$root[[case_i]],
      cases$dataset_id[[case_i]],
      model_id
    )
    stan_data <- fit$stan_data
    inverse_line_map <- setNames(
      names(stan_data$line_map),
      as.character(as.integer(unname(stan_data$line_map)))
    )
    line_ids <- sort(unique(as.integer(stan_data$line_id)))

    # The fitted global coefficient expressed as a physical-parameter log2
    # change for a +1 covariate-unit shift. It does not depend on cell line.
    unit_low <- reconstruct_line_state_parameters(
      fit$draw, model_id, line_ids[[1]], ploidy_metric = 0
    )
    unit_high <- reconstruct_line_state_parameters(
      fit$draw, model_id, line_ids[[1]], ploidy_metric = 1
    )
    for (parameter in focus_parameters) {
      unit_rows[[unit_ptr]] <- data.frame(
        metric = metric,
        metric_label = metric_label,
        model_id = model_id,
        model_alias = model_alias,
        pareto_rank = pareto$pareto_rank[[model_i]],
        parameter = parameter,
        parameter_label = unname(parameter_labels[[parameter]]),
        log2_effect_per_metric_unit = log2(unit_high[[parameter]] / unit_low[[parameter]]),
        fit_id = fit$fit_id,
        best_start = fit$best_start,
        stringsAsFactors = FALSE
      )
      unit_ptr <- unit_ptr + 1L
    }

    for (line_id in line_ids) {
      line_name <- unname(inverse_line_map[[as.character(line_id)]])
      states <- get_line_states(stan_data, line_id, line_name)
      baseline_metric <- states$metric_value[states$pair_state == "baseline"]
      elevated_metric <- states$metric_value[states$pair_state == "elevated"]
      baseline <- reconstruct_line_state_parameters(
        fit$draw, model_id, line_id, baseline_metric
      )
      elevated <- reconstruct_line_state_parameters(
        fit$draw, model_id, line_id, elevated_metric
      )
      for (parameter in focus_parameters) {
        effect_rows[[effect_ptr]] <- data.frame(
          metric = metric,
          metric_label = metric_label,
          model_id = model_id,
          model_alias = model_alias,
          pareto_rank = pareto$pareto_rank[[model_i]],
          line_id = line_id,
          cell_line = line_name,
          parameter = parameter,
          parameter_label = unname(parameter_labels[[parameter]]),
          baseline_metric = baseline_metric,
          elevated_metric = elevated_metric,
          metric_delta = elevated_metric - baseline_metric,
          baseline_parameter = baseline[[parameter]],
          elevated_parameter = elevated[[parameter]],
          log2_elevated_vs_baseline = log2(
            elevated[[parameter]] / baseline[[parameter]]
          ),
          fit_id = fit$fit_id,
          best_start = fit$best_start,
          stringsAsFactors = FALSE
        )
        effect_ptr <- effect_ptr + 1L
      }
    }
  }
}

effects <- do.call(rbind, effect_rows)
unit_effects <- do.call(rbind, unit_rows)
effects$metric_label <- factor(
  effects$metric_label,
  levels = cases$metric_label
)
effects$model_alias <- factor(effects$model_alias, levels = pareto$alias)
effects$parameter_label <- factor(
  effects$parameter_label,
  levels = c("v1", "K1", "y1", "m")
)
effects$cell_line <- factor(effects$cell_line, levels = line_levels)
unit_effects$metric_label <- factor(
  unit_effects$metric_label,
  levels = cases$metric_label
)
unit_effects$model_alias <- factor(unit_effects$model_alias, levels = rev(pareto$alias))
unit_effects$parameter_label <- factor(
  unit_effects$parameter_label,
  levels = c("v1", "K1", "y1", "m")
)

p_realized <- ggplot(
  effects,
  aes(
    metric_label,
    log2_elevated_vs_baseline,
    colour = cell_line,
    group = cell_line
  )
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.6, alpha = 0.72) +
  geom_point(size = 1.8) +
  facet_grid(
    rows = vars(parameter_label),
    cols = vars(model_alias),
    scales = "free_y"
  ) +
  scale_colour_manual(values = line_colors, name = "Cell line") +
  labs(
    title = "Do inferred pair effects persist across covariate definitions?",
    subtitle = paste(
      "Best MAP fits for the five Pareto models;",
      "each line connects the same experimental pair"
    ),
    x = "Ploidy-effect covariate",
    y = "log2(elevated-state parameter / baseline-state parameter)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  )

ggsave(
  file.path(output_dir, "06_realized_parameter_effects_by_metric.png"),
  p_realized, width = 13.5, height = 9.0, dpi = 180, bg = "white"
)

p_unit <- ggplot(
  unit_effects,
  aes(
    log2_effect_per_metric_unit,
    model_alias,
    colour = metric_label
  )
) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.4
  ) +
  facet_wrap(~parameter_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = metric_colors, name = "Covariate") +
  labs(
    title = "Fitted global ploidy-effect coefficients",
    subtitle = "Physical-parameter change for a +1 unit shift in the fitted covariate",
    x = "log2 parameter change per +1 covariate unit",
    y = "Pareto model"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(
  file.path(output_dir, "07_per_unit_parameter_effect_coefficients.png"),
  p_unit, width = 12.5, height = 4.8, dpi = 180, bg = "white"
)

write.csv(
  effects,
  file.path(output_dir, "realized_parameter_effects_by_metric.csv"),
  row.names = FALSE
)
write.csv(
  unit_effects,
  file.path(output_dir, "per_unit_parameter_effect_coefficients.csv"),
  row.names = FALSE
)

unit_reference <- unit_effects[
  unit_effects$metric == "ploidy",
  c("model_id", "parameter", "log2_effect_per_metric_unit")
]
names(unit_reference)[[3]] <- "ploidy_effect"
unit_concordance <- merge(
  unit_effects[unit_effects$metric != "ploidy", , drop = FALSE],
  unit_reference,
  by = c("model_id", "parameter"),
  all.x = TRUE,
  sort = FALSE
)
unit_concordance$sign_matches_ploidy <-
  sign(unit_concordance$log2_effect_per_metric_unit) ==
  sign(unit_concordance$ploidy_effect)
write.csv(
  unit_concordance,
  file.path(output_dir, "per_unit_parameter_effect_sign_concordance.csv"),
  row.names = FALSE
)
unit_concordance_summary <- aggregate(
  sign_matches_ploidy ~ metric_label + parameter_label,
  data = unit_concordance,
  FUN = function(x) sprintf("%d/%d", sum(x), length(x))
)
write.csv(
  unit_concordance_summary,
  file.path(output_dir, "per_unit_parameter_effect_sign_concordance_summary.csv"),
  row.names = FALSE
)

cat(sprintf("Wrote parameter-effect plots and tables to %s\n", output_dir))
