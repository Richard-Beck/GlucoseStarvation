#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
source("R/optim_utils.R")
source("R/model_analysis_utils.R")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
package_root <- dirname(dirname(script_path))
output_dir <- file.path(package_root, "source_tables")
current_root <- "data/modelling/gpath_v1/morphology_metrics_a30_20260729"
reference_root <- "data/modelling/gpath_v1/red_a30_counts_20260722"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

current_assessment <- readRDS(
  file.path(current_root, "derived/optimization/assessment.Rds")
)$fit_summary
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
current_plan <- read.delim(
  file.path(current_root, "manifests/fit_plan.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
reference_plan <- read.delim(
  file.path(reference_root, "manifests/fit_plan.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

read_best_draw <- function(assessment, plan, dataset_id, model_id) {
  fit <- assessment[
    assessment$dataset_id == dataset_id & assessment$model_id == model_id,
    , drop = FALSE
  ]
  if (nrow(fit) != 1L) {
    stop(sprintf("Expected one assessment row for %s / %s", dataset_id, model_id))
  }
  plan_row <- plan[plan$fit_id == fit$fit_id, , drop = FALSE]
  if (nrow(plan_row) != 1L) stop(sprintf("Missing fit plan row for %s", fit$fit_id))
  draws <- readRDS(file.path(plan_row$optim_dir, "optim_draws_all.Rds"))
  draw <- extract_draw_vector(draws[[fit$best_start]])
  list(draw = draw, stan_data_path = plan_row$stan_data_path, fit_id = fit$fit_id)
}

extract_ll_well <- function(draw, n_wells) {
  expected <- sprintf("ll_well[%d]", seq_len(n_wells))
  if (!all(expected %in% names(draw))) stop("Best draw is missing one or more ll_well entries")
  as.numeric(draw[expected])
}

reference_data <- readRDS(file.path(reference_root, "datasets/all_lines/stan_data.Rds"))
inverse_line_map <- setNames(
  names(reference_data$line_map),
  as.character(as.integer(unname(reference_data$line_map)))
)
well_metadata <- data.frame(
  well_idx = seq_len(reference_data$N_wells),
  cell_line = unname(inverse_line_map[as.character(reference_data$line_id)]),
  glucose_mM = as.numeric(reference_data$G0_per_well),
  stringsAsFactors = FALSE
)
if (anyNA(well_metadata$cell_line)) stop("Could not map every well to a cell line")

metric_datasets <- c(
  nuclear_area = "all_lines_nuclear_area",
  cell_area = "all_lines_cell_area"
)
condition_rows <- list()
ptr <- 1L
for (model_i in seq_len(nrow(pareto))) {
  model_id <- pareto$model_id[[model_i]]
  reference_fit <- read_best_draw(
    reference_assessment, reference_plan, "all_lines", model_id
  )
  reference_ll <- extract_ll_well(reference_fit$draw, reference_data$N_wells)

  for (metric in names(metric_datasets)) {
    current_fit <- read_best_draw(
      current_assessment, current_plan, metric_datasets[[metric]], model_id
    )
    current_data <- readRDS(current_fit$stan_data_path)
    if (current_data$N_wells != reference_data$N_wells) {
      stop(sprintf("Well-count mismatch for %s / %s", metric, model_id))
    }
    current_ll <- extract_ll_well(current_fit$draw, current_data$N_wells)
    detail <- well_metadata
    detail$model_id <- model_id
    detail$model_alias <- pareto$alias[[model_i]]
    detail$pareto_rank <- pareto$pareto_rank[[model_i]]
    detail$metric <- metric
    detail$delta_log_lik <- current_ll - reference_ll

    total_from_wells <- sum(detail$delta_log_lik)
    total_from_scalar <- as.numeric(current_fit$draw[["log_lik"]]) -
      as.numeric(reference_fit$draw[["log_lik"]])
    if (!isTRUE(all.equal(total_from_wells, total_from_scalar, tolerance = 1e-4))) {
      stop(sprintf(
        "ll_well does not sum to log_lik for %s / %s: wells=%g scalar=%g diff=%g",
        metric, model_id, total_from_wells, total_from_scalar,
        total_from_wells - total_from_scalar
      ))
    }
    condition_rows[[ptr]] <- detail
    ptr <- ptr + 1L
  }
}
well_contributions <- do.call(rbind, condition_rows)
condition_contributions <- aggregate(
  delta_log_lik ~ model_id + model_alias + pareto_rank + metric +
    cell_line + glucose_mM,
  data = well_contributions,
  FUN = sum
)
well_counts <- aggregate(
  well_idx ~ model_id + model_alias + pareto_rank + metric +
    cell_line + glucose_mM,
  data = well_contributions,
  FUN = length
)
names(well_counts)[names(well_counts) == "well_idx"] <- "n_wells"
condition_contributions <- merge(
  condition_contributions,
  well_counts,
  by = c(
    "model_id", "model_alias", "pareto_rank", "metric",
    "cell_line", "glucose_mM"
  ),
  sort = FALSE
)
across_models <- aggregate(
  delta_log_lik ~ metric + cell_line + glucose_mM,
  data = condition_contributions,
  FUN = function(x) c(mean = mean(x), sum = sum(x), min = min(x), max = max(x))
)
across_values <- as.data.frame(across_models$delta_log_lik)
across_models$delta_log_lik <- NULL
across_models <- cbind(across_models, across_values)
names(across_models)[(ncol(across_models) - 3L):ncol(across_models)] <-
  c("mean_delta_log_lik", "sum_delta_log_lik", "min_delta_log_lik", "max_delta_log_lik")

line_levels <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
metric_levels <- c("nuclear_area", "cell_area")
metric_labels <- c(nuclear_area = "Nuclear area", cell_area = "Cell area")
condition_contributions$cell_line <- factor(
  condition_contributions$cell_line,
  levels = rev(line_levels)
)
condition_contributions$model_alias <- factor(
  condition_contributions$model_alias,
  levels = pareto$alias
)
condition_contributions$metric <- factor(
  condition_contributions$metric,
  levels = metric_levels
)
glucose_levels <- sort(unique(condition_contributions$glucose_mM))
condition_contributions$glucose_label <- factor(
  format(glucose_levels[match(condition_contributions$glucose_mM, glucose_levels)],
         trim = TRUE, scientific = FALSE),
  levels = format(glucose_levels, trim = TRUE, scientific = FALSE)
)

fill_limit <- max(abs(condition_contributions$delta_log_lik), na.rm = TRUE)
p_conditions <- ggplot(
  condition_contributions,
  aes(glucose_label, cell_line, fill = delta_log_lik)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_grid(
    rows = vars(model_alias),
    cols = vars(metric),
    labeller = labeller(metric = metric_labels)
  ) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-fill_limit, fill_limit),
    name = "Delta log\nlikelihood"
  ) +
  labs(
    title = "Where area metrics changed fit quality",
    subtitle = paste(
      "Sum of well-level log-likelihood changes versus ploidy;",
      "blue favours the area metric"
    ),
    x = "Initial glucose (mM)", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = grid::unit(0.55, "lines"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "04_condition_delta_loglik_heatmap.png"),
  p_conditions, width = 11.5, height = 12.0, dpi = 180, bg = "white"
)

# Extract the exact covariate values used by each of the three fitted datasets.
metric_sources <- c(
  ploidy = file.path(reference_root, "datasets/all_lines/stan_data.Rds"),
  cell_area = file.path(current_root, "datasets/all_lines_cell_area/stan_data.Rds"),
  nuclear_area = file.path(current_root, "datasets/all_lines_nuclear_area/stan_data.Rds")
)
pair_map <- read.delim(
  "config/parameter_estimation/morphology_metrics_a30_20260729/pair_map.tsv",
  stringsAsFactors = FALSE
)
metric_value_rows <- lapply(names(metric_sources), function(metric) {
  x <- readRDS(metric_sources[[metric]])
  inverse <- setNames(
    names(x$line_map),
    as.character(as.integer(unname(x$line_map)))
  )
  values <- unique(data.frame(
    cellLine = unname(inverse[as.character(x$line_id)]),
    ploidy_abs = as.numeric(x$ploidy_abs),
    ploidy_effect_metric = as.numeric(x$ploidy_metric),
    stringsAsFactors = FALSE
  ))
  values <- merge(
    values, pair_map,
    by = c("cellLine", "ploidy_abs"),
    all.x = TRUE, sort = FALSE
  )
  values$metric <- metric
  values
})
metric_values <- do.call(rbind, metric_value_rows)
if (anyNA(metric_values$pair_state)) stop("Metric values did not fully match the pair map")
metric_values$cellLine <- factor(metric_values$cellLine, levels = rev(line_levels))
metric_values$metric <- factor(
  metric_values$metric,
  levels = c("ploidy", "cell_area", "nuclear_area")
)

elevated <- metric_values[metric_values$pair_state == "elevated", , drop = FALSE]
metric_plot_labels <- c(
  ploidy = "Ploidy",
  cell_area = "Cell area",
  nuclear_area = "Nuclear area"
)
metric_limit <- max(abs(elevated$ploidy_effect_metric), na.rm = TRUE)
p_metric_values <- ggplot(
  elevated,
  aes(metric, cellLine, fill = ploidy_effect_metric)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.4f", ploidy_effect_metric)),
    size = 4
  ) +
  scale_x_discrete(labels = metric_plot_labels, expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-metric_limit, metric_limit),
    name = "Metric value"
  ) +
  labs(
    title = "Covariate values used for the elevated member of each pair",
    subtitle = "Within-line baseline values are zero in all three definitions",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_dir, "05_ploidy_effect_metric_values.png"),
  p_metric_values, width = 7.8, height = 5.3, dpi = 180, bg = "white"
)

write.csv(
  well_contributions,
  file.path(output_dir, "condition_delta_loglik_by_well.csv"),
  row.names = FALSE
)
write.csv(
  condition_contributions,
  file.path(output_dir, "condition_delta_loglik_summary.csv"),
  row.names = FALSE
)
write.csv(
  across_models,
  file.path(output_dir, "condition_delta_loglik_across_models.csv"),
  row.names = FALSE
)
write.csv(
  metric_values,
  file.path(output_dir, "ploidy_effect_metric_values.csv"),
  row.names = FALSE
)

cat(sprintf("Wrote condition and metric-value plots to %s\n", output_dir))
