#!/usr/bin/env Rscript

# FG6 selection/competition figures, regenerated from the canonical
# red_a30_counts_20260722 fits. All simulations and reductions remain in memory.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(qs)
  library(scales)
  library(tidyr)
})

source("R/selection_strategy_utils.R")
source("R/gpath_posterior_strategy_utils.R")

total_started <- proc.time()[["elapsed"]]
package_id <- "FG6_selection_simulation_revision"
release_root <- file.path("data", "modelling", "gpath_v1", "red_a30_counts_20260722")
target_root <- file.path("manuscript_figures", "20260722T141236_resegmentation_fix")
final_dir <- file.path(target_root, "final_images")
timing_dir <- file.path(target_root, "timings")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
cores <- suppressWarnings(as.integer(Sys.getenv("FG6_CORES", "8")))
if (!is.finite(cores) || cores < 1L) cores <- 1L
cores <- min(cores, parallel::detectCores(logical = FALSE))

models <- c(
  "1R_1P_0W_C0_M1",
  "1R_1P_1W_C0_M0",
  "2R_1P_0W_C0_M1",
  "2R_2P_0W_C0_M1",
  "2R_2P_1W_C0_M0"
)
primary_model <- "2R_2P_0W_C0_M1"
model_alias <- c(
  "1R_1P_0W_C0_M1" = "1R",
  "1R_1P_1W_C0_M0" = "1R,W",
  "2R_1P_0W_C0_M1" = "2R(a)",
  "2R_2P_0W_C0_M1" = "2R(f)",
  "2R_2P_1W_C0_M0" = "2R(f),W"
)
model_role <- c(
  "1R_1P_0W_C0_M1" = "legacy audit",
  "1R_1P_1W_C0_M0" = "robustness",
  "2R_1P_0W_C0_M1" = "robustness",
  "2R_2P_0W_C0_M1" = "primary",
  "2R_2P_1W_C0_M0" = "robustness"
)
model_labels <- paste0(unname(model_alias[models]), "\n", unname(model_role[models]))
names(model_labels) <- models

datasets <- c(no_SUM = "loo_exclude_sum_159_fuse", `SUM incl.` = "all_lines")
objective_levels <- c("Resource starvation\n(favor high ploidy)", "Resource abundance\n(favor low ploidy)")
objective_palette <- c(
  "Resource starvation\n(favor high ploidy)" = "#B23A48",
  "Resource abundance\n(favor low ploidy)" = "#2A6FBB"
)
line_palette <- c(
  "MCF10A" = "#0072B2",
  "MDA-MB-231" = "#D55E00",
  "SNU668" = "#009E73",
  "SUM-159-chem" = "#CC79A7",
  "SUM-159-fuse" = "#E69F00"
)

objective_label <- function(x) {
  recode(
    as.character(x),
    for_high = objective_levels[[1]],
    against_high = objective_levels[[2]]
  )
}

mode_string <- function(x) {
  x <- as.character(x)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

theme_fg6 <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#202020"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E4E4DF", linewidth = 0.22),
      panel.spacing = unit(4, "pt"),
      strip.background = element_rect(fill = "#ECEDE7", color = NA),
      strip.text = element_text(face = "bold", lineheight = 0.92),
      axis.text = element_text(color = "#252525"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

timings <- tibble(
  package_id = character(),
  operation = character(),
  elapsed_seconds = numeric(),
  details = character()
)
record_timing <- function(operation, started, details) {
  timings <<- bind_rows(
    timings,
    tibble(
      package_id = package_id,
      operation = operation,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      details = details
    )
  )
}

load_started <- proc.time()[["elapsed"]]
strategy_root <- file.path(release_root, "derived", "posterior", "strategy_simulations")
strategy_validation <- jsonlite::fromJSON(
  file.path(strategy_root, "validation.json"),
  simplifyVector = TRUE
)
if (!isTRUE(strategy_validation$validation_passed) ||
    as.integer(strategy_validation$observed$valid_tasks) != 70L) {
  stop("Canonical posterior strategy simulations have not passed validation.", call. = FALSE)
}
posterior_contexts <- list(
  no_SUM = readRDS(file.path(strategy_root, "endpoints", "loo_exclude_sum_159_fuse.Rds")),
  `SUM incl.` = readRDS(file.path(strategy_root, "endpoints", "all_lines.Rds"))
)
if (!identical(unname(vapply(posterior_contexts, function(x) x$metadata$source_tasks, integer(1))), c(20L, 25L))) {
  stop("Unexpected endpoint-reduction task coverage.", call. = FALSE)
}
fit_plan <- read.delim(
  file.path(release_root, "manifests", "fit_plan.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assessment <- readRDS(file.path(release_root, "derived", "optimization", "assessment.Rds"))
fit_summary <- assessment$fit_summary
fit_plan <- fit_plan |>
  filter(dataset_id %in% unname(datasets), model_id %in% models, run_optim == 1L)
if (nrow(fit_plan) != length(datasets) * length(models)) {
  stop("Canonical fit plan does not contain the expected ten dataset/model fits.", call. = FALSE)
}
stan_cache <- lapply(unname(datasets), function(dataset_id) {
  x <- readRDS(file.path(release_root, "datasets", dataset_id, "stan_data.Rds"))
  if (is.null(x$group_id)) x <- add_group_structure(x)
  x
})
names(stan_cache) <- unname(datasets)
record_timing(
  "load",
  load_started,
  "Verified 70/70 full-strategy validation; loaded two compact posterior endpoint shards plus canonical optimization/Stan inputs retained for unchanged initial-condition panels."
)

load_best_draw <- function(dataset_id, model_id) {
  row <- fit_plan |>
    filter(.data$dataset_id == .env$dataset_id, .data$model_id == .env$model_id)
  if (nrow(row) != 1L) stop("Missing unique canonical fit: ", dataset_id, " / ", model_id, call. = FALSE)
  draws <- readRDS(file.path(row$optim_dir[[1]], "optim_draws_all.Rds"))
  lp <- readRDS(file.path(row$optim_dir[[1]], "optim_lp_all.Rds"))
  extract_best_draw_from_optim_outputs(draws_all = draws, lp_all = lp)$draw_vec
}

prepare_context <- function(dataset_id, model_id, line_id) {
  stan <- stan_cache[[dataset_id]]
  draw <- load_best_draw(dataset_id, model_id)
  pair <- get_line_ploidy_pair(stan, line_id)
  comp <- reconstruct_mixture_components(
    draw_vec = draw,
    model_id = model_id,
    line_id = line_id,
    low_ploidy = pair$low,
    high_ploidy = pair$high
  )
  seed <- gps_estimate_mixture_seed_total_n(
    draw_vec = draw,
    stan_data = stan,
    line_id = line_id,
    low_ploidy = pair$low,
    high_ploidy = pair$high
  )
  if (!is.null(stan$line_map) && length(stan$line_map)) {
    line_name <- names(stan$line_map)[match(as.integer(line_id), as.integer(stan$line_map))]
  } else {
    labels <- get_line_label_map(stan)
    line_name <- labels$cell_line[match(line_id, labels$line_id)]
  }
  if (length(line_name) != 1L || is.na(line_name) || !nzchar(line_name)) {
    stop("Could not resolve canonical cell-line label for ", dataset_id, " line ", line_id, call. = FALSE)
  }
  list(
    parms = comp$parms_mix,
    seed = seed,
    pair = pair,
    line_name = line_name,
    glucose = resolve_glucose_values(stan, list(glucose_values = NULL))
  )
}

simulate_one <- function(strategy_row, context, seed_total_n = context$seed,
                         initial_high_fraction = 0.5, detailed = FALSE) {
  simulate_strategy(
    strategy_row = strategy_row,
    parms_mix = context$parms,
    seed_total_n = seed_total_n,
    initial_high_fraction = initial_high_fraction,
    interval_hours = 48,
    time_step_hours = 1,
    reset_total_n_on_refresh = FALSE,
    detailed = detailed
  )
}

posterior_started <- proc.time()[["elapsed"]]
strategy_grid <- posterior_contexts$no_SUM$schedule_index |>
  mutate(
    day0_action = paste0(format(g0_day0, trim = TRUE), " mM"),
    day2_action = if_else(is.na(g0_day2), "carry", paste0(format(g0_day2, trim = TRUE), " mM")),
    day4_action = if_else(is.na(g0_day4), "carry", paste0(format(g0_day4, trim = TRUE), " mM"))
  )
action_levels_day0 <- paste0(format(c(0, 0.1, 0.25, 0.5, 1, 5, 25), trim = TRUE), " mM")
action_levels_refresh <- c("carry", action_levels_day0)

rank_endpoint_context <- function(artifact, scope_name) {
  metric_names <- dimnames(artifact$endpoint)$metric
  metric <- function(name) match(name, metric_names)
  schedule_codes <- as.character(artifact$schedule_index$strategy_code)
  rows <- vector("list", nrow(artifact$task_index) * 2L)
  out_i <- 1L
  for (task_i in seq_len(nrow(artifact$task_index))) {
    task <- artifact$task_index[task_i, , drop = FALSE]
    log_ratio <- artifact$endpoint[task_i, , , metric("log_ratio_high_low")]
    total_live <- artifact$endpoint[task_i, , , metric("total_live")]
    folds <- list(
      for_high = artifact$endpoint[task_i, , , metric("high_fold_change")],
      against_high = artifact$endpoint[task_i, , , metric("low_fold_change")]
    )
    for (objective in c("for_high", "against_high")) {
      score <- if (objective == "for_high") log_ratio else -log_ratio
      score[folds[[objective]] < 1 | !is.finite(folds[[objective]])] <- NA_real_
      chosen <- vapply(seq_len(nrow(score)), function(draw_i) {
        viable <- which(is.finite(score[draw_i, ]))
        if (!length(viable)) {
          stop("No viable posterior schedule for ", scope_name, " task ", task_i,
               " draw ", draw_i, " objective ", objective, call. = FALSE)
        }
        viable[order(score[draw_i, viable], total_live[draw_i, viable],
                     decreasing = TRUE)[[1L]]]
      }, integer(1))
      draw_rows <- seq_len(nrow(score))
      chosen_log_ratio <- log_ratio[cbind(draw_rows, chosen)]
      chosen_objective_score <- if (objective == "for_high") chosen_log_ratio else -chosen_log_ratio
      rows[[out_i]] <- tibble(
        scope = scope_name,
        dataset_id = task$dataset_id,
        model_id = task$model_id,
        line_id = as.integer(task$line_id),
        cell_line = task$line_name,
        draw_id = artifact$draw_index[[task_i]]$draw_id,
        objective,
        strategy_code = schedule_codes[chosen],
        log_ratio_final = chosen_log_ratio,
        high_fraction_final = artifact$endpoint[task_i, , , metric("high_fraction")][cbind(draw_rows, chosen)],
        total_live = total_live[cbind(draw_rows, chosen)],
        low_fold_change = artifact$endpoint[task_i, , , metric("low_fold_change")][cbind(draw_rows, chosen)],
        high_fold_change = artifact$endpoint[task_i, , , metric("high_fold_change")][cbind(draw_rows, chosen)],
        seed_total_n = artifact$seed_total_n[task_i, ],
        objective_score = chosen_objective_score
      )
      out_i <- out_i + 1L
    }
  }
  bind_rows(rows)
}

best <- bind_rows(
  rank_endpoint_context(posterior_contexts$no_SUM, "no_SUM"),
  rank_endpoint_context(posterior_contexts[["SUM incl."]], "SUM incl.")
)
if (nrow(best) != 9000L) stop("Unexpected posterior best-row coverage.", call. = FALSE)

primary_modal <- best |>
  filter(scope == "no_SUM", model_id == primary_model) |>
  group_by(objective) |>
  summarize(strategy_code = mode_string(strategy_code), .groups = "drop")

aggregate_draws <- best |>
  group_by(scope, model_id, objective, draw_id) |>
  summarize(
    mean_best_log_ratio = mean(log_ratio_final),
    objective_success = if_else(
      first(objective) == "for_high",
      mean_best_log_ratio > 0,
      mean_best_log_ratio < 0
    ),
    .groups = "drop"
  )
aggregate_df <- aggregate_draws |>
  group_by(scope, model_id, objective) |>
  summarize(
    mean_best_log_ratio = median(mean_best_log_ratio),
    q10_best_log_ratio = quantile(mean_best_log_ratio, 0.1),
    q90_best_log_ratio = quantile(mean_best_log_ratio, 0.9),
    success_fraction = mean(objective_success),
    .groups = "drop"
  ) |>
  mutate(
    scope_display = factor(scope, levels = c("no_SUM", "SUM incl."), labels = c("no-SUM", "SUM incl.")),
    objective_panel = factor(objective_label(objective), levels = objective_levels),
    model_label = factor(model_labels[model_id], levels = model_labels)
  )

no_sum_log_ratio <- posterior_contexts$no_SUM$endpoint[
  , , , match("log_ratio_high_low", dimnames(posterior_contexts$no_SUM$endpoint)$metric)
]
s11a_source <- tibble(
  strategy_code = as.character(posterior_contexts$no_SUM$schedule_index$strategy_code),
  posterior_median_log_ratio = apply(no_sum_log_ratio, 3L, median)
) |>
  left_join(strategy_grid |> select(strategy_code, day0_action, day2_action, day4_action), by = "strategy_code") |>
  mutate(
    day0_action = factor(day0_action, levels = action_levels_day0),
    day2_action = factor(day2_action, levels = rev(action_levels_refresh)),
    day4_action = factor(day4_action, levels = action_levels_refresh)
  )
s11_highlights <- primary_modal |>
  left_join(strategy_grid |> select(strategy_code, day0_action, day2_action, day4_action), by = "strategy_code") |>
  mutate(
    objective_panel = factor(objective_label(objective), levels = objective_levels),
    day0_action = factor(day0_action, levels = action_levels_day0),
    day2_action = factor(day2_action, levels = rev(action_levels_refresh)),
    day4_action = factor(day4_action, levels = action_levels_refresh)
  )
record_timing(
  "posterior_endpoint_ranking",
  posterior_started,
  "Applied low/high fold>=1 viability scoring per posterior draw and line; selected best schedules, modal primary schedules, endpoint distributions, and the 448-schedule posterior landscape."
)

dense_started <- proc.time()[["elapsed"]]
dense_parts <- list()
dense_i <- 1L
primary_sources <- posterior_contexts$no_SUM$task_index |>
  filter(model_id == primary_model)
for (source_i in seq_len(nrow(primary_sources))) {
  source_row <- primary_sources[source_i, , drop = FALSE]
  x <- qs::qread(source_row$output_path[[1]], nthreads = min(cores, 8L))
  for (objective in c("for_high", "against_high")) {
    strategy_code <- primary_modal$strategy_code[primary_modal$objective == objective]
    strategy_i <- match(strategy_code, as.character(x$schedule_index$strategy_code))
    if (!is.finite(strategy_i)) stop("Selected modal schedule is absent from full trajectory shard.", call. = FALSE)
    ntime <- length(x$time_hours)
    g1_i <- match("G1", x$state_names)
    low_i <- match("N_low", x$state_names)
    high_i <- match("N_high", x$state_names)
    g1 <- vapply(x$state, function(a) a[strategy_i, , g1_i], numeric(ntime))
    n_low <- vapply(x$state, function(a) a[strategy_i, , low_i], numeric(ntime))
    n_high <- vapply(x$state, function(a) a[strategy_i, , high_i], numeric(ntime))
    high_fraction <- n_high / pmax(n_low + n_high, 1e-12)
    high_fraction_median <- apply(high_fraction, 1L, median)
    high_fraction_q10 <- apply(high_fraction, 1L, quantile, probs = 0.1)
    high_fraction_q90 <- apply(high_fraction, 1L, quantile, probs = 0.9)
    dense_parts[[dense_i]] <- tibble(
      objective,
      strategy_code,
      cell_line = source_row$line_name,
      day = x$time_hours / 24,
      G1 = apply(g1, 1L, median),
      G1_q10 = apply(g1, 1L, quantile, probs = 0.1),
      G1_q90 = apply(g1, 1L, quantile, probs = 0.9),
      high_fraction = high_fraction_median,
      high_fraction_q10 = high_fraction_q10,
      high_fraction_q90 = high_fraction_q90,
      trajectory_id = paste(objective, source_row$line_name, sep = "::")
    )
    dense_i <- dense_i + 1L
  }
  rm(x)
  gc(verbose = FALSE)
}
dense <- bind_rows(dense_parts) |>
  mutate(objective_panel = factor(objective_label(objective), levels = objective_levels))
record_timing(
  "selected_posterior_trajectory_load",
  dense_started,
  "Loaded only four primary-model no-SUM full .qs shards and reduced the two modal schedules to posterior median and 10-90% ribbons."
)

deterministic_started <- proc.time()[["elapsed"]]
det_context_grid <- expand_grid(
  scope = "no_SUM",
  dataset_id = "loo_exclude_sum_159_fuse",
  model_id = models,
  line_id = sort(unique(as.integer(stan_cache[["loo_exclude_sum_159_fuse"]]$line_id)))
)
run_deterministic_context <- function(i) {
  task <- det_context_grid[i, , drop = FALSE]
  context <- prepare_context(task$dataset_id, task$model_id, task$line_id)
  rows <- vector("list", nrow(strategy_grid))
  for (j in seq_len(nrow(strategy_grid))) {
    rows[[j]] <- simulate_one(strategy_grid[j, , drop = FALSE], context, detailed = FALSE)$final_summary
  }
  bind_rows(rows) |>
    transmute(
      scope = task$scope, dataset_id = task$dataset_id, model_id = task$model_id,
      line_id = as.integer(task$line_id), cell_line = context$line_name,
      strategy_code = as.character(strategy_code), log_ratio_final,
      high_fraction_final, total_live, low_fold_change, high_fold_change,
      for_high_viable = if_else(high_fold_change >= 1, log_ratio_final, NA_real_),
      against_high_viable = if_else(low_fold_change >= 1, -log_ratio_final, NA_real_)
    )
}
det_parts <- if (.Platform$OS.type != "windows" && cores > 1L) {
  parallel::mclapply(seq_len(nrow(det_context_grid)), run_deterministic_context, mc.cores = cores)
} else {
  lapply(seq_len(nrow(det_context_grid)), run_deterministic_context)
}
det_errors <- vapply(det_parts, inherits, logical(1), what = "try-error")
if (any(det_errors)) stop("Deterministic IC baseline worker failure: ", det_parts[[which(det_errors)[1L]]], call. = FALSE)
det_baseline <- bind_rows(det_parts)
det_best <- bind_rows(
  det_baseline |>
    group_by(dataset_id, model_id, line_id, cell_line) |>
    arrange(is.na(for_high_viable), desc(for_high_viable), desc(total_live), .by_group = TRUE) |>
    slice(1L) |>
    ungroup() |>
    mutate(objective = "for_high"),
  det_baseline |>
    group_by(dataset_id, model_id, line_id, cell_line) |>
    arrange(is.na(against_high_viable), desc(against_high_viable), desc(total_live), .by_group = TRUE) |>
    slice(1L) |>
    ungroup() |>
    mutate(objective = "against_high")
)
record_timing(
  "unchanged_initial_condition_baseline",
  deterministic_started,
  sprintf("Recomputed the prior deterministic no-SUM baseline schedules solely to preserve F5d/S11b identity (%d cores).", cores)
)

ic_started <- proc.time()[["elapsed"]]
ic_grid <- bind_rows(
  expand_grid(initial_cell_num = c(50, 100, 200, 400, 800, 1600, 3200),
              initial_hi_ploidy_fraction = seq(0.1, 0.9, 0.1))
) |>
  distinct() |>
  mutate(is_baseline_anchor = FALSE)
ic_tasks <- det_best |>
  select(dataset_id, model_id, line_id, cell_line, objective, strategy_code)

run_ic_context <- function(i) {
  task <- ic_tasks[i, , drop = FALSE]
  context <- prepare_context(task$dataset_id, task$model_id, task$line_id)
  strategy <- strategy_grid |> filter(.data$strategy_code == task$strategy_code[[1]])
  task_grid <- bind_rows(
    ic_grid,
    tibble(
      initial_cell_num = context$seed,
      initial_hi_ploidy_fraction = 0.5,
      is_baseline_anchor = TRUE
    )
  )
  rows <- vector("list", nrow(task_grid))
  for (j in seq_len(nrow(task_grid))) {
    sim <- simulate_one(
      strategy,
      context,
      seed_total_n = task_grid$initial_cell_num[[j]],
      initial_high_fraction = task_grid$initial_hi_ploidy_fraction[[j]],
      detailed = FALSE
    )$final_summary
    rows[[j]] <- tibble(
      model_id = task$model_id,
      line_id = task$line_id,
      cell_line = task$cell_line,
      objective = task$objective,
      strategy_code = task$strategy_code,
      initial_cell_num = task_grid$initial_cell_num[[j]],
      initial_hi_ploidy_fraction = task_grid$initial_hi_ploidy_fraction[[j]],
      is_baseline_anchor = task_grid$is_baseline_anchor[[j]],
      log_ratio_final = sim$log_ratio_final,
      high_fraction_final = sim$high_fraction_final,
      total_live = sim$total_live
    )
  }
  bind_rows(rows)
}
ic_parts <- if (.Platform$OS.type != "windows" && cores > 1L) {
  parallel::mclapply(seq_len(nrow(ic_tasks)), run_ic_context, mc.cores = cores)
} else {
  lapply(seq_len(nrow(ic_tasks)), run_ic_context)
}
ic_errors <- vapply(ic_parts, inherits, logical(1), what = "try-error")
if (any(ic_errors)) {
  stop("Initial-condition worker failure: ", as.character(ic_parts[[which(ic_errors)[1L]]]), call. = FALSE)
}
ic <- bind_rows(ic_parts)
ic_baseline <- ic |>
  filter(is_baseline_anchor) |>
  select(model_id, line_id, objective, baseline_high_fraction = high_fraction_final)
ic_summary <- ic |>
  filter(!is_baseline_anchor) |>
  left_join(ic_baseline, by = c("model_id", "line_id", "objective")) |>
  mutate(delta_high_fraction = high_fraction_final - baseline_high_fraction) |>
  group_by(model_id, objective, initial_cell_num, initial_hi_ploidy_fraction) |>
  summarize(
    mean_delta_high_fraction = mean(delta_high_fraction),
    direction_flip_fraction = mean(
      if_else(objective == "for_high", log_ratio_final <= 0, log_ratio_final >= 0)
    ),
    .groups = "drop"
  ) |>
  mutate(
    objective_panel = factor(objective_label(objective), levels = objective_levels),
    model_label = factor(model_labels[model_id], levels = model_labels)
  )
record_timing(
  "initial_condition_simulation",
  ic_started,
  sprintf("Reintegrated 40 selected schedules over 63 explicit sensitivity conditions plus each context's fitted-N0 baseline, without refresh reseeding (%d simulations; %d cores).", nrow(ic), cores)
)

loo_started <- proc.time()[["elapsed"]]
make_score_rows <- function(scope_name) {
  x <- best |> filter(scope == scope_name)
  all_rows <- x |>
    group_by(model_id, objective, draw_id) |>
    summarize(
      objective_score = mean(objective_score),
      score_source = "all lines",
      .groups = "drop"
    )
  omitted <- bind_rows(lapply(sort(unique(x$cell_line)), function(omit) {
    x |>
      filter(cell_line != omit) |>
      group_by(model_id, objective, draw_id) |>
      summarize(
        objective_score = mean(objective_score),
        score_source = paste("omit", omit),
        .groups = "drop"
      )
  }))
  bind_rows(all_rows, omitted) |>
    mutate(context_label = if_else(scope_name == "no_SUM", "no-SUM", "SUM incl."))
}
score_draw_rows <- bind_rows(make_score_rows("no_SUM"), make_score_rows("SUM incl."))
score_rows <- score_draw_rows |>
  group_by(context_label, model_id, objective, score_source) |>
  summarize(
    posterior_q10 = quantile(objective_score, 0.1),
    posterior_q90 = quantile(objective_score, 0.9),
    success_fraction = mean(objective_score > 0),
    objective_score = median(objective_score),
    .groups = "drop"
  ) |>
  mutate(
    context_label = factor(context_label, levels = c("no-SUM", "SUM incl.")),
    objective_panel = factor(objective_label(objective), levels = objective_levels),
    model_label = factor(model_labels[model_id], levels = rev(model_labels))
  )
score_ranges <- score_rows |>
  group_by(context_label, objective_panel, model_label) |>
  summarize(min_score = min(objective_score), max_score = max(objective_score), .groups = "drop")
score_sources <- unique(score_rows$score_source)
score_palette <- c(
  "all lines" = "#222222",
  "omit MCF10A" = line_palette[["MCF10A"]],
  "omit MDA-MB-231" = line_palette[["MDA-MB-231"]],
  "omit SNU668" = line_palette[["SNU668"]],
  "omit SUM-159-chem" = line_palette[["SUM-159-chem"]],
  "omit SUM-159-fuse" = line_palette[["SUM-159-fuse"]]
)
record_timing(
  "leave_one_line_transforms",
  loo_started,
  "Recomputed posterior all-line and omit-one-line objective-score distributions in memory from the two endpoint arrays."
)

render_started <- proc.time()[["elapsed"]]
p5a <- ggplot(dense, aes(day, G1, color = cell_line, group = trajectory_id)) +
  geom_vline(xintercept = c(2, 4), color = "#D2D2CB", linewidth = 0.32) +
  geom_ribbon(aes(ymin = G1_q10, ymax = G1_q90, fill = cell_line),
              color = NA, alpha = 0.16) +
  geom_line(linewidth = 0.58, alpha = 0.92) +
  facet_wrap(~objective_panel, nrow = 1) +
  scale_color_manual(values = line_palette, name = "Cell line") +
  scale_fill_manual(values = line_palette, guide = "none") +
  scale_x_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.05), breaks = c(0, 0.1, 0.5, 1, 5, 25)) +
  labs(x = NULL, y = "Residual simulated\nmeasured G1 (mM)") +
  theme_fg6(7.2) +
  theme(legend.position = "none")

p5b <- ggplot(dense, aes(day, high_fraction, color = cell_line, group = trajectory_id)) +
  geom_vline(xintercept = c(2, 4), color = "#D2D2CB", linewidth = 0.32) +
  geom_hline(yintercept = 0.5, color = "#5C5C5C", linewidth = 0.3) +
  geom_ribbon(aes(ymin = high_fraction_q10, ymax = high_fraction_q90, fill = cell_line),
              color = NA, alpha = 0.16) +
  geom_line(linewidth = 0.62, alpha = 0.94) +
  facet_wrap(~objective_panel, nrow = 1) +
  scale_color_manual(values = line_palette, name = "Cell line") +
  scale_fill_manual(values = line_palette, guide = "none") +
  scale_x_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Time (days)", y = "High-ploidy\nfraction") +
  theme_fg6(7.2)

aggregate_limit <- max(1.25, ceiling((max(abs(c(aggregate_df$q10_best_log_ratio, aggregate_df$q90_best_log_ratio))) + 0.15) * 10) / 10)
p5c <- ggplot(aggregate_df, aes(model_label, mean_best_log_ratio, color = scope_display)) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.32) +
  geom_linerange(aes(ymin = q10_best_log_ratio, ymax = q90_best_log_ratio),
                 position = position_dodge(width = 0.48), linewidth = 0.62, alpha = 0.82) +
  geom_point(aes(size = success_fraction), position = position_dodge(width = 0.48), alpha = 0.96) +
  facet_wrap(~objective_panel, nrow = 1) +
  scale_color_manual(values = c("no-SUM" = "#1F6F8B", "SUM incl." = "#9C5A35"), name = "Scope") +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(1.4, 4), limits = c(0, 1), name = "Posterior\nsuccess") +
  scale_y_continuous(limits = c(-aggregate_limit, aggregate_limit), breaks = pretty_breaks(n = 5)) +
  labs(x = "Pareto model", y = "Posterior median day-6\nmean log(high / low)") +
  theme_fg6(7) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 6.2), legend.box = "horizontal")

p5d_source <- ic_summary |>
  filter(
    model_id != "1R_1P_0W_C0_M1",
    abs(initial_hi_ploidy_fraction - 0.5) < 1e-8
  )
p5d <- ggplot(p5d_source, aes(initial_cell_num, mean_delta_high_fraction, color = objective_panel, group = objective_panel)) +
  geom_hline(yintercept = 0, color = "#5B5B5B", linewidth = 0.32) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 1.5) +
  facet_wrap(~model_label, nrow = 2) +
  scale_color_manual(values = objective_palette, name = "Objective") +
  scale_x_log10(breaks = c(50, 100, 200, 400, 800, 1600, 3200), labels = label_number()) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Initial total cells (50:50 mixture)", y = "Change in high-ploidy\nfraction") +
  theme_fg6(6.35) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.2), strip.text = element_text(size = 5.9))

figure5 <- plot_grid(
  p5a, p5b, p5c, p5d,
  ncol = 1,
  rel_heights = c(1, 1.05, 1.22, 1.34),
  labels = c("a", "b", "c", "d"),
  label_size = 11,
  label_fontface = "bold",
  label_x = 0.006,
  label_y = 0.992,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr"
)

fill_limit <- max(abs(s11a_source$posterior_median_log_ratio))
p11a <- ggplot(s11a_source, aes(day0_action, 1)) +
  geom_tile(aes(fill = posterior_median_log_ratio), color = "white", linewidth = 0.08, height = 0.94, width = 0.94) +
  geom_point(
    data = s11_highlights,
    aes(x = day0_action, y = 1, shape = objective_panel),
    inherit.aes = FALSE,
    fill = "white",
    color = "#101010",
    size = 2,
    stroke = 0.35
  ) +
  facet_grid(day2_action ~ day4_action, switch = "y",
             labeller = labeller(day2_action = function(x) paste0("day 2: ", x),
                                 day4_action = function(x) paste0("day 4: ", x))) +
  scale_fill_gradient2(low = "#2A6FBB", mid = "#FAF7EF", high = "#B23A48",
                       midpoint = 0, limits = c(-fill_limit, fill_limit),
                       name = "Posterior median final\nlog(high / low)") +
  scale_shape_manual(values = c(24, 25), name = "Figure 5A/B\nexample objective") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "Day 0 glucose", y = NULL) +
  theme_fg6(5.2) +
  theme(
    panel.spacing = unit(1.2, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 3.6),
    strip.placement = "outside",
    strip.text.x = element_text(size = 4),
    strip.text.y.left = element_text(size = 4, angle = 0),
    legend.box = "vertical",
    legend.text = element_text(size = 4.7),
    legend.title = element_text(size = 4.8)
  )

p11b <- ggplot(
  ic_summary,
  aes(factor(initial_cell_num), factor(initial_hi_ploidy_fraction))
) +
  geom_tile(aes(fill = mean_delta_high_fraction), color = "white", linewidth = 0.12) +
  geom_point(
    data = ic_summary |> filter(direction_flip_fraction >= 0.5),
    shape = 4, size = 0.58, stroke = 0.32
  ) +
  facet_grid(objective_panel ~ model_label, drop = FALSE) +
  scale_fill_gradient2(low = "#2A6FBB", mid = "#FAF7EF", high = "#B23A48",
                       midpoint = 0, name = "change in\nhigh-ploidy\nfraction") +
  labs(x = "Initial total cells", y = "Initial high-ploidy fraction") +
  theme_fg6(5.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4.4),
    axis.text.y = element_text(size = 4.4),
    strip.text = element_text(size = 4.8),
    legend.position = "right"
  )

figure_s11 <- plot_grid(
  p11a, p11b,
  ncol = 1,
  rel_heights = c(4.96, 3.25),
  labels = c("a", "b"),
  label_size = 11,
  label_fontface = "bold",
  label_x = 0.006,
  label_y = 0.994,
  hjust = 0,
  vjust = 1
)

score_limit <- max(0.6, ceiling((max(abs(score_rows$objective_score)) + 0.1) * 10) / 10)
p12 <- ggplot(score_rows, aes(model_label, objective_score)) +
  geom_hline(yintercept = 0, color = "#575757", linewidth = 0.32) +
  geom_linerange(
    data = score_ranges,
    aes(x = model_label, ymin = min_score, ymax = max_score),
    inherit.aes = FALSE,
    color = "grey72",
    linewidth = 0.48
  ) +
  geom_point(
    aes(color = score_source, size = success_fraction, group = score_source),
    position = position_dodge(width = 0.62),
    alpha = 0.94
  ) +
  facet_grid(context_label ~ objective_panel) +
  scale_color_manual(values = score_palette, drop = FALSE, name = "Score source") +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(1.1, 3.8),
                        limits = c(0, 1), name = "Objective\nsuccess") +
  scale_y_continuous(limits = c(-score_limit, score_limit), breaks = pretty_breaks(n = 5)) +
  labs(x = "Model", y = "Objective-favoring day-6 log-ratio score") +
  theme_fg6(6.8) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1, size = 5.8),
    strip.text = element_text(size = 6.6),
    legend.box = "vertical"
  )
figure_s12 <- ggdraw() +
  draw_plot(p12, x = 0.018, y = 0.01, width = 0.976, height = 0.98) +
  draw_plot_label("c", x = 0.004, y = 0.996, hjust = 0, vjust = 1, size = 12, fontface = "bold")

ggsave(file.path(final_dir, "figure_5.png"), figure5, width = 7, height = 9.14, units = "in", dpi = 360, bg = "white")
ggsave(file.path(final_dir, "figure_s11.png"), figure_s11, width = 7, height = 8.29, units = "in", dpi = 360, bg = "white")
ggsave(file.path(final_dir, "figure_s12.png"), figure_s12, width = 7, height = 5, units = "in", dpi = 360, bg = "white")
record_timing(
  "render",
  render_started,
  "Rendered Figure 5 (7x9.14 in), Figure S11 (7x8.29 in), and Figure S12 (7x5 in) at 360 dpi."
)

timings <- bind_rows(
  timings,
  tibble(
    package_id = package_id,
    operation = "total",
    elapsed_seconds = proc.time()[["elapsed"]] - total_started,
    details = sprintf("Canonical in-memory rebuild using %d forked cores; persisted only script, three PNGs, and this timing table.", cores)
  )
)
write.table(
  timings,
  file.path(timing_dir, paste0(package_id, ".tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat(sprintf("FG6 complete in %.2f seconds\n", timings$elapsed_seconds[timings$operation == "total"]))
