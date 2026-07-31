#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}

repo <- normalizePath(arg_value("--project-root", "."), mustWork = TRUE)
release_root <- normalizePath(
  arg_value(
    "--release-root",
    file.path(repo, "data/modelling/gpath_v1/red_a30_counts_20260722")
  ),
  mustWork = TRUE
)
output_root <- arg_value(
  "--output-root",
  file.path(repo, "manuscript_figures/20260722T141236_resegmentation_fix")
)
output_root <- normalizePath(output_root, mustWork = TRUE)
dir.create(file.path(output_root, "final_images"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "timings"), recursive = TRUE, showWarnings = FALSE)

package_id <- "FG7_prediction_transfer"
timings <- list()
clock <- function() unname(proc.time()[["elapsed"]])
total_start <- clock()
record_timing <- function(operation, started, details) {
  timings[[length(timings) + 1L]] <<- data.frame(
    package_id = package_id,
    operation = operation,
    elapsed_seconds = round(clock() - started, 3),
    details = details,
    stringsAsFactors = FALSE
  )
}
assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

old_wd <- setwd(repo)
on.exit(setwd(old_wd), add = TRUE)
source("R/selection_strategy_utils.R")

assessment_path <- file.path(release_root, "derived/optimization/assessment.Rds")
dataset_manifest_path <- file.path(release_root, "manifests/datasets.tsv")
fit_plan_path <- file.path(release_root, "manifests/fit_plan.tsv")

t <- clock()
assessment <- readRDS(assessment_path)
dataset_manifest <- read.delim(
  dataset_manifest_path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE
)
fit_plan <- read.delim(
  fit_plan_path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(identical(assessment$metadata$release_id, "red_a30_counts_20260722"),
       "Unexpected optimization-assessment release id")
assert(all(c("fit_summary", "start_summary", "transfer_vs_null") %in% names(assessment)),
       "Optimization assessment does not have the expected schema")
record_timing(
  "load_assessment_and_manifests",
  t,
  paste0("Loaded ", assessment_path, ", ", dataset_manifest_path, ", and ", fit_plan_path)
)

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
line_ids <- 1:5
line_labels <- c(
  "1" = "MCF10A",
  "2" = "MDA-MB-231",
  "3" = "SNU668",
  "4" = "SUM-159-chem",
  "5" = "SUM-159-fuse"
)
directions <- c("low_to_high", "high_to_low")
search_modes <- c("common_r1", "model_native")
objectives <- c("favor_high", "favor_low")
r1_grid <- sort(unique(c(
  seq(0, 0.25, by = 0.005),
  seq(0.25, 1, by = 0.025),
  seq(1, 5, by = 0.1),
  seq(5, 25, by = 0.5)
)))
r2_grid <- seq(0, 1, by = 0.01)
common_r2 <- 1
viability_floor <- 0
robust_rate_margin <- 0.05 / 144

direction_rows <- dataset_manifest[
  dataset_manifest$parent_dataset == "all_lines" &
    dataset_manifest$transform == "ploidy_holdout" &
    dataset_manifest$direction %in% directions &
    dataset_manifest$line_name %in% unname(line_labels),
  , drop = FALSE
]
assert(nrow(direction_rows) == 10L,
       sprintf("Expected 10 all-line directional datasets, found %d", nrow(direction_rows)))
direction_rows$line_id <- as.integer(unname(
  setNames(as.character(line_ids), unname(line_labels))[direction_rows$line_name]
))
direction_rows <- direction_rows[
  order(match(direction_rows$line_id, line_ids), match(direction_rows$direction, directions)),
  , drop = FALSE
]

t <- clock()
wave2 <- assessment$transfer_vs_null
wave2 <- wave2[
  wave2$parent_dataset == "all_lines" &
    wave2$model_id %in% model_ids &
    wave2$direction %in% directions &
    wave2$line_name %in% unname(line_labels),
  , drop = FALSE
]
assert(nrow(wave2) == 50L && all(wave2$complete),
       "Wave-2 real/null assessment is incomplete for the selected 50 directional fits")
full_stan_row <- dataset_manifest[dataset_manifest$dataset_id == "all_lines", , drop = FALSE]
assert(nrow(full_stan_row) == 1L, "Could not resolve canonical all-lines Stan data")
full_stan_path <- full_stan_row$stan_data_path[[1]]
if (!grepl("^/", full_stan_path)) full_stan_path <- file.path(repo, full_stan_path)
stan_data <- readRDS(full_stan_path)
if (is.null(stan_data$group_id)) stan_data <- add_group_structure(stan_data)
assert(identical(sort(unique(as.integer(stan_data$line_id))), line_ids),
       "Canonical all-lines Stan data do not contain the expected five lines")
for (i in seq_len(nrow(direction_rows))) {
  path <- direction_rows$stan_data_path[[i]]
  if (!grepl("^/", path)) path <- file.path(repo, path)
  split_data <- readRDS(path)
  assert(!is.null(split_data$is_train), "Directional Stan data are missing is_train")
  target_id <- direction_rows$line_id[[i]]
  pair <- get_line_ploidy_pair(stan_data, target_id)
  held_out_value <- if (direction_rows$direction[[i]] == "low_to_high") pair$high else pair$low
  held_out_wells <- which(
    split_data$line_id == target_id &
      abs(split_data$ploidy_metric - held_out_value) < 1e-12
  )
  assert(length(held_out_wells) > 0L && all(split_data$is_train[held_out_wells] == 0L),
         "A directional Stan-data mask does not exclude its target state")
}
record_timing(
  "validate_wave2_directional_inputs",
  t,
  paste0(
    "Validated 50 complete selected-model wave-2 real/null assessment rows and ",
    "10 canonical directional Stan-data holdout masks"
  )
)

selected_fit_ids <- c(
  fit_plan$fit_id[
    fit_plan$dataset_id == "all_lines" & fit_plan$model_id %in% model_ids
  ],
  fit_plan$fit_id[
    fit_plan$dataset_id %in% direction_rows$dataset_id & fit_plan$model_id %in% model_ids
  ]
)
assert(length(selected_fit_ids) == 55L && !anyDuplicated(selected_fit_ids),
       "Expected 55 unique global/directional fit families")

load_best_rc0_draw <- function(fit_id) {
  plan_row <- fit_plan[fit_plan$fit_id == fit_id, , drop = FALSE]
  assert(nrow(plan_row) == 1L, paste("Could not resolve fit plan row", fit_id))
  starts <- assessment$start_summary[
    assessment$start_summary$fit_id == fit_id &
      assessment$start_summary$finite_lp &
      assessment$start_summary$rc == 0,
    , drop = FALSE
  ]
  assert(nrow(starts) > 0L, paste("No finite rc-zero starts for", fit_id))
  best_start <- starts$start_id[[which.max(starts$lp)]]
  draw_path <- file.path(plan_row$optim_dir[[1]], "optim_draws_all.Rds")
  if (!grepl("^/", draw_path)) draw_path <- file.path(repo, draw_path)
  draws <- readRDS(draw_path)
  assert(best_start <= length(draws), paste("Best start exceeds draw list for", fit_id))
  list(
    draw_vec = extract_draw_vector(draws[[best_start]]),
    best_start = best_start,
    best_lp = max(starts$lp)
  )
}

t <- clock()
fit_cache <- setNames(vector("list", length(selected_fit_ids)), selected_fit_ids)
for (fit_id in selected_fit_ids) fit_cache[[fit_id]] <- load_best_rc0_draw(fit_id)
record_timing(
  "select_and_load_rc0_fits",
  t,
  paste0(
    "Selected maximum-LP finite rc=0 starts from assessment$start_summary and ",
    "loaded the corresponding draw from each of 55 canonical optimizer archives"
  )
)

fit_for <- function(dataset_id, model_id) {
  row <- fit_plan[
    fit_plan$dataset_id == dataset_id & fit_plan$model_id == model_id,
    , drop = FALSE
  ]
  assert(nrow(row) == 1L, paste("Could not resolve fit", dataset_id, model_id))
  fit_cache[[row$fit_id[[1]]]]
}
rate_sign <- function(x, margin = 0) {
  ifelse(x > margin, "high", ifelse(x < -margin, "low", "neutral"))
}
requested_sign <- function(objective) {
  ifelse(objective == "favor_high", "high", "low")
}
make_environment_grid <- function(dims, search_mode) {
  if (search_mode == "common_r1" || dims$R == 1L) {
    return(data.frame(R1 = r1_grid, R2 = common_r2))
  }
  expand.grid(
    R1 = r1_grid, R2 = r2_grid,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
}
evaluate_growth_grid <- function(parms, environment_grid) {
  gpath_eval_instantaneous_net_growth_grid(
    parms = parms,
    env_grid = data.frame(
      glucose = environment_grid$R1,
      resource2 = environment_grid$R2
    )
  )
}
learn_environment <- function(
  environment_grid, complete_growth, target_observed_growth, objective
) {
  all_growth <- cbind(complete_growth$low, complete_growth$high)
  feasible <- apply(all_growth, 1L, min) >= viability_floor &
    target_observed_growth >= viability_floor
  delta <- complete_growth$high - complete_growth$low
  signed_delta <- if (objective == "favor_high") delta else -delta
  worst_effect <- apply(signed_delta, 1L, min)
  mean_effect <- rowMeans(signed_delta)
  min_growth <- pmin(apply(all_growth, 1L, min), target_observed_growth)
  eligible <- which(feasible)
  if (!length(eligible)) {
    return(list(status = "no_feasible_environment", best_index = NA_integer_))
  }
  ordering <- order(
    -worst_effect[eligible],
    -mean_effect[eligible],
    -min_growth[eligible],
    environment_grid$R1[eligible],
    environment_grid$R2[eligible]
  )
  best <- eligible[ordering[[1]]]
  list(
    status = "learned",
    best_index = best,
    worst_effect = worst_effect[[best]]
  )
}

t <- clock()
result_rows <- vector("list", 200L)
row_ptr <- 1L
for (model_id in model_ids) {
  dims <- parse_run_id(model_id)
  global_fit <- fit_for("all_lines", model_id)
  for (target_line in line_ids) {
    target_pair <- get_line_ploidy_pair(stan_data, target_line)
    global_target <- reconstruct_mixture_components(
      draw_vec = global_fit$draw_vec,
      model_id = model_id,
      line_id = target_line,
      low_ploidy = target_pair$low,
      high_ploidy = target_pair$high
    )
    for (direction in directions) {
      dataset_id <- direction_rows$dataset_id[
        direction_rows$line_id == target_line &
          direction_rows$direction == direction
      ]
      assert(length(dataset_id) == 1L, "Could not resolve directional transfer dataset")
      transfer_fit <- fit_for(dataset_id, model_id)
      transfer_components <- lapply(line_ids, function(line_id) {
        pair <- get_line_ploidy_pair(stan_data, line_id)
        reconstruct_mixture_components(
          draw_vec = transfer_fit$draw_vec,
          model_id = model_id,
          line_id = line_id,
          low_ploidy = pair$low,
          high_ploidy = pair$high
        )
      })
      names(transfer_components) <- as.character(line_ids)
      complete_lines <- setdiff(line_ids, target_line)
      observed_state <- if (direction == "low_to_high") "low" else "high"
      for (search_mode in search_modes) {
        environment_grid <- make_environment_grid(dims, search_mode)
        low_growth <- high_growth <- matrix(
          NA_real_, nrow = nrow(environment_grid), ncol = length(complete_lines)
        )
        for (j in seq_along(complete_lines)) {
          components <- transfer_components[[as.character(complete_lines[[j]])]]
          low_growth[, j] <- evaluate_growth_grid(
            components$parms_low, environment_grid
          )
          high_growth[, j] <- evaluate_growth_grid(
            components$parms_high, environment_grid
          )
        }
        target_transfer <- transfer_components[[as.character(target_line)]]
        target_observed_growth <- evaluate_growth_grid(
          if (observed_state == "low") {
            target_transfer$parms_low
          } else {
            target_transfer$parms_high
          },
          environment_grid
        )
        for (objective in objectives) {
          learned <- learn_environment(
            environment_grid,
            list(low = low_growth, high = high_growth),
            target_observed_growth,
            objective
          )
          assert(learned$status == "learned", "No feasible learned environment")
          best_env <- environment_grid[learned$best_index, , drop = FALSE]
          global_low <- evaluate_growth_grid(
            global_target$parms_low, best_env
          )[[1]]
          global_high <- evaluate_growth_grid(
            global_target$parms_high, best_env
          )[[1]]
          global_delta <- global_high - global_low
          wanted <- requested_sign(objective)
          result_rows[[row_ptr]] <- data.frame(
            model_id = model_id,
            model_alias = unname(model_alias[[model_id]]),
            n_resources = dims$R,
            target_line = target_line,
            target_label = unname(line_labels[[as.character(target_line)]]),
            direction = direction,
            search_mode = search_mode,
            objective = objective,
            training_worst_intended_delta = learned$worst_effect,
            global_delta = global_delta,
            end_to_end_objective_success =
              learned$worst_effect > 0 && rate_sign(global_delta) == wanted,
            end_to_end_robust_objective_success =
              learned$worst_effect > robust_rate_margin &&
              rate_sign(global_delta, robust_rate_margin) == wanted,
            stringsAsFactors = FALSE
          )
          row_ptr <- row_ptr + 1L
        }
      }
    }
  }
}
results <- do.call(rbind, result_rows)
assert(
  nrow(results) == 200L &&
    !anyDuplicated(results[c(
      "model_id", "target_line", "direction", "search_mode", "objective"
    )]),
  "Exact learned-environment recomputation did not yield 200 unique rows"
)
record_timing(
  "recompute_learned_environments",
  t,
  paste0(
    "Reconstructed 50 directional 9/10 fits, searched common-R1/model-native ",
    "R1/R2 grids in memory, and evaluated 200 frozen environments with 10/10 fits"
  )
)

summarise_groups <- function(dat, groups) {
  key <- interaction(dat[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(dat)), key), function(idx) {
    out <- dat[idx[[1]], groups, drop = FALSE]
    out$n_splits <- length(idx)
    out$end_to_end_n <- sum(dat$end_to_end_objective_success[idx])
    out$end_to_end_rate <- mean(dat$end_to_end_objective_success[idx])
    out$robust_end_to_end_n <- sum(dat$end_to_end_robust_objective_success[idx])
    out$robust_end_to_end_rate <- mean(dat$end_to_end_robust_objective_success[idx])
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

t <- clock()
headline <- summarise_groups(results, c("search_mode", "objective"))
models <- summarise_groups(
  results, c("search_mode", "objective", "model_id", "model_alias", "n_resources")
)
lines <- summarise_groups(
  results, c("search_mode", "objective", "target_line", "target_label")
)
assert(nrow(headline) == 4L && all(headline$n_splits == 50L),
       "Headline summaries do not have the expected four 50-split rows")
assert(nrow(models) == 20L && all(models$n_splits == 10L),
       "Model summaries do not have the expected 20 ten-split rows")
assert(nrow(lines) == 20L && all(lines$n_splits == 10L),
       "Line summaries do not have the expected 20 ten-split rows")
record_timing(
  "summarize_exact_and_robust_rates",
  t,
  "Computed headline, model-family, and target-line exact/robust requested-sign fractions"
)

objective_labels <- c(
  favor_high = "Favor high ploidy",
  favor_low = "Favor low ploidy"
)
mode_labels <- c(
  common_r1 = "Common R1",
  model_native = "Model-native"
)
objective_colors <- c(
  favor_high = "#D55E00",
  favor_low = "#0072B2"
)
base_theme <- theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "#D9D9D9"),
    strip.text = element_text(face = "bold", size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

make_design <- function() {
  boxes <- data.frame(
    x = c(1.25, 3.85, 6.45),
    y = 1,
    label = c(
      "Select\n9/10 fitted states",
      "Freeze\nR1/R2 environment",
      "Assess\n10/10 fitted reference"
    ),
    fill = c("#E6E6E6", "#F2F2F2", "#E8F1F8"),
    stringsAsFactors = FALSE
  )
  ggplot(boxes, aes(x, y)) +
    geom_tile(
      aes(fill = fill), width = 2.15, height = 0.94,
      color = "#333333", linewidth = 0.45
    ) +
    scale_fill_identity() +
    geom_text(aes(label = label), size = 2.6, lineheight = 0.92) +
    annotate(
      "segment", x = 2.36, xend = 2.70, y = 1, yend = 1,
      arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.55
    ) +
    annotate(
      "segment", x = 4.96, xend = 5.30, y = 1, yend = 1,
      arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.55
    ) +
    annotate(
      "text", x = 3.85, y = 0.33,
      label = paste0(
        "Held-out state excluded from selection; requested sign compared with ",
        "target-inclusive 10/10 fitted reference"
      ),
      size = 2.35, color = "#444444"
    ) +
    coord_cartesian(xlim = c(0.05, 7.65), ylim = c(0.08, 1.58), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(8, 8, 8, 8))
}

make_rates <- function() {
  exact <- data.frame(
    headline[c("search_mode", "objective", "n_splits", "end_to_end_n", "end_to_end_rate")],
    criterion = "Exact sign",
    stringsAsFactors = FALSE
  )
  names(exact)[names(exact) == "end_to_end_n"] <- "n"
  names(exact)[names(exact) == "end_to_end_rate"] <- "rate"
  robust <- data.frame(
    headline[c(
      "search_mode", "objective", "n_splits",
      "robust_end_to_end_n", "robust_end_to_end_rate"
    )],
    criterion = "Robust sign",
    stringsAsFactors = FALSE
  )
  names(robust)[names(robust) == "robust_end_to_end_n"] <- "n"
  names(robust)[names(robust) == "robust_end_to_end_rate"] <- "rate"
  dat <- rbind(exact, robust)
  dat$objective_label <- unname(objective_labels[dat$objective])
  dat$search_label <- unname(mode_labels[dat$search_mode])
  dat$row <- factor(
    paste(dat$search_label, dat$objective_label, sep = " \u00b7 "),
    levels = rev(c(
      "Common R1 \u00b7 Favor high ploidy",
      "Common R1 \u00b7 Favor low ploidy",
      "Model-native \u00b7 Favor high ploidy",
      "Model-native \u00b7 Favor low ploidy"
    ))
  )
  dat$label <- sprintf("%d/%d", dat$n, dat$n_splits)
  for (key in unique(paste(dat$search_mode, dat$objective))) {
    idx <- which(paste(dat$search_mode, dat$objective) == key)
    if (length(idx) == 2L && length(unique(dat$rate[idx])) == 1L) {
      dat$label[idx[dat$criterion[idx] == "Robust sign"]] <- ""
    }
  }
  ggplot(dat, aes(rate, row, color = objective, shape = criterion)) +
    geom_point(
      aes(group = criterion), size = 3.0, stroke = 0.85, fill = "white",
      position = position_dodge(width = 0.42)
    ) +
    geom_text(
      aes(label = label, group = criterion), hjust = -0.18, size = 2.55,
      position = position_dodge(width = 0.42), show.legend = FALSE
    ) +
    scale_color_manual(values = objective_colors, guide = "none") +
    scale_shape_manual(values = c("Exact sign" = 16, "Robust sign" = 1)) +
    scale_x_continuous(
      limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
      labels = function(x) paste0(100 * x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Requested-sign fraction across 50 fitted splits",
      y = NULL, shape = NULL
    ) +
    base_theme +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5, 18, 5, 5)
    )
}

make_dependence <- function() {
  model_dat <- data.frame(
    source_type = "Model family",
    category = models$model_alias,
    search_mode = models$search_mode,
    objective = models$objective,
    rate = models$end_to_end_rate,
    n = models$end_to_end_n,
    denom = models$n_splits,
    stringsAsFactors = FALSE
  )
  line_dat <- data.frame(
    source_type = "Target line",
    category = lines$target_label,
    search_mode = lines$search_mode,
    objective = lines$objective,
    rate = lines$end_to_end_rate,
    n = lines$end_to_end_n,
    denom = lines$n_splits,
    stringsAsFactors = FALSE
  )
  dat <- rbind(model_dat, line_dat)
  dat$search_mode <- factor(
    unname(mode_labels[dat$search_mode]),
    levels = unname(mode_labels)
  )
  dat$objective <- factor(
    unname(objective_labels[dat$objective]),
    levels = unname(objective_labels)
  )
  dat$facet_label <- factor(
    paste(
      ifelse(dat$source_type == "Model family", "Models", "Lines"),
      dat$objective,
      sep = " \u00b7 "
    ),
    levels = c(
      "Models \u00b7 Favor high ploidy",
      "Models \u00b7 Favor low ploidy",
      "Lines \u00b7 Favor high ploidy",
      "Lines \u00b7 Favor low ploidy"
    )
  )
  dat$category <- factor(
    dat$category,
    levels = rev(c(
      "1R", "1R,W(m)", "2R(a)", "2R(f)", "2R(f),W(m)",
      "MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse"
    ))
  )
  dat$count <- sprintf("%d/%d", dat$n, dat$denom)
  ggplot(dat, aes(1, category, fill = rate)) +
    geom_tile(width = 0.96, color = "white", linewidth = 0.45) +
    geom_text(aes(label = count), size = 2.55, color = "#111111") +
    facet_grid(facet_label ~ search_mode, scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = c("#F7FBFF", "#6BAED6", "#08306B"),
      limits = c(0, 1), breaks = c(0, 0.5, 1),
      labels = c("0%", "50%", "100%")
    ) +
    scale_x_continuous(limits = c(0.5, 1.5), expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL, fill = "Requested-sign\nfraction") +
    base_theme +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      strip.text.y = element_text(angle = 0, size = 7.2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      plot.margin = margin(4, 4, 4, 4)
    )
}

t <- clock()
plots <- list(a = make_design(), b = make_rates(), c = make_dependence())
record_timing(
  "construct_panels",
  t,
  "Constructed the accepted S14 schematic, exact/robust summary plot, and model/line heatmap in memory"
)

t <- clock()
final_path <- file.path(output_root, "final_images", "figure_s14.png")
temporary_png <- tempfile("figure_s14_", tmpdir = tempdir(), fileext = ".png")
ragg::agg_png(
  temporary_png, width = 7, height = 8.76, units = "in", res = 300,
  background = "white"
)
grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))
layout <- data.frame(
  panel = c("a", "b", "c"),
  x = c(0, 0, 0),
  y = c(0.805936073059361, 0.534246575342466, 0),
  width = c(1, 1, 1),
  height = c(0.194063926940639, 0.262557077625571, 0.525114155251142)
)
for (i in seq_len(nrow(layout))) {
  panel <- layout$panel[[i]]
  pushViewport(viewport(
    x = layout$x[[i]] + layout$width[[i]] / 2,
    y = layout$y[[i]] + layout$height[[i]] / 2,
    width = layout$width[[i]],
    height = layout$height[[i]],
    just = "center"
  ))
  grid.draw(ggplotGrob(plots[[panel]]))
  popViewport()
  grid.text(
    panel,
    x = unit(layout$x[[i]] + 0.006, "npc"),
    y = unit(layout$y[[i]] + layout$height[[i]] - 0.006, "npc"),
    just = c("left", "top"),
    gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12)
  )
}
invisible(dev.off())
assert(file.copy(temporary_png, final_path, overwrite = TRUE),
       "Could not promote rendered Figure S14")
unlink(temporary_png)
record_timing(
  "render",
  t,
  paste0("Rendered ", final_path, " at 7 x 8.76 inches and 300 dpi")
)

timings[[length(timings) + 1L]] <- data.frame(
  package_id = package_id,
  operation = "total",
  elapsed_seconds = round(clock() - total_start, 3),
  details = paste0(
    "Exact 9/10 learned-environment requested-sign recomputation from canonical ",
    "rc-zero optimizer draws; no legacy transfer tables or persisted intermediates"
  ),
  stringsAsFactors = FALSE
)
timing_path <- file.path(output_root, "timings", paste0(package_id, ".tsv"))
write.table(
  do.call(rbind, timings), timing_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf("Wrote %s\n", final_path))
cat(sprintf("Wrote %s\n", timing_path))
