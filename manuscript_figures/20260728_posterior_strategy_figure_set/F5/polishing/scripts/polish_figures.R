#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(ggplot2)
  library(posterior)
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
library_only <- tolower(Sys.getenv("F5_POLISH_LIBRARY_ONLY", "false")) %in%
  c("true", "t", "1", "yes")
if (!library_only && !phase %in% c("subpanels", "final")) {
  stop("--phase must be subpanels or final", call. = FALSE)
}

project_root <- normalizePath(".", mustWork = TRUE)
polish_root <- file.path(
  project_root, "manuscript_figures",
  "20260728_posterior_strategy_figure_set", "F5", "polishing"
)
subpanel_dir <- file.path(polish_root, "subpanels")
layout_dir <- file.path(polish_root, "layout")
final_dir <- file.path(polish_root, "final_images")
derived_dir <- file.path(polish_root, "derived")
for (path in c(subpanel_dir, layout_dir, final_dir, derived_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

release_root <- file.path(
  project_root, "data", "modelling", "gpath_v1",
  "red_a30_counts_20260722"
)
simulation_run_root <- file.path(
  release_root, "derived", "posterior", "strategy_simulations", "runs",
  "seed3density_7g_xxx_xcc_axax_20260727"
)
simulation_root <- file.path(
  simulation_run_root, "loo_exclude_sum_159_fuse"
)
prediction_root <- file.path(
  release_root, "derived", "posterior", "predictions",
  "loo_exclude_sum_159_fuse"
)
cache_path <- file.path(derived_dir, "figure_5_panel_data.qs")

source(file.path(project_root, "scripts", "parameter_estimation", "common.R"))
source(file.path(project_root, "R", "project_paths.R"))
release_cfg <- pe_read_config(file.path(release_root, "release.json"))
source(get_model_r_path(release_cfg$model_name, release_cfg$model_version))
source(file.path(project_root, "R", "gpath_run_utils.R"))
source(file.path(project_root, "R", "gpath_derived_utils.R"))

models <- c(
  "1R_1P_0W_C0_M1" = "1R–1P",
  "1R_1P_1W_C0_M0" = "1R–1P–W",
  "2R_1P_0W_C0_M1" = "2R–1P",
  "2R_2P_0W_C0_M1" = "2R–2P",
  "2R_2P_1W_C0_M0" = "2R–2P–W"
)
lines <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem")
line_slugs <- c(
  "MCF10A" = "mcf10a",
  "MDA-MB-231" = "mda_mb_231",
  "SNU668" = "snu668",
  "SUM-159-chem" = "sum_159_chem"
)
strategies <- c(
  "0CC" = "0,C,C",
  "5CC" = "5,C,C",
  "111" = "1,1,1"
)
omitted_support_strategy <- "0,A0,A0"
density_labels <- c("0.5" = "0.5×", "1" = "1×", "2" = "2×")
density_colors <- c(
  "0.5×" = "#332288",
  "1×" = "#EE7733",
  "2×" = "#CCBB44"
)
strategy_colors <- c(
  "0CC" = "#CC3311",
  "5CC" = "#0077BB",
  "111" = "#009988"
)
ploidy_linetypes <- c(
  "Low ploidy" = "solid",
  "High ploidy" = "22"
)
day_occasions <- c(
  "Day 2" = "day2_pre_action",
  "Day 4" = "day4_pre_action",
  "Day 6" = "day6_end"
)
day_colors <- c(
  "Day 2" = "#AA3377",
  "Day 4" = "#33BBEE",
  "Day 6" = "#666666"
)
model_colors <- c(
  "1R–1P" = "#332288",
  "1R–1P–W" = "#AA3377",
  "2R–1P" = "#EE7733",
  "2R–2P" = "#CCBB44",
  "2R–2P–W" = "#33BBEE"
)
line_colors <- c(
  "MCF10A" = "#332288",
  "MDA-MB-231" = "#AA3377",
  "SNU668" = "#EE7733",
  "SUM-159-chem" = "#CCBB44"
)

surface_x0 <- 1e-4
surface_glucose_max <- 25
surface_glucose_points <- 120L
surface_resource2_points <- 40L
surface_draws <- 100L
surface_cores <- 4L
surface_glucose_pad <- -1e-4
surface_resource2_pad <- -0.025
surface_transform <- function(x) asinh(x / surface_x0)

validate_inputs <- function() {
  validation_path <- file.path(simulation_run_root, "validation.json")
  validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
  if (!isTRUE(validation$validation_passed)) {
    stop("The compact posterior strategy run is not validated.", call. = FALSE)
  }
  prediction_tasks <- fread(
    file.path(release_root, "derived", "plans", "prediction_tasks.tsv")
  )[
    dataset_id == "loo_exclude_sum_159_fuse" &
      model_id %in% names(models)
  ]
  if (
    nrow(prediction_tasks) != length(models) ||
    length(unique(prediction_tasks$stan_data_path)) != 1L
  ) {
    stop("Could not resolve the five four-line posterior fits.", call. = FALSE)
  }
  prediction_tasks
}

compact_strategy_label <- function(x) {
  first <- sub(",.*$", "", x)
  family <- fifelse(
    grepl(",C,C$", x), "CC",
    fifelse(grepl(",A[^,]+,A[^,]+$", x), "AA", "XX")
  )
  paste0(first, family)
}

build_support_summaries <- function() {
  tasks <- fread(file.path(simulation_run_root, "plans", "tasks.tsv"))[
    dataset_id == "loo_exclude_sum_159_fuse"
  ]
  schedules <- fread(file.path(simulation_run_root, "schedule_index.tsv"))[
    strategy_code != omitted_support_strategy
  ]
  if (nrow(tasks) != 20L || nrow(schedules) != 20L) {
    stop("Unexpected four-line support input coverage.", call. = FALSE)
  }

  strategy_meta <- copy(schedules)
  strategy_meta[, family_order := fcase(
    strategy_family == "XXX", 1L,
    strategy_family == "XA_xA_x", 2L,
    strategy_family == "XCC", 3L,
    default = NA_integer_
  )]
  strategy_meta[, family_label := fcase(
    strategy_family == "XXX", "XX",
    strategy_family == "XA_xA_x", "AA",
    strategy_family == "XCC", "CC",
    default = NA_character_
  )]
  if (anyNA(strategy_meta$family_order)) {
    stop("Unknown strategy family in schedule index.", call. = FALSE)
  }
  strategy_meta[, glucose_label := fifelse(
    abs(g0_day0 - round(g0_day0)) < 1e-10,
    sprintf("%.0f", g0_day0),
    sub("0+$", "", sub("\\.$", "", sprintf("%.2f", g0_day0)))
  )]
  strategy_meta[, strategy_label := paste0(glucose_label, family_label)]
  setorder(strategy_meta, g0_day0, family_order)
  strategy_codes <- strategy_meta$strategy_code
  strategy_order_top_to_bottom <- strategy_meta$strategy_label
  strategy_label_map <- setNames(
    strategy_meta$strategy_label, strategy_meta$strategy_code
  )

  parts <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    sim <- qread(
      file.path(project_root, tasks$output_path[[i]]),
      nthreads = 2L
    )
    density_i <- match(
      c(0.5, 1, 2), as.numeric(sim$density_multipliers)
    )
    day_i <- match(
      unname(day_occasions), sim$observation_occasions
    )
    low_i <- match("N_low", sim$state_names)
    high_i <- match("N_high", sim$state_names)
    strategy_i <- match(
      strategy_codes, sim$schedule_index$strategy_code
    )
    if (anyNA(c(density_i, day_i, low_i, high_i, strategy_i))) {
      stop("Required support indices are missing.", call. = FALSE)
    }
    low <- sim$observation_state[
      , density_i, strategy_i, day_i, low_i, drop = FALSE
    ]
    high <- sim$observation_state[
      , density_i, strategy_i, day_i, high_i, drop = FALSE
    ]
    dim(low) <- dim(high) <- c(
      dim(sim$observation_state)[[1L]],
      length(density_i),
      length(strategy_i),
      length(day_i)
    )
    grid <- as.data.table(expand.grid(
      draw = seq_len(dim(low)[[1L]]),
      density_multiplier = c(0.5, 1, 2),
      strategy_code = strategy_codes,
      day = names(day_occasions),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
    grid[, `:=`(
      supports_high = as.vector(high) > as.vector(low),
      model_id = tasks$model_id[[i]],
      line_name = tasks$line_name[[i]]
    )]
    parts[[i]] <- grid
  }
  draws <- rbindlist(parts)

  support_by_density <- draws[day == "Day 6", .(
    supporting_draws = sum(supports_high),
    total_draws = .N
  ), by = .(strategy_code, density_multiplier)]
  support_by_day <- draws[density_multiplier == 1, .(
    supporting_draws = sum(supports_high),
    total_draws = .N
  ), by = .(strategy_code, day)]
  support_by_model <- draws[
    density_multiplier == 1 & day == "Day 6", .(
      supporting_draws = sum(supports_high),
      total_draws = .N
    ), by = .(strategy_code, model_id)
  ]
  support_by_line <- draws[
    density_multiplier == 1 & day == "Day 6", .(
      supporting_draws = sum(supports_high),
      total_draws = .N
    ), by = .(strategy_code, line_name)
  ]

  stopifnot(
    nrow(support_by_density) == 60L,
    nrow(support_by_day) == 60L,
    nrow(support_by_model) == 100L,
    nrow(support_by_line) == 80L,
    all(support_by_density$total_draws == 2000L),
    all(support_by_day$total_draws == 2000L),
    all(support_by_model$total_draws == 400L),
    all(support_by_line$total_draws == 500L)
  )

  add_display_fields <- function(x, group_column, group_levels) {
    out <- copy(x)
    out[, strategy_label := unname(
      strategy_label_map[strategy_code]
    )]
    out[, strategy_factor := factor(
      strategy_label,
      levels = rev(strategy_order_top_to_bottom)
    )]
    out[, group := factor(
      get(group_column), levels = group_levels
    )]
    offsets <- setNames(
      seq(-0.24, 0.24, length.out = length(group_levels)),
      group_levels
    )
    out[, strategy_y := as.numeric(strategy_factor) +
      unname(offsets[as.character(group)])]
    out
  }

  support_by_density[, density := unname(
    density_labels[as.character(density_multiplier)]
  )]
  support_by_density <- add_display_fields(
    support_by_density, "density", unname(density_labels)
  )
  support_by_day <- add_display_fields(
    support_by_day, "day", names(day_occasions)
  )
  support_by_model[, model := unname(models[model_id])]
  support_by_model <- add_display_fields(
    support_by_model, "model", unname(models)
  )
  support_by_line[, cell_line := line_name]
  support_by_line <- add_display_fields(
    support_by_line, "cell_line", lines
  )

  finish_summary <- function(points) {
    list(
      points = points,
      ranges = points[, .(
        minimum_support = min(supporting_draws),
        maximum_support = max(supporting_draws),
        strategy_y = unique(as.numeric(strategy_factor))
      ), by = strategy_factor],
      labels = rev(strategy_order_top_to_bottom)
    )
  }

  list(
    density = finish_summary(support_by_density),
    day = finish_summary(support_by_day),
    model = finish_summary(support_by_model),
    line = finish_summary(support_by_line),
    strategy_order_top_to_bottom = strategy_order_top_to_bottom
  )
}

make_dense_environment <- function(
  sim, density_i, schedule_i, strategy_labels,
  glucose_i, r2_i, time_hours
) {
  nstrategy <- length(schedule_i)
  ntime <- length(time_hours)
  q50_i <- match("q50", dimnames(sim$trajectory_quantiles)$quantile)
  glucose <- sim$trajectory_quantiles[
    q50_i, density_i, schedule_i, , glucose_i, drop = FALSE
  ]
  dim(glucose) <- c(nstrategy, ntime)
  resource2 <- if (!is.na(r2_i)) {
    out <- sim$trajectory_quantiles[
      q50_i, density_i, schedule_i, , r2_i, drop = FALSE
    ]
    dim(out) <- c(nstrategy, ntime)
    out
  } else {
    matrix(1, nrow = nstrategy, ncol = ntime)
  }

  dense <- as.data.table(expand.grid(
    strategy_i = seq_len(nstrategy),
    time_i = seq_len(ntime),
    KEEP.OUT.ATTRS = FALSE
  ))
  dense[, `:=`(
    strategy = strategy_labels[strategy_i],
    time_hours = time_hours[time_i],
    point_order = 1L,
    glucose = as.vector(glucose),
    resource2 = as.vector(resource2)
  )]

  post_occasions <- c(day2_post_action = 48, day4_post_action = 96)
  post <- rbindlist(lapply(names(post_occasions), function(occasion_name) {
    occasion_i <- match(occasion_name, sim$observation_occasions)
    time_value <- unname(post_occasions[[occasion_name]])
    rbindlist(lapply(seq_len(nstrategy), function(strategy_j) {
      data.table(
        strategy_i = strategy_j,
        time_i = NA_integer_,
        strategy = strategy_labels[[strategy_j]],
        time_hours = time_value,
        point_order = 2L,
        glucose = median(sim$observation_state[
          , density_i, schedule_i[[strategy_j]], occasion_i, glucose_i
        ]),
        resource2 = if (!is.na(r2_i)) {
          median(sim$observation_state[
            , density_i, schedule_i[[strategy_j]], occasion_i, r2_i
          ])
        } else {
          1
        }
      )
    }))
  }))
  rbind(dense, post, use.names = TRUE)[
    order(strategy_i, time_hours, point_order)
  ]
}

build_panel_data <- function() {
  prediction_tasks <- validate_inputs()
  stan_data <- add_group_structure(
    readRDS(unique(prediction_tasks$stan_data_path)[[1L]])
  )
  line_table <- as.data.table(gpd_line_table(stan_data))

  surface_parts <- list()
  path_parts <- list()
  action_parts <- list()
  marker_parts <- list()
  count_parts <- list()
  growth_parts <- list()

  for (model_id in names(models)) {
    current_model_id <- model_id
    model_label <- unname(models[[model_id]])
    pred <- readRDS(file.path(prediction_root, paste0(model_id, ".Rds")))
    prediction_lines <- unique(as.character(pred$well_metadata$line_name))

    fit_task <- prediction_tasks[model_id == current_model_id]
    surface_draw_info <- gpd_read_fit_draws(
      fit_task$nuts_dir[[1L]],
      fit_task$nuts_chains[[1L]],
      max_draws = 400L
    )
    surface_idx <- unique(as.integer(round(seq(
      1,
      nrow(surface_draw_info$draws),
      length.out = min(surface_draws, nrow(surface_draw_info$draws))
    ))))
    selected_surface_index <- surface_draw_info$index[
      surface_idx, , drop = FALSE
    ]
    production_surface_index <- pred$growth_surface$draw_index
    index_columns <- intersect(
      c("chain_id", "iteration", "draw_id"),
      intersect(
        names(selected_surface_index),
        names(production_surface_index)
      )
    )
    if (
      !length(index_columns) ||
      !identical(
        as.data.frame(
          selected_surface_index[, index_columns, drop = FALSE]
        ),
        as.data.frame(
          production_surface_index[, index_columns, drop = FALSE]
        )
      )
    ) {
      stop(
        "Regenerated surface draw selection differs from production for ",
        model_id,
        call. = FALSE
      )
    }

    glucose_axis <- seq(
      0,
      surface_transform(surface_glucose_max),
      length.out = surface_glucose_points
    )
    glucose_grid <- surface_x0 * sinh(glucose_axis)
    resource2_grid <- seq(
      0, 1, length.out = surface_resource2_points
    )
    evaluate_surface_draw <- function(i) {
      gpd_growth_delta_surface(
        draw_vec = surface_draw_info$draws[
          surface_idx[[i]], , drop = TRUE
        ],
        stan_data = stan_data,
        model_id = model_id,
        glucose_grid = glucose_grid,
        resource2_grid = resource2_grid,
        extra_resource_value = 1
      )
    }
    surface_parts_by_draw <- if (surface_cores > 1L) {
      parallel::mclapply(
        seq_along(surface_idx),
        evaluate_surface_draw,
        mc.cores = surface_cores
      )
    } else {
      lapply(seq_along(surface_idx), evaluate_surface_draw)
    }
    surface_array <- gpd_bind_first_dimension(surface_parts_by_draw)
    surface_median <- apply(
      surface_array, c(2, 3, 4), median, na.rm = TRUE
    )
    surface_lines <- as.character(line_table$line_name)
    if (!setequal(surface_lines, prediction_lines)) {
      stop(
        "Regenerated surface line coverage differs from production for ",
        model_id,
        call. = FALSE
      )
    }

    draw_info <- gpd_read_fit_draws(
      fit_task$nuts_dir[[1L]],
      fit_task$nuts_chains[[1L]],
      max_draws = 100L
    )

    for (line_name in lines) {
      current_line_name <- line_name
      line_i <- match(line_name, surface_lines)
      line_id <- line_table[line_name == current_line_name, line_id]
      if (is.na(line_i) || length(line_id) != 1L) {
        stop(
          "Could not resolve ", line_name, " for ", model_id,
          call. = FALSE
        )
      }
      surface_parts[[length(surface_parts) + 1L]] <- as.data.table(
        expand.grid(
          glucose = glucose_grid,
          R2 = resource2_grid,
          KEEP.OUT.ATTRS = FALSE
        )
      )[, `:=`(
        delta_growth = as.vector(surface_median[line_i, , ]),
        model = model_label,
        line_name = line_name
      )]

      simulation_path <- file.path(
        simulation_root, model_id, paste0(line_slugs[[line_name]], ".qs")
      )
      sim <- qread(simulation_path, nthreads = 2L)
      if (
        !identical(
          as.integer(sim$draw_index$chain_id),
          as.integer(draw_info$index$chain_id)
        ) ||
        !identical(
          as.integer(sim$draw_index$iteration),
          as.integer(draw_info$index$iteration)
        )
      ) {
        stop(
          "Posterior draw mismatch for ", model_id, " / ", line_name,
          call. = FALSE
        )
      }

      density_i <- match(1, as.numeric(sim$density_multipliers))
      q50_i <- match("q50", dimnames(sim$trajectory_quantiles)$quantile)
      schedule_i <- match(
        unname(strategies), sim$schedule_index$strategy_code
      )
      glucose_i <- match("G1", sim$state_names)
      r2_i <- match("R2", sim$state_names)
      low_i <- match("N_low", sim$state_names)
      high_i <- match("N_high", sim$state_names)
      if (anyNA(c(
        density_i, q50_i, schedule_i, glucose_i, low_i, high_i
      ))) {
        stop("Required simulation fields are absent: ", simulation_path,
             call. = FALSE)
      }

      time_hours <- as.numeric(sim$time_hours)
      ntime <- length(time_hours)
      nstrategy <- length(strategies)
      low_n <- sim$trajectory_quantiles[
        q50_i, density_i, schedule_i, , low_i, drop = FALSE
      ]
      high_n <- sim$trajectory_quantiles[
        q50_i, density_i, schedule_i, , high_i, drop = FALSE
      ]
      dim(low_n) <- dim(high_n) <- c(nstrategy, ntime)

      count_grid <- as.data.table(expand.grid(
        strategy_i = seq_len(nstrategy),
        time_i = seq_len(ntime),
        KEEP.OUT.ATTRS = FALSE
      ))
      count_grid[, `:=`(
        strategy = names(strategies)[strategy_i],
        time_hours = time_hours[time_i],
        low_ploidy = as.vector(low_n),
        high_ploidy = as.vector(high_n),
        model = model_label,
        line_name = line_name
      )]
      count_parts[[length(count_parts) + 1L]] <- melt(
        count_grid,
        id.vars = c("strategy", "time_hours", "model", "line_name"),
        measure.vars = c("low_ploidy", "high_ploidy"),
        variable.name = "ploidy",
        value.name = "live_cells"
      )[, ploidy := fifelse(
        ploidy == "low_ploidy", "Low ploidy", "High ploidy"
      )]

      has_r2 <- !is.na(r2_i)
      for (strategy_j in seq_along(strategies)) {
        strategy_label <- names(strategies)[[strategy_j]]
        strategy_index <- schedule_i[[strategy_j]]
        reference_y <- strategy_j / (length(strategies) + 1)
        glucose <- as.numeric(sim$trajectory_quantiles[
          q50_i, density_i, strategy_index, , glucose_i
        ])
        r2 <- if (has_r2) {
          as.numeric(sim$trajectory_quantiles[
            q50_i, density_i, strategy_index, , r2_i
          ])
        } else {
          rep(reference_y, ntime)
        }
        dense_path <- data.table(
          time_hours = time_hours,
          glucose = glucose,
          R2 = r2,
          segment = cut(
            time_hours,
            breaks = c(-Inf, 48, 96, Inf),
            labels = c("0–2 d", "2–4 d", "4–6 d"),
            right = TRUE
          ),
          model = model_label,
          line_name = line_name,
          strategy = strategy_label
        )

        occasions <- c(
          day0_post_seed = 0,
          day2_pre_action = 48,
          day2_post_action = 48,
          day4_pre_action = 96,
          day4_post_action = 96,
          day6_end = 144
        )
        snapshots <- rbindlist(lapply(names(occasions), function(occasion) {
          occasion_i <- match(occasion, sim$observation_occasions)
          data.table(
            occasion = occasion,
            time_hours = unname(occasions[[occasion]]),
            glucose = median(sim$observation_state[
              , density_i, strategy_index, occasion_i, glucose_i
            ]),
            R2 = if (has_r2) {
              median(sim$observation_state[
                , density_i, strategy_index, occasion_i, r2_i
              ])
            } else {
              reference_y
            },
            model = model_label,
            line_name = line_name,
            strategy = strategy_label
          )
        }))
        post_path <- rbind(
          snapshots[occasion == "day2_post_action"][
            , segment := factor("2–4 d", levels = levels(dense_path$segment))
          ],
          snapshots[occasion == "day4_post_action"][
            , segment := factor("4–6 d", levels = levels(dense_path$segment))
          ],
          fill = TRUE
        )[, .(
          time_hours, glucose, R2, segment,
          model, line_name, strategy
        )]
        path_parts[[length(path_parts) + 1L]] <- rbind(
          dense_path, post_path, use.names = TRUE
        )[order(segment, time_hours)]

        action_parts[[length(action_parts) + 1L]] <- rbind(
          snapshots[
            occasion %in% c("day2_pre_action", "day2_post_action")
          ][, action := "day2"],
          snapshots[
            occasion %in% c("day4_pre_action", "day4_post_action")
          ][, action := "day4"],
          use.names = TRUE
        )
        marker_parts[[length(marker_parts) + 1L]] <- snapshots
      }

      environment <- make_dense_environment(
        sim = sim,
        density_i = density_i,
        schedule_i = schedule_i,
        strategy_labels = names(strategies),
        glucose_i = glucose_i,
        r2_i = r2_i,
        time_hours = time_hours
      )
      ploidy_states <- get_line_ploidy_states(stan_data, line_id)
      environment_frame <- data.frame(
        glucose = pmax(environment$glucose, 0),
        resource2 = pmax(environment$resource2, 0)
      )
      low_growth <- matrix(
        NA_real_, nrow = nrow(environment), ncol = nrow(draw_info$draws)
      )
      high_growth <- low_growth
      for (draw_i in seq_len(nrow(draw_info$draws))) {
        draw_vec <- draw_info$draws[draw_i, , drop = TRUE]
        low_state <- gpath_reconstruct_state_from_draw(
          draw_vec, model_id, line_id, ploidy_states$low_value
        )
        high_state <- gpath_reconstruct_state_from_draw(
          draw_vec, model_id, line_id, ploidy_states$high_value
        )
        low_growth[, draw_i] <- gpath_eval_instantaneous_net_growth_grid(
          low_state, environment_frame
        )
        high_growth[, draw_i] <- gpath_eval_instantaneous_net_growth_grid(
          high_state, environment_frame
        )
      }
      environment[, `:=`(
        low_growth = apply(low_growth, 1L, median, na.rm = TRUE),
        high_growth = apply(high_growth, 1L, median, na.rm = TRUE),
        model = model_label,
        line_name = line_name
      )]
      growth_parts[[length(growth_parts) + 1L]] <- melt(
        environment,
        id.vars = c(
          "strategy", "time_hours", "point_order", "model", "line_name"
        ),
        measure.vars = c("low_growth", "high_growth"),
        variable.name = "ploidy",
        value.name = "net_growth"
      )[, ploidy := fifelse(
        ploidy == "low_growth", "Low ploidy", "High ploidy"
      )]
    }
  }

  bundle <- list(
    support = build_support_summaries(),
    surface = rbindlist(surface_parts),
    paths = rbindlist(path_parts),
    actions = rbindlist(action_parts, fill = TRUE),
    markers = rbindlist(marker_parts, fill = TRUE),
    counts = rbindlist(count_parts),
    growth = rbindlist(growth_parts),
    metadata = list(
      dataset_id = "loo_exclude_sum_159_fuse",
      density_multiplier = 1,
      strategies = strategies,
      models = models,
      lines = lines,
      posterior_draws = 100L,
      surface_transform = "asinh(glucose / 1e-4)",
      surface_glucose_points = surface_glucose_points,
      surface_resource2_points = surface_resource2_points,
      surface_glucose_pad = surface_glucose_pad,
      surface_resource2_pad = surface_resource2_pad
    )
  )
  qsave(bundle, cache_path, preset = "high", nthreads = 2L)
  bundle
}

set_plot_factors <- function(bundle) {
  for (name in c("surface", "paths", "actions", "markers", "counts", "growth")) {
    x <- bundle[[name]]
    x[, `:=`(
      model = factor(model, levels = unname(models)),
      line_name = factor(line_name, levels = lines)
    )]
    if ("strategy" %in% names(x)) {
      x[, strategy := factor(strategy, levels = names(strategies))]
    }
    if ("ploidy" %in% names(x)) {
      x[, ploidy := factor(ploidy, levels = names(ploidy_linetypes))]
    }
    if ("time_hours" %in% names(x)) x[, day := time_hours / 24]
    bundle[[name]] <- x
  }
  bundle$surface[, glucose_axis := surface_transform(glucose)]
  bundle$paths[, glucose_axis := surface_transform(glucose)]
  bundle$actions[, glucose_axis := surface_transform(glucose)]
  bundle$markers[, glucose_axis := surface_transform(glucose)]

  centres_to_edges <- function(x, lower, upper) {
    x <- sort(unique(x))
    inner <- (x[-1L] + x[-length(x)]) / 2
    c(lower, inner, upper)
  }
  glucose_centres <- sort(unique(bundle$surface$glucose_axis))
  resource_centres <- sort(unique(bundle$surface$R2))
  if (
    min(bundle$surface$glucose) != 0 ||
    min(bundle$surface$R2) != 0
  ) {
    stop("The regenerated surface is not zero-inclusive.", call. = FALSE)
  }
  glucose_edges <- centres_to_edges(
    glucose_centres,
    lower = surface_transform(surface_glucose_pad),
    upper = surface_transform(surface_glucose_max)
  )
  resource_edges <- centres_to_edges(
    resource_centres,
    lower = surface_resource2_pad,
    upper = 1
  )
  bundle$surface[, `:=`(
    tile_xmin = glucose_edges[
      match(glucose_axis, glucose_centres)
    ],
    tile_xmax = glucose_edges[
      match(glucose_axis, glucose_centres) + 1L
    ],
    tile_ymin = resource_edges[match(R2, resource_centres)],
    tile_ymax = resource_edges[match(R2, resource_centres) + 1L]
  )]
  bundle$surface[, `:=`(
    tile_x = (tile_xmin + tile_xmax) / 2,
    tile_y = (tile_ymin + tile_ymax) / 2,
    tile_width = tile_xmax - tile_xmin,
    tile_height = tile_ymax - tile_ymin
  )]
  bundle
}

theme_f5 <- function() {
  theme_bw(base_size = 8, base_family = "sans") +
    theme(
      text = element_text(colour = "#202020"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "#E7E7E3", linewidth = 0.20
      ),
      panel.spacing = unit(0.36, "mm"),
      strip.background = element_rect(
        fill = "#ECEDE7", colour = "#AFAFA9", linewidth = 0.25
      ),
      strip.text = element_text(face = "bold", size = 6.1),
      strip.text.y = element_text(angle = 0, size = 5.8),
      axis.text = element_text(size = 5.9, colour = "#252525"),
      axis.title = element_text(size = 7.1),
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 6.1),
      legend.title = element_text(size = 6.3, face = "bold"),
      legend.key.height = unit(2.6, "mm"),
      legend.spacing.y = unit(0.8, "mm"),
      plot.margin = margin(7, 3, 2, 4)
    )
}

make_support_panel <- function(
  support, maximum, reference, colors, legend_title,
  legend_labels = waiver(), percent_axis = FALSE
) {
  axis_labels <- support$labels
  axis_labels[axis_labels == "1XX"] <- "111"
  x_breaks <- c(0, reference, maximum)
  x_labels <- if (percent_axis) {
    label_percent(accuracy = 1)(x_breaks / maximum)
  } else {
    label_number(
      scale_cut = cut_short_scale(),
      accuracy = 1
    )(x_breaks)
  }
  x_title <- if (percent_axis) {
    paste0(
      "% favor high ploidy\n(N = ",
      format(maximum, big.mark = ","),
      " draws)"
    )
  } else {
    paste0(
      "Supporting draws\n(out of ",
      format(maximum, big.mark = ","),
      ")"
    )
  }

  ggplot() +
    geom_vline(
      xintercept = reference,
      colour = "grey55", linewidth = 0.22, linetype = 2
    ) +
    geom_segment(
      data = support$ranges,
      aes(
        x = minimum_support, xend = maximum_support,
        y = strategy_y, yend = strategy_y
      ),
      colour = "grey78", linewidth = 0.35
    ) +
    geom_point(
      data = support$points,
      aes(
        x = supporting_draws, y = strategy_y,
        colour = group
      ),
      size = 0.82
    ) +
    scale_colour_manual(
      values = colors,
      name = legend_title,
      labels = legend_labels,
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, maximum),
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      breaks = seq_along(support$labels),
      labels = axis_labels,
      limits = c(0.55, length(support$labels) + 0.45),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = x_title,
      y = NULL
    ) +
    theme_minimal(base_size = 8, base_family = "sans") +
    theme(
      text = element_text(colour = "#202020"),
      panel.grid.major.y = element_line(
        colour = "grey92", linewidth = 0.20
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 4.6),
      axis.text.y = element_text(size = 4.35, colour = "#4D4D4D"),
      axis.title.x = element_text(size = 5.0, lineheight = 0.9),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 4.7),
      legend.text = element_text(size = 4.25),
      legend.key.height = unit(1.25, "mm"),
      legend.key.width = unit(1.25, "mm"),
      legend.spacing.y = unit(0.25, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.5, "mm"),
      plot.margin = margin(6, 1, 1, 2)
    ) +
    guides(
      colour = guide_legend(
        ncol = 1, byrow = TRUE,
        title.position = "top",
        override.aes = list(size = 0.8)
      )
    )
}

color_selected_strategy_axis_labels <- function(plot) {
  recolor_text <- function(grob) {
    if (inherits(grob, "text") && length(grob$label)) {
      labels <- as.character(grob$label)
      selected_i <- match(labels, names(strategy_colors))
      if (any(!is.na(selected_i))) {
        current_colors <- grob$gp$col
        if (is.null(current_colors)) current_colors <- "#4D4D4D"
        current_colors <- rep(
          current_colors, length.out = length(labels)
        )
        current_colors[!is.na(selected_i)] <-
          unname(strategy_colors[selected_i[!is.na(selected_i)]])
        grob$gp$col <- current_colors
      }
    }
    if (!is.null(grob$children)) {
      for (i in seq_along(grob$children)) {
        grob$children[[i]] <- recolor_text(grob$children[[i]])
      }
    }
    if (!is.null(grob$grobs)) {
      for (i in seq_along(grob$grobs)) {
        grob$grobs[[i]] <- recolor_text(grob$grobs[[i]])
      }
    }
    grob
  }
  recolor_text(ggplotGrob(plot))
}

make_panel_a <- function(bundle) {
  make_support_panel(
    bundle$support$density,
    2000, 1000, density_colors, "Seeding",
    percent_axis = TRUE
  )
}

make_panel_b <- function(bundle) {
  make_support_panel(
    bundle$support$day,
    2000, 1000, day_colors, "Day",
    legend_labels = c("2", "4", "6"),
    percent_axis = TRUE
  )
}

make_panel_c <- function(bundle) {
  make_support_panel(
    bundle$support$model,
    400, 200, model_colors, "Model"
  )
}

make_panel_d <- function(bundle) {
  make_support_panel(
    bundle$support$line,
    500, 250, line_colors, "Cell line"
  )
}

make_panel_e <- function(bundle) {
  ggplot() +
    geom_tile(
      data = bundle$surface,
      aes(
        x = tile_x,
        y = tile_y,
        width = tile_width,
        height = tile_height,
        fill = delta_growth
      )
    ) +
    geom_path(
      data = bundle$paths,
      aes(
        x = glucose_axis, y = R2, colour = strategy,
        group = interaction(line_name, model, strategy, segment)
      ),
      inherit.aes = FALSE,
      linewidth = 0.55, lineend = "round"
    ) +
    geom_path(
      data = bundle$actions,
      aes(
        x = glucose_axis, y = R2, colour = strategy,
        group = interaction(line_name, model, strategy, action)
      ),
      inherit.aes = FALSE,
      linewidth = 0.48, linetype = 2, lineend = "round"
    ) +
    geom_point(
      data = bundle$markers,
      aes(x = glucose_axis, y = R2, colour = strategy),
      inherit.aes = FALSE,
      size = 0.75, stroke = 0.35
    ) +
    facet_grid(
      rows = vars(line_name),
      cols = vars(model),
      labeller = label_value
    ) +
    scale_fill_gradient2(
      low = "#3E6FB0", mid = "white", high = "#BE4D47",
      midpoint = 0, limits = c(-0.05, 0.05), oob = squish,
      breaks = c(-0.05, 0.05),
      name = "\u0394 growth"
    ) +
    scale_colour_manual(values = strategy_colors, name = "Strategy") +
    scale_x_continuous(
      breaks = surface_transform(c(0, 0.01, 1)),
      labels = c("0", "0.01", "1")
    ) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75)) +
    coord_cartesian(
      xlim = c(
        surface_transform(surface_glucose_pad),
        surface_transform(surface_glucose_max)
      ),
      ylim = c(surface_resource2_pad, 1),
      expand = FALSE, clip = "on"
    ) +
    labs(
      x = "Glucose (mM)",
      y = "Latent resource R2"
    ) +
    theme_f5() +
    theme(panel.grid = element_blank()) +
    guides(
      fill = guide_colourbar(
        order = 1,
        direction = "vertical",
        title.position = "top",
        barwidth = unit(2.8, "mm"),
        barheight = unit(30, "mm")
      ),
      colour = guide_legend(
        order = 2, ncol = 1, title.position = "top"
      )
    )
}

make_panel_f <- function(bundle) {
  ggplot(
    bundle$counts,
    aes(
      x = day, y = live_cells,
      colour = strategy, linetype = ploidy,
      group = interaction(strategy, ploidy)
    )
  ) +
    geom_line(linewidth = 0.55) +
    facet_grid(
      rows = vars(line_name),
      cols = vars(model),
      scales = "free_y",
      labeller = label_value
    ) +
    scale_colour_manual(values = strategy_colors, name = "Strategy") +
    scale_linetype_manual(values = ploidy_linetypes, name = "Ploidy") +
    scale_x_continuous(
      breaks = c(2, 4), limits = c(0, 6),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_log10(labels = label_number()) +
    labs(
      x = "Day",
      y = "Live cells"
    ) +
    theme_f5() +
    guides(
      colour = guide_legend(
        order = 1, ncol = 1, title.position = "top"
      ),
      linetype = guide_legend(
        order = 2, ncol = 1, title.position = "top"
      )
    )
}

make_panel_g <- function(bundle) {
  ggplot(
    bundle$growth,
    aes(
      x = day, y = net_growth,
      colour = strategy, linetype = ploidy,
      group = interaction(strategy, ploidy)
    )
  ) +
    geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.25) +
    geom_line(linewidth = 0.52) +
    facet_grid(
      rows = vars(line_name),
      cols = vars(model),
      scales = "free_y",
      labeller = label_value
    ) +
    scale_colour_manual(values = strategy_colors, name = "Strategy") +
    scale_linetype_manual(values = ploidy_linetypes, name = "Ploidy") +
    scale_x_continuous(
      breaks = c(2, 4), limits = c(0, 6),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Day",
      y = expression("Growth (h"^-1*")")
    ) +
    theme_f5() +
    guides(
      colour = guide_legend(
        order = 1, ncol = 1, title.position = "top"
      ),
      linetype = guide_legend(
        order = 2, ncol = 1, title.position = "top"
      )
    )
}

make_panels <- function(bundle) {
  bundle <- set_plot_factors(bundle)
  list(
    a = color_selected_strategy_axis_labels(make_panel_a(bundle)),
    b = color_selected_strategy_axis_labels(make_panel_b(bundle)),
    c = color_selected_strategy_axis_labels(make_panel_c(bundle)),
    d = color_selected_strategy_axis_labels(make_panel_d(bundle)),
    e = make_panel_e(bundle),
    f = make_panel_f(bundle),
    g = make_panel_g(bundle)
  )
}

write_subpanels <- function(panels) {
  dimensions <- data.table(
    figure = "Figure 5",
    panel = letters[1:7],
    subpanel_png = file.path(
      "manuscript_figures", "20260728_posterior_strategy_figure_set",
      "F5", "polishing", "subpanels",
      paste0("figure_5", letters[1:7], ".png")
    ),
    width_in = c(rep(1.75, 4), rep(5.25, 3)),
    height_in = c(rep(9.25 / 4, 4), 3.55, 2.85, 2.85),
    order = 1:7
  )
  dpi <- 600
  for (i in seq_len(nrow(dimensions))) {
    ggsave(
      filename = file.path(project_root, dimensions$subpanel_png[[i]]),
      plot = panels[[dimensions$panel[[i]]]],
      width = dimensions$width_in[[i]],
      height = dimensions$height_in[[i]],
      dpi = dpi,
      bg = "white"
    )
  }
  dimensions[, `:=`(
    width_px = as.integer(width_in * dpi),
    height_px = as.integer(height_in * dpi)
  )]
  setcolorder(dimensions, c(
    "figure", "panel", "subpanel_png",
    "width_px", "height_px", "width_in", "height_in", "order"
  ))
  fwrite(
    dimensions,
    file.path(layout_dir, "subpanel_dimensions.csv")
  )
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

write_final_metadata <- function(final_path) {
  relative_final <- sub(
    paste0("^", project_root, "/"), "", normalizePath(final_path)
  )
  script_path <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "F5", "polishing", "scripts", "polish_figures.R"
  )
  layout_path <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "F5", "polishing", "layout", "layout_plan.csv"
  )
  subpanels <- file.path(
    "manuscript_figures", "20260728_posterior_strategy_figure_set",
    "F5", "polishing", "subpanels",
    paste0("figure_5", letters[1:7], ".png")
  )
  source_scripts <- c(
    "tmp/f5_panel_a_draft_options_20260729/build_draft_options.R",
    "tmp/regenerate_figure5b_single_model_asinh_surface.R",
    "tmp/figure5b_single_model_playground.R",
    "tmp/posterior_strategy_success_20260727/heatmap_survey/scripts/build_fourline_allmodels_alllines_strategy_overlay.R",
    "tmp/posterior_strategy_success_20260727/heatmap_survey/scripts/build_fourline_selected_strategy_timecourses.R"
  )
  direct_inputs <- c(
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/posterior/strategy_simulations/runs/seed3density_7g_xxx_xcc_axax_20260727/manifest.tsv",
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/posterior/strategy_simulations/runs/seed3density_7g_xxx_xcc_axax_20260727/validation.json",
    "data/modelling/gpath_v1/red_a30_counts_20260722/derived/plans/prediction_tasks.tsv",
    "data/modelling/gpath_v1/red_a30_counts_20260722/datasets/loo_exclude_sum_159_fuse/stan_data.Rds"
  )
  fwrite(
    data.table(
      figure = "Figure 5",
      source_local_name = relative_final,
      source_root = dirname(relative_final),
      panel_labels = paste(letters[1:7], collapse = ";"),
      polished_final_image = relative_final,
      polished_image_sha256 = sha256_file(final_path),
      source_files = paste(source_scripts, collapse = ";"),
      data_inputs = paste(direct_inputs, collapse = ";"),
      subpanel_audit_paths = paste(subpanels, collapse = ";"),
      polishing_script = script_path
    ),
    file.path(polish_root, "manifest.csv")
  )
  panel_notes <- c(
    "Day-6 percentage of 2,000 posterior draws favoring high ploidy across three seeding densities, five models, and four included lines; 0AA omitted and strategies fixed in ascending initial-glucose then XX/AA/CC order.",
    "Percentage of 2,000 posterior draws favoring high ploidy at days 2, 4, and 6 at 1x seeding, pooled across five models and four included lines.",
    "Day-6 support at 1x seeding separated by five Pareto models; each model pools four included lines.",
    "Day-6 support at 1x seeding separated by four included lines; each line pools five Pareto models.",
    "Regenerated from 100 matched draws per model on exact zero-inclusive 120x40 surfaces uniform in asinh(glucose/1e-4), with zero-boundary tiles extended into display-only negative padding.",
    "Regenerated q50 low/high live-cell trajectories directly from the 20 compact simulation shards.",
    "Regenerated by evaluating 100 matched posterior draws along each median environmental trajectory."
  )
  data_by_panel <- c(
    paste(direct_inputs[c(1, 2)], collapse = ";"),
    paste(direct_inputs[c(1, 2)], collapse = ";"),
    paste(direct_inputs[c(1, 2)], collapse = ";"),
    paste(direct_inputs[c(1, 2)], collapse = ";"),
    paste(direct_inputs, collapse = ";"),
    paste(direct_inputs[c(1, 2)], collapse = ";"),
    paste(direct_inputs, collapse = ";")
  )
  fwrite(
    data.table(
      figure = "Figure 5",
      panel = letters[1:7],
      subpanel_image = subpanels,
      generator = script_path,
      command = paste(
        "scripts/agentRrunner.sh", script_path, "--phase final"
      ),
      data_inputs = data_by_panel,
      layout_plan = layout_path,
      output_image = relative_final,
      notes = panel_notes
    ),
    file.path(polish_root, "provenance.csv")
  )
  fwrite(
    data.table(
      figure = "Figure 5",
      stage = "polishing",
      output_image = relative_final,
      output_sha256 = sha256_file(final_path),
      generator = script_path,
      command = paste(
        "scripts/agentRrunner.sh", script_path, "--phase final"
      ),
      data_inputs = paste(direct_inputs, collapse = ";"),
      layout_plan = layout_path
    ),
    file.path(polish_root, "figure_rebuild_manifest.tsv"),
    sep = "\t"
  )
  writeLines(
    c(
      "# Figure 5 polishing notes",
      "",
      "- Panels a-g were regenerated from validated posterior outputs; no draft PNG was embedded.",
      "- Panels a-d report support counts for 20 strategies; 0AA is omitted and the y axis is fixed by ascending initial glucose followed by XX, AA, and CC.",
      "- Panel a compares 0.5x, 1x, and 2x seeding at day 6; panels b-d use 1x seeding and separate evaluation day, model, or cell line.",
      "- Panels a-b display support as a percentage of 2,000 posterior draws; panels c-d retain their natural per-model and per-line supporting-draw counts.",
      "- In panels a-d, the 0CC, 5CC, and 111 y-axis labels use the strategy colors from panels e-g; support-group points use a separate qualitative palette.",
      "- Panels e-g preserve the four-line fit, four included cell lines, five Pareto models, 1x seeding, and strategies 0CC, 5CC, and 111.",
      "- Panel e uses exact zero-inclusive 120x40 growth surfaces uniform in asinh(glucose / 1e-4); evaluated zero-boundary tiles extend into display-only negative glucose and R2 padding.",
      "- Strategy colors and low/high ploidy linetypes match the approved exploratory plots.",
      "- Panels e-g each display all guides relevant to their encodings.",
      "- Panel f reports modeled live-cell abundance on the mean-per-field scale; it is not a whole-culture cell total.",
      "- The composite uses the full 7 x 9.25 inch envelope; all guides are right-side and panel-specific.",
      "- The final composite contains no figure-level title, subtitle, caption, or semantic interpretation.",
      "- Typography and axis-break density were reduced after visual QC of the compressed layout.",
      "- Project-map decision: no project_map update is needed while this remains a review-stage potential replacement rather than a canonical manuscript figure."
    ),
    file.path(polish_root, "notes.md")
  )
}

load_or_build_panel_data <- function() {
  if (file.exists(cache_path)) {
    bundle <- qread(cache_path, nthreads = 2L)
    cache_current <- identical(
      bundle$metadata$surface_transform,
      "asinh(glucose / 1e-4)"
    ) && identical(
      sort(names(bundle$support)),
      sort(c(
        "density", "day", "model", "line",
        "strategy_order_top_to_bottom"
      ))
    )
    if (isTRUE(cache_current)) return(bundle)
  }
  build_panel_data()
}

if (!library_only && phase == "subpanels") {
  bundle <- load_or_build_panel_data()
  panels <- make_panels(bundle)
  write_subpanels(panels)
  cat("Wrote regenerated Figure 5 audit subpanels and panel-data cache.\n")
}

if (!library_only && phase == "final") {
  bundle <- load_or_build_panel_data()
  panels <- make_panels(bundle)
  layout <- fread(file.path(layout_dir, "layout_plan.csv"))
  layout <- layout[figure == "Figure 5"][match(letters[1:7], panel)]
  if (nrow(layout) != 7L || anyNA(layout$panel)) {
    stop("Layout plan does not contain Figure 5 panels a-g.", call. = FALSE)
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
        x = row$x_npc + 0.002,
        y = row$y_npc + row$height_npc - 0.002,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = 8.5,
        fontfamily = "sans"
      )
  }
  final_path <- file.path(final_dir, "figure_5.png")
  ggsave(
    final_path,
    canvas,
    width = unique(layout$layout_width_in),
    height = unique(layout$layout_height_in),
    dpi = 600,
    bg = "white"
  )
  write_final_metadata(final_path)
  cat("Wrote polished Figure 5 composite: ", final_path, "\n", sep = "")
}
