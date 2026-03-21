source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/parameter_transfer_utils.R")

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

selection_default_config <- function() {
  list(
    project_root = getwd(),
    model_name = "gpath",
    model_version = "v1",
    dataset_label = "gstarvation_v1",
    stan_data_path = NULL,
    assessment_run_id_path = file.path("models", "gpath", "v1", "assessment_run_ids.txt"),
    optim_root = file.path("data", "runs", "gpath", "v1", "optim"),
    transfer_output_root = file.path("data", "gpath_transfer_cv"),
    output_root = file.path("data", "selection_strategy", "gpath", "v1"),
    include_global = TRUE,
    include_transfer = TRUE,
    transfer_fit_types = c("transfer"),
    transfer_directions = c("low_to_high", "high_to_low"),
    model_ids = NULL,
    line_ids = NULL,
    glucose_values = NULL,
    interval_hours = 48,
    total_hours = 144,
    time_step_hours = 1,
    initial_high_fraction = 0.5,
    refresh_resets_total_n = TRUE,
    reset_total_n_mode = "estimated_mean",
    detailed_default = FALSE,
    detailed_time_step_hours = 1,
    representative_model_ids = NULL,
    representative_fit_contexts = c("global"),
    workers = max(1L, parallel::detectCores(logical = TRUE) - 1L)
  )
}

merge_selection_config <- function(defaults, overrides = list()) {
  utils::modifyList(defaults, overrides)
}

load_selection_config <- function(config_path = NULL, overrides = list()) {
  cfg <- selection_default_config()
  if (!is.null(config_path) && nzchar(config_path)) {
    loaded <- source(config_path, local = new.env(parent = globalenv()))$value
    if (!is.list(loaded)) {
      stop(sprintf("Config at '%s' did not evaluate to a list", config_path))
    }
    cfg <- merge_selection_config(cfg, loaded)
  }
  merge_selection_config(cfg, overrides)
}

resolve_selection_output_root <- function(config) {
  file.path(config$output_root, config$dataset_label)
}

resolve_selection_stan_data <- function(config) {
  resolve_stan_data_path(config$stan_data_path)
}

read_assessment_run_ids <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  lines <- trimws(readLines(path, warn = FALSE))
  lines[nzchar(lines)]
}

discover_global_model_ids <- function(config) {
  if (!is.null(config$model_ids) && length(config$model_ids)) {
    return(unique(as.character(config$model_ids)))
  }

  dataset_root <- file.path(config$optim_root, config$dataset_label)
  if (dir.exists(dataset_root)) {
    dirs <- list.dirs(dataset_root, recursive = FALSE, full.names = FALSE)
    dirs <- dirs[nzchar(dirs)]
    if (length(dirs)) {
      return(sort(unique(dirs)))
    }
  }

  assessment_ids <- read_assessment_run_ids(config$assessment_run_id_path)
  if (length(assessment_ids)) {
    return(sort(unique(assessment_ids)))
  }

  stop("Could not resolve any gpath model IDs from config$model_ids, local runs, or assessment_run_ids.txt")
}

resolve_selection_model_ids <- function(config) {
  discover_global_model_ids(config)
}

resolve_selection_stan <- function(config) {
  stan_data <- readRDS(resolve_selection_stan_data(config))
  if (is.null(stan_data$group_id)) {
    stan_data <- add_group_structure(stan_data)
  }
  stan_data
}

resolve_glucose_values <- function(stan_data, config) {
  if (!is.null(config$glucose_values) && length(config$glucose_values)) {
    return(sort(unique(as.numeric(config$glucose_values))))
  }
  vals <- sort(unique(as.numeric(stan_data$G0_per_well)))
  vals[is.finite(vals)]
}

resolve_line_ids <- function(stan_data, config) {
  line_ids <- sort(unique(as.integer(stan_data$line_id)))
  line_ids <- line_ids[is.finite(line_ids)]
  if (!is.null(config$line_ids) && length(config$line_ids)) {
    line_ids <- intersect(line_ids, as.integer(config$line_ids))
  }
  line_ids
}

get_line_label_map <- function(stan_data) {
  line_ids <- as.integer(stan_data$line_id)
  labels <- as.character(stan_data$cellLine_per_well)
  keep <- is.finite(line_ids) & !is.na(labels)
  pairs <- unique(data.frame(line_id = line_ids[keep], cell_line = labels[keep], stringsAsFactors = FALSE))
  pairs[order(pairs$line_id), , drop = FALSE]
}

get_line_ploidy_pair <- function(stan_data, line_id, tol = 1e-12) {
  states <- get_line_ploidy_states(stan_data, line_id = line_id, tol = tol)
  list(
    low = states$low_value,
    high = states$high_value
  )
}

get_group_id_for_state <- function(stan_data, line_id, ploidy_metric, tol = 1e-12) {
  idx <- which(
    as.integer(stan_data$line_id) == as.integer(line_id) &
      abs(as.numeric(stan_data$ploidy_metric) - as.numeric(ploidy_metric)) < tol
  )
  group_ids <- sort(unique(as.integer(stan_data$group_id[idx])))
  group_ids <- group_ids[is.finite(group_ids)]

  if (length(group_ids) != 1L) {
    stop(sprintf(
      "Expected exactly one group_id for line=%d ploidy_metric=%s; found %d",
      as.integer(line_id),
      format(ploidy_metric, digits = 6),
      length(group_ids)
    ))
  }

  group_ids[[1]]
}

extract_best_draw_from_optim_outputs <- function(draws_all, lp_all) {
  lp_vec <- as.numeric(lp_all)
  valid_idx <- which(is.finite(lp_vec))
  if (!length(valid_idx)) {
    stop("No finite optimization scores were found")
  }

  best_idx <- valid_idx[which.max(lp_vec[valid_idx])]
  draw_obj <- draws_all[[best_idx]]
  draw_vec <- extract_draw_vector(draw_obj)

  list(
    draw_vec = draw_vec,
    best_idx = best_idx,
    best_lp = lp_vec[[best_idx]]
  )
}

resolve_global_run_dir <- function(config, model_id, run_label = NULL) {
  model_root <- file.path(config$optim_root, config$dataset_label, model_id)
  if (!dir.exists(model_root)) {
    stop(sprintf("Global optimization directory not found for model '%s': %s", model_id, model_root))
  }

  run_dirs <- list.dirs(model_root, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[file.exists(file.path(run_dirs, "optim_draws_all.Rds")) & file.exists(file.path(run_dirs, "optim_lp_all.Rds"))]

  if (!length(run_dirs)) {
    stop(sprintf("No completed optimization run directories found for model '%s' under %s", model_id, model_root))
  }

  if (!is.null(run_label) && nzchar(run_label)) {
    hit <- run_dirs[basename(run_dirs) == run_label]
    if (!length(hit)) {
      stop(sprintf("Requested run_label '%s' not found for model '%s'", run_label, model_id))
    }
    return(hit[[1]])
  }

  sort(run_dirs)[[length(run_dirs)]]
}

load_best_global_fit <- function(config, model_id, run_label = NULL) {
  run_dir <- resolve_global_run_dir(config, model_id = model_id, run_label = run_label)
  draws_all <- readRDS(file.path(run_dir, "optim_draws_all.Rds"))
  lp_all <- readRDS(file.path(run_dir, "optim_lp_all.Rds"))
  best <- extract_best_draw_from_optim_outputs(draws_all = draws_all, lp_all = lp_all)

  list(
    draw_vec = best$draw_vec,
    best_lp = best$best_lp,
    best_idx = best$best_idx,
    run_dir = run_dir,
    run_label = basename(run_dir),
    ploidy_effect_mask = NULL
  )
}

load_best_transfer_fit_context <- function(config, model_id, line_id, direction, fit_type) {
  fit <- load_best_transfer_fit(
    model_id = model_id,
    line_id = line_id,
    direction = direction,
    fit_type = fit_type,
    output_root = config$transfer_output_root,
    model_name = config$model_name
  )

  list(
    draw_vec = extract_draw_vector(fit$draws),
    best_lp = fit$summary$lp__[[1]],
    run_dir = file.path(config$transfer_output_root, config$model_name, model_id, direction, fit_type),
    run_label = fit$run_tag,
    ploidy_effect_mask = fit$split_meta$ploidy_effect_mask %||% NULL,
    split_meta = fit$split_meta
  )
}

global_fit_available <- function(config, model_id) {
  model_root <- file.path(config$optim_root, config$dataset_label, model_id)
  if (!dir.exists(model_root)) {
    return(FALSE)
  }
  run_dirs <- list.dirs(model_root, recursive = FALSE, full.names = TRUE)
  any(file.exists(file.path(run_dirs, "optim_draws_all.Rds")) & file.exists(file.path(run_dirs, "optim_lp_all.Rds")))
}

transfer_fit_available <- function(config, model_id, line_id, direction, fit_type) {
  best_path <- file.path(config$transfer_output_root, config$model_name, model_id, "transfer_best_start_summary.Rds")
  if (!file.exists(best_path)) {
    return(FALSE)
  }

  best_df <- readRDS(best_path)
  direction <- normalize_transfer_direction(direction)
  fit_type <- normalize_fit_type(fit_type)

  any(
    best_df$line_id == as.integer(line_id) &
      best_df$direction == direction &
      best_df$fit_type == fit_type
  )
}

estimate_state_n0 <- function(draw_vec, stan_data, line_id, ploidy_metric, tol = 1e-12) {
  group_id <- get_group_id_for_state(stan_data, line_id = line_id, ploidy_metric = ploidy_metric, tol = tol)
  raw_name <- sprintf("raw_N0[%d]", group_id)
  raw_val <- draw_vec[[raw_name]]
  if (!length(raw_val) || !is.finite(raw_val)) {
    stop(sprintf("Could not read '%s' from fitted parameter vector", raw_name))
  }
  exp(log(500) + as.numeric(raw_val))
}

estimate_mixture_seed_total_n <- function(draw_vec, stan_data, line_id, low_ploidy, high_ploidy, tol = 1e-12) {
  low_n0 <- estimate_state_n0(draw_vec, stan_data, line_id = line_id, ploidy_metric = low_ploidy, tol = tol)
  high_n0 <- estimate_state_n0(draw_vec, stan_data, line_id = line_id, ploidy_metric = high_ploidy, tol = tol)
  mean(c(low_n0, high_n0))
}

reconstruct_mixture_components <- function(draw_vec, model_id, line_id, low_ploidy, high_ploidy, ploidy_effect_mask = NULL) {
  dims <- parse_run_id(model_id)
  L <- 3 * dims$R + (dims$P - 1) * dims$R + dims$W * dims$R + 1
  raw_theta_line <- draw_vec[sprintf("raw_theta_line[%d,%d]", 1:L, line_id)]
  raw_theta_ploidy <- draw_vec[sprintf("raw_theta_ploidy[%d]", 1:L)]

  parms_low <- reconstruct_parms(
    R = dims$R,
    P = dims$P,
    W = dims$W,
    strict_spec = (dims$C == 1L),
    M = dims$M,
    base_priors = base_priors,
    raw_theta_line = raw_theta_line,
    raw_theta_ploidy = raw_theta_ploidy,
    ploidy_metric = low_ploidy,
    ploidy_effect_mask = ploidy_effect_mask
  )

  parms_high <- reconstruct_parms(
    R = dims$R,
    P = dims$P,
    W = dims$W,
    strict_spec = (dims$C == 1L),
    M = dims$M,
    base_priors = base_priors,
    raw_theta_line = raw_theta_line,
    raw_theta_ploidy = raw_theta_ploidy,
    ploidy_metric = high_ploidy,
    ploidy_effect_mask = ploidy_effect_mask
  )

  list(
    parms_low = parms_low,
    parms_high = parms_high,
    parms_mix = combine_parms(parms_low, parms_high),
    dims = dims
  )
}

build_strategy_grid <- function(glucose_values) {
  glucose_values <- as.numeric(glucose_values)
  refresh_idx <- 0:length(glucose_values)
  grid <- expand.grid(
    day0_idx = seq_along(glucose_values),
    day2_idx = refresh_idx,
    day4_idx = refresh_idx,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  map_refresh_value <- function(idx) {
    ifelse(idx == 0L, NA_real_, glucose_values[idx])
  }

  grid$strategy_code <- sprintf("%d%d%d", grid$day0_idx, grid$day2_idx, grid$day4_idx)
  grid$g0_day0 <- glucose_values[grid$day0_idx]
  grid$g0_day2 <- map_refresh_value(grid$day2_idx)
  grid$g0_day4 <- map_refresh_value(grid$day4_idx)
  grid$day2_action <- ifelse(grid$day2_idx == 0L, "carry", "refresh")
  grid$day4_action <- ifelse(grid$day4_idx == 0L, "carry", "refresh")
  grid
}

make_initial_mix_state <- function(parms_mix, g0, seed_total_n, initial_high_fraction = 0.5) {
  n_high <- seed_total_n * initial_high_fraction
  n_low <- seed_total_n - n_high

  y0 <- c(n_low, n_high, g0)
  if (parms_mix$R > 1L) {
    y0 <- c(y0, rep(1.0, parms_mix$R - 1L))
  }
  if (parms_mix$W > 0L) {
    y0 <- c(y0, rep(0.0, parms_mix$W))
  }
  y0
}

state_column_names <- function(parms_mix) {
  cols <- c("N_low", "N_high", "G1")
  if (parms_mix$R > 1L) {
    cols <- c(cols, sprintf("R%d", 2:parms_mix$R))
  }
  if (parms_mix$W > 0L) {
    cols <- c(cols, sprintf("W%d", 1:parms_mix$W))
  }
  cols
}

simulate_mixture_segment <- function(
  parms_mix,
  state0,
  hours,
  time_step_hours = 1,
  detailed = FALSE,
  method = "lsoda"
) {
  times <- seq(0, hours, by = time_step_hours)
  if (tail(times, 1) < hours) {
    times <- c(times, hours)
  }
  if (length(times) == 1L) {
    times <- c(0, hours)
  }

  out <- deSolve::ode(
    y = state0,
    times = times,
    func = rhs_mix,
    parms = parms_mix,
    method = method
  )

  out <- as.data.frame(out)
  names(out) <- c("time", state_column_names(parms_mix))
  out$total_live <- out$N_low + out$N_high
  out$high_fraction <- ifelse(out$total_live > 0, out$N_high / out$total_live, NA_real_)

  if (isTRUE(detailed)) {
    out
  } else {
    out[nrow(out), , drop = FALSE]
  }
}

reset_state_for_refresh <- function(
  state,
  parms_mix,
  g0,
  seed_total_n,
  reset_total_n = TRUE
) {
  new_state <- as.numeric(state)

  total_live <- sum(new_state[1:2])
  high_fraction <- if (is.finite(total_live) && total_live > 0) {
    new_state[2] / total_live
  } else {
    0.5
  }

  if (isTRUE(reset_total_n)) {
    new_state[2] <- seed_total_n * high_fraction
    new_state[1] <- seed_total_n - new_state[2]
  }

  new_state[3] <- g0
  if (parms_mix$R > 1L) {
    new_state[4:(2 + parms_mix$R)] <- 1.0
  }
  if (parms_mix$W > 0L) {
    waste_idx <- (3 + parms_mix$R):(2 + parms_mix$R + parms_mix$W)
    new_state[waste_idx] <- 0.0
  }

  new_state
}

simulate_strategy <- function(
  strategy_row,
  parms_mix,
  seed_total_n,
  initial_high_fraction = 0.5,
  interval_hours = 48,
  time_step_hours = 1,
  reset_total_n_on_refresh = TRUE,
  detailed = FALSE
) {
  segment_g0 <- c(
    strategy_row$g0_day0,
    strategy_row$g0_day2,
    strategy_row$g0_day4
  )
  segment_action <- c("seed", strategy_row$day2_action, strategy_row$day4_action)

  state <- make_initial_mix_state(
    parms_mix = parms_mix,
    g0 = segment_g0[[1]],
    seed_total_n = seed_total_n,
    initial_high_fraction = initial_high_fraction
  )

  summary_rows <- vector("list", 3L)
  detail_rows <- vector("list", 3L)
  time_offset <- 0

  for (segment_idx in 1:3) {
    if (segment_idx > 1L) {
      if (segment_action[[segment_idx]] == "refresh") {
        state <- reset_state_for_refresh(
          state = state,
          parms_mix = parms_mix,
          g0 = segment_g0[[segment_idx]],
          seed_total_n = seed_total_n,
          reset_total_n = reset_total_n_on_refresh
        )
      }
    }

    seg_out <- simulate_mixture_segment(
      parms_mix = parms_mix,
      state0 = state,
      hours = interval_hours,
      time_step_hours = time_step_hours,
      detailed = TRUE
    )

    seg_out$segment_idx <- segment_idx
    seg_out$segment_start_hours <- time_offset
    seg_out$time_global <- seg_out$time + time_offset

    end_row <- seg_out[nrow(seg_out), , drop = FALSE]
    end_row$strategy_code <- strategy_row$strategy_code
    end_row$day0_idx <- strategy_row$day0_idx
    end_row$day2_idx <- strategy_row$day2_idx
    end_row$day4_idx <- strategy_row$day4_idx
    end_row$segment_end_hours <- time_offset + interval_hours
    end_row$segment_action <- segment_action[[segment_idx]]
    end_row$refresh_glucose <- segment_g0[[segment_idx]]
    summary_rows[[segment_idx]] <- end_row

    if (isTRUE(detailed)) {
      if (segment_idx > 1L) {
        seg_out <- seg_out[-1, , drop = FALSE]
      }
      seg_out$strategy_code <- strategy_row$strategy_code
      seg_out$day0_idx <- strategy_row$day0_idx
      seg_out$day2_idx <- strategy_row$day2_idx
      seg_out$day4_idx <- strategy_row$day4_idx
      seg_out$segment_action <- segment_action[[segment_idx]]
      seg_out$refresh_glucose <- segment_g0[[segment_idx]]
      detail_rows[[segment_idx]] <- seg_out
    }

    state <- unname(as.numeric(end_row[, state_column_names(parms_mix), drop = TRUE]))
    time_offset <- time_offset + interval_hours
  }

  summary_df <- do.call(rbind, summary_rows)
  final_row <- summary_df[nrow(summary_df), , drop = FALSE]
  final_row$seed_total_n <- seed_total_n
  final_row$log_ratio_final <- log((final_row$N_high + 1e-8) / (final_row$N_low + 1e-8))
  final_row$high_fraction_final <- ifelse(final_row$total_live > 0, final_row$N_high / final_row$total_live, NA_real_)
  final_row$low_fold_change <- final_row$N_low / max(seed_total_n * (1 - initial_high_fraction), 1e-8)
  final_row$high_fold_change <- final_row$N_high / max(seed_total_n * initial_high_fraction, 1e-8)

  list(
    interval_summary = summary_df,
    final_summary = final_row,
    detail = if (isTRUE(detailed)) do.call(rbind, detail_rows) else NULL
  )
}

score_strategy_table <- function(final_df) {
  out <- final_df
  out$select_against_high_score <- -out$log_ratio_final
  out$select_for_high_score <- out$log_ratio_final
  out$viability_score <- out$total_live
  out$against_high_viable <- out$select_against_high_score
  out$against_high_viable[out$low_fold_change < 1] <- NA_real_
  out$for_high_viable <- out$select_for_high_score
  out$for_high_viable[out$high_fold_change < 1] <- NA_real_
  out
}

rank_strategy_table <- function(final_df, score_col = "against_high_viable", decreasing = TRUE) {
  if (!(score_col %in% names(final_df))) {
    stop(sprintf("score_col '%s' not found in strategy table", score_col))
  }
  ranked <- final_df[order(final_df[[score_col]], final_df$viability_score, decreasing = decreasing, na.last = TRUE), , drop = FALSE]
  ranked$rank <- seq_len(nrow(ranked))
  ranked
}

select_top_strategies <- function(final_df, n_top = 5L) {
  scored <- score_strategy_table(final_df)
  bind_top <- function(score_col, label) {
    ranked <- rank_strategy_table(scored, score_col = score_col, decreasing = TRUE)
    ranked <- ranked[seq_len(min(n_top, nrow(ranked))), , drop = FALSE]
    ranked$selection_objective <- label
    ranked$selection_score_col <- score_col
    ranked
  }

  do.call(rbind, list(
    bind_top("against_high_viable", "select_against_high_viable"),
    bind_top("select_against_high_score", "select_against_high"),
    bind_top("for_high_viable", "select_for_high_viable"),
    bind_top("select_for_high_score", "select_for_high")
  ))
}

build_selection_tasks <- function(config, stan_data) {
  model_ids <- resolve_selection_model_ids(config)
  line_ids <- resolve_line_ids(stan_data, config)
  glucose_values <- resolve_glucose_values(stan_data, config)

  tasks <- list()
  idx <- 1L

  if (isTRUE(config$include_global)) {
    for (model_id in model_ids) {
      for (line_id in line_ids) {
        tasks[[idx]] <- data.frame(
          model_id = model_id,
          line_id = as.integer(line_id),
          fit_context = "global",
          direction = NA_character_,
          fit_type = NA_character_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (isTRUE(config$include_transfer)) {
    for (model_id in model_ids) {
      for (line_id in line_ids) {
        for (direction in config$transfer_directions) {
          for (fit_type in config$transfer_fit_types) {
            tasks[[idx]] <- data.frame(
              model_id = model_id,
              line_id = as.integer(line_id),
              fit_context = "transfer",
              direction = normalize_transfer_direction(direction),
              fit_type = normalize_fit_type(fit_type),
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }

  out <- do.call(rbind, tasks)
  if (!nrow(out)) {
    stop("No selection-strategy tasks were constructed from the current config")
  }

  keep <- vapply(seq_len(nrow(out)), function(i) {
    row <- out[i, , drop = FALSE]
    if (identical(row$fit_context[[1]], "global")) {
      global_fit_available(config, row$model_id[[1]])
    } else {
      transfer_fit_available(
        config = config,
        model_id = row$model_id[[1]],
        line_id = row$line_id[[1]],
        direction = row$direction[[1]],
        fit_type = row$fit_type[[1]]
      )
    }
  }, logical(1))

  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) {
    stop("No selection-strategy tasks matched available fit outputs")
  }

  attr(out, "glucose_values") <- glucose_values
  out
}

task_output_dir <- function(config, task_row) {
  fit_label <- if (task_row$fit_context == "global") {
    "global"
  } else {
    sprintf("%s_%s", task_row$direction, task_row$fit_type)
  }

  file.path(
    resolve_selection_output_root(config),
    task_row$model_id,
    fit_label,
    sprintf("line_%02d", as.integer(task_row$line_id))
  )
}

save_task_results <- function(config, task_row, result) {
  out_dir <- task_output_dir(config, task_row)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result$interval_summary, file.path(out_dir, "interval_summary.Rds"))
  saveRDS(result$final_summary, file.path(out_dir, "final_summary.Rds"))
  saveRDS(result$task_meta, file.path(out_dir, "task_meta.Rds"))
  if (!is.null(result$detail)) {
    saveRDS(result$detail, file.path(out_dir, "detail.Rds"))
  }
  invisible(out_dir)
}

simulate_selection_task <- function(task_row, config, stan_data, detailed = FALSE) {
  task_row <- as.list(task_row)
  line_pair <- get_line_ploidy_pair(stan_data, line_id = task_row$line_id)

  fit_obj <- if (identical(task_row$fit_context, "global")) {
    load_best_global_fit(config = config, model_id = task_row$model_id)
  } else {
    load_best_transfer_fit_context(
      config = config,
      model_id = task_row$model_id,
      line_id = task_row$line_id,
      direction = task_row$direction,
      fit_type = task_row$fit_type
    )
  }

  comp <- reconstruct_mixture_components(
    draw_vec = fit_obj$draw_vec,
    model_id = task_row$model_id,
    line_id = task_row$line_id,
    low_ploidy = line_pair$low,
    high_ploidy = line_pair$high,
    ploidy_effect_mask = fit_obj$ploidy_effect_mask
  )

  seed_total_n <- estimate_mixture_seed_total_n(
    draw_vec = fit_obj$draw_vec,
    stan_data = stan_data,
    line_id = task_row$line_id,
    low_ploidy = line_pair$low,
    high_ploidy = line_pair$high
  )

  glucose_values <- resolve_glucose_values(stan_data, config)
  strategy_grid <- build_strategy_grid(glucose_values)

  all_interval <- vector("list", nrow(strategy_grid))
  all_final <- vector("list", nrow(strategy_grid))
  all_detail <- vector("list", nrow(strategy_grid))

  for (i in seq_len(nrow(strategy_grid))) {
    sim <- simulate_strategy(
      strategy_row = strategy_grid[i, , drop = FALSE],
      parms_mix = comp$parms_mix,
      seed_total_n = seed_total_n,
      initial_high_fraction = config$initial_high_fraction,
      interval_hours = config$interval_hours,
      time_step_hours = if (isTRUE(detailed)) config$detailed_time_step_hours else config$time_step_hours,
      reset_total_n_on_refresh = isTRUE(config$refresh_resets_total_n),
      detailed = detailed
    )

    interval_df <- sim$interval_summary
    interval_df$model_id <- task_row$model_id
    interval_df$line_id <- task_row$line_id
    interval_df$fit_context <- task_row$fit_context
    interval_df$direction <- task_row$direction %||% NA_character_
    interval_df$fit_type <- task_row$fit_type %||% NA_character_
    interval_df$low_ploidy_metric <- line_pair$low
    interval_df$high_ploidy_metric <- line_pair$high
    all_interval[[i]] <- interval_df

    final_df <- sim$final_summary
    final_df$model_id <- task_row$model_id
    final_df$line_id <- task_row$line_id
    final_df$fit_context <- task_row$fit_context
    final_df$direction <- task_row$direction %||% NA_character_
    final_df$fit_type <- task_row$fit_type %||% NA_character_
    final_df$low_ploidy_metric <- line_pair$low
    final_df$high_ploidy_metric <- line_pair$high
    final_df$seed_total_n <- seed_total_n
    all_final[[i]] <- final_df

    if (isTRUE(detailed)) {
      detail_df <- sim$detail
      detail_df$model_id <- task_row$model_id
      detail_df$line_id <- task_row$line_id
      detail_df$fit_context <- task_row$fit_context
      detail_df$direction <- task_row$direction %||% NA_character_
      detail_df$fit_type <- task_row$fit_type %||% NA_character_
      detail_df$low_ploidy_metric <- line_pair$low
      detail_df$high_ploidy_metric <- line_pair$high
      all_detail[[i]] <- detail_df
    }
  }

  final_summary <- score_strategy_table(do.call(rbind, all_final))

  list(
    interval_summary = do.call(rbind, all_interval),
    final_summary = final_summary,
    detail = if (isTRUE(detailed)) do.call(rbind, all_detail) else NULL,
    task_meta = data.frame(
      model_id = task_row$model_id,
      line_id = task_row$line_id,
      fit_context = task_row$fit_context,
      direction = task_row$direction %||% NA_character_,
      fit_type = task_row$fit_type %||% NA_character_,
      run_label = fit_obj$run_label %||% NA_character_,
      run_dir = fit_obj$run_dir %||% NA_character_,
      seed_total_n = seed_total_n,
      low_ploidy_metric = line_pair$low,
      high_ploidy_metric = line_pair$high,
      stringsAsFactors = FALSE
    )
  )
}

run_selection_tasks <- function(config, tasks = NULL, detailed = NULL) {
  stan_data <- resolve_selection_stan(config)
  tasks <- tasks %||% build_selection_tasks(config, stan_data)
  detailed <- isTRUE(detailed %||% config$detailed_default)
  root_dir <- resolve_selection_output_root(config)
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

  task_rows <- split(tasks, seq_len(nrow(tasks)))
  worker_n <- max(1L, min(as.integer(config$workers), length(task_rows)))

  run_one <- function(task_df) {
    simulate_selection_task(task_row = task_df[1, , drop = FALSE], config = config, stan_data = stan_data, detailed = detailed)
  }

  if (worker_n <= 1L || length(task_rows) <= 1L) {
    results <- lapply(task_rows, run_one)
  } else {
    cl <- parallel::makePSOCKcluster(worker_n)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c("config", "stan_data", "task_rows", "detailed"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      setwd(config$project_root)
      library(deSolve)
      source("R/selection_strategy_utils.R")
      NULL
    })

    results <- parallel::parLapply(cl, task_rows, function(task_df) {
      simulate_selection_task(task_row = task_df[1, , drop = FALSE], config = config, stan_data = stan_data, detailed = detailed)
    })
  }

  for (i in seq_along(results)) {
    save_task_results(config = config, task_row = task_rows[[i]][1, , drop = FALSE], result = results[[i]])
  }

  saveRDS(tasks, file.path(root_dir, "task_index.Rds"))
  saveRDS(config, file.path(root_dir, "selection_strategy_config.Rds"))

  invisible(results)
}

collect_selection_outputs <- function(config) {
  root_dir <- resolve_selection_output_root(config)
  if (!dir.exists(root_dir)) {
    stop(sprintf("Selection strategy output root not found: %s", root_dir))
  }

  bind_align <- function(paths) bind_rows_align(paths, reader = readRDS)

  final_files <- list.files(root_dir, pattern = "^final_summary\\.Rds$", recursive = TRUE, full.names = TRUE)
  interval_files <- list.files(root_dir, pattern = "^interval_summary\\.Rds$", recursive = TRUE, full.names = TRUE)
  detail_files <- list.files(root_dir, pattern = "^detail\\.Rds$", recursive = TRUE, full.names = TRUE)
  meta_files <- list.files(root_dir, pattern = "^task_meta\\.Rds$", recursive = TRUE, full.names = TRUE)

  list(
    final_summary = bind_align(final_files),
    interval_summary = bind_align(interval_files),
    detail = bind_align(detail_files),
    task_meta = bind_align(meta_files)
  )
}

bind_rows_align <- function(objects, reader = identity) {
  if (!length(objects)) {
    return(data.frame())
  }

  dfs <- lapply(objects, reader)
  all_names <- unique(unlist(lapply(dfs, names)))

  dfs <- lapply(dfs, function(df) {
    missing <- setdiff(all_names, names(df))
    if (length(missing)) {
      for (nm in missing) {
        df[[nm]] <- NA
      }
    }
    df[, all_names, drop = FALSE]
  })

  do.call(rbind, dfs)
}

rerun_selected_strategies_detailed <- function(config, selected_rows, time_step_hours = NULL) {
  if (!nrow(selected_rows)) {
    return(data.frame())
  }

  stan_data <- resolve_selection_stan(config)
  out <- vector("list", nrow(selected_rows))
  saved_step <- config$detailed_time_step_hours
  if (!is.null(time_step_hours) && is.finite(time_step_hours)) {
    config$detailed_time_step_hours <- time_step_hours
  }

  for (i in seq_len(nrow(selected_rows))) {
    row <- selected_rows[i, , drop = FALSE]
    task_row <- data.frame(
      model_id = row$model_id,
      line_id = row$line_id,
      fit_context = row$fit_context,
      direction = row$direction,
      fit_type = row$fit_type,
      stringsAsFactors = FALSE
    )
    full_res <- simulate_selection_task(task_row = task_row, config = config, stan_data = stan_data, detailed = TRUE)
    detail_df <- full_res$detail
    detail_df <- detail_df[detail_df$strategy_code == row$strategy_code, , drop = FALSE]
    detail_df$selection_objective <- row$selection_objective %||% NA_character_
    out[[i]] <- detail_df
  }

  config$detailed_time_step_hours <- saved_step
  bind_rows_align(out)
}
