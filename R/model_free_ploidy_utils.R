source("R/project_paths.R")
source("R/gpath_run_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

safe_trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

calc_rate_extrema <- function(hours, value, log1p_scale = FALSE) {
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (length(hours) < 2L) {
    return(list(
      max_rate = NA_real_,
      min_rate = NA_real_,
      time_at_max_rate = NA_real_,
      time_at_min_rate = NA_real_
    ))
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  x <- if (log1p_scale) log1p(pmax(value, 0)) else value
  dt <- diff(hours)
  dy <- diff(x)

  keep <- is.finite(dt) & dt > 0 & is.finite(dy)
  dt <- dt[keep]
  dy <- dy[keep]

  if (!length(dt)) {
    return(list(
      max_rate = NA_real_,
      min_rate = NA_real_,
      time_at_max_rate = NA_real_,
      time_at_min_rate = NA_real_
    ))
  }

  rate <- dy / dt
  mids <- (hours[-1L] + hours[-length(hours)]) / 2
  mids <- mids[keep]

  idx_max <- which.max(rate)
  idx_min <- which.min(rate)

  list(
    max_rate = rate[idx_max],
    min_rate = rate[idx_min],
    time_at_max_rate = mids[idx_max],
    time_at_min_rate = mids[idx_min]
  )
}

time_to_threshold <- function(hours, value, threshold, direction = c("below", "above")) {
  direction <- match.arg(direction)
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (!length(hours)) {
    return(NA_real_)
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  hit_idx <- if (direction == "below") {
    which(value <= threshold)
  } else {
    which(value >= threshold)
  }

  if (!length(hit_idx)) {
    return(NA_real_)
  }

  hours[min(hit_idx)]
}

build_model_free_tables <- function(stan_data_path = NULL) {
  stan_data_path <- resolve_stan_data_path(stan_data_path)
  stan_data <- readRDS(stan_data_path)

  line_names <- names(stan_data$line_map)[match(stan_data$line_id, unname(stan_data$line_map))]

  well_meta <- tibble(
    well_idx = seq_len(stan_data$N_wells),
    line_id = as.integer(stan_data$line_id),
    cellLine = line_names,
    ploidy_metric = as.numeric(stan_data$ploidy_metric),
    ploidy_abs = as.numeric(stan_data$ploidy_abs),
    G0 = as.numeric(stan_data$G0_per_well),
    exp_id = as.integer(stan_data$exp_id),
    has_starvation = as.integer(stan_data$has_starvation)
  )

  count_obs <- tibble(
    well_idx = as.integer(stan_data$well_idx_count),
    hours = as.numeric(stan_data$t_grid[stan_data$grid_idx_count]),
    rep_id = as.character(stan_data$rep_id_count),
    total_cells = as.numeric(stan_data$N_obs),
    dead_cells = as.numeric(stan_data$D_obs)
  ) %>%
    mutate(
      live_cells = pmax(total_cells - dead_cells, 0),
      dead_fraction = if_else(total_cells > 0, dead_cells / total_cells, NA_real_)
    )

  count_summary <- count_obs %>%
    group_by(well_idx, hours) %>%
    summarise(
      n_count_reps = n(),
      live_cells = median(live_cells, na.rm = TRUE),
      total_cells = median(total_cells, na.rm = TRUE),
      dead_cells = median(dead_cells, na.rm = TRUE),
      dead_fraction = median(dead_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(well_meta, by = "well_idx")

  gluc_obs <- tibble(
    well_idx = as.integer(stan_data$well_idx_gluc),
    hours = as.numeric(stan_data$t_grid[stan_data$grid_idx_gluc]),
    lum_obs = as.numeric(stan_data$lum_obs),
    dilution = as.numeric(stan_data$dilution),
    is_censored = as.integer(stan_data$is_censored)
  ) %>%
    mutate(
      exp_id = well_meta$exp_id[well_idx],
      calib_a = as.numeric(stan_data$calib_a_fixed[exp_id]),
      calib_b = as.numeric(stan_data$calib_b_fixed[exp_id]),
      glucose_hat = pmax(0, (lum_obs - calib_b) / (calib_a * pmax(dilution, 1e-12)))
    )

  gluc_summary <- gluc_obs %>%
    group_by(well_idx, hours) %>%
    summarise(
      n_glucose_reps = n(),
      glucose_hat = {
        uncens <- glucose_hat[is_censored == 0L]
        if (length(uncens)) {
          median(uncens, na.rm = TRUE)
        } else {
          median(glucose_hat, na.rm = TRUE)
        }
      },
      glucose_hat_lo = quantile(glucose_hat, probs = 0.25, na.rm = TRUE, names = FALSE),
      glucose_hat_hi = quantile(glucose_hat, probs = 0.75, na.rm = TRUE, names = FALSE),
      any_glucose_censored = as.integer(any(is_censored == 1L)),
      .groups = "drop"
    ) %>%
    left_join(well_meta, by = "well_idx")

  list(
    stan_data_path = stan_data_path,
    stan_data = stan_data,
    well_meta = well_meta,
    count_obs = count_obs,
    count_summary = count_summary,
    glucose_obs = gluc_obs,
    glucose_summary = gluc_summary
  )
}

summarize_one_well_features <- function(well_meta_row, count_df, glucose_df, glucose_floor = 0.1) {
  count_df <- count_df %>% arrange(hours)
  glucose_df <- glucose_df %>% arrange(hours)

  live_rate <- calc_rate_extrema(count_df$hours, count_df$live_cells, log1p_scale = TRUE)
  dead_rate <- calc_rate_extrema(count_df$hours, count_df$dead_cells, log1p_scale = TRUE)
  dead_frac_rate <- calc_rate_extrema(count_df$hours, count_df$dead_fraction, log1p_scale = FALSE)
  glucose_rate <- calc_rate_extrema(glucose_df$hours, glucose_df$glucose_hat, log1p_scale = FALSE)

  peak_live_idx <- if (nrow(count_df)) which.max(count_df$live_cells) else integer(0)
  nadir_live_idx <- if (nrow(count_df)) which.min(count_df$live_cells) else integer(0)
  peak_dead_idx <- if (nrow(count_df)) which.max(count_df$dead_cells) else integer(0)
  min_glucose_idx <- if (nrow(glucose_df)) which.min(glucose_df$glucose_hat) else integer(0)

  out <- tibble(
    well_idx = well_meta_row$well_idx,
    cellLine = well_meta_row$cellLine,
    line_id = well_meta_row$line_id,
    ploidy_metric = well_meta_row$ploidy_metric,
    ploidy_abs = well_meta_row$ploidy_abs,
    G0 = well_meta_row$G0,
    exp_id = well_meta_row$exp_id,
    has_starvation = well_meta_row$has_starvation,
    n_count_times = nrow(count_df),
    n_glucose_times = nrow(glucose_df),
    live_initial = if (nrow(count_df)) count_df$live_cells[1] else NA_real_,
    live_final = if (nrow(count_df)) dplyr::last(count_df$live_cells) else NA_real_,
    live_peak = if (length(peak_live_idx)) count_df$live_cells[peak_live_idx] else NA_real_,
    live_nadir = if (length(nadir_live_idx)) count_df$live_cells[nadir_live_idx] else NA_real_,
    live_peak_time = if (length(peak_live_idx)) count_df$hours[peak_live_idx] else NA_real_,
    live_nadir_time = if (length(nadir_live_idx)) count_df$hours[nadir_live_idx] else NA_real_,
    live_net_change = if (nrow(count_df)) dplyr::last(count_df$live_cells) - count_df$live_cells[1] else NA_real_,
    live_fold_change = if (nrow(count_df)) (dplyr::last(count_df$live_cells) + 1) / (count_df$live_cells[1] + 1) else NA_real_,
    live_auc = safe_trapz(count_df$hours, count_df$live_cells),
    max_growth_rate = live_rate$max_rate,
    max_decline_rate = live_rate$min_rate,
    time_max_growth_rate = live_rate$time_at_max_rate,
    time_max_decline_rate = live_rate$time_at_min_rate,
    dead_initial = if (nrow(count_df)) count_df$dead_cells[1] else NA_real_,
    dead_final = if (nrow(count_df)) dplyr::last(count_df$dead_cells) else NA_real_,
    dead_peak = if (length(peak_dead_idx)) count_df$dead_cells[peak_dead_idx] else NA_real_,
    dead_peak_time = if (length(peak_dead_idx)) count_df$hours[peak_dead_idx] else NA_real_,
    dead_net_change = if (nrow(count_df)) dplyr::last(count_df$dead_cells) - count_df$dead_cells[1] else NA_real_,
    dead_auc = safe_trapz(count_df$hours, count_df$dead_cells),
    max_death_rate = dead_rate$max_rate,
    dead_fraction_initial = if (nrow(count_df)) count_df$dead_fraction[1] else NA_real_,
    dead_fraction_final = if (nrow(count_df)) dplyr::last(count_df$dead_fraction) else NA_real_,
    dead_fraction_peak = if (nrow(count_df)) max(count_df$dead_fraction, na.rm = TRUE) else NA_real_,
    max_dead_fraction_rate = dead_frac_rate$max_rate,
    glucose_initial = if (nrow(glucose_df)) glucose_df$glucose_hat[1] else NA_real_,
    glucose_final = if (nrow(glucose_df)) dplyr::last(glucose_df$glucose_hat) else NA_real_,
    glucose_min = if (length(min_glucose_idx)) glucose_df$glucose_hat[min_glucose_idx] else NA_real_,
    glucose_min_time = if (length(min_glucose_idx)) glucose_df$hours[min_glucose_idx] else NA_real_,
    glucose_drawdown = if (nrow(glucose_df)) glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat) else NA_real_,
    glucose_drawdown_fraction = if (nrow(glucose_df) && glucose_df$glucose_hat[1] > 0) {
      (glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat)) / glucose_df$glucose_hat[1]
    } else {
      NA_real_
    },
    glucose_auc = safe_trapz(glucose_df$hours, glucose_df$glucose_hat),
    max_glucose_drawdown_rate = if (is.finite(glucose_rate$min_rate)) -glucose_rate$min_rate else NA_real_,
    time_max_glucose_drawdown = glucose_rate$time_at_min_rate,
    time_to_glucose_floor = time_to_threshold(glucose_df$hours, glucose_df$glucose_hat, glucose_floor, direction = "below"),
    live_gain_per_glucose = if (nrow(count_df) && nrow(glucose_df)) {
      (dplyr::last(count_df$live_cells) - count_df$live_cells[1]) / max(glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat), 1e-8)
    } else {
      NA_real_
    },
    dead_gain_per_glucose = if (nrow(count_df) && nrow(glucose_df)) {
      (dplyr::last(count_df$dead_cells) - count_df$dead_cells[1]) / max(glucose_df$glucose_hat[1] - dplyr::last(glucose_df$glucose_hat), 1e-8)
    } else {
      NA_real_
    }
  )

  out
}

build_feature_panel <- function(stan_data_path = NULL, glucose_floor = 0.1) {
  tables <- build_model_free_tables(stan_data_path = stan_data_path)

  feature_panel <- tables$well_meta %>%
    split(.$well_idx) %>%
    map_dfr(function(meta_row) {
      well_idx <- meta_row$well_idx[[1]]
      count_df <- tables$count_summary %>% filter(well_idx == !!well_idx)
      glucose_df <- tables$glucose_summary %>% filter(well_idx == !!well_idx)
      summarize_one_well_features(meta_row[1, ], count_df, glucose_df, glucose_floor = glucose_floor)
    }) %>%
    arrange(cellLine, ploidy_metric, G0)

  list(
    stan_data_path = tables$stan_data_path,
    well_meta = tables$well_meta,
    count_summary = tables$count_summary,
    glucose_summary = tables$glucose_summary,
    feature_panel = feature_panel
  )
}

compute_empirical_effects <- function(feature_panel) {
  feature_cols <- names(feature_panel)[vapply(feature_panel, is.numeric, logical(1))]
  feature_cols <- setdiff(
    feature_cols,
    c("well_idx", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id", "has_starvation", "n_count_times", "n_glucose_times")
  )

  paired <- feature_panel %>%
    group_by(cellLine, G0) %>%
    filter(n() == 2L) %>%
    arrange(ploidy_metric, .by_group = TRUE) %>%
    summarise(
      observed_ploidy = first(ploidy_metric),
      holdout_ploidy = last(ploidy_metric),
      delta_ploidy = last(ploidy_metric) - first(ploidy_metric),
      across(all_of(feature_cols), ~ (dplyr::last(.) - dplyr::first(.)), .names = "{.col}_delta"),
      .groups = "drop"
    )

  long <- paired %>%
    pivot_longer(
      cols = ends_with("_delta"),
      names_to = "feature",
      values_to = "delta_value"
    ) %>%
    mutate(
      feature = sub("_delta$", "", feature),
      effect_per_ploidy = delta_value / delta_ploidy
    )

  effect_matrix <- long %>%
    select(cellLine, G0, feature, effect_per_ploidy) %>%
    tidyr::pivot_wider(names_from = feature, values_from = effect_per_ploidy) %>%
    arrange(cellLine, G0)

  pca <- NULL
  score_df <- tibble()
  loading_df <- tibble()

  mat <- effect_matrix %>%
    ungroup() %>%
    select(-cellLine, -G0) %>%
    as.matrix()

  keep_cols <- apply(mat, 2, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  if (nrow(mat) >= 3L && sum(keep_cols) >= 2L) {
    mat_sub <- mat[, keep_cols, drop = FALSE]
    pca <- prcomp(mat_sub, center = TRUE, scale. = TRUE)
    score_df <- bind_cols(
      effect_matrix %>% select(cellLine, G0),
      as_tibble(pca$x, .name_repair = "minimal")
    )
    loading_df <- tibble(
      feature = rownames(pca$rotation),
      PC1 = pca$rotation[, 1],
      PC2 = if (ncol(pca$rotation) >= 2L) pca$rotation[, 2] else NA_real_
    )
  }

  list(
    feature_cols = feature_cols,
    paired_effects_wide = paired,
    paired_effects_long = long,
    effect_matrix = effect_matrix,
    pca = pca,
    pca_scores = score_df,
    pca_loadings = loading_df
  )
}

evaluate_transfer_predictions <- function(feature_panel) {
  feature_cols <- names(feature_panel)[vapply(feature_panel, is.numeric, logical(1))]
  feature_cols <- setdiff(
    feature_cols,
    c("well_idx", "line_id", "ploidy_metric", "ploidy_abs", "G0", "exp_id", "has_starvation", "n_count_times", "n_glucose_times")
  )

  rows <- list()
  idx <- 1L

  for (line_name in sort(unique(feature_panel$cellLine))) {
    line_df <- feature_panel %>% filter(cellLine == line_name)
    states <- sort(unique(line_df$ploidy_metric))

    if (length(states) != 2L) {
      next
    }

    for (direction in c("low_to_high", "high_to_low")) {
      observed_ploidy <- if (direction == "low_to_high") states[1] else states[2]
      holdout_ploidy <- if (direction == "low_to_high") states[2] else states[1]
      delta_ploidy <- holdout_ploidy - observed_ploidy

      obs_df <- feature_panel %>%
        filter(cellLine == line_name, abs(ploidy_metric - observed_ploidy) < 1e-12) %>%
        arrange(G0)

      hold_df <- feature_panel %>%
        filter(cellLine == line_name, abs(ploidy_metric - holdout_ploidy) < 1e-12) %>%
        arrange(G0)

      train_pairs <- feature_panel %>%
        filter(cellLine != line_name) %>%
        group_by(cellLine, G0) %>%
        filter(n() == 2L) %>%
        summarise(
          low_ploidy = min(ploidy_metric),
          high_ploidy = max(ploidy_metric),
          delta_ploidy_train = high_ploidy - low_ploidy,
          across(all_of(feature_cols), ~ .[which.min(ploidy_metric)], .names = "low__{.col}"),
          across(all_of(feature_cols), ~ .[which.max(ploidy_metric)], .names = "high__{.col}"),
          .groups = "drop"
        )

      for (g0 in sort(unique(obs_df$G0))) {
        obs_row <- obs_df %>% filter(abs(G0 - g0) < 1e-12)
        hold_row <- hold_df %>% filter(abs(G0 - g0) < 1e-12)
        train_row <- train_pairs %>% filter(abs(G0 - g0) < 1e-12)

        if (!nrow(obs_row) || !nrow(hold_row) || !nrow(train_row)) {
          next
        }

        for (feature in feature_cols) {
          obs_val <- obs_row[[feature]][1]
          true_val <- hold_row[[feature]][1]

          low_col <- paste0("low__", feature)
          high_col <- paste0("high__", feature)
          low_vals <- train_row[[low_col]]
          high_vals <- train_row[[high_col]]
          slope_vals <- (high_vals - low_vals) / train_row$delta_ploidy_train

          effect_hat <- median(slope_vals[is.finite(slope_vals)], na.rm = TRUE)
          if (!is.finite(effect_hat)) {
            next
          }

          pred_null <- obs_val
          pred_transfer <- obs_val + effect_hat * delta_ploidy
          true_effect <- (true_val - obs_val) / delta_ploidy

          rows[[idx]] <- tibble(
            cellLine = line_name,
            G0 = g0,
            direction = direction,
            feature = feature,
            observed_ploidy = observed_ploidy,
            holdout_ploidy = holdout_ploidy,
            delta_ploidy = delta_ploidy,
            observed_value = obs_val,
            true_value = true_val,
            true_effect = true_effect,
            transfer_effect = effect_hat,
            pred_null = pred_null,
            pred_transfer = pred_transfer,
            abs_err_null = abs(pred_null - true_val),
            abs_err_transfer = abs(pred_transfer - true_val),
            sq_err_null = (pred_null - true_val) ^ 2,
            sq_err_transfer = (pred_transfer - true_val) ^ 2,
            abs_err_improvement = abs(pred_null - true_val) - abs(pred_transfer - true_val),
            sign_match = as.integer(sign(true_effect) == sign(effect_hat)),
            n_training_lines = sum(is.finite(slope_vals))
          )
          idx <- idx + 1L
        }
      }
    }
  }

  predictions <- bind_rows(rows)

  feature_summary <- predictions %>%
    group_by(feature) %>%
    summarise(
      n_cases = n(),
      mae_null = mean(abs_err_null, na.rm = TRUE),
      mae_transfer = mean(abs_err_transfer, na.rm = TRUE),
      rmse_null = sqrt(mean(sq_err_null, na.rm = TRUE)),
      rmse_transfer = sqrt(mean(sq_err_transfer, na.rm = TRUE)),
      mean_abs_err_improvement = mean(abs_err_improvement, na.rm = TRUE),
      median_abs_err_improvement = median(abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_err_improvement), mae_transfer)

  case_summary <- predictions %>%
    group_by(cellLine, direction, G0) %>%
    summarise(
      n_features = n(),
      mean_abs_err_improvement = mean(abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cellLine, direction, G0)

  list(
    feature_cols = feature_cols,
    predictions = predictions,
    feature_summary = feature_summary,
    case_summary = case_summary
  )
}
