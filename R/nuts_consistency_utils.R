if (!exists("%>%", mode = "function")) {
  `%>%` <- magrittr::`%>%`
}

oracle_nuts_dir <- function(model_id, transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath") {
  file.path(transfer_nuts_root, model_id, "low_to_high", "oracle")
}

load_oracle_nuts_draws <- function(
  model_id,
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath",
  pattern = "^nuts_draws_.*\\.Rds$"
) {
  dir_path <- oracle_nuts_dir(model_id, transfer_nuts_root = transfer_nuts_root)
  posterior::as_draws_matrix(
    bind_rds_draws(
      list_matching_rds(dir_path, pattern = pattern),
      along = "chain"
    )
  )
}

thin_draws_matrix <- function(draws_matrix, max_draws = 1000L) {
  n <- nrow(draws_matrix)
  if (!is.finite(max_draws) || max_draws <= 0L || n <= max_draws) {
    return(draws_matrix)
  }

  idx <- unique(round(seq(1, n, length.out = max_draws)))
  draws_matrix[idx, , drop = FALSE]
}

get_model_latent_dim <- function(model_id) {
  dims <- parse_run_id(model_id)
  3L * dims$R + (dims$P - 1L) * dims$R + dims$W * dims$R + 1L
}

reconstruct_draw_parameters <- function(draw_vec, model_id, line_id, ploidy_metric) {
  dims <- parse_run_id(model_id)
  L <- get_model_latent_dim(model_id)

  reconstruct_parms(
    R = dims$R,
    P = dims$P,
    W = dims$W,
    strict_spec = (dims$C == 1L),
    M = dims$M,
    base_priors = base_priors,
    raw_theta_line = draw_vec[sprintf("raw_theta_line[%d,%d]", seq_len(L), line_id)],
    raw_theta_ploidy = draw_vec[sprintf("raw_theta_ploidy[%d]", seq_len(L))],
    ploidy_metric = ploidy_metric
  )
}

collect_oracle_nuts_qc <- function(
  model_ids,
  focus_parameters = c("ae[1]", "ah[1]", "Y_R[1]", "m"),
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath"
) {
  summary_rows <- list()
  parameter_rows <- list()

  for (model_id in model_ids) {
    dir_path <- oracle_nuts_dir(model_id, transfer_nuts_root = transfer_nuts_root)
    if (!dir.exists(dir_path)) {
      next
    }

    draw_files <- list.files(dir_path, pattern = "^nuts_draws_.*\\.Rds$", full.names = TRUE)
    if (!length(draw_files)) {
      next
    }

    draws <- load_oracle_nuts_draws(model_id, transfer_nuts_root = transfer_nuts_root)
    qc_df <- posterior::summarise_draws(
      draws,
      mean = mean,
      sd = sd,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )
    qc_df <- tibble::as_tibble(qc_df)

    summary_rows[[length(summary_rows) + 1L]] <- tibble::tibble(
      model_id = model_id,
      oracle_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE),
      n_chains = length(draw_files),
      n_draws = nrow(draws),
      n_parameters = nrow(qc_df),
      bad_rhat = sum(qc_df$rhat > 1.01, na.rm = TRUE),
      bad_ess_bulk = sum(qc_df$ess_bulk < 400, na.rm = TRUE),
      bad_ess_tail = sum(qc_df$ess_tail < 400, na.rm = TRUE),
      max_rhat = max(qc_df$rhat, na.rm = TRUE),
      min_ess_bulk = min(qc_df$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(qc_df$ess_tail, na.rm = TRUE)
    )

    parameter_rows[[length(parameter_rows) + 1L]] <- qc_df %>%
      dplyr::filter(.data$variable %in% focus_parameters) %>%
      dplyr::mutate(model_id = model_id) %>%
      dplyr::rename(parameter = "variable") %>%
      dplyr::select("model_id", "parameter", dplyr::everything())
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    parameter = dplyr::bind_rows(parameter_rows)
  )
}

collect_oracle_focus_parameter_draws <- function(
  model_ids,
  stan_data,
  focus_parameters = c("ae[1]", "ah[1]", "Y_R[1]", "m"),
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath",
  max_draws_per_model = 1000L
) {
  line_ids <- sort(unique(as.integer(stan_data$line_id)))
  rows <- list()

  for (model_id in model_ids) {
    dir_path <- oracle_nuts_dir(model_id, transfer_nuts_root = transfer_nuts_root)
    if (!dir.exists(dir_path)) {
      next
    }

    draws <- tryCatch(
      load_oracle_nuts_draws(model_id, transfer_nuts_root = transfer_nuts_root),
      error = function(e) NULL
    )
    if (is.null(draws)) {
      next
    }

    draws <- thin_draws_matrix(draws, max_draws = max_draws_per_model)
    for (line_id in line_ids) {
      states <- tryCatch(
        get_line_ploidy_states(stan_data, line_id = line_id),
        error = function(e) NULL
      )
      if (is.null(states)) {
        next
      }

      state_tbl <- tibble::tribble(
        ~ploidy_label, ~ploidy_metric,
        "low_ploidy", states$low_value,
        "high_ploidy", states$high_value
      )

      for (draw_idx in seq_len(nrow(draws))) {
        draw_vec <- draws[draw_idx, , drop = TRUE]
        names(draw_vec) <- colnames(draws)

        for (state_i in seq_len(nrow(state_tbl))) {
          parms <- reconstruct_draw_parameters(
            draw_vec = draw_vec,
            model_id = model_id,
            line_id = line_id,
            ploidy_metric = state_tbl$ploidy_metric[[state_i]]
          )

          value_map <- c(
            "ae[1]" = unname(parms$ae[[1]]),
            "ah[1]" = unname(parms$ah[[1]]),
            "Y_R[1]" = unname(parms$Y_R[[1]]),
            "m" = unname(parms$m)
          )
          keep <- intersect(names(value_map), focus_parameters)

          rows[[length(rows) + 1L]] <- tibble::tibble(
            model_id = model_id,
            line_id = line_id,
            ploidy_label = state_tbl$ploidy_label[[state_i]],
            ploidy_metric = state_tbl$ploidy_metric[[state_i]],
            .draw = draw_idx,
            parameter = keep,
            value = as.numeric(value_map[keep])
          )
        }
      }
    }
  }

  dplyr::bind_rows(rows)
}

collect_oracle_focus_raw_theta_ploidy_draws <- function(
  model_ids,
  focus_parameters = c("ae[1]", "ah[1]", "Y_R[1]", "m"),
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath",
  max_draws_per_model = 1000L
) {
  rows <- list()

  for (model_id in model_ids) {
    dir_path <- oracle_nuts_dir(model_id, transfer_nuts_root = transfer_nuts_root)
    if (!dir.exists(dir_path)) {
      next
    }

    draws <- tryCatch(
      load_oracle_nuts_draws(model_id, transfer_nuts_root = transfer_nuts_root),
      error = function(e) NULL
    )
    if (is.null(draws)) {
      next
    }

    dims <- parse_run_id(model_id)
    param_labels <- get_param_names(dims$R, dims$P, dims$W, dims$C, dims$M)
    keep_idx <- which(param_labels %in% focus_parameters)
    if (!length(keep_idx)) {
      next
    }

    draws <- thin_draws_matrix(draws, max_draws = max_draws_per_model)
    coef_names <- sprintf("raw_theta_ploidy[%d]", keep_idx)

    for (draw_idx in seq_len(nrow(draws))) {
      draw_vec <- draws[draw_idx, coef_names, drop = TRUE]
      names(draw_vec) <- param_labels[keep_idx]

      rows[[length(rows) + 1L]] <- tibble::tibble(
        model_id = model_id,
        .draw = draw_idx,
        parameter = names(draw_vec),
        value = unname(as.numeric(draw_vec))
      )
    }
  }

  dplyr::bind_rows(rows)
}

summarize_focus_parameter_posteriors <- function(focus_draw_df) {
  if (is.null(focus_draw_df) || !nrow(focus_draw_df)) {
    return(tibble::tibble())
  }

  focus_draw_df %>%
    dplyr::group_by(.data$model_id, .data$line_id, .data$ploidy_label, .data$parameter) %>%
    dplyr::summarise(
      n_draws = dplyr::n(),
      median = stats::median(.data$value, na.rm = TRUE),
      mean = mean(.data$value, na.rm = TRUE),
      sd = stats::sd(.data$value, na.rm = TRUE),
      q25 = stats::quantile(.data$value, probs = 0.25, na.rm = TRUE, names = FALSE),
      q75 = stats::quantile(.data$value, probs = 0.75, na.rm = TRUE, names = FALSE),
      q025 = stats::quantile(.data$value, probs = 0.025, na.rm = TRUE, names = FALSE),
      q975 = stats::quantile(.data$value, probs = 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

compute_instantaneous_net_growth <- function(parms, glucose, aux_resource_value = 1.0) {
  R_vec <- c(glucose, rep(aux_resource_value, max(parms$R - 1L, 0L)))
  R_vec <- R_vec[seq_len(parms$R)]
  u <- (parms$ae * R_vec) / (parms$ah + R_vec)
  Phi <- parms$A_mat %*% (parms$Y_R * u)
  mu_growth <- smooth_min_vec(Phi)
  as.numeric(mu_growth - parms$m)
}

compute_instantaneous_net_growth_with_resource2 <- function(
  parms,
  glucose,
  resource2_value,
  extra_resource_value = 1.0
) {
  if (parms$R <= 1L) {
    R_vec <- glucose
  } else {
    R_vec <- c(glucose, resource2_value, rep(extra_resource_value, max(parms$R - 2L, 0L)))
  }
  R_vec <- R_vec[seq_len(parms$R)]
  u <- (parms$ae * R_vec) / (parms$ah + R_vec)
  Phi <- parms$A_mat %*% (parms$Y_R * u)
  mu_growth <- smooth_min_vec(Phi)
  as.numeric(mu_growth - parms$m)
}

collect_oracle_growth_curve_draws <- function(
  model_ids,
  stan_data,
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath",
  glucose_grid = NULL,
  max_draws_per_model = 1000L,
  aux_resource_value = 1.0
) {
  if (is.null(glucose_grid)) {
    glucose_grid <- sort(unique(as.numeric(stan_data$G0_per_well)))
  }
  glucose_grid <- glucose_grid[is.finite(glucose_grid)]

  line_ids <- sort(unique(as.integer(stan_data$line_id)))
  rows <- list()
  idx <- 1L

  for (model_id in model_ids) {
    draws <- tryCatch(
      load_oracle_nuts_draws(model_id, transfer_nuts_root = transfer_nuts_root),
      error = function(e) NULL
    )
    if (is.null(draws)) {
      next
    }

    draws <- thin_draws_matrix(draws, max_draws = max_draws_per_model)

    for (line_id in line_ids) {
      states <- tryCatch(
        get_line_ploidy_states(stan_data, line_id = line_id),
        error = function(e) NULL
      )
      if (is.null(states)) {
        next
      }

      state_tbl <- tibble::tribble(
        ~ploidy_label, ~ploidy_metric,
        "low_ploidy", states$low_value,
        "high_ploidy", states$high_value
      )

      for (draw_idx in seq_len(nrow(draws))) {
        draw_vec <- draws[draw_idx, , drop = TRUE]
        names(draw_vec) <- colnames(draws)

        for (state_i in seq_len(nrow(state_tbl))) {
          parms <- reconstruct_draw_parameters(
            draw_vec = draw_vec,
            model_id = model_id,
            line_id = line_id,
            ploidy_metric = state_tbl$ploidy_metric[[state_i]]
          )

          growth_vals <- vapply(
            glucose_grid,
            function(g) compute_instantaneous_net_growth(parms, glucose = g, aux_resource_value = aux_resource_value),
            numeric(1)
          )

          rows[[idx]] <- tibble::tibble(
            model_id = model_id,
            line_id = line_id,
            ploidy_label = state_tbl$ploidy_label[[state_i]],
            ploidy_metric = state_tbl$ploidy_metric[[state_i]],
            .draw = draw_idx,
            glucose = glucose_grid,
            growth_rate = growth_vals
          )
          idx <- idx + 1L
        }
      }
    }
  }

  dplyr::bind_rows(rows)
}

summarize_model_growth_curves <- function(curve_draw_df) {
  if (is.null(curve_draw_df) || !nrow(curve_draw_df)) {
    return(tibble::tibble())
  }

  curve_draw_df %>%
    dplyr::group_by(.data$model_id, .data$line_id, .data$ploidy_label, .data$glucose) %>%
    dplyr::summarise(
      median_growth = stats::median(.data$growth_rate, na.rm = TRUE),
      q25_growth = stats::quantile(.data$growth_rate, probs = 0.25, na.rm = TRUE, names = FALSE),
      q75_growth = stats::quantile(.data$growth_rate, probs = 0.75, na.rm = TRUE, names = FALSE),
      q025_growth = stats::quantile(.data$growth_rate, probs = 0.025, na.rm = TRUE, names = FALSE),
      q975_growth = stats::quantile(.data$growth_rate, probs = 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

summarize_pooled_growth_curves <- function(model_curve_df) {
  if (is.null(model_curve_df) || !nrow(model_curve_df)) {
    return(tibble::tibble())
  }

  model_curve_df %>%
    dplyr::group_by(.data$line_id, .data$ploidy_label, .data$glucose) %>%
    dplyr::summarise(
      n_models = dplyr::n_distinct(.data$model_id),
      pooled_median_growth = stats::median(.data$median_growth, na.rm = TRUE),
      pooled_q25_growth = stats::quantile(.data$median_growth, probs = 0.25, na.rm = TRUE, names = FALSE),
      pooled_q75_growth = stats::quantile(.data$median_growth, probs = 0.75, na.rm = TRUE, names = FALSE),
      pooled_q025_growth = stats::quantile(.data$median_growth, probs = 0.025, na.rm = TRUE, names = FALSE),
      pooled_q975_growth = stats::quantile(.data$median_growth, probs = 0.975, na.rm = TRUE, names = FALSE),
      mean_within_model_band_width = mean(.data$q975_growth - .data$q025_growth, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_growth_curve_consistency <- function(model_curve_df) {
  if (is.null(model_curve_df) || !nrow(model_curve_df)) {
    return(tibble::tibble())
  }

  contrast_df <- model_curve_df %>%
    dplyr::select("model_id", "line_id", "ploidy_label", "glucose", "median_growth") %>%
    tidyr::pivot_wider(names_from = "ploidy_label", values_from = "median_growth")

  base_df <- model_curve_df %>%
    dplyr::group_by(.data$line_id, .data$ploidy_label, .data$glucose) %>%
    dplyr::summarise(
      n_models = dplyr::n_distinct(.data$model_id),
      between_model_sd = stats::sd(.data$median_growth, na.rm = TRUE),
      model_median_range = diff(range(.data$median_growth, na.rm = TRUE)),
      .groups = "drop"
    )

  contrast_summary <- contrast_df %>%
    dplyr::mutate(high_minus_low = .data$high_ploidy - .data$low_ploidy) %>%
    dplyr::group_by(.data$line_id, .data$glucose) %>%
    dplyr::summarise(
      ploidy_sign_agreement = max(mean(.data$high_minus_low >= 0, na.rm = TRUE), mean(.data$high_minus_low <= 0, na.rm = TRUE)),
      median_high_minus_low = stats::median(.data$high_minus_low, na.rm = TRUE),
      .groups = "drop"
    )

  base_df %>%
    dplyr::left_join(contrast_summary, by = c("line_id", "glucose"))
}

summarize_oracle_growth_delta_surface <- function(
  model_ids,
  stan_data,
  transfer_nuts_root = "data/gpath_transfer_cv_nuts/gpath",
  glucose_grid = NULL,
  resource2_grid = NULL,
  max_draws_per_model = 200L,
  extra_resource_value = 1.0
) {
  if (is.null(glucose_grid)) {
    glucose_grid <- sort(unique(as.numeric(stan_data$G0_per_well)))
  }
  glucose_grid <- glucose_grid[is.finite(glucose_grid)]
  if (is.null(resource2_grid)) {
    resource2_grid <- seq(0, 1, length.out = length(glucose_grid))
  }
  resource2_grid <- resource2_grid[is.finite(resource2_grid)]

  grid_df <- expand.grid(
    glucose = glucose_grid,
    resource2 = resource2_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  line_ids <- sort(unique(as.integer(stan_data$line_id)))
  rows <- list()
  idx <- 1L

  for (model_id in model_ids) {
    draws <- tryCatch(
      load_oracle_nuts_draws(model_id, transfer_nuts_root = transfer_nuts_root),
      error = function(e) NULL
    )
    if (is.null(draws)) {
      next
    }

    draws <- thin_draws_matrix(draws, max_draws = max_draws_per_model)

    for (line_id in line_ids) {
      states <- tryCatch(
        get_line_ploidy_states(stan_data, line_id = line_id),
        error = function(e) NULL
      )
      if (is.null(states)) {
        next
      }

      delta_store <- matrix(NA_real_, nrow = nrow(draws), ncol = nrow(grid_df))
      for (draw_idx in seq_len(nrow(draws))) {
        draw_vec <- draws[draw_idx, , drop = TRUE]
        names(draw_vec) <- colnames(draws)

        parms_low <- reconstruct_draw_parameters(
          draw_vec = draw_vec,
          model_id = model_id,
          line_id = line_id,
          ploidy_metric = states$low_value
        )
        parms_high <- reconstruct_draw_parameters(
          draw_vec = draw_vec,
          model_id = model_id,
          line_id = line_id,
          ploidy_metric = states$high_value
        )

        delta_store[draw_idx, ] <- vapply(
          seq_len(nrow(grid_df)),
          function(i) {
            growth_high <- compute_instantaneous_net_growth_with_resource2(
              parms_high,
              glucose = grid_df$glucose[[i]],
              resource2_value = grid_df$resource2[[i]],
              extra_resource_value = extra_resource_value
            )
            growth_low <- compute_instantaneous_net_growth_with_resource2(
              parms_low,
              glucose = grid_df$glucose[[i]],
              resource2_value = grid_df$resource2[[i]],
              extra_resource_value = extra_resource_value
            )
            growth_high - growth_low
          },
          numeric(1)
        )
      }

      rows[[idx]] <- tibble::tibble(
        model_id = model_id,
        line_id = line_id,
        glucose = grid_df$glucose,
        resource2 = grid_df$resource2,
        n_draws = nrow(draws),
        median_delta_growth = apply(delta_store, 2L, stats::median, na.rm = TRUE),
        mean_delta_growth = apply(delta_store, 2L, mean, na.rm = TRUE),
        q025_delta_growth = apply(delta_store, 2L, stats::quantile, probs = 0.025, na.rm = TRUE, names = FALSE),
        q25_delta_growth = apply(delta_store, 2L, stats::quantile, probs = 0.25, na.rm = TRUE, names = FALSE),
        q75_delta_growth = apply(delta_store, 2L, stats::quantile, probs = 0.75, na.rm = TRUE, names = FALSE),
        q975_delta_growth = apply(delta_store, 2L, stats::quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(rows)
}
