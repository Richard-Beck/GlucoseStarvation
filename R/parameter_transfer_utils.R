source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/elpd_transfer_utils.R")
source("R/optim_utils.R")

flatten_numeric_object <- function(x, prefix = NULL) {
  out <- numeric()

  if (is.null(x)) {
    return(out)
  }

  if (is.numeric(x) && is.null(dim(x)) && length(x) == 1L) {
    nm <- prefix %||% ""
    out[nm] <- as.numeric(x)
    return(out)
  }

  if (is.vector(x) && is.numeric(x)) {
    idx <- seq_along(x)
    nm <- sprintf("%s[%d]", prefix, idx)
    out[nm] <- as.numeric(x)
    return(out)
  }

  if (is.matrix(x)) {
    for (i in seq_len(nrow(x))) {
      for (j in seq_len(ncol(x))) {
        nm <- sprintf("%s[%d,%d]", prefix, i, j)
        out[nm] <- as.numeric(x[i, j])
      }
    }
    return(out)
  }

  if (is.list(x)) {
    for (nm in names(x)) {
      out <- c(out, flatten_numeric_object(x[[nm]], prefix = nm))
    }
  }

  out
}

flatten_gpath_reconstructed_parms <- function(parms) {
  out <- numeric()

  if (is.null(parms) || !is.list(parms)) {
    return(out)
  }

  if (!is.null(parms$R) && length(parms$R) == 1L) out["R"] <- as.numeric(parms$R)
  if (!is.null(parms$P) && length(parms$P) == 1L) out["P"] <- as.numeric(parms$P)
  if (!is.null(parms$W) && length(parms$W) == 1L) out["W"] <- as.numeric(parms$W)

  if (!is.null(parms$ae)) {
    out[sprintf("ae[%d]", seq_along(parms$ae))] <- as.numeric(parms$ae)
  }
  if (!is.null(parms$ah)) {
    out[sprintf("ah[%d]", seq_along(parms$ah))] <- as.numeric(parms$ah)
  }
  if (!is.null(parms$Y_R)) {
    out[sprintf("Y_R[%d]", seq_along(parms$Y_R))] <- as.numeric(parms$Y_R)
  }

  if (!is.null(parms$A_mat) && is.matrix(parms$A_mat)) {
    for (i in seq_len(nrow(parms$A_mat))) {
      for (j in seq_len(ncol(parms$A_mat))) {
        out[sprintf("A_mat[%d,%d]", i, j)] <- as.numeric(parms$A_mat[i, j])
      }
    }
  }

  if (!is.null(parms$K_mat) && is.matrix(parms$K_mat) && nrow(parms$K_mat) > 0L && ncol(parms$K_mat) > 0L) {
    for (i in seq_len(nrow(parms$K_mat))) {
      for (j in seq_len(ncol(parms$K_mat))) {
        out[sprintf("K_mat[%d,%d]", i, j)] <- as.numeric(parms$K_mat[i, j])
      }
    }
  }

  if (!is.null(parms$m) && length(parms$m) == 1L) {
    out["m"] <- as.numeric(parms$m)
  }

  if (!is.null(parms$waste_mech)) {
    out[sprintf("waste_mech[%d]", seq_along(parms$waste_mech))] <- as.numeric(parms$waste_mech)
  }

  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

load_best_transfer_fit <- function(
  model_id,
  line_id,
  direction,
  fit_type,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  best_path <- file.path(output_root, model_name, model_id, "transfer_best_start_summary.Rds")
  if (!file.exists(best_path)) {
    stop(sprintf("Best-start summary not found: %s", best_path))
  }

  best_df <- readRDS(best_path)
  direction <- normalize_transfer_direction(direction)
  fit_type <- normalize_fit_type(fit_type)

  hit <- best_df[
    best_df$line_id == line_id &
      best_df$direction == direction &
      best_df$fit_type == fit_type,
  ]

  if (nrow(hit) != 1L) {
    stop(sprintf(
      "Expected exactly one best-start row for model=%s line=%d direction=%s fit=%s",
      model_id,
      line_id,
      direction,
      fit_type
    ))
  }

  run_tag <- hit$run_tag[[1]]
  fit_dir <- file.path(output_root, model_name, model_id, direction, fit_type)
  draws_path <- file.path(fit_dir, sprintf("optim_draws_%s.Rds", run_tag))
  meta_path <- file.path(fit_dir, sprintf("split_meta_%s.Rds", run_tag))

  list(
    summary = hit,
    draws = readRDS(draws_path),
    split_meta = readRDS(meta_path),
    run_tag = run_tag
  )
}

reconstruct_line_state_parameters <- function(draw_vec, model_id, line_id, ploidy_metric, ploidy_effect_mask = NULL) {
  dims <- parse_run_id(model_id)
  R <- dims$R
  P <- dims$P
  W <- dims$W
  strict_spec <- (dims$C == 1)
  M <- dims$M
  L <- 3 * R + (P - 1) * R + W * R + 1

  raw_theta_line <- draw_vec[sprintf("raw_theta_line[%d,%d]", 1:L, line_id)]
  raw_theta_ploidy <- draw_vec[sprintf("raw_theta_ploidy[%d]", 1:L)]

  parms <- reconstruct_parms(
    R = R,
    P = P,
    W = W,
    strict_spec = strict_spec,
    M = M,
    base_priors = base_priors,
    raw_theta_line = raw_theta_line,
    raw_theta_ploidy = raw_theta_ploidy,
    ploidy_metric = ploidy_metric,
    ploidy_effect_mask = ploidy_effect_mask
  )

  flatten_gpath_reconstructed_parms(parms)
}

build_parameter_transfer_tables <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  stan_data_path = "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds",
  model_name = "gpath"
) {
  best_path <- file.path(output_root, model_name, model_id, "transfer_best_start_summary.Rds")
  if (!file.exists(best_path)) {
    stop(sprintf("Best-start summary not found: %s", best_path))
  }

  best_df <- readRDS(best_path)
  best_df <- best_df[best_df$fit_type %in% c("transfer", "oracle"), ]
  if (!nrow(best_df)) {
    stop("No transfer/oracle rows found in best-start summary")
  }

  stan_data <- readRDS(resolve_stan_data_path(stan_data_path))

  state_rows <- list()
  shift_rows <- list()
  compare_rows <- list()
  summary_rows <- list()
  idx_state <- 1L
  idx_shift <- 1L
  idx_compare <- 1L
  idx_summary <- 1L

  keys <- unique(best_df[, c("line_id", "direction")])

  for (k in seq_len(nrow(keys))) {
    line_id <- keys$line_id[k]
    direction <- keys$direction[k]
    split_meta <- get_directional_transfer_split(stan_data, line_id, direction)

    fits <- lapply(c("transfer", "oracle"), function(ft) {
      load_best_transfer_fit(
        model_id = model_id,
        line_id = line_id,
        direction = direction,
        fit_type = ft,
        output_root = output_root,
        model_name = model_name
      )
    })
    names(fits) <- c("transfer", "oracle")

    state_map <- list(
      observed = split_meta$observed_value,
      holdout = split_meta$holdout_value
    )

    fit_state_vals <- list()
    for (fit_type in names(fits)) {
      draw_vec <- extract_draw_vector(fits[[fit_type]]$draws)
      fit_mask <- fits[[fit_type]]$split_meta$ploidy_effect_mask %||% NULL
      fit_state_vals[[fit_type]] <- list()

      for (state_name in names(state_map)) {
        vals <- reconstruct_line_state_parameters(
          draw_vec = draw_vec,
          model_id = model_id,
          line_id = line_id,
          ploidy_metric = state_map[[state_name]],
          ploidy_effect_mask = fit_mask
        )
        fit_state_vals[[fit_type]][[state_name]] <- vals

        state_rows[[idx_state]] <- data.frame(
          line_id = line_id,
          direction = direction,
          fit_type = fit_type,
          state = state_name,
          ploidy_metric = state_map[[state_name]],
          parameter = names(vals),
          value = as.numeric(vals),
          stringsAsFactors = FALSE
        )
        idx_state <- idx_state + 1L
      }

      common <- intersect(
        names(fit_state_vals[[fit_type]]$observed),
        names(fit_state_vals[[fit_type]]$holdout)
      )
      shift_vals <- fit_state_vals[[fit_type]]$holdout[common] - fit_state_vals[[fit_type]]$observed[common]

      shift_rows[[idx_shift]] <- data.frame(
        line_id = line_id,
        direction = direction,
        fit_type = fit_type,
        parameter = common,
        observed_value = as.numeric(fit_state_vals[[fit_type]]$observed[common]),
        holdout_value = as.numeric(fit_state_vals[[fit_type]]$holdout[common]),
        shift = as.numeric(shift_vals),
        stringsAsFactors = FALSE
      )
      idx_shift <- idx_shift + 1L
    }

    common_params <- Reduce(intersect, lapply(fit_state_vals, function(x) names(x$holdout)))
    missing_transfer <- fit_state_vals$transfer$holdout[common_params]
    missing_oracle <- fit_state_vals$oracle$holdout[common_params]
    observed_transfer <- fit_state_vals$transfer$observed[common_params]
    observed_oracle <- fit_state_vals$oracle$observed[common_params]
    null_holdout <- observed_transfer
    shift_transfer <- missing_transfer - observed_transfer
    shift_oracle <- missing_oracle - observed_oracle

    compare_rows[[idx_compare]] <- data.frame(
      line_id = line_id,
      direction = direction,
      parameter = common_params,
      null_holdout = as.numeric(null_holdout),
      transfer_holdout = as.numeric(missing_transfer),
      oracle_holdout = as.numeric(missing_oracle),
      null_vs_oracle_diff = as.numeric(null_holdout - missing_oracle),
      holdout_diff = as.numeric(missing_transfer - missing_oracle),
      transfer_improvement_over_null = as.numeric(abs(null_holdout - missing_oracle) - abs(missing_transfer - missing_oracle)),
      null_shift = 0,
      transfer_shift = as.numeric(shift_transfer),
      oracle_shift = as.numeric(shift_oracle),
      shift_diff = as.numeric(shift_transfer - shift_oracle),
      stringsAsFactors = FALSE
    )
    idx_compare <- idx_compare + 1L

    positive_mask <- (missing_transfer > 0) & (missing_oracle > 0)
    log_err <- rep(NA_real_, length(common_params))
    log_err[positive_mask] <- log(missing_transfer[positive_mask]) - log(missing_oracle[positive_mask])
    null_positive_mask <- (null_holdout > 0) & (missing_oracle > 0)
    null_log_err <- rep(NA_real_, length(common_params))
    null_log_err[null_positive_mask] <- log(null_holdout[null_positive_mask]) - log(missing_oracle[null_positive_mask])

    summary_rows[[idx_summary]] <- data.frame(
      line_id = line_id,
      direction = direction,
      n_parameters = length(common_params),
      mean_abs_null_vs_oracle_diff = mean(abs(null_holdout - missing_oracle), na.rm = TRUE),
      mean_abs_holdout_diff = mean(abs(missing_transfer - missing_oracle), na.rm = TRUE),
      mean_transfer_improvement_over_null = mean(abs(null_holdout - missing_oracle) - abs(missing_transfer - missing_oracle), na.rm = TRUE),
      prop_parameters_transfer_better_than_null = mean(abs(missing_transfer - missing_oracle) < abs(null_holdout - missing_oracle), na.rm = TRUE),
      rmse_null_vs_oracle_diff = sqrt(mean((null_holdout - missing_oracle)^2, na.rm = TRUE)),
      rmse_holdout_diff = sqrt(mean((missing_transfer - missing_oracle)^2, na.rm = TRUE)),
      mean_abs_shift_diff = mean(abs(shift_transfer - shift_oracle), na.rm = TRUE),
      rmse_shift_diff = sqrt(mean((shift_transfer - shift_oracle)^2, na.rm = TRUE)),
      mean_abs_log_null_holdout_ratio = mean(abs(null_log_err), na.rm = TRUE),
      mean_abs_log_holdout_ratio = mean(abs(log_err), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    idx_summary <- idx_summary + 1L
  }

  list(
    parameter_states = do.call(rbind, state_rows),
    parameter_shifts = do.call(rbind, shift_rows),
    parameter_comparison = do.call(rbind, compare_rows),
    comparison_summary = do.call(rbind, summary_rows)
  )
}
