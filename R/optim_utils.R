# Functions shared across workflows for setting up or summarizing parameter optimization.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

load_optim_outputs_from_dir <- function(path) {
  list(
    lp = readRDS(file.path(path, "optim_lp_all.Rds")),
    draws = readRDS(file.path(path, "optim_draws_all.Rds"))
  )
}

summarize_lp_neighborhood <- function(lp, near_cut = 5) {
  finite_mask <- is.finite(lp)
  n_total <- length(lp)
  n_finite <- sum(finite_mask)

  if (!n_finite) {
    return(list(
      best_lp = NA_real_,
      lp_diff = rep(NA_real_, n_total),
      near_idx = integer(0),
      finite_mask = finite_mask
    ))
  }

  best_lp <- max(lp[finite_mask])
  lp_diff <- rep(NA_real_, n_total)
  lp_diff[finite_mask] <- best_lp - lp[finite_mask]
  near_idx <- which(finite_mask & lp_diff < near_cut)

  list(
    best_lp = best_lp,
    lp_diff = lp_diff,
    near_idx = near_idx,
    finite_mask = finite_mask
  )
}

extract_draw_vector <- function(draws) {
  if (is.matrix(draws)) {
    return(draws[1, , drop = TRUE])
  }

  if (is.data.frame(draws)) {
    return(as.matrix(draws)[1, , drop = TRUE])
  }

  if (is.numeric(draws) && !is.null(names(draws))) {
    return(draws)
  }

  stop("Unsupported draws format for parameter extraction")
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

read_optim_spec <- function(spec_path) {
  # NOTE: This mirrors spec-definition logic that currently lives in
  # `scripts/make_gpath_optim_spec.R` and `slurm/submit_gpath_optim_pipeline.sh`.
  # That duplication should ideally be removed in future by centralizing spec
  # schema and path resolution in one place.
  spec_df <- utils::read.delim(
    spec_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "enabled", "model_name", "model_version", "run_id", "n_starts",
    "num_threads", "stan_data_path", "dataset_label", "run_label",
    "array_cpus", "array_mem_gb", "array_time", "combine_mem_gb",
    "combine_time", "qos", "delete_after"
  )
  missing_cols <- setdiff(required_cols, names(spec_df))
  if (length(missing_cols)) {
    stop(sprintf(
      "Spec is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  spec_df$spec_row_id <- seq_len(nrow(spec_df))
  batch_label <- basename(normalizePath(dirname(spec_path), winslash = "/", mustWork = FALSE))
  run_label_resolved <- trimws(spec_df$run_label)
  run_label_resolved[!nzchar(run_label_resolved)] <- batch_label
  spec_df$optim_output_dir <- file.path(
    "data",
    "runs",
    spec_df$model_name,
    spec_df$model_version,
    "optim",
    spec_df$dataset_label,
    spec_df$run_id,
    run_label_resolved
  )
  spec_df$optim_lp_path <- file.path(spec_df$optim_output_dir, "optim_lp_all.Rds")
  spec_df$optim_draws_path <- file.path(spec_df$optim_output_dir, "optim_draws_all.Rds")
  spec_df$optim_rc_path <- file.path(spec_df$optim_output_dir, "optim_rc_all.Rds")
  spec_df$optim_lp_exists <- file.exists(spec_df$optim_lp_path)
  spec_df$optim_draws_exists <- file.exists(spec_df$optim_draws_path)
  spec_df$optim_complete <- spec_df$optim_lp_exists & spec_df$optim_draws_exists
  spec_df
}

extract_draw_metrics <- function(draw_vec) {
  get_metric <- function(name) {
    if (is.null(draw_vec) || is.null(names(draw_vec)) || !(name %in% names(draw_vec))) {
      return(NA_real_)
    }
    as.numeric(draw_vec[[name]])
  }

  list(
    log_lik_train = get_metric("log_lik_train"),
    log_lik_holdout = get_metric("log_lik_holdout")
  )
}

compute_optim_k <- function(spec_row) {
  if (is.data.frame(spec_row)) {
    if (nrow(spec_row) != 1L) {
      stop("spec_row data.frame must have exactly one row")
    }
    spec_row <- as.list(spec_row[1, , drop = FALSE])
  }

  source("R/project_paths.R")
  source(get_model_r_path(
    as.character(spec_row$model_name),
    as.character(spec_row$model_version)
  ))
  source("R/gpath_run_utils.R")

  stan_data <- readRDS(resolve_stan_data_path(as.character(spec_row$stan_data_path)))
  get_hierarchical_k(as.character(spec_row$run_id), N_lines = stan_data$N_lines)
}

summarize_optim_row <- function(spec_row, include_k = TRUE) {
  if (is.data.frame(spec_row)) {
    if (nrow(spec_row) != 1L) {
      stop("spec_row data.frame must have exactly one row")
    }
    spec_row <- as.list(spec_row[1, , drop = FALSE])
  }

  output_dir <- as.character(spec_row$optim_output_dir %||% "")
  if (!nzchar(output_dir)) {
    stop("spec_row must contain optim_output_dir")
  }

  opt <- load_optim_outputs_from_dir(output_dir)
  nbhd <- summarize_lp_neighborhood(opt$lp, near_cut = 5)
  best <- extract_best_draw_from_optim_outputs(
    draws_all = opt$draws,
    lp_all = opt$lp
  )
  draw_metrics <- extract_draw_metrics(best$draw_vec)

  out <- data.frame(
    spec_row_id = as.integer(spec_row$spec_row_id %||% NA_integer_),
    n_starts_total = length(opt$lp),
    n_starts_valid = sum(nbhd$finite_mask),
    best_draw_index = as.integer(best$best_idx),
    best_lp = as.numeric(best$best_lp),
    n_within_1 = sum(nbhd$finite_mask & nbhd$lp_diff < 1),
    n_within_5 = sum(nbhd$finite_mask & nbhd$lp_diff < 5),
    log_lik_train = draw_metrics$log_lik_train,
    log_lik_holdout = draw_metrics$log_lik_holdout,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_k)) {
    out$k <- compute_optim_k(spec_row)
  }

  out
}

extract_best_draw_from_optim_row <- function(spec_row) {
  if (is.data.frame(spec_row)) {
    if (nrow(spec_row) != 1L) {
      stop("spec_row data.frame must have exactly one row")
    }
    spec_row <- as.list(spec_row[1, , drop = FALSE])
  }

  output_dir <- as.character(spec_row$optim_output_dir %||% "")
  if (!nzchar(output_dir)) {
    stop("spec_row must contain optim_output_dir")
  }

  opt <- load_optim_outputs_from_dir(output_dir)
  extract_best_draw_from_optim_outputs(
    draws_all = opt$draws,
    lp_all = opt$lp
  )$draw_vec
}

build_optim_objects_from_spec <- function(
  spec_path,
  include_disabled = FALSE,
  include_k = TRUE,
  keep_incomplete_summary = FALSE
) {
  spec_df <- read_optim_spec(spec_path)
  if (!include_disabled) {
    spec_df <- spec_df[as.integer(spec_df$enabled) != 0L, , drop = FALSE]
  }

  if (!nrow(spec_df)) {
    return(list(
      summary_df = spec_df,
      best_draws = list()
    ))
  }

  complete_idx <- which(spec_df$optim_complete)
  incomplete_idx <- which(!spec_df$optim_complete)

  if (length(incomplete_idx)) {
    warning(sprintf(
      "Skipping %d incomplete optimization row(s) with missing optim_lp_all.Rds and/or optim_draws_all.Rds.",
      length(incomplete_idx)
    ))
  }

  if (!length(complete_idx)) {
    return(list(
      summary_df = spec_df[0, , drop = FALSE],
      best_draws = list()
    ))
  }

  spec_df_complete <- spec_df[complete_idx, , drop = FALSE]

  summary_rows <- lapply(seq_len(nrow(spec_df_complete)), function(i) {
    summarize_optim_row(spec_row = spec_df_complete[i, , drop = FALSE], include_k = include_k)
  })
  best_draws <- lapply(seq_len(nrow(spec_df_complete)), function(i) {
    extract_best_draw_from_optim_row(spec_row = spec_df_complete[i, , drop = FALSE])
  })

  summary_df <- cbind(
    spec_df_complete,
    do.call(rbind, summary_rows)[, setdiff(names(summary_rows[[1]]), "spec_row_id"), drop = FALSE]
  )

  if (isTRUE(keep_incomplete_summary) && length(incomplete_idx)) {
    incomplete_df <- spec_df[incomplete_idx, , drop = FALSE]
    metric_cols <- setdiff(names(summary_df), names(spec_df_complete))
    for (nm in metric_cols) {
      incomplete_df[[nm]] <- NA
    }
    summary_df <- rbind(summary_df, incomplete_df[, names(summary_df), drop = FALSE])
    summary_df <- summary_df[order(summary_df$spec_row_id), , drop = FALSE]
  }

  list(
    summary_df = summary_df,
    best_draws = best_draws
  )
}
