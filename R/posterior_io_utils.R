list_matching_rds <- function(path, pattern) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (!length(files)) {
    stop(sprintf("No files found in '%s' for pattern '%s'", path, pattern))
  }
  sort(files)
}

bind_rds_draws <- function(files, along = "chain") {
  draws_list <- lapply(files, readRDS)
  do.call(posterior::bind_draws, c(draws_list, along = along))
}

load_draws_matrix_from_dir <- function(path, pattern = "^nuts_draws_[0-9]+\\.Rds$") {
  files <- list_matching_rds(path, pattern)
  posterior::as_draws_matrix(bind_rds_draws(files, along = "chain"))
}

load_gpath_nuts_matrix <- function(model_id, base_path = "ecology/model_selection/data/gpath") {
  nuts_dir <- file.path(base_path, model_id, "hier", "nuts")
  load_draws_matrix_from_dir(nuts_dir, pattern = "^nuts_draws_[0-9]+\\.Rds$")
}

load_optim_outputs <- function(model_id, base_path = "ecology/model_selection/data/gpath") {
  fit_dir <- file.path(base_path, model_id, "hier")
  list(
    lp = readRDS(file.path(fit_dir, "optim_lp_all.Rds")),
    draws = readRDS(file.path(fit_dir, "optim_draws_all.Rds"))
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

