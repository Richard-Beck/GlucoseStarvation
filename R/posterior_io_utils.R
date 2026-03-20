list_matching_rds <- function(path, pattern) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (!length(files)) {
    stop(sprintf("No files found in '%s' for pattern '%s'", path, pattern))
  }
  sort(files)
}

log_mean_exp <- function(x) {
  if (!length(x)) {
    return(NA_real_)
  }

  xmax <- max(x)
  xmax + log(mean(exp(x - xmax)))
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

clamp_constrained_init <- function(init, maxG0) {
  if (!is.null(init$raw_sigma_base)) {
    init$raw_sigma_base <- as.numeric(init$raw_sigma_base)
  }

  if (!is.null(init$raw_sigma_rel)) {
    init$raw_sigma_rel <- as.numeric(init$raw_sigma_rel)
  }

  if (!is.null(init$nu_N)) {
    init$nu_N <- max(as.numeric(init$nu_N), 2 + 1e-12)
  }

  if (!is.null(init$G1_0)) {
    init$G1_0 <- pmax(init$G1_0, 1e-12)
    init$G1_0 <- pmin(init$G1_0, (2.0 * maxG0) - 1e-12)
  }

  init
}

flat_to_init_list <- function(
  x_named,
  param_names = c(
    "raw_theta_line",
    "raw_theta_ploidy",
    "raw_N0",
    "raw_sigma_N",
    "raw_sigma_base",
    "raw_sigma_rel",
    "nu_N",
    "G1_0"
  ),
  maxG0
) {
  keep <- !grepl("__$", names(x_named)) & names(x_named) != "lp__"
  x_named <- x_named[keep]
  x_named <- x_named[is.finite(x_named)]

  nms <- names(x_named)
  if (is.null(nms) || !length(nms)) {
    stop("init vector has no usable names")
  }

  out <- list()

  scalars <- !grepl("\\[", nms)
  if (any(scalars)) {
    for (nm in nms[scalars]) {
      if (nm %in% param_names) {
        out[[nm]] <- unname(x_named[[nm]])
      }
    }
  }

  idx_nms <- nms[!scalars]
  if (length(idx_nms)) {
    base <- sub("\\[.*$", "", idx_nms)

    for (bn in unique(base)) {
      if (!(bn %in% param_names)) {
        next
      }

      these <- idx_nms[base == bn]
      idx_list <- lapply(these, function(s) {
        as.integer(strsplit(gsub("^.*\\[|\\]$", "", s), ",")[[1]])
      })

      dims <- pmax(1L, Reduce(pmax, idx_list))
      arr <- array(0.0, dim = dims)

      for (k in seq_along(these)) {
        ii <- idx_list[[k]]
        arr[matrix(ii, nrow = 1)] <- unname(x_named[[these[[k]]]])
      }

      if (length(dims) == 1L) {
        out[[bn]] <- as.numeric(arr)
      } else if (length(dims) == 2L) {
        out[[bn]] <- matrix(arr, nrow = dims[1], ncol = dims[2])
      } else {
        out[[bn]] <- arr
      }
    }
  }

  clamp_constrained_init(out, maxG0 = maxG0)
}

sample_posterior_init <- function(draws, seed, maxG0, param_names = NULL) {
  dm <- posterior::as_draws_matrix(draws)
  if (nrow(dm) < 1) {
    stop("Posterior draws matrix is empty")
  }

  set.seed(seed)
  pick <- sample.int(nrow(dm), 1)
  x_named <- dm[pick, , drop = TRUE]
  x_named <- x_named[is.finite(x_named)]
  x_named <- x_named[!is.na(x_named)]

  if (is.null(param_names)) {
    return(flat_to_init_list(x_named = x_named, maxG0 = maxG0))
  }

  flat_to_init_list(
    x_named = x_named,
    param_names = param_names,
    maxG0 = maxG0
  )
}

load_posterior_init_from_dir <- function(
  path,
  seed,
  maxG0,
  pattern = "^nuts_draws_[0-9]+\\.Rds$",
  param_names = NULL
) {
  draws <- bind_rds_draws(list_matching_rds(path, pattern), along = "chain")
  sample_posterior_init(
    draws = draws,
    seed = seed,
    maxG0 = maxG0,
    param_names = param_names
  )
}

load_optim_init_from_dir <- function(
  path,
  chain_id = 1L,
  maxG0,
  param_names = NULL
) {
  opt <- load_optim_outputs_from_dir(path)
  lp <- opt$lp
  draws_list <- opt$draws

  valid_idx <- which(is.finite(lp) & !vapply(draws_list, is.null, logical(1)))
  if (!length(valid_idx)) {
    stop(sprintf("No valid optimization starts found in '%s'", path))
  }

  sorted_idx <- valid_idx[order(lp[valid_idx], decreasing = TRUE)]
  pick <- sorted_idx[(as.integer(chain_id) - 1L) %% length(sorted_idx) + 1L]
  draw_obj <- draws_list[[pick]]

  if (is.matrix(draw_obj)) {
    x_named <- draw_obj[1, , drop = TRUE]
    if (!is.null(colnames(draw_obj))) {
      names(x_named) <- colnames(draw_obj)
    }
  } else {
    x_named <- draw_obj
  }

  if (is.null(names(x_named)) || !length(names(x_named))) {
    stop(sprintf("Chosen optimization start %d in '%s' has no parameter names", pick, path))
  }
  init <- flat_to_init_list(
    x_named = x_named,
    param_names = param_names,
    maxG0 = maxG0
  )

  if (!length(init)) {
    stop(sprintf(
      "Chosen optimization start %d in '%s' produced an empty init list; parameter names may be incompatible with the current Stan model",
      pick,
      path
    ))
  }

  init
}

get_well_loglik_draws <- function(draws, well_idx) {
  dm <- posterior::as_draws_matrix(draws)
  cols <- sprintf("ll_well[%d]", as.integer(well_idx))
  available <- colnames(dm)
  missing <- setdiff(cols, available)

  if (length(missing)) {
    alt_cols <- make.names(cols)
    if (all(alt_cols %in% available)) {
      cols <- alt_cols
      missing <- character(0)
    }
  }

  if (length(missing)) {
    stop(sprintf(
      "Missing ll_well columns: %s",
      paste(utils::head(missing, 5), collapse = ", ")
    ))
  }

  rowSums(dm[, cols, drop = FALSE])
}

compute_well_subset_elpd <- function(draws, well_idx) {
  ll_draws <- get_well_loglik_draws(draws, well_idx)
  list(
    elpd = log_mean_exp(ll_draws),
    ll_draws = ll_draws
  )
}

