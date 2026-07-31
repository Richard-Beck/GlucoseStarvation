gpd_atomic_save_rds <- function(x, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = compress)
  if (!file.rename(tmp, path)) stop(sprintf("Could not atomically install %s", path), call. = FALSE)
  invisible(path)
}

gpd_expected_chain_paths <- function(nuts_dir, n_chains) {
  file.path(nuts_dir, sprintf("nuts_draws_chain%02d.Rds", seq_len(as.integer(n_chains))))
}

gpd_expected_diagnostic_paths <- function(nuts_dir, n_chains) {
  file.path(nuts_dir, sprintf("nuts_diagnostics_chain%02d.Rds", seq_len(as.integer(n_chains))))
}

gpd_require_files <- function(paths, label = "files") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(sprintf("Missing %s (%d): %s", label, length(missing), paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(paths)
}

gpd_chain_balanced_indices <- function(n_per_chain, max_draws = Inf) {
  n_per_chain <- as.integer(n_per_chain)
  if (!length(n_per_chain) || any(n_per_chain < 1L)) stop("Each chain must contain at least one draw", call. = FALSE)
  total <- sum(n_per_chain)
  if (!is.finite(max_draws) || max_draws >= total) return(lapply(n_per_chain, seq_len))
  max_draws <- as.integer(max_draws)
  if (max_draws < length(n_per_chain)) stop("max_draws must be at least the number of chains", call. = FALSE)
  allocation <- rep(max_draws %/% length(n_per_chain), length(n_per_chain))
  allocation[seq_len(max_draws %% length(n_per_chain))] <- allocation[seq_len(max_draws %% length(n_per_chain))] + 1L
  allocation <- pmin(allocation, n_per_chain)
  while (sum(allocation) < max_draws && any(allocation < n_per_chain)) {
    available <- which(allocation < n_per_chain)
    take <- available[seq_len(min(length(available), max_draws - sum(allocation)))]
    allocation[take] <- allocation[take] + 1L
  }
  Map(function(n, take) unique(as.integer(round(seq(1, n, length.out = take)))), n_per_chain, allocation)
}

gpd_read_fit_draws <- function(nuts_dir, n_chains, max_draws = Inf) {
  paths <- gpd_expected_chain_paths(nuts_dir, n_chains)
  gpd_require_files(paths, "NUTS chain files")
  chains <- lapply(paths, function(path) posterior::as_draws_matrix(readRDS(path)))
  variables <- colnames(chains[[1]])
  if (any(vapply(chains, function(x) !identical(colnames(x), variables), logical(1)))) {
    stop(sprintf("NUTS variable mismatch across chains in %s", nuts_dir), call. = FALSE)
  }
  selected <- gpd_chain_balanced_indices(vapply(chains, nrow, integer(1)), max_draws = max_draws)
  matrices <- Map(function(x, idx) x[idx, , drop = FALSE], chains, selected)
  index <- do.call(rbind, Map(function(chain_id, idx) {
    data.frame(chain_id = as.integer(chain_id), iteration = as.integer(idx), stringsAsFactors = FALSE)
  }, seq_along(selected), selected))
  index$draw_id <- seq_len(nrow(index))
  list(
    draws = do.call(rbind, matrices),
    index = index[, c("draw_id", "chain_id", "iteration")],
    chain_paths = paths
  )
}

gpd_bind_fit_draws <- function(nuts_dir, n_chains) {
  paths <- gpd_expected_chain_paths(nuts_dir, n_chains)
  gpd_require_files(paths, "NUTS chain files")
  draws <- lapply(paths, readRDS)
  do.call(posterior::bind_draws, c(draws, along = "chain"))
}

gpd_model_k <- function(model_id, stan_data) {
  get_hierarchical_k(model_id, N_lines = stan_data$N_lines)
}

gpd_pareto_membership <- function(k, deviance) {
  valid <- is.finite(k) & is.finite(deviance)
  out <- rep(FALSE, length(k))
  for (i in which(valid)) {
    dominated <- valid & k <= k[[i]] & deviance <= deviance[[i]] & (k < k[[i]] | deviance < deviance[[i]])
    out[[i]] <- !any(dominated)
  }
  out
}

gpd_line_table <- function(stan_data) {
  ids <- sort(unique(as.integer(stan_data$line_id)))
  names_by_id <- setNames(rep(NA_character_, length(ids)), ids)
  if (!is.null(stan_data$line_map) && length(stan_data$line_map)) {
    mapped <- as.integer(unname(stan_data$line_map))
    names(mapped) <- names(stan_data$line_map)
    for (id in ids) {
      hit <- names(mapped)[mapped == id]
      if (length(hit)) names_by_id[[as.character(id)]] <- hit[[1]]
    }
  }
  data.frame(line_id = ids, line_name = unname(names_by_id[as.character(ids)]), stringsAsFactors = FALSE)
}

gpd_reconstruct_parameter_table <- function(draw_info, stan_data, dataset_id, model_id) {
  lines <- gpd_line_table(stan_data)
  rows <- vector("list", nrow(lines) * 2L * nrow(draw_info$draws))
  ptr <- 1L
  for (line_i in seq_len(nrow(lines))) {
    states <- get_line_ploidy_states(stan_data, lines$line_id[[line_i]])
    state_values <- c(low = states$low_value, high = states$high_value)
    for (state_name in names(state_values)) {
      ploidy_metric <- state_values[[state_name]]
      for (draw_i in seq_len(nrow(draw_info$draws))) {
        draw_vec <- draw_info$draws[draw_i, , drop = TRUE]
        values <- reconstruct_line_state_parameters(
          draw_vec = draw_vec,
          model_id = model_id,
          line_id = lines$line_id[[line_i]],
          ploidy_metric = ploidy_metric
        )
        rows[[ptr]] <- data.frame(
          dataset_id = dataset_id,
          model_id = model_id,
          draw_id = draw_info$index$draw_id[[draw_i]],
          chain_id = draw_info$index$chain_id[[draw_i]],
          iteration = draw_info$index$iteration[[draw_i]],
          line_id = lines$line_id[[line_i]],
          line_name = lines$line_name[[line_i]],
          ploidy = state_name,
          ploidy_metric = ploidy_metric,
          parameter = names(values),
          value = as.numeric(values),
          stringsAsFactors = FALSE
        )
        ptr <- ptr + 1L
      }
    }
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

gpd_well_metadata <- function(stan_data) {
  line_names <- gpd_line_table(stan_data)
  out <- data.frame(
    well_idx = seq_len(stan_data$N_wells),
    line_id = as.integer(stan_data$line_id),
    ploidy_metric = as.numeric(stan_data$ploidy_metric),
    ploidy_abs = as.numeric(stan_data$ploidy_abs),
    exp_id = as.integer(stan_data$exp_id),
    G0 = as.numeric(stan_data$G0_per_well),
    group_id = as.integer(stan_data$group_id),
    g1_id = as.integer(stan_data$g1_id),
    stringsAsFactors = FALSE
  )
  merge(out, line_names, by = "line_id", all.x = TRUE, sort = FALSE)
}

gpd_simulate_draw <- function(draw_vec, stan_data, model_id, times) {
  dims <- parse_run_id(model_id)
  state_names <- c("N", paste0("R", seq_len(dims$R)), if (dims$W > 0L) paste0("W", seq_len(dims$W)))
  states <- array(NA_real_, dim = c(stan_data$N_wells, length(times), length(state_names)),
                  dimnames = list(well = seq_len(stan_data$N_wells), time = as.character(times), state = state_names))
  prior_N0_center <- as.numeric(stan_data$prior_N0_center)
  if (length(prior_N0_center) != 1L || !is.finite(prior_N0_center) || prior_N0_center <= 0) {
    stop("stan_data$prior_N0_center must be positive and finite", call. = FALSE)
  }
  for (well_idx in seq_len(stan_data$N_wells)) {
    parms <- gpath_reconstruct_state_from_draw(
      draw_vec = draw_vec,
      model_id = model_id,
      line_id = stan_data$line_id[[well_idx]],
      ploidy_metric = stan_data$ploidy_metric[[well_idx]]
    )
    n0_name <- sprintf("raw_N0[%d]", stan_data$group_id[[well_idx]])
    g1_name <- sprintf("G1_0[%d]", stan_data$g1_id[[well_idx]])
    if (!all(c(n0_name, g1_name) %in% names(draw_vec))) stop("Posterior draw lacks initial-condition parameters", call. = FALSE)
    y0 <- c(N = exp(log(prior_N0_center) + as.numeric(draw_vec[[n0_name]])),
            R1 = as.numeric(draw_vec[[g1_name]]))
    if (dims$R > 1L) y0 <- c(y0, setNames(rep(1, dims$R - 1L), paste0("R", 2:dims$R)))
    if (dims$W > 0L) y0 <- c(y0, setNames(rep(0, dims$W), paste0("W", seq_len(dims$W))))
    sim <- as.matrix(deSolve::ode(y = y0, times = times, func = rhs, parms = parms, method = "lsoda"))
    states[well_idx, , ] <- sim[, state_names, drop = FALSE]
  }
  states
}

gpd_growth_delta_surface <- function(draw_vec, stan_data, model_id, glucose_grid, resource2_grid, extra_resource_value = 1) {
  lines <- gpd_line_table(stan_data)
  grid <- expand.grid(glucose = glucose_grid, resource2 = resource2_grid, KEEP.OUT.ATTRS = FALSE)
  out <- array(NA_real_, dim = c(nrow(lines), length(glucose_grid), length(resource2_grid)),
               dimnames = list(line = lines$line_name, glucose = as.character(glucose_grid), resource2 = as.character(resource2_grid)))
  for (line_i in seq_len(nrow(lines))) {
    states <- get_line_ploidy_states(stan_data, lines$line_id[[line_i]])
    low <- gpath_reconstruct_state_from_draw(draw_vec, model_id, lines$line_id[[line_i]], states$low_value)
    high <- gpath_reconstruct_state_from_draw(draw_vec, model_id, lines$line_id[[line_i]], states$high_value)
    low_growth <- gpath_eval_instantaneous_net_growth_grid(low, grid, extra_resource_value = extra_resource_value)
    high_growth <- gpath_eval_instantaneous_net_growth_grid(high, grid, extra_resource_value = extra_resource_value)
    out[line_i, , ] <- matrix(high_growth - low_growth, nrow = length(glucose_grid), ncol = length(resource2_grid))
  }
  out
}

gpd_bind_first_dimension <- function(x) {
  if (!length(x)) return(NULL)
  d <- dim(x[[1]])
  out <- array(unlist(x, use.names = FALSE), dim = c(d, length(x)))
  out <- aperm(out, c(length(d) + 1L, seq_along(d)))
  dimnames(out) <- c(list(NULL), dimnames(x[[1]]))
  names(dimnames(out)) <- c("draw", names(dimnames(x[[1]])))
  out
}
