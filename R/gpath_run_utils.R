resolve_stan_data_path <- function(
  stan_data_path = NULL,
  candidates = NULL
) {
  if (is.null(candidates)) {
    if (!exists("default_stan_data_candidates")) {
      source("R/project_paths.R")
    }
    candidates <- default_stan_data_candidates()
  }

  paths <- unique(c(stan_data_path, candidates))
  paths <- paths[!is.na(paths) & nzchar(paths)]

  for (path in paths) {
    if (file.exists(path)) {
      return(path)
    }
  }

  stop(sprintf(
    "Could not find stan_ready_data.Rds. Tried: %s",
    paste(paths, collapse = ", ")
  ))
}

add_group_structure <- function(stan_data) {
  grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
  unique_grps <- unique(grp_keys)

  stan_data$N_groups <- length(unique_grps)
  stan_data$group_id <- match(grp_keys, unique_grps)
  stan_data
}

subset_stan_data_to_line <- function(stan_data, target_line_id) {
  target_line_id <- as.integer(target_line_id)[1]
  if (is.na(target_line_id)) {
    stop("target_line_id must be a finite integer")
  }

  old_line_id <- as.integer(stan_data$line_id)
  keep_wells <- which(old_line_id == target_line_id)
  if (!length(keep_wells)) {
    stop(sprintf("No wells found for line_id == %d", target_line_id))
  }

  remap_index <- function(idx, keep) {
    out <- integer(length(idx))
    out[keep] <- seq_along(keep)
    out
  }

  subset_by_well <- function(x, keep) {
    if (is.null(x)) {
      return(NULL)
    }

    if (is.matrix(x)) {
      return(x[keep, , drop = FALSE])
    }

    if (is.array(x) && length(dim(x)) > 1L) {
      return(x[keep, , drop = FALSE])
    }

    x[keep]
  }

  well_map <- remap_index(seq_len(stan_data$N_wells), keep_wells)

  keep_count <- stan_data$well_idx_count %in% keep_wells
  keep_gluc <- stan_data$well_idx_gluc %in% keep_wells

  out <- stan_data

  well_fields <- c(
    "line_id", "ploidy_metric", "ploidy_abs", "has_starvation", "exp_id",
    "G0_per_well", "g1_id", "group_id", "R_init_base", "is_train"
  )
  for (nm in intersect(well_fields, names(out))) {
    out[[nm]] <- subset_by_well(out[[nm]], keep_wells)
  }

  count_fields <- c("well_idx_count", "grid_idx_count", "rep_id_count", "N_obs", "D_obs")
  for (nm in intersect(count_fields, names(out))) {
    out[[nm]] <- out[[nm]][keep_count]
  }
  out$well_idx_count <- well_map[out$well_idx_count]

  gluc_fields <- c("well_idx_gluc", "grid_idx_gluc", "lum_obs", "dilution", "is_censored")
  for (nm in intersect(gluc_fields, names(out))) {
    out[[nm]] <- out[[nm]][keep_gluc]
  }
  out$well_idx_gluc <- well_map[out$well_idx_gluc]

  exp_keep <- sort(unique(out$exp_id))
  exp_map <- remap_index(seq_len(stan_data$N_exps), exp_keep)
  out$exp_id <- exp_map[out$exp_id]
  out$N_exps <- length(exp_keep)

  calib_fields <- c("calib_a_fixed", "calib_b_fixed", "calib_sigma_fixed")
  for (nm in intersect(calib_fields, names(out))) {
    out[[nm]] <- out[[nm]][exp_keep]
  }

  if ("calib_exp_idx" %in% names(out)) {
    keep_calib <- out$calib_exp_idx %in% exp_keep
    out$calib_exp_idx <- exp_map[out$calib_exp_idx[keep_calib]]

    for (nm in c("calib_G", "calib_Lum")) {
      if (nm %in% names(out)) {
        out[[nm]] <- out[[nm]][keep_calib]
      }
    }

    if ("N_obs_calib" %in% names(out)) {
      out$N_obs_calib <- length(out$calib_exp_idx)
    }
  }

  line_name <- NULL
  if ("line_map" %in% names(out) && length(out$line_map)) {
    line_matches <- names(out$line_map)[as.integer(out$line_map) == target_line_id]
    if (length(line_matches)) {
      line_name <- line_matches[1]
      out$line_map <- setNames(1L, line_name)
    } else {
      out$line_map <- setNames(1L, as.character(target_line_id))
    }
  }

  out$line_id <- rep.int(1L, length(keep_wells))
  out$N_wells <- length(keep_wells)
  out$N_lines <- 1L
  out$N_obs_count <- length(out$well_idx_count)
  out$N_obs_gluc <- length(out$well_idx_gluc)

  out <- add_group_structure(out)

  g1_key <- interaction(out$line_id, out$exp_id, out$G0_per_well, drop = TRUE)
  out$g1_id <- as.integer(g1_key)
  out$N_G1 <- max(out$g1_id)
  out$g1_ref_well <- match(seq_len(out$N_G1), out$g1_id)

  out$subset_source_line_id <- target_line_id
  if (!is.null(line_name)) {
    out$subset_source_line_name <- line_name
  }

  out
}

build_R_init_base <- function(stan_data, R_val) {
  R_init_base <- matrix(0.0, nrow = stan_data$N_wells, ncol = R_val)
  if (R_val > 1) {
    for (c in 2:R_val) {
      R_init_base[, c] <- 1.0
    }
  }
  R_init_base
}

apply_gpath_run_config <- function(
  stan_data,
  R_val,
  P_val,
  W_val,
  constraint_flag,
  waste_mech_flag,
  base_priors,
  ploidy_effect_mask_spec = NULL,
  drop_character = TRUE
) {
  config <- generate_stan_config(
    R = R_val,
    P = P_val,
    W = W_val,
    strict_spec = (constraint_flag == 1L),
    M = waste_mech_flag,
    base_priors = base_priors
  )

  for (nm in names(config)) {
    stan_data[[nm]] <- config[[nm]]
  }

  mask_info <- build_ploidy_effect_mask(
    R = R_val,
    P = P_val,
    W = W_val,
    C = constraint_flag,
    M = waste_mech_flag,
    target_spec = ploidy_effect_mask_spec,
    priors = base_priors
  )
  stan_data$ploidy_effect_mask <- as.numeric(mask_info$mask)
  attr(stan_data, "ploidy_effect_mask_info") <- mask_info

  stan_data$waste_mech <- if (W_val > 0) {
    rep(as.numeric(waste_mech_flag), W_val)
  } else {
    numeric(0)
  }

  stan_data$R_init_base <- build_R_init_base(stan_data, R_val)

  if (drop_character) {
    stan_data <- stan_data[!sapply(stan_data, is.character)]
  }

  stan_data
}

set_high_ploidy_holdout <- function(stan_data) {
  is_train <- rep(1L, stan_data$N_wells)
  is_train[stan_data$ploidy_metric > 0] <- 0L
  stan_data$is_train <- is_train
  stan_data
}

set_cv_fold_holdout <- function(stan_data, fold_id, tol = 1e-12) {
  line_vec <- as.integer(stan_data$line_id)
  ploidy_vec <- as.numeric(stan_data$ploidy_metric)

  idx_line <- which(line_vec == fold_id)
  if (!length(idx_line)) {
    stop(sprintf("No wells found for line_id == %d", fold_id))
  }

  hi_ploidy <- max(ploidy_vec[idx_line], na.rm = TRUE)
  idx_holdout <- which(line_vec == fold_id & abs(ploidy_vec - hi_ploidy) < tol)
  if (!length(idx_holdout)) {
    stop(sprintf("No heldout wells found for line=%d", fold_id))
  }

  stan_data$is_train <- rep(1L, stan_data$N_wells)
  stan_data$is_train[idx_holdout] <- 0L
  stan_data
}

normalize_transfer_direction <- function(direction) {
  if (!length(direction) || is.na(direction)) {
    stop("direction must be provided")
  }

  key <- gsub("\\s+", "", tolower(direction))

  if (key %in% c("low_to_high", "low->high", "l2h", "lowhigh")) {
    return("low_to_high")
  }

  if (key %in% c("high_to_low", "high->low", "h2l", "highlow")) {
    return("high_to_low")
  }

  stop(sprintf("Unsupported transfer direction '%s'", direction))
}

get_line_ploidy_states <- function(stan_data, line_id, tol = 1e-12) {
  line_vec <- as.integer(stan_data$line_id)
  ploidy_vec <- as.numeric(stan_data$ploidy_metric)

  idx_line <- which(line_vec == line_id)
  if (!length(idx_line)) {
    stop(sprintf("No wells found for line_id == %d", line_id))
  }

  values <- sort(unique(ploidy_vec[idx_line]))
  values <- values[is.finite(values)]

  if (length(values) != 2L) {
    stop(sprintf(
      "Directional transfer requires exactly 2 ploidy states for line %d; found %d",
      line_id,
      length(values)
    ))
  }

  low_value <- values[1]
  high_value <- values[length(values)]

  list(
    line_id = as.integer(line_id),
    low_value = low_value,
    high_value = high_value,
    tol = tol
  )
}

get_directional_transfer_split <- function(stan_data, line_id, direction, tol = 1e-12) {
  direction <- normalize_transfer_direction(direction)
  states <- get_line_ploidy_states(stan_data, line_id = line_id, tol = tol)

  observed_value <- if (direction == "low_to_high") states$low_value else states$high_value
  holdout_value <- if (direction == "low_to_high") states$high_value else states$low_value

  line_vec <- as.integer(stan_data$line_id)
  ploidy_vec <- as.numeric(stan_data$ploidy_metric)

  observed_wells <- which(line_vec == line_id & abs(ploidy_vec - observed_value) < tol)
  holdout_wells <- which(line_vec == line_id & abs(ploidy_vec - holdout_value) < tol)

  if (!length(observed_wells) || !length(holdout_wells)) {
    stop(sprintf(
      "Could not build transfer split for line %d (%s)",
      line_id,
      direction
    ))
  }

  count_obs <- sum(stan_data$well_idx_count %in% holdout_wells)
  gluc_obs <- sum(stan_data$well_idx_gluc %in% holdout_wells)
  holdout_times <- sort(unique(stan_data$t_grid[
    unique(c(
      stan_data$grid_idx_count[stan_data$well_idx_count %in% holdout_wells],
      stan_data$grid_idx_gluc[stan_data$well_idx_gluc %in% holdout_wells]
    ))
  ]))

  list(
    line_id = as.integer(line_id),
    direction = direction,
    observed_value = observed_value,
    holdout_value = holdout_value,
    observed_wells = observed_wells,
    holdout_wells = holdout_wells,
    holdout_count_obs = as.integer(count_obs),
    holdout_gluc_obs = as.integer(gluc_obs),
    holdout_total_obs = as.integer(count_obs + gluc_obs),
    holdout_times = holdout_times
  )
}

apply_directional_transfer_holdout <- function(stan_data, line_id, direction, tol = 1e-12) {
  split <- get_directional_transfer_split(
    stan_data = stan_data,
    line_id = line_id,
    direction = direction,
    tol = tol
  )

  stan_data$is_train <- rep(1L, stan_data$N_wells)
  stan_data$is_train[split$holdout_wells] <- 0L
  attr(stan_data, "transfer_split") <- split
  stan_data
}

count_holdout_observations <- function(stan_data, holdout_wells) {
  count_obs <- sum(stan_data$well_idx_count %in% holdout_wells)
  gluc_obs <- sum(stan_data$well_idx_gluc %in% holdout_wells)

  list(
    count = as.integer(count_obs),
    glucose = as.integer(gluc_obs),
    total = as.integer(count_obs + gluc_obs)
  )
}

prepare_gpath_stan_data <- function(
  stan_data_path = NULL,
  R_val,
  P_val,
  W_val,
  constraint_flag,
  waste_mech_flag,
  base_priors,
  ploidy_effect_mask_spec = NULL,
  holdout_high_ploidy = FALSE,
  cv_fold_id = NULL,
  transfer_line_id = NULL,
  transfer_direction = NULL,
  drop_character = TRUE
) {
  stan_data_path <- resolve_stan_data_path(stan_data_path)
  stan_data <- readRDS(stan_data_path)
  stan_data <- add_group_structure(stan_data)
  stan_data <- apply_gpath_run_config(
    stan_data = stan_data,
    R_val = R_val,
    P_val = P_val,
    W_val = W_val,
    constraint_flag = constraint_flag,
    waste_mech_flag = waste_mech_flag,
    base_priors = base_priors,
    ploidy_effect_mask_spec = ploidy_effect_mask_spec,
    drop_character = drop_character
  )

  stan_data$is_train <- rep(1L, stan_data$N_wells)

  if (isTRUE(holdout_high_ploidy)) {
    stan_data <- set_high_ploidy_holdout(stan_data)
  }

  if (!is.null(cv_fold_id)) {
    stan_data <- set_cv_fold_holdout(stan_data, cv_fold_id)
  }

  if (!is.null(transfer_line_id)) {
    stan_data <- apply_directional_transfer_holdout(
      stan_data = stan_data,
      line_id = transfer_line_id,
      direction = transfer_direction
    )
  }

  stan_data
}
