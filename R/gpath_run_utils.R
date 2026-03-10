add_group_structure <- function(stan_data) {
  grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
  unique_grps <- unique(grp_keys)

  stan_data$N_groups <- length(unique_grps)
  stan_data$group_id <- match(grp_keys, unique_grps)
  stan_data
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

prepare_gpath_stan_data <- function(
  stan_data_path = "data/stan_ready_data.Rds",
  R_val,
  P_val,
  W_val,
  constraint_flag,
  waste_mech_flag,
  base_priors,
  holdout_high_ploidy = FALSE,
  cv_fold_id = NULL,
  drop_character = TRUE
) {
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
    drop_character = drop_character
  )

  stan_data$is_train <- rep(1L, stan_data$N_wells)

  if (isTRUE(holdout_high_ploidy)) {
    stan_data <- set_high_ploidy_holdout(stan_data)
  }

  if (!is.null(cv_fold_id)) {
    stan_data <- set_cv_fold_holdout(stan_data, cv_fold_id)
  }

  stan_data
}
