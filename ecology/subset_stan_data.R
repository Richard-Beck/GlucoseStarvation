# ---- subset function ----
subset_stan_data <- function(stan_data, target_line, ploidy_level = "hi") {
  
  if (is.character(target_line)) {
    if (!target_line %in% names(stan_data$line_map)) stop("Cell line name not found in stan_data$line_map")
    line_idx <- stan_data$line_map[[target_line]]
  } else {
    line_idx <- target_line
  }
  
  wells_in_line <- which(stan_data$line_id == line_idx)
  if (!length(wells_in_line)) stop("No wells found for this line ID.")
  
  current_ploidies <- stan_data$ploidy_metric[wells_in_line]
  unique_ploidies <- sort(unique(current_ploidies))
  
  if (length(unique_ploidies) == 1) {
    target_ploidy <- unique_ploidies[1]
  } else {
    if (ploidy_level == "hi") target_ploidy <- max(unique_ploidies)
    else if (ploidy_level == "lo") target_ploidy <- min(unique_ploidies)
    else stop("ploidy_level must be 'hi' or 'lo'")
  }
  
  target_wells_orig_idx <- which(stan_data$line_id == line_idx & stan_data$ploidy_metric == target_ploidy)
  if (!length(target_wells_orig_idx)) stop("No wells match this Line + Ploidy combo.")
  
  well_map <- rep(NA_integer_, stan_data$N_wells)
  well_map[target_wells_orig_idx] <- seq_along(target_wells_orig_idx)
  
  new_data <- list()
  
  # Globals
  new_data$N_exps    <- stan_data$N_exps
  new_data$N_grid    <- stan_data$N_grid
  new_data$t_grid    <- stan_data$t_grid
  
  # Calibration
  new_data$calib_a_fixed     <- stan_data$calib_a_fixed
  new_data$calib_b_fixed     <- stan_data$calib_b_fixed
  new_data$calib_sigma_fixed <- stan_data$calib_sigma_fixed
  
  # Well-level
  new_data$N_wells     <- length(target_wells_orig_idx)
  new_data$G0_per_well <- stan_data$G0_per_well[target_wells_orig_idx]
  new_data$exp_id      <- stan_data$exp_id[target_wells_orig_idx]
  
  # Counts obs
  keep_count_idx <- which(stan_data$well_idx_count %in% target_wells_orig_idx)
  new_data$N_obs_count    <- length(keep_count_idx)
  new_data$well_idx_count <- well_map[stan_data$well_idx_count[keep_count_idx]]
  new_data$grid_idx_count <- stan_data$grid_idx_count[keep_count_idx]
  new_data$N_obs          <- stan_data$N_obs[keep_count_idx]
  new_data$D_obs          <- stan_data$D_obs[keep_count_idx]
  
  # Glucose/lum obs
  keep_gluc_idx <- which(stan_data$well_idx_gluc %in% target_wells_orig_idx)
  new_data$N_obs_gluc    <- length(keep_gluc_idx)
  new_data$well_idx_gluc <- well_map[stan_data$well_idx_gluc[keep_gluc_idx]]
  new_data$grid_idx_gluc <- stan_data$grid_idx_gluc[keep_gluc_idx]
  new_data$lum_obs       <- stan_data$lum_obs[keep_gluc_idx]
  new_data$dilution      <- stan_data$dilution[keep_gluc_idx]
  new_data$is_censored   <- stan_data$is_censored[keep_gluc_idx]
  
  new_data
}