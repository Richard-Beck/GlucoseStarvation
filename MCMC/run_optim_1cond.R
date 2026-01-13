library(cmdstanr)
library(dplyr)
library(posterior)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME <- if (length(args) >= 1) args[1] else "model_B1cond"
cellLine <- if (length(args) >= 2) args[2] else "MDA-MB-231"
ploidy <- if (length(args) >= 3) args[3] else "lo"
NUM_PATHS <- if (length(args) >= 4) as.integer(args[4]) else 4
N_CORES <- if (length(args) >= 5) as.integer(args[5]) else 4
# Usually valid to simulate data during pathfinder to check fit visually
CALC_SIM   <- 1 
NUM_THREADS <- 1
cat(sprintf(">>> RUNNING PATHFINDER (Approx) FOR: %s\n", MODEL_NAME))


# 1. LOAD & PREP (Identical logic to run_fit.R)
stan_data <- readRDS("data/stan_ready_data.Rds")
config    <- read_json(paste0("MCMC/",MODEL_NAME, ".json"), simplifyVector = TRUE)

if (MODEL_NAME == "model_B1cond") {
  max_N <- max(stan_data$N_obs, na.rm = TRUE)
  config$prior_means[1] <- log(max_N * 1.5)
  config$prior_sds[1]   <- 0.5
  stan_data$lower_b <- config$lower
  stan_data$upper_b <- config$upper
  stan_data$N_params <- length(config$lower)
  stan_data$has_starvation <- stan_data$has_starvation*0
}

stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0 
stan_data$calc_sim       <- CALC_SIM
stan_data$structure_mode <- 0
# 2. COMPILE
stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
#mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE),force_recompile=T)
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

subset_stan_data <- function(stan_data, target_line, ploidy_level = "hi") {
  # target_line: Can be integer ID (e.g., 3) or string name (e.g., "SNU668")
  # ploidy_level: "hi" or "lo"
  
  # 1. Resolve Line ID
  if (is.character(target_line)) {
    if (!target_line %in% names(stan_data$line_map)) {
      stop("Cell line name not found in stan_data$line_map")
    }
    line_idx <- stan_data$line_map[[target_line]]
  } else {
    line_idx <- target_line
  }
  
  # 2. Identify wells belonging to this line
  # We look at the original indices of wells that match the line_id
  wells_in_line <- which(stan_data$line_id == line_idx)
  
  if (length(wells_in_line) == 0) stop("No wells found for this line ID.")
  
  # 3. Resolve Ploidy Logic
  # Get the ploidy metrics for just these wells
  current_ploidies <- stan_data$ploidy_metric[wells_in_line]
  unique_ploidies <- sort(unique(current_ploidies))
  
  if (length(unique_ploidies) < 1) {
    stop("No ploidy data found for this line.")
  } else if (length(unique_ploidies) == 1) {
    message(paste("Only one ploidy level found:", unique_ploidies[1], "- Using it."))
    target_ploidy <- unique_ploidies[1]
  } else {
    # If there are multiple, pick based on hi/lo arg
    if (ploidy_level == "hi") {
      target_ploidy <- max(unique_ploidies)
    } else if (ploidy_level == "lo") {
      target_ploidy <- min(unique_ploidies)
    } else {
      stop("ploidy_level must be 'hi' or 'lo'")
    }
  }
  
  # 4. Final Well Selection
  # Intersect line wells with ploidy wells
  # We use the original indices to compare against the full vectors
  target_wells_orig_idx <- which(stan_data$line_id == line_idx & 
                                   stan_data$ploidy_metric == target_ploidy)
  
  if (length(target_wells_orig_idx) == 0) stop("No wells match this Line + Ploidy combo.")
  
  # 5. Create Mapping: Old Index -> New Index
  # This map helps us re-index the observation vectors
  # Initialize with NA; only target wells get a new index 1:N
  well_map <- rep(NA, stan_data$N_wells)
  well_map[target_wells_orig_idx] <- 1:length(target_wells_orig_idx)
  
  # 6. Construct New Data List
  new_data <- list()
  
  # --- A. Pass through global constants ---
  new_data$N_params  <- stan_data$N_params
  new_data$lower_b   <- stan_data$lower_b
  new_data$upper_b   <- stan_data$upper_b
  new_data$N_exps    <- stan_data$N_exps
  new_data$N_grid    <- stan_data$N_grid
  new_data$t_grid    <- stan_data$t_grid
  new_data$grainsize <- stan_data$grainsize
  new_data$calc_sim  <- stan_data$calc_sim # Keep simulation flag
  
  # --- B. Priors (Pass through) ---
  new_data$prior_mu_N0_mean <- stan_data$prior_mu_N0_mean
  new_data$prior_mu_N0_sd   <- stan_data$prior_mu_N0_sd
  new_data$prior_mu_D0_mean <- stan_data$prior_mu_D0_mean
  new_data$prior_mu_D0_sd   <- stan_data$prior_mu_D0_sd
  new_data$prior_ode_mean   <- stan_data$prior_ode_mean
  new_data$prior_ode_sd     <- stan_data$prior_ode_sd
  
  # --- C. Calibration Params ---
  new_data$calib_a_fixed     <- stan_data$calib_a_fixed
  new_data$calib_b_fixed     <- stan_data$calib_b_fixed
  new_data$calib_sigma_fixed <- stan_data$calib_sigma_fixed
  
  # --- D. Subset Well-Level Data ---
  new_data$N_wells        <- length(target_wells_orig_idx)
  new_data$G0_per_well    <- stan_data$G0_per_well[target_wells_orig_idx]
  new_data$has_starvation <- stan_data$has_starvation[target_wells_orig_idx]
  new_data$exp_id         <- stan_data$exp_id[target_wells_orig_idx]
  
  # --- E. Subset & Re-index Count Observations ---
  # Keep only obs where the well index is in our target set
  keep_count_idx <- which(stan_data$well_idx_count %in% target_wells_orig_idx)
  
  new_data$N_obs_count    <- length(keep_count_idx)
  new_data$well_idx_count <- well_map[stan_data$well_idx_count[keep_count_idx]] # Map to new 1:N
  new_data$grid_idx_count <- stan_data$grid_idx_count[keep_count_idx]
  new_data$N_obs          <- stan_data$N_obs[keep_count_idx]
  new_data$D_obs          <- stan_data$D_obs[keep_count_idx]
  
  # --- F. Subset & Re-index Glucose/Lum Observations ---
  keep_gluc_idx <- which(stan_data$well_idx_gluc %in% target_wells_orig_idx)
  
  new_data$N_obs_gluc    <- length(keep_gluc_idx)
  new_data$well_idx_gluc <- well_map[stan_data$well_idx_gluc[keep_gluc_idx]] # Map to new 1:N
  new_data$grid_idx_gluc <- stan_data$grid_idx_gluc[keep_gluc_idx]
  new_data$lum_obs       <- stan_data$lum_obs[keep_gluc_idx]
  new_data$dilution      <- stan_data$dilution[keep_gluc_idx]
  new_data$is_censored   <- stan_data$is_censored[keep_gluc_idx]
  
  # --- G. Clean up obsolete variables ---
  # Explicitly NOT copying: line_id, ploidy_metric, N_lines, structure_mode, line_map
  
  message(paste("Subset complete."))
  message(paste("Line:", ifelse(is.character(target_line), target_line, line_idx)))
  message(paste("Ploidy Metric:", target_ploidy, "(", ploidy_level, ")"))
  message(paste("Wells retained:", new_data$N_wells))
  message(paste("Observations retained:", new_data$N_obs_count + new_data$N_obs_gluc))
  
  return(new_data)
}

stan_data <- subset_stan_data(stan_data,cellLine,ploidy)

run_opt <- function(NUM_THREADS,stan_data,init=2,algorithm="bfgs"){
  # 4. RUN PATHFINDER
  # runs 20 paths, returns the best approximation
  # ---- 3. Run Optimize ----
  opt <- mod$optimize(
    data = stan_data,
    algorithm = "bfgs",
    init = 2, 
    #iter = 100,
    refresh = 50,
    threads = NUM_THREADS
  )
  
  rc <- opt$return_codes()[1]; ok <- rc == 0L
  if (ok) {
    m <- opt$draws(format = "matrix")
    keep <- grepl("^logp\\[|^mu_IC_raw\\[|^phi_total$|^phi_frac$", colnames(m))
    pars <- setNames(as.numeric(m[1, keep]), colnames(m)[keep])
    lp <- opt$lp()
  } else {
    pars <- NaN
    lp <- NaN
  }
  list(retcode = rc, pars = pars, lp = lp)
  
}

library(parallel)

N_CORES  
N_WORKERS <- min(NUM_PATHS,N_CORES)

cl <- makeCluster(N_WORKERS)
clusterExport(cl, "mod", envir = environment())
clusterEvalQ(cl, library(cmdstanr))

opt_smm <- parLapplyLB(
  cl,
  rep(NUM_THREADS,NUM_PATHS),
  run_opt, stan_data = stan_data
)

stopCluster(cl)

outDir <- file.path("data/pathfinder_1cond",MODEL_NAME,cellLine,ploidy)
if(!dir.exists(outDir)) dir.create(outDir,recursive = T)

save_path <- file.path(outDir, paste0("stan_data", ".Rds"))
saveRDS(stan_data,save_path)

save_path <- file.path(outDir, paste0("optim_", MODEL_NAME, ".Rds"))
saveRDS(opt_smm,save_path)

cat(sprintf("\nOptim complete. Saved to %s\n", save_path))
