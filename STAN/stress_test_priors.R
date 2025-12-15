library(cmdstanr)
library(dplyr)
library(posterior)

# ==============================================================================
# 1. SETUP & SUBSETTING (The "Single Reactor" simplification)
# ==============================================================================
# We don't need to crash 1000 wells to prove a point. 
# We subset to the first 12 wells to make the crash fast and obvious.

stan_data <- readRDS("data/stan_ready_data.Rds")

# Keep N_lines/N_exps intact (so param dimensions match), but slice the wells
N_sub <- 12 
cat(sprintf("Subsetting data to first %d wells for stress testing...\n", N_sub))
stan_data$calib_a_fixed <- rep(stan_data$prior_calib_a_mean, stan_data$N_exps)
stan_data$calib_b_fixed <- rep(stan_data$prior_calib_b_mean, stan_data$N_exps)
stan_data$N_wells <- N_sub
stan_data$line_id <- stan_data$line_id[1:N_sub]
stan_data$is_high_ploidy <- stan_data$is_high_ploidy[1:N_sub]
stan_data$exp_id <- stan_data$exp_id[1:N_sub]
stan_data$G0_per_well <- stan_data$G0_per_well[1:N_sub]

# Filter observations to match these 12 wells
obs_mask <- stan_data$well_idx_count <= N_sub
stan_data$well_idx_count <- stan_data$well_idx_count[obs_mask]
stan_data$grid_idx_count <- stan_data$grid_idx_count[obs_mask]
stan_data$N_obs <- stan_data$N_obs[obs_mask]
stan_data$D_obs <- stan_data$D_obs[obs_mask]
stan_data$N_obs_count <- sum(obs_mask)

gluc_mask <- stan_data$well_idx_gluc <= N_sub
stan_data$well_idx_gluc <- stan_data$well_idx_gluc[gluc_mask]
stan_data$grid_idx_gluc <- stan_data$grid_idx_gluc[gluc_mask]
stan_data$lum_obs <- stan_data$lum_obs[gluc_mask]
stan_data$dilution <- stan_data$dilution[gluc_mask]
stan_data$N_obs_gluc <- sum(gluc_mask)

# ==============================================================================
# 2. STEP 1: GENERATE VALID PARAMETERS (No ODE)
# ==============================================================================
cat("\nStep 1: Sampling 100 parameter sets from priors (skipping ODE)...\n")

# [cite_start]Ensure prior_only=1 so the model block skips likelihood and GQ skips ODE [cite: 53]
stan_data$prior_only <- 1 

model_path <- "STAN/model_B.stan"
mod <- cmdstan_model(model_path)

# This will succeed because it's just drawing from Normal/Exponential distributions
fit_prior <- mod$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 100, 
  fixed_param = FALSE,
  parallel_chains = 1,
  refresh = 0,
  show_messages = FALSE
)

cat("  -> Draws generated successfully.\n")

# ==============================================================================
# 3. STEP 2: THE "SWITCHEROO" STRESS TEST
# ==============================================================================
cat("\nStep 2: Feeding these priors into the ODE solver...\n")

# MAGIC TRICK: We feed the valid parameters back in, but tell the data 
# that prior_only=0. This forces the 'generated quantities' block 
# [cite_start]to enter the 'else' branch and run 'ode_bdf_tol'[cite: 54].
stan_data$prior_only <- 0 

# We wrap in try() because we expect it to explode
gq_task <- try(mod$generate_quantities(
  fitted_params = fit_prior$output_files(),
  data = stan_data,
  parallel_chains = 1
), silent = TRUE)

# ==============================================================================
# 4. DIAGNOSIS
# ==============================================================================
if (inherits(gq_task, "try-error")) {
  cat("\n[CRASH CONFIRMED] 💥 The ODE solver failed on your priors.\n")
  cat("------------------------------------------------------------\n")
  
  # Extract the error message
  err_msg <- as.character(gq_task)
  
  # Look for the smoking gun
  if (grepl("CVode", err_msg) || grepl("error flag -4", err_msg)) {
    cat("Detected: 'CVode failed with error flag -4'\n")
    cat("Meaning:  Your priors permit parameters that create singularities.\n")
    cat("          (Likely Glucose -> 0 or Rates -> Infinity)\n")
  } else {
    cat("Error output:\n")
    cat(substr(err_msg, 1, 600))
    cat("...\n")
  }
  
  cat("------------------------------------------------------------\n")
  cat("ACTION: You MUST tighten priors or rescale the ODE before fitting.\n")
  cat("        See 'Direction 1' (Rescaling/Logging) in the chat.\n")
  
} else {
  cat("\n[PASSED] ✅ The ODE solver survived 100 prior draws!\n")
  cat("This implies the crashing during fitting is due to the specific posterior\n")
  cat("geometry (funnels) or data conflict, not the fundamental parameterization.\n")
}