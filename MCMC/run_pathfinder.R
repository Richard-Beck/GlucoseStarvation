library(cmdstanr)
library(dplyr)
library(posterior)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME <- if (length(args) >= 1) args[1] else "model_B"
# Usually valid to simulate data during pathfinder to check fit visually
CALC_SIM   <- 1 

cat(sprintf(">>> RUNNING PATHFINDER (Approx) FOR: %s\n", MODEL_NAME))

# 1. LOAD & PREP (Identical logic to run_fit.R)
stan_data <- readRDS("data/stan_ready_data.Rds")
config    <- read_json(paste0("MCMC/",MODEL_NAME, ".json"), simplifyVector = TRUE)

if (MODEL_NAME == "model_B") {
  max_N <- max(stan_data$N_obs, na.rm = TRUE)
  config$prior_means[1] <- log(max_N * 1.5)
  config$prior_sds[1]   <- 0.5
  stan_data$lower_b <- config$lower
  stan_data$upper_b <- config$upper
  stan_data$N_params <- length(config$lower)
}

stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0 
stan_data$calc_sim       <- CALC_SIM


# 2. COMPILE
stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# 3. INIT (Single point init for pathfinder is often sufficient, or use function)
# ==============================================================================
# 3. INITIALIZATION (Full Prior Sampling)
# ==============================================================================
init_fun <- function() {
  # 1. ODE Parameters (mu_global): Sample from JSON Config Priors
  #    mu_global ~ normal(prior_mean, prior_sd)
  jit_mu <- rnorm(config$N_params, 
                  mean = stan_data$prior_ode_mean, 
                  sd = stan_data$prior_ode_sd)
  
  # 2. Sigmas: Sample from Exponential(1) 
  #    sigma_line ~ exponential(1)
  jit_sigma_line <- rexp(config$N_params, rate = 1)
  
  # 3. Random Effects (z_line): Sample from Standard Normal
  #    to_vector(z_line) ~ std_normal()
  jit_z_line <- matrix(rnorm(config$N_params * stan_data$N_lines, 0, 1), 
                       nrow = config$N_params, ncol = stan_data$N_lines)
  
  # 4. WGD Effect (beta): Sample from Normal(0, 1)
  #    beta_high ~ normal(0, 1)
  jit_beta_high <- rnorm(config$N_params, 0, 1)
  
  # 5. ICs: Sample from Data-Derived Priors
  #    mu_IC ~ normal(prior_mu_N0/D0, ...)
  jit_mu_IC <- c(
    rnorm(1, stan_data$prior_mu_N0_mean, stan_data$prior_mu_N0_sd),
    rnorm(1, stan_data$prior_mu_D0_mean, stan_data$prior_mu_D0_sd)
  )
  
  # 6. IC Variance: Sample from Exponential(1)
  #    sigma_IC ~ exponential(1)
  jit_sigma_IC <- rexp(2, rate = 1)
  
  # 7. IC Random Effects: Sample from Standard Normal
  jit_z_IC <- matrix(rnorm(2 * stan_data$N_lines, 0, 1), nrow = 2)
  
  # 8. IC WGD Effect: Sample from Normal(0, 1)
  jit_beta_IC <- rnorm(2, 0, 1)
  
  # 9. Nuisance / Calibration
  #    calib_sigma ~ exponential(1)
  #    phi ~ exponential(0.1)
  jit_calib_sigma <- rexp(stan_data$N_exps, rate = 1)
  jit_phi_N       <- rexp(1, rate = 0.1)
  jit_phi_D       <- rexp(1, rate = 0.1)
  
  list(
    mu_global   = jit_mu,
    sigma_line  = jit_sigma_line,
    beta_high   = jit_beta_high,
    z_line      = jit_z_line,
    mu_IC       = jit_mu_IC,
    sigma_IC    = jit_sigma_IC,
    beta_IC     = jit_beta_IC,
    z_IC        = jit_z_IC,
    calib_sigma = jit_calib_sigma,
    phi_N       = jit_phi_N,
    phi_D       = jit_phi_D
  )
}

# 4. RUN PATHFINDER
# runs 20 paths, returns the best approximation
fit_pf <- mod$pathfinder(
  data = stan_data,
#  init = init_fun,
  seed = 123,
  num_threads = 60,
  num_paths = 36
)

# 5. SAVE
if (!dir.exists("results")) dir.create("results")
save_path <- file.path("results", paste0("pathfinder_", MODEL_NAME, ".Rds"))
fit_pf$save_object(file = save_path)

cat(sprintf("\nPathfinder complete. Saved to %s\n", save_path))
print(fit_pf$summary(variables = c("mu_global")))