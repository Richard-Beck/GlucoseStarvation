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
  stan_data$has_starvation <- stan_data$has_starvation*0
}

stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0 
stan_data$calc_sim       <- CALC_SIM
stan_data$structure_mode <- 0
# 2. COMPILE
stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))


gen_inits <- function(d) {
  # Calculate log bounds locally 
  log_lo <- log(d$lower_b)
  log_up <- log(d$upper_b)
  
  list(
    # FIX: Do NOT exponentiate. mu_global is defined in log-space in the Stan file.
    mu_global = log_lo + (log_up - log_lo) * runif(length(log_lo), 0.25, 0.75),
    
    # Variances (Small start)
    sigma_line    = rep(0.1, 10),
    sigma_beta    = rep(0.1, 10),
    sigma_IC      = rep(0.1, 2),
    sigma_beta_IC = rep(0.1, 2),
    
    # Random effects (Zeros)
    z_line    = matrix(0, 10, d$N_lines),
    z_beta    = matrix(0, 10, d$N_lines),
    z_IC      = matrix(0, 2,  d$N_lines),
    z_beta_IC = matrix(0, 2,  d$N_lines),
    
    # Slopes and Raw ICs
    beta_high = rep(0, 10),
    beta_IC   = rep(0, 2),
    mu_IC_raw = rep(0, 2),
    
    # Noise
    phi_total = 1.0,
    phi_frac  = 10.0
  )
}

# Generate
init_list <- replicate(60, gen_inits(stan_data), simplify = FALSE)

# 4. RUN PATHFINDER
# runs 20 paths, returns the best approximation
fit_pf <- mod$pathfinder(
  data = stan_data,
  seed = 123,
  init = 0.2, 
  num_threads = 60,
  num_paths = 60
)

# 5. SAVE
if (!dir.exists("results")) dir.create("results")
save_path <- file.path("results", paste0("pathfinder_", MODEL_NAME, ".Rds"))
fit_pf$save_object(file = save_path)

cat(sprintf("\nPathfinder complete. Saved to %s\n", save_path))
print(fit_pf$summary(variables = c("mu_global")))