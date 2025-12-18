library(cmdstanr)
library(dplyr)
library(posterior)
library(jsonlite)

# ==============================================================================
# USER SETTINGS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
# Default to model_B if no arg provided
MODEL_NAME <- if (length(args) >= 1) args[1] else "model_B"
# Default to NOT calculating sim (y_sim) unless flag says otherwise
CALC_SIM   <- if (length(args) >= 2) as.integer(args[2]) else 0 

cat(sprintf(">>> SETTING UP RUN FOR: %s (Sim Output: %s)\n", MODEL_NAME, ifelse(CALC_SIM, "ON", "OFF")))

# ==============================================================================
# 1. LOAD DATA & CONFIG
# ==============================================================================
stan_data <- readRDS("data/stan_ready_data.Rds")
config    <- read_json(paste0("MCMC/",MODEL_NAME, ".json"), simplifyVector = TRUE)

# ==============================================================================
# 2. INJECT PRIORS
# ==============================================================================
# Handle Model B specific dynamic prior for Theta (Capacity)
if (MODEL_NAME == "model_B") {
  max_N <- max(stan_data$N_obs, na.rm = TRUE)
  # Update index 1 (theta) to be log(1.5 * max_obs)
  config$prior_means[1] <- log(max_N * 1.5)
  config$prior_sds[1]   <- 0.5
  cat(sprintf("   -> Updated Theta prior based on data: %.2f (SD: %.2f)\n", config$prior_means[1], config$prior_sds[1]))
}

stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0  # 0 = Sampling Mode
stan_data$calc_sim       <- CALC_SIM

# ==============================================================================
# 3. COMPILATION
# ==============================================================================
stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
if (!file.exists(stan_file)) stop(paste("Stan file not found:", stan_file))

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ==============================================================================
# 4. INITIALIZATION FUNCTION
# ==============================================================================
init_fun <- function() {
  list(
    mu_global   = stan_data$prior_ode_mean,
    sigma_line  = rep(0.05, config$N_params),
    beta_high   = rep(0, config$N_params),
    z_line      = matrix(0, config$N_params, stan_data$N_lines),
    
    mu_IC       = c(stan_data$prior_mu_N0_mean, stan_data$prior_mu_D0_mean),
    sigma_IC    = rep(0.01, 2),
    beta_IC     = rep(0, 2),
    z_IC        = matrix(0, 2, stan_data$N_lines),
    
    calib_sigma = rep(0.1, stan_data$N_exps),
    phi_N       = 10,
    phi_D       = 10
  )
}

# ==============================================================================
# 5. SAMPLE
# ==============================================================================
output_dir <- path.expand(file.path("~/tmp", MODEL_NAME))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists("results")) dir.create("results")

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 15,
  init = init_fun,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 1,
  adapt_delta = 0.90,
  max_treedepth = 12,
  output_dir = output_dir,
  save_warmup=T
)

# ==============================================================================
# 6. SAVE
# ==============================================================================
save_path <- file.path("results", paste0("fit_", MODEL_NAME, ".Rds"))
fit$save_object(file = save_path)
cat(sprintf("\nDone. Saved object to %s\n", save_path))

print(fit$summary(variables = c("mu_global", "phi_N", "phi_D", "calib_sigma")))