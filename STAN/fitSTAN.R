library(cmdstanr)
library(dplyr)
library(posterior)
library(bayesplot)

color_scheme_set("mix-blue-pink")

cat("Loading Data...\n")
stan_data <- readRDS("data/stan_ready_data.Rds")
if (is.null(stan_data$prior_ode_mean)) {
  stop("stan_data is missing 'prior_ode_mean'. Please re-run prepare_stan_data.R!")
}
# ----------------------------------------------------------------
# Compile
# ----------------------------------------------------------------
cat("Compiling Model...\n")
model_path <- "STAN/model_B.stan"
mod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

# ----------------------------------------------------------------
# Init
# ----------------------------------------------------------------
init_fun <- function() {
  list(
    # [ACTION] Point the ODE parameters to the safe, calculated means
    mu_global = stan_data$prior_ode_mean, 
    
    # [ACTION] Initialize variability to be small (start chains close together)
    sigma_line = rep(0.01, 10),
    beta_high  = rep(0, 10),
    z_line     = matrix(0, 10, stan_data$N_lines),
    
    # [ACTION] Point ICs to their calculated means
    mu_IC    = c(stan_data$prior_mu_N0_mean, stan_data$prior_mu_D0_mean),
    sigma_IC = rep(0.01, 2),
    beta_IC  = rep(0, 2),
    z_IC     = matrix(0, 2, stan_data$N_lines),
    
    # Keep the rest
    calib_sigma = rep(0.1, stan_data$N_exps),
    phi_N = 10,
    phi_D = 10
  )
}

if (!dir.exists("results")) dir.create("results")


# ----------------------------------------------------------------
# Sample (real fit)
# ----------------------------------------------------------------
stan_data$mode <- 0

N_CHAINS <- 1
N_THREADS_PER_CHAIN <- 3
cat(sprintf("\nRunning %d chains with %d threads each.\n", N_CHAINS, N_THREADS_PER_CHAIN))

output_dir <- path.expand("~/tmp")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fit <- mod$sample(
  data = stan_data,
  seed = 42,
  chains = N_CHAINS,
  parallel_chains = N_CHAINS,
  threads_per_chain = N_THREADS_PER_CHAIN,
  init = init_fun,
  iter_warmup = 1500,
  iter_sampling = 1000,
  refresh = 1,
  adapt_delta = 0.9,
  max_treedepth = 12,
  output_dir = output_dir,
  save_warmup = TRUE
)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------
cat("Saving fit object...\n")
fit$save_object(file = "results/fit_model_B.Rds")

cat("Done. Summary:\n")
print(fit$summary(variables = c("mu_global", "beta_high", "phi_N", "phi_D", "calib_sigma")))
