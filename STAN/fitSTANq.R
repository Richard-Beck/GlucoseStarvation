library(cmdstanr)
library(dplyr)
library(posterior)
library(bayesplot)

color_scheme_set("mix-blue-pink")

cat("Loading Data...\n")
stan_data <- readRDS("data/stan_ready_data.Rds")

# ----------------------------------------------------------------
# Define Quota Model Priors (Overwriting generic priors)
# ----------------------------------------------------------------
logit <- function(p) log(p/(1-p))

prior_means <- c(
  20.0,            # 1 log_theta
  -2.91,           # 2 log_kp
  -4.34,           # 3 log_kd
  -9.96,           # 4 log_kd2
  0.0,             # 5 log_nu   (nu ~ 1 to start)
  logit(0.10),     # 6 rho_logit (alpha_sum = nu * rho)
  logit(0.10),     # 7 r_logit   (alpha_c share of alpha_sum)
  7.6,             # 8 log_sigma_G
  0.0,             # 9 log_KG
  logit(0.95),     # 10 q50g_frac_logit (q50gN = q50g_frac * q_star)
  logit(0.05),     # 11 q50d_frac_logit (q50dN = q50d_frac * q_star)
  2.3,             # 12 log_na
  2.3              # 13 log_nd
)

prior_sds <- c(
  2.0,   # 1 log_theta
  1.0,   # 2 log_kp
  1.0,   # 3 log_kd
  2.0,   # 4 log_kd2
  1.0,   # 5 log_nu
  1.0,   # 6 rho_logit (broad on logit scale)
  1.0,   # 7 r_logit   (broad on logit scale)
  1.0,   # 8 log_sigma_G
  1.0,   # 9 log_KG
  1.0,   # 10 q50g_frac_logit
  1.0,   # 11 q50d_frac_logit
  0.5,   # 12 log_na
  0.5    # 13 log_nd
)


# Inject into stan_data
stan_data$prior_ode_mean <- prior_means
stan_data$prior_ode_sd   <- prior_sds

# ----------------------------------------------------------------
# Compile
# ----------------------------------------------------------------
cat("Compiling Model...\n")
model_path <- "STAN/model_q.stan"
mod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))

# ----------------------------------------------------------------
# Init
# ----------------------------------------------------------------
N_PARAMS <- 13

init_fun <- function() {
  list(
    # Global parameters
    mu_global = stan_data$prior_ode_mean, 
    
    # Variability (start small)
    sigma_line = rep(0.05, N_PARAMS),
    beta_high  = rep(0, N_PARAMS),
    z_line     = matrix(0, N_PARAMS, stan_data$N_lines),
    
    # ICs (u0, v0, q0 is fixed to 1.0 in code usually, but we fit u0/v0)
    mu_IC    = c(stan_data$prior_mu_N0_mean, stan_data$prior_mu_D0_mean),
    sigma_IC = rep(0.01, 2),
    beta_IC  = rep(0, 2),
    z_IC     = matrix(0, 2, stan_data$N_lines),
    
    # Noise/Overdispersion
    calib_sigma = rep(0.1, stan_data$N_exps),
    phi_N = 10,
    phi_D = 10
  )
}

if (!dir.exists("results")) dir.create("results")

# ----------------------------------------------------------------
# Sample
# ----------------------------------------------------------------
stan_data$mode <- 0

N_CHAINS <- 1
N_THREADS_PER_CHAIN <- 3
cat(sprintf("\nRunning %d chains with %d threads each.\n", N_CHAINS, N_THREADS_PER_CHAIN))

output_dir <- path.expand("~/tmp_quota")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fit <- mod$sample(
  data = stan_data,
  seed = 42,
  chains = N_CHAINS,
  parallel_chains = N_CHAINS,
  threads_per_chain = N_THREADS_PER_CHAIN,
  init = init_fun,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 1,
  adapt_delta = 0.95,       # Higher delta for stiff quota dynamics
  max_treedepth = 10,
  output_dir = output_dir,
  save_warmup = TRUE
)

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------
cat("Saving fit object...\n")
fit$save_object(file = "results/fit_model_q.Rds")

cat("Done. Summary:\n")
print(fit$summary(variables = c("mu_global", "phi_N", "phi_D")))