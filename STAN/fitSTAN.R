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

# ----------------------------------------------------------------
# NEW: Inspect init_fun() output (and derived per-well ODE params)
# ----------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

init0 <- init_fun()
saveRDS(init0, "results/init_fun_example.Rds")

# Derive the 10 kinetic params used in the ODE for a few representative cases:
# raw = mu_global + sigma_line*z + beta_high*h; here z=0 in init, so raw = mu_global (+0)
# ODE uses exp(raw) (and for na/nd you may later enforce >=1 etc)
derive_theta <- function(mu_global, sigma_line, z_line, beta_high, h, l) {
  raw <- mu_global + sigma_line * z_line[, l] + beta_high * h
  exp(raw)
}

# Pick up to 5 lines to print
Lpick <- seq_len(min(5, stan_data$N_lines))
for (h in 0:1) {
  cat(sprintf("\nInit-derived exp(raw) kinetics (h=%d) for first %d lines:\n", h, length(Lpick)))
  for (l in Lpick) {
    th <- derive_theta(init0$mu_global, init0$sigma_line, init0$z_line, init0$beta_high, h, l)
    cat(sprintf("  line %d: min=%g max=%g any_inf=%s\n", l, min(th), max(th), any(!is.finite(th))))
  }
}
cat("\nInit mu_global:\n"); print(init0$mu_global)
cat("\nInit sigma_line:\n"); print(init0$sigma_line)
cat("\nInit phi_N, phi_D:\n"); print(c(init0$phi_N, init0$phi_D))

# ----------------------------------------------------------------
# NEW: Prior-only sampling + auto PDF of priors
# Requires model_B.stan to have data int prior_only and likelihood guarded.
# ----------------------------------------------------------------
cat("\nRunning prior-only sampling and saving prior PDF...\n")
stan_data_prior <- stan_data
stan_data_prior$mode <- 1

fit_prior <- mod$sample(
  data = stan_data_prior,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 200,
  adapt_delta = 0.9,
  max_treedepth = 10
)

draws_prior <- fit_prior$draws(format = "draws_df")
# only plot the parameter blocks that define the prior geometry
keep <- grep("^(mu_global|sigma_line|beta_high|mu_IC|sigma_IC|beta_IC|calib_sigma|phi_N|phi_D)(\\[|$)",
             names(draws_prior), value = TRUE)

pdf("results/priors_model_B.pdf", width = 10, height = 7)
chunk_size <- 12
for (i in seq(1, length(keep), by = chunk_size)) {
  p <- keep[i:min(i + chunk_size - 1, length(keep))]
  print(bayesplot::mcmc_dens(draws_prior, pars = p))
}
dev.off()
cat("Wrote results/priors_model_B.pdf\n")


# ----------------------------------------------------------------
# Sample (real fit)
# ----------------------------------------------------------------
stan_data$mode <- 0

N_CHAINS <- 4
N_THREADS_PER_CHAIN <- 15
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
  adapt_delta = 0.99,
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
