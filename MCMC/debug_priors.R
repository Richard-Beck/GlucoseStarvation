library(cmdstanr)
library(jsonlite)

library(cmdstanr)

prior_eval_return_codes <- function(mod, stan_data, n_draws = 500, seed = 1,
                                    threads_per_chain = 1,
                                    capture_params = TRUE,
                                    quiet = TRUE) {
  stan_data$mode <- 2
  stan_data$calc_sim <- 1
  
  rc <- integer(n_draws)
  has_exception <- logical(n_draws)
  ok <- logical(n_draws)
  
  # store params for successful draws only
  params_list <- if (capture_params) vector("list", n_draws) else NULL
  
  for (i in seq_len(n_draws)) {
    fit_i <- mod$sample(
      data = stan_data,
      fixed_param = TRUE,
      chains = 1, parallel_chains = 1,
      threads_per_chain = threads_per_chain,
      iter_warmup = 0,
      iter_sampling = 1,
      seed = seed + i,
      refresh = if (quiet) 0 else 1
    )
    
    rci <- fit_i$return_codes()
    rc[i] <- if (length(rci) == 1 && is.finite(rci)) as.integer(rci) else NA_integer_
    
    out <- paste(fit_i$output(), collapse = "\n")
    has_exception[i] <- grepl("Exception:", out, fixed = TRUE) ||
      grepl("Error evaluating model log probability", out, fixed = TRUE)
    
    ok[i] <- (rc[i] == 0L) && !has_exception[i]
    
    if (capture_params && ok[i]) {
      # 1 iteration, 1 chain -> pull that single draw
      # variables = NULL returns all parameters + gq; you can restrict if you want
      params_list[[i]] <- as.data.frame(fit_i$draws(format = "df"))
      # That df will include lp__, and also generated quantities if any; trim later.
    }
  }
  
  params_df <- if (capture_params) do.call(rbind, params_list[ok]) else NULL
  
  list(rc = rc, has_exception = has_exception, ok = ok, params = params_df)
}

# Re-run only the failing prior draws (same seed convention: seed + i)
# Assumes `res$ok` is length n_draws from the original run.
rerun_failing_prior_draws <- function(mod, stan_data, res, seed = 9000,
                                      threads_per_chain = 1,
                                      keep_fits = FALSE,
                                      refresh = 1) {
  stopifnot(inherits(mod, "CmdStanModel"), is.list(stan_data), is.list(res), !is.null(res$ok))
  fail_idx <- which(!res$ok)
  if (!length(fail_idx)) {
    return(list(fail_idx = integer(0), seeds = integer(0), fits = NULL, output = character(0), return_codes = integer(0)))
  }
  
  # Ensure we are in the same “prior+sim” regime you used
  stan_data$mode <- 2
  stan_data$calc_sim <- 1
  # If you added debug flags in Stan data, enable them here (safe if unused)
  if (!is.null(stan_data$debug_print)) stan_data$debug_print <- 1
  
  fits <- if (keep_fits) vector("list", length(fail_idx)) else NULL
  out  <- character(length(fail_idx))
  rc   <- integer(length(fail_idx))
  seeds <- seed + fail_idx
  
  for (j in seq_along(fail_idx)) {
    i <- fail_idx[j]
    fit_j <- mod$sample(
      data = stan_data,
      fixed_param = TRUE,
      chains = 1, parallel_chains = 1,
      threads_per_chain = threads_per_chain,
      iter_warmup = 0,
      iter_sampling = 1,
      seed = seed + i,
      refresh = refresh
    )
    
    rcj <- fit_j$return_codes()
    rc[j] <- if (length(rcj) == 1 && is.finite(rcj)) as.integer(rcj) else NA_integer_
    out[j] <- paste(fit_j$output(), collapse = "\n")
    if (keep_fits) fits[[j]] <- fit_j
  }
  
  list(
    fail_idx = fail_idx,
    seeds = seeds,
    return_codes = rc,
    output = out,
    fits = fits
  )
}


# ----------------------------
# Example usage with your setup
# ----------------------------
MODEL_NAME <- "model_B"
stan_data <- readRDS("data/stan_ready_data.Rds")
config    <- read_json(paste0("MCMC/", MODEL_NAME, ".json"), simplifyVector = TRUE)
stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds

stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

## run the model on the priors
res <- prior_eval_return_codes(mod, stan_data, n_draws = 500, seed = 9000, threads_per_chain = 1)
table(res$rc)
mean(!res$ok)                 # failure rate
sum(res$has_exception)        # how many had Exception text

fail_runs <- rerun_failing_prior_draws(mod, stan_data, res, seed = 9000, threads_per_chain = 1, keep_fits = FALSE)
fail_runs$fail_idx
fail_runs$seeds
cat(fail_runs$output[1])
