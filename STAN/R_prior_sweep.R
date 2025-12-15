suppressPackageStartupMessages({
  library(data.table)
})

# ---- Inputs ----
N <- 100L
data_rds <- "data/stan_ready_data.Rds"  # adjust if needed
out_csv  <- sprintf("results/prior_ode_sweep_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
out_rds  <- sub("\\.csv$", ".fail_params.Rds", out_csv)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

# ---- Load data & helper ----
stan_data <- readRDS(data_rds)
source("STAN/model_B.R")

# sanity: ensure fixed calibration present (required by log_lik_model)
if (is.null(stan_data$calib_a_fixed) || is.null(stan_data$calib_b_fixed)) {
  stop("stan_data$calib_a_fixed / calib_b_fixed missing. Add them during preprocessing.")
}

set.seed(42)

res <- vector("list", N)
fail_params <- list()

for (i in seq_len(N)) {
  status <- "ok"
  ll <- NA_real_
  msg <- NA_character_
  
  # 1) prior draw
  p <- draw_from_prior(stan_data)
  
  # 2) simulate + 3) loglik
  yhat <- NULL
  tryCatch({
    yhat <- simulate_full_model(p, stan_data)
    ll   <- log_lik_model(p, stan_data, yhat)
  }, error = function(e){
    status <<- "ode_fail"
    msg    <<- conditionMessage(e)
  }, warning = function(w){
    # record warnings but still take ll if available
    if (is.na(msg)) msg <<- conditionMessage(w)
  })
  
  # If not ODE failure but something else went wrong
  if (status == "ok" && !is.numeric(ll) || is.na(ll)) {
    status <- "other_fail"
    if (is.na(msg)) msg <- "non-numeric or NA ll"
  }
  
  # store results
  res[[i]] <- data.frame(iter = i, status = status, log_lik = ll, note = msg, stringsAsFactors = FALSE)
  if (status != "ok") {
    fail_params[[length(fail_params) + 1L]] <- list(iter = i, params = p, note = msg)
  }
}

cat("Wrote: ", out_csv, "\n", sep="")
cat("Wrote: ", out_rds, "\n", sep="")

# quick console summary
cat("\nSummary:\n")
print(tab[, .(n = .N, frac = .N / .N[1L]), by = status])
if (any(tab$status == "ok")) {
  ok <- tab[status == "ok", log_lik]
  cat(sprintf("OK draws: %d, mean ll = %.3f, median ll = %.3f\n",
              length(ok), mean(ok), median(ok)))
}