library(cmdstanr)
library(jsonlite)
library(data.table)

# ---- Settings ----
MODEL_NAME <- "model_B"
ALG <- "bfgs"
SEED <- 123
THREADS <- 3
ITERS_PRE  <- 1000 
ITERS_FINAL <- 2000

# ---- ROBUST PATH HANDLING ----
# 1. Try a local relative path first
OUTDIR <- "results/opt_chain"

# 2. Try to create it. If this fails (permissions), we catch it.
tryCatch({
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
}, warning = function(w) {
  message("Warning: Could not create local results directory. Switching to temporary directory.")
  OUTDIR <<- file.path(tempdir(), "opt_chain")
  dir.create(OUTDIR, recursive = TRUE)
}, error = function(e) {
  message("Error: Could not create local results directory. Switching to temporary directory.")
  OUTDIR <<- file.path(tempdir(), "opt_chain")
  dir.create(OUTDIR, recursive = TRUE)
})

# 3. Final check: Can we actually write here?
if (file.access(OUTDIR, mode = 2) != 0) {
  message("Warning: No write access to ", OUTDIR, ". Switching to temporary directory.")
  OUTDIR <- file.path(tempdir(), "opt_chain")
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
}

# 4. Convert to absolute path for CmdStan to be happy
OUTDIR <- normalizePath(OUTDIR)
print(paste("Writing output to:", OUTDIR))

# ---- Data Prep ----
stan_data <- readRDS("data/stan_ready_data.Rds")
config <- read_json(paste0("MCMC/", MODEL_NAME, ".json"), simplifyVector = TRUE)

# Config patches
stan_data$lower_b  <- config$lower
stan_data$upper_b  <- config$upper
stan_data$N_params <- length(config$lower)
stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$calc_sim       <- 0 
stan_data$log_lower <- log(stan_data$lower_b)
stan_data$log_upper <- log(stan_data$upper_b)

# ---- Compile ----
mod <- cmdstan_model(paste0("MCMC/", MODEL_NAME, ".stan"), 
                     cpp_options = list(stan_threads = TRUE))

# ---- Helper: Full Init Generator ----
gen_base_inits <- function(d) {
  log_lo <- d$log_lower
  log_up <- d$log_upper
  
  list(
    mu_global = log_lo + (log_up - log_lo) * runif(length(log_lo), 0.25, 0.75),
    sigma_line    = rep(0.1, 10),
    sigma_beta    = rep(0.1, 10),
    sigma_IC      = rep(0.1, 2),
    sigma_beta_IC = rep(0.1, 2),
    z_line    = matrix(0, 10, d$N_lines),
    z_beta    = matrix(0, 10, d$N_lines),
    z_IC      = matrix(0, 2,  d$N_lines),
    z_beta_IC = matrix(0, 2,  d$N_lines),
    beta_high = rep(0, 10),
    beta_IC   = rep(0, 2),
    mu_IC_raw = rep(0, 2),
    phi_total = 1.0,
    phi_frac  = 10.0
  )
}

# ---- Robust Parameter Extractor ----
# Now accepts manual paths in case the fit object crashes
get_best_params <- function(fit_obj, out_dir = NULL, base_name = NULL) {
  
  csv_files <- fit_obj$output_files()
  
  # FALLBACK: If optimization crashed, cmdstanr loses the file handle.
  # We manually reconstruct the filename: <basename>-1.csv
  if (length(csv_files) == 0) {
    if (!is.null(out_dir) && !is.null(base_name)) {
      manual_file <- file.path(out_dir, paste0(base_name, "-1.csv"))
      if (file.exists(manual_file)) {
        message("Optimization crashed, but found CSV manually: ", manual_file)
        csv_files <- manual_file
      } else {
        stop("Optimization failed and no CSV found at: ", manual_file)
      }
    } else {
      stop("Optimization failed and no CSV output attached to object.")
    }
  }
  
  # Read CSV
  headers <- names(fread(csv_files[1], skip = "lp__", nrows = 0))
  drop_cols <- grep("^y_sim|__$", headers, value = TRUE)
  dt <- fread(csv_files[1], skip = "lp__", drop = drop_cols, fill = TRUE)
  
  # Sort by lp__ descending to get the best iteration (even if crash happened at end)
  dt <- dt[order(-lp__)]
  best_row <- dt[1]
  
  # Parse Parameters
  out_list <- list()
  param_names <- names(best_row)
  param_names <- param_names[!grepl("^theta_", param_names)] # remove tuple dupes
  
  bases <- unique(sapply(strsplit(param_names, "\\."), `[`, 1))
  bases <- setdiff(bases, "lp__")
  
  for (b in bases) {
    cols <- grep(paste0("^", b, "(\\.|$)"), param_names, value = TRUE)
    
    if (length(cols) == 1 && !grepl("\\.", cols)) {
      out_list[[b]] <- as.numeric(best_row[[b]])
    } else {
      meta <- lapply(cols, function(c) {
        parts <- as.numeric(unlist(strsplit(sub(paste0(b, "\\."), "", c), "\\.")))
        list(idx = parts, val = as.numeric(best_row[[c]]))
      })
      
      indices <- do.call(rbind, lapply(meta, function(x) x$idx))
      
      if (ncol(indices) == 1) {
        arr <- numeric(max(indices))
        for(m in meta) arr[m$idx] <- m$val
        out_list[[b]] <- arr
      } else {
        max_dims <- apply(indices, 2, max)
        arr <- array(0, dim = max_dims)
        for(m in meta) arr[matrix(m$idx, nrow=1)] <- m$val
        out_list[[b]] <- arr
      }
    }
  }
  return(out_list)
}

# ---- Re-Define Init Generator (Just in case) ----
gen_base_inits <- function(d) {
  log_lo <- d$log_lower
  log_up <- d$log_upper
  
  list(
    mu_global     = log_lo + (log_up - log_lo) * runif(length(log_lo), 0.25, 0.75),
    sigma_line    = rep(0.1, 10),
    sigma_beta    = rep(0.1, 10),
    sigma_IC      = rep(0.1, 2),
    sigma_beta_IC = rep(0.1, 2),
    z_line        = matrix(0, 10, d$N_lines),
    z_beta        = matrix(0, 10, d$N_lines),
    z_IC          = matrix(0, 2,  d$N_lines),
    z_beta_IC     = matrix(0, 2,  d$N_lines),
    beta_high     = rep(0, 10),
    beta_IC       = rep(0, 2),
    mu_IC_raw     = rep(0, 2),
    phi_total     = 1.0,
    phi_frac      = 10.0
  )
}

# ==============================================================================
# STEP 1: Run Mode 2 (Rigid Slope)
# ==============================================================================
print("=== STARTING MODE 2 (Rigid Slope) ===")
stan_data$structure_mode <- 2
stan_data$mode <- 0

init_m2 <- gen_base_inits(stan_data)
bn_2 <- paste0(MODEL_NAME, "_step1_mode2")

# We ignore warnings here because we expect LS failure might happen, 
# but we want to salvage the partial result.
fit_2 <- mod$optimize(
  data = stan_data, algorithm = ALG, seed = SEED, 
  init = list(init_m2), iter = ITERS_PRE, threads = THREADS,
  init_alpha = 0.0001,  # <--- Smaller step size to prevent LS failure
  output_dir = OUTDIR, output_basename = bn_2
)

# Pass the manual paths so we can find the CSV even if fit_2 crashed
params_2 <- get_best_params(fit_2, out_dir = OUTDIR, base_name = bn_2)

# If we extracted parameters, we are good to go, even if LP is bad.
print(paste("Mode 2 Params Extracted. Proceeding to Step 2."))

# ==============================================================================
# STEP 2: Map Mode 2 -> Mode 0 (Hierarchical)
# ==============================================================================
print("=== STARTING MODE 0 (Hierarchical) ===")
init_m0 <- params_2

# Reset hierarchical variances to small (allow growth)
init_m0$sigma_beta    <- rep(0.01, 10)
init_m0$sigma_beta_IC <- rep(0.01, 2)
init_m0$z_beta        <- matrix(0, 10, stan_data$N_lines)
init_m0$z_beta_IC     <- matrix(0, 2, stan_data$N_lines)

stan_data$structure_mode <- 0
bn_0 <- paste0(MODEL_NAME, "_step2_mode0")

fit_0 <- mod$optimize(
  data = stan_data, algorithm = ALG, seed = SEED, 
  init = list(init_m0), iter = ITERS_PRE, threads = THREADS,
  init_alpha = 0.001, # Back to standard step size
  output_dir = OUTDIR, output_basename = bn_0
)

params_0 <- get_best_params(fit_0, out_dir = OUTDIR, base_name = bn_0)
print(paste("Mode 0 Finished. Proceeding to Step 3."))

# ==============================================================================
# STEP 3: Map Mode 0 -> Mode 1 (Fully Independent)
# ==============================================================================
print("=== STARTING MODE 1 (Fully Free) ===")
init_m1 <- params_0

# Transform: z_1 = (sigma * z_0) / 5.0
init_m1$z_line <- (params_0$sigma_line * params_0$z_line) / 5.0
init_m1$z_IC   <- (params_0$sigma_IC * params_0$z_IC) / 5.0

stan_data$structure_mode <- 1
bn_1 <- paste0(MODEL_NAME, "_step3_mode1")

fit_1 <- mod$optimize(
  data = stan_data, algorithm = ALG, seed = SEED, 
  init = list(init_m1), iter = ITERS_FINAL, threads = THREADS,
  output_dir = OUTDIR, output_basename = bn_1
)

print("=== CHAIN COMPLETE ===")
# Check final result
final_params <- get_best_params(fit_1, out_dir = OUTDIR, base_name = bn_1)
print("Final extraction successful.")