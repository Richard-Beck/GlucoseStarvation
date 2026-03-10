library(cmdstanr)
library(posterior)
library(jsonlite)
library(parallel)

source("ecology/subset_stan_data.R")

args <- commandArgs(trailingOnly = TRUE)

# ---- Argument Parsing ----
# 1. Model Name
MODEL_NAME      <- if (length(args) >= 1) args[1] else "cobbdoug"
# 2. Constraint Name (New)
CONSTRAINT_NAME <- if (length(args) >= 2) args[2] else "unconstrained" 
# 3. Cell Line
cellLine        <- if (length(args) >= 3) args[3] else "SUM-159-fuse"
# 4. Ploidy
ploidy          <- if (length(args) >= 4) args[4] else "hi"
# 5. Starts
N_STARTS        <- if (length(args) >= 5) as.integer(args[5]) else 4L
# 6. Cores
N_CORES         <- if (length(args) >= 6) as.integer(args[6]) else 1L
# 7. Threads (New)
NUM_THREADS     <- if (length(args) >= 7) as.integer(args[7]) else 4L

cat(sprintf(">>> OPTIM: %s | %s | %s | %s | starts=%d | cores=%d | threads=%d\n",
            MODEL_NAME, CONSTRAINT_NAME, cellLine, ploidy, N_STARTS, N_CORES, NUM_THREADS))

# ---- Define Constraints ----
source(file.path("ecology/model_selection/models",MODEL_NAME,"constraints.R"))

# Validate constraint choice
if (!CONSTRAINT_NAME %in% names(constraints_list)) {
  stop(sprintf("Constraint '%s' not found in constraints_list.", CONSTRAINT_NAME))
}
active_constraint <- constraints_list[[CONSTRAINT_NAME]]


# ---- Load Prepared Data ----
stan_data <- readRDS("data/stan_ready_data.Rds")
data <- subset_stan_data(stan_data, cellLine, ploidy)
data$R2_0 <- 1.0
# Assume these priors are hardcoded in the model file
source(file.path("ecology/model_selection/models",MODEL_NAME,paste0(MODEL_NAME,".R")))
data$prior_centers <- prior_centers
data$prior_scales <- prior_scales

# 4. Apply Constraint Map & Mask
data$param_map  <- active_constraint$map
data$param_mask <- active_constraint$mask

# ---- Compile ----
stan_file <-file.path("ecology/model_selection/models",MODEL_NAME,paste0(MODEL_NAME,".stan"))
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- Define Parallel Wrapper ----
# ---- Define Parallel Wrapper ----
run_one_start <- function(seed_idx) {
  
  # --- Timeout Configuration ---
  TIMEOUT_SEC <- 60  # e.g., 10 minutes per start
  
  # Set the limit. 'transient = TRUE' ensures it applies only to the current calculation
  # and not the rest of the R session if an error occurs.
  setTimeLimit(elapsed = TIMEOUT_SEC, transient = TRUE)
  
  # Ensure the limit is strictly unset when the function exits (success or fail)
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)
  
  # --- Optimization Block ---
  tryCatch({
    
    # Optimize with random init
    opt <- mod$optimize(
      data = data,
      algorithm = "bfgs",
      init = 3, 
      refresh = 0,
      threads = NUM_THREADS 
    )
    
    rc <- opt$return_codes()[1]
    
    if (rc == 0L) {
      return(list(rc = rc, lp = opt$lp(), draws = opt$draws(format = "matrix")))
    } else {
      return(list(rc = rc, lp = NA_real_, draws = NULL))
    }
    
  }, error = function(e) {
    # This block runs if the timeout is hit OR if Stan crashes
    # We return a failure code (e.g., 99) so the pipeline continues
    return(list(rc = 99L, lp = NA_real_, draws = NULL, error_msg = e$message))
  })
}

message(sprintf("Starting %d optim starts across %d cores (Stan threads: %d)...", 
                N_STARTS, N_CORES, NUM_THREADS))

# ---- Run Parallel ----
results <- mclapply(seq_len(N_STARTS), run_one_start, mc.cores = N_CORES)

# ---- Unpack Results ----
rc_vec <- vapply(results, function(x) x$rc, integer(1))
lp_vec <- vapply(results, function(x) if(is.null(x$lp)) NA_real_ else x$lp, numeric(1))
opt_list <- lapply(results, function(x) x$draws)

# ---- Save Output ----
# Path now includes CONSTRAINT_NAME
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, CONSTRAINT_NAME, cellLine, ploidy)
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

saveRDS(data,       file.path(outDir, "stan_data.Rds"))
saveRDS(opt_list,   file.path(outDir, "optim_draws_all.Rds"))
saveRDS(lp_vec,     file.path(outDir, "optim_lp_all.Rds"))
saveRDS(rc_vec,     file.path(outDir, "optim_rc_all.Rds"))


