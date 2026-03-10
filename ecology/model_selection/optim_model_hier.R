library(cmdstanr)
library(posterior)
library(jsonlite)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

# ---- Argument Parsing ----
# 1. Model Name
MODEL_NAME      <- if (length(args) >= 1) args[1] else "pathways"
# 2. Constraint Name
CONSTRAINT_NAME <- if (length(args) >= 2) args[2] else "glucose_glutamine_limiting" 
# 3. Starts
N_STARTS        <- if (length(args) >= 3) as.integer(args[3]) else 4L
# 4. Cores
N_CORES         <- if (length(args) >= 4) as.integer(args[4]) else 1L
# 5. Threads
NUM_THREADS     <- if (length(args) >= 5) as.integer(args[5]) else 4L

cat(sprintf(">>> OPTIM HIERARCHICAL: %s | %s | starts=%d | cores=%d | threads=%d\n",
            MODEL_NAME, CONSTRAINT_NAME, N_STARTS, N_CORES, NUM_THREADS))

# ---- Define Constraints ----
source(file.path("ecology/model_selection/models", MODEL_NAME, "constraints.R"))

if (!CONSTRAINT_NAME %in% names(constraints_list)) {
  stop(sprintf("Constraint '%s' not found in constraints_list.", CONSTRAINT_NAME))
}
active_constraint <- constraints_list[[CONSTRAINT_NAME]]

# ---- Load Data ----
stan_data <- readRDS("data/stan_ready_data.Rds")

# ---- Modify Data on the Fly ----
# 1. Create Group ID for N0 (Unique Line x Ploidy combo)
grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
unique_grps <- unique(grp_keys)

stan_data$N_groups <- length(unique_grps)
stan_data$group_id <- match(grp_keys, unique_grps)

# 2. Add other constants
stan_data$R2_0 <- 1.0

# 3. Load Model-Specific Priors
source(file.path("ecology/model_selection/models", MODEL_NAME, paste0(MODEL_NAME, ".R")))
stan_data$prior_centers <- prior_centers
stan_data$prior_scales  <- prior_scales

# 4. Apply Constraint Map & Mask
stan_data$param_map  <- active_constraint$map
stan_data$param_mask <- active_constraint$mask

stan_data <- stan_data[!sapply(stan_data,is.character)]

# ---- Compile ----
stan_file <- file.path("ecology/model_selection/models", MODEL_NAME, paste0(MODEL_NAME, "_hier.stan"))
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- Define Parallel Wrapper (Size Diagnostic) ----
run_one_start <- function(seed_idx) {
  
  pfx <- sprintf("[Worker %d]:", seed_idx)
  
  # 1. Log Start
  cat(sprintf("%s STARTING wrapper. \n", pfx), file = stderr())
  
  TIMEOUT_SEC <- 600  
  setTimeLimit(elapsed = TIMEOUT_SEC, transient = TRUE)
  
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)
  
  tryCatch({
    
    cat(sprintf("%s Calling mod$optimize...\n", pfx), file = stderr())
    
    # Run Optimization
    opt <- mod$optimize(
      data = stan_data,
      algorithm = "bfgs", 
      init = 2, 
      refresh = 0,       
      threads = NUM_THREADS 
    )
    
    rc <- opt$return_codes()[1]
    
    if (rc == 0L) {
      cat(sprintf("%s RC is 0. Extracting...\n", pfx), file = stderr())
      
      val_lp <- opt$lp()
      val_draws <- opt$draws(format = "matrix")
      
      # --- DIAGNOSTIC: CHECK SIZE ---
      size_draws <- object.size(val_draws)
      cat(sprintf("%s 'draws' matrix size: %s\n", pfx, format(size_draws, units="auto")), file = stderr())
      
      res <- list(rc = rc, lp = val_lp, draws = val_draws)
      size_list <- object.size(res)
      cat(sprintf("%s Full return list size: %s\n", pfx, format(size_list, units="auto")), file = stderr())
      
      # If it's huge, warn explicitly
      if (size_list > 2e9) { # > 2GB
        cat(sprintf("%s CRITICAL WARNING: Object is >2GB. R serialization will likely fail.\n", pfx), file = stderr())
      }
      
      return(res)
      
    } else {
      return(list(rc = rc, lp = NA_real_, draws = NULL))
    }
    
  }, error = function(e) {
    cat(sprintf("%s ERROR: %s\n", pfx, e$message), file = stderr())
    return(list(rc = 99L, lp = NA_real_, draws = NULL, error_msg = e$message))
  })
}

message(sprintf("Starting %d optim starts across %d cores (Stan threads: %d)...", 
                N_STARTS, N_CORES, NUM_THREADS))

# ---- Run Parallel ----
results <- mclapply(seq_len(N_STARTS), run_one_start, mc.cores = N_CORES)

print("otpimization complete")

# ---- Unpack Results ----
rc_vec <- vapply(results, function(x) x$rc, integer(1))
lp_vec <- vapply(results, function(x) if(is.null(x$lp) || length(x$lp)==0) NA_real_ else x$lp[1], numeric(1))
opt_list <- lapply(results, function(x) x$draws)

# ---- Save Output ----
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, CONSTRAINT_NAME, "hier")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("Saving results to %s...", 
                outDir))

saveRDS(stan_data,  file.path(outDir, "stan_data.Rds"))
saveRDS(opt_list,   file.path(outDir, "optim_draws_all.Rds"))
saveRDS(lp_vec,     file.path(outDir, "optim_lp_all.Rds"))
saveRDS(rc_vec,     file.path(outDir, "optim_rc_all.Rds"))

cat(">>> Done.\n")