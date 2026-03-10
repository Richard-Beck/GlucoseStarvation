library(cmdstanr)
library(posterior)
library(jsonlite)
library(parallel)

source("ecology/subset_stan_data.R")

args <- commandArgs(trailingOnly = TRUE)

# ---- Argument Parsing ----
# 1. Model Name
MODEL_NAME      <- if (length(args) >= 1) args[1] else "pathways"
# 2. Constraint Name
CONSTRAINT_NAME <- if (length(args) >= 2) args[2] else "glucose_simple"
# 3. Cell Line
cellLine        <- if (length(args) >= 3) args[3] else "MDA-MB-231"
# 4. Ploidy
ploidy          <- if (length(args) >= 4) args[4] else "lo"

# 5. Chains (replaces N_STARTS)
N_CHAINS        <- if (length(args) >= 5) as.integer(args[5]) else 4L
# 6. Cores (for parallel chains)
N_CORES         <- if (length(args) >= 6) as.integer(args[6]) else 4L
# 7. Threads per chain (Stan intra-chain threading)
THREADS_PER_CHAIN <- if (length(args) >= 7) as.integer(args[7]) else 1L

cat(sprintf(">>> NUTS: %s | %s | %s | %s | chains=%d | cores=%d | threads/chain=%d\n",
            MODEL_NAME, CONSTRAINT_NAME, cellLine, ploidy, N_CHAINS, N_CORES, THREADS_PER_CHAIN))

# ---- Define Constraints ----
source(file.path("ecology/model_selection/models",MODEL_NAME,"constraints.R"))

if (!CONSTRAINT_NAME %in% names(constraints_list)) {
  stop(sprintf("Constraint '%s' not found.", CONSTRAINT_NAME))
}
active_constraint <- constraints_list[[CONSTRAINT_NAME]]

# ---- Load Prepared Data ----
stan_data <- readRDS("data/stan_ready_data.Rds")
data <- subset_stan_data(stan_data, cellLine, ploidy)

# ---- Add Model Constants & Priors ----
data$R2_0 <- 1.0
source(file.path("ecology/model_selection/models",MODEL_NAME,paste0(MODEL_NAME,".R")))
data$prior_centers <- prior_centers
data$prior_scales <- prior_scales
# 4. Apply Constraint Map & Mask
data$param_map  <- active_constraint$map
data$param_mask <- active_constraint$mask

# ---- Compile ----
stan_file <-file.path("ecology/model_selection/models",MODEL_NAME,paste0(MODEL_NAME,".stan"))
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- Prepare Inits from Optimization ----
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, CONSTRAINT_NAME, cellLine, ploidy)

optim_draws_file <- file.path(outDir, "optim_draws_all.Rds")
optim_lp_file    <- file.path(outDir, "optim_lp_all.Rds")

if (!file.exists(optim_draws_file) || !file.exists(optim_lp_file)) {
  stop(sprintf("Optimization results not found in %s. Please run the optimization script first.", outDir))
}

# Load results
opt_draws_all <- readRDS(optim_draws_file)
opt_lp_all    <- readRDS(optim_lp_file)

# Filter for valid runs
valid_indices <- which(!is.na(opt_lp_all))
if (length(valid_indices) == 0) stop("No valid optimization runs found to use as inits.")

valid_lps   <- opt_lp_all[valid_indices]
valid_draws <- opt_draws_all[valid_indices]

# Sort by best LP (descending)
sorted_idx <- order(valid_lps, decreasing = TRUE)
best_indices <- sorted_idx[1:min(length(sorted_idx), N_CHAINS)]

cat(sprintf("Loaded inits from %d best optimization runs.\n", length(best_indices)))

# ---- Helper Function to Unflatten Parameters ----
# Transforms single-row DF with cols like "raw_theta[1]", "raw_theta[2]" 
# into list(raw_theta = c(...))
unflatten_optim_draw <- function(draw_matrix) {
  # Convert to DF to access names easily
  d <- as.data.frame(draw_matrix)
  # Remove Stan internal cols (ending in __)
  d <- d[, !grepl("__$", names(d)), drop = FALSE]
  
  col_names <- names(d)
  # Identify unique base names (e.g., "raw_theta" from "raw_theta[1]")
  base_names <- unique(gsub("\\[.*", "", col_names))
  
  init_out <- list()
  
  for (p in base_names) {
    # Find all columns matching this parameter
    # Pattern: ^name$ OR ^name[...
    pattern <- paste0("^", p, "(\\[|$)")
    matching_cols <- grep(pattern, col_names, value = TRUE)
    
    vals <- as.numeric(d[1, matching_cols])
    
    # If it is a vector (has brackets), we MUST ensure it is sorted by index
    # (e.g. raw_theta[1], raw_theta[2], ... not raw_theta[10], raw_theta[1])
    if (length(matching_cols) > 1 && any(grepl("\\[", matching_cols))) {
      # Extract index number
      indices <- as.integer(sub(paste0("^", p, "\\[([0-9]+)\\].*"), "\\1", matching_cols))
      if (!any(is.na(indices))) {
        vals <- vals[order(indices)]
      }
    }
    
    init_out[[p]] <- vals
  }
  return(init_out)
}

# Create init list of length N_CHAINS
init_list <- list()
for (i in seq_len(N_CHAINS)) {
  # Recycle best indices if needed
  idx_to_use <- best_indices[(i - 1) %% length(best_indices) + 1]
  init_list[[i]] <- unflatten_optim_draw(valid_draws[[idx_to_use]])
}

# ---- Run NUTS ----

fit <- mod$sample(
  data = data,
  seed = 12345,
  chains = N_CHAINS,
  parallel_chains = N_CORES,       
  threads_per_chain = THREADS_PER_CHAIN, 
  init = init_list,
  iter_warmup = 500,              
  iter_sampling = 1000,            
  refresh = 20
)

# ---- Check & Save ----

if (is.null(fit) || fit$return_codes()[1] != 0) {
  stop("NUTS sampling failed.")
}

# Unpack results into list of matrices
draws_array <- fit$draws()
nuts_draws_list <- lapply(seq_len(N_CHAINS), function(i) {
  posterior::as_draws_matrix(subset_draws(draws_array, chain = i))
})

# Save Output
saveRDS(data,            file.path(outDir, "stan_data.Rds")) 
saveRDS(nuts_draws_list, file.path(outDir, "nuts_draws_list.Rds"))
saveRDS(fit$summary(),   file.path(outDir, "nuts_summary.Rds"))

cat(">>> NUTS Completed and saved.\n")


