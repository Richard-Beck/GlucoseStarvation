library(cmdstanr)
library(posterior)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

# ---- Argument Parsing ----
MODEL_NAME      <- if (length(args) >= 1) args[1] else "cobbdoug"
CONSTRAINT_NAME <- if (length(args) >= 2) args[2] else "unconstrained"
CHAIN_ID        <- if (length(args) >= 3) as.integer(args[3]) else 1L
THREADS         <- if (length(args) >= 4) as.integer(args[4]) else 32L

cat(sprintf(">>> NUTS ARRAY JOB: %s | %s | Chain=%d | Threads=%d\n",
            MODEL_NAME, CONSTRAINT_NAME, CHAIN_ID, THREADS))

# ---- Define Paths ----
base_dir <- "ecology/model_selection"
# Update output directory to include a 'nuts' subfolder
outDir <- file.path(base_dir, "data", MODEL_NAME, CONSTRAINT_NAME, "hier", "nuts")

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
}

# ---- Load Data & Constraints ----
source(file.path(base_dir, "models", MODEL_NAME, "constraints.R"))
if (!CONSTRAINT_NAME %in% names(constraints_list)) {
  stop(sprintf("Constraint '%s' not found.", CONSTRAINT_NAME))
}
active_constraint <- constraints_list[[CONSTRAINT_NAME]]

stan_data_raw <- readRDS("data/stan_ready_data.Rds")

# ---- Prepare Hierarchical Data ----
data <- stan_data_raw
grp_keys <- paste(data$line_id, data$ploidy_metric, sep = "_")
unique_grps <- unique(grp_keys)

data$N_groups <- length(unique_grps)
data$group_id <- match(grp_keys, unique_grps)
data$R2_0     <- 1.0

source(file.path(base_dir, "models", MODEL_NAME, paste0(MODEL_NAME, ".R")))
data$prior_centers <- prior_centers
data$prior_scales  <- prior_scales
data$param_map  <- active_constraint$map
data$param_mask <- active_constraint$mask
data <- data[!sapply(data,is.character)]

# ---- Compile ----
stan_file <- file.path(base_dir, "models", MODEL_NAME, paste0(MODEL_NAME, "_hier.stan"))
# We compile once; if multiple array jobs hit this simultaneously, file locking usually handles it,
# but ideally, compilation is done beforehand.
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- Prepare Single Init from Optimization ----
# Parent dir for optimization results
optimDir <- dirname(outDir) 
optim_draws_file <- file.path(optimDir, "optim_draws_all.Rds")
optim_lp_file    <- file.path(optimDir, "optim_lp_all.Rds")
optim_rc_file    <- file.path(optimDir, "optim_rc_all.Rds")

if (!file.exists(optim_draws_file)) stop("Optimization results not found.")

opt_draws_list <- readRDS(optim_draws_file)
opt_lp_vec     <- readRDS(optim_lp_file)
opt_rc_vec     <- readRDS(optim_rc_file)

valid_idx <- which((opt_rc_vec == 0L | opt_rc_vec == 70L) & is.finite(opt_lp_vec))
if (length(valid_idx) == 0) stop("No valid optimization runs found.")

# Sort by best LP
sorted_idx <- valid_idx[order(opt_lp_vec[valid_idx], decreasing = TRUE)]

# Select specific index for this chain (cycling if CHAIN_ID > available inits)
# We use modulo arithmetic to ensure we always get an init even if we run 100 chains
selected_idx <- sorted_idx[(CHAIN_ID - 1) %% length(sorted_idx) + 1]

cat(sprintf("Chain %d using optimization run index %d (LP = %f)\n", 
            CHAIN_ID, selected_idx, opt_lp_vec[selected_idx]))

# Helper to unflatten (kept compact for brevity)
unflatten_hier_draw <- function(flat_vec, d) {
  out <- list()
  get_scalar <- function(n) if(n %in% names(flat_vec)) as.numeric(flat_vec[n]) else 0
  
  out$raw_sigma_N <- get_scalar("raw_sigma_N")
  out$sigma_base  <- get_scalar("sigma_base")
  out$sigma_rel   <- get_scalar("sigma_rel")
  out$raw_N0_mean <- get_scalar("raw_N0_mean")
  
  idx <- 1:13
  out$raw_theta_ploidy <- as.numeric(flat_vec[paste0("raw_theta_ploidy[", idx, "]")])
  
  idx <- 1:d$N_groups
  out$raw_N0 <- as.numeric(flat_vec[paste0("raw_N0[", idx, "]")])
  
  idx <- 1:d$N_wells
  out$G1_0_raw <- as.numeric(flat_vec[paste0("G1_0_raw[", idx, "]")])
  
  mat <- matrix(0, nrow = 13, ncol = d$N_lines)
  for(j in 1:d$N_lines) {
    for(i in 1:13) {
      key <- paste0("raw_theta_line[", i, ",", j, "]")
      if(key %in% names(flat_vec)) mat[i,j] <- as.numeric(flat_vec[key])
    }
  }
  out$raw_theta_line <- mat
  return(out)
}

src_mat <- opt_draws_list[[selected_idx]]
flat_vec <- setNames(as.numeric(src_mat[1,]), colnames(src_mat))
init_list <- list(unflatten_hier_draw(flat_vec, data))

# ---- Run NUTS (Single Chain) ----
fit <- mod$sample(
  data = data,
  seed = 12345 + CHAIN_ID, # Unique seed per chain
  chains = 1,
  parallel_chains = 1,        
  threads_per_chain = THREADS, 
  init = init_list,
  iter_warmup = 500,               
  iter_sampling = 1000,              
  refresh = 10,
  save_warmup = FALSE,
  adapt_delta = 0.6
)

# ---- Check & Save ----
if (is.null(fit) || fit$return_codes()[1] != 0) stop("NUTS sampling failed.")

# Save Individual Chain Results
cat(sprintf("Saving results for chain %d...\n", CHAIN_ID))

saveRDS(fit$summary(), file.path(outDir, sprintf("nuts_summary_%d.Rds", CHAIN_ID)))
saveRDS(fit$draws(),   file.path(outDir, sprintf("nuts_draws_%d.Rds", CHAIN_ID)))

try({
  saveRDS(fit$sampler_diagnostics(), file.path(outDir, sprintf("nuts_diagnostics_%d.Rds", CHAIN_ID)))
})

cat(">>> Chain Completed.\n")


