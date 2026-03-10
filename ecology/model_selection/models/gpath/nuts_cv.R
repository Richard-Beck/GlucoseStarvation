#!/usr/bin/env Rscript

library(cmdstanr)
library(posterior)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

# ---- 1) Args ----
MODEL_NAME      <- if (length(args) >= 1)  args[1] else "gpath"
R_val           <- if (length(args) >= 2)  as.integer(args[2]) else 2L
P_val           <- if (length(args) >= 3)  as.integer(args[3]) else 2L
W_val           <- if (length(args) >= 4)  as.integer(args[4]) else 0L
CONSTRAINT_FLAG <- if (length(args) >= 5)  as.integer(args[5]) else 0L
WASTE_MECH_FLAG <- if (length(args) >= 6)  as.integer(args[6]) else 1L
NUM_THREADS     <- if (length(args) >= 7)  as.integer(args[7]) else 4L

FOLD_ID         <- if (length(args) >= 8)  as.integer(args[8]) else 1L  # 1..N_lines
CHAIN_ID        <- if (length(args) >= 9)  as.integer(args[9]) else 1L

ITER_WARMUP     <- if (length(args) >= 10) as.integer(args[10]) else 500L
ITER_SAMPLING   <- if (length(args) >= 11) as.integer(args[11]) else 1000L
ADAPT_DELTA     <- if (length(args) >= 12) as.numeric(args[12]) else 0.9
MAX_TREED       <- if (length(args) >= 13) as.integer(args[13]) else 12L

run_id <- sprintf("%dR_%dP_%dW_C%d_M%d", R_val, P_val, W_val, CONSTRAINT_FLAG, WASTE_MECH_FLAG)

cat(sprintf(
  ">>> CV NUTS | %s | %s | fold(line)=%d | chain=%d | warmup=%d | sampling=%d\n",
  MODEL_NAME, run_id, FOLD_ID, CHAIN_ID, ITER_WARMUP, ITER_SAMPLING
))

# ---- 2) Source shared utilities (expects generate_stan_config, base_priors, etc.) ----
source("ecology/model_selection/models/gpath/gpath.R")

if (!exists("generate_stan_config")) stop("generate_stan_config() not found after sourcing gpath.R")

# If base_priors is not defined in gpath.R, define it here
if (!exists("base_priors")) {
  base_priors <- list(
    ae_c = 6.0e-5, ae_s = 0.5,
    ah_c = 1.0,    ah_s = 0.5,
    Y_R_c = 2500.0, Y_R_s = 0.5,
    A_c  = 1.0,    A_s  = 0.5,
    A_gap_c = 1.0, A_gap_s = 0.5,
    K_add_c = 0.01, K_add_s = 1.5,
    K_mult_c = .01, K_mult_s = 1.5,
    m_c  = 0.05,   m_s  = 0.2
  )
}

# ---- 3) Load base data ----
stan_data <- readRDS("data/stan_ready_data.Rds")

# Rebuild group_id exactly as in your existing pipeline
grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
unique_grps <- unique(grp_keys)
stan_data$N_groups <- length(unique_grps)
stan_data$group_id <- match(grp_keys, unique_grps)

# ---- 4) Apply configuration ----
config <- generate_stan_config(
  R = R_val,
  P = P_val,
  W = W_val,
  strict_spec = (CONSTRAINT_FLAG == 1L),
  M = WASTE_MECH_FLAG,
  base_priors = base_priors
)

for (nm in names(config)) stan_data[[nm]] <- config[[nm]]

if (W_val > 0) stan_data$waste_mech <- rep(as.numeric(WASTE_MECH_FLAG), W_val) else stan_data$waste_mech <- numeric(0)

# Build R_init_base consistent with your current script
R_init_base <- matrix(0.0, nrow = stan_data$N_wells, ncol = R_val)
if (R_val > 1) for (c in 2:R_val) R_init_base[, c] <- 1.0
stan_data$R_init_base <- R_init_base

# Drop accidental characters
stan_data <- stan_data[!sapply(stan_data, is.character)]

# ---- 5) Create CV mask is_train (hold out high ploidy wells for this line) ----
if (!("N_lines" %in% names(stan_data))) stop("stan_data$N_lines missing")
if (FOLD_ID < 1L || FOLD_ID > stan_data$N_lines) stop("FOLD_ID out of range")

line_vec   <- as.integer(stan_data$line_id)
ploidy_vec <- as.numeric(stan_data$ploidy_metric)

idx_line <- which(line_vec == FOLD_ID)
if (!length(idx_line)) stop(sprintf("No wells found for line_id == %d", FOLD_ID))

# Define "high ploidy" as max(ploidy_metric) within this line (robust to 0/1 or -0.5/+0.5)
hi_ploidy <- max(ploidy_vec[idx_line], na.rm = TRUE)
tol <- 1e-12
idx_holdout <- which(line_vec == FOLD_ID & abs(ploidy_vec - hi_ploidy) < tol)

if (!length(idx_holdout)) {
  stop(sprintf("No heldout wells found for line=%d with ploidy_metric == max(ploidy_metric) == %g", FOLD_ID, hi_ploidy))
}

is_train <- rep(1L, stan_data$N_wells)
is_train[idx_holdout] <- 0L
stan_data$is_train <- as.integer(is_train)
#stan_data$ploidy_metric <- log2(stan_data$ploidy_abs) - 1
cat(sprintf(">>> CV split: holding out %d wells (line=%d, ploidy_metric=%g)\n", length(idx_holdout), FOLD_ID, hi_ploidy))

# ---- 6) Compile CV Stan model ----
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, "hier")

# Prefer a dedicated CV Stan file
stan_file_cv <- file.path("ecology/model_selection/models", MODEL_NAME, paste0(MODEL_NAME, "_hier.stan"))
if (!file.exists(stan_file_cv)) stop(sprintf("CV Stan file not found: %s", stan_file_cv))

mod <- cmdstan_model(stan_file_cv, cpp_options = list(stan_threads = TRUE))

# ---- 7) Output paths ----
nutsDir <- file.path(outDir, "nuts_cv")
dir.create(nutsDir, recursive = TRUE, showWarnings = FALSE)

tag <- sprintf("fold%02d_chain%02d", FOLD_ID, CHAIN_ID)

summary_path <- file.path(nutsDir, sprintf("nuts_summary_%s.Rds", tag))
draws_path   <- file.path(nutsDir, sprintf("nuts_draws_%s.Rds", tag))
diag_path    <- file.path(nutsDir, sprintf("nuts_diagnostics_%s.Rds", tag))

# Skip completed
#if (file.exists(draws_path)) {
 # cat(">>> Output exists; skipping.\n")
  #quit(save = "no", status = 0)
#}

# ---- 8) Build init from FULL posterior draws (same run_id), seeded by fold+chain ----
PARAM_NAMES <- c("raw_theta_line","raw_theta_ploidy","raw_N0","raw_sigma_N","sigma_base","sigma_rel","G1_0")

clamp_constrained <- function(init, maxG0) {
  if (!is.null(init$sigma_base)) init$sigma_base <- max(init$sigma_base, 1e-12)
  if (!is.null(init$sigma_rel))  init$sigma_rel  <- max(init$sigma_rel,  1e-12)
  if (!is.null(init$G1_0)) {
    init$G1_0 <- pmax(init$G1_0, 1e-12)
    init$G1_0 <- pmin(init$G1_0, (2.0 * maxG0) - 1e-12)
  }
  init
}

flat_to_init_list <- function(x_named, param_names = PARAM_NAMES, maxG0) {
  keep <- !grepl("__$", names(x_named)) & names(x_named) != "lp__"
  x_named <- x_named[keep]
  x_named <- x_named[is.finite(x_named)]
  
  nms <- names(x_named)
  if (is.null(nms) || !length(nms)) stop("init vector has no usable names")
  
  out <- list()
  
  scal <- !grepl("\\[", nms)
  if (any(scal)) {
    for (nm in nms[scal]) if (nm %in% param_names) out[[nm]] <- unname(x_named[[nm]])
  }
  
  idx_nms <- nms[!scal]
  if (length(idx_nms)) {
    base <- sub("\\[.*$", "", idx_nms)
    
    for (bn in unique(base)) {
      if (!(bn %in% param_names)) next
      
      these <- idx_nms[base == bn]
      
      idx_list <- lapply(these, function(s) {
        # strip "...[", trailing "]", then split by comma
        as.integer(strsplit(gsub("^.*\\[|\\]$", "", s), ",")[[1]])
      })
      
      dims <- pmax(1L, Reduce(pmax, idx_list))
      arr <- array(0.0, dim = dims)
      
      for (k in seq_along(these)) {
        ii <- idx_list[[k]]
        arr[matrix(ii, nrow = 1)] <- unname(x_named[[these[[k]]]])
      }
      
      if (length(dims) == 1L)      out[[bn]] <- as.numeric(arr)
      else if (length(dims) == 2L) out[[bn]] <- matrix(arr, nrow = dims[1], ncol = dims[2])
      else                         out[[bn]] <- arr
    }
  }
  
  clamp_constrained(out, maxG0 = maxG0)
}

# Load full posterior draws from the *non-CV* directory
nutsFullDir <- file.path(outDir, "nuts")
full_draw_files <- list.files(nutsFullDir, pattern = "^nuts_draws_[0-9]+\\.Rds$", full.names = TRUE)

init_arg <- 2  # fallback

if (length(full_draw_files) > 0) {
  cat(sprintf(">>> Loading %d full-fit draw files from %s\n", length(full_draw_files), nutsFullDir))
  
  dr_list <- lapply(full_draw_files, readRDS)
  # bind along chain dimension (works for draws_array; otherwise convert)
  dr <- tryCatch({
    Reduce(function(a,b) posterior::bind_draws(a, b, along = "chain"), dr_list)
  }, error = function(e) NULL)
  
  if (!is.null(dr)) {
    dm <- posterior::as_draws_matrix(dr)
    if (nrow(dm) < 1) stop("Full-fit draws matrix is empty")
    
    seed_init <- 1000000L + 1000L * FOLD_ID + CHAIN_ID
    set.seed(seed_init)
    
    pick <- sample.int(nrow(dm), 1)
    x_named <- dm[pick, , drop = TRUE]               # now a named numeric vector (names = colnames)
    x_named <- x_named[is.finite(x_named)]           # keeps names
    x_named <- x_named[!is.na(x_named)]
    
    maxG0 <- max(as.numeric(stan_data$G0_per_well), na.rm = TRUE)
    init_list <- tryCatch(flat_to_init_list(x_named, maxG0 = maxG0), error = function(e) NULL)
    
    if (!is.null(init_list)) {
      init_arg <- list(init_list)  # one chain
      cat(sprintf(">>> Init from full posterior: picked draw row=%d (seed_init=%d)\n", pick, seed_init))
    } else {
      cat(">>> Init conversion failed; using init=2\n")
    }
  } else {
    cat(">>> Failed to bind full-fit draws; using init=2\n")
  }
} else {
  cat(sprintf(">>> No full-fit draw files found in %s; using init=2\n", nutsFullDir))
}

# ---- 9) Run NUTS (single chain) ----
seed_fit <- 2000000L + 1000L * FOLD_ID + CHAIN_ID

fit <- mod$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1,
  threads_per_chain = NUM_THREADS,
  seed = seed_fit,
  iter_warmup = ITER_WARMUP,
  iter_sampling = ITER_SAMPLING,
  adapt_delta = ADAPT_DELTA,
  max_treedepth = MAX_TREED,
  save_warmup = TRUE,
  init = init_arg,
  metric = "dense_e",
  refresh = 10
)

cat(sprintf(">>> Saving CV results: %s\n", tag))

saveRDS(fit$summary(), summary_path)
saveRDS(fit$draws(),   draws_path)

# sampler diagnostics are optional
diag <- tryCatch(fit$sampler_diagnostics(), error = function(e) NULL)
if (!is.null(diag)) saveRDS(diag, diag_path)

cat(">>> CV chain completed.\n")
