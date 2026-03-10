library(cmdstanr)
library(posterior)
library(jsonlite)

# ---- 1. Argument Parsing ----
args <- commandArgs(trailingOnly = TRUE)

MODEL_NAME      <- if (length(args) >= 1) args[1] else "gpath"
R_val           <- if (length(args) >= 2) as.integer(args[2]) else 2L
P_val           <- if (length(args) >= 3) as.integer(args[3]) else 2L
W_val           <- if (length(args) >= 4) as.integer(args[4]) else 0L
CONSTRAINT_FLAG <- if (length(args) >= 5) as.integer(args[5]) else 0L
WASTE_MECH_FLAG <- if (length(args) >= 6) as.integer(args[6]) else 1L
NUM_THREADS     <- if (length(args) >= 7) as.integer(args[7]) else 4L

# ---- Pathfinder controls ----
SEED             <- if (length(args) >= 8)  as.integer(args[8]) else 1L
PF_DRAWS         <- if (length(args) >= 9)  as.integer(args[9]) else 1000L
PF_SINGLE_DRAWS  <- if (length(args) >= 10) as.integer(args[10]) else 1000L
PF_MAX_LBFGS     <- if (length(args) >= 11) as.integer(args[11]) else 1000L
PF_INIT_ALPHA    <- if (length(args) >= 12) as.numeric(args[12]) else 1e-2
PF_NUM_PATHS     <- if (length(args) >= 13) as.integer(args[13]) else 4L

HOLDOUT <- FALSE
tmp <- strsplit(MODEL_NAME, split = "_")[[1]]
if (length(tmp) > 1) {
  if (tmp[[2]] == "holdout") {
    HOLDOUT <- TRUE
    MODEL_NAME <- tmp[[1]]
  }
}

run_id <- sprintf("%dR_%dP_%dW_C%d_M%d", R_val, P_val, W_val, CONSTRAINT_FLAG, WASTE_MECH_FLAG)

source("ecology/model_selection/models/gpath/gpath.R")

cat(sprintf(">>> RUNNING PATHFINDER | %s | %s | num_paths=%d | draws=%d\n",
            MODEL_NAME, run_id, PF_NUM_PATHS, PF_DRAWS))

# ---- 2. Load Base Data ----
stan_data <- readRDS("data/stan_ready_data.Rds")

grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
unique_grps <- unique(grp_keys)
stan_data$N_groups <- length(unique_grps)
stan_data$group_id <- match(grp_keys, unique_grps)

# ---- 3. Apply Configuration ----
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

R_init_base <- matrix(0.0, nrow = stan_data$N_wells, ncol = R_val)
if (R_val > 1) for (c in 2:R_val) R_init_base[, c] <- 1.0
stan_data$R_init_base <- R_init_base
stan_data$is_train <- rep(1L, stan_data$N_wells)
if (HOLDOUT) stan_data$is_train[stan_data$ploidy_metric > 0] <- 0
stan_data <- stan_data[!sapply(stan_data, is.character)]

# ---- 4. Compile Model ----
stan_file <- file.path("ecology/model_selection/models", MODEL_NAME, paste0(MODEL_NAME, "_hier.stan"))
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- 5. Output directory ----
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, ifelse(HOLDOUT, "holdout", "hier"))
pfDir <- file.path(outDir, "pathfinder")
dir.create(pfDir, recursive = TRUE, showWarnings = FALSE)

outfile_base <- file.path(pfDir, "pf")

# ---- 5b. Helper: convert flat 'name[index,...]' vector to nested init list ----
PARAM_NAMES <- c(
  "raw_theta_line",
  "raw_theta_ploidy",
  "raw_N0",
  "raw_sigma_N",
  "raw_sigma_base",
  "raw_sigma_rel",
  "G1_0"
)

clamp_constrained <- function(init) {
  if (!is.null(init$G1_0)) init$G1_0 <- pmax(init$G1_0, 1e-12)
  init
}

flat_to_init_list <- function(x_named, param_names = PARAM_NAMES) {
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
      idx_list <- lapply(these, function(s) as.integer(strsplit(gsub("^.*\\[|\\]$", "", s), ",")[[1]]))
      dims <- pmax(1L, Reduce(pmax, idx_list))
      arr <- array(0.0, dim = dims)
      
      for (k in seq_along(these)) {
        ii <- idx_list[[k]]
        arr[matrix(ii, nrow = 1)] <- unname(x_named[[these[[k]]]])
      }
      
      if (length(dims) == 1L) out[[bn]] <- as.numeric(arr)
      else if (length(dims) == 2L) out[[bn]] <- matrix(arr, nrow = dims[1], ncol = dims[2])
      else out[[bn]] <- arr
    }
  }
  
  clamp_constrained(out)
}

# ---- 5c. Init from best optimization fit ----
init_arg <- 2
opt_lp_path <- file.path(outDir, "optim_lp_all.Rds")
opt_dr_path <- file.path(outDir, "optim_draws_all.Rds")

if (file.exists(opt_lp_path) && file.exists(opt_dr_path)) {
  lp_vec <- readRDS(opt_lp_path)
  dr_lst <- readRDS(opt_dr_path)
  
  ok <- which(is.finite(lp_vec))
  if (length(ok) >= 1) {
    ord <- ok[order(lp_vec[ok], decreasing = TRUE)]
    pick <- ord[[1]]
    
    dm <- dr_lst[[pick]]
    if (!is.null(dm)) {
      row <- if (is.matrix(dm)) dm[1, , drop = TRUE] else dm
      if (is.matrix(dm)) names(row) <- colnames(dm)
      init_list <- tryCatch(flat_to_init_list(row), error = function(e) NULL)
      
      if (!is.null(init_list)) {
        init_arg <- list(init_list)   # one sublist, reused for pathfinder
        cat(sprintf(">>> Init from optimization: rank=1 (start=%d), lp=%g\n", pick, lp_vec[[pick]]))
      } else {
        cat(">>> Init conversion failed; using init=2\n")
      }
    } else {
      cat(">>> Picked optimization draw is NULL; using init=2\n")
    }
  } else {
    cat(">>> No finite optimization fits; using init=2\n")
  }
} else {
  cat(">>> No optimization init files found; using init=2\n")
}

# ---- 6. Run multi-path Pathfinder from same init ----
fit_pf <- mod$pathfinder(
  data = stan_data,
  seed = SEED,
  init = 2,
  num_paths = PF_NUM_PATHS,
  draws = PF_DRAWS,
  single_path_draws = PF_SINGLE_DRAWS,
  max_lbfgs_iters = PF_MAX_LBFGS,
  init_alpha = PF_INIT_ALPHA,
  num_threads = NUM_THREADS,
  save_single_paths = TRUE,
  psis_resample = TRUE,
  calculate_lp = TRUE,
  refresh = 10
#  output_dir = pfDir,
#  output_basename = "pf"
)

cat(">>> Pathfinder completed.\n")

saveRDS(fit_pf, file.path(pfDir, "pathfinder_fit.Rds"))
saveRDS(fit_pf$draws(), file.path(pfDir, "pathfinder_draws.Rds"))

# optional: if available in your cmdstanr version
try(saveRDS(fit_pf$summary(), file.path(pfDir, "pathfinder_summary.Rds")), silent = TRUE)
