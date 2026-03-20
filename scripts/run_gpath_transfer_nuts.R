#!/usr/bin/env Rscript

library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)

MODEL_NAME <- if (length(args) >= 1) args[1] else "gpath"
R_VAL <- if (length(args) >= 2) as.integer(args[2]) else 1L
P_VAL <- if (length(args) >= 3) as.integer(args[3]) else 1L
W_VAL <- if (length(args) >= 4) as.integer(args[4]) else 0L
CONSTRAINT_FLAG <- if (length(args) >= 5) as.integer(args[5]) else 0L
WASTE_MECH_FLAG <- if (length(args) >= 6) as.integer(args[6]) else 1L
LINE_ID <- if (length(args) >= 7) as.integer(args[7]) else 1L
DIRECTION <- if (length(args) >= 8) args[8] else "low_to_high"
FIT_TYPE <- if (length(args) >= 9) args[9] else "transfer"
CHAIN_ID <- if (length(args) >= 10) as.integer(args[10]) else 1L
ITER_WARMUP <- if (length(args) >= 11) as.integer(args[11]) else 500L
ITER_SAMPLING <- if (length(args) >= 12) as.integer(args[12]) else 1000L
ADAPT_DELTA <- if (length(args) >= 13) as.numeric(args[13]) else 0.99
MAX_TREED <- if (length(args) >= 14) as.integer(args[14]) else 12L
NUM_THREADS <- if (length(args) >= 15) as.integer(args[15]) else 4L
STAN_DATA_PATH <- if (length(args) >= 16) args[16] else ""
OUTPUT_ROOT <- if (length(args) >= 17) args[17] else "data/gpath_transfer_cv_nuts"
PLOIDY_EFFECT_MASK_SPEC <- if (length(args) >= 18) args[18] else ""
INIT_MODE <- if (length(args) >= 19) args[19] else "auto"
INIT_SOURCE <- if (length(args) >= 20) args[20] else ""

run_id <- sprintf("%dR_%dP_%dW_C%d_M%d", R_VAL, P_VAL, W_VAL, CONSTRAINT_FLAG, WASTE_MECH_FLAG)

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/posterior_io_utils.R")
source("R/elpd_transfer_utils.R")

DIRECTION <- normalize_transfer_direction(DIRECTION)
FIT_TYPE <- normalize_fit_type(FIT_TYPE)
STAN_DATA_PATH <- resolve_stan_data_path(STAN_DATA_PATH)
INIT_MODE <- tolower(trimws(INIT_MODE))

if (!(INIT_MODE %in% c("auto", "posterior", "optim", "random"))) {
  stop(sprintf("Unsupported init mode '%s'", INIT_MODE))
}

cat(sprintf(">>> run_gpath_transfer_nuts.R cwd: %s\n", getwd()))
cat(sprintf(">>> Transfer CV NUTS | %s | %s | line=%d | direction=%s | fit=%s | chain=%d\n",
            MODEL_NAME, run_id, LINE_ID, DIRECTION, FIT_TYPE, CHAIN_ID))
cat(sprintf(">>> Stan data: %s\n", STAN_DATA_PATH))
cat(sprintf(">>> Init mode: %s\n", INIT_MODE))
if (nzchar(trimws(INIT_SOURCE))) {
  cat(sprintf(">>> Init source: %s\n", INIT_SOURCE))
}

stan_data <- prepare_gpath_stan_data(
  stan_data_path = STAN_DATA_PATH,
  R_val = R_VAL,
  P_val = P_VAL,
  W_val = W_VAL,
  constraint_flag = CONSTRAINT_FLAG,
  waste_mech_flag = WASTE_MECH_FLAG,
  base_priors = base_priors,
  ploidy_effect_mask_spec = PLOIDY_EFFECT_MASK_SPEC
)
mask_info <- attr(stan_data, "ploidy_effect_mask_info")

split_meta <- get_directional_transfer_split(stan_data = stan_data, line_id = LINE_ID, direction = DIRECTION)

if (FIT_TYPE %in% c("transfer", "null")) {
  stan_data <- apply_directional_transfer_holdout(stan_data = stan_data, line_id = LINE_ID, direction = DIRECTION)
} else {
  stan_data$is_train <- rep(1L, stan_data$N_wells)
}

if (FIT_TYPE == "null") {
  stan_data$ploidy_metric[] <- 0.0
}

split_meta$fit_type <- FIT_TYPE
split_meta$model_name <- MODEL_NAME
split_meta$run_id <- run_id
split_meta$stan_data_path <- STAN_DATA_PATH
split_meta$ploidy_effect_mask_spec <- if (nzchar(trimws(PLOIDY_EFFECT_MASK_SPEC))) PLOIDY_EFFECT_MASK_SPEC else "all"
split_meta$ploidy_effect_mask_label <- if (!is.null(mask_info$label)) mask_info$label else "all"
split_meta$ploidy_effect_mask <- stan_data$ploidy_effect_mask

stan_file <- get_model_stan_path(MODEL_NAME, "v1")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

posterior_dir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, "hier", "nuts")
seed_init <- 3000000L + 1000L * LINE_ID + CHAIN_ID
seed_fit <- 4000000L + 1000L * LINE_ID + CHAIN_ID
maxG0 <- max(as.numeric(stan_data$G0_per_well), na.rm = TRUE)
cat(sprintf(">>> NUTS config: warmup=%d sampling=%d adapt_delta=%.3f max_treedepth=%d seed=%d\n",
            ITER_WARMUP, ITER_SAMPLING, ADAPT_DELTA, MAX_TREED, seed_fit))

init_arg <- 2
if (INIT_MODE == "random") {
  cat(">>> Using randomized CmdStan init=2\n")
} else if (INIT_MODE == "optim") {
  if (!nzchar(trimws(INIT_SOURCE))) {
    stop("Optimization init requested but init source directory was not provided")
  }
  if (!dir.exists(INIT_SOURCE)) {
    stop(sprintf("Optimization init requested but directory is missing: %s", INIT_SOURCE))
  }

  init_list <- load_optim_init_from_dir(
    path = INIT_SOURCE,
    chain_id = CHAIN_ID,
    maxG0 = maxG0
  )
  init_arg <- list(init_list)
  cat(sprintf(">>> Init from optimization outputs in %s\n", INIT_SOURCE))
} else if (INIT_MODE == "posterior") {
  posterior_source <- if (nzchar(trimws(INIT_SOURCE))) INIT_SOURCE else posterior_dir
  if (!dir.exists(posterior_source)) {
    stop(sprintf("Posterior init requested but directory is missing: %s", posterior_source))
  }

  init_list <- load_posterior_init_from_dir(
    path = posterior_source,
    seed = seed_init,
    maxG0 = maxG0
  )
  init_arg <- list(init_list)
  cat(sprintf(">>> Init from posterior draws in %s\n", posterior_source))
} else if (dir.exists(posterior_dir)) {
  init_list <- tryCatch(
    load_posterior_init_from_dir(path = posterior_dir, seed = seed_init, maxG0 = maxG0),
    error = function(e) NULL
  )
  if (!is.null(init_list)) {
    init_arg <- list(init_list)
    cat(sprintf(">>> Init from posterior draws in %s\n", posterior_dir))
  } else if (nzchar(trimws(INIT_SOURCE)) && dir.exists(INIT_SOURCE)) {
    init_list <- tryCatch(
      load_optim_init_from_dir(path = INIT_SOURCE, chain_id = CHAIN_ID, maxG0 = maxG0),
      error = function(e) NULL
    )
    if (!is.null(init_list)) {
      init_arg <- list(init_list)
      cat(sprintf(">>> Posterior init unavailable; fell back to optimization outputs in %s\n", INIT_SOURCE))
    } else {
      cat(">>> Posterior and optimization init unavailable; using init=2\n")
    }
  } else {
    cat(">>> Posterior init unavailable; using init=2\n")
  }
} else if (nzchar(trimws(INIT_SOURCE)) && dir.exists(INIT_SOURCE)) {
  init_list <- tryCatch(
    load_optim_init_from_dir(path = INIT_SOURCE, chain_id = CHAIN_ID, maxG0 = maxG0),
    error = function(e) NULL
  )
  if (!is.null(init_list)) {
    init_arg <- list(init_list)
    cat(sprintf(">>> Init from optimization outputs in %s\n", INIT_SOURCE))
  } else {
    cat(">>> Optimization init unavailable; using init=2\n")
  }
} else {
  cat(">>> No posterior/optimization init source available; using init=2\n")
}

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

run_tag <- sprintf(
  "line%02d_%s_%s_chain%02d",
  LINE_ID,
  DIRECTION,
  FIT_TYPE,
  CHAIN_ID
)
output_dir <- file.path(OUTPUT_ROOT, MODEL_NAME, run_id, DIRECTION, FIT_TYPE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(output_dir, sprintf("nuts_summary_%s.Rds", run_tag))
draws_path <- file.path(output_dir, sprintf("nuts_draws_%s.Rds", run_tag))
diag_path <- file.path(output_dir, sprintf("nuts_diagnostics_%s.Rds", run_tag))
meta_path <- file.path(output_dir, sprintf("split_meta_%s.Rds", run_tag))
chain_summary_path <- file.path(output_dir, sprintf("chain_summary_%s.Rds", run_tag))

saveRDS(fit$summary(), summary_path)
saveRDS(fit$draws(), draws_path)
diag <- tryCatch(fit$sampler_diagnostics(), error = function(e) NULL)
if (!is.null(diag)) saveRDS(diag, diag_path)
saveRDS(split_meta, meta_path)

chain_summary <- summarize_transfer_draws(fit$draws(), split_meta)
chain_summary$run_tag <- run_tag
chain_summary$chain_id <- CHAIN_ID
saveRDS(chain_summary, chain_summary_path)

cat(sprintf(">>> Saved NUTS draws to %s\n", draws_path))
cat(sprintf(">>> Chain ELPD on held-out wells: %.3f\n", chain_summary$elpd))
