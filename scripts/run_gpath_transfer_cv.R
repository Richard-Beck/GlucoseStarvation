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
START_ID <- if (length(args) >= 10) as.integer(args[10]) else 1L
ITER_WARMUP <- if (length(args) >= 11) as.integer(args[11]) else 500L
ITER_SAMPLING <- if (length(args) >= 12) as.integer(args[12]) else 1000L
ADAPT_DELTA <- if (length(args) >= 13) as.numeric(args[13]) else 0.9
MAX_TREED <- if (length(args) >= 14) as.integer(args[14]) else 12L
NUM_THREADS <- if (length(args) >= 15) as.integer(args[15]) else 4L
STAN_DATA_PATH <- if (length(args) >= 16) args[16] else ""
OUTPUT_ROOT <- if (length(args) >= 17) args[17] else "data/gpath_transfer_cv"
INIT_MODE <- if (length(args) >= 18) args[18] else "auto"
PLOIDY_EFFECT_MASK_SPEC <- if (length(args) >= 19) args[19] else ""

run_id <- sprintf(
  "%dR_%dP_%dW_C%d_M%d",
  R_VAL,
  P_VAL,
  W_VAL,
  CONSTRAINT_FLAG,
  WASTE_MECH_FLAG
)

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/posterior_io_utils.R")
source("R/elpd_transfer_utils.R")

DIRECTION <- normalize_transfer_direction(DIRECTION)
FIT_TYPE <- normalize_fit_type(FIT_TYPE)
STAN_DATA_PATH <- resolve_stan_data_path(STAN_DATA_PATH)
INIT_MODE <- tolower(trimws(INIT_MODE))

if (!(INIT_MODE %in% c("auto", "posterior", "random"))) {
  stop(sprintf("Unsupported init mode '%s'", INIT_MODE))
}

cat(sprintf(">>> run_gpath_transfer_cv.R cwd: %s\n", getwd()))

cat(sprintf(
  ">>> Transfer CV optimize | %s | %s | line=%d | direction=%s | fit=%s | start=%d\n",
  MODEL_NAME,
  run_id,
  LINE_ID,
  DIRECTION,
  FIT_TYPE,
  START_ID
))
cat(sprintf(">>> Stan data: %s\n", STAN_DATA_PATH))
cat(sprintf(">>> Init mode: %s\n", INIT_MODE))

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

split_meta <- get_directional_transfer_split(
  stan_data = stan_data,
  line_id = LINE_ID,
  direction = DIRECTION
)

if (FIT_TYPE %in% c("transfer", "null")) {
  stan_data <- apply_directional_transfer_holdout(
    stan_data = stan_data,
    line_id = LINE_ID,
    direction = DIRECTION
  )
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
if (!file.exists(stan_file)) {
  stop(sprintf("Stan file not found: %s", stan_file))
}

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

posterior_dir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, "hier", "nuts")
seed_init <- 1000000L + 1000L * LINE_ID + START_ID
seed_fit <- 2000000L + 1000L * LINE_ID + START_ID
maxG0 <- max(as.numeric(stan_data$G0_per_well), na.rm = TRUE)
cat(sprintf(">>> Optimization config: algorithm=lbfgs iter=%d seed=%d\n", ITER_SAMPLING, seed_fit))

init_arg <- 2
if (INIT_MODE == "posterior" && !dir.exists(posterior_dir)) {
  stop(sprintf("Posterior init requested but directory is missing: %s", posterior_dir))
}

if (INIT_MODE == "random") {
  cat(">>> Using randomized CmdStan init=2\n")
} else if (dir.exists(posterior_dir)) {
  init_list <- tryCatch(
    load_posterior_init_from_dir(
      path = posterior_dir,
      seed = seed_init,
      maxG0 = maxG0
    ),
    error = function(e) NULL
  )

  if (!is.null(init_list)) {
    init_arg <- list(init_list)
    cat(sprintf(">>> Init from posterior draws in %s\n", posterior_dir))
  } else {
    if (INIT_MODE == "posterior") {
      stop("Posterior init requested but no usable posterior init could be loaded")
    }
    cat(">>> Posterior init unavailable; using init=2\n")
  }
} else {
  cat(sprintf(">>> Posterior directory missing (%s); using init=2\n", posterior_dir))
}

opt_res <- tryCatch({
  opt <- mod$optimize(
    data = stan_data,
    algorithm = "lbfgs",
    init = init_arg,
    refresh = 10,
    threads = NUM_THREADS,
    save_latent_dynamics = TRUE,
    jacobian = TRUE,
    show_messages = TRUE,
    seed = seed_fit,
    iter = ITER_SAMPLING
  )

  rc <- opt$return_codes()[1]

  if (rc == 0L) {
    list(
      rc = rc,
      lp = opt$lp(),
      draws = opt$draws(format = "matrix"),
      summary = list(
        return_code = rc,
        lp = opt$lp(),
        seed = seed_fit,
        algorithm = "lbfgs",
        iter = ITER_SAMPLING
      ),
      diagnostics = tryCatch(opt$output(), error = function(e) NULL)
    )
  } else {
    txt <- capture.output(opt$output())
    csv_line <- grep("^\\s*file = ", txt, value = TRUE)
    csv_line <- tail(csv_line, 1)
    csv <- sub("^\\s*file =\\s*", "", csv_line)

    if (!nzchar(csv) || !file.exists(csv)) {
      stop(sprintf("Optimization failed with rc=%s and no readable CSV fallback", rc))
    }

    csv_fit <- read_cmdstan_csv(csv)
    point_est <- csv_fit$point_estimates
    lp_val <- if ("lp__" %in% colnames(point_est)) point_est[, "lp__"] else NA_real_

    list(
      rc = rc,
      lp = lp_val,
      draws = point_est,
      summary = list(
        return_code = rc,
        lp = lp_val,
        seed = seed_fit,
        algorithm = "lbfgs",
        iter = ITER_SAMPLING,
        csv_fallback = csv
      ),
      diagnostics = txt
    )
  }
}, error = function(e) {
  stop(sprintf("Transfer optimization failed for line=%d direction=%s fit=%s start=%d: %s",
               LINE_ID, DIRECTION, FIT_TYPE, START_ID, conditionMessage(e)))
})

run_tag <- build_transfer_run_id(
  line_id = LINE_ID,
  direction = DIRECTION,
  fit_type = FIT_TYPE,
  start_id = START_ID
)
output_dir <- build_transfer_output_dir(
  output_root = OUTPUT_ROOT,
  model_name = MODEL_NAME,
  run_id = run_id,
  direction = DIRECTION,
  fit_type = FIT_TYPE
)

paths <- save_transfer_run_outputs(
  draws = opt_res$draws,
  output_dir = output_dir,
  run_tag = run_tag,
  split_meta = split_meta,
  summary_obj = opt_res$summary,
  diag_obj = opt_res$diagnostics
)

start_summary <- summarize_transfer_draws(opt_res$draws, split_meta)
start_summary$run_tag <- run_tag
start_summary$start_id <- START_ID
start_summary$return_code <- opt_res$rc
saveRDS(
  start_summary,
  file.path(output_dir, sprintf("start_summary_%s.Rds", run_tag))
)

cat(sprintf(">>> Saved optimize draws to %s\n", paths$draws_path))
cat(sprintf(">>> Optimize return code: %s\n", opt_res$rc))
cat(sprintf(">>> Start ELPD on held-out wells: %.3f\n", start_summary$elpd))
