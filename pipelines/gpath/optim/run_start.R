#!/usr/bin/env Rscript

library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)

MODEL_NAME <- if (length(args) >= 1) args[1] else "gpath"
MODEL_VERSION <- if (length(args) >= 2) args[2] else "v1"
RUN_ID <- if (length(args) >= 3) args[3] else stop("run_id required")
NUM_THREADS <- if (length(args) >= 4) as.integer(args[4]) else 4L
TASK_ID <- if (length(args) >= 5) as.integer(args[5]) else 1L
STAN_DATA_PATH <- if (length(args) >= 6) args[6] else ""
DATASET_LABEL <- if (length(args) >= 7) args[7] else "gstarvation_v1"
OUTPUT_DIR <- if (length(args) >= 8) args[8] else stop("output_dir required")
INIT_RADIUS <- suppressWarnings(as.numeric(Sys.getenv("GPATH_OPTIM_INIT_RADIUS", "2")))
if (length(INIT_RADIUS) != 1L || !is.finite(INIT_RADIUS) || INIT_RADIUS <= 0) {
  stop("GPATH_OPTIM_INIT_RADIUS must be a finite number greater than zero")
}
INIT_SOURCE_ROOT <- Sys.getenv("GPATH_OPTIM_INIT_SOURCE_ROOT", "")
INIT_OFFSET <- suppressWarnings(as.integer(Sys.getenv("GPATH_OPTIM_INIT_OFFSET", "0")))
if (length(INIT_OFFSET) != 1L || is.na(INIT_OFFSET) || INIT_OFFSET < 0L) {
  stop("GPATH_OPTIM_INIT_OFFSET must be a non-negative integer")
}

source("R/project_paths.R")
source(get_model_r_path(MODEL_NAME, MODEL_VERSION))
source("R/gpath_run_utils.R")
source("R/run_manifest_utils.R")

run_parts <- parse_run_id(RUN_ID)
R_VAL <- run_parts$R
P_VAL <- run_parts$P
W_VAL <- run_parts$W
CONSTRAINT_FLAG <- run_parts$C
WASTE_MECH_FLAG <- run_parts$M
stan_data_path <- resolve_stan_data_path(STAN_DATA_PATH)
stan_file <- get_model_stan_path(MODEL_NAME, MODEL_VERSION)

cat(sprintf(">>> run_start.R cwd: %s\n", getwd()))
cat(sprintf(">>> OPTIM START | model=%s version=%s run_id=%s task=%d\n", MODEL_NAME, MODEL_VERSION, RUN_ID, TASK_ID))
cat(sprintf(">>> Stan data: %s\n", stan_data_path))
cat(sprintf(">>> Stan file: %s\n", stan_file))
cat(sprintf(">>> Output dir: %s\n", OUTPUT_DIR))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

manifest <- build_run_manifest(
  model_name = MODEL_NAME,
  model_version = MODEL_VERSION,
  pipeline_name = "optim",
  run_id = RUN_ID,
  dataset_label = DATASET_LABEL,
  stan_data_path = stan_data_path,
  stan_file_path = stan_file,
  script_path = "pipelines/gpath/optim/run_start.R",
  command_args = list(
    R = R_VAL,
    P = P_VAL,
    W = W_VAL,
    constraint_flag = CONSTRAINT_FLAG,
    waste_mech_flag = WASTE_MECH_FLAG,
    num_threads = NUM_THREADS,
    task_id = TASK_ID,
    init_radius = INIT_RADIUS
  ),
  output_dir = OUTPUT_DIR
)
write_run_manifest(OUTPUT_DIR, manifest)

stan_data <- readRDS(stan_data_path)
stan_data <- add_group_structure(stan_data)
if (is.null(stan_data$is_train)) {
  stan_data$is_train <- rep(1L, stan_data$N_wells)
}
stan_data <- apply_gpath_run_config(
  stan_data = stan_data,
  R_val = R_VAL,
  P_val = P_VAL,
  W_val = W_VAL,
  constraint_flag = CONSTRAINT_FLAG,
  waste_mech_flag = WASTE_MECH_FLAG,
  base_priors = base_priors
)

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

init_spec <- INIT_RADIUS
init_source_dir <- NULL
if (nzchar(INIT_SOURCE_ROOT)) {
  source("R/posterior_io_utils.R")
  init_source_dir <- file.path(
    INIT_SOURCE_ROOT,
    RUN_ID,
    sprintf("all_lines__%s", tolower(RUN_ID))
  )
  init_list <- load_optim_init_from_dir(
    path = init_source_dir,
    chain_id = TASK_ID + INIT_OFFSET,
    maxG0 = max(stan_data$G0_per_well)
  )
  init_spec <- function(chain_id = 1L) init_list
}

manifest$command_args$init_source_dir <- init_source_dir
manifest$command_args$init_source_rank <- TASK_ID + INIT_OFFSET
write_run_manifest(OUTPUT_DIR, manifest)

res <- tryCatch({
  opt <- mod$optimize(
    data = stan_data,
    algorithm = "bfgs",
    init = init_spec,
    refresh = 0,
    threads = NUM_THREADS,
    save_latent_dynamics = TRUE,
    jacobian = TRUE,
    show_messages = TRUE
  )

  rc <- opt$return_codes()[1]
  if (rc == 0L) {
    list(
      rc = rc,
      lp = opt$lp(),
      draws = opt$draws(format = "matrix")
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
    list(
      rc = rc,
      lp = c(csv_fit$point_estimates[, "lp__"]),
      draws = csv_fit$point_estimates
    )
  }
}, error = function(e) {
  stop(sprintf("Optimization start failed for task %d: %s", TASK_ID, conditionMessage(e)))
})

saveRDS(res$draws, file.path(OUTPUT_DIR, sprintf("optim_draws_%d.Rds", TASK_ID)))
saveRDS(res$lp, file.path(OUTPUT_DIR, sprintf("optim_lp_%d.Rds", TASK_ID)))
saveRDS(res$rc, file.path(OUTPUT_DIR, sprintf("optim_rc_%d.Rds", TASK_ID)))

cat(">>> Optimization start completed.\n")
