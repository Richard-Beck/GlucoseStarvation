#!/usr/bin/env Rscript

# Reduce validated full posterior strategy trajectories to one compact endpoint
# shard per modelling context. The full .qs sources remain canonical and are not
# duplicated: each output retains only endpoint arrays, fitted seed sizes, draw
# indices, schedule metadata, and source identities.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("usage: reduce_endpoints.R RELEASE_CONFIG [--overwrite]", call. = FALSE)
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(qs)
})
source("scripts/parameter_estimation/common.R")
source("R/gpath_posterior_strategy_utils.R")

cfg <- pe_read_config(args[[1]])
overwrite <- "--overwrite" %in% args[-1L]
strategy_cfg <- gps_read_strategy_config(cfg)
if (identical(as.character(strategy_cfg$storage$mode), "compact_observation_quantiles")) {
  pe_fail("Endpoint reduction is unnecessary for schema-v2 compact strategy artifacts")
}
root <- gps_output_root(cfg, strategy_cfg)
validation_path <- file.path(root, "validation.json")
if (!file.exists(validation_path)) {
  pe_fail("Strategy validation is absent: %s", validation_path)
}
validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
if (!isTRUE(validation$validation_passed) ||
    as.integer(validation$observed$valid_tasks) != as.integer(validation$expected$tasks)) {
  pe_fail("Full posterior strategy simulations have not passed validation")
}

tasks <- pe_read_tsv(
  file.path(root, "plans", "tasks.tsv"),
  c("task_id", "dataset_id", "model_id", "line_id", "line_name", "output_path")
)
manifest <- pe_read_tsv(
  file.path(root, "manifest.tsv"),
  c("task_id", "output_path", "valid", "sha256")
)
manifest$valid <- tolower(as.character(manifest$valid)) %in% c("true", "t", "1")
tasks <- merge(
  tasks,
  manifest[, c("task_id", "valid", "sha256")],
  by = "task_id",
  all.x = TRUE,
  sort = FALSE
)
tasks <- tasks[order(as.integer(tasks$task_id)), , drop = FALSE]
if (nrow(tasks) != 70L || any(!tasks$valid) || any(!file.exists(tasks$output_path))) {
  pe_fail("Expected 70 present, validator-passing posterior strategy sources")
}

expected_draws <- as.integer(strategy_cfg$posterior_draws)
expected_schedules <- 448L
expected_metrics <- c(
  "N_low", "N_high", "total_live", "high_fraction",
  "log_ratio_high_low", "low_fold_change", "high_fold_change"
)
contexts <- unique(as.character(tasks$dataset_id))
output_dir <- file.path(root, "endpoints")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
nthreads <- max(1L, min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 8L))

reduce_context <- function(dataset_id) {
  context_tasks <- tasks[tasks$dataset_id == dataset_id, , drop = FALSE]
  context_tasks <- context_tasks[order(
    match(context_tasks$model_id, unique(tasks$model_id)),
    as.integer(context_tasks$line_id)
  ), , drop = FALSE]
  output_path <- file.path(output_dir, paste0(dataset_id, ".Rds"))
  if (file.exists(output_path) && !overwrite) {
    pe_fail("Refusing to overwrite endpoint shard: %s", output_path)
  }

  endpoint <- array(
    NA_real_,
    dim = c(nrow(context_tasks), expected_draws, expected_schedules, length(expected_metrics)),
    dimnames = list(
      task = paste(context_tasks$model_id, context_tasks$line_name, sep = "::"),
      draw = NULL,
      schedule = NULL,
      metric = expected_metrics
    )
  )
  seed_total_n <- matrix(
    NA_real_,
    nrow = nrow(context_tasks),
    ncol = expected_draws,
    dimnames = list(task = dimnames(endpoint)$task, draw = NULL)
  )
  draw_index <- vector("list", nrow(context_tasks))
  schedule_index <- NULL
  source_metadata <- vector("list", nrow(context_tasks))
  started <- proc.time()[["elapsed"]]

  for (i in seq_len(nrow(context_tasks))) {
    task <- context_tasks[i, , drop = FALSE]
    x <- qs::qread(task$output_path[[1]], nthreads = nthreads)
    if (!identical(dim(x$endpoint), c(expected_draws, expected_schedules, length(expected_metrics))) ||
        !identical(dimnames(x$endpoint)[[3]], expected_metrics) ||
        length(x$seed_total_n) != expected_draws ||
        any(!is.finite(x$seed_total_n))) {
      pe_fail("Endpoint schema mismatch in %s", task$output_path[[1]])
    }
    if (is.null(schedule_index)) {
      schedule_index <- x$schedule_index
      dimnames(endpoint)$draw <- as.character(x$draw_index$draw_id)
      dimnames(endpoint)$schedule <- as.character(x$schedule_index$strategy_code)
      dimnames(seed_total_n)$draw <- as.character(x$draw_index$draw_id)
    } else if (!identical(
      as.character(schedule_index$strategy_code),
      as.character(x$schedule_index$strategy_code)
    )) {
      pe_fail("Schedule index mismatch in %s", task$output_path[[1]])
    }
    endpoint[i, , , ] <- x$endpoint
    seed_total_n[i, ] <- x$seed_total_n
    draw_index[[i]] <- x$draw_index
    source_metadata[[i]] <- x$metadata
    rm(x)
    gc(verbose = FALSE)
    cat(sprintf("context=%s reduced=%d/%d\n", dataset_id, i, nrow(context_tasks)))
  }

  names(draw_index) <- dimnames(endpoint)$task
  names(source_metadata) <- dimnames(endpoint)$task
  source_index <- context_tasks[, c(
    "task_id", "dataset_id", "model_id", "line_id", "line_name",
    "output_path", "sha256"
  )]
  artifact <- list(
    schema_version = 1L,
    metadata = list(
      release_id = cfg$release_id,
      dataset_id = dataset_id,
      source_validation_path = validation_path,
      source_validation_generated_at = validation$generated_at,
      source_tasks = nrow(context_tasks),
      posterior_draws = expected_draws,
      schedules = expected_schedules,
      metrics = expected_metrics,
      reduction = "endpoints_only_no_full_trajectory_duplication",
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      generated_at = format(Sys.time(), usetz = TRUE)
    ),
    task_index = source_index,
    draw_index = draw_index,
    schedule_index = schedule_index,
    endpoint = endpoint,
    seed_total_n = seed_total_n,
    source_metadata = source_metadata
  )

  tmp_path <- tempfile(pattern = paste0(".", basename(output_path), "."),
                       tmpdir = output_dir)
  on.exit(unlink(tmp_path), add = TRUE)
  saveRDS(artifact, tmp_path, compress = "gzip")
  check <- readRDS(tmp_path)
  if (!identical(dim(check$endpoint), dim(endpoint)) ||
      !identical(check$task_index[, c("dataset_id", "model_id", "line_id", "line_name")],
                 source_index[, c("dataset_id", "model_id", "line_id", "line_name")]) ||
      any(!is.finite(check$seed_total_n))) {
    pe_fail("Endpoint shard failed immediate readback validation: %s", output_path)
  }
  rm(check)
  if (file.exists(output_path)) unlink(output_path)
  if (!file.rename(tmp_path, output_path)) pe_fail("Could not install %s", output_path)
  cat(sprintf("Wrote %s (%d tasks x %d draws x %d schedules)\n",
              output_path, nrow(context_tasks), expected_draws, expected_schedules))
  invisible(output_path)
}

outputs <- vapply(contexts, reduce_context, character(1))
if (length(outputs) != 7L || any(!file.exists(outputs))) {
  pe_fail("Expected seven completed context-level endpoint shards")
}
cat(sprintf("Endpoint reduction complete: %d context shards\n", length(outputs)))
