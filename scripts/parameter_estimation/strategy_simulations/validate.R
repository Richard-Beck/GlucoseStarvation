#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: validate.R RELEASE_CONFIG", call. = FALSE)
suppressPackageStartupMessages({library(jsonlite); library(qs)})
source("scripts/parameter_estimation/common.R")
source("R/gpath_posterior_strategy_utils.R")

cfg <- pe_read_config(args[[1]])
strategy_cfg <- gps_read_strategy_config(cfg)
root <- gps_output_root(cfg, strategy_cfg)
tasks <- pe_read_tsv(
  file.path(root, "plans", "tasks.tsv"),
  c("task_id", "dataset_id", "model_id", "line_id", "line_name", "output_path")
)
expected_time <- length(seq(
  0,
  as.numeric(strategy_cfg$interval_hours) * as.integer(strategy_cfg$segments),
  by = as.numeric(strategy_cfg$time_step_hours)
))
schedule <- gps_build_strategy_panel(strategy_cfg)
expected_schedules <- nrow(schedule)
compact_mode <- identical(
  as.character(strategy_cfg$storage$mode),
  "compact_observation_quantiles"
)
expected_draws <- as.integer(strategy_cfg$posterior_draws)
expected_densities <- if (compact_mode) {
  length(strategy_cfg$seed_total_multipliers)
} else {
  1L
}
expected_quantiles <- if (compact_mode) length(strategy_cfg$storage$quantiles) else 0L
approximately_equal <- function(x, y, tolerance = 1e-9) {
  scale <- 1 + max(abs(c(x, y)), na.rm = TRUE)
  all(is.finite(c(x, y))) && max(abs(x - y)) <= tolerance * scale
}
records <- vector("list", nrow(tasks))

for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  path <- task$output_path[[1]]
  valid <- FALSE
  detail <- "missing"
  dims <- NA_character_
  if (file.exists(path)) {
    check <- tryCatch({
      x <- qs::qread(path, nthreads = 2L)
      if (compact_mode) {
        state_count <- length(x$state_names)
        observation_dim <- c(
          expected_draws, expected_densities, expected_schedules, 6L, state_count
        )
        quantile_dim <- c(
          expected_quantiles, expected_densities, expected_schedules,
          expected_time, state_count
        )
        valid <- identical(as.integer(x$schema_version), 2L) &&
          identical(dim(x$observation_state), observation_dim) &&
          identical(dim(x$trajectory_quantiles), quantile_dim) &&
          all(is.finite(x$observation_state)) &&
          all(is.finite(x$trajectory_quantiles)) &&
          length(x$fitted_seed_total_n) == expected_draws &&
          identical(dim(x$actual_seed_total_n), c(expected_draws, expected_densities)) &&
          all(is.finite(x$fitted_seed_total_n)) &&
          all(is.finite(x$actual_seed_total_n))
        if (valid) {
          expected_seed <- outer(
            x$fitted_seed_total_n,
            as.numeric(strategy_cfg$seed_total_multipliers),
            `*`
          )
          valid <- approximately_equal(x$actual_seed_total_n, expected_seed)
        }
        if (valid) {
          low0 <- x$observation_state[, , , "day0_post_seed", "N_low"]
          high0 <- x$observation_state[, , , "day0_post_seed", "N_high"]
          valid <- approximately_equal(
            high0 / (low0 + high0),
            rep(as.numeric(strategy_cfg$initial_high_fraction), length(high0))
          )
        }
        if (valid && all(c("0,C,C", "0,A0,A0") %in% x$schedule_index$strategy_code)) {
          carry_i <- match("0,C,C", x$schedule_index$strategy_code)
          add_i <- match("0,A0,A0", x$schedule_index$strategy_code)
          valid <- approximately_equal(
            x$observation_state[, , carry_i, , ],
            x$observation_state[, , add_i, , ]
          ) && approximately_equal(
            x$trajectory_quantiles[, , carry_i, , ],
            x$trajectory_quantiles[, , add_i, , ]
          )
        }
        if (valid) {
          add_i <- which(x$schedule_index$day2_action == "add_glucose")
          for (j in add_i) {
            increment <- as.numeric(x$schedule_index$g0_day2[[j]])
            before <- x$observation_state[, , j, "day2_pre_action", , drop = FALSE]
            after <- x$observation_state[, , j, "day2_post_action", , drop = FALSE]
            g_idx <- match("G1", x$state_names)
            if (!approximately_equal(after[, , , , g_idx] - before[, , , , g_idx], increment) ||
                !approximately_equal(
                  after[, , , , -g_idx, drop = FALSE],
                  before[, , , , -g_idx, drop = FALSE]
                )) {
              valid <- FALSE
              break
            }
          }
        }
        state_dims <- paste(dim(x$observation_state), collapse = "x")
      } else {
        endpoint_dim <- dim(x$endpoint)
        legacy_dims <- unique(vapply(
          x$state, function(y) paste(dim(y), collapse = "x"), character(1)
        ))
        valid <- length(x$state) == expected_draws &&
          identical(as.integer(endpoint_dim[1:2]), c(expected_draws, expected_schedules)) &&
          length(legacy_dims) == 1L &&
          identical(
            as.integer(strsplit(legacy_dims, "x", fixed = TRUE)[[1]][1:2]),
            c(expected_schedules, expected_time)
          ) &&
          length(x$seed_total_n) == expected_draws &&
          all(is.finite(x$seed_total_n))
        state_dims <- paste(length(x$state), legacy_dims[[1]], sep = "x")
      }
      list(
        valid = valid,
        detail = if (valid) "complete" else "dimension or completion mismatch",
        dims = state_dims
      )
    }, error = function(e) {
      list(valid = FALSE, detail = conditionMessage(e), dims = NA_character_)
    })
    valid <- check$valid
    detail <- check$detail
    dims <- check$dims
  }
  records[[i]] <- data.frame(
    task_id = task$task_id,
    dataset_id = task$dataset_id,
    model_id = task$model_id,
    line_id = task$line_id,
    line_name = task$line_name,
    output_path = path,
    exists = file.exists(path),
    valid = valid,
    bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
    state_dimensions = dims,
    detail = detail,
    sha256 = if (valid) pe_sha256(path) else NA_character_,
    stringsAsFactors = FALSE
  )
  cat(sprintf("validated=%d/%d valid=%s\n", i, nrow(tasks), valid))
}
manifest <- do.call(rbind, records)
utils::write.table(
  manifest, file.path(root, "manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
passed <- nrow(manifest) == 70L && all(manifest$valid)
report <- list(
  schema_version = if (compact_mode) 2 else 1,
  validation_passed = passed,
  generated_at = format(Sys.time(), usetz = TRUE),
  expected = list(
    tasks = 70L, draws = expected_draws, densities = expected_densities,
    schedules = expected_schedules, time_points = expected_time,
    quantiles = expected_quantiles
  ),
  observed = list(tasks = nrow(manifest), valid_tasks = sum(manifest$valid),
                  total_bytes = sum(manifest$bytes, na.rm = TRUE)),
  errors = as.list(manifest$detail[!manifest$valid])
)
pe_write_json(report, file.path(root, "validation.json"))
if (!passed) pe_fail("Strategy validation failed: %d/%d valid", sum(manifest$valid), nrow(manifest))
cat(sprintf("Strategy validation passed: %d files, %.3f GB\n",
            nrow(manifest), sum(manifest$bytes) / 1024^3))
