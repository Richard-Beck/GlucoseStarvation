#!/usr/bin/env Rscript

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L) {
      stop(sprintf("Arguments must be provided as key=value, got '%s'", arg))
    }
    out[[pieces[[1]]]] <- pieces[[2]]
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

coerce_flag <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
CONFIG_PATH <- args$config %||% file.path("scripts", "selection_strategy_config.R")
MANIFEST_PATH <- args$manifest %||% stop("manifest=<path> is required")
TASK_ID <- as.integer(args$task_id %||% stop("task_id=<int> is required"))
DETAILED <- coerce_flag(args$detailed, default = FALSE)

source("R/selection_strategy_utils.R")

cfg <- load_selection_config(CONFIG_PATH)
cfg$project_root <- getwd()

manifest <- read.delim(
  MANIFEST_PATH,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

task_row <- manifest[manifest$task_id == TASK_ID, , drop = FALSE]
if (nrow(task_row) != 1L) {
  stop(sprintf("Expected exactly one manifest row for task_id=%d", TASK_ID))
}

stan_data <- resolve_selection_stan(cfg)
result <- simulate_selection_task(
  task_row = task_row[, c("model_id", "line_id", "fit_context", "direction", "fit_type"), drop = FALSE],
  config = cfg,
  stan_data = stan_data,
  detailed = DETAILED
)

save_task_results(
  config = cfg,
  task_row = task_row[, c("model_id", "line_id", "fit_context", "direction", "fit_type"), drop = FALSE],
  result = result
)

cat(sprintf(
  "Completed selection-strategy task %d | model=%s line=%s context=%s\n",
  TASK_ID,
  task_row$model_id[[1]],
  task_row$line_id[[1]],
  task_row$fit_context[[1]]
))
