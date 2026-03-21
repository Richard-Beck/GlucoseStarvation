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

config_path <- args$config %||% file.path("scripts", "selection_strategy_config.R")
overrides <- list()
if (!is.null(args$workers)) {
  overrides$workers <- as.integer(args$workers)
}
if (!is.null(args$output_root)) {
  overrides$output_root <- args$output_root
}
if (!is.null(args$dataset_label)) {
  overrides$dataset_label <- args$dataset_label
}
if (!is.null(args$model_ids) && nzchar(args$model_ids)) {
  overrides$model_ids <- trimws(strsplit(args$model_ids, ",", fixed = TRUE)[[1]])
}
if (!is.null(args$line_ids) && nzchar(args$line_ids)) {
  overrides$line_ids <- as.integer(trimws(strsplit(args$line_ids, ",", fixed = TRUE)[[1]]))
}
if (!is.null(args$include_global)) {
  overrides$include_global <- coerce_flag(args$include_global, default = TRUE)
}
if (!is.null(args$include_transfer)) {
  overrides$include_transfer <- coerce_flag(args$include_transfer, default = TRUE)
}

source("R/selection_strategy_utils.R")

cfg <- load_selection_config(config_path = config_path, overrides = overrides)
cfg$project_root <- getwd()
detailed <- coerce_flag(args$detailed, default = FALSE)

cat(sprintf(">>> run_gpath_selection_strategy_batch.R cwd: %s\n", getwd()))
cat(sprintf(">>> config: %s\n", config_path))
cat(sprintf(">>> output root: %s\n", resolve_selection_output_root(cfg)))
cat(sprintf(">>> detailed trajectories: %s\n", detailed))

tasks <- build_selection_tasks(cfg, resolve_selection_stan(cfg))
cat(sprintf(">>> queued tasks: %d\n", nrow(tasks)))

run_selection_tasks(
  config = cfg,
  tasks = tasks,
  detailed = detailed
)

cat(">>> selection strategy batch complete.\n")
