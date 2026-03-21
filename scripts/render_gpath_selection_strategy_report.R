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

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
config_path <- args$config %||% file.path("scripts", "selection_strategy_config.R")
output_file <- args$output_file %||% file.path("workflow", "selection_strategy_analysis.html")

source("R/selection_strategy_utils.R")
cfg <- load_selection_config(config_path = config_path)

rmarkdown::render(
  input = file.path("workflow", "selection_strategy_analysis.Rmd"),
  output_file = output_file,
  params = list(
    config_path = config_path,
    output_root = resolve_selection_output_root(cfg)
  ),
  envir = new.env(parent = globalenv())
)
