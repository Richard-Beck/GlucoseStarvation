#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input <- if (length(args) >= 1L) {
  args[[1]]
} else {
  file.path("workflow", "model_analysis.Rmd")
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the notebook.")
}

rmarkdown::render(
  input = input,
  output_format = "html_document",
  envir = new.env(parent = globalenv())
)
