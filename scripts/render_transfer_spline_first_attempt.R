args <- commandArgs(trailingOnly = TRUE)

input <- if (length(args) >= 1L) {
  args[[1]]
} else {
  file.path("workflow", "transfer_spline_first_attempt.Rmd")
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the notebook.")
}

rmarkdown::render(
  input = input,
  output_format = "html_document",
  envir = new.env(parent = globalenv())
)
