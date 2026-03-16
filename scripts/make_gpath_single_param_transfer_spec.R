#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

output_path <- if (length(args) >= 1) args[1] else stop("output_path required")
run_id_arg <- if (length(args) >= 2) args[2] else file.path("models", "gpath", "v1", "assessment_run_ids.txt")
start_spec <- if (length(args) >= 3) args[3] else "1:50"
num_threads <- if (length(args) >= 4) as.integer(args[4]) else 16L
stan_data_path <- if (length(args) >= 5) args[5] else file.path("data", "inputs", "stan", "gstarvation_v1", "stan_ready_data.Rds")
output_root <- if (length(args) >= 6) args[6] else "data/gpath_transfer_cv_single_param"
iter_sampling <- if (length(args) >= 7) as.integer(args[7]) else 1000L
max_treedepth <- if (length(args) >= 8) as.integer(args[8]) else 12L
init_mode <- if (length(args) >= 9) args[9] else "random"
qos <- if (length(args) >= 10) args[10] else "medium"
target_regex <- if (length(args) >= 11) args[11] else ".*"

source("models/gpath/v1/model.R")

read_run_ids <- function(path_or_csv) {
  clean_vals <- function(vals) {
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    vals <- vals[!grepl("^#", vals)]
    vals <- gsub('^"(.*)"$', "\\1", vals)
    vals <- gsub("^'(.*)'$", "\\1", vals)
    vals[nzchar(vals)]
  }

  if (file.exists(path_or_csv)) {
    clean_vals(readLines(path_or_csv, warn = FALSE))
  } else {
    clean_vals(strsplit(path_or_csv, ",", fixed = TRUE)[[1]])
  }
}

run_ids <- read_run_ids(run_id_arg)
if (!length(run_ids)) {
  stop("No run_ids supplied")
}

rows <- list()
idx <- 1L

for (run_id in run_ids) {
  dims <- parse_run_id(run_id)
  labels <- get_param_names(dims$R, dims$P, dims$W, dims$C, dims$M)
  labels <- labels[grepl(target_regex, labels)]

  if (!length(labels)) {
    next
  }

  for (label in labels) {
    mask_info <- build_ploidy_effect_mask(
      R = dims$R,
      P = dims$P,
      W = dims$W,
      C = dims$C,
      M = dims$M,
      target_spec = label
    )

    rows[[idx]] <- data.frame(
      enabled = 1L,
      model_name = "gpath",
      run_id = run_id,
      line_spec = "all",
      direction_spec = "low_to_high,high_to_low",
      fit_spec = "null,transfer,oracle",
      start_spec = start_spec,
      iter_warmup = 500L,
      iter_sampling = iter_sampling,
      adapt_delta = 0.99,
      max_treedepth = max_treedepth,
      num_threads = num_threads,
      stan_data_path = stan_data_path,
      output_root = file.path(output_root, mask_info$label),
      ploidy_effect_mask_spec = label,
      run_prefit = 0L,
      prefit_chains = 4L,
      init_mode = init_mode,
      qos = qos,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

spec <- if (length(rows)) {
  do.call(rbind, rows)
} else {
  data.frame(
    enabled = integer(),
    model_name = character(),
    run_id = character(),
    line_spec = character(),
    direction_spec = character(),
    fit_spec = character(),
    start_spec = character(),
    iter_warmup = integer(),
    iter_sampling = integer(),
    adapt_delta = numeric(),
    max_treedepth = integer(),
    num_threads = integer(),
    stan_data_path = character(),
    output_root = character(),
    ploidy_effect_mask_spec = character(),
    run_prefit = integer(),
    prefit_chains = integer(),
    init_mode = character(),
    qos = character(),
    stringsAsFactors = FALSE
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(spec, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %d single-parameter transfer screening jobs to %s\n", nrow(spec), output_path))
