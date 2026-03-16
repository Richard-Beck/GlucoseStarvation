#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

output_path <- if (length(args) >= 1) args[1] else stop("output_path required")
run_id_arg <- if (length(args) >= 2) args[2] else file.path("models", "gpath", "v1", "assessment_run_ids.txt")
start_spec <- if (length(args) >= 3) args[3] else "1:50"
num_threads <- if (length(args) >= 4) as.integer(args[4]) else 16L
stan_data_path <- if (length(args) >= 5) args[5] else file.path("data", "inputs", "stan", "gstarvation_v1", "stan_ready_data.Rds")
output_root <- if (length(args) >= 6) args[6] else "data/gpath_transfer_cv"
iter_sampling <- if (length(args) >= 7) as.integer(args[7]) else 1000L
max_treedepth <- if (length(args) >= 8) as.integer(args[8]) else 12L
init_mode <- if (length(args) >= 9) args[9] else "random"
qos <- if (length(args) >= 10) args[10] else "medium"
ploidy_effect_mask_spec <- if (length(args) >= 11) args[11] else "all"

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

spec <- data.frame(
  enabled = 1L,
  model_name = "gpath",
  run_id = run_ids,
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
  output_root = output_root,
  ploidy_effect_mask_spec = ploidy_effect_mask_spec,
  run_prefit = 0L,
  prefit_chains = 4L,
  init_mode = init_mode,
  qos = qos,
  stringsAsFactors = FALSE
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(spec, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote %d transfer screening jobs to %s\n", nrow(spec), output_path))
