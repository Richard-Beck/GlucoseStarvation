#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

output_path <- if (length(args) >= 1) args[1] else stop("output_path required")
run_id_arg <- if (length(args) >= 2) args[2] else stop("run_ids required")
n_starts <- if (length(args) >= 3) as.integer(args[3]) else 250L
num_threads <- if (length(args) >= 4) as.integer(args[4]) else 8L
stan_data_path <- if (length(args) >= 5) args[5] else file.path("data", "inputs", "stan", "gstarvation_v1", "stan_ready_data.Rds")
dataset_label <- if (length(args) >= 6) args[6] else "gstarvation_v1"
model_name <- if (length(args) >= 7) args[7] else "gpath"
model_version <- if (length(args) >= 8) args[8] else "v1"
run_label <- if (length(args) >= 9) args[9] else ""
array_mem_gb <- if (length(args) >= 10) as.integer(args[10]) else 8L
array_time <- if (length(args) >= 11) args[11] else "00:15:00"
combine_mem_gb <- if (length(args) >= 12) as.integer(args[12]) else 4L
combine_time <- if (length(args) >= 13) args[13] else "00:15:00"
qos <- if (length(args) >= 14) args[14] else "normal"
delete_after <- if (length(args) >= 15) as.integer(args[15]) else 1L

read_run_ids <- function(x) {
  if (file.exists(x)) {
    vals <- readLines(x, warn = FALSE)
    vals <- trimws(vals)
    vals[nzchar(vals)]
  } else {
    vals <- strsplit(x, ",", fixed = TRUE)[[1]]
    vals <- trimws(vals)
    vals[nzchar(vals)]
  }
}

run_ids <- read_run_ids(run_id_arg)

if (!length(run_ids)) {
  stop("No run_ids supplied")
}

spec <- data.frame(
  enabled = 1L,
  model_name = model_name,
  model_version = model_version,
  run_id = run_ids,
  n_starts = n_starts,
  num_threads = num_threads,
  stan_data_path = stan_data_path,
  dataset_label = dataset_label,
  run_label = run_label,
  array_cpus = num_threads,
  array_mem_gb = array_mem_gb,
  array_time = array_time,
  combine_mem_gb = combine_mem_gb,
  combine_time = combine_time,
  qos = qos,
  delete_after = delete_after,
  stringsAsFactors = FALSE
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(spec, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote %d optimisation jobs to %s\n", nrow(spec), output_path))
