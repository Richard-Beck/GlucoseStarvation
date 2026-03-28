#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_PATH <- if (length(args) >= 1) args[1] else stop("output_path required")
DATASET_MANIFEST_PATH <- if (length(args) >= 2) args[2] else stop("dataset_manifest required")
RUN_ID_ARG <- if (length(args) >= 3) args[3] else stop("run_ids required")
N_STARTS <- if (length(args) >= 4) as.integer(args[4]) else 250L
NUM_THREADS <- if (length(args) >= 5) as.integer(args[5]) else 8L
MODEL_NAME <- if (length(args) >= 6) args[6] else "gpath"
MODEL_VERSION <- if (length(args) >= 7) args[7] else "v1"
RUN_LABEL <- if (length(args) >= 8) args[8] else ""
ARRAY_MEM_GB <- if (length(args) >= 9) as.integer(args[9]) else 8L
ARRAY_TIME <- if (length(args) >= 10) args[10] else "00:15:00"
COMBINE_MEM_GB <- if (length(args) >= 11) as.integer(args[11]) else 4L
COMBINE_TIME <- if (length(args) >= 12) args[12] else "00:15:00"
QOS <- if (length(args) >= 13) args[13] else "medium"
DELETE_AFTER <- if (length(args) >= 14) as.integer(args[14]) else 1L

read_run_ids <- function(x) {
  clean_vals <- function(vals) {
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    vals <- vals[!grepl("^#", vals)]
    vals <- gsub('^"(.*)"$', "\\1", vals)
    vals <- gsub("^'(.*)'$", "\\1", vals)
    vals[nzchar(vals)]
  }

  if (file.exists(x)) {
    vals <- readLines(x, warn = FALSE)
    clean_vals(vals)
  } else {
    vals <- strsplit(x, ",", fixed = TRUE)[[1]]
    clean_vals(vals)
  }
}

if (!file.exists(DATASET_MANIFEST_PATH)) {
  stop(sprintf("Dataset manifest not found: %s", DATASET_MANIFEST_PATH))
}

run_ids <- read_run_ids(RUN_ID_ARG)
if (!length(run_ids)) {
  stop("No run_ids supplied")
}

datasets <- read.delim(DATASET_MANIFEST_PATH, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("dataset_label", "stan_data_path")
missing_cols <- setdiff(required_cols, names(datasets))
if (length(missing_cols)) {
  stop(sprintf("Dataset manifest missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

datasets <- datasets[nzchar(datasets$dataset_label) & nzchar(datasets$stan_data_path), , drop = FALSE]
if (!nrow(datasets)) {
  stop("Dataset manifest contains no usable rows")
}

spec_rows <- vector("list", length(run_ids) * nrow(datasets))
idx <- 1L
for (dataset_i in seq_len(nrow(datasets))) {
  for (run_id in run_ids) {
    spec_rows[[idx]] <- data.frame(
      enabled = 1L,
      model_name = MODEL_NAME,
      model_version = MODEL_VERSION,
      run_id = run_id,
      n_starts = N_STARTS,
      num_threads = NUM_THREADS,
      stan_data_path = datasets$stan_data_path[[dataset_i]],
      dataset_label = datasets$dataset_label[[dataset_i]],
      run_label = RUN_LABEL,
      array_cpus = NUM_THREADS,
      array_mem_gb = ARRAY_MEM_GB,
      array_time = ARRAY_TIME,
      combine_mem_gb = COMBINE_MEM_GB,
      combine_time = COMBINE_TIME,
      qos = QOS,
      delete_after = DELETE_AFTER,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

spec <- do.call(rbind, spec_rows)
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
write.table(spec, file = OUTPUT_PATH, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(
  "Wrote %d optimisation rows (%d datasets x %d run_ids) to %s\n",
  nrow(spec),
  nrow(datasets),
  length(run_ids),
  OUTPUT_PATH
))
