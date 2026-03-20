#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_PATH <- if (length(args) >= 1) args[1] else stop("output manifest path required")
RUN_ID_PATH <- if (length(args) >= 2) args[2] else file.path("models", "gpath", "v1", "assessment_run_ids.txt")
CHAIN_ARG <- if (length(args) >= 3) args[3] else "1,2,3,4"
LINE_ID <- if (length(args) >= 4) as.integer(args[4]) else 1L
DIRECTION <- if (length(args) >= 5) args[5] else "low_to_high"
STAN_DATA_PATH <- if (length(args) >= 6) args[6] else ""
INIT_MODE <- if (length(args) >= 7) args[7] else "optim"
OPTIM_DATASET_LABEL <- if (length(args) >= 8) args[8] else ""
OPTIM_RUN_LABEL <- if (length(args) >= 9) args[9] else ""
OPTIM_ROOT_OVERRIDE <- if (length(args) >= 10) args[10] else ""

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/elpd_transfer_utils.R")

parse_csv_arg <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",")))
  vals <- vals[nzchar(vals)]

  expand_token <- function(tok) {
    if (grepl("^[0-9]+:[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      return(as.character(seq(parts[1], parts[2])))
    }
    if (grepl("^[0-9]+-[0-9]+$", tok)) {
      parts <- as.integer(strsplit(tok, "-", fixed = TRUE)[[1]])
      return(as.character(seq(parts[1], parts[2])))
    }
    tok
  }

  unlist(lapply(vals, expand_token), use.names = FALSE)
}

resolve_init_source <- function(run_id) {
  if (nzchar(trimws(OPTIM_ROOT_OVERRIDE))) {
    return(file.path(OPTIM_ROOT_OVERRIDE, run_id))
  }

  if (nzchar(trimws(OPTIM_DATASET_LABEL)) && nzchar(trimws(OPTIM_RUN_LABEL))) {
    return(get_run_output_dir(
      model_name = "gpath",
      model_version = "v1",
      pipeline_name = "optim",
      dataset_label = OPTIM_DATASET_LABEL,
      run_id = run_id,
      run_label = OPTIM_RUN_LABEL
    ))
  }

  ""
}

run_ids <- readLines(RUN_ID_PATH, warn = FALSE)
run_ids <- trimws(run_ids)
run_ids <- run_ids[nzchar(run_ids)]

chain_ids <- as.integer(parse_csv_arg(CHAIN_ARG))
DIRECTION <- normalize_transfer_direction(DIRECTION)
STAN_DATA_PATH <- resolve_stan_data_path(STAN_DATA_PATH)
INIT_MODE <- tolower(trimws(INIT_MODE))

rows <- lapply(seq_along(run_ids), function(i) {
  run_id <- run_ids[[i]]
  do.call(rbind, lapply(chain_ids, function(chain_id) {
    data.frame(
      task_id = NA_integer_,
      run_id = run_id,
      chain_id = as.integer(chain_id),
      line_id = as.integer(LINE_ID),
      direction = DIRECTION,
      stan_data_path = STAN_DATA_PATH,
      init_mode = INIT_MODE,
      init_source = resolve_init_source(run_id),
      stringsAsFactors = FALSE
    )
  }))
})

manifest <- do.call(rbind, rows)
manifest$task_id <- seq_len(nrow(manifest))

dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
utils::write.table(manifest, file = OUTPUT_PATH, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Wrote %d tasks to %s\n", nrow(manifest), OUTPUT_PATH))
