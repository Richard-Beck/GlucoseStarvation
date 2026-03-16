#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_PATH <- if (length(args) >= 1) args[1] else stop("output manifest path required")
STAN_DATA_PATH <- if (length(args) >= 2) args[2] else ""
LINE_ARG <- if (length(args) >= 3) args[3] else "all"
DIRECTION_ARG <- if (length(args) >= 4) args[4] else "low_to_high,high_to_low"
FIT_ARG <- if (length(args) >= 5) args[5] else "null,transfer,oracle"
START_ARG <- if (length(args) >= 6) args[6] else "1,2,3,4"
PLOIDY_EFFECT_MASK_SPEC <- if (length(args) >= 7) args[7] else ""

source("R/gpath_run_utils.R")
source("R/elpd_transfer_utils.R")

cat(sprintf(">>> make_gpath_transfer_manifest.R cwd: %s\n", getwd()))

stan_data_path <- resolve_stan_data_path(STAN_DATA_PATH)
stan_data <- readRDS(stan_data_path)

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

line_ids <- if (tolower(trimws(LINE_ARG)) == "all") {
  sort(unique(as.integer(stan_data$line_id)))
} else {
  as.integer(parse_csv_arg(LINE_ARG))
}

directions <- vapply(parse_csv_arg(DIRECTION_ARG), normalize_transfer_direction, character(1))
fit_types <- vapply(parse_csv_arg(FIT_ARG), normalize_fit_type, character(1))
start_ids <- as.integer(parse_csv_arg(START_ARG))

rows <- list()
idx <- 1L

for (line_id in line_ids) {
  for (direction in directions) {
    split_meta <- get_directional_transfer_split(
      stan_data = stan_data,
      line_id = line_id,
      direction = direction
    )

    for (fit_type in fit_types) {
      for (start_id in start_ids) {
        rows[[idx]] <- data.frame(
          task_id = idx,
          line_id = as.integer(line_id),
          direction = direction,
          fit_type = fit_type,
          start_id = as.integer(start_id),
          observed_value = split_meta$observed_value,
          holdout_value = split_meta$holdout_value,
          holdout_wells = length(split_meta$holdout_wells),
          holdout_obs = split_meta$holdout_total_obs,
          stan_data_path = stan_data_path,
          ploidy_effect_mask_spec = if (nzchar(trimws(PLOIDY_EFFECT_MASK_SPEC))) PLOIDY_EFFECT_MASK_SPEC else "all",
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
}

manifest <- do.call(rbind, rows)
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  manifest,
  file = OUTPUT_PATH,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Wrote %d tasks to %s\n", nrow(manifest), OUTPUT_PATH))
