#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_PATH <- if (length(args) >= 1) args[1] else stop("output manifest path required")
STAN_DATA_PATH <- if (length(args) >= 2) args[2] else ""
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else file.path("data", "inputs", "stan", "gpath_optim_datasets")
LINE_ARG <- if (length(args) >= 4) args[4] else "all"
DIRECTION_ARG <- if (length(args) >= 5) args[5] else "low_to_high,high_to_low"
FIT_ARG <- if (length(args) >= 6) args[6] else "full,transfer,null"
DATASET_PREFIX <- if (length(args) >= 7) args[7] else "gstarvation_v1"

source("R/gpath_run_utils.R")

cat(sprintf(">>> make_gpath_optim_dataset_manifest.R cwd: %s\n", getwd()))

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

normalize_dataset_fit_type <- function(x) {
  key <- tolower(trimws(x))
  valid <- c("full", "transfer", "null")
  if (!(key %in% valid)) {
    stop(sprintf("Unsupported fit_type '%s'; expected one of: %s", x, paste(valid, collapse = ", ")))
  }
  key
}

build_dataset_label <- function(prefix, fit_type, line_id = NULL, direction = NULL) {
  if (fit_type == "full") {
    return(prefix)
  }

  sprintf(
    "%s_%s_line%02d_%s",
    prefix,
    fit_type,
    as.integer(line_id),
    normalize_transfer_direction(direction)
  )
}

write_dataset_meta <- function(path, meta) {
  saveRDS(meta, path)
}

stan_data_path <- resolve_stan_data_path(STAN_DATA_PATH)
base_stan_data <- readRDS(stan_data_path)

fit_types <- vapply(parse_csv_arg(FIT_ARG), normalize_dataset_fit_type, character(1))
fit_types <- unique(fit_types)

directions <- vapply(parse_csv_arg(DIRECTION_ARG), normalize_transfer_direction, character(1))
directions <- unique(directions)

line_ids <- if (tolower(trimws(LINE_ARG)) == "all") {
  sort(unique(as.integer(base_stan_data$line_id)))
} else {
  as.integer(parse_csv_arg(LINE_ARG))
}

rows <- list()
idx <- 1L

if ("full" %in% fit_types) {
  rows[[idx]] <- data.frame(
    dataset_label = build_dataset_label(DATASET_PREFIX, "full"),
    stan_data_path = stan_data_path,
    fit_type = "full",
    line_id = NA_integer_,
    direction = NA_character_,
    observed_value = NA_real_,
    holdout_value = NA_real_,
    holdout_wells = NA_integer_,
    holdout_obs = NA_integer_,
    source_stan_data_path = stan_data_path,
    stringsAsFactors = FALSE
  )
  idx <- idx + 1L
}

derived_fit_types <- setdiff(fit_types, "full")
if (length(derived_fit_types)) {
  dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
}

for (line_id in line_ids) {
  for (direction in directions) {
    split_meta <- get_directional_transfer_split(
      stan_data = base_stan_data,
      line_id = line_id,
      direction = direction
    )

    for (fit_type in derived_fit_types) {
      dataset_label <- build_dataset_label(
        prefix = DATASET_PREFIX,
        fit_type = fit_type,
        line_id = line_id,
        direction = direction
      )

      out_dir <- file.path(OUTPUT_ROOT, dataset_label)
      out_path <- file.path(out_dir, "stan_ready_data.Rds")
      meta_path <- file.path(out_dir, "dataset_meta.Rds")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      derived_stan_data <- apply_directional_transfer_holdout(
        stan_data = base_stan_data,
        line_id = line_id,
        direction = direction
      )

      if (fit_type == "null") {
        derived_stan_data$ploidy_metric[] <- 0.0
      }

      saveRDS(derived_stan_data, out_path)
      write_dataset_meta(meta_path, list(
        dataset_label = dataset_label,
        fit_type = fit_type,
        line_id = as.integer(line_id),
        direction = direction,
        observed_value = split_meta$observed_value,
        holdout_value = split_meta$holdout_value,
        holdout_wells = split_meta$holdout_wells,
        holdout_obs = split_meta$holdout_total_obs,
        source_stan_data_path = stan_data_path,
        derived_stan_data_path = out_path
      ))

      rows[[idx]] <- data.frame(
        dataset_label = dataset_label,
        stan_data_path = out_path,
        fit_type = fit_type,
        line_id = as.integer(line_id),
        direction = direction,
        observed_value = split_meta$observed_value,
        holdout_value = split_meta$holdout_value,
        holdout_wells = length(split_meta$holdout_wells),
        holdout_obs = split_meta$holdout_total_obs,
        source_stan_data_path = stan_data_path,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
}

manifest <- do.call(rbind, rows)
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
write.table(manifest, file = OUTPUT_PATH, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote %d dataset rows to %s\n", nrow(manifest), OUTPUT_PATH))
