#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_ROOT <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "inputs", "stan", "wp4b_no_sum159_fuse")
}
DATASET_PREFIX <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  "gstarvation_v1_no_sum159_fuse"
}
SOURCE_STAN_DATA_PATH <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("data", "stan_ready_data.Rds")
}
RUN_ID_SOURCE <- if (length(args) >= 4 && nzchar(args[4])) {
  args[4]
} else {
  file.path("models", "gpath", "v1", "assessment_run_ids.txt")
}
SMOKE_RUN_IDS <- if (length(args) >= 5 && nzchar(args[5])) {
  args[5]
} else {
  "1R_1P_0W_C0_M1,2R_2P_1W_C0_M0"
}
SMOKE_LINES <- if (length(args) >= 6 && nzchar(args[6])) {
  args[6]
} else {
  "1"
}
N_STARTS_SMOKE <- if (length(args) >= 7 && nzchar(args[7])) as.integer(args[7]) else 5L
N_STARTS_ASSESSMENT <- if (length(args) >= 8 && nzchar(args[8])) as.integer(args[8]) else 500L

EXCLUDED_LINE <- "SUM-159-fuse"
DATASET_MANIFEST_PATH <- file.path("data", "specs", "wp4b_no_sum159_fuse_datasets.tsv")
SMOKE_SPEC_PATH <- file.path("data", "specs", "wp4b_no_sum159_fuse_smoke.tsv")
ASSESSMENT_SPEC_PATH <- file.path("data", "specs", "wp4b_no_sum159_fuse_assessment.tsv")
SUMMARY_PATH <- file.path("data", "report_exports", "wp4b_no_sum159_fuse_gpath", "wp4b_input_summary.txt")

source("R/project_paths.R")
source("R/gpath_run_utils.R")

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
    return(clean_vals(readLines(x, warn = FALSE)))
  }

  clean_vals(strsplit(x, ",", fixed = TRUE)[[1]])
}

build_dataset_label <- function(prefix, fit_type, line_id = NULL, direction = NULL) {
  if (identical(fit_type, "full")) {
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

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_optim_spec <- function(
  datasets,
  run_ids,
  output_path,
  n_starts,
  run_label,
  num_threads = 8L,
  array_mem_gb = 8L,
  array_time = "00:15:00",
  combine_mem_gb = 4L,
  combine_time = "00:15:00",
  qos = "medium",
  delete_after = 1L
) {
  if (!nrow(datasets)) {
    stop("No datasets were supplied for spec generation")
  }
  if (!length(run_ids)) {
    stop("No run_ids were supplied for spec generation")
  }

  rows <- vector("list", nrow(datasets) * length(run_ids))
  idx <- 1L
  for (dataset_i in seq_len(nrow(datasets))) {
    for (run_id in run_ids) {
      rows[[idx]] <- data.frame(
        enabled = 1L,
        model_name = "gpath",
        model_version = "v1",
        run_id = run_id,
        n_starts = as.integer(n_starts),
        num_threads = as.integer(num_threads),
        stan_data_path = datasets$stan_data_path[[dataset_i]],
        dataset_label = datasets$dataset_label[[dataset_i]],
        run_label = run_label,
        array_cpus = as.integer(num_threads),
        array_mem_gb = as.integer(array_mem_gb),
        array_time = array_time,
        combine_mem_gb = as.integer(combine_mem_gb),
        combine_time = combine_time,
        qos = qos,
        delete_after = as.integer(delete_after),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  write_tsv(do.call(rbind, rows), output_path)
}

source_stan_data_path <- resolve_stan_data_path(SOURCE_STAN_DATA_PATH)
source_stan_data <- readRDS(source_stan_data_path)

if (!("line_map" %in% names(source_stan_data)) || is.null(source_stan_data$line_map)) {
  stop("Source Stan data must contain line_map so SUM-159-fuse can be excluded by name")
}
if (!(EXCLUDED_LINE %in% names(source_stan_data$line_map))) {
  stop(sprintf("Could not find excluded line '%s' in line_map", EXCLUDED_LINE))
}

old_line_map <- as.integer(unname(source_stan_data$line_map))
names(old_line_map) <- names(source_stan_data$line_map)
excluded_line_id <- unname(old_line_map[[EXCLUDED_LINE]])
keep_line_ids <- sort(setdiff(unique(as.integer(source_stan_data$line_id)), excluded_line_id))

no_sum_stan_data <- subset_stan_data_to_lines(source_stan_data, keep_line_ids)
no_sum_stan_data$subset_excluded_line_id <- as.integer(excluded_line_id)
no_sum_stan_data$subset_excluded_line_name <- EXCLUDED_LINE

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
no_sum_path <- file.path(OUTPUT_ROOT, "stan_ready_data.Rds")
saveRDS(no_sum_stan_data, no_sum_path)

line_lookup <- data.frame(
  line_id = as.integer(no_sum_stan_data$line_map),
  line_name = names(no_sum_stan_data$line_map),
  stringsAsFactors = FALSE
)

dataset_rows <- list()
idx <- 1L
dataset_rows[[idx]] <- data.frame(
  dataset_label = build_dataset_label(DATASET_PREFIX, "full"),
  stan_data_path = no_sum_path,
  fit_type = "full",
  line_id = NA_integer_,
  line_name = NA_character_,
  direction = NA_character_,
  observed_value = NA_real_,
  holdout_value = NA_real_,
  holdout_wells = NA_integer_,
  holdout_obs = NA_integer_,
  source_stan_data_path = source_stan_data_path,
  stringsAsFactors = FALSE
)
idx <- idx + 1L

dataset_root <- file.path(OUTPUT_ROOT, "datasets")
dir.create(dataset_root, recursive = TRUE, showWarnings = FALSE)
directions <- c("low_to_high", "high_to_low")

for (line_id in sort(unique(as.integer(no_sum_stan_data$line_id)))) {
  line_name <- line_lookup$line_name[match(line_id, line_lookup$line_id)]
  for (direction in directions) {
    split_meta <- get_directional_transfer_split(
      stan_data = no_sum_stan_data,
      line_id = line_id,
      direction = direction
    )

    for (fit_type in c("transfer", "null")) {
      dataset_label <- build_dataset_label(
        prefix = DATASET_PREFIX,
        fit_type = fit_type,
        line_id = line_id,
        direction = direction
      )
      out_dir <- file.path(dataset_root, dataset_label)
      out_path <- file.path(out_dir, "stan_ready_data.Rds")
      meta_path <- file.path(out_dir, "dataset_meta.Rds")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      derived_stan_data <- apply_directional_transfer_holdout(
        stan_data = no_sum_stan_data,
        line_id = line_id,
        direction = direction
      )
      if (identical(fit_type, "null")) {
        derived_stan_data$ploidy_metric[] <- 0.0
      }

      saveRDS(derived_stan_data, out_path)
      write_dataset_meta(meta_path, list(
        dataset_label = dataset_label,
        fit_type = fit_type,
        line_id = as.integer(line_id),
        line_name = line_name,
        direction = direction,
        observed_value = split_meta$observed_value,
        holdout_value = split_meta$holdout_value,
        holdout_wells = split_meta$holdout_wells,
        holdout_obs = split_meta$holdout_total_obs,
        source_stan_data_path = source_stan_data_path,
        no_sum_stan_data_path = no_sum_path,
        derived_stan_data_path = out_path,
        excluded_line_id = as.integer(excluded_line_id),
        excluded_line_name = EXCLUDED_LINE
      ))

      dataset_rows[[idx]] <- data.frame(
        dataset_label = dataset_label,
        stan_data_path = out_path,
        fit_type = fit_type,
        line_id = as.integer(line_id),
        line_name = line_name,
        direction = direction,
        observed_value = split_meta$observed_value,
        holdout_value = split_meta$holdout_value,
        holdout_wells = length(split_meta$holdout_wells),
        holdout_obs = split_meta$holdout_total_obs,
        source_stan_data_path = source_stan_data_path,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
}

dataset_manifest <- do.call(rbind, dataset_rows)
write_tsv(dataset_manifest, DATASET_MANIFEST_PATH)

assessment_run_ids <- read_run_ids(RUN_ID_SOURCE)
smoke_run_ids <- read_run_ids(SMOKE_RUN_IDS)
smoke_lines <- as.integer(parse_csv_arg(SMOKE_LINES))
smoke_datasets <- dataset_manifest[
  dataset_manifest$fit_type == "full" |
    (!is.na(dataset_manifest$line_id) & dataset_manifest$line_id %in% smoke_lines),
  ,
  drop = FALSE
]

write_optim_spec(
  datasets = smoke_datasets,
  run_ids = smoke_run_ids,
  output_path = SMOKE_SPEC_PATH,
  n_starts = N_STARTS_SMOKE,
  run_label = "wp4b_no_sum159_fuse_smoke"
)

write_optim_spec(
  datasets = dataset_manifest,
  run_ids = assessment_run_ids,
  output_path = ASSESSMENT_SPEC_PATH,
  n_starts = N_STARTS_ASSESSMENT,
  run_label = "wp4b_no_sum159_fuse_assessment"
)

dir.create(dirname(SUMMARY_PATH), recursive = TRUE, showWarnings = FALSE)
summary_lines <- c(
  "WP4B_NO_SUM159_FUSE_GPATH_INPUTS",
  sprintf("generated\t%s", Sys.time()),
  sprintf("source_stan_data\t%s", source_stan_data_path),
  sprintf("no_sum_stan_data\t%s", no_sum_path),
  sprintf("excluded_line\t%s", EXCLUDED_LINE),
  sprintf("excluded_line_id\t%d", as.integer(excluded_line_id)),
  sprintf("N_lines\t%d", no_sum_stan_data$N_lines),
  sprintf("N_wells\t%d", no_sum_stan_data$N_wells),
  sprintf("N_obs_count\t%d", no_sum_stan_data$N_obs_count),
  sprintf("N_obs_gluc\t%d", no_sum_stan_data$N_obs_gluc),
  sprintf("dataset_manifest\t%s", DATASET_MANIFEST_PATH),
  sprintf("n_datasets\t%d", nrow(dataset_manifest)),
  sprintf("smoke_spec\t%s", SMOKE_SPEC_PATH),
  sprintf("n_smoke_rows\t%d", nrow(smoke_datasets) * length(smoke_run_ids)),
  sprintf("assessment_spec\t%s", ASSESSMENT_SPEC_PATH),
  sprintf("n_assessment_rows\t%d", nrow(dataset_manifest) * length(assessment_run_ids)),
  "",
  "Suggested execution order:",
  sprintf("1. bash slurm/submit_gpath_optim_pipeline.sh %s", SMOKE_SPEC_PATH),
  sprintf("2. After smoke results pass, bash slurm/submit_gpath_optim_pipeline.sh %s", ASSESSMENT_SPEC_PATH),
  "3. Summarize completed results with scripts/agentRrunner.sh scripts/make_wp4b_no_sum_gpath_report.R"
)
writeLines(summary_lines, SUMMARY_PATH, useBytes = TRUE)

cat(sprintf("Wrote no-SUM Stan data to %s\n", no_sum_path))
cat(sprintf("Wrote dataset manifest to %s\n", DATASET_MANIFEST_PATH))
cat(sprintf("Wrote smoke spec to %s\n", SMOKE_SPEC_PATH))
cat(sprintf("Wrote assessment spec to %s\n", ASSESSMENT_SPEC_PATH))
cat(sprintf("Wrote input summary to %s\n", SUMMARY_PATH))
