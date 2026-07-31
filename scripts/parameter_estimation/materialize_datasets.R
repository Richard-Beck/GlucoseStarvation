#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")
source("R/gpath_run_utils.R")
source("R/morphology_metric_utils.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
spec <- pe_read_tsv(
  cfg$datasets_tsv,
  c("dataset_id", "parent_dataset", "transform", "line_name", "direction", "notes")
)
if (!nrow(spec)) pe_fail("Dataset specification is empty")
if (anyDuplicated(spec$dataset_id)) pe_fail("Dataset specification contains duplicate dataset_id values")
invisible(lapply(spec$dataset_id, pe_safe_id, label = "dataset_id"))

manifest_rows <- vector("list", nrow(spec))
available <- character()

morphology_transforms <- c("morphology_cell_area", "morphology_nuclear_area")
morphology_table <- NULL
morphology_cfg <- cfg$morphology_metrics
if (any(tolower(spec$transform) %in% morphology_transforms)) {
  if (is.null(morphology_cfg)) pe_fail("Release config must declare morphology_metrics")
  morphology_source <- as.character(morphology_cfg$summary_path %||% "")
  morphology_pair_map <- as.character(morphology_cfg$pair_map_tsv %||% "")
  expected_sha <- as.character(morphology_cfg$expected_sha256 %||% "")
  if (!file.exists(morphology_source)) pe_fail("Morphology summary does not exist: %s", morphology_source)
  if (nzchar(expected_sha) && !identical(pe_sha256(morphology_source), expected_sha)) {
    pe_fail("Morphology summary SHA256 does not match configured corrected-run source")
  }
  morphology_table <- build_morphology_metric_table(
    summary_path = morphology_source,
    pair_map_path = morphology_pair_map,
    min_hours = as.numeric(morphology_cfg$min_hours %||% 24),
    max_hours = as.numeric(morphology_cfg$max_hours %||% 48),
    object_scope = as.character(morphology_cfg$object_scope %||% "alive")
  )
  morphology_input_dir <- file.path(root, "inputs", "morphology")
  dir.create(morphology_input_dir, recursive = TRUE, showWarnings = FALSE)
  morphology_table_path <- file.path(morphology_input_dir, "metric_values.tsv")
  utils::write.table(
    morphology_table, morphology_table_path, sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
  pe_write_json(list(
    summary_path = morphology_source,
    summary_sha256 = pe_sha256(morphology_source),
    pair_map_path = morphology_pair_map,
    pair_map_sha256 = pe_sha256(morphology_pair_map),
    object_scope = as.character(morphology_cfg$object_scope %||% "alive"),
    min_hours = as.numeric(morphology_cfg$min_hours %||% 24),
    max_hours = as.numeric(morphology_cfg$max_hours %||% 48),
    glucose_filter = NULL,
    aggregation = "median of corrected-run replicate-time q50 summaries",
    metric_table_path = morphology_table_path,
    metric_table_sha256 = pe_sha256(morphology_table_path),
    generated_at = format(Sys.time(), usetz = TRUE)
  ), file.path(morphology_input_dir, "provenance.json"))
}

resolve_line_id <- function(stan_data, line_name) {
  if (is.na(line_name) || !nzchar(line_name)) pe_fail("A line_name is required for this transform")
  if (is.null(stan_data$line_map) || !(line_name %in% names(stan_data$line_map))) {
    pe_fail("Line '%s' is absent from parent dataset", line_name)
  }
  as.integer(stan_data$line_map[[line_name]])
}

set_prior_N0_center <- function(stan_data) {
  is_t0 <- stan_data$t_grid[stan_data$grid_idx_count] == 0
  center <- mean(as.numeric(stan_data$N_obs[is_t0]))
  if (!is.finite(center) || center <= 0) pe_fail("Dataset has no positive finite t=0 live-count mean")
  stan_data$prior_N0_center <- center
  stan_data
}

for (i in seq_len(nrow(spec))) {
  dataset_id <- spec$dataset_id[[i]]
  transform <- tolower(spec$transform[[i]])
  parent_id <- spec$parent_dataset[[i]]
  out_dir <- file.path(root, "datasets", dataset_id)
  out_path <- file.path(out_dir, "stan_data.Rds")

  if (transform == "identity") {
    if (dataset_id != "all_lines") pe_fail("Only all_lines may use identity transform")
    if (!file.exists(out_path)) pe_fail("all_lines data does not exist: %s", out_path)
    parent_path <- NA_character_
    line_id <- NA_integer_
    direction <- NA_character_
    stan_data <- readRDS(out_path)
  } else {
    if (is.na(parent_id) || !nzchar(parent_id) || !(parent_id %in% available)) {
      pe_fail("Dataset %s refers to unavailable/forward parent %s", dataset_id, parent_id)
    }
    parent_path <- file.path(root, "datasets", parent_id, "stan_data.Rds")
    parent <- readRDS(parent_path)
    line_id <- if (transform %in% morphology_transforms) NA_integer_ else resolve_line_id(parent, spec$line_name[[i]])
    direction <- if (is.na(spec$direction[[i]])) "" else spec$direction[[i]]

    if (transform == "exclude_line") {
      keep <- setdiff(sort(unique(as.integer(parent$line_id))), line_id)
      stan_data <- subset_stan_data_to_lines(parent, keep)
    } else if (transform == "single_line") {
      stan_data <- subset_stan_data_to_line(parent, line_id)
    } else if (transform %in% c("ploidy_holdout", "ploidy_null")) {
      stan_data <- apply_directional_transfer_holdout(parent, line_id, direction)
      if (transform == "ploidy_null") stan_data$ploidy_metric[] <- 0.0
    } else if (transform == "morphology_cell_area") {
      stan_data <- apply_morphology_metric(parent, morphology_table, "cell_area_metric")
    } else if (transform == "morphology_nuclear_area") {
      stan_data <- apply_morphology_metric(parent, morphology_table, "nuclear_area_metric")
    } else {
      pe_fail("Unsupported dataset transform '%s' for %s", transform, dataset_id)
    }

    stan_data <- set_prior_N0_center(stan_data)

    if (file.exists(out_path)) pe_fail("Refusing to overwrite derived dataset: %s", out_path)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stan_data$parameter_estimation_dataset_id <- dataset_id
    stan_data$parameter_estimation_parent_dataset <- parent_id
    stan_data$parameter_estimation_transform <- transform
    saveRDS(stan_data, out_path)
    pe_write_json(list(
      dataset_id = dataset_id,
      parent_dataset = parent_id,
      parent_path = parent_path,
      parent_sha256 = pe_sha256(parent_path),
      transform = transform,
      line_name = spec$line_name[[i]],
      parent_line_id = line_id,
      direction = if (nzchar(direction)) direction else NULL,
      morphology_source_path = if (transform %in% morphology_transforms) as.character(morphology_cfg$summary_path) else NULL,
      morphology_source_sha256 = if (transform %in% morphology_transforms) pe_sha256(as.character(morphology_cfg$summary_path)) else NULL,
      morphology_metric_column = if (transform == "morphology_cell_area") "cell_area_metric" else if (transform == "morphology_nuclear_area") "nuclear_area_metric" else NULL,
      stan_data_path = out_path,
      stan_data_sha256 = pe_sha256(out_path),
      generated_at = format(Sys.time(), usetz = TRUE)
    ), file.path(out_dir, "provenance.json"))
  }

  manifest_rows[[i]] <- data.frame(
    dataset_id = dataset_id,
    parent_dataset = if (is.na(parent_id)) "" else parent_id,
    transform = transform,
    line_name = if (is.na(spec$line_name[[i]])) "" else spec$line_name[[i]],
    parent_line_id = if (is.na(line_id)) "" else as.character(line_id),
    direction = if (is.na(direction)) "" else direction,
    stan_data_path = out_path,
    stan_data_sha256 = pe_sha256(out_path),
    N_lines = as.integer(stan_data$N_lines),
    N_wells = as.integer(stan_data$N_wells),
    N_obs_count = as.integer(stan_data$N_obs_count),
    N_obs_gluc = as.integer(stan_data$N_obs_gluc),
    notes = if (is.na(spec$notes[[i]])) "" else spec$notes[[i]],
    stringsAsFactors = FALSE
  )
  available <- c(available, dataset_id)
}

manifest <- do.call(rbind, manifest_rows)
manifest_path <- file.path(root, "manifests", "datasets.tsv")
utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("Materialized/registered %d datasets in %s\n", nrow(manifest), manifest_path))
