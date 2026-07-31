#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
counts_path <- pe_counts_path(cfg, root)
stan_dir <- file.path(root, "datasets", "all_lines")
stan_path <- file.path(stan_dir, "stan_data.Rds")
area_path <- file.path(stan_dir, "stan_data_with_area.Rds")
area_source <- as.character(cfg$area_objects_csv %||% "data/counts/batch_results_TEMP_MEASURE_SIZE.csv")
build_area <- isTRUE(cfg$build_area_stan_data %||% FALSE)
base_stan_cfg <- cfg$base_stan_data

if (!file.exists(counts_path)) pe_fail("Release counts are missing; run build_counts first")
if (build_area && !file.exists(area_source)) pe_fail("Area-object source does not exist: %s", area_source)
if (file.exists(stan_path) || (build_area && file.exists(area_path))) pe_fail("Refusing to overwrite existing all-lines Stan data")
dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)

base_stan_source <- NULL
if (!is.null(base_stan_cfg)) {
  if (build_area) pe_fail("base_stan_data cannot be combined with build_area_stan_data")
  base_stan_source <- as.character(base_stan_cfg$input_path %||% "")
  expected_base_sha <- as.character(base_stan_cfg$expected_sha256 %||% "")
  if (!file.exists(base_stan_source)) pe_fail("Configured base Stan data is missing: %s", base_stan_source)
  if (nzchar(expected_base_sha) && !identical(pe_sha256(base_stan_source), expected_base_sha)) {
    pe_fail("Configured base Stan-data SHA256 does not match")
  }
  if (!file.copy(base_stan_source, stan_path, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)) {
    pe_fail("Failed to copy configured base Stan data to %s", stan_path)
  }
} else {
  status <- system2(
    "Rscript",
    c("R/prepare_data.R", stan_path, if (build_area) area_path else "NONE", counts_path, area_source),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) pe_fail("R/prepare_data.R failed with status %s", status)
}

stan <- readRDS(stan_path)
required <- c("N_lines", "N_wells", "N_obs_count", "N_obs_gluc", "line_map", "N_obs", "D_obs")
missing <- setdiff(required, names(stan))
if (length(missing)) pe_fail("All-lines Stan data lacks: %s", paste(missing, collapse = ", "))

# Regression gate for the July 2026 SUM-159-fuse Ploidy/ploidy bug. The glucose
# inputs are unchanged by count reprocessing, so the corrected mapping must
# retain the authenticated 912-row / 60+60 / 44-censored structure.
if ("SUM-159-fuse" %in% names(stan$line_map)) {
  affected_exp <- 3L
  affected <- which(stan$exp_id[stan$well_idx_gluc] == affected_exp)
  affected_ploidy <- stan$ploidy_abs[stan$well_idx_gluc[affected]]
  mapped_counts <- as.integer(table(factor(affected_ploidy, levels = c(2, 4))))
  if (stan$N_obs_gluc != 912L || length(affected) != 120L ||
      !identical(mapped_counts, c(60L, 60L)) ||
      sum(stan$is_censored[affected]) != 44L) {
    pe_fail(
      "SUM-159 glucose mapping regression gate failed: N_obs_gluc=%d, affected=%d, split=%s, censored=%d",
      stan$N_obs_gluc, length(affected), paste(mapped_counts, collapse = "/"),
      sum(stan$is_censored[affected])
    )
  }
}

provenance <- list(
  dataset_id = "all_lines",
  parent_dataset = NULL,
  transform = "identity",
  stan_data_path = stan_path,
  stan_data_sha256 = pe_sha256(stan_path),
  area_stan_data_path = if (build_area) area_path else NULL,
  area_stan_data_sha256 = if (build_area) pe_sha256(area_path) else NULL,
  counts_path = counts_path,
  counts_sha256 = pe_sha256(counts_path),
  base_stan_data_path = base_stan_source,
  base_stan_data_sha256 = if (!is.null(base_stan_source)) pe_sha256(base_stan_source) else NULL,
  area_objects_csv = if (build_area) area_source else NULL,
  area_objects_sha256 = if (build_area) pe_sha256(area_source) else NULL,
  builder_path = if (is.null(base_stan_source)) "R/prepare_data.R" else NULL,
  builder_sha256 = if (is.null(base_stan_source)) pe_sha256("R/prepare_data.R") else NULL,
  builder_mtime = if (is.null(base_stan_source)) format(file.info("R/prepare_data.R")$mtime, usetz = TRUE) else NULL,
  sum159_glucose_mapping_gate = "PASS: 912 total; affected batch 120 (60/60), 44 censored",
  dimensions = list(N_lines = stan$N_lines, N_wells = stan$N_wells, N_obs_count = stan$N_obs_count, N_obs_gluc = stan$N_obs_gluc),
  generated_at = format(Sys.time(), usetz = TRUE)
)
pe_write_json(provenance, file.path(stan_dir, "provenance.json"))
cat(sprintf("Built all-lines Stan data at %s\n", stan_path))
