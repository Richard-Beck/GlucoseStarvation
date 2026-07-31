#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
if (dir.exists(root)) pe_fail("Refusing to overwrite existing release: %s", root)
root <- pe_copy_release_inputs(cfg)
counts_cfg <- cfg$counts
mode <- tolower(as.character(counts_cfg$mode %||% "existing_rds"))
materialize <- isTRUE(counts_cfg$materialize %||% TRUE)
count_definition <- as.character(counts_cfg$count_definition %||% "mean_per_imaged_field_including_zero-count_fields")
out_dir <- file.path(root, "inputs", "counts")
out_rds <- file.path(out_dir, "uncorrected.Rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!materialize && mode != "existing_rds") {
  pe_fail("counts.materialize=false is supported only for existing_rds")
}
should_materialize <- mode != "existing_rds" || materialize
if (should_materialize && file.exists(out_rds)) pe_fail("Refusing to overwrite existing release counts: %s", out_rds)

source_files <- character()
field_audit <- NULL

if (mode == "existing_rds") {
  input <- as.character(counts_cfg$input_path %||% pe_fail("counts.input_path is required"))
  if (!file.exists(input)) pe_fail("Count RDS does not exist: %s", input)
  x <- readRDS(input)
  required <- c("cellLine", "experiment", "plateID", "ploidy", "glucose", "hours")
  if (!all(required %in% names(x))) pe_fail("Existing count RDS lacks required metadata columns")
  if (!all(c("alive", "dead") %in% names(x)) && !all(c("N", "D") %in% names(x))) {
    pe_fail("Existing count RDS must contain alive/dead or N/D")
  }
  expected_sha <- as.character(counts_cfg$expected_sha256 %||% "")
  observed_sha <- pe_sha256(input)
  if (nzchar(expected_sha) && !identical(observed_sha, expected_sha)) {
    pe_fail("Count RDS SHA256 mismatch: expected %s, observed %s", expected_sha, observed_sha)
  }
  source_files <- input
} else if (mode == "image_status_csv") {
  images_path <- as.character(counts_cfg$images_csv %||% pe_fail("counts.images_csv is required"))
  validation_path <- as.character(counts_cfg$validation_json %||% "")
  if (!file.exists(images_path)) pe_fail("Image-status CSV does not exist: %s", images_path)
  if (nzchar(validation_path)) {
    validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
    if (!isTRUE(validation$validation_passed)) pe_fail("Image-processing validation did not pass: %s", validation_path)
    source_files <- c(source_files, validation_path)
  }
  images <- data.table::fread(images_path)
  required <- c("image_key", "processing_status", "n_alive", "n_dead")
  missing <- setdiff(required, names(images))
  if (length(missing)) pe_fail("Image-status CSV lacks: %s", paste(missing, collapse = ", "))
  if (anyDuplicated(images$image_key)) pe_fail("Image-status CSV contains duplicate image_key rows")
  bad <- images[!processing_status %in% c("ok", "zero_objects")]
  if (nrow(bad)) pe_fail("Image-status CSV contains %d failed/nonfinal fields", nrow(bad))
  zero_bad <- images[processing_status == "zero_objects" & (n_alive != 0 | n_dead != 0)]
  if (nrow(zero_bad)) pe_fail("zero_objects fields contain nonzero classified counts")

  source("scripts/define_plate_maps.R")
  meta <- data.table::rbindlist(lapply(images$image_key, function(key) {
    value <- get_meta(key)
    value$image_key <- key
    value
  }), fill = TRUE)
  joined <- merge(images, meta, by = "image_key", all.x = TRUE, sort = FALSE)
  metadata_cols <- c("cellLine", "experiment", "plateID", "ploidy", "glucose", "hours")
  if (any(!stats::complete.cases(joined[, ..metadata_cols]))) pe_fail("Some image keys could not be mapped to complete experiment metadata")
  field_audit <- joined[, .(
    n_fields = .N,
    n_zero_object_fields = sum(processing_status == "zero_objects"),
    alive = mean(as.numeric(n_alive)),
    dead = mean(as.numeric(n_dead))
  ), by = metadata_cols]
  x <- field_audit[, c(metadata_cols, "alive", "dead"), with = FALSE]
  source_files <- c(source_files, images_path)
} else {
  pe_fail("Unsupported counts.mode: %s", mode)
}

if (!nrow(x)) pe_fail("Count construction produced no rows")
artifact <- if (should_materialize) out_rds else as.character(counts_cfg$input_path)
if (should_materialize) {
  saveRDS(x, out_rds)
  utils::write.csv(as.data.frame(x), file.path(out_dir, "counts.csv"), row.names = FALSE, na = "")
  if (!is.null(field_audit)) utils::write.csv(as.data.frame(field_audit), file.path(out_dir, "field_audit.csv"), row.names = FALSE, na = "")
}

provenance <- list(
  artifact = artifact,
  sha256 = pe_sha256(artifact),
  mode = mode,
  materialized = should_materialize,
  count_definition = count_definition,
  zero_object_policy = if (mode == "image_status_csv") "included_as_zero" else "inherited_from_existing_rds",
  failed_field_policy = if (mode == "image_status_csv") "release_build_fails" else "unknown_in_legacy_input",
  source_files = lapply(source_files, function(path) list(path = path, sha256 = pe_sha256(path))),
  n_rows = nrow(x),
  generated_at = format(Sys.time(), usetz = TRUE),
  git_commit = pe_git_value(c("rev-parse", "HEAD"))
)
pe_write_json(provenance, file.path(out_dir, "provenance.json"))
cat(sprintf("%s %d count rows at %s\n", if (should_materialize) "Built" else "Registered", nrow(x), artifact))
