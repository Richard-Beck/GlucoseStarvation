#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
param_plan <- pe_read_tsv(file.path(derived_root, "plans", "parameter_tasks.tsv"), c("dataset_id", "output_path"))
pred_plan <- pe_read_tsv(file.path(derived_root, "plans", "prediction_tasks.tsv"), c("dataset_id", "model_id", "output_path"))

records <- list()
errors <- character()
add_record <- function(kind, dataset_id, model_id, path, valid, detail = "") {
  records[[length(records) + 1L]] <<- data.frame(
    artifact_type = kind, dataset_id = dataset_id, model_id = model_id,
    path = path, exists = file.exists(path), valid = valid,
    bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
    sha256 = if (file.exists(path)) pe_sha256(path) else NA_character_,
    detail = detail, stringsAsFactors = FALSE
  )
  if (!valid) errors <<- c(errors, sprintf("%s: %s", path, detail))
}

assessment_path <- file.path(derived_root, "optimization", "assessment.Rds")
assessment_ok <- FALSE
assessment_detail <- "missing"
if (file.exists(assessment_path)) {
  x <- tryCatch(readRDS(assessment_path), error = identity)
  assessment_ok <- !inherits(x, "error") && all(c("fit_summary", "start_summary", "pareto", "transfer_vs_null") %in% names(x))
  assessment_detail <- if (assessment_ok) "schema ok" else "invalid assessment schema"
}
add_record("optimization_assessment", "", "", assessment_path, assessment_ok, assessment_detail)

qc_path <- file.path(derived_root, "posterior", "qc.Rds")
qc_ok <- FALSE
qc_detail <- "missing"
if (file.exists(qc_path)) {
  x <- tryCatch(readRDS(qc_path), error = identity)
  qc_ok <- !inherits(x, "error") && all(c("fit_summary", "chain_summary", "parameter_summary") %in% names(x))
  qc_detail <- if (qc_ok) "schema ok" else "invalid QC schema"
}
add_record("nuts_qc", "", "", qc_path, qc_ok, qc_detail)

for (i in seq_len(nrow(param_plan))) {
  path <- param_plan$output_path[[i]]
  ok <- FALSE; detail <- "missing"
  if (file.exists(path)) {
    x <- tryCatch(readRDS(path), error = identity)
    required <- c("dataset_id", "model_id", "draw_id", "chain_id", "iteration", "line_id", "ploidy", "parameter", "value")
    ok <- !inherits(x, "error") && is.list(x) && is.data.frame(x$draws) && all(required %in% names(x$draws))
    detail <- if (ok) sprintf("%d rows", nrow(x$draws)) else "invalid parameter schema"
  }
  add_record("posterior_parameters", param_plan$dataset_id[[i]], "", path, ok, detail)
}

for (i in seq_len(nrow(pred_plan))) {
  path <- pred_plan$output_path[[i]]
  ok <- FALSE; detail <- "missing"
  if (file.exists(path)) {
    x <- tryCatch(readRDS(path), error = identity)
    ok <- !inherits(x, "error") && is.list(x) &&
      all(c("draw_index", "well_metadata", "time", "state", "glucose_observation", "growth_surface") %in% names(x)) &&
      length(dim(x$state)) == 4L && length(dim(x$growth_surface$high_minus_low_growth)) == 4L &&
      dim(x$state)[[1]] == nrow(x$draw_index)
    detail <- if (ok) paste(dim(x$state), collapse = "x") else "invalid prediction schema"
  }
  add_record("posterior_predictions", pred_plan$dataset_id[[i]], pred_plan$model_id[[i]], path, ok, detail)
}

manifest <- do.call(rbind, records)
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
utils::write.table(manifest, file.path(derived_root, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
report <- list(
  validation_passed = !length(errors), release_id = cfg$release_id,
  checked_at = format(Sys.time(), usetz = TRUE), n_artifacts = nrow(manifest), errors = as.list(errors)
)
pe_write_json(report, file.path(derived_root, "validation.json"))
if (length(errors)) pe_fail("Derived validation failed for %d artifact(s); see %s", length(errors), file.path(derived_root, "manifest.tsv"))
cat(sprintf("PASS: %d derived artifacts\n", nrow(manifest)))
