#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: build_cell_nuclear_area_summaries.R OBJECTS_CSV OUTPUT_RDS",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(data.table))

objects_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- args[[2]]
if (file.exists(output_path)) {
  stop("Refusing to overwrite existing output: ", output_path, call. = FALSE)
}

plate_map_path <- normalizePath("scripts/define_plate_maps.R", mustWork = TRUE)
source(plate_map_path)

quantile_probs <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)
quantile_labels <- c("q0", "q01", "q05", "q10", "q25", "q50", "q75", "q90", "q95", "q99", "q100")

required_columns <- c(
  "image_key",
  "predicted_label_name",
  "segmented_area_px",
  "cell_area_pass",
  "nuclear_area_px",
  "nuclear_to_segmented_area_ratio"
)

named_quantiles <- function(x, prefix) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  values <- if (length(x)) {
    as.numeric(stats::quantile(x, probs = quantile_probs, names = FALSE, type = 7))
  } else {
    rep(NA_real_, length(quantile_probs))
  }
  stats::setNames(as.list(values), paste0(prefix, "_", quantile_labels))
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, call. = FALSE)
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_args)) {
  normalizePath(sub("^--file=", "", script_args[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

message("Reading selected object columns from ", objects_path)
objects <- fread(objects_path, select = required_columns, showProgress = TRUE)
missing_columns <- setdiff(required_columns, names(objects))
if (length(missing_columns)) {
  stop("Object table lacks required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (!nrow(objects)) {
  stop("Object table is empty.", call. = FALSE)
}

objects[, predicted_label_name := as.character(predicted_label_name)]
objects[, cell_area_pass := as.integer(cell_area_pass)]
numeric_columns <- c("segmented_area_px", "nuclear_area_px", "nuclear_to_segmented_area_ratio")
objects[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]

if (anyNA(objects$image_key) || any(!nzchar(objects$image_key))) {
  stop("Object table contains missing image keys.", call. = FALSE)
}
if (any(!objects$predicted_label_name %chin% c("alive", "dead", "junk"))) {
  stop("Object table contains unexpected predicted labels.", call. = FALSE)
}
if (any(!objects$cell_area_pass %in% c(0L, 1L))) {
  stop("cell_area_pass contains values other than 0 and 1.", call. = FALSE)
}
if (any(!is.finite(objects$segmented_area_px) | objects$segmented_area_px <= 0)) {
  stop("segmented_area_px contains nonpositive or nonfinite values.", call. = FALSE)
}
if (any(!is.finite(objects$nuclear_area_px) | objects$nuclear_area_px < 0)) {
  stop("nuclear_area_px contains negative or nonfinite values.", call. = FALSE)
}
area_exceeds <- which(objects$nuclear_area_px > objects$segmented_area_px)
if (length(area_exceeds)) {
  message(
    sprintf(
      "%d nuclear areas exceed segmented area; retaining without clipping (maximum ratio %.4f)",
      length(area_exceeds),
      max(objects$nuclear_area_px[area_exceeds] / objects$segmented_area_px[area_exceeds])
    )
  )
}
if (any(objects$cell_area_pass != as.integer(objects$segmented_area_px >= 50))) {
  stop("cell_area_pass is inconsistent with the configured 50-pixel threshold.", call. = FALSE)
}

image_keys <- unique(objects$image_key)
message("Parsing metadata for ", length(image_keys), " image keys")
metadata <- rbindlist(lapply(image_keys, function(image_key) {
  row <- get_meta(image_key)
  row$image_key <- image_key
  row
}), use.names = TRUE, fill = TRUE)

metadata_columns <- c("cellLine", "experiment", "plateID", "pos", "ploidy", "glucose", "hours")
if (nrow(metadata) != length(image_keys) || anyDuplicated(metadata$image_key)) {
  stop("Image metadata parsing did not produce one unique row per image key.", call. = FALSE)
}
if (any(!stats::complete.cases(metadata[, ..metadata_columns]))) {
  stop("Image metadata contains missing values.", call. = FALSE)
}

setkey(metadata, image_key)
objects[metadata, on = "image_key", `:=`(
  cellLine = i.cellLine,
  experiment = i.experiment,
  plateID = i.plateID,
  pos = i.pos,
  ploidy = i.ploidy,
  glucose = i.glucose,
  hours = i.hours
)]
if (any(!stats::complete.cases(objects[, ..metadata_columns]))) {
  stop("Metadata join left object rows unmatched.", call. = FALSE)
}

group_columns <- c("cellLine", "experiment", "plateID", "ploidy", "glucose", "hours")

summarize_replicate_time <- function(x, scope) {
  x[, {
    detected <- is.finite(nuclear_area_px) & nuclear_area_px > 0
    c(
      list(
        object_scope = scope,
        n_fields = uniqueN(image_key),
        n_objects = .N,
        n_cell_area_pass = sum(cell_area_pass == 1L),
        cell_area_pass_fraction = mean(cell_area_pass == 1L),
        n_with_nucleus = sum(detected),
        nuclear_detection_fraction = mean(detected),
        n_nuclear_area_exceeds_segmented = sum(nuclear_area_px > segmented_area_px),
        nuclear_area_exceeds_segmented_fraction = mean(nuclear_area_px > segmented_area_px),
        segmented_area_px_sum = sum(segmented_area_px),
        nuclear_area_px_sum = sum(nuclear_area_px)
      ),
      named_quantiles(segmented_area_px, "segmented_area_px"),
      named_quantiles(nuclear_area_px[detected], "nuclear_area_detected_px"),
      named_quantiles(
        nuclear_to_segmented_area_ratio[detected],
        "nuclear_to_segmented_area_ratio_detected"
      )
    )
  }, by = group_columns]
}

summarize_field_qc <- function(x, scope) {
  x[, {
    detected <- is.finite(nuclear_area_px) & nuclear_area_px > 0
    list(
      object_scope = scope,
      n_objects = .N,
      n_cell_area_pass = sum(cell_area_pass == 1L),
      cell_area_pass_fraction = mean(cell_area_pass == 1L),
      n_with_nucleus = sum(detected),
      nuclear_detection_fraction = mean(detected),
      n_nuclear_area_exceeds_segmented = sum(nuclear_area_px > segmented_area_px),
      nuclear_area_exceeds_segmented_fraction = mean(nuclear_area_px > segmented_area_px),
      segmented_area_px_sum = sum(segmented_area_px),
      nuclear_area_px_sum = sum(nuclear_area_px)
    )
  }, by = "image_key"]
}

scopes <- list(
  all = objects,
  alive = objects[predicted_label_name == "alive"],
  alive_area_pass = objects[predicted_label_name == "alive" & cell_area_pass == 1L]
)

message("Computing replicate-by-time summaries")
replicate_time <- rbindlist(
  Map(summarize_replicate_time, scopes, names(scopes)),
  use.names = TRUE,
  fill = TRUE
)
setcolorder(replicate_time, c(group_columns, "object_scope", setdiff(names(replicate_time), c(group_columns, "object_scope"))))
setorder(replicate_time, cellLine, experiment, plateID, ploidy, glucose, hours, object_scope)

message("Computing field-level QC summaries")
field_qc <- rbindlist(
  Map(summarize_field_qc, scopes, names(scopes)),
  use.names = TRUE,
  fill = TRUE
)
field_qc[metadata, on = "image_key", `:=`(
  cellLine = i.cellLine,
  experiment = i.experiment,
  plateID = i.plateID,
  pos = i.pos,
  ploidy = i.ploidy,
  glucose = i.glucose,
  hours = i.hours
)]
setcolorder(field_qc, c("image_key", metadata_columns, "object_scope", setdiff(names(field_qc), c("image_key", metadata_columns, "object_scope"))))
setorder(field_qc, cellLine, experiment, plateID, ploidy, glucose, hours, pos, object_scope)

git_commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]],
  error = function(e) NA_character_
)

settings <- list(
  schema_version = 1L,
  generated_at = format(Sys.time(), usetz = TRUE),
  quantile_probs = quantile_probs,
  quantile_labels = quantile_labels,
  replicate_keys = group_columns,
  object_scopes = list(
    all = "all segmented and classified objects",
    alive = "objects classified as alive",
    alive_area_pass = "objects classified as alive with segmented_area_px >= 50"
  ),
  cell_area_pass_definition = "segmented_area_px >= 50; inherited from the segmentation run",
  nuclear_detection_definition = "nuclear_area_px > 0",
  nuclear_quantile_population = "objects with nuclear_area_px > 0 within each object scope",
  nuclear_area_overshoot_policy = paste0(
    "Retained without clipping and counted in QC. Nuclear-mask morphological closing/fill was not ",
    "re-intersected with the cell mask in the source segmentation implementation."
  ),
  n_nuclear_area_exceeds_segmented = length(area_exceeds),
  aggregation_definition = "objects pooled across image fields within each experimental replicate and timepoint",
  provenance = list(
    objects_path = objects_path,
    objects_sha256 = sha256_file(objects_path),
    objects_n_rows = nrow(objects),
    objects_n_images = length(image_keys),
    plate_map_path = plate_map_path,
    plate_map_sha256 = sha256_file(plate_map_path),
    builder_path = script_path,
    builder_sha256 = if (!is.na(script_path)) sha256_file(script_path) else NA_character_,
    git_commit = git_commit
  )
)

result <- list(
  replicate_time = as.data.frame(replicate_time),
  field_qc = as.data.frame(field_qc),
  settings = settings
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_path, compress = "xz")

round_trip <- readRDS(output_path)
if (!identical(names(round_trip), c("replicate_time", "field_qc", "settings"))) {
  stop("Output RDS failed list-schema validation.", call. = FALSE)
}
if (nrow(round_trip$replicate_time) != nrow(replicate_time) || nrow(round_trip$field_qc) != nrow(field_qc)) {
  stop("Output RDS failed row-count validation.", call. = FALSE)
}

message("Wrote ", output_path)
message("  replicate_time rows: ", nrow(replicate_time))
message("  field_qc rows: ", nrow(field_qc))
message("  source object rows: ", nrow(objects))
