#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: build_segmentation_counts.R IMAGES_CSV OUTPUT_RDS", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

images_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- args[[2]]
if (file.exists(output_path)) {
  stop("Refusing to overwrite existing output: ", output_path, call. = FALSE)
}

plate_map_path <- normalizePath("scripts/define_plate_maps.R", mustWork = TRUE)
source(plate_map_path)

validation_path <- file.path(dirname(images_path), "validation.json")
if (!file.exists(validation_path)) {
  stop("Expected sibling validation file: ", validation_path, call. = FALSE)
}
validation <- jsonlite::fromJSON(validation_path, simplifyVector = TRUE)
if (!isTRUE(validation$validation_passed)) {
  stop("Segmentation validation did not pass: ", validation_path, call. = FALSE)
}

required_columns <- c(
  "image_key",
  "processing_status",
  "n_alive",
  "n_dead"
)
images <- fread(images_path)
missing_columns <- setdiff(required_columns, names(images))
if (length(missing_columns)) {
  stop("Image table lacks required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (!nrow(images)) {
  stop("Image table is empty.", call. = FALSE)
}
if (anyNA(images$image_key) || any(!nzchar(images$image_key)) || anyDuplicated(images$image_key)) {
  stop("Image keys must be present and unique.", call. = FALSE)
}
if (any(!images$processing_status %chin% c("ok", "zero_objects"))) {
  stop("Image table contains failed or nonfinal processing statuses.", call. = FALSE)
}

count_columns <- c("n_alive", "n_dead")
if (any(!stats::complete.cases(images[, ..count_columns])) ||
    any(unlist(images[, ..count_columns], use.names = FALSE) < 0)) {
  stop("Image counts contain missing or negative values.", call. = FALSE)
}
if (any(images$processing_status == "zero_objects" &
        (images$n_alive != 0 | images$n_dead != 0))) {
  stop("zero_objects rows contain nonzero object counts.", call. = FALSE)
}

message("Parsing metadata for ", nrow(images), " image keys")
metadata <- rbindlist(lapply(images$image_key, function(image_key) {
  value <- get_meta(image_key)
  value$image_key <- image_key
  value
}), fill = TRUE)

metadata_columns <- c("cellLine", "experiment", "plateID", "ploidy", "glucose", "hours")
if (nrow(metadata) != nrow(images) || anyDuplicated(metadata$image_key)) {
  stop("Metadata parsing did not produce one unique row per image key.", call. = FALSE)
}
if (any(!stats::complete.cases(metadata[, ..metadata_columns]))) {
  stop("Parsed image metadata contains missing values.", call. = FALSE)
}

joined <- merge(images, metadata, by = "image_key", all.x = TRUE, sort = FALSE)
if (any(!stats::complete.cases(joined[, ..metadata_columns]))) {
  stop("Metadata join left image rows unmatched.", call. = FALSE)
}

field_audit <- joined[, .(
  n_fields = .N,
  n_zero_object_fields = sum(processing_status == "zero_objects"),
  alive = mean(as.numeric(n_alive)),
  dead = mean(as.numeric(n_dead))
), by = metadata_columns]
counts <- field_audit[, c(metadata_columns, "alive", "dead"), with = FALSE]

if (sum(field_audit$n_fields) != nrow(images)) {
  stop("Grouped field count does not equal the image-table row count.", call. = FALSE)
}
if (!isTRUE(all.equal(sum(field_audit$alive * field_audit$n_fields),
                      sum(as.numeric(images$n_alive)))) ||
    !isTRUE(all.equal(sum(field_audit$dead * field_audit$n_fields),
                      sum(as.numeric(images$n_dead))))) {
  stop("Grouped field means do not reproduce image-table totals.", call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(counts, output_path)

round_trip <- readRDS(output_path)
required_output_columns <- c(metadata_columns, "alive", "dead")
if (!identical(names(round_trip), required_output_columns) ||
    !inherits(round_trip, "data.table") ||
    nrow(round_trip) != nrow(counts)) {
  stop("Output RDS failed schema or row-count validation.", call. = FALSE)
}

message("Wrote ", output_path)
message("  count rows: ", nrow(counts))
message("  image rows: ", nrow(images))
message("  count definition: arithmetic mean per imaged field (zero-object fields included)")
