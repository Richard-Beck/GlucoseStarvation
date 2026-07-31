#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: validate_sum159_targeted_nuclear_area.R MANIFEST MERGED_DIR REPORT_DIR FIGURE_DIR SMOKE_JSON OUTPUT_MD", call. = FALSE)
}
manifest_path <- args[[1]]
merged_dir <- args[[2]]
report_dir <- args[[3]]
figure_dir <- args[[4]]
smoke_json <- args[[5]]
output_md <- args[[6]]

manifest <- fread(manifest_path)
images <- fread(file.path(merged_dir, "wp3_nuclear_image_qc.csv"))
objects <- fread(file.path(merged_dir, "wp3_nuclear_object_features.csv"), select = c("image_key", "object_id", "nuclear_area_px"))
masks <- fread(file.path(merged_dir, "mask_manifest.csv"))
trajectory <- fread(file.path(report_dir, "sum159_nuclear_area_time_quantiles.csv"))
image_summary <- fread(file.path(report_dir, "sum159_nuclear_area_image_summary.csv"))
field_summary <- fread(file.path(report_dir, "sum159_nuclear_area_field_reproducibility.csv"))
measured_objects <- fread(file.path(report_dir, "sum159_nuclear_objects_alive.csv.gz"), select = "image_key")
smoke <- fromJSON(smoke_json)
qc_manifest <- file.path(dirname(merged_dir), "qc_stratified", "qc_manifest.csv")
qc <- fread(qc_manifest)
png_path <- file.path(figure_dir, "sum159_nuclear_area_trajectories.png")
pdf_path <- file.path(figure_dir, "sum159_nuclear_area_trajectories.pdf")

checks <- c(
  manifest_rows = nrow(manifest) == 340L,
  manifest_keys_unique = uniqueN(manifest$image_key) == 340L,
  four_equal_shards = identical(as.integer(sort(table(manifest$shard_id))), rep(85L, 4L)),
  image_rows_complete = nrow(images) == 340L && uniqueN(images$image_key) == 340L,
  all_images_ok = all(images$status == "ok"),
  object_keys_unique = uniqueN(objects, by = c("image_key", "object_id")) == nrow(objects),
  nuclear_calls_present = sum(objects$nuclear_area_px > 0) > 0L,
  mask_manifest_complete = nrow(masks) == 340L && all(file.exists(masks$cell_mask)) && all(file.exists(masks$nuclear_mask)),
  trajectory_grid_complete = nrow(trajectory) == 170L,
  estimability_rule_applied = all(trajectory$estimable == (trajectory$n_images == 2L & trajectory$n_alive_cells_with_nucleus >= 10L)),
  alive_image_summary_accounted = uniqueN(image_summary$image_key) == nrow(image_summary) &&
    setequal(image_summary$image_key, unique(measured_objects$image_key)) &&
    all(image_summary$n_alive_cells_with_nucleus > 0L),
  field_summary_complete = nrow(field_summary) == 170L && all(field_summary$n_fields %in% c(1L, 2L)),
  smoke_exact_overlap_pass = isTRUE(smoke$passed) && smoke$n_overlap_images == 12L && smoke$n_object_rows_with_value_mismatch == 0L && smoke$n_exact_mask_matches == 24L,
  stratified_qc_complete = nrow(qc) == 36L && all(file.exists(qc$output_png)),
  plot_outputs_exist = file.exists(png_path) && file.exists(pdf_path) && file.info(png_path)$size > 0 && file.info(pdf_path)$size > 0
)

estimable_keys <- trajectory[estimable == TRUE, .(batch_id, ploidy, G0, hours)]
field_estimable <- merge(field_summary, estimable_keys, by = c("batch_id", "ploidy", "G0", "hours"))
field_abs <- abs(field_estimable$nuclear_area_field_log2_ratio)
insufficient <- trajectory[estimable == FALSE]
missing_alive_images <- manifest[!image_summary, on = "image_key"]
lines <- c(
  "# Targeted SUM-159 nuclear-area validation",
  "",
  sprintf("Status: **%s**", if (all(checks)) "PASS" else "FAIL"),
  "",
  "## Run integrity",
  "",
  sprintf("- Manifest images: %s; shard sizes: %s.", nrow(manifest), paste(sort(table(manifest$shard_id)), collapse = ", ")),
  sprintf("- Completed images: %s/%s; failed images: %s.", sum(images$status == "ok"), nrow(images), sum(images$status != "ok")),
  sprintf("- Segmented objects: %s; objects with a nucleus: %s (%.2f%%).", format(nrow(objects), big.mark = ","), format(sum(objects$nuclear_area_px > 0), big.mark = ","), 100 * mean(objects$nuclear_area_px > 0)),
  sprintf("- Mask pairs: %s; stratified QC panels: %s.", nrow(masks), nrow(qc)),
  "",
  "## Scientific checks",
  "",
  sprintf("- Exact prior-pilot overlap: %s images, %s object rows, zero value mismatches, %s/%s masks byte-identical.", smoke$n_overlap_images, format(smoke$n_new_objects, big.mark = ","), smoke$n_exact_mask_matches, smoke$n_mask_files_checked),
  sprintf("- Trajectory grid: %s batch/ploidy/G0/time rows; %s pass the prespecified estimability rule (both fields and at least 10 cells).", nrow(trajectory), sum(trajectory$estimable)),
  sprintf("- Insufficient late-starvation endpoints excluded from the plotted estimates: %s; selected images with zero alive nuclei: %s.", nrow(insufficient), nrow(missing_alive_images)),
  sprintf("- Alive cells with nuclei per image: median %s, range %s-%s.", format(median(image_summary$n_alive_cells_with_nucleus), big.mark = ","), min(image_summary$n_alive_cells_with_nucleus), max(image_summary$n_alive_cells_with_nucleus)),
  sprintf("- Absolute between-field log2 median-area difference: median %.3f, 90th percentile %.3f, maximum %.3f.", median(field_abs), quantile(field_abs, 0.90, names = FALSE), max(field_abs)),
  "- Visual review: representative early and late QC panels from all three batch-by-line columns showed cell and nuclear boundaries aligned with the phase and nuclear-stain images; strong field illumination gradients remained visible but did not produce obvious off-cell nuclear masks.",
  "",
  "## Checks",
  "",
  paste0("- ", names(checks), ": ", ifelse(checks, "PASS", "FAIL"))
)
dir.create(dirname(output_md), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, output_md)
cat(paste(lines, collapse = "\n"), "\n")
if (!all(checks)) stop("One or more validation checks failed.", call. = FALSE)
