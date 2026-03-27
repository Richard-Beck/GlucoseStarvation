#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: make_prelim_counts_from_object_features.R <input_csv> [output_csv] [meta_json]", call. = FALSE)
}

input_csv <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
input_dir <- dirname(input_csv)
output_csv <- if (length(args) >= 2) normalizePath(args[[2]], winslash = "/", mustWork = FALSE) else file.path(input_dir, "prelim_counts_by_image.csv")
meta_json <- if (length(args) >= 3) normalizePath(args[[3]], winslash = "/", mustWork = FALSE) else file.path(input_dir, "prelim_counts_meta.json")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

repo_root <- tryCatch(normalizePath(".", winslash = "/", mustWork = TRUE), error = function(e) ".")
plate_map_script <- file.path(repo_root, "scripts", "define_plate_maps.R")
if (!file.exists(plate_map_script)) {
  stop("Could not find scripts/define_plate_maps.R from current working directory.", call. = FALSE)
}
source(plate_map_script)
if (!exists("get_meta", mode = "function")) {
  stop("get_meta() was not loaded from scripts/define_plate_maps.R", call. = FALSE)
}

dt <- data.table::fread(
  input_csv,
  select = c("image_key", "area_px", "live_bg_z", "dead_bg_z")
)

n_rows_total <- nrow(dt)

dt[, area_px := as.numeric(area_px)]
dt[, live_bg_z := as.numeric(live_bg_z)]
dt[, dead_bg_z := as.numeric(dead_bg_z)]

valid_mask <- is.finite(dt$area_px) & is.finite(dt$live_bg_z) & is.finite(dt$dead_bg_z) & dt$area_px > 0

dt_valid <- dt[valid_mask]
rm(dt)

n_rows_valid <- nrow(dt_valid)
if (n_rows_valid == 0) {
  stop("No valid rows found after applying basic validity checks.", call. = FALSE)
}

median_area_px <- as.numeric(stats::median(dt_valid$area_px))
area_min_px <- 0.05 * median_area_px

agg_all <- dt_valid[, .(
  total_count_valid_all_objects = .N,
  total_area_px_all_objects = sum(area_px)
), by = .(image_key)]

dt_kept <- dt_valid[area_px >= area_min_px]
dt_kept[, score := live_bg_z - dead_bg_z]

agg_kept <- dt_kept[, .(
  alive_count = sum(score >= 0),
  dead_count = sum(score < 0),
  total_count_kept = .N
), by = .(image_key)]

out <- merge(agg_all, agg_kept, by = "image_key", all.x = TRUE)
out[is.na(alive_count), alive_count := 0L]
out[is.na(dead_count), dead_count := 0L]
out[is.na(total_count_kept), total_count_kept := 0L]
out[, filtered_small_count := total_count_valid_all_objects - total_count_kept]

out[, time := sub("^.*_([0-9]{2}d[0-9]{2}h[0-9]{2}m)$", "\\1", image_key)]
out[, base_key := sub("_[0-9]{2}d[0-9]{2}h[0-9]{2}m$", "", image_key)]
out[, days := as.integer(substr(time, 1, 2))]
out[, hours_part := as.integer(substr(time, 4, 5))]
out[, mins := as.integer(substr(time, 7, 8))]
out[, hours := days * 24 + hours_part + mins / 60]
out[, c("days", "hours_part", "mins") := NULL]

meta_input <- paste0(unique(out$base_key), "_00d00h00m")
meta_df <- data.table::as.data.table(do.call(rbind, lapply(meta_input, get_meta)))
meta_df[, base_key := unique(out$base_key)]
meta_df <- meta_df[, .(base_key, cellLine, experiment, plateID, pos, ploidy, glucose)]
out <- merge(out, meta_df, by = "base_key", all.x = TRUE)
out[, G0 := glucose]

data.table::setcolorder(
  out,
  c(
    "cellLine",
    "experiment",
    "plateID",
    "pos",
    "ploidy",
    "glucose",
    "G0",
    "hours",
    "image_key",
    "base_key",
    "time",
    "alive_count",
    "dead_count",
    "total_count_kept",
    "total_count_valid_all_objects",
    "filtered_small_count",
    "total_area_px_all_objects"
  )
)

data.table::setorder(out, image_key)
data.table::fwrite(out, output_csv)

sha_line <- tryCatch(
  system2("sha256sum", input_csv, stdout = TRUE, stderr = TRUE),
  error = function(e) character(0)
)
input_csv_sha256 <- if (length(sha_line) > 0 && nzchar(sha_line[1])) {
  strsplit(sha_line[1], "\\s+")[[1]][1]
} else {
  NA_character_
}

git_commit <- tryCatch(
  {
    x <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)
    if (length(x) > 0 && nzchar(x[1])) x[1] else NA_character_
  },
  error = function(e) NA_character_
)

meta <- list(
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_csv = input_csv,
  input_csv_sha256 = input_csv_sha256,
  git_commit = git_commit,
  rules = list(
    valid_object = "finite(area_px, live_bg_z, dead_bg_z) and area_px > 0",
    size_filter = "area_px >= 0.05 * median_area_px_over_valid_objects",
    classification = "alive if (live_bg_z - dead_bg_z) >= 0 else dead",
    confluence_proxy = "sum(area_px) over all valid objects per image"
  ),
  summary = list(
    n_rows_total = as.integer(n_rows_total),
    n_rows_valid = as.integer(n_rows_valid),
    median_area_px = median_area_px,
    area_min_px = area_min_px,
    n_images_out = as.integer(nrow(out)),
    alive_count_total = as.integer(sum(out$alive_count)),
    dead_count_total = as.integer(sum(out$dead_count)),
    kept_count_total = as.integer(sum(out$total_count_kept))
  ),
  outputs = list(
    prelim_counts_by_image_csv = normalizePath(output_csv, winslash = "/", mustWork = FALSE),
    prelim_counts_meta_json = normalizePath(meta_json, winslash = "/", mustWork = FALSE)
  )
)

jsonlite::write_json(meta, meta_json, auto_unbox = TRUE, pretty = TRUE)

cat("Wrote:", output_csv, "\n")
cat("Wrote:", meta_json, "\n")
