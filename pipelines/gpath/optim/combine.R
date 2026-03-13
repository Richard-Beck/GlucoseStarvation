#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_DIR <- if (length(args) >= 1) args[1] else stop("output_dir required")
DELETE_AFTER <- if (length(args) >= 2) as.logical(as.integer(args[2])) else TRUE

cat(sprintf(">>> combine.R cwd: %s\n", getwd()))
cat(sprintf(">>> Combining optimisation starts in %s\n", OUTPUT_DIR))

draw_files <- list.files(OUTPUT_DIR, pattern = "^optim_draws_[0-9]+\\.Rds$", full.names = TRUE)
lp_files <- list.files(OUTPUT_DIR, pattern = "^optim_lp_[0-9]+\\.Rds$", full.names = TRUE)
rc_files <- list.files(OUTPUT_DIR, pattern = "^optim_rc_[0-9]+\\.Rds$", full.names = TRUE)

get_id <- function(x, prefix) as.integer(sub(sprintf("^%s_([0-9]+)\\.Rds$", prefix), "\\1", basename(x)))

draw_ids <- get_id(draw_files, "optim_draws")
lp_ids <- get_id(lp_files, "optim_lp")
rc_ids <- get_id(rc_files, "optim_rc")
ids <- sort(Reduce(intersect, list(draw_ids, lp_ids, rc_ids)))

if (!length(ids)) {
  stop(sprintf("No matching optim_draws/optim_lp/optim_rc files found in %s", OUTPUT_DIR))
}

draw_map <- setNames(draw_files, draw_ids)
lp_map <- setNames(lp_files, lp_ids)
rc_map <- setNames(rc_files, rc_ids)

opt_list <- lapply(ids, function(i) readRDS(draw_map[as.character(i)]))
lp_vec <- vapply(ids, function(i) readRDS(lp_map[as.character(i)]), numeric(1))
rc_vec <- vapply(ids, function(i) readRDS(rc_map[as.character(i)]), numeric(1))

saveRDS(opt_list, file.path(OUTPUT_DIR, "optim_draws_all.Rds"))
saveRDS(lp_vec, file.path(OUTPUT_DIR, "optim_lp_all.Rds"))
saveRDS(rc_vec, file.path(OUTPUT_DIR, "optim_rc_all.Rds"))

if (DELETE_AFTER) {
  file.remove(unname(draw_map[as.character(ids)]), unname(lp_map[as.character(ids)]), unname(rc_map[as.character(ids)]))
}

cat(sprintf(">>> Consolidation complete. Combined %d starts.\n", length(ids)))
