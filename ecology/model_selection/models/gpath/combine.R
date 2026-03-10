args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME      <- if (length(args) >= 1) args[1] else "gpath"
R_val           <- if (length(args) >= 2) as.integer(args[2]) else 2L
P_val           <- if (length(args) >= 3) as.integer(args[3]) else 2L
W_val           <- if (length(args) >= 4) as.integer(args[4]) else 0L
CONSTRAINT_FLAG <- if (length(args) >= 5) as.integer(args[5]) else 0L
WASTE_MECH_FLAG <- if (length(args) >= 6) as.integer(args[6]) else 1L
DELETE_AFTER    <- TRUE

HOLDOUT <- FALSE
tmp <- strsplit(MODEL_NAME, "_")[[1]]
if (length(tmp) > 1 && tmp[[2]] == "holdout") {
  HOLDOUT <- TRUE
  MODEL_NAME <- tmp[[1]]
}

run_id <- sprintf("%dR_%dP_%dW_C%d_M%d", R_val, P_val, W_val, CONSTRAINT_FLAG, WASTE_MECH_FLAG)
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, ifelse(HOLDOUT, "holdout", "hier"))

draw_files <- list.files(outDir, pattern = "^optim_draws_[0-9]+\\.Rds$", full.names = TRUE)
lp_files   <- list.files(outDir, pattern = "^optim_lp_[0-9]+\\.Rds$", full.names = TRUE)

get_id <- function(x, prefix) as.integer(sub(sprintf("^%s_([0-9]+)\\.Rds$", prefix), "\\1", basename(x)))

draw_ids <- get_id(draw_files, "optim_draws")
lp_ids   <- get_id(lp_files, "optim_lp")
ids <- sort(intersect(draw_ids, lp_ids))

if (!length(ids)) stop("No matching optim_draws_*.Rds / optim_lp_*.Rds files found in: ", outDir)

draw_map <- setNames(draw_files, draw_ids)
lp_map   <- setNames(lp_files, lp_ids)

opt_list <- lapply(ids, function(i) readRDS(draw_map[as.character(i)]))
lp_vec   <- vapply(ids, function(i) readRDS(lp_map[as.character(i)]), numeric(1))

saveRDS(opt_list, file.path(outDir, "optim_draws_all.Rds"))
saveRDS(lp_vec,   file.path(outDir, "optim_lp_all.Rds"))

if (DELETE_AFTER) {
  file.remove(unname(draw_map[as.character(ids)]), unname(lp_map[as.character(ids)]))
}

cat(sprintf(">>> Consolidation complete. Combined %d starts.\n", length(ids)))

if (length(setdiff(draw_ids, ids)) || length(setdiff(lp_ids, ids))) {
  cat(">>> Warning: unmatched files were ignored.\n")
}
