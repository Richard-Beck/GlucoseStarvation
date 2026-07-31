#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")
cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
plan <- pe_read_tsv(file.path(root, "manifests", "fit_plan.tsv"))

rows <- lapply(seq_len(nrow(plan)), function(i) {
  fit <- plan[i, , drop = FALSE]
  shard_ids <- if (dir.exists(fit$optim_dir[[1]])) {
    files <- list.files(fit$optim_dir[[1]], pattern = "^optim_draws_[0-9]+\\.Rds$")
    as.integer(sub("^optim_draws_([0-9]+)\\.Rds$", "\\1", files))
  } else integer()
  nuts_draws <- if (dir.exists(fit$nuts_dir[[1]])) {
    length(list.files(fit$nuts_dir[[1]], pattern = "^nuts_draws_chain[0-9]+\\.Rds$"))
  } else 0L
  data.frame(
    fit_id = fit$fit_id, dataset_id = fit$dataset_id, model_id = fit$model_id,
    optim_shards = length(unique(shard_ids)),
    optim_expected = if (fit$run_optim[[1]] == 1L) fit$optim_starts else 0L,
    optim_combined = as.integer(file.exists(file.path(fit$optim_dir, "optim_draws_all.Rds"))),
    nuts_chains = nuts_draws,
    nuts_expected = if (fit$run_nuts[[1]] == 1L) fit$nuts_chains else 0L,
    stringsAsFactors = FALSE
  )
})
status <- do.call(rbind, rows)
out <- file.path(root, "status.tsv")
utils::write.table(status, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(status, row.names = FALSE)
cat(sprintf("status_file=%s\n", out))
