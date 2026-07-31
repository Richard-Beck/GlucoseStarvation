#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
source("scripts/parameter_estimation/common.R")
source("R/project_paths.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
source(get_model_r_path(cfg$model_name, cfg$model_version))
source("R/gpath_run_utils.R")
source("R/optim_utils.R")
source("R/gpath_derived_utils.R")

plan <- pe_read_tsv(file.path(root, "manifests", "fit_plan.tsv"), c(
  "fit_id", "dataset_id", "model_id", "stan_data_path", "optim_starts", "optim_dir"
))
datasets <- pe_read_tsv(file.path(root, "manifests", "datasets.tsv"), c(
  "dataset_id", "parent_dataset", "transform", "line_name", "direction"
))

stan_cache <- new.env(parent = emptyenv())
get_stan <- function(path) {
  if (!exists(path, envir = stan_cache, inherits = FALSE)) assign(path, readRDS(path), envir = stan_cache)
  get(path, envir = stan_cache, inherits = FALSE)
}

fit_rows <- vector("list", nrow(plan))
start_rows <- list()
start_ptr <- 1L
for (i in seq_len(nrow(plan))) {
  row <- plan[i, , drop = FALSE]
  required <- file.path(row$optim_dir[[1]], c("optim_draws_all.Rds", "optim_lp_all.Rds", "optim_rc_all.Rds"))
  complete <- all(file.exists(required))
  fit_row <- data.frame(
    fit_id = row$fit_id[[1]], dataset_id = row$dataset_id[[1]], model_id = row$model_id[[1]],
    complete = complete, n_starts_expected = as.integer(row$optim_starts[[1]]),
    n_starts = 0L, n_finite_lp = 0L, n_rc_zero = 0L,
    best_start = NA_integer_, best_lp = NA_real_, log_lik = NA_real_,
    log_lik_train = NA_real_, log_lik_holdout = NA_real_,
    k = gpd_model_k(row$model_id[[1]], get_stan(row$stan_data_path[[1]])),
    n_within_1 = NA_integer_, n_within_5 = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (complete) {
    lp <- as.numeric(readRDS(required[[2]]))
    rc <- as.numeric(unlist(readRDS(required[[3]]), use.names = FALSE))
    draws <- readRDS(required[[1]])
    if (length(lp) != length(draws) || length(rc) != length(lp)) pe_fail("Optimization length mismatch for %s", row$fit_id[[1]])
    finite <- is.finite(lp)
    delta <- if (any(finite)) max(lp[finite]) - lp else rep(NA_real_, length(lp))
    start_rows[[start_ptr]] <- data.frame(
      fit_id = row$fit_id[[1]], dataset_id = row$dataset_id[[1]], model_id = row$model_id[[1]],
      start_id = seq_along(lp), lp = lp, rc = rc, finite_lp = finite,
      delta_lp = delta, within_1 = finite & delta < 1, within_5 = finite & delta < 5,
      stringsAsFactors = FALSE
    )
    start_ptr <- start_ptr + 1L
    fit_row$n_starts <- length(lp)
    fit_row$n_finite_lp <- sum(finite)
    fit_row$n_rc_zero <- sum(rc == 0, na.rm = TRUE)
    fit_row$n_within_1 <- sum(finite & delta < 1)
    fit_row$n_within_5 <- sum(finite & delta < 5)
    if (any(finite)) {
      best <- which.max(replace(lp, !finite, -Inf))
      draw_vec <- extract_draw_vector(draws[[best]])
      metric <- function(name) if (name %in% names(draw_vec)) as.numeric(draw_vec[[name]]) else NA_real_
      fit_row$best_start <- best
      fit_row$best_lp <- lp[[best]]
      fit_row$log_lik <- metric("log_lik")
      fit_row$log_lik_train <- metric("log_lik_train")
      fit_row$log_lik_holdout <- metric("log_lik_holdout")
    }
  }
  fit_rows[[i]] <- fit_row
}

fit_summary <- rbindlist(fit_rows, fill = TRUE)
start_summary <- if (length(start_rows)) rbindlist(start_rows, fill = TRUE) else data.table()
fit_summary[, `:=`(
  deviance = ifelse(is.finite(best_lp), -2 * best_lp, NA_real_),
  AIC = ifelse(is.finite(best_lp), -2 * best_lp + 2 * k, NA_real_)
)]
fit_summary[, delta_AIC := if (any(is.finite(AIC))) AIC - min(AIC, na.rm = TRUE) else NA_real_, by = dataset_id]
fit_summary[, pareto_member := gpd_pareto_membership(k, deviance), by = dataset_id]
fit_summary[, pareto_rank := ifelse(pareto_member, frank(k, ties.method = "dense"), NA_integer_), by = dataset_id]

transfer_rows <- list()
ptr <- 1L
holdouts <- datasets[datasets$transform == "ploidy_holdout", , drop = FALSE]
nulls <- datasets[datasets$transform == "ploidy_null", , drop = FALSE]
for (i in seq_len(nrow(holdouts))) {
  h <- holdouts[i, , drop = FALSE]
  match_null <- nulls[
    nulls$parent_dataset == h$parent_dataset & nulls$line_name == h$line_name & nulls$direction == h$direction,
    , drop = FALSE
  ]
  if (nrow(match_null) != 1L) next
  real <- fit_summary[dataset_id == h$dataset_id]
  null <- fit_summary[dataset_id == match_null$dataset_id[[1]]]
  merged <- merge(real, null, by = "model_id", suffixes = c("_real", "_null"))
  if (!nrow(merged)) next
  transfer_rows[[ptr]] <- data.frame(
    dataset_id = h$dataset_id, null_dataset_id = match_null$dataset_id[[1]],
    parent_dataset = h$parent_dataset, line_name = h$line_name, direction = h$direction,
    model_id = merged$model_id,
    complete = merged$complete_real & merged$complete_null,
    real_log_lik_holdout = merged$log_lik_holdout_real,
    null_log_lik_holdout = merged$log_lik_holdout_null,
    transfer_minus_null = merged$log_lik_holdout_real - merged$log_lik_holdout_null,
    stringsAsFactors = FALSE
  )
  ptr <- ptr + 1L
}
transfer_vs_null <- if (length(transfer_rows)) rbindlist(transfer_rows, fill = TRUE) else data.table()
pareto <- fit_summary[pareto_member %in% TRUE & is.finite(deviance)]

artifact <- list(
  schema_version = 1L,
  metadata = list(release_id = cfg$release_id, generated_at = format(Sys.time(), usetz = TRUE)),
  fit_summary = as.data.frame(fit_summary),
  start_summary = as.data.frame(start_summary),
  pareto = as.data.frame(pareto),
  transfer_vs_null = as.data.frame(transfer_vs_null)
)
out <- file.path(derived_root, "optimization", "assessment.Rds")
gpd_atomic_save_rds(artifact, out, compress = "gzip")
cat(sprintf("Wrote %s (%d/%d complete optimization families)\n", out, sum(fit_summary$complete), nrow(fit_summary)))
