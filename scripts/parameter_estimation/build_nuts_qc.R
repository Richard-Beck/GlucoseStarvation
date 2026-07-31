#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
suppressPackageStartupMessages(library(posterior))
source("scripts/parameter_estimation/common.R")
source("R/gpath_derived_utils.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
fits <- pe_read_tsv(file.path(root, "manifests", "fit_plan.tsv"), c(
  "fit_id", "dataset_id", "model_id", "run_nuts", "nuts_chains", "nuts_dir"
))
fits <- fits[as.integer(fits$run_nuts) == 1L, , drop = FALSE]

fit_rows <- list()
chain_rows <- list()
parameter_rows <- list()
for (i in seq_len(nrow(fits))) {
  fit <- fits[i, , drop = FALSE]
  draw_paths <- gpd_expected_chain_paths(fit$nuts_dir[[1]], fit$nuts_chains[[1]])
  diag_paths <- gpd_expected_diagnostic_paths(fit$nuts_dir[[1]], fit$nuts_chains[[1]])
  gpd_require_files(c(draw_paths, diag_paths), sprintf("NUTS outputs for %s", fit$fit_id[[1]]))
  bound <- gpd_bind_fit_draws(fit$nuts_dir[[1]], fit$nuts_chains[[1]])
  summary <- as.data.frame(posterior::summarise_draws(
    bound, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail
  ))
  names(summary)[names(summary) == "variable"] <- "parameter"
  parameter_rows[[i]] <- cbind(
    data.frame(fit_id = fit$fit_id[[1]], dataset_id = fit$dataset_id[[1]], model_id = fit$model_id[[1]], stringsAsFactors = FALSE),
    summary
  )
  per_chain <- lapply(seq_along(diag_paths), function(chain_id) {
    x <- posterior::as_draws_matrix(readRDS(diag_paths[[chain_id]]))
    get <- function(name) if (name %in% colnames(x)) as.numeric(x[, name]) else numeric()
    energy <- get("energy__")
    ebfmi <- if (length(energy) > 1L && stats::var(energy) > 0) mean(diff(energy)^2) / stats::var(energy) else NA_real_
    data.frame(
      fit_id = fit$fit_id[[1]], dataset_id = fit$dataset_id[[1]], model_id = fit$model_id[[1]], chain_id = chain_id,
      n_draws = nrow(x), divergences = sum(get("divergent__") > 0),
      max_treedepth_observed = if (length(get("treedepth__"))) max(get("treedepth__")) else NA_real_,
      ebfmi = ebfmi, stringsAsFactors = FALSE
    )
  })
  chain_df <- do.call(rbind, per_chain)
  chain_rows[[i]] <- chain_df
  fit_rows[[i]] <- data.frame(
    fit_id = fit$fit_id[[1]], dataset_id = fit$dataset_id[[1]], model_id = fit$model_id[[1]],
    n_chains = nrow(chain_df), n_draws = sum(chain_df$n_draws),
    divergences = sum(chain_df$divergences), max_treedepth_observed = max(chain_df$max_treedepth_observed, na.rm = TRUE),
    min_ebfmi = min(chain_df$ebfmi, na.rm = TRUE), max_rhat = max(summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE), min_ess_tail = min(summary$ess_tail, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

artifact <- list(
  schema_version = 1L,
  metadata = list(release_id = cfg$release_id, generated_at = format(Sys.time(), usetz = TRUE)),
  fit_summary = do.call(rbind, fit_rows),
  chain_summary = do.call(rbind, chain_rows),
  parameter_summary = do.call(rbind, parameter_rows)
)
out <- file.path(derived_root, "posterior", "qc.Rds")
gpd_atomic_save_rds(artifact, out, compress = "gzip")
cat(sprintf("Wrote %s (%d fits)\n", out, nrow(artifact$fit_summary)))
