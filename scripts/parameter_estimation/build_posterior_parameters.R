#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: build_posterior_parameters.R CONFIG TASK_ID", call. = FALSE)
suppressPackageStartupMessages({library(data.table); library(posterior)})
source("scripts/parameter_estimation/common.R")
source("R/project_paths.R")

cfg <- pe_read_config(args[[1]])
task_id <- as.integer(args[[2]])
root <- pe_release_root(cfg)
derived_root <- pe_derived_root(cfg, root)
tasks <- pe_read_tsv(file.path(derived_root, "plans", "parameter_tasks.tsv"), c(
  "task_id", "dataset_id", "stan_data_path", "output_path"
))
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) pe_fail("Unknown parameter task id: %s", args[[2]])
if (file.exists(task$output_path[[1]])) pe_fail("Refusing to overwrite parameter shard: %s", task$output_path[[1]])

source(get_model_r_path(cfg$model_name, cfg$model_version))
source("R/gpath_run_utils.R")
source("R/parameter_transfer_utils.R")
source("R/gpath_derived_utils.R")

fits <- pe_read_tsv(file.path(root, "manifests", "fit_plan.tsv"), c(
  "fit_id", "dataset_id", "model_id", "run_nuts", "nuts_chains", "nuts_dir"
))
fits <- fits[fits$dataset_id == task$dataset_id[[1]] & as.integer(fits$run_nuts) == 1L, , drop = FALSE]
if (!nrow(fits)) pe_fail("No NUTS fits for dataset %s", task$dataset_id[[1]])
stan_data <- readRDS(task$stan_data_path[[1]])
cores <- max(1L, min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), nrow(fits)))

build_one <- function(i) {
  fit <- fits[i, , drop = FALSE]
  draw_info <- gpd_read_fit_draws(fit$nuts_dir[[1]], fit$nuts_chains[[1]], max_draws = Inf)
  gpd_reconstruct_parameter_table(draw_info, stan_data, task$dataset_id[[1]], fit$model_id[[1]])
}
parts <- if (cores > 1L) parallel::mclapply(seq_len(nrow(fits)), build_one, mc.cores = cores) else lapply(seq_len(nrow(fits)), build_one)
draws <- rbindlist(parts, use.names = TRUE, fill = TRUE)
chain_paths <- unlist(lapply(seq_len(nrow(fits)), function(i) {
  gpd_expected_chain_paths(fits$nuts_dir[[i]], fits$nuts_chains[[i]])
}), use.names = FALSE)

artifact <- list(
  schema_version = 1L,
  metadata = list(
    release_id = cfg$release_id,
    dataset_id = task$dataset_id[[1]],
    stan_data_path = task$stan_data_path[[1]],
    stan_data_sha256 = pe_sha256(task$stan_data_path[[1]]),
    models = as.list(fits$model_id),
    source_chain_sha256 = setNames(lapply(chain_paths, pe_sha256), chain_paths),
    draw_policy = "all sampling draws",
    generated_at = format(Sys.time(), usetz = TRUE)
  ),
  draws = as.data.frame(draws)
)
gpd_atomic_save_rds(artifact, task$output_path[[1]], compress = "gzip")
cat(sprintf("Wrote %s (%d reconstructed parameter rows)\n", task$output_path[[1]], nrow(draws)))
