#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
datasets <- pe_read_tsv(file.path(root, "manifests", "datasets.tsv"), c("dataset_id", "stan_data_path", "stan_data_sha256"))
fits <- pe_read_tsv(cfg$fits_tsv, c(
  "fit_id", "dataset_id", "model_id", "run_optim", "optim_starts", "optim_threads",
  "optim_mem_gb", "optim_time", "optim_qos", "run_nuts", "nuts_chains",
  "nuts_warmup", "nuts_sampling", "adapt_delta", "max_treedepth",
  "nuts_threads", "nuts_mem_gb", "nuts_time", "nuts_qos"
))
if (!nrow(fits)) pe_fail("Fit specification is empty")
if (anyDuplicated(fits$fit_id)) pe_fail("fit_id values must be unique")
invisible(lapply(fits$fit_id, pe_safe_id, label = "fit_id"))
unknown <- setdiff(fits$dataset_id, datasets$dataset_id)
if (length(unknown)) pe_fail("Fits refer to unknown datasets: %s", paste(unique(unknown), collapse = ", "))

plan_rows <- vector("list", nrow(fits))
optim_rows <- list()
nuts_rows <- list()
nuts_task <- 1L

for (i in seq_len(nrow(fits))) {
  row <- fits[i, , drop = FALSE]
  fit_id <- row$fit_id[[1]]
  dataset_id <- row$dataset_id[[1]]
  model_id <- row$model_id[[1]]
  run_optim <- pe_as_flag(row$run_optim[[1]], sprintf("%s.run_optim", fit_id))
  run_nuts <- pe_as_flag(row$run_nuts[[1]], sprintf("%s.run_nuts", fit_id))
  if (run_nuts == 1L && run_optim != 1L) pe_fail("%s requests NUTS without an optimization initializer", fit_id)
  drow <- datasets[datasets$dataset_id == dataset_id, , drop = FALSE]
  stan_path <- drow$stan_data_path[[1]]
  if (!file.exists(stan_path)) pe_fail("Missing Stan data for %s", fit_id)
  optim_dir <- file.path(root, "datasets", dataset_id, "optim", model_id, fit_id)
  nuts_dir <- file.path(root, "datasets", dataset_id, "nuts", model_id, fit_id)

  plan_rows[[i]] <- data.frame(
    fit_id = fit_id, dataset_id = dataset_id, model_id = model_id,
    stan_data_path = stan_path, stan_data_sha256 = drow$stan_data_sha256[[1]],
    run_optim = run_optim, optim_starts = as.integer(row$optim_starts[[1]]), optim_dir = optim_dir,
    run_nuts = run_nuts, nuts_chains = as.integer(row$nuts_chains[[1]]), nuts_dir = nuts_dir,
    stringsAsFactors = FALSE
  )

  if (run_optim == 1L) {
    optim_rows[[length(optim_rows) + 1L]] <- data.frame(
      fit_id = fit_id, dataset_id = dataset_id, model_id = model_id,
      n_starts = pe_as_positive_int(row$optim_starts[[1]], sprintf("%s.optim_starts", fit_id)),
      threads = pe_as_positive_int(row$optim_threads[[1]], sprintf("%s.optim_threads", fit_id)),
      stan_data_path = stan_path, output_dir = optim_dir,
      cpus = pe_as_positive_int(row$optim_threads[[1]], sprintf("%s.optim_threads", fit_id)),
      mem_gb = pe_as_positive_int(row$optim_mem_gb[[1]], sprintf("%s.optim_mem_gb", fit_id)),
      time = row$optim_time[[1]], qos = row$optim_qos[[1]], stringsAsFactors = FALSE
    )
  }

  if (run_nuts == 1L) {
    chains <- pe_as_positive_int(row$nuts_chains[[1]], sprintf("%s.nuts_chains", fit_id))
    for (chain_id in seq_len(chains)) {
      config_dir <- file.path(root, "manifests", "nuts_configs", fit_id)
      config_path <- file.path(config_dir, sprintf("chain_%02d.json", chain_id))
      run_tag <- sprintf("chain%02d", chain_id)
      task_cfg <- list(
        model_name = cfg$model_name,
        model_version = cfg$model_version,
        stan_data_path = stan_path,
        init = list(mode = "optim", source_path = optim_dir, seed = 8100000L + nuts_task),
        sampling = list(
          chains = 1L, parallel_chains = 1L,
          threads_per_chain = pe_as_positive_int(row$nuts_threads[[1]], sprintf("%s.nuts_threads", fit_id)),
          chain_id = chain_id, seed = 8200000L + nuts_task,
          iter_warmup = pe_as_positive_int(row$nuts_warmup[[1]], sprintf("%s.nuts_warmup", fit_id), allow_zero = TRUE),
          iter_sampling = pe_as_positive_int(row$nuts_sampling[[1]], sprintf("%s.nuts_sampling", fit_id)),
          adapt_delta = as.numeric(row$adapt_delta[[1]]),
          max_treedepth = pe_as_positive_int(row$max_treedepth[[1]], sprintf("%s.max_treedepth", fit_id)),
          save_warmup = TRUE, metric = "dense_e", refresh = 10L
        ),
        output = list(output_root = nuts_dir, output_dir = nuts_dir, run_tag = run_tag),
        model_specific_params = list(
          ploidy_effect_mask_spec = "all", run_id = model_id,
          posterior_dir = nuts_dir, optim_init_dir = optim_dir,
          auto_init_source_paths = list(optim_dir, nuts_dir)
        ),
        metadata = list(
          release_id = cfg$release_id, fit_id = fit_id, dataset_id = dataset_id,
          model_id = model_id, task_id = nuts_task, chain_id = chain_id,
          config_path_hint = config_path
        )
      )
      pe_write_json(task_cfg, config_path)
      nuts_rows[[length(nuts_rows) + 1L]] <- data.frame(
        task_id = nuts_task, fit_id = fit_id, run_id = model_id, chain_id = chain_id,
        config_path = config_path, output_dir = nuts_dir, run_tag = run_tag,
        cpus = as.integer(row$nuts_threads[[1]]), mem_gb = as.integer(row$nuts_mem_gb[[1]]),
        time = row$nuts_time[[1]], qos = row$nuts_qos[[1]], stringsAsFactors = FALSE
      )
      nuts_task <- nuts_task + 1L
    }
  }
}

plan <- do.call(rbind, plan_rows)
optim <- if (length(optim_rows)) do.call(rbind, optim_rows) else data.frame(
  fit_id = character(), dataset_id = character(), model_id = character(), n_starts = integer(),
  threads = integer(), stan_data_path = character(), output_dir = character(), cpus = integer(),
  mem_gb = integer(), time = character(), qos = character(), stringsAsFactors = FALSE
)
nuts <- if (length(nuts_rows)) do.call(rbind, nuts_rows) else data.frame(
  task_id = integer(), fit_id = character(), run_id = character(), chain_id = integer(),
  config_path = character(), output_dir = character(), run_tag = character(), cpus = integer(),
  mem_gb = integer(), time = character(), qos = character(), stringsAsFactors = FALSE
)
utils::write.table(plan, file.path(root, "manifests", "fit_plan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
utils::write.table(optim, file.path(root, "manifests", "optim_plan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
utils::write.table(nuts, file.path(root, "manifests", "nuts_plan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("Planned %d fits, %d optimization families and %d NUTS chain tasks\n", nrow(plan), nrow(optim), nrow(nuts)))
