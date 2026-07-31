#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: make_optim_fit_spec.R DATASETS.tsv MODELS.txt RESOURCES.tsv OUTPUT.tsv", call. = FALSE)
}

datasets <- utils::read.delim(args[[1]], sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
models <- trimws(readLines(args[[2]], warn = FALSE))
models <- models[nzchar(models) & !startsWith(models, "#")]
resources <- utils::read.delim(args[[3]], sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_datasets <- c("dataset_id", "transform")
required_resources <- c("transform", "optim_starts", "optim_threads", "optim_mem_gb", "optim_time", "optim_qos")
if (length(missing <- setdiff(required_datasets, names(datasets)))) stop("dataset spec lacks: ", paste(missing, collapse = ", "), call. = FALSE)
if (length(missing <- setdiff(required_resources, names(resources)))) stop("resource spec lacks: ", paste(missing, collapse = ", "), call. = FALSE)
if (anyDuplicated(datasets$dataset_id)) stop("duplicate dataset_id", call. = FALSE)
if (!length(models) || anyDuplicated(models)) stop("model list is empty or duplicated", call. = FALSE)
if (anyDuplicated(resources$transform)) stop("duplicate transform resource rule", call. = FALSE)

spec <- merge(datasets[c("dataset_id", "transform")], resources, by = "transform", all.x = TRUE, sort = FALSE)
if (any(!stats::complete.cases(spec[required_resources]))) stop("some dataset transforms lack resource rules", call. = FALSE)
spec <- spec[match(datasets$dataset_id, spec$dataset_id), , drop = FALSE]

rows <- lapply(seq_len(nrow(spec)), function(i) {
  data.frame(
    fit_id = paste0(spec$dataset_id[[i]], "__", tolower(models)),
    dataset_id = spec$dataset_id[[i]],
    model_id = models,
    run_optim = 1L,
    optim_starts = as.integer(spec$optim_starts[[i]]),
    optim_threads = as.integer(spec$optim_threads[[i]]),
    optim_mem_gb = as.integer(spec$optim_mem_gb[[i]]),
    optim_time = spec$optim_time[[i]],
    optim_qos = spec$optim_qos[[i]],
    run_nuts = 0L,
    nuts_chains = 4L,
    nuts_warmup = 500L,
    nuts_sampling = 1000L,
    adapt_delta = 0.99,
    max_treedepth = 12L,
    nuts_threads = 16L,
    nuts_mem_gb = 32L,
    nuts_time = "12:00:00",
    nuts_qos = "normal",
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
dir.create(dirname(args[[4]]), recursive = TRUE, showWarnings = FALSE)
utils::write.table(out, args[[4]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("Wrote %d dataset-model fits and %d starts to %s\n", nrow(out), sum(out$optim_starts), args[[4]]))
