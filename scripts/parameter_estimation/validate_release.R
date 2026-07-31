#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")

cfg <- pe_read_config(args[[1]])
phase <- if (length(args) >= 2L) tolower(args[[2]]) else "prepare"
if (!phase %in% c("prepare", "complete")) pe_fail("phase must be prepare or complete")
root <- pe_release_root(cfg)
errors <- character()
checks <- list()

check <- function(label, condition, detail = "") {
  ok <- isTRUE(condition)
  checks[[length(checks) + 1L]] <<- list(label = label, status = if (ok) "PASS" else "FAIL", detail = detail)
  if (!ok) errors <<- c(errors, sprintf("%s%s", label, if (nzchar(detail)) paste0(": ", detail) else ""))
  invisible(ok)
}

check("release root exists", dir.exists(root), root)
builder <- "R/prepare_data.R"
check("Stan-data builder exists", file.exists(builder), builder)
if (file.exists(builder)) {
  age_days <- as.numeric(difftime(Sys.time(), file.info(builder)$mtime, units = "days"))
  check("Stan-data builder modified within 14 days", is.finite(age_days) && age_days <= 14, sprintf("mtime=%s", file.info(builder)$mtime))
  builder_text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  check(
    "Stan-data builder contains case-normalization fix",
    grepl("normalize_glucose_ploidy_column", builder_text, fixed = TRUE),
    "guards Ploidy/ploidy input variants"
  )
}

counts_path <- pe_counts_path(cfg, root)
check("configured counts exist", file.exists(counts_path), counts_path)
expected_counts_sha <- as.character(cfg$counts$expected_sha256 %||% "")
if (file.exists(counts_path) && nzchar(expected_counts_sha)) {
  check(
    "configured counts SHA256 matches",
    identical(pe_sha256(counts_path), expected_counts_sha),
    sprintf("expected=%s observed=%s", expected_counts_sha, pe_sha256(counts_path))
  )
}

dataset_manifest_path <- file.path(root, "manifests", "datasets.tsv")
fit_plan_path <- file.path(root, "manifests", "fit_plan.tsv")
optim_plan_path <- file.path(root, "manifests", "optim_plan.tsv")
nuts_plan_path <- file.path(root, "manifests", "nuts_plan.tsv")
check("dataset manifest exists", file.exists(dataset_manifest_path), dataset_manifest_path)
check("fit plan exists", file.exists(fit_plan_path), fit_plan_path)
check("optimization plan exists", file.exists(optim_plan_path), optim_plan_path)
check("NUTS plan exists", file.exists(nuts_plan_path), nuts_plan_path)

datasets <- NULL
if (file.exists(dataset_manifest_path)) {
  datasets <- tryCatch(pe_read_tsv(dataset_manifest_path, c(
    "dataset_id", "parent_dataset", "transform", "line_name", "stan_data_path", "stan_data_sha256",
    "N_lines", "N_wells", "N_obs_count", "N_obs_gluc"
  )), error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
}

validate_stan <- function(path, expected_sha = NULL) {
  if (!file.exists(path)) return(sprintf("missing %s", path))
  x <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(x, "error")) return(conditionMessage(x))
  required <- c("N_lines", "N_wells", "N_obs_count", "N_obs_gluc", "line_id",
                "well_idx_count", "well_idx_gluc", "N_obs", "D_obs")
  missing <- setdiff(required, names(x))
  if (length(missing)) return(sprintf("missing fields: %s", paste(missing, collapse = ", ")))
  problems <- character()
  if (length(x$line_id) != x$N_wells) problems <- c(problems, "line_id/N_wells mismatch")
  if (length(x$well_idx_count) != x$N_obs_count || length(x$N_obs) != x$N_obs_count || length(x$D_obs) != x$N_obs_count) {
    problems <- c(problems, "count-vector length mismatch")
  }
  if (length(x$well_idx_gluc) != x$N_obs_gluc) problems <- c(problems, "glucose-vector length mismatch")
  if (length(x$well_idx_count) && any(x$well_idx_count < 1L | x$well_idx_count > x$N_wells)) problems <- c(problems, "invalid count well index")
  if (length(x$well_idx_gluc) && any(x$well_idx_gluc < 1L | x$well_idx_gluc > x$N_wells)) problems <- c(problems, "invalid glucose well index")
  if (!is.null(expected_sha) && !identical(pe_sha256(path), expected_sha)) problems <- c(problems, "SHA256 mismatch")
  paste(problems, collapse = "; ")
}

if (!is.null(datasets)) {
  check("dataset ids are unique", !anyDuplicated(datasets$dataset_id))
  for (i in seq_len(nrow(datasets))) {
    detail <- validate_stan(datasets$stan_data_path[[i]], datasets$stan_data_sha256[[i]])
    check(sprintf("dataset %s is internally valid", datasets$dataset_id[[i]]), !nzchar(detail), detail)
    transform <- datasets$transform[[i]]
    if (transform != "identity" && file.exists(datasets$stan_data_path[[i]])) {
      x <- readRDS(datasets$stan_data_path[[i]])
      parent_row <- datasets[datasets$dataset_id == datasets$parent_dataset[[i]], , drop = FALSE]
      semantic <- character()
      if (nrow(parent_row) != 1L || !file.exists(parent_row$stan_data_path[[1]])) {
        semantic <- c(semantic, "missing parent dataset")
      } else {
        parent <- readRDS(parent_row$stan_data_path[[1]])
        line_name <- datasets$line_name[[i]]
        if (transform == "exclude_line") {
          if (x$N_lines != parent$N_lines - 1L || line_name %in% names(x$line_map)) semantic <- c(semantic, "line exclusion mismatch")
        } else if (transform == "single_line") {
          if (x$N_lines != 1L || !identical(names(x$line_map), line_name)) semantic <- c(semantic, "single-line mapping mismatch")
        } else if (transform %in% c("ploidy_holdout", "ploidy_null")) {
          if (x$N_wells != parent$N_wells || is.null(x$is_train) || !any(x$is_train == 0L) || !any(x$is_train == 1L)) {
            semantic <- c(semantic, "directional holdout mismatch")
          }
          if (transform == "ploidy_null" && any(x$ploidy_metric != 0)) semantic <- c(semantic, "ploidy-null metric is nonzero")
        } else if (transform %in% c("morphology_cell_area", "morphology_nuclear_area")) {
          same_shape <- x$N_lines == parent$N_lines && x$N_wells == parent$N_wells &&
            x$N_obs_count == parent$N_obs_count && x$N_obs_gluc == parent$N_obs_gluc
          if (!same_shape) semantic <- c(semantic, "morphology transform changed dataset dimensions")
          if (length(x$ploidy_metric) != x$N_wells || any(!is.finite(x$ploidy_metric))) {
            semantic <- c(semantic, "morphology metric is incomplete or non-finite")
          }
          if (!identical(as.numeric(x$ploidy_abs), as.numeric(parent$ploidy_abs))) {
            semantic <- c(semantic, "morphology transform changed ploidy_abs pair identity")
          }
          baseline_ok <- vapply(sort(unique(x$line_id)), function(line_id) {
            idx <- which(x$line_id == line_id)
            base_abs <- min(x$ploidy_abs[idx])
            all(abs(x$ploidy_metric[idx][x$ploidy_abs[idx] == base_abs]) < 1e-10)
          }, logical(1))
          if (!all(baseline_ok)) semantic <- c(semantic, "morphology metric baselines are not zero")
        }
      }
      check(sprintf("dataset %s transform is correct", datasets$dataset_id[[i]]), !length(semantic), paste(semantic, collapse = "; "))
    }
  }
  all_path <- datasets$stan_data_path[datasets$dataset_id == "all_lines"]
  if (length(all_path) == 1L && file.exists(all_path)) {
    x <- readRDS(all_path)
    if ("SUM-159-fuse" %in% names(x$line_map)) {
      affected <- which(x$exp_id[x$well_idx_gluc] == 3L)
      split <- as.integer(table(factor(x$ploidy_abs[x$well_idx_gluc[affected]], levels = c(2, 4))))
      gate_ok <- x$N_obs_gluc == 912L && length(affected) == 120L &&
        identical(split, c(60L, 60L)) && sum(x$is_censored[affected]) == 44L
      check("SUM-159 glucose mapping regression", gate_ok,
            sprintf("N=%d affected=%d split=%s censored=%d", x$N_obs_gluc, length(affected), paste(split, collapse = "/"), sum(x$is_censored[affected])))
    }
  }
}

fits <- NULL
if (file.exists(fit_plan_path)) {
  fits <- tryCatch(pe_read_tsv(fit_plan_path, c(
    "fit_id", "dataset_id", "model_id", "stan_data_path", "stan_data_sha256",
    "run_optim", "optim_starts", "optim_dir", "run_nuts", "nuts_chains", "nuts_dir"
  )), error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
}
if (!is.null(fits)) {
  check("fit ids are unique", !anyDuplicated(fits$fit_id))
  if (!is.null(datasets)) check("fit datasets exist", all(fits$dataset_id %in% datasets$dataset_id))
  check("fit Stan-data hashes match", all(vapply(seq_len(nrow(fits)), function(i) {
    file.exists(fits$stan_data_path[[i]]) && identical(pe_sha256(fits$stan_data_path[[i]]), fits$stan_data_sha256[[i]])
  }, logical(1))))
}

if (file.exists(nuts_plan_path)) {
  nuts <- pe_read_tsv(nuts_plan_path)
  if (nrow(nuts)) {
    required <- c("task_id", "fit_id", "config_path", "output_dir", "run_tag")
    check("NUTS plan schema", all(required %in% names(nuts)), paste(setdiff(required, names(nuts)), collapse = ", "))
    if (all(required %in% names(nuts))) {
      config_ok <- vapply(nuts$config_path, function(path) {
        file.exists(path) && !inherits(try(jsonlite::fromJSON(path, simplifyVector = FALSE), silent = TRUE), "try-error")
      }, logical(1))
      check("all NUTS configs parse", all(config_ok), sprintf("%d/%d", sum(config_ok), length(config_ok)))
    }
  }
}

if (phase == "complete" && !is.null(fits)) {
  for (i in seq_len(nrow(fits))) {
    fit <- fits[i, , drop = FALSE]
    if (fit$run_optim[[1]] == 1L) {
      ok <- all(file.exists(file.path(fit$optim_dir[[1]], c("optim_draws_all.Rds", "optim_lp_all.Rds", "optim_rc_all.Rds"))))
      check(sprintf("fit %s optimization complete", fit$fit_id[[1]]), ok, fit$optim_dir[[1]])
    }
    if (fit$run_nuts[[1]] == 1L) {
      expected <- sprintf("nuts_draws_chain%02d.Rds", seq_len(fit$nuts_chains[[1]]))
      ok <- all(file.exists(file.path(fit$nuts_dir[[1]], expected)))
      check(sprintf("fit %s NUTS complete", fit$fit_id[[1]]), ok, fit$nuts_dir[[1]])
    }
  }
}

out_dir <- file.path(root, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
report <- list(
  validation_passed = !length(errors), phase = phase, release_root = root,
  checked_at = format(Sys.time(), usetz = TRUE), checks = checks, errors = as.list(errors)
)
pe_write_json(report, file.path(out_dir, sprintf("%s_validation.json", phase)))
writeLines(c(
  sprintf("validation=%s", if (!length(errors)) "PASS" else "FAIL"),
  sprintf("phase=%s", phase),
  sprintf("checks=%d", length(checks)),
  sprintf("errors=%d", length(errors)),
  if (length(errors)) paste0("- ", errors)
), file.path(out_dir, sprintf("%s_validation.txt", phase)))

if (length(errors)) pe_fail("Release validation failed with %d error(s); see %s", length(errors), out_dir)
cat(sprintf("PASS: %d checks (%s phase)\n", length(checks), phase))
