normalize_fit_type <- function(fit_type) {
  if (!length(fit_type) || is.na(fit_type)) {
    stop("fit_type must be provided")
  }

  key <- tolower(trimws(fit_type))
  valid <- c("null", "transfer", "oracle")

  if (!(key %in% valid)) {
    stop(sprintf("Unsupported fit_type '%s'", fit_type))
  }

  key
}

build_transfer_run_id <- function(
  line_id,
  direction,
  fit_type,
  start_id
) {
  sprintf(
    "line%02d_%s_%s_start%02d",
    as.integer(line_id),
    normalize_transfer_direction(direction),
    normalize_fit_type(fit_type),
    as.integer(start_id)
  )
}

build_transfer_output_dir <- function(
  output_root,
  model_name,
  run_id,
  direction,
  fit_type
) {
  file.path(
    output_root,
    model_name,
    run_id,
    normalize_transfer_direction(direction),
    normalize_fit_type(fit_type)
  )
}

save_transfer_run_outputs <- function(
  fit,
  output_dir,
  run_tag,
  split_meta
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summary_path <- file.path(output_dir, sprintf("optim_summary_%s.Rds", run_tag))
  draws_path <- file.path(output_dir, sprintf("optim_draws_%s.Rds", run_tag))
  diag_path <- file.path(output_dir, sprintf("optim_diagnostics_%s.Rds", run_tag))
  meta_path <- file.path(output_dir, sprintf("split_meta_%s.Rds", run_tag))

  saveRDS(fit$summary(), summary_path)
  saveRDS(fit$draws(), draws_path)

  diag <- tryCatch(fit$diagnostic_summary(), error = function(e) NULL)
  if (!is.null(diag)) {
    saveRDS(diag, diag_path)
  }

  saveRDS(split_meta, meta_path)

  invisible(list(
    summary_path = summary_path,
    draws_path = draws_path,
    diag_path = diag_path,
    meta_path = meta_path
  ))
}

summarize_transfer_draws <- function(draws, split_meta) {
  ll <- compute_well_subset_elpd(draws, split_meta$holdout_wells)
  dm <- posterior::as_draws_matrix(draws)
  lp_col <- if ("lp__" %in% colnames(dm)) "lp__" else NULL
  lp_value <- if (is.null(lp_col)) NA_real_ else as.numeric(dm[1, lp_col])

  data.frame(
    line_id = split_meta$line_id,
    direction = split_meta$direction,
    fit_type = split_meta$fit_type,
    holdout_value = split_meta$holdout_value,
    observed_value = split_meta$observed_value,
    holdout_wells = length(split_meta$holdout_wells),
    holdout_obs = split_meta$holdout_total_obs,
    lp__ = lp_value,
    elpd = ll$elpd,
    stringsAsFactors = FALSE
  )
}

collect_transfer_result_files <- function(output_dir) {
  draw_files <- list.files(output_dir, pattern = "^optim_draws_.*\\.Rds$", full.names = TRUE)
  meta_files <- list.files(output_dir, pattern = "^split_meta_.*\\.Rds$", full.names = TRUE)

  list(
    draw_files = sort(draw_files),
    meta_files = sort(meta_files)
  )
}

summarize_transfer_directory <- function(output_dir) {
  files <- collect_transfer_result_files(output_dir)
  if (!length(files$draw_files) || !length(files$meta_files)) {
    stop(sprintf("No transfer result files found in '%s'", output_dir))
  }

  meta_by_tag <- stats::setNames(
    lapply(files$meta_files, readRDS),
    sub("^split_meta_(.*)\\.Rds$", "\\1", basename(files$meta_files))
  )

  rows <- lapply(files$draw_files, function(draw_file) {
    tag <- sub("^optim_draws_(.*)\\.Rds$", "\\1", basename(draw_file))
    split_meta <- meta_by_tag[[tag]]
    if (is.null(split_meta)) {
      stop(sprintf("Missing split metadata for '%s'", tag))
    }

    draws <- readRDS(draw_file)
    out <- summarize_transfer_draws(draws, split_meta)
    out$start_tag <- tag
    out
  })

  do.call(rbind, rows)
}
