#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

MODEL_SPEC <- if (length(args) >= 1) args[1] else "all"
LINE_SPEC <- if (length(args) >= 2) args[2] else "all"
TRANSFER_ROOT <- if (length(args) >= 3) args[3] else "data/gpath_transfer_cv"
OUTPUT_DIR <- if (length(args) >= 4) args[4] else "data/gpath_transfer_cv_overlays"
STAN_DATA_PATH <- if (length(args) >= 5) args[5] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"
MODEL_NAME <- if (length(args) >= 6) args[6] else "gpath"

library(ggplot2)

source("R/transfer_cv_plot_utils.R")

cat(sprintf(">>> export_gpath_transfer_overlays.R cwd: %s\n", getwd()))

parse_csv_arg <- function(x) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

discover_transfer_models <- function(transfer_root, model_name) {
  base_dir <- file.path(transfer_root, model_name)
  if (!dir.exists(base_dir)) {
    stop(sprintf("Transfer model directory not found: %s", base_dir))
  }

  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  dirs[file.exists(file.path(base_dir, dirs, "transfer_best_start_summary.Rds"))]
}

resolve_line_ids <- function(model_id, line_spec, transfer_root, model_name, fallback_lines) {
  if (!identical(line_spec, "all")) {
    return(as.integer(parse_csv_arg(line_spec)))
  }

  best_path <- file.path(transfer_root, model_name, model_id, "transfer_best_start_summary.Rds")
  if (!file.exists(best_path)) {
    return(fallback_lines)
  }

  sort(unique(readRDS(best_path)$line_id))
}

stan_data <- readRDS(resolve_stan_data_path(STAN_DATA_PATH))
all_line_ids <- sort(unique(stan_data$line_id))

model_ids <- if (identical(MODEL_SPEC, "all")) {
  discover_transfer_models(TRANSFER_ROOT, MODEL_NAME)
} else {
  parse_csv_arg(MODEL_SPEC)
}

if (!length(model_ids)) {
  stop("No model IDs found for overlay export")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

manifest_rows <- list()
idx <- 1L

for (model_id in model_ids) {
  line_ids <- resolve_line_ids(
    model_id = model_id,
    line_spec = LINE_SPEC,
    transfer_root = TRANSFER_ROOT,
    model_name = MODEL_NAME,
    fallback_lines = all_line_ids
  )

  if (!length(line_ids)) {
    cat(sprintf(">>> Skipping %s: no line IDs resolved\n", model_id))
    next
  }

  out_dir_model <- file.path(OUTPUT_DIR, MODEL_NAME, model_id)
  dir.create(out_dir_model, recursive = TRUE, showWarnings = FALSE)

  for (line_id in line_ids) {
    cat(sprintf(">>> Rendering overlay | model=%s | line=%d\n", model_id, line_id))

    overlay_data <- tryCatch(
      generate_transfer_overlay_data(
        model_id = model_id,
        line_id = line_id,
        output_root = TRANSFER_ROOT,
        stan_data_path = STAN_DATA_PATH
      ),
      error = function(e) {
        cat(sprintf(">>> Skip model=%s line=%d: %s\n", model_id, line_id, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(overlay_data)) {
      next
    }

    score_text <- tryCatch(
      format_transfer_overlay_scores(
        model_id = model_id,
        line_id = line_id,
        output_root = TRANSFER_ROOT,
        model_name = MODEL_NAME
      ),
      error = function(e) NULL
    )

    plt <- plot_transfer_line_trajectories(
      transfer_data = overlay_data,
      line_id = line_id,
      model_id = model_id,
      score_text = score_text
    )

    out_path <- file.path(out_dir_model, sprintf("line_%02d_overlay.png", as.integer(line_id)))
    ggsave(
      filename = out_path,
      plot = plt,
      width = 16,
      height = 8,
      units = "in",
      dpi = 200
    )

    manifest_rows[[idx]] <- data.frame(
      model_id = model_id,
      line_id = as.integer(line_id),
      output_path = out_path,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

manifest <- if (length(manifest_rows)) {
  do.call(rbind, manifest_rows)
} else {
  data.frame(
    model_id = character(),
    line_id = integer(),
    output_path = character(),
    stringsAsFactors = FALSE
  )
}
manifest_path <- file.path(OUTPUT_DIR, MODEL_NAME, "overlay_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
cat(sprintf(">>> Wrote %d overlays to %s\n", nrow(manifest), file.path(OUTPUT_DIR, MODEL_NAME)))
cat(sprintf(">>> Manifest: %s\n", manifest_path))
