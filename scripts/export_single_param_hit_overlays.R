#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

SCREEN_ROOT <- if (length(args) >= 1) args[1] else "data/gpath_transfer_cv_single_param"
MASK_SPEC <- if (length(args) >= 2) args[2] else "ae-1,ah-1"
OUTPUT_DIR <- if (length(args) >= 3) args[3] else file.path(SCREEN_ROOT, "hit_overlays")
STAN_DATA_PATH <- if (length(args) >= 4) args[4] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"
MODEL_NAME <- if (length(args) >= 5) args[5] else "gpath"

library(ggplot2)

source("R/project_paths.R")
source("R/transfer_cv_plot_utils.R")

cat(sprintf(">>> export_single_param_hit_overlays.R cwd: %s\n", getwd()))

parse_csv_arg <- function(x) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

mask_labels <- parse_csv_arg(MASK_SPEC)
if (!length(mask_labels)) {
  stop("No mask labels supplied")
}

stan_data <- readRDS(resolve_stan_data_path(STAN_DATA_PATH))
line_ids_all <- sort(unique(stan_data$line_id))

manifest_rows <- list()
idx <- 1L

for (mask_label in mask_labels) {
  mask_root <- file.path(SCREEN_ROOT, mask_label, MODEL_NAME)
  if (!dir.exists(mask_root)) {
    cat(sprintf(">>> Skip mask %s: directory not found %s\n", mask_label, mask_root))
    next
  }

  run_ids <- list.dirs(mask_root, recursive = FALSE, full.names = FALSE)
  run_ids <- run_ids[file.exists(file.path(mask_root, run_ids, "transfer_best_start_summary.Rds"))]

  if (!length(run_ids)) {
    cat(sprintf(">>> Skip mask %s: no summarized run_ids found\n", mask_label))
    next
  }

  for (run_id in run_ids) {
    best_path <- file.path(mask_root, run_id, "transfer_best_start_summary.Rds")
    best_df <- readRDS(best_path)
    line_ids <- sort(unique(best_df$line_id))
    if (!length(line_ids)) {
      line_ids <- line_ids_all
    }

    out_dir_model <- file.path(OUTPUT_DIR, mask_label, MODEL_NAME, run_id)
    dir.create(out_dir_model, recursive = TRUE, showWarnings = FALSE)

    for (line_id in line_ids) {
      cat(sprintf(">>> Rendering overlay | mask=%s | model=%s | line=%d\n", mask_label, run_id, line_id))

      overlay_data <- tryCatch(
        generate_transfer_overlay_data(
          model_id = run_id,
          line_id = line_id,
          output_root = file.path(SCREEN_ROOT, mask_label),
          stan_data_path = STAN_DATA_PATH
        ),
        error = function(e) {
          cat(sprintf(">>> Skip mask=%s model=%s line=%d: %s\n", mask_label, run_id, line_id, conditionMessage(e)))
          NULL
        }
      )

      if (is.null(overlay_data)) {
        next
      }

      score_text <- tryCatch(
        format_transfer_overlay_scores(
          model_id = run_id,
          line_id = line_id,
          output_root = file.path(SCREEN_ROOT, mask_label),
          model_name = MODEL_NAME
        ),
        error = function(e) NULL
      )

      plt <- plot_transfer_line_trajectories(
        transfer_data = overlay_data,
        line_id = line_id,
        model_id = sprintf("%s | %s", run_id, mask_label),
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
        mask_label = mask_label,
        model_id = run_id,
        line_id = as.integer(line_id),
        output_path = out_path,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
}

manifest <- if (length(manifest_rows)) {
  do.call(rbind, manifest_rows)
} else {
  data.frame(
    mask_label = character(),
    model_id = character(),
    line_id = integer(),
    output_path = character(),
    stringsAsFactors = FALSE
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(OUTPUT_DIR, "overlay_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

cat(sprintf(">>> Wrote %d overlays to %s\n", nrow(manifest), OUTPUT_DIR))
cat(sprintf(">>> Manifest: %s\n", manifest_path))
