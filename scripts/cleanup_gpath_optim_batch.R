#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

batch_label <- if (length(args) >= 1) args[1] else stop("batch_label required")
delete_metadata <- if (length(args) >= 2) as.logical(as.integer(args[2])) else TRUE
dry_run <- if (length(args) >= 3) as.logical(as.integer(args[3])) else TRUE

batch_dir <- file.path("slurm", "runs", "gpath_optim_batch", batch_label)
job_table <- file.path(batch_dir, "job_ids.tsv")

cat(sprintf(">>> cleanup_gpath_optim_batch.R cwd: %s\n", getwd()))
cat(sprintf(">>> Batch dir: %s\n", batch_dir))
cat(sprintf(">>> Dry run: %s\n", dry_run))

if (!file.exists(job_table)) {
  stop(sprintf("Job table not found: %s", job_table))
}

jobs <- read.delim(job_table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!nrow(jobs)) {
  cat("No jobs found in job table.\n")
} else {
  output_dirs <- unique(jobs$output_dir)
  output_dirs <- output_dirs[nzchar(output_dirs)]

  cat("Output directories targeted for cleanup:\n")
  cat(paste0(" - ", output_dirs), sep = "\n")
  cat("\n")

  if (!dry_run) {
    for (out_dir in output_dirs) {
      if (dir.exists(out_dir)) {
        unlink(out_dir, recursive = TRUE, force = TRUE)
        cat(sprintf("Removed output dir: %s\n", out_dir))
      } else {
        cat(sprintf("Output dir already absent: %s\n", out_dir))
      }
    }
  }
}

if (delete_metadata) {
  if (dry_run) {
    cat(sprintf("Would remove batch metadata dir: %s\n", batch_dir))
  } else if (dir.exists(batch_dir)) {
    unlink(batch_dir, recursive = TRUE, force = TRUE)
    cat(sprintf("Removed batch metadata dir: %s\n", batch_dir))
  }
}
