#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_DIR <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp4b_no_sum159_fuse_gpath")
}
FIGURE_DIR <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp4b_no_sum159_fuse_gpath")
}
NO_SUM_SPEC_PATH <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("data", "specs", "wp4b_no_sum159_fuse_assessment.tsv")
}
ALL_LINE_SPEC_PATH <- if (length(args) >= 4 && nzchar(args[4])) {
  args[4]
} else {
  file.path("data", "specs", "gpath_optim_from_datasets.tsv")
}
NO_SUM_STAN_DATA_PATH <- if (length(args) >= 5 && nzchar(args[5])) {
  args[5]
} else {
  file.path("data", "inputs", "stan", "wp4b_no_sum159_fuse", "stan_ready_data.Rds")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/optim_utils.R")

wp4b_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey72"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

save_plot_pair <- function(plot, basename, width, height, dpi = 450) {
  pdf_path <- file.path(FIGURE_DIR, paste0(basename, ".pdf"))
  png_path <- file.path(FIGURE_DIR, paste0(basename, ".png"))

  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)

  invisible(c(pdf = pdf_path, png = png_path))
}

parse_dataset_label <- function(dataset_label) {
  out <- data.frame(
    fit_type = "full",
    line_id = NA_integer_,
    direction = NA_character_,
    stringsAsFactors = FALSE
  )

  m <- regexec("_(transfer|null)_line([0-9]+)_(low_to_high|high_to_low)$", dataset_label)
  p <- regmatches(dataset_label, m)[[1]]
  if (length(p) == 4L) {
    out$fit_type <- p[[2]]
    out$line_id <- as.integer(p[[3]])
    out$direction <- p[[4]]
  }

  out
}

summarize_spec_scope <- function(spec_path, scope) {
  if (!file.exists(spec_path)) {
    return(tibble(
      scope = character(),
      spec_path = character(),
      spec_row_id = integer(),
      run_id = character(),
      dataset_label = character(),
      run_label = character(),
      optim_output_dir = character(),
      optim_complete = logical(),
      fit_type = character(),
      line_id = integer(),
      direction = character(),
      n_starts_total = integer(),
      n_starts_valid = integer(),
      best_lp = numeric(),
      log_lik_train = numeric(),
      log_lik_holdout = numeric()
    ))
  }

  spec_df <- read_optim_spec(spec_path)
  spec_df <- spec_df[as.integer(spec_df$enabled) != 0L, , drop = FALSE]
  if (!nrow(spec_df)) {
    return(tibble())
  }

  parsed <- bind_rows(lapply(spec_df$dataset_label, parse_dataset_label))
  spec_df <- bind_cols(as_tibble(spec_df), parsed)

  metric_rows <- lapply(seq_len(nrow(spec_df)), function(i) {
    if (!isTRUE(spec_df$optim_complete[[i]])) {
      return(tibble(
        spec_row_id = spec_df$spec_row_id[[i]],
        n_starts_total = NA_integer_,
        n_starts_valid = NA_integer_,
        best_lp = NA_real_,
        log_lik_train = NA_real_,
        log_lik_holdout = NA_real_
      ))
    }

    x <- tryCatch(
      summarize_optim_row(spec_df[i, , drop = FALSE], include_k = FALSE),
      error = function(e) NULL
    )
    if (is.null(x)) {
      return(tibble(
        spec_row_id = spec_df$spec_row_id[[i]],
        n_starts_total = NA_integer_,
        n_starts_valid = NA_integer_,
        best_lp = NA_real_,
        log_lik_train = NA_real_,
        log_lik_holdout = NA_real_
      ))
    }
    as_tibble(x) %>%
      select(spec_row_id, n_starts_total, n_starts_valid, best_lp, log_lik_train, log_lik_holdout)
  })

  spec_df %>%
    left_join(bind_rows(metric_rows), by = "spec_row_id") %>%
    mutate(scope = scope, spec_path = spec_path, .before = 1) %>%
    select(
      scope, spec_path, spec_row_id, run_id, dataset_label, run_label,
      optim_output_dir, optim_complete, fit_type, line_id, direction,
      n_starts_total, n_starts_valid, best_lp, log_lik_train, log_lik_holdout
    )
}

line_label_map <- function(path) {
  if (!file.exists(path)) {
    return(setNames(character(0), character(0)))
  }
  x <- readRDS(path)
  if (!("line_map" %in% names(x)) || is.null(x$line_map)) {
    return(setNames(character(0), character(0)))
  }
  setNames(names(x$line_map), as.character(as.integer(x$line_map)))
}

label_lines <- function(line_id, labels) {
  out <- unname(labels[as.character(as.integer(line_id))])
  fallback <- sprintf("line%02d", as.integer(line_id))
  out[is.na(out) | !nzchar(out)] <- fallback[is.na(out) | !nzchar(out)]
  out
}

build_transfer_gain_table <- function(status_df) {
  x <- status_df %>%
    filter(
      optim_complete,
      fit_type %in% c("transfer", "null"),
      !is.na(line_id),
      is.finite(log_lik_holdout)
    ) %>%
    select(
      scope, run_id, line_id, direction, fit_type,
      log_lik_holdout, log_lik_train, best_lp, n_starts_valid
    )

  if (!nrow(x)) {
    return(tibble())
  }

  x %>%
    pivot_wider(
      names_from = fit_type,
      values_from = c(log_lik_holdout, log_lik_train, best_lp, n_starts_valid)
    ) %>%
    mutate(
      transfer_minus_null_holdout = log_lik_holdout_transfer - log_lik_holdout_null,
      transfer_minus_null_train = log_lik_train_transfer - log_lik_train_null,
      transfer_better_than_null = transfer_minus_null_holdout > 0
    )
}

build_gain_comparison <- function(gain_df) {
  if (!nrow(gain_df)) {
    return(tibble())
  }

  all_line <- gain_df %>%
    filter(scope == "all_lines", line_id <= 4L) %>%
    select(
      run_id, line_id, direction,
      all_line_transfer_minus_null_holdout = transfer_minus_null_holdout,
      all_line_transfer_minus_null_train = transfer_minus_null_train,
      all_line_transfer_better_than_null = transfer_better_than_null
    )

  no_sum <- gain_df %>%
    filter(scope == "no_sum159_fuse") %>%
    select(
      run_id, line_id, direction,
      no_sum_transfer_minus_null_holdout = transfer_minus_null_holdout,
      no_sum_transfer_minus_null_train = transfer_minus_null_train,
      no_sum_transfer_better_than_null = transfer_better_than_null
    )

  inner_join(no_sum, all_line, by = c("run_id", "line_id", "direction")) %>%
    mutate(
      delta_transfer_minus_null_holdout =
        no_sum_transfer_minus_null_holdout - all_line_transfer_minus_null_holdout,
      delta_transfer_minus_null_train =
        no_sum_transfer_minus_null_train - all_line_transfer_minus_null_train,
      no_sum_improved_over_all_line = delta_transfer_minus_null_holdout > 0
    )
}

write_decision_note <- function(path, status_df, gain_df, comparison_df) {
  no_sum_status <- status_df %>% filter(scope == "no_sum159_fuse")
  all_line_status <- status_df %>% filter(scope == "all_lines")

  no_sum_complete <- sum(no_sum_status$optim_complete, na.rm = TRUE)
  no_sum_total <- nrow(no_sum_status)
  all_line_complete <- sum(all_line_status$optim_complete, na.rm = TRUE)
  all_line_total <- nrow(all_line_status)

  no_sum_transfer_complete <- no_sum_status %>%
    filter(fit_type %in% c("transfer", "null")) %>%
    summarise(n = sum(optim_complete, na.rm = TRUE), total = n(), .groups = "drop")

  result_sentence <- if (nrow(comparison_df)) {
    sprintf(
      paste(
        "Completed matched comparisons: %d model-line-direction cases.",
        "Mean no-SUM minus all-line transfer/null holdout gain is %.3f;",
        "%.1f%% of matched cases improve."
      ),
      nrow(comparison_df),
      mean(comparison_df$delta_transfer_minus_null_holdout, na.rm = TRUE),
      100 * mean(comparison_df$no_sum_improved_over_all_line, na.rm = TRUE)
    )
  } else {
    "No completed matched no-SUM versus all-line transfer/null comparisons were available when this report was generated."
  }

  nuts_sentence <- if (nrow(comparison_df)) {
    paste(
      "NUTS remains intentionally skipped by this WP4b report.",
      "Run posterior sampling only after selecting representative no-SUM models from completed optimisation results."
    )
  } else {
    paste(
      "NUTS and visual-fit overlays are intentionally skipped because the no-SUM optimisation benchmark is not complete enough",
      "to justify posterior sampling or representative-model plotting."
    )
  }

  lines <- c(
    "# WP4b no-SUM-159-fuse gpath decision note",
    "",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Rationale",
    "",
    paste(
      "SUM-159-fuse is treated as an exception candidate before this gpath rerun because it is",
      "the only cell-line dataset with two batches and differing protocol details, it has weak or",
      "reversed high-versus-low ploidy volume separation, and the transferable mechanistic conclusion",
      "is directional rather than generic holdout fit: ploidy-effect directions on y1, v1, K1, and m",
      "are conserved across tested models, move further from zero, and show stronger inter-model",
      "agreement after SUM-159-fuse is excluded."
    ),
    "",
    "## Optimisation status",
    "",
    sprintf("- No-SUM spec rows complete: %d/%d", no_sum_complete, no_sum_total),
    sprintf("- No-SUM transfer/null rows complete: %d/%d", no_sum_transfer_complete$n, no_sum_transfer_complete$total),
    sprintf("- All-line baseline spec rows complete: %d/%d", all_line_complete, all_line_total),
    "",
    "## Current comparison",
    "",
    result_sentence,
    "",
    "## Expensive jobs",
    "",
    nuts_sentence,
    "",
    "## Files",
    "",
    sprintf("- Status table: %s", file.path(OUTPUT_DIR, "wp4b_optim_status_by_scope.csv")),
    sprintf("- Transfer gain table: %s", file.path(OUTPUT_DIR, "wp4b_transfer_gain_by_scope.csv")),
    sprintf("- Matched comparison table: %s", file.path(OUTPUT_DIR, "wp4b_transfer_gain_comparison.csv"))
  )

  writeLines(lines, path, useBytes = TRUE)
}

status_df <- bind_rows(
  summarize_spec_scope(ALL_LINE_SPEC_PATH, "all_lines"),
  summarize_spec_scope(NO_SUM_SPEC_PATH, "no_sum159_fuse")
)

labels <- line_label_map(NO_SUM_STAN_DATA_PATH)
status_df <- status_df %>%
  mutate(
    line_label = if_else(is.na(line_id), NA_character_, label_lines(line_id, labels))
  )

gain_df <- build_transfer_gain_table(status_df) %>%
  mutate(
    line_label = label_lines(line_id, labels),
    line_direction = paste(line_label, direction, sep = " / ")
  )
comparison_df <- build_gain_comparison(gain_df) %>%
  mutate(
    line_label = label_lines(line_id, labels),
    line_direction = paste(line_label, direction, sep = " / ")
  )

gain_summary_df <- if (nrow(gain_df)) {
  gain_df %>%
    group_by(scope, run_id) %>%
    summarise(
      n_cases = n(),
      mean_transfer_minus_null_holdout = mean(transfer_minus_null_holdout, na.rm = TRUE),
      median_transfer_minus_null_holdout = median(transfer_minus_null_holdout, na.rm = TRUE),
      win_rate = mean(transfer_better_than_null, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scope, desc(mean_transfer_minus_null_holdout))
} else {
  tibble()
}

write.csv(status_df, file.path(OUTPUT_DIR, "wp4b_optim_status_by_scope.csv"), row.names = FALSE)
write.csv(gain_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_by_scope.csv"), row.names = FALSE)
write.csv(comparison_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_comparison.csv"), row.names = FALSE)
write.csv(gain_summary_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_summary.csv"), row.names = FALSE)
saveRDS(status_df, file.path(OUTPUT_DIR, "wp4b_optim_status_by_scope.Rds"))
saveRDS(gain_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_by_scope.Rds"))
saveRDS(comparison_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_comparison.Rds"))
saveRDS(gain_summary_df, file.path(OUTPUT_DIR, "wp4b_transfer_gain_summary.Rds"))

if (nrow(gain_df)) {
  gain_plot <- gain_df %>%
    filter(scope %in% c("all_lines", "no_sum159_fuse"), line_id <= 4L) %>%
    mutate(
      scope = recode(scope, all_lines = "All lines", no_sum159_fuse = "SUM-159-fuse excluded"),
      run_id = factor(run_id, levels = unique(run_id))
    ) %>%
    ggplot(aes(run_id, line_direction, fill = transfer_minus_null_holdout)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~scope, ncol = 1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      oob = squish,
      name = "transfer - null\nholdout log lik."
    ) +
    labs(
      title = "WP4b gpath transfer/null holdout gain by scope",
      x = NULL,
      y = NULL
    ) +
    wp4b_theme(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  save_plot_pair(gain_plot, "wp4b_transfer_gain_by_scope", width = 11.5, height = 7.5)
}

if (nrow(comparison_df)) {
  comparison_plot <- comparison_df %>%
    mutate(run_id = factor(run_id, levels = unique(run_id))) %>%
    ggplot(aes(run_id, line_direction, fill = delta_transfer_minus_null_holdout)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      oob = squish,
      name = "no-SUM - all-line\nholdout gain"
    ) +
    labs(
      title = "WP4b no-SUM versus all-line gpath transfer gain",
      x = NULL,
      y = NULL,
      caption = "Positive values mean excluding SUM-159-fuse improved transfer-minus-null holdout log likelihood on matched non-SUM holdouts."
    ) +
    wp4b_theme(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  save_plot_pair(comparison_plot, "wp4b_no_sum_vs_all_line_transfer_gain", width = 11.5, height = 5.2)
}

write_decision_note(
  path = file.path(OUTPUT_DIR, "wp4b_no_sum159_fuse_gpath_decision_note.md"),
  status_df = status_df,
  gain_df = gain_df,
  comparison_df = comparison_df
)

cat(sprintf("Wrote WP4b report tables to %s\n", OUTPUT_DIR))
cat(sprintf("Wrote WP4b figures, where data were available, to %s\n", FIGURE_DIR))
