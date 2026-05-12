#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_DIR <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp5_mechanistic_representative_model")
}
FIGURE_DIR <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp5_mechanistic_representative_model")
}
SPEC_PATH <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("data", "specs", "gpath_optim_from_datasets.tsv")
}
STAN_DATA_PATH <- if (length(args) >= 4 && nzchar(args[4])) {
  args[4]
} else {
  file.path("data", "stan_ready_data.Rds")
}
NUTS_JOB_CONFIG_PATH <- if (length(args) >= 5 && nzchar(args[5])) {
  args[5]
} else {
  file.path("config", "jobs", "gpath_full_oracle_assessment.json")
}
MODEL_ANALYSIS_EXPORT_DIR <- if (length(args) >= 6 && nzchar(args[6])) {
  args[6]
} else {
  file.path("data", "report_exports", "model_analysis")
}
MODEL_ANALYSIS_PARETO_OVERLAY_DIR <- if (length(args) >= 7 && nzchar(args[7])) {
  args[7]
} else {
  file.path("dev-figs", "model_analysis", "pareto_global_fit_overlays")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/optim_utils.R")
source("R/model_analysis_utils.R")
source("R/sim_utils.R")
source("R/plot_utils.R")
source("R/job_config_utils.R")
source(get_model_r_path("gpath", "v1"))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

wp5_theme <- function(base_size = 8) {
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

write_table_pair <- function(x, stem) {
  csv_path <- file.path(OUTPUT_DIR, paste0(stem, ".csv"))
  rds_path <- file.path(OUTPUT_DIR, paste0(stem, ".Rds"))
  utils::write.csv(as.data.frame(x), csv_path, row.names = FALSE)
  saveRDS(x, rds_path)
  invisible(c(csv = csv_path, rds = rds_path))
}

parse_dataset_label <- function(dataset_label) {
  out <- tibble(
    fit_type = "full",
    line_id = NA_integer_,
    direction = NA_character_
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

safe_metric_summary <- function(spec_row) {
  if (!isTRUE(spec_row$optim_complete[[1]])) {
    return(tibble(
      spec_row_id = spec_row$spec_row_id[[1]],
      n_starts_total = NA_integer_,
      n_starts_valid = NA_integer_,
      best_draw_index = NA_integer_,
      best_lp = NA_real_,
      n_within_1 = NA_integer_,
      n_within_5 = NA_integer_,
      log_lik_train = NA_real_,
      log_lik_holdout = NA_real_
    ))
  }

  tryCatch(
    as_tibble(summarize_optim_row(spec_row, include_k = FALSE)),
    error = function(e) {
      warning(sprintf(
        "Could not summarize optimisation row %s (%s / %s): %s",
        spec_row$spec_row_id[[1]], spec_row$dataset_label[[1]], spec_row$run_id[[1]],
        conditionMessage(e)
      ))
      tibble(
        spec_row_id = spec_row$spec_row_id[[1]],
        n_starts_total = NA_integer_,
        n_starts_valid = NA_integer_,
        best_draw_index = NA_integer_,
        best_lp = NA_real_,
        n_within_1 = NA_integer_,
        n_within_5 = NA_integer_,
        log_lik_train = NA_real_,
        log_lik_holdout = NA_real_
      )
    }
  )
}

summarize_spec <- function(spec_path, stan_data_ref) {
  spec_df <- read_optim_spec(spec_path)
  spec_df <- spec_df[as.integer(spec_df$enabled) != 0L, , drop = FALSE]
  if (!nrow(spec_df)) {
    return(tibble())
  }

  parsed <- bind_rows(lapply(spec_df$dataset_label, parse_dataset_label))
  metric_df <- bind_rows(lapply(seq_len(nrow(spec_df)), function(i) {
    safe_metric_summary(spec_df[i, , drop = FALSE])
  }))
  k_map <- tibble(
    run_id = unique(as.character(spec_df$run_id)),
    k = vapply(
      unique(as.character(spec_df$run_id)),
      function(run_id) get_hierarchical_k(run_id, N_lines = stan_data_ref$N_lines),
      numeric(1)
    )
  )

  bind_cols(as_tibble(spec_df), parsed) %>%
    left_join(metric_df, by = "spec_row_id") %>%
    left_join(k_map, by = "run_id") %>%
    mutate(
      model_id = run_id,
      alias = vapply(model_id, build_model_alias, character(1), format = "text")
    )
}

build_transfer_gain_table <- function(spec_summary_df) {
  x <- spec_summary_df %>%
    filter(
      optim_complete,
      fit_type %in% c("transfer", "null"),
      !is.na(line_id),
      is.finite(log_lik_holdout)
    ) %>%
    select(
      model_id, alias, line_id, direction, fit_type,
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

build_nuts_qc <- function(job_config_path, rhat_threshold = 1.01, ess_threshold = 400) {
  if (!file.exists(job_config_path)) {
    return(tibble())
  }

  job_info <- read_job_info(job_config_path)
  bind_rows(lapply(job_info$run_ids, function(model_id) {
    dir_path <- get_job_output_dir(job_info$job_cfg, model_id)
    diag_path <- file.path(dir_path, "nuts_chain_diagnostics.Rds")
    if (!file.exists(diag_path)) {
      return(tibble(
        model_id = model_id,
        nuts_available = FALSE,
        n_parameters = NA_integer_,
        n_bad_rhat = NA_integer_,
        n_bad_ess_bulk = NA_integer_,
        n_bad_ess_tail = NA_integer_,
        max_rhat = NA_real_,
        min_ess_bulk = NA_real_,
        min_ess_tail = NA_real_,
        nuts_qc_pass = NA
      ))
    }

    chain_diag_df <- readRDS(diag_path)
    tibble(
      model_id = model_id,
      nuts_available = TRUE,
      n_parameters = nrow(chain_diag_df),
      n_bad_rhat = sum(chain_diag_df$rhat > rhat_threshold, na.rm = TRUE),
      n_bad_ess_bulk = sum(chain_diag_df$ess_bulk < ess_threshold, na.rm = TRUE),
      n_bad_ess_tail = sum(chain_diag_df$ess_tail < ess_threshold, na.rm = TRUE),
      max_rhat = max(chain_diag_df$rhat, na.rm = TRUE),
      min_ess_bulk = min(chain_diag_df$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(chain_diag_df$ess_tail, na.rm = TRUE),
      nuts_qc_pass =
        sum(chain_diag_df$rhat > rhat_threshold, na.rm = TRUE) == 0L &&
        sum(chain_diag_df$ess_bulk < ess_threshold, na.rm = TRUE) == 0L &&
        sum(chain_diag_df$ess_tail < ess_threshold, na.rm = TRUE) == 0L
    )
  }))
}

build_nuts_qc_from_model_analysis <- function(export_dir) {
  path <- file.path(export_dir, "nuts_qc_summary_df.csv")
  if (!file.exists(path)) {
    return(tibble())
  }

  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("model_id", "n_parameters", "bad_rhat", "bad_ess_bulk", "bad_ess_tail", "max_rhat", "min_ess_bulk", "min_ess_tail")
  if (!all(required_cols %in% names(x))) {
    warning(sprintf("Ignoring %s because it does not contain the expected columns.", path))
    return(tibble())
  }

  as_tibble(x) %>%
    transmute(
      model_id,
      nuts_available = TRUE,
      n_parameters,
      n_bad_rhat = bad_rhat,
      n_bad_ess_bulk = bad_ess_bulk,
      n_bad_ess_tail = bad_ess_tail,
      max_rhat,
      min_ess_bulk,
      min_ess_tail,
      nuts_qc_pass = bad_rhat == 0L & bad_ess_bulk == 0L & bad_ess_tail == 0L,
      nuts_qc_source = "model_analysis/nuts_qc_summary_df.csv"
    )
}

build_model_analysis_screen_summary <- function(export_dir) {
  path <- file.path(export_dir, "model_screen_df.csv")
  if (!file.exists(path)) {
    return(tibble())
  }

  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE) %>%
    as_tibble() %>%
    select(
      model_id,
      model_analysis_mean_rank = mean_rank,
      model_analysis_directional_agreement = prop_directional_agreement_model,
      model_analysis_parameter_log_error_improvement = mean_parameter_log_error_improvement,
      model_analysis_parameter_win_rate = mean_prop_parameters_better
    )
}

build_pareto_front <- function(ranking_df) {
  if (!nrow(ranking_df)) {
    return(tibble())
  }

  ranking_df %>%
    arrange(k, deviance) %>%
    filter(deviance == cummin(deviance)) %>%
    group_by(deviance) %>%
    slice_min(order_by = k, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(k) %>%
    mutate(pareto_rank = row_number())
}

build_existing_pareto_overlay_manifest <- function(overlay_dir, stan_data, pareto_model_ids) {
  line_ids <- sort(unique(as.integer(stan_data$line_id)))
  line_label_map <- make_line_label_map(stan_data$line_map)
  bind_rows(lapply(line_ids, function(line_id) {
    line_name <- unname(line_label_map[as.character(line_id)])
    if (is.na(line_name) || !nzchar(line_name)) {
      line_name <- sprintf("line_%02d", line_id)
    }
    out_path <- file.path(overlay_dir, sprintf("line_%02d_pareto_overlay.png", line_id))
    tibble(
      line_id = line_id,
      line_name = line_name,
      output_path = out_path,
      file_exists = file.exists(out_path),
      source = "workflow/model_analysis.Rmd:selected-global-fit-overlay",
      pareto_model_ids = paste(pareto_model_ids, collapse = ",")
    )
  }))
}

build_visual_review_table <- function(pareto_front_df) {
  if (!nrow(pareto_front_df)) {
    return(tibble())
  }

  pareto_front_df %>%
    transmute(
      model_id,
      alias,
      pareto_rank,
      visual_fit_status = if_else(
        model_id == "1R_1P_0W_C0_M1",
        "rejected_by_human_review",
        "acceptable_by_human_review"
      ),
      visual_review_source = "user_review_2026-05-12_existing_model_analysis_pareto_overlays",
      visual_review_note = if_else(
        model_id == "1R_1P_0W_C0_M1",
        "Legacy 1R overlay judged visually unacceptable for representative-model use.",
        "Existing model_analysis Pareto-front overlay judged visually acceptable."
      )
    )
}

extract_full_fit_draw <- function(full_fit_row) {
  tryCatch(
    extract_best_draw_from_optim_row(full_fit_row),
    error = function(e) {
      warning(sprintf(
        "Could not extract best draw for %s: %s",
        full_fit_row$model_id[[1]],
        conditionMessage(e)
      ))
      NULL
    }
  )
}

build_fit_metric_for_model <- function(fit_row, draw_vec, stan_data, obs_df, sim_times) {
  if (is.null(draw_vec)) {
    return(tibble())
  }

  sim_df <- bind_rows(lapply(sort(unique(as.integer(stan_data$line_id))), function(line_id) {
    line_name <- names(stan_data$line_map)[as.integer(unname(stan_data$line_map)) == line_id]
    if (!length(line_name)) {
      line_name <- sprintf("line %d", line_id)
    }
    simulate_line_from_draw(
      model_id = fit_row$model_id[[1]],
      model_r_path = get_model_r_path("gpath", "v1"),
      draw_vec = draw_vec,
      stan_data = stan_data,
      line_id = line_id,
      times = sim_times,
      line_name = line_name[[1]]
    )
  })) %>%
    filter(variable %in% c("N", "R1"))

  cmp_df <- obs_df %>%
    filter(variable %in% c("N", "R1")) %>%
    inner_join(
      sim_df %>%
        select(well_idx, time, variable, pred_value = value),
      by = c("well_idx", "time", "variable")
    )

  if (!nrow(cmp_df)) {
    return(tibble())
  }

  metric_by_var <- cmp_df %>%
    group_by(variable) %>%
    summarise(
      n_obs = n(),
      rmse_log1p = sqrt(mean((log1p(pmax(pred_value, 0)) - log1p(pmax(value, 0)))^2, na.rm = TRUE)),
      rmse_raw = sqrt(mean((pred_value - value)^2, na.rm = TRUE)),
      obs_sd = stats::sd(value, na.rm = TRUE),
      scaled_rmse = rmse_raw / pmax(obs_sd, 1e-8),
      .groups = "drop"
    )

  n_metric <- metric_by_var %>% filter(variable == "N")
  r_metric <- metric_by_var %>% filter(variable == "R1")

  tibble(
    model_id = fit_row$model_id[[1]],
    alias = fit_row$alias[[1]],
    n_obs_compared = nrow(cmp_df),
    live_rmse_log1p = if (nrow(n_metric)) n_metric$rmse_log1p[[1]] else NA_real_,
    glucose_rmse_scaled = if (nrow(r_metric)) r_metric$scaled_rmse[[1]] else NA_real_,
    fit_quality_score = live_rmse_log1p + glucose_rmse_scaled
  )
}

make_overlay_for_line <- function(fit_rows, draw_map, stan_data, line_id, line_name, sim_times, out_path) {
  sim_df <- bind_rows(lapply(seq_len(nrow(fit_rows)), function(i) {
    draw_vec <- draw_map[[fit_rows$model_id[[i]]]]
    if (is.null(draw_vec)) {
      return(tibble())
    }
    simulate_line_from_draw(
      model_id = fit_rows$model_id[[i]],
      model_r_path = get_model_r_path("gpath", "v1"),
      draw_vec = draw_vec,
      stan_data = stan_data,
      line_id = line_id,
      times = sim_times,
      line_name = line_name,
      extra_cols = list(model_alias = fit_rows$alias[[i]])
    )
  })) %>%
    filter(variable %in% c("N", "R1"))

  obs_df <- build_obs_df_from_stan_data(stan_data, line_id = line_id, line_name = line_name) %>%
    filter(variable %in% c("N", "R1"))

  palette_vals <- setNames(
    grDevices::hcl.colors(nrow(fit_rows), palette = "Dark 3"),
    fit_rows$alias
  )

  plt <- plot_fit_overlays(
    sim_df = sim_df,
    obs_df = obs_df,
    color_by = "model_alias",
    color_values = palette_vals,
    color_label = "Model",
    title = line_name,
    subtitle = "Observed points with optimized full-data trajectories. Panels are alive cells (N) and glucose (R1).",
    print_plot = FALSE
  )

  ggsave(out_path, plot = plt, width = 12, height = 6, units = "in", dpi = 220)
  invisible(out_path)
}

write_decision_note <- function(path, decision_df, selected_model_id, overlay_manifest_df, existing_pareto_overlay_manifest_df) {
  selected_row <- decision_df %>% filter(model_id == selected_model_id)
  if (!nrow(selected_row)) {
    selected_sentence <- "No provisional representative model was selected from the available complete full-data fits."
  } else {
    selected_sentence <- sprintf(
      paste(
        "Provisional representative candidate: `%s` (`%s`).",
        "It has delta AIC %.3f, automated fit-quality score %.3f, mean transfer-minus-null holdout gain %.3f,",
        "NUTS QC status `%s`, and visual status `%s`."
      ),
      selected_row$model_id[[1]],
      selected_row$alias[[1]],
      selected_row$delta_AIC[[1]],
      selected_row$fit_quality_score[[1]],
      selected_row$mean_transfer_minus_null_holdout[[1]],
      selected_row$nuts_qc_status[[1]],
      selected_row$visual_fit_status[[1]]
    )
  }

  legacy_row <- decision_df %>% filter(model_id == "1R_1P_0W_C0_M1")
  legacy_sentence <- if (nrow(legacy_row)) {
    sprintf(
      "`1R_1P_0W_C0_M1` is no longer treated as the default representative in this WP5 output; its status is `%s`.",
      legacy_row$decision_status[[1]]
    )
  } else {
    "`1R_1P_0W_C0_M1` was not present among complete full-data fits in the supplied spec."
  }

  existing_overlay_sentence <- if (nrow(existing_pareto_overlay_manifest_df)) {
    sprintf(
      "Registered %d existing Pareto-front overlay files from `%s`; %d currently exist on disk.",
      nrow(existing_pareto_overlay_manifest_df),
      MODEL_ANALYSIS_PARETO_OVERLAY_DIR,
      sum(existing_pareto_overlay_manifest_df$file_exists, na.rm = TRUE)
    )
  } else {
    "No existing model_analysis Pareto overlay manifest was available."
  }

  overlay_sentence <- if (nrow(overlay_manifest_df)) {
    sprintf("Also generated %d WP5 candidate overlay files under `%s`.", nrow(overlay_manifest_df), FIGURE_DIR)
  } else {
    "No additional WP5 overlay files were generated."
  }

  lines <- c(
    "# WP5 mechanistic representative model decision note",
    "",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Scope",
    "",
    paste(
      "This report fixes the representative-model decision layer without changing the upstream gpath model code.",
      "It combines full-data optimisation AIC, transfer-minus-null holdout likelihood, optimisation stability,",
      "available oracle NUTS diagnostics, and automated alive/glucose fit-quality scores."
    ),
    "",
    "Visual credibility is a human-review criterion. This script records the current human review of the existing",
    "model_analysis Pareto-front overlays and still marks the selected model as provisional until the manuscript-level",
    "interpretation is finalized.",
    "",
    "## Provisional decision",
    "",
    selected_sentence,
    "",
    legacy_sentence,
    "",
    "## Generated review artifacts",
    "",
    existing_overlay_sentence,
    overlay_sentence,
    sprintf("- Decision table: `%s`", file.path(OUTPUT_DIR, "wp5_decision_table.csv")),
    sprintf("- Visual review table: `%s`", file.path(OUTPUT_DIR, "wp5_visual_review_table.csv")),
    sprintf("- Existing Pareto overlay manifest: `%s`", file.path(OUTPUT_DIR, "wp5_existing_pareto_overlay_manifest.csv")),
    sprintf("- Ranking/fit-quality tradeoff table: `%s`", file.path(OUTPUT_DIR, "wp5_model_ranking_tradeoff.csv")),
    sprintf("- Transfer gain table: `%s`", file.path(OUTPUT_DIR, "wp5_transfer_gain_by_case.csv")),
    sprintf("- Ranking/fit-quality figure: `%s`", file.path(FIGURE_DIR, "wp5_model_ranking_fit_quality_tradeoff.png")),
    "",
    "## Decision rule",
    "",
    paste(
      "A model is eligible only if it has a complete full-data optimisation row, finite automated fit metrics,",
      "acceptable human visual review from the existing model_analysis Pareto-front overlays,",
      "no obvious NUTS pathology in the preferred model_analysis oracle diagnostics, and is not the legacy 1R representative.",
      "It must also have at least one of: delta AIC <= 1000, positive mean transfer-minus-null holdout gain, or top-three",
      "automated fit quality. Among eligible models, the provisional selection minimizes a rank sum across delta AIC,",
      "automated fit quality, and transfer-minus-null holdout gain. Models without passing NUTS QC remain review candidates",
      "but are not selected."
    ),
    "",
    "## Caveats",
    "",
    paste(
      "The current canonical overlay helper exposes alive-cell burden and glucose trajectories, not a separate observed dead-cell",
      "trajectory. Dead counts are present in the Stan data object, but the active gpath likelihood and plotting helper do not",
      "currently produce a fitted dead-cell state for overlay. Strong WP5 claims should therefore be phrased around the fitted",
      "alive/glucose dynamics unless the model is extended to fit dead-cell observations explicitly."
    )
  )

  writeLines(lines, path, useBytes = TRUE)
}

stan_data <- add_group_structure(readRDS(resolve_stan_data_path(STAN_DATA_PATH)))
line_label_map <- make_line_label_map(stan_data$line_map)
obs_df_all <- build_obs_df_from_stan_data(stan_data)
sim_times <- sort(unique(c(stan_data$t_grid, obs_df_all$time)))

spec_summary_df <- summarize_spec(SPEC_PATH, stan_data)
write_table_pair(spec_summary_df, "wp5_optim_status_by_spec_row")

transfer_gain_df <- build_transfer_gain_table(spec_summary_df)
write_table_pair(transfer_gain_df, "wp5_transfer_gain_by_case")

transfer_summary_df <- if (nrow(transfer_gain_df)) {
  transfer_gain_df %>%
    group_by(model_id) %>%
    summarise(
      n_transfer_cases = n(),
      mean_transfer_minus_null_holdout = mean(transfer_minus_null_holdout, na.rm = TRUE),
      median_transfer_minus_null_holdout = median(transfer_minus_null_holdout, na.rm = TRUE),
      transfer_win_rate = mean(transfer_better_than_null, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  tibble(
    model_id = character(),
    n_transfer_cases = integer(),
    mean_transfer_minus_null_holdout = numeric(),
    median_transfer_minus_null_holdout = numeric(),
    transfer_win_rate = numeric()
  )
}
write_table_pair(transfer_summary_df, "wp5_transfer_gain_summary_by_model")

nuts_qc_df <- build_nuts_qc_from_model_analysis(MODEL_ANALYSIS_EXPORT_DIR)
if (!nrow(nuts_qc_df)) {
  nuts_qc_df <- build_nuts_qc(NUTS_JOB_CONFIG_PATH) %>%
    mutate(nuts_qc_source = NUTS_JOB_CONFIG_PATH)
}
write_table_pair(nuts_qc_df, "wp5_nuts_qc_summary")

model_analysis_screen_df <- build_model_analysis_screen_summary(MODEL_ANALYSIS_EXPORT_DIR)
write_table_pair(model_analysis_screen_df, "wp5_model_analysis_screen_summary")

full_fit_df <- spec_summary_df %>%
  filter(fit_type == "full", dataset_label == "gstarvation_v1", optim_complete, is.finite(best_lp)) %>%
  group_by(model_id) %>%
  slice_max(order_by = run_label, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    deviance = -2 * best_lp,
    AIC = 2 * k - 2 * best_lp,
    delta_AIC = AIC - min(AIC, na.rm = TRUE),
    prop_within_1 = n_within_1 / pmax(n_starts_total, 1),
    prop_within_5 = n_within_5 / pmax(n_starts_total, 1)
  ) %>%
  arrange(delta_AIC)

if (!nrow(full_fit_df)) {
  stop("No complete full-data optimisation rows were found for dataset_label == 'gstarvation_v1'.")
}

draw_map <- setNames(vector("list", nrow(full_fit_df)), full_fit_df$model_id)
for (i in seq_len(nrow(full_fit_df))) {
  draw_map[[full_fit_df$model_id[[i]]]] <- extract_full_fit_draw(full_fit_df[i, , drop = FALSE])
}

fit_quality_df <- bind_rows(lapply(seq_len(nrow(full_fit_df)), function(i) {
  build_fit_metric_for_model(
    fit_row = full_fit_df[i, , drop = FALSE],
    draw_vec = draw_map[[full_fit_df$model_id[[i]]]],
    stan_data = stan_data,
    obs_df = obs_df_all,
    sim_times = sim_times
  )
}))
write_table_pair(fit_quality_df, "wp5_fit_quality_summary_by_model")

pareto_front_df <- build_pareto_front(full_fit_df)
write_table_pair(
  pareto_front_df %>% select(model_id, alias, pareto_rank, k, deviance, AIC, delta_AIC),
  "wp5_model_analysis_style_pareto_front"
)

existing_pareto_overlay_manifest_df <- build_existing_pareto_overlay_manifest(
  overlay_dir = MODEL_ANALYSIS_PARETO_OVERLAY_DIR,
  stan_data = stan_data,
  pareto_model_ids = pareto_front_df$model_id
)
write_table_pair(existing_pareto_overlay_manifest_df, "wp5_existing_pareto_overlay_manifest")

visual_review_df <- build_visual_review_table(pareto_front_df)
write_table_pair(visual_review_df, "wp5_visual_review_table")

ranking_df <- full_fit_df %>%
  select(
    model_id, alias, run_label, k, n_starts_total, n_starts_valid,
    n_within_1, n_within_5, prop_within_1, prop_within_5,
    best_lp, deviance, AIC, delta_AIC
  ) %>%
  left_join(transfer_summary_df, by = "model_id") %>%
  left_join(nuts_qc_df, by = "model_id") %>%
  left_join(fit_quality_df, by = c("model_id", "alias")) %>%
  left_join(model_analysis_screen_df, by = "model_id") %>%
  left_join(visual_review_df, by = c("model_id", "alias")) %>%
  mutate(
    visual_fit_status = if_else(is.na(visual_fit_status), "not_reviewed_in_existing_pareto_overlays", visual_fit_status),
    visual_review_source = if_else(is.na(visual_review_source), NA_character_, visual_review_source),
    visual_review_note = if_else(is.na(visual_review_note), NA_character_, visual_review_note),
    on_model_analysis_pareto_front = model_id %in% pareto_front_df$model_id,
    nuts_qc_status = case_when(
      nuts_available %in% TRUE & nuts_qc_pass %in% TRUE ~ "pass",
      nuts_available %in% TRUE & !is.na(nuts_qc_pass) & !nuts_qc_pass ~ "pathology",
      TRUE ~ "not_available"
    ),
    aic_rank = rank(delta_AIC, ties.method = "min", na.last = "keep"),
    fit_quality_rank = rank(fit_quality_score, ties.method = "min", na.last = "keep"),
    transfer_rank = rank(-mean_transfer_minus_null_holdout, ties.method = "min", na.last = "keep"),
    rank_sum = aic_rank + fit_quality_rank + transfer_rank
  ) %>%
  arrange(rank_sum, delta_AIC)

write_table_pair(ranking_df, "wp5_model_ranking_tradeoff")

eligible_df <- ranking_df %>%
  filter(
    model_id != "1R_1P_0W_C0_M1",
    visual_fit_status == "acceptable_by_human_review",
    on_model_analysis_pareto_front,
    is.finite(rank_sum),
    is.finite(fit_quality_score),
    nuts_qc_status == "pass",
    delta_AIC <= 1000 |
      mean_transfer_minus_null_holdout > 0 |
      fit_quality_rank <= 3
  ) %>%
  arrange(rank_sum, delta_AIC)

selected_model_id <- if (nrow(eligible_df)) eligible_df$model_id[[1]] else NA_character_

candidate_ids <- unique_preserve(c(
  selected_model_id,
  pareto_front_df$model_id,
  ranking_df %>% arrange(fit_quality_score) %>% pull(model_id) %>% head(1),
  ranking_df %>% arrange(desc(mean_transfer_minus_null_holdout)) %>% pull(model_id) %>% head(1)
))
candidate_ids <- candidate_ids[!is.na(candidate_ids) & candidate_ids %in% full_fit_df$model_id]
candidate_ids <- unique_preserve(candidate_ids)

overlay_fit_rows <- ranking_df %>%
  filter(model_id %in% candidate_ids) %>%
  arrange(match(model_id, candidate_ids))

overlay_manifest_df <- bind_rows(lapply(sort(unique(as.integer(stan_data$line_id))), function(line_id) {
  line_name <- unname(line_label_map[as.character(line_id)])
  if (is.na(line_name) || !nzchar(line_name)) {
    line_name <- sprintf("line_%02d", line_id)
  }
  safe_name <- gsub("[^A-Za-z0-9_+-]+", "_", line_name)
  out_path <- file.path(FIGURE_DIR, sprintf("wp5_candidate_fit_overlay_line%02d_%s.png", line_id, safe_name))
  make_overlay_for_line(
    fit_rows = overlay_fit_rows,
    draw_map = draw_map,
    stan_data = stan_data,
    line_id = line_id,
    line_name = line_name,
    sim_times = sim_times,
    out_path = out_path
  )
  tibble(
    line_id = line_id,
    line_name = line_name,
    output_path = out_path
  )
}))
write_table_pair(overlay_manifest_df, "wp5_candidate_overlay_manifest")

decision_df <- ranking_df %>%
  mutate(
    is_candidate_overlay = model_id %in% candidate_ids,
    decision_status = case_when(
      model_id == selected_model_id ~ "provisional_representative_visual_review_passed",
      model_id == "1R_1P_0W_C0_M1" ~ "legacy_default_rejected_for_wp5_review",
      on_model_analysis_pareto_front & visual_fit_status == "acceptable_by_human_review" ~ "candidate_accepted_visual_review",
      on_model_analysis_pareto_front & visual_fit_status != "acceptable_by_human_review" ~ "rejected_visual_review",
      nuts_qc_status == "pathology" ~ "rejected_nuts_pathology",
      !is.finite(fit_quality_score) ~ "rejected_missing_fit_quality",
      model_id %in% candidate_ids ~ "candidate_visual_review",
      TRUE ~ "not_shortlisted"
    ),
    rejection_reason = case_when(
      model_id == selected_model_id ~ NA_character_,
      model_id == "1R_1P_0W_C0_M1" ~ "Legacy 1R representative is not allowed to remain the default without visual-fit support.",
      on_model_analysis_pareto_front & visual_fit_status == "acceptable_by_human_review" ~ NA_character_,
      on_model_analysis_pareto_front & visual_fit_status != "acceptable_by_human_review" ~ "Existing model_analysis Pareto-front overlays did not pass human visual review.",
      nuts_qc_status == "pathology" ~ "Oracle NUTS diagnostics have R-hat or ESS pathologies under WP5 thresholds.",
      !is.finite(fit_quality_score) ~ "Could not compute automated alive/glucose fit-quality metrics.",
      model_id %in% candidate_ids ~ NA_character_,
      TRUE ~ "Lower joint ranking across AIC, fit quality, transfer gain, and QC than shortlisted candidates."
    )
  ) %>%
  arrange(match(decision_status, c(
    "provisional_representative_visual_review_passed",
    "candidate_accepted_visual_review",
    "candidate_visual_review",
    "legacy_default_rejected_for_wp5_review",
    "rejected_visual_review",
    "rejected_nuts_pathology",
    "rejected_missing_fit_quality",
    "not_shortlisted"
  )), rank_sum)
write_table_pair(decision_df, "wp5_decision_table")

plot_df <- ranking_df %>%
  mutate(
    alias = factor(alias, levels = rev(unique_preserve(ranking_df$alias[order(ranking_df$delta_AIC)]))),
    decision_group = case_when(
      model_id == selected_model_id ~ "Provisional representative",
      model_id == "1R_1P_0W_C0_M1" ~ "Legacy default",
      model_id %in% candidate_ids ~ "Candidate overlay",
      TRUE ~ "Other assessment model"
    )
  )

tradeoff_plot <- ggplot(plot_df, aes(delta_AIC, fit_quality_score)) +
  geom_point(aes(size = pmax(n_within_1, 0), color = decision_group), alpha = 0.9) +
  ggrepel::geom_label_repel(
    aes(label = alias),
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    size = 2.5
  ) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_size_continuous(range = c(1.8, 6), name = "starts within 1 log unit") +
  scale_color_manual(
    values = c(
      "Provisional representative" = "#B2182B",
      "Candidate overlay" = "#2166AC",
      "Legacy default" = "#666666",
      "Other assessment model" = "#999999"
    )
  ) +
  labs(
    title = "WP5 model ranking versus automated fit quality",
    x = expression(Delta ~ AIC),
    y = "alive/glucose fit-quality score (lower is better)",
    color = NULL,
    caption = "Fit quality is an automated score from optimized full-data trajectories at observed alive-cell and glucose time points; final WP5 eligibility still requires human visual review of overlays."
  ) +
  wp5_theme(base_size = 8)
save_plot_pair(tradeoff_plot, "wp5_model_ranking_fit_quality_tradeoff", width = 7.2, height = 5.2)

if (nrow(transfer_summary_df)) {
  transfer_plot <- transfer_summary_df %>%
    left_join(full_fit_df %>% select(model_id, alias, delta_AIC), by = "model_id") %>%
    mutate(alias = factor(alias, levels = alias[order(mean_transfer_minus_null_holdout)])) %>%
    ggplot(aes(alias, mean_transfer_minus_null_holdout, fill = mean_transfer_minus_null_holdout > 0)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
    labs(
      title = "WP5 transfer-minus-null holdout gain by model",
      x = NULL,
      y = "mean holdout log-likelihood gain"
    ) +
    wp5_theme(base_size = 8)
  save_plot_pair(transfer_plot, "wp5_transfer_gain_by_model", width = 6.8, height = 4.8)
}

write_decision_note(
  path = file.path(OUTPUT_DIR, "wp5_mechanistic_representative_decision_note.md"),
  decision_df = decision_df,
  selected_model_id = selected_model_id,
  overlay_manifest_df = overlay_manifest_df,
  existing_pareto_overlay_manifest_df = existing_pareto_overlay_manifest_df
)

cat(sprintf("Wrote WP5 tables to %s\n", OUTPUT_DIR))
cat(sprintf("Wrote WP5 figures to %s\n", FIGURE_DIR))
if (!is.na(selected_model_id)) {
  cat(sprintf("Provisional WP5 representative candidate: %s\n", selected_model_id))
} else {
  cat("No provisional WP5 representative candidate passed the automated eligibility filter.\n")
}
