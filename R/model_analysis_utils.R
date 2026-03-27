# Helper utilities for workflow/model_analysis.Rmd.
# NOTE: This file is currently too bloated and should be streamlined in a refactor ASAP.

safe_read_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

announce_missing <- function(...) {
  cat("SKIP:", sprintf(...), "\n")
}

ensure_report_export_root <- function(path = project_paths$report_export_root) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_audit_artifact <- function(x, stem, dir_path = project_paths$report_export_root) {
  dir_path <- ensure_report_export_root(dir_path)
  rds_path <- file.path(dir_path, sprintf("%s.Rds", stem))
  csv_path <- file.path(dir_path, sprintf("%s.csv", stem))
  saveRDS(x, rds_path)

  if (is.data.frame(x) || tibble::is_tibble(x)) {
    utils::write.csv(as.data.frame(x), csv_path, row.names = FALSE)
  } else {
    csv_path <- NA_character_
  }

  tibble(
    artifact = stem,
    rds_path = normalizePath(rds_path, winslash = "/", mustWork = FALSE),
    csv_path = if (is.na(csv_path)) NA_character_ else normalizePath(csv_path, winslash = "/", mustWork = FALSE)
  )
}

preview_table <- function(x, n = 20L) {
  if (is.null(x) || !nrow(x)) {
    return(x)
  }
  head(as.data.frame(x), n)
}

unique_preserve <- function(x) {
  x[!duplicated(x)]
}

# Editable model-ID alias conventions for plotting/reporting.
# Users can modify these labels without changing downstream alias logic.
model_alias_conventions <- list(
  separator = ",",
  family_alias = list(
    `1R_1P` = "1R",
    `2R_1P` = "2R(a)",
    `2R_2P_C0` = "2R(f)",
    `2R_2P_C1` = "2R(l)"
  ),
  waste_alias = list(
    `0W` = "",
    `1W_M0` = "W(m)",
    `1W_M1` = "W(a)"
  )
)

parse_run_id_for_alias <- function(run_id) {
  if (exists("parse_run_id", mode = "function")) {
    return(parse_run_id(run_id))
  }

  m <- regexec("^(\\d+)R_(\\d+)P_(\\d+)W_C(\\d+)_M(\\d+)$", run_id)
  p <- regmatches(run_id, m)[[1]]
  if (length(p) != 6) {
    stop(sprintf("Could not parse run_id: %s", run_id))
  }
  list(
    R = as.integer(p[2]),
    P = as.integer(p[3]),
    W = as.integer(p[4]),
    C = as.integer(p[5]),
    M = as.integer(p[6])
  )
}

build_model_alias <- function(run_id, conventions = model_alias_conventions, format = c("text", "plotmath")) {
  format <- match.arg(format)
  dims <- tryCatch(parse_run_id_for_alias(run_id), error = function(e) NULL)
  if (is.null(dims)) {
    return(run_id)
  }

  family_key <- if (dims$R == 1 && dims$P == 1) {
    "1R_1P"
  } else if (dims$R == 2 && dims$P == 1) {
    "2R_1P"
  } else if (dims$R == 2 && dims$P == 2 && dims$C == 0) {
    "2R_2P_C0"
  } else if (dims$R == 2 && dims$P == 2 && dims$C == 1) {
    "2R_2P_C1"
  } else {
    NA_character_
  }

  waste_key <- if (dims$W == 0) {
    "0W"
  } else if (dims$W == 1 && dims$M == 0) {
    "1W_M0"
  } else if (dims$W == 1 && dims$M == 1) {
    "1W_M1"
  } else {
    NA_character_
  }

  family_alias <- conventions$family_alias[[family_key]]
  waste_alias <- conventions$waste_alias[[waste_key]]

  if (is.null(family_alias) || !nzchar(family_alias)) {
    family_alias <- run_id
  }
  if (is.null(waste_alias)) {
    waste_alias <- ""
  }

  sep <- conventions$separator
  if (is.null(sep) || !nzchar(sep)) {
    sep <- ","
  }

  alias_text <- if (nzchar(waste_alias)) {
    paste(family_alias, waste_alias, sep = sep)
  } else {
    family_alias
  }

  if (identical(format, "text")) {
    return(alias_text)
  }

  alias_text
}

ranking_model_levels <- function(ranking_df, reverse = FALSE) {
  levels <- unique_preserve(ranking_df$model_id)
  if (reverse) rev(levels) else levels
}

ranking_model_order <- function(x, ranking_df) {
  match(x, ranking_model_levels(ranking_df))
}

make_line_label_map <- function(line_map) {
  if (is.null(line_map) || !length(line_map)) {
    return(setNames(character(0), character(0)))
  }
  ids <- as.integer(unname(line_map))
  labels <- names(line_map)
  out <- setNames(labels, as.character(ids))
  out[order(as.integer(names(out)))]
}

label_line_ids <- function(line_ids, line_label_map, as_factor = TRUE) {
  line_ids <- as.integer(line_ids)
  labels <- unname(line_label_map[as.character(line_ids)])
  fallback <- sprintf("Line %d", line_ids)
  labels[is.na(labels) | !nzchar(labels)] <- fallback[is.na(labels) | !nzchar(labels)]
  if (!as_factor) {
    return(labels)
  }

  known_ids <- suppressWarnings(as.integer(names(line_label_map)))
  known_ids <- known_ids[is.finite(known_ids)]
  level_ids <- sort(unique(c(known_ids, line_ids)))
  level_labels <- unname(line_label_map[as.character(level_ids)])
  level_fallback <- sprintf("Line %d", level_ids)
  level_labels[is.na(level_labels) | !nzchar(level_labels)] <- level_fallback[is.na(level_labels) | !nzchar(level_labels)]
  factor(labels, levels = unique_preserve(level_labels))
}

build_density_ridgeline_df <- function(
  draw_df,
  y_var,
  y_levels = NULL,
  facet_vars = c("line_label", "parameter"),
  ploidy_offsets = c("Low ploidy" = -0.18, "High ploidy" = 0.18),
  ridge_height = 0.34,
  density_n = 256L
) {
  if (is.null(draw_df) || !nrow(draw_df)) {
    return(tibble())
  }

  req_cols <- c(y_var, facet_vars, "ploidy_label", "value")
  missing_cols <- setdiff(req_cols, names(draw_df))
  if (length(missing_cols)) {
    stop(sprintf("Missing ridgeline columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (is.null(y_levels)) {
    y_levels <- unique_preserve(as.character(draw_df[[y_var]]))
  }
  y_map <- tibble(
    y_value = y_levels,
    y_index = seq_along(y_levels)
  )

  facet_bounds_df <- draw_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
    dplyr::summarise(
      x_min = min(.data$value, na.rm = TRUE),
      x_max = max(.data$value, na.rm = TRUE),
      .groups = "drop"
    )

  density_df <- draw_df %>%
    dplyr::mutate(.y_value = as.character(.data[[y_var]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, ".y_value", "ploidy_label")))) %>%
    dplyr::group_modify(function(.x, .y) {
      bounds <- facet_bounds_df
      for (nm in facet_vars) {
        bounds <- bounds[bounds[[nm]] == .y[[nm]][[1]], , drop = FALSE]
      }
      vals <- .x$value[is.finite(.x$value)]
      if (!length(vals) || !nrow(bounds)) {
        return(tibble())
      }
      dens <- stats::density(
        vals,
        from = bounds$x_min[[1]],
        to = bounds$x_max[[1]],
        n = density_n
      )
      tibble(value = dens$x, density = dens$y)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(y_map, by = c(".y_value" = "y_value")) %>%
    dplyr::mutate(
      ploidy_offset = unname(ploidy_offsets[as.character(.data$ploidy_label)]),
      ploidy_offset = dplyr::if_else(is.na(.data$ploidy_offset), 0, .data$ploidy_offset)
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, ".y_value", "ploidy_label")))) %>%
    dplyr::mutate(
      density_scaled = ridge_height * .data$density / max(.data$density, na.rm = TRUE),
      y_base = .data$y_index + .data$ploidy_offset,
      y_top = .data$y_base + .data$density_scaled
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(y_label = ".y_value")

  density_df
}

artifact_csv_path <- function(stem, dir_path = project_paths$report_export_root) {
  file.path(ensure_report_export_root(dir_path), sprintf("%s.csv", stem))
}

emit_artifact_note <- function(stem, label = "Audit CSV") {
  knitr::asis_output(sprintf("%s: `%s`\n", label, artifact_csv_path(stem)))
}

write_agent_text_report <- function(
  path = project_paths$agent_text_path,
  ranking_df,
  support_set_df = NULL,
  shortlist_df = NULL,
  provenance_df = NULL,
  metric_glossary_df = NULL,
  selected_model_summary_df = NULL,
  raw_theta_agreement_df = NULL,
  focused_hit_df = NULL,
  report_export_root = project_paths$report_export_root
) {
  lines <- c(
    "Model Analysis Summary",
    "",
    sprintf("Generated: %s", Sys.time()),
    sprintf("Stan data: %s", project_paths$stan_data),
    sprintf("Audit export root: %s", ensure_report_export_root(report_export_root)),
    ""
  )

  if (!is.null(selected_model_summary_df) && nrow(selected_model_summary_df)) {
    x <- selected_model_summary_df[1, , drop = FALSE]
    lines <- c(
      lines,
      "Representative Model",
      sprintf("model_id: %s", x$model_id[[1]]),
      sprintf("alias: %s", x$alias[[1]]),
      sprintf("selection_basis: best mean rank within the aggregated shortlist; used for single-model plots only"),
      sprintf("delta_AIC: %.3f", x$delta_AIC[[1]]),
      sprintf("prop_within_1: %.3f", x$prop_within_1[[1]]),
      sprintf("prop_within_5: %.3f", x$prop_within_5[[1]]),
      sprintf("directional_agreement: %s", format(round(x$prop_directional_agreement_model[[1]], 3), nsmall = 3)),
      sprintf("mean_transfer_gain: %s", format(round(x$mean_transfer_gain[[1]], 3), nsmall = 3)),
      sprintf("mean_parameter_log_error_improvement: %s", format(round(x$mean_parameter_log_error_improvement[[1]], 3), nsmall = 3)),
      sprintf("parameter_win_rate: %s", format(round(x$mean_prop_parameters_better[[1]], 3), nsmall = 3)),
      sprintf("mean_rank: %s", format(round(x$mean_rank[[1]], 3), nsmall = 3)),
      ""
    )
  }

  if (!is.null(shortlist_df) && nrow(shortlist_df)) {
    lines <- c(
      lines,
      "Aggregated Shortlist",
      capture.output(print(as.data.frame(shortlist_df), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(support_set_df) && nrow(support_set_df)) {
    lines <- c(
      lines,
      "AIC Support Set",
      capture.output(print(as.data.frame(support_set_df), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(provenance_df) && nrow(provenance_df)) {
    lines <- c(
      lines,
      "Current Provenance",
      capture.output(print(as.data.frame(provenance_df), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(ranking_df) && nrow(ranking_df)) {
    lines <- c(
      lines,
      "Model Ranking",
      capture.output(print(as.data.frame(ranking_df %>% select(model_id, alias, best_lp, k, AIC, delta_AIC)), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(raw_theta_agreement_df) && nrow(raw_theta_agreement_df)) {
    lines <- c(
      lines,
      "Directional Agreement Of raw_theta_ploidy With Global Fit",
      capture.output(print(as.data.frame(raw_theta_agreement_df), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(focused_hit_df) && nrow(focused_hit_df)) {
    focused_summary_df <- focused_hit_df %>%
      group_by(parameter) %>%
      summarise(
        mean_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
        mean_win_rate = mean(prop_transfer_better_than_null, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_log_error_improvement), desc(mean_win_rate))
    lines <- c(
      lines,
      "Focused Parameter Summary",
      capture.output(print(as.data.frame(focused_summary_df), row.names = FALSE)),
      ""
    )
  }

  if (!is.null(metric_glossary_df) && nrow(metric_glossary_df)) {
    lines <- c(
      lines,
      "Metric Glossary",
      capture.output(print(as.data.frame(metric_glossary_df), row.names = FALSE)),
      ""
    )
  }

  lines <- c(
    lines,
    "Interpretation Verdict",
    if (!is.null(shortlist_df) && nrow(shortlist_df)) {
      pos_pred <- shortlist_df %>% filter(!is.na(mean_transfer_gain), mean_transfer_gain > 0)
      pos_param <- shortlist_df %>% filter(!is.na(mean_parameter_log_error_improvement), mean_parameter_log_error_improvement > 0)
      c(
        sprintf("shortlist_models: %s", paste(shortlist_df$model_id, collapse = ", ")),
        if (!nrow(pos_pred) && !nrow(pos_param)) {
          "verdict: no shortlisted model is positive on average for either predictive transfer or parameter transfer."
        } else if (nrow(pos_pred) && !nrow(pos_param)) {
          "verdict: some shortlisted models improve predictive transfer, but none are positive on average for parameter transfer."
        } else if (!nrow(pos_pred) && nrow(pos_param)) {
          "verdict: some shortlisted models improve parameter transfer, but none are positive on average for predictive transfer."
        } else {
          "verdict: at least some shortlisted models are positive on both predictive and parameter-transfer summaries."
        }
      )
    } else {
      "verdict: no shortlist available."
    },
    "",
    "Key Audit Artifacts",
    sprintf("- %s", artifact_csv_path("raw_theta_ploidy_df", report_export_root)),
    sprintf("- %s", artifact_csv_path("raw_theta_ploidy_transfer_df", report_export_root)),
    sprintf("- %s", artifact_csv_path("focused_parameter_df", report_export_root)),
    sprintf("- %s", artifact_csv_path("focused_parameter_hit_df", report_export_root))
  )

  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

build_stability_summary <- function(ranking_df) {
  if (!nrow(ranking_df)) {
    return(tibble())
  }

  purrr::map_dfr(seq_len(nrow(ranking_df)), function(i) {
    run_dir <- ranking_df$run_dir[i]
    lp <- readRDS(file.path(run_dir, "optim_lp_all.Rds"))
    nbhd <- summarize_lp_neighborhood(lp, near_cut = 5)

    tibble(
      model_id = ranking_df$model_id[i],
      alias = ranking_df$alias[i],
      run_label = ranking_df$run_label[i],
      n_total = length(lp),
      n_finite = sum(nbhd$finite_mask),
      best_lp = nbhd$best_lp,
      n_within_1 = sum(nbhd$finite_mask & nbhd$lp_diff < 1),
      n_within_5 = sum(nbhd$finite_mask & nbhd$lp_diff < 5),
      n_within_10 = sum(nbhd$finite_mask & nbhd$lp_diff < 10),
      prop_within_1 = mean(nbhd$finite_mask & nbhd$lp_diff < 1),
      prop_within_5 = mean(nbhd$finite_mask & nbhd$lp_diff < 5),
      prop_within_10 = mean(nbhd$finite_mask & nbhd$lp_diff < 10),
      delta_AIC = ranking_df$delta_AIC[i]
    )
  }) %>%
    arrange(desc(prop_within_5), desc(prop_within_1), delta_AIC)
}

build_support_set <- function(ranking_df, stability_df, delta_aic_cutoff = 10) {
  if (!nrow(ranking_df)) {
    return(tibble())
  }

  ranking_df %>%
    left_join(
      stability_df %>% select(model_id, prop_within_1, prop_within_5),
      by = "model_id"
    ) %>%
    filter(delta_AIC <= delta_aic_cutoff) %>%
    select(model_id, alias, delta_AIC, k, prop_within_1, prop_within_5) %>%
    arrange(delta_AIC, desc(prop_within_1), desc(prop_within_5))
}

build_aggregated_shortlist <- function(model_screen_df, top_n_mean_rank = 3L, min_prop_within_1 = 0.05) {
  if (!nrow(model_screen_df)) {
    return(tibble())
  }

  stable_df <- model_screen_df %>%
    filter(!is.na(prop_within_1), prop_within_1 >= min_prop_within_1)

  if (!nrow(stable_df)) {
    stable_df <- model_screen_df
  }

  stable_df %>%
    arrange(mean_rank, delta_AIC) %>%
    slice_head(n = top_n_mean_rank) %>%
    mutate(shortlist_source = sprintf("top_mean_rank_after_prop_within_1_filter>=%.2f", min_prop_within_1)) %>%
    arrange(mean_rank, delta_AIC)
}

build_metric_glossary <- function() {
  tibble::tribble(
    ~metric, ~interpretation,
    "AIC", "Lower is better among trusted optimisation fits; combines fit and parameter count.",
    "delta_AIC", "AIC minus the best AIC in the ranked set; values near zero indicate comparable support.",
    "prop_within_1", "Fraction of optimisation starts within 1 log-posterior unit of the best start; higher means a sharper and more trustworthy optimum neighborhood.",
    "prop_within_5", "Fraction of optimisation starts within 5 log-posterior units of the best start; higher means more stable optimisation.",
    "prop_directional_agreement_model", "Across the common raw_theta_ploidy parameters, fraction of transfer fits whose sign matches the global all-data fit for that model.",
    "transfer_gain", "Held-out predictive improvement of transfer over null; positive favors transfer.",
    "transfer_regret", "Held-out predictive gap between oracle and transfer; smaller is better.",
    "mean_log_error_improvement", "For reconstructed held-out parameters, positive means transfer is closer to oracle than null on average.",
    "prop_parameters_transfer_better_than_null", "Fraction of reconstructed held-out parameters for which transfer beats null; above 0.5 favors transfer."
  )
}

build_current_provenance_table <- function(
  ranking_df,
  selected_model_id,
  project_paths,
  stan_data_command = "Rscript R/prepare_data.R"
) {
  if (!nrow(ranking_df) || !"model_id" %in% names(ranking_df)) {
    selected_row <- tibble()
  } else {
    selected_row <- ranking_df %>% filter(model_id == selected_model_id) %>% slice_head(n = 1)
  }
  selected_run_dir <- if (nrow(selected_row)) selected_row$run_dir[[1]] else NA_character_
  selected_manifest <- if (nrow(selected_row)) selected_row$manifest_path[[1]] else NA_character_

  tibble::tribble(
    ~artifact_group, ~path, ~provenance, ~status,
    "Stan data", project_paths$stan_data,
    stan_data_command,
    "Confirmed on 2026-03-20 by byte-identical rerun.",
    "Assessment model set", project_paths$assessment_run_ids,
    "Current ranking is restricted to these model IDs.",
    "Concrete input file.",
    "Selected optimisation run", selected_run_dir,
    if (!is.na(selected_manifest)) sprintf("Manifest: %s", selected_manifest) else "No selected optimisation run available.",
    if (!is.na(selected_run_dir)) "Concrete current source for selected-model ranking and best-fit coefficients." else "Missing.",
    "Selected optimisation transfer summary", file.path(project_paths$transfer_optim_root, selected_model_id, "transfer_comparison_summary.Rds"),
    "Loaded by the predictive transfer summary blocks for the selected model.",
    "Concrete current source when present.",
    "Selected parameter transfer summary", file.path(project_paths$transfer_optim_root, selected_model_id, "parameter_transfer_summary.Rds"),
    "Loaded by the parameter-transfer interpretation blocks for the selected model.",
    "Concrete current source when present.",
    "Selected NUTS transfer summary", file.path(project_paths$transfer_nuts_root, selected_model_id, "transfer_comparison_summary.Rds"),
    "Optional comparison source for posterior-averaged transfer checks.",
    "Concrete comparison path when present."
  )
}

summarize_model_screen <- function(ranking_df, stability_df, directional_agreement_model_df, screen_transfer_df, screen_parameter_df) {
  base <- ranking_df %>%
    select(model_id, alias, delta_AIC, k, best_lp)

  if (nrow(stability_df)) {
    base <- base %>%
      left_join(
        stability_df %>% select(model_id, prop_within_1, prop_within_5, n_within_1, n_within_5),
        by = "model_id"
      )
  }

  if (nrow(directional_agreement_model_df)) {
    base <- base %>%
      left_join(
        directional_agreement_model_df %>%
          select(model_id, prop_directional_agreement_model = prop_directional_agreement, mean_abs_shift_from_global),
        by = "model_id"
      )
  }

  if (nrow(screen_transfer_df)) {
    transfer_summary <- screen_transfer_df %>%
      group_by(model_id) %>%
      summarise(
        mean_transfer_gain = mean(transfer_gain, na.rm = TRUE),
        median_transfer_gain = median(transfer_gain, na.rm = TRUE),
        mean_transfer_regret = mean(transfer_regret, na.rm = TRUE),
        prop_transfer_better_than_null = mean(transfer > null, na.rm = TRUE),
        n_transfer_cases = dplyr::n(),
        .groups = "drop"
      )
    base <- base %>% left_join(transfer_summary, by = "model_id")
  }

  if (nrow(screen_parameter_df)) {
    parameter_summary <- screen_parameter_df %>%
      mutate(
        mean_log_error_improvement = mean_abs_log_null_holdout_ratio - mean_abs_log_holdout_ratio
      ) %>%
      group_by(model_id) %>%
      summarise(
        mean_parameter_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
        median_parameter_log_error_improvement = median(mean_log_error_improvement, na.rm = TRUE),
        mean_prop_parameters_better = mean(prop_parameters_transfer_better_than_null, na.rm = TRUE),
        n_parameter_cases = dplyr::n(),
        .groups = "drop"
      )
    base <- base %>% left_join(parameter_summary, by = "model_id")
  }

  base %>%
    mutate(
      rank_delta_AIC = rank(delta_AIC, ties.method = "min", na.last = "keep"),
      rank_stability = rank(-prop_within_1, ties.method = "min", na.last = "keep"),
      rank_directional_agreement = rank(-prop_directional_agreement_model, ties.method = "min", na.last = "keep"),
      rank_predictive_transfer = rank(-mean_transfer_gain, ties.method = "min", na.last = "keep"),
      rank_parameter_transfer = rank(-mean_parameter_log_error_improvement, ties.method = "min", na.last = "keep")
    ) %>%
    rowwise() %>%
    mutate(mean_rank = mean(c_across(starts_with("rank_")), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(mean_rank, delta_AIC, desc(prop_within_1), desc(mean_parameter_log_error_improvement))
}

build_selected_model_takeaway <- function(selected_model_summary_df) {
  if (is.null(selected_model_summary_df) || !nrow(selected_model_summary_df)) {
    return("No selected model summary could be built.")
  }

  x <- selected_model_summary_df[1, , drop = FALSE]
  predictive_phrase <- dplyr::case_when(
    is.na(x$mean_transfer_gain[[1]]) ~ "predictive transfer is unavailable",
    x$mean_transfer_gain[[1]] > 0 ~ "predictive transfer beats the null baseline on average",
    TRUE ~ "predictive transfer does not beat the null baseline on average"
  )
  parameter_phrase <- dplyr::case_when(
    is.na(x$mean_parameter_log_error_improvement[[1]]) ~ "parameter-transfer evidence is unavailable",
    x$mean_parameter_log_error_improvement[[1]] > 0 ~ "parameter transfer moves held-out parameters closer to oracle on average",
    TRUE ~ "parameter transfer does not move held-out parameters closer to oracle on average"
  )

  sprintf(
    "Selected model %s (%s) has delta_AIC %.3f, prop_within_1 %.3f, prop_within_5 %.3f, directional agreement %.3f, and suggests that %s while %s.",
    x$model_id[[1]],
    x$alias[[1]],
    x$delta_AIC[[1]],
    x$prop_within_1[[1]],
    x$prop_within_5[[1]],
    x$prop_directional_agreement_model[[1]],
    predictive_phrase,
    parameter_phrase
  )
}

read_run_id_file <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }

  vals <- readLines(path, warn = FALSE)
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  vals <- vals[!grepl("^#", vals)]
  vals
}

discover_model_runs <- function(base_path) {
  if (!dir.exists(base_path)) {
    return(character(0))
  }
  run_dirs <- list.dirs(base_path, recursive = TRUE, full.names = TRUE)
  run_dirs[
    file.exists(file.path(run_dirs, "optim_lp_all.Rds")) &
    file.exists(file.path(run_dirs, "optim_draws_all.Rds")) &
    file.exists(file.path(run_dirs, "run_manifest.json"))
  ]
}

build_model_ranking <- function(
  base_path = project_paths$optim_root,
  stan_data_path = project_paths$stan_data,
  assessment_ids = read_run_id_file(project_paths$assessment_run_ids),
  dataset_label_filter = project_paths$report_dataset_label
) {
  run_dirs <- discover_model_runs(base_path)
  if (!length(run_dirs)) {
    return(tibble())
  }

  stan_data <- readRDS(stan_data_path)

  rows <- lapply(run_dirs, function(run_dir) {
    parts <- strsplit(normalizePath(run_dir, winslash = "/"), "/", fixed = TRUE)[[1]]
    n <- length(parts)
    run_label <- parts[n]
    run_id <- parts[n - 1]
    dataset_label <- parts[n - 2]

    lp <- readRDS(file.path(run_dir, "optim_lp_all.Rds"))
    lp <- lp[is.finite(lp)]
    if (!length(lp)) {
      return(NULL)
    }

    manifest_path <- file.path(run_dir, "run_manifest.json")
    manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)

    tibble(
      model_id = run_id,
      dataset_label = dataset_label,
      run_label = run_label,
      run_dir = run_dir,
      manifest_path = manifest_path,
      best_lp = max(lp),
      k = get_hierarchical_k(run_id, N_lines = stan_data$N_lines),
      manifest_stan_data = manifest$stan_data_path,
      manifest_stan_file = manifest$stan_file_path,
      alias = build_model_alias(run_id, format = "text"),
      alias_plot = build_model_alias(run_id, format = "plotmath")
    )
  })

  out <- bind_rows(rows) %>%
    filter(model_id %in% assessment_ids)

  if (!is.null(dataset_label_filter)) {
    out <- out %>% filter(dataset_label %in% dataset_label_filter)
  }

  out %>%
    group_by(model_id, alias, alias_plot, dataset_label, k) %>%
    slice_max(order_by = run_label, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      deviance = -2 * best_lp,
      AIC = 2 * k - 2 * best_lp,
      delta_AIC = AIC - min(AIC)
    ) %>%
    arrange(delta_AIC)
}

summarize_manifest_status <- function(ranking_df) {
  if (!nrow(ranking_df)) {
    return(tibble())
  }

  ranking_df %>%
    transmute(
      model_id,
      run_label,
      run_dir,
      manifest_present = file.exists(manifest_path),
      optim_complete = file.exists(file.path(run_dir, "optim_lp_all.Rds")) &
        file.exists(file.path(run_dir, "optim_draws_all.Rds")),
      stan_data_matches_current = normalizePath(manifest_stan_data, winslash = "/", mustWork = FALSE) ==
        normalizePath(project_paths$stan_data, winslash = "/", mustWork = FALSE)
    ) %>%
    arrange(model_id)
}

find_selected_model <- function(ranking_df) {
  if (!nrow(ranking_df)) {
    return(NA_character_)
  }

  ranking_df$model_id[which.min(ranking_df$delta_AIC)]
}

extract_best_draw_vector <- function(run_dir) {
  lp <- readRDS(file.path(run_dir, "optim_lp_all.Rds"))
  draws <- readRDS(file.path(run_dir, "optim_draws_all.Rds"))
  rc_path <- file.path(run_dir, "optim_rc_all.Rds")
  rc <- if (file.exists(rc_path)) readRDS(rc_path) else rep(0L, length(lp))

  ok <- which(is.finite(lp) & rc == 0L)
  if (!length(ok)) {
    ok <- which(is.finite(lp))
  }
  if (!length(ok)) {
    return(NULL)
  }

  best_i <- ok[which.max(lp[ok])]
  best_draw <- draws[[best_i]]

  if (is.data.frame(best_draw)) {
    best_draw <- as.matrix(best_draw)
  }

  if (!is.matrix(best_draw) || !nrow(best_draw)) {
    stop(sprintf("Unsupported best-draw format in %s", run_dir))
  }

  out <- best_draw[1, , drop = TRUE]
  names(out) <- colnames(best_draw)
  out
}

collect_best_raw_theta_ploidy_global <- function(ranking_df) {
  if (!nrow(ranking_df)) {
    return(tibble())
  }

  rows <- lapply(seq_len(nrow(ranking_df)), function(i) {
    model_id <- ranking_df$model_id[[i]]
    run_dir <- ranking_df$run_dir[[i]]
    dims <- parse_run_id(model_id)
    param_labels <- get_param_names(
      R = dims$R,
      P = dims$P,
      W = dims$W,
      C = dims$C,
      M = dims$M,
      priors = base_priors
    )

    draw_vec <- tryCatch(extract_best_draw_vector(run_dir), error = function(e) NULL)
    if (is.null(draw_vec)) {
      return(NULL)
    }

    coef_names <- sprintf("raw_theta_ploidy[%d]", seq_along(param_labels))
    values <- unname(draw_vec[coef_names])

    tibble(
      source = "global_optim",
      model_id = model_id,
      alias = ranking_df$alias[[i]],
      run_label = ranking_df$run_label[[i]],
      line_id = NA_integer_,
      direction = NA_character_,
      fit_type = "all_data_best",
      run_tag = NA_character_,
      coefficient = coef_names,
      parameter = param_labels,
      raw_theta_ploidy = values,
      missing = is.na(values),
      return_code = NA_real_,
      lp__ = unname(draw_vec["lp__"]),
      elpd = NA_real_,
      ll_transfer_train = NA_real_,
      ll_transfer_holdout = NA_real_,
      holdout_wells = NA_integer_,
      holdout_obs = NA_integer_
    )
  })

  bind_rows(rows)
}

collect_best_raw_theta_ploidy_transfer <- function(
  ranking_df,
  output_root = project_paths$transfer_optim_root,
  stan_data_path = project_paths$stan_data
) {
  model_ids <- ranking_df$model_id
  if (!length(model_ids)) {
    return(tibble())
  }

  stan_data <- readRDS(resolve_stan_data_path(stan_data_path))
  all_wells <- seq_len(stan_data$N_wells)
  rows <- list()
  idx <- 1L

  for (model_id in model_ids) {
    base_dir <- file.path(output_root, model_id)
    best_path <- file.path(base_dir, "transfer_best_start_summary.Rds")
    if (!file.exists(best_path)) {
      next
    }

    best_df <- readRDS(best_path)
    if (is.null(best_df) || !nrow(best_df)) {
      next
    }

    dims <- parse_run_id(model_id)
    param_labels <- get_param_names(
      R = dims$R,
      P = dims$P,
      W = dims$W,
      C = dims$C,
      M = dims$M,
      priors = base_priors
    )
    coef_names <- sprintf("raw_theta_ploidy[%d]", seq_along(param_labels))

    for (j in seq_len(nrow(best_df))) {
      fit_row <- best_df[j, , drop = FALSE]
      fit_obj <- tryCatch(
        load_best_transfer_fit(
          model_id = model_id,
          line_id = fit_row$line_id[[1]],
          direction = fit_row$direction[[1]],
          fit_type = fit_row$fit_type[[1]],
          output_root = dirname(output_root),
          model_name = "gpath"
        ),
        error = function(e) NULL
      )

      if (is.null(fit_obj)) {
        next
      }

      draw_vec <- extract_draw_vector(fit_obj$draws)
      values <- unname(draw_vec[coef_names])
      holdout_wells <- fit_obj$split_meta$holdout_wells
      train_wells <- setdiff(all_wells, holdout_wells)

      rows[[idx]] <- tibble(
        source = "transfer_cv_best",
        model_id = model_id,
        alias = ranking_df$alias[match(model_id, ranking_df$model_id)],
        run_label = ranking_df$run_label[match(model_id, ranking_df$model_id)],
        line_id = as.integer(fit_row$line_id[[1]]),
        direction = as.character(fit_row$direction[[1]]),
        fit_type = as.character(fit_row$fit_type[[1]]),
        run_tag = as.character(fit_row$run_tag[[1]]),
        coefficient = coef_names,
        parameter = param_labels,
        raw_theta_ploidy = values,
        missing = is.na(values),
        return_code = as.numeric(fit_row$return_code[[1]]),
        lp__ = as.numeric(fit_row$lp__[[1]]),
        elpd = as.numeric(fit_row$elpd[[1]]),
        ll_transfer_train = compute_point_loglik_subset(fit_obj$draws, train_wells),
        ll_transfer_holdout = compute_point_loglik_subset(fit_obj$draws, holdout_wells),
        holdout_wells = as.integer(length(holdout_wells)),
        holdout_obs = as.integer(fit_row$holdout_obs[[1]])
      )
      idx <- idx + 1L
    }
  }

  bind_rows(rows)
}

load_transfer_outputs <- function(model_id, output_root) {
  base_dir <- file.path(output_root, model_id)
  list(
    transfer_summary = safe_read_rds(file.path(base_dir, "transfer_comparison_summary.Rds")),
    parameter_states = safe_read_rds(file.path(base_dir, "parameter_transfer_states.Rds")),
    parameter_summary = safe_read_rds(file.path(base_dir, "parameter_transfer_summary.Rds")),
    parameter_comparison = safe_read_rds(file.path(base_dir, "parameter_transfer_comparison.Rds")),
    best_start = safe_read_rds(file.path(base_dir, "transfer_best_start_summary.Rds"))
  )
}

collect_transfer_screening <- function(model_ids, output_root = project_paths$transfer_optim_root) {
  rows <- lapply(model_ids, function(model_id) {
    x <- load_transfer_outputs(model_id, output_root = output_root)
    ts <- x$transfer_summary

    if (is.null(ts) || !nrow(ts)) {
      return(NULL)
    }

    ts$model_id <- model_id
    ts
  })

  bind_rows(rows)
}

collect_parameter_screening <- function(model_ids, output_root = project_paths$transfer_optim_root) {
  rows <- lapply(model_ids, function(model_id) {
    x <- load_transfer_outputs(model_id, output_root = output_root)
    ps <- x$parameter_summary

    if (is.null(ps) || !nrow(ps)) {
      return(NULL)
    }

    ps$model_id <- model_id
    ps
  })

  bind_rows(rows)
}

collect_parameter_family_transfer <- function(model_ids, output_root = project_paths$transfer_optim_root) {
  rows <- lapply(model_ids, function(model_id) {
    x <- load_transfer_outputs(model_id, output_root = output_root)
    pc <- x$parameter_comparison

    if (is.null(pc) || !nrow(pc)) {
      return(NULL)
    }

    pc %>%
      mutate(
        model_id = model_id,
        family = vapply(parameter, get_parameter_family, character(1))
      ) %>%
      filter(!is.na(family)) %>%
      mutate(
        abs_null_err = pmax(abs(null_vs_oracle_diff), 1e-12),
        abs_transfer_err = pmax(abs(holdout_diff), 1e-12),
        log10_error_ratio = log10(abs_transfer_err / abs_null_err),
        log_error_improvement = log(abs_null_err) - log(abs_transfer_err),
        transfer_better_than_null = abs_transfer_err < abs_null_err
      )
  })

  bind_rows(rows)
}

summarize_predictive_screening <- function(screen_df) {
  if (!nrow(screen_df)) {
    return(tibble())
  }

  screen_df %>%
    group_by(model_id, direction) %>%
    summarise(
      mean_transfer_gain = mean(transfer_gain, na.rm = TRUE),
      median_transfer_gain = median(transfer_gain, na.rm = TRUE),
      mean_transfer_regret = mean(transfer_regret, na.rm = TRUE),
      prop_transfer_beats_null = mean(transfer > null, na.rm = TRUE),
      n_cases = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = direction,
      values_from = c(
        mean_transfer_gain,
        median_transfer_gain,
        mean_transfer_regret,
        prop_transfer_beats_null,
        n_cases
      )
    ) %>%
    arrange(desc(mean_transfer_gain_low_to_high), desc(mean_transfer_gain_high_to_low))
}

summarize_parameter_screening <- function(screen_df) {
  if (!nrow(screen_df)) {
    return(tibble())
  }

  screen_df %>%
    mutate(
      mean_log_error_improvement = mean_abs_log_null_holdout_ratio - mean_abs_log_holdout_ratio
    ) %>%
    group_by(model_id, direction) %>%
    summarise(
      mean_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
      median_log_error_improvement = median(mean_log_error_improvement, na.rm = TRUE),
      mean_prop_parameters_better = mean(prop_parameters_transfer_better_than_null, na.rm = TRUE),
      n_cases = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = direction,
      values_from = c(
        mean_log_error_improvement,
        median_log_error_improvement,
        mean_prop_parameters_better,
        n_cases
      )
    ) %>%
    arrange(desc(mean_log_error_improvement_low_to_high), desc(mean_log_error_improvement_high_to_low))
}

summarize_parameter_family_transfer <- function(family_df) {
  if (!nrow(family_df)) {
    return(tibble())
  }

  family_df %>%
    group_by(model_id, line_id, direction, family) %>%
    summarise(
      median_log10_error_ratio = median(log10_error_ratio, na.rm = TRUE),
      mean_log_error_improvement = mean(log_error_improvement, na.rm = TRUE),
      prop_transfer_better_than_null = mean(transfer_better_than_null, na.rm = TRUE),
      n_parameters = dplyr::n(),
      .groups = "drop"
    )
}

collect_focused_parameter_transfer <- function(
  model_ids,
  focus_parameters = focus_common_parameters,
  output_root = project_paths$transfer_optim_root
) {
  rows <- lapply(model_ids, function(model_id) {
    x <- load_transfer_outputs(model_id, output_root = output_root)
    pc <- x$parameter_comparison

    if (is.null(pc) || !nrow(pc)) {
      return(NULL)
    }

    pc %>%
      mutate(parameter_focus = normalize_focus_parameter(parameter)) %>%
      filter(parameter_focus %in% focus_parameters) %>%
      mutate(
        model_id = model_id,
        parameter = parameter_focus,
        abs_null_err = abs(null_vs_oracle_diff),
        abs_transfer_err = abs(holdout_diff),
        abs_error_improvement = abs_null_err - abs_transfer_err,
        positive_triplet = null_holdout > 0 & transfer_holdout > 0 & oracle_holdout > 0,
        log_null_ratio = ifelse(
          positive_triplet,
          log(null_holdout / oracle_holdout),
          NA_real_
        ),
        log_transfer_ratio = ifelse(
          positive_triplet,
          log(transfer_holdout / oracle_holdout),
          NA_real_
        ),
        abs_log_null_ratio = abs(log_null_ratio),
        abs_log_transfer_ratio = abs(log_transfer_ratio),
        ratio_error_improvement = abs(log_null_ratio) - abs(log_transfer_ratio),
        transfer_beats_null_abs = abs_transfer_err < abs_null_err,
        transfer_beats_null_ratio = abs(log_transfer_ratio) < abs(log_null_ratio),
        line_direction = sprintf("line %d | %s", line_id, direction)
      )
  })

  bind_rows(rows)
}

summarize_parameter_family_across_models <- function(family_summary_df) {
  if (!nrow(family_summary_df)) {
    return(tibble())
  }

  family_summary_df %>%
    group_by(model_id, direction, family) %>%
    summarise(
      mean_family_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
      median_family_log10_error_ratio = median(median_log10_error_ratio, na.rm = TRUE),
      mean_family_win_rate = mean(prop_transfer_better_than_null, na.rm = TRUE),
      n_line_cases = dplyr::n(),
      .groups = "drop"
    )
}

rank_parameter_families <- function(model_family_df) {
  if (!nrow(model_family_df)) {
    return(tibble())
  }

  model_family_df %>%
    group_by(direction, family) %>%
    summarise(
      mean_log_error_improvement = mean(mean_family_log_error_improvement, na.rm = TRUE),
      median_log10_error_ratio = median(median_family_log10_error_ratio, na.rm = TRUE),
      mean_win_rate = mean(mean_family_win_rate, na.rm = TRUE),
      prop_models_positive = mean(mean_family_log_error_improvement > 0, na.rm = TRUE),
      n_models = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(direction, desc(mean_log_error_improvement), desc(mean_win_rate))
}

derive_glucose_transfer_summaries <- function(state_df) {
  if (is.null(state_df) || !nrow(state_df)) {
    return(tibble(
      line_id = integer(),
      direction = character(),
      fit_type = character(),
      state = character(),
      summary_name = character(),
      summary_value = numeric()
    ))
  }

  glucose_name_map <- c(
    "ae[1]" = "ae1",
    "ae" = "ae1",
    "ah[1]" = "ah1",
    "ah" = "ah1",
    "Y_R[1]" = "YR1",
    "Y_R" = "YR1",
    "m" = "m"
  )

  wide <- state_df %>%
    filter(parameter %in% names(glucose_name_map), fit_type %in% c("transfer", "oracle")) %>%
    mutate(parameter_key = unname(glucose_name_map[parameter])) %>%
    select(line_id, direction, fit_type, state, parameter_key, value) %>%
    tidyr::pivot_wider(names_from = parameter_key, values_from = value)

  needed <- c("ae1", "ah1", "YR1", "m")
  if (!all(needed %in% names(wide))) {
    return(tibble(
      line_id = integer(),
      direction = character(),
      fit_type = character(),
      state = character(),
      summary_name = character(),
      summary_value = numeric()
    ))
  }

  wide %>%
    transmute(
      line_id,
      direction,
      fit_type,
      state,
      Umax_1 = ae1,
      Khalf_1 = ah1,
      Efflow_1 = ae1 / ah1,
      YieldCap_1 = YR1 * ae1,
      NetGrow_1 = (YR1 * ae1) - m
    ) %>%
    tidyr::pivot_longer(
      cols = c(Umax_1, Khalf_1, Efflow_1, YieldCap_1, NetGrow_1),
      names_to = "summary_name",
      values_to = "summary_value"
    )
}

collect_derived_transfer_summaries <- function(model_ids, output_root = project_paths$transfer_optim_root) {
  rows <- lapply(model_ids, function(model_id) {
    x <- load_transfer_outputs(model_id, output_root = output_root)
    ds <- derive_glucose_transfer_summaries(x$parameter_states)

    if (is.null(ds) || !nrow(ds)) {
      return(NULL)
    }

    ds$model_id <- model_id
    ds
  })

  out <- bind_rows(rows)
  if (!nrow(out)) {
    return(tibble(
      line_id = integer(),
      direction = character(),
      fit_type = character(),
      state = character(),
      summary_name = character(),
      summary_value = numeric(),
      model_id = character()
    ))
  }
  out
}

build_derived_transfer_comparison <- function(derived_df) {
  if (!nrow(derived_df)) {
    return(tibble(
      model_id = character(),
      line_id = integer(),
      direction = character(),
      summary_name = character(),
      observed_transfer = numeric(),
      transfer_holdout = numeric(),
      oracle_holdout = numeric(),
      observed_oracle = numeric(),
      null_holdout = numeric(),
      null_vs_oracle_diff = numeric(),
      holdout_diff = numeric(),
      abs_null_err = numeric(),
      abs_transfer_err = numeric(),
      log10_error_ratio = numeric(),
      log_error_improvement = numeric(),
      transfer_better_than_null = logical(),
      oracle_shift = numeric(),
      transfer_shift = numeric(),
      shift_diff = numeric()
    ))
  }

  transfer_obs <- derived_df %>%
    filter(fit_type == "transfer", state == "observed") %>%
    select(model_id, line_id, direction, summary_name, observed_transfer = summary_value)

  transfer_holdout <- derived_df %>%
    filter(fit_type == "transfer", state == "holdout") %>%
    select(model_id, line_id, direction, summary_name, transfer_holdout = summary_value)

  oracle_holdout <- derived_df %>%
    filter(fit_type == "oracle", state == "holdout") %>%
    select(model_id, line_id, direction, summary_name, oracle_holdout = summary_value)

  oracle_obs <- derived_df %>%
    filter(fit_type == "oracle", state == "observed") %>%
    select(model_id, line_id, direction, summary_name, observed_oracle = summary_value)

  transfer_obs %>%
    left_join(transfer_holdout, by = c("model_id", "line_id", "direction", "summary_name")) %>%
    left_join(oracle_holdout, by = c("model_id", "line_id", "direction", "summary_name")) %>%
    left_join(oracle_obs, by = c("model_id", "line_id", "direction", "summary_name")) %>%
    mutate(
      null_holdout = observed_transfer,
      null_vs_oracle_diff = null_holdout - oracle_holdout,
      holdout_diff = transfer_holdout - oracle_holdout,
      abs_null_err = pmax(abs(null_vs_oracle_diff), 1e-12),
      abs_transfer_err = pmax(abs(holdout_diff), 1e-12),
      log10_error_ratio = log10(abs_transfer_err / abs_null_err),
      log_error_improvement = log(abs_null_err) - log(abs_transfer_err),
      transfer_better_than_null = abs_transfer_err < abs_null_err,
      oracle_shift = oracle_holdout - observed_oracle,
      transfer_shift = transfer_holdout - observed_transfer,
      shift_diff = transfer_shift - oracle_shift
    )
}

summarize_derived_transfer <- function(derived_cmp_df) {
  if (!nrow(derived_cmp_df)) {
    return(tibble(
      model_id = character(),
      direction = character(),
      summary_name = character(),
      mean_log_error_improvement = numeric(),
      median_log10_error_ratio = numeric(),
      mean_win_rate = numeric(),
      prop_lines_positive = numeric(),
      mean_abs_shift_diff = numeric(),
      n_line_cases = integer()
    ))
  }

  derived_cmp_df %>%
    group_by(model_id, direction, summary_name) %>%
    summarise(
      mean_log_error_improvement = mean(log_error_improvement, na.rm = TRUE),
      median_log10_error_ratio = median(log10_error_ratio, na.rm = TRUE),
      mean_win_rate = mean(transfer_better_than_null, na.rm = TRUE),
      prop_lines_positive = mean(log_error_improvement > 0, na.rm = TRUE),
      mean_abs_shift_diff = mean(abs(shift_diff), na.rm = TRUE),
      n_line_cases = dplyr::n(),
      .groups = "drop"
    )
}

rank_derived_transfer <- function(derived_summary_df) {
  if (!nrow(derived_summary_df)) {
    return(tibble(
      direction = character(),
      summary_name = character(),
      mean_log_error_improvement = numeric(),
      mean_win_rate = numeric(),
      prop_models_positive = numeric(),
      mean_abs_shift_diff = numeric(),
      n_models = integer()
    ))
  }

  derived_summary_df %>%
    group_by(direction, summary_name) %>%
    summarise(
      mean_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
      mean_win_rate = mean(mean_win_rate, na.rm = TRUE),
      prop_models_positive = mean(mean_log_error_improvement > 0, na.rm = TRUE),
      mean_abs_shift_diff = mean(mean_abs_shift_diff, na.rm = TRUE),
      n_models = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(direction, desc(mean_log_error_improvement), desc(mean_win_rate))
}

compute_point_loglik_subset <- function(draws, well_idx) {
  ll <- get_well_loglik_draws(draws, well_idx)
  if (!length(ll)) {
    return(NA_real_)
  }
  as.numeric(ll[1])
}

build_transfer_subset_diagnostics <- function(
  model_ids,
  output_root = project_paths$transfer_optim_root,
  stan_data_path = project_paths$stan_data
) {
  if (!length(model_ids)) {
    return(tibble(
      model_id = character(),
      line_id = integer(),
      direction = character(),
      fit_type = character(),
      return_code = numeric(),
      lp__ = numeric(),
      ll_transfer_train = numeric(),
      ll_transfer_holdout = numeric(),
      n_transfer_train_wells = integer(),
      n_transfer_holdout_wells = integer()
    ))
  }

  stan_data <- readRDS(resolve_stan_data_path(stan_data_path))
  all_wells <- seq_len(stan_data$N_wells)
  rows <- list()
  idx <- 1L

  for (model_id in model_ids) {
    best_path <- file.path(output_root, model_id, "transfer_best_start_summary.Rds")
    if (!file.exists(best_path)) {
      next
    }

    best_df <- readRDS(best_path)
    keys <- unique(best_df[, c("line_id", "direction")])

    for (k in seq_len(nrow(keys))) {
      line_id <- keys$line_id[k]
      direction <- keys$direction[k]
      split_meta <- get_directional_transfer_split(stan_data, line_id, direction)
      transfer_train_wells <- setdiff(all_wells, split_meta$holdout_wells)
      transfer_holdout_wells <- split_meta$holdout_wells

      for (fit_type in c("null", "transfer", "oracle")) {
        fit_obj <- tryCatch(
          load_transfer_best_fit(
            model_id = model_id,
            line_id = line_id,
            direction = direction,
            fit_type = fit_type,
            output_root = dirname(output_root),
            model_name = "gpath"
          ),
          error = function(e) NULL
        )

        if (is.null(fit_obj)) {
          next
        }

        rows[[idx]] <- data.frame(
          model_id = model_id,
          line_id = as.integer(line_id),
          direction = direction,
          fit_type = fit_type,
          return_code = as.numeric(fit_obj$summary$return_code[[1]]),
          lp__ = as.numeric(fit_obj$summary$lp__[[1]]),
          ll_transfer_train = compute_point_loglik_subset(fit_obj$draws, transfer_train_wells),
          ll_transfer_holdout = compute_point_loglik_subset(fit_obj$draws, transfer_holdout_wells),
          n_transfer_train_wells = length(transfer_train_wells),
          n_transfer_holdout_wells = length(transfer_holdout_wells),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (!length(rows)) {
    return(tibble(
      model_id = character(),
      line_id = integer(),
      direction = character(),
      fit_type = character(),
      return_code = numeric(),
      lp__ = numeric(),
      ll_transfer_train = numeric(),
      ll_transfer_holdout = numeric(),
      n_transfer_train_wells = integer(),
      n_transfer_holdout_wells = integer()
    ))
  }

  bind_rows(rows)
}

summarize_transfer_subset_diagnostics <- function(diag_df) {
  if (!nrow(diag_df)) {
    return(tibble(
      model_id = character(),
      line_id = integer(),
      direction = character(),
      transfer_minus_oracle_train = numeric(),
      transfer_minus_null_train = numeric(),
      transfer_minus_oracle_holdout = numeric(),
      transfer_minus_null_holdout = numeric()
    ))
  }

  wide <- diag_df %>%
    select(model_id, line_id, direction, fit_type, ll_transfer_train, ll_transfer_holdout) %>%
    tidyr::pivot_wider(
      names_from = fit_type,
      values_from = c(ll_transfer_train, ll_transfer_holdout)
    )

  wide %>%
    transmute(
      model_id,
      line_id,
      direction,
      transfer_minus_oracle_train = ll_transfer_train_transfer - ll_transfer_train_oracle,
      transfer_minus_null_train = ll_transfer_train_transfer - ll_transfer_train_null,
      transfer_minus_oracle_holdout = ll_transfer_holdout_transfer - ll_transfer_holdout_oracle,
      transfer_minus_null_holdout = ll_transfer_holdout_transfer - ll_transfer_holdout_null
    )
}
