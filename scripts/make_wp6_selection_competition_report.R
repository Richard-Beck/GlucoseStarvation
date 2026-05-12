#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

OUTPUT_DIR <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp6_selection_competition")
}
FIGURE_DIR <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp6_selection_competition")
}
WP5_DECISION_PATH <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("data", "report_exports", "wp5_mechanistic_representative_model", "wp5_decision_table.csv")
}
ALL_LINE_DATASET <- if (length(args) >= 4 && nzchar(args[4])) {
  args[4]
} else {
  "gstarvation_v1"
}
NO_SUM_DATASET <- if (length(args) >= 5 && nzchar(args[5])) {
  args[5]
} else {
  "gstarvation_v1_no_sum159_fuse"
}
COMPETITION_SCORE_PATH <- if (length(args) >= 6 && nzchar(args[6])) {
  args[6]
} else {
  file.path("data", "selection_strategy", "gpath", "v1", "gstarvation_v1_single_line_competition_validation", "model_scores.csv")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("R/selection_strategy_utils.R")

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

wp6_theme <- function(base_size = 8) {
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

write_table_pair <- function(x, stem) {
  csv_path <- file.path(OUTPUT_DIR, paste0(stem, ".csv"))
  rds_path <- file.path(OUTPUT_DIR, paste0(stem, ".Rds"))
  utils::write.csv(as.data.frame(x), csv_path, row.names = FALSE)
  saveRDS(x, rds_path)
  invisible(c(csv = csv_path, rds = rds_path))
}

save_plot_pair <- function(plot, stem, width, height, dpi = 450) {
  pdf_path <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURE_DIR, paste0(stem, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf = pdf_path, png = png_path))
}

mode_string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

load_wp5_model_info <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("WP5 decision table not found: %s", path))
  }

  read.csv(path, stringsAsFactors = FALSE) %>%
    mutate(
      wp6_role = case_when(
        decision_status == "provisional_representative_visual_review_passed" ~ "primary",
        decision_status == "candidate_accepted_visual_review" ~ "robustness",
        decision_status == "candidate_visual_review" & nuts_qc_status == "pass" ~ "sensitivity_pending_visual_review",
        decision_status == "legacy_default_rejected_for_wp5_review" ~ "legacy_rejected",
        TRUE ~ "excluded"
      ),
      wp6_main = wp6_role %in% c("primary", "robustness"),
      alias = if_else(is.na(alias) | !nzchar(alias), model_id, alias)
    )
}

load_selection_scope <- function(scope, dataset_label, stan_data_path = NULL) {
  cfg <- selection_default_config()
  cfg$dataset_label <- dataset_label
  cfg$stan_data_path <- stan_data_path

  out <- tryCatch(
    collect_selection_outputs(cfg),
    error = function(e) {
      warning(sprintf("Could not collect selection outputs for %s (%s): %s", scope, dataset_label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(out) || !nrow(out$final_summary)) {
    return(tibble())
  }

  as_tibble(out$final_summary) %>%
    mutate(
      scope = scope,
      dataset_label = dataset_label,
      context_label = if_else(
        fit_context == "global",
        "global",
        paste(direction, fit_type, sep = "_")
      ),
      scope_context = paste(scope, context_label, sep = " / "),
      direction = if_else(is.na(direction), NA_character_, as.character(direction)),
      fit_type = if_else(is.na(fit_type), NA_character_, as.character(fit_type))
    )
}

rank_selection_contexts <- function(final_df) {
  final_df %>%
    mutate(
      for_rank_score = if_else(is.finite(for_high_viable), for_high_viable, select_for_high_score),
      against_rank_score = if_else(is.finite(against_high_viable), against_high_viable, select_against_high_score)
    ) %>%
    group_by(scope, dataset_label, model_id, line_id, context_label) %>%
    mutate(
      rank_for_high = min_rank(desc(for_rank_score)),
      rank_against_high = min_rank(desc(against_rank_score))
    ) %>%
    ungroup()
}

build_strategy_summary <- function(ranked_df) {
  bind_rows(
    ranked_df %>%
      group_by(scope, scope_context, strategy_code) %>%
      summarise(
        objective = "for_high",
        mean_rank = mean(rank_for_high, na.rm = TRUE),
        top5_fraction = mean(rank_for_high <= 5, na.rm = TRUE),
        favorable_fraction = mean(log_ratio_final > 0, na.rm = TRUE),
        mean_log_ratio = mean(log_ratio_final, na.rm = TRUE),
        q10_log_ratio = quantile(log_ratio_final, 0.10, na.rm = TRUE),
        q90_log_ratio = quantile(log_ratio_final, 0.90, na.rm = TRUE),
        n_contexts = n(),
        .groups = "drop"
      ),
    ranked_df %>%
      group_by(scope, scope_context, strategy_code) %>%
      summarise(
        objective = "against_high",
        mean_rank = mean(rank_against_high, na.rm = TRUE),
        top5_fraction = mean(rank_against_high <= 5, na.rm = TRUE),
        favorable_fraction = mean(log_ratio_final < 0, na.rm = TRUE),
        mean_log_ratio = mean(log_ratio_final, na.rm = TRUE),
        q10_log_ratio = quantile(log_ratio_final, 0.10, na.rm = TRUE),
        q90_log_ratio = quantile(log_ratio_final, 0.90, na.rm = TRUE),
        n_contexts = n(),
        .groups = "drop"
      )
  ) %>%
    arrange(scope_context, objective, mean_rank)
}

build_best_context_summary <- function(ranked_df) {
  key_cols <- c("scope", "scope_context", "dataset_label", "model_id", "line_id", "context_label")

  for_best <- ranked_df %>%
    group_by(across(all_of(key_cols))) %>%
    slice_min(order_by = rank_for_high, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(objective = "for_high", objective_rank = rank_for_high)

  against_best <- ranked_df %>%
    group_by(across(all_of(key_cols))) %>%
    slice_min(order_by = rank_against_high, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(objective = "against_high", objective_rank = rank_against_high)

  bind_rows(for_best, against_best) %>%
    group_by(scope, scope_context, dataset_label, model_id, objective) %>%
    summarise(
      n_contexts = n(),
      modal_strategy_code = mode_string(strategy_code),
      mean_best_log_ratio = mean(log_ratio_final, na.rm = TRUE),
      median_best_log_ratio = median(log_ratio_final, na.rm = TRUE),
      q10_best_log_ratio = quantile(log_ratio_final, 0.10, na.rm = TRUE),
      q90_best_log_ratio = quantile(log_ratio_final, 0.90, na.rm = TRUE),
      success_fraction = if_else(
        objective[[1]] == "for_high",
        mean(log_ratio_final > 0, na.rm = TRUE),
        mean(log_ratio_final < 0, na.rm = TRUE)
      ),
      mean_high_fraction_final = mean(high_fraction_final, na.rm = TRUE),
      .groups = "drop"
    )
}

build_balance_summary <- function(best_summary) {
  best_summary %>%
    select(scope, scope_context, dataset_label, model_id, objective, mean_best_log_ratio, success_fraction, modal_strategy_code) %>%
    pivot_wider(
      names_from = objective,
      values_from = c(mean_best_log_ratio, success_fraction, modal_strategy_code),
      names_sep = "__"
    ) %>%
    mutate(
      for_abs_log_ratio = abs(mean_best_log_ratio__for_high),
      against_abs_log_ratio = abs(mean_best_log_ratio__against_high),
      for_minus_against_abs = for_abs_log_ratio - against_abs_log_ratio,
      for_success_minus_against = success_fraction__for_high - success_fraction__against_high,
      selection_balance = case_when(
        for_minus_against_abs > 0 & for_success_minus_against >= 0 ~ "for_high_stronger",
        for_minus_against_abs < 0 & for_success_minus_against <= 0 ~ "against_high_stronger",
        TRUE ~ "mixed"
      )
    )
}

write_decision_note <- function(path, model_info, balance_summary, competition_scores, no_sum_available) {
  main_models <- model_info %>%
    filter(wp6_main) %>%
    arrange(match(wp6_role, c("primary", "robustness")), model_id)

  primary_model <- main_models$model_id[main_models$wp6_role == "primary"][[1]] %||% NA_character_
  scope_lines <- balance_summary %>%
    select(-any_of(c("alias", "wp6_role"))) %>%
    left_join(model_info %>% select(model_id, alias, wp6_role), by = "model_id") %>%
    filter(wp6_role %in% c("primary", "robustness")) %>%
    group_by(scope_context) %>%
    summarise(
      n_models = n_distinct(model_id),
      mean_for_abs_log_ratio = mean(for_abs_log_ratio, na.rm = TRUE),
      mean_against_abs_log_ratio = mean(against_abs_log_ratio, na.rm = TRUE),
      mean_for_success = mean(success_fraction__for_high, na.rm = TRUE),
      mean_against_success = mean(success_fraction__against_high, na.rm = TRUE),
      n_for_high_stronger = sum(selection_balance == "for_high_stronger", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      line = sprintf(
        "- `%s`: %d models; mean |for| %.3f vs |against| %.3f; success %.2f vs %.2f; for-high stronger in %d models.",
        scope_context, n_models, mean_for_abs_log_ratio, mean_against_abs_log_ratio,
        mean_for_success, mean_against_success, n_for_high_stronger
      )
    ) %>%
    pull(line)

  comp_lines <- character(0)
  if (nrow(competition_scores)) {
    comp_lines <- competition_scores %>%
      select(-any_of(c("alias", "wp6_role", "wp6_main"))) %>%
      left_join(model_info %>% select(model_id, alias, wp6_role), by = "model_id") %>%
      filter(wp6_role %in% c("primary", "robustness")) %>%
      arrange(validation_score) %>%
      transmute(
        line = sprintf(
          "- `%s` (%s): sign match %.2f, high-fraction MAE %.3f, validation score %.3f.",
          model_id, alias, glucose_sign_match, mae_high_fraction, validation_score
        )
      ) %>%
      pull(line)
  }

  no_sum_sentence <- if (isTRUE(no_sum_available)) {
    "No-SUM-159-fuse global selection simulations were available and included."
  } else {
    "No-SUM-159-fuse global selection simulations were not available; run `scripts/run_gpath_selection_strategy_batch.R` with `scripts/wp6_no_sum_selection_config.R`."
  }

  lines <- c(
    "# WP6 selection and competition decision note",
    "",
    sprintf("Generated: %s", as.character(Sys.time())),
    "",
    "## Model set",
    "",
    sprintf("Primary WP6 model: `%s`.", primary_model),
    "",
    main_models %>%
      transmute(line = sprintf("- `%s` (`%s`): %s.", model_id, alias, wp6_role)) %>%
      pull(line),
    "",
    "`1R_1P_0W_C0_M1` remains excluded from main WP6 claims because WP5 rejected its visual fit.",
    "",
    "## Selection result",
    "",
    scope_lines,
    "",
    no_sum_sentence,
    "",
    "All-line global and transfer simulations are reused from `data/selection_strategy/gpath/v1/gstarvation_v1/`. The no-SUM check currently uses global no-SUM fits from the maintained optimisation workflow; no-SUM transfer simulations are intentionally not mixed in because the selection transfer loader still targets the older `data/gpath_transfer_cv` layout.",
    "",
    "## Competition check",
    "",
    if (length(comp_lines)) comp_lines else "No single-line competition validation score table was available.",
    "",
    "Competition validation is treated as a check on direction and scale, not as the primary model-selection criterion."
  )

  writeLines(lines, path)
  invisible(path)
}

model_info <- load_wp5_model_info(WP5_DECISION_PATH)
write_table_pair(model_info, "wp6_model_set_from_wp5")

main_model_ids <- model_info %>%
  filter(wp6_main) %>%
  pull(model_id) %>%
  unique()

if (!length(main_model_ids)) {
  stop("No WP6 main models were selected from the WP5 decision table")
}

all_line_final <- load_selection_scope("all_lines", ALL_LINE_DATASET)
no_sum_final <- load_selection_scope(
  "no_sum159_fuse",
  NO_SUM_DATASET,
  stan_data_path = file.path("data", "inputs", "stan", "wp4b_no_sum159_fuse", "stan_ready_data.Rds")
)
no_sum_available <- nrow(no_sum_final) > 0

final_df <- bind_rows(all_line_final, no_sum_final) %>%
  filter(model_id %in% main_model_ids) %>%
  left_join(model_info %>% select(model_id, alias, wp6_role), by = "model_id")

if (!nrow(final_df)) {
  stop("No selection-strategy rows matched the WP6 main model set")
}

ranked_df <- rank_selection_contexts(final_df)
write_table_pair(ranked_df, "wp6_selection_final_ranked")

strategy_summary <- build_strategy_summary(ranked_df)
write_table_pair(strategy_summary, "wp6_strategy_summary_by_scope")

best_context_summary <- build_best_context_summary(ranked_df) %>%
  left_join(model_info %>% select(model_id, alias, wp6_role), by = "model_id")
write_table_pair(best_context_summary, "wp6_best_context_summary_by_model")

balance_summary <- build_balance_summary(best_context_summary) %>%
  left_join(model_info %>% select(model_id, alias, wp6_role), by = "model_id")
write_table_pair(balance_summary, "wp6_selection_balance_by_model")

competition_scores <- if (file.exists(COMPETITION_SCORE_PATH)) {
  read.csv(COMPETITION_SCORE_PATH, stringsAsFactors = FALSE) %>%
    left_join(model_info %>% select(model_id, alias, wp6_role, wp6_main), by = "model_id")
} else {
  tibble()
}
write_table_pair(competition_scores, "wp6_single_line_competition_scores")

top_strategy_plot_df <- strategy_summary %>%
  group_by(scope_context, objective) %>%
  slice_min(order_by = mean_rank, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    objective = factor(objective, levels = c("against_high", "for_high")),
    strategy_label = factor(strategy_code, levels = rev(unique(strategy_code[order(scope_context, objective, mean_rank)])))
  )

strategy_plot <- ggplot(top_strategy_plot_df, aes(x = reorder(strategy_code, mean_rank), y = mean_log_ratio, fill = favorable_fraction)) +
  geom_col(width = 0.75, color = "grey25", linewidth = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey35") +
  coord_flip() +
  facet_grid(scope_context ~ objective, scales = "free_y") +
  scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title = "WP6 strategy landscape after WP5 model filtering",
    x = "strategy code",
    y = "mean final log( high / low )",
    fill = "directional\nfraction"
  ) +
  wp6_theme()
save_plot_pair(strategy_plot, "wp6_strategy_landscape", width = 11, height = 7.5)

balance_plot_df <- best_context_summary %>%
  filter(wp6_role %in% c("primary", "robustness")) %>%
  mutate(
    model_label = paste0(alias, "\n", model_id),
    objective = factor(objective, levels = c("against_high", "for_high"))
  )

balance_plot <- ggplot(balance_plot_df, aes(x = model_label, y = mean_best_log_ratio, fill = objective)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, color = "grey25", linewidth = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey35") +
  facet_wrap(~ scope_context, scales = "free_y") +
  scale_fill_manual(values = c(against_high = "#4777B3", for_high = "#D55E00")) +
  labs(
    title = "Best attainable selection direction by accepted model",
    x = NULL,
    y = "mean best final log( high / low )",
    fill = "objective"
  ) +
  wp6_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
save_plot_pair(balance_plot, "wp6_model_selection_balance", width = 11, height = 6.5)

if (nrow(competition_scores)) {
  competition_plot_df <- competition_scores %>%
    filter(wp6_role %in% c("primary", "robustness")) %>%
    mutate(model_label = paste0(alias, "\n", model_id))

  competition_plot <- ggplot(competition_plot_df, aes(x = reorder(model_label, validation_score), y = validation_score, fill = glucose_sign_match)) +
    geom_col(width = 0.7, color = "grey25", linewidth = 0.15) +
    coord_flip() +
    scale_fill_viridis_c(option = "D", limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(
      title = "Single-line competition validation check",
      x = NULL,
      y = "validation score (lower is better)",
      fill = "glucose sign\nmatch"
    ) +
    wp6_theme()
  save_plot_pair(competition_plot, "wp6_single_line_competition_validation", width = 8.2, height = 4.8)
}

decision_note_path <- file.path(OUTPUT_DIR, "wp6_selection_competition_decision_note.md")
write_decision_note(
  path = decision_note_path,
  model_info = model_info,
  balance_summary = balance_summary,
  competition_scores = competition_scores,
  no_sum_available = no_sum_available
)

summary_lines <- c(
  "WP6_SELECTION_COMPETITION",
  sprintf("generated\t%s", as.character(Sys.time())),
  sprintf("wp5_decision_path\t%s", WP5_DECISION_PATH),
  sprintf("all_line_dataset\t%s", ALL_LINE_DATASET),
  sprintf("no_sum_dataset\t%s", NO_SUM_DATASET),
  sprintf("n_main_models\t%d", length(main_model_ids)),
  sprintf("n_selection_rows\t%d", nrow(final_df)),
  sprintf("no_sum_available\t%s", no_sum_available),
  sprintf("decision_note\t%s", decision_note_path)
)
writeLines(summary_lines, file.path(OUTPUT_DIR, "wp6_selection_competition_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
