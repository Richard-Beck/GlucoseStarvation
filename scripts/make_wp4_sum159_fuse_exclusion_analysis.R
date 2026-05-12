args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp4_sum159_fuse_exclusion")
}

figure_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp4_sum159_fuse_exclusion")
}

stan_data_path_arg <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("R/project_paths.R")
source("R/model_free_ploidy_utils.R")

sum_fuse_line <- "SUM-159-fuse"
sign_tol <- 1e-8

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

wp4_theme <- function(base_size = 8) {
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
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))

  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)

  invisible(c(pdf = pdf_path, png = png_path))
}

wrap_label <- function(x, width = 22) {
  vapply(
    x,
    function(one) paste(strwrap(one, width = width), collapse = "\n"),
    character(1)
  )
}

feature_catalog_for_wp4 <- function() {
  get_model_free_feature_catalog() %>%
    mutate(
      display_label = case_when(
        feature == "growth_lowG_median" ~ "Low-G growth",
        feature == "growth_highG_median" ~ "High-G exp. growth",
        feature == "death_lowG_median" ~ "Low-G death",
        feature == "death_highG_median" ~ "High-G death",
        feature == "yield_alive_auc_intercept" ~ "Alive AUC baseline",
        feature == "yield_alive_auc_slope" ~ "Alive AUC glucose response",
        feature == "peak_total_yield_intercept" ~ "Yield baseline",
        feature == "peak_total_yield_slope" ~ "Yield glucose response",
        TRUE ~ short_label
      ),
      display_label = factor(display_label, levels = display_label),
      category = factor(category, levels = c("Growth", "Death", "Yield"))
    )
}

build_effect_scope_tables <- function(signature_panel) {
  catalog <- feature_catalog_for_wp4()
  signature_no_sum <- signature_panel %>% filter(cellLine != sum_fuse_line)

  effects_all <- compute_empirical_effects(signature_panel)
  effects_no_sum <- compute_empirical_effects(signature_no_sum)

  paired_effects <- bind_rows(
    effects_all$paired_effects_long %>% mutate(scope = "All lines"),
    effects_no_sum$paired_effects_long %>% mutate(scope = "SUM-159-fuse excluded")
  ) %>%
    left_join(catalog, by = "feature") %>%
    mutate(
      scope = factor(scope, levels = c("All lines", "SUM-159-fuse excluded")),
      effect_direction = case_when(
        effect_per_ploidy > sign_tol ~ "higher in high ploidy",
        effect_per_ploidy < -sign_tol ~ "lower in high ploidy",
        TRUE ~ "near zero"
      )
    )

  effect_scale <- paired_effects %>%
    filter(scope == "All lines") %>%
    group_by(feature) %>%
    summarise(
      center = median(effect_per_ploidy, na.rm = TRUE),
      spread = safe_mad(effect_per_ploidy),
      .groups = "drop"
    )

  paired_effects <- paired_effects %>%
    left_join(effect_scale, by = "feature") %>%
    mutate(
      effect_z = if_else(
        is.finite(effect_per_ploidy) & is.finite(center) & is.finite(spread) & spread > 0,
        (effect_per_ploidy - center) / spread,
        NA_real_
      )
    )

  consistency <- paired_effects %>%
    group_by(scope, feature, display_label, category) %>%
    summarise(
      n_lines = sum(is.finite(effect_per_ploidy)),
      n_positive = sum(effect_per_ploidy > sign_tol, na.rm = TRUE),
      n_negative = sum(effect_per_ploidy < -sign_tol, na.rm = TRUE),
      n_near_zero = sum(abs(effect_per_ploidy) <= sign_tol | !is.finite(effect_per_ploidy), na.rm = TRUE),
      consensus_n = max(n_positive, n_negative),
      consensus_fraction = if_else(n_lines > 0, consensus_n / n_lines, NA_real_),
      consensus_direction = case_when(
        n_positive > n_negative ~ "positive",
        n_negative > n_positive ~ "negative",
        TRUE ~ "tie"
      ),
      mean_effect_per_ploidy = mean(effect_per_ploidy, na.rm = TRUE),
      median_effect_per_ploidy = median(effect_per_ploidy, na.rm = TRUE),
      median_abs_effect_per_ploidy = median(abs(effect_per_ploidy), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scope, category, display_label)

  consistency_comparison <- consistency %>%
    mutate(
      scope_key = recode(
        as.character(scope),
        "All lines" = "all_lines",
        "SUM-159-fuse excluded" = "no_sum159_fuse"
      )
    ) %>%
    select(
      scope_key, feature, display_label, category,
      n_lines, consensus_n, consensus_fraction, consensus_direction,
      mean_effect_per_ploidy, median_abs_effect_per_ploidy
    ) %>%
    pivot_wider(
      names_from = scope_key,
      values_from = c(
        n_lines, consensus_n, consensus_fraction, consensus_direction,
        mean_effect_per_ploidy, median_abs_effect_per_ploidy
      ),
      names_glue = "{.value}__{scope_key}"
    ) %>%
    mutate(
      delta_consensus_fraction = consensus_fraction__no_sum159_fuse - consensus_fraction__all_lines,
      delta_consensus_n = consensus_n__no_sum159_fuse - consensus_n__all_lines,
      consensus_direction_changed = consensus_direction__no_sum159_fuse != consensus_direction__all_lines,
      exclusion_consistency_call = case_when(
        delta_consensus_fraction > 0 ~ "improved",
        delta_consensus_fraction < 0 ~ "worse",
        TRUE ~ "unchanged"
      )
    ) %>%
    arrange(category, display_label)

  effect_matrices <- bind_rows(
    effects_all$effect_matrix %>% mutate(scope = "All lines", .before = 1),
    effects_no_sum$effect_matrix %>% mutate(scope = "SUM-159-fuse excluded", .before = 1)
  )

  list(
    paired_effects = paired_effects,
    consistency = consistency,
    consistency_comparison = consistency_comparison,
    effect_matrices = effect_matrices,
    pca_scores_all = effects_all$pca_scores,
    pca_scores_no_sum = effects_no_sum$pca_scores
  )
}

build_transfer_scope_tables <- function(signature_panel) {
  signature_no_sum <- signature_panel %>% filter(cellLine != sum_fuse_line)

  transfer_all <- evaluate_transfer_predictions(signature_panel)
  transfer_no_sum <- evaluate_transfer_predictions(signature_no_sum)

  predictions <- bind_rows(
    transfer_all$predictions %>% mutate(scope = "All lines"),
    transfer_all$predictions %>%
      filter(cellLine != sum_fuse_line) %>%
      mutate(scope = "All-line estimator on non-SUM holdouts"),
    transfer_no_sum$predictions %>% mutate(scope = "SUM-159-fuse excluded")
  ) %>%
    mutate(
      scope = factor(
        scope,
        levels = c("All lines", "All-line estimator on non-SUM holdouts", "SUM-159-fuse excluded")
      )
    )

  overall <- predictions %>%
    group_by(scope) %>%
    summarise(
      n_predictions = n(),
      n_cell_lines = n_distinct(cellLine),
      mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
      median_scaled_err_improvement = median(scaled_abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      mean_abs_err_improvement = mean(abs_err_improvement, na.rm = TRUE),
      median_abs_err_improvement = median(abs_err_improvement, na.rm = TRUE),
      .groups = "drop"
    )

  feature_summary <- predictions %>%
    group_by(scope, feature) %>%
    summarise(
      n_cases = n(),
      mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
      median_scaled_err_improvement = median(scaled_abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(feature_catalog_for_wp4() %>% select(feature, display_label, category), by = "feature") %>%
    arrange(scope, category, display_label)

  case_summary <- predictions %>%
    group_by(scope, cellLine, direction) %>%
    summarise(
      n_features = n(),
      mean_scaled_err_improvement = mean(scaled_abs_err_improvement, na.rm = TRUE),
      transfer_win_rate = mean(abs_err_transfer < abs_err_null, na.rm = TRUE),
      sign_accuracy = mean(sign_match, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scope, cellLine, direction)

  matched_comparison <- overall %>%
    filter(scope %in% c("All-line estimator on non-SUM holdouts", "SUM-159-fuse excluded")) %>%
    mutate(scope_key = recode(
      as.character(scope),
      "All-line estimator on non-SUM holdouts" = "all_line_estimator_non_sum_holdouts",
      "SUM-159-fuse excluded" = "no_sum159_fuse"
    )) %>%
    select(
      scope_key, n_predictions, n_cell_lines, mean_scaled_err_improvement,
      median_scaled_err_improvement, transfer_win_rate, sign_accuracy
    ) %>%
    pivot_wider(
      names_from = scope_key,
      values_from = c(
        n_predictions, n_cell_lines, mean_scaled_err_improvement,
        median_scaled_err_improvement, transfer_win_rate, sign_accuracy
      ),
      names_glue = "{.value}__{scope_key}"
    ) %>%
    mutate(
      delta_mean_scaled_err_improvement =
        mean_scaled_err_improvement__no_sum159_fuse -
        mean_scaled_err_improvement__all_line_estimator_non_sum_holdouts,
      delta_transfer_win_rate =
        transfer_win_rate__no_sum159_fuse -
        transfer_win_rate__all_line_estimator_non_sum_holdouts,
      delta_sign_accuracy =
        sign_accuracy__no_sum159_fuse -
        sign_accuracy__all_line_estimator_non_sum_holdouts
    )

  list(
    predictions = predictions,
    overall = overall,
    feature_summary = feature_summary,
    case_summary = case_summary,
    matched_comparison = matched_comparison,
    transfer_all = transfer_all,
    transfer_no_sum = transfer_no_sum
  )
}

load_morphology_context <- function() {
  path <- file.path("data", "report_exports", "wp2_morphology_volume", "wp2_sum159_fuse_size_test.csv")
  if (!file.exists(path)) {
    return(tibble())
  }

  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(df) || !"cellLine" %in% names(df)) {
    return(tibble())
  }

  df %>% filter(cellLine == sum_fuse_line)
}

make_effect_consistency_plot <- function(consistency) {
  consistency %>%
    mutate(
      display_label_wrapped = wrap_label(as.character(display_label), width = 16),
      display_label_wrapped = factor(display_label_wrapped, levels = unique(display_label_wrapped)),
      consensus_label = paste0(consensus_n, "/", n_lines)
    ) %>%
    ggplot(aes(display_label_wrapped, consensus_fraction, fill = consensus_direction)) +
    geom_col(width = 0.72, color = "grey20", linewidth = 0.18) +
    geom_text(aes(label = consensus_label), vjust = -0.25, size = 2.35) +
    facet_wrap(~scope, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.03))) +
    scale_fill_manual(
      values = c(positive = "#B2182B", negative = "#2166AC", tie = "grey70"),
      name = "consensus\nsign"
    ) +
    labs(
      title = "A. Effect sign consistency before and after SUM-159-fuse exclusion",
      x = NULL,
      y = "same-sign line fraction"
    ) +
    wp4_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
}

make_effect_heatmap <- function(paired_effects) {
  paired_effects %>%
    mutate(
      display_label_wrapped = wrap_label(as.character(display_label), width = 16),
      display_label_wrapped = factor(display_label_wrapped, levels = unique(display_label_wrapped))
    ) %>%
    ggplot(aes(display_label_wrapped, cellLine, fill = effect_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_point(
      aes(size = abs(effect_per_ploidy), shape = effect_direction),
      color = "grey12",
      fill = "white",
      stroke = 0.35,
      alpha = 0.85
    ) +
    facet_wrap(~scope, nrow = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "all-line\nrobust z"
    ) +
    scale_size_continuous(range = c(1.1, 4.1), guide = "none") +
    scale_shape_manual(
      values = c("higher in high ploidy" = 24, "lower in high ploidy" = 25, "near zero" = 21),
      name = "raw effect sign"
    ) +
    labs(
      title = "B. Paired high-minus-low ploidy contrasts",
      x = NULL,
      y = NULL
    ) +
    wp4_theme(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
}

make_transfer_summary_plot <- function(overall) {
  overall %>%
    select(scope, mean_scaled_err_improvement, transfer_win_rate, sign_accuracy) %>%
    pivot_longer(cols = -scope, names_to = "metric", values_to = "value") %>%
    mutate(
      metric = recode(
        metric,
        mean_scaled_err_improvement = "Mean scaled error improvement",
        transfer_win_rate = "Transfer win rate",
        sign_accuracy = "Sign accuracy"
      ),
      metric = factor(
        metric,
        levels = c("Mean scaled error improvement", "Transfer win rate", "Sign accuracy")
      )
    ) %>%
    ggplot(aes(scope, value, fill = scope)) +
    geom_hline(
      data = tibble(
        metric = factor(
          c("Mean scaled error improvement", "Transfer win rate", "Sign accuracy"),
          levels = c("Mean scaled error improvement", "Transfer win rate", "Sign accuracy")
        ),
        yint = c(0, 0.5, 0.5)
      ),
      aes(yintercept = yint),
      color = "grey55",
      linewidth = 0.25,
      linetype = "dashed"
    ) +
    geom_col(width = 0.68, color = "grey20", linewidth = 0.18) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(
      values = c(
        "All lines" = "#4D4D4D",
        "All-line estimator on non-SUM holdouts" = "#2C7BB6",
        "SUM-159-fuse excluded" = "#D7191C"
      ),
      guide = "none"
    ) +
    labs(
      title = "Transfer/generalization benchmark",
      x = NULL,
      y = NULL
    ) +
    wp4_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
}

make_transfer_feature_plot <- function(feature_summary) {
  feature_summary %>%
    filter(scope %in% c("All-line estimator on non-SUM holdouts", "SUM-159-fuse excluded")) %>%
    mutate(
      display_label_wrapped = wrap_label(as.character(display_label), width = 16),
      display_label_wrapped = factor(display_label_wrapped, levels = unique(display_label_wrapped))
    ) %>%
    ggplot(aes(display_label_wrapped, mean_scaled_err_improvement, fill = scope)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.25, linetype = "dashed") +
    geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "grey20", linewidth = 0.15) +
    scale_fill_manual(
      values = c(
        "All-line estimator on non-SUM holdouts" = "#2C7BB6",
        "SUM-159-fuse excluded" = "#D7191C"
      ),
      name = NULL
    ) +
    labs(
      title = "Feature-level transfer improvement on non-SUM holdouts",
      x = NULL,
      y = "mean scaled error improvement"
    ) +
    wp4_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
}

format_num <- function(x, digits = 3) {
  if (!length(x) || !is.finite(x[1])) {
    return("NA")
  }
  formatC(x[1], format = "f", digits = digits)
}

write_wp4_text_report <- function(
  path,
  stan_data_path,
  figure_dir,
  output_dir,
  consistency_comparison,
  transfer_overall,
  transfer_matched_comparison,
  morphology_context
) {
  mean_all <- mean(consistency_comparison$consensus_fraction__all_lines, na.rm = TRUE)
  mean_no_sum <- mean(consistency_comparison$consensus_fraction__no_sum159_fuse, na.rm = TRUE)
  n_improved <- sum(consistency_comparison$delta_consensus_fraction > 0, na.rm = TRUE)
  n_worse <- sum(consistency_comparison$delta_consensus_fraction < 0, na.rm = TRUE)

  matched <- transfer_matched_comparison
  delta_transfer <- if (nrow(matched)) matched$delta_mean_scaled_err_improvement[1] else NA_real_
  delta_win <- if (nrow(matched)) matched$delta_transfer_win_rate[1] else NA_real_
  delta_sign <- if (nrow(matched)) matched$delta_sign_accuracy[1] else NA_real_

  morphology_sentence <- if (nrow(morphology_context)) {
    sprintf(
      paste(
        "WP2 morphology context: SUM-159-fuse has log2 attached-area ratio %s",
        "and log2 latent-volume ratio %s; relative to non-SUM lines, z-scores are %s",
        "for attached area and %s for latent volume."
      ),
      format_num(morphology_context$log2_attached_area_ratio),
      format_num(morphology_context$log2_latent_volume_ratio),
      format_num(morphology_context$z_vs_non_sum_q50),
      format_num(morphology_context$z_vs_non_sum_volume)
    )
  } else {
    "WP2 morphology context was not available in data/report_exports/wp2_morphology_volume/wp2_sum159_fuse_size_test.csv."
  }

  rationale_sentence <- paste(
    "The rationale for treating SUM-159-fuse as an exception candidate is not performance-only:",
    "it is the only cell-line dataset with two batches and differing protocol details,",
    "it is a volume outlier, and the transferable mechanistic conclusion",
    "(the direction of ploidy effects on y1, v1, K1, and m) is conserved across tested models",
    "but moves further from zero and shows stronger inter-model agreement after excluding SUM-159-fuse."
  )

  exception_call <- case_when(
    is.finite(delta_transfer) && delta_transfer > 0 && n_improved > n_worse ~
      "The model-free evidence supports treating SUM-159-fuse as an exception candidate.",
    is.finite(delta_transfer) && delta_transfer > 0 ~
      "The transfer benchmark improves after exclusion, but sign-consistency support is mixed.",
    n_improved > n_worse ~
      "Sign consistency improves after exclusion, but transfer support is limited.",
    TRUE ~
      "The model-free evidence does not support excluding SUM-159-fuse on performance grounds alone."
  )

  paragraph <- paste(
    sprintf(
      paste(
        "Excluding SUM-159-fuse changes the mean same-sign effect-consistency fraction from %s to %s",
        "across the curated model-free features (%d improved, %d worsened)."
      ),
      format_num(mean_all),
      format_num(mean_no_sum),
      n_improved,
      n_worse
    ),
    sprintf(
      paste(
        "On matched non-SUM holdouts, the no-SUM estimator changes mean scaled transfer improvement by %s,",
        "transfer win rate by %s, and sign accuracy by %s relative to the all-line estimator."
      ),
      format_num(delta_transfer),
      format_num(delta_win),
      format_num(delta_sign)
    ),
    morphology_sentence,
    rationale_sentence,
    exception_call
  )

  lines <- c(
    "WP4_SUM159_FUSE_EXCLUSION_ANALYSIS",
    sprintf("generated\t%s", Sys.time()),
    sprintf("stan_data\t%s", stan_data_path),
    sprintf("figure_dir\t%s", figure_dir),
    sprintf("table_dir\t%s", output_dir),
    "",
    "SECTION\tSHORT_EXCEPTION_ANALYSIS_PARAGRAPH",
    paragraph,
    "",
    "SECTION\tEFFECT_CONSISTENCY_COMPARISON",
    capture.output(write.table(
      as.data.frame(consistency_comparison),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "SECTION\tTRANSFER_OVERALL",
    capture.output(write.table(
      as.data.frame(transfer_overall),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "SECTION\tMATCHED_TRANSFER_COMPARISON",
    capture.output(write.table(
      as.data.frame(transfer_matched_comparison),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "SECTION\tMORPHOLOGY_CONTEXT",
    if (nrow(morphology_context)) {
      capture.output(write.table(as.data.frame(morphology_context), sep = "\t", row.names = FALSE, quote = FALSE))
    } else {
      "No morphology context table available."
    }
  )

  writeLines(lines, path, useBytes = TRUE)
  invisible(paragraph)
}

write_wp4_markdown_report <- function(
  path,
  paragraph,
  consistency_comparison,
  transfer_overall,
  transfer_matched_comparison
) {
  top_consistency <- consistency_comparison %>%
    arrange(desc(abs(delta_consensus_fraction)), feature) %>%
    select(
      feature, display_label,
      consensus_fraction__all_lines,
      consensus_fraction__no_sum159_fuse,
      delta_consensus_fraction,
      consensus_direction__all_lines,
      consensus_direction__no_sum159_fuse
    )

  lines <- c(
    "# WP4 SUM-159-fuse exclusion analysis",
    "",
    "## Exception analysis paragraph",
    "",
    paragraph,
    "",
    "## Transfer metrics",
    "",
    capture.output(write.table(
      as.data.frame(transfer_overall),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "## Matched non-SUM transfer comparison",
    "",
    capture.output(write.table(
      as.data.frame(transfer_matched_comparison),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )),
    "",
    "## Effect consistency shifts",
    "",
    capture.output(write.table(
      as.data.frame(top_consistency),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    ))
  )

  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

stan_data_path <- resolve_stan_data_path(stan_data_path_arg)
message("Building WP4 model-free tables from ", stan_data_path)

tables <- build_feature_panel(stan_data_path = stan_data_path)
feature_panel_all <- tables$feature_panel
feature_panel_no_sum <- feature_panel_all %>% filter(cellLine != sum_fuse_line)
signature_panel_all <- summarize_glucose_signatures(feature_panel_all)
signature_panel_no_sum <- summarize_glucose_signatures(feature_panel_no_sum)

effect_tables <- build_effect_scope_tables(signature_panel_all)
transfer_tables <- build_transfer_scope_tables(signature_panel_all)
morphology_context <- load_morphology_context()

write.csv(feature_panel_all, file.path(output_dir, "wp4_feature_panel_all_lines.csv"), row.names = FALSE)
write.csv(feature_panel_no_sum, file.path(output_dir, "wp4_feature_panel_no_sum159_fuse.csv"), row.names = FALSE)
write.csv(signature_panel_all, file.path(output_dir, "wp4_signature_panel_all_lines.csv"), row.names = FALSE)
write.csv(signature_panel_no_sum, file.path(output_dir, "wp4_signature_panel_no_sum159_fuse.csv"), row.names = FALSE)
write.csv(effect_tables$paired_effects, file.path(output_dir, "wp4_paired_effects_by_scope.csv"), row.names = FALSE)
write.csv(effect_tables$effect_matrices, file.path(output_dir, "wp4_effect_matrix_by_scope.csv"), row.names = FALSE)
write.csv(effect_tables$consistency, file.path(output_dir, "wp4_effect_consistency_by_scope.csv"), row.names = FALSE)
write.csv(effect_tables$consistency_comparison, file.path(output_dir, "wp4_effect_consistency_comparison.csv"), row.names = FALSE)
write.csv(transfer_tables$predictions, file.path(output_dir, "wp4_transfer_predictions_by_scope.csv"), row.names = FALSE)
write.csv(transfer_tables$overall, file.path(output_dir, "wp4_transfer_overall_by_scope.csv"), row.names = FALSE)
write.csv(transfer_tables$feature_summary, file.path(output_dir, "wp4_transfer_feature_summary_by_scope.csv"), row.names = FALSE)
write.csv(transfer_tables$case_summary, file.path(output_dir, "wp4_transfer_case_summary_by_scope.csv"), row.names = FALSE)
write.csv(transfer_tables$matched_comparison, file.path(output_dir, "wp4_transfer_matched_comparison.csv"), row.names = FALSE)
write.csv(morphology_context, file.path(output_dir, "wp4_sum159_fuse_morphology_context.csv"), row.names = FALSE)

saveRDS(feature_panel_all, file.path(output_dir, "wp4_feature_panel_all_lines.Rds"))
saveRDS(feature_panel_no_sum, file.path(output_dir, "wp4_feature_panel_no_sum159_fuse.Rds"))
saveRDS(signature_panel_all, file.path(output_dir, "wp4_signature_panel_all_lines.Rds"))
saveRDS(signature_panel_no_sum, file.path(output_dir, "wp4_signature_panel_no_sum159_fuse.Rds"))
saveRDS(effect_tables, file.path(output_dir, "wp4_effect_tables.Rds"))
saveRDS(transfer_tables, file.path(output_dir, "wp4_transfer_tables.Rds"))
saveRDS(morphology_context, file.path(output_dir, "wp4_sum159_fuse_morphology_context.Rds"))

effect_figure <- make_effect_consistency_plot(effect_tables$consistency) /
  make_effect_heatmap(effect_tables$paired_effects) +
  plot_layout(heights = c(1, 1.55)) +
  plot_annotation(
    title = "WP4. SUM-159-fuse exclusion: effect consistency",
    caption = "Effects are high-ploidy minus low-ploidy contrasts per ploidy-metric unit. The excluded scope recomputes paired effects after dropping SUM-159-fuse."
  )
save_plot_pair(effect_figure, "wp4_sum159_fuse_effect_consistency", width = 11.8, height = 9)

transfer_figure <- make_transfer_summary_plot(transfer_tables$overall) /
  make_transfer_feature_plot(transfer_tables$feature_summary) +
  plot_layout(heights = c(1, 1.35)) +
  plot_annotation(
    title = "WP4. SUM-159-fuse exclusion: transfer/generalization",
    caption = "The matched comparison contrasts non-SUM holdouts under an all-line estimator versus an estimator trained and evaluated after excluding SUM-159-fuse."
  )
save_plot_pair(transfer_figure, "wp4_sum159_fuse_transfer_generalization", width = 11.8, height = 8.2)

paragraph <- write_wp4_text_report(
  path = file.path(output_dir, "wp4_sum159_fuse_exception_analysis.txt"),
  stan_data_path = stan_data_path,
  figure_dir = figure_dir,
  output_dir = output_dir,
  consistency_comparison = effect_tables$consistency_comparison,
  transfer_overall = transfer_tables$overall,
  transfer_matched_comparison = transfer_tables$matched_comparison,
  morphology_context = morphology_context
)

write_wp4_markdown_report(
  path = file.path(output_dir, "wp4_sum159_fuse_exception_analysis.md"),
  paragraph = paragraph,
  consistency_comparison = effect_tables$consistency_comparison,
  transfer_overall = transfer_tables$overall,
  transfer_matched_comparison = transfer_tables$matched_comparison
)

cat(sprintf("Wrote WP4 figures to %s\n", figure_dir))
cat(sprintf("Wrote WP4 audit tables to %s\n", output_dir))
