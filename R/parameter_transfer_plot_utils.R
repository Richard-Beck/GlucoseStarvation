library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

get_parameter_family <- function(parameter) {
  base <- sub("\\[.*$", "", parameter)

  if (base %in% c("R", "P", "W", "waste_mech")) {
    return(NA_character_)
  }

  if (base %in% c("ae", "ah", "Y_R", "m")) {
    return(base)
  }

  if (grepl("^A_mat", base)) {
    return("A_mat")
  }

  base
}

load_parameter_transfer_outputs <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  base_dir <- file.path(output_root, model_name, model_id)

  comparison_path <- file.path(base_dir, "parameter_transfer_comparison.Rds")
  summary_path <- file.path(base_dir, "parameter_transfer_summary.Rds")

  if (!file.exists(comparison_path)) {
    stop(sprintf("Parameter transfer comparison not found: %s", comparison_path))
  }

  if (!file.exists(summary_path)) {
    stop(sprintf("Parameter transfer summary not found: %s", summary_path))
  }

  list(
    comparison = readRDS(comparison_path),
    summary = readRDS(summary_path)
  )
}

plot_parameter_transfer_improvement <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  dat <- load_parameter_transfer_outputs(
    model_id = model_id,
    output_root = output_root,
    model_name = model_name
  )

  df <- dat$comparison
  df$family <- vapply(df$parameter, get_parameter_family, character(1))
  df <- df[!is.na(df$family), , drop = FALSE]
  df$direction <- factor(df$direction, levels = c("low_to_high", "high_to_low"))
  eps <- 1e-12
  df$abs_null_err <- pmax(abs(df$null_vs_oracle_diff), eps)
  df$abs_transfer_err <- pmax(abs(df$holdout_diff), eps)
  df$log10_error_ratio <- log10(df$abs_transfer_err / df$abs_null_err)

  ggplot(df, aes(x = family, y = log10_error_ratio, color = direction)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.65, size = 1.6) +
    facet_wrap(~ line_id, scales = "free_x") +
    scale_color_manual(values = c(low_to_high = "#1b9e77", high_to_low = "#d95f02")) +
    theme_minimal() +
    labs(
      title = sprintf("Parameter Transfer Improvement Over Null | %s", model_id),
      x = "Parameter family",
      y = "log10(|transfer-oracle| / |null-oracle|)",
      color = "Direction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey85", fill = NA)
    )
}

prepare_parameter_transfer_comparison <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  dat <- load_parameter_transfer_outputs(
    model_id = model_id,
    output_root = output_root,
    model_name = model_name
  )

  df <- dat$comparison
  df$family <- vapply(df$parameter, get_parameter_family, character(1))
  df <- df[!is.na(df$family), , drop = FALSE]
  df$direction <- factor(df$direction, levels = c("low_to_high", "high_to_low"))
  eps <- 1e-12
  df$abs_null_err <- pmax(abs(df$null_vs_oracle_diff), eps)
  df$abs_transfer_err <- pmax(abs(df$holdout_diff), eps)
  df$log10_error_ratio <- log10(df$abs_transfer_err / df$abs_null_err)
  df$log_error_improvement <- log(df$abs_null_err) - log(df$abs_transfer_err)
  df$transfer_better_than_null <- df$abs_transfer_err < df$abs_null_err
  df
}

plot_parameter_transfer_summary <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  dat <- load_parameter_transfer_outputs(
    model_id = model_id,
    output_root = output_root,
    model_name = model_name
  )

  df <- dat$summary
  df$direction <- factor(df$direction, levels = c("low_to_high", "high_to_low"))
  df$line_direction <- sprintf("line %d | %s", df$line_id, df$direction)
  df$line_direction <- factor(df$line_direction, levels = df$line_direction[order(df$direction, df$line_id)])
  df$mean_log_error_improvement <- df$mean_abs_log_null_holdout_ratio - df$mean_abs_log_holdout_ratio

  p1 <- ggplot(df, aes(x = line_direction, y = mean_log_error_improvement, fill = direction)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c(low_to_high = "#1b9e77", high_to_low = "#d95f02")) +
    theme_minimal() +
    labs(
      title = sprintf("Parameter Transfer Log-Error Improvement | %s", model_id),
      x = "",
      y = "Mean log-error improvement over null",
      fill = "Direction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey85", fill = NA)
    )

  p2 <- ggplot(df, aes(x = line_direction, y = prop_parameters_transfer_better_than_null, fill = direction)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c(low_to_high = "#1b9e77", high_to_low = "#d95f02")) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal() +
    labs(
      title = sprintf("Parameter-Wise Win Rate | %s", model_id),
      x = "",
      y = "Proportion transfer better than null",
      fill = "Direction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey85", fill = NA)
    )

  p1 / p2
}

plot_parameter_transfer_family_summary <- function(
  model_id,
  output_root = "data/gpath_transfer_cv",
  model_name = "gpath"
) {
  df <- prepare_parameter_transfer_comparison(
    model_id = model_id,
    output_root = output_root,
    model_name = model_name
  )

  family_df <- df %>%
    group_by(line_id, direction, family) %>%
    summarise(
      median_log10_error_ratio = median(log10_error_ratio, na.rm = TRUE),
      prop_transfer_better_than_null = mean(transfer_better_than_null, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      line_direction = sprintf("line %d | %s", line_id, direction),
      line_direction = factor(line_direction, levels = unique(line_direction[order(direction, line_id)]))
    )

  p1 <- ggplot(family_df, aes(x = family, y = line_direction, fill = median_log10_error_ratio)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#1b9e77",
      mid = "white",
      high = "#d95f02",
      midpoint = 0
    ) +
    theme_minimal() +
    labs(
      title = sprintf("Family-Level Relative Error Shift | %s", model_id),
      x = "Parameter family",
      y = "",
      fill = "Median log10\n(transfer/null error)"
    ) +
    theme(panel.border = element_rect(color = "grey85", fill = NA))

  p2 <- ggplot(family_df, aes(x = family, y = line_direction, fill = prop_transfer_better_than_null)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#d95f02",
      mid = "white",
      high = "#1b9e77",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    theme_minimal() +
    labs(
      title = sprintf("Family-Level Win Rate | %s", model_id),
      x = "Parameter family",
      y = "",
      fill = "Prop.\ntransfer better"
    ) +
    theme(panel.border = element_rect(color = "grey85", fill = NA))

  p1 / p2
}
