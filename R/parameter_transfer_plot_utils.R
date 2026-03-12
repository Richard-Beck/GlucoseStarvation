library(ggplot2)
library(dplyr)

get_parameter_family <- function(parameter) {
  base <- sub("\\[.*$", "", parameter)

  if (base %in% c("ae", "ah", "Y_R", "m", "R", "P", "W")) {
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
  df$direction <- factor(df$direction, levels = c("low_to_high", "high_to_low"))

  ggplot(df, aes(x = family, y = transfer_improvement_over_null, color = direction)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.65, size = 1.6) +
    facet_wrap(~ line_id, scales = "free_x") +
    scale_color_manual(values = c(low_to_high = "#1b9e77", high_to_low = "#d95f02")) +
    theme_minimal() +
    labs(
      title = sprintf("Parameter Transfer Improvement Over Null | %s", model_id),
      x = "Parameter family",
      y = "Improvement over null\n(abs(null-oracle) - abs(transfer-oracle))",
      color = "Direction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey85", fill = NA)
    )
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
  scale_ref <- max(abs(df$mean_transfer_improvement_over_null), na.rm = TRUE)
  if (!is.finite(scale_ref) || scale_ref <= 0) {
    scale_ref <- 1
  }

  ggplot(df, aes(x = line_direction, y = mean_transfer_improvement_over_null, fill = direction)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_col(width = 0.75) +
    geom_point(aes(y = prop_parameters_transfer_better_than_null * scale_ref),
               color = "black", size = 2) +
    scale_fill_manual(values = c(low_to_high = "#1b9e77", high_to_low = "#d95f02")) +
    theme_minimal() +
    labs(
      title = sprintf("Parameter Transfer Summary | %s", model_id),
      x = "",
      y = "Mean improvement over null",
      fill = "Direction",
      caption = "Black points are scaled from prop_parameters_transfer_better_than_null."
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey85", fill = NA)
    )
}
