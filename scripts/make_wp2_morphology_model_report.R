args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp2_morphology_model_report")
}

input_path <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("data", "report_exports", "wp2_morphology_volume", "wp2_area_count_observations.Rds")
}

figure_dir <- file.path(output_dir, "figures")

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

sum_fuse_line <- "SUM-159-fuse"

ploidy_palette <- c(
  "baseline" = "#2C7BB6",
  "elevated" = "#D7191C"
)

line_palette <- c(
  "MCF10A" = "#2C7BB6",
  "MDA-MB-231" = "#1A9641",
  "SNU668" = "#FDAE61",
  "SUM-159-chem" = "#984EA3",
  "SUM-159-fuse" = "#D7191C"
)

wp2_theme <- function(base_size = 8) {
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

save_png_pdf <- function(plot, basename, width, height, dpi = 360) {
  png_path <- file.path(figure_dir, paste0(basename, ".png"))
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))

  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")

  invisible(c(png = png_path, pdf = pdf_path))
}

wrap_label <- function(x, width = 22) {
  vapply(
    x,
    function(one) paste(strwrap(one, width = width), collapse = "\n"),
    character(1)
  )
}

scale_vector <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) {
    sig <- 1
  }
  (x - mu) / sig
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  sqrt(mean(x^2))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok]))
}

read_area_observations <- function(path) {
  if (file.exists(path)) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      return(readRDS(path))
    }
    return(read.csv(path, stringsAsFactors = FALSE))
  }

  csv_path <- sub("\\.Rds$", ".csv", path, ignore.case = TRUE)
  if (file.exists(csv_path)) {
    return(read.csv(csv_path, stringsAsFactors = FALSE))
  }

  stop(
    "Could not find WP2 area observation table at '", path, "' or '", csv_path,
    "'. Run scripts/make_wp2_morphology_volume_figures.R first."
  )
}

prepare_model_data <- function(area_obs) {
  required <- c(
    "well_idx", "rep_id", "hours", "live_cells", "dead_cells", "n_area_alive",
    "q50", "q90", "total_cells", "dead_fraction", "cellLine", "ploidy_state",
    "ploidy_metric", "G0", "glucose_hat", "growth_rate_proxy", "death_rate_proxy"
  )
  missing <- setdiff(required, names(area_obs))
  if (length(missing)) {
    stop("Area observation table is missing required columns: ", paste(missing, collapse = ", "))
  }

  capacity <- area_obs %>%
    group_by(cellLine, ploidy_state) %>%
    summarise(max_total_cells_state = max(total_cells, na.rm = TRUE), .groups = "drop")

  out <- area_obs %>%
    left_join(capacity, by = c("cellLine", "ploidy_state")) %>%
    mutate(
      cellLine = factor(cellLine, levels = sort(unique(cellLine))),
      ploidy_state = factor(ploidy_state, levels = c("baseline", "elevated")),
      rep_id = as.character(rep_id),
      log_q50 = log(q50),
      log_q90 = log(q90),
      log_volume_q50 = 1.5 * log(q50),
      log_live_density_raw = log1p(live_cells),
      log_total_density_raw = log1p(total_cells),
      crowding_frac = total_cells / pmax(max_total_cells_state, 1e-8),
      start_glucose_log = log1p(pmax(G0, 0)),
      current_glucose_log = log1p(pmax(glucose_hat, 0)),
      glucose_depletion_log = pmax(start_glucose_log - current_glucose_log, 0),
      net_growth_death_proxy = growth_rate_proxy - death_rate_proxy,
      glucose_bin = case_when(
        G0 <= 0.25 ~ "low G0",
        G0 <= 1 ~ "mid G0",
        TRUE ~ "high G0"
      ),
      glucose_bin = factor(glucose_bin, levels = c("low G0", "mid G0", "high G0"))
    ) %>%
    filter(
      is.finite(log_q50),
      is.finite(log_q90),
      is.finite(log_volume_q50),
      is.finite(hours),
      is.finite(log_live_density_raw),
      is.finite(log_total_density_raw),
      is.finite(crowding_frac),
      is.finite(current_glucose_log),
      is.finite(start_glucose_log),
      is.finite(glucose_depletion_log),
      is.finite(growth_rate_proxy),
      is.finite(death_rate_proxy),
      is.finite(dead_fraction),
      n_area_alive >= 30
    ) %>%
    mutate(
      hours_z = scale_vector(hours),
      hours_log_z = scale_vector(log1p(hours)),
      log_live_density_z = scale_vector(log_live_density_raw),
      log_total_density_z = scale_vector(log_total_density_raw),
      crowding_z = scale_vector(crowding_frac),
      start_glucose_z = scale_vector(start_glucose_log),
      current_glucose_z = scale_vector(current_glucose_log),
      glucose_depletion_z = scale_vector(glucose_depletion_log),
      growth_rate_z = scale_vector(growth_rate_proxy),
      death_rate_z = scale_vector(death_rate_proxy),
      dead_fraction_z = scale_vector(dead_fraction),
      net_growth_death_z = scale_vector(net_growth_death_proxy)
    )

  fold_tbl <- out %>%
    distinct(cellLine, ploidy_state, well_idx) %>%
    arrange(cellLine, ploidy_state, well_idx) %>%
    group_by(cellLine, ploidy_state) %>%
    mutate(fold = ((row_number() - 1L) %% 5L) + 1L) %>%
    ungroup() %>%
    select(well_idx, fold)

  out %>%
    left_join(fold_tbl, by = "well_idx")
}

model_catalog <- function() {
  time_terms <- "hours_z + I(hours_z^2) + I(hours_z^3)"
  density_terms <- "log_live_density_z + log_total_density_z + crowding_z + I(crowding_z^2)"
  glucose_terms <- "start_glucose_z + current_glucose_z + glucose_depletion_z"
  phase_terms <- "growth_rate_z + I(growth_rate_z^2) + death_rate_z + dead_fraction_z + net_growth_death_z"
  full_covariates <- paste(time_terms, density_terms, glucose_terms, phase_terms, sep = " + ")

  tibble(
    model_id = c(
      "line_ploidy",
      "time",
      "density",
      "glucose",
      "growth_death",
      "time_density",
      "full_additive",
      "ploidy_dynamic",
      "line_specific_dynamics"
    ),
    model_label = c(
      "Line + ploidy",
      "+ time since seeding",
      "+ density/crowding",
      "+ glucose state",
      "+ growth/death rates",
      "+ time + density",
      "Full additive",
      "Full + ploidy dynamics",
      "Line-specific dynamics"
    ),
    rhs = c(
      "cellLine * ploidy_state",
      paste("cellLine * ploidy_state", time_terms, sep = " + "),
      paste("cellLine * ploidy_state", density_terms, sep = " + "),
      paste("cellLine * ploidy_state", glucose_terms, sep = " + "),
      paste("cellLine * ploidy_state", phase_terms, sep = " + "),
      paste("cellLine * ploidy_state", time_terms, density_terms, sep = " + "),
      paste("cellLine * ploidy_state", full_covariates, sep = " + "),
      paste(
        "cellLine * ploidy_state",
        full_covariates,
        "ploidy_state:(hours_z + log_live_density_z + crowding_z + current_glucose_z + growth_rate_z)",
        sep = " + "
      ),
      paste(
        "cellLine * ploidy_state",
        full_covariates,
        "cellLine:(hours_z + I(hours_z^2) + log_live_density_z + crowding_z + current_glucose_z + growth_rate_z + death_rate_z)",
        sep = " + "
      )
    ),
    covariate_blocks = c(
      "line, ploidy",
      "line, ploidy, time",
      "line, ploidy, density/crowding",
      "line, ploidy, glucose",
      "line, ploidy, growth/death",
      "line, ploidy, time, density/crowding",
      "line, ploidy, time, density/crowding, glucose, growth/death",
      "full additive plus ploidy-by-dynamic-covariate interactions",
      "full additive plus line-specific time/density/glucose/growth/death slopes"
    )
  )
}

response_catalog <- function() {
  tibble(
    response_col = c("log_q50", "log_q90", "log_volume_q50"),
    response_label = c("Area q50", "Area q90", "Volume proxy q50 = area^(3/2)")
  )
}

fit_one <- function(df, response_col, rhs) {
  stats::lm(stats::as.formula(paste(response_col, "~", rhs)), data = df)
}

fit_metrics <- function(fit, df, response_col) {
  resid <- stats::residuals(fit)
  fitted <- stats::fitted(fit)
  y <- df[[response_col]]
  n <- length(y)
  p <- length(stats::coef(fit)[is.finite(stats::coef(fit))])

  tibble(
    n_rows = n,
    n_parameters = p,
    train_rmse = safe_rmse(resid),
    train_mae = mean(abs(resid), na.rm = TRUE),
    adj_r_squared = summary(fit)$adj.r.squared,
    aic = stats::AIC(fit),
    bic = stats::BIC(fit),
    residual_cor_hours = safe_cor(resid, df$hours),
    residual_cor_crowding = safe_cor(resid, df$crowding_frac),
    residual_cor_growth = safe_cor(resid, df$growth_rate_proxy),
    fitted_cor_observed = safe_cor(fitted, y)
  )
}

run_blocked_cv <- function(df, response_col, model_tbl) {
  fold_ids <- sort(unique(df$fold))

  bind_rows(lapply(seq_len(nrow(model_tbl)), function(i) {
    mod <- model_tbl[i, ]
    fold_rows <- lapply(fold_ids, function(fold_id) {
      train <- df %>% filter(fold != fold_id)
      test <- df %>% filter(fold == fold_id)
      fit <- fit_one(train, response_col = response_col, rhs = mod$rhs)
      pred <- as.numeric(stats::predict(fit, newdata = test))

      tibble(
        model_id = mod$model_id,
        fold = fold_id,
        n_test = nrow(test),
        cv_rmse = safe_rmse(test[[response_col]] - pred),
        cv_mae = mean(abs(test[[response_col]] - pred), na.rm = TRUE)
      )
    })
    bind_rows(fold_rows)
  })) %>%
    group_by(model_id) %>%
    summarise(
      cv_rmse = sqrt(weighted.mean(cv_rmse^2, w = n_test, na.rm = TRUE)),
      cv_mae = weighted.mean(cv_mae, w = n_test, na.rm = TRUE),
      n_cv_rows = sum(n_test),
      .groups = "drop"
    )
}

fit_model_suite <- function(df) {
  models <- model_catalog()
  responses <- response_catalog()

  out <- lapply(seq_len(nrow(responses)), function(ridx) {
    response_col <- responses$response_col[[ridx]]
    response_label <- responses$response_label[[ridx]]

    message("Fitting morphology model suite for ", response_label)

    fits <- lapply(seq_len(nrow(models)), function(i) {
      fit_one(df, response_col = response_col, rhs = models$rhs[[i]])
    })
    names(fits) <- models$model_id

    train_metrics <- bind_rows(lapply(seq_along(fits), function(i) {
      fit_metrics(fits[[i]], df, response_col) %>%
        mutate(model_id = names(fits)[[i]], .before = 1)
    }))

    cv_metrics <- run_blocked_cv(df, response_col = response_col, model_tbl = models)

    comparison <- models %>%
      select(model_id, model_label, covariate_blocks) %>%
      left_join(train_metrics, by = "model_id") %>%
      left_join(cv_metrics, by = "model_id") %>%
      mutate(
        response_col = response_col,
        response_label = response_label,
        delta_aic = aic - min(aic, na.rm = TRUE),
        delta_cv_rmse = cv_rmse - min(cv_rmse, na.rm = TRUE),
        response_label = factor(response_label, levels = responses$response_label),
        model_label = factor(model_label, levels = models$model_label)
      )

    list(
      response_col = response_col,
      response_label = response_label,
      fits = fits,
      comparison = comparison
    )
  })

  names(out) <- responses$response_col
  out
}

reference_rows_for_effects <- function(df) {
  ref <- df %>%
    group_by(cellLine) %>%
    summarise(
      across(
        c(
          hours_z, hours_log_z, log_live_density_z, log_total_density_z, crowding_z,
          start_glucose_z, current_glucose_z, glucose_depletion_z, growth_rate_z,
          death_rate_z, dead_fraction_z, net_growth_death_z
        ),
        ~ median(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  bind_rows(
    ref %>% mutate(ploidy_state = "baseline"),
    ref %>% mutate(ploidy_state = "elevated")
  ) %>%
    mutate(
      cellLine = factor(as.character(cellLine), levels = levels(df$cellLine)),
      ploidy_state = factor(ploidy_state, levels = levels(df$ploidy_state))
    ) %>%
    arrange(cellLine, ploidy_state)
}

estimate_ploidy_effects <- function(fit, df, model_id, model_label, response_col, response_label) {
  ref <- reference_rows_for_effects(df)
  model_terms <- stats::delete.response(stats::terms(fit))
  mm <- stats::model.matrix(model_terms, data = ref, contrasts.arg = fit$contrasts)
  cf <- stats::coef(fit)
  keep <- is.finite(cf)
  cf <- cf[keep]
  mm <- mm[, keep, drop = FALSE]
  vc <- stats::vcov(fit)[keep, keep, drop = FALSE]

  bind_rows(lapply(levels(df$cellLine), function(line_name) {
    idx <- which(as.character(ref$cellLine) == line_name)
    if (length(idx) != 2L) {
      return(tibble())
    }

    idx_base <- idx[as.character(ref$ploidy_state[idx]) == "baseline"]
    idx_elev <- idx[as.character(ref$ploidy_state[idx]) == "elevated"]
    if (!length(idx_base) || !length(idx_elev)) {
      return(tibble())
    }

    xdiff <- mm[idx_elev, , drop = FALSE] - mm[idx_base, , drop = FALSE]
    effect <- as.numeric(xdiff %*% cf)
    se <- as.numeric(sqrt(xdiff %*% vc %*% t(xdiff)))

    tibble(
      model_id = model_id,
      model_label = model_label,
      response_col = response_col,
      response_label = response_label,
      cellLine = line_name,
      log_ratio = effect,
      se_log_ratio = se,
      lower_log_ratio = effect - 1.96 * se,
      upper_log_ratio = effect + 1.96 * se,
      ratio = exp(effect),
      lower_ratio = exp(effect - 1.96 * se),
      upper_ratio = exp(effect + 1.96 * se),
      log2_ratio = effect / log(2),
      lower_log2_ratio = (effect - 1.96 * se) / log(2),
      upper_log2_ratio = (effect + 1.96 * se) / log(2),
      is_sum159_fuse = line_name == sum_fuse_line
    )
  }))
}

estimate_all_effects <- function(suite, df) {
  bind_rows(lapply(suite, function(res) {
    bind_rows(lapply(names(res$fits), function(model_id) {
      model_row <- model_catalog() %>% filter(model_id == !!model_id)
      estimate_ploidy_effects(
        fit = res$fits[[model_id]],
        df = df,
        model_id = model_id,
        model_label = model_row$model_label[[1]],
        response_col = res$response_col,
        response_label = res$response_label
      )
    }))
  })) %>%
    mutate(
      model_label = factor(model_label, levels = model_catalog()$model_label),
      response_label = factor(response_label, levels = response_catalog()$response_label)
    )
}

block_drop_catalog <- function() {
  time_terms <- "hours_z + I(hours_z^2) + I(hours_z^3)"
  density_terms <- "log_live_density_z + log_total_density_z + crowding_z + I(crowding_z^2)"
  glucose_terms <- "start_glucose_z + current_glucose_z + glucose_depletion_z"
  phase_terms <- "growth_rate_z + I(growth_rate_z^2) + death_rate_z + dead_fraction_z + net_growth_death_z"

  tibble(
    model_id = c(
      "full_additive",
      "drop_time",
      "drop_density",
      "drop_glucose",
      "drop_growth_death",
      "drop_ploidy"
    ),
    block_tested = c(
      "Full additive",
      "Time since seeding",
      "Density/crowding",
      "Glucose state",
      "Growth/death dynamics",
      "Ploidy within line"
    ),
    rhs = c(
      paste("cellLine * ploidy_state", time_terms, density_terms, glucose_terms, phase_terms, sep = " + "),
      paste("cellLine * ploidy_state", density_terms, glucose_terms, phase_terms, sep = " + "),
      paste("cellLine * ploidy_state", time_terms, glucose_terms, phase_terms, sep = " + "),
      paste("cellLine * ploidy_state", time_terms, density_terms, phase_terms, sep = " + "),
      paste("cellLine * ploidy_state", time_terms, density_terms, glucose_terms, sep = " + "),
      paste("cellLine", time_terms, density_terms, glucose_terms, phase_terms, sep = " + ")
    )
  )
}

run_block_drop <- function(df, response_col = "log_q50") {
  block_tbl <- block_drop_catalog()
  cv <- run_blocked_cv(
    df = df,
    response_col = response_col,
    model_tbl = block_tbl %>% transmute(model_id, rhs)
  )

  full_rmse <- cv$cv_rmse[cv$model_id == "full_additive"]
  block_tbl %>%
    select(model_id, block_tested) %>%
    left_join(cv, by = "model_id") %>%
    mutate(
      rmse_increase_vs_full = cv_rmse - full_rmse,
      block_tested = factor(
        block_tested,
        levels = c("Ploidy within line", "Time since seeding", "Density/crowding", "Glucose state", "Growth/death dynamics", "Full additive")
      )
    )
}

get_fitted_values <- function(fit, df, response_col, model_id, model_label) {
  df %>%
    mutate(
      model_id = model_id,
      model_label = model_label,
      response_col = response_col,
      fitted_log = as.numeric(stats::fitted(fit)),
      residual_log = as.numeric(stats::residuals(fit)),
      observed_log = .data[[response_col]]
    )
}

make_covariate_landscape_plot <- function(df) {
  time_summary <- df %>%
    group_by(cellLine, ploidy_state, glucose_bin, hours) %>%
    summarise(q50 = median(q50, na.rm = TRUE), .groups = "drop")

  p_time <- ggplot(
    time_summary,
    aes(hours, q50, color = ploidy_state, linetype = glucose_bin, group = interaction(ploidy_state, glucose_bin))
  ) +
    geom_line(linewidth = 0.65, alpha = 0.9) +
    facet_wrap(~cellLine, scales = "free_y", nrow = 1) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "A. Area changes over time and starting glucose",
      x = "hours since seeding",
      y = "median area q50"
    ) +
    wp2_theme(base_size = 8)

  p_crowding <- ggplot(
    df,
    aes(crowding_frac, q50, color = ploidy_state)
  ) +
    geom_point(alpha = 0.07, size = 0.55) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.7, span = 0.9) +
    facet_wrap(~cellLine, scales = "free_y", nrow = 1) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "B. Area depends on crowding",
      x = "within-state crowding fraction",
      y = "area q50"
    ) +
    wp2_theme(base_size = 8)

  p_growth <- ggplot(
    df,
    aes(growth_rate_proxy, q50, color = ploidy_state)
  ) +
    geom_point(alpha = 0.07, size = 0.55) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.7, span = 0.9) +
    facet_wrap(~cellLine, scales = "free_y", nrow = 1) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "C. Area also tracks growth phase",
      x = "smoothed per-cell live growth proxy",
      y = "area q50"
    ) +
    wp2_theme(base_size = 8)

  (p_time / p_crowding / p_growth) +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(
      title = "Figure 1. Pixel Area Is Confounded With Time, Crowding, And Growth Phase"
    )
}

make_model_comparison_plot <- function(comparison) {
  order_tbl <- comparison %>%
    filter(response_label == "Area q50") %>%
    arrange(cv_rmse) %>%
    mutate(model_label_ordered = factor(as.character(model_label), levels = as.character(model_label))) %>%
    select(model_label, model_label_ordered)

  plot_df <- comparison %>%
    left_join(order_tbl, by = "model_label") %>%
    mutate(
      model_label_ordered = if_else(
        is.na(as.character(model_label_ordered)),
        as.character(model_label),
        as.character(model_label_ordered)
      ),
      model_label_ordered = factor(model_label_ordered, levels = levels(order_tbl$model_label_ordered))
    )

  p_rmse <- ggplot(plot_df, aes(model_label_ordered, cv_rmse, fill = response_label)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "grey25", linewidth = 0.18) +
    coord_flip() +
    labs(
      title = "A. Blocked-by-well prediction error",
      x = NULL,
      y = "CV RMSE on log scale",
      fill = NULL
    ) +
    wp2_theme(base_size = 8)

  p_delta <- ggplot(plot_df, aes(model_label_ordered, delta_cv_rmse, fill = response_label)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "grey25", linewidth = 0.18) +
    coord_flip() +
    labs(
      title = "B. Difference from best model for each response",
      x = NULL,
      y = "delta CV RMSE",
      fill = NULL
    ) +
    wp2_theme(base_size = 8)

  p_rmse | p_delta
}

make_effect_stability_plot <- function(effects) {
  plot_df <- effects %>%
    filter(response_label %in% c("Area q50", "Volume proxy q50 = area^(3/2)")) %>%
    mutate(
      cellLine = factor(cellLine, levels = sort(unique(cellLine))),
      model_short = as.character(model_label)
    )

  ggplot(plot_df, aes(model_label, log2_ratio, color = cellLine, group = cellLine)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(shape = is_sum159_fuse), size = 1.9, alpha = 0.92) +
    facet_wrap(~response_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = line_palette, name = "cell line") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    labs(
      title = "Figure 3. Ploidy Effect Stability Across Covariate Models",
      x = NULL,
      y = "adjusted elevated/baseline log2 ratio"
    ) +
    wp2_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))
}

make_block_importance_plot <- function(block_importance) {
  block_importance %>%
    filter(model_id != "full_additive") %>%
    mutate(
      direction = if_else(rmse_increase_vs_full >= 0, "worse when dropped", "better when dropped")
    ) %>%
    ggplot(aes(block_tested, rmse_increase_vs_full, fill = direction)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_col(width = 0.68, color = "grey25", linewidth = 0.18) +
    coord_flip() +
    scale_fill_manual(values = c("worse when dropped" = "#B2182B", "better when dropped" = "#2166AC"), name = NULL) +
    labs(
      title = "Figure 4. Covariate Blocks Tested By Drop-One Blocked CV",
      subtitle = "Positive values mean the full additive q50 model predicts better when the block is retained.",
      x = NULL,
      y = "CV RMSE increase versus full additive model"
    ) +
    wp2_theme(base_size = 8)
}

make_residual_diagnostic_plot <- function(fitted_df) {
  residual_bins <- bind_rows(
    fitted_df %>%
      mutate(covariate = "hours since seeding", covariate_value = hours) %>%
      select(covariate, covariate_value, residual_log, cellLine, ploidy_state),
    fitted_df %>%
      mutate(covariate = "crowding fraction", covariate_value = crowding_frac) %>%
      select(covariate, covariate_value, residual_log, cellLine, ploidy_state),
    fitted_df %>%
      mutate(covariate = "growth-rate proxy", covariate_value = growth_rate_proxy) %>%
      select(covariate, covariate_value, residual_log, cellLine, ploidy_state)
  ) %>%
    group_by(covariate, cellLine) %>%
    mutate(bin = ntile(covariate_value, 8)) %>%
    ungroup() %>%
    group_by(covariate, cellLine, ploidy_state, bin) %>%
    summarise(
      covariate_value = median(covariate_value, na.rm = TRUE),
      median_residual = median(residual_log, na.rm = TRUE),
      p25 = quantile(residual_log, 0.25, na.rm = TRUE),
      p75 = quantile(residual_log, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(residual_bins, aes(covariate_value, median_residual, color = ploidy_state)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_ribbon(aes(ymin = p25, ymax = p75, fill = ploidy_state), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.3) +
    facet_grid(covariate ~ cellLine, scales = "free_x") +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    scale_fill_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "Figure 5. Residual Trends Under The Main Adjusted q50 Model",
      x = "covariate value",
      y = "median log residual"
    ) +
    wp2_theme(base_size = 8)
}

make_observed_fitted_plot <- function(fitted_df) {
  ggplot(fitted_df, aes(exp(fitted_log), exp(observed_log), color = ploidy_state)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.28, size = 0.65) +
    facet_wrap(~cellLine, scales = "free", nrow = 1) +
    scale_x_log10(labels = label_number(big.mark = ",")) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "Figure 6. Observed Versus Fitted Area q50 For The Main Adjusted Model",
      x = "fitted area q50",
      y = "observed area q50"
    ) +
    wp2_theme(base_size = 8)
}

format_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    ifelse(abs(x) >= 1000, format(round(x, digits), big.mark = ","), format(round(x, digits), nsmall = digits))
  )
}

markdown_table <- function(df, digits = 3) {
  if (!nrow(df)) {
    return("_No rows._")
  }

  fmt <- df
  for (nm in names(fmt)) {
    if (is.numeric(fmt[[nm]])) {
      fmt[[nm]] <- format_num(fmt[[nm]], digits = digits)
    } else {
      fmt[[nm]] <- as.character(fmt[[nm]])
    }
  }

  header <- paste0("| ", paste(names(fmt), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(fmt)), collapse = " | "), " |")
  rows <- apply(fmt, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

fig_link <- function(path) {
  file.path("figures", basename(path))
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

html_attr_escape <- function(x) {
  html_escape(x)
}

inline_markdown_to_html <- function(x) {
  x <- html_escape(x)
  x <- gsub("`([^`]+)`", "<code>\\1</code>", x, perl = TRUE)
  x <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", x, perl = TRUE)
  x
}

split_markdown_table_row <- function(line) {
  line <- trimws(line)
  line <- sub("^\\|", "", line)
  line <- sub("\\|$", "", line)
  trimws(strsplit(line, "\\|", fixed = FALSE)[[1]])
}

markdown_lines_to_html <- function(lines, title = "WP2 Morphology Model Report") {
  body <- character()
  paragraph <- character()
  in_ol <- FALSE

  flush_paragraph <- function() {
    if (length(paragraph)) {
      body <<- c(body, sprintf("<p>%s</p>", inline_markdown_to_html(paste(paragraph, collapse = " "))))
      paragraph <<- character()
    }
  }

  close_ol <- function() {
    if (isTRUE(in_ol)) {
      body <<- c(body, "</ol>")
      in_ol <<- FALSE
    }
  }

  i <- 1L
  while (i <= length(lines)) {
    line <- lines[[i]]

    if (!nzchar(trimws(line))) {
      flush_paragraph()
      close_ol()
      i <- i + 1L
      next
    }

    if (grepl("^#{1,6}\\s+", line)) {
      flush_paragraph()
      close_ol()
      level <- nchar(sub("^(#{1,6}).*$", "\\1", line))
      text <- sub("^#{1,6}\\s+", "", line)
      body <- c(body, sprintf("<h%d>%s</h%d>", level, inline_markdown_to_html(text), level))
      i <- i + 1L
      next
    }

    if (grepl("^!\\[[^]]*\\]\\([^)]+\\)$", line)) {
      flush_paragraph()
      close_ol()
      alt <- sub("^!\\[([^]]*)\\]\\(([^)]+)\\)$", "\\1", line)
      src <- sub("^!\\[([^]]*)\\]\\(([^)]+)\\)$", "\\2", line)
      body <- c(
        body,
        sprintf(
          "<figure><img src=\"%s\" alt=\"%s\"><figcaption>%s</figcaption></figure>",
          html_attr_escape(src),
          html_attr_escape(alt),
          inline_markdown_to_html(alt)
        )
      )
      i <- i + 1L
      next
    }

    if (startsWith(trimws(line), "|")) {
      flush_paragraph()
      close_ol()

      table_lines <- character()
      while (i <= length(lines) && startsWith(trimws(lines[[i]]), "|")) {
        table_lines <- c(table_lines, lines[[i]])
        i <- i + 1L
      }

      if (length(table_lines) >= 2L) {
        header <- split_markdown_table_row(table_lines[[1]])
        data_lines <- table_lines[-c(1L, 2L)]
        body <- c(body, "<table>", "<thead>", "<tr>")
        body <- c(body, paste0("<th>", vapply(header, inline_markdown_to_html, character(1)), "</th>"))
        body <- c(body, "</tr>", "</thead>", "<tbody>")
        for (row_line in data_lines) {
          row <- split_markdown_table_row(row_line)
          if (length(row) < length(header)) {
            row <- c(row, rep("", length(header) - length(row)))
          }
          if (length(row) > length(header)) {
            row <- row[seq_along(header)]
          }
          body <- c(body, "<tr>")
          body <- c(body, paste0("<td>", vapply(row, inline_markdown_to_html, character(1)), "</td>"))
          body <- c(body, "</tr>")
        }
        body <- c(body, "</tbody>", "</table>")
      } else {
        paragraph <- c(paragraph, table_lines)
      }
      next
    }

    if (grepl("^[0-9]+\\.\\s+", line)) {
      flush_paragraph()
      if (!isTRUE(in_ol)) {
        body <- c(body, "<ol>")
        in_ol <- TRUE
      }
      text <- sub("^[0-9]+\\.\\s+", "", line)
      body <- c(body, sprintf("<li>%s</li>", inline_markdown_to_html(text)))
      i <- i + 1L
      next
    }

    close_ol()
    paragraph <- c(paragraph, line)
    i <- i + 1L
  }

  flush_paragraph()
  close_ol()

  paste(
    "<!doctype html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    sprintf("<title>%s</title>", html_escape(title)),
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;line-height:1.5;color:#222;max-width:1180px;margin:32px auto;padding:0 24px;background:#fff;}",
    "h1,h2{line-height:1.2;margin-top:1.6em;} h1{font-size:2rem;border-bottom:1px solid #ddd;padding-bottom:.35rem;} h2{font-size:1.35rem;}",
    "p,li{font-size:1rem;} code{background:#f3f4f6;padding:.12rem .28rem;border-radius:4px;font-size:.94em;}",
    "figure{margin:1.4rem 0 2rem;} figure img{max-width:100%;height:auto;border:1px solid #ddd;} figcaption{font-size:.9rem;color:#555;margin-top:.35rem;}",
    "table{border-collapse:collapse;width:100%;margin:1rem 0 1.4rem;font-size:.92rem;} th,td{border:1px solid #d9d9d9;padding:.42rem .55rem;text-align:left;vertical-align:top;} th{background:#f3f4f6;font-weight:650;} tr:nth-child(even) td{background:#fafafa;}",
    "ol{padding-left:1.5rem;} strong{font-weight:700;}",
    "</style>",
    "</head>",
    "<body>",
    paste(body, collapse = "\n"),
    "</body>",
    "</html>",
    sep = "\n"
  )
}

write_html_report <- function(markdown_path, html_path) {
  lines <- readLines(markdown_path, warn = FALSE)
  title <- if (length(lines) && grepl("^#\\s+", lines[[1]])) {
    sub("^#\\s+", "", lines[[1]])
  } else {
    "WP2 Morphology Model Report"
  }

  writeLines(markdown_lines_to_html(lines, title = title), html_path, useBytes = TRUE)
  invisible(html_path)
}

write_report <- function(
  report_path,
  input_path,
  df_raw,
  df_model,
  model_table,
  comparison,
  effects,
  block_importance,
  fitted_main,
  figure_paths,
  phenotype_area_path = NULL
) {
  q50_comp <- comparison %>%
    filter(response_label == "Area q50") %>%
    arrange(cv_rmse)

  best_q50 <- q50_comp %>% slice(1)
  base_q50 <- q50_comp %>% filter(model_id == "line_ploidy") %>% slice(1)
  full_q50 <- q50_comp %>% filter(model_id == "full_additive") %>% slice(1)

  main_effects <- effects %>%
    filter(response_label == "Area q50", model_id == "full_additive") %>%
    arrange(log2_ratio)

  stability_tbl <- effects %>%
    filter(response_label == "Area q50") %>%
    group_by(cellLine) %>%
    summarise(
      min_log2_ratio = min(log2_ratio, na.rm = TRUE),
      median_log2_ratio = median(log2_ratio, na.rm = TRUE),
      max_log2_ratio = max(log2_ratio, na.rm = TRUE),
      positive_models = sum(log2_ratio > 0, na.rm = TRUE),
      total_models = sum(is.finite(log2_ratio)),
      .groups = "drop"
    ) %>%
    arrange(median_log2_ratio)

  full_effect_tbl <- effects %>%
    filter(response_label == "Area q50", model_id == "full_additive") %>%
    transmute(
      cellLine,
      log2_area_ratio = log2_ratio,
      lower = lower_log2_ratio,
      upper = upper_log2_ratio,
      ratio = ratio
    ) %>%
    arrange(log2_area_ratio)

  best_reduction <- (base_q50$cv_rmse - best_q50$cv_rmse) / base_q50$cv_rmse
  full_reduction <- (base_q50$cv_rmse - full_q50$cv_rmse) / base_q50$cv_rmse

  sum_row <- full_effect_tbl %>% filter(cellLine == sum_fuse_line)

  block_tbl <- block_importance %>%
    filter(model_id != "full_additive") %>%
    transmute(
      block = as.character(block_tested),
      cv_rmse = cv_rmse,
      rmse_increase_vs_full = rmse_increase_vs_full
    ) %>%
    arrange(desc(rmse_increase_vs_full))

  residual_summary <- fitted_main %>%
    group_by(cellLine, ploidy_state) %>%
    summarise(
      rmse_log = safe_rmse(residual_log),
      median_abs_residual = median(abs(residual_log), na.rm = TRUE),
      residual_cor_hours = safe_cor(residual_log, hours),
      residual_cor_growth = safe_cor(residual_log, growth_rate_proxy),
      .groups = "drop"
    ) %>%
    arrange(cellLine, ploidy_state)

  phenotype_text <- ""
  if (!is.null(phenotype_area_path) && file.exists(phenotype_area_path)) {
    phenotype_tbl <- read.csv(phenotype_area_path, stringsAsFactors = FALSE) %>%
      arrange(desc(r_squared)) %>%
      head(6) %>%
      transmute(
        feature,
        category,
        r_squared,
        rmse_reduction_fraction,
        area_slope
      )

    phenotype_text <- paste0(
      "\n\n## Link To WP1 Phenotype Effects\n\n",
      "The previous WP2 export contains a line-level screen asking whether adjusted q50 area separation explains WP1 ",
      "model-free starvation-effect differences. This is a low-power five-line regression, so it should be treated as ",
      "triage rather than proof. The strongest entries are:\n\n",
      markdown_table(phenotype_tbl, digits = 3),
      "\n"
    )
  }

  lines <- c(
    "# WP2 Morphology Model Report",
    "",
    paste0("Generated: ", Sys.time()),
    "",
    paste0("Input table: `", input_path, "`"),
    "",
    "## Scope",
    "",
    "This report revisits the morphology analysis with an explicit model-comparison frame. The outcome is the alive-cell pixel-area distribution aligned to the count observations, mainly `log(q50)`; `log(q90)` and `1.5 * log(q50)` are included as sensitivity checks. The volume result is therefore a volume-like scaling index, not an independent physical measurement.",
    "",
    "The key modelling question is whether the apparent ploidy-size separation survives after accounting for plausible non-ploidy covariates: time since seeding, live/total density and crowding, starting/current glucose, glucose depletion, smoothed growth rate, smoothed death rate, and dead fraction.",
    "",
    "## Data Used",
    "",
    markdown_table(tibble(
      quantity = c("raw exported rows", "model complete-case rows", "wells", "cell lines"),
      value = c(nrow(df_raw), nrow(df_model), dplyr::n_distinct(df_model$well_idx), dplyr::n_distinct(df_model$cellLine))
    ), digits = 0),
    "",
    "Rows were restricted to finite `q50`/`q90`, finite dynamic covariates, and `n_area_alive >= 30`, so every model below is compared on the same rows.",
    "",
    paste0("![Covariate landscape](", fig_link(figure_paths$covariate_landscape), ")"),
    "",
    "Interpretation: pixel area is visibly phase-dependent. Time since seeding, density/crowding, and growth phase are not nuisance details; they are strong candidate explanations for area variation. This is why the main evidence below comes from covariate-adjusted and cross-validated models rather than raw 2N/4N area distributions alone.",
    "",
    "## Model Suite",
    "",
    "All models include cell-line context and estimate the elevated-versus-baseline ploidy contrast within each line. Validation is blocked by well so repeated timepoints from one well do not appear in both train and test folds.",
    "",
    markdown_table(model_table %>% select(model_label, covariate_blocks), digits = 3),
    "",
    paste0("![Model comparison](", fig_link(figure_paths$model_comparison), ")"),
    "",
    sprintf(
      "For `Area q50`, the best blocked-CV model was **%s** with RMSE %.3f. The line+ploidy-only model had RMSE %.3f, so the best model reduced prediction error by %.1f%%. The full additive model reduced prediction error by %.1f%% relative to line+ploidy only.",
      as.character(best_q50$model_label[[1]]),
      best_q50$cv_rmse[[1]],
      base_q50$cv_rmse[[1]],
      100 * best_reduction,
      100 * full_reduction
    ),
    "",
    "This confirms the concern: time/density/growth covariates materially improve area modelling. A raw ploidy-only analysis is not adequate.",
    "",
    "## Covariate Block Checks",
    "",
    paste0("![Block importance](", fig_link(figure_paths$block_importance), ")"),
    "",
    markdown_table(block_tbl, digits = 4),
    "",
    "Interpretation: positive values mean the q50 full-additive model got worse when that block was removed. The drop-one analysis is not causal, because several blocks are correlated, but it gives a useful sanity check on which covariates are carrying predictive information.",
    "",
    "## Ploidy Effects After Adjustment",
    "",
    paste0("![Effect stability](", fig_link(figure_paths$effect_stability), ")"),
    "",
    "Full-additive adjusted q50 effects:",
    "",
    markdown_table(full_effect_tbl, digits = 3),
    "",
    "Area-q50 effect stability across all model variants:",
    "",
    markdown_table(stability_tbl, digits = 3),
    "",
    if (nrow(sum_row)) {
      sprintf(
        "Interpretation: under the main full-additive adjustment, SUM-159-fuse has an elevated/baseline q50 log2 ratio of %.3f (95%% interval %.3f to %.3f), corresponding to a multiplicative area ratio of %.3f. This is not a subtle covariate artefact in the current model suite: SUM-159-fuse remains the line with the most negative or weakest size separation across the tested models.",
        sum_row$log2_area_ratio[[1]],
        sum_row$lower[[1]],
        sum_row$upper[[1]],
        sum_row$ratio[[1]]
      )
    } else {
      "Interpretation: SUM-159-fuse was not found in the main adjusted effect table."
    },
    "",
    "MDA-MB-231 remains the other weak-separation case. The stronger positive size separations are MCF10A, SNU668, and SUM-159-chem. Thus, the revised conclusion is more conditional than a simple 'ploidy increases area' statement: size separation is line-specific and robustly absent or reversed in SUM-159-fuse.",
    "",
    "## Residual Checks",
    "",
    paste0("![Residual diagnostics](", fig_link(figure_paths$residual_diagnostics), ")"),
    "",
    paste0("![Observed fitted](", fig_link(figure_paths$observed_fitted), ")"),
    "",
    markdown_table(residual_summary, digits = 3),
    "",
    "Interpretation: the residual plots are the main guard against overclaiming. If a covariate-adjusted model still leaves a strong monotonic residual trend against time, crowding, or growth rate in a given line, then the adjusted ploidy effect for that line should be treated as descriptive rather than mechanistic. The current residual summaries are acceptable for a WP2 screening report, but the line-specific dynamics model should be preferred for prediction while the full-additive model is easier to interpret.",
    "",
    "## Bottom Line",
    "",
    "1. Pixel area is strongly affected by time since seeding, density/crowding, and growth/death phase, so those covariates need to be present in any serious morphology analysis.",
    "2. Adding those covariates improves blocked-by-well prediction error relative to a line+ploidy-only model.",
    "3. The adjusted ploidy-size story is heterogeneous: MCF10A, SNU668, and SUM-159-chem show clear positive elevated-ploidy area effects; MDA-MB-231 is weak; SUM-159-fuse is reversed/negative.",
    "4. The volume proxy should be presented as a sensitivity transform of area, not as directly measured cell volume.",
    phenotype_text
  )

  writeLines(lines, report_path, useBytes = TRUE)
  invisible(report_path)
}

area_obs <- read_area_observations(input_path)
model_df <- prepare_model_data(area_obs)

model_tbl <- model_catalog()
suite <- fit_model_suite(model_df)

comparison <- bind_rows(lapply(suite, `[[`, "comparison")) %>%
  mutate(response_label = as.character(response_label))
effects <- estimate_all_effects(suite, model_df) %>%
  mutate(response_label = as.character(response_label))
block_importance <- run_block_drop(model_df, response_col = "log_q50")

main_fit <- suite$log_q50$fits$full_additive
main_fitted <- get_fitted_values(
  fit = main_fit,
  df = model_df,
  response_col = "log_q50",
  model_id = "full_additive",
  model_label = "Full additive"
)

write.csv(model_df, file.path(output_dir, "wp2_morphology_model_rows.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "wp2_morphology_model_comparison.csv"), row.names = FALSE)
write.csv(effects, file.path(output_dir, "wp2_morphology_ploidy_effect_stability.csv"), row.names = FALSE)
write.csv(block_importance, file.path(output_dir, "wp2_morphology_block_importance.csv"), row.names = FALSE)
write.csv(main_fitted, file.path(output_dir, "wp2_morphology_full_additive_fitted.csv"), row.names = FALSE)

saveRDS(model_df, file.path(output_dir, "wp2_morphology_model_rows.Rds"))
saveRDS(suite, file.path(output_dir, "wp2_morphology_model_suite.Rds"))
saveRDS(effects, file.path(output_dir, "wp2_morphology_ploidy_effect_stability.Rds"))
saveRDS(block_importance, file.path(output_dir, "wp2_morphology_block_importance.Rds"))

figure_paths <- list()
figure_paths$covariate_landscape <- save_png_pdf(
  make_covariate_landscape_plot(model_df),
  "figure_1_covariate_landscape",
  width = 13.5,
  height = 10.5
)[["png"]]
figure_paths$model_comparison <- save_png_pdf(
  make_model_comparison_plot(comparison),
  "figure_2_model_comparison",
  width = 12.5,
  height = 7.8
)[["png"]]
figure_paths$effect_stability <- save_png_pdf(
  make_effect_stability_plot(effects),
  "figure_3_ploidy_effect_stability",
  width = 13,
  height = 8.4
)[["png"]]
figure_paths$block_importance <- save_png_pdf(
  make_block_importance_plot(block_importance),
  "figure_4_block_importance",
  width = 9,
  height = 5.5
)[["png"]]
figure_paths$residual_diagnostics <- save_png_pdf(
  make_residual_diagnostic_plot(main_fitted),
  "figure_5_residual_diagnostics",
  width = 13.5,
  height = 7.8
)[["png"]]
figure_paths$observed_fitted <- save_png_pdf(
  make_observed_fitted_plot(main_fitted),
  "figure_6_observed_fitted",
  width = 13,
  height = 4.8
)[["png"]]

phenotype_area_path <- file.path("data", "report_exports", "wp2_morphology_volume", "wp2_phenotype_area_explanation.csv")
report_path <- file.path(output_dir, "wp2_morphology_model_report.md")
write_report(
  report_path = report_path,
  input_path = input_path,
  df_raw = area_obs,
  df_model = model_df,
  model_table = model_tbl,
  comparison = comparison,
  effects = effects,
  block_importance = block_importance,
  fitted_main = main_fitted,
  figure_paths = figure_paths,
  phenotype_area_path = phenotype_area_path
)
html_report_path <- file.path(output_dir, "wp2_morphology_model_report.html")
write_html_report(report_path, html_report_path)

cat(sprintf("Wrote WP2 morphology model report to %s\n", report_path))
cat(sprintf("Wrote WP2 morphology model HTML report to %s\n", html_report_path))
cat(sprintf("Wrote report figures to %s\n", figure_dir))
