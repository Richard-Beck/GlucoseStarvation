args <- commandArgs(trailingOnly = TRUE)

run_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "image_processing_runs", "wp3_nuclear_size_pilot")
}

output_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("data", "report_exports", "wp3_nuclear_size")
}

figure_dir <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("figures", "wp3_nuclear_size")
}

alive_feature_csv <- if (length(args) >= 4 && nzchar(args[4])) {
  args[4]
} else {
  file.path("data", "image_processing_runs", "run_20260324_233122", "object_features.csv")
}

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
})

wp3_theme <- function(base_size = 8) {
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

resolve_run_dir <- function(path) {
  if (
    file.exists(file.path(path, "wp3_nuclear_object_features.csv")) &&
      file.exists(file.path(path, "summary.json"))
  ) {
    return(path)
  }

  candidates <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  completed <- candidates[
    file.exists(file.path(candidates, "wp3_nuclear_object_features.csv")) &
      file.exists(file.path(candidates, "summary.json"))
  ]
  if (length(completed)) {
    return(completed[order(file.info(completed)$mtime, decreasing = TRUE)][[1]])
  }

  candidates <- candidates[file.exists(file.path(candidates, "wp3_nuclear_object_features.csv"))]
  if (!length(candidates)) {
    stop("Could not find wp3_nuclear_object_features.csv in run_dir or its immediate children: ", path, call. = FALSE)
  }
  warning("No completed WP3 run with summary.json found; using the most recently modified partial run.", call. = FALSE)
  candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][[1]]
}

ploidy_state_from_label <- function(x) {
  dplyr::case_when(
    x %in% c("2N", "parental", "low") ~ "baseline",
    x %in% c("3N", "4N", "high") ~ "elevated",
    TRUE ~ as.character(x)
  )
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  sqrt(mean(x^2))
}

safe_mae <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(abs(x))
}

safe_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 3L) return(NA_real_)
  denom <- sum((obs[ok] - mean(obs[ok]))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  1 - sum((obs[ok] - pred[ok])^2) / denom
}

scale_vector <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) sig <- 1
  (x - mu) / sig
}

load_alive_scores <- function(path, image_keys) {
  if (!file.exists(path)) {
    warning("Alive-score source not found: ", path, ". Falling back to normalized live-minus-dead channel means.", call. = FALSE)
    return(NULL)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    warning("Package data.table is unavailable. Falling back to normalized live-minus-dead channel means.", call. = FALSE)
    return(NULL)
  }

  message("Reading alive/dead z-scores from ", path)
  dt <- data.table::fread(
    path,
    select = c("image_key", "object_id", "live_bg_z", "dead_bg_z"),
    showProgress = FALSE
  )
  dt <- dt[image_key %in% image_keys]
  dt[, object_id := as.integer(object_id)]
  dt[, live_bg_z := as.numeric(live_bg_z)]
  dt[, dead_bg_z := as.numeric(dead_bg_z)]
  as_tibble(dt)
}

build_model_rows <- function(objects_for_measure) {
  image_covariates <- objects_for_measure %>%
    group_by(image_key) %>%
    summarise(
      segmented_cells_image = n(),
      alive_cells_image = sum(is_alive, na.rm = TRUE),
      alive_area_image = sum(cell_area_px[is_alive], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      objects_for_measure %>%
        distinct(image_key, cellLine) %>%
        mutate(cellLine = as.character(cellLine)),
      by = "image_key"
    ) %>%
    group_by(cellLine) %>%
    mutate(
      confluency_rank = percent_rank(alive_area_image),
      confluency_bin = factor(
        c("low", "mid", "high")[pmax(1L, pmin(3L, dplyr::ntile(alive_area_image, 3L)))],
        levels = c("low", "mid", "high")
      )
    ) %>%
    ungroup() %>%
    select(-cellLine)

  out <- objects_for_measure %>%
    filter(is_alive, has_nucleus, is.finite(nuclear_area_px), nuclear_area_px > 0, is.finite(cell_area_px), cell_area_px > 0) %>%
    left_join(image_covariates, by = "image_key") %>%
    mutate(
      cellLine = factor(cellLine, levels = sort(unique(cellLine))),
      ploidy_state = factor(ploidy_state, levels = c("baseline", "elevated")),
      ploidy_display = factor(ploidy_display),
      glucose_bin = factor(glucose_bin, levels = c("low", "mid", "high")),
      log_nuclear_area = log(nuclear_area_px),
      log_cell_area = log(cell_area_px),
      log_alive_cells_image = log1p(alive_cells_image),
      log_segmented_cells_image = log1p(segmented_cells_image),
      log_alive_area_image = log1p(alive_area_image),
      confluency_bin = factor(confluency_bin, levels = c("low", "mid", "high")),
      G0_log = log1p(pmax(G0, 0)),
      hours_z = scale_vector(hours),
      hours_log_z = scale_vector(log1p(hours)),
      log_cell_area_z = scale_vector(log_cell_area),
      log_alive_cells_image_z = scale_vector(log_alive_cells_image),
      log_segmented_cells_image_z = scale_vector(log_segmented_cells_image),
      log_alive_area_image_z = scale_vector(log_alive_area_image),
      G0_log_z = scale_vector(G0_log)
    )

  fold_tbl <- out %>%
    distinct(cellLine, ploidy_state, image_key) %>%
    arrange(cellLine, ploidy_state, image_key) %>%
    group_by(cellLine, ploidy_state) %>%
    mutate(fold = ((row_number() - 1L) %% 5L) + 1L) %>%
    ungroup() %>%
    select(image_key, fold)

  out %>%
    left_join(fold_tbl, by = "image_key")
}

model_catalog <- function() {
  time_terms <- "hours_z + I(hours_z^2)"
  glucose_terms <- "G0_log_z + glucose_bin"
  density_terms <- "log_alive_cells_image_z + log_segmented_cells_image_z + log_alive_area_image_z"
  cell_area_terms <- "log_cell_area_z"
  full_terms <- paste(cell_area_terms, time_terms, glucose_terms, density_terms, sep = " + ")

  tibble(
    model_id = c(
      "line_ploidy",
      "cell_area",
      "time",
      "glucose",
      "image_density",
      "cell_area_time_glucose",
      "full_additive",
      "ploidy_cell_area_interaction"
    ),
    model_label = c(
      "Line + ploidy",
      "+ cell area",
      "+ time",
      "+ glucose",
      "+ image density",
      "+ cell area + time + glucose",
      "Full additive",
      "Full + ploidy x cell area"
    ),
    rhs = c(
      "cellLine * ploidy_state",
      paste("cellLine * ploidy_state", cell_area_terms, sep = " + "),
      paste("cellLine * ploidy_state", time_terms, sep = " + "),
      paste("cellLine * ploidy_state", glucose_terms, sep = " + "),
      paste("cellLine * ploidy_state", density_terms, sep = " + "),
      paste("cellLine * ploidy_state", cell_area_terms, time_terms, glucose_terms, sep = " + "),
      paste("cellLine * ploidy_state", full_terms, sep = " + "),
      paste("cellLine * ploidy_state", full_terms, "ploidy_state:log_cell_area_z", sep = " + ")
    ),
    covariate_blocks = c(
      "line, ploidy",
      "line, ploidy, cell area",
      "line, ploidy, time",
      "line, ploidy, glucose",
      "line, ploidy, image density",
      "line, ploidy, cell area, time, glucose",
      "line, ploidy, cell area, time, glucose, image density",
      "full additive plus ploidy-specific cell-area scaling"
    )
  )
}

fit_one_model <- function(df, response_col, rhs) {
  stats::lm(stats::as.formula(paste(response_col, "~", rhs)), data = df)
}

run_blocked_cv <- function(df, response_col, model_tbl) {
  folds <- sort(unique(df$fold))
  rows <- vector("list", length = nrow(model_tbl) * length(folds))
  idx <- 1L

  for (i in seq_len(nrow(model_tbl))) {
    for (fold in folds) {
      train <- df %>% filter(fold != !!fold)
      test <- df %>% filter(fold == !!fold)
      fit <- try(fit_one_model(train, response_col, model_tbl$rhs[[i]]), silent = TRUE)
      if (inherits(fit, "try-error")) {
        rows[[idx]] <- tibble(
          model_id = model_tbl$model_id[[i]],
          fold = fold,
          n_test = nrow(test),
          rmse = NA_real_,
          mae = NA_real_,
          r2 = NA_real_
        )
      } else {
        pred <- as.numeric(stats::predict(fit, newdata = test))
        obs <- test[[response_col]]
        rows[[idx]] <- tibble(
          model_id = model_tbl$model_id[[i]],
          fold = fold,
          n_test = nrow(test),
          rmse = safe_rmse(obs - pred),
          mae = safe_mae(obs - pred),
          r2 = safe_r2(obs, pred)
        )
      }
      idx <- idx + 1L
    }
  }

  bind_rows(rows) %>%
    group_by(model_id) %>%
    summarise(
      cv_rmse = weighted.mean(rmse, w = n_test, na.rm = TRUE),
      cv_mae = weighted.mean(mae, w = n_test, na.rm = TRUE),
      cv_r2 = weighted.mean(r2, w = n_test, na.rm = TRUE),
      n_test = sum(n_test, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(model_tbl, by = "model_id")
}

run_model_suite <- function(df, response_col = "log_nuclear_area") {
  models <- model_catalog()
  fits <- lapply(models$rhs, function(rhs) fit_one_model(df, response_col, rhs))
  names(fits) <- models$model_id

  fit_metrics <- bind_rows(lapply(seq_len(nrow(models)), function(i) {
    fit <- fits[[models$model_id[[i]]]]
    pred <- as.numeric(stats::fitted(fit))
    obs <- df[[response_col]]
    tibble(
      model_id = models$model_id[[i]],
      train_rmse = safe_rmse(obs - pred),
      train_mae = safe_mae(obs - pred),
      train_r2 = safe_r2(obs, pred),
      aic = stats::AIC(fit),
      n = stats::nobs(fit),
      df_model = length(stats::coef(fit))
    )
  }))

  cv_metrics <- run_blocked_cv(df, response_col, models)
  comparison <- cv_metrics %>%
    left_join(fit_metrics, by = "model_id") %>%
    mutate(
      response = response_col,
      delta_cv_rmse_vs_line_ploidy = cv_rmse - cv_rmse[model_id == "line_ploidy"][[1]],
      pct_cv_rmse_reduction_vs_line_ploidy = 100 * (cv_rmse[model_id == "line_ploidy"][[1]] - cv_rmse) / cv_rmse[model_id == "line_ploidy"][[1]]
    ) %>%
    arrange(cv_rmse)

  list(models = models, fits = fits, comparison = comparison)
}

block_drop_catalog <- function() {
  full_rhs <- model_catalog() %>% filter(model_id == "full_additive") %>% pull(rhs)
  tibble(
    model_id = c("full_additive", "drop_cell_area", "drop_time", "drop_glucose", "drop_image_density"),
    block_tested = c("none", "cell area", "time", "glucose", "image density"),
    rhs = c(
      full_rhs,
      "cellLine * ploidy_state + hours_z + I(hours_z^2) + G0_log_z + glucose_bin + log_alive_cells_image_z + log_segmented_cells_image_z + log_alive_area_image_z",
      "cellLine * ploidy_state + log_cell_area_z + G0_log_z + glucose_bin + log_alive_cells_image_z + log_segmented_cells_image_z + log_alive_area_image_z",
      "cellLine * ploidy_state + log_cell_area_z + hours_z + I(hours_z^2) + log_alive_cells_image_z + log_segmented_cells_image_z + log_alive_area_image_z",
      "cellLine * ploidy_state + log_cell_area_z + hours_z + I(hours_z^2) + G0_log_z + glucose_bin"
    )
  )
}

run_block_drop <- function(df, response_col = "log_nuclear_area") {
  block_tbl <- block_drop_catalog()
  cv <- run_blocked_cv(
    df = df,
    response_col = response_col,
    model_tbl = block_tbl %>% transmute(model_id, rhs)
  )
  full_rmse <- cv$cv_rmse[cv$model_id == "full_additive"][[1]]
  block_tbl %>%
    select(model_id, block_tested) %>%
    left_join(cv, by = "model_id") %>%
    mutate(
      rmse_increase_vs_full = cv_rmse - full_rmse,
      direction = if_else(rmse_increase_vs_full >= 0, "worse without block", "better without block"),
      block_tested = factor(block_tested, levels = c("none", "cell area", "time", "glucose", "image density"))
    )
}

make_model_comparison_plot <- function(comparison) {
  comparison %>%
    mutate(model_label = factor(model_label, levels = rev(model_label[order(cv_rmse)]))) %>%
    ggplot(aes(cv_rmse, model_label)) +
    geom_col(fill = "#386CB0", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", cv_rmse)), hjust = -0.08, size = 2.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "WP3 alive-cell nuclear-area model comparison",
      subtitle = "Five-fold cross-validation blocked by image.",
      x = "CV RMSE on log(nuclear area)",
      y = NULL
    ) +
    wp3_theme()
}

make_block_importance_plot <- function(block_importance) {
  block_importance %>%
    filter(model_id != "full_additive") %>%
    mutate(block_tested = factor(block_tested, levels = rev(levels(block_tested)))) %>%
    ggplot(aes(rmse_increase_vs_full, block_tested, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    scale_fill_manual(values = c("worse without block" = "#D95F02", "better without block" = "#1B9E77")) +
    labs(
      title = "WP3 nuclear-area covariate block importance",
      subtitle = "Positive values mean the full model predicts better when the block is retained.",
      x = "CV RMSE increase versus full additive model",
      y = NULL,
      fill = NULL
    ) +
    wp3_theme()
}

run_dir <- resolve_run_dir(run_dir)
object_csv <- file.path(run_dir, "wp3_nuclear_object_features.csv")
image_csv <- file.path(run_dir, "wp3_nuclear_image_qc.csv")

objects <- read.csv(object_csv, stringsAsFactors = FALSE, check.names = FALSE)
images <- read.csv(image_csv, stringsAsFactors = FALSE, check.names = FALSE)

objects <- objects %>%
  mutate(
    object_id = as.integer(object_id),
    G0 = as.numeric(G0),
    hours = as.numeric(hours),
    cell_area_px = as.numeric(cell_area_px),
    nuclear_area_px = as.numeric(nuclear_area_px),
    largest_nuclear_area_px = as.numeric(largest_nuclear_area_px),
    nuclear_to_cell_area_ratio = as.numeric(nuclear_to_cell_area_ratio),
    nuclear_mean_intensity = as.numeric(nuclear_mean_intensity),
    nuclear_integrated_intensity = as.numeric(nuclear_integrated_intensity),
    cell_mean_nuclear_channel = as.numeric(cell_mean_nuclear_channel),
    cell_mean_alive_channel = as.numeric(cell_mean_alive_channel),
    cell_mean_dead_channel = as.numeric(cell_mean_dead_channel),
    cell_area_pass = as.integer(cell_area_pass),
    has_nucleus = nuclear_area_px > 0,
    ploidy_state = ploidy_state_from_label(ploidy),
    ploidy_display = paste(ploidy_state, ploidy, sep = ": ")
  )

alive_scores <- load_alive_scores(alive_feature_csv, unique(objects$image_key))
if (!is.null(alive_scores)) {
  objects <- objects %>%
    left_join(alive_scores, by = c("image_key", "object_id")) %>%
    mutate(
      alive_score = live_bg_z - dead_bg_z,
      is_alive = is.finite(alive_score) & alive_score >= 0,
      alive_gate_source = "live_bg_z - dead_bg_z >= 0 from object_features.csv"
    )
} else {
  objects <- objects %>%
    mutate(
      live_bg_z = NA_real_,
      dead_bg_z = NA_real_,
      alive_score = cell_mean_alive_channel - cell_mean_dead_channel,
      is_alive = is.finite(alive_score) & alive_score >= 0,
      alive_gate_source = "cell_mean_alive_channel - cell_mean_dead_channel >= 0 fallback"
    )
}

images <- images %>%
  mutate(
    G0 = as.numeric(G0),
    hours = as.numeric(hours),
    n_cells = as.integer(n_cells),
    n_cells_with_nucleus = as.integer(n_cells_with_nucleus),
    nuclear_success_fraction = as.numeric(nuclear_success_fraction),
    ploidy_state = ploidy_state_from_label(ploidy),
    ploidy_display = paste(ploidy_state, ploidy, sep = ": ")
  )

objects_for_measure <- objects %>%
  filter(cell_area_pass == 1L)

condition_summary <- objects_for_measure %>%
  group_by(cellLine, ploidy, ploidy_state, glucose_bin, time_bin, G0) %>%
  summarise(
    n_segmented_cells = n(),
    n_alive_cells = sum(is_alive, na.rm = TRUE),
    n_alive_cells_with_nucleus = sum(is_alive & has_nucleus, na.rm = TRUE),
    alive_nuclear_success_fraction = n_alive_cells_with_nucleus / pmax(n_alive_cells, 1),
    median_cell_area_px = safe_median(cell_area_px[is_alive]),
    median_nuclear_area_px = safe_median(nuclear_area_px[is_alive & has_nucleus]),
    median_nuclear_to_cell_area_ratio = safe_median(nuclear_to_cell_area_ratio[is_alive & has_nucleus]),
    median_nuclear_mean_intensity = safe_median(nuclear_mean_intensity[is_alive & has_nucleus]),
    .groups = "drop"
  )

image_qc_summary <- objects_for_measure %>%
  group_by(image_key, cellLine, ploidy, ploidy_state, glucose_bin, time_bin, G0, hours) %>%
  summarise(
    n_segmented_cells = n(),
    n_alive_cells = sum(is_alive, na.rm = TRUE),
    n_alive_cells_with_nucleus = sum(is_alive & has_nucleus, na.rm = TRUE),
    alive_nuclear_success_fraction = n_alive_cells_with_nucleus / pmax(n_alive_cells, 1),
    .groups = "drop"
  )

ploidy_effects <- condition_summary %>%
  filter(ploidy_state %in% c("baseline", "elevated")) %>%
  group_by(cellLine, glucose_bin, time_bin) %>%
  filter(n_distinct(ploidy_state) == 2L) %>%
  ungroup() %>%
  select(cellLine, glucose_bin, time_bin, ploidy_state, median_cell_area_px, median_nuclear_area_px, median_nuclear_to_cell_area_ratio) %>%
  pivot_wider(
    names_from = ploidy_state,
    values_from = c(median_cell_area_px, median_nuclear_area_px, median_nuclear_to_cell_area_ratio)
  ) %>%
  mutate(
    log2_nuclear_area_ratio_elevated_vs_baseline = log2(median_nuclear_area_px_elevated / median_nuclear_area_px_baseline),
    log2_cell_area_ratio_elevated_vs_baseline = log2(median_cell_area_px_elevated / median_cell_area_px_baseline),
    delta_nuclear_to_cell_area_ratio = median_nuclear_to_cell_area_ratio_elevated - median_nuclear_to_cell_area_ratio_baseline
  )

high_low_effects <- condition_summary %>%
  filter(ploidy_state %in% c("baseline", "elevated")) %>%
  select(cellLine, glucose_bin, time_bin, ploidy_state, ploidy, n_alive_cells_with_nucleus, median_cell_area_px, median_nuclear_area_px, median_nuclear_to_cell_area_ratio) %>%
  pivot_wider(
    names_from = ploidy_state,
    values_from = c(ploidy, n_alive_cells_with_nucleus, median_cell_area_px, median_nuclear_area_px, median_nuclear_to_cell_area_ratio)
  ) %>%
  filter(
    is.finite(median_nuclear_area_px_baseline),
    is.finite(median_nuclear_area_px_elevated),
    n_alive_cells_with_nucleus_baseline >= 10,
    n_alive_cells_with_nucleus_elevated >= 10
  ) %>%
  mutate(
    ploidy_comparison = paste0(ploidy_elevated, " / ", ploidy_baseline),
    nuclear_area_ratio_elevated_vs_baseline = median_nuclear_area_px_elevated / median_nuclear_area_px_baseline,
    log2_nuclear_area_ratio_elevated_vs_baseline = log2(nuclear_area_ratio_elevated_vs_baseline),
    cell_area_ratio_elevated_vs_baseline = median_cell_area_px_elevated / median_cell_area_px_baseline,
    nuclear_to_cell_ratio_delta_elevated_minus_baseline = median_nuclear_to_cell_area_ratio_elevated - median_nuclear_to_cell_area_ratio_baseline
  )

model_rows <- build_model_rows(objects_for_measure)
model_suite <- run_model_suite(model_rows, response_col = "log_nuclear_area")
model_comparison <- model_suite$comparison
block_importance <- run_block_drop(model_rows, response_col = "log_nuclear_area")
full_fit <- model_suite$fits$full_additive
full_coefficients <- tibble(
  term = names(stats::coef(full_fit)),
  estimate = as.numeric(stats::coef(full_fit))
)

write.csv(condition_summary, file.path(output_dir, "wp3_nuclear_condition_summary.csv"), row.names = FALSE)
write.csv(image_qc_summary, file.path(output_dir, "wp3_nuclear_image_qc_summary.csv"), row.names = FALSE)
write.csv(ploidy_effects, file.path(output_dir, "wp3_nuclear_ploidy_effects.csv"), row.names = FALSE)
write.csv(high_low_effects, file.path(output_dir, "wp3_nuclear_high_vs_low_effects.csv"), row.names = FALSE)
write.csv(model_rows, file.path(output_dir, "wp3_nuclear_model_rows_alive.csv"), row.names = FALSE)
write.csv(model_comparison, file.path(output_dir, "wp3_nuclear_model_comparison.csv"), row.names = FALSE)
write.csv(block_importance, file.path(output_dir, "wp3_nuclear_model_block_importance.csv"), row.names = FALSE)
write.csv(full_coefficients, file.path(output_dir, "wp3_nuclear_full_additive_coefficients.csv"), row.names = FALSE)

saveRDS(condition_summary, file.path(output_dir, "wp3_nuclear_condition_summary.Rds"))
saveRDS(image_qc_summary, file.path(output_dir, "wp3_nuclear_image_qc_summary.Rds"))
saveRDS(ploidy_effects, file.path(output_dir, "wp3_nuclear_ploidy_effects.Rds"))
saveRDS(high_low_effects, file.path(output_dir, "wp3_nuclear_high_vs_low_effects.Rds"))
saveRDS(model_rows, file.path(output_dir, "wp3_nuclear_model_rows_alive.Rds"))
saveRDS(model_suite, file.path(output_dir, "wp3_nuclear_model_suite.Rds"))
saveRDS(block_importance, file.path(output_dir, "wp3_nuclear_model_block_importance.Rds"))

set.seed(1)
plot_objects <- objects_for_measure %>%
  filter(is_alive, has_nucleus, is.finite(cell_area_px), is.finite(nuclear_area_px), cell_area_px > 0, nuclear_area_px > 0)
if (nrow(plot_objects) > 12000) {
  plot_objects <- plot_objects[sample.int(nrow(plot_objects), 12000), , drop = FALSE]
}

fig5 <- ggplot(plot_objects, aes(x = cell_area_px, y = nuclear_area_px, color = ploidy_display)) +
  geom_point(alpha = 0.18, size = 0.35) +
  stat_ellipse(aes(group = ploidy_display), linewidth = 0.35, alpha = 0.65, show.legend = FALSE) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  facet_wrap(~ cellLine, scales = "free") +
  labs(
    title = "WP3 pilot: alive-cell area versus nuclear area",
    subtitle = "Alive gate: live_bg_z - dead_bg_z >= 0.",
    x = "Cell area (px, CPSAM mask)",
    y = "Nuclear area (px, thresholded stain within cell)",
    color = "Ploidy"
  ) +
  wp3_theme()
save_plot_pair(fig5, "wp3_fig5_cell_vs_nuclear_area", width = 9, height = 6.2)

ratio_plot <- ggplot(plot_objects, aes(x = ploidy_display, y = nuclear_to_cell_area_ratio, fill = ploidy_display)) +
  geom_violin(scale = "width", linewidth = 0.25, alpha = 0.75, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.8) +
  facet_wrap(~ cellLine, scales = "free_x") +
  coord_cartesian(ylim = quantile(plot_objects$nuclear_to_cell_area_ratio, c(0.01, 0.99), na.rm = TRUE)) +
  labs(
    title = "WP3 pilot: alive-cell nuclear-to-cell area ratio",
    x = NULL,
    y = "Nuclear:cell area ratio",
    fill = "Ploidy"
  ) +
  wp3_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot_pair(ratio_plot, "wp3_nuclear_to_cell_ratio_by_line", width = 9, height = 5.8)

success_plot <- ggplot(image_qc_summary, aes(x = glucose_bin, y = time_bin, fill = alive_nuclear_success_fraction)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = percent(alive_nuclear_success_fraction, accuracy = 1)), size = 2.3) +
  scale_fill_viridis_c(labels = percent, limits = c(0, 1), option = "C") +
  facet_grid(cellLine ~ ploidy) +
  labs(
    title = "WP3 pilot: alive-cell nuclear segmentation success",
    x = "Glucose bin",
    y = "Time bin",
    fill = "Alive cells\nwith nucleus"
  ) +
  wp3_theme(base_size = 7)
save_plot_pair(success_plot, "wp3_nuclear_segmentation_success", width = 9, height = 8)

model_plot_rows <- model_rows
if (nrow(model_plot_rows) > 24000) {
  set.seed(2)
  model_plot_rows <- model_plot_rows[sample.int(nrow(model_plot_rows), 24000), , drop = FALSE]
}

time_violin_plot <- ggplot(
  model_plot_rows,
  aes(x = ploidy_state, y = nuclear_area_px, fill = ploidy_state)
) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.2, alpha = 0.75) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_manual(values = c("baseline" = "#2C7BB6", "elevated" = "#D7191C")) +
  facet_grid(cellLine ~ time_bin, scales = "free_y") +
  labs(
    title = "WP3 pilot: alive-cell nuclear area by time bin",
    subtitle = "Baseline/elevated mapping is line-specific; y-axis is log scaled.",
    x = "Ploidy state",
    y = "Nuclear area (px)",
    fill = "Ploidy state"
  ) +
  wp3_theme(base_size = 7)
save_plot_pair(time_violin_plot, "wp3_high_low_nuclear_area_by_time_violin", width = 9.5, height = 8.5)

confluency_violin_plot <- ggplot(
  model_plot_rows,
  aes(x = ploidy_state, y = nuclear_area_px, fill = ploidy_state)
) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.2, alpha = 0.75) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_manual(values = c("baseline" = "#2C7BB6", "elevated" = "#D7191C")) +
  facet_grid(cellLine ~ confluency_bin, scales = "free_y") +
  labs(
    title = "WP3 pilot: alive-cell nuclear area by confluency bin",
    subtitle = "Confluency is approximated by image-level summed alive-cell area, binned within cell line.",
    x = "Ploidy state",
    y = "Nuclear area (px)",
    fill = "Ploidy state"
  ) +
  wp3_theme(base_size = 7)
save_plot_pair(confluency_violin_plot, "wp3_high_low_nuclear_area_by_confluency_violin", width = 9.5, height = 8.5)

save_plot_pair(
  make_model_comparison_plot(model_comparison),
  "wp3_nuclear_model_comparison",
  width = 8.5,
  height = 4.8,
  dpi = 360
)
save_plot_pair(
  make_block_importance_plot(block_importance),
  "wp3_nuclear_model_block_importance",
  width = 7.5,
  height = 3.8,
  dpi = 360
)

best_model <- model_comparison %>% arrange(cv_rmse) %>% slice(1)
line_model <- model_comparison %>% filter(model_id == "line_ploidy") %>% slice(1)
full_model <- model_comparison %>% filter(model_id == "full_additive") %>% slice(1)

summary_lines <- c(
  "WP3 nuclear-size pilot summary",
  "",
  paste("Run directory:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)),
  paste("Alive gate:", unique(objects$alive_gate_source)[[1]]),
  paste("Object rows:", nrow(objects)),
  paste("Objects passing cell-area gate:", nrow(objects_for_measure)),
  paste("Alive objects passing cell-area gate:", sum(objects_for_measure$is_alive, na.rm = TRUE)),
  paste("Alive objects passing gate with nucleus:", nrow(model_rows)),
  paste("Images:", nrow(images)),
  paste("Images with status ok:", sum(images$status == "ok", na.rm = TRUE)),
  paste("Alive-cell nuclear success fraction:", percent(nrow(model_rows) / max(sum(objects_for_measure$is_alive, na.rm = TRUE), 1), accuracy = 0.1)),
  paste("Best nuclear-area CV model:", best_model$model_label[[1]], sprintf("(RMSE %.3f)", best_model$cv_rmse[[1]])),
  paste("Line+ploidy CV RMSE:", sprintf("%.3f", line_model$cv_rmse[[1]])),
  paste("Full additive CV RMSE:", sprintf("%.3f", full_model$cv_rmse[[1]])),
  "",
  "Primary outputs:",
  "- wp3_fig5_cell_vs_nuclear_area.{png,pdf}",
  "- wp3_high_low_nuclear_area_by_time_violin.{png,pdf}",
  "- wp3_high_low_nuclear_area_by_confluency_violin.{png,pdf}",
  "- wp3_nuclear_to_cell_ratio_by_line.{png,pdf}",
  "- wp3_nuclear_segmentation_success.{png,pdf}",
  "- wp3_nuclear_model_comparison.{png,pdf}",
  "- wp3_nuclear_model_block_importance.{png,pdf}",
  "- wp3_nuclear_high_vs_low_effects.csv",
  "- wp3_nuclear_model_comparison.csv",
  "- wp3_nuclear_model_block_importance.csv"
)
writeLines(summary_lines, file.path(output_dir, "wp3_nuclear_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
