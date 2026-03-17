#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

SCREEN_ROOT <- if (length(args) >= 1) args[1] else "data/gpath_transfer_cv_single_param"
OUTPUT_DIR <- if (length(args) >= 2) args[2] else file.path(SCREEN_ROOT, "report")
STAN_DATA_PATH <- if (length(args) >= 3) args[3] else "data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"
MODEL_NAME <- if (length(args) >= 4) args[4] else "gpath"

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(posterior)

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/posterior_io_utils.R")
source("R/transfer_cv_plot_utils.R")

cat(sprintf(">>> report_gpath_single_param_screen.R cwd: %s\n", getwd()))
cat(sprintf(">>> Screen root: %s\n", SCREEN_ROOT))

safe_read_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

discover_mask_runs <- function(screen_root, model_name) {
  mask_dirs <- list.dirs(screen_root, recursive = FALSE, full.names = TRUE)
  rows <- list()
  idx <- 1L

  for (mask_root in mask_dirs) {
    model_dir <- file.path(mask_root, model_name)
    if (!dir.exists(model_dir)) {
      next
    }

    run_ids <- list.dirs(model_dir, recursive = FALSE, full.names = FALSE)
    for (run_id in run_ids) {
      base_dir <- file.path(model_dir, run_id)
      if (!file.exists(file.path(base_dir, "transfer_best_start_summary.Rds"))) {
        next
      }

      rows[[idx]] <- data.frame(
        mask_label = basename(mask_root),
        mask_root = mask_root,
        run_id = run_id,
        base_dir = base_dir,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  if (!length(rows)) {
    return(tibble(
      mask_label = character(),
      mask_root = character(),
      run_id = character(),
      base_dir = character()
    ))
  }

  bind_rows(rows)
}

collect_predictive_summaries <- function(run_df) {
  rows <- lapply(seq_len(nrow(run_df)), function(i) {
    x <- safe_read_rds(file.path(run_df$base_dir[i], "transfer_comparison_summary.Rds"))
    if (is.null(x) || !nrow(x)) {
      return(NULL)
    }
    x$mask_label <- run_df$mask_label[i]
    x$run_id <- run_df$run_id[i]
    x
  })
  bind_rows(rows)
}

collect_parameter_summaries <- function(run_df) {
  rows <- lapply(seq_len(nrow(run_df)), function(i) {
    x <- safe_read_rds(file.path(run_df$base_dir[i], "parameter_transfer_summary.Rds"))
    if (is.null(x) || !nrow(x)) {
      return(NULL)
    }
    x$mask_label <- run_df$mask_label[i]
    x$run_id <- run_df$run_id[i]
    x$mean_log_error_improvement <- x$mean_abs_log_null_holdout_ratio - x$mean_abs_log_holdout_ratio
    x
  })
  bind_rows(rows)
}

collect_parameter_states <- function(run_df) {
  rows <- lapply(seq_len(nrow(run_df)), function(i) {
    x <- safe_read_rds(file.path(run_df$base_dir[i], "parameter_transfer_states.Rds"))
    if (is.null(x) || !nrow(x)) {
      return(NULL)
    }
    x$mask_label <- run_df$mask_label[i]
    x$run_id <- run_df$run_id[i]
    x
  })
  bind_rows(rows)
}

derive_glucose_transfer_summaries <- function(state_df) {
  if (is.null(state_df) || !nrow(state_df)) {
    return(tibble())
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
    select(mask_label, run_id, line_id, direction, fit_type, state, parameter_key, value) %>%
    tidyr::pivot_wider(names_from = parameter_key, values_from = value)

  needed <- c("ae1", "ah1", "YR1", "m")
  if (!all(needed %in% names(wide))) {
    return(tibble())
  }

  wide %>%
    transmute(
      mask_label,
      run_id,
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

build_derived_transfer_comparison <- function(derived_df) {
  if (!nrow(derived_df)) {
    return(tibble())
  }

  transfer_obs <- derived_df %>%
    filter(fit_type == "transfer", state == "observed") %>%
    select(mask_label, run_id, line_id, direction, summary_name, observed_transfer = summary_value)

  transfer_holdout <- derived_df %>%
    filter(fit_type == "transfer", state == "holdout") %>%
    select(mask_label, run_id, line_id, direction, summary_name, transfer_holdout = summary_value)

  oracle_holdout <- derived_df %>%
    filter(fit_type == "oracle", state == "holdout") %>%
    select(mask_label, run_id, line_id, direction, summary_name, oracle_holdout = summary_value)

  oracle_obs <- derived_df %>%
    filter(fit_type == "oracle", state == "observed") %>%
    select(mask_label, run_id, line_id, direction, summary_name, observed_oracle = summary_value)

  transfer_obs %>%
    left_join(transfer_holdout, by = c("mask_label", "run_id", "line_id", "direction", "summary_name")) %>%
    left_join(oracle_holdout, by = c("mask_label", "run_id", "line_id", "direction", "summary_name")) %>%
    left_join(oracle_obs, by = c("mask_label", "run_id", "line_id", "direction", "summary_name")) %>%
    mutate(
      null_holdout = observed_transfer,
      abs_null_err = pmax(abs(null_holdout - oracle_holdout), 1e-12),
      abs_transfer_err = pmax(abs(transfer_holdout - oracle_holdout), 1e-12),
      log_error_improvement = log(abs_null_err) - log(abs_transfer_err),
      transfer_better_than_null = abs_transfer_err < abs_null_err
    )
}

compute_point_loglik_subset <- function(draws, well_idx) {
  ll <- get_well_loglik_draws(draws, well_idx)
  if (!length(ll)) {
    return(NA_real_)
  }
  as.numeric(ll[1])
}

build_transfer_subset_diagnostics <- function(run_df, stan_data_path) {
  if (!nrow(run_df)) {
    return(tibble())
  }

  stan_data <- readRDS(resolve_stan_data_path(stan_data_path))
  all_wells <- seq_len(stan_data$N_wells)
  rows <- list()
  idx <- 1L

  for (i in seq_len(nrow(run_df))) {
    best_path <- file.path(run_df$base_dir[i], "transfer_best_start_summary.Rds")
    best_df <- safe_read_rds(best_path)
    if (is.null(best_df) || !nrow(best_df)) {
      next
    }

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
            model_id = run_df$run_id[i],
            line_id = line_id,
            direction = direction,
            fit_type = fit_type,
            output_root = run_df$mask_root[i],
            model_name = MODEL_NAME
          ),
          error = function(e) NULL
        )

        if (is.null(fit_obj)) {
          next
        }

        rows[[idx]] <- data.frame(
          mask_label = run_df$mask_label[i],
          run_id = run_df$run_id[i],
          line_id = as.integer(line_id),
          direction = direction,
          fit_type = fit_type,
          ll_transfer_train = compute_point_loglik_subset(fit_obj$draws, transfer_train_wells),
          ll_transfer_holdout = compute_point_loglik_subset(fit_obj$draws, transfer_holdout_wells),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  bind_rows(rows)
}

summarize_transfer_subset_diagnostics <- function(diag_df) {
  if (!nrow(diag_df)) {
    return(tibble())
  }

  wide <- diag_df %>%
    tidyr::pivot_wider(
      names_from = fit_type,
      values_from = c(ll_transfer_train, ll_transfer_holdout)
    )

  wide %>%
    transmute(
      mask_label,
      run_id,
      line_id,
      direction,
      transfer_minus_oracle_train = ll_transfer_train_transfer - ll_transfer_train_oracle,
      transfer_minus_null_train = ll_transfer_train_transfer - ll_transfer_train_null,
      transfer_minus_oracle_holdout = ll_transfer_holdout_transfer - ll_transfer_holdout_oracle,
      transfer_minus_null_holdout = ll_transfer_holdout_transfer - ll_transfer_holdout_null
    )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

run_df <- discover_mask_runs(SCREEN_ROOT, MODEL_NAME)
if (!nrow(run_df)) {
  stop(sprintf("No summarized masked transfer runs found under %s", SCREEN_ROOT))
}

predictive_df <- collect_predictive_summaries(run_df)
parameter_df <- collect_parameter_summaries(run_df)
state_df <- collect_parameter_states(run_df)
derived_df <- derive_glucose_transfer_summaries(state_df)
derived_cmp_df <- build_derived_transfer_comparison(derived_df)
subset_diag_df <- build_transfer_subset_diagnostics(run_df, stan_data_path = STAN_DATA_PATH)
subset_summary_df <- summarize_transfer_subset_diagnostics(subset_diag_df)

predictive_mask_summary <- predictive_df %>%
  group_by(mask_label, direction) %>%
  summarise(
    mean_transfer_gain = mean(transfer_gain, na.rm = TRUE),
    sum_transfer_gain = sum(transfer_gain, na.rm = TRUE),
    mean_transfer_regret = mean(transfer_regret, na.rm = TRUE),
    prop_transfer_beats_null = mean(transfer > null, na.rm = TRUE),
    n_cases = dplyr::n(),
    .groups = "drop"
  )

parameter_mask_summary <- parameter_df %>%
  group_by(mask_label, direction) %>%
  summarise(
    mean_log_error_improvement = mean(mean_log_error_improvement, na.rm = TRUE),
    mean_prop_parameters_better = mean(prop_parameters_transfer_better_than_null, na.rm = TRUE),
    n_cases = dplyr::n(),
    .groups = "drop"
  )

derived_mask_summary <- derived_cmp_df %>%
  group_by(mask_label, direction, summary_name) %>%
  summarise(
    mean_log_error_improvement = mean(log_error_improvement, na.rm = TRUE),
    mean_win_rate = mean(transfer_better_than_null, na.rm = TRUE),
    n_cases = dplyr::n(),
    .groups = "drop"
  )

subset_mask_summary <- subset_summary_df %>%
  group_by(mask_label, direction) %>%
  summarise(
    mean_transfer_minus_oracle_train = mean(transfer_minus_oracle_train, na.rm = TRUE),
    mean_transfer_minus_null_train = mean(transfer_minus_null_train, na.rm = TRUE),
    mean_transfer_minus_oracle_holdout = mean(transfer_minus_oracle_holdout, na.rm = TRUE),
    mean_transfer_minus_null_holdout = mean(transfer_minus_null_holdout, na.rm = TRUE),
    n_cases = dplyr::n(),
    .groups = "drop"
  )

write.csv(predictive_mask_summary, file.path(OUTPUT_DIR, "predictive_mask_summary.csv"), row.names = FALSE)
write.csv(parameter_mask_summary, file.path(OUTPUT_DIR, "parameter_mask_summary.csv"), row.names = FALSE)
write.csv(derived_mask_summary, file.path(OUTPUT_DIR, "derived_mask_summary.csv"), row.names = FALSE)
write.csv(subset_mask_summary, file.path(OUTPUT_DIR, "subset_mask_summary.csv"), row.names = FALSE)

plot_predictive <- ggplot(
  predictive_mask_summary,
  aes(x = reorder(mask_label, mean_transfer_gain), y = mean_transfer_gain, color = direction)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5) +
  coord_flip() +
  facet_wrap(~ direction, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Mean transfer gain")

plot_parameter <- ggplot(
  parameter_mask_summary,
  aes(x = reorder(mask_label, mean_log_error_improvement), y = mean_log_error_improvement, color = direction)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5) +
  coord_flip() +
  facet_wrap(~ direction, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Mean log-error improvement")

plot_subset <- ggplot(
  subset_mask_summary,
  aes(x = reorder(mask_label, mean_transfer_minus_null_holdout), y = mean_transfer_minus_null_holdout, color = direction)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5) +
  coord_flip() +
  facet_wrap(~ direction, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Transfer - Null on transfer-holdout wells")

plot_derived <- ggplot(
  derived_mask_summary,
  aes(x = summary_name, y = mask_label, fill = mean_log_error_improvement)
) +
  geom_tile(color = "white") +
  facet_wrap(~ direction) +
  scale_fill_gradient2(low = "#d95f02", mid = "white", high = "#1b9e77", midpoint = 0) +
  theme_minimal() +
  labs(x = "Derived glucose summary", y = "Mask", fill = "Mean log-error\nimprovement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "predictive_mask_ranking.png"), plot_predictive, width = 11, height = 8, dpi = 200)
ggsave(file.path(OUTPUT_DIR, "parameter_mask_ranking.png"), plot_parameter, width = 11, height = 8, dpi = 200)
ggsave(file.path(OUTPUT_DIR, "subset_mask_ranking.png"), plot_subset, width = 11, height = 8, dpi = 200)
ggsave(file.path(OUTPUT_DIR, "derived_mask_heatmap.png"), plot_derived, width = 12, height = 9, dpi = 200)

combined <- (plot_predictive | plot_parameter) / (plot_subset | plot_derived)
ggsave(file.path(OUTPUT_DIR, "single_param_screen_dashboard.png"), combined, width = 16, height = 12, dpi = 200)

cat(sprintf(">>> Wrote report outputs to %s\n", OUTPUT_DIR))
cat(">>> Top predictive masks by direction:\n")
print(predictive_mask_summary %>% arrange(direction, desc(mean_transfer_gain)))
