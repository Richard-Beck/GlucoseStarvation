args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "report_exports", "wp2_morphology_volume")
}

figure_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("figures", "wp2_morphology_volume")
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
  library(patchwork)
})

source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/model_free_ploidy_utils.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

sum_fuse_line <- "SUM-159-fuse"

ploidy_palette <- c(
  "baseline" = "#2C7BB6",
  "elevated" = "#D7191C"
)

metric_palette <- c(
  "Area q50" = "#2C7BB6",
  "Area q90" = "#1A9641",
  "Volume q50 = area^(3/2)" = "#D7191C"
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

format_ploidy_abs <- function(x) {
  out <- ifelse(abs(x - round(x)) < 1e-8, sprintf("%.0fN", x), sprintf("%.1fN", x))
  out[!is.finite(x)] <- NA_character_
  out
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_quantile <- function(x, prob = 0.9) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 8))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  sqrt(mean(x^2))
}

scale_vector <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) {
    sig <- 1
  }
  (x - mu) / sig
}

resolve_area_stan_data_path <- function() {
  candidates <- c(
    file.path("data", "stan_ready_data_with_area.Rds"),
    file.path("dev_data", "stan_ready_data_with_area.Rds")
  )

  for (path in candidates) {
    if (file.exists(path)) {
      return(path)
    }
  }

  stop(
    "Could not find an area-aware stan data object. Tried: ",
    paste(candidates, collapse = ", ")
  )
}

smooth_derivative <- function(hours, value, spar = 0.55) {
  ok <- is.finite(hours) & is.finite(value)
  hours <- hours[ok]
  value <- value[ok]

  if (length(hours) < 4L || length(unique(hours)) < 4L) {
    return(tibble(
      hours = hours,
      value_smooth = value,
      deriv = c(NA_real_, diff(value) / diff(hours))
    ))
  }

  ord <- order(hours)
  hours <- hours[ord]
  value <- value[ord]

  fit <- try(stats::smooth.spline(hours, value, spar = spar), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(
      hours = hours,
      value_smooth = value,
      deriv = c(NA_real_, diff(value) / diff(hours))
    ))
  }

  tibble(
    hours = hours,
    value_smooth = as.numeric(stats::predict(fit, x = hours, deriv = 0)$y),
    deriv = as.numeric(stats::predict(fit, x = hours, deriv = 1)$y)
  )
}

build_wp2_tables <- function(stan_path) {
  stan <- readRDS(stan_path)

  if (is.null(stan$area_alive_quantiles)) {
    stop("Loaded stan data does not contain area_alive_quantiles.")
  }

  line_names <- names(stan$line_map)[match(stan$line_id, unname(stan$line_map))]
  well_meta <- tibble(
    well_idx = seq_len(stan$N_wells),
    cellLine = line_names,
    line_id = as.integer(stan$line_id),
    ploidy_metric = as.numeric(stan$ploidy_metric),
    ploidy_abs = as.numeric(stan$ploidy_abs),
    G0 = as.numeric(stan$G0_per_well),
    exp_id = as.integer(stan$exp_id),
    has_starvation = as.integer(stan$has_starvation)
  ) %>%
    group_by(cellLine) %>%
    mutate(
      ploidy_state = if_else(ploidy_metric == min(ploidy_metric, na.rm = TRUE), "baseline", "elevated"),
      ploidy_abs_label = format_ploidy_abs(ploidy_abs),
      ploidy_display = paste0(ploidy_state, ": ", ploidy_abs_label)
    ) %>%
    ungroup()

  area_q <- as.data.frame(stan$area_alive_quantiles, check.names = FALSE)

  count_obs <- tibble(
    well_idx = as.integer(stan$well_idx_count),
    rep_id = as.character(stan$rep_id_count),
    hours = as.numeric(stan$t_grid[stan$grid_idx_count]),
    live_cells = as.numeric(stan$N_obs),
    dead_cells = as.numeric(stan$D_obs),
    n_area_alive = as.numeric(stan$n_area_alive %||% rep(NA_real_, length(stan$N_obs))),
    q0 = area_q[["q0"]],
    q01 = area_q[["q0.01"]],
    q05 = area_q[["q0.05"]],
    q10 = area_q[["q0.1"]],
    q25 = area_q[["q0.25"]],
    q50 = area_q[["q0.5"]],
    q75 = area_q[["q0.75"]],
    q90 = area_q[["q0.9"]],
    q95 = area_q[["q0.95"]],
    q99 = area_q[["q0.99"]],
    q100 = area_q[["q1"]]
  ) %>%
    mutate(
      total_cells = live_cells + dead_cells,
      dead_fraction = if_else(total_cells > 0, dead_cells / total_cells, NA_real_),
      tail_ratio = q99 / pmax(q50, 1e-8),
      volume_q50 = q50^1.5,
      log_live_density = log1p(live_cells)
    ) %>%
    left_join(well_meta, by = "well_idx")

  glucose_obs <- tibble(
    well_idx = as.integer(stan$well_idx_gluc),
    hours = as.numeric(stan$t_grid[stan$grid_idx_gluc]),
    lum = as.numeric(stan$lum_obs),
    dilution = as.numeric(stan$dilution),
    is_censored = as.integer(stan$is_censored)
  ) %>%
    left_join(well_meta %>% select(well_idx, exp_id), by = "well_idx") %>%
    mutate(
      calib_a = as.numeric(stan$calib_a_fixed[exp_id]),
      calib_b = as.numeric(stan$calib_b_fixed[exp_id]),
      glucose_hat = pmax(0, (lum - calib_b) / (calib_a * pmax(dilution, 1e-12)))
    ) %>%
    group_by(well_idx, hours) %>%
    summarise(
      glucose_hat = {
        uncens <- glucose_hat[is_censored == 0L]
        if (length(uncens)) median(uncens, na.rm = TRUE) else NA_real_
      },
      any_glucose_censored = any(is_censored == 1L, na.rm = TRUE),
      n_glucose_uncensored = sum(is_censored == 0L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(well_meta, by = "well_idx")

  phase_obs <- count_obs %>%
    group_by(well_idx, rep_id, cellLine, ploidy_metric, ploidy_state, G0) %>%
    group_modify(function(df, keys) {
      live_fit <- smooth_derivative(df$hours, df$live_cells)
      death_fit <- smooth_derivative(df$hours, df$dead_cells)
      total_ref <- df %>%
        select(hours, live_cells, dead_cells, total_cells) %>%
        distinct()

      full_join(
        live_fit %>% rename(live_cells_smooth = value_smooth, live_deriv = deriv),
        death_fit %>% rename(dead_cells_smooth = value_smooth, death_deriv = deriv),
        by = "hours"
      ) %>%
        left_join(total_ref, by = "hours") %>%
        transmute(
          hours = hours,
          live_cells_smooth = live_cells_smooth,
          dead_cells_smooth = dead_cells_smooth,
          live_deriv = live_deriv,
          death_deriv = death_deriv,
          growth_rate_proxy = pmax(live_deriv, 0) / pmax(total_cells, 1e-8),
          death_rate_proxy = pmax(death_deriv, 0) / pmax(total_cells, 1e-8)
        )
    }) %>%
    ungroup()

  peak_phase_ref <- phase_obs %>%
    group_by(well_idx, rep_id) %>%
    summarise(
      peak_growth_rate_proxy = safe_quantile(growth_rate_proxy[growth_rate_proxy > 0], prob = 0.75),
      .groups = "drop"
    )

  count_glucose_aligned <- count_obs %>%
    group_by(well_idx, rep_id) %>%
    group_modify(function(df, keys) {
      gluc_df <- glucose_obs %>%
        filter(well_idx == keys$well_idx[[1]]) %>%
        arrange(hours)

      if (!nrow(gluc_df)) {
        return(df %>%
          mutate(glucose_hat = G0))
      }

      ok <- is.finite(gluc_df$hours) & is.finite(gluc_df$glucose_hat)
      if (sum(ok) < 2L) {
        return(df %>%
          mutate(glucose_hat = G0))
      }

      df %>%
        arrange(hours) %>%
        mutate(
          glucose_hat = approx(
            x = gluc_df$hours[ok],
            y = gluc_df$glucose_hat[ok],
            xout = hours,
            rule = 2
          )$y
        )
    }) %>%
    ungroup() %>%
    left_join(
      phase_obs %>%
        select(well_idx, rep_id, hours, growth_rate_proxy, death_rate_proxy),
      by = c("well_idx", "rep_id", "hours")
    )

  cellline_ploidy_capacity <- count_glucose_aligned %>%
    group_by(cellLine, ploidy_metric, ploidy_state) %>%
    summarise(max_total_cells = max(total_cells, na.rm = TRUE), .groups = "drop")

  volume_proxy_rows <- count_glucose_aligned %>%
    left_join(peak_phase_ref, by = c("well_idx", "rep_id")) %>%
    left_join(
      cellline_ploidy_capacity,
      by = c("cellLine", "ploidy_metric", "ploidy_state")
    ) %>%
    mutate(
      crowding_frac = total_cells / pmax(max_total_cells, 1e-8),
      core_area = sqrt(pmax(q25, 1e-8) * pmax(q50, 1e-8)),
      likely_attached = hours >= 16,
      growth_phase = growth_rate_proxy > 0.25 * pmax(peak_growth_rate_proxy, 1e-8),
      no_death_phase = death_rate_proxy <= pmax(0.002, 0.2 * pmax(peak_growth_rate_proxy, 1e-8)),
      usable_for_volume = is.finite(core_area) &
        n_area_alive >= 30 &
        likely_attached &
        crowding_frac <= 0.7 &
        growth_phase &
        no_death_phase
    )

  crowding_slopes <- volume_proxy_rows %>%
    filter(usable_for_volume, is.finite(crowding_frac)) %>%
    group_by(cellLine) %>%
    group_modify(function(df, keys) {
      if (nrow(df) < 4L || dplyr::n_distinct(df$crowding_frac) < 2L) {
        return(tibble(beta_crowding = 0, n_fit = nrow(df)))
      }

      fit <- try(stats::lm(log(core_area) ~ crowding_frac, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) {
        tibble(beta_crowding = 0, n_fit = nrow(df))
      } else {
        tibble(
          beta_crowding = unname(stats::coef(fit)[["crowding_frac"]] %||% 0),
          n_fit = nrow(df)
        )
      }
    }) %>%
    ungroup()

  volume_proxy_rows <- volume_proxy_rows %>%
    left_join(crowding_slopes, by = "cellLine") %>%
    mutate(
      beta_crowding = if_else(is.finite(beta_crowding), beta_crowding, 0),
      log_area_adj = if_else(
        usable_for_volume,
        log(core_area) - beta_crowding * (crowding_frac - 0.2),
        NA_real_
      ),
      adj_area_ref = exp(log_area_adj)
    )

  well_volume_proxy <- volume_proxy_rows %>%
    group_by(well_idx, rep_id, cellLine, line_id, ploidy_metric, ploidy_abs, ploidy_state, ploidy_display, G0) %>%
    summarise(
      n_time = sum(is.finite(adj_area_ref)),
      n_usable = sum(usable_for_volume, na.rm = TRUE),
      attached_area_hat = safe_quantile(adj_area_ref[usable_for_volume], prob = 0.9),
      latent_volume_index = attached_area_hat^1.5,
      naive_q50_24h = safe_median(q50[hours >= 20 & hours <= 28]),
      naive_volume_24h = naive_q50_24h^1.5,
      .groups = "drop"
    ) %>%
    mutate(
      latent_vs_naive_ratio = latent_volume_index / naive_volume_24h
    )

  condition_volume_proxy <- well_volume_proxy %>%
    group_by(cellLine, line_id, ploidy_metric, ploidy_abs, ploidy_state, ploidy_display, G0) %>%
    summarise(
      attached_area_hat = median(attached_area_hat, na.rm = TRUE),
      latent_volume_index = median(latent_volume_index, na.rm = TRUE),
      naive_volume_24h = median(naive_volume_24h, na.rm = TRUE),
      latent_vs_naive_ratio = median(latent_vs_naive_ratio, na.rm = TRUE),
      n_reps = sum(is.finite(attached_area_hat)),
      median_usable_timepoints = median(n_usable, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    stan = stan,
    well_meta = well_meta,
    count_obs = count_obs,
    glucose_obs = glucose_obs,
    phase_obs = phase_obs,
    count_glucose_aligned = count_glucose_aligned,
    volume_proxy_rows = volume_proxy_rows,
    crowding_slopes = crowding_slopes,
    well_volume_proxy = well_volume_proxy,
    condition_volume_proxy = condition_volume_proxy
  )
}

fit_adjusted_morphology_model <- function(model_source, response_col, metric_label) {
  model_df <- model_source %>%
    mutate(
      response = .data[[response_col]],
      log_response = log(response),
      glucose_model = log1p(pmax(glucose_hat, 0)),
      hours_z = scale_vector(hours),
      log_live_density_z = scale_vector(log_live_density),
      glucose_z = scale_vector(glucose_model),
      growth_rate_z = scale_vector(growth_rate_proxy),
      cellLine = factor(cellLine, levels = sort(unique(cellLine)))
    ) %>%
    filter(
      is.finite(log_response),
      is.finite(ploidy_metric),
      is.finite(hours_z),
      is.finite(log_live_density_z),
      is.finite(glucose_z),
      is.finite(growth_rate_z),
      n_area_alive >= 30
    )

  fit <- stats::lm(
    log_response ~ ploidy_metric * cellLine +
      log_live_density_z + glucose_z + hours_z + I(hours_z^2) + growth_rate_z,
    data = model_df
  )

  model_terms <- stats::delete.response(stats::terms(fit))
  vc <- stats::vcov(fit)
  cf <- stats::coef(fit)

  effect_rows <- model_df %>%
    group_by(cellLine) %>%
    group_modify(function(df, keys) {
      ploidy_vals <- sort(unique(df$ploidy_metric[is.finite(df$ploidy_metric)]))
      if (length(ploidy_vals) < 2L) {
        return(tibble())
      }

      low <- ploidy_vals[1]
      high <- ploidy_vals[length(ploidy_vals)]
      ref <- tibble(
        ploidy_metric = c(low, high),
        cellLine = factor(rep(as.character(keys$cellLine[[1]]), 2), levels = levels(model_df$cellLine)),
        log_live_density_z = median(df$log_live_density_z, na.rm = TRUE),
        glucose_z = median(df$glucose_z, na.rm = TRUE),
        hours_z = median(df$hours_z, na.rm = TRUE),
        growth_rate_z = median(df$growth_rate_z, na.rm = TRUE)
      )

      mm <- stats::model.matrix(model_terms, ref)
      xdiff <- mm[2, , drop = FALSE] - mm[1, , drop = FALSE]
      effect <- as.numeric(xdiff %*% cf)
      se <- as.numeric(sqrt(xdiff %*% vc %*% t(xdiff)))

      tibble(
        ploidy_low = low,
        ploidy_high = high,
        ploidy_delta = high - low,
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
        n_model_rows = nrow(df)
      )
    }) %>%
    ungroup() %>%
    mutate(
      metric = metric_label,
      is_sum159_fuse = cellLine == sum_fuse_line
    )

  coeff_tbl <- as.data.frame(summary(fit)$coefficients, check.names = FALSE) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `t value`,
      p_value = `Pr(>|t|)`
    ) %>%
    mutate(metric = metric_label, response_col = response_col) %>%
    select(metric, response_col, term, estimate, std_error, statistic, p_value)

  fitted_tbl <- model_df %>%
    mutate(
      metric = metric_label,
      response_col = response_col,
      fitted_log = as.numeric(stats::fitted(fit)),
      residual_log = as.numeric(stats::residuals(fit)),
      fitted = exp(fitted_log)
    ) %>%
    select(
      metric, response_col, well_idx, rep_id, cellLine, ploidy_metric, ploidy_state, G0,
      hours, response, log_response, fitted, fitted_log, residual_log,
      live_cells, dead_cells, total_cells, log_live_density, glucose_hat,
      growth_rate_proxy, n_area_alive
    )

  residual_diag <- fitted_tbl %>%
    group_by(metric, cellLine, ploidy_state) %>%
    summarise(
      n_rows = n(),
      rmse_log = safe_rmse(residual_log),
      median_abs_resid_log = median(abs(residual_log), na.rm = TRUE),
      residual_sd_log = safe_sd(residual_log),
      .groups = "drop"
    )

  list(
    fit = fit,
    model_data = model_df,
    effects = effect_rows,
    coefficients = coeff_tbl,
    fitted = fitted_tbl,
    residual_diag = residual_diag
  )
}

summarize_line_volume <- function(well_volume_proxy, condition_volume_proxy) {
  well_summary <- well_volume_proxy %>%
    group_by(cellLine, line_id, ploidy_metric, ploidy_abs, ploidy_state) %>%
    summarise(
      attached_area_hat = median(attached_area_hat, na.rm = TRUE),
      latent_volume_index = median(latent_volume_index, na.rm = TRUE),
      naive_volume_24h = median(naive_volume_24h, na.rm = TRUE),
      n_wells_with_volume = sum(is.finite(latent_volume_index)),
      median_usable_timepoints = median(n_usable, na.rm = TRUE),
      .groups = "drop"
    )

  ratio_summary <- well_summary %>%
    group_by(cellLine) %>%
    arrange(ploidy_metric, .by_group = TRUE) %>%
    summarise(
      ploidy_low = first(ploidy_metric),
      ploidy_high = last(ploidy_metric),
      attached_area_low = first(attached_area_hat),
      attached_area_high = last(attached_area_hat),
      latent_volume_low = first(latent_volume_index),
      latent_volume_high = last(latent_volume_index),
      attached_area_ratio = attached_area_high / attached_area_low,
      latent_volume_ratio = latent_volume_high / latent_volume_low,
      log2_attached_area_ratio = log(attached_area_ratio) / log(2),
      log2_latent_volume_ratio = log(latent_volume_ratio) / log(2),
      min_n_wells_with_volume = min(n_wells_with_volume),
      min_median_usable_timepoints = min(median_usable_timepoints, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(is_sum159_fuse = cellLine == sum_fuse_line)

  condition_summary <- condition_volume_proxy %>%
    group_by(cellLine, line_id, ploidy_metric, ploidy_abs, ploidy_state) %>%
    summarise(
      attached_area_hat = median(attached_area_hat, na.rm = TRUE),
      latent_volume_index = median(latent_volume_index, na.rm = TRUE),
      n_conditions_with_volume = sum(is.finite(latent_volume_index)),
      .groups = "drop"
    )

  list(
    well_summary = well_summary,
    ratio_summary = ratio_summary,
    condition_summary = condition_summary
  )
}

summarize_sum_fuse_size_test <- function(adjusted_effects, volume_ratios) {
  q50_effects <- adjusted_effects %>%
    filter(metric == "Area q50") %>%
    select(cellLine, log2_ratio, lower_log2_ratio, upper_log2_ratio, ratio)

  out <- q50_effects %>%
    left_join(
      volume_ratios %>%
        select(cellLine, log2_attached_area_ratio, log2_latent_volume_ratio, latent_volume_ratio),
      by = "cellLine"
    ) %>%
    mutate(
      adjusted_q50_rank_smallest = rank(log2_ratio, ties.method = "min"),
      adjusted_q50_rank_abs_smallest = rank(abs(log2_ratio), ties.method = "min"),
      latent_volume_rank_smallest = rank(log2_latent_volume_ratio, ties.method = "min")
    )

  non_sum <- out %>%
    filter(cellLine != sum_fuse_line)

  out %>%
    mutate(
      non_sum_mean_log2_q50_ratio = mean(non_sum$log2_ratio, na.rm = TRUE),
      non_sum_sd_log2_q50_ratio = safe_sd(non_sum$log2_ratio),
      z_vs_non_sum_q50 = (log2_ratio - non_sum_mean_log2_q50_ratio) /
        pmax(non_sum_sd_log2_q50_ratio, 1e-8),
      non_sum_mean_log2_volume_ratio = mean(non_sum$log2_latent_volume_ratio, na.rm = TRUE),
      non_sum_sd_log2_volume_ratio = safe_sd(non_sum$log2_latent_volume_ratio),
      z_vs_non_sum_volume = (log2_latent_volume_ratio - non_sum_mean_log2_volume_ratio) /
        pmax(non_sum_sd_log2_volume_ratio, 1e-8),
      is_sum159_fuse = cellLine == sum_fuse_line
    ) %>%
    arrange(adjusted_q50_rank_smallest, cellLine)
}

summarize_phenotype_area_explanation <- function(base_stan_data_path, q50_effects) {
  if (is.null(base_stan_data_path) || !file.exists(base_stan_data_path)) {
    return(tibble())
  }

  feature_tables <- try(build_feature_panel(stan_data_path = base_stan_data_path), silent = TRUE)
  if (inherits(feature_tables, "try-error")) {
    warning("Could not build WP1 feature panel for phenotype-area explanation: ", conditionMessage(attr(feature_tables, "condition")))
    return(tibble())
  }

  signature_panel <- summarize_glucose_signatures(feature_tables$feature_panel)
  phenotype_effects <- compute_empirical_effects(signature_panel)$paired_effects_long

  joined <- phenotype_effects %>%
    inner_join(
      q50_effects %>%
        filter(metric == "Area q50") %>%
        select(cellLine, area_log2_ratio = log2_ratio),
      by = "cellLine"
    ) %>%
    filter(is.finite(effect_per_ploidy), is.finite(area_log2_ratio))

  if (!nrow(joined)) {
    return(tibble())
  }

  catalog <- get_model_free_feature_catalog() %>%
    select(feature, category, short_label, interpretation)

  joined %>%
    group_by(feature) %>%
    group_modify(function(df, keys) {
      if (nrow(df) < 4L || dplyr::n_distinct(df$area_log2_ratio) < 2L) {
        return(tibble(
          n_lines = nrow(df),
          intercept_only_rmse = NA_real_,
          area_model_rmse = NA_real_,
          rmse_reduction_fraction = NA_real_,
          r_squared = NA_real_,
          area_slope = NA_real_,
          area_slope_p_value = NA_real_
        ))
      }

      intercept_resid <- df$effect_per_ploidy - mean(df$effect_per_ploidy, na.rm = TRUE)
      area_fit <- stats::lm(effect_per_ploidy ~ area_log2_ratio, data = df)
      coef_tbl <- summary(area_fit)$coefficients
      intercept_rmse <- safe_rmse(intercept_resid)
      area_rmse <- safe_rmse(stats::residuals(area_fit))

      tibble(
        n_lines = nrow(df),
        intercept_only_rmse = intercept_rmse,
        area_model_rmse = area_rmse,
        rmse_reduction_fraction = (intercept_rmse - area_rmse) / intercept_rmse,
        r_squared = summary(area_fit)$r.squared,
        area_slope = unname(stats::coef(area_fit)[["area_log2_ratio"]]),
        area_slope_p_value = coef_tbl["area_log2_ratio", "Pr(>|t|)"]
      )
    }) %>%
    ungroup() %>%
    left_join(catalog, by = "feature") %>%
    arrange(desc(r_squared), feature)
}

make_distribution_figure <- function(count_glucose_aligned, well_volume_proxy) {
  area_rows <- count_glucose_aligned %>%
    filter(is.finite(q50), n_area_alive >= 30) %>%
    mutate(
      cellLine_wrapped = wrap_label(cellLine, width = 13),
      row_metric = "Observed area q50"
    )

  area_plot <- ggplot(
    area_rows,
    aes(ploidy_state, q50, fill = ploidy_state)
  ) +
    geom_violin(width = 0.84, alpha = 0.55, color = NA, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.9, linewidth = 0.25) +
    facet_wrap(~cellLine_wrapped, scales = "free_y", nrow = 1) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_fill_manual(values = ploidy_palette, guide = "none") +
    labs(
      title = "A. Alive-cell area distributions",
      x = NULL,
      y = "area q50 per image row"
    ) +
    wp2_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  volume_rows <- well_volume_proxy %>%
    filter(is.finite(attached_area_hat), is.finite(latent_volume_index)) %>%
    select(cellLine, ploidy_state, G0, well_idx, attached_area_hat, latent_volume_index) %>%
    pivot_longer(
      cols = c(attached_area_hat, latent_volume_index),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        attached_area_hat = "Attached area proxy",
        latent_volume_index = "Latent volume index"
      ),
      cellLine_wrapped = wrap_label(cellLine, width = 13)
    )

  volume_plot <- ggplot(
    volume_rows,
    aes(ploidy_state, value, color = ploidy_state)
  ) +
    geom_point(
      aes(group = interaction(G0, well_idx)),
      position = position_jitter(width = 0.08, height = 0),
      size = 1.15,
      alpha = 0.58
    ) +
    stat_summary(fun = median, geom = "crossbar", width = 0.55, color = "grey8", linewidth = 0.25) +
    facet_grid(metric ~ cellLine_wrapped, scales = "free_y") +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, guide = "none") +
    labs(
      title = "B. Crowding-adjusted area and volume proxies",
      x = NULL,
      y = "proxy value"
    ) +
    wp2_theme(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.spacing.x = grid::unit(0.45, "lines")
    )

  area_plot / volume_plot +
    plot_layout(heights = c(1, 1.35)) +
    plot_annotation(
      title = "Figure 3. Morphology And Volume Separate Most Ploidy Pairs",
      caption = "Area q50 rows are alive-cell median areas aligned to the count observations. Volume is a sensitivity index, area^(3/2), after growth-phase and crowding adjustment."
    )
}

make_adjusted_effect_figure <- function(adjusted_effects, sum_fuse_test, phenotype_area_explanation) {
  line_order <- adjusted_effects %>%
    filter(metric == "Area q50") %>%
    arrange(log2_ratio) %>%
    pull(cellLine)

  effect_plot <- adjusted_effects %>%
    mutate(
      cellLine = factor(cellLine, levels = line_order),
      cellLine_label = if_else(is_sum159_fuse, paste0(cellLine, " *"), cellLine)
    ) %>%
    ggplot(aes(cellLine, log2_ratio, color = metric)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(
      aes(ymin = lower_log2_ratio, ymax = upper_log2_ratio),
      width = 0.16,
      linewidth = 0.35,
      position = position_dodge(width = 0.48)
    ) +
    geom_point(
      aes(shape = is_sum159_fuse),
      size = 2.35,
      position = position_dodge(width = 0.48)
    ) +
    scale_color_manual(values = metric_palette, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    labs(
      title = "A. Adjusted high-versus-low ploidy morphology effect",
      x = NULL,
      y = "adjusted log2 ratio"
    ) +
    wp2_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))

  sum_plot <- sum_fuse_test %>%
    mutate(
      cellLine = factor(cellLine, levels = line_order),
      line_status = if_else(is_sum159_fuse, "SUM-159-fuse", "other lines")
    ) %>%
    select(cellLine, line_status, `Adjusted area q50` = log2_ratio, `Latent volume proxy` = log2_latent_volume_ratio) %>%
    pivot_longer(
      cols = c(`Adjusted area q50`, `Latent volume proxy`),
      names_to = "metric",
      values_to = "log2_ratio"
    ) %>%
    ggplot(aes(cellLine, log2_ratio, fill = line_status)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.18) +
    facet_wrap(~metric, nrow = 1) +
    scale_fill_manual(values = c("SUM-159-fuse" = "#D7191C", "other lines" = "grey70"), name = NULL) +
    labs(
      title = "B. SUM-159-fuse size-separation check",
      x = NULL,
      y = "high/low log2 ratio"
    ) +
    wp2_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))

  if (nrow(phenotype_area_explanation)) {
    explanation_plot <- phenotype_area_explanation %>%
      filter(is.finite(r_squared), n_lines >= 4) %>%
      slice_max(order_by = r_squared, n = 8, with_ties = FALSE) %>%
      mutate(
        display = if_else(!is.na(short_label), short_label, feature),
        display = factor(wrap_label(display, width = 18), levels = rev(wrap_label(display, width = 18)))
      ) %>%
      ggplot(aes(display, r_squared, fill = category)) +
      geom_col(width = 0.72, color = "grey25", linewidth = 0.18) +
      coord_flip() +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
      labs(
        title = "C. Starvation-effect variance explained by adjusted area separation",
        x = NULL,
        y = "line-level R^2"
      ) +
      wp2_theme(base_size = 8) +
      theme(legend.title = element_blank())
  } else {
    explanation_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "Phenotype explanation table was not available.", size = 3) +
      theme_void(base_size = 8) +
      labs(title = "C. Starvation-effect variance explained by adjusted area separation")
  }

  effect_plot / (sum_plot | explanation_plot) +
    plot_layout(heights = c(1.08, 1)) +
    plot_annotation(
      title = "Figure 4. Adjusted Morphology Effects And Phenotype Alignment",
      caption = "Effects come from log(area) models adjusted for cell line, starting/current glucose, time, live-cell density, and growth-rate proxy. SUM-159-fuse is highlighted with a triangle or red bar."
    )
}

make_diagnostic_figure <- function(fitted_values, residual_diagnostics) {
  q50_fit <- fitted_values %>%
    filter(metric == "Area q50")

  observed_fitted <- ggplot(
    q50_fit,
    aes(fitted, response, color = ploidy_state)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.28, size = 0.7) +
    facet_wrap(~cellLine, scales = "free", nrow = 1) +
    scale_x_log10(labels = label_number(big.mark = ",")) +
    scale_y_log10(labels = label_number(big.mark = ",")) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "A. Area q50 observed versus fitted",
      x = "fitted area q50",
      y = "observed area q50"
    ) +
    wp2_theme(base_size = 7)

  residual_plot <- ggplot(
    q50_fit,
    aes(fitted_log, residual_log, color = ploidy_state)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.28, size = 0.7) +
    facet_wrap(~cellLine, nrow = 1) +
    scale_color_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "B. Area q50 log-residuals",
      x = "fitted log(area q50)",
      y = "residual"
    ) +
    wp2_theme(base_size = 7)

  diag_plot <- residual_diagnostics %>%
    mutate(cellLine = factor(cellLine, levels = sort(unique(cellLine)))) %>%
    ggplot(aes(cellLine, rmse_log, fill = ploidy_state)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "grey25", linewidth = 0.18) +
    facet_wrap(~metric, nrow = 1) +
    scale_fill_manual(values = ploidy_palette, name = "ploidy") +
    labs(
      title = "C. Residual RMSE by line and ploidy state",
      x = NULL,
      y = "log-scale RMSE"
    ) +
    wp2_theme(base_size = 8) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))

  (observed_fitted / residual_plot / diag_plot) +
    plot_layout(heights = c(1, 1, 0.95)) +
    plot_annotation(
      title = "Supplement. Morphology Model Residual Diagnostics"
    )
}

write_wp2_text_report <- function(
  path,
  stan_path,
  base_stan_data_path,
  figure_dir,
  output_dir,
  coverage_tbl,
  adjusted_effects,
  sum_fuse_test,
  line_volume_ratios,
  phenotype_area_explanation,
  residual_diagnostics
) {
  q50_effects <- adjusted_effects %>%
    filter(metric == "Area q50") %>%
    arrange(desc(log2_ratio))

  n_positive_q50 <- sum(q50_effects$lower_log2_ratio > 0, na.rm = TRUE)
  n_effect_rows <- sum(is.finite(q50_effects$log2_ratio))

  sum_row <- sum_fuse_test %>%
    filter(cellLine == sum_fuse_line)

  top_explained <- phenotype_area_explanation %>%
    filter(is.finite(r_squared), n_lines >= 4) %>%
    arrange(desc(r_squared)) %>%
    head(10)

  lines <- c(
    "WP2_MORPHOLOGY_VOLUME_ANALYSIS",
    sprintf("generated\t%s", Sys.time()),
    sprintf("area_stan_data\t%s", stan_path),
    sprintf("base_stan_data_for_wp1_features\t%s", base_stan_data_path %||% NA_character_),
    sprintf("figure_dir\t%s", figure_dir),
    sprintf("table_dir\t%s", output_dir),
    "",
    "SECTION\tCOVERAGE",
    capture.output(write.table(as.data.frame(coverage_tbl), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tPASS_CRITERION_QUICK_CHECK",
    sprintf("adjusted_q50_effects_with_95pct_interval_above_zero\t%d/%d", n_positive_q50, n_effect_rows),
    if (nrow(sum_row)) {
      c(
        sprintf("SUM159_fuse_adjusted_q50_log2_ratio\t%.4f", sum_row$log2_ratio[[1]]),
        sprintf("SUM159_fuse_adjusted_q50_rank_smallest\t%d", as.integer(sum_row$adjusted_q50_rank_smallest[[1]])),
        sprintf("SUM159_fuse_adjusted_q50_rank_abs_smallest\t%d", as.integer(sum_row$adjusted_q50_rank_abs_smallest[[1]])),
        sprintf("SUM159_fuse_latent_volume_log2_ratio\t%.4f", sum_row$log2_latent_volume_ratio[[1]]),
        sprintf("SUM159_fuse_z_vs_non_sum_q50\t%.4f", sum_row$z_vs_non_sum_q50[[1]])
      )
    } else {
      "SUM-159-fuse was not found in the adjusted effect table."
    },
    "",
    "SECTION\tADJUSTED_PLOIDY_EFFECTS",
    capture.output(write.table(as.data.frame(adjusted_effects), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tLATENT_VOLUME_RATIOS",
    capture.output(write.table(as.data.frame(line_volume_ratios), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tSUM159_FUSE_SIZE_SEPARATION_TEST",
    capture.output(write.table(as.data.frame(sum_fuse_test), sep = "\t", row.names = FALSE, quote = FALSE)),
    "",
    "SECTION\tTOP_PHENOTYPE_EFFECTS_EXPLAINED_BY_AREA",
    if (nrow(top_explained)) {
      capture.output(write.table(as.data.frame(top_explained), sep = "\t", row.names = FALSE, quote = FALSE))
    } else {
      "Phenotype-area explanation table was unavailable or empty."
    },
    "",
    "SECTION\tRESIDUAL_DIAGNOSTICS",
    capture.output(write.table(as.data.frame(residual_diagnostics), sep = "\t", row.names = FALSE, quote = FALSE))
  )

  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

stan_path <- resolve_area_stan_data_path()
base_stan_data_path <- tryCatch(resolve_stan_data_path(), error = function(e) NA_character_)

message("Building WP2 morphology tables from ", stan_path)
tables <- build_wp2_tables(stan_path)

coverage_tbl <- tibble(
  metric = c(
    "count rows",
    "rows with any area summary",
    "rows with finite q50",
    "rows with n_area_alive >= 30",
    "wells with latent volume proxy"
  ),
  value = c(
    nrow(tables$count_obs),
    sum(rowSums(is.na(tables$stan$area_alive_quantiles)) < ncol(tables$stan$area_alive_quantiles)),
    sum(is.finite(tables$count_obs$q50)),
    sum(is.finite(tables$count_obs$q50) & tables$count_obs$n_area_alive >= 30, na.rm = TRUE),
    sum(is.finite(tables$well_volume_proxy$latent_volume_index))
  )
)

model_source <- tables$count_glucose_aligned %>%
  mutate(volume_q50 = q50^1.5)

model_specs <- tibble(
  response_col = c("q50", "q90", "volume_q50"),
  metric = c("Area q50", "Area q90", "Volume q50 = area^(3/2)")
)

model_results <- lapply(seq_len(nrow(model_specs)), function(i) {
  fit_adjusted_morphology_model(
    model_source = model_source,
    response_col = model_specs$response_col[[i]],
    metric_label = model_specs$metric[[i]]
  )
})

adjusted_effects <- bind_rows(lapply(model_results, `[[`, "effects"))
coefficient_table <- bind_rows(lapply(model_results, `[[`, "coefficients"))
fitted_values <- bind_rows(lapply(model_results, `[[`, "fitted"))
residual_diagnostics <- bind_rows(lapply(model_results, `[[`, "residual_diag"))

volume_summaries <- summarize_line_volume(
  well_volume_proxy = tables$well_volume_proxy,
  condition_volume_proxy = tables$condition_volume_proxy
)

sum_fuse_test <- summarize_sum_fuse_size_test(
  adjusted_effects = adjusted_effects,
  volume_ratios = volume_summaries$ratio_summary
)

phenotype_area_explanation <- summarize_phenotype_area_explanation(
  base_stan_data_path = base_stan_data_path,
  q50_effects = adjusted_effects
)

write.csv(coverage_tbl, file.path(output_dir, "wp2_area_coverage.csv"), row.names = FALSE)
write.csv(tables$count_glucose_aligned, file.path(output_dir, "wp2_area_count_observations.csv"), row.names = FALSE)
write.csv(tables$well_volume_proxy, file.path(output_dir, "wp2_well_volume_proxy.csv"), row.names = FALSE)
write.csv(tables$condition_volume_proxy, file.path(output_dir, "wp2_condition_volume_proxy.csv"), row.names = FALSE)
write.csv(volume_summaries$well_summary, file.path(output_dir, "wp2_line_ploidy_volume_summary.csv"), row.names = FALSE)
write.csv(volume_summaries$ratio_summary, file.path(output_dir, "wp2_line_volume_ratios.csv"), row.names = FALSE)
write.csv(adjusted_effects, file.path(output_dir, "wp2_adjusted_ploidy_effects.csv"), row.names = FALSE)
write.csv(coefficient_table, file.path(output_dir, "wp2_model_coefficients.csv"), row.names = FALSE)
write.csv(fitted_values, file.path(output_dir, "wp2_model_fitted_values.csv"), row.names = FALSE)
write.csv(residual_diagnostics, file.path(output_dir, "wp2_residual_diagnostics.csv"), row.names = FALSE)
write.csv(sum_fuse_test, file.path(output_dir, "wp2_sum159_fuse_size_test.csv"), row.names = FALSE)
write.csv(phenotype_area_explanation, file.path(output_dir, "wp2_phenotype_area_explanation.csv"), row.names = FALSE)
write.csv(tables$crowding_slopes, file.path(output_dir, "wp2_crowding_slopes.csv"), row.names = FALSE)

saveRDS(tables$count_glucose_aligned, file.path(output_dir, "wp2_area_count_observations.Rds"))
saveRDS(tables$well_volume_proxy, file.path(output_dir, "wp2_well_volume_proxy.Rds"))
saveRDS(tables$condition_volume_proxy, file.path(output_dir, "wp2_condition_volume_proxy.Rds"))
saveRDS(volume_summaries, file.path(output_dir, "wp2_volume_summaries.Rds"))
saveRDS(adjusted_effects, file.path(output_dir, "wp2_adjusted_ploidy_effects.Rds"))
saveRDS(coefficient_table, file.path(output_dir, "wp2_model_coefficients.Rds"))
saveRDS(residual_diagnostics, file.path(output_dir, "wp2_residual_diagnostics.Rds"))
saveRDS(sum_fuse_test, file.path(output_dir, "wp2_sum159_fuse_size_test.Rds"))
saveRDS(phenotype_area_explanation, file.path(output_dir, "wp2_phenotype_area_explanation.Rds"))

figure_3 <- make_distribution_figure(
  count_glucose_aligned = tables$count_glucose_aligned,
  well_volume_proxy = tables$well_volume_proxy
)
save_plot_pair(figure_3, "figure_3_area_volume_distributions", width = 12.8, height = 8.6)

figure_4 <- make_adjusted_effect_figure(
  adjusted_effects = adjusted_effects,
  sum_fuse_test = sum_fuse_test,
  phenotype_area_explanation = phenotype_area_explanation
)
save_plot_pair(figure_4, "figure_4_adjusted_area_volume_effects", width = 12.8, height = 8.8)

diagnostic_figure <- make_diagnostic_figure(
  fitted_values = fitted_values,
  residual_diagnostics = residual_diagnostics
)
save_plot_pair(diagnostic_figure, "supplement_morphology_model_diagnostics", width = 13.5, height = 10.2)

write_wp2_text_report(
  path = file.path(output_dir, "wp2_summary.txt"),
  stan_path = stan_path,
  base_stan_data_path = base_stan_data_path,
  figure_dir = figure_dir,
  output_dir = output_dir,
  coverage_tbl = coverage_tbl,
  adjusted_effects = adjusted_effects,
  sum_fuse_test = sum_fuse_test,
  line_volume_ratios = volume_summaries$ratio_summary,
  phenotype_area_explanation = phenotype_area_explanation,
  residual_diagnostics = residual_diagnostics
)

cat(sprintf("Wrote WP2 figures to %s\n", figure_dir))
cat(sprintf("Wrote WP2 audit tables to %s\n", output_dir))
