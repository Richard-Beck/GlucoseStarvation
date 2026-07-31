#!/usr/bin/env Rscript

# FG3 mechanistic diagnostics: canonical red_a30_counts_20260722 rebuild.
# All analytical transformations and panel objects remain in memory. The only
# outputs are the six final PNGs and one operation-timing TSV.

options(warn = 1)
suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(scales)
  library(tidyr)
  library(tibble)
})
source("R/project_paths.R")
source("R/gpath_run_utils.R")
source("R/sim_utils.R")

total_start <- proc.time()[["elapsed"]]
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
release_root <- file.path(project_root, "data", "modelling", "gpath_v1", "red_a30_counts_20260722")
package_root <- file.path(project_root, "manuscript_figures", "20260722T141236_resegmentation_fix")
final_dir <- file.path(package_root, "final_images")
timing_dir <- file.path(package_root, "timings")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)

raster_root <- file.path(project_root, "figures", "user-approved-raster-figures")
raster_candidates <- list.files(raster_root, pattern = "\\.png$", full.names = TRUE)
if (length(raster_candidates) != 1L) {
  stop("Expected exactly one canonical immutable raster in ", raster_root, call. = FALSE)
}
f3a_path <- raster_candidates[[1]]
f3a_sha <- strsplit(system2("sha256sum", f3a_path, stdout = TRUE), "[[:space:]]+")[[1]][[1]]
f3a_expected_sha <- "8f222b659b718ecd74e27df08ca6e4abd66f83564a7f74892d3facc38bc60c37"
if (
  basename(f3a_path) != "model_family_schematic_placeholder.png" ||
  !identical(tolower(f3a_sha), f3a_expected_sha)
) {
  stop("Canonical raster does not match the approved F3A asset name and SHA-256.", call. = FALSE)
}

timings <- list()
time_operation <- function(operation, details, expression) {
  start <- proc.time()[["elapsed"]]
  value <- force(expression)
  elapsed <- proc.time()[["elapsed"]] - start
  timings[[length(timings) + 1L]] <<- tibble(
    package_id = "FG3_mechanistic_diagnostics",
    operation = operation,
    elapsed_seconds = round(elapsed, 3),
    details = details
  )
  value
}

model_ids <- c(
  "1R_1P_0W_C0_M1", "1R_1P_1W_C0_M0", "2R_1P_0W_C0_M1",
  "2R_2P_0W_C0_M1", "2R_2P_1W_C0_M0"
)
model_alias <- c(
  "1R_1P_0W_C0_M1" = "1R",
  "1R_1P_1W_C0_M0" = "1R,W(m)",
  "2R_1P_0W_C0_M1" = "2R(a)",
  "2R_2P_0W_C0_M1" = "2R(f)",
  "2R_2P_1W_C0_M0" = "2R(f),W(m)"
)
model_colors <- c(
  "1R" = "#DF6C8B", "1R,W(m)" = "#B79A00", "2R(a)" = "#18A765",
  "2R(f)" = "#08A9C4", "2R(f),W(m)" = "#B36FE0"
)
parse_model <- function(ids) {
  bind_rows(lapply(ids, function(id) {
    hit <- regmatches(id, regexec("^(\\d+)R_(\\d+)P_(\\d+)W_C(\\d+)_M(\\d+)$", id))[[1]]
    if (length(hit) != 6L) stop("Unparseable model id: ", id, call. = FALSE)
    tibble(
      model_id = id, R = as.integer(hit[2]), P = as.integer(hit[3]),
      W = as.integer(hit[4]), C = as.integer(hit[5]), M = as.integer(hit[6])
    )
  }))
}
alias_all <- function(ids) {
  dims <- parse_model(ids)
  base <- case_when(
    dims$R == 1L ~ "1R",
    dims$R == 2L & dims$P == 1L ~ "2R(a)",
    dims$R == 2L & dims$P == 2L & dims$C == 0L ~ "2R(f)",
    dims$R == 2L & dims$P == 2L & dims$C == 1L ~ "2R(l)",
    TRUE ~ dims$model_id
  )
  waste <- case_when(
    dims$W == 0L ~ "",
    dims$W == 1L & dims$M == 0L ~ ",W(m)",
    dims$W == 1L & dims$M == 1L ~ ",W(a)",
    TRUE ~ ""
  )
  setNames(paste0(base, waste), dims$model_id)
}
line_labels <- c(
  "MCF10A" = "MCF10A", "MDA-MB-231" = "MDA-MB-231", "SNU668" = "SNU668",
  "SUM-159-chem" = "SUM-159-chem", "SUM-159-fuse" = "SUM-159-fuse"
)
single_dataset_ids <- c(
  "MCF10A" = "single_mcf10a",
  "MDA-MB-231" = "single_mda_mb_231",
  "SNU668" = "single_snu668",
  "SUM-159-chem" = "single_sum_159_chem",
  "SUM-159-fuse" = "single_sum_159_fuse"
)

theme_fg3 <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.22, colour = "#E6E6E6"),
      strip.background = element_rect(fill = "#EEEEEE", colour = "#999999", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = rel(0.92)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(colour = "#333333"),
      legend.key.height = grid::unit(3, "mm"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

loaded <- time_operation(
  "data_loading",
  paste(
    "optimization/assessment.Rds; posterior/qc.Rds; all-lines Stan input;",
    "five all-lines posterior prediction objects"
  ),
  {
    assessment <- readRDS(file.path(release_root, "derived", "optimization", "assessment.Rds"))
    posterior_qc <- readRDS(file.path(release_root, "derived", "posterior", "qc.Rds"))
    all_predictions <- setNames(lapply(model_ids, function(model_id) {
      readRDS(file.path(release_root, "derived", "posterior", "predictions", "all_lines", paste0(model_id, ".Rds")))
    }), model_ids)
    stan_all <- readRDS(file.path(release_root, "datasets", "all_lines", "stan_data.Rds"))
    list(
      assessment = assessment, posterior_qc = posterior_qc,
      all_predictions = all_predictions, stan_all = stan_all
    )
  }
)
assessment <- loaded$assessment
posterior_qc <- loaded$posterior_qc
all_predictions <- loaded$all_predictions
stan_all <- loaded$stan_all
rm(loaded)

f3_s5 <- time_operation(
  "transformations_f3_s5",
  "Pareto diagnostics, no-SUM-159-fuse per-well likelihood loss, optimization stability, AIC and complexity decomposition",
  {
    fit_all <- assessment$fit_summary %>%
      filter(.data$dataset_id == "all_lines") %>%
      mutate(
        alias = unname(model_alias[.data$model_id]),
        selected = .data$model_id %in% model_ids
      )

    no_sum_stan <- readRDS(file.path(release_root, "datasets", "loo_exclude_sum_159_fuse", "stan_data.Rds"))
    line_map <- tibble(line_name = names(no_sum_stan$line_map), line_id = as.integer(no_sum_stan$line_map))
    well_meta <- tibble(
      well_idx = seq_len(no_sum_stan$N_wells),
      line_id = as.integer(no_sum_stan$line_id),
      ploidy = as.numeric(no_sum_stan$ploidy_abs),
      G0 = as.numeric(no_sum_stan$G0_per_well)
    ) %>%
      left_join(line_map, by = "line_id") %>%
      mutate(condition = paste0(.data$line_name, "\n", format(.data$ploidy, trim = TRUE), "N"))

    ll_rows <- bind_rows(lapply(model_ids, function(model_id) {
      fit_row <- assessment$fit_summary %>%
        filter(.data$dataset_id == "loo_exclude_sum_159_fuse", .data$model_id == .env$model_id) %>%
        slice(1)
      fit_id <- fit_row$fit_id[[1]]
      best_start <- fit_row$best_start[[1]]
      draw_path <- file.path(
        release_root, "datasets", "loo_exclude_sum_159_fuse", "optim", model_id,
        fit_id, paste0("optim_draws_", best_start, ".Rds")
      )
      draw <- readRDS(draw_path)
      values <- as.numeric(draw[1, ])
      names(values) <- colnames(draw)
      ll_names <- sprintf("ll_well[%d]", seq_len(no_sum_stan$N_wells))
      if (!all(ll_names %in% names(values))) {
        stop("Canonical best optimization draw lacks ll_well values for ", model_id, call. = FALSE)
      }
      tibble(model_id = model_id, well_idx = seq_len(no_sum_stan$N_wells), ll = values[ll_names])
    })) %>%
      left_join(well_meta, by = "well_idx") %>%
      group_by(.data$model_id, .data$line_name, .data$ploidy, .data$G0, .data$condition) %>%
      summarise(ll = mean(.data$ll), .groups = "drop") %>%
      group_by(.data$line_name, .data$ploidy, .data$G0, .data$condition) %>%
      mutate(loss = pmax(0, -2 * (.data$ll - max(.data$ll)))) %>%
      ungroup() %>%
      mutate(
        alias = factor(unname(model_alias[.data$model_id]), levels = unname(model_alias)),
        G0_label = factor(format(.data$G0, trim = TRUE), levels = format(sort(unique(.data$G0)), trim = TRUE)),
        condition = factor(.data$condition, levels = rev(unique(well_meta$condition)))
      )
    list(fit_all = fit_all, ll_rows = ll_rows)
  }
)
fit_all <- f3_s5$fit_all
ll_rows <- f3_s5$ll_rows

panel_f3a <- ggdraw() + draw_image(f3a_path, x = 0, y = 0, width = 1, height = 1)
panel_f3b <- ggplot(fit_all, aes(x = .data$k, y = .data$deviance)) +
  geom_point(colour = "#CFCFCF", size = 1.6) +
  geom_point(data = filter(fit_all, .data$selected), shape = 21, fill = "#EAF5FC", colour = "#0072B2", size = 3.2, stroke = 0.9) +
  geom_text(
    data = filter(fit_all, .data$selected),
    aes(label = .data$alias), size = 2.8, nudge_y = -250, check_overlap = TRUE
  ) +
  labs(x = "effective parameter count k", y = "deviance\n(lower is better)") +
  theme_fg3(9) + theme(legend.position = "none")

loss_cap <- quantile(ll_rows$loss, 0.95, na.rm = TRUE)
panel_f3c <- ggplot(ll_rows, aes(x = .data$G0_label, y = .data$condition, fill = pmin(.data$loss, loss_cap))) +
  geom_tile(colour = "white", linewidth = 0.18) +
  facet_wrap(~alias, nrow = 1) +
  scale_fill_gradient(low = "#F6F6F2", high = "#B94D49", oob = squish, name = "Deviance loss\nper well") +
  labs(x = "starting glucose G0", y = NULL) +
  theme_fg3(7) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position = "bottom")

all_alias <- alias_all(fit_all$model_id)
s5_fit <- fit_all %>%
  mutate(
    alias12 = unname(all_alias[.data$model_id]),
    candidate = factor(
      if_else(.data$pareto_member, "Pareto-front candidate", "Other assessed model"),
      levels = c("Pareto-front candidate", "Other assessed model")
    ),
    stability = .data$n_within_1 / .data$n_starts_expected,
    stability_label = paste0(.data$n_within_1, "/", .data$n_starts_expected),
    stability_name = paste0(.data$alias12, "  (k=", .data$k, ")")
  )
s5a_data <- s5_fit %>%
  arrange(.data$stability, .data$k, .data$model_id) %>%
  mutate(stability_name = factor(.data$stability_name, levels = .data$stability_name))
panel_s5a <- ggplot(s5a_data, aes(x = .data$stability, y = .data$stability_name, fill = .data$candidate)) +
  geom_col(width = 0.62, colour = "grey45", linewidth = 0.12) +
  geom_text(aes(label = .data$stability_label), hjust = -0.1, size = 2.15, colour = "grey35") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(s5a_data$stability) * 1.24), expand = expansion(mult = c(0, 0.035))) +
  scale_fill_manual(values = c("Pareto-front candidate" = "#007C89", "Other assessed model" = "#B8B8B8"), name = NULL) +
  labs(x = "Valid starts within 1", y = NULL) +
  theme_fg3(7) +
  theme(legend.position = "bottom", legend.box = "vertical", axis.text.y = element_text(size = 6.2))

s5b_data <- s5_fit %>%
  arrange(.data$delta_AIC, .data$k, .data$model_id) %>%
  mutate(alias12 = factor(.data$alias12, levels = rev(unique(.data$alias12))))
panel_s5b <- ggplot(s5b_data, aes(x = .data$delta_AIC, y = .data$alias12, fill = .data$candidate)) +
  geom_col(width = 0.66, colour = "grey34", linewidth = 0.18, show.legend = FALSE) +
  scale_fill_manual(values = c("Pareto-front candidate" = "#0072B2", "Other assessed model" = "#BDBDBD"), guide = "none") +
  scale_x_continuous(labels = label_number(big.mark = ","), breaks = seq(0, 4000, by = 1000), expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = c(0, max(s5b_data$delta_AIC) + 80), clip = "on") +
  labs(x = "Delta AIC", y = NULL) +
  theme_fg3(7) + theme(axis.text.y = element_text(size = 6.4))

dims12 <- parse_model(s5_fit$model_id) %>%
  left_join(s5_fit %>% select("model_id", "k"), by = "model_id") %>%
  mutate(
    alias12 = unname(all_alias[.data$model_id]),
    resource_kinetics = 3L * .data$R,
    resource_allocation = if_else(.data$C == 0L, (.data$P - 1L) * .data$R, 0L),
    waste_kinetics = .data$W * .data$R,
    maintenance = 1L,
    model_label = paste0(.data$alias12, "\nk=", .data$k)
  ) %>%
  arrange(.data$k, .data$R, .data$P, .data$C, .data$W, .data$M) %>%
  mutate(model_label = factor(.data$model_label, levels = rev(.data$model_label)))
complexity <- dims12 %>%
  select("model_id", "model_label", "resource_kinetics", "resource_allocation", "waste_kinetics", "maintenance") %>%
  pivot_longer(
    c("resource_kinetics", "resource_allocation", "waste_kinetics", "maintenance"),
    names_to = "component", values_to = "base_count"
  ) %>%
  mutate(
    component = factor(
      .data$component,
      levels = c("resource_kinetics", "resource_allocation", "waste_kinetics", "maintenance"),
      labels = c("Resource kinetics", "Resource allocation", "Waste kinetics", "Maintenance")
    )
  )
complexity_line <- complexity %>% mutate(value = .data$base_count * stan_all$N_lines)
complexity_ploidy <- complexity %>% mutate(value = .data$base_count)
complexity_colors <- c(
  "Resource kinetics" = "#0072B2", "Resource allocation" = "#009E73",
  "Waste kinetics" = "#D55E00", "Maintenance" = "#7F7F7F"
)
p_s5c_line <- ggplot(complexity_line, aes(x = .data$value, y = .data$model_label, fill = .data$component)) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = complexity_colors,
    name = NULL,
    guide = guide_legend(nrow = 2, byrow = TRUE, direction = "horizontal", reverse = TRUE)
  ) +
  labs(x = "Contribution to k from effective structural slots", y = NULL) +
  ggtitle("Line effects") + theme_fg3(6.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7), legend.position = "bottom", axis.text.y = element_text(size = 5.4))
p_s5c_ploidy <- ggplot(complexity_ploidy, aes(x = .data$value, y = .data$model_label, fill = .data$component)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = complexity_colors, guide = "none") +
  labs(x = NULL, y = NULL) + ggtitle("Ploidy effects") + theme_fg3(6.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7), axis.text.y = element_blank(), axis.ticks.y = element_blank())
panel_s5c <- (p_s5c_line | p_s5c_ploidy) + plot_layout(widths = c(5, 1), guides = "collect")

s6_panels <- time_operation(
  "transformations_s6",
  "Posterior median N/R1/R2/W1 trajectories and observed-cell overlays from the all-lines fit for five lines and five model families",
  {
    lapply(names(single_dataset_ids), function(line_name) {
      stan <- stan_all
      predictions <- all_predictions
      this_line_id <- as.integer(stan$line_map[[line_name]])
      pred_long <- bind_rows(lapply(model_ids, function(model_id) {
        pred <- predictions[[model_id]]
        med <- apply(pred$state, c(2, 3, 4), median, na.rm = TRUE)
        meta <- pred$well_metadata
        bind_rows(lapply(seq_len(dim(med)[3]), function(state_idx) {
          state_table <- as_tibble(as.data.frame.table(med[, , state_idx], responseName = "value"))
          names(state_table) <- c("well_idx", "time_idx", "value")
          state_table %>%
            transmute(
              well_idx = as.integer(.data$well_idx),
              time = pred$time[as.integer(.data$time_idx)],
              state = dimnames(med)[[3]][state_idx],
              value = .data$value
            )
        })) %>%
          left_join(meta %>% select(.data$well_idx, .data$ploidy_abs, .data$G0), by = "well_idx") %>%
          left_join(meta %>% select(.data$well_idx, .data$line_name), by = "well_idx") %>%
          filter(.data$line_name == .env$line_name) %>%
          mutate(model = factor(unname(model_alias[model_id]), levels = unname(model_alias)))
      }))
      observed <- build_obs_df_from_stan_data(stan, line_id = this_line_id, line_name = line_name) %>%
        transmute(
          well_idx = .data$well_idx, time = .data$time, state = .data$variable,
          value = .data$value, ploidy_abs = .data$ploidy, G0 = .data$G0,
          line_name = .data$line_name
        )
      display_limits <- observed %>%
        group_by(.data$state) %>%
        summarise(display_y_max = 1.5 * max(.data$value, na.rm = TRUE), .groups = "drop")
      pred_long <- pred_long %>%
        left_join(display_limits, by = "state") %>%
        mutate(value = if_else(!is.na(.data$display_y_max) & .data$value > .data$display_y_max, NA_real_, .data$value))
      state_order <- intersect(c("N", "R1", "W1", "R2"), unique(pred_long$state))
      state_plots <- lapply(seq_along(state_order), function(state_idx) {
        state_name <- state_order[[state_idx]]
        one <- filter(pred_long, .data$state == .env$state_name)
        p <- ggplot(one, aes(x = .data$time, y = .data$value, colour = .data$model, group = interaction(.data$model, .data$well_idx))) +
          geom_line(linewidth = 0.35, alpha = 0.9)
        if (state_name == "N") {
          p <- p + geom_point(data = observed, aes(x = .data$time, y = .data$value), inherit.aes = FALSE, size = 0.35, colour = "black", alpha = 0.65)
        }
        p +
          facet_grid(
            rows = vars(G0), cols = vars(ploidy_abs), scales = "free_y",
            labeller = labeller(
              G0 = function(x) paste0("G0=", x),
              ploidy_abs = function(x) paste0(x, "N")
            )
          ) +
          scale_colour_manual(values = model_colors, drop = FALSE) +
          labs(
            x = "time (h)",
            y = if (state_idx == 1L) line_name else NULL,
            colour = "Model", title = state_name
          ) +
          theme_fg3(5.4) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
            legend.position = "none", legend.text = element_text(size = 5),
            legend.title = element_text(size = 5.5), axis.text = element_text(size = 4.5),
            axis.title.y = element_text(face = "bold", size = 7),
            strip.text = element_text(size = 4.7), panel.spacing = grid::unit(0.45, "mm")
          )
      })
      wrap_elements(full = wrap_plots(state_plots, nrow = 1))
    })
  }
)

s6_legend <- cowplot::get_legend(
  ggplot(
    tibble(model = factor(unname(model_alias), levels = unname(model_alias)), x = seq_along(model_ids), y = 1),
    aes(x = .data$x, y = .data$y, colour = .data$model)
  ) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = model_colors, drop = FALSE) +
    labs(colour = "Model") +
    theme_void(base_family = "sans", base_size = 6.5) +
    theme(legend.position = "right")
)
s6_figure <- cowplot::plot_grid(
  wrap_plots(s6_panels, ncol = 1) + plot_annotation(tag_levels = "a"),
  s6_legend,
  ncol = 2,
  rel_widths = c(0.88, 0.12)
)

s7 <- time_operation(
  "transformations_s7",
  "Prior-design transfer schematic, labelled heatmap, and canonical null/transfer/full-data case overlays",
  {
    transfer <- assessment$transfer_vs_null %>%
      filter(.data$parent_dataset == "loo_exclude_sum_159_fuse") %>%
      mutate(
        line_name = factor(.data$line_name, levels = names(line_labels)),
        direction_label = recode(.data$direction, low_to_high = "low → high", high_to_low = "high → low"),
        alias = unname(alias_all(unique(assessment$transfer_vs_null$model_id))[.data$model_id])
      )
    transfer_dims <- parse_model(unique(transfer$model_id)) %>%
      arrange(.data$R, .data$P, .data$C, .data$W, .data$M, .data$model_id)
    transfer <- transfer %>%
      mutate(
        alias = factor(.data$alias, levels = unname(alias_all(transfer_dims$model_id))),
        line_name = factor(.data$line_name, levels = rev(c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem"))),
        gain_label = case_when(
          abs(.data$transfer_minus_null) >= 1000 ~ sprintf("%+.1fk", .data$transfer_minus_null / 1000),
          TRUE ~ sprintf("%+.0f", .data$transfer_minus_null)
        ),
        clipped_gain = squish(.data$transfer_minus_null, c(-1500, 1500))
      )
    step_df <- tibble(
      x = c(0.9, 2.45, 4.0), y = 1.25,
      label = c(
        "Train source\nother lines + one\nploidy state",
        "Hold out target\nopposite ploidy in\nthe test line",
        "Score prediction\nholdout log-lik:\ntransfer - null"
      )
    )
    fit_df <- tibble(
      x = c(0.9, 2.45, 4.0), y = 0.45,
      label = c("Null\nno ploidy effect", "Transfer\npredict held-out ploidy", "Oracle\nsees all data")
    )
    schematic <- ggplot() +
      geom_rect(data = step_df, aes(xmin = .data$x - 0.55, xmax = .data$x + 0.55, ymin = .data$y - 0.32, ymax = .data$y + 0.32), fill = "grey97", colour = "grey30", linewidth = 0.35) +
      geom_segment(data = tibble(x = c(1.43, 2.98), xend = c(1.82, 3.37), y = 1.25), aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y), arrow = arrow(length = unit(0.07, "in")), linewidth = 0.35) +
      geom_text(data = step_df, aes(.data$x, .data$y, label = .data$label), size = 2.35, lineheight = 0.92) +
      geom_rect(data = fit_df, aes(xmin = .data$x - 0.55, xmax = .data$x + 0.55, ymin = .data$y - 0.22, ymax = .data$y + 0.22), fill = "#F6FBFD", colour = "#5B7F93", linewidth = 0.3) +
      geom_text(data = fit_df, aes(.data$x, .data$y, label = .data$label), size = 2.25, lineheight = 0.92) +
      coord_cartesian(xlim = c(0.15, 4.75), ylim = c(0.12, 1.62), clip = "off") +
      theme_void(base_size = 8) + theme(plot.margin = margin(3, 3, 3, 3))

    case_marks <- tibble(
      line_name = factor(c("MCF10A", "SNU668"), levels = levels(transfer$line_name)),
      direction = c("low_to_high", "high_to_low"),
      model_id = c("2R_2P_1W_C0_M0", "2R_2P_0W_C0_M1"),
      case = c("C", "D")
    )
    heat_marks <- transfer %>% inner_join(case_marks, by = c("line_name", "direction", "model_id"))
    heatmap <- ggplot(transfer, aes(x = .data$alias, y = .data$line_name, fill = .data$transfer_minus_null)) +
      geom_tile(colour = "white", linewidth = 0.18) +
      geom_text(aes(label = .data$gain_label), size = 1.72, colour = "grey12") +
      geom_tile(data = heat_marks, aes(x = .data$alias, y = .data$line_name), inherit.aes = FALSE, fill = NA, colour = "grey8", linewidth = 0.72) +
      geom_label(data = heat_marks, aes(x = .data$alias, y = .data$line_name, label = .data$case), inherit.aes = FALSE, nudge_x = 0.32, nudge_y = 0.32, size = 2.05, label.size = 0.18, label.padding = unit(0.035, "lines"), fill = "white") +
      facet_grid(rows = vars(direction_label)) +
      scale_fill_gradient2(
        low = "#355F9C", mid = "white", high = "#BD4A45", midpoint = 0,
        limits = c(-1500, 1500), oob = squish,
        name = "holdout log-lik\ntransfer − null"
      ) +
      labs(x = "model", y = "cell line") +
      guides(fill = guide_colorbar(barwidth = unit(1.9, "in"), barheight = unit(0.12, "in"))) +
      theme_fg3(7.3) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.7), axis.text.y = element_text(size = 6.3), legend.position = "bottom")

    best_draw <- function(dataset_id, model_id) {
      row <- assessment$fit_summary %>%
        filter(.data$dataset_id == .env$dataset_id, .data$model_id == .env$model_id) %>%
        slice(1)
      path <- file.path(
        release_root, "datasets", dataset_id, "optim", model_id, row$fit_id[[1]],
        paste0("optim_draws_", row$best_start[[1]], ".Rds")
      )
      draw <- readRDS(path)
      out <- as.numeric(draw[1, ])
      names(out) <- colnames(draw)
      out
    }
    no_sum_stan_case <- add_group_structure(readRDS(file.path(release_root, "datasets", "loo_exclude_sum_159_fuse", "stan_data.Rds")))
    make_case_overlay <- function(case, line_name, direction, model_id, transfer_id, null_id) {
      line_id <- as.integer(no_sum_stan_case$line_map[[line_name]])
      split <- get_directional_transfer_split(no_sum_stan_case, line_id = line_id, direction = direction)
      fit_specs <- list(
        Null = list(dataset = null_id, stan = add_group_structure(readRDS(file.path(release_root, "datasets", null_id, "stan_data.Rds")))),
        Transfer = list(dataset = transfer_id, stan = add_group_structure(readRDS(file.path(release_root, "datasets", transfer_id, "stan_data.Rds")))),
        `Full-data reference` = list(dataset = "loo_exclude_sum_159_fuse", stan = no_sum_stan_case)
      )
      sim <- bind_rows(lapply(names(fit_specs), function(fit_name) {
        spec <- fit_specs[[fit_name]]
        draw <- best_draw(spec$dataset, model_id)
        bind_rows(lapply(split$holdout_wells, function(well_idx) {
          simulate_well_from_index(
            model_id = model_id, model_r_path = get_model_r_path("gpath", "v1"),
            draw_vec = draw, stan_data = spec$stan, well_idx = well_idx,
            times = seq(0, max(no_sum_stan_case$t_grid), by = 1), line_name = line_name
          ) %>%
            filter(.data$variable %in% c("N", "R1")) %>%
            mutate(fit = fit_name, value = pmax(.data$value, 0))
        }))
      }))
      observed <- build_obs_df_from_stan_data(no_sum_stan_case, line_id = line_id, line_name = line_name) %>%
        filter(.data$well_idx %in% split$holdout_wells, .data$variable %in% c("N", "R1"))
      limits <- observed %>%
        group_by(.data$variable) %>%
        summarise(ymax = 1.6 * max(.data$value, na.rm = TRUE), .groups = "drop")
      sim <- sim %>% left_join(limits, by = "variable") %>% mutate(value = if_else(.data$value > .data$ymax, NA_real_, .data$value))
      score <- transfer %>%
        filter(as.character(.data$line_name) == .env$line_name, .data$direction == .env$direction, .data$model_id == .env$model_id) %>%
        slice(1) %>% pull(.data$transfer_minus_null)
      fit_cols <- c("Null" = "#6F6F6F", "Transfer" = "#2B6CB0", "Full-data reference" = "#C85247")
      fit_types <- c("Null" = "22", "Transfer" = "solid", "Full-data reference" = "42")
      state_plots <- lapply(c("N", "R1"), function(state_name) {
        ggplot(filter(sim, .data$variable == .env$state_name), aes(x = .data$time, y = .data$value, colour = .data$fit, linetype = .data$fit, group = interaction(.data$fit, .data$well_idx))) +
          geom_line(linewidth = 0.38) +
          geom_point(data = filter(observed, .data$variable == .env$state_name), aes(x = .data$time, y = .data$value), inherit.aes = FALSE, size = 0.42, alpha = 0.6) +
          facet_grid(rows = vars(G0), cols = vars(ploidy), scales = "free") +
          scale_colour_manual(values = fit_cols, name = "fit") +
          scale_linetype_manual(values = fit_types, name = "fit") +
          labs(x = "time (h)", y = NULL, title = if (state_name == "N") "N (alive cells)" else "R1 (glucose)") +
          theme_fg3(5.8) +
          theme(plot.title = element_text(face = "bold", size = 7), legend.position = "bottom", strip.text = element_text(size = 4.7), axis.text = element_text(size = 4.2))
      })
      (state_plots[[1]] | state_plots[[2]]) +
        plot_layout(guides = "collect") +
        plot_annotation(
          title = sprintf(
            "   Case %s: %s, %s, %s",
            case, line_name,
            ifelse(direction == "low_to_high", "low -> high", "high -> low"),
            unname(all_alias[[model_id]])
          ),
          subtitle = sprintf("   No-SUM-159-fuse fits; transfer-null LL=%+.0f.", score),
          theme = theme(plot.title = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 6.3))
        ) &
        theme(legend.position = "bottom")
    }
    list(
      schematic = schematic, heatmap = heatmap,
      case_c = make_case_overlay("C", "MCF10A", "low_to_high", "2R_2P_1W_C0_M0", "nosum_l2h_mcf10a", "nosum_l2h_null_mcf10a"),
      case_d = make_case_overlay("D", "SNU668", "high_to_low", "2R_2P_0W_C0_M1", "nosum_h2l_snu668", "nosum_h2l_null_snu668")
    )
  }
)

s8_s9 <- time_operation(
  "transformations_s8_s9",
  "Posterior growth-difference surfaces, fraction-positive map, and R1/R2 phase trajectories from all-lines predictions",
  {
    surface_rows <- list()
    phase_rows <- list()
    for (model_id in model_ids) {
      pred <- all_predictions[[model_id]]
      gs <- pred$growth_surface
      surface_array <- gs$high_minus_low_growth
      surface_med <- apply(surface_array, c(2, 3, 4), median, na.rm = TRUE)
      surface_pos <- apply(surface_array > 0, c(2, 3, 4), mean, na.rm = TRUE)
      line_names <- unique(pred$well_metadata$line_name)
      for (line_idx in seq_along(line_names)) {
        surface_rows[[length(surface_rows) + 1L]] <- expand_grid(
          glucose = gs$glucose, latent_R2 = gs$resource2
        ) %>%
          mutate(
            delta_growth = as.vector(surface_med[line_idx, , ]),
            fraction_positive = as.vector(surface_pos[line_idx, , ]),
            line_name = line_names[[line_idx]], model = unname(model_alias[model_id])
          )
      }
      if (model_id %in% model_ids[3:5]) {
        med <- apply(pred$state, c(2, 3, 4), median, na.rm = TRUE)
        meta <- pred$well_metadata
        phase_rows[[length(phase_rows) + 1L]] <- bind_rows(lapply(seq_len(dim(med)[1]), function(well_idx) {
          tibble(
            time = pred$time,
            R1 = med[well_idx, , "R1"],
            R2 = med[well_idx, , "R2"],
            line_name = meta$line_name[[well_idx]],
            ploidy = meta$ploidy_abs[[well_idx]],
            G0 = meta$G0[[well_idx]],
            model = unname(model_alias[model_id]),
            well_idx = well_idx
          )
        }))
      }
    }
    list(surface = bind_rows(surface_rows), phase = bind_rows(phase_rows))
  }
)
surface <- s8_s9$surface %>%
  mutate(
    model = factor(.data$model, levels = unname(model_alias)),
    line_name = factor(.data$line_name, levels = names(line_labels)),
    log10_glucose = log10(pmax(.data$glucose, 1e-4))
  )
phase <- s8_s9$phase %>%
  filter(.data$G0 %in% c(0.25, 1, 5), .data$R1 > 0) %>%
  mutate(
    model = factor(.data$model, levels = unname(model_alias[model_ids[3:5]])),
    line_name = factor(.data$line_name, levels = names(line_labels)),
    G0 = factor(.data$G0)
  ) %>%
  group_by(.data$line_name) %>%
  mutate(ploidy_label = if_else(.data$ploidy == max(.data$ploidy), "high ploidy", "low ploidy")) %>%
  ungroup() %>%
  mutate(ploidy_label = factor(.data$ploidy_label, levels = c("high ploidy", "low ploidy")))
panel_s8a <- ggplot(surface, aes(x = .data$log10_glucose, y = .data$latent_R2, fill = .data$delta_growth)) +
  geom_raster() +
  facet_grid(rows = vars(model), cols = vars(line_name)) +
  scale_fill_gradient2(low = "#3E6FB0", mid = "white", high = "#BE4D47", midpoint = 0, limits = c(-0.05, 0.05), oob = squish, name = "Δ growth") +
  labs(x = "log10 measured glucose", y = "latent R2") +
  theme_fg3(6.7) +
  theme(panel.spacing = unit(0.5, "mm"), strip.text.y = element_text(angle = 0), legend.position = "bottom")
surface_summary <- surface %>%
  group_by(.data$line_name, .data$log10_glucose, .data$latent_R2) %>%
  summarise(
    fraction_positive = mean(.data$delta_growth > 0, na.rm = TRUE),
    n_pareto_models = n_distinct(.data$model),
    .groups = "drop"
  )
surface_ann <- surface_summary %>%
  group_by(.data$line_name) %>%
  summarise(n = min(.data$n_pareto_models), .groups = "drop") %>%
  mutate(
    log10_glucose = min(surface_summary$log10_glucose) + 0.08 * diff(range(surface_summary$log10_glucose)),
    latent_R2 = 0.88,
    label = paste0("n=", .data$n)
  )
panel_s8b <- ggplot(surface_summary, aes(x = .data$log10_glucose, y = .data$latent_R2, fill = .data$fraction_positive)) +
  geom_raster() +
  geom_text(data = surface_ann, aes(x = .data$log10_glucose, y = .data$latent_R2, label = .data$label), inherit.aes = FALSE, size = 1.8, colour = "grey15") +
  facet_wrap(~line_name, nrow = 1) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  scale_fill_gradient(low = "#3666A8", high = "#B54542", limits = c(0, 1), name = "fraction\npositive") +
  labs(x = "log10 measured glucose", y = "latent R2") +
  theme_fg3(6.4) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5), panel.spacing = unit(0.6, "mm"), legend.position = "bottom")
panel_s9 <- ggplot(phase, aes(x = .data$R1, y = .data$R2, colour = .data$G0, group = interaction(.data$well_idx, .data$model))) +
  geom_path(linewidth = 0.45, arrow = arrow(type = "closed", length = unit(1.1, "mm")), alpha = 0.85) +
  facet_grid(rows = vars(model, ploidy_label), cols = vars(line_name), scales = "free") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  scale_colour_manual(values = c(`0.25` = "#E97C74", `1` = "#2FB271", `5` = "#6FA4E8"), name = "G0 (mM)") +
  labs(x = "measured glucose R1 (log scale)", y = "latent R2") +
  theme_fg3(6.5) +
  theme(panel.spacing = unit(0.55, "mm"), strip.text.y = element_text(angle = 0), legend.position = "bottom")

figures <- list(
  figure_3 = (panel_f3a | panel_f3b) / panel_f3c + plot_layout(heights = c(0.72, 1.28), widths = c(1, 1)) + plot_annotation(tag_levels = "a"),
  figure_s5 = ((panel_s5a | panel_s5b) / wrap_elements(full = panel_s5c)) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "a"),
  figure_s6 = s6_figure,
  figure_s7 = (
    s7$schematic /
      (
        s7$heatmap |
          (
            wrap_elements(full = s7$case_c) /
              wrap_elements(full = s7$case_d)
          )
      )
  ) + plot_layout(heights = c(0.42, 1.58), widths = c(1, 1)) + plot_annotation(tag_levels = "a"),
  figure_s8 = panel_s8a / panel_s8b + plot_layout(heights = c(1.65, 0.55)) + plot_annotation(tag_levels = "a"),
  figure_s9 = panel_s9 + plot_annotation(tag_levels = "a")
)
dimensions <- list(
  figure_3 = c(7.0, 6.72),
  figure_s5 = c(7.0, 5.85),
  figure_s6 = c(7.0, 9.12),
  figure_s7 = c(7.0, 7.10),
  figure_s8 = c(7.0, 9.00),
  figure_s9 = c(7.0, 7.45)
)

time_operation(
  "rendering",
  "Six final composites rendered directly from in-memory plot objects at 360 dpi; F3A read directly from canonical immutable raster",
  {
    for (name in names(figures)) {
      ggsave(
        file.path(final_dir, paste0(name, ".png")),
        figures[[name]],
        width = dimensions[[name]][[1]], height = dimensions[[name]][[2]],
        units = "in", dpi = 360, bg = "white", limitsize = FALSE
      )
    }
    invisible(NULL)
  }
)

total_elapsed <- proc.time()[["elapsed"]] - total_start
timings[[length(timings) + 1L]] <- tibble(
  package_id = "FG3_mechanistic_diagnostics",
  operation = "total",
  elapsed_seconds = round(total_elapsed, 3),
  details = paste0(
    "Canonical modelling release red_a30_counts_20260722; F3A canonical raster basename=",
    basename(f3a_path), "; raster_sha256=", f3a_sha
  )
)
write.table(
  bind_rows(timings),
  file.path(timing_dir, "FG3_mechanistic_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
message(sprintf("FG3 complete in %.3f seconds", total_elapsed))
