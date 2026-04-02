#!/usr/bin/env Rscript

library(dplyr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

TRANSFER_NUTS_ROOT <- if (length(args) >= 1) args[1] else file.path("data", "gpath_transfer_cv_nuts", "gpath")
STAN_DATA_PATH <- if (length(args) >= 2) args[2] else "data/stan_ready_data.Rds"
OUTPUT_ROOT <- if (length(args) >= 3) args[3] else file.path("data", "gpath_transfer_cv_nuts", "report")
MAX_PARAM_DRAWS <- if (length(args) >= 4) as.integer(args[4]) else 1000L
MAX_GROWTH_DRAWS <- if (length(args) >= 5) as.integer(args[5]) else 400L
N_CORES <- if (length(args) >= 6) as.integer(args[6]) else 4L
MAX_SURFACE_DRAWS <- if (length(args) >= 7) as.integer(args[7]) else min(MAX_GROWTH_DRAWS, 100L)

source("R/project_paths.R")
source(get_model_r_path("gpath", "v1"))
source("R/gpath_run_utils.R")
source("R/posterior_io_utils.R")
source("R/nuts_consistency_utils.R")

cat(sprintf(">>> summarize_gpath_nuts_consistency.R cwd: %s\n", getwd()))
cat(sprintf(">>> Using %d core(s)\n", N_CORES))

stan_data_path <- resolve_stan_data_path(STAN_DATA_PATH)
stan_data <- readRDS(stan_data_path)

model_dirs <- list.dirs(TRANSFER_NUTS_ROOT, recursive = FALSE, full.names = TRUE)
model_ids <- basename(model_dirs)
model_ids <- model_ids[dir.exists(file.path(TRANSFER_NUTS_ROOT, model_ids, "low_to_high", "oracle"))]

if (!length(model_ids)) {
  stop(sprintf("No oracle NUTS model directories found under %s", TRANSFER_NUTS_ROOT))
}

cat(sprintf(">>> Found %d oracle NUTS model directories\n", length(model_ids)))

focus_parameters <- c("ae[1]", "ah[1]", "Y_R[1]", "m")
glucose_grid <- exp(seq(log(0.0001), log(25), length.out = 40L))
resource2_grid <- seq(0, 1, length.out = length(glucose_grid))

summarize_one_model_id <- function(model_id) {
  cat(sprintf(">>> [%s] starting\n", model_id))

  qc <- collect_oracle_nuts_qc(
    model_ids = model_id,
    focus_parameters = focus_parameters,
    transfer_nuts_root = TRANSFER_NUTS_ROOT
  )
  focus_draw_df <- collect_oracle_focus_parameter_draws(
    model_ids = model_id,
    stan_data = stan_data,
    focus_parameters = focus_parameters,
    transfer_nuts_root = TRANSFER_NUTS_ROOT,
    max_draws_per_model = MAX_PARAM_DRAWS
  )
  focus_raw_theta_draw_df <- collect_oracle_focus_raw_theta_ploidy_draws(
    model_ids = model_id,
    focus_parameters = focus_parameters,
    transfer_nuts_root = TRANSFER_NUTS_ROOT,
    max_draws_per_model = MAX_PARAM_DRAWS
  )
  focus_summary_df <- summarize_focus_parameter_posteriors(focus_draw_df)
  growth_draw_df <- collect_oracle_growth_curve_draws(
    model_ids = model_id,
    stan_data = stan_data,
    transfer_nuts_root = TRANSFER_NUTS_ROOT,
    glucose_grid = glucose_grid,
    max_draws_per_model = MAX_GROWTH_DRAWS,
    aux_resource_value = 1.0
  )
  growth_model_df <- summarize_model_growth_curves(growth_draw_df)
  growth_delta_surface_df <- summarize_oracle_growth_delta_surface(
    model_ids = model_id,
    stan_data = stan_data,
    transfer_nuts_root = TRANSFER_NUTS_ROOT,
    glucose_grid = glucose_grid,
    resource2_grid = resource2_grid,
    max_draws_per_model = MAX_SURFACE_DRAWS,
    extra_resource_value = 1.0
  )

  cat(sprintf(">>> [%s] done | qc=%d focus_draws=%d raw_theta=%d focus=%d growth=%d surface=%d\n",
              model_id, nrow(qc$summary), nrow(focus_draw_df), nrow(focus_raw_theta_draw_df), nrow(focus_summary_df), nrow(growth_model_df), nrow(growth_delta_surface_df)))

  list(
    qc_summary = qc$summary,
    focus_draws = focus_draw_df,
    focus_raw_theta_draws = focus_raw_theta_draw_df,
    focus_summary = focus_summary_df,
    growth_model = growth_model_df,
    growth_delta_surface = growth_delta_surface_df
  )
}

per_model_results <- if (N_CORES > 1L) {
  parallel::mclapply(model_ids, summarize_one_model_id, mc.cores = N_CORES)
} else {
  lapply(model_ids, summarize_one_model_id)
}

qc_summary_df <- bind_rows(lapply(per_model_results, `[[`, "qc_summary"))
focus_draws_df <- bind_rows(lapply(per_model_results, `[[`, "focus_draws"))
focus_raw_theta_draws_df <- bind_rows(lapply(per_model_results, `[[`, "focus_raw_theta_draws"))
focus_summary_df <- bind_rows(lapply(per_model_results, `[[`, "focus_summary"))
growth_model_df <- bind_rows(lapply(per_model_results, `[[`, "growth_model"))
growth_delta_surface_df <- bind_rows(lapply(per_model_results, `[[`, "growth_delta_surface"))
growth_pooled_df <- summarize_pooled_growth_curves(growth_model_df)
growth_consistency_df <- summarize_growth_curve_consistency(growth_model_df)

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

saveRDS(qc_summary_df, file.path(OUTPUT_ROOT, "nuts_qc_summary.Rds"))
saveRDS(focus_draws_df, file.path(OUTPUT_ROOT, "nuts_focus_draws.Rds"))
saveRDS(focus_raw_theta_draws_df, file.path(OUTPUT_ROOT, "nuts_focus_raw_theta_ploidy_draws.Rds"))
saveRDS(focus_summary_df, file.path(OUTPUT_ROOT, "nuts_focus_summary.Rds"))
saveRDS(growth_model_df, file.path(OUTPUT_ROOT, "nuts_growth_model_summary.Rds"))
saveRDS(growth_delta_surface_df, file.path(OUTPUT_ROOT, "nuts_growth_delta_surface_summary.Rds"))
saveRDS(growth_pooled_df, file.path(OUTPUT_ROOT, "nuts_growth_pooled_summary.Rds"))
saveRDS(growth_consistency_df, file.path(OUTPUT_ROOT, "nuts_growth_consistency_summary.Rds"))

cat(sprintf(">>> Wrote NUTS consistency summaries to %s\n", normalizePath(OUTPUT_ROOT, winslash = "/", mustWork = FALSE)))
cat(sprintf(">>> QC rows: %d | focus draws: %d | raw-theta draws: %d | focus rows: %d | growth rows: %d | surface rows: %d\n",
            nrow(qc_summary_df), nrow(focus_draws_df), nrow(focus_raw_theta_draws_df), nrow(focus_summary_df), nrow(growth_model_df), nrow(growth_delta_surface_df)))
