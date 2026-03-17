args <- commandArgs(trailingOnly = TRUE)

stan_data_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else NULL
output_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path("data", "model_free_ploidy")
glucose_floor <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 0.1

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

source("R/model_free_ploidy_utils.R")

res <- build_feature_panel(
  stan_data_path = stan_data_path,
  glucose_floor = glucose_floor
)

effects <- compute_empirical_effects(res$feature_panel)

write.csv(res$feature_panel, file.path(output_dir, "feature_panel.csv"), row.names = FALSE)
write.csv(res$count_summary, file.path(output_dir, "count_summary.csv"), row.names = FALSE)
write.csv(res$glucose_summary, file.path(output_dir, "glucose_summary.csv"), row.names = FALSE)
write.csv(effects$paired_effects_long, file.path(output_dir, "paired_effects_long.csv"), row.names = FALSE)
write.csv(effects$effect_matrix, file.path(output_dir, "effect_matrix.csv"), row.names = FALSE)

saveRDS(res$feature_panel, file.path(output_dir, "feature_panel.Rds"))
saveRDS(res$count_summary, file.path(output_dir, "count_summary.Rds"))
saveRDS(res$glucose_summary, file.path(output_dir, "glucose_summary.Rds"))
saveRDS(effects$paired_effects_wide, file.path(output_dir, "paired_effects_wide.Rds"))
saveRDS(effects$paired_effects_long, file.path(output_dir, "paired_effects_long.Rds"))
saveRDS(effects$effect_matrix, file.path(output_dir, "effect_matrix.Rds"))
saveRDS(effects$pca_scores, file.path(output_dir, "effect_pca_scores.Rds"))
saveRDS(effects$pca_loadings, file.path(output_dir, "effect_pca_loadings.Rds"))

cat(sprintf("Wrote model-free feature outputs to %s\n", output_dir))
