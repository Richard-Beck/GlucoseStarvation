args <- commandArgs(trailingOnly = TRUE)

signature_panel_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else file.path("data", "model_free_ploidy", "signature_panel.Rds")
output_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path("data", "model_free_ploidy")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

source("R/model_free_ploidy_utils.R")

signature_panel <- readRDS(signature_panel_path)
transfer <- evaluate_transfer_predictions(signature_panel)
robustness <- analyze_model_free_outlier_robustness(signature_panel)

write.csv(transfer$predictions, file.path(output_dir, "transfer_predictions.csv"), row.names = FALSE)
write.csv(transfer$feature_summary, file.path(output_dir, "transfer_feature_summary.csv"), row.names = FALSE)
write.csv(transfer$case_summary, file.path(output_dir, "transfer_case_summary.csv"), row.names = FALSE)
write.csv(robustness$outlier_scores, file.path(output_dir, "outlier_scores.csv"), row.names = FALSE)
write.csv(robustness$exclusion_summary, file.path(output_dir, "outlier_exclusion_summary.csv"), row.names = FALSE)
write.csv(robustness$feature_exclusion_summary, file.path(output_dir, "outlier_feature_exclusion_summary.csv"), row.names = FALSE)
write.csv(robustness$recommended_exclusion, file.path(output_dir, "outlier_recommended_exclusion.csv"), row.names = FALSE)

saveRDS(transfer$predictions, file.path(output_dir, "transfer_predictions.Rds"))
saveRDS(transfer$feature_summary, file.path(output_dir, "transfer_feature_summary.Rds"))
saveRDS(transfer$case_summary, file.path(output_dir, "transfer_case_summary.Rds"))
saveRDS(robustness$outlier_scores, file.path(output_dir, "outlier_scores.Rds"))
saveRDS(robustness$outlier_feature_stats, file.path(output_dir, "outlier_feature_stats.Rds"))
saveRDS(robustness$baseline_overall, file.path(output_dir, "outlier_baseline_overall.Rds"))
saveRDS(robustness$exclusion_summary, file.path(output_dir, "outlier_exclusion_summary.Rds"))
saveRDS(robustness$feature_exclusion_summary, file.path(output_dir, "outlier_feature_exclusion_summary.Rds"))
saveRDS(robustness$recommended_exclusion, file.path(output_dir, "outlier_recommended_exclusion.Rds"))

cat(sprintf("Wrote model-free transfer outputs to %s\n", output_dir))
