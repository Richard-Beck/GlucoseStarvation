args <- commandArgs(trailingOnly = TRUE)

feature_panel_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else file.path("data", "model_free_ploidy", "feature_panel.Rds")
output_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path("data", "model_free_ploidy")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

source("R/model_free_ploidy_utils.R")

feature_panel <- readRDS(feature_panel_path)
transfer <- evaluate_transfer_predictions(feature_panel)

write.csv(transfer$predictions, file.path(output_dir, "transfer_predictions.csv"), row.names = FALSE)
write.csv(transfer$feature_summary, file.path(output_dir, "transfer_feature_summary.csv"), row.names = FALSE)
write.csv(transfer$case_summary, file.path(output_dir, "transfer_case_summary.csv"), row.names = FALSE)

saveRDS(transfer$predictions, file.path(output_dir, "transfer_predictions.Rds"))
saveRDS(transfer$feature_summary, file.path(output_dir, "transfer_feature_summary.Rds"))
saveRDS(transfer$case_summary, file.path(output_dir, "transfer_case_summary.Rds"))

cat(sprintf("Wrote model-free transfer outputs to %s\n", output_dir))
