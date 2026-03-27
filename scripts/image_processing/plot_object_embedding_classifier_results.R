#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: plot_object_embedding_classifier_results.R <results_dir> [top_n]", call. = FALSE)
}

results_dir <- normalizePath(args[[1]], mustWork = TRUE)
top_n <- if (length(args) >= 2) as.integer(args[[2]]) else 20L
if (is.na(top_n) || top_n <= 0L) {
  stop("top_n must be a positive integer.", call. = FALSE)
}

metrics_path <- file.path(results_dir, "metrics.json")
pred_path <- file.path(results_dir, "heldout_predictions.csv")
coef_path <- file.path(results_dir, "coefficients.csv")

if (!file.exists(metrics_path)) stop("Missing metrics.json", call. = FALSE)
if (!file.exists(pred_path)) stop("Missing heldout_predictions.csv", call. = FALSE)
if (!file.exists(coef_path)) stop("Missing coefficients.csv", call. = FALSE)

metrics <- fromJSON(metrics_path, simplifyVector = TRUE)
pred <- read.csv(pred_path, stringsAsFactors = FALSE)
coef <- read.csv(coef_path, stringsAsFactors = FALSE)

test_pred <- pred[pred$split == "test", , drop = FALSE]
if (nrow(test_pred) == 0L) {
  stop("No test rows found in heldout_predictions.csv", call. = FALSE)
}

test_metrics <- metrics$test_metrics
cm <- matrix(
  c(test_metrics$tn, test_metrics$fp, test_metrics$fn, test_metrics$tp),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("true_0", "true_1"), c("pred_0", "pred_1"))
)

cm_df <- expand.grid(
  truth = rownames(cm),
  pred = colnames(cm),
  stringsAsFactors = FALSE
)
cm_df$count <- as.vector(cm)

score_grid <- data.frame(threshold = seq(0, 1, by = 0.01))
score_grid$tp <- vapply(score_grid$threshold, function(th) sum(test_pred$y_true == 1L & test_pred$y_prob >= th), numeric(1))
score_grid$fp <- vapply(score_grid$threshold, function(th) sum(test_pred$y_true == 0L & test_pred$y_prob >= th), numeric(1))
score_grid$tn <- vapply(score_grid$threshold, function(th) sum(test_pred$y_true == 0L & test_pred$y_prob < th), numeric(1))
score_grid$fn <- vapply(score_grid$threshold, function(th) sum(test_pred$y_true == 1L & test_pred$y_prob < th), numeric(1))

score_grid$tpr <- ifelse((score_grid$tp + score_grid$fn) > 0, score_grid$tp / (score_grid$tp + score_grid$fn), NA_real_)
score_grid$fpr <- ifelse((score_grid$fp + score_grid$tn) > 0, score_grid$fp / (score_grid$fp + score_grid$tn), NA_real_)
score_grid$precision <- ifelse((score_grid$tp + score_grid$fp) > 0, score_grid$tp / (score_grid$tp + score_grid$fp), 1)
score_grid$recall <- score_grid$tpr

roc_df <- unique(score_grid[, c("fpr", "tpr")])
roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), , drop = FALSE]

pr_df <- unique(score_grid[, c("recall", "precision")])
pr_df <- pr_df[order(pr_df$recall, pr_df$precision), , drop = FALSE]

top_coef <- coef[seq_len(min(nrow(coef), top_n)), , drop = FALSE]
top_coef$feature <- factor(top_coef$feature, levels = rev(top_coef$feature))
top_coef$direction <- ifelse(top_coef$coefficient >= 0, "positive", "negative")

p_cm <- ggplot(cm_df, aes(x = pred, y = truth, fill = count)) +
  geom_tile() +
  geom_text(aes(label = count), size = 5) +
  scale_fill_gradient(low = "white", high = "#2166ac") +
  labs(title = "Held-Out Confusion Matrix", x = NULL, y = NULL, fill = "count") +
  theme_minimal(base_size = 12)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_path(color = "#1b9e77", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Held-Out ROC",
    x = "False positive rate",
    y = "True positive rate",
    subtitle = paste0("ROC AUC: ", format(round(test_metrics$roc_auc, 4), nsmall = 4))
  ) +
  theme_minimal(base_size = 12)

p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_path(color = "#d95f02", linewidth = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Held-Out Precision-Recall",
    x = "Recall",
    y = "Precision",
    subtitle = paste0("Average precision: ", format(round(test_metrics$average_precision, 4), nsmall = 4))
  ) +
  theme_minimal(base_size = 12)

p_coef <- ggplot(top_coef, aes(x = feature, y = coefficient, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(negative = "#d95f02", positive = "#1b9e77")) +
  labs(title = paste0("Top ", nrow(top_coef), " Coefficients"), x = NULL, y = "coefficient", fill = NULL) +
  theme_minimal(base_size = 12)

png(file.path(results_dir, "summary.png"), width = 1800, height = 1400, res = 180)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
print(p_cm, vp = grid::viewport())
print(p_roc, vp = grid::viewport())
print(p_pr, vp = grid::viewport())
print(p_coef, vp = grid::viewport())
dev.off()

summary_lines <- c(
  "Object Embedding Classifier Summary",
  paste0("Input CSV: ", metrics$input_csv),
  paste0("Split type: ", metrics$split_type),
  paste0("Rows used: ", metrics$n_rows_used),
  paste0("Features: ", metrics$n_features),
  paste0("Held-out n: ", test_metrics$n),
  paste0("Held-out prevalence: ", format(round(test_metrics$prevalence, 4), nsmall = 4)),
  paste0("Held-out accuracy: ", format(round(test_metrics$accuracy, 4), nsmall = 4)),
  paste0("Held-out balanced_accuracy: ", format(round(test_metrics$balanced_accuracy, 4), nsmall = 4)),
  paste0("Held-out precision: ", format(round(test_metrics$precision, 4), nsmall = 4)),
  paste0("Held-out recall: ", format(round(test_metrics$recall, 4), nsmall = 4)),
  paste0("Held-out f1: ", format(round(test_metrics$f1, 4), nsmall = 4)),
  paste0("Held-out roc_auc: ", format(round(test_metrics$roc_auc, 4), nsmall = 4)),
  paste0("Held-out average_precision: ", format(round(test_metrics$average_precision, 4), nsmall = 4)),
  "",
  "Top coefficients by absolute value:"
)
for (i in seq_len(nrow(top_coef))) {
  summary_lines <- c(
    summary_lines,
    paste0(
      top_coef$feature[[i]],
      ": coef=",
      format(round(top_coef$coefficient[[i]], 6), nsmall = 6),
      ", abs=",
      format(round(top_coef$abs_coefficient[[i]], 6), nsmall = 6)
    )
  )
}
writeLines(summary_lines, con = file.path(results_dir, "summary.txt"))

cat("Wrote figure: ", normalizePath(file.path(results_dir, "summary.png")), "\n", sep = "")
cat("Wrote summary: ", normalizePath(file.path(results_dir, "summary.txt")), "\n", sep = "")
