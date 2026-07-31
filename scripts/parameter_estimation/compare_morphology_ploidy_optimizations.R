#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("morphology release config path required", call. = FALSE)
source("scripts/parameter_estimation/common.R")
source("R/optim_utils.R")

best_converged_diagnostics <- function(summary, assessment_path) {
  release_root <- dirname(dirname(dirname(assessment_path)))
  plan <- utils::read.delim(
    file.path(release_root, "manifests", "fit_plan.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rows <- lapply(seq_len(nrow(summary)), function(i) {
    fit <- summary[i, , drop = FALSE]
    optim_dir <- plan$optim_dir[match(fit$fit_id, plan$fit_id)]
    lp <- as.numeric(readRDS(file.path(optim_dir, "optim_lp_all.Rds")))
    rc <- as.numeric(unlist(readRDS(file.path(optim_dir, "optim_rc_all.Rds")), use.names = FALSE))
    candidates <- which(is.finite(lp) & rc == 0)
    if (!length(candidates)) {
      return(data.frame(
        fit_id = fit$fit_id, rc0_best_start = NA_integer_,
        rc0_best_lp = NA_real_, rc0_log_lik = NA_real_
      ))
    }
    best <- candidates[[which.max(lp[candidates])]]
    draws <- readRDS(file.path(optim_dir, "optim_draws_all.Rds"))
    draw <- extract_draw_vector(draws[[best]])
    data.frame(
      fit_id = fit$fit_id,
      rc0_best_start = best,
      rc0_best_lp = lp[[best]],
      rc0_log_lik = if ("log_lik" %in% names(draw)) as.numeric(draw[["log_lik"]]) else NA_real_
    )
  })
  do.call(rbind, rows)
}

cfg <- pe_read_config(args[[1]])
root <- pe_release_root(cfg)
reference_path <- if (length(args) >= 2L) args[[2]] else
  "data/modelling/gpath_v1/red_a30_counts_20260722/derived/optimization/assessment.Rds"
output_dir <- if (length(args) >= 3L) args[[3]] else
  "data/report_exports/morphology_vs_ploidy_five_line_20260729"

current_path <- file.path(root, "derived", "optimization", "assessment.Rds")
if (!file.exists(current_path)) {
  pe_fail("Morphology optimization assessment is missing; run derive-optim after fits complete: %s", current_path)
}
if (!file.exists(reference_path)) pe_fail("Reference ploidy assessment is missing: %s", reference_path)

current <- readRDS(current_path)$fit_summary
reference <- readRDS(reference_path)$fit_summary
current <- current[current$dataset_id %in% c("all_lines_cell_area", "all_lines_nuclear_area"), , drop = FALSE]
reference <- reference[reference$dataset_id == "all_lines", , drop = FALSE]
if (!nrow(current) || !nrow(reference)) pe_fail("Current or reference assessment has no requested rows")
current <- merge(
  current, best_converged_diagnostics(current, current_path),
  by = "fit_id", all.x = TRUE, sort = FALSE
)
reference <- merge(
  reference, best_converged_diagnostics(reference, reference_path),
  by = "fit_id", all.x = TRUE, sort = FALSE
)

keep <- c(
  "model_id", "complete", "best_start", "best_lp", "log_lik", "k",
  "n_starts", "n_finite_lp", "n_rc_zero",
  "rc0_best_start", "rc0_best_lp", "rc0_log_lik"
)
reference <- reference[, intersect(keep, names(reference)), drop = FALSE]
names(reference)[names(reference) != "model_id"] <- paste0("ploidy_", names(reference)[names(reference) != "model_id"])
joined <- merge(current, reference, by = "model_id", all.x = TRUE, sort = TRUE)
joined$metric <- ifelse(
  joined$dataset_id == "all_lines_cell_area", "cell_area",
  ifelse(joined$dataset_id == "all_lines_nuclear_area", "nuclear_area", joined$dataset_id)
)
joined$comparable <- joined$complete %in% TRUE & joined$ploidy_complete %in% TRUE &
  is.finite(joined$log_lik) & is.finite(joined$ploidy_log_lik)
joined$delta_log_lik_vs_ploidy <- joined$log_lik - joined$ploidy_log_lik
joined$likelihood_AIC <- -2 * joined$log_lik + 2 * joined$k
joined$ploidy_likelihood_AIC <- -2 * joined$ploidy_log_lik + 2 * joined$ploidy_k
joined$delta_likelihood_AIC_vs_ploidy <- joined$likelihood_AIC - joined$ploidy_likelihood_AIC
joined$rc0_likelihood_AIC <- -2 * joined$rc0_log_lik + 2 * joined$k
joined$ploidy_rc0_likelihood_AIC <- -2 * joined$ploidy_rc0_log_lik + 2 * joined$ploidy_k
joined$rc0_delta_likelihood_AIC_vs_ploidy <-
  joined$rc0_likelihood_AIC - joined$ploidy_rc0_likelihood_AIC

rank_input <- rbind(
  data.frame(
    model_id = reference$model_id, metric = "ploidy",
    log_lik = reference$ploidy_log_lik, k = reference$ploidy_k,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model_id = joined$model_id, metric = joined$metric,
    log_lik = joined$log_lik, k = joined$k,
    stringsAsFactors = FALSE
  )
)
rank_input$likelihood_AIC <- -2 * rank_input$log_lik + 2 * rank_input$k
rank_input <- rank_input[order(rank_input$model_id, rank_input$likelihood_AIC), , drop = FALSE]
rank_input$rank_within_model <- ave(
  rank_input$likelihood_AIC, rank_input$model_id,
  FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
)

rc0_rank_input <- rbind(
  data.frame(
    model_id = reference$model_id, metric = "ploidy",
    log_lik = reference$ploidy_rc0_log_lik, k = reference$ploidy_k,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model_id = joined$model_id, metric = joined$metric,
    log_lik = joined$rc0_log_lik, k = joined$k,
    stringsAsFactors = FALSE
  )
)
rc0_rank_input$likelihood_AIC <- -2 * rc0_rank_input$log_lik + 2 * rc0_rank_input$k
rc0_rank_input <- rc0_rank_input[
  order(rc0_rank_input$model_id, rc0_rank_input$likelihood_AIC),
  , drop = FALSE
]
rc0_rank_input$rank_within_model <- ave(
  rc0_rank_input$likelihood_AIC, rc0_rank_input$model_id,
  FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(joined, file.path(output_dir, "morphology_vs_ploidy_by_model.csv"), row.names = FALSE)
utils::write.csv(rank_input, file.path(output_dir, "metric_ranking_by_model.csv"), row.names = FALSE)
utils::write.csv(
  rc0_rank_input,
  file.path(output_dir, "metric_ranking_by_model_rc0_only.csv"),
  row.names = FALSE
)
pe_write_json(list(
  morphology_assessment_path = current_path,
  morphology_assessment_sha256 = pe_sha256(current_path),
  ploidy_assessment_path = reference_path,
  ploidy_assessment_sha256 = pe_sha256(reference_path),
  comparison_basis = "full-data generated log_lik with matched model_id and parameter count",
  warning = "Optimized lp__ is retained for diagnostics but is not used as the primary cross-metric comparison because covariate scaling changes effective prior regularization.",
  outputs = list(
    by_model = file.path(output_dir, "morphology_vs_ploidy_by_model.csv"),
    ranking = file.path(output_dir, "metric_ranking_by_model.csv"),
    ranking_rc0_only = file.path(output_dir, "metric_ranking_by_model_rc0_only.csv")
  ),
  generated_at = format(Sys.time(), usetz = TRUE)
), file.path(output_dir, "provenance.json"))

cat(sprintf("Wrote morphology/ploidy comparison outputs to %s\n", output_dir))
