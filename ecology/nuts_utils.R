library(posterior)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

extract_theta_ploidy <- function(nutsDir, run_id, chains = 1:4) {
  
  # 1. Parse dimensions using your sourced helper
  run_params <- parse_run_id(run_id)
  R <- run_params$R
  P <- run_params$P
  W <- run_params$W
  C <- run_params$C
  M <- run_params$M
  # 2. Get exact prior scales from your Stan config generator
  config <- generate_stan_config(
    R = R, 
    P = P, 
    W = W, 
    strict_spec = (C == 1L),
    M=M,
    base_priors = base_priors
  )
  p_scales <- config$prior_scales
  
  param_labels <- get_param_names(R,P,W,C,M)
  
  # 4. Load draws and subset
  draw_files <- file.path(nutsDir, sprintf("nuts_draws_%d.Rds", chains))
  draws_raw <- reduce(map(draw_files, readRDS), bind_draws, along="chain")
  draws_sub <- subset_draws(draws_raw, variable = c("raw_theta_ploidy", "lp__"))
  
  # 5. Filter out stuck chains
  lp_df <- as_draws_df(subset_draws(draws_sub, variable = "lp__")) %>%
    group_by(.chain) %>%
    summarise(q10 = quantile(lp__, 0.10), q90 = quantile(lp__, 0.90), .groups="drop")
  
  best_q10 <- max(lp_df$q10)
  stuck_chains <- filter(lp_df, q90 < best_q10)$.chain
  good_chains <- setdiff(chains, stuck_chains)
  
  if (length(good_chains) == 0) stop(sprintf("[%s] All chains failed lp__ overlap check.", run_id))
  
  # 6. Extract clean draws, merge with scales, and calculate fold change
  clean_ploidy_draws <- as_draws_df(subset_draws(draws_sub, variable = "raw_theta_ploidy", chain = good_chains))
  
  df_long <- clean_ploidy_draws %>%
    pivot_longer(
      cols = starts_with("raw_theta_ploidy["),
      names_to = "raw_name",
      values_to = "parval"
    ) %>%
    mutate(
      model_id = run_id,
      chain_id = .chain,
      idx = as.integer(str_extract(raw_name, "(?<=\\[)\\d+(?=\\])")),
      parname = param_labels[idx],
      prior_scale = p_scales[idx],
      fold_change = exp(parval * prior_scale),
      scaling_exponent = log2(fold_change)
    ) %>%
    select(model_id, chain_id, parname, parval, fold_change,scaling_exponent)
  
  return(df_long)
}

analyze_nuts_run <- function(nutsDir, run_id, chains=1:4, cor_threshold=0.8) {
  
  parts <- strsplit(run_id, "_")[[1]]
  R <- as.integer(sub("R", "", parts[1]))
  P <- as.integer(sub("P", "", parts[2]))
  W <- as.integer(sub("W", "", parts[3]))
  
  param_labels <- c(sprintf("ae_%d", 1:R), sprintf("ah_%d", 1:R))
  for(c in 1:R) for(r in 1:P) param_labels <- c(param_labels, sprintf("Y_mat_%d_%d", r, c))
  if (W > 0) {
    for(c in 1:R) for(r in 1:W) param_labels <- c(param_labels, sprintf("K_mat_%d_%d", r, c))
    param_labels <- c(param_labels, sprintf("tox_%d", 1:W))
  }
  param_labels <- c(param_labels, "m")
  
  draw_files <- file.path(nutsDir, sprintf("nuts_draws_%d.Rds", chains))
  diag_files <- file.path(nutsDir, sprintf("nuts_diagnostics_%d.Rds", chains))
  
  draws_raw <- reduce(map(draw_files, readRDS), bind_draws, along="chain")
  diags_raw <- reduce(map(diag_files, readRDS), bind_draws, along="chain")
  
  vars <- variables(draws_raw)
  
  base_vars <- c("raw","sigma_base","sigma_rel","G1_0","lp__")
  
  target_vars <- vars[
    Reduce(`|`, lapply(base_vars, function(b) startsWith(vars, b)))
  ]
  draws_combined <- subset_draws(draws_raw, variable = target_vars)
  
  lp_df <- as_draws_df(subset_draws(draws_combined, variable = "lp__")) %>%
    group_by(.chain) %>%
    summarise(q10 = quantile(lp__, 0.10), q90 = quantile(lp__, 0.90), med = median(lp__), .groups="drop")
  
  best_chain_q10 <- max(lp_df$q10)
  stuck_chains <- filter(lp_df, q90 < best_chain_q10)$.chain
  good_chains <- setdiff(chains, stuck_chains)
  
  #if (length(good_chains) <=1) stop("All chains failed lp__ overlap check. Severe convergence failure.")
  
  draws_clean <- subset_draws(draws_combined, chain = good_chains)
  diags_clean <- subset_draws(diags_raw, chain = good_chains)
  
  conv_stats <- summarise_draws(draws_clean, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  flagged_params <- filter(conv_stats, rhat > 1.01 | ess_bulk < 400) %>% arrange(ess_bulk)
  
  diag_stats <- summarise_draws(diags_clean, total=~sum(.x), max=~max(.x), mean=~mean(.x))
  key_diags <- filter(diag_stats, grepl("divergent__|treedepth__|accept_stat__", variable))
  
  annotate_param <- function(p) {
    if (grepl("\\[\\d+,\\s*\\d+\\]", p)) {
      base <- sub("\\[.*", "", p)
      idx1 <- as.integer(sub(".*\\[(\\d+),.*", "\\1", p))
      idx2 <- as.integer(sub(".*,\\s*(\\d+)\\].*", "\\1", p))
      if (base == "raw_theta_line") return(sprintf("%s (Line %d)", param_labels[idx1], idx2))
    }
    if (grepl("\\[\\d+\\]", p)) {
      base <- sub("\\[.*", "", p)
      idx1 <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", p))
      if (base == "raw_theta_ploidy") return(sprintf("%s (Ploidy Effect)", param_labels[idx1]))
      if (base == "raw_N0") return(sprintf("N0 (Group %d)", idx1))
    }
    return(p)
  }
  
  cor_df <- as.data.frame(as.table(cor(as_draws_matrix(draws_clean)))) %>%
    rename(param1 = Var1, param2 = Var2, correlation = Freq) %>%
    mutate(param1 = as.character(param1), param2 = as.character(param2), abs_cor = abs(correlation)) %>%
    filter(param1 < param2, abs_cor > cor_threshold)
  
  if (nrow(cor_df) > 0) {
    cor_summary <- cor_df %>%
      rowwise() %>%
      mutate(
        p1_desc = annotate_param(param1),
        p2_desc = annotate_param(param2),
        base_p1 = trimws(str_remove(p1_desc, "\\(.*\\)")),
        base_p2 = trimws(str_remove(p2_desc, "\\(.*\\)")),
        pair_name = paste(sort(c(base_p1, base_p2)), collapse = " - "),
        is_l1 = str_detect(p1_desc, "Line"),
        is_l2 = str_detect(p2_desc, "Line"),
        lid1 = str_extract(p1_desc, "(?<=Line )\\d+"),
        lid2 = str_extract(p2_desc, "(?<=Line )\\d+"),
        rel_type = case_when(
          is_l1 & is_l2 & (lid1 == lid2) ~ "Within-Line",
          is_l1 & is_l2 & (lid1 != lid2) ~ "Between-Line",
          TRUE ~ "Global/Other"
        )
      ) %>%
      ungroup() %>%
      group_by(pair_name, rel_type) %>%
      summarize(count = n(), avg_abs_cor = round(mean(abs_cor), 3), .groups = "drop") %>%
      arrange(desc(count), desc(avg_abs_cor))
  } else {
    cor_summary <- data.frame()
  }
  
  cat("\n=== NUTS DIAGNOSTIC REPORT ===\n")
  cat("\n[1] CHAIN CONVERGENCE (lp__ overlap)\n")
  print(lp_df)
  if (length(stuck_chains) > 0) {
    cat(sprintf("-> WARNING: Dropped chains %s (stuck in suboptimal modes).\n", paste(stuck_chains, collapse=", ")))
  } else {
    cat("-> PASS: All chains overlap.\n")
  }
  
  cat(sprintf("\n[2] SAMPLER EFFICIENCY (Chains %s)\n", paste(good_chains, collapse=", ")))
  print(select(key_diags, variable, total, max, mean))
  
  cat(sprintf("\n[3] PARAMETER CONVERGENCE (Rhat > 1.01 | Bulk ESS < 400): %d / %d failed\n", nrow(flagged_params), nrow(conv_stats)))
  if (nrow(flagged_params) > 0) print(head(select(flagged_params, variable, rhat, ess_bulk, ess_tail), 10))
  
  cat(sprintf("\n[4] STRUCTURAL COLLINEARITY (|r| > %s)\n", cor_threshold))
  if (nrow(cor_summary) > 0) print(head(cor_summary, 15)) else cat("-> PASS: No highly correlated parameter groups.\n")
  
  invisible(list(
    clean_draws = draws_clean,
    stuck_chains = stuck_chains, 
    lp_stats = lp_df, 
    diagnostics = key_diags, 
    flagged = flagged_params, 
    correlations = cor_summary
  ))
}