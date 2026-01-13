library(cmdstanr)
library(jsonlite)
library(data.table)

plot_latent_like <- function(fit, stan_data, MODEL_NAME = "model_B", draw = 1, ploidy_cut = 0.01) {
  suppressPackageStartupMessages({
    library(posterior)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
  })
  
  # ---- line id -> name map (best-effort) ----
  id_map <- NULL
  if (!is.null(stan_data$line_map)) id_map <- setNames(names(stan_data$line_map), stan_data$line_map)
  if (is.null(id_map)) {
    lids <- sort(unique(stan_data$line_id))
    id_map <- setNames(as.character(lids), as.character(lids))
  }
  
  # ---- extract y_sim as matrix, pick ONE draw (row) ----
  mu_mat <- posterior::as_draws_matrix(fit$draws("y_sim"))
  if (draw < 1 || draw > nrow(mu_mat)) stop("draw out of range (1..", nrow(mu_mat), ")")
  
  y_row <- as.numeric(mu_mat[draw, ])
  var_names <- colnames(mu_mat)
  
  # ---- parse indices from y_sim[w,t,s] colnames ----
  idx_df <- tibble(variable = var_names, value = y_row) %>%
    mutate(clean = gsub("^y_sim\\[|\\]$", "", variable)) %>%
    separate(clean, c("w","t","s"), sep = ",", convert = TRUE) %>%
    mutate(
      well_idx = as.integer(w),
      t_idx    = as.integer(t),
      s_idx    = as.integer(s),
      time     = stan_data$t_grid[t_idx],
      type     = case_when(s_idx == 1 ~ "NL", s_idx == 2 ~ "ND", s_idx == 3 ~ "G", TRUE ~ NA_character_)
    ) %>%
    filter(!is.na(type)) %>%
    select(well_idx, time, type, value)
  
  # ---- metadata ----
  meta_df <- data.frame(
    well_idx  = seq_len(stan_data$N_wells),
    line_id   = stan_data$line_id,
    line_name = id_map[as.character(stan_data$line_id)],
    metric    = stan_data$ploidy_metric,
    G0        = stan_data$G0_per_well,
    exp_id    = stan_data$exp_id
  )
  meta_df$ploidy_lbl <- ifelse(meta_df$metric > ploidy_cut, "High Ploidy", "Baseline Ploidy")
  meta_df$G0_lbl <- factor(paste0(meta_df$G0, " mM"), levels = paste0(sort(unique(meta_df$G0)), " mM"))
  
  sim_all <- idx_df %>% left_join(meta_df, by = "well_idx")
  
  # ---- observed ----
  obs_counts <- data.frame(
    well_idx = stan_data$well_idx_count,
    time     = stan_data$t_grid[stan_data$grid_idx_count],
    value    = stan_data$N_obs,
    type     = "NL"
  ) %>%
    bind_rows(data.frame(
      well_idx = stan_data$well_idx_count,
      time     = stan_data$t_grid[stan_data$grid_idx_count],
      value    = stan_data$D_obs,
      type     = "ND"
    )) %>%
    left_join(meta_df %>% select(well_idx, exp_id), by = "well_idx")
  
  obs_gluc <- data.frame(
    well_idx = stan_data$well_idx_gluc,
    time     = stan_data$t_grid[stan_data$grid_idx_gluc],
    lum      = stan_data$lum_obs,
    dilution = stan_data$dilution
  ) %>%
    left_join(meta_df %>% select(well_idx, exp_id), by = "well_idx") %>%
    mutate(
      a = stan_data$calib_a_fixed[exp_id],
      b = stan_data$calib_b_fixed[exp_id],
      value = pmax(0, (lum - b) / (a * dilution)),
      type  = "G"
    ) %>%
    select(well_idx, time, value, type, exp_id)
  
  obs_all <- bind_rows(obs_counts, obs_gluc) %>% left_join(meta_df, by = c("well_idx", "exp_id"))
  
  # ---- plotting ----
  cols <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")
  
  unique_lines <- sort(unique(meta_df$line_id))
  for (lid in unique_lines) {
    lname <- id_map[as.character(lid)]
    d_sim <- sim_all %>% filter(line_id == lid)
    d_obs <- obs_all %>% filter(line_id == lid)
    if (nrow(d_sim) == 0) next
    
    p1 <- ggplot() +
      geom_line(
        data = d_sim %>% filter(type != "G"),
        aes(x = time, y = value, color = type, group = interaction(well_idx, type)),
        linewidth = 0.8
      ) +
      geom_point(
        data = d_obs %>% filter(type != "G"),
        aes(x = time, y = value, color = type, shape = factor(exp_id)),
        size = 2, alpha = 0.7
      ) +
      facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
      scale_y_continuous(labels = scales::comma) +
      scale_color_manual(values = cols) +
      labs(title = paste0("Latent Cell Counts | Line ", lname), x = "Time (h)", y = "Count", shape = "Exp ID") +
      theme_bw() +
      theme(legend.position = "bottom", strip.background = element_rect(fill = "#f0f0f0"))
    
    p2 <- ggplot() +
      geom_line(
        data = d_sim %>% filter(type == "G"),
        aes(x = time, y = value, color = type, group = interaction(well_idx, type)),
        linewidth = 0.8
      ) +
      geom_point(
        data = d_obs %>% filter(type == "G"),
        aes(x = time, y = value, color = type, shape = factor(exp_id)),
        size = 2, alpha = 0.7
      ) +
      facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
      scale_color_manual(values = cols) +
      labs(title = paste0("Latent Glucose | Line ", lname), x = "Time (h)", y = "Conc (mM)", shape = "Exp ID") +
      theme_bw() +
      theme(legend.position = "bottom", strip.background = element_rect(fill = "#f0f0f0"))
    
    print((p1 | p2) +
            plot_layout(guides = "collect") +
            plot_annotation(
              title = paste0("Latent Trajectories: Cell Line ", lid, " (", MODEL_NAME, ")"),
              subtitle = paste0("Single draw (row): ", draw),
              theme = theme(plot.title = element_text(size = 16, face = "bold"))
            )
    )
  }
  
  invisible(TRUE)
}

# ---- 1. Compile ----
library(cmdstanr)
library(jsonlite)

# ---- settings ----
MODEL_NAME <- "model_B"   # or "model_B0"
ALG <- "bfgs"             # "lbfgs" or "bfgs"
SEED <- 123
ITERS <- 2000
THREADS <- 3

OUTDIR <- file.path("results", "opt_csv")
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
BASE <- sprintf("opt_%s_%s_seed%d", MODEL_NAME, ALG, SEED)
CSV  <- file.path(OUTDIR, paste0(BASE, ".csv"))

# ---- data prep (mirrors your runner) ----
stan_data <- readRDS("data/stan_ready_data.Rds")
config <- read_json(paste0("MCMC/", MODEL_NAME, ".json"), simplifyVector = TRUE)

max_N <- max(stan_data$N_obs, na.rm = TRUE)
config$prior_means[1] <- log(max_N * 1.5)
config$prior_sds[1]   <- 0.5

# Patch in config data
if (MODEL_NAME == "model_B") {
  stan_data$lower_b  <- config$lower
  stan_data$upper_b  <- config$upper
  stan_data$N_params <- length(config$lower)
}

# Add derived data required for Inits or Model
stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0
stan_data$calc_sim       <- 1
stan_data$structure_mode <- 0

# [CRITICAL FIX] Calculate log bounds in R for the init generator
stan_data$log_lower <- log(stan_data$lower_b)
stan_data$log_upper <- log(stan_data$upper_b)

# ---- compile ----
stan_file <- paste0("MCMC/", MODEL_NAME, ".stan")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# ---- 2. Robust Init Generator (CORRECTED) ----
gen_inits <- function(d) {
  # Calculate log bounds locally 
  log_lo <- log(d$lower_b)
  log_up <- log(d$upper_b)
  
  list(
    # FIX: Do NOT exponentiate. mu_global is defined in log-space in the Stan file.
    mu_global = log_lo + (log_up - log_lo) * runif(length(log_lo), 0.25, 0.75),
    
    # Variances (Small start)
    sigma_line    = rep(0.1, 10),
    sigma_beta    = rep(0.1, 10),
    sigma_IC      = rep(0.1, 2),
    sigma_beta_IC = rep(0.1, 2),
    
    # Random effects (Zeros)
    z_line    = matrix(0, 10, d$N_lines),
    z_beta    = matrix(0, 10, d$N_lines),
    z_IC      = matrix(0, 2,  d$N_lines),
    z_beta_IC = matrix(0, 2,  d$N_lines),
    
    # Slopes and Raw ICs
    beta_high = rep(0, 10),
    beta_IC   = rep(0, 2),
    mu_IC_raw = rep(0, 2),
    
    # Noise
    phi_total = 1.0,
    phi_frac  = 10.0
  )
}

# Generate
init_list <- gen_inits(stan_data)

fit_sim <- mod$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE,
  init = list(init_list),   # or a function returning your parameter list
  seed = 123,
  threads_per_chain = 3,
  refresh = 1
)

# y_sim will be generated from exactly init_list's parameters
draws <- fit_sim$draws("y_sim")
stan_data_sim <- stan_data
stan_data_sim$mode <- 2
stan_data_sim$calc_sim <- 1
plot_latent_like(fit_sim, stan_data_sim)

# ---- 3. Run Optimize ----
fit_opt <- mod$optimize(
  data = stan_data,
  algorithm = "bfgs",
  init = list(init_list), 
  iter = 100,
  refresh = 50,
  threads = 3
)
# ---- 4. Check Output ----
if (!is.null(fit_opt$return_codes()) && fit_opt$return_codes() == 0) {
  print("Success!")
  print(head(fit_opt$mle()))
} else {
  print("Optimization failed. Checking CSV for final state...")
  csv_files <- fit_opt$output_files()
  if (length(csv_files) > 0 && file.exists(csv_files)) {
    headers <- names(fread(csv_files, skip = "lp__", nrows = 0))
    drop_cols <- grep("^y_sim|__$", headers, value = TRUE)
    res_dt <- fread(csv_files, skip = "lp__", drop = drop_cols, fill = TRUE)
    
    # Clean up duplicates
    real_params <- res_dt[, !grep("^theta_", names(res_dt)), with = FALSE]
    
    print(paste("Final LP:", tail(res_dt$lp__, 1)))
    print(tail(real_params, 1))
  }
}

# ==============================================================================
# 5. Generate latent trajectories from the OPTIMIZED parameter set (MLE)
# ==============================================================================

# Helper: convert cmdstanr matrix outputs to base types (important for JSON/init)
.as_numeric_vec <- function(x) as.numeric(x)
.as_numeric_mat <- function(x) { x <- as.matrix(x); storage.mode(x) <- "double"; x }

if (!is.null(fit_opt$return_codes()) && fit_opt$return_codes() == 0) {
  
  # --- A) Build init list from MLE (must match Stan 'parameters' block exactly) ---
  mle <- fit_opt$mle()
  
  init_mle <- list(
    mu_global      = .as_numeric_vec(mle$mu_global),
    sigma_line     = .as_numeric_vec(mle$sigma_line),
    beta_high      = .as_numeric_vec(mle$beta_high),
    z_line         = .as_numeric_mat(mle$z_line),
    sigma_beta     = .as_numeric_vec(mle$sigma_beta),
    z_beta         = .as_numeric_mat(mle$z_beta),
    
    mu_IC_raw      = .as_numeric_vec(mle$mu_IC_raw),
    sigma_IC       = .as_numeric_vec(mle$sigma_IC),
    beta_IC        = .as_numeric_vec(mle$beta_IC),
    z_IC           = .as_numeric_mat(mle$z_IC),
    sigma_beta_IC  = .as_numeric_vec(mle$sigma_beta_IC),
    z_beta_IC      = .as_numeric_mat(mle$z_beta_IC),
    
    phi_total      = as.numeric(mle$phi_total),
    phi_frac       = as.numeric(mle$phi_frac)
  )
  
  
  # --- C) One-shot forward run using fixed parameters (exactly the MLE) ---
  fit_sim_mle <- mod$sample(
    data = stan_data_sim,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 0,
    iter_sampling = 1,
    fixed_param = TRUE,
    init = list(init_mle),
    seed = SEED + 1,
    threads_per_chain = THREADS,
    refresh = 1
  )
  
  # --- D) Plot latent trajectories from the optimized params ---
  cat("\n>>> Plotting latent trajectories using OPTIMIZED (MLE) parameters...\n")
  plot_latent_like(fit_sim_mle, stan_data_sim, MODEL_NAME = paste0(MODEL_NAME, " (MLE)"), draw = 1)
  
} else {
  cat("\n>>> Skipping MLE latent simulation: optimize did not succeed.\n")
}

