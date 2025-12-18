library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
library(jsonlite)

# ==============================================================================
# 1. SETUP & LOAD
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME <- if (length(args) >= 1) args[1] else "model_B"

# Files
res_file  <- file.path("results", paste0("fit_", MODEL_NAME, ".Rds"))
data_file <- "data/stan_ready_data.Rds" 
json_file <- file.path("MCMC", paste0(MODEL_NAME, ".json"))
stan_file <- file.path("MCMC", paste0(MODEL_NAME, ".stan"))

if (!file.exists(res_file)) stop("Result file not found: ", res_file)
if (!file.exists(data_file)) stop("Data file not found: ", data_file)

cat(sprintf("\n>>> Loading MCMC Results for: %s\n", MODEL_NAME))
fit <- readRDS(res_file)
stan_data <- readRDS(data_file)
config <- read_json(json_file, simplifyVector = TRUE)
param_names <- config$param_names

# ==============================================================================
# 2. MCMC DIAGNOSTICS
# ==============================================================================
cat("\n[1/6] Running MCMC Diagnostics...\n")

# 2a. NUTS Statistics
sampler_diagnostics <- fit$sampler_diagnostics()
divergences <- sum(sampler_diagnostics[,,"divergent__"])
max_treedepths <- sum(sampler_diagnostics[,,"treedepth__"] >= 12)

cat(sprintf("  Total Divergences: %d  %s\n", divergences, ifelse(divergences > 0, "<-- WARNING!", "(Good)")))
cat(sprintf("  Max Treedepth Hits: %d\n", max_treedepths))

# 2b. Convergence (Rhat & ESS)
draws_mu <- fit$draws("mu_global")
summ_stats <- summarise_draws(draws_mu, rhat, ess_bulk, ess_tail)

pdf("results/mcmc_diagnostics.pdf", width = 12, height = 8)

p_rhat <- ggplot(summ_stats, aes(x = rhat)) +
  geom_histogram(bins = 30, fill = "#377EB8", color = "white") +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "R-hat Distribution (Should be < 1.01)", x = "R-hat", y = "Count")

p_ess <- ggplot(summ_stats, aes(x = ess_bulk)) +
  geom_histogram(bins = 30, fill = "#4DAF4A", color = "white") +
  geom_vline(xintercept = 400, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "ESS Bulk (Should be > 400)", x = "Effective Sample Size", y = "Count")

# Trace Plots
p_trace <- mcmc_trace(draws_mu, pars = variables(draws_mu)[1:min(6, config$N_params)]) +
  ggtitle("Trace Plots: First 6 Global Parameters")

layout <- (p_rhat | p_ess) / p_trace
print(layout)

p_rank <- mcmc_rank_overlay(draws_mu, pars = variables(draws_mu)[1:min(6, config$N_params)]) +
  ggtitle("Rank Overlay Plots (Should be uniform)")
print(p_rank)

dev.off()
cat("  Diagnostics saved to 'results/mcmc_diagnostics.pdf'\n")

if (divergences > 0) cat("  WARNING: Divergences detected. Check 'mcmc_diagnostics.pdf'.\n")
if (any(summ_stats$rhat > 1.05, na.rm=T)) cat("  WARNING: Rhat > 1.05 detected.\n")


# ==============================================================================
# 3. GENERATE QUANTITIES (Smart Detection & Data Injection)
# ==============================================================================
cat("\n[2/6] Preparing Posterior Predictions...\n")

# --- CRITICAL FIX: Replicate run_fit.R logic to complete stan_data ---
# The .stan file expects these variables, but they aren't in the .Rds file.
if (MODEL_NAME == "model_B") {
  max_N <- max(stan_data$N_obs, na.rm = TRUE)
  # Re-calculate the dynamic Theta prior used during fitting
  config$prior_means[1] <- log(max_N * 1.5)
  config$prior_sds[1]   <- 0.5
}
stan_data$prior_ode_mean <- config$prior_means
stan_data$prior_ode_sd   <- config$prior_sds
stan_data$mode           <- 0  # 0 = Sampling Mode
# ---------------------------------------------------------------------

# Now check for y_sim
all_vars <- variables(fit$draws())
ysim_vars <- all_vars[grepl("^y_sim", all_vars)]
need_regen <- TRUE

if (length(ysim_vars) > 0) {
  cat("  'y_sim' variables found in output. Checking validity of first element...\n")
  first_var <- ysim_vars[1]
  check_subset <- fit$draws(variables = first_var)
  
  if (all(is.na(check_subset))) {
    cat("  'y_sim' is filled with NaNs. Regeneration REQUIRED.\n")
    need_regen <- TRUE
  } else {
    cat("  'y_sim' appears valid. Extracting full posterior...\n")
    draws_sim <- fit$draws(variables = "y_sim")
    need_regen <- FALSE
  }
} else {
  cat("  'y_sim' NOT found in fit object. Regeneration REQUIRED.\n")
  need_regen <- TRUE
}

if (need_regen) {
  cat("  Running standalone Generate Quantities...\n")
  
  mod_gq <- cmdstan_model(stan_file, quiet = TRUE, force_recompile = FALSE)
  
  # Ensure calc_sim is ON for this pass
  stan_data_gq <- stan_data
  stan_data_gq$calc_sim <- 1
  
  gq_fit <- mod_gq$generate_quantities(
    fit,
    data = stan_data_gq,
    parallel_chains = 4,
    seed = 123
  )
  
  draws_sim <- gq_fit$draws("y_sim")
  cat("  Generation complete.\n")
}

# ==============================================================================
# 4. POSTERIOR PREDICTIVE CHECK (PPC)
# ==============================================================================
cat("\n[3/6] Generating PPC Plots...\n")

# Robust Summary (na.rm = TRUE handles rare solver failures)
summ_sim <- summarise_draws(
  draws_sim, 
  median = function(x) quantile(x, 0.5, names = FALSE, na.rm = TRUE),
  lo     = function(x) quantile(x, 0.05, names = FALSE, na.rm = TRUE),
  hi     = function(x) quantile(x, 0.95, names = FALSE, na.rm = TRUE)
)

summ_clean <- summ_sim %>%
  mutate(clean = gsub("y_sim\\[|\\]", "", variable)) %>%
  separate(clean, c("w", "t", "s"), sep = ",", convert = TRUE) %>%
  mutate(
    well_idx = w,
    time     = stan_data$t_grid[t],
    type     = case_when(s==1 ~ "NL", s==2 ~ "ND", s==3 ~ "G", TRUE ~ "Other")
  ) %>%
  filter(type %in% c("NL", "ND", "G")) %>%
  select(well_idx, time, type, median, lo, hi)

# Rebuild Metadata
meta_df <- data.frame(
  well_idx = 1:stan_data$N_wells,
  line_id  = stan_data$line_id,
  is_high  = stan_data$is_high_ploidy,
  G0       = stan_data$G0_per_well,
  exp_id   = stan_data$exp_id
)
meta_df$ploidy_lbl <- ifelse(meta_df$is_high == 1, "High Ploidy", "Low Ploidy")
meta_df$G0_lbl     <- factor(paste0(meta_df$G0, " mM"), levels = paste0(sort(unique(meta_df$G0)), " mM"))

# Observed Data
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
  ))

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
    value = pmax(0, (lum - b)/(a*dilution)),
    type  = "G"
  ) %>%
  select(well_idx, time, value, type)

obs_all <- bind_rows(obs_counts, obs_gluc) %>% left_join(meta_df, by = "well_idx")
sim_all <- summ_clean %>% left_join(meta_df, by = "well_idx")

# PLOTTING
pdf("results/final_ppc_check_mcmc.pdf", width = 16, height = 12)

unique_lines <- sort(unique(meta_df$line_id))
cols  <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")
fills <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")

for(lid in unique_lines) {
  d_sim <- sim_all %>% filter(line_id == lid)
  d_obs <- obs_all %>% filter(line_id == lid)
  if(nrow(d_sim) == 0) next
  
  # Panel 1: Cells
  p1 <- ggplot() +
    geom_ribbon(data = d_sim %>% filter(type != "G"), 
                aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
    geom_line(data = d_sim %>% filter(type != "G"), 
              aes(x=time, y=median, color=type, group=interaction(well_idx, type)), linewidth=0.8) +
    geom_point(data = d_obs %>% filter(type != "G"), 
               aes(x=time, y=value, color=type, shape=type), size=1.5, alpha=0.7) +
    facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
    scale_y_continuous(labels = scales::comma) + 
    scale_color_manual(values = cols) +
    scale_fill_manual(values = fills) +
    scale_shape_manual(values = c("NL"=16, "ND"=17)) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill="#f0f0f0")) +
    labs(title = paste0("Cell Counts | Line ", lid), x="Time (h)", y="Count")
  
  # Panel 2: Glucose
  p2 <- ggplot() +
    geom_ribbon(data = d_sim %>% filter(type == "G"), 
                aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
    geom_line(data = d_sim %>% filter(type == "G"), 
              aes(x=time, y=median, color=type, group=interaction(well_idx, type)), linewidth=0.8) +
    geom_point(data = d_obs %>% filter(type == "G"), 
               aes(x=time, y=value, color=type), size=1.5, alpha=0.7) +
    facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
    scale_y_continuous() +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = fills) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill="#f0f0f0")) +
    labs(title = paste0("Glucose | Line ", lid), x="Time (h)", y="Conc (mM)")
  
  combined <- (p1 | p2) + plot_annotation(title = paste0("MCMC Fit: Line ", lid), theme=theme(plot.title=element_text(size=16, face="bold")))
  print(combined)
  cat("  Plot generated for Line", lid, "\n")
}
dev.off()


# ==============================================================================
# 5. RECONSTRUCT TRANSFORMED PARAMETERS
# ==============================================================================
cat("\n[4/6] Reconstructing Parameters...\n")

vars_needed <- c("mu_global", "sigma_line", "z_line", "beta_high")
draws_df <- as_draws_df(fit$draws(variables = vars_needed)) 

softcap <- function(x, cap) cap - log1p(exp(cap - x))
inv_logit <- function(x) 1 / (1 + exp(-x))

groups <- tibble(
  line_id = stan_data$line_id, 
  high = stan_data$is_high_ploidy
) %>% distinct(line_id, high) %>% arrange(line_id, high) %>%
  mutate(group_lbl = paste0("Line ", line_id, " | ", ifelse(high==1, "High", "Low")))

reconstructed_list <- list()
N_draws <- nrow(draws_df)
n_params <- length(param_names)

for(g in 1:nrow(groups)) {
  l <- groups$line_id[g]; h <- groups$high[g]; lbl <- groups$group_lbl[g]
  p_mat <- matrix(NA_real_, nrow = N_draws, ncol = n_params)
  colnames(p_mat) <- param_names
  
  for(pp in 1:n_params) {
    mu    <- draws_df[[paste0("mu_global[", pp, "]")]]
    sigma <- draws_df[[paste0("sigma_line[", pp, "]")]]
    z     <- draws_df[[paste0("z_line[", pp, ",", l, "]")]]
    beta  <- draws_df[[paste0("beta_high[", pp, "]")]]
    raw <- mu + sigma * z + beta * h
    
    if (MODEL_NAME == "model_B") {
      cap_main <- 40.0; cap_hill <- 6.0
      if (pp %in% c(6, 8)) { p_mat[, pp] <- 1.0 + exp(softcap(raw, cap_hill)) } 
      else { p_mat[, pp] <- exp(softcap(raw, cap_main)) }
    } else if (MODEL_NAME == "model_q") {
      cap_main <- 15.0; cap_hill <- 2.7
      if (pp %in% c(6, 7, 10, 11)) { p_mat[, pp] <- inv_logit(raw) }
      else if (pp %in% c(12, 13)) { p_mat[, pp] <- exp(softcap(raw, cap_hill)) }
      else if (pp == 1) { p_mat[, pp] <- exp(raw) }
      else { p_mat[, pp] <- exp(softcap(raw, cap_main)) }
    }
  }
  df_g <- as.data.frame(p_mat)
  df_g$group <- lbl
  if(".chain" %in% names(draws_df)) df_g$chain <- as.factor(draws_df$.chain)
  reconstructed_list[[g]] <- pivot_longer(df_g, cols = all_of(param_names), names_to = "param", values_to = "value")
}

reconstructed_df <- bind_rows(reconstructed_list)
reconstructed_df$group <- factor(reconstructed_df$group, levels = groups$group_lbl)
reconstructed_df$param <- factor(reconstructed_df$param, levels = param_names)

pdf("results/final_posterior_parameters_mcmc.pdf", width = 20, height = max(8, length(unique(reconstructed_df$group)) * 2.5))
p_trans <- ggplot(reconstructed_df, aes(x = value)) +
  geom_histogram(bins = 50, fill = "#377EB8", alpha = 0.6, color = NA) +
  facet_grid(group ~ param, scales = "free") +
  scale_x_log10() +
  theme_bw() +
  labs(title = paste0("Posterior Parameter Distributions (Transformed) - ", MODEL_NAME), x = "Value (Log Scale)", y = "Count")
print(p_trans)
dev.off()


# ==============================================================================
# 6. WGD EFFECT ANALYSIS
# ==============================================================================
cat("\n[5/6] Analyzing WGD Effects...\n")

draws_beta <- fit$draws("beta_high")
beta_df <- as_draws_df(draws_beta) %>% pivot_longer(cols = starts_with("beta_high"), names_to = "param_idx", values_to = "value")
beta_df$idx <- as.integer(str_extract(beta_df$param_idx, "[0-9]+"))
beta_df$param_name <- factor(param_names[beta_df$idx], levels = param_names)

pdf("results/wgd_effect_posterior_mcmc.pdf", width = 12, height = 8)
p_beta <- ggplot(beta_df, aes(x = value, fill = stat(x > 0))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_manual(values = c("TRUE"="#E41A1C", "FALSE"="#377EB8"), guide="none") +
  theme_bw() +
  labs(title = "WGD Effect (Raw Beta Parameters)", x = "Raw Beta Value", y = "Density")
print(p_beta)

p_fc <- ggplot(beta_df %>% mutate(fc=exp(value)), aes(x = fc)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_density(fill = "gray40", alpha = 0.5) +
  facet_wrap(~ param_name, scales = "free") +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 10), labels = c("0.1x", "0.5x", "1x", "2x", "10x")) +
  theme_bw() +
  labs(title = "Implied WGD Fold Change", x = "Fold Change", y = "Density")
print(p_fc)
dev.off()


# ==============================================================================
# 7. RAW PARAMETERS
# ==============================================================================
cat("\n[6/6] Plotting Raw Parameters...\n")
vars_raw <- c("mu_global", "beta_high", "sigma_line")
draws_all_raw <- fit$draws(vars_raw)
all_raw_df <- as_draws_df(draws_all_raw) %>%
  pivot_longer(cols = everything(), names_to = "full_name", values_to = "value") %>%
  filter(!full_name %in% c(".chain", ".iteration", ".draw")) %>%
  mutate(
    type = str_extract(full_name, "^[a-z_]+"),
    idx = as.integer(str_extract(full_name, "[0-9]+")),
    param_name = factor(param_names[idx], levels = param_names)
  )

pdf("results/all_raw_parameters_mcmc.pdf", width = 14, height = 10)
p_all <- ggplot(all_raw_df, aes(x = value, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  labs(title = "All Hierarchical Components (Raw Scale)", x = "Raw Value", y = "Count")
print(p_all)
dev.off()

cat("\n>>> MCMC Assessment Complete.\n")