library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)

# ==============================================================================
# 1. SETUP & LOAD
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME <- if (length(args) >= 1) args[1] else "model_q"

# Files
res_file  <- file.path("results", paste0("pathfinder_", MODEL_NAME, ".Rds"))
data_file <- "data/stan_ready_data.Rds" 
json_file <- file.path("MCMC", paste0(MODEL_NAME, ".json"))

if (!file.exists(res_file)) stop("Result file not found: ", res_file)
if (!file.exists(data_file)) stop("Data file not found: ", data_file)
if (!file.exists(json_file)) stop("JSON config not found: ", json_file)

cat(sprintf(">>> Loading Pathfinder Results for: %s\n", MODEL_NAME))
fit <- readRDS(res_file)
stan_data <- readRDS(data_file)
config <- jsonlite::read_json(json_file, simplifyVector = TRUE)
param_names <- config$param_names

# ==============================================================================
# 2. DATA PREPARATION FOR PPC (Posterior Predictive Intervals)
# ==============================================================================
cat("\nExtracting y_sim and noise parameters...\n")

# 1. Extract Latent Trajectories (Draws x Variables)
mu_mat <- fit$draws("y_sim", inc_warmup = FALSE, format = "matrix") 
var_names <- colnames(mu_mat)

# 2. Extract Noise Parameters (Draws) - UPDATED NAMES
phi_draws <- fit$draws(c("phi_total", "phi_frac"), inc_warmup = FALSE, format = "df")
phi_total_vec <- phi_draws$phi_total
phi_frac_vec  <- phi_draws$phi_frac

# 3. Helper: Map Dilution to Well (needed for Glucose calc)
dil_by_well <- rep(NA_real_, stan_data$N_wells)
for (w in 1:stan_data$N_wells) {
  ii <- which(stan_data$well_idx_gluc == w)
  if (length(ii) > 0) {
    u <- unique(stan_data$dilution[ii])
    dil_by_well[w] <- u[1] 
  } else {
    dil_by_well[w] <- 1.0 
  }
}

# 4. Parse Indices
idx_df <- tibble(variable = var_names) %>%
  mutate(clean = gsub("y_sim\\[|\\]", "", variable)) %>%
  separate(clean, c("w","t","s"), sep = ",", convert = TRUE) %>%
  mutate(
    well_idx = w,
    time     = stan_data$t_grid[t],
    type     = case_when(s==1 ~ "NL", s==2 ~ "ND", s==3 ~ "G", TRUE ~ "Other"),
    col_idx  = 1:n()
  ) %>%
  filter(type %in% c("NL","ND","G"))

cat("Simulating posterior predictive replicates (Total/Fraction Model)...\n")

# Helper for quantiles
q_summ <- function(x) {
  c(median = quantile(x, 0.5, names=F), lo = quantile(x, 0.05, names=F), hi = quantile(x, 0.95, names=F))
}

res_list <- list()

# --- Separate Indices for Pairing ---
idx_NL <- idx_df %>% filter(type == "NL")
idx_ND <- idx_df %>% filter(type == "ND")
idx_G  <- idx_df %>% filter(type == "G")

# Join NL and ND to ensure we simulate them as a pair
# Assumes the structure of y_sim is consistent (it is, coming from Stan)
idx_pairs <- inner_join(idx_NL, idx_ND, by = c("well_idx", "time"), suffix = c("_NL", "_ND"))

# --- A. Simulate Counts (NL + ND) ---
for (i in 1:nrow(idx_pairs)) {
  c_nl <- idx_pairs$col_idx_NL[i]
  c_nd <- idx_pairs$col_idx_ND[i]
  
  # Latent means
  mu_nl_vec <- pmax(0, as.numeric(mu_mat[, c_nl]))
  mu_nd_vec <- pmax(0, as.numeric(mu_mat[, c_nd]))
  
  # 1. Calculate Latent Total and Fraction
  mu_tot <- mu_nl_vec + mu_nd_vec
  p_hat  <- mu_nl_vec / (mu_tot + 1e-12)
  p_hat  <- pmax(1e-6, pmin(1 - 1e-6, p_hat))
  
  # 2. Simulate Total Count ~ NegBinomial(mu_tot, phi_total)
  rep_tot <- rnbinom(n = length(mu_tot), size = phi_total_vec, mu = mu_tot)
  
  # 3. Simulate Alive Count ~ BetaBinomial(rep_tot, p_hat, phi_frac)
  alpha <- p_hat * phi_frac_vec
  beta  <- (1 - p_hat) * phi_frac_vec
  
  p_star <- rbeta(n = length(mu_tot), shape1 = alpha, shape2 = beta)
  rep_nl <- rbinom(n = length(mu_tot), size = rep_tot, prob = p_star)
  rep_nd <- rep_tot - rep_nl
  
  # 4. Summarize
  qs_nl <- q_summ(rep_nl)
  qs_nd <- q_summ(rep_nd)
  
  res_list[[length(res_list) + 1]] <- data.frame(
    well_idx = idx_pairs$well_idx[i], time = idx_pairs$time[i], type = "NL",
    median = qs_nl[1], lo = qs_nl[2], hi = qs_nl[3]
  )
  res_list[[length(res_list) + 1]] <- data.frame(
    well_idx = idx_pairs$well_idx[i], time = idx_pairs$time[i], type = "ND",
    median = qs_nd[1], lo = qs_nd[2], hi = qs_nd[3]
  )
}

# --- B. Simulate Glucose (Standard Lognormal) ---
for (i in 1:nrow(idx_G)) {
  c_g <- idx_G$col_idx[i]
  mu_g <- pmax(0, as.numeric(mu_mat[, c_g]))
  
  w <- idx_G$well_idx[i]
  e <- stan_data$exp_id[w]
  d <- dil_by_well[w]
  
  a <- stan_data$calib_a_fixed[e]
  b <- stan_data$calib_b_fixed[e]
  sig <- stan_data$calib_sigma_fixed[e] 
  
  mu_lum <- a * mu_g * d + b
  yrep_lum <- exp(rnorm(n = length(mu_g), mean = log(mu_lum + 1e-12), sd = sig))
  yrep <- pmax(0, (yrep_lum - b) / (a * d))
  
  qs <- q_summ(yrep)
  res_list[[length(res_list) + 1]] <- data.frame(
    well_idx = w, time = idx_G$time[i], type = "G",
    median = qs[1], lo = qs[2], hi = qs[3]
  )
}

cat("Formatting Simulation Data...\n")
summ_clean <- bind_rows(res_list)

# Metadata Map
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

# ==============================================================================
# 3. GENERATE PPC PLOTS (Free Scales)
# ==============================================================================
if (!dir.exists("results")) dir.create("results")
cat("Generating PPC Plots (results/final_ppc_check.pdf)...\n")

pdf("results/final_ppc_check.pdf", width = 16, height = 12)

unique_lines <- sort(unique(meta_df$line_id))
cols  <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")
fills <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")

for(lid in unique_lines) {
  d_sim <- sim_all %>% filter(line_id == lid)
  d_obs <- obs_all %>% filter(line_id == lid)
  if(nrow(d_sim) == 0) next
  
  # --- LEFT PANEL: Live/Dead Cells ---
  p1 <- ggplot() +
    # PREDICTION INTERVAL (Shaded Region)
    geom_ribbon(data = d_sim %>% filter(type != "G"), 
                aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
    # MEDIAN TRAJECTORY
    geom_line(data = d_sim %>% filter(type != "G"), 
              aes(x=time, y=median, color=type, group=interaction(well_idx, type)), linewidth=0.8) +
    # OBSERVED DATA
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
  
  # --- RIGHT PANEL: Glucose ---
  p2 <- ggplot() +
    # PREDICTION INTERVAL (Shaded Region)
    geom_ribbon(data = d_sim %>% filter(type == "G"), 
                aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
    # MEDIAN TRAJECTORY
    geom_line(data = d_sim %>% filter(type == "G"), 
              aes(x=time, y=median, color=type, group=interaction(well_idx, type)), linewidth=0.8) +
    # OBSERVED DATA
    geom_point(data = d_obs %>% filter(type == "G"), 
               aes(x=time, y=value, color=type), size=1.5, alpha=0.7) +
    facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
    scale_y_continuous() +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = fills) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill="#f0f0f0")) +
    labs(title = paste0("Glucose | Line ", lid), x="Time (h)", y="Conc (mM)")
  
  combined_plot <- (p1 | p2) + 
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = paste0("Posterior Check: Cell Line ", lid, " (", MODEL_NAME, ")"),
      subtitle = "Shaded regions represent 90% Posterior Predictive Intervals (simulated data)",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  print(combined_plot)
  cat("  Plot generated for Line", lid, "\n")
}
dev.off()

# ==============================================================================
# 4. RECONSTRUCT TRANSFORMED PARAMETERS
# ==============================================================================
cat("\nReconstructing effective ODE parameters (Transformed)...\n")
vars_needed <- c("mu_global", "sigma_line", "z_line", "beta_high")
draws_df <- as_draws_df(fit$draws(variables = vars_needed, inc_warmup = FALSE))

softcap <- function(x, cap) cap - log1p(exp(cap - x))
inv_logit <- function(x) 1 / (1 + exp(-x))

groups <- tibble(
  line_id = stan_data$line_id, 
  high = stan_data$is_high_ploidy
) %>%
  distinct(line_id, high) %>%
  arrange(line_id, high) %>%
  mutate(group_lbl = paste0("Line ", line_id, " | ", ifelse(high==1, "High", "Low")))

reconstructed_list <- list()
N_draws <- nrow(draws_df)
n_params <- length(param_names)

cat(sprintf("  Processing %d groups...\n", nrow(groups)))

for(g in 1:nrow(groups)) {
  l <- groups$line_id[g]
  h <- groups$high[g]
  lbl <- groups$group_lbl[g]
  
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
      if (pp %in% c(6, 8)) {
        p_mat[, pp] <- 1.0 + exp(softcap(raw, cap_hill))
      } else {
        p_mat[, pp] <- exp(softcap(raw, cap_main))
      }
    } else if (MODEL_NAME == "model_q") {
      cap_main <- 15.0; cap_hill <- 2.7
      if (pp %in% c(6, 7, 10, 11)) {
        p_mat[, pp] <- inv_logit(raw)
      } else if (pp %in% c(12, 13)) {
        p_mat[, pp] <- exp(softcap(raw, cap_hill))
      } else if (pp == 1) {
        p_mat[, pp] <- exp(raw)
      } else {
        p_mat[, pp] <- exp(softcap(raw, cap_main))
      }
    }
  }
  
  df_g <- as.data.frame(p_mat)
  df_g$group <- lbl
  df_g$chain <- if(".chain" %in% names(draws_df)) as.factor(draws_df$.chain) else as.factor(1)
  reconstructed_list[[g]] <- pivot_longer(df_g, cols = all_of(param_names), names_to = "param", values_to = "value")
}

reconstructed_df <- bind_rows(reconstructed_list)
reconstructed_df$group <- factor(reconstructed_df$group, levels = groups$group_lbl)
reconstructed_df$param <- factor(reconstructed_df$param, levels = param_names)

cat("Generating Transformed Parameter Plot (results/final_posterior_parameters.pdf)...\n")
plot_height <- max(8, length(unique(reconstructed_df$group)) * 2.5)

pdf("results/final_posterior_parameters.pdf", width = 20, height = plot_height)
p_trans <- ggplot(reconstructed_df, aes(x = value)) +
  geom_histogram(bins = 50, fill = "#377EB8", alpha = 0.6, color = NA) +
  facet_grid(group ~ param, scales = "free") +
  scale_x_log10() +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 9),
    strip.text.x = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  ) +
  labs(
    title = paste0("Posterior Parameter Distributions (Transformed) - ", MODEL_NAME),
    subtitle = "Effective parameters per cell line (including random effects & transforms)",
    x = "Value (Log Scale)", y = "Count"
  )
print(p_trans)
dev.off()

# ==============================================================================
# 5. WGD EFFECT (beta_high) ANALYSIS
# ==============================================================================
cat("\nAnalyzing Whole Genome Duplication (WGD) Effects...\n")

draws_beta <- fit$draws("beta_high", inc_warmup = FALSE)
beta_df <- as_draws_df(draws_beta) %>%
  pivot_longer(cols = starts_with("beta_high"), names_to = "param_idx", values_to = "value")

beta_df$idx <- as.integer(str_extract(beta_df$param_idx, "[0-9]+"))
beta_df$param_name <- factor(param_names[beta_df$idx], levels = param_names)

pdf("results/wgd_effect_posterior.pdf", width = 12, height = 8)

p_beta_raw <- ggplot(beta_df, aes(x = value, fill = stat(x > 0))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_histogram(bins = 60, alpha = 0.7, color = "white", size = 0.1) +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_manual(values = c("TRUE"="#E41A1C", "FALSE"="#377EB8"), guide="none") +
  theme_bw() +
  labs(title = paste0("WGD Effect (Raw Beta Parameters) - ", MODEL_NAME),
       subtitle = "Value > 0 (Red) = WGD increases raw value. Value < 0 (Blue) = WGD decreases it.",
       x = "Raw Beta Value", y = "Density")
print(p_beta_raw)

fc_df <- beta_df %>% mutate(fold_change = exp(value))
p_fc <- ggplot(fc_df, aes(x = fold_change)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_density(fill = "gray40", alpha = 0.5) +
  facet_wrap(~ param_name, scales = "free") +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 10), labels = c("0.1x", "0.5x", "1x", "2x", "10x")) +
  theme_bw() +
  labs(title = "Implied WGD Fold Change", x = "Fold Change (High vs Low)", y = "Density")
print(p_fc)
dev.off()

# ==============================================================================
# 6. ALL RAW PARAMETERS (mu, sigma, beta)
# ==============================================================================
cat("\nGenerating Comprehensive Raw Parameter Plots (results/all_raw_parameters.pdf)...\n")

vars_raw <- c("mu_global", "beta_high", "sigma_line")
draws_all_raw <- fit$draws(vars_raw, inc_warmup = FALSE)
all_raw_df <- as_draws_df(draws_all_raw) %>%
  pivot_longer(cols = everything(), names_to = "full_name", values_to = "value") %>%
  filter(!full_name %in% c(".chain", ".iteration", ".draw"))

all_raw_df <- all_raw_df %>%
  mutate(
    type = str_extract(full_name, "^[a-z_]+"),
    idx_s = str_extract(full_name, "[0-9]+"),
    idx = as.integer(idx_s),
    param_name = factor(param_names[idx], levels = param_names)
  )

pdf("results/all_raw_parameters.pdf", width = 14, height = 10)
p_all <- ggplot(all_raw_df, aes(x = value, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  labs(title = paste0("All Hierarchical Components (Raw) - ", MODEL_NAME), x = "Raw Value", y = "Count")
print(p_all)
dev.off()

cat("\n>>> Assessment Complete.\n")

# ==============================================================================
# 7. DIAGNOSTICS
# ==============================================================================
cat("\nGenerating Diagnostics...\n")

draws_trace <- fit$draws("mu_global", inc_warmup = FALSE)
trace_df <- as_draws_df(draws_trace) %>%
  mutate(.draw = 1:n()) %>% 
  pivot_longer(cols = starts_with("mu_global"), names_to = "param_idx", values_to = "value")

trace_df$idx <- as.integer(str_extract(trace_df$param_idx, "[0-9]+"))
trace_df$param_name <- factor(param_names[trace_df$idx], levels = param_names)

pdf("results/pseudo_trace_check.pdf", width = 12, height = 8)
p_trace <- ggplot(trace_df, aes(x = .draw, y = value)) +
  geom_point(alpha = 0.4, size = 0.8, color = "#2c3e50") +
  facet_wrap(~ param_name, scales = "free_y") +
  theme_bw() +
  labs(title = "Pseudo-Trace of Pathfinder Draws", x = "Draw Index", y = "Raw Parameter Value")
print(p_trace)
dev.off()

codes <- fit$return_codes()
n_fail <- sum(codes != 0)
cat(sprintf("  Total Paths Run: %d\n", length(codes)))
cat(sprintf("  Failed (Crashed): %d\n", n_fail))

draws_lp <- fit$draws("lp__", inc_warmup = FALSE)
lp_df <- as_draws_df(draws_lp)
pdf("results/diagnostic_lp_check.pdf", width = 8, height = 6)
p_lp <- ggplot(lp_df, aes(x = lp__)) +
  geom_histogram(bins = 50, fill = "purple", alpha = 0.7, color = "black") +
  theme_bw() +
  labs(title = "Log-Probability (lp__) Distribution", x = "Log Probability", y = "Count of Draws")
print(p_lp)
dev.off()