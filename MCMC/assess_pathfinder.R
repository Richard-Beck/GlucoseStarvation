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
id_map <- setNames(names(stan_data$line_map), stan_data$line_map)
config <- jsonlite::read_json(json_file, simplifyVector = TRUE)
param_names <- config$param_names

# ==============================================================================
# 2. DATA PREPARATION FOR PPC (Posterior Predictive Intervals)
# ==============================================================================
cat("\nExtracting y_sim and noise parameters...\n")

# 1. Extract Latent Trajectories (Draws x Variables)
mu_mat <- fit$draws("y_sim", inc_warmup = FALSE, format = "matrix") 
var_names <- colnames(mu_mat)

# 2. Extract Noise Parameters (Draws)
phi_draws <- fit$draws(c("phi_total", "phi_frac"), inc_warmup = FALSE, format = "df")
phi_total_vec <- phi_draws$phi_total
phi_frac_vec  <- phi_draws$phi_frac

# 3. Helper: Map Dilution to Well
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

q_summ <- function(x) {
  c(median = quantile(x, 0.5, names=F), lo = quantile(x, 0.05, names=F), hi = quantile(x, 0.95, names=F))
}

res_list <- list()

# --- Separate Indices for Pairing ---
idx_NL <- idx_df %>% filter(type == "NL")
idx_ND <- idx_df %>% filter(type == "ND")
idx_G  <- idx_df %>% filter(type == "G")

idx_pairs <- inner_join(idx_NL, idx_ND, by = c("well_idx", "time"), suffix = c("_NL", "_ND"))

# --- A. Simulate Counts (NL + ND) ---
for (i in 1:nrow(idx_pairs)) {
  c_nl <- idx_pairs$col_idx_NL[i]
  c_nd <- idx_pairs$col_idx_ND[i]
  
  mu_nl_vec <- pmax(0, as.numeric(mu_mat[, c_nl]))
  mu_nd_vec <- pmax(0, as.numeric(mu_mat[, c_nd]))
  
  mu_tot <- mu_nl_vec + mu_nd_vec
  p_hat  <- mu_nl_vec / (mu_tot + 1e-12)
  p_hat  <- pmax(1e-6, pmin(1 - 1e-6, p_hat))
  
  rep_tot <- rnbinom(n = length(mu_tot), size = phi_total_vec, mu = mu_tot)
  
  alpha <- p_hat * phi_frac_vec
  beta  <- (1 - p_hat) * phi_frac_vec
  
  p_star <- rbeta(n = length(mu_tot), shape1 = alpha, shape2 = beta)
  rep_nl <- rbinom(n = length(mu_tot), size = rep_tot, prob = p_star)
  rep_nd <- rep_tot - rep_nl
  
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

# --- B. Simulate Glucose ---
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

# CHANGED: Metadata Map now uses ploidy_metric
meta_df <- data.frame(
  well_idx = 1:stan_data$N_wells,
  line_id  = stan_data$line_id,
  line_name = id_map[as.character(stan_data$line_id)],
  metric   = stan_data$ploidy_metric, # New continuous variable
  G0       = stan_data$G0_per_well,
  exp_id   = stan_data$exp_id
)

# Create synthetic labels for visualization: 0 is Baseline, >0 is High
meta_df$ploidy_lbl <- ifelse(meta_df$metric > 0.01, "High Ploidy", "Baseline Ploidy")
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
# 3. GENERATE PPC PLOTS
# ==============================================================================
if (!dir.exists("results")) dir.create("results")
cat("Generating PPC Plots (results/final_ppc_check.pdf)...\n")

pdf("results/final_ppc_check.pdf", width = 16, height = 12)

unique_lines <- sort(unique(meta_df$line_id))
cols  <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")
fills <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#0072B2")

for(lid in unique_lines) {
  lname <- id_map[as.character(lid)]
  d_sim <- sim_all %>% filter(line_id == lid)
  d_obs <- obs_all %>% filter(line_id == lid)
  if(nrow(d_sim) == 0) next
  
  # --- LEFT PANEL: Live/Dead Cells ---
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
    labs(title = paste0("Cell Counts | Line ", lname), x="Time (h)", y="Count")
  
  # --- RIGHT PANEL: Glucose ---
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
    labs(title = paste0("Glucose | Line ", lname), x="Time (h)", y="Conc (mM)")
  
  combined_plot <- (p1 | p2) + 
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = paste0("Posterior Check: Cell Line ", lid, " (", MODEL_NAME, ")"),
      subtitle = "Shaded regions represent 90% Posterior Predictive Intervals",
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

# CHANGED: Group by metric instead of binary flag
groups <- tibble(
  line_id = stan_data$line_id, 
  line_name = id_map[as.character(stan_data$line_id)],
  metric  = stan_data$ploidy_metric
) %>%
  distinct(line_id,line_name, metric) %>%
  arrange(line_id, metric) %>%
  mutate(
    # Create label: "Line 1 | Baseline" or "Line 1 | +2.0N"
    group_lbl = paste0(line_name, " | ", 
                       ifelse(metric < 0.01, "Baseline", paste0("+", round(metric,1), "N")))
  )

reconstructed_list <- list()
N_draws <- nrow(draws_df)
n_params <- length(param_names)

cat(sprintf("  Processing %d groups...\n", nrow(groups)))

for(g in 1:nrow(groups)) {
  l <- groups$line_id[g]
  p_met <- groups$metric[g] # Continuous value
  lbl <- groups$group_lbl[g]
  
  p_mat <- matrix(NA_real_, nrow = N_draws, ncol = n_params)
  colnames(p_mat) <- param_names
  
  for(pp in 1:n_params) {
    mu    <- draws_df[[paste0("mu_global[", pp, "]")]]
    sigma <- draws_df[[paste0("sigma_line[", pp, "]")]]
    z     <- draws_df[[paste0("z_line[", pp, ",", l, "]")]]
    beta  <- draws_df[[paste0("beta_high[", pp, "]")]]
    
    # CHANGED: Linear predictor uses p_met instead of h
    raw <- mu + sigma * z + beta * p_met
    
    # Transform back to constrained scale
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
    subtitle = "Effective parameters per cell line (incorporating ploidy offset)",
    x = "Value (Log Scale)", y = "Count"
  )
print(p_trans)
dev.off()

# ==============================================================================
# 5. PLOIDY EFFECT (beta_high) ANALYSIS
# ==============================================================================
cat("\nAnalyzing Continuous Ploidy Effects...\n")

draws_beta <- fit$draws("beta_high", inc_warmup = FALSE)
beta_df <- as_draws_df(draws_beta) %>%
  pivot_longer(cols = starts_with("beta_high"), names_to = "param_idx", values_to = "value")

beta_df$idx <- as.integer(str_extract(beta_df$param_idx, "[0-9]+"))
beta_df$param_name <- factor(param_names[beta_df$idx], levels = param_names)

pdf("results/ploidy_effect_posterior.pdf", width = 12, height = 8)

# CHANGED: Interpret beta as "Change per unit ploidy"
p_beta_raw <- ggplot(beta_df, aes(x = value, fill = stat(x > 0))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_histogram(bins = 60, alpha = 0.7, color = "white", size = 0.1) +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_manual(values = c("TRUE"="#E41A1C", "FALSE"="#377EB8"), guide="none") +
  theme_bw() +
  labs(title = paste0("Ploidy Effect (Beta Coefficient) - ", MODEL_NAME),
       subtitle = "Change in linear predictor per unit increase in Ploidy Metric\nRed (>0) = Increases with ploidy; Blue (<0) = Decreases with ploidy",
       x = "Beta Value", y = "Density")
print(p_beta_raw)

# Fold Change per Unit
fc_df <- beta_df %>% mutate(fold_change = exp(value))
p_fc <- ggplot(fc_df, aes(x = fold_change)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_density(fill = "gray40", alpha = 0.5) +
  facet_wrap(~ param_name, scales = "free") +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 10), labels = c("0.1x", "0.5x", "1x", "2x", "10x")) +
  theme_bw() +
  labs(title = "Implied Multiplier per Unit of Ploidy", 
       subtitle = "How much the parameter scales for every +1.0N difference in ploidy",
       x = "Multiplier (per unit N)", y = "Density")
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


library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# 8. PHYSIOLOGICAL PROFILING (Robust Math Fix)
# ==============================================================================

# ------------------------------------------------------------------------------
# A. Robust Math Helpers
# ------------------------------------------------------------------------------

# Safe log1p_exp: log(1 + exp(x))
# If x is large (> 50), log(1 + exp(x)) approx x.
# This prevents exp(700) -> Inf crashes.
log1p_exp_r <- function(x) {
  if (is.na(x)) return(NA)
  if (x > 50) return(x) 
  log1p(exp(x))
}

softcap_r <- function(x, cap) {
  if (is.na(x)) return(NA)
  if (x > cap) return(cap) 
  # cap - log(1 + exp(cap - x))
  # Use our robust helper here too for safety
  cap - log1p_exp_r(cap - x)
}

inv_logit_r <- function(x) {
  if (is.na(x)) return(NA)
  if (x > 100) return(1.0)
  if (x < -100) return(0.0)
  1 / (1 + exp(-x))
}

# ------------------------------------------------------------------------------
# B. ODE Derivative Function
# ------------------------------------------------------------------------------
ode_derivs <- function(t, y, p) {
  # 1. Unpack
  theta   <- p[1]; kp      <- p[2]; kd      <- p[3]
  nu      <- p[5]; rho     <- p[6]; r       <- p[7]
  sigma_G <- p[8]; KG      <- p[9]
  q50g_f  <- p[10]; q50d_f  <- p[11]
  na      <- p[12]; nd      <- p[13]
  
  fixed_G <- if ("fixed_G_val" %in% names(p)) p["fixed_G_val"] else NA
  
  # 2. State
  logNL   <- y[1]
  logND   <- y[2]
  G_state <- y[3]
  q_raw   <- y[4]
  
  G_use <- if (!is.na(fixed_G)) fixed_G else G_state
  
  # 3. Calculations
  NL <- exp(logNL)
  
  k_home_cap <- 4.0
  k_home <- exp(softcap_r(log(nu * kp + 1e-12), log(k_home_cap)))
  
  s <- nu * rho
  alpha_c <- s * r
  alpha_m <- s * (1.0 - r)
  
  eps_qstar <- 1e-6
  q_star <- (nu * (1.0 - rho)) / (nu + 1.0)
  q_star <- min(max(q_star, eps_qstar), 1.0 - eps_qstar)
  
  q50gN <- q50g_f * q_star
  q50dN <- q50d_f * q_star
  
  # --- FIX: ROBUST G SMOOTHING ---
  k_smooth_G <- 100.0
  # G_smooth = log(1 + exp(100*G)) / 100
  # Uses robust log1p_exp_r to avoid overflow when G > 7.1
  G_smooth <- log1p_exp_r(k_smooth_G * G_use) / k_smooth_G
  
  # --- FIX: ROBUST Q SMOOTHING ---
  k_smooth_q <- 50.0
  qN_low <- log1p_exp_r(k_smooth_q * q_raw) / k_smooth_q
  # For the second term, we do: 1.0 - (log1p_exp(...))
  term_high <- log1p_exp_r(k_smooth_q * (1.0 - qN_low)) / k_smooth_q
  qN <- 1.0 - term_high
  
  # Fluxes
  drive_term <- log1p_exp_r(k_smooth_q * (1.0 - qN)) / k_smooth_q
  gate  <- G_smooth / (KG + G_smooth + 1e-12)
  Jin   <- k_home * gate * drive_term
  
  # Regulation
  log_qN <- log(qN + 1e-12)
  reg_growth  <- inv_logit_r(na * (log_qN - log(q50gN + 1e-12)))
  term_d_hill <- inv_logit_r(nd * (log_qN - log(q50dN + 1e-12)))
  
  mu    <- kp * reg_growth
  delta <- kd * (1.0 - term_d_hill)
  
  b <- mu * (1.0 - NL/theta)
  
  # Derivatives
  du <- b - delta 
  dv <- delta * exp(logNL - logND)
  
  dG <- if (!is.na(fixed_G)) 0 else -NL * Jin / sigma_G
  
  dq_raw <- Jin - (alpha_m * kp) - (qN * b) - (alpha_c * b)
  
  list(c(du, dv, dG, dq_raw))
}

# ------------------------------------------------------------------------------
# C. Post-Processing Helper (Updated with Robust Math)
# ------------------------------------------------------------------------------
compute_derived_rates <- function(out_matrix, p) {
  theta   <- p[1]; kp      <- p[2]; kd      <- p[3]
  nu      <- p[5]; rho     <- p[6]; r       <- p[7]
  sigma_G <- p[8]; KG      <- p[9]
  q50g_f  <- p[10]; q50d_f  <- p[11]
  na      <- p[12]; nd      <- p[13]
  
  fixed_G <- if ("fixed_G_val" %in% names(p)) p["fixed_G_val"] else NA
  
  n <- nrow(out_matrix)
  mu_vec <- numeric(n); delta_vec <- numeric(n); qN_vec <- numeric(n); cons_vec <- numeric(n)
  
  k_home_cap <- 4.0
  k_home <- exp(softcap_r(log(nu * kp + 1e-12), log(k_home_cap)))
  
  q_star <- (nu * (1.0 - rho)) / (nu + 1.0)
  q_star <- min(max(q_star, 1e-6), 1.0 - 1e-6)
  
  q50gN <- q50g_f * q_star
  q50dN <- q50d_f * q_star
  
  k_smooth_G <- 100.0
  k_smooth_q <- 50.0
  
  for(i in 1:n) {
    logNL <- out_matrix[i, "logNL"]
    G_state <- out_matrix[i, "G"]
    q_raw <- out_matrix[i, "q_raw"]
    
    G_use <- if (!is.na(fixed_G)) fixed_G else G_state
    
    # Robust smoothing
    G_smooth <- log1p_exp_r(k_smooth_G * G_use) / k_smooth_G
    
    qN_low <- log1p_exp_r(k_smooth_q * q_raw) / k_smooth_q
    qN <- 1.0 - (log1p_exp_r(k_smooth_q * (1.0 - qN_low)) / k_smooth_q)
    
    drive <- log1p_exp_r(k_smooth_q * (1.0 - qN)) / k_smooth_q
    gate  <- G_smooth / (KG + G_smooth + 1e-12)
    Jin   <- k_home * gate * drive
    
    log_qN <- log(qN + 1e-12)
    reg_growth  <- inv_logit_r(na * (log_qN - log(q50gN + 1e-12)))
    term_d_hill <- inv_logit_r(nd * (log_qN - log(q50dN + 1e-12)))
    
    mu_vec[i]    <- kp * reg_growth
    delta_vec[i] <- kd * (1.0 - term_d_hill)
    qN_vec[i]    <- qN
    cons_vec[i]  <- Jin / sigma_G
  }
  
  return(data.frame(
    time = out_matrix[, "time"],
    mu = mu_vec,
    delta = delta_vec,
    q = qN_vec,
    consumption = cons_vec
  ))
}

# The rest of the plotting logic (E and F) remains exactly the same as the previous block.
# Re-run sections D, E, and F with these updated functions defined.

# ------------------------------------------------------------------------------
# D. Parameter Preparation
# ------------------------------------------------------------------------------
param_wide <- reconstructed_df %>%
  group_by(group, param) %>%
  summarise(val = median(value, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = val)

get_p_vec <- function(row_idx, df, p_names) {
  vec <- unlist(df[row_idx, p_names])
  return(as.numeric(vec))
}

# ------------------------------------------------------------------------------
# E. Plot 1: Steady State Glucose Consumption
# ------------------------------------------------------------------------------
cat("Calculating Steady State Curves (Separated Logic)...\n")
g_grid <- seq(0, 15, length.out = 50) 
ss_results <- list()

for(i in 1:nrow(param_wide)) {
  grp_name <- param_wide$group[i]
  p_vec <- get_p_vec(i, param_wide, param_names)
  if(any(is.na(p_vec))) next
  
  cons_vec <- numeric(length(g_grid))
  q_vec    <- numeric(length(g_grid))
  
  for(j in seq_along(g_grid)) {
    g_val <- g_grid[j]
    p_run <- c(p_vec, 1.0); names(p_run) <- c(param_names, "fixed_G_val")
    p_run["fixed_G_val"] <- g_val
    
    y0 <- c(logNL = log(1e-6), logND = log(1e-6), G = g_val, q_raw = 0.5)
    
    # 1. Solve ODE
    out <- suppressWarnings(
      ode(y = y0, times = c(0, 100), func = ode_derivs, parms = p_run, method = "lsoda")
    )
    
    # 2. Compute Rates on the output
    rates_df <- compute_derived_rates(out, p_run)
    
    # 3. Take final row
    final_rates <- rates_df[nrow(rates_df), ]
    
    cons_vec[j] <- final_rates$consumption
    q_vec[j]    <- final_rates$q
  }
  
  ss_results[[length(ss_results)+1]] <- data.frame(
    group = grp_name,
    Glucose = g_grid,
    Consumption = cons_vec,
    q_ss = q_vec
  )
}


# ------------------------------------------------------------------------------
# F. Plot 2: Step Drop Response
# ------------------------------------------------------------------------------
cat("Calculating Step Drop Dynamics (Separated Logic)...\n")
step_results <- list()

for(i in 1:nrow(param_wide)) {
  print(i)
  grp_name <- param_wide$group[i]
  p_vec <- get_p_vec(i, param_wide, param_names)
  if(any(is.na(p_vec))) next
  
  # --- Phase 1: Warmup at 10mM ---
  p_phase1 <- c(p_vec, 10.0); names(p_phase1) <- c(param_names, "fixed_G_val")
  y0 <- c(logNL = log(1e-6), logND = log(1e-6), G = 10.0, q_raw = 0.5)
  
  out_warmup <- suppressWarnings(
    ode(y = y0, times = seq(0,48,.1), func = ode_derivs, parms = p_phase1, method = "lsoda")
  )
  state_warmup <- out_warmup[nrow(out_warmup), -1] # Exclude time column
  
  # --- Phase 2: Drop to 0mM ---
  state_drop <- state_warmup
  state_drop["G"] <- 0.0
  
  p_phase2 <- c(p_vec, 0.0); names(p_phase2) <- c(param_names, "fixed_G_val")
  times_drop <- seq(0, 48, by = 0.1)
  
  # 1. Solve ODE
  out_drop <- suppressWarnings(
    ode(y = state_drop, times = times_drop, func = ode_derivs, parms = p_phase2, method = "lsoda")
  )
  
  # 2. Compute Rates
  rates_df <- compute_derived_rates(out_drop, p_phase2)
  rates_df$group <- grp_name
  
  step_results[[length(step_results)+1]] <- rates_df
}

# ==============================================================================
# G. Updated Visualization (2 Colors + Linetypes)
# ==============================================================================
library(stringr)

# 1. Process Steady State Data
# ------------------------------------------------------------------------------
if(length(ss_results) > 0) {
  df_ss <- bind_rows(ss_results) %>%
    mutate(
      # Extract "Line X"
      Line_ID = str_split_fixed(group, " \\| ", 2)[,1],
      # Determine Status
      Status = ifelse(str_detect(group, "Baseline"), "Baseline", "Above Baseline")
    )
  
  # Plot
  p_ss <- ggplot(df_ss, aes(x = Glucose, y = Consumption, color = Status, linetype = Line_ID)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c("Baseline" = "#377EB8", "Above Baseline" = "#E41A1C")) +
    theme_bw() +
    labs(
      title = "Per Capita Glucose Consumption (Steady State)",
      subtitle = "Consumption = Jin / sigma_G",
      x = "Extracellular Glucose (mM)",
      y = "Consumption Rate (mM/cell/h)",
      color = "Ploidy Status",
      linetype = "Cell Line"
    )
  print(p_ss)
}

# 2. Process Step Drop Data (Mu and Delta Only)
# ------------------------------------------------------------------------------
if(length(step_results) > 0) {
  df_step <- bind_rows(step_results) %>%
    mutate(
      Line_ID = str_split_fixed(group, " \\| ", 2)[,1],
      Status = ifelse(str_detect(group, "Baseline"), "Baseline", "Above Baseline")
    ) %>%
    pivot_longer(cols = c("mu", "delta", "q"), names_to = "metric", values_to = "value") %>%
    # FILTER: Remove 'q'
    filter(metric %in% c("mu", "delta"))
  
  # Define clean labels for the facets
  metric_labs <- c(
    "mu" = "Growth Rate (mu)", 
    "delta" = "Death Rate (delta)"
  )
  
  # Plot
  p_step <- ggplot(df_step, aes(x = time, y = value, color = Status, linetype = Line_ID)) +
    geom_line(linewidth = 1.0) +
    facet_wrap(metric~Line_ID, scales = "free_y", ncol = 5, labeller = as_labeller(metric_labs)) +
    scale_color_manual(values = c("Baseline" = "#377EB8", "Above Baseline" = "#E41A1C")) +
    theme_bw() +
    scale_x_log10()+
    labs(
      title = "Response to Glucose Drop (10mM -> 0mM)",
      subtitle = "Immediate physiological adaptation",
      x = "Time since Drop (h)",
      y = "Rate (1/h)",
      color = "Ploidy Status",
      linetype = "Cell Line"
    )
  print(p_step)
}

# Save Steady State Consumption Plot
ggsave("results/physio_steady_state_uptake.pdf", plot = p_ss, width = 8, height = 6)

# Save Step Drop Response Plot
ggsave("results/physio_starvation_step_drop.pdf", plot = p_step, width = 10, height = 8)
cat(">>> Physiological plots generated.\n")