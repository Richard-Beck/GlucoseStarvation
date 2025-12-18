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
data_file <- "data/stan_ready_data.Rds" # Assuming data stays in data/
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
# 2. DATA PREPARATION FOR PPC
# ==============================================================================
cat("\nExtracting y_sim for PPC...\n")
draws_sim <- fit$draws("y_sim", inc_warmup = FALSE)

cat("Calculating Credible Intervals (Summary)...\n")
summ_sim <- summarize_draws(
  draws_sim, 
  median = function(x) quantile(x, 0.5, names = FALSE),
  lo     = function(x) quantile(x, 0.05, names = FALSE),
  hi     = function(x) quantile(x, 0.95, names = FALSE)
)

cat("Formatting Simulation Data...\n")
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
    geom_ribbon(data = d_sim %>% filter(type != "G"), 
                aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
    geom_line(data = d_sim %>% filter(type != "G"), 
              aes(x=time, y=median, color=type, group=interaction(well_idx, type)), linewidth=0.8) +
    geom_point(data = d_obs %>% filter(type != "G"), 
               aes(x=time, y=value, color=type, shape=type), size=1.5, alpha=0.7) +
    # FACET: Rows=Glucose, Cols=Ploidy. Scales=FREE.
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
  
  combined_plot <- (p1 | p2) + 
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = paste0("Posterior Check: Cell Line ", lid, " (", MODEL_NAME, ")"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  print(combined_plot)
  cat("  Plot generated for Line", lid, "\n")
}
dev.off()

# ==============================================================================
# 4. RECONSTRUCT TRANSFORMED PARAMETERS (RESTORED)
# ==============================================================================
cat("\nReconstructing effective ODE parameters (Transformed)...\n")

vars_needed <- c("mu_global", "sigma_line", "z_line", "beta_high")
draws_df <- as_draws_df(fit$draws(variables = vars_needed, inc_warmup = FALSE))

# Functions used in Stan
softcap <- function(x, cap) cap - log1p(exp(cap - x))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Define Groups
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
    # Extract hierarchical components
    mu    <- draws_df[[paste0("mu_global[", pp, "]")]]
    sigma <- draws_df[[paste0("sigma_line[", pp, "]")]]
    z     <- draws_df[[paste0("z_line[", pp, ",", l, "]")]]
    beta  <- draws_df[[paste0("beta_high[", pp, "]")]]
    
    # Linear Predictor
    raw <- mu + sigma * z + beta * h
    
    # --- MODEL SPECIFIC TRANSFORMS ---
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
  if(".chain" %in% names(draws_df)) {
    df_g$chain <- as.factor(draws_df$.chain)
  } else {
    df_g$chain <- as.factor(1)
  }
  
  reconstructed_list[[g]] <- pivot_longer(df_g, cols = all_of(param_names), names_to = "param", values_to = "value")
}

reconstructed_df <- bind_rows(reconstructed_list)
reconstructed_df$group <- factor(reconstructed_df$group, levels = groups$group_lbl)
reconstructed_df$param <- factor(reconstructed_df$param, levels = param_names)

# --- PLOT TRANSFORMED ---
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

# Plot 1: Raw Beta
p_beta_raw <- ggplot(beta_df, aes(x = value, fill = stat(x > 0))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_histogram(bins = 60, alpha = 0.7, color = "white", size = 0.1) +
  facet_wrap(~ param_name, scales = "free") +
  scale_fill_manual(values = c("TRUE"="#E41A1C", "FALSE"="#377EB8"), guide="none") +
  theme_bw() +
  labs(
    title = paste0("WGD Effect (Raw Beta Parameters) - ", MODEL_NAME),
    subtitle = "Value > 0 (Red) = WGD increases raw value. Value < 0 (Blue) = WGD decreases it.",
    x = "Raw Beta Value", y = "Density"
  )
print(p_beta_raw)

# Plot 2: Fold Change
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
# 7. PSEUDO-TRACE DIAGNOSTIC (Draw Diversity Check)
# ==============================================================================
cat("\nGenerating Pseudo-Trace Plots (results/pseudo_trace_check.pdf)...\n")

# Re-extract mu_global WITH draw indices (previous df might have filtered them)
draws_trace <- fit$draws("mu_global", inc_warmup = FALSE)
trace_df <- as_draws_df(draws_trace) %>%
  mutate(.draw = 1:n()) %>% # Ensure explicit draw index
  pivot_longer(cols = starts_with("mu_global"), names_to = "param_idx", values_to = "value")

# Map Index to Name
trace_df$idx <- as.integer(str_extract(trace_df$param_idx, "[0-9]+"))
trace_df$param_name <- factor(param_names[trace_df$idx], levels = param_names)

pdf("results/pseudo_trace_check.pdf", width = 12, height = 8)

p_trace <- ggplot(trace_df, aes(x = .draw, y = value)) +
  # Use jitter to see density if draws are identical
  geom_point(alpha = 0.4, size = 0.8, color = "#2c3e50") +
  facet_wrap(~ param_name, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Pseudo-Trace of Pathfinder Draws", 
    subtitle = "Diagnostic: 'Fuzzy Bar' = Good Approximation. 'Discrete Lines' = Disjoint/Fragile Modes.",
    x = "Draw Index (Random Sample)", 
    y = "Raw Parameter Value"
  )

print(p_trace)
dev.off()

cat("  Saved pseudo-trace to 'results/pseudo_trace_check.pdf'\n")

# ==============================================================================
# 8. PATHFINDER DIAGNOSTICS (Success/Fail Check)
# ==============================================================================
cat("\nChecking Pathfinder Convergence Diagnostics...\n")

# 1. Return Codes (Did the processes crash?)
# 0 = Success, Non-Zero = Crash/Failure
codes <- fit$return_codes()
n_fail <- sum(codes != 0)
cat(sprintf("  Total Paths Run: %d\n", length(codes)))
cat(sprintf("  Failed (Crashed): %d\n", n_fail))
if (n_fail > 0) {
  cat("  WARNING: Some paths crashed. Check prior initialization range.\n")
}

# 2. Log-Probability (lp__) Distribution
# This tells you if the 'successful' paths found the same mode.
# If you see a bimodal histogram, some paths got stuck in a bad local optima.
draws_lp <- fit$draws("lp__", inc_warmup = FALSE)
lp_df <- as_draws_df(draws_lp)

pdf("results/diagnostic_lp_check.pdf", width = 8, height = 6)

p_lp <- ggplot(lp_df, aes(x = lp__)) +
  geom_histogram(bins = 50, fill = "purple", alpha = 0.7, color = "black") +
  theme_bw() +
  labs(
    title = "Log-Probability (lp__) Distribution",
    subtitle = "Left-tail outliers indicate paths stuck in poor local optima.",
    x = "Log Probability",
    y = "Count of Draws"
  )

print(p_lp)
dev.off()

cat("  Saved lp__ diagnostic to 'results/diagnostic_lp_check.pdf'\n")