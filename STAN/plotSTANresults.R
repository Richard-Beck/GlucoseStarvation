library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# ==============================================================================
# 1. SETUP & LOAD
# ==============================================================================
cat("Loading Data and Fit...\n")
stan_data <- readRDS("data/stan_ready_data.Rds")
fit_file  <- "results/fit_model_B.Rds"

if (!file.exists(fit_file)) stop("Fit file not found!")
fit <- readRDS(fit_file)

if (!dir.exists("results")) dir.create("results")

# ==============================================================================
# 2. POSTERIOR PARAMETER DISTRIBUTIONS (Single Page)
# ==============================================================================
cat("Extracting hierarchical parameters for reconstruction...\n")

vars_needed <- c("mu_global", "sigma_line", "z_line", "beta_high")
draws_df <- as_draws_df(fit$draws(variables = vars_needed, inc_warmup = FALSE))

softcap <- function(x, cap) cap - log1p(exp(cap - x))
cap_log_main <- 40.0
cap_log_hill <- 6.0
par10_names  <- c("theta","kp","kd","kd2","g50a","na","g50d","nd","v1","v2")

groups <- tibble(
  line_id = stan_data$line_id, 
  high = stan_data$is_high_ploidy
) %>%
  distinct(line_id, high) %>%
  arrange(line_id, high) %>%
  mutate(group_lbl = paste0("Line ", line_id, " | ", ifelse(high==1, "High", "Low")))

cat("Reconstructing effective ODE parameters...\n")
reconstructed_list <- list()
N_draws <- nrow(draws_df)

for(g in 1:nrow(groups)) {
  l <- groups$line_id[g]
  h <- groups$high[g]
  lbl <- groups$group_lbl[g]
  
  p_mat <- matrix(NA_real_, nrow = N_draws, ncol = 10)
  colnames(p_mat) <- par10_names
  
  for(pp in 1:10) {
    mu    <- draws_df[[paste0("mu_global[", pp, "]")]]
    sigma <- draws_df[[paste0("sigma_line[", pp, "]")]]
    z     <- draws_df[[paste0("z_line[", pp, ",", l, "]")]]
    beta  <- draws_df[[paste0("beta_high[", pp, "]")]]
    
    raw <- mu + sigma * z + beta * h
    
    if(pp %in% c(6, 8)) {
      p_mat[, pp] <- 1.0 + exp(softcap(raw, cap_log_hill))
    } else {
      p_mat[, pp] <- exp(softcap(raw, cap_log_main))
    }
  }
  
  df_g <- as.data.frame(p_mat)
  df_g$chain <- as.factor(draws_df$.chain)
  df_g$group <- lbl
  
  reconstructed_list[[g]] <- pivot_longer(df_g, cols = all_of(par10_names), names_to = "param", values_to = "value")
}

reconstructed_df <- bind_rows(reconstructed_list)
reconstructed_df$group <- factor(reconstructed_df$group, levels = groups$group_lbl)
reconstructed_df$param <- factor(reconstructed_df$param, levels = par10_names)

cat("Generating Parameter Distribution Plot...\n")
plot_height <- max(8, length(unique(reconstructed_df$group)) * 2.5)

pdf("results/final_posterior_parameters.pdf", width = 20, height = plot_height)

p <- ggplot(reconstructed_df, aes(x = value, fill = chain)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.5, color = NA) +
  facet_grid(group ~ param, scales = "free") +
  scale_x_log10() +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.text.x = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  ) +
  labs(
    title = "Posterior Parameter Distributions by Cell Line",
    subtitle = "Reconstructed Effective Parameters (Sampling Phase) | Colored by Chain",
    x = "Value (Log Scale)", 
    y = "Count",
    fill = "Chain"
  )

print(p)
dev.off()


# ==============================================================================
# 3. POSTERIOR PREDICTIVE CHECKS (PPC)
# ==============================================================================
cat("\nExtracting y_sim for PPC...\n")
draws_sim <- fit$draws("y_sim", inc_warmup = FALSE)

cat("Calculating Credible Intervals...\n")
summ_sim <- summarize_draws(
  draws_sim, 
  median = function(x) quantile(x, 0.5, names = FALSE),
  lo     = function(x) quantile(x, 0.025, names = FALSE),
  hi     = function(x) quantile(x, 0.975, names = FALSE)
)

summ_clean <- summ_sim %>%
  mutate(clean = gsub("y_sim\\[|\\]", "", variable)) %>%
  separate(clean, c("w", "t", "s"), sep = ",", convert = TRUE) %>%
  mutate(
    well_idx = w,
    time     = stan_data$t_grid[t],
    type     = case_when(s==1 ~ "NL", s==2 ~ "ND", s==3 ~ "G")
  ) %>%
  select(well_idx, time, type, median, lo, hi)

meta_df <- data.frame(
  well_idx = 1:stan_data$N_wells,
  line_id  = stan_data$line_id,
  is_high  = stan_data$is_high_ploidy,
  G0       = stan_data$G0_per_well,
  exp_id   = stan_data$exp_id
)
meta_df$ploidy_lbl <- ifelse(meta_df$is_high == 1, "High Ploidy", "Low Ploidy")
meta_df$G0_lbl     <- factor(paste0(meta_df$G0, " mM"), levels = paste0(sort(unique(meta_df$G0)), " mM"))

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

cat("Generating PPC Plot...\n")
pdf("results/final_ppc_check.pdf", width = 12, height = 10)

scale_factor <- 1000 
sim_scaled <- sim_all %>%
  mutate(
    median = ifelse(type == "G", median * scale_factor, median),
    lo     = ifelse(type == "G", lo * scale_factor, lo),
    hi     = ifelse(type == "G", hi * scale_factor, hi)
  )
obs_scaled <- obs_all %>%
  mutate(
    value = ifelse(type == "G", value * scale_factor, value)
  )

unique_lines <- sort(unique(meta_df$line_id))
for(lid in unique_lines) {
  p_sim <- sim_scaled %>% filter(line_id == lid)
  p_obs <- obs_scaled %>% filter(line_id == lid)
  if(nrow(p_sim) == 0) next
  
  cols <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#E69F00")
  fills <- c("NL" = "#009E73", "ND" = "#D55E00", "G" = "#E69F00")
  
  p <- ggplot() +
    geom_ribbon(data = p_sim, aes(x = time, ymin = lo, ymax = hi, fill = type, group = interaction(well_idx, type)), alpha = 0.2) +
    geom_line(data = p_sim, aes(x = time, y = median, color = type, group = interaction(well_idx, type)), linewidth = 0.8) +
    geom_point(data = p_obs, aes(x = time, y = value, color = type, shape = type), size = 1.8, alpha = 0.6) +
    facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
    scale_y_continuous(name = "Cell Count (N/D)", sec.axis = sec_axis(~ . / scale_factor, name = "Glucose (mM)")) +
    scale_color_manual(values = cols, labels = c("G" = "Glucose", "ND" = "Dead", "NL" = "Live")) +
    scale_fill_manual(values = fills, labels = c("G" = "Glucose", "ND" = "Dead", "NL" = "Live")) +
    scale_shape_manual(values = c("NL" = 16, "ND" = 17, "G" = 15), labels = c("G" = "Glucose", "ND" = "Dead", "NL" = "Live")) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"),
          axis.title.y.right = element_text(color = "#E69F00"), axis.text.y.right = element_text(color = "#E69F00")) +
    labs(title = paste("Posterior Predictive Check: Cell Line", lid), subtitle = "Ribbons = 95% CI | Lines = Median | Points = Data", x = "Time (hours)")
  print(p)
  cat("  Plot generated for Line", lid, "\n")
}
dev.off()


# ==============================================================================
# 4. CHAIN DIAGNOSTICS (NEW)
# ==============================================================================
cat("\nExtracting Sampler Diagnostics...\n")
sampler_diag <- fit$sampler_diagnostics(format = "df")

# Summarize metrics per chain
chain_summary <- sampler_diag %>%
  group_by(.chain) %>%
  summarise(
    divergences = sum(divergent__),
    div_rate    = mean(divergent__),
    max_tree_hit = sum(treedepth__ >= 12), # Assuming max_treedepth=12 from fitSTAN
    max_tree_rate = mean(treedepth__ >= 12),
    mean_stepsize = mean(stepsize__),
    mean_energy = mean(energy__),
    .groups = "drop"
  )

print(chain_summary)

cat("Generating Chain Diagnostics Plot...\n")
pdf("results/final_chain_diagnostics.pdf", width = 12, height = 8)

# 1. Divergence Rate (Bar Plot)
p_div <- ggplot(chain_summary, aes(x = factor(.chain), y = div_rate, fill = factor(.chain))) +
  geom_col() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Divergence Rate per Chain", 
       subtitle = "Target: 0%. High divergence indicates pathology.",
       x = "Chain", y = "Divergence Rate", fill = "Chain")

# 2. Max Treedepth Hits (Bar Plot)
p_tree <- ggplot(chain_summary, aes(x = factor(.chain), y = max_tree_rate, fill = factor(.chain))) +
  geom_col() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Max Treedepth Hit Rate per Chain", 
       subtitle = "Target: 0%. Hitting limit implies inefficiency (step size too small).",
       x = "Chain", y = "Hit Rate (Max Depth)", fill = "Chain")

# 3. Step Size Distribution (Violin Plot)
p_step <- ggplot(sampler_diag, aes(x = factor(.chain), y = stepsize__, fill = factor(.chain))) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Step Size Distribution", 
       subtitle = "Should be consistent across chains. Tiny steps = stiff curvature.",
       x = "Chain", y = "Step Size", fill = "Chain")

# 4. Energy Distribution (Violin Plot)
# Overlapping energy distributions indicate chains explored same physics
p_energy <- ggplot(sampler_diag, aes(x = factor(.chain), y = energy__, fill = factor(.chain))) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Hamiltonian Energy Distribution", 
       subtitle = "Distributions should overlap significantly.",
       x = "Chain", y = "Energy", fill = "Chain")

# 5. Acceptance Rate vs Treedepth (Scatter - Diagnostics in action)
p_accept <- ggplot(sampler_diag, aes(x = factor(treedepth__), y = accept_stat__, fill = factor(.chain))) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Acceptance Rate vs Tree Depth",
       subtitle = "High depth usually correlates with lower acceptance if struggling.",
       x = "Tree Depth", y = "Acceptance Stat", fill = "Chain")

print(p_div)
print(p_tree)
print(p_step)
print(p_energy)
print(p_accept)

dev.off()
cat("Done. Check 'results/final_chain_diagnostics.pdf'.\n")
