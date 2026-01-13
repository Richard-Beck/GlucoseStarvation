suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 1) SETUP
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
TARGET_LINE <- if (length(args) >= 1) args[1] else "MCF10A"
MODEL_NAME  <- if (length(args) >= 2) args[2] else "model_B1cond"

base_dir <- file.path("data", "pathfinder_1cond", MODEL_NAME, TARGET_LINE)
if (!dir.exists(base_dir)) stop("Base dir not found: ", base_dir)

out_dir <- base_dir
cat(sprintf(">>> Comparison plots for %s (%s)\n", TARGET_LINE, MODEL_NAME))
cat(sprintf("    Base dir: %s\n", base_dir))

# Parameter names in the *p_phys input* order (i.e., what logp represents)
# 1=theta, 2=kp, 3=kd, 4=kd2, 5=g50a, 6=na, 7=gamma, 8=nd, 9=yield_inv (or yield proxy), 10=m
param_names <- c("theta","kp","kd","kd2","g50a","na","gamma","nd","yield_inv","m")

state_names <- c("NL","ND","G")

# ==============================================================================
# 2) HELPERS
# ==============================================================================
safe_read_fit <- function(path) {
  if (!file.exists(path)) stop("Missing fit file: ", path)
  readRDS(path)
}

safe_read_data <- function(path) {
  if (!file.exists(path)) stop("Missing stan data file: ", path)
  readRDS(path)
}

# Extract parameter draws (logp -> p_phys = exp(logp))
extract_params <- function(fit, cond_label) {
  logp_mat <- fit$draws("logp", format = "matrix")
  if (ncol(logp_mat) != length(param_names)) {
    stop("logp has ", ncol(logp_mat), " cols; expected ", length(param_names),
         ". Check N_params / param order.")
  }
  p_phys <- exp(logp_mat)
  colnames(p_phys) <- param_names
  as_tibble(p_phys) %>% mutate(Condition = cond_label)
}

extract_sims <- function(fit, stan_data, cond_label, probs = c(0.05, 0.5, 0.95)) {
  # If y_sim wasn't saved, return empty
  if (!("y_sim" %in% fit$metadata()$stan_variables)) {
    message("    (No y_sim in fit for ", cond_label, ")")
    return(tibble())
  }
  
  mu_mat <- fit$draws("y_sim", inc_warmup = FALSE, format = "matrix")
  var_names <- colnames(mu_mat)
  if (is.null(var_names) || length(var_names) == 0) stop("No y_sim columns found.")
  
  # Quantiles column-wise: returns length(probs) x n_vars
  qs <- apply(mu_mat, 2, quantile, probs = probs, names = FALSE)
  # ensure orientation [q x var]
  if (ncol(qs) != length(var_names)) qs <- t(qs)
  
  # Parse "y_sim[w,g,s]" from names
  # Works whether cmdstanr uses "y_sim[1,2,3]" or "y_sim.1.2.3" (rare)
  idx_df <- tibble(variable = var_names) %>%
    mutate(clean = gsub("^y_sim\\[|\\]$", "", variable)) %>%
    tidyr::separate(clean, into = c("w", "g", "s"), sep = ",", convert = TRUE, remove = FALSE)
  
  # Basic sanity checks
  if (anyNA(idx_df$w) || anyNA(idx_df$g) || anyNA(idx_df$s)) {
    stop("Failed to parse some y_sim indices from names. Example names:\n",
         paste(head(var_names, 10), collapse = "\n"))
  }
  
  # Map indices -> data (SUBSET LOCAL INDEXING)
  out <- idx_df %>%
    mutate(
      lo     = as.numeric(qs[1, ]),
      median = as.numeric(qs[2, ]),
      hi     = as.numeric(qs[3, ]),
      well_idx = w,
      time     = as.numeric(stan_data$t_grid[g]),
      type     = dplyr::case_when(
        s == 1 ~ "NL",
        s == 2 ~ "ND",
        s == 3 ~ "G",
        TRUE ~ NA_character_
      ),
      Condition = cond_label,
      exp_id = as.integer(stan_data$exp_id[well_idx]),
      G0 = as.numeric(stan_data$G0_per_well[well_idx]),
      G0_lbl = factor(paste0(G0, " mM"), levels = paste0(sort(unique(stan_data$G0_per_well)), " mM")),
      well_key = paste0(cond_label, "_", well_idx)
    ) %>%
    filter(!is.na(type)) %>%
    select(well_idx, well_key, exp_id, G0_lbl, time, type, lo, median, hi, Condition)
  
  out
}


# Observed data (counts + glucose) for a subset stan_data
extract_observed <- function(stan_data, cond_label) {
  meta <- tibble(
    well_idx = seq_len(stan_data$N_wells),
    G0       = as.numeric(stan_data$G0_per_well),
    exp_id   = as.integer(stan_data$exp_id)
  ) %>%
    mutate(
      Condition = cond_label,
      G0_lbl = factor(paste0(G0, " mM"), levels = paste0(sort(unique(G0)), " mM")),
      well_key = paste0(cond_label, "_", well_idx)
    )
  
  # Counts
  df_counts <- tibble(
    well_idx = as.integer(stan_data$well_idx_count),
    time     = as.numeric(stan_data$t_grid[stan_data$grid_idx_count]),
    NL       = as.numeric(stan_data$N_obs),
    ND       = as.numeric(stan_data$D_obs)
  ) %>%
    pivot_longer(cols = c("NL","ND"), names_to = "type", values_to = "value") %>%
    left_join(meta, by = "well_idx")
  
  # Glucose: convert lum->G using fixed calibration (same as your old plotting logic)
  df_gluc <- tibble(
    well_idx = as.integer(stan_data$well_idx_gluc),
    time     = as.numeric(stan_data$t_grid[stan_data$grid_idx_gluc]),
    lum      = as.numeric(stan_data$lum_obs),
    dil      = as.numeric(stan_data$dilution)
  ) %>%
    left_join(meta %>% select(well_idx, exp_id, Condition, G0_lbl, well_key), by = "well_idx") %>%
    mutate(
      a = as.numeric(stan_data$calib_a_fixed[exp_id]),
      b = as.numeric(stan_data$calib_b_fixed[exp_id]),
      value = pmax(0, (lum - b) / (a * dil)),
      type  = "G"
    ) %>%
    select(well_idx, time, type, value, exp_id, Condition, G0_lbl, well_key)
  
  bind_rows(df_counts, df_gluc)
}

# Hill curve from median params
calc_hill <- function(df_median, G_seq) {
  kp   <- df_median$kp
  kd   <- df_median$kd
  g50a <- df_median$g50a
  na   <- df_median$na
  nd   <- df_median$nd
  gamma <- df_median$gamma
  
  g50d <- g50a * gamma
  log_G <- log(G_seq + 1e-9)
  
  actA   <- 1 / (1 + exp(-na * (log_G - log(g50a))))
  term_d <- 1 / (1 + exp(-nd * (log_G - log(g50d))))
  inhD   <- 1 - term_d
  
  tibble(
    G = G_seq,
    Growth = kp * actA,
    Death  = kd * inhD
  )
}

calc_G_threshold_lines_medianparams <- function(params_df, qs = c(0.05, 0.5, 0.95)) {
  qlab <- paste0(qs * 100, "%")
  logit_q  <- qlogis(qs)
  logit_1q <- qlogis(1 - qs)
  
  med <- params_df %>%
    group_by(Condition) %>%
    summarise(across(all_of(param_names), median), .groups = "drop") %>%
    mutate(g50d = g50a * gamma)
  
  out <- bind_rows(lapply(seq_len(nrow(med)), function(i) {
    row <- med[i, ]
    # Growth: actA(G)=q -> G = g50a * exp(logit(q)/na)
    G_grow <- row$g50a * exp(logit_q / row$na)
    # Death: inhD(G)=q -> inv_logit(...)=1-q -> G = g50d * exp(logit(1-q)/nd)
    G_dead <- row$g50d * exp(logit_1q / row$nd)
    
    tibble(
      Condition = row$Condition,
      Curve = rep(c("Growth", "Death"), each = length(qs)),
      q = rep(qlab, times = 2),
      G_line = c(G_grow, G_dead)
    )
  }))
  
  out
}




# ==============================================================================
# 3) LOAD HI/LO RESULTS USING SUBSETTED DATA FILES
# ==============================================================================
conditions <- c("hi","lo")
labels     <- c("High Ploidy","Low Ploidy")

all_params <- list()
all_sims   <- list()
all_obs    <- list()

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  lbl  <- labels[i]
  
  fit_path  <- file.path(base_dir, cond, paste0("pathfinder_", MODEL_NAME, ".Rds"))
  data_path <- file.path(base_dir, cond, "stan_data.Rds")
  
  if (!file.exists(fit_path)) {
    message("!!! Missing fit for ", cond, ": ", fit_path)
    next
  }
  if (!file.exists(data_path)) {
    message("!!! Missing stan_data for ", cond, ": ", data_path)
    next
  }
  
  cat(sprintf("... Loading %s\n", cond))
  fit_i  <- safe_read_fit(fit_path)
  data_i <- safe_read_data(data_path)
  ref_stan_data <- data_i
  
  all_params[[cond]] <- extract_params(fit_i, lbl)
  all_sims[[cond]]   <- extract_sims(fit_i, data_i, lbl)
  all_obs[[cond]]    <- extract_observed(data_i, lbl)
}

params_df <- bind_rows(all_params)
sims_df   <- bind_rows(all_sims)
obs_df    <- bind_rows(all_obs)

if (nrow(params_df) == 0) stop("No parameter draws loaded.")
if (nrow(obs_df) == 0) stop("No observed data loaded.")
if (nrow(sims_df) == 0) message("Note: no simulations loaded (calc_sim may have been 0 or y_sim not saved).")

# ==============================================================================
# 4) PLOTS
# ==============================================================================

# --- Plot 1: PPC (counts + glucose + hill-quantile bands) ---
cat("Writing comparison_ppc.pdf ...\n")
pdf(file.path(out_dir, "comparison_ppc.pdf"), width = 21, height = 10)  # wider


# no manual colors per your preference; keep defaults
# Counts panel
p_counts <- ggplot() +
  geom_ribbon(
    data = sims_df %>% filter(type != "G"),
    aes(x = time, ymin = lo, ymax = hi, group = interaction(well_key, type), fill = type),
    alpha = 0.2
  ) +
  geom_line(
    data = sims_df %>% filter(type != "G"),
    aes(x = time, y = median, group = interaction(well_key, type), color = type),
    linewidth = 0.6
  ) +
  geom_point(
    data = obs_df %>% filter(type != "G"),
    aes(x = time, y = value, shape = factor(exp_id), color = type),
    size = 1.8, alpha = 0.6
  ) +
  facet_grid(G0_lbl ~ Condition,scales="free") +
  theme_bw() +
  labs(title = paste0(TARGET_LINE, ": Cell counts (NL/ND)"), shape = "Exp")

# Glucose panel
p_gluc <- ggplot() +
  geom_ribbon(
    data = sims_df %>% filter(type == "G"),
    aes(x = time, ymin = lo, ymax = hi, group = interaction(well_key, type), fill = type),
    alpha = 0.2
  ) +
  geom_line(
    data = sims_df %>% filter(type == "G"),
    aes(x = time, y = median, group = interaction(well_key, type), color = type),
    linewidth = 0.6
  ) +
  geom_point(
    data = obs_df %>% filter(type == "G"),
    aes(x = time, y = value, shape = factor(exp_id), color = type),
    size = 1.8, alpha = 0.6
  ) +
  facet_grid(G0_lbl ~ Condition,scales="free") +
  theme_bw() +
  labs(title = paste0(TARGET_LINE, ": Glucose (inferred from lum)"), shape = "Exp")

# Threshold lines: glucose values where growth/death terms are at 5/50/95% of max
thr_df <- calc_G_threshold_lines_medianparams(params_df, qs = c(0.05, 0.5, 0.95))

p_thresh <- ggplot() +
  geom_line(
    data = sims_df %>% filter(type == "G"),
    aes(x = time, y = median, group = well_key),
    linewidth = 0.5, alpha = 0.7
  ) +
  geom_hline(
    data = thr_df,
    aes(yintercept = G_line, color = Curve, linetype = q),
    linewidth = 0.6
  ) +
  scale_y_log10() +
  facet_grid(G0_lbl ~ Condition, scales = "free_x") +
  theme_bw() +
  labs(
    title = paste0(TARGET_LINE, ": Glucose (log) with Hill thresholds (median params)"),
    y = "Glucose (mM, log scale)",
    color = "Term",
    linetype = "Target %"
  )



print(p_counts + p_gluc + p_thresh)
dev.off()

# --- Plot 2: Parameter densities ---
cat("Writing comparison_params.pdf ...\n")
pdf(file.path(out_dir, "comparison_params.pdf"), width = 12, height = 8)

long_p <- params_df %>%
  pivot_longer(cols = all_of(param_names), names_to = "Param", values_to = "Val")

print(
  ggplot(long_p, aes(x = Val, fill = Condition)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Param, scales = "free") +
    theme_bw() +
    labs(title = paste0("Parameter comparison: ", TARGET_LINE))
)
dev.off()

# --- Plot 2b: Parameter posteriors on log scale (boxplot) + bounds + prior mean ---
cat("Writing comparison_params_logbox.pdf ...\n")
pdf(file.path(out_dir, "comparison_params_logbox.pdf"), width = 12, height = 8)

bounds_df <- tibble(
  Param = param_names,
  log_lower = log(as.numeric(ref_stan_data$lower_b)),
  log_upper = log(as.numeric(ref_stan_data$upper_b)),
  log_prior_mean = as.numeric(ref_stan_data$prior_ode_mean),
  log_prior_sd   = as.numeric(ref_stan_data$prior_ode_sd)
)

long_log <- params_df %>%
  pivot_longer(cols = all_of(param_names), names_to = "Param", values_to = "Val") %>%
  mutate(logVal = log(Val)) %>%
  left_join(bounds_df, by = "Param")

p_logbox <- ggplot(long_log, aes(x = Condition, y = logVal, fill = Condition)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
  facet_wrap(~Param, scales = "free_y") +
  geom_hline(aes(yintercept = log_lower), linetype = 2, linewidth = 0.3) +
  geom_hline(aes(yintercept = log_upper), linetype = 2, linewidth = 0.3) +
  geom_hline(aes(yintercept = log_prior_mean), linetype = 3, linewidth = 0.4) +
  theme_bw() +
  labs(
    title = paste0("Posterior (log scale) + bounds + prior mean: ", TARGET_LINE),
    x = NULL, y = "log(parameter)"
  )

print(p_logbox)
dev.off()


# --- Plot 3: Hill curves (growth vs death terms) ---
cat("Writing comparison_hill.pdf ...\n")
pdf(file.path(out_dir, "comparison_hill.pdf"), width = 10, height = 6)

medians <- params_df %>%
  group_by(Condition) %>%
  summarise(across(all_of(param_names), median), .groups = "drop")

G_range <- exp(seq(log(0.001), log(25), length.out = 200))

hill_df <- bind_rows(lapply(seq_len(nrow(medians)), function(i) {
  h <- calc_hill(medians[i,], G_range)
  h$Condition <- medians$Condition[i]
  h
})) %>%
  pivot_longer(cols = c("Growth","Death"), names_to = "RateType", values_to = "Rate")

print(
  ggplot(hill_df, aes(x = G, y = Rate, color = Condition)) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~RateType, scales = "free") +
    scale_x_log10() +
    theme_bw() +
    labs(title = paste0("Physiology curves: ", TARGET_LINE),
         x = "Glucose (mM)", y = "Rate (1/h)")
)
dev.off()

cat(">>> Done.\n")


