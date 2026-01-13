library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)

# ==============================================================================
# 1. SETUP
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
TARGET_LINE <- if (length(args) >= 1) args[1] else "MCF10A"
MODEL_NAME  <- "model_B1cond"

# Input/Output Paths
base_dir  <- file.path("data", "pathfinder_1cond", MODEL_NAME, TARGET_LINE)
data_file <- "data/stan_ready_data.Rds" 
out_dir   <- base_dir 

cat(sprintf(">>> Generating Comparison Plots for Line: %s\n", TARGET_LINE))

# Load Master Data
if (!file.exists(data_file)) stop("Master data file not found.")
stan_data <- readRDS(data_file)

# Parameter Names in p_phys order (Input Parameterization)
# 1=theta, 2=kp, 3=kd, 4=kd2, 5=g50a, 6=na, 7=gamma, 8=nd, 9=yield, 10=m
param_names <- c("theta", "kp", "kd", "kd2", "g50a", "na", "gamma", "nd", "yield", "m")

# ==============================================================================
# 2. HELPER: INDICES & TRANSFORMS
# ==============================================================================

# A. Robust Data Subsetter (Returns GLOBAL indices for observations)
get_global_indices <- function(m_data, t_line, t_ploidy) {
  line_idx <- m_data$line_map[[t_line]]
  
  # Find all wells for this line
  wells_in_line <- which(m_data$line_id == line_idx)
  
  # Determine specific ploidy value
  p_vals <- m_data$ploidy_metric[wells_in_line]
  u_p    <- sort(unique(p_vals))
  
  if(length(u_p) == 1) {
    target_p <- u_p[1]
  } else {
    target_p <- if(t_ploidy == "hi") max(u_p) else min(u_p)
  }
  
  # Return global indices
  which(m_data$line_id == line_idx & m_data$ploidy_metric == target_p)
}

# B. Parameter Transformer (Raw Unconstrained -> Physical)
transform_params <- function(draws_mat, lower_vec, upper_vec) {
  lb <- log(as.numeric(lower_vec))
  ub <- log(as.numeric(upper_vec))
  
  n_d <- nrow(draws_mat)
  n_p <- ncol(draws_mat)
  
  phys <- matrix(NA, nrow=n_d, ncol=n_p)
  
  for(j in 1:n_p) {
    inv_logit_x <- 1 / (1 + exp(-draws_mat[,j]))
    phys[,j] <- exp(lb[j] + (ub[j] - lb[j]) * inv_logit_x)
  }
  return(phys)
}

# ==============================================================================
# 3. EXTRACTION LOOP (HI & LO)
# ==============================================================================
conditions <- c("hi", "lo")
labels     <- c("High Ploidy", "Low Ploidy")

res_sims   <- list()
res_params <- list()

for(i in 1:2) {
  cond <- conditions[i]
  lbl  <- labels[i]
  fpath <- file.path(base_dir, cond, paste0("pathfinder_", MODEL_NAME, ".Rds"))
  
  if(!file.exists(fpath)) {
    cat(sprintf("!!! Missing fit file for %s. Skipping.\n", cond))
    next
  }
  
  cat(sprintf("... Processing %s ...\n", cond))
  fit <- readRDS(fpath)
  
  # --- 1. Extract & Transform Parameters ---
  N_p <- length(stan_data$lower_b)
  p_raw <- fit$draws("p_raw", format="matrix")[, 1:N_p, drop=FALSE]
  
  p_phys_mat <- transform_params(p_raw, stan_data$lower_b, stan_data$upper_b)
  colnames(p_phys_mat) <- param_names[1:N_p]
  
  df_p <- as.data.frame(p_phys_mat)
  df_p$Condition <- lbl
  res_params[[i]] <- df_p
  
  # --- 2. Extract Simulations (Robust String Parsing) ---
  if("y_sim" %in% fit$metadata()$stan_variables) {
    
    # Extract as flat matrix [draws x variables]
    mu_mat <- fit$draws("y_sim", inc_warmup = FALSE, format = "matrix")
    var_names <- colnames(mu_mat)
    
    cat("    Calculating quantiles for", length(var_names), "variables...\n")
    
    # Calculate quantiles column-wise (returns 3 x n_vars matrix)
    stats <- apply(mu_mat, 2, quantile, probs = c(0.05, 0.5, 0.95))
    
    # Create Index Map from names (e.g., "y_sim[1,10,2]")
    idx_df <- tibble(variable = var_names) %>%
      mutate(clean = gsub("y_sim\\[|\\]", "", variable)) %>%
      separate(clean, c("w_local", "t_idx", "s_idx"), sep = ",", convert = TRUE)
    
    # Map Global Indices
    global_idx_map <- get_global_indices(stan_data, TARGET_LINE, cond)
    
    # Combine
    sim_df <- idx_df %>%
      mutate(
        lo     = stats[1, ],
        median = stats[2, ],
        hi     = stats[3, ],
        
        # Map indices
        well_idx = global_idx_map[w_local],
        time     = stan_data$t_grid[t_idx],
        type     = case_when(s_idx==1 ~ "NL", s_idx==2 ~ "ND", s_idx==3 ~ "G", TRUE ~ "Other"),
        Condition = lbl
      ) %>%
      filter(type %in% c("NL", "ND", "G")) %>%
      select(well_idx, time, type, lo, median, hi, Condition)
    
    res_sims[[i]] <- sim_df
    
  } else {
    cat("    (No y_sim found in fit object)\n")
  }
}

all_params <- bind_rows(res_params)
all_sims   <- bind_rows(res_sims)

if(nrow(all_params) == 0) stop("No data loaded.")

# ==============================================================================
# 4. PREPARE OBSERVED DATA
# ==============================================================================
# Get global indices for ALL wells in this line
line_idx_glob <- stan_data$line_map[[TARGET_LINE]]
all_wells <- which(stan_data$line_id == line_idx_glob)

# Create Condition Map
cond_map <- rep(NA, stan_data$N_wells)
cond_map[get_global_indices(stan_data, TARGET_LINE, "hi")] <- "High Ploidy"
cond_map[get_global_indices(stan_data, TARGET_LINE, "lo")] <- "Low Ploidy"

# Obs: Metadata
obs_meta <- data.frame(
  well_idx = all_wells,
  G0       = stan_data$G0_per_well[all_wells],
  exp_id   = stan_data$exp_id[all_wells],
  Condition = cond_map[all_wells]
) %>% 
  filter(!is.na(Condition)) %>%
  mutate(G0_lbl = factor(paste0(G0, " mM"), levels=paste0(sort(unique(G0)), " mM")))

# Obs: Counts (NL, ND)
df_counts <- data.frame(
  well_idx = stan_data$well_idx_count,
  time     = stan_data$t_grid[stan_data$grid_idx_count],
  NL       = stan_data$N_obs,
  ND       = stan_data$D_obs
) %>%
  filter(well_idx %in% obs_meta$well_idx) %>%
  pivot_longer(cols=c("NL","ND"), names_to="type", values_to="value")

# Obs: Glucose
df_gluc <- data.frame(
  well_idx = stan_data$well_idx_gluc,
  time     = stan_data$t_grid[stan_data$grid_idx_gluc],
  lum      = stan_data$lum_obs,
  dil      = stan_data$dilution
) %>%
  filter(well_idx %in% obs_meta$well_idx) %>%
  left_join(obs_meta %>% select(well_idx, exp_id), by="well_idx") %>%
  mutate(
    a = stan_data$calib_a_fixed[exp_id],
    b = stan_data$calib_b_fixed[exp_id],
    value = pmax(0, (lum - b)/(a*dil)),
    type  = "G"
  ) %>%
  select(well_idx, time, value, type)

obs_all <- bind_rows(df_counts, df_gluc) %>%
  left_join(obs_meta, by="well_idx")

# Join Metadata to Sims
all_sims <- all_sims %>%
  left_join(obs_meta %>% select(well_idx, G0_lbl, exp_id), by="well_idx")


# ==============================================================================
# 5. PLOT GENERATION
# ==============================================================================

# --- Plot 1: PPC (Trajectories) ---
cat("Generating PPC Plot...\n")
pdf(file.path(out_dir, "comparison_ppc.pdf"), width=16, height=12)

cols <- c("NL"="#009E73", "ND"="#D55E00", "G"="#0072B2")
fills <- c("NL"="#009E73", "ND"="#D55E00", "G"="#0072B2")

# Cell Counts
p1 <- ggplot() +
  geom_ribbon(data=all_sims %>% filter(type!="G"), 
              aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
  geom_line(data=all_sims %>% filter(type!="G"), 
            aes(x=time, y=median, color=type, group=interaction(well_idx, type)), size=0.8) +
  geom_point(data=obs_all %>% filter(type!="G"), 
             aes(x=time, y=value, color=type, shape=factor(exp_id)), size=2, alpha=0.6) +
  facet_grid(G0_lbl ~ Condition) +
  scale_color_manual(values=cols) + scale_fill_manual(values=fills) +
  labs(title=paste0(TARGET_LINE, ": Cell Counts"), shape="Exp") + theme_bw() +
  theme(legend.position="bottom")

# Glucose
p2 <- ggplot() +
  geom_ribbon(data=all_sims %>% filter(type=="G"), 
              aes(x=time, ymin=lo, ymax=hi, fill=type, group=interaction(well_idx, type)), alpha=0.2) +
  geom_line(data=all_sims %>% filter(type=="G"), 
            aes(x=time, y=median, color=type, group=interaction(well_idx, type)), size=0.8) +
  geom_point(data=obs_all %>% filter(type=="G"), 
             aes(x=time, y=value, color=type, shape=factor(exp_id)), size=2, alpha=0.6) +
  facet_grid(G0_lbl ~ Condition) +
  scale_color_manual(values=cols) + scale_fill_manual(values=fills) +
  labs(title=paste0(TARGET_LINE, ": Glucose"), shape="Exp") + theme_bw() +
  theme(legend.position="bottom")

print(p1 | p2)
dev.off()


# --- Plot 2: Parameter Densities ---
cat("Generating Parameter Plot...\n")
pdf(file.path(out_dir, "comparison_params.pdf"), width=12, height=8)

long_p <- all_params %>%
  pivot_longer(cols=-Condition, names_to="Param", values_to="Val")

print(
  ggplot(long_p, aes(x=Val, fill=Condition)) +
    geom_density(alpha=0.5) +
    facet_wrap(~Param, scales="free") +
    scale_fill_manual(values=c("Low Ploidy"="#377EB8", "High Ploidy"="#E41A1C")) +
    theme_bw() +
    labs(title=paste0("Parameter Comparison: ", TARGET_LINE))
)
dev.off()


# --- Plot 3: Hill Functions (Physiology) ---
cat("Generating Physiology Plot...\n")
pdf(file.path(out_dir, "comparison_hill.pdf"), width=10, height=6)

# Helper to calc rate curves
calc_hill <- function(df_median, G_seq) {
  # Map variables
  kp   <- df_median$kp
  kd   <- df_median$kd
  g50a <- df_median$g50a
  na   <- df_median$na
  nd   <- df_median$nd
  g50d <- g50a * df_median$gamma # derived
  
  log_G <- log(G_seq + 1e-9)
  
  actA   <- 1 / (1 + exp(-na * (log_G - log(g50a))))
  term_d <- 1 / (1 + exp(-nd * (log_G - log(g50d))))
  inhD   <- 1 - term_d
  
  data.frame(
    G = G_seq,
    Growth = kp * actA,
    Death  = kd * inhD
  )
}

# Calculate Medians
medians <- all_params %>%
  group_by(Condition) %>%
  summarise(across(everything(), median))

G_range <- exp(seq(log(0.01), log(25), length.out=100))
hill_list <- list()

for(i in 1:nrow(medians)) {
  res <- calc_hill(medians[i,], G_range)
  res$Condition <- medians$Condition[i]
  hill_list[[i]] <- res
}

hill_df <- bind_rows(hill_list) %>%
  pivot_longer(cols=c("Growth","Death"), names_to="RateType", values_to="Rate")

print(
  ggplot(hill_df, aes(x=G, y=Rate, color=Condition)) +
    geom_line(linewidth=1.2) +
    facet_wrap(~RateType, scales="free") +
    scale_x_log10() +
    scale_color_manual(values=c("Low Ploidy"="#377EB8", "High Ploidy"="#E41A1C")) +
    theme_bw() +
    labs(title=paste0("Physiological Functions: ", TARGET_LINE), x="Glucose (mM)", y="Rate (1/h)")
)
dev.off()

cat(">>> Done.\n")
