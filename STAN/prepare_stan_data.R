library(dplyr)
library(data.table)
library(stringr)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
glu_map <- c(
  'MCF10A-ctrl.A00-IncucyteRawDataLiveDead-varyGlucose-241015' = NULL,
  'MCF10A-hras.A00-IncucyteRawDataLiveDead-varyGlucose-241015' = NULL,
  'MCF10A.A00b-IncucyteRawDataLiveDead-varyGlucose-241015'     = "MCF10A.MDA-MB-231",
  'MCF10A.A00c-IncucyteRawDataLiveDead-varyGlucose-241015'     = "MCF10A.MDA-MB-231",
  'MDA-MB-231.B00-IncucyteRawDataLiveDead-varyGlucose-250213'  = "MCF10A.MDA-MB-231",
  'SNU668.A00-IncucyteRawDataLiveDead-varyGlucose-250324'      = "SUM-159-chem.SNU668",
  'SUM-159-chem.M00b-IncucyteRawDataLiveDead-varyGlucose'      = "SUM-159-chem.SNU668",
  'SUM-159-fuse.C00-IncucyteRawDataLiveDead-varyGlucose'       = NULL,
  'SUM-159-fuse.I00-IncucyteRawDataLiveDead-varyGlucose'       = "SUM-159-fuse"
)

ploidy_map <- list(
  'SUM-159-chem' = c("2N"=0, "4N"=1),
  'SNU668'       = c("high"=1, "low"=0),
  'MCF10A'       = c("4N"=1, "2N"=0),
  'MDA-MB-231'   = c("parental"=1, "3N"=0),
  'SUM-159-fuse' = c("2N"=0, "4N"=1)
)

# ==============================================================================
# 2. LOAD & CLEAN COUNT DATA
# ==============================================================================
counts_raw <- readRDS("data/counts/corrected.Rds")

counts <- counts_raw %>%
  mutate(exp_key = paste0(cellLine, ".", experiment)) %>%
  filter(exp_key %in% names(glu_map)) %>%
  filter(!is.null(glu_map[exp_key]))

counts$is_high_ploidy <- NA_integer_
unique_lines <- unique(counts$cellLine)
for (cline in unique_lines) {
  map <- ploidy_map[[cline]]
  idx <- counts$cellLine == cline
  counts$is_high_ploidy[idx] <- map[counts$ploidy[idx]]
}

counts$N_live <- counts$N - counts$D

counts <- counts %>%
  rename(G0 = glucose) %>%
  select(cellLine, experiment, is_high_ploidy, G0, hours, N_live, D) %>%
  arrange(cellLine, experiment, G0, hours)

# ==============================================================================
# 3. LOAD GLUCOSE DATA
# ==============================================================================
message("Loading Glucose Data...")

unique_folders <- unique(glu_map[!sapply(glu_map, is.null)])
calib_store <- list()
gluc_store  <- list()

for (folder in unique_folders) {
  path_base <- file.path("data/glucose/processed", folder)
  
  calib_df <- fread(file.path(path_base, "calibration.csv"))
  calib_df$calib_batch <- folder
  calib_store[[folder]] <- calib_df[, c("G", "Lum", "calib_batch")]
  
  xg <- fread(file.path(path_base, "data.csv"))
  xg$is_high_ploidy <- ifelse(xg$ploidy == "high", 1, ifelse(xg$ploidy == "low", 0, NA_integer_))
  xg_long <- reshape2::melt(xg, measure.vars = paste0("R", 1:3), variable.name = "Replicate", value.name = "lum")
  xg_long$hours <- xg_long$Day * 24
  xg_long$calib_batch <- folder
  xg_long$dilution <- 1000 / xg_long$`Dilution Factor`
  
  gluc_store[[folder]] <- xg_long %>%
    select(CellLine, is_high_ploidy, G0, hours, lum, dilution, calib_batch)
}

all_glucose_raw <- do.call(rbind, gluc_store)
all_calib_raw   <- do.call(rbind, calib_store)

# ==============================================================================
# 4. CONSTRUCT STAN STRUCTURE
# ==============================================================================
conditions <- counts %>%
  distinct(cellLine, experiment, is_high_ploidy, G0) %>%
  mutate(well_id = row_number())

N_wells <- nrow(conditions)
G0_vector <- numeric(N_wells)
line_ids  <- integer(N_wells)
exp_ids   <- integer(N_wells)
high_p    <- integer(N_wells)

u_lines <- unique(conditions$cellLine)
line_map_int <- setNames(seq_along(u_lines), u_lines)

u_calibs <- unique(all_calib_raw$calib_batch)
calib_map_int <- setNames(seq_along(u_calibs), u_calibs)

max_time <- max(counts$hours, na.rm = TRUE)
t_grid <- seq(0, max_time + 8, by = 8)
get_grid_idx <- function(t) { round(t / 8) + 1 }

obs_counts_list <- list()
obs_gluc_list   <- list()

for (i in 1:N_wells) {
  cond <- conditions[i, ]
  
  line_ids[i]  <- line_map_int[[cond$cellLine]]
  high_p[i]    <- cond$is_high_ploidy
  exp_ids[i]   <- calib_map_int[[glu_map[[paste0(cond$cellLine, ".", cond$experiment)]]]]
  G0_vector[i] <- as.numeric(cond$G0)
  
  sub_c <- counts %>%
    filter(cellLine == cond$cellLine, experiment == cond$experiment,
           is_high_ploidy == cond$is_high_ploidy, G0 == cond$G0)
  sub_c$well_idx <- i
  sub_c$grid_idx <- get_grid_idx(sub_c$hours)
  obs_counts_list[[i]] <- sub_c %>% select(well_idx, grid_idx, N_live, D)
  
  sub_g <- all_glucose_raw %>%
    filter(calib_batch == glu_map[[paste0(cond$cellLine, ".", cond$experiment)]],
           CellLine == cond$cellLine, G0 == cond$G0) %>%
    filter(is_high_ploidy == cond$is_high_ploidy | is.na(is_high_ploidy))
  
  if (nrow(sub_g) > 0) {
    sub_g$well_idx <- i
    sub_g$grid_idx <- get_grid_idx(sub_g$hours)
    obs_gluc_list[[i]] <- sub_g %>% select(well_idx, grid_idx, lum, dilution)
  }
}

long_counts <- do.call(rbind, obs_counts_list)
long_gluc   <- do.call(rbind, obs_gluc_list)

all_calib_raw$exp_idx <- calib_map_int[all_calib_raw$calib_batch]

# ==============================================================================
# 5. CALCULATE DATA-DRIVEN PRIORS
# ==============================================================================
idx_t0 <- which(long_counts$grid_idx == 1)
prior_mu_N0 <- mean(log(long_counts$N_live[idx_t0] + 1e-9))
prior_mu_D0 <- mean(log(long_counts$D[idx_t0] + 1e-9))

df_calib <- data.frame(G = all_calib_raw$G, Lum = all_calib_raw$Lum)
min_G <- min(df_calib$G); b_est <- median(df_calib$Lum[df_calib$G == min_G])
max_G <- max(df_calib$G); max_Lum <- median(df_calib$Lum[df_calib$G == max_G])
a_est <- (max_Lum - b_est) / (max_G - min_G)

prior_calib_a_mean <- max(10.0, as.numeric(a_est))
prior_calib_b_mean <- max(10.0, as.numeric(b_est))

cat("Calculated Priors:\n")
cat(sprintf("  N0: %.2f, D0: %.2f\n", prior_mu_N0, prior_mu_D0))
cat(sprintf("  Calib a_mean: %.2f, b_mean: %.2f\n", prior_calib_a_mean, prior_calib_b_mean))

# ==============================================================================
# 6. FIXED CALIBRATION PER EXP (log-fit)  -> calib_a_fixed / calib_b_fixed
# ==============================================================================
cat("Estimating fixed calibration (log-fit) per experiment...\n")

N_exps <- length(u_calibs)
a_fix <- rep(NA_real_, N_exps)
b_fix <- rep(NA_real_, N_exps)
sdlog_fix <- rep(NA_real_, N_exps)

cal_df <- data.frame(
  e   = all_calib_raw$exp_idx,
  G   = all_calib_raw$G,
  Lum = all_calib_raw$Lum
)

costf <- function(pars, x) {
  a <- abs(pars[1]); b <- abs(pars[2])
  mu <- x$G * a + b
  if (any(!is.finite(mu)) || any(mu <= 0) || any(!is.finite(x$Lum)) || any(x$Lum <= 0)) return(1e12)
  err <- sum((log(x$Lum) - log(mu))^2)
  if (!is.finite(err)) 1e12 else err
}

for (e in 1:N_exps) {
  sub <- cal_df[cal_df$e == e & is.finite(cal_df$G) & is.finite(cal_df$Lum) & cal_df$G >= 0 & cal_df$Lum > 0, , drop = FALSE]
  
  if (nrow(sub) < 2) {
    a_fix[e] <- max(1e-6, prior_calib_a_mean)
    b_fix[e] <- max(1e-6, prior_calib_b_mean)
    sdlog_fix[e] <- 0.1
    next
  }
  
  fit0 <- try(lm(Lum ~ G, data = sub), silent = TRUE)
  a0 <- if (!inherits(fit0, "try-error")) max(1e-6, as.numeric(coef(fit0)["G"])) else 1.0
  b0 <- max(1e-6, min(sub$Lum, na.rm = TRUE))
  
  opt <- optim(c(a0, b0), costf, x = sub, method = "Nelder-Mead", control = list(maxit = 5000))
  a_hat <- max(1e-6, abs(opt$par[1]))
  b_hat <- max(1e-6, abs(opt$par[2]))
  
  mu_hat <- sub$G * a_hat + b_hat
  sd_hat <- sd(log(sub$Lum) - log(mu_hat))
  if (!is.finite(sd_hat) || sd_hat <= 0) sd_hat <- 0.1
  
  a_fix[e] <- a_hat
  b_fix[e] <- b_hat
  sdlog_fix[e] <- sd_hat
}

# ==============================================================================
# 6b. PRINCIPLED ODE PRIORS (Scale Estimation)
# ==============================================================================
# 1. theta (Carrying Capacity): Should be slightly above max observed count
max_N_obs <- max(long_counts$N_live, na.rm=TRUE)
prior_theta <- log(max_N_obs * 1.5) 

# 2. kp (Proliferation): Doubling time ~24h => ln(2)/24 ~ 0.029
prior_kp <- log(0.03)

# 3. kd (Death): Usually 10-20% of growth rate?
prior_kd <- log(0.005)

# 4. kd2 (Density dependent): Hard to guess, keep small relative to theta
prior_kd2 <- log(0.001)

# 5. g50a (Glucose sensitivity): mM range. ~0.5 mM is typical Km
prior_g50a <- log(0.5)

# 6. na (Hill Coeff): Prior is on raw parameter for 1 + exp(raw).
# We want na ~ 2 or 3. 1 + exp(0.7) ~ 3. Center at 0.
prior_na <- 0.0

# 7. g50d (Death sensitivity): Similar to g50a
prior_g50d <- log(0.5)

# 8. nd (Hill):
prior_nd <- 0.0

# 9 & 10. v1, v2 (Consumption): 
# Estimate: (Max Glucose) / (Avg Cells * Duration)
# 25 mM / (5000 cells * 100 hours) = 5e-5
prior_v <- log(5e-5)

# Construct the vectors
# Map: 1:theta, 2:kp, 3:kd, 4:kd2, 5:g50a, 6:na, 7:g50d, 8:nd, 9:v1, 10:v2
prior_ode_mean <- c(
  prior_theta, 
  prior_kp, 
  prior_kd, 
  prior_kd2, 
  prior_g50a, 
  prior_na, 
  prior_g50d, 
  prior_nd, 
  prior_v, 
  prior_v
)

# Set widths (SDs). 
# Use 1.0 for most (allows +/- 2.7x fold change at 1 sigma).
# Use 0.5 for tight biological constraints (like max count).
prior_ode_sd <- c(
  0.5, # theta: fairly certain about scale (it's N_obs)
  1.0, # kp
  1.0, # kd
  2.0, # kd2: very uncertain
  1.5, # g50: uncertain
  1.0, # na
  1.5, # g50d
  1.0, # nd
  1.5, # v1: uncertain
  1.5  # v2
)

cat("Computed Principled ODE Priors:\n")
print(data.frame(param=c("theta","kp","kd","kd2","g50a","na","g50d","nd","v1","v2"), 
                 log_mean=round(prior_ode_mean,2), 
                 natural_scale=round(exp(prior_ode_mean), 6)))

# ==============================================================================
# 7. SAVE COMPLETE OBJECT
# ==============================================================================
stan_data <- list(
  N_wells     = N_wells,
  N_lines     = length(u_lines),
  N_exps      = N_exps,
  
  N_obs_count = nrow(long_counts),
  N_obs_gluc  = nrow(long_gluc),
  N_obs_calib = nrow(all_calib_raw),
  
  N_grid      = length(t_grid),
  t_grid      = t_grid,
  
  line_id     = line_ids,
  is_high_ploidy = high_p,
  exp_id      = exp_ids,
  
  # Observations
  well_idx_count = long_counts$well_idx,
  grid_idx_count = long_counts$grid_idx,
  N_obs          = long_counts$N_live,
  D_obs          = long_counts$D,
  
  well_idx_gluc  = long_gluc$well_idx,
  grid_idx_gluc  = long_gluc$grid_idx,
  lum_obs        = long_gluc$lum,
  dilution       = long_gluc$dilution,
  
  calib_exp_idx  = all_calib_raw$exp_idx,
  calib_G        = all_calib_raw$G,
  calib_Lum      = all_calib_raw$Lum,
  
  G0_per_well    = G0_vector,
  grainsize      = 1,
  
  # Priors used in Stan (IC priors)
  prior_mu_N0_mean = prior_mu_N0,
  prior_mu_N0_sd   = 0.5,
  prior_mu_D0_mean = prior_mu_D0,
  prior_mu_D0_sd   = 1.0,
  
  prior_ode_mean   = prior_ode_mean,
  prior_ode_sd     = prior_ode_sd,
  
  # Keep these if your Stan still expects them (even if not used now)
  prior_calib_a_mean = prior_calib_a_mean,
  prior_calib_a_sd   = 200,
  prior_calib_b_mean = prior_calib_b_mean,
  prior_calib_b_sd   = 100,
  
  # Fixed calibration per experiment (required by your Stan data block)
  calib_a_fixed = a_fix,
  calib_b_fixed = b_fix
  

)

saveRDS(stan_data, "data/stan_ready_data.Rds")
cat("Data saved to 'data/stan_ready_data.Rds'.\n")
cat(sprintf("Fixed calibration saved: calib_a_fixed (%d), calib_b_fixed (%d)\n", length(a_fix), length(b_fix)))
