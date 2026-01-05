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
  'SUM-159-fuse.C00-IncucyteRawDataLiveDead-varyGlucose'       = "SUM-159-fuse-gluccell",
  'SUM-159-fuse.I00-IncucyteRawDataLiveDead-varyGlucose'       = "SUM-159-fuse"
)

# CHANGED: Updated to numeric values for continuous variable
# Based on user notes:
# SUM159: ~2N vs ~4N
# MCF10A: ~2N vs 4N
# SNU-668: ~2.6N vs ~5.2N
# MDA-MB-231: ~3.5N vs ~4.4N
ploidy_map <- list(
  'SUM-159-chem' = c("2N"=2.0, "4N"=4.0),
  'SNU668'       = c("high"=5.2, "low"=2.6),
  'MCF10A'       = c("4N"=4.0, "2N"=2.0),
  'MDA-MB-231'   = c("parental"=4.4, "3N"=3.5), # Assuming 'parental' is the higher ploidy ~4.4N based on context
  'SUM-159-fuse' = c("2N"=2.0, "4N"=4.0)
)

if (!dir.exists("data")) dir.create("data")

# ==============================================================================
# 2. LOAD & CLEAN COUNT DATA
# ==============================================================================
as_corrected_counts <- function(x){
  is_dt <- inherits(x, "data.table")
  y <- if (is_dt && requireNamespace("data.table", quietly=TRUE)) data.table::copy(x) else x
  
  if (all(c("time_numeric","N","D") %in% names(y))) {
    if (!("hours" %in% names(y))) y$hours <- y$time_numeric * 24
    y$time_numeric <- as.numeric(y$time_numeric); y$N <- as.numeric(y$N); y$D <- as.numeric(y$D); y$hours <- as.numeric(y$hours)
  } else {
    if (!("hours" %in% names(y))) {
      if ("time" %in% names(y)) {
        m <- regexec("^(\\d+)d(\\d+)h(\\d+)m$", as.character(y$time))
        g <- regmatches(as.character(y$time), m)
        y$hours <- vapply(g, function(z) if (length(z)==4) as.numeric(z[2])*24 + as.numeric(z[3]) + as.numeric(z[4])/60 else NA_real_, 0)
      } else stop("Need either `hours` or `time` to compute time_numeric.")
    }
    if (!("alive" %in% names(y))) stop("Uncorrected input must have `alive`.")
    if (!("dead"  %in% names(y))) stop("Uncorrected input must have `dead`.")
    y$time_numeric <- as.numeric(y$hours) / 24
    y$N <- as.numeric(y$alive)+as.numeric(y$dead)
    y$D <- as.numeric(y$dead)
  }
  
  for (nm in c("cellLine","experiment","plateID","ploidy","glucose")) if (nm %in% names(y)) y[[nm]] <- as.character(y[[nm]])
  for (nm in c("time_numeric","N","D","hours"))                 if (nm %in% names(y)) y[[nm]] <- as.numeric(y[[nm]])
  
  want <- c("time_numeric","N","D","cellLine","experiment","plateID","ploidy","glucose","hours")
  for (nm in setdiff(want, names(y))) y[[nm]] <- if (nm %in% c("time_numeric","N","D","hours")) NA_real_ else NA_character_
  
  if (is_dt) y <- y[, ..want] else y <- y[, want, drop=FALSE]
  y
}

cat("Processing Count Data...\n")
counts_raw <- as_corrected_counts(readRDS("data/counts/uncorrected.Rds"))

counts <- counts_raw %>%
  mutate(exp_key = paste0(cellLine, ".", experiment)) %>%
  filter(exp_key %in% names(glu_map)) %>%
  filter(!is.null(glu_map[exp_key]))

# CHANGED: Calculate continuous ploidy metric (Difference from baseline)
counts$ploidy_val <- NA_real_
counts$ploidy_metric <- NA_real_

unique_lines <- unique(counts$cellLine)
for (cline in unique_lines) {
  map <- ploidy_map[[cline]]
  idx <- counts$cellLine == cline
  
  # Assign absolute numeric values
  counts$ploidy_val[idx] <- map[counts$ploidy[idx]]
  
  # Calculate difference from baseline (minimum ploidy for that line)
  baseline <- min(map)
  counts$ploidy_metric[idx] <- counts$ploidy_val[idx] - baseline
}

counts$N_live <- counts$N - counts$D

counts <- counts %>%
  rename(G0 = glucose) %>%
  select(cellLine, experiment, ploidy_metric, G0, hours, N_live, D) %>%
  arrange(cellLine, experiment, G0, hours)

# ==============================================================================
# 3. LOAD GLUCOSE DATA (Robust "LO" and NA handling)
# ==============================================================================
cat("Processing Glucose Data...\n")

unique_folders <- unique(glu_map[!sapply(glu_map, is.null)])
calib_store <- list()
gluc_store  <- list()

for (folder in unique_folders) {
  path_base <- file.path("data/glucose/processed", folder)
  
  # 1. Load Calibration
  calib_df <- fread(file.path(path_base, "calibration.csv"))
  calib_df$calib_batch <- folder
  calib_store[[folder]] <- calib_df[, c("G", "Lum", "calib_batch")]
  
  # 2. Load Data
  xg <- fread(file.path(path_base, "data.csv"))
  
  # Identify valid replicate columns (R1:4 or R1:3 depending on file)
  measure_cols <- intersect(paste0("R", 1:4), names(xg))
  
  LOD_value <- 20.0 
  xg_flags  <- xg # Copy for censorship flags
  
  # Clean "LO" and NAs
  for(col in measure_cols) {
    # Check for specific strings/NAs
    is_lo <- xg[[col]] == "LO"
    is_missing <- xg[[col]] == "" | is.na(xg[[col]]) | xg[[col]] == "NA"
    
    if(is.character(xg[[col]])) {
      xg[[col]][is_lo] <- as.character(LOD_value)
      xg[[col]] <- as.numeric(xg[[col]])
      
      # Flag: 1=Censored(LO), 0=Observed, NA=Missing
      xg_flags[[col]] <- ifelse(is_lo, 1, 0)
      xg_flags[[col]][is_missing] <- NA
    } else {
      xg_flags[[col]] <- 0
      xg_flags[[col]][is_missing] <- NA
    }
  }
  
  # Map ploidy strings to numeric metric (Continuous)
  xg$ploidy_metric <- NA_real_
  for(ln in unique(xg$CellLine)) {
    if(ln %in% names(ploidy_map)) {
      map <- ploidy_map[[ln]]
      base_val <- min(map)
      vals <- map[xg$ploidy[xg$CellLine == ln]]
      xg$ploidy_metric[xg$CellLine == ln] <- vals - base_val
    }
  }
  
  # Melt Values and Flags
  xg_long <- reshape2::melt(xg, measure.vars = measure_cols, variable.name = "Replicate", value.name = "lum")
  xg_flags_long <- reshape2::melt(xg_flags, measure.vars = measure_cols, variable.name = "Replicate", value.name = "is_censored")
  
  xg_long$is_censored <- xg_flags_long$is_censored
  
  # Drop rows where data was missing (NA)
  xg_long <- xg_long %>% filter(!is.na(lum))
  
  xg_long$hours <- xg_long$Day * 24
  xg_long$calib_batch <- folder
  
  # Fix dilution: 1 means undiluted (1.0), otherwise use 1000/X legacy logic
  xg_long$dilution <- ifelse(xg_long$`Dilution Factor` == 1, 
                             1.0, 
                             1000 / xg_long$`Dilution Factor`)
  
  gluc_store[[folder]] <- xg_long %>%
    select(CellLine, ploidy_metric, G0, hours, lum, dilution, calib_batch, is_censored)
}

all_glucose_raw <- do.call(rbind, gluc_store)
if(any(is.na(all_glucose_raw$lum))) stop("Error: NAs remaining in glucose luminescence data.")

all_calib_raw   <- do.call(rbind, calib_store)
# ==============================================================================
# 4. CONSTRUCT STAN STRUCTURE
# ==============================================================================
# CHANGED: Group by ploidy_metric instead of is_high_ploidy
conditions <- counts %>%
  distinct(cellLine, experiment, ploidy_metric, G0) %>%
  mutate(well_id = row_number())

N_wells <- nrow(conditions)
G0_vector <- numeric(N_wells)
line_ids  <- integer(N_wells)
exp_ids   <- integer(N_wells)
ploidy_vec <- numeric(N_wells) # CHANGED: numeric vector
has_starvation <- integer(N_wells) # CHANGED: New flag for protocol

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
  ploidy_vec[i] <- cond$ploidy_metric
  exp_ids[i]   <- calib_map_int[[glu_map[[paste0(cond$cellLine, ".", cond$experiment)]]]]
  G0_vector[i] <- as.numeric(cond$G0)
  
  # CHANGED: Set starvation flag. 
  # Logic: MCF10A = 0 (No starve), All others = 1 (Starve -6 to 0)
  if (grepl("MCF10A", cond$cellLine)) {
    has_starvation[i] <- 0
  } else {
    has_starvation[i] <- 1
  }
  
  sub_c <- counts %>%
    filter(cellLine == cond$cellLine, experiment == cond$experiment,
           ploidy_metric == cond$ploidy_metric, G0 == cond$G0)
  sub_c$well_idx <- i
  sub_c$grid_idx <- get_grid_idx(sub_c$hours)
  obs_counts_list[[i]] <- sub_c %>% select(well_idx, grid_idx, N_live, D)
  
  # For glucose data matching, we might have NAs in ploidy if not mapped perfectly,
  # but main mapping is done.
  sub_g <- all_glucose_raw %>%
    filter(calib_batch == glu_map[[paste0(cond$cellLine, ".", cond$experiment)]],
           CellLine == cond$cellLine, G0 == cond$G0)
  
  # Try to match ploidy if available, else take all (often glucose data is pooled or scarce)
  sub_g_p <- sub_g %>% filter(abs(ploidy_metric - cond$ploidy_metric) < 1e-3 | is.na(ploidy_metric))
  
  if (nrow(sub_g_p) > 0) {
    sub_g_p$well_idx <- i
    sub_g_p$grid_idx <- get_grid_idx(sub_g_p$hours)
    obs_gluc_list[[i]] <- sub_g_p %>% select(well_idx, grid_idx, lum, dilution, is_censored)
  }
}

long_counts <- do.call(rbind, obs_counts_list)
long_gluc   <- do.call(rbind, obs_gluc_list)

all_calib_raw$exp_idx <- calib_map_int[all_calib_raw$calib_batch]

# ==============================================================================
# 5. CALCULATE DATA-SPECIFIC PRIORS
# ==============================================================================
idx_t0 <- which(long_counts$grid_idx == 1)
prior_mu_N0 <- mean(log(long_counts$N_live[idx_t0] + 1e-9))
prior_mu_D0 <- mean(log(long_counts$D[idx_t0] + 1e-9))

# Estimate Fixed Calibration per Experiment
cat("Estimating fixed calibration (log-fit) per experiment...\n")
N_exps <- length(u_calibs)
a_fix <- rep(NA_real_, N_exps)
b_fix <- rep(NA_real_, N_exps)
sigma_fix <- rep(NA_real_, N_exps)

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
    a_fix[e] <- 10.0
    b_fix[e] <- 10.0
    sigma_fix[e] <- 0.5
    next
  }
  
  fit0 <- try(lm(Lum ~ G, data = sub), silent = TRUE)
  a0 <- if (!inherits(fit0, "try-error")) max(1e-6, as.numeric(coef(fit0)["G"])) else 1.0
  b0 <- max(1e-6, min(sub$Lum, na.rm = TRUE))
  
  opt <- optim(c(a0, b0), costf, x = sub, method = "Nelder-Mead")
  a_est <- max(1e-6, abs(opt$par[1]))
  b_est <- max(1e-6, abs(opt$par[2]))
  
  a_fix[e] <- a_est
  b_fix[e] <- b_est
  
  mu_est <- sub$G * a_est + b_est
  valid_idx <- mu_est > 0 & sub$Lum > 0
  if (sum(valid_idx) > 2) {
    resid <- log(sub$Lum[valid_idx]) - log(mu_est[valid_idx])
    sigma_fix[e] <- max(0.1, sd(resid)) 
  } else {
    sigma_fix[e] <- 0.1 
  }
}

# ==============================================================================
# 6. SAVE
# ==============================================================================
stan_data <- list(
  
  line_map = line_map_int,
  
  N_wells     = N_wells,
  N_lines     = length(u_lines),
  N_exps      = N_exps,
  N_obs_count = nrow(long_counts),
  N_obs_gluc  = nrow(long_gluc),
  N_obs_calib = nrow(all_calib_raw),
  
  N_grid      = length(t_grid),
  t_grid      = t_grid,
  
  line_id     = line_ids,
  ploidy_metric = ploidy_vec, # CHANGED: Continuous variable
  has_starvation = has_starvation, # CHANGED: Protocol flag
  
  exp_id      = exp_ids,
  
  well_idx_count = long_counts$well_idx,
  grid_idx_count = long_counts$grid_idx,
  N_obs          = long_counts$N_live,
  D_obs          = long_counts$D,
  
  well_idx_gluc  = long_gluc$well_idx,
  grid_idx_gluc  = long_gluc$grid_idx,
  lum_obs        = long_gluc$lum,
  dilution       = long_gluc$dilution,
  is_censored    = as.integer(long_gluc$is_censored),
  
  calib_exp_idx  = all_calib_raw$exp_idx,
  calib_G        = all_calib_raw$G,
  calib_Lum      = all_calib_raw$Lum,
  
  G0_per_well    = G0_vector,
  grainsize      = 1,
  
  prior_mu_N0_mean = prior_mu_N0,
  prior_mu_N0_sd   = 0.5,
  prior_mu_D0_mean = prior_mu_D0,
  prior_mu_D0_sd   = 1.0,
  
  calib_a_fixed = a_fix,
  calib_b_fixed = b_fix,
  calib_sigma_fixed = sigma_fix
)

saveRDS(stan_data, "data/stan_ready_data.Rds")
cat("Data saved to 'data/stan_ready_data.Rds'.\n")