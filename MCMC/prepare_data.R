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


# usage:
# counts <- as_corrected_counts(readRDS("data/counts/uncorrected.Rds"))
# counts <- as_corrected_counts(readRDS("data/counts/corrected.Rds"))



cat("Processing Count Data...\n")
counts_raw <- as_corrected_counts(readRDS("data/counts/uncorrected.Rds"))

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
cat("Processing Glucose Data...\n")

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
# 5. CALCULATE DATA-SPECIFIC PRIORS (ICs and Calib)
# ==============================================================================
idx_t0 <- which(long_counts$grid_idx == 1)
prior_mu_N0 <- mean(log(long_counts$N_live[idx_t0] + 1e-9))
prior_mu_D0 <- mean(log(long_counts$D[idx_t0] + 1e-9))

# Estimate Fixed Calibration per Experiment
cat("Estimating fixed calibration (log-fit) per experiment...\n")
N_exps <- length(u_calibs)
a_fix <- rep(NA_real_, N_exps)
b_fix <- rep(NA_real_, N_exps)

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
    next
  }
  
  fit0 <- try(lm(Lum ~ G, data = sub), silent = TRUE)
  a0 <- if (!inherits(fit0, "try-error")) max(1e-6, as.numeric(coef(fit0)["G"])) else 1.0
  b0 <- max(1e-6, min(sub$Lum, na.rm = TRUE))
  
  opt <- optim(c(a0, b0), costf, x = sub, method = "Nelder-Mead")
  a_fix[e] <- max(1e-6, abs(opt$par[1]))
  b_fix[e] <- max(1e-6, abs(opt$par[2]))
}

# ==============================================================================
# 6. SAVE
# ==============================================================================
# Note: ODE priors (means/sds) are now removed. They will be injected by the runner.
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
  
  # Data-Derived IC Priors
  prior_mu_N0_mean = prior_mu_N0,
  prior_mu_N0_sd   = 0.5,
  prior_mu_D0_mean = prior_mu_D0,
  prior_mu_D0_sd   = 1.0,
  
  # Fixed Calibration
  calib_a_fixed = a_fix,
  calib_b_fixed = b_fix
)

saveRDS(stan_data, "data/stan_ready_data.Rds")
cat("Data saved to 'data/stan_ready_data.Rds'.\n")