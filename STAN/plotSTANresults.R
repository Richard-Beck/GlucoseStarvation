# ============================================================
# Pair-plot + chain-overlap diagnostics for model_B fit (FIXED CALIB VERSION)
# - Works whether calib_a/calib_b are parameters (old) or fixed as data (new)
# ============================================================

library(cmdstanr)
library(bayesplot)
library(posterior)
library(ggplot2)

color_scheme_set("mix-blue-pink")

fit <- readRDS("results/fit_model_B.Rds")

# ------------------------------------------------------------
# 0) Helpers (robust variable detection)
# ------------------------------------------------------------
safe_draws_matrix <- function(vars){
  tryCatch(
    fit$draws(variables = vars, format = "matrix"),
    error = function(e) NULL
  )
}

# Prefer old calibration axis if present, otherwise use calib_sigma[1] if present,
# otherwise fall back to something always-present to keep the pipeline running.
calib_axis <- NA_character_

# Try old-style calib_a[1]
tmp <- safe_draws_matrix("calib_a")
if (!is.null(tmp) && any(grepl("^calib_a\\[1\\]$", colnames(tmp)))) {
  calib_axis <- "calib_a[1]"
}

# Try calib_sigma[1] (new-style fixed calib still estimates calib_sigma)
if (is.na(calib_axis)) {
  tmp <- safe_draws_matrix("calib_sigma")
  if (!is.null(tmp)) {
    # pick first indexed element available (usually [1])
    cs <- colnames(tmp)
    cs1 <- cs[grepl("^calib_sigma\\[", cs)][1]
    if (!is.na(cs1)) calib_axis <- cs1
  }
}

# Last-resort fallback so the rest of the script runs
if (is.na(calib_axis)) {
  tmp <- safe_draws_matrix(c("phi_N","phi_D","mu_IC"))
  if (!is.null(tmp)) {
    # pick a scalar param if present
    if ("phi_N" %in% colnames(tmp)) calib_axis <- "phi_N"
    else if ("phi_D" %in% colnames(tmp)) calib_axis <- "phi_D"
    else if (any(grepl("^mu_IC\\[1\\]$", colnames(tmp)))) calib_axis <- "mu_IC[1]"
  }
}

if (is.na(calib_axis)) stop("Could not find any usable axis (calib_a[1] / calib_sigma[*] / phi_N / phi_D / mu_IC[1]). Check fit object.")

cat("Using axis:", calib_axis, "\n")

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}


# ------------------------------------------------------------
# 1) Extract only parameters we care about (avoid huge y_sim)
# ------------------------------------------------------------
vars_wanted <- c("mu_global","mu_IC","beta_high")

# include calibration pieces only if they exist
if (!is.null(safe_draws_matrix("calib_a"))) vars_wanted <- c(vars_wanted, "calib_a","calib_b")
if (!is.null(safe_draws_matrix("calib_sigma"))) vars_wanted <- c(vars_wanted, "calib_sigma")

m <- fit$draws(variables = unique(vars_wanted), format = "matrix")
pn <- colnames(m)
cat("draw matrix dim:", paste(dim(m), collapse=" x "), "\n")


# ------------------------------------------------------------
# 2) Identify candidate "uptake-like", "growth-like", "gsens-like" via correlations
# ------------------------------------------------------------
mg <- paste0("mu_global[", 1:10, "]")

# "uptake-like": the mu_global most correlated with calib axis
cors_cal <- sapply(mg, \(v) safe_cor(m[, calib_axis], m[, v]))
cat("\nCorrelation with ", calib_axis, ":\n", sep = "")
print(sort(cors_cal))

idx_upt <- which.max(abs(cors_cal))
uptake_par <- sprintf("mu_global[%d]", idx_upt)
cat("\nChosen uptake_par (max |cor| with ", calib_axis, "): ", uptake_par, "\n", sep = "")

# "growth-like": the mu_global most correlated with N0 proxy
cors_N0 <- sapply(mg, \(v) safe_cor(m[, "mu_IC[1]"], m[, v]))
cat("\nCorrelation with mu_IC[1] (N0 proxy):\n")
print(sort(cors_N0))

idx_grow <- which.max(abs(cors_N0))
grow_par <- sprintf("mu_global[%d]", idx_grow)
cat("\nChosen grow_par (max |cor| with mu_IC[1]):", grow_par, "\n")

# "gsens-like": strong partner with uptake_par excluding itself
cors_upt <- sapply(mg, \(v) safe_cor(m[, uptake_par], m[, v]))
cat("\nCorrelation with uptake_par:\n")
print(sort(cors_upt))

exclude <- unique(c(uptake_par))
cand <- setdiff(mg, exclude)
cors_partner <- sapply(cand, \(v) abs(safe_cor(m[, uptake_par], m[, v])))
gsens_par <- cand[which.max(cors_partner)]
cat("\nChosen gsens_par (max |cor| with uptake_par among others):", gsens_par, "\n")

# ------------------------------------------------------------
# 3) Bayesplot pairs
# ------------------------------------------------------------
cat("\nRendering bayesplot mcmc_pairs...\n")

p1 <- mcmc_pairs(m, pars = c(calib_axis, uptake_par))
p2 <- mcmc_pairs(m, pars = c("mu_IC[1]",   grow_par))
p3 <- mcmc_pairs(m, pars = c(uptake_par,   gsens_par))
print(p1); print(p2); print(p3)

# Optional: explicitly split chains 1-2 vs 3-4 (requires draws_array)
a_vars <- unique(c("mu_global","mu_IC","calib_sigma"))
a <- fit$draws(variables = a_vars, format = "draws_array")

cond_12_34 <- pairs_condition(chains = list(c(1,2), c(3,4)))
print(mcmc_pairs(a, pars = c("mu_IC[1]", grow_par, calib_axis, uptake_par), condition = cond_12_34))

# ------------------------------------------------------------
# 4) Better "do chains overlap?" visuals: ggplot colored by chain
# ------------------------------------------------------------
d_vars <- unique(c("mu_global","mu_IC","calib_sigma"))

d <- as_draws_df(fit$draws(variables = d_vars))
d$.chain <- factor(d$.chain)

make_scatter <- function(x, y, title=NULL) {
  ggplot(d, aes(x = .data[[x]], y = .data[[y]], color = .chain)) +
    geom_point(alpha = 0.15, size = 0.6) +
    theme_bw() +
    labs(x = x, y = y, color = "chain", title = title)
}

make_contour <- function(x, y, title=NULL) {
  ggplot(d, aes(x = .data[[x]], y = .data[[y]], color = .chain)) +
    geom_density_2d(linewidth = 0.8) +
    theme_bw() +
    labs(x = x, y = y, color = "chain", title = title)
}

# Core three diagnostics, colored by chain:
print(make_scatter(calib_axis, uptake_par, paste0("Scatter by chain: ", calib_axis, " vs uptake_par")))
print(make_scatter("mu_IC[1]",  grow_par,   "Scatter by chain: mu_IC[1] vs grow_par"))
print(make_scatter(uptake_par,  gsens_par,  "Scatter by chain: uptake_par vs gsens_par"))

# Often clearer than scatter:
print(make_contour(calib_axis, uptake_par, paste0("Contours by chain: ", calib_axis, " vs uptake_par")))
print(make_contour("mu_IC[1]",  grow_par,   "Contours by chain: mu_IC[1] vs grow_par"))
print(make_contour(uptake_par,  gsens_par,  "Contours by chain: uptake_par vs gsens_par"))

# ------------------------------------------------------------
# 5) Extra quick sanity diagnostics (cheap, useful)
# ------------------------------------------------------------

# A) Rhat / ESS summary on key blocks (will be garbage if solver fails, but still informative)
cat("\nSummary (mu_global, mu_IC, calib_sigma):\n")
print(fit$summary(variables = c("mu_global","mu_IC","calib_sigma")))

# B) Divergence and treedepth rates (if present)
diag <- fit$sampler_diagnostics(format = "draws_df")
if ("divergent__" %in% names(diag)) {
  div_rate <- diag %>% group_by(.chain) %>% summarise(div_rate = mean(divergent__ > 0), .groups = "drop")
  cat("\nDivergence rate by chain:\n"); print(div_rate)
}
if ("treedepth__" %in% names(diag)) {
  td_summ <- diag %>% group_by(.chain) %>% summarise(p_hit_max = mean(treedepth__ >= max(treedepth__)), .groups = "drop")
  cat("\nTreedepth quick summary by chain (rough):\n"); print(td_summ)
}

