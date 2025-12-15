# stan_prior_debug.R
suppressPackageStartupMessages({
  library(cmdstanr)
  library(data.table)
})
source("STAN/model_B.R")

# ---- Load data and ensure fixed calibration present ----
stan_data <- readRDS("data/stan_ready_data.Rds")

# Backfill calib_a_fixed/b_fixed if missing (simple log-fit per experiment)
ensure_fixed_calibration <- function(stan_data){
  if (!is.null(stan_data$calib_a_fixed) && !is.null(stan_data$calib_b_fixed)) return(stan_data)
  message("calib_a_fixed / calib_b_fixed missing; estimating quickly from calibration points...")
  e <- stan_data$calib_exp_idx; G <- stan_data$calib_G; Lum <- stan_data$calib_Lum
  N_exps <- stan_data$N_exps
  a_fix <- b_fix <- sdlog_fix <- rep(NA_real_, N_exps)
  for (i in seq_len(N_exps)) {
    idx <- which(e == i & is.finite(G) & is.finite(Lum) & G >= 0 & Lum > 0)
    if (length(idx) < 2) {
      a_fix[i] <- max(1e-6, stan_data$prior_calib_a_mean %||% 1.0)
      b_fix[i] <- max(1e-6, stan_data$prior_calib_b_mean %||% 1.0)
      sdlog_fix[i] <- 0.1
      next
    }
    sub <- data.frame(G=G[idx], Lum=Lum[idx])
    fit0 <- try(lm(Lum ~ G, data = sub), silent = TRUE)
    a0 <- if (!inherits(fit0, "try-error")) max(1e-6, as.numeric(coef(fit0)["G"])) else 1.0
    b0 <- max(1e-6, min(sub$Lum, na.rm = TRUE))
    costf <- function(pars){
      a <- abs(pars[1]); b <- abs(pars[2])
      mu <- sub$G * a + b
      if (any(mu <= 0)) return(1e12)
      sum((log(sub$Lum) - log(mu))^2)
    }
    opt <- optim(c(a0, b0), costf, method = "Nelder-Mead", control = list(maxit = 5000))
    a_hat <- max(1e-6, abs(opt$par[1])); b_hat <- max(1e-6, abs(opt$par[2]))
    mu_hat <- sub$G * a_hat + b_hat
    sdlog_hat <- sd(log(sub$Lum) - log(mu_hat)); if (!is.finite(sdlog_hat) || sdlog_hat <= 0) sdlog_hat <- 0.1
    a_fix[i] <- a_hat; b_fix[i] <- b_hat
  }
  stan_data$calib_a_fixed <- a_fix
  stan_data$calib_b_fixed <- b_fix
  stan_data
}
`%||%` <- function(a,b) if(!is.null(a)) a else b
stan_data <- ensure_fixed_calibration(stan_data)

# ---- Compile Stan model ----
model_path <- "STAN/model_B.stan"
mod <- cmdstan_model(model_path, quiet = TRUE)

stan_data$prior_only <- 1L
logfile <- file("logs/stan_console.log", open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")
fit <- mod$sample(
  data = stan_data,
  chains = 1, seed = 42,
  iter_warmup = 50, iter_sampling = 100,
  init = 0, refresh = 0
)
sink(type = "message")
sink(type = "output")
close(logfile)
## ---- Extract offending p_w from log and replay in R ----

## Robust bracketed-vector parser for "Testing params: [ ... ]"
strip_ansi <- function(x) gsub("\\x1b\\[[0-9;]*[A-Za-z]", "", x)
parse_pw <- function(line) {
  if (length(line) != 1L || is.na(line)) return(NULL)
  s <- strip_ansi(line)
  # normalize fullwidth brackets → ASCII
  s <- chartr("\uFF3B\uFF3D", "[]", s)
  # find first '[' and last ']'
  i1 <- regexpr("\\[", s, perl = TRUE)[1]
  if (i1 < 0) return(NULL)
  i2 <- regexpr("\\][^\\]]*$", s, perl = TRUE)[1]  # last closing bracket
  if (i2 < 0) return(NULL)
  inner <- substr(s, i1 + 1L, i2 - 1L)
  inner <- gsub("\\s+", "", inner)
  parts <- strsplit(inner, ",", fixed = TRUE)[[1]]
  nums <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(nums))) return(NULL)
  nums
}
log_path <- "logs/stan_console.log"
stopifnot(file.exists(log_path))
logs <- readLines(log_path, warn = FALSE)

p_idx <- which(grepl("^Chain\\s*\\d+\\s+Testing params:\\s*\\[", logs, perl = TRUE))
e_idx <- which(grepl("CVode\\(cvodes_mem.*error flag -4", logs))
trial_p <- e_idx - 1L

# robust parser you already fixed (use yours)
offenders <- lapply(logs[p_idx[p_idx %in% trial_p]], parse_pw)
offenders <- Filter(function(x) is.numeric(x) && length(x) == 13, offenders)
cat("Captured", length(offenders), "offending vectors\n")

stopifnot(requireNamespace("deSolve", quietly = TRUE))
times <- as.numeric(stan_data$t_grid); if (abs(times[1]) < 1e-14) times[1] <- 1e-8

# neutral IC in transformed space (log NL, log ND, logit s)
u0 <- if (!is.null(stan_data$prior_mu_N0_mean)) stan_data$prior_mu_N0_mean else 0
v0 <- if (!is.null(stan_data$prior_mu_D0_mean)) stan_data$prior_mu_D0_mean else 0
y0 <- c(u0, v0, 5.0)

ode_fun <- function(t, y, parms) .model_b_ode_R(t, y, parms)

# --- CVODE-like settings (BDF, Newton) ---
rtol_cv   <- 1e-8     # match Stan call
atol_cv   <- 1e-8
maxsteps  <- 100000L
hini_cv   <- 1e-10    # tiny first step like your t0 hack
hmax_cv   <- diff(range(times))/1000 
simulate_with_pw <- function(pw) try(
  deSolve::ode(
    y = y0, times = times, func = ode_fun, parms = pw,
    method   = "bdf",        # VODE BDF (stiff), Newton solver; closest to CVODE BDF
    rtol     = rtol_cv,
    atol     = atol_cv,
    maxsteps = maxsteps,
    hini     = hini_cv,
    hmax     = hmax_cv
  ),
  silent = TRUE
)

out <- lapply(seq_along(offenders), function(i) simulate_with_pw(offenders[[i]]))
status <- vapply(out, function(oi) if (inherits(oi, "try-error") || !is.matrix(oi)) "bdf_fail_R" else "ok_R", "")
print(table(status))

param_names <- c("theta","kp","kd","kd2","g50a","na","g50d","nd","v1","v2","G0","G_floor","s_eps")

offenders <- lapply(offenders,function(oi) {
  names(oi) <- param_names
  oi
})

# Optional: tighten/alter if mismatches remain
# rtol_cv <- 1e-7; atol_cv <- 1e-7                 # stricter tolerances
# hmax_cv <- diff(range(times))/1000                # cap step size
# method = "vode" with nlsolver="Newton", maxord=5  # another close variant:
# deSolve::ode(..., method="vode", nlsolver="Newton", maxord=5, rtol=..., atol=..., ...)

# Optional: write them out for later
# df <- as.data.frame(do.call(rbind, offenders)); names(df) <- paste0("p",1:13)
# write.csv(df, "results/offending_pw.csv", row.names = FALSE)
# saveRDS(replays, "results/offending_pw_replays.Rds")



