# model_B.R
# Utilities to debug STAN model_B in pure R:
#  1) Draw one parameter set from the prior
#  2) Simulate the full model (deSolve)
#  3) Compute the full data log-likelihood for a parameter set
#
# Expected data structure: the same list produced by prepare_stan_data.R,
# with fields used below (see ?prepare_stan_data.R).

suppressPackageStartupMessages({
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("Please install 'deSolve' to use simulate_full_model().")
  }
})

# -------------------------------
# Helpers mirroring Stan code
# -------------------------------
.softcap <- function(x, cap) cap - log1p(exp(cap - x))

# R version of model_b_ode working in transformed space (u=log NL, v=log ND, w=logit s)
.model_b_ode_R <- function(t, y, p) {
  # p is numeric vector of length 13:
  # [theta,kp,kd,kd2,g50a,na,g50d,nd,v1,v2,G0,G_floor,s_eps]
  theta <- p[1]; kp <- p[2]; kd <- p[3]; kd2 <- p[4]
  g50a  <- p[5]; na <- p[6]; g50d <- p[7]; nd <- p[8]
  v1    <- p[9]; v2 <- p[10]
  G0    <- p[11]; G_floor <- p[12]; s_eps <- p[13]
  
  NL <- exp(y[1]); ND <- exp(y[2])
  s_raw <- 1/(1 + exp(-y[3]))
  s     <- s_eps + (1 - 2*s_eps) * s_raw
  G     <- G_floor + (G0 - G_floor) * s
  
  eps <- 1e-12
  actA   <- 1 / (1 + (g50a / (G + eps))^na)
  term_d <- 1 / (1 + (g50d / (G + eps))^nd)
  inhD   <- 1 - term_d
  confl  <- kd2 * (NL * NL) / (theta + eps)
  
  dNL <- kp * NL * (1 - NL / (theta + eps)) * actA - kd * NL * inhD - confl
  dND <- kd * NL * inhD + confl
  
  taper <- (s - s_eps) / (1 - 2*s_eps + eps)
  dG    <- -NL * (v1 * actA + v2 * term_d) / 2 * taper
  
  denom <- (G0 - G_floor) * (1 - 2*s_eps) * s_raw * (1 - s_raw) + eps
  
  list(c(dNL / (NL + eps),
         dND / (ND + eps),
         dG  / denom))
}

# Build per-well ODE parameter vector from top-level parameters
.build_well_params <- function(params, stan_data, l, h, G0) {
  # Smooth caps
  cap_log_main <- 40.0
  cap_log_hill <- 6.0
  
  p10 <- numeric(10)
  raw <- params$mu_global + params$sigma_line * params$z_line[, l] + params$beta_high * h
  for (pp in 1:10) {
    if (pp %in% c(6,8)) {
      p10[pp] <- 1 + exp(.softcap(raw[pp], cap_log_hill))  # na, nd >= 1 + ...
    } else {
      p10[pp] <- exp(.softcap(raw[pp], cap_log_main))
    }
  }
  G_floor <- min(1e-12, 0.5 * G0)
  s_eps   <- 1e-6
  c(p10, G0, G_floor, s_eps)
}

# Per-well initial condition (in u,v,w space)
.build_y0 <- function(params, l, h) {
  u0 <- params$mu_IC[1] + params$sigma_IC[1] * params$z_IC[1, l] + params$beta_IC[1] * h
  v0 <- params$mu_IC[2] + params$sigma_IC[2] * params$z_IC[2, l] + params$beta_IC[2] * h
  w0 <- 5.0
  c(u0, v0, w0)
}

# -------------------------------
# (1) Draw from priors
# -------------------------------
draw_from_prior <- function(stan_data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  N_lines <- stan_data$N_lines
  N_exps  <- stan_data$N_exps
  
  list(
    mu_global  = rnorm(10, 0, 2),
    sigma_line = rexp(10, rate = 1),
    beta_high  = rnorm(10, 0, 1),
    z_line     = matrix(rnorm(10 * N_lines), nrow = 10, ncol = N_lines),
    
    mu_IC      = c(rnorm(1, stan_data$prior_mu_N0_mean, stan_data$prior_mu_N0_sd),
                   rnorm(1, stan_data$prior_mu_D0_mean, stan_data$prior_mu_D0_sd)),
    sigma_IC   = rexp(2, rate = 1),
    beta_IC    = rnorm(2, 0, 1),
    z_IC       = matrix(rnorm(2 * N_lines), nrow = 2, ncol = N_lines),
    
    calib_sigma = rexp(N_exps, rate = 1),
    
    phi_N = rexp(1, rate = 0.1),
    phi_D = rexp(1, rate = 0.1)
  )
}

# -------------------------------
# (2) Simulate full model
# -------------------------------
simulate_full_model <- function(params, stan_data, method = "lsoda", rtol = 1e-6, atol = 1e-6, maxsteps = 1e5) {
  times <- as.numeric(stan_data$t_grid)
  if (abs(times[1]) < 1e-14) times[1] <- 1e-8
  
  Nw <- stan_data$N_wells
  Ng <- length(times)
  
  # storage: 3 states per grid point per well
  y_hat <- array(NA_real_, dim = c(Nw, Ng, 3))
  
  for (w in seq_len(Nw)) {
    l <- stan_data$line_id[w]
    h <- stan_data$is_high_ploidy[w]
    G0 <- stan_data$G0_per_well[w]
    
    p_w  <- .build_well_params(params, stan_data, l, h, G0)
    y0_w <- .build_y0(params, l, h)
    
    # deSolve expects a function(y, t, parms) but also supports events
    ode_fun <- function(t, y, parms) .model_b_ode_R(t, y, parms)
    
    sol <- try(
      deSolve::ode(y = y0_w, times = times, func = ode_fun, parms = p_w,
                   method = method, rtol = rtol, atol = atol, maxsteps = maxsteps),
      silent = TRUE
    )
    if (inherits(sol, "try-error") || !is.matrix(sol)) {
      stop(sprintf("deSolve failed for well %d (line=%d, high=%d, G0=%.4g): %s",
                   w, l, h, G0, as.character(sol)))
    }
    
    # columns: time, u, v, w
    u <- sol[, "1"]; v <- sol[, "2"]; wlogit <- sol[, "3"]
    # Transform back to (NL, ND, G)
    NL <- exp(u); ND <- exp(v)
    s_raw <- 1/(1 + exp(-wlogit))
    s_eps <- 1e-6
    s     <- s_eps + (1 - 2*s_eps) * s_raw
    G     <- (min(1e-12, 0.5 * G0)) + (G0 - min(1e-12, 0.5 * G0)) * s
    
    y_hat[w, , 1] <- NL
    y_hat[w, , 2] <- ND
    y_hat[w, , 3] <- G
  }
  
  y_hat
}

# -------------------------------
# (3) Log-likelihood
# -------------------------------
log_lik_model <- function(params, stan_data, y_hat = NULL) {
  if (is.null(y_hat)) {
    y_hat <- simulate_full_model(params, stan_data)
  }
  
  # Ensure fixed calibration is present
  if (is.null(stan_data$calib_a_fixed) || is.null(stan_data$calib_b_fixed)) {
    stop("stan_data must contain calib_a_fixed and calib_b_fixed (see fitSTAN.R / prepare_stan_data.R).")
  }
  
  # Count likelihoods (NegBin2 with mean=NL_hat/ND_hat and size=phi)
  ll <- 0
  
  # Pre-extract for speed
  t_grid <- stan_data$t_grid
  for (n in seq_len(stan_data$N_obs_count)) {
    w <- stan_data$well_idx_count[n]
    g <- stan_data$grid_idx_count[n]
    NL_hat <- y_hat[w, g, 1] + 1e-12
    ND_hat <- y_hat[w, g, 2] + 1e-12
    ll <- ll + dnbinom(stan_data$N_obs[n], size = params$phi_N, mu = NL_hat, log = TRUE)
    ll <- ll + dnbinom(stan_data$D_obs[n], size = params$phi_D, mu = ND_hat, log = TRUE)
  }
  
  # Glucose lum (lognormal on mu = a_e*G_hat*dilution + b_e)
  for (n in seq_len(stan_data$N_obs_gluc)) {
    w <- stan_data$well_idx_gluc[n]
    g <- stan_data$grid_idx_gluc[n]
    G_hat <- y_hat[w, g, 3]
    e <- stan_data$exp_id[w]
    mu <- stan_data$calib_a_fixed[e] * G_hat * stan_data$dilution[n] + stan_data$calib_b_fixed[e] + 1e-12
    ll <- ll + dlnorm(stan_data$lum_obs[n], meanlog = log(mu), sdlog = params$calib_sigma[e], log = TRUE)
  }
  
  as.numeric(ll)
}

# -------------------------------
# Convenience: convert a single CmdStan draw (data.frame) to params list
# -------------------------------
stan_draw_to_params <- function(draw_df, stan_data) {
  # draw_df should have scalar columns like "mu_global[1]",..., consistent with CmdStan CSV
  N_lines <- stan_data$N_lines
  N_exps  <- stan_data$N_exps
  
  get_vec <- function(prefix, n) as.numeric(draw_df[1, paste0(prefix, "[", seq_len(n), "]"), drop = TRUE])
  get_mat <- function(prefix, r, c) {
    m <- matrix(NA_real_, nrow = r, ncol = c)
    for (i in 1:r) for (j in 1:c) {
      m[i, j] <- as.numeric(draw_df[1, paste0(prefix, "[", i, ",", j, "]")])
    }
    m
  }
  
  list(
    mu_global   = get_vec("mu_global", 10),
    sigma_line  = get_vec("sigma_line", 10),
    beta_high   = get_vec("beta_high", 10),
    z_line      = get_mat("z_line", 10, N_lines),
    mu_IC       = get_vec("mu_IC", 2),
    sigma_IC    = get_vec("sigma_IC", 2),
    beta_IC     = get_vec("beta_IC", 2),
    z_IC        = get_mat("z_IC", 2, N_lines),
    calib_sigma = get_vec("calib_sigma", N_exps),
    phi_N       = as.numeric(draw_df[1, "phi_N"]),
    phi_D       = as.numeric(draw_df[1, "phi_D"])
  )
}

# -------------------------------
# End of file
# -------------------------------