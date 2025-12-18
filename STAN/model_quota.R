library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------
# 1. ODE Model (q is normalized quota fraction qN in [0,1])
# ----------------------------
ode_quota_full_qN <- function(t, y, pars) {
  with(as.list(pars), {
    u <- y[1]; v <- y[2]; G_raw <- y[3]; qN_raw <- y[4]
    NL <- exp(u); ND <- exp(v)
    
    # bounds
    G  <- pmax(G_raw, 0)
    qN <- pmin(pmax(qN_raw, 0), 1)
    
    # uptake gate + drive (dimensionless)
    gate  <- G / (KG + G + 1e-12)
    drive <- pmax(1 - qN, 0)
    
    # influx in "quota-fraction units per time" (q_set == 1 permanently)
    Jin <- k_home * gate * drive
    
    # growth/death use normalized half-sats directly in (0,1)
    mu    <- kp * (qN^na) / (qN^na + q50gN^na + 1e-18)
    delta <- kd * (q50dN^nd) / (q50dN^nd + qN^nd + 1e-18)
    
    b <- mu * (1 - NL/theta)
    d <- delta + kd2 * NL/theta
    
    du <- b - d
    dv <- (delta * NL + kd2 * (NL*NL)/theta) / ND
    
    # glucose mass balance (sigma_G is a true yield/stoichiometry parameter)
    dG  <- -NL * Jin / sigma_G
    
    # quota-fraction balance
    dqN <- Jin - (alpha_m * kp) - qN*b - (alpha_c * b)
    
    list(c(du, dv, dG, dqN))
  })
}

# ----------------------------
# 2. Parameter Constructor
#   - positive params: log -> exp()
#   - (0,1) params: raw in (0,1)
#   - NEW: alpha_m/alpha_c reparam via rho + r, enforcing alpha_m+alpha_c < nu
#   - OPTIONAL: q50gN/q50dN as fractions of abundant-glucose QSS q_star
# ----------------------------
make_pars_full_qN <- function(p, G0, use_qstar_scaling=TRUE) {
  
  get_log <- function(name) {
    val <- p[name]
    if (is.na(val)) stop(paste("Missing parameter:", name))
    exp(val)
  }
  
  get_raw01 <- function(name) {
    val <- p[name]
    if (is.na(val)) stop(paste("Missing parameter:", name))
    if (!is.finite(val) || val <= 0 || val >= 1) stop(paste(name, "must be in (0,1)"))
    as.numeric(val)
  }
  
  theta   <- get_log("log_theta")
  kp      <- get_log("log_kp")
  kd      <- get_log("log_kd")
  kd2     <- get_log("log_kd2")
  
  nu      <- get_log("log_nu")
  sigma_G <- get_log("log_sigma_G")
  KG      <- get_log("log_KG")
  
  # ---- NEW: alpha_m/alpha_c via (rho, r) ----
  # rho in (0,1): total treadmill burden as a fraction of nu
  # r   in (0,1): split of that burden going to alpha_c vs alpha_m
  rho <- get_raw01("rho")   # (alpha_m + alpha_c)/nu
  r   <- get_raw01("r")     # alpha_c/(alpha_m + alpha_c)
  
  s <- nu * rho
  alpha_c <- s * r
  alpha_m <- s * (1 - r)
  
  # ---- Half-sats ----
  # Option A (recommended for ridge-killing): express as fractions of q_star
  # Option B: use raw half-sats directly (set use_qstar_scaling=FALSE)
  if (use_qstar_scaling) {
    q_star <- (nu - (alpha_m + alpha_c)) / (nu + 1)    # always > 0 by construction
    q_star <- pmin(pmax(q_star, 1e-6), 1-1e-6)
    
    q50g_frac <- get_raw01("q50g_frac")
    q50d_frac <- get_raw01("q50d_frac")
    
    q50gN <- q50g_frac * q_star
    q50dN <- q50d_frac * q_star
  } else {
    q50gN <- get_raw01("q50gN")
    q50dN <- get_raw01("q50dN")
    q_star <- NA_real_
  }
  
  na <- get_log("log_na")
  nd <- get_log("log_nd")
  
  k_home <- nu * kp
  
  list(
    theta=theta, kp=kp, kd=kd, kd2=kd2,
    nu=nu, rho=rho, r=r,
    alpha_m=alpha_m, alpha_c=alpha_c,
    sigma_G=sigma_G, KG=KG,
    q50gN=q50gN, q50dN=q50dN, na=na, nd=nd,
    k_home=k_home,
    q_star=q_star,
    G0=G0
  )
}

# ----------------------------
# 3. Simulation Wrapper
# ----------------------------
simulate_full_qN <- function(G0, p, times, use_qstar_scaling=TRUE) {
  pars <- make_pars_full_qN(p, G0, use_qstar_scaling=use_qstar_scaling)
  y0 <- c(u=log(500), v=log(1), G=G0, qN=1.0)
  
  out <- as.data.frame(ode(y=y0, times=times, func=ode_quota_full_qN, parms=pars, method="lsoda"))
  colnames(out) <- c("time", "u", "v", "G", "qN")
  
  out %>%
    mutate(
      NL = exp(u),
      ND = exp(v),
      G  = pmax(G, 0),
      qN = pmin(pmax(qN, 0), 1),
      q  = qN,
      G0 = G0
    )
}

# ----------------------------
# 4. Example parameter set (matches your last Stan-derived values where possible)
#    NOTE: alpha_m/alpha_c are now defined by (rho, r) instead of logs.
# ----------------------------
inv_logit <- function(x) 1/(1+exp(-x))

p_all_qN <- c(
  # main (log-scale)
  log_theta   = 28.954064,
  log_kp      = -2.4937429,
  log_kd      = -4.4905072,
  log_kd2     = -10.481447,
  log_nu      =  4.0433448,
  log_sigma_G = 10.302311,
  log_KG      = -2.7379615,
  
  # NEW alpha parametrization (0,1)
  # choose something plausible; tune as needed
  rho         = 0.30,   # (alpha_m+alpha_c) = 0.10 * nu
  r           = 0.99,   # 10% of total goes to alpha_c, 90% to alpha_m
  
  # half-sat fractions (0,1); you can reuse your old logits here
  q50g_frac   = 0.99,
  q50d_frac   = 0.05,
  
  # hills (log-scale)
  log_na      =  1.6554994,
  log_nd      =  2.5897
)

# ----------------------------
# 5. Run & Plot
# ----------------------------
G0s <- c(0, 0.25, 1, 5, 25)
times <- seq(0, 140, by=0.25)

df_full <- bind_rows(lapply(G0s, function(g) simulate_full_qN(g, p_all_qN, times, use_qstar_scaling=TRUE))) %>%
  mutate(G0 = factor(G0, levels=G0s)) %>%
  pivot_longer(c(NL, ND, G, q), names_to="metric") %>%
  mutate(metric = factor(metric, levels=c("NL","ND","G","q")))

ggplot(df_full, aes(time, value, color=G0)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  scale_color_viridis_d(option="magma", end=0.9) +
  labs(
    title="Quota Model (alpha reparam: rho+r; optional q* scaling for half-sats)",
    subtitle="alpha_m+alpha_c constrained by construction; q50* optionally expressed as fractions of abundant-glucose q*",
    x="Time (h)", y=NULL
  )

