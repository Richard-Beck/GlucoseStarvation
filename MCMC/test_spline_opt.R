library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(patchwork); library(deSolve)

# ---- 1. Fixed Data Extraction (Adds 'Condition') ----
extract_stan_data <- function(d, cond_label="Condition") {
  # Metadata
  meta <- tibble(
    well_idx = seq_len(d$N_wells), 
    G0 = as.numeric(d$G0_per_well), 
    exp_id = as.integer(d$exp_id)
  ) %>%
    mutate(
      Condition = cond_label, # <--- Added this to fix the error
      G0_lbl = factor(paste0(G0, " mM"), levels = paste0(sort(unique(G0)), " mM"))
    )
  
  # Counts
  df_N <- tibble(
    well_idx = as.integer(d$well_idx_count), 
    time = as.numeric(d$t_grid[d$grid_idx_count]),
    NL = as.numeric(d$N_obs), 
    ND = as.numeric(d$D_obs)
  ) %>% pivot_longer(c(NL, ND), names_to="type", values_to="value")
  
  # Glucose (Back-calc)
  df_G <- tibble(
    well_idx = as.integer(d$well_idx_gluc), 
    time = as.numeric(d$t_grid[d$grid_idx_gluc]),
    lum = as.numeric(d$lum_obs), 
    dil = as.numeric(d$dilution)
  ) %>%
    left_join(meta, by="well_idx") %>%
    mutate(
      a = as.numeric(d$calib_a_fixed[exp_id]), 
      b = as.numeric(d$calib_b_fixed[exp_id]),
      # Clamp to small epsilon to avoid log(0) issues later
      value = pmax(1e-6, (lum - b) / (a * dil)), 
      type = "G"
    ) %>% select(well_idx, time, type, value)
  
  bind_rows(df_N, df_G) %>% left_join(meta, by="well_idx")
}

# ---- 2. Fixed Splines (Safe Select & Grouping) ----
get_splines <- function(obs_df, n_dense=200, spar=NULL) {
  
  fit_one <- function(curr_df, keys) {
    if(nrow(curr_df) < 4) return(tibble())
    
    # Aggregate duplicates
    d0 <- curr_df %>% group_by(time) %>% summarise(value=mean(value), .groups="drop")
    t_rng <- seq(min(d0$time), max(d0$time), length.out=n_dense)
    
    # Check 'type' from the grouping keys (safe access)
    is_gluc <- keys$type == "G"
    
    # Fit on log-scale for Glucose to prevent negative values/oscillations
    y_fit <- if(is_gluc) log(d0$value) else d0$value
    
    sp <- smooth.spline(d0$time, y_fit, spar=spar)
    
    pred0 <- predict(sp, t_rng, deriv=0)$y
    pred1 <- predict(sp, t_rng, deriv=1)$y
    
    # Transform back: if y = log(G), then G = exp(y) and dG/dt = G * dy/dt
    if(is_gluc) {
      val <- exp(pred0)
      grad <- val * pred1 
    } else {
      val <- pred0
      grad <- pred1
    }
    
    tibble(time=t_rng, value_hat=val, d_dt=grad)
  }
  
  obs_df %>% 
    group_by(well_idx, type) %>% 
    group_modify(fit_one) %>% 
    ungroup() %>%
    # Safe join: uses any_of so it won't crash if a column is missing
    left_join(
      obs_df %>% select(well_idx, any_of(c("G0_lbl", "Condition", "exp_id"))) %>% distinct(), 
      by="well_idx"
    )
}
# ---- ODE RHS (Robust) ----
rhs_optim <- function(NL, ND, G, p) {
  # Unpack parameters (p is vector)
  theta <- p[1]; kp <- p[2]; kd <- p[3]; kd2 <- 0; g50a <- p[5]; na <- p[6]
  gamma <- p[7]; nd <- p[8]; yield <- p[9]; v2 <- p[10]; g50g_mult <- p[11]; ng <- p[12]
  
  # Safety Clamps
  G    <- pmax(G, 1e-12) 
  g50d <- pmax(g50a * gamma, 1e-12)
  g50g <- pmax(g50d * g50g_mult, 1e-12) # Prevents log(0) if mult is 0
  
  # Hill Functions (on log scale for stability)
  logG <- log(G)
  actA <- plogis(na * (logG - log(g50a)))
  inhD <- 1 - plogis(nd * (logG - log(g50d)))
  term_g <- plogis(ng * (logG - log(g50g)))
  
  # Equations
  v1 <- 0#kp / yield
  growth <- kp * (1 - NL/theta) * actA - kd * inhD - kd2 * NL/theta
  
  dNL <- growth * NL
  dND <- kd * NL * inhD + kd2 * (NL^2)/theta
  dG  <- -NL * (v1 * actA + v2 * term_g) / 2
  
  c(dNL, dND, dG)
}

# ---- Optimization Wrapper ----
fit_ode_params <- function(spline_df, stan_data, n_starts=20) {
  
  # Pivot splines to wide format (Time x State)
  wide <- spline_df %>%
    select(well_idx, time, type, value_hat, d_dt) %>%
    pivot_wider(names_from = type, values_from = c(value_hat, d_dt)) %>%
    drop_na()
  
  # Objective Function
  obj_fn <- function(p) {
    if(any(p < 0)) return(1e20) # Hard lower bound
    
    # Calc Derivatives based on params
    res <- mapply(function(n, d, g) rhs_optim(n, d, g, p), 
                  wide$value_hat_NL, wide$value_hat_ND, wide$value_hat_G)
    
    # Errors
    e_nl <- wide$d_dt_NL - res[1,]
    e_nd <- wide$d_dt_ND - res[2,]
    e_g  <- wide$d_dt_G  - res[3,]
    
    # Sum of Squares (handle NA)
    sse <- sum(c(e_nl^2, e_nd^2, e_g^2), na.rm = TRUE)
    if(!is.finite(sse)) return(1e20)
    return(sse)
  }
  
  # Multistart setup
  lb <- as.numeric(stan_data$lower_b)
  ub <- as.numeric(stan_data$upper_b)
  
  best_fit <- list(value = Inf)
  
  for(i in 1:n_starts) {
    # Log-uniform random start
    start_p <- exp(runif(length(lb), log(pmax(lb, 1e-10)), log(ub)))
    
    # Try/Catch for optimizer safety
    try({
      opt <- optim(start_p, obj_fn, method="L-BFGS-B", lower=lb, upper=ub, control=list(maxit=1000))
      if(opt$value < best_fit$value) best_fit <- opt
    }, silent=TRUE)
    print(opt$value)
    print(opt$par)
  }
  
  return(best_fit$par)
}

# ---- Simulator ----
run_sim <- function(p_vec, stan_data, obs_df) {
  
  # Determine initial conditions per well
  inits <- obs_df %>% 
    group_by(well_idx) %>%
    summarise(
      NL = mean(value[type=="NL" & time==min(time[type=="NL"])]),
      ND = 0, # Assume 0 dead at start
      G  = unique(G0[1]), # Use theoretical G0
      .groups="drop"
    )
  
  t_vec <- sort(unique(obs_df$time))
  
  # Solve ODE per well
  run_one_well <- function(w) {
    y0 <- inits %>% filter(well_idx == w) %>% select(NL, ND, G) %>% unlist()
    out <- ode(y=y0, times=t_vec, func=function(t, y, p) list(rhs_optim(y[1], y[2], y[3], p)), parms=p_vec)
    as.data.frame(out) %>% mutate(well_idx = w)
  }
  
  lapply(unique(obs_df$well_idx), run_one_well) %>% bind_rows() %>% as_tibble()
}

# ---- Plotter ----
plot_results <- function(obs_df, sim_df, title="") {
  
  # Combine Data for ggplot
  sim_long <- sim_df %>% 
    pivot_longer(c(NL, ND, G), names_to="type", values_to="pred") %>%
    left_join(obs_df %>% select(well_idx, G0_lbl, Condition) %>% distinct(), by="well_idx")
  
  p_cnt <- ggplot() +
    geom_line(data=sim_long %>% filter(type!="G"), aes(time, pred, color=type, group=interaction(well_idx, type))) +
    geom_point(data=obs_df %>% filter(type!="G"), aes(time, value, color=type), alpha=0.5) +
    facet_wrap(~G0_lbl, scales="free") + theme_bw() + labs(title="Cell Counts")
  
  p_gluc <- ggplot() +
    geom_line(data=sim_long %>% filter(type=="G"), aes(time, pred, group=well_idx)) +
    geom_point(data=obs_df %>% filter(type=="G"), aes(time, value), alpha=0.5) +
    scale_y_log10() + facet_wrap(~G0_lbl, scales="free") + theme_bw() + labs(title="Glucose")
  
  (p_cnt | p_gluc) + plot_annotation(title=title)
}

# ---- Run It ----
stan_data <- readRDS("data/pathfinder_1cond/model_B1cond/MDA-MB-231/hi/stan_data.Rds")

# 1. Clean Data
obs <- extract_stan_data(stan_data)

# 2. Fit Splines (Log-G ensures stability)
splines <- get_splines(obs, spar=0.7) # spar=0.7 smoothes noise better

# 3. Fit ODE Parameters
best_params <- fit_ode_params(splines, stan_data, n_starts=3)

# 4. Simulate & Plot
sims <- run_sim(best_params, stan_data, obs)
plot_results(obs, sims, title="Cleaned Pipeline Fit")