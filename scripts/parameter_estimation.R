## ------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------

library(data.table)
library(rxode2) # Replaces deSolve for high-performance ODE solving
library(future)
library(future.apply)
library(progressr)
library(ggplot2)

## ------------------------------------------------------------
## 1. Model Definitions (using rxode2)
## ------------------------------------------------------------
# This section contains the two model definitions. You can pass either
# 'yang_model' or 'hsieh_model' to the fitting functions.

# -- Model 1: The original model, renamed to 'yang_model'
yang_model <- rxode2({
  # Define model parameters with dummy initial values.
  # These will be overridden by the optimizer during fitting.
  kp    <- 0.1
  Kp    <- 100
  kdmax <- 0.1
  Kd50  <- 50
  nd    <- 4
  kbys  <- 0
  vglc  <- 0.01
  
  # Model equations
  Gp    = 0.5 * (G + sqrt(G^2 + 1e-9^2)) # ~max(G,0)
  mu    = kp * Gp / (Kp + Gp)
  kS    = kdmax * Kd50^nd / (Gp^nd + Kd50^nd)
  kB    = kbys * ND
  gate  = 1/(1 + exp(-G/1e-3)) # ~Heaviside(G)
  
  # Differential equations for the state variables
  d/dt(NL) = (mu - kS - kB) * NL
  d/dt(ND) = (kS + kB) * NL
  d/dt(G)  = -vglc * NL * gate
})


# -- Model 2: The new Hsieh et al. (2022) model 🔬
hsieh_model <- rxode2({
  # Define model parameters
  Vmax <- 0.1  # Max substrate uptake rate
  Km   <- 100  # Michaelis-Menten constant for uptake
  m    <- 0.01 # Maintenance energy rate (mass substrate per mass biomass per time)
  Yxs  <- 0.5  # True growth yield (mass biomass per mass substrate)
  kmax <- 0.1  # Max death rate under starvation
  
  # Intermediate calculations
  Gp = 0.5 * (G + sqrt(G^2 + 1e-9^2)) # Differentiable version of max(G, 0)
  qs = Vmax * Gp / (Km + Gp)         # Specific substrate uptake rate
  
  # Core of Hsieh model: Net growth rate depends on maintenance surplus/deficit.
  # This uses a C-style ternary operator: (condition) ? (value_if_true) : (value_if_false)
  # rxode2 compiles this efficiently.
  mu_net = (qs >= m) * (Yxs * (qs - m)) + (qs < m) * (kmax * (qs - m) / m)
  
  # The death rate for ND accumulation is the positive part of the net death rate
  # k_death = -min(0, mu_net). This is 0 if mu_net > 0 (growth) and > 0 if mu_net < 0 (death).
  k_death = -min(0, mu_net)
  
  # Differential equations
  d/dt(NL) = mu_net * NL       # Live cells change based on the net growth rate
  d/dt(ND) = k_death * NL      # Dead cells accumulate when live cells die
  d/dt(G)  = -qs * NL         # Substrate is consumed for both growth and maintenance
})


model_A <- rxode2({
  ## parameters (defaults for fitting)
  kp    <- 4.227693e-02
  theta <- 4.630158e+03
  na    <- 7.429587e+00
  g50a  <- 6.670658e-01
  kd    <- 5.935262e-01
  nd    <- 6.114026e+00
  g50d  <- 1.237991e-03
  v     <- 1.799280e-05
  m     <- 2.321256e+00
  g50c  <- 2.270614e-02
  g50y  <- 3.341612e+03
  ky    <- 3.731199e-03
  
  ## growth term
  growth = kp * NL * (1 - NL/theta) * G^na / (G^na + g50a^na)
  ## death term (starvation + bystander)
  death  = kd * NL * (ky/kd * Y/(Y + g50y) + 1 - G^nd / (G^nd + g50d^nd))
  ## substrate consumption
  cons   = v  * NL * G^m  / (G^m  + g50c^m)
  
  ## ODEs
  d/dt(NL) = growth - death
  d/dt(ND) = death
  d/dt(G)  = -cons
  d/dt(Y)  = NL
})


# -- Model B: constrained glucose consumption & confluency‑driven death --

model_B <- rxode2({
  ## parameters (defaults for fitting)
  theta <- 1.586453e+04   # carrying capacity
  kp    <- 1.622025e-02   # proliferation rate
  kd    <- 7.569006e-02   # death rate
  kd2   <- 4.858519e-02   # confluency‑driven death
  g50a  <- 4.570815e-02   # half‑max for growth gating
  na    <- 5.912126e+00   # Hill coefficient for growth
  g50d  <- 3.475508e-03   # half‑max for death gating
  nd    <- 4.167656e+00   # Hill coefficient for death
  v1    <- 6.851104e-05   # consumption rate component A
  v2    <- 1.319974e-06   # consumption rate component D
  epsilon <- 1e-9 ## small constant enables 0 glucose simulation
  
  ## gating functions
  actA  = 1 / (1 + (g50a / (G+epsilon))^na)
  inhD  = 1 - 1 / (1 + (g50d / (G+epsilon))^nd)
  confl = kd2 * NL^2 / theta
  
  ## ODE system
  d/dt(NL) = kp * NL * (1 - NL/theta) * actA
  - kd * NL * inhD
  - confl
  
  d/dt(ND) = kd * NL * inhD
  + confl
  
  d/dt(G)  = - NL * (v1 * actA + v2 * (1 / (1 + (g50d / (G+epsilon))^nd))) / 2
})

## intracellular glucose store model
intra_model <- rxode2({
  ## -----------------------------------------------------------------------
  ## Model Parameters
  ## -----------------------------------------------------------------------
  # Glucose Uptake
  Vmax_uptake <- 0.0001    # mM/(#*hr)      ; max glucose uptake rate per cell
  Km_uptake   <- 0.1       # mM              ; half-saturation for uptake
  
  # Growth & Maintenance from Internal Store (R)
  kp          <- 1.6e-2    # 1/hr            ; max growth rate
  Y_xr        <- 500       # #/mM            ; yield of cells from internal resource R
  m_r         <- 1e-5      # mM/(#*hr)       ; maintenance cost per cell
  
  # Death Rates
  kdStarv     <- 8e-2      # 1/hr            ; max death rate from starvation
  kw          <- 5e-8      # 1/(mM*hr)       ; toxicity constant for waste
  
  # Waste Dynamics
  deltaW      <- 0         # 1/hr            ; waste removal/decay rate
  
  # Growth & Death Switches (based on internal store R)
  r_half_g    <- 0.001     # mM/#            ; per-cell store R for half-maximal growth
  nr_g        <- 4         # unitless        ; hill coef for growth switch
  r_half_d    <- 0.00002   # mM/#            ; per-cell store R for half-maximal death
  nr_d        <- 4         # unitless        ; hill coef for death switch
  
  ## -----------------------------------------------------------------------
  ## Intermediate Calculations & Differential Equations
  ## -----------------------------------------------------------------------
  # Per-cell internal store (average)
  r_cell <- max(0,R / (NL + 1e-9))
  
  # Key rates
  glucose_uptake_rate <- Vmax_uptake * (G / (Km_uptake + G)) 
  growth_rate <- kp * (r_cell^nr_g / (r_half_g^nr_g + r_cell^nr_g)) 
  starvation_death_rate <- kdStarv * (r_half_d^nr_d / (r_half_d^nr_d + r_cell^nr_d)) 
  waste_death_rate <- kw * W  # SIMPLIFIED: Linear toxicity
  
  # State variables
  d/dt(G) = -glucose_uptake_rate*NL
  
  d/dt(NL) = (growth_rate - starvation_death_rate - waste_death_rate)*NL
  
  d/dt(ND) = (starvation_death_rate + waste_death_rate)*NL
  
  d/dt(R) = glucose_uptake_rate - NL*((1/Y_xr) * growth_rate + m_r + r_cell * (starvation_death_rate + waste_death_rate))
  
  d/dt(W) = (glucose_uptake_rate - deltaW * W)*NL # SIMPLIFIED: Waste integrates glucose use
  
})

## ------------------------------------------------------------
## 2. Objective function builder (NEW GENERIC VERSION)
## ------------------------------------------------------------
# This function is now extremely fast to create. It returns a generic
# objective function that takes two arguments: the parameter vector
# and the data subset to fit against.
make_obj_fun <- function(model, weight_dead = 1, eps = 1, penalty = 1e12) {
  
  function(par_vec, data_subset) {
    
    # 1. Get model parameters and set up initial conditions
    # Use the correct accessor for default parameters
    model_pars <- model$.mv$ini
    pars_to_update <- intersect(names(par_vec), names(model_pars))
    model_pars[pars_to_update] <- par_vec[pars_to_update]
    
    # Create default initial conditions (all zeros)
    state_names <- model$state
    model_inits <- setNames(rep(0.0, length(state_names)), state_names)
    
    # Update initial conditions if they are provided in par_vec (e.g., NL0)
    ic_param_names <- paste0(names(model_inits), "0")
    ics_to_update <- intersect(names(par_vec), ic_param_names)
    
    if (length(ics_to_update) > 0) {
      inits_from_pars <- par_vec[ics_to_update]
      names(inits_from_pars) <- sub("0$", "", names(inits_from_pars))
      model_inits[names(inits_from_pars)] <- inits_from_pars
    }
    
    # 2. Run simulations for each unique initial glucose concentration
    all_sims <- rbindlist(lapply(unique(data_subset$glucose), function(g_init) {
      times_for_g <- sort(unique(data_subset[glucose == g_init, hours]))
      if (length(times_for_g) == 0) return(NULL)
      
      et <- rxode2::eventTable()$add.sampling(times_for_g)
      inits <- model_inits
      inits["G"] <- as.numeric(g_init)
      
      sim <- try(model$solve(params = model_pars, events = et, inits = inits, cores = 1), silent = TRUE)
      if (inherits(sim, "try-error")) return(NULL)
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      sim_dt[, glucose := g_init]
      return(sim_dt)
    }))
    
    if (is.null(all_sims) || nrow(all_sims) == 0) return(penalty*nrow(data_subset))
    
    # 3. Merge simulations with data and calculate sum of squared errors
    mm <- merge(data_subset, all_sims, by = c("hours", "glucose"))
    if (nrow(mm) != nrow(data_subset)) return(penalty)
    
    r_alive <- log(mm$alive + eps) - log(mm$NL + eps)
    r_dead  <- log(mm$dead + eps) - log(mm$ND + eps)
    SSE     <- r_alive^2 + weight_dead * r_dead^2
    SSE[!is.finite(SSE)] <- penalty
    SSE <- sum(SSE)
    print(SSE)
    attr(SSE, "n_obs") <- nrow(data_subset) * 2L
    return(SSE)
  }
}

## ------------------------------------------------------------
## 3. Multi-start optimizer (UPDATED with progress bar)
## ------------------------------------------------------------
# This function now uses the 'future' and 'progressr' packages to provide
# a real-time progress bar for the multi-start optimization.
multi_start_fit <- function(obj_fun, par_bounds, n_starts = 100,
                            cores = max(1, parallel::detectCores() - 1)) {
  par_names <- par_bounds$par
  lower     <- par_bounds$lower
  upper     <- par_bounds$upper
  
  starts <- replicate(n_starts, runif(length(par_names), lower, upper))
  colnames(starts) <- paste0("start", seq_len(n_starts))
  rownames(starts) <- par_names
  
  run_one <- function(init_par) {
    tryCatch(
      optim(par = init_par, fn  = obj_fun, method = "L-BFGS-B",
            lower  = lower, upper  = upper),
      error = function(e) list(error = e$message)
    )
  }
  
  # This block uses the future framework for parallel processing and
  # wraps the call in with_progress to show a progress bar.
  with_progress({
    p <- progressor(steps = n_starts)
    out_list <- future_lapply(seq_len(n_starts), function(i) {
      res <- run_one(starts[, i])
      p() # This updates the progress bar after each start is complete
      res
    }, future.seed = TRUE) # future.seed = TRUE for reproducible results
  })
  
  ok    <- vapply(out_list, function(x) is.null(x$error), logical(1))
  fails <- out_list[!ok]
  sucs  <- out_list[ok]
  
  if (!length(sucs)) {
    if (length(fails)) cat("All optimizations failed. Example error:\n", fails[[1]]$error, "\n")
    stop("All optimizations failed.")
  }
  
  best_idx <- which.min(vapply(sucs, `[[`, numeric(1), "value"))
  best     <- sucs[[best_idx]]
  
  list(best = best, all = sucs, fails = fails, starts = starts)
}

## ------------------------------------------------------------
## 4. Fit wrapper (UPDATED with log-space search and auto-bounds)
## ------------------------------------------------------------
# This is the main user-facing function. It orchestrates fitting the model
# to each group defined by 'group_cols' (e.g., each cellLine/ploidy combo).
# It returns a list of fit results, one for each group.
fit_model <- function(dt, model, group_cols, par_bounds, n_starts = 200,
                      weight_dead = 1, eps = 1,
                      cores = max(1, parallel::detectCores() - 1)) {
  
  setDT(dt)
  
  # --- Setup for progress reporting and parallel processing ---
  plan(multisession, workers = cores)
  handlers("progress")
  
  # 1. Create the generic objective function ONCE. This is fast.
  obj_fun_generic <- make_obj_fun(model, weight_dead, eps)
  
  # 2. Split the data into the main groups to be fitted
  groups_to_fit <- split(dt, by = group_cols, keep.by = FALSE)
  
  # 3. Fit each group and collect the results
  all_results <- lapply(names(groups_to_fit), function(group_name) {
    
    cat("Fitting group:", group_name, "\n")
    
    # Ensure the data subset is a data.table.
    data_subset <- groups_to_fit[[group_name]]
    if (!is.data.table(data_subset)) setDT(data_subset)
    
    # --- Automatic bounds for initial conditions ---
    group_par_bounds <- copy(par_bounds)
    
    # If NL0 bounds are not provided, calculate them from the data
    if (!"NL0" %in% group_par_bounds$par) {
      r_alive <- data_subset[hours == 0, alive]
      mean_a <- mean(r_alive, na.rm = TRUE)
      sd_a <- if (length(r_alive) > 1) sd(r_alive, na.rm = TRUE) else 0
      sd_a <- ifelse(is.na(sd_a), 0, sd_a)
      nl0_bounds <- data.table(par = "NL0", lower = max(0, mean_a - 3*sd_a), upper = mean_a + 3*sd_a)
      group_par_bounds <- rbindlist(list(group_par_bounds, nl0_bounds))
    }
    
    # If ND0 bounds are not provided, calculate them from the data
    if (!"ND0" %in% group_par_bounds$par) {
      r_dead <- data_subset[hours == 0, dead]
      mean_d <- mean(r_dead, na.rm = TRUE)
      sd_d <- if (length(r_dead) > 1) sd(r_dead, na.rm = TRUE) else 0
      sd_d <- ifelse(is.na(sd_d), 0, sd_d)
      nd0_bounds <- data.table(par = "ND0", lower = max(0, mean_d - 3*sd_d), upper = mean_d + 3*sd_d)
      group_par_bounds <- rbindlist(list(group_par_bounds, nd0_bounds))
    }
    
    # --- Log-space transformation for optimization ---
    log_par_bounds <- copy(group_par_bounds)
    # Use pmax to avoid log(0)
    log_par_bounds[, `:=`(lower = log(pmax(1e-9, lower)), upper = log(pmax(1e-9, upper)))]
    
    # Create a wrapper for optim. This is a "closure" that "remembers" the data_subset.
    obj_fun_wrapped <- function(p_log) {
      p_natural <- exp(p_log)
      p_named <- setNames(p_natural, log_par_bounds$par)
      obj_fun_generic(p_named, data_subset)
    }
    
    # Run the optimizer in log-space
    opt_res <- multi_start_fit(obj_fun_wrapped, log_par_bounds, n_starts, cores)
    
    # Transform best parameters back to natural space for reporting
    best_par_log <- opt_res$best$par
    best_par <- setNames(exp(best_par_log), log_par_bounds$par)
    
    # Calculate fit statistics
    RSS   <- opt_res$best$value
    n_obs <- attr(RSS, "n_obs")
    if(is.null(n_obs)) n_obs <- nrow(data_subset) * 2L # Fallback
    
    k       <- length(best_par)
    sigma2  <- RSS / n_obs
    logLik  <- -0.5 * n_obs * (log(2 * pi * sigma2) + 1)
    AIC     <- 2 * k - 2 * logLik
    BIC     <- log(n_obs) * k - 2 * logLik
    RMSE    <- sqrt(RSS / n_obs)
    
    # Simulate the final best-fit trajectory for this group
    traj_best <- rbindlist(lapply(unique(data_subset$glucose), function(g_init) {
      times <- sort(unique(data_subset[glucose == g_init, hours]))
      et <- rxode2::eventTable()
      et$add.sampling(times)
      inits <- c(NL = best_par[["NL0"]], ND = best_par[["ND0"]], G = g_init)
      model_pars <- best_par[model$params]
      sim <- as.data.table(model$solve(params = model_pars, events = et, inits = inits, cores = 1))
      setnames(sim, "time", "hours")
      sim[, glucose := g_init]
      # UPDATE: Add the grouping columns to the trajectory output
      for(col in group_cols) {
        sim[, (col) := data_subset[[col]][1]]
      }
      return(sim)
    }))
    
    # Return a structured list for this group
    list(
      group = group_name,
      par_best = best_par,
      fit_stats = data.table(RSS=RSS, RMSE=RMSE, logLik=logLik, AIC=AIC, BIC=BIC, n_obs=n_obs, k=k),
      traj_best = traj_best,
      optim_res = opt_res
    )
  })
  
  names(all_results) <- names(groups_to_fit)
  return(all_results)
}

create_fit_data <- function(opt_results, data_subset, model) {
  
  # Loop over the list of optimization results using base R
  results_list <- lapply(names(opt_results), function(name) {
    
    res <- opt_results[[name]]
    
    # Skip any optimizations that failed
    if (inherits(res, "try-error")) {
      warning(paste("Skipping failed optimization for:", name))
      return(NULL)
    }
    
    # 1. Extract optimized parameters and the corresponding original data
    optimized_pars <- res$par
    original_data <- data_subset[[name]]
    if (is.null(original_data)) return(NULL)
    
    # 2. Set up parameters and initial states using the optimized values
    model_pars <- model$.mv$ini
    pars_to_update <- intersect(names(optimized_pars), names(model_pars))
    model_pars[pars_to_update] <- optimized_pars[pars_to_update]
    
    state_names <- model$state
    model_inits <- setNames(rep(0.0, length(state_names)), state_names)
    
    ic_param_names <- paste0(names(model_inits), "0")
    ics_to_update <- intersect(names(optimized_pars), ic_param_names)
    
    if (length(ics_to_update) > 0) {
      inits_from_pars <- optimized_pars[ics_to_update]
      names(inits_from_pars) <- sub("0$", "", names(inits_from_pars))
      model_inits[names(inits_from_pars)] <- inits_from_pars
    }
    
    # 3. Define a smooth time sequence for plotting
    min_hr <- min(original_data$hours)
    max_hr <- max(original_data$hours)
    smooth_times <- seq(min_hr, max_hr, length.out = 200)
    
    # 4. Run a simulation for each unique glucose condition
    sim_data <- rbindlist(lapply(unique(original_data$glucose), function(g_init) {
      
      et <- rxode2::eventTable()$add.sampling(smooth_times)
      inits <- model_inits
      inits["G"] <- as.numeric(g_init)
      
      sim <- model$solve(params = model_pars, events = et, inits = inits, cores = 1)
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      sim_dt[, glucose := g_init]
      
      return(sim_dt)
    }))
    
    # 5. Add identifying columns by parsing the experiment's name
    id_parts <- strsplit(name, " ")[[1]]
    sim_data[, `:=`(cellLine = id_parts[1], ploidy = id_parts[2])]
    
    return(sim_data)
  })
  
  # Combine all individual data.tables into one and return
  all_fit_data <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  
  return(all_fit_data)
}

#' @title Objective Function for DEoptim (Log-Space Version)
#' @description Accepts log-transformed parameters and converts them to natural space.

deoptim_objective_function <- function(global_par_vec_log, model, data_subsets, par_names, weight_dead = 1, eps = 1, penalty = 1e12) {
  
  # --- NEW: Convert parameters from log-space to natural space ---
  global_par_vec <- exp(global_par_vec_log)
  # -------------------------------------------------------------
  
  # Name the parameter vector for use in the model
  names(global_par_vec) <- par_names
  
  total_sse <- 0
  
  for (d_sub in data_subsets) {
    
    # --- Inner Optimization for Initial Conditions (NL0, ND0) ---
    inner_obj_fun <- function(ic_par_vec) {
      times_for_g <- sort(unique(d_sub$hours))
      et <- rxode2::eventTable()$add.sampling(times_for_g)
      inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = d_sub$glucose[1])
      
      sim <- try(model$solve(params = global_par_vec, events = et, inits = inits, cores = 1), silent = TRUE)
      
      if (inherits(sim, "try-error")) return(penalty)
      
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      mm <- merge(d_sub, sim_dt, by = "hours")
      if (nrow(mm) != nrow(d_sub)) return(penalty)
      
      r_alive <- log(mm$alive + eps) - log(mm$NL + eps)
      r_dead  <- log(mm$dead + eps) - log(mm$ND + eps)
      sse     <- sum(r_alive^2 + weight_dead * r_dead^2)
      
      return(ifelse(!is.finite(sse), penalty, sse))
    }
    
    ic_start_guess <- c(d_sub[hours == 0, mean(alive)], d_sub[hours == 0, mean(dead)])
    ic_start_guess[is.na(ic_start_guess)] <- 0
    
    inner_opt <- optim(
      par = ic_start_guess,
      fn = inner_obj_fun,
      method = "Nelder-Mead",
      control = list(maxit = 50, reltol = 1e-4)
    )
    
    total_sse <- total_sse + inner_opt$value
  }
  
  return(total_sse)
}
library(DEoptim)
library(ggplot2)

#' @title Main Fitting Function using DEoptim (Log-Space Version)
#' @description Orchestrates model fitting in log-space using Differential Evolution.

fit_model_deoptim <- function(dt, model, group_cols, par_bounds,
                              n_pop = 10 * length(par_bounds$par),
                              n_gen_chunk = 20,
                              total_gens = 200,
                              cores = parallel::detectCores()) {
  
  setDT(dt)
  plan(multisession, workers = cores)
  
  groups_to_fit <- split(dt, by = group_cols, keep.by = TRUE)
  
  # --- NEW: Log-transform the parameter bounds for the optimizer ---
  # A small minimum value avoids log(0) --> -Inf
  log_lower_bounds <- log(pmax(par_bounds$lower, 1e-12))
  log_upper_bounds <- log(pmax(par_bounds$upper, 1e-12))
  # ---------------------------------------------------------------
  
  all_results <- lapply(names(groups_to_fit), function(group_name) {
    
    cat("\n--- Fitting Group:", group_name, "---\n")
    data_for_group <- groups_to_fit[[group_name]]
    data_subsets <- split(data_for_group, by = "glucose")
    
    de_optim_result <- NULL
    n_chunks <- ceiling(total_gens / n_gen_chunk)
    
    for(i in 1:n_chunks) {
      cat(sprintf("Running chunk %d of %d (Generations %d to %d)...\n",
                  i, n_chunks, (i-1)*n_gen_chunk + 1, i*n_gen_chunk))
                  
      de_optim_result <- DEoptim(
        fn = deoptim_objective_function,
        lower = log_lower_bounds, # Use log-transformed bounds
        upper = log_upper_bounds, # Use log-transformed bounds
        model = model,
        data_subsets = data_subsets,
        par_names = par_bounds$par,
        control = DEoptim.control(
          NP = n_pop,
          itermax = n_gen_chunk,
          initialpop = if(i > 1) de_optim_result$member$pop else NULL,
          steptol = 50,
          strategy = 2,
          parallelType = "auto",
          trace = TRUE
        )
      )

      cat("Generating intermediate plot...\n")
      
      # --- NEW: Convert best parameters back to natural space for plotting ---
      best_pars_log <- de_optim_result$optim$bestmem
      best_pars_natural <- setNames(exp(best_pars_log), par_bounds$par)
      # --------------------------------------------------------------------
      
      plot_title <- sprintf("Fit for %s - Gen %d\nSSE: %.2f", 
                            group_name, i * n_gen_chunk, de_optim_result$optim$bestval)
      
      tryCatch({
        # The plotting helper receives parameters in natural space
        fit_traj <- create_fit_data(best_pars_natural, model, data_subsets)
        
        p <- ggplot() +
          geom_line(data = fit_traj, aes(x = hours, y = NL, color = factor(glucose)), linewidth = 1, alpha = 0.8) +
          geom_point(data = data_for_group, aes(x = hours, y = alive, color = factor(glucose))) +
          geom_line(data = fit_traj, aes(x = hours, y = ND, color = factor(glucose)), linetype = "dashed", linewidth = 1, alpha = 0.8) +
          geom_point(data = data_for_group, aes(x = hours, y = dead, color = factor(glucose)), shape = 2) +
          facet_wrap(~glucose, scales = "free_y") +
          labs(title = plot_title, x = "Hours", y = "Cell Count (Alive/Dead)", color = "Initial Glucose") +
          theme_bw()
        
        ggsave(filename = sprintf("tmp/fit_progress_%s_gen%04d.png", gsub("[ .]", "_", group_name), i * n_gen_chunk),
               plot = p, width = 10, height = 6)
      },error=function(e) cat("plot generation failed!"))
      
      
    }
    
    # --- NEW: Also convert the final returned object to natural space ---
    de_optim_result$optim$bestmem <- setNames(exp(de_optim_result$optim$bestmem), par_bounds$par)
    de_optim_result$member$pop <- exp(de_optim_result$member$pop)
    # --------------------------------------------------------------------
    
    return(de_optim_result)
  })
  
  names(all_results) <- names(groups_to_fit)
  return(all_results)
}

#' @title Helper to create smooth trajectories for plotting (Corrected)
#' @description Adapted to work with the new direct-argument structure.

create_fit_data <- function(best_pars_vec, model, data_subsets) {
  
  all_sims <- rbindlist(lapply(data_subsets, function(d_sub) {
    # Re-run the quick inner optimization to find the specific initial conditions
    inner_obj_fun <- function(ic_par_vec) {
      times <- sort(unique(d_sub$hours))
      et <- rxode2::eventTable()$add.sampling(times)
      inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = d_sub$glucose[1])
      sim <- try(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1), silent = TRUE)
      if (inherits(sim, "try-error")) return(1e12)
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      mm <- merge(d_sub, sim_dt, by = "hours")
      if (nrow(mm) != nrow(d_sub)) return(1e12)
      sse <- sum((log(mm$alive + 1) - log(mm$NL + 1))^2) + sum((log(mm$dead + 1) - log(mm$ND + 1))^2)
      return(sse)
    }
    
    ic_start_guess <- c(d_sub[hours == 0, mean(alive)], d_sub[hours == 0, mean(dead)])
    ic_start_guess[is.na(ic_start_guess)] <- 0
    inner_opt <- optim(par = ic_start_guess, fn = inner_obj_fun, method = "Nelder-Mead", control = list(maxit = 50))
    best_ics <- inner_opt$par
    
    # Run a final, smooth simulation
    smooth_times <- seq(min(d_sub$hours), max(d_sub$hours), length.out = 200)
    et <- rxode2::eventTable()$add.sampling(smooth_times)
    inits <- c(NL = best_ics[1], ND = best_ics[2], G = d_sub$glucose[1])
    sim <- as.data.table(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1))
    setnames(sim, "time", "hours")
    sim[, glucose := d_sub$glucose[1]]
    return(sim)
  }))
  return(all_sims)
}
library(DEoptim)
library(ggplot2)

#' @title Main Fitting Function (Final Optimized Version)
#' @description Implements data slimming, pre-computation, and a robust
#' closure-based parallelization strategy for maximum efficiency.

fit_model_deoptim_lean <- function(dt, model, group_cols, par_bounds,
                              n_pop = 10 * length(par_bounds$par),
                              n_gen_chunk = 20,
                              total_gens = 200,
                              cores = parallel::detectCores()) {
  
  # --- Create a subfolder for intermediate plots ---
  dir.create("tmp", showWarnings = FALSE)
  
  # 1. SLIM DATA
  required_cols <- c(group_cols, "glucose", "hours", "alive", "dead")
  dt_slim <- dt[, ..required_cols]
  
  setDT(dt_slim)
  plan(multisession, workers = cores)
  
  groups_to_fit <- split(dt_slim, by = group_cols, keep.by = TRUE)
  
  log_lower_bounds <- log(pmax(par_bounds$lower, 1e-12))
  log_upper_bounds <- log(pmax(par_bounds$upper, 1e-12))
  
  all_results <- lapply(names(groups_to_fit), function(group_name) {
    
    cat("\n--- Fitting Group:", group_name, "---\n")
    data_for_group <- groups_to_fit[[group_name]]
    data_subsets_df <- split(data_for_group, by = "glucose")
    
    # 2. PRE-COMPUTE: Create the lean data structure
    lean_data_subsets <- lapply(data_subsets_df, function(d) {
      d <- d[order(hours)] # Ensure alignment
      t0_data <- d[hours == 0]
      list(
        hours = d$hours,
        alive = d$alive,
        dead = d$dead,
        glucose = d$glucose[1],
        ic_alive = if(nrow(t0_data) > 0) mean(t0_data$alive, na.rm=TRUE) else 0,
        ic_dead = if(nrow(t0_data) > 0) mean(t0_data$dead, na.rm=TRUE) else 0
      )
    })
    
    # --- 3. CREATE CLOSURE: Define the objective function here ---
    # It automatically captures 'model' and 'lean_data_subsets' from this environment.
    # The 'future' backend will handle sending them to the workers.
    objective_function_closure <- function(global_par_vec_log, par_names = par_bounds$par, weight_dead = 1, eps = 1, penalty = 1e12) {
      # Convert parameters from log-space to natural space
      global_par_vec <- exp(global_par_vec_log)
      names(global_par_vec) <- par_names
      total_sse <- 0
      
      for (d_sub in lean_data_subsets) {
        inner_obj_fun <- function(ic_par_vec) {
          et <- rxode2::eventTable()$add.sampling(d_sub$hours)
          inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = d_sub$glucose)
          sim <- try(model$solve(params = global_par_vec, events = et, inits = inits, cores = 1), silent = TRUE)
          
          if (inherits(sim, "try-error") || nrow(sim) != length(d_sub$hours)) return(penalty)
          
          r_alive <- log(d_sub$alive + eps) - log(sim$NL + eps)
          r_dead  <- log(d_sub$dead + eps) - log(sim$ND + eps)
          sse     <- sum(r_alive^2 + weight_dead * r_dead^2)
          return(ifelse(!is.finite(sse), penalty, sse))
        }
        ic_start_guess <- c(d_sub$ic_alive, d_sub$ic_dead)
        inner_opt <- optim(
          par = ic_start_guess, fn = inner_obj_fun,
          method = "Nelder-Mead", control = list(maxit = 50, reltol = 1e-4)
        )
        total_sse <- total_sse + inner_opt$value
      }
      return(total_sse)
    }
    
    
    # --- 4. OPTIMIZE ---
    de_optim_result <- NULL
    n_chunks <- ceiling(total_gens / n_gen_chunk)
    
    for(i in 1:n_chunks) {
      cat(sprintf("Running chunk %d of %d (Generations %d to %d)...\n",
                  i, n_chunks, (i-1)*n_gen_chunk + 1, i*n_gen_chunk))
      
      de_optim_result <- DEoptim(
        fn = objective_function_closure, # Use the closure
        lower = log_lower_bounds,
        upper = log_upper_bounds,
        control = DEoptim.control(
          NP = n_pop, itermax = n_gen_chunk,
          initialpop = if(i > 1) de_optim_result$member$pop else NULL,
          steptol = 50, strategy = 2, parallelType = "auto", trace = TRUE
        )
      )
      
      # --- Plotting ---
      cat("Generating intermediate plot...\n")
      best_pars_log <- de_optim_result$optim$bestmem
      best_pars_natural <- setNames(exp(best_pars_log), par_bounds$par)
      plot_title <- sprintf("Fit for %s - Gen %d\nSSE: %.2f",
                            group_name, i * n_gen_chunk, de_optim_result$optim$bestval)
      
      fit_traj <- create_fit_data_lean(best_pars_natural, model, lean_data_subsets)
      
      p <- ggplot() +
        geom_line(data = fit_traj, aes(x = hours, y = NL, color = factor(glucose)), linewidth = 1, alpha = 0.8) +
        geom_point(data = data_for_group, aes(x = hours, y = alive, color = factor(glucose))) +
        geom_line(data = fit_traj, aes(x = hours, y = ND, color = factor(glucose)), linetype = "dashed", linewidth = 1, alpha = 0.8) +
        geom_point(data = data_for_group, aes(x = hours, y = dead, color = factor(glucose)), shape = 2) +
        facet_wrap(~glucose, scales = "free_y") +
        labs(title = plot_title, x = "Hours", y = "Cell Count", color = "Initial Glucose") +
        theme_bw()
      
      # Save plot to the 'tmp/' subfolder
      plot_filename <- sprintf("tmp/fit_progress_%s_gen%04d.png", gsub("[ .]", "_", group_name), i * n_gen_chunk)
      ggsave(filename = plot_filename, plot = p, width = 10, height = 6)
    }
    
    de_optim_result$optim$bestmem <- setNames(exp(de_optim_result$optim$bestmem), par_bounds$par)
    de_optim_result$member$pop <- exp(de_optim_result$member$pop)
    return(de_optim_result)
  })
  
  names(all_results) <- names(groups_to_fit)
  return(all_results)
}


#' @title Lean Plotting Helper
#' @description Creates smooth trajectories using the lean data structure.
create_fit_data_lean <- function(best_pars_vec, model, lean_data_subsets) {
  all_sims <- rbindlist(lapply(lean_data_subsets, function(d_sub) {
    inner_obj_fun <- function(ic_par_vec) {
      et <- rxode2::eventTable()$add.sampling(d_sub$hours)
      inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = d_sub$glucose)
      sim <- try(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1), silent = TRUE)
      if (inherits(sim, "try-error") || nrow(sim) != length(d_sub$hours)) return(1e12)
      sse <- sum((log(d_sub$alive + 1) - log(sim$NL + 1))^2) + sum((log(d_sub$dead + 1) - log(sim$ND + 1))^2)
      return(sse)
    }
    ic_start_guess <- c(d_sub$ic_alive, d_sub$ic_dead)
    inner_opt <- optim(par = ic_start_guess, fn = inner_obj_fun, method = "Nelder-Mead", control = list(maxit = 50))
    best_ics <- inner_opt$par
    
    smooth_times <- seq(min(d_sub$hours), max(d_sub$hours), length.out = 200)
    et <- rxode2::eventTable()$add.sampling(smooth_times)
    inits <- c(NL = best_ics[1], ND = best_ics[2], G = d_sub$glucose)
    sim <- as.data.table(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1))
    setnames(sim, "time", "hours")
    sim[, glucose := d_sub$glucose]
    return(sim)
  }))
  return(all_sims)
}