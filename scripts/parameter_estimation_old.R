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
  
  d/dt(R) = NL*glucose_uptake_rate - NL*((1/Y_xr) * growth_rate + m_r + r_cell * (starvation_death_rate + waste_death_rate))
  
  d/dt(W) = (glucose_uptake_rate - deltaW * W)*NL # SIMPLIFIED: Waste integrates glucose use
  
})


deoptim_objective_function <- function(global_par_vec_log, model, data_subsets, par_names, weight_dead = 1, weight_glc = 1, eps = 1, penalty = 1e12) {
  
  global_par_vec <- exp(global_par_vec_log)
  names(global_par_vec) <- par_names
  
  total_sse <- 0
  
  for (d_sub in data_subsets) {
    
    inner_obj_fun <- function(ic_par_vec) {
      ic_par_vec <- abs(ic_par_vec)##prevent negative
      times_for_g <- sort(unique(d_sub$hours))
      et <- rxode2::eventTable()$add.sampling(times_for_g)
      # MODIFIED: G is now the 3rd parameter to fit (ic_par_vec[3])
      inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = ic_par_vec[3])
      
      if("r_cell"%in%names(global_par_vec)) inits["R"] <- inits["NL"]*global_par_vec["r_cell"]
      
      sim <- try(model$solve(params = global_par_vec, events = et, inits = inits, cores = 1), silent = TRUE)
      
      if (inherits(sim, "try-error")) return(penalty)
      
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      mm <- merge(d_sub, sim_dt, by = "hours")
      if (nrow(mm) != nrow(d_sub)) return(penalty)
      
      if (any(!is.finite(mm$NL)) || any(!is.finite(mm$ND)) || any(!is.finite(mm$G))) {
        return(penalty)
      }
      
      r_alive <- log(mm$alive + eps) - log(mm$NL + eps)
      r_dead  <- log(mm$dead + eps) - log(mm$ND + eps)
      r_glc   <- log(mm$glucose_measured + eps) - log(mm$G + eps)
      
      sse <- sum(r_alive^2, na.rm = TRUE) +
        weight_dead * sum(r_dead^2, na.rm = TRUE) +
        weight_glc * sum(r_glc^2, na.rm = TRUE)
      
      return(ifelse(!is.finite(sse), penalty, sse))
    }
    
    # MODIFIED: Add initial glucose from data as the 3rd start guess
    ic_start_guess <- c(d_sub[hours == 0, mean(alive)], d_sub[hours == 0, mean(dead)], d_sub$glucose[1])
    ic_start_guess[is.na(ic_start_guess)] <- 0
    
    inner_opt <- optim(
      par = ic_start_guess,
      fn = inner_obj_fun,
      method = "Nelder-Mead",
      control = list(maxit = 100, reltol = 1e-4) # Increased maxit for 3 params
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
                              weight_dead = 1, weight_glc = 1,
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
        weight_dead = weight_dead, # Pass weight
        weight_glc = weight_glc,   # Pass weight
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
        fit_traj <- create_fit_data(best_pars_natural, model, data_subsets, weight_dead = weight_dead, weight_glc = weight_glc)
        
        # --- Reshape data for faceted plotting ---
        sim_long <- melt(fit_traj, id.vars = c("hours", "glucose"), measure.vars = c("NL", "ND", "G"), variable.name = "type", value.name = "sim_value")
        sim_long[type %in% c("NL", "ND"), panel := "Cell Count"]
        sim_long[type == "G", panel := "Glucose (mM)"]
        
        exp_long <- melt(data_for_group, id.vars = c("hours", "glucose"), measure.vars = c("alive", "dead", "glucose_measured"), variable.name = "type", value.name = "exp_value")
        exp_long[type %in% c("alive", "dead"), panel := "Cell Count"]
        exp_long[type == "glucose_measured", panel := "Glucose (mM)"]
        
        p <- ggplot() +
          geom_line(data = sim_long, aes(x = hours, y = sim_value, color = factor(glucose), linetype = type), linewidth = 1) +
          geom_point(data = exp_long, aes(x = hours, y = exp_value, color = factor(glucose), shape = type), size = 2) +
          facet_wrap(vars(panel, glucose), scales = "free") +
          scale_linetype_manual(name = "Simulated", values = c("NL" = "solid", "ND" = "dashed", "G" = "solid")) +
          scale_shape_manual(name = "Measured", values = c("alive" = 16, "dead" = 17, "glucose_measured" = 15)) +
          labs(title = plot_title, x = "Hours", y = "Value", color = "Initial Glc") +
          theme_bw() +
          theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")
        
        ggsave(filename = sprintf("tmp/fit_progress_%s_gen%04d.png", gsub("[ .]", "_", group_name), i * n_gen_chunk),
               plot = p, width = 12, height = 7)
      }, error = function(e) cat(paste("Plot generation failed:", e$message, "\n")))
      
      
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

#' @title Create Fit Data and Optionally Return SSE
#' @description Solves the model and generates smooth trajectories for plotting.
#' @param best_pars_vec A named vector of model parameters.
#' @param model The compiled rxode2 model object.
#' @param data_subsets A list of data.tables, split by condition.
#' @param return_list A logical flag. If TRUE, the function returns a list
#'   containing the simulation data and the total SSE. Defaults to FALSE for
#'   backward compatibility.
#' @return A data.table of simulation results or a list containing the
#'   simulation data and the total SSE.
create_fit_data <- function(best_pars_vec, model, data_subsets, return_list = FALSE, weight_dead = 1, weight_glc = 1, eps = 1) {
  
  sim_and_sse_list <- lapply(data_subsets, function(d_sub) {
    
    inner_obj_fun <- function(ic_par_vec) {
      ic_par_vec <- abs(ic_par_vec)
      times <- sort(unique(d_sub$hours))
      et <- rxode2::eventTable()$add.sampling(times)
      # MODIFIED: G is now the 3rd parameter to fit (ic_par_vec[3])
      inits <- c(NL = ic_par_vec[1], ND = ic_par_vec[2], G = ic_par_vec[3])
      sim <- try(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1), silent = TRUE)
      
      if (inherits(sim, "try-error")) return(1e12)
      
      sim_dt <- as.data.table(sim)
      setnames(sim_dt, "time", "hours")
      mm <- merge(d_sub, sim_dt, by = "hours")
      
      if (nrow(mm) != nrow(d_sub)) return(1e12)
      
      if (any(!is.finite(mm$NL)) || any(!is.finite(mm$ND)) || any(!is.finite(mm$G))) {
        return(1e12)
      }
      
      r_alive <- log(mm$alive + eps) - log(mm$NL + eps)
      r_dead  <- log(mm$dead + eps) - log(mm$ND + eps)
      r_glc   <- log(mm$glucose_measured + eps) - log(mm$G + eps)
      
      sse <- sum(r_alive^2, na.rm = TRUE) +
        weight_dead * sum(r_dead^2, na.rm = TRUE) +
        weight_glc * sum(r_glc^2, na.rm = TRUE)
      
      return(ifelse(!is.finite(sse), 1e12, sse))
    }
    
    # MODIFIED: Add initial glucose from data as the 3rd start guess
    ic_start_guess <- c(d_sub[hours == 0, mean(alive)], d_sub[hours == 0, mean(dead)], d_sub$glucose[1])
    ic_start_guess[is.na(ic_start_guess)] <- 0
    inner_opt <- optim(par = ic_start_guess, fn = inner_obj_fun, method = "Nelder-Mead", control = list(maxit = 100))
    
    sse_for_subset <- inner_opt$value
    best_ics <- abs(inner_opt$par)
    
    smooth_times <- seq(min(d_sub$hours), max(d_sub$hours), length.out = 200)
    et <- rxode2::eventTable()$add.sampling(smooth_times)
    # MODIFIED: Use the optimized initial conditions for the smooth plot
    inits <- c(NL = best_ics[1], ND = best_ics[2], G = best_ics[3])
    sim <- as.data.table(model$solve(params = best_pars_vec, events = et, inits = inits, cores = 1))
    setnames(sim, "time", "hours")
    sim[, glucose := d_sub$glucose[1]] # Keep original glucose for labeling
    
    return(list(sim_data = sim, sse = sse_for_subset))
  })
  
  all_sims <- rbindlist(lapply(sim_and_sse_list, `[[`, "sim_data"))
  
  if (return_list) {
    total_sse <- sum(sapply(sim_and_sse_list, `[[`, "sse"))
    return(list(sim_data = all_sims, sse = total_sse))
  } else {
    return(all_sims)
  }
}

#' @title Plot Model Fit Against Experimental Data (with SSE)
#' @description Generates a ggplot object showing fits for cell counts and glucose.

plot_model_fit <- function(par_vec, data_subset, model, title = "Model Fit vs. Data",
                           weight_dead = 1, weight_glc = 1, eps = 1) {
  
  setDT(data_subset)
  
  has_glucose_data <- "glucose_measured" %in% names(data_subset) && 
    any(!is.na(data_subset$glucose_measured))
  
  data_subsets <- split(data_subset, by = "glucose")
  
  fit_results <- create_fit_data(
    best_pars_vec = par_vec, model = model, data_subsets = data_subsets,
    return_list = TRUE, weight_dead = weight_dead, weight_glc = weight_glc, eps = eps
  )
  fit_trajectory <- fit_results$sim_data
  sse_value <- fit_results$sse
  plot_title <- sprintf("%s\nTotal Weighted SSE: %.3f", title, sse_value)
  
  if (has_glucose_data) {
    sim_long <- melt(fit_trajectory, id.vars = c("hours", "glucose"), measure.vars = c("NL", "ND", "G"),
                     variable.name = "type", value.name = "sim_value")
    sim_long[type %in% c("NL", "ND"), panel := "Cell Count"]
    sim_long[type == "G", panel := "Glucose (mM)"]
    
    exp_long <- melt(data_subset, id.vars = c("hours", "glucose"), measure.vars = c("alive", "dead", "glucose_measured"),
                     variable.name = "type", value.name = "exp_value")
    exp_long[type %in% c("alive", "dead"), panel := "Cell Count"]
    exp_long[type == "glucose_measured", panel := "Glucose (mM)"]
    
    p <- ggplot() +
      geom_line(data = sim_long, aes(x = hours, y = sim_value, color = factor(glucose), linetype = type), linewidth = 1) +
      geom_point(data = exp_long, aes(x = hours, y = exp_value, color = factor(glucose), shape = type), size = 2.5, alpha = 0.8) +
      # MODIFIED: Use facet_wrap for fully independent panel scales
      facet_wrap(vars(panel, glucose), scales = "free") +
      scale_linetype_manual(name = "Simulated", values = c("NL" = "solid", "ND" = "dashed", "G" = "solid")) +
      scale_shape_manual(name = "Measured", values = c("alive" = 16, "dead" = 17, "glucose_measured" = 15)) +
      labs(title = plot_title, x = "Hours", y = "Value", color = "Initial Glc") +
      theme_bw() + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
    
  } else {
    # Fallback plot for cell counts only
    p <- ggplot() +
      geom_line(data = fit_trajectory, aes(x = hours, y = NL, color = factor(glucose)), linewidth = 1) +
      geom_point(data = data_subset, aes(x = hours, y = alive, color = factor(glucose))) +
      geom_line(data = fit_trajectory, aes(x = hours, y = ND, color = factor(glucose)), linetype = "dashed") +
      geom_point(data = data_subset, aes(x = hours, y = dead, color = factor(glucose)), shape = 2) +
      # MODIFIED: Use facet_wrap here as well for consistency
      facet_wrap(~glucose, scales = "free") +
      labs(title = plot_title, x = "Hours", y = "Cell Count (Alive/Dead)", color = "Initial Glucose") +
      theme_bw()
  }
  
  return(p)
}
