library(data.table)
library(rxode2) # Replaces deSolve for high-performance ODE solving
library(future)
library(future.apply)
library(DEoptim)

#' @title DEoptim objective (Poisson counts + fluorescence via dt$ll_lum)
#' @description
#' Joint negative log-likelihood for an ODE with states `NL` (alive), `ND` (dead), `G` (glucose).
#' One ODE solve per glucose arm (`G0`). Initial conditions `N0` and `D0` are **shared across all arms**
#' and taken from rows with `hours == 0` in `dt$cells`. Initial glucose is fixed to the arm's design `G0`.
#' Counts are modeled independently as Poisson: `N ~ Pois(muN)`, `D ~ Pois(muD)`.
#' Fluorescence likelihood is supplied by `dt$ll_lum(G, dat)` where `dat` is `dt$glucose` filtered
#' to the arm/time; its parameters (e.g., `a, b, sdlog`) are already captured in the closure.
#'
#' @param global_par_vec_log Named numeric vector of log-parameters (same order as `par_names`).
#' @param model RxODE-like model with `$solve(params, events, inits, cores)` returning columns `time, NL, ND, G`.
#'        If `"r_cell"` is in the parameter vector, an extra state `R` is initialized as `R = N0 * r_cell`.
#' @param dt A **list** with:
#'   - `cells`: data.frame/data.table with columns
#'       * `G0` (character or numeric arm ID; if character it should coerce to numeric),
#'       * `hours` (numeric; includes 0),
#'       * `N`, `D` (numeric counts; optional—if absent, counts term is skipped).
#'   - `glucose`: data.frame/data.table with columns
#'       * `G0` (matching the arms in `cells`), `hours` (numeric),
#'       * `lum` (numeric fluorescence), and may include `Dilution Factor` (used by your `ll_lum`).
#'   - `ll_lum`: function `function(G, dat) -> density vector p(lum | G)`; called with `dat` filtered
#'       to a single arm/time. We take `log(p)` inside; non-finite logs incur a per-observation penalty.
#' @param par_names Character vector naming the parameters (must match `names(global_par_vec_log)` order).
#' @param eps Small positive number for numerical guards (default `1e-12`).
#' @param penalty Penalty added **per issue** (per failed row / non-finite term / failed timepoint) (default `1e12`).
#'
#' @return Scalar numeric: negative log-likelihood plus accumulated penalties (for minimization).
deoptim_objective_function <- function(
    global_par_vec_log, model, dt, par_names,
    eps = 1e-12, penalty = 1e12
){
  # unpack params (positivity)
  theta <- stats::setNames(exp(global_par_vec_log), par_names)
  
  # coerce tables
  cells <- data.table::as.data.table(dt$cells)
  glu   <- data.table::as.data.table(dt$glucose)
  ll_lum <- dt$ll_lum
  
  # shared ICs from t=0 (assumed exactly 0)
  N0 <- stats::median(cells[hours == 0]$N, na.rm = TRUE); if (!is.finite(N0)) N0 <- 0
  D0 <- stats::median(cells[hours == 0]$D, na.rm = TRUE); if (!is.finite(D0)) D0 <- 0
  
  # define arms and helper to coerce G0 to numeric for the ODE initial G
  arms <- sort(unique(c(cells$G0, glu$G0)))
  as_num <- function(x) suppressWarnings(as.numeric(x))
  
  total_ll <- 0.0
  penalty_count <- 0L
  
  for (g0 in arms) {
    g0_num <- as_num(g0)
    if (!is.finite(g0_num)) { penalty_count <- penalty_count + 1L; next }
    
    # union of sampling times from cells & glucose for this arm
    t_cells <- cells[G0 == g0, sort(unique(hours))]
    t_glu   <- glu[G0 == g0,   sort(unique(hours))]
    times   <- sort(unique(c(t_cells, t_glu)))
    if (!length(times)) next
    
    # ODE solve
    et <- rxode2::eventTable(); et$add.sampling(times)
    inits <- c(NL = N0, ND = D0, G = g0_num)
    if ("r_cell" %in% names(theta)) inits["R"] <- N0 * theta["r_cell"]
    
    sim <- try(model$solve(params = theta, events = et, inits = inits, cores = 1), silent = TRUE)
    if (inherits(sim, "try-error")) { penalty_count <- penalty_count + length(times); next }
    
    sim_dt <- data.table::as.data.table(sim)
    data.table::setnames(sim_dt, "time", "hours")
    
    ## --- counts term (Poisson for N and D) ---
    if (all(c("N","D","hours") %in% names(cells))) {
      mmC <- merge(cells[G0 == g0], sim_dt, by = "hours", all.x = TRUE, sort = FALSE)
      idx <- is.finite(mmC$N) & is.finite(mmC$D) & is.finite(mmC$NL) & is.finite(mmC$ND)
      if (any(idx)) {
        muN <- pmax(mmC$NL[idx], eps); muD <- pmax(mmC$ND[idx], eps)
        llN <- dpois(mmC$N[idx], lambda = muN, log = TRUE)
        llD <- dpois(mmC$D[idx], lambda = muD, log = TRUE)
        bad <- !is.finite(llN) | !is.finite(llD)
        penalty_count <- penalty_count + sum(bad)
        total_ll <- total_ll + sum(llN[!bad]) + sum(llD[!bad])
      }
    }
    
    ## --- fluorescence term (provided by dt$ll_lum) ---
    if (is.function(ll_lum) && all(c("hours","lum") %in% names(glu))) {
      # evaluate by time so each call uses the correct predicted G
      mmG <- merge(glu[G0 == g0], sim_dt, by = "hours", all.x = TRUE, sort = FALSE)
      # split indices by time to pass per-time data frame to ll_lum
      if (nrow(mmG)) {
        sp <- split(seq_len(nrow(mmG)), mmG$hours)
        for (ix in sp) {
          if (!all(is.finite(mmG$G[ix]))) { penalty_count <- penalty_count + length(ix); next }
          G_here <- unique(mmG$G[ix])
          if (length(G_here) != 1L || !is.finite(G_here)) { penalty_count <- penalty_count + length(ix); next }
          dat_here <- mmG[ix]
          # ll_lum returns densities p(lum|G); convert to log and penalize non-finites
          p_vec <- try(ll_lum(G_here, dat_here), silent = TRUE)
          if (inherits(p_vec, "try-error")) {
            penalty_count <- penalty_count + length(ix)
          } else {
            lp <- log(p_vec)
            bad <- !is.finite(lp)
            penalty_count <- penalty_count + sum(bad)
            total_ll <- total_ll + sum(lp[!bad])
          }
        }
      }
    }
  }
  
  -(total_ll) + penalty * penalty_count
}

## basically shadows the objective function defined above but returns data and/or (TBd) plots
plot_objective_function <- function(
    global_par_vec_log, model, dt, par_names,
    sample_freq=1 
){
  # unpack params (positivity)
  theta <- stats::setNames(exp(global_par_vec_log), par_names)
  
  # coerce tables
  cells <- data.table::as.data.table(dt$cells)
  glu   <- data.table::as.data.table(dt$glucose)
  ll_lum <- dt$ll_lum
  
  # shared ICs from t=0 (assumed exactly 0)
  N0 <- stats::median(cells[hours == 0]$N, na.rm = TRUE); if (!is.finite(N0)) N0 <- 0
  D0 <- stats::median(cells[hours == 0]$D, na.rm = TRUE); if (!is.finite(D0)) D0 <- 0
  
  # define arms and helper to coerce G0 to numeric for the ODE initial G
  arms <- sort(unique(c(cells$G0, glu$G0)))
  as_num <- function(x) suppressWarnings(as.numeric(x))
  
  
  
  do.call(rbind,lapply(arms,function(g0){
    g0_num <- as_num(g0)
    
    # union of sampling times from cells & glucose for this arm
    t_cells <- cells[G0 == g0, sort(unique(hours))]
    t_glu   <- glu[G0 == g0,   sort(unique(hours))]
    times   <- unique(c(t_cells, t_glu))
    times <- seq(0,max(times),by=sample_freq)
    
    if (!length(times)) next
    
    # ODE solve
    et <- rxode2::eventTable(); et$add.sampling(times)
    inits <- c(NL = N0, ND = D0, G = g0_num)
    if ("r_cell" %in% names(theta)) inits["R"] <- N0 * theta["r_cell"]
    
    sim <- try(model$solve(params = theta, events = et, inits = inits, cores = 1), silent = TRUE)
    if (inherits(sim, "try-error")) { penalty_count <- penalty_count + length(times); next }
    
    sim_dt <- data.table::as.data.table(sim)
    data.table::setnames(sim_dt, "time", "hours")
    sim_dt$G0 <- g0_num
    
    return(sim_dt)
  }))
}

#' @title Main Fitting Function using DEoptim (Log-Space Version)
#' @description Orchestrates model fitting in log-space using Differential Evolution.

fit_model_deoptim <- function(data, model,par_info,
                              n_pop = 10 * length(par_info$par),
                              total_gens = 200,
                              cores = parallel::detectCores()) {
  
  plan(multisession, workers = cores)
  
  de_optim_result <- DEoptim(
    fn = deoptim_objective_function,
    lower = log(par_info$lower), # Use log-transformed bounds
    upper = log(par_info$upper), # Use log-transformed bounds
    model = model,
    dt = data,
    par_names = par_info$par,
    control = DEoptim.control(
      NP = n_pop,
      itermax = total_gens,
      initialpop = NULL,
      steptol = 50,
      strategy = 2,
      parallelType = "auto",
      trace = TRUE
    )
  )
  plan(sequential)
  return(de_optim_result)
}


