library(data.table)
optimize_counts <- function(dt, method="coordinate_descent",
                            lambda_N_mono=1, lambda_D_mono=1, lambda_create=1,
                            max_iter=10, tol=0.01, grid_radius=5, ...) {
  if (!"time_numeric" %in% names(dt)) {
    dt[, time_numeric := {
      m <- str_match(image_key, "_(\\d{2})d(\\d{2})h\\d{2}m$")
      pmax(as.numeric(m[,2]),0) + pmax(as.numeric(m[,3]),0)/24
    }]
  }
  setorder(dt, time_numeric)
  times <- unique(dt$time_numeric); T <- length(times)
  init <- dt[, .(N = sum(predicted_label_name %in% c("alive","dead"), na.rm=TRUE),
                 D = sum(predicted_label_name=="dead", na.rm=TRUE)),
             by=time_numeric][order(time_numeric)]
  ref <- preprocess_dt(dt$time_numeric, dt$prob_alive, dt$prob_dead, dt$prob_junk)
  cost_fn <- function(p) {
    p <- round(p); N <- p[1:T]; D <- p[(T+1):(2*T)]
    D <- pmin(pmax(D,0), N); N <- pmax(N,0)
    calculate_total_cost(ref, times, N, D,
                         lambda_N_mono, lambda_D_mono, lambda_create)
  }
  current <- c(init$N, init$D); cost <- cost_fn(current)
  if (tolower(method) == "coordinate_descent") {
    message("Using coordinate descent")
    for (it in 1:max_iter) {
      prev <- cost
      for (i in sample(T)) {
        best_cost <- cost; best_par <- current
        ni <- i; di <- T + i
        for (dN in -grid_radius:grid_radius) for (dD in -grid_radius:grid_radius) {
          if (dN==0 && dD==0) next
          tmp <- current
          tmp[ni] <- tmp[ni] + dN; tmp[di] <- tmp[di] + dD
          if (tmp[ni]<0 || tmp[di]<0 || tmp[di]>tmp[ni]) next
          cst <- cost_fn(tmp)
          if (cst < best_cost) { best_cost <- cst; best_par <- tmp }
        }
        current <- best_par; cost <- best_cost
      }
      message(sprintf("Iter %d cost: %.2f", it, cost))
      if ((prev - cost) < tol) break
    }
    best_par <- round(current)
  } else {
    message(sprintf("Using optim(method='%s')", method))
    best_par <- round(optim(par=current, fn=cost_fn, method=method, ...)$par)
  }
  N_best <- best_par[1:T]
  D_best <- pmin(best_par[(T+1):(2*T)], N_best)
  list(final=data.table(time_numeric=times, N=N_best, D=D_best),
       init=init)
}