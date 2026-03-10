library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source("ecology/model_selection/models/gpath/gpath.R")

# 1. Generate Data (Simulations +  Observations)
generate_line_data <- function(model_id) {
  # Parse Dimensions
  dims <- parse_run_id(model_id)
  R <- dims$R; P <- dims$P; W <- dims$W; C <- dims$C; M <- dims$M
  strict_spec <- (C == 1)
  L <- 3 * R + (P - 1) * R + W * R + 1
  
  # Load Data and map groups
  stan_data <- readRDS("data/stan_ready_data.Rds")
  grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
  unique_grps <- unique(grp_keys)
  stan_data$N_groups <- length(unique_grps)
  stan_data$group_id <- match(grp_keys, unique_grps)
  
  times <- seq(0, max(stan_data$t_grid), 0.5)
  # Load fits
  
  chain_dir <- file.path("ecology/model_selection/data/gpath", model_id, "hier/optim_draws_all.Rds")
  fits <- readRDS(chain_dir)
  fits <- do.call(rbind,fits)
  
  pars <- fits[which.max(fits[,"log_lik"]),]
  
  # Helper to simulate a single well
  run_sim <- function(draw_vals, well_idx) {
    group_id <- stan_data$group_id[well_idx]
    p_metric <- stan_data$ploidy_metric[well_idx]
    
    target_line_id <- stan_data$line_id[well_idx]
    
    N0_raw <- draw_vals[sprintf("raw_N0[%d]", group_id)]
    N0 <- exp(log(500) + N0_raw * 1.0) 
    G1_0 <- draw_vals[sprintf("G1_0[%d]", stan_data$g1_id[well_idx])]
    
    raw_theta_line <- draw_vals[sprintf("raw_theta_line[%d,%d]", 1:L, target_line_id)]
    raw_theta_ploidy <- draw_vals[sprintf("raw_theta_ploidy[%d]", 1:L)]
    
    parms <- reconstruct_parms(
      R = R, P = P, W = W, strict_spec = strict_spec, M = M, 
      base_priors = base_priors, raw_theta_line = raw_theta_line, 
      raw_theta_ploidy = raw_theta_ploidy, ploidy_metric = p_metric
    )
    
    y0 <- c(N = max(as.numeric(N0), 1e-6), R1 = as.numeric(G1_0))
    if (R > 1) y0 <- c(y0, setNames(pmax(rep(1.0, R-1), 0), paste0("R", 2:R)))
    if (W > 0) y0 <- c(y0, setNames(rep(0.0, W), paste0("W", 1:W)))
    
    out <- as.data.frame(deSolve::ode(y = y0, times = times, func = rhs, parms = parms, method = "lsoda"))
    
    out_long <- tidyr::pivot_longer(out, cols = -time, names_to = "Variable", values_to = "Value")
    out_long$well_idx <- well_idx
    out_long$ploidy_abs <- stan_data$ploidy_abs[well_idx]
    out_long$G0 <- stan_data$G0_per_well[well_idx]
    out_long$line_id <- stan_data$line_id[well_idx]
    return(out_long)
  }
  
  
  sim_df <- do.call(rbind,pbapply::pblapply(1:stan_data$N_wells,function(i){
    run_sim(pars, i)
  }))
  
  
  
  # --- 1C. Experimental Observations ---
  obsN <- data.frame(
    well_idx = stan_data$well_idx_count,
    time = stan_data$t_grid[stan_data$grid_idx_count],
    Value = as.numeric(stan_data$N_obs),
    Variable = "N"
  )
  
  exps_G <- stan_data$exp_id[stan_data$well_idx_gluc]
  
  lum <- as.numeric(stan_data$lum_obs)
  dil <- as.numeric(stan_data$dilution)
  obsG <- data.frame(
    well_idx = stan_data$well_idx_gluc,
    time = stan_data$t_grid[stan_data$grid_idx_gluc],
    Value = pmax(0, (lum - stan_data$calib_b_fixed[exps_G]) / (stan_data$calib_a_fixed[exps_G] * dil + 1e-12)),
    Variable = "R1"
  )
  
  df_obs <- rbind(obsN, obsG)
  df_obs$ploidy_abs <- stan_data$ploidy_abs[df_obs$well_idx]
  df_obs$G0 <- stan_data$G0_per_well[df_obs$well_idx]
  df_obs$line_id <- stan_data$line_id[df_obs$well_idx]
  
  return(list(sim_df=sim_df,obs_df=df_obs))
}

# 2. Plot Data
plot_line_trajectories <- function(line_data,line_id,model_id) {
  stan_data <- readRDS("data/stan_ready_data.Rds")
  line_name <- names(stan_data$line_map)[line_id]
  sim_df <- line_data$sim_df[line_data$sim_df$line_id==line_id,]
  obs_df <- line_data$obs_df[line_data$obs_df$line_id==line_id,]
  
  
  # Format G0 levels for faceting
  g0_levels <- sort(unique(sim_df$G0))
  g0_labels <- paste0("G[0]=", g0_levels)
  
  sim_df$G0_label <- factor(paste0("G[0]=", sim_df$G0), levels = g0_labels)
  obs_df$G0_label <- factor(paste0("G[0]=", obs_df$G0), levels = g0_labels)
  
  vars_to_plot <- unique(sim_df$Variable)
  
  plot_list <- lapply(vars_to_plot, function(v) {
    sub_mean <- sim_df[sim_df$Variable == v, ]
    sub_obs <- obs_df[obs_df$Variable == v, ]
    
    title_str <- v
    if (v == "N") title_str <- "N (alive cells)"
    if (v == "R1") title_str <- "R1: Glucose"
    
    p <- ggplot() +
      geom_line(data = sub_mean, aes(x = time, y = Value), linewidth = 1, color = "dodgerblue4")
    
    if (nrow(sub_obs) > 0) {
      p <- p + geom_point(data = sub_obs, aes(x = time, y = Value), color = "black", size = 1.2, alpha = 0.8)
    }
    
    p + facet_grid(rows = vars(G0_label), cols = vars(ploidy_abs), scales = "free") +
      theme_minimal() +
      labs(title = title_str, x = "Time", y = "") +
      theme(
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
  })
  
  combined_plot <- wrap_plots(plot_list, nrow = 1) + 
    plot_annotation(
      title = sprintf("Model: %s | Cell Line: %s", model_id, line_name),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  print(combined_plot)
  return(invisible(combined_plot))
}

prep_line_diag_exact <- function(line_data, line_id, model_id) {
  stan_data <- readRDS("data/stan_ready_data.Rds")
  fits <- readRDS(file.path("ecology/model_selection/data/gpath", model_id, "hier/optim_draws_all.Rds"))
  fits <- do.call(rbind, fits)
  pars <- fits[which.max(fits[, "log_lik"]), ]
  
  line_name <- names(stan_data$line_map)[line_id]
  sigma_N <- exp(log(0.2) + as.numeric(pars["raw_sigma_N"]) * 0.5)
  
  sim_df <- line_data$sim_df[line_data$sim_df$line_id == line_id, ]
  pred_df <- sim_df[, c("well_idx", "time", "Variable", "Value", "ploidy_abs", "G0", "line_id")]
  colnames(pred_df)[colnames(pred_df) == "Value"] <- "Pred"
  
  obsN <- data.frame(
    well_idx = stan_data$well_idx_count,
    time = stan_data$t_grid[stan_data$grid_idx_count],
    Variable = "N",
    ObsRaw = as.numeric(stan_data$N_obs),
    DisplayValue = as.numeric(stan_data$N_obs),
    obs_type = "count",
    exp_id = NA_integer_,
    dilution = NA_real_,
    is_censored = 0L,
    sigma = sigma_N
  )
  
  exps_G <- stan_data$exp_id[stan_data$well_idx_gluc]
  lum <- as.numeric(stan_data$lum_obs)
  dil <- as.numeric(stan_data$dilution)
  a <- stan_data$calib_a_fixed[exps_G]
  b <- stan_data$calib_b_fixed[exps_G]
  sigG <- stan_data$calib_sigma_fixed[exps_G]
  
  obsG <- data.frame(
    well_idx = stan_data$well_idx_gluc,
    time = stan_data$t_grid[stan_data$grid_idx_gluc],
    Variable = "R1",
    ObsRaw = lum,
    DisplayValue = pmax(0, (lum - b) / (a * dil + 1e-12)),
    obs_type = "glucose",
    exp_id = exps_G,
    dilution = dil,
    is_censored = as.integer(stan_data$is_censored),
    sigma = sigG
  )
  
  obs_df <- dplyr::bind_rows(obsN, obsG)
  obs_df$ploidy_abs <- stan_data$ploidy_abs[obs_df$well_idx]
  obs_df$G0 <- stan_data$G0_per_well[obs_df$well_idx]
  obs_df$line_id <- stan_data$line_id[obs_df$well_idx]
  obs_df <- obs_df[obs_df$line_id == line_id, ]
  
  dd <- dplyr::left_join(obs_df, pred_df, by = c("well_idx", "time", "Variable", "ploidy_abs", "G0", "line_id"))
  
  dd$mu_lum <- NA_real_
  dd$ll_i <- NA_real_
  dd$z_i <- NA_real_
  
  iN <- dd$Variable == "N"
  Nhat <- pmax(dd$Pred[iN], 1e-6)
  dd$ll_i[iN] <- dlnorm(dd$ObsRaw[iN] + 1.0, meanlog = log(Nhat + 1.0), sdlog = sigma_N, log = TRUE)
  dd$z_i[iN] <- (log(dd$ObsRaw[iN] + 1.0) - log(Nhat + 1.0)) / sigma_N
  
  iG <- dd$Variable == "R1"
  mu <- stan_data$calib_a_fixed[dd$exp_id[iG]] * pmax(dd$Pred[iG], 0.0) * dd$dilution[iG] + stan_data$calib_b_fixed[dd$exp_id[iG]]
  mu <- pmax(mu, 1e-12)
  dd$mu_lum[iG] <- mu
  
  iG_unc <- iG & dd$is_censored == 0
  iG_cen <- iG & dd$is_censored == 1
  
  dd$ll_i[iG_unc] <- dlnorm(dd$ObsRaw[iG_unc], meanlog = log(dd$mu_lum[iG_unc]), sdlog = dd$sigma[iG_unc], log = TRUE)
  dd$ll_i[iG_cen] <- plnorm(dd$ObsRaw[iG_cen], meanlog = log(dd$mu_lum[iG_cen]), sdlog = dd$sigma[iG_cen], log.p = TRUE)
  dd$z_i[iG_unc] <- (log(pmax(dd$ObsRaw[iG_unc], 1e-12)) - log(dd$mu_lum[iG_unc])) / dd$sigma[iG_unc]
  
  g0_levels <- sort(unique(sim_df$G0))
  g0_labels <- paste0("G[0]=", g0_levels)
  sim_df$G0_label <- factor(paste0("G[0]=", sim_df$G0), levels = g0_labels)
  dd$G0_label <- factor(paste0("G[0]=", dd$G0), levels = g0_labels)
  
  facet_ll <- dd %>%
    dplyr::group_by(G0_label, ploidy_abs) %>%
    dplyr::summarize(facet_ll = sum(ll_i, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(label = sprintf("LL = %d", round(facet_ll)))
  
  list(
    sim_df = sim_df,
    obs_df = dd,
    facet_ll = facet_ll,
    line_ll = sum(dd$ll_i, na.rm = TRUE),
    global_ll = as.numeric(pars["log_lik"]),
    line_name = line_name
  )
}
plot_line_trajectories_llcolor <- function(line_data, line_id, model_id) {
  dd <- prep_line_diag_exact(line_data, line_id, model_id)
  sim_df <- dd$sim_df
  obs_df <- dd$obs_df
  facet_ll <- dd$facet_ll
  
  vars_to_plot <- unique(sim_df$Variable)
  
  plot_list <- lapply(vars_to_plot, function(v) {
    sub_sim <- sim_df[sim_df$Variable == v, ]
    sub_obs <- obs_df[obs_df$Variable == v, ]
    title_str <- if (v == "N") "N (alive cells)" else if (v == "R1") "R1: Glucose" else v
    
    ggplot() +
      geom_line(data = sub_sim, aes(time, Value), linewidth = 1, color = "dodgerblue4") +
      geom_point(data = sub_obs, aes(time, DisplayValue, color = ll_i, shape = factor(is_censored)), size = 1.8, alpha = 0.9) +
      geom_text(data = facet_ll, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3.4) +
      facet_grid(rows = vars(G0_label), cols = vars(ploidy_abs), scales = "free") +
      scale_shape_manual(values = c(`0` = 16, `1` = 1), guide = "none") +
      scale_color_viridis_c(name = "LL contrib", option = "C") +
      theme_minimal() +
      labs(title = title_str, x = "Time", y = "") +
      theme(strip.background = element_rect(fill = "grey90", color = NA),
            strip.text = element_text(face = "bold"),
            panel.border = element_rect(color = "grey80", fill = NA))
  })
  
  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = sprintf("Model: %s | Cell Line: %s | line LL = %d | global LL = %d", model_id, dd$line_name, round(dd$line_ll), round(dd$global_ll)),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  print(combined_plot)
  invisible(combined_plot)
}

plot_line_residuals_z <- function(line_data, line_id, model_id) {
  dd <- prep_line_diag_exact(line_data, line_id, model_id)
  obs_df <- dd$obs_df
  facet_ll <- dd$facet_ll
  
  vars_to_plot <- unique(obs_df$Variable)
  
  plot_list <- lapply(vars_to_plot, function(v) {
    sub_obs <- obs_df[obs_df$Variable == v & !is.na(obs_df$z_i), ]
    title_str <- if (v == "N") "N residual z-scores" else if (v == "R1") "R1 residual z-scores" else paste(v, "z-scores")
    
    ggplot(sub_obs, aes(time, z_i)) +
      geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey40") +
      geom_hline(yintercept = c(-2, 2), linewidth = 0.4, linetype = 3, color = "grey60") +
      geom_point(aes(color = ll_i), size = 1.8, alpha = 0.9) +
      geom_text(data = facet_ll, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3.4) +
      facet_grid(rows = vars(G0_label), cols = vars(ploidy_abs)) +
      scale_color_viridis_c(name = "LL contrib", option = "C") +
      theme_minimal() +
      labs(title = title_str, x = "Time", y = "z score") +
      theme(strip.background = element_rect(fill = "grey90", color = NA),
            strip.text = element_text(face = "bold"),
            panel.border = element_rect(color = "grey80", fill = NA))
  })
  
  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = sprintf("Model: %s | Cell Line: %s | line LL = %d | global LL = %d", model_id, dd$line_name, round(dd$line_ll), round(dd$global_ll)),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  print(combined_plot)
  invisible(combined_plot)
}


model_id <- c(  "1R_1P_0W_C0_M1","1R_1P_1W_C0_M1","2R_1P_0W_C0_M1",
                 "2R_2P_0W_C0_M1","2R_2P_1W_C0_M1")[4]
line_data <- generate_line_data(model_id)

p <- lapply(1:5,function(line_id){
#  plot_line_trajectories_llcolor(line_data, line_id,model_id)
  plot_line_residuals_z(line_data, line_id,model_id)
})
