library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)

# --- 1. Robust Log-Spline Helper (Unchanged) ---
get_log_spline_derivs <- function(df, time_res = 0.5, sm_par = 0.6) {
  if(nrow(df) < 5) return(NULL)
  
  epsilon <- 1e-4
  log_y <- log(df$val + epsilon)
  
  fit <- smooth.spline(x = df$hours, y = log_y, spar = sm_par)
  t_seq <- seq(min(df$hours), max(df$hours), by = time_res)
  
  pred <- predict(fit, t_seq)       
  pred_deriv <- predict(fit, t_seq, deriv = 1) 
  
  y_est <- exp(pred$y)
  dy_dt <- y_est * pred_deriv$y
  
  data.frame(hours = t_seq, X = y_est, dX = dy_dt)
}

# --- 2. Process File (Unchanged) ---
process_ploidy_file <- function(file_path, ploidy_label, sm_par) {
  
  dat <- readRDS(file_path)
  env <- environment(dat$ll_lum)
  a_est <- env$a; b_est <- env$b
  
  # Prepare Data
  raw_G <- dat$glucose %>%
    mutate(
      dilution_multiplier = 1000 / `Dilution Factor`,
      val = (lum - b_est) / (a_est * dilution_multiplier),
      val = pmax(1e-3, val), 
      type = "G"
    ) %>% select(G0, hours, val, type)
  
  raw_cells <- dat$cells %>%
    mutate(A = N - D, G0 = as.numeric(G0)) %>%
    pivot_longer(cols = c(A, D), names_to = "type", values_to = "val") %>%
    mutate(val = pmax(1, val)) %>% 
    select(G0, hours, val, type)
  
  all_raw <- bind_rows(raw_G, raw_cells)
  
  # Fit Splines
  fitted <- all_raw %>%
    group_by(G0, type) %>%
    group_modify(~ get_log_spline_derivs(.x, time_res = 0.5, sm_par = sm_par)) %>%
    ungroup() %>%
    mutate(ploidy = ploidy_label)
  
  # Calculate Metrics (Wide format for rates)
  fitted_wide <- fitted %>%
    pivot_wider(names_from = type, values_from = c(X, dX), names_glue = "{type}_{.value}") %>%
    rename(A = A_X, dA = A_dX, D = D_X, dD = D_dX, G = G_X, dG = G_dX) %>%
    mutate(
      birth_rate = (dA + dD) / pmax(1, A), 
      death_rate = dD / pmax(1, A),
      spec_consumption = -dG / pmax(1, A) 
    )
  
  return(list(raw = all_raw %>% mutate(ploidy = ploidy_label), fit = fitted_wide))
}

# --- 3. Main Dashboard Generator (Updated) ---
generate_dashboards <- function(cell_line, sm_par = 0.6) {
  
  # Locate and Process Files
  files <- list.files("data/model_inputs", pattern = paste0("^", cell_line, "\\."), full.names = T)
  if(length(files) < 2) {
    message(paste("Skipping", cell_line, "- missing pair."))
    return(NULL)
  }
  
  file_high <- files[grep("\\.high\\.|\\.4N\\.", files)]
  file_low  <- files[grep("\\.low\\.|\\.2N\\.|\\.parental\\.", files)]
  
  res_high <- process_ploidy_file(file_high, "High Ploidy", sm_par)
  res_low  <- process_ploidy_file(file_low, "Low Ploidy", sm_par)
  
  # Combine Data
  dat_fit <- bind_rows(res_high$fit, res_low$fit)
  dat_raw_all <- bind_rows(res_high$raw, res_low$raw)
  
  # Convert fits to long format for time-course plotting
  dat_fit_long <- dat_fit %>%
    select(hours, G0, ploidy, A, D, G) %>%
    pivot_longer(cols = c(A, D, G), names_to = "type", values_to = "val")
  
  # Arrow data for phase planes
  arrow_dat <- dat_fit %>% group_by(ploidy, G0) %>% slice(seq(1, n(), by = 10))
  
  # Theme
  base_theme <- theme_bw() + 
    theme(legend.position = "bottom",
          strip.background = element_rect(fill="grey95"),
          plot.title = element_text(size = 11, face="bold"),
          axis.title = element_text(size = 9))
  
  # ---------------------------------------------------------
  # DASHBOARD 1: TIME COURSE FITS (Updated Layout)
  # ---------------------------------------------------------
  
  # 1. Glucose Plot (Left Column)
  p_gluc <- ggplot() +
    geom_point(data = filter(dat_raw_all, type == "G"), 
               aes(x = hours, y = val), color = "blue", alpha = 0.3) +
    geom_line(data = filter(dat_fit_long, type == "G"), 
              aes(x = hours, y = val), color = "blue", size = 1) +
    # Rows = G0, Cols = Ploidy
    facet_grid(G0 ~ ploidy, scales = "free", labeller = label_both) +
    labs(title = "Glucose Fits", 
         y = "[Glucose] (mM)", x = "Time (h)") +
    base_theme
  
  # 2. Cells Plot (Right Column)
  p_cells <- ggplot() +
    geom_point(data = filter(dat_raw_all, type %in% c("A", "D")), 
               aes(x = hours, y = val, color = type), alpha = 0.3) +
    geom_line(data = filter(dat_fit_long, type %in% c("A", "D")), 
              aes(x = hours, y = val, color = type), size = 1) +
    # Rows = G0, Cols = Ploidy
    facet_grid(G0 ~ ploidy, scales = "free", labeller = label_both) +
    scale_color_manual(values = c("A" = "#2E8B57", "D" = "#CD5C5C"), 
                       labels = c("A" = "Alive", "D" = "Dead")) +
    labs(title = "Cell Count Fits", 
         y = "Count", x = "Time (h)") +
    base_theme
  
  # Combine into one grid (2 Columns)
  dash_fits <- arrangeGrob(p_gluc, p_cells, ncol = 2,
                           top = textGrob(paste(cell_line, "- Spline Fits (spar =", sm_par, ")"), 
                                          gp=gpar(fontsize=14)))
  
  # ---------------------------------------------------------
  # DASHBOARD 2: RATES (Phase Planes)
  # ---------------------------------------------------------
  
  make_rate_plot <- function(data, y_metric, title, y_lab, global_limits) {
    ggplot(data, aes(x=G, y=.data[[y_metric]], group=G0, color=hours)) +
      geom_path(alpha=0.7, size=0.8, na.rm=TRUE) + 
      geom_path(data=filter(arrow_dat, ploidy==data$ploidy[1]), 
                aes(y=.data[[y_metric]]),
                arrow = arrow(length=unit(0.15,"cm"), type="closed"), size=0.5, na.rm=TRUE) +
      scale_color_viridis_c() + 
      scale_x_log10() +
      coord_cartesian(ylim = global_limits) +
      labs(title=title, y=y_lab, x="Log [G] (mM)") + base_theme +
      theme(legend.position = "none") # Hide legend for cleanliness in grid
  }
  
  # Limits
  get_limits <- function(metric) {
    q <- quantile(dat_fit[[metric]], c(0.05, 0.95), na.rm=TRUE)
    range <- diff(q)
    c(q[1] - 0.1*range, q[2] + 0.1*range)
  }
  
  lim_birth <- get_limits("birth_rate")
  lim_death <- get_limits("death_rate")
  lim_cons  <- get_limits("spec_consumption")
  
  # Row 1: Low Ploidy
  low_dat <- filter(dat_fit, ploidy=="Low Ploidy")
  r1 <- make_rate_plot(low_dat, "birth_rate", "Low: Birth Rate", "Birth (1/h)", lim_birth)
  r2 <- make_rate_plot(low_dat, "death_rate", "Low: Death Rate", "Death (1/h)", lim_death)
  r3 <- make_rate_plot(low_dat, "spec_consumption", "Low: Spec. Consump", "mM/cell/h", lim_cons)
  
  # Row 2: High Ploidy
  high_dat <- filter(dat_fit, ploidy=="High Ploidy")
  r4 <- make_rate_plot(high_dat, "birth_rate", "High: Birth Rate", "Birth (1/h)", lim_birth)
  r5 <- make_rate_plot(high_dat, "death_rate", "High: Death Rate", "Death (1/h)", lim_death)
  r6 <- make_rate_plot(high_dat, "spec_consumption", "High: Spec. Consump", "mM/cell/h", lim_cons)
  
  dash_rates <- arrangeGrob(r1, r2, r3, r4, r5, r6, nrow = 2, 
                            top = textGrob(paste(cell_line, "- Kinetic Phase Planes"), gp=gpar(fontsize=14)))
  
  return(list(rates = dash_rates, fits = dash_fits))
}


# --- Usage ---
lines <- c("MCF10A","MDA-MB-231","SNU668","SUM-159-chem","SUM-159-fuse")

for(line in lines){
  plots <- generate_dashboards(line, sm_par = 0.5)
  
  if(!is.null(plots)) {
    # 1. Phase Planes
    grid.newpage()
    grid.draw(plots$rates)
    
    # 2. Spline Fits
    grid.newpage()
    grid.draw(plots$fits)
  }
}

