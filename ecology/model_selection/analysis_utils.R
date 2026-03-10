library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggh4x)

# ==============================================================================
# 1. SETUP & HELPER FUNCTIONS
# ==============================================================================

# Ensure ModelPathways is available. 
# If it is defined in "ecology/pathways_model/pathways.R", source it here:
if(file.exists("ecology/pathways_model/pathways.R")) {
  source("ecology/pathways_model/pathways.R")
} else {
  stop("Please ensure ModelPathways is defined or sourced.")
}

# Helper to parse Line and Ploidy from directory path
# Assumes structure: .../LineName/PloidyLabel (e.g. .../MCF10A/hi)
parse_metadata <- function(path) {
  parts <- unlist(strsplit(path, "/"))
  n <- length(parts)
  list(
    line_name = parts[n-1],
    ploidy_lbl = parts[n]
  )
}

# ==============================================================================
# 2. DATA PROCESSING ENGINE
# ==============================================================================

process_directory <- function(fit_dir, model_obj) {
  
  # --- A. Validation ---
  d_path  <- file.path(fit_dir, "stan_data.Rds")
  lp_path <- file.path(fit_dir, "optim_lp_all.Rds")
  
  if(!file.exists(d_path) || !file.exists(lp_path)) return(NULL)
  
  # --- B. Load Data & Metadata ---
  meta <- parse_metadata(fit_dir)
  d    <- readRDS(d_path)
  
  # --- C. Load Best Fit Parameters ---
  lp    <- readRDS(lp_path)
  rc    <- readRDS(file.path(fit_dir, "optim_rc_all.Rds"))
  draws <- readRDS(file.path(fit_dir, "optim_draws_all.Rds"))
  
  # Find best successful optimization
  ok <- which(rc == 0L & is.finite(lp))
  if (!length(ok)) {
    warning(paste("No successful fits in:", fit_dir))
    return(NULL)
  }
  best_i <- ok[which.max(lp[ok])]
  
  best_draw_mat <- draws[[best_i]]
  # Extract the single best row
  row_vals_raw <- as.numeric(best_draw_mat[1, ]) 
  names(row_vals_raw) <- colnames(best_draw_mat)
  
  # --- Map 'theta' Array to Named Parameters ---
  p_list <- list()
  
  # 1. Biological Parameters
  for(i in seq_along(model_obj$param_names)) {
    p_name <- model_obj$param_names[i]
    col_key <- paste0("theta[", i, "]") 
    
    if(col_key %in% names(row_vals_raw)) {
      p_list[[p_name]] <- row_vals_raw[[col_key]]
    } else {
      if(p_name %in% names(row_vals_raw)) {
        p_list[[p_name]] <- row_vals_raw[[p_name]]
      } else {
        warning(paste("Parameter missing:", p_name))
        p_list[[p_name]] <- NA
      }
    }
  }
  
  # 2. Add N0 (Global Parameter)
  if("N0" %in% names(row_vals_raw)) {
    p_list[["N0"]] <- row_vals_raw[["N0"]]
  }
  
  final_params_vec <- unlist(p_list)
  
  master_params <- tibble(
    param = names(final_params_vec),
    value = as.numeric(final_params_vec),
    line_name = meta$line_name,
    ploidy = meta$ploidy_lbl
  )
  
  # --- D. SIMULATION (Vectorized over Wells) ---
  sim_results <- list()
  
  # Prepare params list for ODE solver
  p_sim <- as.list(final_params_vec)
  
  if(!is.null(model_obj$data_names)) {
    for(dn in model_obj$data_names) {
      if(!is.null(d[[dn]])) p_sim[[dn]] <- d[[dn]]
    }
  }
  
  for(w in 1:d$N_wells) {
    # N0 handling (Global)
    n0_val <- if("N0" %in% names(p_sim)) p_sim[["N0"]] else 0.1
    
    # State Init
    y0 <- setNames(rep(0, length(model_obj$state_names)), model_obj$state_names)
    
    #    Apply Model-Specific Fixed Defaults 
    #    e.g., if model_obj$init_defaults is list(x1=1, x2=1)
    if (!is.null(model_obj$state_defaults)) {
      # Intersection ensures we only update states that actually exist in this model
      overlap <- intersect(names(model_obj$state_defaults), names(y0))
      y0[overlap] <- unlist(model_obj$state_defaults[overlap])
    }
    
    # R2_0 is passed in stan_data
    r2_val <- if(!is.null(d$R2_0)) d$R2_0 else 1.0
    
    # --- CHANGED: Use Estimated G1_0 for Initial Condition ---
    # Construct key for this well, e.g., "G1_0[1]"
    g_key <- paste0("G1_0[", w, "]")
    
    # Check if the estimated value exists in the draw
    if(g_key %in% names(row_vals_raw)) {
      est_g1_0 <- row_vals_raw[[g_key]]
    } else {
      # Fallback to nominal if missing (shouldn't happen with correct Stan output)
      est_g1_0 <- d$G0_per_well[w]
    }
    
    y0["N"]  <- n0_val
    y0["G1"] <- est_g1_0  # <--- Using estimated value here
    if("G2"%in%model_obj$state_names) y0["G2"] <- r2_val
    
    # Run ODE
    out <- ode(y = y0, times = d$t_grid, func = model_obj$rhs, parms = p_sim, method = "lsoda")
    
    # Format
    df_sim <- as.data.frame(out)
    df_sim$well <- w
    
    # Keep Nominal G0 for Labeling consistency
    # This ensures your plots still say "G0=0" or "G0=5" even if the simulation used 0.3 or 4.8
    df_sim$G0_val <- d$G0_per_well[w] 
    df_sim$G1_0_est <- est_g1_0       # Store estimated value for reference/QC
    
    sim_results[[w]] <- df_sim
  }
  
  master_sim <- bind_rows(sim_results) %>%
    mutate(
      line_name = meta$line_name,
      ploidy = meta$ploidy_lbl,
      # Labels derived from NOMINAL value
      G0_lbl = paste0("G0=", signif(G0_val, 3)) 
    )
  
  # --- E. OBSERVATIONS PREP ---
  obsN <- tibble(
    well = d$well_idx_count,
    time = d$t_grid[d$grid_idx_count],
    value = as.numeric(d$N_obs),
    type = "N"
  )
  
  obsG <- tibble(
    well = d$well_idx_gluc,
    time = d$t_grid[d$grid_idx_gluc],
    lum = as.numeric(d$lum_obs),
    dil = as.numeric(d$dilution)
  ) %>%
    mutate(
      exp = d$exp_id[well],
      value = pmax(0, (lum - d$calib_b_fixed[exp]) / (d$calib_a_fixed[exp] * dil + 1e-12)),
      type = "G1"
    ) %>%
    select(well, time, value, type)
  
  master_obs <- bind_rows(obsN, obsG) %>%
    mutate(
      line_name = meta$line_name,
      ploidy = meta$ploidy_lbl,
      G0_val = d$G0_per_well[well],
      G0_lbl = paste0("G0=", signif(G0_val, 3))
    )
  
  return(list(sim = master_sim, obs = master_obs, params = master_params))
}

# ==============================================================================
# 3. PLOTTER (Directly adapted from your reference)
# ==============================================================================


plot_results_grid <- function(obs, sim, target_line, 
                              vars_to_plot = c(N=TRUE, G1=TRUE, G2=TRUE, Lac=TRUE, Amm=TRUE),
                              y_buffer = 1.1) {
  
  # --- 1. Filter Data ---
  sim_sub <- sim %>% filter(line_name == target_line)
  obs_sub <- obs %>% filter(line_name == target_line)
  
  if(nrow(sim_sub) == 0) return(NULL)
  
  # --- 2. Formatting (Sort G0 High -> Low) ---
  g0_levels <- sim_sub %>% 
    select(G0_lbl, G0_val) %>% 
    distinct() %>% 
    arrange(desc(G0_val)) %>% 
    pull(G0_lbl)
  
  sim_sub$G0_lbl <- factor(sim_sub$G0_lbl, levels = g0_levels)
  if(nrow(obs_sub) > 0) obs_sub$G0_lbl <- factor(obs_sub$G0_lbl, levels = g0_levels)
  
  # --- 3. Configuration ---
  obs_mapping <- list(N = "N", G1 = "G1", G2 = NA, Lac = NA, Amm = NA)
  has_models <- "model" %in% colnames(sim_sub)
  
  # --- 4. Panel Constructor ---
  make_var_plot <- function(var_name, include_var) {
    
    if(!include_var) return(NULL)
    
    # Identify Obs Data
    obs_type <- if(var_name %in% names(obs_mapping)) obs_mapping[[var_name]] else NA
    dat_obs <- if(!is.na(obs_type)) filter(obs_sub, type == obs_type) else obs_sub[0,]
    
    p <- ggplot()
    
    # --- LAYER 1: Observations (Grey, Background) ---
    if(nrow(dat_obs) > 0) {
      p <- p + geom_point(data = dat_obs, aes(x = time, y = value), 
                          color = "grey50", size = 1.5, alpha = 0.6)
    }
    
    # --- LAYER 2: Simulation Lines (Colored, Foreground) ---
    if (has_models) {
      p <- p + 
        geom_line(data = sim_sub, 
                  aes(x = time, y = .data[[var_name]],group=well, color = model), 
                  linewidth = 0.8) +
        scale_color_brewer(palette = "Dark2") 
    } else {
      p <- p + 
        geom_line(data = sim_sub, 
                  aes(x = time, y = .data[[var_name]],group=well), 
                  color = "black", linewidth = 0.8)
    }
    
    # --- Formatting ---
    p <- p +
      facet_grid(rows = vars(G0_lbl), cols = vars(ploidy), scales = "free_y") +
      # Global font size 8pt
      theme_bw(base_size = 6) + 
      # Remove axis labels (titles), but keep ticks (numbers)
      labs(title = var_name, y = NULL, x = NULL, color = NULL) + 
      theme(
        # Center title, bold, ensure size 8
        plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
        # Ensure strip text (facet labels) is 8pt
        strip.text = element_text(size = 6),
        strip.background = element_rect(fill = "grey95"),
        # Ensure axis tick text is 8pt
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6)
      )
    
    # --- ggh4x SCALES (Clipping Logic) ---
    if(nrow(dat_obs) > 0) {
      row_limits <- dat_obs %>%
        group_by(G0_lbl) %>%
        summarise(max_val = max(value, na.rm = TRUE) * y_buffer, .groups = "drop")
      
      my_scales <- list()
      for(lvl in levels(sim_sub$G0_lbl)) {
        limit_val <- row_limits$max_val[row_limits$G0_lbl == lvl]
        
        if(length(limit_val) > 0 && !is.na(limit_val)) {
          # Keep oob=function(x) x to allow lines to run off chart
          my_scales[[lvl]] <- scale_y_continuous(
            limits = c(0, limit_val), 
            oob = function(x, ...) x 
          )
        } else {
          my_scales[[lvl]] <- scale_y_continuous() 
        }
      }
      p <- p + facetted_pos_scales(y = my_scales)
    }
    
    return(p)
  }
  
  # --- 5. Assemble ---
  plot_list <- mapply(make_var_plot, names(vars_to_plot), vars_to_plot, SIMPLIFY = FALSE)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if(length(plot_list) == 0) return(NULL)
  
  wrap_plots(plot_list, nrow = 1, guides = "collect") + 
    plot_annotation(
      title = paste("Individual Model Fit:", target_line),
      # Main title slightly larger (10pt), or 8pt if strictly preferred
      theme = theme(plot.title = element_text(size = 6, face = "bold"))
    ) & theme(legend.position = "bottom", legend.box = "horizontal")
}
# ==============================================================================
# HELPER: PARAMETER PLOT FUNCTION
# ==============================================================================
plot_param_comparison <- function(param_df) {
  # Define order: Rates -> Yields -> Saturation -> Tox -> Maintenance
  desired_order <- c(
    "ae1", "ah1", "ae2", "ah2",       
    "Y_11", "Y_21", "Y_12", "Y_22",   
    "k_lac", "k_amm",                 
    "tox_lac", "tox_amm",             
    "m"                               
  )
  desired_order <- unique(param_df$param) 
  # Filter and factorize
  plot_df <- param_df %>%
    filter(param %in% desired_order) %>%
    mutate(param = factor(param, levels = rev(desired_order))) 
  
  ggplot(plot_df, aes(x = value, y = param, color = ploidy)) +
    geom_vline(xintercept = c(1e-4, 1e-2, 1, 100), color = "grey95") +
    geom_point(size = 3, alpha = 0.8) +
    scale_x_log10(labels = scales::label_log()) +
    facet_wrap(~line_name, ncol = 3) +
    scale_color_manual(values = c("hi" = "#E41A1C", "lo" = "#377EB8")) +
    theme_bw() +
    labs(
      title = "Parameter Estimates by Cell Line",
      x = "Parameter Value (Log Scale)", y = NULL, color = "Ploidy"
    ) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# ==============================================================================
# HELPER: PARAMETER RATIO PLOT
# ==============================================================================
plot_param_ratios <- function(param_df) {
  
  # 1. Define Order
  desired_order <- c(
    "ae1", "ah1", "ae2", "ah2",       
    "Y_11", "Y_21", "Y_12", "Y_22",   
    "k_lac", "k_amm",                 
    "tox_lac", "tox_amm",             
    "m"                               
  )
  desired_order <- unique(param_df$param) 
  
  # 2. Reshape & Calculate Ratio
  # We need to pivot so we have 'hi' and 'lo' values in the same row
  plot_df <- param_df %>%
    filter(param %in% desired_order) %>%
    # Pivot wider to get columns 'hi' and 'lo'
    pivot_wider(names_from = ploidy, values_from = value) %>%
    # Filter out lines that don't have both hi and lo fits (incomplete pairs)
    filter(!is.na(hi) & !is.na(lo)) %>%
    mutate(
      ratio = hi / lo,
      log_ratio = log2(ratio),
      param = factor(param, levels = rev(desired_order))
    )
  
  # 3. Plot
  ggplot(plot_df, aes(x = log_ratio, y = param, color = line_name)) +
    # Center line at 0 (Ratio = 1)
    geom_vline(xintercept = 0, size = 0.8, color = "grey40") +
    # Reference lines for 2x and 0.5x
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey80") +
    
    # Points
    geom_point(size = 4, alpha = 0.8) +
    
    # Scales
    # We plot on Log2 scale for symmetry, but label with the actual Ratio
    scale_x_continuous(
      breaks = c(-2, -1, 0, 1, 2, 3), 
      labels = c("0.25x", "0.5x", "1x", "2x", "4x", "8x"),
      limits = c(min(plot_df$log_ratio, -1.5), max(plot_df$log_ratio, 1.5)) # Ensure view is wide enough
    ) +
    
    scale_color_brewer(palette = "Dark2") +
    
    theme_bw() +
    labs(
      title = "Shift in Parameters (High vs Low Ploidy)",
      subtitle = "Points to the right mean High Ploidy has a higher value",
      x = "Ratio (High / Low)",
      y = NULL,
      color = "Cell Line"
    ) +
    theme(
      axis.text.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    )
}

# ==============================================================================
# LOAD FULL POSTERIORS FOR GIVEN MODEL
# ==============================================================================

load_post <- function(MODEL_NAME,CONSTRAINT_NAME){
  # Define your experimental design
  cell_lines <- c("MCF10A", "MDA-MB-231", "SNU668", "SUM-159-chem", "SUM-159-fuse")
  ploidies   <- c("hi", "lo")
  
  cat(">>> Aggregating posteriors across conditions...\n")
  
  all_draws <- list()
  
  for (cl in cell_lines) {
    for (pl in ploidies) {
      
      # Construct path
      fpath <- file.path("ecology/model_selection/data", MODEL_NAME, CONSTRAINT_NAME, cl, pl, "nuts_draws_list.Rds")
      
      if (file.exists(fpath)) {
        cat(sprintf("  Loading: %s %s\n", cl, pl))
        draws_list <- readRDS(fpath)
        
        # Iterate over chains
        for (chain_idx in seq_along(draws_list)) {
          # Convert chain matrix to DF
          d <- as.data.frame(draws_list[[chain_idx]])

          # Keep only the 'theta' columns - and log likelihood
          theta_cols <- c(grep("^theta\\[", names(d), value = TRUE),"log_lik")
          d <- d[, theta_cols, drop = FALSE]
          
          # Add metadata columns
          d$Chain     <- as.factor(chain_idx)
          d$CellLine  <- cl
          d$Ploidy    <- pl
          d$Condition <- paste(cl, pl, sep = "\n") # Formatting for plot axis
          
          all_draws[[length(all_draws) + 1]] <- d
        }
      } else {
        warning(sprintf("  MISSING: %s", fpath))
      }
    }
  }
  
  if (length(all_draws) == 0) stop("No data found.")
  
  
  # Combine into one big tidy DataFrame
  cat(">>> Processing dataframe...\n")
  big_df <- bind_rows(all_draws) %>%
    pivot_longer(
      cols = starts_with("theta"), 
      names_to = "ParamRaw", 
      values_to = "Value"
    ) %>%
    # Apply the readable names
    mutate(ParamName = param_map[ParamRaw]) %>%
    # Order parameters so they appear 1-Npar in the plot, not alphabetically
    mutate(ParamName = factor(ParamName, levels = param_map))
}

# ==============================================================================
# GENERATE A PROBABILITY MAP THAT LO BEATS HI PLOIDY
# ==============================================================================

generate_probability_map <- function(input_df, target_cell_line, model_obj = Model, n_draws = 1000) {
  
  message(sprintf("Processing Cell Line: %s", target_cell_line))
  
  # --- 1. Data Prep Helper ---
  prep_data <- function(p_ploidy) {
    d <- input_df %>%
      filter(CellLine == target_cell_line, Ploidy == p_ploidy)
    
    if(nrow(d) == 0) return(NULL)
    
    # Pivot Wide using ParamName
    # We construct a DrawID based on the order of appearance per parameter
    d_wide <- d %>%
      group_by(ParamName) %>%
      mutate(DrawID = row_number()) %>%
      ungroup() %>%
      select(DrawID, ParamName, Value) %>%
      pivot_wider(names_from = ParamName, values_from = Value) %>%
      select(-DrawID) # Remove ID after pivoting
    
    # Filter to keep only parameters relevant to the model
    # (Prevents errors if input_df has extra tracking vars)
    valid_cols <- intersect(names(d_wide), model_obj$param_names)
    d_wide <- d_wide[, valid_cols, drop = FALSE]
    
    # Subsample if necessary
    if(nrow(d_wide) > n_draws) {
      d_wide <- d_wide[sample(nrow(d_wide), n_draws), ]
    }
    
    return(d_wide)
  }
  
  # Prepare Draws
  draws_lo <- prep_data("lo")
  draws_hi <- prep_data("hi")
  
  if(is.null(draws_lo) || is.null(draws_hi)) {
    warning("Missing data for one or both ploidies.")
    return(NULL)
  }
  
  # Ensure equal lengths for pairwise comparison
  n_min <- min(nrow(draws_lo), nrow(draws_hi))
  draws_lo <- draws_lo[1:n_min, ]
  draws_hi <- draws_hi[1:n_min, ]
  
  # --- 2. Define Grid ---
  # (Using your specific log-scale logic)
  grid_G2 <- c(0, exp(seq(log(0.001), log(1), length.out = 15)))
  grid_G1 <- c(0, exp(seq(log(0.001), log(10), length.out = 20)))
  grid_df <- expand.grid(G1 = grid_G1, G2 = grid_G2)
  
  # --- 3. Compute Probabilities via get_derivs ---
  
  # Helper: Calculates growth rates for a whole data.frame of parameters
  # Returns a numeric vector of rates
  calc_batch_rates <- function(param_df, g1, g2) {
    # We set N=1. Since dN = u * N, if N=1, dN = u (Growth Rate).
    state_in <- c(N = 1, G1 = g1, G2 = g2)
    
    # Apply row-wise (vapply is faster/safer than apply here)
    # We iterate by index to extract the parameter row as a named vector
    vapply(seq_len(nrow(param_df)), function(i) {
      pars <- unlist(param_df[i, ])
      # get_derivs handles defaults (e.g. Lac=0, x1=1) automatically
      derivs <- get_derivs(model_obj, pars, state_in)
      return(derivs[["N"]]) 
    }, FUN.VALUE = numeric(1))
  }
  
  # Loop over grid
  probs <- numeric(nrow(grid_df))
  
  for(i in seq_len(nrow(grid_df))) {
    g1 <- grid_df$G1[i]
    g2 <- grid_df$G2[i]
    
    r_lo <- calc_batch_rates(draws_lo, g1, g2)
    r_hi <- calc_batch_rates(draws_hi, g1, g2)
    
    probs[i] <- mean(r_lo > r_hi)
  }
  
  # --- 4. Final Formatting ---
  grid_df$Prob_Lo_Gt_Hi <- probs
  grid_df$CellLine <- target_cell_line
  
  return(grid_df)
}