library(cmdstanr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
# Try to load gridExtra, handle if missing
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  warning("gridExtra not installed. Diagnostics layout might be simple.")
  has_grid <- FALSE
} else {
  library(gridExtra)
  has_grid <- TRUE
}

# ============================================================
# 1) CONFIG & DATA
# ============================================================
output_dir <- "~/tmp" 

# [HARDCODED SETTINGS] 
WARMUP_ITERS <- 1500
SAVE_WARMUP  <- TRUE 

cat("Loading Data...\n")
stan_data <- readRDS("data/stan_ready_data.Rds")
stan_data$mode <- 0 

if(is.null(stan_data$calib_a_fixed)) stop("stan_data missing 'calib_a_fixed'.")

# Compile model
mod_gen <- cmdstan_model("STAN/model_B.stan", force_recompile = FALSE)

# ============================================================
# 2) PROCESS CHAIN CSVs
# ============================================================
csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("temp_best_chain_", basename(csv_files))]
csv_files <- csv_files[order(file.info(csv_files)$mtime, decreasing = TRUE)]
csv_files <- na.omit(csv_files[1:min(4, length(csv_files))])

if(length(csv_files) == 0) stop("No chain CSV files found in ", output_dir)

cat("Found", length(csv_files), "chain logs. Processing...\n")

# Helper to copy comments
get_csv_comments <- function(csv_file) {
  grep("^#", readLines(csv_file, n = 50), value = TRUE)
}

# Regex for params (Loose match to grab everything first)
param_patterns <- "mu_global|sigma_line|beta_high|z_line|mu_IC|phi_|calib_sigma|lp__"

all_chains_sim <- list()
trace_list     <- list()
diag_list      <- list()
chain_stats    <- rep(NA_character_, length(csv_files))

for(i in seq_along(csv_files)) {
  chain_file <- csv_files[i]
  
  # A. Read Data
  raw_dt <- suppressWarnings(fread(chain_file, skip = "lp__"))
  if(nrow(raw_dt) < 1) { cat("  Chain", i, ": Empty.\n"); next }
  
  # B. Assign Phases
  if (SAVE_WARMUP) {
    raw_dt$iter  <- 1:nrow(raw_dt)
    raw_dt$phase <- ifelse(raw_dt$iter <= WARMUP_ITERS, "Warmup", "Sampling")
  } else {
    raw_dt$iter  <- (WARMUP_ITERS + 1):(WARMUP_ITERS + nrow(raw_dt))
    raw_dt$phase <- "Sampling"
  }
  raw_dt$chain_id <- as.factor(i)
  
  # C. Extract Diagnostics
  diag_cols <- grep("__$", names(raw_dt), value = TRUE)
  if(length(diag_cols) > 0) {
    diag_list[[i]] <- raw_dt[, c("iter", "phase", "chain_id", diag_cols), with = FALSE]
  }
  
  # D. Extract Traces (Parameters)
  # Grab any column that looks like our parameters
  keep_cols <- grep(param_patterns, names(raw_dt), value = TRUE)
  if(length(keep_cols) > 0) {
    trace_list[[i]] <- raw_dt[, c("iter", "phase", "chain_id", keep_cols), with = FALSE]
  }
  
  # E. Find Best Draw & Run GQ
  post_dt <- raw_dt[raw_dt$phase == "Sampling", ]
  if (nrow(post_dt) > 0) {
    status_label <- "Sampling"
    best_idx_local <- which.max(post_dt$lp__)
    best_row <- post_dt[best_idx_local, ]
  } else {
    status_label <- "Warmup"
    best_idx_local <- which.max(raw_dt$lp__)
    best_row <- raw_dt[best_idx_local, ]
  }
  
  best_lp <- best_row$lp__
  chain_stats[i] <- sprintf("Ch%d (%s): lp=%.1f", i, status_label, best_lp)
  
  # Write temp CSV for GQ
  temp_csv <- file.path(output_dir, paste0("temp_best_chain_", i, ".csv"))
  comments <- get_csv_comments(chain_file)
  writeLines(comments, temp_csv)
  
  # Use write.table with quote=FALSE
  cols_to_write <- setdiff(names(best_row), c("iter", "phase", "chain_id"))
  write.table(best_row[, ..cols_to_write], file = temp_csv, sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
  
  gq <- try(mod_gen$generate_quantities(
    fitted_params = temp_csv, data = stan_data, parallel_chains = 1
  ), silent = TRUE)
  
  if(inherits(gq, "try-error")) { 
    cat("  Chain", i, ": GQ execution failed.\n") 
  } else {
    d_sim <- suppressWarnings(as_draws_df(gq$draws()))
    d_sim <- as.data.frame(d_sim)
    
    # Check y_sim exists
    if (!any(grepl("^y_sim", names(d_sim)))) {
      cat("  Chain", i, ": Output missing 'y_sim'. Skipping trajectories.\n")
    } else {
      # Extract y_sim cleanly
      d_sim_long <- d_sim %>%
        select(starts_with("y_sim")) %>%
        pivot_longer(everything(), names_to = "param", values_to = "value") %>%
        mutate(clean = gsub("y_sim\\[|\\]", "", param)) %>%
        separate(clean, c("w", "t", "s"), sep = ",", convert = TRUE) %>%
        mutate(
          well_idx = w, time = stan_data$t_grid[t],
          type = case_when(s==1 ~ "NL", s==2 ~ "ND", s==3 ~ "G"),
          chain_id = as.factor(i)
        ) %>% select(well_idx, time, type, value, chain_id)
      
      all_chains_sim[[i]] <- d_sim_long
    }
  }
}

if(length(diag_list) == 0 && length(trace_list) == 0) stop("No valid data found in CSVs.")

master_diag_df  <- bind_rows(diag_list)
master_trace_df <- bind_rows(trace_list)
master_sim_df   <- if(length(all_chains_sim) > 0) bind_rows(all_chains_sim) %>% pivot_wider(names_from = type, values_from = value) else NULL

# ============================================================
# 3) RECONSTRUCT EFFECTIVE ODE PARAMETERS (ROBUST)
# ============================================================
cat("Reconstructing effective parameters...\n")

softcap <- function(x, cap) cap - log1p(exp(cap - x))
cap_log_main <- 40.0
cap_log_hill <- 2.5 
par10_names  <- c("theta","kp","kd","kd2","g50a","na","g50d","nd","v1","v2")

groups <- tibble(line_id = stan_data$line_id, high = stan_data$is_high_ploidy) %>%
  distinct(line_id, high) %>%
  arrange(line_id, high) %>%
  mutate(group_lbl = paste0("Line ", line_id, " | ", ifelse(high==1, "High", "Low")))

reconstructed_list <- list()
draws_to_use <- master_trace_df %>% filter(phase == "Sampling")
plot_phase_label <- "Sampling Phase"

if(nrow(draws_to_use) == 0) {
  cat("  NOTE: No sampling draws found. Using Warmup draws.\n")
  draws_to_use <- master_trace_df
  plot_phase_label <- "Warmup Phase"
}

# --- ROBUST COLUMN FINDER ---
# Explicitly checks for name.1, name[1], name.1.1, name[1,1]
find_col <- function(df, base, idx1, idx2=NULL) {
  possibles <- if (is.null(idx2)) {
    c(paste0(base, ".", idx1), paste0(base, "[", idx1, "]"))
  } else {
    c(paste0(base, ".", idx1, ".", idx2), paste0(base, "[", idx1, ",", idx2, "]"))
  }
  match <- intersect(possibles, names(df))
  if(length(match) > 0) return(match[1])
  return(NULL)
}

if(nrow(draws_to_use) > 0) {
  for(g in 1:nrow(groups)) {
    l <- groups$line_id[g]
    h <- groups$high[g]
    lbl <- groups$group_lbl[g]
    
    group_df <- data.frame(chain_id = draws_to_use$chain_id)
    
    # Debug print for first group only
    if(g == 1) cat("Debug Group 1 (Line", l, "):\n")
    
    for(pp in 1:10) {
      c_mu   <- find_col(draws_to_use, "mu_global", pp)
      c_sig  <- find_col(draws_to_use, "sigma_line", pp)
      c_z    <- find_col(draws_to_use, "z_line", pp, l)
      c_beta <- find_col(draws_to_use, "beta_high", pp)
      
      # Debug missing cols for first param of first group
      if(g == 1 && pp == 1) {
        if(is.null(c_mu)) cat("  MISSING: mu_global[1]\n")
        if(is.null(c_z))  cat("  MISSING: z_line[1,", l, "]\n")
      }
      
      # Defaults if missing (so we don't crash)
      val_mu   <- if(!is.null(c_mu)) draws_to_use[[c_mu]] else 0
      val_sig  <- if(!is.null(c_sig)) draws_to_use[[c_sig]] else 0
      val_z    <- if(!is.null(c_z)) draws_to_use[[c_z]] else 0
      val_beta <- if(!is.null(c_beta)) draws_to_use[[c_beta]] else 0
      
      # Only calculate if we have at least the global mean (bare minimum)
      if(!is.null(c_mu)) {
        raw <- val_mu + val_sig * val_z + val_beta * h
        val <- if(pp %in% c(6, 8)) 1.0 + exp(softcap(raw, cap_log_hill)) else exp(softcap(raw, cap_log_main))
        group_df[[par10_names[pp]]] <- val
      }
    }
    
    if(ncol(group_df) > 1) {
      reconstructed_list[[g]] <- group_df %>%
        mutate(group = lbl) %>%
        pivot_longer(cols = -c(chain_id, group), names_to = "param", values_to = "value")
    }
  }
}

reconstructed_df <- bind_rows(reconstructed_list)
cat(sprintf("Reconstructed parameters for %d groups.\n", length(unique(reconstructed_df$group))))

# ============================================================
# 4) DATA PREP
# ============================================================
meta_df <- data.frame(well_idx = 1:stan_data$N_wells, line_id = stan_data$line_id, is_high = stan_data$is_high_ploidy, G0 = stan_data$G0_per_well, exp_id = stan_data$exp_id)
meta_df$ploidy_lbl <- ifelse(meta_df$is_high == 1, "High Ploidy", "Low Ploidy")
meta_df$G0_lbl     <- factor(paste0(meta_df$G0, " mM"), levels = paste0(sort(unique(meta_df$G0)), " mM"))

obs_merged <- full_join(
  data.frame(well_idx = stan_data$well_idx_count, time = stan_data$t_grid[stan_data$grid_idx_count], N_obs = stan_data$N_obs, D_obs = stan_data$D_obs),
  data.frame(well_idx = stan_data$well_idx_gluc, time = stan_data$t_grid[stan_data$grid_idx_gluc], lum = stan_data$lum_obs, dilution = stan_data$dilution) %>%
    left_join(meta_df %>% select(well_idx, exp_id), by = "well_idx") %>%
    mutate(a = stan_data$calib_a_fixed[exp_id], b = stan_data$calib_b_fixed[exp_id], G_obs = pmax(0, (lum - b)/(a*dilution))) %>%
    select(well_idx, time, G_obs),
  by = c("well_idx", "time")
) %>% left_join(meta_df, by = "well_idx")

if(!is.null(master_sim_df)) sim_merged <- master_sim_df %>% left_join(meta_df, by = "well_idx")

# ============================================================
# 5) PLOTTING (Multi-Page PDF)
# ============================================================
pdf("chain_comparison_by_line.pdf", width = 14, height = 10)

# PAGE 1: Likelihood
cat("Plotting Likelihood...\n")
print(ggplot(master_trace_df, aes(x = iter, y = lp__, color = chain_id)) +
        geom_line(alpha = 0.8, size = 0.3) + facet_wrap(~phase, scales = "free_x") + theme_bw() +
        labs(title = "Evolution of Log-Posterior (lp__)", y = "Log Probability"))

# PAGE 2: Diagnostics
if(!is.null(master_diag_df) && has_grid) {
  cat("Plotting Diagnostics...\n")
  chain_summ <- master_diag_df %>% group_by(chain_id, phase) %>%
    summarise(div_rate = mean(divergent__ > 0), max_tree_rate = mean(treedepth__ >= 12), .groups="drop")
  
  p1 <- ggplot(chain_summ, aes(x = chain_id, y = div_rate, fill = phase)) + geom_col(position = "dodge") + labs(title="Divergences")
  p2 <- ggplot(chain_summ, aes(x = chain_id, y = max_tree_rate, fill = phase)) + geom_col(position = "dodge") + labs(title="Max Tree Hit")
  p3 <- ggplot(master_diag_df, aes(x = chain_id, y = stepsize__, fill = chain_id)) + geom_violin() + facet_wrap(~phase) + labs(title="Step Size")
  p4 <- ggplot(master_diag_df, aes(x = chain_id, y = energy__, fill = chain_id)) + geom_violin() + facet_wrap(~phase) + labs(title="Energy")
  
  grid.arrange(p1, p2, p3, p4, nrow=2, top = "Sampler Diagnostics")
}

# PAGE 3+: Posteriors
if(nrow(reconstructed_df) > 0) {
  cat("Plotting Distributions...\n")
  for(grp in unique(reconstructed_df$group)) {
    sub_df <- reconstructed_df %>% filter(group == grp)
    print(ggplot(sub_df, aes(x = value, fill = chain_id)) +
            geom_histogram(bins = 50, position = "identity", alpha = 0.5, color=NA) +
            facet_wrap(~param, scales = "free") + scale_x_log10() + scale_fill_brewer(palette = "Set1") +
            theme_bw() + labs(title = paste("Posteriors:", grp), subtitle = plot_phase_label))
  }
}

# PAGE X+: Trajectories
if(!is.null(master_sim_df)) {
  cat("Plotting Trajectories...\n")
  scale_factor <- 1000
  for(lid in sort(unique(meta_df$line_id))) {
    p_sim <- sim_merged %>% filter(line_id == lid)
    p_obs <- obs_merged %>% filter(line_id == lid)
    if(nrow(p_sim) == 0) next
    
    print(ggplot() +
            geom_point(data = p_obs, aes(x = time, y = N_obs), color = "#009E73", size=1.5, alpha=0.4) +
            geom_point(data = p_obs, aes(x = time, y = D_obs), color = "#D55E00", size=1.5, alpha=0.4) +
            geom_point(data = p_obs, aes(x = time, y = G_obs * scale_factor), color = "#E69F00", size=1.5, alpha=0.4, shape=15) +
            geom_line(data = p_sim, aes(x = time, y = NL, group = interaction(well_idx, chain_id), linetype = chain_id), color = "#009E73") +
            geom_line(data = p_sim, aes(x = time, y = ND, group = interaction(well_idx, chain_id), linetype = chain_id), color = "#D55E00") +
            geom_line(data = p_sim, aes(x = time, y = G * scale_factor, group = interaction(well_idx, chain_id), linetype = chain_id), color = "#E69F00") +
            facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
            scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor, name = "Glucose (mM)")) +
            theme_bw() + theme(legend.position="bottom") +
            labs(title = paste("Line", lid, "-", paste(na.omit(chain_stats), collapse="|")), linetype="Chain"))
  }
}
dev.off()
cat("Done. Saved 'chain_comparison_by_line.pdf'.\n")