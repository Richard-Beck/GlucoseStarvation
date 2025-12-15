library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

setwd("/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation")

# ============================================================
# 0) Helpers: fixed calibration (log-fit) matching Stan
# ============================================================
compute_fixed_calibration <- function(stan_data){
  N_exps <- stan_data$N_exps
  cal_df <- data.frame(e = stan_data$calib_exp_idx, G = stan_data$calib_G, Lum = stan_data$calib_Lum) %>%
    filter(is.finite(e), is.finite(G), is.finite(Lum), G >= 0, Lum > 0)
  
  costf <- function(pars, x){
    a <- abs(pars[1]); b <- abs(pars[2])
    mu <- x$G * a + b
    if(any(!is.finite(mu)) || any(mu <= 0)) return(1e12)
    err <- sum((log(x$Lum) - log(mu))^2)
    if(!is.finite(err)) 1e12 else err
  }
  
  a_fix <- rep(NA_real_, N_exps); b_fix <- rep(NA_real_, N_exps); sdlog_fix <- rep(NA_real_, N_exps)
  
  for(e in 1:N_exps){
    sub <- cal_df[cal_df$e == e, , drop = FALSE]
    if(nrow(sub) < 2){
      a_fix[e] <- if(!is.null(stan_data$prior_calib_a_mean)) max(1e-6, stan_data$prior_calib_a_mean) else 1.0
      b_fix[e] <- if(!is.null(stan_data$prior_calib_b_mean)) max(1e-6, stan_data$prior_calib_b_mean) else 1.0
      sdlog_fix[e] <- 0.1
      next
    }
    fit0 <- try(lm(Lum ~ G, data = sub), silent = TRUE)
    a0 <- if(!inherits(fit0,"try-error")) max(1e-6, as.numeric(coef(fit0)["G"])) else 1.0
    b0 <- max(1e-6, min(sub$Lum, na.rm = TRUE))
    opt <- optim(c(a0, b0), costf, x = sub, method = "Nelder-Mead", control = list(maxit = 5000))
    a_hat <- max(1e-6, abs(opt$par[1])); b_hat <- max(1e-6, abs(opt$par[2]))
    mu_hat <- sub$G * a_hat + b_hat
    sd_hat <- sd(log(sub$Lum) - log(mu_hat))
    if(!is.finite(sd_hat) || sd_hat <= 0) sd_hat <- 0.1
    a_fix[e] <- a_hat; b_fix[e] <- b_hat; sdlog_fix[e] <- sd_hat
  }
  list(a_fix = a_fix, b_fix = b_fix, sdlog_fix = sdlog_fix)
}

# Parse cmdstan CSV header for warmup info etc
read_cmdstan_header <- function(csv_file, n = 200){
  hdr <- readLines(csv_file, n = n, warn = FALSE)
  cm  <- grep("^#", hdr, value = TRUE)
  
  grab_int <- function(key){
    x <- grep(paste0("^#\\s*", key, "\\s*="), cm, value = TRUE)
    if(length(x) == 0) return(NA_integer_)
    y <- sub(".*=\\s*", "", x[1])
    suppressWarnings(as.integer(y))
  }
  grab_lgl <- function(key){
    x <- grep(paste0("^#\\s*", key, "\\s*="), cm, value = TRUE)
    if(length(x) == 0) return(NA)
    y <- tolower(sub(".*=\\s*", "", x[1]))
    y %in% c("1","true","yes")
  }
  
  list(
    num_warmup    = grab_int("num_warmup"),
    save_warmup   = grab_lgl("save_warmup"),
    comment_lines = cm
  )
}

# ============================================================
# 1) Setup
# ============================================================
stan_data <- readRDS("data/stan_ready_data.Rds")

if(!all(c("calib_a_fixed","calib_b_fixed") %in% names(stan_data))){
  cat("calib_a_fixed / calib_b_fixed missing from stan_data; computing now...\n")
  cf <- compute_fixed_calibration(stan_data)
  stan_data$calib_a_fixed <- cf$a_fix
  stan_data$calib_b_fixed <- cf$b_fix
} else {
  cat("Found calib_a_fixed / calib_b_fixed in stan_data.\n")
}

stan_data$prior_only <- 0

output_dir <- path.expand("~/tmp")
mod_gen <- cmdstan_model("STAN/model_B.stan")

# Find chain CSVs, excluding our own temp files
csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("temp_best_chain_", basename(csv_files))]
csv_files <- csv_files[order(file.info(csv_files)$mtime, decreasing = TRUE)]
csv_files <- na.omit(csv_files[1:min(4, length(csv_files))])
if(length(csv_files) == 0) stop("No chain CSV files found in output_dir!")

cat("Found", length(csv_files), "chains.\n")
chain_stats <- rep(NA_character_, length(csv_files))
best_param_csvs <- rep(NA_character_, length(csv_files))

# ============================================================
# 2) Harvest best *post-warmup* draw + diagnostics
# ============================================================
all_chains_sim <- list()
diag_rows      <- list()

for(i in seq_along(csv_files)) {
  chain_file <- csv_files[i]
  cat("\n--- Reading:", chain_file, "\n")
  
  hinfo <- read_cmdstan_header(chain_file)
  
  # read draws table (skipping comment lines automatically via skip="lp__" anchor)
  raw_dt <- suppressWarnings(fread(chain_file, skip = "lp__"))
  if(nrow(raw_dt) < 1) { cat("No draws yet.\n"); next }
  
  raw_dt$iter <- seq_len(nrow(raw_dt))
  raw_dt$chain_id <- as.factor(i)
  
  # diagnostics table
  diag_cols <- intersect(c("lp__","accept_stat__","stepsize__","treedepth__","n_leapfrog__","divergent__","energy__"), names(raw_dt))
  if(length(diag_cols) > 0) diag_rows[[i]] <- raw_dt[, c(diag_cols, "iter", "chain_id"), with = FALSE]
  
  warmup_kept <- 0L
  if (isTRUE(hinfo$save_warmup) && !is.na(hinfo$num_warmup) && hinfo$num_warmup > 0) {
    warmup_kept <- min(nrow(raw_dt), hinfo$num_warmup)  # thin ignored, but never crashes
  }
  
  post_dt <- if (warmup_kept > 0L && nrow(raw_dt) > warmup_kept) raw_dt[(warmup_kept+1L):nrow(raw_dt)] else raw_dt
  if(nrow(post_dt) < 1) {
    cat("No post-warmup draws yet (still in warmup, or save_warmup=FALSE).\n")
    post_dt <- raw_dt
    chain_stats[i] <- paste0("Ch", i, ": best(lp) from ALL iters (warmup-only so far)")
  }
  
  if(!("lp__" %in% names(post_dt))) stop("lp__ not found in CSV (unexpected).")
  
  best_idx <- which.max(post_dt$lp__)
  best_lp  <- post_dt$lp__[best_idx]
  cat(sprintf("Chain %d: Best lp__ used for GQ = %.2f\n", i, best_lp))
  
  if(is.na(chain_stats[i])) chain_stats[i] <- paste0("Ch", i, ": best lp=", round(best_lp, 0))
  
  # write temp "best draw" CSV with original comment header
  temp_csv <- file.path(output_dir, paste0("temp_best_chain_", i, ".csv"))
  writeLines(hinfo$comment_lines, temp_csv)
  fwrite(post_dt[best_idx, setdiff(names(post_dt), c("iter","chain_id")), with = FALSE],
         file = temp_csv, append = TRUE, col.names = TRUE)
  
  best_param_csvs[i] <- temp_csv
  
  # Generate quantities (trap errors)
  gq <- try(mod_gen$generate_quantities(fitted_params = temp_csv, data = stan_data, parallel_chains = 1), silent = TRUE)
  if(inherits(gq, "try-error")) {
    cat("GQ failed for chain", i, ":\n", as.character(gq), "\n")
    next
  }
  draws <- gq$draws()
  
  # extract y_sim for all wells; IMPORTANT: group by (well,chain) later
  chain_list <- vector("list", stan_data$N_wells)
  for(w in 1:stan_data$N_wells) {
    extract_state <- function(s_idx) {
      vars <- paste0("y_sim[", w, ",", 1:stan_data$N_grid, ",", s_idx, "]")
      as.numeric(subset_draws(draws, variable = vars))
    }
    chain_list[[w]] <- data.frame(
      time     = stan_data$t_grid,
      NL       = extract_state(1),
      ND       = extract_state(2),
      G        = extract_state(3),
      well_idx = w,
      chain_id = as.factor(i)
    )
  }
  all_chains_sim[[i]] <- bind_rows(chain_list)
}

master_sim_df <- bind_rows(all_chains_sim)
diag_df <- bind_rows(diag_rows)

if(nrow(master_sim_df) == 0) stop("No successful GQ draws produced yet (likely still too early or GQ failing).")

# ============================================================
# 3) Enrich with metadata + observed data
# ============================================================
meta_df <- data.frame(
  well_idx = 1:stan_data$N_wells,
  line_id  = stan_data$line_id,
  is_high  = stan_data$is_high_ploidy,
  G0       = stan_data$G0_per_well,
  exp_id   = stan_data$exp_id
)

meta_df$ploidy_lbl <- ifelse(meta_df$is_high == 1, "High Ploidy", "Low Ploidy")
meta_df$G0_lbl <- factor(paste0(meta_df$G0, " mM"), levels = paste0(sort(unique(meta_df$G0)), " mM"))

master_sim_df <- master_sim_df %>% left_join(meta_df, by = "well_idx")

obs_df <- data.frame(
  time     = stan_data$t_grid[stan_data$grid_idx_count],
  N_obs    = stan_data$N_obs,
  D_obs    = stan_data$D_obs,
  well_idx = stan_data$well_idx_count
) %>% left_join(meta_df, by = "well_idx")

# ============================================================
# 4) Grab a few parameters (best-draw per chain) for quick comparisons
# ============================================================
extract_param_from_best <- function(temp_csv, pars){
  dt <- suppressWarnings(fread(temp_csv, skip = "lp__"))
  if(nrow(dt) != 1) return(NULL)
  keep <- intersect(pars, names(dt))
  if(length(keep) == 0) return(NULL)
  out <- as.data.frame(dt[, ..keep])
  out$chain_id <- factor(gsub(".*temp_best_chain_(\\d+)\\.csv$", "\\1", temp_csv))
  out
}

pars_to_track <- c("mu_IC[1]","mu_IC[2]","mu_global[7]","mu_global[5]","mu_global[10]","phi_N","phi_D",
                   "calib_sigma[1]","calib_sigma[2]","calib_sigma[3]")

par_df <- bind_rows(lapply(na.omit(best_param_csvs), extract_param_from_best, pars = pars_to_track))

# ============================================================
# 5) Plotting PDF report (correct grouping!)
# ============================================================
pdf("chain_comparison_by_line.pdf", width = 10, height = 18)
scale_factor <- 1000
unique_lines <- sort(unique(meta_df$line_id))

cat("\nGenerating PDF report...\n")

# ---- Page 0: sampler diagnostics ----
if(nrow(diag_df) > 0){
  p_lp <- ggplot(diag_df, aes(x = iter, y = lp__, color = chain_id)) +
    geom_line(alpha = 0.7, linewidth = 0.4) +
    theme_bw() + labs(title = "Sampler diagnostic: lp__ trace", x = "Iteration (saved draws)", y = "lp__")
  print(p_lp)
  
  if("divergent__" %in% names(diag_df)){
    div_summ <- diag_df %>%
      group_by(chain_id) %>%
      summarise(div_rate = mean(divergent__ > 0), .groups = "drop")
    p_div <- ggplot(div_summ, aes(x = chain_id, y = div_rate)) +
      geom_col() + theme_bw() + labs(title = "Divergence rate by chain (saved draws)", x = "Chain", y = "Fraction divergent")
    print(p_div)
  } else {
    print(ggplot() + theme_void() + labs(title = "No divergent__ column found in CSV"))
  }
  
  if("treedepth__" %in% names(diag_df)){
    p_td <- ggplot(diag_df, aes(x = treedepth__, fill = chain_id)) +
      geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
      theme_bw() + labs(title = "Treedepth distribution by chain (saved draws)", x = "treedepth__", y = "count")
    print(p_td)
  } else {
    print(ggplot() + theme_void() + labs(title = "No treedepth__ column found in CSV"))
  }
}

# ---- Page 1: parameter comparisons ----
if(nrow(par_df) > 0){
  par_long <- par_df %>% pivot_longer(cols = setdiff(names(par_df), "chain_id"), names_to = "param", values_to = "value")
  p_params <- ggplot(par_long, aes(x = chain_id, y = value)) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~ param, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(title = "Best draw per chain (POST-warmup if available): selected parameters", x = "Chain", y = "Value")
  print(p_params)
}

# ---- Main pages: trajectory overlays by cell line ----
for(lid in unique_lines) {
  page_sim <- master_sim_df %>% filter(line_id == lid)
  page_obs <- obs_df %>% filter(line_id == lid)
  
  page_title <- paste0("Cell Line ", lid, " - Chain Comparison\n", paste(na.omit(chain_stats), collapse = " | "))
  
  p <- ggplot() +
    geom_point(data = page_obs, aes(x = time, y = N_obs), color = "#009E73", size = 1.5, alpha = 0.35, shape = 16) +
    geom_point(data = page_obs, aes(x = time, y = D_obs), color = "#D55E00", size = 1.5, alpha = 0.35, shape = 17) +
    
    geom_line(data = page_sim, aes(x = time, y = NL, group = interaction(well_idx, chain_id), linetype = chain_id),
              color = "#009E73", linewidth = 0.7, alpha = 0.75) +
    geom_line(data = page_sim, aes(x = time, y = ND, group = interaction(well_idx, chain_id), linetype = chain_id),
              color = "#D55E00", linewidth = 0.7, alpha = 0.75) +
    geom_line(data = page_sim, aes(x = time, y = G * scale_factor, group = interaction(well_idx, chain_id), linetype = chain_id),
              color = "#E69F00", linewidth = 0.5, alpha = 0.55) +
    
    facet_grid(G0_lbl ~ ploidy_lbl, scales = "free") +
    scale_y_continuous(name = "Cell Count", sec.axis = sec_axis(~ . / scale_factor, name = "Glucose (mM)")) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold", size = 10),
          axis.title.y.right = element_text(color = "#E69F00"),
          legend.position = "bottom") +
    labs(title = page_title,
         subtitle = "Green=Live, Red=Dead, Orange=Glucose (scaled) | Linetype = Chain ID | Each well plotted separately",
         linetype = "Chain ID")
  
  print(p)
  cat("  Page for Cell Line", lid, "generated.\n")
}

dev.off()
cat("Done. Saved to 'chain_comparison_by_line.pdf'\n")
