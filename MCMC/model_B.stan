functions {
  real softcap(real x, real cap) {
    return cap - log1p_exp(cap - x);
  }

  // The ODE function stays "Physical" (receives v1, v2, g50d)
  // We calculate these from Yield/Gamma/m before calling this.
  vector model_b_ode(real t, vector y, array[] real p) {
    real theta = p[1];
    real kp    = p[2];
    real kd    = p[3];
    real kd2   = 0;     // p[4] unused/fixed to 0
    real g50a  = p[5];
    real na    = p[6];
    real g50d  = p[7];  // Physical g50d
    real nd    = p[8];
    real v1    = p[9];  // Physical v1
    real v2    = p[10]; // Physical v2
    
    // [CHANGE 1] Read states directly (Real Space)
    real NL = y[1]; 
    real ND = y[2];
    
    real k_smooth = 100.0;
    real G = log1p_exp(k_smooth * y[3]) / k_smooth;
    real log_G = log(G + 1e-18);
    
    real actA   = inv_logit(na * (log_G - log(g50a)));
    real term_d = inv_logit(nd * (log_G - log(g50d)));
    real inhD   = 1.0 - term_d;
    
    // [CHANGE 2] Derivatives in Real Space
    // du_dt = N * (Growth_Rate - Death_Rate)
    real specific_growth = kp * (1.0 - NL/theta) * actA - kd * inhD - kd2 * NL/theta;
    real du_dt = specific_growth * NL; 
    
    // dv_dt = Influx of Dead Cells (No division by ND)
    real dv_dt = kd * NL * inhD + kd2 * (NL * NL) / theta;
    
    real dG_dt = -NL * (v1 * actA + v2 * term_d) / 2.0;
    
    vector[3] dydt;
    dydt[1] = du_dt;
    dydt[2] = dv_dt;
    dydt[3] = dG_dt;

    return dydt;
  }

  real partial_sum_lpmf(
      array[] int slice_wells,
      int start, int end,
      array[] int line_id, vector ploidy_metric, array[] int exp_id, vector G0_per_well, array[] real t_grid,
      array[] int w_idx_count, array[] int g_idx_count, array[] int N_obs, array[] int D_obs,
      array[] int w_idx_gluc, array[] int g_idx_gluc, array[] real lum_obs, array[] real dilution,
      array[] int is_censored,
      vector calib_a_fixed, vector calib_b_fixed,
      array[] real log_lower, array[] real log_upper,
      vector mu_global, vector sigma_line, 
      vector beta_high, matrix z_line,
      vector sigma_beta, matrix z_beta,
      vector mu_IC, vector sigma_IC, vector beta_IC, matrix z_IC,
      vector sigma_beta_IC, matrix z_beta_IC,
      vector calib_sigma_fixed,
      real phi_total, real phi_frac,
      array[] int has_starvation
  ) {
    real log_lik = 0;
    int N_grid = size(t_grid);
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;

    for (i in 1:size(slice_wells)) {
      int w = slice_wells[i];
      int l = line_id[w];
      real p_met = ploidy_metric[w];

      // 1. Construct Meta-Parameters from Hierarchical Model
      //    Indices: 
      //    7 = Gamma (Gap Ratio)
      //    9 = Yield 
      //    10 = m (Idle Burn)
      array[10] real p_meta;
      for (pp in 1:10) {
        real beta_eff = beta_high[pp] + sigma_beta[pp] * z_beta[pp, l];
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_eff * p_met;

        // raw -> log(parameter) constrained to [log_lower, log_upper]
        real logp = log_lower[pp] + (log_upper[pp] - log_lower[pp]) * inv_logit(raw);

        // back to natural units
        p_meta[pp] = exp(logp);
      }


      // 2. Convert to Physical Parameters for ODE
      //    p_ode must match: [theta, kp, kd, kd2, g50a, na, g50d, nd, v1, v2]
      array[10] real p_ode;
      
      // Direct pass-through
      p_ode[1] = p_meta[1]; // theta
      p_ode[2] = p_meta[2]; // kp
      p_ode[3] = p_meta[3]; // kd
      p_ode[4] = p_meta[4]; // kd2
      p_ode[5] = p_meta[5]; // g50a
      p_ode[6] = p_meta[6]; // na
      p_ode[8] = p_meta[8]; // nd
      
      // Derived / Coupled Parameters
      // g50d = g50a * gamma (p_meta[7])
      p_ode[7] = p_ode[5] * p_meta[7]; 
      
      // v1 = kp / yield (p_meta[9])
      p_ode[9] = p_ode[2] / p_meta[9];
      
      // v2 = m (p_meta[10])
      p_ode[10] = p_meta[10];


      // 3. Initial Conditions
      real beta_eff_IC1 = beta_IC[1] + sigma_beta_IC[1] * z_beta_IC[1, l];
      real beta_eff_IC2 = beta_IC[2] + sigma_beta_IC[2] * z_beta_IC[2, l];

      vector[3] y0_inferred;
      y0_inferred[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_eff_IC1 * p_met;
      y0_inferred[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_eff_IC2 * p_met;
      y0_inferred[3] = 0.0; 

      vector[3] y_start_main;

      // 4. Starvation Protocol
      // [FIX] Clamp to log(theta) + buffer.
      // p_ode[1] is theta. We allow IC to be at most ~2.7x theta (log space +1).
      // This prevents 1,000,000x overshoot which creates infinite gradients.
      real log_theta_cap = log(p_ode[1]) + 1.0;
      real safe_IC_N = fmin(y0_inferred[1], log_theta_cap);
      real safe_IC_D = fmin(y0_inferred[2], log_theta_cap);

      if (has_starvation[w] == 1) {
        vector[3] y0_starve;
        y0_starve[1] = exp(safe_IC_N); 
        y0_starve[2] = exp(safe_IC_D); 
        y0_starve[3] = 0.0; 

        array[1] real t_starve = {0.0};
        array[1] vector[3] y_res_starve;
        y_res_starve = ode_bdf_tol(model_b_ode, y0_starve, -6.0, t_starve, 1e-4, 1e-5, 5000, p_ode);

        y_start_main = y_res_starve[1];
        y_start_main[3] = G0_per_well[w];
      } else {
        y_start_main[1] = exp(safe_IC_N);
        y_start_main[2] = exp(safe_IC_D);
        y_start_main[3] = G0_per_well[w];
      }

      // 5. Main Experiment
      array[N_grid] vector[3] y_hat;
      
      if (is_nan(y_start_main[1]) || is_inf(y_start_main[1]) ||
          is_nan(y_start_main[2]) || is_inf(y_start_main[2]) ||
          is_nan(y_start_main[3]) || is_inf(y_start_main[3])) {
        // Fill with zeros if ICs are not finite
        for (g in 1:N_grid) y_hat[g] = rep_vector(1e-6, 3);
      } else {
        y_hat = ode_bdf_tol(model_b_ode, y_start_main, 0.0, t_eval, 1e-4, 1e-5, 5000, p_ode);
      }


      for (n in 1:size(w_idx_count)) {
        if (w_idx_count[n] == w) {
          int idx = g_idx_count[n];
          real NL_hat = fmax(y_hat[idx, 1], 1e-12); 
          real ND_hat = fmax(y_hat[idx, 2], 1e-12);
          real total_hat = NL_hat + ND_hat;
          real p_hat     = NL_hat / (total_hat + 1e-18);
          p_hat = fmin(fmax(p_hat, 1e-6), 1.0 - 1e-6);
          int total_obs = N_obs[n] + D_obs[n];

          log_lik += neg_binomial_2_lpmf(total_obs | total_hat, phi_total);
          
          real alpha_p = p_hat * phi_frac;
          real beta_p  = (1.0 - p_hat) * phi_frac;
          log_lik += beta_binomial_lpmf(N_obs[n] | total_obs, alpha_p, beta_p);
        }
      }

      for (n in 1:size(w_idx_gluc)) {
        if (w_idx_gluc[n] == w) {
          int idx = g_idx_gluc[n];
          real k_smooth = 100.0;
          real G_hat = log1p_exp(k_smooth * y_hat[idx, 3]) / k_smooth;
          int e = exp_id[w];
          real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e];
          
          if (is_censored[n] == 1) {
             log_lik += lognormal_lcdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
          } else {
             log_lik += lognormal_lpdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
          }
        }
      }
    }
    return log_lik;
  }
}

data {
  
  int<lower=1> N_params;
  array[N_params] real<lower=0> lower_b;
  array[N_params] real<lower=0> upper_b;
  
  int<lower=1> N_wells;
  int<lower=1> N_lines;
  int<lower=1> N_exps;
  int<lower=1> N_obs_count;
  int<lower=1> N_obs_gluc;
  int<lower=1> N_obs_calib;

  int<lower=1> N_grid;
  array[N_grid] real t_grid;
  array[N_wells] int line_id;
  
  vector[N_wells] ploidy_metric;
  array[N_wells] int has_starvation;
  
  array[N_wells] int exp_id;

  array[N_obs_count] int well_idx_count;
  array[N_obs_count] int grid_idx_count;
  array[N_obs_count] int N_obs;
  array[N_obs_count] int D_obs;

  array[N_obs_gluc] int well_idx_gluc;
  array[N_obs_gluc] int grid_idx_gluc;
  array[N_obs_gluc] real lum_obs;
  array[N_obs_gluc] real dilution;
  array[N_obs_gluc] int<lower=0, upper=1> is_censored;

  array[N_obs_calib] int calib_exp_idx;
  array[N_obs_calib] real calib_G;
  array[N_obs_calib] real calib_Lum;

  vector[N_wells] G0_per_well;
  int<lower=1> grainsize;

  real prior_mu_N0_mean;
  real prior_mu_N0_sd;
  real prior_mu_D0_mean;
  real prior_mu_D0_sd;

  vector[10] prior_ode_mean;
  vector[10] prior_ode_sd;

  vector<lower=0>[N_exps] calib_a_fixed;
  vector<lower=0>[N_exps] calib_b_fixed;
  vector<lower=0>[N_exps] calib_sigma_fixed;

  int<lower=0,upper=2> mode;
  int<lower=0,upper=1> calc_sim;
}

transformed data {
  array[N_params] real log_lower;
  array[N_params] real log_upper;
  for (pp in 1:N_params) {
    log_lower[pp] = log(lower_b[pp]);
    log_upper[pp] = log(upper_b[pp]);
  }
}

parameters {
  // Same 10 parameters, but interpretation changes:
  // 7=Gamma, 9=Yield, 10=m
  vector[10] mu_global_raw;
  vector<lower=0>[10] sigma_line;
  vector[10] beta_high;
  matrix[10, N_lines] z_line;

  vector<lower=0>[10] sigma_beta;
  matrix[10, N_lines] z_beta;

  vector[2] mu_IC;
  vector<lower=0>[2] sigma_IC;
  vector[2] beta_IC;
  matrix[2, N_lines] z_IC;

  vector<lower=0>[2] sigma_beta_IC;
  matrix[2, N_lines] z_beta_IC;

  real<lower=1e-4> phi_total;
  real<lower=1e-4> phi_frac;
}

transformed parameters {
  vector[N_params] mu_global;
  for (pp in 1:N_params) {
    mu_global[pp] =
      log_lower[pp] + (log_upper[pp] - log_lower[pp]) * inv_logit(mu_global_raw[pp]);
  }
}

model {
  mu_global_raw ~ logistic(0, 1);
  sigma_line ~ exponential(1);
  to_vector(z_line) ~ std_normal();
  beta_high ~ normal(0, 1);

  sigma_beta ~ exponential(1);
  to_vector(z_beta) ~ std_normal();

  mu_IC[1] ~ normal(prior_mu_N0_mean, prior_mu_N0_sd);
  mu_IC[2] ~ normal(prior_mu_D0_mean, prior_mu_D0_sd);
  sigma_IC ~ normal(0, 0.5);
  to_vector(z_IC) ~ std_normal();
  beta_IC ~ normal(0, 1);

  sigma_beta_IC ~ exponential(1);
  to_vector(z_beta_IC) ~ std_normal();

  phi_total ~ exponential(0.1); 
  phi_frac ~ exponential(0.1); 
  
  if (mode == 0) {
    array[N_wells] int seq_wells;
    for (i in 1:N_wells) seq_wells[i] = i;
    target += reduce_sum(
      partial_sum_lpmf,
      seq_wells,
      grainsize,
      line_id, ploidy_metric, exp_id, G0_per_well, t_grid,
      well_idx_count, grid_idx_count, N_obs, D_obs,
      well_idx_gluc, grid_idx_gluc, lum_obs, dilution,
      is_censored,
      calib_a_fixed, calib_b_fixed,
      log_lower, log_upper,
      mu_global, sigma_line, beta_high, z_line,
      sigma_beta, z_beta,
      mu_IC, sigma_IC, beta_IC, z_IC,
      sigma_beta_IC, z_beta_IC,
      calib_sigma_fixed,
      phi_total, phi_frac,
      has_starvation
    );
  }
}

generated quantities {
  array[N_wells, N_grid] vector[3] y_sim;
  real log_lik = 0;

  // Only calculate if not in prior-predictive mode (mode!=1) 
  // and if simulation flag is on (calc_sim==1)
  if (mode != 1 && calc_sim == 1) {
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;
    
    real cap_log_main = 40.0;
    real cap_log_hill = 1.6;
    
    // ------------------------------------------------------------
    // 1. Simulate Trajectories (y_sim) for ALL Wells
    // ------------------------------------------------------------
    for (w in 1:N_wells) {
      int l = line_id[w];
      real p_met = ploidy_metric[w];
      
      // A. Construct Meta Parameters (Yield, m, Gamma)
      array[10] real p_meta;
      for (pp in 1:10) {
        real beta_eff = beta_high[pp] + sigma_beta[pp] * z_beta[pp, l];
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_eff * p_met;
      
        // Constrain raw -> log-parameter in [log_lower, log_upper]
        real logp = log_lower[pp] + (log_upper[pp] - log_lower[pp]) * inv_logit(raw);
      
        // Convert back to natural units
        p_meta[pp] = exp(logp);
      }

      // B. Convert to Physical Parameters for ODE
      // p_ode: [theta, kp, kd, kd2, g50a, na, g50d, nd, v1, v2]
      array[10] real p_ode;
      p_ode[1] = p_meta[1];
      p_ode[2] = p_meta[2];
      p_ode[3] = p_meta[3];
      p_ode[4] = p_meta[4];
      p_ode[5] = p_meta[5];
      p_ode[6] = p_meta[6];
      p_ode[8] = p_meta[8];
      
      p_ode[7]  = p_ode[5] * p_meta[7];    // g50d = g50a * gamma
      p_ode[9]  = p_ode[2] / p_meta[9];    // v1 = kp / Yield
      p_ode[10] = p_meta[10];              // v2 = m

      // C. Initial Conditions
      real beta_eff_IC1 = beta_IC[1] + sigma_beta_IC[1] * z_beta_IC[1, l];
      real beta_eff_IC2 = beta_IC[2] + sigma_beta_IC[2] * z_beta_IC[2, l];

      vector[3] y0_inferred;
      y0_inferred[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_eff_IC1 * p_met;
      y0_inferred[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_eff_IC2 * p_met;
      y0_inferred[3] = 0.0;
      
      vector[3] y_start_main;
      
      if (mode == 2) {
        print(
          "DEBUG_PRIOR ",
          "w=", w,
          " line=", line_id[w],
          " ploidy=", ploidy_metric[w],
          " starve=", has_starvation[w],
          " G0=", G0_per_well[w],
          " | y0=[", y0_inferred[1], ",", y0_inferred[2], ",", y0_inferred[3], "]",
          " | p_ode=[",
            p_ode[1], ",", p_ode[2], ",", p_ode[3], ",", p_ode[4], ",",
            p_ode[5], ",", p_ode[6], ",", p_ode[7], ",", p_ode[8], ",",
            p_ode[9], ",", p_ode[10],
          "]"
        );
      }


      // D. Starvation Protocol
      // [FIX] Apply same clamping to simulation
      real log_theta_cap = log(p_ode[1]) + 1.0;
      real safe_IC_N = fmin(y0_inferred[1], log_theta_cap);
      real safe_IC_D = fmin(y0_inferred[2], log_theta_cap);

      if (has_starvation[w] == 1) {
        vector[3] y0_starve;
        y0_starve[1] = exp(safe_IC_N);
        y0_starve[2] = exp(safe_IC_D);
        y0_starve[3] = 0.0;

        array[1] real t_starve = {0.0};
        array[1] vector[3] y_res_starve;
        y_res_starve = ode_bdf_tol(model_b_ode, y0_starve, -6.0, t_starve, 1e-4, 1e-5, 50000, p_ode);

        y_start_main = y_res_starve[1];
        y_start_main[3] = G0_per_well[w];
      } else {
        y_start_main[1] = exp(safe_IC_N);
        y_start_main[2] = exp(safe_IC_D);
        y_start_main[3] = G0_per_well[w];
      }
      
      // E. Solve Main ODE
      array[N_grid] vector[3] y_hat;

      if (is_nan(y_start_main[1]) || is_inf(y_start_main[1]) ||
          is_nan(y_start_main[2]) || is_inf(y_start_main[2]) ||
          is_nan(y_start_main[3]) || is_inf(y_start_main[3])) {
        for (g in 1:N_grid) y_hat[g] = rep_vector(1e-6, 3);
      } else {
        y_hat = ode_bdf_tol(model_b_ode, y_start_main, 0.0, t_eval, 1e-4, 1e-5, 50000, p_ode);
      }
        
      // F. Store Results
      for (g in 1:N_grid) {
        // [CHANGE 8] No exp(), just read straight values
        real NL_hat = fmax(y_hat[g, 1], 1e-12);
        real ND_hat = fmax(y_hat[g, 2], 1e-12);
        
        real k_smooth = 100.0; 
        real G_hat = log1p_exp(k_smooth * y_hat[g, 3]) / k_smooth;
        
        y_sim[w, g, 1] = NL_hat;
        y_sim[w, g, 2] = ND_hat;
        y_sim[w, g, 3] = G_hat;
      }
    }
    
    // ------------------------------------------------------------
    // 2. Calculate Likelihoods (Loop over Observations)
    // ------------------------------------------------------------
    if (mode == 0) {
       
       // A. Count Data Likelihood
       for (n in 1:N_obs_count) {
          int w = well_idx_count[n];
          int g = grid_idx_count[n];
          
          real NL_hat = y_sim[w, g, 1];
          real ND_hat = y_sim[w, g, 2];
          
          real total_hat = NL_hat + ND_hat;
          real p_hat     = NL_hat / (total_hat + 1e-18);
          p_hat = fmin(fmax(p_hat, 1e-6), 1.0 - 1e-6); // Clamp prob
          
          int total_obs = N_obs[n] + D_obs[n];
          
          // Total Count (NegBinomial)
          log_lik += neg_binomial_2_lpmf(total_obs | total_hat, phi_total);
          
          // Fraction Live (BetaBinomial)
          real alpha_p = p_hat * phi_frac;
          real beta_p  = (1.0 - p_hat) * phi_frac;
          log_lik += beta_binomial_lpmf(N_obs[n] | total_obs, alpha_p, beta_p);
       }

       // B. Glucose Data Likelihood
       for (n in 1:N_obs_gluc) {
          int w = well_idx_gluc[n];
          int g = grid_idx_gluc[n];
          int e = exp_id[w];
          
          real G_hat = y_sim[w, g, 3];
          real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e];
          
          // Censored (Left-Censored / LO) vs Observed
          if (is_censored[n] == 1) {
             // Observation is "somewhere below this value" -> CDF
             log_lik += lognormal_lcdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
          } else {
             // Observation is exactly this value -> PDF
             log_lik += lognormal_lpdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
          }
       }
    }
  }
}
