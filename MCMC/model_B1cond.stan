functions {
  // --- UTILITY FUNCTIONS ---
  real softcap(real x, real cap) {
    return cap - log1p_exp(cap - x);
  }

  real soft_lower(real x, real limit) {
    return limit + log1p_exp(x - limit);
  }

  // --- ODE SYSTEM (Unchanged) ---
  vector model_b_ode(real t, vector y, array[] real p) {
    real theta = p[1];
    real kp    = p[2];
    real kd    = p[3];
    real kd2   = p[4];
    real g50a  = p[5];
    real na    = p[6];
    real g50d  = p[7];
    real nd    = p[8];
    real v1    = p[9];
    real v2    = p[10];
    
    // control glucose uptake independent of birth/death
    real g50g = p[11]; 
    real ng = p[12];
    
    real NL = y[1];
    real ND = y[2];
    
    real k_smooth = 100.0;
    real G = log1p_exp(k_smooth * y[3]) / k_smooth;
    real log_G = log(G + 1e-18);
    
    real actA   = inv_logit(na * (log_G - log(g50a)));
    real term_d = inv_logit(nd * (log_G - log(g50d)));
    real term_g = inv_logit(ng * (log_G - log(g50g))); // invented term. I could also have subbed nd... this is hacky!
    real inhD   = 1.0 - term_d;
    
    real specific_growth = kp * (1.0 - NL/theta) * actA - kd * inhD - kd2 * NL/theta;
    real du_dt = specific_growth * NL; 
    real dv_dt = kd * NL * inhD + kd2 * (NL * NL) / theta;
    // real dG_dt = -NL * (v1 * actA + v2 * term_d) / 2.0; prior version
    real dG_dt = -NL * (v1 * actA + v2 * term_g) / 2.0;
    
    return [du_dt, dv_dt, dG_dt]';
  }

  // --- PARAMETER MAPPING (Simplified) ---
  // Converts the raw constrained parameters to the physical ODE scales
  // Handles the logic for p[7] and p[9] dependencies
  array[] real get_p_ode(vector p_phys) {
    array[12] real p;
    p[1] = p_phys[1]; // theta
    p[2] = p_phys[2]; // kp
    p[3] = p_phys[3]; // kd
    p[4] = p_phys[4]; // kd2
    p[5] = p_phys[5]; // g50a
    p[6] = p_phys[6]; // na
    p[8] = p_phys[8]; // nd
    p[7]  = p_phys[5] * p_phys[7]; // p[7] is a multiplier in inputs, converted to abs here
    p[9]  = p_phys[2] / p_phys[9]; // p[9] is yield inverse
    p[10] = p_phys[10];       
    p[11] = p_phys[7] * p_phys[11];
    p[12] = p_phys[12];
    return p;
  }

  // --- SOLVER CORE ---
  array[] vector solve_well_trajectory(
    real G0, int starve_flag,
    array[] real t_eval,
    array[] real p_ode, // Now passed directly, no calculation inside needed
    vector y0_params    // [log_N0, log_D0]
  ) {
    
    real log_theta_cap = log(p_ode[1]) + 1.0;
    real safe_IC_N = softcap(y0_params[1], log_theta_cap);
    real safe_IC_D = softcap(y0_params[2], log_theta_cap);

    vector[3] y_start_main;
    
    // Starvation Logic
    if (starve_flag == 1) {
       vector[3] y0_starve = [exp(safe_IC_N), exp(safe_IC_D), 0.0]';
       array[1] real t_starve = {0.0};
       array[1] vector[3] y_res = ode_bdf_tol(model_b_ode, y0_starve, -6.0, t_starve, 1e-4, 1e-5, 100000, p_ode);
       y_start_main = y_res[1];
    } else {
       y_start_main = [exp(safe_IC_N), exp(safe_IC_D), 0.0]';
    }
    y_start_main[3] = G0; 

    // [FIX] Sanitise time grid to avoid t=0 crash
    array[size(t_eval)] real t_safe = t_eval;
    if (abs(t_safe[1]) < 1e-10) {
        t_safe[1] = 1e-8; 
    }

    // Safety check for ICs
    if (is_nan(y_start_main[1]) || is_inf(y_start_main[1])) {
        int N = size(t_safe);
        return rep_array(rep_vector(1e-6, 3), N);
    }
    
    return ode_bdf_tol(model_b_ode, y_start_main, 0.0, t_safe, 1e-4, 1e-5, 100000, p_ode);
  }

  // --- LIKELIHOOD REDUCER ---
  real partial_sum_lpmf(
      array[] int slice_wells,
      int start, int end,
      // Data args
      array[] int exp_id, vector G0_per_well, 
      array[] real t_grid,
      array[] int w_idx_count, array[] int g_idx_count, array[] int N_obs, array[] int D_obs,
      array[] int w_idx_gluc, array[] int g_idx_gluc, array[] real lum_obs, array[] real dilution, 
      array[] int is_censored, array[] int has_starvation,
      vector calib_a_fixed, vector calib_b_fixed, vector calib_sigma_fixed,
      // Parameters (Now simplified)
      array[] real p_ode,
      vector y0_params,
      real phi_total, real phi_frac
  ) {
      real log_lik = 0;
      
      for (i in 1:size(slice_wells)) {
          int w = slice_wells[i];
          
          array[size(t_grid)] vector[3] y_hat;
          
          // No line_id or ploidy_metric passed here anymore
          y_hat = solve_well_trajectory(
              G0_per_well[w], has_starvation[w],
              t_grid, 
              p_ode,
              y0_params
          );

          for (n in 1:size(w_idx_count)) {
             if (w_idx_count[n] == w) {
                 int idx = g_idx_count[n];
                 real NL = soft_lower(y_hat[idx, 1], 1e-12);
                 real ND = soft_lower(y_hat[idx, 2], 1e-12);
                 real total = NL + ND;
                 real p = softcap(soft_lower(NL/(total+1e-18), 1e-6), 1.0-1e-6);
                 
                 log_lik += neg_binomial_2_lpmf(N_obs[n] + D_obs[n] | total, phi_total);
                 log_lik += beta_binomial_lpmf(N_obs[n] | N_obs[n] + D_obs[n], p * phi_frac, (1-p)*phi_frac);
             }
          }

          for (n in 1:size(w_idx_gluc)) {
             if (w_idx_gluc[n] == w) {
                 int idx = g_idx_gluc[n];
                 real G = log1p_exp(100.0 * y_hat[idx, 3]) / 100.0;
                 int e = exp_id[w];
                 real mu = calib_a_fixed[e] * G * dilution[n] + calib_b_fixed[e];
                 
                 if (is_censored[n]) 
                    log_lik += lognormal_lcdf(lum_obs[n] | log(mu+1e-12), calib_sigma_fixed[e]);
                 else 
                    log_lik += lognormal_lpdf(lum_obs[n] | log(mu+1e-12), calib_sigma_fixed[e]);
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
  int<lower=1> N_exps;
  int<lower=1> N_obs_count;
  int<lower=1> N_obs_gluc;

  int<lower=1> N_grid;
  array[N_grid] real t_grid;
  
  // Removed: line_id, ploidy_metric
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

  // Calibration data removed for brevity if not fitting calib, but kept calib params
  // If you are fitting calibration inside here, keep the block from original.
  // Assuming calibration parameters are fixed constants passed in:
  
  vector[N_wells] G0_per_well;
  int<lower=1> grainsize;

  real prior_mu_N0_mean;
  real prior_mu_N0_sd;
  real prior_mu_D0_mean;
  real prior_mu_D0_sd;

  vector[12] prior_ode_mean;
  vector[12] prior_ode_sd;
  
  vector<lower=0>[N_exps] calib_a_fixed;
  vector<lower=0>[N_exps] calib_b_fixed;
  vector<lower=0>[N_exps] calib_sigma_fixed;

  int<lower=0,upper=1> calc_sim;
}

transformed data {
  vector[N_params] log_lower = log(to_vector(lower_b));
  vector[N_params] log_upper = log(to_vector(upper_b));
}

parameters {
  // Single set of global parameters for this Line/Ploidy combo
  vector<lower=log_lower, upper=log_upper>[N_params] logp;
  
  // Single set of Initial Condition parameters (Pooled ICs)
  vector[2] mu_IC_raw; 
  
  real<lower=1e-4> phi_total;
  real<lower=1e-4> phi_frac;
}

transformed parameters {
  // 1. Calculate the single physical parameter set here (Optimization)
  vector[N_params] p_phys = exp(logp);
  array[N_params] real p_ode = get_p_ode(p_phys);


  // 2. Calculate the ICs
  vector[2] mu_IC;
  mu_IC[1] = prior_mu_N0_mean + mu_IC_raw[1] * prior_mu_N0_sd;
  mu_IC[2] = prior_mu_D0_mean + mu_IC_raw[2] * prior_mu_D0_sd;
}

model {
  // Priors
  logp ~ normal(prior_ode_mean, prior_ode_sd);
  mu_IC_raw ~ std_normal();
  phi_total ~ exponential(0.1); 
  phi_frac ~ exponential(0.1); 
  
  array[N_wells] int seq_wells;
  for (i in 1:N_wells) seq_wells[i] = i;
  
  // Arguments unpacked because reduce_sum cannot deep copy tuples
  target += reduce_sum(
    partial_sum_lpmf,
    seq_wells,
    grainsize,
    // Data
    exp_id, G0_per_well, t_grid,
    well_idx_count, grid_idx_count, N_obs, D_obs,
    well_idx_gluc, grid_idx_gluc, lum_obs, dilution,
    is_censored, has_starvation,
    calib_a_fixed, calib_b_fixed, calib_sigma_fixed,
    // Params (Transformed and passed directly)
    p_ode,
    mu_IC, 
    phi_total, phi_frac
  );
}

generated quantities {
  array[N_wells, N_grid] vector[3] y_sim;
  
  if (calc_sim == 1) {
    for (w in 1:N_wells) {
       array[N_grid] vector[3] y_hat;
       
       // Call solver with the single global p_ode and mu_IC
       y_hat = solve_well_trajectory(
          G0_per_well[w], has_starvation[w],
          t_grid, p_ode, mu_IC
       );
       
       for (g in 1:N_grid) {
         y_sim[w, g, 1] = soft_lower(y_hat[g, 1], 1e-12);
         y_sim[w, g, 2] = soft_lower(y_hat[g, 2], 1e-12);
         y_sim[w, g, 3] = log1p_exp(100.0 * y_hat[g, 3]) / 100.0;
       }
    }
  }
}
