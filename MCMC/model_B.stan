functions {
  real softcap(real x, real cap) {
    return cap - log1p_exp(cap - x);
  }

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
    // G0 is no longer needed inside the ODE for normalization

    real NL = exp(y[1]);
    real ND = exp(y[2]);
    
    // --- CHANGED: Absolute Glucose Logic ---
    real k_smooth = 100.0;
    // y[3] is now absolute glucose (mM). Apply smoothing directly.
    real G = log1p_exp(k_smooth * y[3]) / k_smooth;

    real log_G = log(G + 1e-18);
    real actA   = inv_logit(na * (log_G - log(g50a)));
    real term_d = inv_logit(nd * (log_G - log(g50d)));
    real inhD   = 1.0 - term_d;

    real du_dt = kp * (1.0 - NL/theta) * actA - kd * inhD - kd2 * NL/theta;
    real dv_dt = (kd * NL * inhD + kd2 * (NL * NL) / theta) / ND;
    
    // --- CHANGED: Derivative is not divided by G0 ---
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
      array[] int line_id, array[] int high_p, array[] int exp_id, vector G0_per_well, array[] real t_grid,
      array[] int w_idx_count, array[] int g_idx_count, array[] int N_obs, array[] int D_obs,
      array[] int w_idx_gluc, array[] int g_idx_gluc, array[] real lum_obs, array[] real dilution,
      vector calib_a_fixed, vector calib_b_fixed,
      vector mu_global, vector sigma_line, 
      vector beta_high, matrix z_line,
      vector mu_IC, vector sigma_IC, vector beta_IC, matrix z_IC,
      vector calib_sigma_fixed,
      real phi_total, real phi_frac
  ) {
    real log_lik = 0;
    int N_grid = size(t_grid);
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;
    real cap_log_main = 40.0;
    real cap_log_hill = 1.6;

    for (i in 1:size(slice_wells)) {
      int w = slice_wells[i];
      int l = line_id[w];
      int h = high_p[w];

      array[10] real p_w; // CHANGED: Size reduced to 10 (G0 removed)
      for (pp in 1:10) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;
        if (pp == 6 || pp == 8) {
          p_w[pp] = 1.0 + exp(softcap(raw, cap_log_hill));
        } else {
          p_w[pp] = exp(softcap(raw, cap_log_main));
        }
      }

      real G0 = G0_per_well[w];
      // p_w[11] = G0; // REMOVED

      vector[3] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = G0; // CHANGED: Initial Glucose is now G0 (Absolute), not 1.0

      array[N_grid] vector[3] y_hat;
      y_hat = ode_bdf_tol(model_b_ode, y0_w, 0.0, t_eval, 1e-4, 1e-5, 50000, p_w);
      
      for (n in 1:size(w_idx_count)) {
        if (w_idx_count[n] == w) {
          int idx = g_idx_count[n];
          real NL_hat = exp(y_hat[idx, 1]);
          real ND_hat = exp(y_hat[idx, 2]);
          
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
          // CHANGED: Use y_hat[3] directly (absolute)
          real G_hat = log1p_exp(k_smooth * y_hat[idx, 3]) / k_smooth;
          
          int e = exp_id[w];
          real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e];
          log_lik += lognormal_lpdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
        }
      }
    }
    return log_lik;
  }
}

data {
  int<lower=1> N_wells;
  int<lower=1> N_lines;
  int<lower=1> N_exps;
  int<lower=1> N_obs_count;
  int<lower=1> N_obs_gluc;
  int<lower=1> N_obs_calib;

  int<lower=1> N_grid;
  array[N_grid] real t_grid;
  array[N_wells] int line_id;
  array[N_wells] int is_high_ploidy;
  array[N_wells] int exp_id;

  array[N_obs_count] int well_idx_count;
  array[N_obs_count] int grid_idx_count;
  array[N_obs_count] int N_obs;
  array[N_obs_count] int D_obs;

  array[N_obs_gluc] int well_idx_gluc;
  array[N_obs_gluc] int grid_idx_gluc;
  array[N_obs_gluc] real lum_obs;
  array[N_obs_gluc] real dilution;

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

parameters {
  vector[10] mu_global;
  vector<lower=0>[10] sigma_line;
  vector[10] beta_high;
  matrix[10, N_lines] z_line;
  vector[2] mu_IC;
  vector<lower=0>[2] sigma_IC;
  vector[2] beta_IC;
  matrix[2, N_lines] z_IC;

  real<lower=0> phi_total;
  real<lower=0> phi_frac;
}

model {
  mu_global ~ normal(prior_ode_mean, prior_ode_sd);
  sigma_line ~ exponential(1);
  to_vector(z_line) ~ std_normal();
  beta_high ~ normal(0, 1);
  mu_IC[1] ~ normal(prior_mu_N0_mean, prior_mu_N0_sd);
  mu_IC[2] ~ normal(prior_mu_D0_mean, prior_mu_D0_sd);
  sigma_IC ~ exponential(1);
  to_vector(z_IC) ~ std_normal();
  beta_IC ~ normal(0, 1);
  
  phi_total ~ exponential(0.1); 
  phi_frac ~ exponential(0.1); 
  
  if (mode == 0) {
    array[N_wells] int seq_wells;
    for (i in 1:N_wells) seq_wells[i] = i;
    target += reduce_sum(
      partial_sum_lpmf,
      seq_wells,
      grainsize,
      line_id, is_high_ploidy, exp_id, G0_per_well, t_grid,
      well_idx_count, grid_idx_count, N_obs, D_obs,
      well_idx_gluc, grid_idx_gluc, lum_obs, dilution,
      calib_a_fixed, calib_b_fixed,
      mu_global, sigma_line, beta_high, z_line,
      mu_IC, sigma_IC, beta_IC, z_IC,
      calib_sigma_fixed,
      phi_total, phi_frac
    );
  }
}

generated quantities {
  array[N_wells, N_grid] vector[3] y_sim;
  real log_lik = 0;
  if (mode != 1 && calc_sim == 1) {
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;
    real cap_log_main = 40.0;
    real cap_log_hill = 3.0;
    for (w in 1:N_wells) {
      int l = line_id[w];
      int h = is_high_ploidy[w];
      array[10] real p_w; // CHANGED

      for (pp in 1:10) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;
        if (pp == 6 || pp == 8) p_w[pp] = 1.0 + exp(softcap(raw, cap_log_hill));
        else                    p_w[pp] = exp(softcap(raw, cap_log_main));
      }

      real G0 = G0_per_well[w];
      // p_w[11] = G0; // REMOVED

      vector[3] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = G0; // CHANGED: Absolute G0
      
      array[N_grid] vector[3] y_hat =
        ode_bdf_tol(model_b_ode, y0_w, 0.0, t_eval, 1e-4, 1e-5, 50000, p_w);
      for (g in 1:N_grid) {
        real NL_hat = exp(y_hat[g, 1]);
        real ND_hat = exp(y_hat[g, 2]);
        real k_smooth = 100.0; 
        // CHANGED: Absolute smoothing directly on output
        real G_hat = log1p_exp(k_smooth * y_hat[g, 3]) / k_smooth;

        y_sim[w, g, 1] = NL_hat;
        y_sim[w, g, 2] = ND_hat;
        y_sim[w, g, 3] = G_hat;
      }
    }
  }

  if (mode == 0) {
    if (calc_sim == 1) {
       for (n in 1:N_obs_count) {
          int w = well_idx_count[n];
          int g = grid_idx_count[n];
          real NL_hat = y_sim[w, g, 1];
          real ND_hat = y_sim[w, g, 2];
          
          real total_hat = NL_hat + ND_hat;
          real p_hat     = NL_hat / (total_hat + 1e-18);
          p_hat = fmin(fmax(p_hat, 1e-6), 1.0 - 1e-6);
          int total_obs = N_obs[n] + D_obs[n];
          
          real alpha_p = p_hat * phi_frac;
          real beta_p  = (1.0 - p_hat) * phi_frac;

          log_lik += neg_binomial_2_lpmf(total_obs | total_hat, phi_total);
          log_lik += beta_binomial_lpmf(N_obs[n] | total_obs, alpha_p, beta_p);
       }
       for (n in 1:N_obs_gluc) {
          int w = well_idx_gluc[n];
          int g = grid_idx_gluc[n];
          int e = exp_id[w];
          real G_hat = y_sim[w, g, 3];
          real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e];
          log_lik += lognormal_lpdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
       }
    }
  }
}