functions {
  real softcap(real x, real cap) {
    return cap - log1p_exp(cap - x);
  }

  vector model_q_ode(real t, vector y, array[] real p) {
    real theta = p[1];
    real kp    = p[2];
    real kd    = p[3];
    real kd2   = p[4];
    real nu    = p[5];
    real rho   = p[6];
    real r     = p[7];
    real sigma_G = p[8];
    real KG      = p[9];
    real q50g_frac = p[10];
    real q50d_frac = p[11];
    real na = p[12];
    real nd = p[13];

    real k_home_cap = 4.0;
    real k_home = exp(softcap(log(nu * kp + 1e-12), log(k_home_cap)));
    
    real s = nu * rho;
    real alpha_c = s * r;
    real alpha_m = s * (1.0 - r);

    real eps_qstar = 1e-6;
    real q_star = (nu * (1.0 - rho)) / (nu + 1.0);
    q_star = fmin(fmax(q_star, eps_qstar), 1.0 - eps_qstar);
    real q50gN = q50g_frac * q_star;
    real q50dN = q50d_frac * q_star;

    real NL = exp(y[1]);
    real ND = exp(y[2]);
    real k_smooth_G = 100.0; // CHANGED: Consistent smoothing constant
    real G_raw = y[3];
    real G = log1p_exp(k_smooth_G * G_raw) / k_smooth_G;
    
    real k_smooth_q = 50.0;
    real qN_raw = y[4];
    real qN_low = log1p_exp(k_smooth_q * qN_raw) / k_smooth_q;
    real qN = 1.0 - (log1p_exp(k_smooth_q * (1.0 - qN_low)) / k_smooth_q);
    real drive = log1p_exp(k_smooth_q * (1.0 - qN)) / k_smooth_q;
    real gate  = G / (KG + G + 1e-12);
    real Jin = k_home * gate * drive;
    real log_qN = log(qN + 1e-12);
    real reg_growth  = inv_logit(na * (log_qN - log(q50gN + 1e-12)));
    real term_d_hill = inv_logit(nd * (log_qN - log(q50dN + 1e-12)));
    
    real mu    = kp * reg_growth;
    real delta = kd * (1.0 - term_d_hill);

    real b = mu * (1.0 - NL/theta);
    real d = delta + kd2 * NL/theta;

    real du = b - d;
    real dv = (delta * NL + kd2 * (NL*NL)/theta) / ND;
    real dG = -NL * Jin / sigma_G;
    real dqN = Jin - (alpha_m * kp) - (qN * b) - (alpha_c * b);
    return [du, dv, dG, dqN]';
  }

  real partial_sum_lpmf(
      array[] int slice_wells,
      int start, int end,
      array[] int line_id, array[] int high_p, array[] int exp_id, vector G0_per_well, array[] real t_grid,
      array[] int w_idx_count, array[] int g_idx_count, array[] int N_obs, array[] int D_obs,
      array[] int w_idx_gluc, array[] int g_idx_gluc, array[] real lum_obs, array[] real dilution,
      vector calib_a_fixed, vector calib_b_fixed,
      vector mu_global, vector sigma_line, vector beta_high, matrix z_line,
      vector mu_IC, vector sigma_IC, vector beta_IC, matrix z_IC,
      vector calib_sigma_fixed,
      real phi_total, real phi_frac
  ) {
    real log_lik = 0;
    int N_grid = size(t_grid);
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;

    real cap_log_main = 15.0;
    real cap_log_hill = 2.7;

    for (i in 1:size(slice_wells)) {
      int w = slice_wells[i];
      int l = line_id[w];
      int h = high_p[w];

      array[13] real p_w;
      for (pp in 1:13) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;
        if (pp == 6 || pp == 7 || pp == 10 || pp == 11) {
          p_w[pp] = inv_logit(raw);
        } else if (pp == 12 || pp == 13) {
          p_w[pp] = exp(softcap(raw, cap_log_hill));
        } else if (pp == 1) {
          p_w[pp] = exp(raw);
        } else {
          p_w[pp] = exp(softcap(raw, cap_log_main));
        }
      }

      real G0 = G0_per_well[w];
      vector[4] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = G0;
      y0_w[4] = 1.0;
      array[N_grid] vector[4] y_hat;
      y_hat = ode_bdf_tol(model_q_ode, y0_w, 0.0, t_eval, 1e-3, 1e-4, 5000, p_w);
      
      for (n in 1:size(w_idx_count)) if (w_idx_count[n] == w) {
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

      for (n in 1:size(w_idx_gluc)) if (w_idx_gluc[n] == w) {
        int idx = g_idx_gluc[n];
        real k_smooth_G = 100.0; // CHANGED: Consistent smoothing constant
        real G_hat = log1p_exp(k_smooth_G * y_hat[idx, 3]) / k_smooth_G;
        int e = exp_id[w];
        real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e];
        log_lik += lognormal_lpdf(lum_obs[n] | log(mu + 1e-12), calib_sigma_fixed[e]);
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
  vector[13] prior_ode_mean;
  vector[13] prior_ode_sd;
  vector<lower=0>[N_exps] calib_a_fixed;
  vector<lower=0>[N_exps] calib_b_fixed;
  vector<lower=0>[N_exps] calib_sigma_fixed;
  int<lower=0,upper=2> mode;
  int<lower=0,upper=1> calc_sim;
}

parameters {
  vector[13] mu_global;
  vector<lower=0>[13] sigma_line;
  vector[13] beta_high;
  matrix[13, N_lines] z_line;
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
  array[N_wells, N_grid] vector[4] y_sim;
  real log_lik = 0;
  if (mode != 1 && calc_sim == 1) {
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;
    real cap_log_main = 15.0;
    real cap_log_hill = 2.7;
    for (w in 1:N_wells) {
      int l = line_id[w];
      int h = is_high_ploidy[w];
      array[13] real p_w;
      for (pp in 1:13) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;
        if (pp == 6 || pp == 7 || pp == 10 || pp == 11) {
          p_w[pp] = inv_logit(raw);
        } else if (pp == 12 || pp == 13) {
          p_w[pp] = exp(softcap(raw, cap_log_hill));
        } else if (pp == 1) {
          p_w[pp] = exp(raw);
        } else {
          p_w[pp] = exp(softcap(raw, cap_log_main));
        }
      }

      real G0 = G0_per_well[w];
      vector[4] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = G0;
      y0_w[4] = 1.0;
      array[N_grid] vector[4] y_hat =
        ode_bdf_tol(model_q_ode, y0_w, 0.0, t_eval, 1e-3, 1e-4, 5000, p_w);
      for (g in 1:N_grid) {
        y_sim[w, g, 1] = exp(y_hat[g, 1]);
        y_sim[w, g, 2] = exp(y_hat[g, 2]);
        real k_smooth_G = 100.0;
        y_sim[w, g, 3] = log1p_exp(k_smooth_G * y_hat[g, 3]) / k_smooth_G;
        real k_smooth_q = 50.0;
        real qN_low = log1p_exp(k_smooth_q * y_hat[g, 4]) / k_smooth_q;
        y_sim[w, g, 4] = 1.0 - (log1p_exp(k_smooth_q * (1.0 - qN_low)) / k_smooth_q);
      }
    }
    
    if (mode == 0) {
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