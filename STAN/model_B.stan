functions {
  // Smooth, differentiable upper cap for log-parameters.
  // For x << cap: ~ x
  // For x >> cap: ~ cap
  real softcap(real x, real cap) {
    return cap - log1p_exp(cap - x);
  }

  // ============================================================
  // ODE in transformed state space
  // y = [u, v, w]
  //   NL = exp(u)
  //   ND = exp(v)
  //   s_raw = inv_logit(w) in (0,1)
  //   s = s_eps + (1 - 2*s_eps) * s_raw in (s_eps, 1-s_eps)
  //   G = G_floor + (G0 - G_floor) * s   in (G_floor, G0)
  //
  // p = [theta,kp,kd,kd2,g50a,na,g50d,nd,v1,v2, G0, G_floor, s_eps]
  // ============================================================
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

    real G0      = p[11];
    real G_floor = p[12];
    real s_eps   = p[13];

    real NL = exp(y[1]);
    real ND = exp(y[2]);

    real s_raw = inv_logit(y[3]);
    real s     = s_eps + (1.0 - 2.0 * s_eps) * s_raw;
    real G     = G_floor + (G0 - G_floor) * s;

    real eps = 1e-12;

    real actA   = 1.0 / (1.0 + pow(g50a / (G + eps), na));
    real term_d = 1.0 / (1.0 + pow(g50d / (G + eps), nd));
    real inhD   = 1.0 - term_d;

    real confl = kd2 * (NL * NL) / (theta + eps);

    real dNL = kp * NL * (1.0 - NL / (theta + eps)) * actA - kd * NL * inhD - confl;
    real dND = kd * NL * inhD + confl;

    // taper -> 0 as G -> G_floor (s -> s_eps)
    real taper = (s - s_eps) / (1.0 - 2.0 * s_eps + eps);

    real dG  = -NL * (v1 * actA + v2 * term_d) / 2.0 * taper;

    // chain rule denom; s_raw*(1-s_raw) never exactly 0
    real denom = (G0 - G_floor) * (1.0 - 2.0 * s_eps) * s_raw * (1.0 - s_raw) + eps;

    vector[3] dydt;
    dydt[1] = dNL / (NL + eps);
    dydt[2] = dND / (ND + eps);
    dydt[3] = dG  / denom;

    return dydt;
  }

  // ============================================================
  // Parallel partial sum
  // ============================================================
  real partial_sum_lpmf(
      array[] int slice_wells,
      int start, int end,

      // Data
      array[] int line_id, array[] int high_p, array[] int exp_id, vector G0_per_well, array[] real t_grid,
      array[] int w_idx_count, array[] int g_idx_count, array[] int N_obs, array[] int D_obs,
      array[] int w_idx_gluc, array[] int g_idx_gluc, array[] real lum_obs, array[] real dilution,

      // Fixed calibration
      vector calib_a_fixed, vector calib_b_fixed,

      // Parameters
      vector mu_global, vector sigma_line, vector beta_high, matrix z_line,
      vector mu_IC, vector sigma_IC, vector beta_IC, matrix z_IC,
      vector calib_sigma,
      real phi_N, real phi_D
  ) {
    real log_lik = 0;

    int N_grid = size(t_grid);
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;

    // Smooth caps on log-params before exp()
    real cap_log_main = 40.0;  // exp(40) ~ 2e17 (finite; already "astronomical")
    real cap_log_hill = 6.0;   // exp(6) ~ 403, so na/nd <= ~404+1

    for (i in 1:size(slice_wells)) {
      int w = slice_wells[i];
      int l = line_id[w];
      int h = high_p[w];

      array[13] real p_w;

      for (pp in 1:10) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;

        if (pp == 6 || pp == 8) {
          p_w[pp] = 1.0 + exp(softcap(raw, cap_log_hill));
        } else {
          p_w[pp] = exp(softcap(raw, cap_log_main));
        }
      }

      real G0 = G0_per_well[w];
      real G_floor = fmin(1e-12, 0.5 * G0);
      real s_eps = 1e-6;

      p_w[11] = G0;
      p_w[12] = G_floor;
      p_w[13] = s_eps;

      // ICs directly in log space (avoids exp->inf->log)
      vector[3] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = 5.0;

      array[N_grid] vector[3] y_hat;
      y_hat = ode_bdf_tol(model_b_ode, y0_w, 0.0, t_eval, 1e-6, 1e-6, 100000, p_w);

      for (n in 1:size(w_idx_count)) {
        if (w_idx_count[n] == w) {
          int idx = g_idx_count[n];
          real NL_hat = exp(y_hat[idx, 1]) + 1e-12;
          real ND_hat = exp(y_hat[idx, 2]) + 1e-12;

          log_lik += neg_binomial_2_lpmf(N_obs[n] | NL_hat, phi_N);
          log_lik += neg_binomial_2_lpmf(D_obs[n] | ND_hat, phi_D);
        }
      }

      for (n in 1:size(w_idx_gluc)) {
        if (w_idx_gluc[n] == w) {
          int idx = g_idx_gluc[n];

          real s_raw = inv_logit(y_hat[idx, 3]);
          real s     = s_eps + (1.0 - 2.0 * s_eps) * s_raw;
          real G_hat = G_floor + (G0 - G_floor) * s;

          int e = exp_id[w];
          real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e] + 1e-12;

          log_lik += lognormal_lpdf(lum_obs[n] | log(mu), calib_sigma[e]);
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

  vector<lower=0>[N_exps] calib_a_fixed;
  vector<lower=0>[N_exps] calib_b_fixed;

  int<lower=0,upper=2> mode;
  // 0 = posterior (likelihood on; ODE in GQ)
  // 1 = prior draws, parameter-only (no likelihood; NO ODE anywhere)
  // 2 = prior predictive (no likelihood; ODE in GQ)
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

  vector<lower=0>[N_exps] calib_sigma;

  real<lower=0> phi_N;
  real<lower=0> phi_D;
}

model {
  // Priors
  mu_global ~ normal(0, 2);
  sigma_line ~ exponential(1);
  to_vector(z_line) ~ std_normal();
  beta_high ~ normal(0, 1);

  mu_IC[1] ~ normal(prior_mu_N0_mean, prior_mu_N0_sd);
  mu_IC[2] ~ normal(prior_mu_D0_mean, prior_mu_D0_sd);
  sigma_IC ~ exponential(1);
  to_vector(z_IC) ~ std_normal();
  beta_IC ~ normal(0, 1);

  phi_N ~ exponential(0.1);
  phi_D ~ exponential(0.1);

  calib_sigma ~ exponential(1);

  // Likelihood
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
      calib_sigma,
      phi_N, phi_D
    );
  }
}

generated quantities {
  array[N_wells, N_grid] vector[3] y_sim;
  real log_lik = 0;

  // ------------------------------------------------------------
  // Always initialize so y_sim is defined even when we skip ODEs
  // ------------------------------------------------------------
  for (w in 1:N_wells)
    for (g in 1:N_grid)
      y_sim[w, g] = rep_vector(0.0, 3);

  // ------------------------------------------------------------
  // Modes:
  //   mode = 0 : posterior  (ODE + log_lik)
  //   mode = 1 : prior only (NO ODE, no log_lik)
  //   mode = 2 : prior predictive (ODE, no log_lik)
  // ------------------------------------------------------------
  if (mode != 1) {
    // Solve ODE for every well to build predictions
    array[N_grid] real t_eval = t_grid;
    if (abs(t_eval[1]) < 1e-14) t_eval[1] = 1e-8;

    real cap_log_main = 40.0;
    real cap_log_hill = 6.0;

    for (w in 1:N_wells) {
      int l = line_id[w];
      int h = is_high_ploidy[w];

      array[13] real p_w;

      // top-10 ODE params (soft-capped, with Hill exponents >= 1)
      for (pp in 1:10) {
        real raw = mu_global[pp] + sigma_line[pp] * z_line[pp, l] + beta_high[pp] * h;
        if (pp == 6 || pp == 8) p_w[pp] = 1.0 + exp(softcap(raw, cap_log_hill));
        else                    p_w[pp] = exp(softcap(raw, cap_log_main));
      }

      // per-well glucose anchors
      real G0 = G0_per_well[w];
      real G_floor = fmin(1e-12, 0.5 * G0);
      real s_eps = 1e-6;
      p_w[11] = G0;
      p_w[12] = G_floor;
      p_w[13] = s_eps;

      // initial conditions in transformed space (u,v,w)
      vector[3] y0_w;
      y0_w[1] = mu_IC[1] + sigma_IC[1] * z_IC[1, l] + beta_IC[1] * h;
      y0_w[2] = mu_IC[2] + sigma_IC[2] * z_IC[2, l] + beta_IC[2] * h;
      y0_w[3] = 5.0;

      // integrate
      array[N_grid] vector[3] y_hat =
        ode_bdf_tol(model_b_ode, y0_w, 0.0, t_eval, 1e-6, 1e-6, 100000, p_w);

      // back-transform and store NL, ND, G
      for (g in 1:N_grid) {
        real NL_hat = exp(y_hat[g, 1]);
        real ND_hat = exp(y_hat[g, 2]);

        real s_raw = inv_logit(y_hat[g, 3]);
        real s     = s_eps + (1.0 - 2.0 * s_eps) * s_raw;
        real G_hat = G_floor + (G0 - G_floor) * s;

        y_sim[w, g, 1] = NL_hat;
        y_sim[w, g, 2] = ND_hat;
        y_sim[w, g, 3] = G_hat;
      }
    }
  }

  // Only compute full-data log-likelihood in posterior mode
  if (mode == 0) {
    // Count data (NegBin2 on NL, ND)
    for (n in 1:N_obs_count) {
      int w = well_idx_count[n];
      int g = grid_idx_count[n];
      real NL_hat = y_sim[w, g, 1] + 1e-12;
      real ND_hat = y_sim[w, g, 2] + 1e-12;
      log_lik += neg_binomial_2_lpmf(N_obs[n] | NL_hat, phi_N);
      log_lik += neg_binomial_2_lpmf(D_obs[n] | ND_hat, phi_D);
    }

    // Glucose lum (lognormal on a_e*G_hat*dilution + b_e)
    for (n in 1:N_obs_gluc) {
      int w = well_idx_gluc[n];
      int g = grid_idx_gluc[n];
      int e = exp_id[w];
      real G_hat = y_sim[w, g, 3];
      real mu = calib_a_fixed[e] * G_hat * dilution[n] + calib_b_fixed[e] + 1e-12;
      log_lik += lognormal_lpdf(lum_obs[n] | log(mu), calib_sigma[e]);
    }
  }
}




