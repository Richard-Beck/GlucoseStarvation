functions {
  real smooth_min_vec(vector v) {
    real k = 10.0; 
    return -log_sum_exp(-k * v) / k;
  }
  
  // Generalized ODE System
// Generalized ODE System
// Generalized ODE System
  vector ode_general(real t, vector y, 
                     int R, int P, int W,
                     vector ae, vector ah, vector Y_R, matrix A_mat, matrix K_mat, 
                     real m, vector waste_mech) {
    
    real N_cells = fmax(y[1], 1e-12);
    
    vector[R] R_vec;
    for(i in 1:R) R_vec[i] = fmax(y[1+i], 0.0);
    
    // Uptake and Pathway Fluxes
    vector[R] u = (ae .* R_vec) ./ (ah + R_vec);
    vector[P] Phi = A_mat * (Y_R .* u);
    
    real mu_growth = smooth_min_vec(Phi);
    
    real inhibition = 1.0;
    real death = 0.0;
    vector[W] dW = rep_vector(0.0, W);
    
    // Waste Kinetics and Toxicity
    if (W > 0) {
      vector[W] W_vec;
      for(i in 1:W) W_vec[i] = fmax(y[1+R+i], 0.0);
      
      vector[W] mult_mask = 1.0 - waste_mech;
      vector[W] add_mask = waste_mech;
      
      inhibition = 1.0 / (1.0 + dot_product(mult_mask, W_vec));
      death = dot_product(add_mask, W_vec);
      
      dW = (K_mat * u) * N_cells;
    }
    
    real g = (mu_growth * inhibition) - m - death;
    
    vector[1 + R + W] dydt;
    dydt[1] = g * N_cells;
    for(i in 1:R) dydt[1+i] = -u[i] * N_cells;
    for(i in 1:W) dydt[1+R+i] = dW[i];
    
    return dydt;
  }

  // Solver Wrapper
  array[] vector solve_well(real N0, vector R_init,
                            array[] real t_eval, 
                            int R, int P, int W,
                            vector ae, vector ah, vector Y_R, matrix A_mat, matrix K_mat, 
                            real m, vector waste_mech) {
    vector[1 + R + W] y0;
    y0[1] = fmax(N0, 1e-6);
    for(i in 1:R) y0[1+i] = fmax(R_init[i], 0.0);
    for(i in 1:W) y0[1+R+i] = 0.0; 

    array[size(t_eval)] real t_safe = t_eval;
    if (abs(t_safe[1]) < 1e-12) t_safe[1] = 1e-8;

    return ode_bdf_tol(ode_general, y0, 0.0, t_safe, 1e-6, 1e-6, 50000, 
                       R, P, W, ae, ah, Y_R, A_mat, K_mat, m, waste_mech);
  }

// Partial Sum Function
  real partial_sum_likelihood(
      array[] int well_indices_slice, 
      int start, int end,              
      // Hierarchical Base Parameters
      matrix raw_theta_line,
      vector raw_theta_ploidy,
      vector raw_N0,
      real prior_N0_center,
      real nu_N,
      real sigma_N,
      // Mapping and Priors
      array[] int line_id,
      vector ploidy_metric,
      vector ploidy_effect_mask,
      array[] int group_id,
      vector prior_centers,
      vector prior_scales,
      array[] int param_map,
      vector param_mask,
      // Model Dimensions
      int R, int P, int W, int L,
      vector waste_mech,
      // Data
      vector G1_0, array[] int g1_id,
      array[] vector R_init_base,
      array[] int exp_id,
      array[] real t_grid,
      int N_grid,
      int N_obs_count,
      array[] int well_idx_count,
      array[] int grid_idx_count,
      array[] real N_obs,
      int N_obs_gluc,
      array[] int well_idx_gluc,
      array[] int grid_idx_gluc,
      array[] real lum_obs,
      array[] real dilution,
      array[] int is_censored,
      array[] int is_train,
      vector calib_a_fixed,
      vector calib_b_fixed,
      vector calib_sigma_fixed
  ) {
      real lp = 0;
      int size_slice = size(well_indices_slice);

      for (i in 1:size_slice) {
          int w = well_indices_slice[i];
          
          if (is_train[w] == 1) {
          
            real N0_w = exp(log(prior_N0_center) + raw_N0[group_id[w]] * 1.0);
            
            vector[L] raw_w = raw_theta_line[, line_id[w]] + (raw_theta_ploidy .* ploidy_effect_mask) * ploidy_metric[w];
            vector[L] theta_phys = exp(log(prior_centers) + raw_w .* prior_scales);
            
            vector[L] theta_mapped;
            for (k in 1:L) theta_mapped[k] = theta_phys[param_map[k]] * param_mask[k];
            
            int ptr = 1;
            vector[R] ae = theta_mapped[ptr : ptr+R-1]; ptr += R;
            vector[R] ah = theta_mapped[ptr : ptr+R-1]; ptr += R;
            
            // Unpack Y_R
            vector[R] Y_R = theta_mapped[ptr : ptr+R-1]; ptr += R;
            
            // Unpack and normalize A_mat with anchored diagonals
            matrix[P, R] A_mat;
            for(c in 1:R) {
              vector[P] raw_alloc;
              int ref_idx = c > P ? P : c; 
              
              for(r in 1:P) {
                if (r == ref_idx) {
                  raw_alloc[r] = 1.0; // Anchor the diagonal
                } else {
                  if (c == 1) {
                    // Constrain Column 1 off-diagonals to [0, 1) to break symmetry
                    // Preserves exactly 0.0 when strict_spec mask is applied
                    raw_alloc[r] = theta_mapped[ptr] / (1.0 + theta_mapped[ptr]); 
                  } else {
                    // Other columns remain unconstrained (0, infinity)
                    raw_alloc[r] = theta_mapped[ptr];
                  }
                  ptr += 1;
                }
              }
              // Normalize so the column sums to 1.0
              A_mat[:, c] = raw_alloc / sum(raw_alloc); 
            }
            
            matrix[W, R] K_mat;
            for(c in 1:R) {
              for(r in 1:W) {
                 if (W > 0) K_mat[r, c] = theta_mapped[ptr]; 
                 ptr += 1; 
              }
            }
            
            real m = theta_mapped[ptr];
            
            vector[R] R_init = R_init_base[w];
            R_init[1] = G1_0[g1_id[w]]; 
            
            array[N_grid] vector[1+R+W] yhat = solve_well(N0_w, R_init, t_grid, R, P, W, 
                                                          ae, ah, Y_R, A_mat, K_mat, m, waste_mech);
  
            for (n in 1:N_obs_count) {
                if (well_idx_count[n] == w) {
                    int ti = grid_idx_count[n];
                    real Nhat = fmax(yhat[ti,1], 1e-6); 
                    real log_N_obs = log(N_obs[n] + 1.0);
                    //lp += lognormal_lpdf(N_obs[n] + 1.0 | log(Nhat + 1.0), sigma_N);
                    lp += student_t_lpdf(log_N_obs | nu_N, log(Nhat + 1.0), sigma_N) - log_N_obs;
                }
            }
  
            for (n in 1:N_obs_gluc) {
                if (well_idx_gluc[n] == w) {
                    int ti = grid_idx_gluc[n];
                    real G1hat = fmax(yhat[ti,2], 0.0); 
                    
                    int e = exp_id[w];
                    real mu = calib_a_fixed[e] * G1hat * dilution[n] + calib_b_fixed[e];
                    mu = fmax(mu, 1e-12);
  
                    if (is_censored[n] == 1)
                      lp += lognormal_lcdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
                    else
                      lp += lognormal_lpdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
                }
            }
        } else {
            // Minimal evaluation for holdout wells at t=0 (Bypasses ODE solve)
            real N0_w = exp(log(prior_N0_center) + raw_N0[group_id[w]] * 1.0);
            real G1_0_w = G1_0[g1_id[w]];

            // Evaluate Initial Cell Counts
            for (n in 1:N_obs_count) {
                if (well_idx_count[n] == w) {
                    int ti = grid_idx_count[n];
                    if (t_grid[ti] == 0.0) { 
                        real Nhat = fmax(N0_w, 1e-6); 
                        //lp += lognormal_lpdf(N_obs[n] + 1.0 | log(Nhat + 1.0), sigma_N);
                        real log_N_obs = log(N_obs[n] + 1.0);
                        lp += student_t_lpdf(log_N_obs | nu_N, log(Nhat + 1.0), sigma_N) - log_N_obs;
                    }
                }
            }

            // Evaluate Initial Glucose Observations
            for (n in 1:N_obs_gluc) {
                if (well_idx_gluc[n] == w) {
                    int ti = grid_idx_gluc[n];
                    if (t_grid[ti] == 0.0) { 
                        real G1hat = fmax(G1_0_w, 0.0); 
                        
                        int e = exp_id[w];
                        real mu = calib_a_fixed[e] * G1hat * dilution[n] + calib_b_fixed[e];
                        mu = fmax(mu, 1e-12);

                        if (is_censored[n] == 1)
                          lp += lognormal_lcdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
                        else
                          lp += lognormal_lpdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
                    }
                }
            }
          }
      }
      return lp;
  }
}

data {
  // Model Dimensions
  int<lower=1> R; 
  int<lower=1> P; 
  int<lower=0> W;
  int<lower=1> L; // Total parameters: 2R + P*R + W*R + 1
  vector[W] waste_mech; // 1 = additive, 0 = multiplicative
  
  // Existing Data Arrays
  int<lower=1> N_wells;
  int<lower=1> N_exps;
  int<lower=1> N_obs_count;
  int<lower=1> N_obs_gluc;
  int<lower=1> N_grid;
  array[N_grid] real t_grid;

  array[N_wells] int<lower=1, upper=N_exps> exp_id;
  vector[N_wells] G0_per_well;
  array[N_wells] vector[R] R_init_base; // Initial conditions for all resources

  // HIERARCHICAL DATA
  int<lower=1> N_lines;                           
  array[N_wells] int<lower=1, upper=N_lines> line_id; 
  vector[N_wells] ploidy_metric;                  
  
  int<lower=1> N_groups;                          
  array[N_wells] int<lower=1, upper=N_groups> group_id; 

  // Observation Mappings
  array[N_obs_count] int<lower=1, upper=N_wells> well_idx_count;
  array[N_obs_count] int<lower=1, upper=N_grid>  grid_idx_count;
  array[N_obs_count] real<lower=0> N_obs;
  array[N_obs_count] real<lower=0> D_obs;

  array[N_obs_gluc] int<lower=1, upper=N_wells> well_idx_gluc;
  array[N_obs_gluc] int<lower=1, upper=N_grid>  grid_idx_gluc;
  array[N_obs_gluc] real<lower=0> lum_obs;
  array[N_obs_gluc] real<lower=0> dilution;
  array[N_obs_gluc] int<lower=0, upper=1> is_censored;

  vector[N_exps] calib_a_fixed;
  vector[N_exps] calib_b_fixed;
  vector[N_exps] calib_sigma_fixed;
  
  // Priors and Masks
  vector<lower=0>[L] prior_centers; 
  vector<lower=0>[L] prior_scales;
  array[L] int<lower=1, upper=L> param_map; 
  vector[L] param_mask; 
  vector[L] ploidy_effect_mask;
  real<lower=0> prior_N0_center;
  
  //CV masks
  array[N_wells] int<lower=0, upper=1> is_train;
  
  //infrastructure for G0 sharing.
  int<lower=1> N_G1;
  array[N_wells] int<lower=1, upper=N_G1> g1_id;
  array[N_G1] int<lower=1, upper=N_wells> g1_ref_well;
}

transformed data {
  array[N_wells] int well_indices;
  for (i in 1:N_wells) well_indices[i] = i;
}

parameters {
  matrix[L, N_lines] raw_theta_line; 
  vector[L] raw_theta_ploidy;        
  real<lower=2> nu_N;
  vector[N_groups] raw_N0;
  
  real raw_sigma_N; 
  real raw_sigma_base; 
  real raw_sigma_rel;
  
  //real<lower=0> sigma_base; 
  //real<lower=0> sigma_rel;
  //vector[N_wells] G1_0_raw;
  // 1. Shift G1_0 to CP and apply lower bound
  // In the parameters block:
  //vector<lower=0, upper=max(G0_per_well) * 2.0>[N_wells] G1_0;
  vector<lower=0, upper=max(G0_per_well) * 2.0>[N_G1] G1_0;
}

transformed parameters {
  real sigma_N = exp(log(0.2) + raw_sigma_N * 0.5);
  real sigma_base = exp(log(0.01) + 1.0 * raw_sigma_base);
  real sigma_rel = exp(log(0.05) + 0.2 * raw_sigma_rel);
  //vector[N_wells] G1_0;

  vector[N_wells] sigma_vec = sigma_base + sigma_rel * G0_per_well;
  //G1_0 = G0_per_well + sigma_vec .* G1_0_raw;
}

model {
  to_vector(raw_theta_line) ~ std_normal();
  raw_theta_ploidy ~ std_normal();
  raw_N0 ~ std_normal();
  nu_N ~ gamma(2, 0.1);
  raw_sigma_N ~ std_normal();
  //sigma_base ~ lognormal(log(0.01), 1.0); 
  //sigma_rel ~ normal(0.05, 0.01); 
  raw_sigma_base ~ std_normal();
  raw_sigma_rel ~ std_normal();
  //G1_0_raw ~ std_normal();
  // 3. Apply the centered prior
  //G1_0 ~ normal(G0_per_well, sigma_vec);
  for (j in 1:N_G1) G1_0[j] ~ normal(G0_per_well[g1_ref_well[j]], sigma_vec[g1_ref_well[j]]);
  target += reduce_sum(
      partial_sum_likelihood, well_indices, 1, 
      raw_theta_line, raw_theta_ploidy, raw_N0, prior_N0_center, nu_N,sigma_N,
      line_id, ploidy_metric, ploidy_effect_mask, group_id,
      prior_centers, prior_scales, param_map, param_mask,
      R, P, W, L, waste_mech, G1_0,g1_id, R_init_base, exp_id, t_grid, N_grid,
      N_obs_count, well_idx_count, grid_idx_count, N_obs,
      N_obs_gluc, well_idx_gluc, grid_idx_gluc, lum_obs, dilution, is_censored,is_train,
      calib_a_fixed, calib_b_fixed, calib_sigma_fixed
  );
}

generated quantities {
  real log_lik;
  real log_lik_train=0;           // sum over training wells only
  real log_lik_holdout=0;         // sum over held-out wells only
  vector[N_wells] ll_well = rep_vector(0.0, N_wells);

  for (w in 1:N_wells) {
      real N0_w = exp(log(prior_N0_center) + raw_N0[group_id[w]] * 1.0);
          
      vector[L] raw_w = raw_theta_line[, line_id[w]] + (raw_theta_ploidy .* ploidy_effect_mask) * ploidy_metric[w];
      vector[L] theta_phys = exp(log(prior_centers) + raw_w .* prior_scales);
      
      vector[L] theta_mapped;
      for (k in 1:L) theta_mapped[k] = theta_phys[param_map[k]] * param_mask[k];

// (Inside generated quantities...)
      int ptr = 1;
      vector[R] ae = theta_mapped[ptr : ptr+R-1]; ptr += R;
      vector[R] ah = theta_mapped[ptr : ptr+R-1]; ptr += R;
      vector[R] Y_R = theta_mapped[ptr : ptr+R-1]; ptr += R;
      
// Unpack and normalize A_mat
          matrix[P, R] A_mat;
          for(c in 1:R) {
            vector[P] raw_alloc;
            int ref_idx = c > P ? P : c; // For each column, a single target row index (ref_idx) is defined where the 1.0 anchor will be placed.
            
            for(r in 1:P) {
              if (r == ref_idx) {
                raw_alloc[r] = 1.0; // Anchor the diagonal
              } else {
                if (c == 1) {
                  // Constrain Column 1 off-diagonals to [0, 1) to break symmetry
                  // Preserves exactly 0.0 when strict_spec mask is applied
                  raw_alloc[r] = theta_mapped[ptr] / (1.0 + theta_mapped[ptr]); 
                } else {
                  // Other columns remain unconstrained (0, infinity)
                  raw_alloc[r] = theta_mapped[ptr];
                }
                ptr += 1;
              }
            }
            // Normalize so the column sums to 1.0
            A_mat[:, c] = raw_alloc / sum(raw_alloc); 
          } 
          matrix[W, R] K_mat;
          for(c in 1:R) {
            for(r in 1:W) {
               if (W > 0) K_mat[r, c] = theta_mapped[ptr]; 
               ptr += 1; 
            }
          }
          
          real m = theta_mapped[ptr];
      
      vector[R] R_init = R_init_base[w];
      R_init[1] = G1_0[g1_id[w]]; 
      
      array[N_grid] vector[1+R+W] yhat = solve_well(N0_w, R_init, t_grid, R, P, W, 
                                                    ae, ah, Y_R,A_mat, K_mat, m, waste_mech);

      for (n in 1:N_obs_count) {
          if (well_idx_count[n] == w) {
              int ti = grid_idx_count[n];
              real Nhat = fmax(yhat[ti,1], 1e-6); 
              real log_N_obs = log(N_obs[n] + 1.0);
              ll_well[w] += student_t_lpdf(log_N_obs | nu_N, log(Nhat + 1.0), sigma_N) - log_N_obs;
              //ll_well[w] += lognormal_lpdf(N_obs[n] + 1.0 | log(Nhat + 1.0), sigma_N);
          }
      }

      for (n in 1:N_obs_gluc) {
          if (well_idx_gluc[n] == w) {
              int ti = grid_idx_gluc[n];
              real G1hat = fmax(yhat[ti,2], 0.0); 
              int e = exp_id[w];
              real mu = calib_a_fixed[e] * G1hat * dilution[n] + calib_b_fixed[e];
              mu = fmax(mu, 1e-12);

              if (is_censored[n] == 1)
                ll_well[w] += lognormal_lcdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
              else
                ll_well[w] += lognormal_lpdf(lum_obs[n] | log(mu), calib_sigma_fixed[e]);
          }
      }
      if (is_train[w] == 1) log_lik_train += ll_well[w];
      else                 log_lik_holdout += ll_well[w];
  }
  log_lik = sum(ll_well);
}
