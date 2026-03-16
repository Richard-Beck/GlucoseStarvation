parse_run_id <- function(run_id) {
  parts <- strsplit(run_id, "_")[[1]]
  list(
    R = as.integer(sub("R", "", parts[1])),
    P = as.integer(sub("P", "", parts[2])),
    W = as.integer(sub("W", "", parts[3])),
    C = as.integer(sub("C", "", parts[4])),
    M = as.integer(sub("M", "", parts[5]))
  )
}

get_hierarchical_k <- function(run_id, N_lines) {
  matches <- as.integer(regmatches(run_id, gregexpr("[0-9]+", run_id))[[1]])
  
  if (length(matches) < 4) {
    return(NaN)
  }
  
  R <- matches[1]
  P <- matches[2]
  W <- matches[3]
  C <- matches[4]
  
  # Base params: ae(R) + ah(R) + Y_R(R) + K(W*R) + m(1)
  base_params <- 3 * R + W * R + 1
  
  # When C=1 (strict_spec), off-diagonals are masked, meaning no effective parameters in A_mat.
  # If C=0, we estimate the full (P - 1) * R allocation matrix parameters.
  yield_alloc_params <- if (C == 0) ((P - 1) * R) else 0
  
  L_eff <- base_params + yield_alloc_params
  k_hierarchical <- L_eff * (N_lines + 1)
  
  return(k_hierarchical)
}

generate_stan_config <- function(R, P, W, strict_spec, M, base_priors) {
  # L total allocated space: ae(R) + ah(R) + Y_R(R) + A_mat((P-1)*R) + K(W*R) + m(1)
  L <- 3 * R + (P - 1) * R + W * R + 1
  centers <- numeric(L); scales <- numeric(L)
  map <- 1:L; mask <- rep(1.0, L)
  ptr <- 1
  
  centers[ptr:(ptr + R - 1)] <- base_priors$ae_c; scales[ptr:(ptr + R - 1)] <- base_priors$ae_s; ptr <- ptr + R
  centers[ptr:(ptr + R - 1)] <- base_priors$ah_c; scales[ptr:(ptr + R - 1)] <- base_priors$ah_s; ptr <- ptr + R
  centers[ptr:(ptr + R - 1)] <- base_priors$Y_R_c; scales[ptr:(ptr + R - 1)] <- base_priors$Y_R_s; ptr <- ptr + R
  
  for (c in 1:R) {
    ref_idx <- min(c, P)
    for (r in 1:P) {
      if (r == ref_idx) next # Skip the fixed 1.0 diagonal
      
      if (c == 1) {
        centers[ptr] <- base_priors$A_gap_c
        scales[ptr] <- base_priors$A_gap_s
      } else {
        centers[ptr] <- base_priors$A_c
        scales[ptr] <- base_priors$A_s
      }
      
      if (strict_spec && (r != c)) { mask[ptr] <- 0.0; map[ptr] <- 1L }
      ptr <- ptr + 1
    }
  }
  
  if (W > 0) {
    k_center <- if (M == 1) base_priors$K_add_c else base_priors$K_mult_c
    k_scale  <- if (M == 1) base_priors$K_add_s else base_priors$K_mult_s
    
    for (c in 1:R) {
      for (r in 1:W) { 
        centers[ptr] <- k_center
        scales[ptr] <- k_scale
        ptr <- ptr + 1 
      }
    }
  }
  
  centers[ptr] <- base_priors$m_c; scales[ptr] <- base_priors$m_s
  
  list(R = as.integer(R), P = as.integer(P), W = as.integer(W), L = as.integer(L), 
       prior_centers = centers, prior_scales = scales, param_map = as.integer(map), param_mask = mask)
}

smooth_min_vec <- function(v, k = 10.0) {
  -log(sum(exp(-k * v))) / k
}

rhs <- function(t, y, parms) {
  N <- max(y[1], 1e-12)
  R_vec <- pmax(y[2:(1 + parms$R)], 0.0)
  
  u <- (parms$ae * R_vec) / (parms$ah + R_vec)
  Phi <- parms$A_mat %*% (parms$Y_R * u)
  mu_growth <- smooth_min_vec(Phi)
  
  inhibition <- 1.0
  death <- 0.0
  dW <- rep(0.0, parms$W)
  
  if (parms$W > 0) {
    W_vec <- pmax(y[(2 + parms$R):(1 + parms$R + parms$W)], 0.0)
    
    inhibition <- 1.0 / (1.0 + sum((1.0 - parms$waste_mech) * W_vec))
    death <- sum(parms$waste_mech * W_vec)
    
    dW <- as.vector((parms$K_mat %*% u) * N)
  }
  
  g <- (mu_growth * inhibition) - parms$m - death
  
  list(c(g * N, -u * N, dW))
}

base_priors <- list(
  ae_c = 6.0e-5, ae_s = 0.5,
  ah_c = 1.0,    ah_s = 0.5,
  Y_R_c = 2500.0, Y_R_s = 0.5,
  A_c  = 1.0,    A_s  = 0.5,
  A_gap_c = 1.0, A_gap_s = 0.5,
  K_add_c = 0.01, K_add_s = 1.5,
  K_mult_c = .01, K_mult_s = 1.5,
  m_c  = 0.05,   m_s  = 0.2
)

get_param_names <- function(R, P, W, C, M, priors = base_priors) {
  # Extract prefixes dynamically from the priors list keys (stripping "_c")
  # This creates a named vector, e.g., base_names["ae"] yields "ae"
  keys <- grep("_c$", names(priors), value = TRUE)
  base_names <- setNames(sub("_c$", "", keys), sub("_c$", "", keys))
  
  labels <- character()
  
  # 1. Uptake (ae)
  labels <- c(labels, sprintf("%s[%d]", base_names["ae"], 1:R))
  
  # 2. Half-saturation (ah)
  labels <- c(labels, sprintf("%s[%d]", base_names["ah"], 1:R))
  
  # 3. Yield (Y_R)
  labels <- c(labels, sprintf("%s[%d]", base_names["Y_R"], 1:R))
  
  # 4. Allocation Matrix (A_mat) - skipping diagonal anchors
  for (c in 1:R) {
    ref_idx <- min(c, P)
    for (r in 1:P) {
      if (r == ref_idx) next # Skip the fixed 1.0 diagonal
      
      if (c == 1) {
        labels <- c(labels, sprintf("%s[%d,%d]", base_names["A_gap"], r, c))
      } else {
        labels <- c(labels, sprintf("%s[%d,%d]", base_names["A"], r, c))
      }
    }
  }
  
  # 5. Waste Kinetics (K_mat)
  if (W > 0) {
    k_prefix <- if (M == 1) base_names["K_add"] else base_names["K_mult"]
    for (c in 1:R) {
      for (r in 1:W) {
        labels <- c(labels, sprintf("%s[%d,%d]", k_prefix, r, c))
      }
    }
  }
  
  # 6. Maintenance/Death (m)
  labels <- c(labels, base_names["m"])
  
  return(labels)
}

normalize_ploidy_target_label <- function(label) {
  sub("\\[1\\]$", "", label)
}

build_ploidy_effect_mask <- function(
  R,
  P,
  W,
  C,
  M,
  target_spec = NULL,
  priors = base_priors
) {
  labels <- get_param_names(R, P, W, C, M, priors = priors)
  L <- length(labels)

  if (is.null(target_spec) || !length(target_spec) || all(is.na(target_spec)) || all(!nzchar(trimws(target_spec)))) {
    return(list(
      mask = rep(1.0, L),
      selected = labels,
      label = "all"
    ))
  }

  tokens <- unique(trimws(unlist(strsplit(paste(target_spec, collapse = ","), ",", fixed = TRUE))))
  tokens <- tokens[nzchar(tokens)]
  if (!length(tokens)) {
    return(list(
      mask = rep(1.0, L),
      selected = labels,
      label = "all"
    ))
  }

  token_key <- tolower(tokens)
  if (length(token_key) == 1L && token_key %in% c("all", "*")) {
    return(list(
      mask = rep(1.0, L),
      selected = labels,
      label = "all"
    ))
  }
  if (length(token_key) == 1L && token_key %in% c("none", "null")) {
    return(list(
      mask = rep(0.0, L),
      selected = character(0),
      label = "none"
    ))
  }

  norm_labels <- normalize_ploidy_target_label(labels)
  selected_idx <- integer(0)

  for (tok in tokens) {
    if (grepl("^[0-9]+$", tok)) {
      idx <- as.integer(tok)
      if (is.na(idx) || idx < 1L || idx > L) {
        stop(sprintf("Ploidy target index out of range: %s", tok))
      }
      selected_idx <- c(selected_idx, idx)
      next
    }

    exact_idx <- which(labels == tok)
    if (length(exact_idx) == 1L) {
      selected_idx <- c(selected_idx, exact_idx)
      next
    }

    alias_idx <- which(norm_labels == tok)
    if (length(alias_idx) == 1L) {
      selected_idx <- c(selected_idx, alias_idx)
      next
    }

    if (length(alias_idx) > 1L) {
      stop(sprintf("Ambiguous ploidy target '%s'; matches: %s", tok, paste(labels[alias_idx], collapse = ", ")))
    }

    stop(sprintf("Unknown ploidy target '%s'", tok))
  }

  selected_idx <- sort(unique(selected_idx))
  mask <- rep(0.0, L)
  mask[selected_idx] <- 1.0

  safe_label <- gsub("[^A-Za-z0-9]+", "-", paste(labels[selected_idx], collapse = "_"))
  safe_label <- gsub("(^-+|-+$)", "", safe_label)
  if (!nzchar(safe_label)) {
    safe_label <- "none"
  }

  list(
    mask = mask,
    selected = labels[selected_idx],
    label = safe_label
  )
}

reconstruct_parms <- function(R, P, W, strict_spec, M, 
                              base_priors, 
                              raw_theta_line, 
                              raw_theta_ploidy, 
                              ploidy_metric,
                              ploidy_effect_mask = NULL) {
  
  # 1. Generate structural maps, masks, centers, and scales
  config <- generate_stan_config(R, P, W, strict_spec, M, base_priors)
  
  prior_centers <- config$prior_centers
  prior_scales <- config$prior_scales
  param_map <- config$param_map
  param_mask <- config$param_mask
  if (is.null(ploidy_effect_mask)) {
    ploidy_effect_mask <- rep(1.0, length(raw_theta_ploidy))
  }
  
  # 2. Reconstruct physical parameters exactly as in Stan
  # Combine line and ploidy effects
  raw_w <- raw_theta_line + (raw_theta_ploidy * ploidy_effect_mask) * ploidy_metric
  
  # Transform to physical space
  theta_phys <- exp(log(prior_centers) + raw_w * prior_scales)
  
  # Apply mapping and masking logic
  theta_mapped <- theta_phys[param_map] * param_mask
  
  # 3. Unpack parameters sequentially via pointer
  ptr <- 1
  
  ae <- theta_mapped[ptr:(ptr + R - 1)]; ptr <- ptr + R
  ah <- theta_mapped[ptr:(ptr + R - 1)]; ptr <- ptr + R
  Y_R <- theta_mapped[ptr:(ptr + R - 1)]; ptr <- ptr + R
  
  # Unpack and normalize A_mat with anchored diagonals
  A_mat <- matrix(0, nrow = P, ncol = R)
  for (c in 1:R) {
    ref_idx <- min(c, P)
    raw_alloc <- numeric(P)
    
    for (r in 1:P) {
      if (r == ref_idx) {
        raw_alloc[r] <- 1.0 # Anchor the diagonal
      } else {
        if (c == 1) {
          # Symmetry breaking constraint for column 1
          raw_alloc[r] <- theta_mapped[ptr] / (1.0 + theta_mapped[ptr])
        } else {
          raw_alloc[r] <- theta_mapped[ptr]
        }
        ptr <- ptr + 1
      }
    }
    # Normalize column
    A_mat[, c] <- raw_alloc / sum(raw_alloc)
  }
  
  # Unpack K_mat and waste logic
  if (W > 0) {
    K_mat <- matrix(0, nrow = W, ncol = R)
    for (c in 1:R) {
      for (r in 1:W) {
        K_mat[r, c] <- theta_mapped[ptr]
        ptr <- ptr + 1
      }
    }
    waste_mech <- rep(as.numeric(M), W)
  } else {
    K_mat <- matrix(0, nrow = 0, ncol = 0)
    waste_mech <- numeric(0)
  }
  
  # Unpack maintenance/death
  m <- theta_mapped[ptr]
  
  # Return the list structured perfectly for the rhs() function
  list(
    R = R,
    P = P,
    W = W,
    ae = ae,
    ah = ah,
    Y_R = Y_R,
    A_mat = A_mat,
    K_mat = K_mat,
    m = m,
    waste_mech = waste_mech
  )
}

combine_parms <- function(p1, p2) {
  # Ensure the environment dimensions match
  if (p1$R != p2$R) stop("Mismatch in number of resources (R) between cell lines.")
  if (p1$W != p2$W) stop("Mismatch in number of wastes (W) between cell lines.")
  
  list(
    R = p1$R,
    W = p1$W,
    
    # Cell Line 1 specific parameters
    ae_1 = p1$ae,
    ah_1 = p1$ah,
    A_mat_1 = p1$A_mat,
    Y_R_1 = p1$Y_R,
    m_1 = p1$m,
    waste_mech_1 = p1$waste_mech,
    K_mat_1 = p1$K_mat,
    
    # Cell Line 2 specific parameters
    ae_2 = p2$ae,
    ah_2 = p2$ah,
    A_mat_2 = p2$A_mat,
    Y_R_2 = p2$Y_R,
    m_2 = p2$m,
    waste_mech_2 = p2$waste_mech,
    K_mat_2 = p2$K_mat
  )
}

rhs_mix <- function(t, y, parms) {
  # Extract biomass for both cell lines
  N1 <- max(y[1], 1e-12)
  N2 <- max(y[2], 1e-12)
  
  # Extract shared resources
  R_vec <- pmax(y[3:(2 + parms$R)], 0.0)
  
  # Cell Line 1 Uptake and Growth
  u_1 <- (parms$ae_1 * R_vec) / (parms$ah_1 + R_vec)
  Phi_1 <- parms$A_mat_1 %*% (parms$Y_R_1 * u_1)
  mu_growth_1 <- smooth_min_vec(Phi_1)
  
  # Cell Line 2 Uptake and Growth
  u_2 <- (parms$ae_2 * R_vec) / (parms$ah_2 + R_vec)
  Phi_2 <- parms$A_mat_2 %*% (parms$Y_R_2 * u_2)
  mu_growth_2 <- smooth_min_vec(Phi_2)
  
  # Default waste variables
  inhibition_1 <- 1.0
  inhibition_2 <- 1.0
  death_1 <- 0.0
  death_2 <- 0.0
  dW <- rep(0.0, parms$W)
  
  if (parms$W > 0) {
    # Extract shared wastes
    W_vec <- pmax(y[(3 + parms$R):(2 + parms$R + parms$W)], 0.0)
    
    # Cell Line 1 Waste Impact
    inhibition_1 <- 1.0 / (1.0 + sum((1.0 - parms$waste_mech_1) * W_vec))
    death_1 <- sum(parms$waste_mech_1 * W_vec)
    
    # Cell Line 2 Waste Impact
    inhibition_2 <- 1.0 / (1.0 + sum((1.0 - parms$waste_mech_2) * W_vec))
    death_2 <- sum(parms$waste_mech_2 * W_vec)
    
    # Total Waste Production (combined from both cell lines)
    dW <- as.vector((parms$K_mat_1 %*% u_1) * N1 + (parms$K_mat_2 %*% u_2) * N2)
  }
  
  # Net growth rates
  g_1 <- (mu_growth_1 * inhibition_1) - parms$m_1 - death_1
  g_2 <- (mu_growth_2 * inhibition_2) - parms$m_2 - death_2
  
  # Total Resource Depletion
  dR <- -(u_1 * N1 + u_2 * N2)
  
  # Return combined state derivatives
  list(c(g_1 * N1, g_2 * N2, dR, dW))
}

plot_trial <- function(G0, R = 2L, P = 2L, W = 0L, M = 1L) {
  library(deSolve)
  
  ae <- rep(base_priors$ae_c, R)
  ah <- rep(base_priors$ah_c, R)
  Y_R <- rep(base_priors$Y_R_c, R)
  
  # Construct and normalize A_mat (P x R) dynamically matching Stan logic
  A_mat <- matrix(0, nrow = P, ncol = R)
  for (c in 1:R) {
    ref_idx <- min(c, P)
    raw_alloc <- numeric(P)
    
    for (r in 1:P) {
      if (r == ref_idx) {
        raw_alloc[r] <- 1.0 # Anchor the diagonal
      } else {
        if (c == 1) {
          # Apply symmetry breaking constraint
          raw_alloc[r] <- base_priors$A_gap_c / (1.0 + base_priors$A_gap_c)
        } else {
          raw_alloc[r] <- base_priors$A_c
        }
      }
    }
    # Normalize column
    A_mat[, c] <- raw_alloc / sum(raw_alloc)
  }
  
  # Handle W logic dynamically
  if (W > 0) {
    K_val <- if (M == 1) base_priors$K_add_c else base_priors$K_mult_c
    K_mat <- matrix(K_val, nrow = W, ncol = R)
    waste_mech <- rep(as.numeric(M), W)
  } else {
    K_mat <- matrix(0, nrow = 0, ncol = 0)
    waste_mech <- numeric(0)
  }
  
  parms <- list(
    R = R, P = P, W = W,
    ae = ae, ah = ah,
    Y_R = Y_R,
    A_mat = A_mat,
    m = base_priors$m_c,
    waste_mech = waste_mech, 
    K_mat = K_mat
  )
  
  # Dynamic initial conditions: N (1), Resources (R), Waste (W)
  y0 <- c(N = 500, G1 = G0)
  if (R > 1) {
    y0 <- c(y0, setNames(rep(1.0, R - 1), paste0("R", 2:R)))
  }
  if (W > 0) {
    y0 <- c(y0, setNames(rep(0.0, W), paste0("W", 1:W)))
  }
  
  times <- seq(0, 130, by = 0.1)     
  
  out <- as.data.frame(ode(y = y0, times = times, func = rhs, parms = parms, method = "lsoda"))
  
  plot(out$time, out$N, xlab = "hours", ylab = "N", type = "l", 
       main = sprintf("Simulation: R=%d, P=%d, W=%d", R, P, W))
}


#R <- 1
#P <- 1
#W <- 1
#M <- 1
#plot_trial(1,R,P,W,M)
#plot_trial(25,R,P,W,M)
