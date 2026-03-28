get_line_name_from_stan_data <- function(stan_data, well_idx) {
  if (!("line_map" %in% names(stan_data)) || is.null(stan_data$line_map) || !length(stan_data$line_map)) {
    return(sprintf("line %d", as.integer(stan_data$line_id[[well_idx]])))
  }

  hit <- names(stan_data$line_map)[as.integer(unname(stan_data$line_map)) == as.integer(stan_data$line_id[[well_idx]])]
  if (!length(hit)) {
    return(sprintf("line %d", as.integer(stan_data$line_id[[well_idx]])))
  }

  hit[[1]]
}

# TODO: despite the name, this currently does more than pure well-index resolution:
# it also reads fit-specific entries from draw_vec and constructs model-specific
# simulation inputs (y0, parms). A future refactor should split this into a pure
# metadata resolver for stan_data/well_idx and a separate fit-to-simulation mapper.
resolve_well_index <- function(
  model_id,
  model_r_path,
  draw_vec,
  stan_data,
  well_idx,
  line_name = NULL,
  ploidy_effect_mask = NULL
) {
  source(model_r_path)

  well_idx <- as.integer(well_idx)[1]
  if (is.na(well_idx) || well_idx < 1L || well_idx > length(stan_data$line_id)) {
    stop("well_idx must refer to a valid row in stan_data")
  }

  if (is.null(line_name)) {
    line_name <- get_line_name_from_stan_data(stan_data, well_idx)
  }

  dims <- parse_run_id(model_id)
  R <- dims$R
  P <- dims$P
  W <- dims$W
  strict_spec <- (dims$C == 1)
  M <- dims$M
  L <- 3 * R + (P - 1) * R + W * R + 1

  line_id <- as.integer(stan_data$line_id[[well_idx]])
  ploidy_metric <- as.numeric(stan_data$ploidy_metric[[well_idx]])
  ploidy_abs <- as.numeric(stan_data$ploidy_abs[[well_idx]])
  plated_G0 <- as.numeric(stan_data$G0_per_well[[well_idx]])
  group_id <- as.integer(stan_data$group_id[[well_idx]])
  g1_id <- as.integer(stan_data$g1_id[[well_idx]])

  n0_param_name <- sprintf("raw_N0[%d]", group_id)
  g1_param_name <- sprintf("G1_0[%d]", g1_id)

  if (!n0_param_name %in% names(draw_vec)) {
    stop(sprintf("Could not find %s in draw_vec", n0_param_name))
  }
  if (!g1_param_name %in% names(draw_vec)) {
    stop(sprintf("Could not find %s in draw_vec", g1_param_name))
  }

  N0 <- exp(log(500) + as.numeric(draw_vec[[n0_param_name]]) * 1.0)
  fitted_G0 <- as.numeric(draw_vec[[g1_param_name]])

  raw_theta_line <- draw_vec[sprintf("raw_theta_line[%d,%d]", seq_len(L), line_id)]
  raw_theta_ploidy <- draw_vec[sprintf("raw_theta_ploidy[%d]", seq_len(L))]

  if (any(is.na(raw_theta_line))) {
    stop(sprintf("Missing line-level parameters for line_id=%d", line_id))
  }
  if (any(is.na(raw_theta_ploidy))) {
    stop("Missing ploidy-level parameters in draw_vec")
  }

  parms <- reconstruct_parms(
    R = R,
    P = P,
    W = W,
    strict_spec = strict_spec,
    M = M,
    base_priors = base_priors,
    raw_theta_line = raw_theta_line,
    raw_theta_ploidy = raw_theta_ploidy,
    ploidy_metric = ploidy_metric,
    ploidy_effect_mask = ploidy_effect_mask
  )

  y0 <- c(N = max(N0, 1e-6), R1 = fitted_G0)
  if (R > 1) {
    y0 <- c(y0, setNames(rep(1.0, R - 1), paste0("R", 2:R)))
  }
  if (W > 0) {
    y0 <- c(y0, setNames(rep(0.0, W), paste0("W", 1:W)))
  }

  list(
    well_idx = well_idx,
    line_id = line_id,
    line_name = line_name,
    ploidy_metric = ploidy_metric,
    ploidy = ploidy_abs,
    G0 = plated_G0,
    fitted_G0 = fitted_G0,
    y0 = y0,
    parms = parms
  )
}

simulate_well_from_index <- function(
  model_id,
  model_r_path,
  draw_vec,
  stan_data,
  well_idx,
  times = seq(0, 144, by = 0.5),
  line_name = NULL,
  ploidy_effect_mask = NULL
) {
  source(model_r_path)

  well_metadata <- resolve_well_index(
    model_id = model_id,
    model_r_path = model_r_path,
    draw_vec = draw_vec,
    stan_data = stan_data,
    well_idx = well_idx,
    line_name = line_name,
    ploidy_effect_mask = ploidy_effect_mask
  )

  out <- as.data.frame(deSolve::ode(
    y = well_metadata$y0,
    times = times,
    func = rhs,
    parms = well_metadata$parms,
    method = "lsoda"
  ))

  out_long <- tidyr::pivot_longer(out, cols = -time, names_to = "variable", values_to = "value")
  out_long$well_idx <- as.integer(well_metadata$well_idx)
  out_long$line_id <- as.integer(well_metadata$line_id)
  out_long$line_name <- well_metadata$line_name
  out_long$ploidy <- as.numeric(well_metadata$ploidy)
  out_long$G0 <- as.numeric(well_metadata$G0)
  out_long
}

simulate_line_from_draw <- function(
  model_id,
  model_r_path,
  draw_vec,
  stan_data,
  line_id,
  times = seq(0, 144, by = 0.5),
  line_name = NULL,
  ploidy_effect_mask = NULL,
  extra_cols = list()
) {
  line_id <- as.integer(line_id)[1]
  line_wells <- which(as.integer(stan_data$line_id) == line_id)
  if (!length(line_wells)) {
    stop(sprintf("No wells found for line_id == %d", line_id))
  }

  sim_df <- dplyr::bind_rows(lapply(line_wells, function(well_idx) {
    simulate_well_from_index(
      model_id = model_id,
      model_r_path = model_r_path,
      draw_vec = draw_vec,
      stan_data = stan_data,
      well_idx = well_idx,
      times = times,
      line_name = line_name,
      ploidy_effect_mask = ploidy_effect_mask
    )
  }))

  if (length(extra_cols)) {
    for (nm in names(extra_cols)) {
      sim_df[[nm]] <- extra_cols[[nm]]
    }
  }

  sim_df
}

build_obs_df_from_stan_data <- function(
  stan_data,
  line_id = NULL,
  line_name = NULL
) {
  obs_n_df <- tibble::tibble(
    well_idx = as.integer(stan_data$well_idx_count),
    time = stan_data$t_grid[stan_data$grid_idx_count],
    variable = "N",
    value = as.numeric(stan_data$N_obs)
  )

  exps_g <- stan_data$exp_id[stan_data$well_idx_gluc]
  obs_r1_df <- tibble::tibble(
    well_idx = as.integer(stan_data$well_idx_gluc),
    time = stan_data$t_grid[stan_data$grid_idx_gluc],
    variable = "R1",
    value = pmax(
      0,
      (
        as.numeric(stan_data$lum_obs) - stan_data$calib_b_fixed[exps_g]
      ) / (stan_data$calib_a_fixed[exps_g] * as.numeric(stan_data$dilution) + 1e-12)
    )
  )

  obs_df <- dplyr::bind_rows(obs_n_df, obs_r1_df) %>%
    dplyr::mutate(
      line_id = as.integer(stan_data$line_id[well_idx]),
      line_name = vapply(well_idx, function(idx) get_line_name_from_stan_data(stan_data, idx), character(1)),
      ploidy = as.numeric(stan_data$ploidy_abs[well_idx]),
      G0 = as.numeric(stan_data$G0_per_well[well_idx])
    )

  if (!is.null(line_id)) {
    line_id <- as.integer(line_id)[1]
    obs_df <- obs_df %>% dplyr::filter(line_id == !!line_id)

    if (!is.null(line_name)) {
      obs_df$line_name <- line_name
    }
  }

  obs_df
}
