get_line_ploidy_values <- function(stan_data, line_id, tol = 1e-12) {
  line_vec <- as.integer(stan_data$line_id)
  ploidy_vec <- as.numeric(stan_data$ploidy_metric)

  idx_line <- which(line_vec == as.integer(line_id))
  if (!length(idx_line)) {
    stop(sprintf("No wells found for line_id == %d", as.integer(line_id)))
  }

  values <- sort(unique(ploidy_vec[idx_line]))
  values <- values[is.finite(values)]

  if (length(values) != 2L) {
    stop(sprintf(
      "Holdout design requires exactly 2 distinct ploidy values for line_id == %d; found %d",
      as.integer(line_id),
      length(values)
    ))
  }

  list(
    low = values[1],
    high = values[2],
    tol = tol
  )
}

apply_full_data_design <- function(stan_data) {
  stan_data$is_train <- rep(1L, stan_data$N_wells)
  attr(stan_data, "design") <- list(type = "full_data")
  stan_data
}

resolve_holdout_ploidy <- function(stan_data, line_id, ploidy, tol = 1e-12) {
  if (!is.character(ploidy) || length(ploidy) != 1L || is.na(ploidy)) {
    stop("ploidy must be a single string: 'hi' or 'lo'")
  }

  key <- tolower(trimws(ploidy))
  values <- get_line_ploidy_values(stan_data = stan_data, line_id = line_id, tol = tol)

  if (key %in% c("lo", "low")) {
    return(values$low)
  }
  if (key %in% c("hi", "high")) {
    return(values$high)
  }

  stop(sprintf("Unsupported ploidy '%s'; use 'hi' or 'lo'", ploidy))
}

apply_holdout <- function(stan_data, line_id, ploidy, tol = 1e-12) {
  holdout_value <- resolve_holdout_ploidy(
    stan_data = stan_data,
    line_id = line_id,
    ploidy = ploidy,
    tol = tol
  )

  line_vec <- as.integer(stan_data$line_id)
  ploidy_vec <- as.numeric(stan_data$ploidy_metric)

  holdout_wells <- which(
    line_vec == as.integer(line_id) &
      abs(ploidy_vec - holdout_value) < tol
  )

  if (!length(holdout_wells)) {
    stop(sprintf(
      "No wells matched holdout for line_id == %d and ploidy == %s",
      as.integer(line_id),
      ploidy
    ))
  }

  stan_data$is_train <- rep(1L, stan_data$N_wells)
  stan_data$is_train[holdout_wells] <- 0L

  attr(stan_data, "holdout_design") <- list(
    line_id = as.integer(line_id),
    ploidy = tolower(trimws(ploidy)),
    holdout_value = holdout_value,
    holdout_wells = holdout_wells
  )

  stan_data
}

apply_design <- function(stan_data, config = list()) {
  stan_data <- apply_full_data_design(stan_data)

  holdout_cfg <- NULL
  if (is.list(config) && is.list(config$design) && !is.null(config$design$holdout)) {
    holdout_cfg <- config$design$holdout
  }
  if (is.null(holdout_cfg)) {
    return(stan_data)
  }

  apply_holdout(
    stan_data = stan_data,
    line_id = holdout_cfg$line_id,
    ploidy = holdout_cfg$ploidy
  )
}

get_line_state_table <- function(stan_data, line_ids = NULL, low_label = "low_ploidy", high_label = "high_ploidy", tol = 1e-12) {
  if (is.null(line_ids)) {
    line_ids <- sort(unique(as.integer(stan_data$line_id)))
  }
  line_ids <- sort(unique(as.integer(line_ids)))

  rows <- lapply(line_ids, function(line_id) {
    states <- get_line_ploidy_values(stan_data = stan_data, line_id = line_id, tol = tol)
    tibble::tibble(
      line_id = as.integer(line_id),
      ploidy_label = c(low_label, high_label),
      ploidy_metric = c(states$low, states$high)
    )
  })

  dplyr::bind_rows(rows)
}
