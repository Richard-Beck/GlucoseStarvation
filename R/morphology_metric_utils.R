build_morphology_metric_table <- function(
  summary_path,
  pair_map_path,
  min_hours = 24,
  max_hours = 48,
  object_scope = "alive"
) {
  if (!file.exists(summary_path)) stop(sprintf("Morphology summary is missing: %s", summary_path))
  if (!file.exists(pair_map_path)) stop(sprintf("Morphology pair map is missing: %s", pair_map_path))

  source_obj <- readRDS(summary_path)
  if (!is.list(source_obj) || is.null(source_obj$replicate_time)) {
    stop("Morphology summary must contain a replicate_time table")
  }
  x <- as.data.frame(source_obj$replicate_time, stringsAsFactors = FALSE)
  required <- c(
    "cellLine", "ploidy", "glucose", "hours", "object_scope", "n_objects",
    "n_with_nucleus", "segmented_area_px_q50", "nuclear_area_detected_px_q50"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(sprintf("Morphology replicate_time table lacks: %s", paste(missing, collapse = ", ")))

  pair_map <- utils::read.delim(
    pair_map_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  map_required <- c("cellLine", "ploidy_label", "pair_state", "ploidy_abs")
  map_missing <- setdiff(map_required, names(pair_map))
  if (length(map_missing)) stop(sprintf("Morphology pair map lacks: %s", paste(map_missing, collapse = ", ")))
  if (nrow(pair_map) != 10L || anyDuplicated(paste(pair_map$cellLine, pair_map$ploidy_label))) {
    stop("Morphology pair map must contain exactly 10 unique cell-line/ploidy-label rows")
  }
  if (!all(pair_map$pair_state %in% c("baseline", "elevated"))) {
    stop("Morphology pair_state values must be baseline or elevated")
  }
  state_counts <- table(pair_map$cellLine, pair_map$pair_state)
  if (!all(state_counts == 1L)) stop("Each cell line must have exactly one baseline and one elevated pair state")

  keep <- x$object_scope == object_scope &
    is.finite(as.numeric(x$hours)) &
    as.numeric(x$hours) >= min_hours &
    as.numeric(x$hours) <= max_hours &
    x$cellLine %in% pair_map$cellLine
  x <- x[keep, , drop = FALSE]
  if (!nrow(x)) stop("No morphology rows remain after applying the configured filters")

  x$key <- paste(x$cellLine, x$ploidy, sep = "\r")
  pair_map$key <- paste(pair_map$cellLine, pair_map$ploidy_label, sep = "\r")
  unknown <- setdiff(unique(x$key), pair_map$key)
  if (length(unknown)) stop(sprintf("Unmapped morphology states: %s", paste(unknown, collapse = ", ")))

  rows <- lapply(seq_len(nrow(pair_map)), function(i) {
    y <- x[x$key == pair_map$key[[i]], , drop = FALSE]
    if (!nrow(y)) stop(sprintf("No morphology rows for %s / %s", pair_map$cellLine[[i]], pair_map$ploidy_label[[i]]))
    cell_values <- as.numeric(y$segmented_area_px_q50)
    nuclear_values <- as.numeric(y$nuclear_area_detected_px_q50)
    if (!any(is.finite(cell_values)) || !any(is.finite(nuclear_values))) {
      stop(sprintf("Non-finite morphology summaries for %s / %s", pair_map$cellLine[[i]], pair_map$ploidy_label[[i]]))
    }
    n_objects <- sum(as.numeric(y$n_objects), na.rm = TRUE)
    n_with_nucleus <- sum(as.numeric(y$n_with_nucleus), na.rm = TRUE)
    glucose_values <- suppressWarnings(as.numeric(y$glucose))
    data.frame(
      cellLine = pair_map$cellLine[[i]],
      ploidy_label = pair_map$ploidy_label[[i]],
      pair_state = pair_map$pair_state[[i]],
      ploidy_abs = as.numeric(pair_map$ploidy_abs[[i]]),
      cell_area_q50_px = stats::median(cell_values[is.finite(cell_values)]),
      nuclear_area_q50_px = stats::median(nuclear_values[is.finite(nuclear_values)]),
      n_summaries = sum(is.finite(cell_values) & is.finite(nuclear_values)),
      n_glucose_conditions = length(unique(glucose_values[is.finite(glucose_values)])),
      min_glucose = min(glucose_values, na.rm = TRUE),
      max_glucose = max(glucose_values, na.rm = TRUE),
      total_alive_objects = n_objects,
      total_alive_with_nucleus = n_with_nucleus,
      nuclear_detection_fraction = n_with_nucleus / n_objects,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  for (metric in c("cell_area_q50_px", "nuclear_area_q50_px")) {
    metric_out <- numeric(nrow(out))
    for (line in unique(out$cellLine)) {
      idx <- which(out$cellLine == line)
      base_idx <- idx[out$pair_state[idx] == "baseline"]
      base_value <- out[[metric]][base_idx]
      if (length(base_value) != 1L || !is.finite(base_value) || base_value <= 0) {
        stop(sprintf("Invalid %s baseline for %s", metric, line))
      }
      metric_out[idx] <- log2(out[[metric]][idx] / base_value)
    }
    out[[sub("_q50_px$", "_metric", metric)]] <- metric_out
  }

  out <- out[order(out$cellLine, match(out$pair_state, c("baseline", "elevated"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

apply_morphology_metric <- function(stan_data, metric_table, metric_column) {
  if (!metric_column %in% names(metric_table)) stop(sprintf("Unknown morphology metric column: %s", metric_column))
  required <- c("N_wells", "line_id", "line_map", "ploidy_abs", "ploidy_metric")
  missing <- setdiff(required, names(stan_data))
  if (length(missing)) stop(sprintf("Stan data lacks fields needed for morphology mapping: %s", paste(missing, collapse = ", ")))

  inverse_line_map <- setNames(names(stan_data$line_map), as.character(unlist(stan_data$line_map)))
  well_lines <- unname(inverse_line_map[as.character(stan_data$line_id)])
  if (anyNA(well_lines)) stop("Could not invert line_map for all Stan wells")

  mapped <- rep(NA_real_, stan_data$N_wells)
  for (i in seq_len(nrow(metric_table))) {
    idx <- well_lines == metric_table$cellLine[[i]] &
      abs(as.numeric(stan_data$ploidy_abs) - metric_table$ploidy_abs[[i]]) < 1e-8
    mapped[idx] <- as.numeric(metric_table[[metric_column]][[i]])
  }
  if (length(mapped) != stan_data$N_wells || any(!is.finite(mapped))) {
    stop(sprintf("Failed to map %s to every Stan well", metric_column))
  }

  stan_data$ploidy_metric <- mapped
  attr(stan_data, "morphology_metric") <- list(
    metric_column = metric_column,
    baseline_definition = "within-line log2 ratio to baseline pair state",
    original_field_replaced = "ploidy_metric"
  )
  stan_data
}
