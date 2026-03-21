selection_strategy_config <- list(
  dataset_label = "gstarvation_v1",
  include_global = TRUE,
  include_transfer = TRUE,
  transfer_fit_types = c("transfer"),
  transfer_directions = c("low_to_high", "high_to_low"),
  model_ids = NULL,
  line_ids = NULL,
  glucose_values = NULL,
  interval_hours = 48,
  total_hours = 144,
  time_step_hours = 1,
  detailed_time_step_hours = 1,
  initial_high_fraction = 0.5,
  refresh_resets_total_n = TRUE,
  workers = max(1L, parallel::detectCores(logical = TRUE) - 1L),
  representative_model_ids = NULL,
  representative_fit_contexts = c("global")
)

selection_strategy_config
