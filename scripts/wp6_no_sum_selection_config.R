selection_strategy_config <- list(
  dataset_label = "gstarvation_v1_no_sum159_fuse",
  stan_data_path = file.path("data", "inputs", "stan", "wp4b_no_sum159_fuse", "stan_ready_data.Rds"),
  include_global = TRUE,
  include_transfer = FALSE,
  model_ids = c(
    "2R_2P_0W_C0_M1",
    "2R_2P_1W_C0_M0",
    "2R_1P_0W_C0_M1",
    "1R_1P_1W_C0_M0"
  ),
  line_ids = 1:4,
  glucose_values = NULL,
  interval_hours = 48,
  total_hours = 144,
  time_step_hours = 1,
  detailed_time_step_hours = 1,
  initial_high_fraction = 0.5,
  refresh_resets_total_n = TRUE,
  workers = 4L,
  representative_model_ids = c("2R_2P_0W_C0_M1"),
  representative_fit_contexts = c("global")
)

selection_strategy_config
