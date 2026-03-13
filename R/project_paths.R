get_model_root <- function(model_name = "gpath", model_version = "v1") {
  file.path("models", model_name, model_version)
}

get_model_r_path <- function(model_name = "gpath", model_version = "v1") {
  file.path(get_model_root(model_name, model_version), "model.R")
}

get_model_stan_path <- function(model_name = "gpath", model_version = "v1", stan_basename = NULL) {
  if (is.null(stan_basename)) {
    stan_basename <- sprintf("%s_hier.stan", model_name)
  }
  file.path(get_model_root(model_name, model_version), stan_basename)
}

get_run_output_dir <- function(
  model_name = "gpath",
  model_version = "v1",
  pipeline_name,
  dataset_label,
  run_id,
  run_label
) {
  file.path("data", "runs", model_name, model_version, pipeline_name, dataset_label, run_id, run_label)
}

default_stan_data_candidates <- function() {
  c(
    file.path("data", "inputs", "stan", "gstarvation_v1", "stan_ready_data.Rds"),
    file.path("data", "stan_ready_data.Rds"),
    file.path("ecology", "stan_ready_data.Rds")
  )
}
