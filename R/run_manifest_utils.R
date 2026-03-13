library(jsonlite)

file_md5_or_na <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

build_run_manifest <- function(
  model_name,
  model_version,
  pipeline_name,
  run_id,
  dataset_label,
  stan_data_path,
  stan_file_path,
  script_path,
  command_args,
  output_dir
) {
  git_commit <- tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) NA_character_
  )
  if (!length(git_commit)) {
    git_commit <- NA_character_
  }

  list(
    model_name = model_name,
    model_version = model_version,
    pipeline_name = pipeline_name,
    run_id = run_id,
    dataset_label = dataset_label,
    stan_data_path = stan_data_path,
    stan_data_md5 = file_md5_or_na(stan_data_path),
    stan_file_path = stan_file_path,
    stan_file_md5 = file_md5_or_na(stan_file_path),
    script_path = script_path,
    script_md5 = file_md5_or_na(script_path),
    command_args = command_args,
    output_dir = output_dir,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    host = Sys.info()[["nodename"]],
    user = Sys.info()[["user"]],
    git_commit = git_commit
  )
}

write_run_manifest <- function(output_dir, manifest) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(output_dir, "run_manifest.json")
  write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  manifest_path
}
