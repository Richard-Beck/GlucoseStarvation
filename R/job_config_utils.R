`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) {
    y
  } else {
    x
  }
}

runner_read_json <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("JSON file does not exist: %s", path))
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

runner_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    path,
    useBytes = TRUE
  )
}

runner_merge_lists <- function(base, override) {
  if (is.null(override)) {
    return(base)
  }
  if (!is.list(base) || !is.list(override)) {
    return(override)
  }

  out <- base
  for (nm in names(override)) {
    out[[nm]] <- runner_merge_lists(base[[nm]], override[[nm]])
  }
  out
}

runner_read_run_ids <- function(cfg) {
  run_cfg <- cfg$run_ids %||% stop("Job config must define `run_ids`")
  source_type <- tolower(trimws(as.character(run_cfg$source %||% "file")))

  if (source_type != "file") {
    stop(sprintf("Unsupported run_ids source '%s'", source_type))
  }

  path <- as.character(run_cfg$path %||% stop("run_ids.path is required"))
  run_ids <- readLines(path, warn = FALSE)
  run_ids <- trimws(run_ids)
  run_ids[nzchar(run_ids)]
}

runner_get_chain_ids <- function(cfg) {
  chain_ids <- cfg$expand$chain_ids %||% list(1L)
  as.integer(unlist(chain_ids, use.names = FALSE))
}

get_job_output_root <- function(job_cfg) {
  output_root <- job_cfg$shared$output$output_root %||% NULL
  output_root <- as.character(output_root)
  if (is.null(output_root) || !length(output_root) || !nzchar(output_root)) {
    stop("Job config must define shared.output.output_root")
  }
  output_root[[1]]
}

get_job_run_ids <- function(job_cfg) {
  runner_read_run_ids(job_cfg)
}

get_job_output_dir <- function(job_cfg, run_id) {
  if (!is.character(run_id) || length(run_id) != 1L || !nzchar(run_id)) {
    stop("run_id must be a single non-empty string")
  }
  file.path(get_job_output_root(job_cfg), run_id)
}

read_job_info <- function(job_config_path) {
  job_cfg <- runner_read_json(job_config_path)
  list(
    job_config_path = job_config_path,
    job_name = as.character(job_cfg$job_name %||% ""),
    output_root = get_job_output_root(job_cfg),
    run_ids = get_job_run_ids(job_cfg),
    job_cfg = job_cfg
  )
}

runner_apply_template <- function(template, shared) {
  runner_merge_lists(template, shared)
}

runner_render_string <- function(template, values) {
  out <- template
  for (nm in names(values)) {
    out <- gsub(sprintf("\\{%s\\}", nm), as.character(values[[nm]]), out, fixed = FALSE)
  }
  out
}

runner_resolve_init_source <- function(job_cfg, run_id) {
  source("R/project_paths.R")

  init_cfg <- job_cfg$init_source %||% list()
  mode <- tolower(trimws(as.character(init_cfg$mode %||% "")))

  if (!nzchar(mode)) {
    return(NULL)
  }

  if (mode == "optim_run_dir") {
    root_override <- as.character(init_cfg$root_override %||% "")
    if (nzchar(trimws(root_override))) {
      return(file.path(root_override, run_id))
    }

    dataset_label <- as.character(init_cfg$dataset_label %||% stop("init_source.dataset_label is required"))
    run_label <- as.character(init_cfg$run_label %||% stop("init_source.run_label is required"))
    model_name <- as.character(job_cfg$shared$model_name %||% stop("shared.model_name is required"))
    model_version <- as.character(job_cfg$shared$model_version %||% "v1")

    return(get_run_output_dir(
      model_name = model_name,
      model_version = model_version,
      pipeline_name = "optim",
      dataset_label = dataset_label,
      run_id = run_id,
      run_label = run_label
    ))
  }

  stop(sprintf("Unsupported init_source.mode '%s'", mode))
}

runner_get_seed_value <- function(job_cfg, which = c("init", "sampling"), task_id, chain_id) {
  which <- match.arg(which)
  seeds_cfg <- job_cfg$seeds %||% list()
  base <- switch(
    which,
    init = as.integer(seeds_cfg$init_base %||% 3000000L),
    sampling = as.integer(seeds_cfg$sampling_base %||% 4000000L)
  )
  stride <- as.integer(seeds_cfg$task_stride %||% 1000L)

  base + stride * (as.integer(task_id) - 1L) + as.integer(chain_id)
}

runner_build_template_values <- function(job_cfg, cfg, run_id, chain_id, task_id) {
  list(
    job_name = as.character(job_cfg$job_name %||% ""),
    model_name = as.character(cfg$model_name %||% ""),
    model_version = as.character(cfg$model_version %||% ""),
    output_root = as.character(cfg$output$output_root %||% "data/nuts"),
    run_id = as.character(run_id),
    chain_id = as.integer(chain_id),
    chain_id_padded = sprintf("%02d", as.integer(chain_id)),
    task_id = as.integer(task_id),
    task_id_padded = sprintf("%03d", as.integer(task_id))
  )
}

runner_build_task_config <- function(job_cfg, template_cfg, run_id, chain_id, task_id) {
  source("R/gpath_run_utils.R")

  cfg <- runner_apply_template(template_cfg, job_cfg$shared %||% list())
  naming <- job_cfg$naming %||% list()
  init_source <- runner_resolve_init_source(job_cfg, run_id)
  template_values <- runner_build_template_values(
    job_cfg = job_cfg,
    cfg = cfg,
    run_id = run_id,
    chain_id = chain_id,
    task_id = task_id
  )

  output_dir <- runner_render_string(
    as.character(naming$output_path_template %||% "{output_root}/{run_id}"),
    template_values
  )
  run_tag <- runner_render_string(
    as.character(naming$run_tag_template %||% "{run_id}_chain{chain_id_padded}"),
    template_values
  )

  cfg$stan_data_path <- resolve_stan_data_path(as.character(cfg$stan_data_path))
  cfg$init$seed <- runner_get_seed_value(job_cfg, "init", task_id = task_id, chain_id = chain_id)
  if (!is.null(init_source) && nzchar(init_source)) {
    cfg$init$source_path <- init_source
  }

  cfg$sampling$chain_id <- as.integer(chain_id)
  cfg$sampling$seed <- runner_get_seed_value(job_cfg, "sampling", task_id = task_id, chain_id = chain_id)

  cfg$output$output_dir <- output_dir
  cfg$output$run_tag <- run_tag

  cfg$model_specific_params$run_id <- run_id
  cfg$model_specific_params$posterior_dir <- output_dir
  if (!is.null(init_source) && nzchar(init_source)) {
    cfg$model_specific_params$optim_init_dir <- init_source
  }
  cfg$model_specific_params$auto_init_source_paths <- as.list(unique(c(
    if (!is.null(init_source) && nzchar(init_source)) init_source,
    output_dir
  )))

  cfg$metadata$job_name <- as.character(job_cfg$job_name %||% "")
  cfg$metadata$task_id <- as.integer(task_id)
  cfg$metadata$run_id <- as.character(run_id)
  cfg$metadata$chain_id <- as.integer(chain_id)

  cfg
}

prepare_jobs <- function(job_config_path, output_dir = NULL) {
  library(jsonlite)

  job_cfg <- runner_read_json(job_config_path)
  template_path <- as.character(job_cfg$template_path %||% stop("Job config must define `template_path`"))
  template_cfg <- runner_read_json(template_path)

  out_dir <- output_dir %||% file.path("data", "prepared_jobs", tools::file_path_sans_ext(basename(job_config_path)))
  materialized_cfg <- job_cfg$materialized %||% list()
  config_dir <- file.path(out_dir, as.character(materialized_cfg$config_dir %||% "configs"))
  manifest_name <- as.character(materialized_cfg$manifest_name %||% "config_manifest.tsv")

  run_ids <- runner_read_run_ids(job_cfg)
  chain_ids <- runner_get_chain_ids(job_cfg)

  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- vector("list", length(run_ids) * length(chain_ids))
  idx <- 1L
  for (run_id in run_ids) {
    for (chain_id in chain_ids) {
      cfg <- runner_build_task_config(job_cfg, template_cfg, run_id, chain_id, task_id = idx)
      config_name <- sprintf("%s__%s.json", run_id, cfg$output$run_tag)
      config_path <- file.path(config_dir, config_name)
      cfg$metadata$config_path_hint <- config_path
      runner_write_json(cfg, config_path)

      rows[[idx]] <- data.frame(
        task_id = idx,
        run_id = run_id,
        chain_id = as.integer(chain_id),
        config_path = config_path,
        output_dir = cfg$output$output_dir,
        run_tag = cfg$output$run_tag,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  manifest <- do.call(rbind, rows)
  manifest_path <- file.path(out_dir, manifest_name)
  write.table(manifest, file = manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(list(
    job_config_path = job_config_path,
    output_dir = out_dir,
    config_dir = config_dir,
    manifest_path = manifest_path,
    manifest = manifest
  ))
}

run_job <- function(config_path) {
  source("R/nuts_utils.R")
  nuts(config_path)
}
