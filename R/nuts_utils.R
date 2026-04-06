`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) {
    y
  } else {
    x
  }
}

nuts_read_config <- function(config_path) {
  if (!is.character(config_path) || length(config_path) != 1L || !nzchar(config_path)) {
    stop("config_path must be a single non-empty string")
  }
  if (!file.exists(config_path)) {
    stop(sprintf("Config file does not exist: %s", config_path))
  }
  if (tolower(tools::file_ext(config_path)) != "json") {
    stop("NUTS config files must be JSON")
  }

  cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  if (!is.list(cfg)) {
    stop("Parsed NUTS config must be a list")
  }

  cfg$config_path <- normalizePath(config_path, winslash = "/", mustWork = TRUE)
  cfg
}

nuts_get_config_value <- function(cfg, path, default = NULL, required = FALSE) {
  node <- cfg
  for (key in path) {
    if (!is.list(node) || is.null(node[[key]])) {
      if (required) {
        stop(sprintf("Missing required config field: %s", paste(path, collapse = ".")))
      }
      return(default)
    }
    node <- node[[key]]
  }

  if (required && (is.null(node) || !length(node))) {
    stop(sprintf("Missing required config field: %s", paste(path, collapse = ".")))
  }

  node
}

nuts_get_scalar <- function(cfg, path, default = NULL, mode = c("character", "integer", "numeric", "logical"), required = FALSE) {
  mode <- match.arg(mode)
  value <- nuts_get_config_value(cfg, path, default = default, required = required)

  if (is.null(value) || !length(value)) {
    return(value)
  }

  value <- value[[1]]
  switch(
    mode,
    character = as.character(value),
    integer = as.integer(value),
    numeric = as.numeric(value),
    logical = as.logical(value)
  )
}

nuts_load_model_api <- function(cfg) {
  source("R/project_paths.R")

  model_name <- nuts_get_scalar(cfg, c("model_name"), required = TRUE, mode = "character")
  model_version <- nuts_get_scalar(cfg, c("model_version"), default = "v1", mode = "character")
  source(get_model_r_path(model_name, model_version))

  if (!exists("get_nuts_model_api", mode = "function")) {
    stop(sprintf(
      "Model file %s does not define get_nuts_model_api()",
      get_model_r_path(model_name, model_version)
    ))
  }

  api <- get_nuts_model_api()
  if (!is.list(api)) {
    stop("get_nuts_model_api() must return a list")
  }
  if (!is.function(api$prepare_nuts_data)) {
    stop("Model NUTS API must define prepare_nuts_data()")
  }

  api$model_name <- model_name
  api$model_version <- model_version
  api$stan_file <- api$stan_file %||% get_model_stan_path(model_name, model_version)
  api
}

nuts_read_base_stan_data <- function(cfg) {
  source("R/gpath_run_utils.R")
  stan_data_path <- resolve_stan_data_path(
    nuts_get_scalar(cfg, c("stan_data_path"), default = NULL, mode = "character")
  )

  list(
    path = stan_data_path,
    data = readRDS(stan_data_path)
  )
}

nuts_default_output_dir <- function(cfg, model_api) {
  output_dir <- nuts_get_scalar(cfg, c("output", "output_dir"), default = NULL, mode = "character")
  if (!is.null(output_dir) && nzchar(trimws(output_dir))) {
    return(output_dir)
  }

  output_root <- nuts_get_scalar(cfg, c("output", "output_root"), default = file.path("data", "nuts"), mode = "character")
  file.path(output_root, model_api$model_name, model_api$model_version)
}

nuts_default_run_tag <- function(cfg) {
  run_tag <- nuts_get_scalar(cfg, c("output", "run_tag"), default = NULL, mode = "character")
  if (!is.null(run_tag) && nzchar(trimws(run_tag))) {
    return(run_tag)
  }

  chain_id <- nuts_get_scalar(cfg, c("sampling", "chain_id"), default = 1L, mode = "integer")
  sprintf("chain%02d", chain_id)
}

nuts_list_output_files <- function(output_dir, type = c("draws", "summary", "diagnostics", "metadata", "config")) {
  type <- match.arg(type)

  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("output_dir must be a single non-empty string")
  }
  if (!dir.exists(output_dir)) {
    stop(sprintf("NUTS output directory does not exist: %s", output_dir))
  }

  pattern <- switch(
    type,
    draws = "^nuts_draws_.*\\.Rds$",
    summary = "^nuts_summary_.*\\.Rds$",
    diagnostics = "^nuts_diagnostics_.*\\.Rds$",
    metadata = "^nuts_metadata_.*\\.Rds$",
    config = "^nuts_config_.*\\.json$"
  )

  files <- sort(list.files(output_dir, pattern = pattern, full.names = TRUE))
  if (!length(files)) {
    stop(sprintf("No NUTS %s files found in '%s'", type, output_dir))
  }

  files
}

nuts_list_run_tags <- function(output_dir, type = "draws") {
  files <- nuts_list_output_files(output_dir, type = type)
  base <- basename(files)
  run_tags <- sub("^nuts_[^_]+_(.*)$", "\\1", tools::file_path_sans_ext(base))
  sort(unique(sub("_chain[0-9]+$", "", run_tags)))
}

nuts_bind_run <- function(output_dir, run_tag = NULL) {
  draw_files <- nuts_list_output_files(output_dir, type = "draws")
  draw_tags <- sub(
    "^nuts_draws_(.*)$",
    "\\1",
    tools::file_path_sans_ext(basename(draw_files))
  )
  available_run_tags <- sort(unique(sub("_chain[0-9]+$", "", draw_tags)))

  if (is.null(run_tag)) {
    run_tag <- available_run_tags[[1]]
  }
  if (!(run_tag %in% available_run_tags)) {
    stop(sprintf(
      "Run tag '%s' not found in %s; available tags: %s",
      run_tag, output_dir, paste(available_run_tags, collapse = ", ")
    ))
  }

  draw_files <- draw_files[sub("_chain[0-9]+$", "", draw_tags) == run_tag]
  chain_info <- lapply(draw_files, function(path) {
    dm <- posterior::as_draws_matrix(readRDS(path))
    list(
      file = basename(path),
      chain_id = as.integer(sub("^.*_chain([0-9]+)\\.Rds$", "\\1", basename(path))),
      n_iter = nrow(dm),
      vars = colnames(dm),
      has_draw_cols = all(c(".draw", ".iteration", ".chain") %in% colnames(dm))
    )
  })

  chain_ids <- vapply(chain_info, `[[`, integer(1), "chain_id")
  if (anyNA(chain_ids)) {
    bad_files <- vapply(chain_info[is.na(chain_ids)], `[[`, character(1), "file")
    stop(sprintf(
      "Failed to parse chain ids in %s: %s",
      output_dir, paste(bad_files, collapse = ", ")
    ))
  }
  if (anyDuplicated(chain_ids)) {
    dup_ids <- sort(unique(chain_ids[duplicated(chain_ids)]))
    stop(sprintf(
      "Duplicate chain ids in %s: %s",
      output_dir, paste(dup_ids, collapse = ", ")
    ))
  }
  if (!setequal(chain_ids, seq_len(max(chain_ids)))) {
    stop(sprintf(
      "Missing chain ids in %s; observed: %s",
      output_dir, paste(sort(chain_ids), collapse = ", ")
    ))
  }

  n_iter <- vapply(chain_info, `[[`, integer(1), "n_iter")
  if (length(unique(n_iter)) != 1L) {
    stop(sprintf("Iteration-count mismatch across chains in %s", output_dir))
  }

  has_draw_cols <- vapply(chain_info, `[[`, logical(1), "has_draw_cols")
  if (length(unique(has_draw_cols)) != 1L) {
    stop(sprintf("Warmup/sampling convention mismatch across chains in %s", output_dir))
  }

  ref_vars <- chain_info[[1]]$vars
  bad_vars <- Filter(function(x) !setequal(x$vars, ref_vars), chain_info)
  if (length(bad_vars)) {
    msg <- vapply(bad_vars, function(x) {
      sprintf(
        "%s | missing: [%s] | extra: [%s]",
        x$file,
        paste(setdiff(ref_vars, x$vars), collapse = ", "),
        paste(setdiff(x$vars, ref_vars), collapse = ", ")
      )
    }, character(1))
    stop(sprintf("Variable-set mismatch in %s\n%s", output_dir, paste(msg, collapse = "\n")))
  }

  ordered_files <- draw_files[order(chain_ids)]
  list(
    draws = do.call(posterior::bind_draws, c(lapply(ordered_files, readRDS), along = "chain")),
    chain_ids = sort(chain_ids),
    run_tag = run_tag,
    available_run_tags = available_run_tags,
    output_dir = output_dir,
    draw_files = ordered_files
  )
}

nuts_summarize_draw_diagnostics <- function(draws) {
  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    sd = stats::sd,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  out <- tibble::as_tibble(out)
  dplyr::rename(out, variable = variable)
}

nuts_collect_init_sources <- function(cfg, prep) {
  sources <- c(
    nuts_get_scalar(cfg, c("init", "source_path"), default = NULL, mode = "character"),
    unlist(nuts_get_config_value(cfg, c("init", "auto_source_paths"), default = list()), use.names = FALSE),
    unlist(prep$init_context$auto_source_paths %||% list(), use.names = FALSE)
  )
  unique(sources[!is.na(sources) & nzchar(sources)])
}

nuts_build_init <- function(cfg, prep) {
  source("R/posterior_io_utils.R")

  mode <- tolower(trimws(nuts_get_scalar(cfg, c("init", "mode"), default = "random", mode = "character")))
  if (!(mode %in% c("auto", "posterior", "optim", "random"))) {
    stop(sprintf("Unsupported init mode '%s'", mode))
  }

  if (mode == "random") {
    return(list(init = 2, mode = mode, source = NULL))
  }

  init_context <- prep$init_context %||% list()
  maxG0 <- init_context$maxG0 %||% NULL
  if (is.null(maxG0) || !is.finite(maxG0)) {
    stop("Non-random init requires prep$init_context$maxG0")
  }

  chain_id <- nuts_get_scalar(cfg, c("sampling", "chain_id"), default = 1L, mode = "integer")
  seed_init <- nuts_get_scalar(cfg, c("init", "seed"), default = 3000000L + chain_id, mode = "integer")
  source_path <- nuts_get_scalar(cfg, c("init", "source_path"), default = NULL, mode = "character")
  auto_sources <- nuts_collect_init_sources(cfg, prep)

  if (mode == "optim") {
    if (is.null(source_path) || !dir.exists(source_path)) {
      stop("Optimization init requested but `init.source_path` is missing or invalid")
    }
    return(list(
      init = list(load_optim_init_from_dir(path = source_path, chain_id = chain_id, maxG0 = maxG0)),
      mode = mode,
      source = source_path
    ))
  }

  if (mode == "posterior") {
    if (is.null(source_path) || !dir.exists(source_path)) {
      stop("Posterior init requested but `init.source_path` is missing or invalid")
    }
    return(list(
      init = list(load_posterior_init_from_dir(path = source_path, seed = seed_init, maxG0 = maxG0)),
      mode = mode,
      source = source_path
    ))
  }

  for (candidate in auto_sources) {
    if (!dir.exists(candidate)) {
      next
    }

    posterior_init <- tryCatch(
      load_posterior_init_from_dir(path = candidate, seed = seed_init, maxG0 = maxG0),
      error = function(e) NULL
    )
    if (!is.null(posterior_init)) {
      return(list(init = list(posterior_init), mode = mode, source = candidate))
    }

    optim_init <- tryCatch(
      load_optim_init_from_dir(path = candidate, chain_id = chain_id, maxG0 = maxG0),
      error = function(e) NULL
    )
    if (!is.null(optim_init)) {
      return(list(init = list(optim_init), mode = mode, source = candidate))
    }
  }

  list(init = 2, mode = mode, source = NULL)
}

nuts_write_outputs <- function(fit, cfg, model_api, prep) {
  output_dir <- nuts_default_output_dir(cfg, model_api)
  run_tag <- nuts_default_run_tag(cfg)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summary_path <- file.path(output_dir, sprintf("nuts_summary_%s.Rds", run_tag))
  draws_path <- file.path(output_dir, sprintf("nuts_draws_%s.Rds", run_tag))
  diagnostics_path <- file.path(output_dir, sprintf("nuts_diagnostics_%s.Rds", run_tag))
  metadata_path <- file.path(output_dir, sprintf("nuts_metadata_%s.Rds", run_tag))
  config_copy_path <- file.path(output_dir, sprintf("nuts_config_%s.json", run_tag))

  saveRDS(fit$summary(), summary_path)
  saveRDS(fit$draws(), draws_path)
  diagnostics <- tryCatch(fit$sampler_diagnostics(), error = function(e) NULL)
  if (!is.null(diagnostics)) {
    saveRDS(diagnostics, diagnostics_path)
  }
  saveRDS(prep$metadata %||% list(), metadata_path)
  file.copy(cfg$config_path, config_copy_path, overwrite = TRUE)

  extra_paths <- list()
  if (is.function(model_api$collect_nuts_outputs)) {
    extra_outputs <- model_api$collect_nuts_outputs(
      fit = fit,
      prep = prep,
      config = cfg
    )
    if (!is.null(extra_outputs)) {
      if (!is.list(extra_outputs) || is.null(names(extra_outputs)) || any(!nzchar(names(extra_outputs)))) {
        stop("collect_nuts_outputs() must return a named list")
      }

      for (nm in names(extra_outputs)) {
        path <- file.path(output_dir, sprintf("%s_%s.Rds", nm, run_tag))
        saveRDS(extra_outputs[[nm]], path)
        extra_paths[[nm]] <- path
      }
    }
  }

  list(
    output_dir = output_dir,
    run_tag = run_tag,
    summary_path = summary_path,
    draws_path = draws_path,
    diagnostics_path = if (!is.null(diagnostics)) diagnostics_path else NULL,
    metadata_path = metadata_path,
    config_copy_path = config_copy_path,
    extra_paths = extra_paths
  )
}

# Generic JSON-configured NUTS runner.
nuts <- function(config_path) {
  library(cmdstanr)
  library(jsonlite)

  cfg <- nuts_read_config(config_path)
  model_api <- nuts_load_model_api(cfg)
  source("R/stan_design_utils.R")
  base_input <- nuts_read_base_stan_data(cfg)
  designed_input <- apply_design(
    stan_data = base_input$data,
    config = cfg
  )

  prep <- model_api$prepare_nuts_data(
    stan_data = designed_input,
    model_specific_params = nuts_get_config_value(cfg, c("model_specific_params"), default = list())
  )

  if (!is.list(prep) || is.null(prep$stan_data)) {
    stop("prepare_nuts_data() must return a list with at least `stan_data`")
  }

  init_info <- nuts_build_init(cfg, prep)
  mod <- cmdstanr::cmdstan_model(
    model_api$stan_file,
    cpp_options = list(stan_threads = TRUE)
  )

  chains <- nuts_get_scalar(cfg, c("sampling", "chains"), default = 1L, mode = "integer")
  parallel_chains <- nuts_get_scalar(cfg, c("sampling", "parallel_chains"), default = chains, mode = "integer")
  threads_per_chain <- nuts_get_scalar(cfg, c("sampling", "threads_per_chain"), default = 1L, mode = "integer")
  seed <- nuts_get_scalar(cfg, c("sampling", "seed"), default = 4000000L, mode = "integer")

  fit <- mod$sample(
    data = prep$stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    seed = seed,
    iter_warmup = nuts_get_scalar(cfg, c("sampling", "iter_warmup"), default = 500L, mode = "integer"),
    iter_sampling = nuts_get_scalar(cfg, c("sampling", "iter_sampling"), default = 1000L, mode = "integer"),
    adapt_delta = nuts_get_scalar(cfg, c("sampling", "adapt_delta"), default = 0.99, mode = "numeric"),
    max_treedepth = nuts_get_scalar(cfg, c("sampling", "max_treedepth"), default = 12L, mode = "integer"),
    save_warmup = nuts_get_scalar(cfg, c("sampling", "save_warmup"), default = TRUE, mode = "logical"),
    init = init_info$init,
    metric = nuts_get_scalar(cfg, c("sampling", "metric"), default = "dense_e", mode = "character"),
    refresh = nuts_get_scalar(cfg, c("sampling", "refresh"), default = 10L, mode = "integer")
  )

  write_info <- nuts_write_outputs(
    fit = fit,
    cfg = cfg,
    model_api = model_api,
    prep = prep
  )

  invisible(list(
    exit_code = 0L,
    config_path = cfg$config_path,
    stan_data_path = base_input$path,
    init = init_info,
    output = write_info
  ))
}
