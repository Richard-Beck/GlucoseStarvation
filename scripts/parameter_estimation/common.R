`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

pe_fail <- function(...) stop(sprintf(...), call. = FALSE)

pe_read_config <- function(path) {
  if (!file.exists(path)) pe_fail("Release config does not exist: %s", path)
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  required <- c("release_id", "model_name", "model_version", "output_root", "counts", "datasets_tsv", "fits_tsv")
  missing <- required[vapply(required, function(nm) is.null(cfg[[nm]]), logical(1))]
  if (length(missing)) pe_fail("Release config is missing: %s", paste(missing, collapse = ", "))
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", cfg$release_id)) pe_fail("Unsafe release_id: %s", cfg$release_id)
  cfg$config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cfg
}

pe_release_root <- function(cfg) {
  file.path(
    as.character(cfg$output_root),
    sprintf("%s_%s", as.character(cfg$model_name), as.character(cfg$model_version)),
    as.character(cfg$release_id)
  )
}

pe_counts_path <- function(cfg, root = pe_release_root(cfg)) {
  counts_cfg <- cfg$counts
  mode <- tolower(as.character(counts_cfg$mode %||% "existing_rds"))
  materialize <- isTRUE(counts_cfg$materialize %||% TRUE)
  if (mode == "existing_rds" && !materialize) {
    path <- as.character(counts_cfg$input_path %||% pe_fail("counts.input_path is required"))
    if (!nzchar(path)) pe_fail("counts.input_path is empty")
    return(path)
  }
  file.path(root, "inputs", "counts", "uncorrected.Rds")
}

pe_derived_root <- function(cfg, root = pe_release_root(cfg)) {
  file.path(root, "derived")
}

pe_derived_config_path <- function(cfg) {
  path <- as.character(cfg$derived_config %||% "")
  if (!nzchar(path)) pe_fail("release config does not declare derived_config")
  path
}

pe_read_derived_config <- function(cfg) {
  path <- pe_derived_config_path(cfg)
  if (!file.exists(path)) pe_fail("Derived-data config does not exist: %s", path)
  out <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  out$config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out
}

pe_sha256 <- function(path) {
  if (!file.exists(path)) pe_fail("Cannot hash missing file: %s", path)
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) pe_fail("sha256sum failed for %s", path)
  strsplit(out[[1]], "[[:space:]]+")[[1]][[1]]
}

pe_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"),
    path,
    useBytes = TRUE
  )
}

pe_read_tsv <- function(path, required = character()) {
  if (!file.exists(path)) pe_fail("TSV does not exist: %s", path)
  x <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) pe_fail("%s is missing columns: %s", path, paste(missing, collapse = ", "))
  x
}

pe_as_flag <- function(x, label) {
  out <- suppressWarnings(as.integer(x))
  if (length(out) != 1L || is.na(out) || !out %in% c(0L, 1L)) pe_fail("%s must be 0 or 1", label)
  out
}

pe_as_positive_int <- function(x, label, allow_zero = FALSE) {
  out <- suppressWarnings(as.integer(x))
  lower <- if (allow_zero) 0L else 1L
  if (length(out) != 1L || is.na(out) || out < lower) pe_fail("%s must be an integer >= %d", label, lower)
  out
}

pe_safe_id <- function(x, label = "identifier") {
  if (length(x) != 1L || is.na(x) || !grepl("^[a-z0-9][a-z0-9_.-]*$", x)) pe_fail("Unsafe %s: %s", label, x)
  x
}

pe_copy_release_inputs <- function(cfg) {
  root <- pe_release_root(cfg)
  dir.create(file.path(root, "manifests"), recursive = TRUE, showWarnings = FALSE)
  file.copy(cfg$config_path, file.path(root, "release.json"), overwrite = TRUE)
  file.copy(cfg$datasets_tsv, file.path(root, "manifests", "dataset_spec.tsv"), overwrite = TRUE)
  file.copy(cfg$fits_tsv, file.path(root, "manifests", "fits.tsv"), overwrite = TRUE)
  if (!is.null(cfg$derived_config)) {
    derived_config <- pe_derived_config_path(cfg)
    if (!file.exists(derived_config)) pe_fail("Derived-data config does not exist: %s", derived_config)
    file.copy(derived_config, file.path(root, "manifests", "derived_config.json"), overwrite = TRUE)
  }
  invisible(root)
}

pe_git_value <- function(args) {
  out <- tryCatch(system2("git", args, stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (!length(out)) NA_character_ else paste(out, collapse = "\n")
}
