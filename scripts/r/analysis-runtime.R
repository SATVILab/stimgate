.get_qmd_param <- function(nm, default = NULL) {
  if (
    exists("params", inherits = TRUE) &&
      nm %in% names(params) &&
      !is.null(params[[nm]])
  ) {
    params[[nm]]
  } else {
    default
  }
}

.get_qmd_param_env <- function(param_nm, env_nm, default = NULL) {
  env_val <- Sys.getenv(env_nm, unset = NA_character_)

  if (!is.na(env_val) && nzchar(env_val)) {
    env_val
  } else {
    .get_qmd_param(param_nm, default)
  }
}

.as_flag <- function(x) {
  if (is.logical(x)) {
    return(isTRUE(x))
  }

  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "yes", "y", "1")
}

.validate_sim_grid_chunk_settings <- function(sim_grid_chunk_index, sim_grid_n_chunks) {
  sim_grid_chunk_index <- as.integer(sim_grid_chunk_index)
  sim_grid_n_chunks <- as.integer(sim_grid_n_chunks)

  if (
    is.na(sim_grid_chunk_index) ||
      is.na(sim_grid_n_chunks) ||
      sim_grid_n_chunks < 1L ||
      sim_grid_chunk_index < 1L ||
      sim_grid_chunk_index > sim_grid_n_chunks
  ) {
    stop(
      "Invalid sim grid chunk settings: sim_grid_chunk_index must be between 1 and sim_grid_n_chunks."
    )
  }

  invisible(list(
    sim_grid_chunk_index = sim_grid_chunk_index,
    sim_grid_n_chunks = sim_grid_n_chunks
  ))
}

.sim_chunk_label <- function(sim_grid_chunk_index, sim_grid_n_chunks) {
  sprintf(
    "%03d-of-%03d",
    as.integer(sim_grid_chunk_index),
    as.integer(sim_grid_n_chunks)
  )
}

.write_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  tmp <- paste0(
    path,
    ".tmp-",
    Sys.getpid(),
    "-",
    as.integer(Sys.time())
  )

  saveRDS(object, file = tmp)

  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp, force = TRUE)
  }

  invisible(path)
}
