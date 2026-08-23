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

.random_hex <- function(n = 6L) {
  paste0(sample(c(letters[1:6], 0:9), n, replace = TRUE), collapse = "")
}

.sanitize_run_id <- function(x) {
  out <- gsub("[^A-Za-z0-9_.-]+", "-", as.character(x))
  out <- gsub("(^-+|-+$)", "", out)
  out[nzchar(out)]
}

.default_analysis_run_id <- function() {
  paste0(
    format(Sys.time(), "%Y%m%dT%H%M%S"),
    "-p", Sys.getpid(),
    "-",
    .random_hex(6L)
  )
}

.git_sha <- function(path_root = ".") {
  if (!nzchar(path_root)) {
    return(NA_character_)
  }

  out <- suppressWarnings(
    tryCatch(
      system2(
        "git",
        c("-C", path_root, "rev-parse", "--short", "HEAD"),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
  )

  if (length(out) == 0L || !nzchar(out[[1]])) {
    NA_character_
  } else {
    out[[1]]
  }
}

.analysis_manifest_params_for_compat <- function(params) {
  if (is.null(params) || length(params) == 0L) {
    return(list())
  }

  operational_params <- c(
    "run_simulations",
    "run_plots",
    "sim_grid_chunk_index"
  )
  keep <- sort(setdiff(names(params), operational_params))
  params[keep]
}

.analysis_validate_manifest_compat <- function(
    manifest,
    analysis_key,
    run_id,
    expected_n_chunks,
    git_sha,
    params) {
  mismatch <- character()

  if (!identical(as.character(manifest$analysis_key), as.character(analysis_key))) {
    mismatch <- c(mismatch, "analysis_key")
  }

  if (!identical(as.character(manifest$run_id), as.character(run_id))) {
    mismatch <- c(mismatch, "run_id")
  }

  manifest_chunks <- as.integer(manifest$expected_n_chunks)
  if (!identical(manifest_chunks, as.integer(expected_n_chunks))) {
    mismatch <- c(mismatch, "expected_n_chunks")
  }

  manifest_sha <- as.character(manifest$git_sha)
  current_sha <- as.character(git_sha)
  if (
    length(manifest_sha) > 0L &&
      length(current_sha) > 0L &&
      !is.na(manifest_sha) &&
      !is.na(current_sha) &&
      nzchar(manifest_sha) &&
      nzchar(current_sha) &&
      !identical(manifest_sha, current_sha)
  ) {
    mismatch <- c(mismatch, "git_sha")
  }

  manifest_params <- .analysis_manifest_params_for_compat(manifest$params)
  current_params <- .analysis_manifest_params_for_compat(params)
  if (!identical(manifest_params, current_params)) {
    mismatch <- c(mismatch, "params")
  }

  if (length(mismatch) > 0L) {
    stop(
      paste0(
        "Run ID '", run_id,
        "' is incompatible with existing manifest. Mismatched fields: ",
        paste(unique(mismatch), collapse = ", "),
        ". Use a new analysis_run_id."
      )
    )
  }

  invisible(TRUE)
}

.analysis_read_manifest <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

.analysis_current_file <- function(
    run_ctx,
    relative_path,
    required_params = list()) {
  if (
    length(relative_path) == 0L ||
      any(is.na(relative_path)) ||
      any(!nzchar(as.character(relative_path)))
  ) {
    stop("relative_path must contain at least one non-empty path component.")
  }
  if (
    length(required_params) > 0L &&
      (is.null(names(required_params)) || any(!nzchar(names(required_params))))
  ) {
    stop("required_params must be a fully named list.")
  }

  complete_path <- file.path(run_ctx$current_dir, "COMPLETE")
  if (!file.exists(complete_path)) {
    stop(
      "No complete canonical current result is available for analysis key: ",
      paste(run_ctx$analysis_key, collapse = "/"),
      "."
    )
  }

  manifest_path <- file.path(run_ctx$current_dir, "manifest.rds")
  manifest <- .analysis_read_manifest(manifest_path)
  if (is.null(manifest)) {
    stop("Canonical current result has no readable manifest.rds.")
  }

  if (!identical(
    as.character(manifest$analysis_key),
    as.character(run_ctx$analysis_key)
  )) {
    stop("Canonical current manifest does not match the requested analysis key.")
  }

  if (length(required_params) > 0L) {
    manifest_params <- manifest$params
    mismatched_params <- names(required_params)[vapply(
      names(required_params),
      function(nm) {
        is.null(manifest_params) ||
          !nm %in% names(manifest_params) ||
          !isTRUE(all.equal(
            manifest_params[[nm]],
            required_params[[nm]],
            check.attributes = FALSE
          ))
      },
      logical(1)
    )]

    if (length(mismatched_params) > 0L) {
      stop(
        "Canonical current result is incompatible with the required manifest ",
        "parameters: ",
        paste(mismatched_params, collapse = ", "),
        ". Rerun the analysis before using these results."
      )
    }
  }

  path <- do.call(
    file.path,
    c(list(run_ctx$current_dir), as.list(as.character(relative_path)))
  )
  if (!file.exists(path)) {
    stop(
      "Canonical current result is complete but the required file is missing: ",
      paste(as.character(relative_path), collapse = "/"),
      "."
    )
  }

  path
}

.analysis_find_existing_run_manifest <- function(staging_root, run_id) {
  if (!dir.exists(staging_root)) {
    return(NULL)
  }

  date_dir_vec <- list.dirs(staging_root, recursive = FALSE, full.names = TRUE)
  if (length(date_dir_vec) == 0L) {
    return(NULL)
  }

  manifest_path_vec <- file.path(date_dir_vec, run_id, "manifest.rds")
  manifest_path_vec <- manifest_path_vec[file.exists(manifest_path_vec)]

  if (length(manifest_path_vec) == 0L) {
    return(NULL)
  }
  if (length(manifest_path_vec) > 1L) {
    stop(
      "Run ID '", run_id, "' exists under multiple staging dates. ",
      "Use a new analysis_run_id."
    )
  }

  manifest_path <- manifest_path_vec[[1]]
  manifest <- .analysis_read_manifest(manifest_path)
  if (is.null(manifest)) {
    return(NULL)
  }

  list(
    manifest = manifest,
    manifest_path = manifest_path,
    staging_run_dir = dirname(manifest_path),
    run_date = basename(dirname(dirname(manifest_path)))
  )
}

.analysis_lock_path <- function(run_ctx, lock_name) {
  file.path(run_ctx$progress_run_dir, paste0(lock_name, ".lock"))
}

.analysis_acquire_lock <- function(lock_path, timeout_sec = 120) {
  dir.create(dirname(lock_path), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("filelock", quietly = TRUE)) {
    stop("Package 'filelock' is required for transactional analysis locking.")
  }
  filelock::lock(
    lock_path,
    exclusive = TRUE,
    timeout = timeout_sec * 1000
  )
}

.analysis_release_lock <- function(lock) {
  if (!is.null(lock) && inherits(lock, "filelock_lock")) {
    tryCatch(filelock::unlock(lock), error = function(e) NULL)
  }
  invisible(TRUE)
}

.analysis_run_context <- function(
    analysis_key,
    run_id = NULL,
    path_root = NULL,
    params = list(),
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L) {
  if (length(analysis_key) == 0L || !all(nzchar(analysis_key))) {
    stop("analysis_key must be a non-empty character vector.")
  }

  run_id_env <- Sys.getenv("ANALYSIS_RUN_ID", unset = "")
  run_id <- if (!is.null(run_id) && nzchar(as.character(run_id))) {
    as.character(run_id)
  } else if (nzchar(run_id_env)) {
    run_id_env
  } else {
    .default_analysis_run_id()
  }

  run_id <- .sanitize_run_id(run_id)
  if (length(run_id) == 0L || !nzchar(run_id[[1]])) {
    stop("run_id cannot be empty after sanitization.")
  }
  run_id <- run_id[[1]]

  start_time <- Sys.time()
  run_date_now <- format(start_time, "%Y-%m-%d")
  run_time <- format(start_time, "%H%M%S")

  if (requireNamespace("projr", quietly = TRUE)) {
    sim_root <- do.call(projr::projr_path_get_dir, c(list("cache"), as.list(analysis_key)))
    log_root <- do.call(
      projr::projr_path_get_dir,
      c(list("cache", "log", "analysis"), as.list(analysis_key))
    )
  } else {
    root_local <- normalizePath(
      if (is.null(path_root) || !nzchar(path_root)) "." else path_root,
      mustWork = FALSE
    )
    sim_root <- do.call(file.path, c(list(root_local, "cache"), as.list(analysis_key)))
    log_root <- do.call(
      file.path,
      c(list(root_local, "cache", "log", "analysis"), as.list(analysis_key))
    )
  }

  staging_root <- file.path(sim_root, "staging")
  current_dir <- file.path(sim_root, "current")
  existing_run <- .analysis_find_existing_run_manifest(staging_root, run_id)
  if (is.null(existing_run)) {
    run_date <- run_date_now
    staging_run_dir <- file.path(staging_root, run_date, run_id)
    progress_run_dir <- file.path(log_root, run_date, run_id)
  } else {
    run_date <- existing_run$run_date
    staging_run_dir <- existing_run$staging_run_dir

    progress_run_dir <- if (
      !is.null(existing_run$manifest$path_log_run) &&
        nzchar(as.character(existing_run$manifest$path_log_run))
    ) {
      as.character(existing_run$manifest$path_log_run)
    } else {
      file.path(log_root, run_date, run_id)
    }
  }

  chunk_label <- .sim_chunk_label(sim_grid_chunk_index, sim_grid_n_chunks)
  chunk_dir <- file.path(staging_run_dir, "chunks", chunk_label)
  chunk_output_dir <- file.path(chunk_dir, "output")
  chunk_jobs_dir <- file.path(progress_run_dir, "jobs", chunk_label)
  progress_file <- file.path(progress_run_dir, "progress.txt")

  dir.create(staging_run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(staging_run_dir, "output"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(staging_run_dir, "collated"), recursive = TRUE, showWarnings = FALSE)
  dir.create(chunk_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(chunk_jobs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(progress_run_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_path <- file.path(staging_run_dir, "manifest.rds")
  status_path <- file.path(progress_run_dir, "status.rds")
  chunk_status_dir <- file.path(progress_run_dir, "chunks")
  dir.create(chunk_status_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(manifest_path)) {
    manifest_git_sha <- .git_sha(if (is.null(path_root)) "." else path_root)
    manifest <- list(
      analysis_key = analysis_key,
      run_id = run_id,
      run_date = run_date,
      run_time = run_time,
      started_at = start_time,
      status = "running",
      git_sha = manifest_git_sha,
      params = params,
      expected_n_chunks = as.integer(sim_grid_n_chunks),
      path_staging_run = staging_run_dir,
      path_current = current_dir,
      path_log_run = progress_run_dir
    )
    .write_rds_atomic(manifest, manifest_path)
    .write_rds_atomic(manifest, file.path(progress_run_dir, "manifest.rds"))
  } else {
    manifest <- readRDS(manifest_path)
    .analysis_validate_manifest_compat(
      manifest = manifest,
      analysis_key = analysis_key,
      run_id = run_id,
      expected_n_chunks = sim_grid_n_chunks,
      git_sha = .git_sha(if (is.null(path_root)) "." else path_root),
      params = params
    )
    .write_rds_atomic(manifest, file.path(progress_run_dir, "manifest.rds"))
  }

  if (!file.exists(status_path)) {
    status <- list(
      run_id = run_id,
      analysis_key = analysis_key,
      started_at = start_time,
      ended_at = as.POSIXct(NA),
      status = "running",
      expected_n_chunks = as.integer(sim_grid_n_chunks),
      chunks = list(),
      n_completed = 0L,
      n_failed = 0L,
      n_outstanding = NA_integer_,
      promotion_done = FALSE
    )
    .write_rds_atomic(status, status_path)
  }

  list(
    analysis_key = analysis_key,
    run_id = run_id,
    run_date = run_date,
    run_time = run_time,
    sim_root = sim_root,
    log_root = log_root,
    current_dir = current_dir,
    staging_root = staging_root,
    staging_run_dir = staging_run_dir,
    staging_output_dir = file.path(staging_run_dir, "output"),
    staging_collated_dir = file.path(staging_run_dir, "collated"),
    chunk_label = chunk_label,
    chunk_dir = chunk_dir,
    chunk_output_dir = chunk_output_dir,
    progress_run_dir = progress_run_dir,
    progress_file = progress_file,
    chunk_jobs_dir = chunk_jobs_dir,
    chunk_status_dir = chunk_status_dir,
    manifest_path = manifest_path,
    status_path = status_path
  )
}

.analysis_read_status <- function(run_ctx) {
  if (!file.exists(run_ctx$status_path)) {
    return(NULL)
  }
  readRDS(run_ctx$status_path)
}

.analysis_chunk_status_path <- function(run_ctx, chunk_label = run_ctx$chunk_label) {
  file.path(run_ctx$chunk_status_dir, paste0(chunk_label, ".rds"))
}

.analysis_read_chunk_statuses <- function(run_ctx) {
  if (!dir.exists(run_ctx$chunk_status_dir)) {
    return(list())
  }

  path_vec <- list.files(
    run_ctx$chunk_status_dir,
    pattern = "[.]rds$",
    full.names = TRUE
  )

  if (length(path_vec) == 0L) {
    return(list())
  }

  out <- lapply(path_vec, function(path) {
    tryCatch(readRDS(path), error = function(e) NULL)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(list())
  }

  nm <- vapply(
    out,
    function(x) {
      if (!is.null(x$chunk_label) && nzchar(as.character(x$chunk_label))) {
        as.character(x$chunk_label)
      } else {
        NA_character_
      }
    },
    character(1)
  )
  keep <- !is.na(nm) & nzchar(nm)
  out <- out[keep]
  nm <- nm[keep]
  if (length(out) == 0L) {
    return(list())
  }
  names(out) <- nm
  out
}

.analysis_output_has_error <- function(existing_output, error_col = "error_message") {
  if (is.null(existing_output) || !error_col %in% names(existing_output)) {
    return(FALSE)
  }
  any(
    !is.na(existing_output[[error_col]]) &
      nzchar(as.character(existing_output[[error_col]]))
  )
}

.analysis_reconcile_resume_markers <- function(
    existing_output,
    file_completed,
    file_error,
    file_running = NULL,
    error_col = "error_message") {
  existing_is_error <- .analysis_output_has_error(
    existing_output = existing_output,
    error_col = error_col
  )

  if (isTRUE(existing_is_error)) {
    file.create(file_error)
    if (file.exists(file_completed)) {
      file.remove(file_completed)
    }
  } else {
    file.create(file_completed)
    if (file.exists(file_error)) {
      file.remove(file_error)
    }
  }

  if (!is.null(file_running) && file.exists(file_running)) {
    file.remove(file_running)
  }

  invisible(existing_is_error)
}

.analysis_update_status <- function(
    run_ctx,
    status = NULL,
    chunk_status = NULL,
    n_completed = NULL,
    n_failed = NULL,
    n_outstanding = NULL,
    ended_at = NULL,
    promotion_done = NULL) {
  lock_path <- .analysis_lock_path(run_ctx, "status-update")
  lock <- .analysis_acquire_lock(lock_path, timeout_sec = 300)
  if (is.null(lock)) {
    stop("Timed out acquiring status update lock for run_id: ", run_ctx$run_id)
  }
  on.exit(.analysis_release_lock(lock), add = TRUE)

  st <- .analysis_read_status(run_ctx)
  if (is.null(st)) {
    st <- list(
      run_id = run_ctx$run_id,
      analysis_key = run_ctx$analysis_key,
      started_at = Sys.time(),
      ended_at = as.POSIXct(NA),
      status = "running",
      expected_n_chunks = 1L,
      chunks = list(),
      promotion_done = FALSE
    )
  }

  if (!is.null(chunk_status) && !is.null(chunk_status$chunk_label)) {
    .write_rds_atomic(
      chunk_status,
      .analysis_chunk_status_path(run_ctx, chunk_status$chunk_label)
    )
  }

  chunk_list <- .analysis_read_chunk_statuses(run_ctx)
  st$chunks <- chunk_list

  as_nonneg_int <- function(x) {
    out <- suppressWarnings(as.integer(x))
    if (length(out) == 0L || is.na(out)) {
      return(0L)
    }
    max(0L, out)
  }

  if (length(chunk_list) > 0L) {
    st$n_completed <- sum(vapply(chunk_list, function(cs) {
      as_nonneg_int(cs$completed_sims)
    }, integer(1)))
    st$n_failed <- sum(vapply(chunk_list, function(cs) {
      as_nonneg_int(cs$failed_sims)
    }, integer(1)))
    st$n_outstanding <- sum(vapply(chunk_list, function(cs) {
      total <- as_nonneg_int(cs$total_sims)
      done <- as_nonneg_int(cs$completed_sims) + as_nonneg_int(cs$failed_sims)
      as.integer(max(0L, total - done))
    }, integer(1)))
  } else {
    if (!is.null(n_completed)) st$n_completed <- as.integer(n_completed)
    if (!is.null(n_failed)) st$n_failed <- as.integer(n_failed)
    if (!is.null(n_outstanding)) st$n_outstanding <- as.integer(n_outstanding)
  }

  if (!is.null(status)) {
    if (!identical(st$status, "completed") || identical(status, "completed")) {
      st$status <- status
    }
  }
  if (!is.null(ended_at)) st$ended_at <- ended_at
  if (!is.null(promotion_done)) st$promotion_done <- isTRUE(promotion_done)

  .write_rds_atomic(st, run_ctx$status_path)
  invisible(st)
}

.analysis_mark_chunk <- function(
    run_ctx,
    total_sims,
    completed_sims,
    failed_sims = 0L,
    collate_ok = NULL,
    validation_ok = NULL,
    error_message = NULL) {
  total_sims <- as.integer(total_sims)
  completed_sims <- as.integer(completed_sims)
  failed_sims <- as.integer(failed_sims)

  total_sims <- max(0L, total_sims)
  completed_sims <- max(0L, completed_sims)
  failed_sims <- max(0L, failed_sims)

  done_sims <- as.integer(completed_sims + failed_sims)
  outstanding <- as.integer(max(0L, total_sims - done_sims))

  cs_path <- .analysis_chunk_status_path(run_ctx)
  existing_cs <- if (file.exists(cs_path)) {
    tryCatch(
      readRDS(cs_path),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  if (is.null(collate_ok) && !is.null(existing_cs$collate_ok)) {
    collate_ok <- existing_cs$collate_ok
  }
  if (is.null(validation_ok) && !is.null(existing_cs$validation_ok)) {
    validation_ok <- existing_cs$validation_ok
  }

  chunk_status <- list(
    chunk_label = run_ctx$chunk_label,
    total_sims = total_sims,
    completed_sims = completed_sims,
    failed_sims = failed_sims,
    collate_ok = collate_ok,
    validation_ok = validation_ok,
    error_message = error_message,
    updated_at = Sys.time()
  )

  chunk_status$is_complete <- isTRUE(done_sims >= total_sims)

  .analysis_update_status(
    run_ctx = run_ctx,
    chunk_status = chunk_status,
    n_completed = completed_sims,
    n_failed = failed_sims,
    n_outstanding = outstanding,
    status = if (
      (!is.null(error_message) && nzchar(error_message)) ||
        failed_sims > 0L
    ) {
      "failed"
    } else {
      "running"
    }
  )
}

.analysis_can_promote <- function(run_ctx) {
  manifest <- .analysis_read_manifest(run_ctx$manifest_path)
  if (is.null(manifest)) {
    return(FALSE)
  }

  expected <- as.integer(manifest$expected_n_chunks)
  if (is.na(expected) || expected < 1L) {
    expected <- 1L
  }

  chunk_list <- .analysis_read_chunk_statuses(run_ctx)
  chunk_labels_expected <- vapply(
    seq_len(expected),
    function(i) .sim_chunk_label(i, expected),
    character(1)
  )

  if (!all(chunk_labels_expected %in% names(chunk_list))) {
    return(FALSE)
  }

  for (label in chunk_labels_expected) {
    cs <- chunk_list[[label]]
    if (!isTRUE(cs$is_complete)) {
      return(FALSE)
    }
    if (!isTRUE(cs$collate_ok) || !isTRUE(cs$validation_ok)) {
      return(FALSE)
    }
    if (!is.null(cs$failed_sims) && as.integer(cs$failed_sims) > 0L) {
      return(FALSE)
    }
    if (!is.null(cs$error_message) && nzchar(cs$error_message)) {
      return(FALSE)
    }
  }

  TRUE
}

.copy_dir_contents <- function(from, to) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(from, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(files) == 0L) {
    return(invisible(TRUE))
  }
  ok <- file.copy(files, to = to, recursive = TRUE, copy.mode = TRUE)
  invisible(all(ok))
}

.analysis_promote_run <- function(run_ctx) {
  lock_path <- .analysis_lock_path(run_ctx, "promotion")
  lock <- .analysis_acquire_lock(lock_path, timeout_sec = 300)
  if (is.null(lock)) {
    return(invisible(FALSE))
  }
  on.exit(.analysis_release_lock(lock), add = TRUE)

  if (!.analysis_can_promote(run_ctx)) {
    return(invisible(FALSE))
  }

  current_manifest <- .analysis_read_manifest(file.path(run_ctx$current_dir, "manifest.rds"))
  if (
    !is.null(current_manifest) &&
      identical(as.character(current_manifest$run_id), as.character(run_ctx$run_id))
  ) {
    .analysis_update_status(
      run_ctx = run_ctx,
      status = "completed",
      ended_at = Sys.time(),
      promotion_done = TRUE
    )
    return(invisible(TRUE))
  }

  complete_path <- file.path(run_ctx$staging_run_dir, "COMPLETE")
  cat(
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = complete_path
  )

  current_parent <- dirname(run_ctx$current_dir)
  dir.create(current_parent, recursive = TRUE, showWarnings = FALSE)

  tmp_current <- file.path(
    current_parent,
    paste0(".current-next-", run_ctx$run_id, "-", .random_hex(4L))
  )
  backup_current <- file.path(
    current_parent,
    paste0(".current-prev-", run_ctx$run_id, "-", .random_hex(4L))
  )

  unlink(tmp_current, recursive = TRUE, force = TRUE)
  unlink(backup_current, recursive = TRUE, force = TRUE)

  copied <- .copy_dir_contents(run_ctx$staging_run_dir, tmp_current)
  if (!isTRUE(copied)) {
    stop("Failed to stage promoted output set for current run.")
  }

  had_current <- dir.exists(run_ctx$current_dir)
  if (had_current) {
    if (!file.rename(run_ctx$current_dir, backup_current)) {
      stop("Failed to move existing current directory before promotion.")
    }
  }

  promoted <- file.rename(tmp_current, run_ctx$current_dir)
  if (!isTRUE(promoted)) {
    if (had_current && dir.exists(backup_current)) {
      file.rename(backup_current, run_ctx$current_dir)
    }
    stop("Failed to promote staged run into canonical current location.")
  }

  if (dir.exists(backup_current)) {
    unlink(backup_current, recursive = TRUE, force = TRUE)
  }

  .analysis_update_status(
    run_ctx = run_ctx,
    status = "completed",
    ended_at = Sys.time(),
    promotion_done = TRUE
  )

  invisible(TRUE)
}
