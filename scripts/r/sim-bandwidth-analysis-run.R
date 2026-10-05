# Shared run helpers for the bandwidth-simulation analyses (QMDs 2-6).
#
# Source after analysis-runtime.R, sim-misc.R, sim-bandwidth.R and
# sim-bandwidth-analysis-io.R.
#
# Each analysis supplies:
# - a scenario function `scenario_fn(row, settings)` that runs one row of
#   `sim_grid` (a one-row tibble) with the fixed `settings` list and returns
#   a data frame of results. It must not seed the RNG itself.
# - optionally a `validate_fn(tbl)` returning problem strings (character(0)
#   when valid), and a `collate_fn(tbl)` returning a named list of objects
#   saved as `collated/<name>.rds` at promotion.
#
# `.simBandwidthRunRow()` is the single code path for one row: the furrr
# workers call it, and so does an interactive rerun of one `sim_id`, so both
# use exactly the same random numbers.

# RNG kinds fixed for every simulation row, independent of the caller's
# (or furrr's L'Ecuyer-CMRG) RNG state.
.simBandwidthRngKind <- c(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

#' Run one simulation-grid row with its own seed
#'
#' Seeds with `row$sim_seed` under fixed RNG kinds, runs `scenario_fn`,
#' restores the caller's RNG state, and prefixes the grid-row columns to the
#' result. The output is what the parallel run stores for that `sim_id`.
#'
#' @param row data.frame One row of `sim_grid`, with `sim_id` and `sim_seed`.
#' @param scenario_fn function `function(row, settings)` returning a data
#'   frame.
#' @param settings list Fixed analysis settings passed to `scenario_fn`.
#' @param error_col character Name of the error column (NA on success).
#' @return tibble Grid-row columns, scenario results and `error_col`.
.simBandwidthRunRow <- function(
    row,
    scenario_fn,
    settings = list(),
    error_col = "error_message") {
  if (nrow(row) != 1L || !all(c("sim_id", "sim_seed") %in% names(row))) {
    stop("row must be one sim_grid row with sim_id and sim_seed columns.")
  }
  res <- .analysis_with_seed(row$sim_seed[[1]], scenario_fn(row, settings))
  res <- tibble::as_tibble(res)
  if (nrow(res) == 0L) {
    stop("scenario_fn must return at least one result row.")
  }
  res <- res[, setdiff(names(res), c(names(row), error_col)), drop = FALSE]
  out <- dplyr::bind_cols(
    tibble::as_tibble(row)[rep(1L, nrow(res)), , drop = FALSE],
    res
  )
  out[[error_col]] <- rep(NA_character_, nrow(out))
  out
}

#' Error row for a failed simulation-grid row
#'
#' Grid-row columns plus the error message only. Result columns are absent,
#' so `dplyr::bind_rows()` / `purrr::list_rbind()` fill them with typed NA
#' from the success rows instead of clashing with hand-typed NA columns.
#'
#' @param row data.frame One row of `sim_grid`.
#' @param message character Error message.
#' @param error_col character Name of the error column.
#' @return tibble One row.
.simBandwidthErrorRow <- function(row, message, error_col = "error_message") {
  out <- tibble::as_tibble(row)
  out[[error_col]] <- as.character(message)
  out
}

#' Worker body for one simulation-grid row
#'
#' Handles resume, running/completed/error markers, atomic output writing
#' and error rows around `.simBandwidthRunRow()`. Existing outputs are reused
#' unless the output or marker recorded an error and `retry_errors` is TRUE,
#' in which case the row is run again.
#'
#' @param row data.frame One row of `sim_grid`.
#' @param scenario_fn,settings,error_col See `.simBandwidthRunRow()`.
#' @param run_ctx list Run context from `.analysis_run_context()`.
#' @param total_sims integer Number of rows in this chunk (for progress).
#' @param sim_grid_chunk_index,sim_grid_n_chunks integer Chunk settings.
#' @param retry_errors logical Rerun rows whose saved output has an error.
#' @param heading character Progress dashboard heading.
#' @param path_root character or NULL Checkout to load in workers.
#' @param p function or NULL progressr progressor.
#' @return tibble Saved or newly computed output for the row.
.simBandwidthRunRowResumable <- function(
    row,
    scenario_fn,
    settings,
    run_ctx,
    total_sims,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L,
    retry_errors = TRUE,
    error_col = "error_message",
    heading = "BANDWIDTH SIMULATION PROGRESS",
    path_root = NULL,
    p = NULL) {
  if (isTRUE(run_ctx$read_only)) {
    stop("Cannot run simulations in a read-only results context.")
  }
  if (!is.null(path_root)) {
    .simBandwidthEnsureCurrentCheckout(path_root)
  }
  sim_id <- row$sim_id[[1]]
  dir_jobs <- run_ctx$chunk_jobs_dir
  file_running <- file.path(dir_jobs, paste0("running-", sim_id))
  file_completed <- file.path(dir_jobs, paste0("completed-", sim_id))
  file_error <- file.path(dir_jobs, paste0("error-", sim_id))
  file_output <- .path_sim_output(
    sim_id,
    dir_output = run_ctx$chunk_output_dir,
    sim_grid_chunk_index = sim_grid_chunk_index,
    sim_grid_n_chunks = sim_grid_n_chunks
  )
  update_progress <- function() {
    .update_progress_summary(
      path_progress_file = run_ctx$progress_file,
      dir_jobs_chunk = dir_jobs,
      total_sims = total_sims,
      sim_grid_chunk_index = sim_grid_chunk_index,
      sim_grid_n_chunks = sim_grid_n_chunks,
      dir_output = run_ctx$chunk_output_dir,
      heading = heading
    )
  }
  tick <- function(msg) {
    update_progress()
    if (!is.null(p)) p(sprintf("%s sim_id: %s", msg, sim_id))
  }

  existing <- if (file.exists(file_output)) {
    tryCatch(readRDS(file_output), error = function(e) NULL)
  }
  if (
    !is.null(existing) &&
      !(isTRUE(retry_errors) && (
        file.exists(file_error) ||
          .analysis_output_has_error(existing, error_col)
      ))
  ) {
    .analysis_reconcile_resume_markers(
      existing_output = existing,
      file_completed = file_completed,
      file_error = file_error,
      file_running = file_running,
      error_col = error_col
    )
    tick("Skipped existing")
    return(existing)
  }

  unlink(c(file_error, file_completed))
  file.create(file_running)
  update_progress()
  out <- tryCatch(
    .simBandwidthRunRow(row, scenario_fn, settings, error_col),
    error = function(e) {
      .simBandwidthErrorRow(row, conditionMessage(e), error_col)
    }
  )
  .write_rds_atomic(out, file_output)
  failed <- .analysis_output_has_error(out, error_col)
  file.create(if (failed) file_error else file_completed)
  unlink(file_running)
  tick(if (failed) "ERROR on" else "Completed")
  out
}

#' Run this chunk's simulation-grid rows in parallel
#'
#' @param sim_grid data.frame Rows assigned to this chunk.
#' @param workers integer Number of multisession workers (1 = sequential).
#' @param ... Passed to `.simBandwidthRunRowResumable()` (`scenario_fn`,
#'   `settings`, `run_ctx`, chunk settings, `retry_errors`, `error_col`,
#'   `heading`, `path_root`).
#' @return list Per-row outputs, in `sim_grid` order.
.simBandwidthRunGrid <- function(sim_grid, ..., workers = .simGetCores()) {
  if (nrow(sim_grid) == 0L) {
    return(list())
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  workers <- max(1L, as.integer(workers))
  if (workers > 1L) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  row_list <- lapply(seq_len(nrow(sim_grid)), function(i) {
    sim_grid[i, , drop = FALSE]
  })
  progressr::with_progress({
    p <- progressr::progressor(steps = length(row_list))
    # Arguments go through furrr's `...` (not a closure) so that functions
    # such as `scenario_fn` ship to workers with their globals. Seeds come
    # from each row's sim_seed (see .simBandwidthRunRow); furrr's
    # seed = TRUE only silences RNG warnings and does not affect results.
    do.call(furrr::future_map, c(
      list(
        .x = row_list,
        .f = .simBandwidthRunRowResumable,
        total_sims = length(row_list),
        p = p
      ),
      list(...),
      list(.options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
    ))
  })
}

#' Locate the saved output file for one sim_id
#'
#' @param dir character Directory searched recursively (a staged run, a
#'   chunk directory or `current/`).
#' @param sim_id integer Simulation ID.
#' @return character Path, or `character(0)` if not found.
.simBandwidthFindSimOutput <- function(dir, sim_id) {
  paths <- .find_bw_list_output_files(
    output_dir = dir,
    allow_cache_fallback = FALSE
  )
  paths[grepl(sprintf("-sim_id_%06d[.]rds$", as.integer(sim_id)), paths)]
}

#' Read saved simulation outputs into one tibble
#'
#' @param path_vec character Output file paths.
#' @return tibble Row-bound outputs (empty tibble for no paths).
.simBandwidthReadOutputs <- function(path_vec) {
  out <- lapply(path_vec, function(path) {
    tryCatch(
      readRDS(path),
      error = function(e) {
        stop(
          "Could not read saved simulation result: ", path, ". ",
          conditionMessage(e)
        )
      }
    )
  })
  if (length(out) == 0L) {
    return(tibble::tibble())
  }
  purrr::list_rbind(out)
}

#' Validate collated outputs against the expected grid
#'
#' Checks that exactly the expected `sim_id` set is present, that no output
#' recorded an error, that `sim_seed` matches the grid, and any analysis
#' checks from `validate_fn`.
#'
#' @param tbl data.frame Collated outputs.
#' @param expected_grid data.frame Grid rows expected in `tbl`.
#' @param validate_fn function or NULL `function(tbl)` returning problems.
#' @param error_col character Name of the error column.
#' @return list `ids_ok`, `problems` (character) and `validation_ok`.
.simBandwidthValidateOutputs <- function(
    tbl,
    expected_grid,
    validate_fn = NULL,
    error_col = "error_message") {
  expected <- dplyr::distinct(
    tibble::tibble(
      sim_id = as.integer(expected_grid$sim_id),
      sim_seed = as.integer(expected_grid$sim_seed)
    )
  ) |>
    dplyr::arrange(.data$sim_id)
  observed <- if (all(c("sim_id", "sim_seed") %in% names(tbl))) {
    tibble::tibble(
      sim_id = as.integer(tbl$sim_id),
      sim_seed = as.integer(tbl$sim_seed)
    ) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$sim_id)
  } else {
    expected[0, ]
  }
  ids_ok <- identical(unique(observed$sim_id), expected$sim_id)
  problems <- character()
  if (!ids_ok) {
    problems <- "outputs do not contain exactly the expected sim_id set"
  } else if (!identical(observed, expected)) {
    problems <- "output sim_seed values do not match the simulation grid"
  }
  if (error_col %in% names(tbl)) {
    err <- !is.na(tbl[[error_col]]) & nzchar(as.character(tbl[[error_col]]))
    error_ids <- sort(unique(as.integer(tbl$sim_id[err])))
    if (length(error_ids) > 0L) {
      problems <- c(problems, paste0(
        "simulation errors for sim_id: ", paste(error_ids, collapse = ", ")
      ))
    }
  }
  if (length(problems) == 0L && !is.null(validate_fn)) {
    problems <- as.character(validate_fn(tbl))
  }
  list(
    ids_ok = ids_ok,
    problems = problems,
    validation_ok = length(problems) == 0L
  )
}

#' Count completed and failed markers for this chunk
#'
#' @param run_ctx list Run context.
#' @return list `completed` and `failed` counts.
.simBandwidthChunkMarkerCounts <- function(run_ctx) {
  files <- list.files(run_ctx$chunk_jobs_dir)
  ids <- function(prefix) {
    unique(sub(prefix, "", files[grepl(prefix, files)]))
  }
  failed <- ids("^error-")
  completed <- ids("^completed-")
  list(
    completed = length(setdiff(completed, failed)),
    failed = length(failed)
  )
}

#' Promote a staged run once every chunk is complete and valid
#'
#' Under a collation lock, checks that the staged outputs cover exactly the
#' full grid, validates them, writes `collate_fn(tbl)` to `collated/` and
#' promotes the run. Returns FALSE (without error) while other chunks are
#' outstanding.
#'
#' @param run_ctx list Run context.
#' @param sim_grid_all data.frame Full grid (all chunks).
#' @param collate_fn function `function(tbl)` returning a named list.
#' @param validate_fn,error_col See `.simBandwidthValidateOutputs()`.
#' @param label character Analysis label for error messages.
#' @param counts list Chunk counts `total`, `completed`, `failed`.
#' @return logical Whether the run was promoted.
.simBandwidthPromoteIfReady <- function(
    run_ctx,
    sim_grid_all,
    collate_fn,
    validate_fn = NULL,
    error_col = "error_message",
    label = "analysis",
    counts = list(total = 0L, completed = 0L, failed = 0L)) {
  if (isTRUE(run_ctx$read_only)) {
    return(invisible(FALSE))
  }
  if (!.analysis_can_promote(run_ctx)) {
    return(invisible(FALSE))
  }
  lock <- .analysis_acquire_lock(
    .analysis_lock_path(run_ctx, "collation"),
    timeout_sec = 300
  )
  if (is.null(lock)) {
    return(invisible(FALSE))
  }
  on.exit(.analysis_release_lock(lock), add = TRUE)
  if (!.analysis_can_promote(run_ctx)) {
    return(invisible(FALSE))
  }

  paths <- .find_bw_list_output_files(
    output_dir = run_ctx$staging_run_dir,
    allow_cache_fallback = FALSE
  )
  tbl <- .simBandwidthReadOutputs(paths)
  validation <- .simBandwidthValidateOutputs(
    tbl,
    expected_grid = sim_grid_all,
    validate_fn = validate_fn,
    error_col = error_col
  )
  if (length(paths) != nrow(sim_grid_all)) {
    validation$problems <- c(
      "staged output files do not match the full simulation grid",
      validation$problems
    )
  }
  if (length(validation$problems) > 0L) {
    error_message <- paste0(
      "Refusing to promote ", label, ": ",
      paste(validation$problems, collapse = "; "), "."
    )
    .analysis_mark_chunk(
      run_ctx = run_ctx,
      total_sims = counts$total,
      completed_sims = counts$completed,
      failed_sims = counts$failed,
      collate_ok = FALSE,
      validation_ok = FALSE,
      error_message = error_message
    )
    stop(error_message)
  }

  collated <- collate_fn(tbl)
  if (
    !is.list(collated) || length(collated) == 0L ||
      is.null(names(collated)) || anyNA(names(collated)) ||
      any(!nzchar(names(collated))) || anyDuplicated(names(collated))
  ) {
    stop("collate_fn must return a non-empty list with unique object names.")
  }
  for (nm in names(collated)) {
    .write_rds_atomic(
      collated[[nm]],
      file.path(run_ctx$staging_collated_dir, paste0(nm, ".rds"))
    )
  }
  .analysis_promote_run(run_ctx)
}

#' Validate this chunk, record its status and promote if all chunks are done
#'
#' Reads this chunk's outputs, validates them against the chunk grid, marks
#' the chunk status, stops on errors and then calls
#' `.simBandwidthPromoteIfReady()`. An empty chunk is marked complete.
#' Writes no per-chunk collated files.
#'
#' @param run_ctx list Run context.
#' @param sim_grid data.frame Rows assigned to this chunk.
#' @param sim_grid_all data.frame Full grid.
#' @param ... Passed to `.simBandwidthPromoteIfReady()` (`collate_fn`,
#'   `validate_fn`, `error_col`, `label`).
#' @return logical Whether the run was promoted.
.simBandwidthFinishChunk <- function(
    run_ctx,
    sim_grid,
    sim_grid_all,
    validate_fn = NULL,
    error_col = "error_message",
    label = "analysis",
    ...) {
  if (isTRUE(run_ctx$read_only)) {
    return(invisible(FALSE))
  }
  counts <- .simBandwidthChunkMarkerCounts(run_ctx)
  counts$total <- nrow(sim_grid)
  validation <- if (nrow(sim_grid) == 0L) {
    list(ids_ok = TRUE, problems = character(), validation_ok = TRUE)
  } else {
    .simBandwidthValidateOutputs(
      .simBandwidthReadOutputs(.find_bw_list_output_files(
        output_dir = run_ctx$chunk_dir,
        allow_cache_fallback = FALSE
      )),
      expected_grid = sim_grid,
      validate_fn = validate_fn,
      error_col = error_col
    )
  }
  error_message <- if (!validation$validation_ok) {
    paste0(
      "Chunk ", run_ctx$chunk_label, " of ", label, ": ",
      paste(validation$problems, collapse = "; "), "."
    )
  }
  .analysis_mark_chunk(
    run_ctx = run_ctx,
    total_sims = counts$total,
    completed_sims = counts$completed,
    failed_sims = counts$failed,
    collate_ok = validation$ids_ok,
    validation_ok = validation$validation_ok,
    error_message = error_message
  )
  if (!is.null(error_message)) {
    stop(error_message)
  }
  .simBandwidthPromoteIfReady(
    run_ctx = run_ctx,
    sim_grid_all = sim_grid_all,
    validate_fn = validate_fn,
    error_col = error_col,
    label = label,
    counts = counts,
    ...
  )
}

# ---------------------------------------------------------------------------
# Analysis 2: fixed-bandwidth background-subtracted frequency (global)
# ---------------------------------------------------------------------------

#' Analysis 2 scenario: one fixed-bandwidth frequency simulation
#'
#' @param row data.frame One row of the analysis 2 `sim_grid`.
#' @param settings list Fixed `.simBandwidthBsFreq()` arguments (e.g.
#'   `nSample`, `nIter`, `covEvMin`, `clusterGates`).
#' @return tibble `.simBandwidthBsFreq()` output.
.simBandwidthFreqBsGlobalScenario <- function(row, settings) {
  do.call(.simBandwidthBsFreq, c(settings, list(
    biasUns = row$bias_uns[[1]],
    bw = row$bw[[1]],
    bwFallback = row$bw[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]],
    samplePerturbationSd = row$sample_perturbation_sd[[1]],
    conditionPerturbationSd = row$condition_perturbation_sd[[1]],
    clusterPerturbationSd = row$cluster_perturbation_sd[[1]],
    backgroundRelativeToResponse = row$background_relative_to_response[[1]],
    ncellUnsRelativeToStim = row$n_cell_uns_relative_to_stim[[1]]
  )))
}

#' Analysis 2 final sample-level results and scenario summary
#'
#' @param tbl data.frame Collated analysis 2 outputs.
#' @param grid_cols character Grid column names.
#' @return list `bw_tbl_results_raw` and `bw_tbl_results_summary`.
.simBandwidthFreqBsGlobalCollate <- function(tbl, grid_cols) {
  if (all(is.na(tbl$threshold))) {
    stop(
      "No valid threshold results were collated. ",
      "Check the progress log for simulation-level errors."
    )
  }
  results_raw <- tbl |>
    dplyr::filter(
      .data$method == "loc_sample",
      is.finite(.data$threshold),
      is.finite(.data$propRespTruth),
      is.finite(.data$propRespEst)
    ) |>
    dplyr::select(
      dplyr::any_of(grid_cols),
      "iter", "sample", "ind", "method",
      "propRespTruth", "propRespEst", "threshold"
    )
  if (anyDuplicated(results_raw[c("sim_id", "iter", "ind")]) > 0L) {
    stop(
      "Expected exactly one final loc_sample result per sim_id/iter/ind, ",
      "but duplicate result keys were found."
    )
  }
  results_summary <- results_raw |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(grid_cols))) |>
    dplyr::summarise(
      dplyr::across(
        c("threshold", "propRespTruth", "propRespEst"),
        list(
          min = ~ min(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE),
          mean = ~ mean(.x, na.rm = TRUE),
          median = ~ stats::median(.x, na.rm = TRUE)
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      propRespEst_median_diff =
        .data$propRespEst_median - .data$propRespTruth_median,
      propRespEst_median_diff_rel =
        .data$propRespEst_median_diff / .data$propRespTruth_median
    )
  list(
    bw_tbl_results_raw = results_raw,
    bw_tbl_results_summary = results_summary
  )
}

# ---------------------------------------------------------------------------
# Analysis 2b: biasUns tuning across fixed bandwidths and batch mismatch
# ---------------------------------------------------------------------------

.simBandwidthBiasUnsScenario <- function(row, settings) {
  cluster_arg <- function(x) {
    x <- as.character(x)[1]
    if (is.na(x) || !nzchar(x)) NULL else x
  }

  bias_basis <- row$bias_uns_basis[[1]]
  bias_multiplier <- row$bias_uns_multiplier[[1]]
  if (!bias_basis %in% c("bandwidth", "negative_width")) {
    stop("Unknown bias_uns_basis: ", bias_basis)
  }

  do.call(.simBandwidthBsFreq, c(settings, list(
    biasUns = if (bias_basis == "bandwidth") {
      bias_multiplier * row$bw[[1]]
    } else {
      0
    },
    biasUnsWidthMultiplier = if (bias_basis == "negative_width") {
      bias_multiplier
    } else {
      NULL
    },
    bw = row$bw[[1]],
    bwFallback = row$bw[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]],
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = row$background_relative_to_response[[1]],
    ncellUnsRelativeToStim = row$n_cell_uns_relative_to_stim[[1]],
    stimMeanShift = row$stim_mean_shift[[1]],
    stimSdMultiplier = row$stim_sd_multiplier[[1]],
    stimMeanShiftClusters = cluster_arg(row$stim_mean_shift_clusters[[1]]),
    stimSdMultiplierClusters = cluster_arg(row$stim_sd_multiplier_clusters[[1]])
  )))
}

.simBandwidthBiasUnsCollate <- function(tbl, grid_cols, n_sample_expected = NULL) {
  results_raw <- tbl |>
    dplyr::filter(.data$method == "loc_sample") |>
    dplyr::select(
      dplyr::any_of(grid_cols),
      "iter", "sample", "ind", "method",
      "propRespTruth", "propRespEst", "threshold",
      "nCellStim", "nCellUns", "nPosStim", "nPosUns",
      "propStim", "propUns",
      "thresholdOrigin", "gateReturnPoint",
      "locGenerated", "locGeneratedDirect", "locSource", "locReason",
      "biasUns", "biasUnsNegativeWidth"
    ) |>
    dplyr::mutate(
      valid_estimate = is.finite(.data$threshold) &
        is.finite(.data$propRespTruth) & .data$propRespTruth > 0 &
        is.finite(.data$propRespEst),
      error = dplyr::if_else(
        .data$valid_estimate,
        .data$propRespEst - .data$propRespTruth,
        NA_real_
      ),
      rel_error = .data$error / .data$propRespTruth,
      abs_rel_error = abs(.data$rel_error)
    )

  if (anyDuplicated(results_raw[c("sim_id", "iter", "ind")]) > 0L) {
    stop(
      "Expected exactly one final loc_sample result per sim_id/iter/ind, ",
      "but duplicate result keys were found."
    )
  }

  if (!setequal(unique(results_raw$sim_id), unique(tbl$sim_id))) {
    stop("Missing final loc_sample results for one or more simulation IDs.")
  }
  if (!is.null(n_sample_expected)) {
    counts <- dplyr::count(results_raw, .data$sim_id, .data$iter)
    if (any(counts$n != n_sample_expected)) {
      stop("Expected ", n_sample_expected, " final sample results per sim_id/iter.")
    }
  }

  results_summary <- results_raw |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(grid_cols))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_valid = sum(.data$valid_estimate),
      n_failed = sum(!.data$valid_estimate),
      failure_fraction = mean(!.data$valid_estimate),
      propRespTruth = stats::median(.data$propRespTruth, na.rm = TRUE),
      propRespEst_median = stats::median(
        .data$propRespEst[.data$valid_estimate], na.rm = TRUE
      ),
      propRespEst_mean = mean(.data$propRespEst[.data$valid_estimate], na.rm = TRUE),
      median_rel_error = stats::median(.data$rel_error, na.rm = TRUE),
      median_abs_rel_error = stats::median(.data$abs_rel_error, na.rm = TRUE),
      q90_abs_rel_error = stats::quantile(
        .data$abs_rel_error,
        probs = 0.9,
        na.rm = TRUE,
        names = FALSE
      ),
      bias_uns_realised = stats::median(.data$biasUns, na.rm = TRUE),
      negative_width = if (any(is.finite(.data$biasUnsNegativeWidth))) {
        stats::median(.data$biasUnsNegativeWidth, na.rm = TRUE)
      } else {
        NA_real_
      },
      threshold_median = stats::median(
        .data$threshold[.data$valid_estimate], na.rm = TRUE
      ),
      prop_stim_median = stats::median(.data$propStim, na.rm = TRUE),
      prop_uns_median = stats::median(.data$propUns, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    bias_uns_results_raw = results_raw,
    bias_uns_results_summary = results_summary
  )
}

# ---------------------------------------------------------------------------
# Analysis 3: base bandwidth estimators
# ---------------------------------------------------------------------------

# One scenario, without seeding; both workers and reruns use RunRow's RNG.
.simBandwidthEstBaseScenario <- function(row, settings) {
  do.call(.simBandwidthEstBwDirect, c(settings, list(
    biasUns = row$bias_uns[[1]],
    bwMtd = row$bw_mtd[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]]
  )))
}

# Non-finite estimates are scientific outcomes, not runtime errors.
# Require the full sample/iteration count even when estimates are non-finite.
.simBandwidthEstBaseValidate <- function(tbl, settings) {
  expected_rows_per_sim <- as.integer(settings$nSample * settings$nIter)
  bad_ids <- tbl |>
    dplyr::count(.data$sim_id, name = "n_rows") |>
    dplyr::filter(.data$n_rows != expected_rows_per_sim) |>
    dplyr::pull(.data$sim_id)
  if (length(bad_ids) == 0L) {
    character()
  } else {
    paste0(
      "unexpected sample-row counts for sim_id: ",
      paste(sort(bad_ids), collapse = ", ")
    )
  }
}

.simBandwidthEstBaseSummary <- function(.data, grid_cols) {
  .data |>
    dplyr::group_by(
      dplyr::pick(dplyr::any_of(grid_cols))
    ) |>
    dplyr::summarise(
      n_bw_total = dplyr::n(),
      n_bw_stim_finite = sum(is.finite(.data$bw_stim)),
      n_bw_uns_finite = sum(is.finite(.data$bw_uns)),
      n_bw_finite = sum(
        is.finite(.data$bw_stim) & is.finite(.data$bw_uns)
      ),
      prop_bw_finite = mean(
        is.finite(.data$bw_stim) & is.finite(.data$bw_uns)
      ),
      mean_bw_stim = .simBandwidthFiniteMean(.data$bw_stim),
      mean_bw_uns = .simBandwidthFiniteMean(.data$bw_uns),
      mean_bw = .simBandwidthFiniteMean(
        dplyr::if_else(
          is.finite(.data$bw_stim) & is.finite(.data$bw_uns),
          pmin(.data$bw_stim, .data$bw_uns),
          NA_real_
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::select(
      dplyr::any_of(grid_cols),
      n_bw_total,
      n_bw_stim_finite,
      n_bw_uns_finite,
      n_bw_finite,
      prop_bw_finite,
      mean_bw_stim,
      mean_bw_uns,
      mean_bw
    )
}


.simBandwidthEstBaseCollate <- function(tbl, grid_cols) {
  list(
    bw_list_raw_mtd = tbl,
    bw_tbl_results = .simBandwidthEstBaseSummary(tbl, grid_cols)
  )
}

# ---------------------------------------------------------------------------
# Analysis 4: ordinary vs normalised bandwidth estimators
# ---------------------------------------------------------------------------

#' Analysis 4 scenario: one bandwidth-estimator simulation
#'
#' @param row data.frame One row of the analysis 4 `sim_grid`.
#' @param settings list Fixed `.simBandwidthEstBwDirect()` arguments.
#' @return tibble `.simBandwidthEstBwDirect()` output (one row per sample
#'   and iteration).
.simBandwidthEstNormScenario <- function(row, settings) {
  do.call(.simBandwidthEstBwDirect, c(settings, list(
    biasUns = row$bias_uns[[1]],
    bwMtd = row$bw_mtd[[1]],
    bwNcellMax = row$bw_ncell_upper[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]]
  )))
}

#' Analysis 4 validation of collated outputs
#'
#' @param n_rows_per_sim integer Expected result rows per `sim_id`.
#' @return function `function(tbl)` returning problem strings: `sim_id`s with
#'   an unexpected row count, or with no finite stim/unstim bandwidth pair.
.simBandwidthEstNormValidator <- function(n_rows_per_sim) {
  function(tbl) {
    if (nrow(tbl) == 0L) {
      return(character())
    }
    per_sim <- tbl |>
      dplyr::group_by(.data$sim_id) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_pair = sum(is.finite(.data$bw_stim) & is.finite(.data$bw_uns)),
        .groups = "drop"
      )
    ids <- function(x) paste(sort(as.integer(x)), collapse = ", ")
    c(
      if (any(per_sim$n_rows != n_rows_per_sim)) {
        paste0(
          "unexpected row counts for sim_id: ",
          ids(per_sim$sim_id[per_sim$n_rows != n_rows_per_sim])
        )
      },
      if (any(per_sim$n_pair == 0L)) {
        paste0(
          "no finite stim/unstim bandwidth pair for sim_id: ",
          ids(per_sim$sim_id[per_sim$n_pair == 0L])
        )
      }
    )
  }
}

#' Analysis 4 collation: raw rows and per-scenario summary
#'
#' `mean_bw` averages `pmin(bw_stim, bw_uns)` over complete finite
#' stim/unstim pairs.
#'
#' @param tbl data.frame Collated analysis 4 outputs.
#' @param grid_cols character Grid column names.
#' @return list `bw_list_raw_mtd` (raw rows) and `bw_tbl_results`.
.simBandwidthEstNormCollate <- function(tbl, grid_cols) {
  results <- tbl |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(grid_cols))) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_bw_stim_finite = sum(is.finite(.data$bw_stim)),
      n_bw_uns_finite = sum(is.finite(.data$bw_uns)),
      n_est = sum(is.finite(.data$bw_stim) & is.finite(.data$bw_uns)),
      n_norm_fallback = sum(
        is.finite(.data$bw_stim) &
          is.finite(.data$bw_uns) &
          .data$bw_norm_fallback %in% TRUE
      ),
      mean_bw_stim = .simBandwidthFiniteMean(.data$bw_stim),
      mean_bw_uns = .simBandwidthFiniteMean(.data$bw_uns),
      mean_bw = .simBandwidthFiniteMean(
        dplyr::if_else(
          is.finite(.data$bw_stim) & is.finite(.data$bw_uns),
          pmin(.data$bw_stim, .data$bw_uns),
          NA_real_
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      estimate_rate = .data$n_est / .data$n_total,
      norm_fallback_rate = dplyr::if_else(
        .data$n_est > 0L,
        .data$n_norm_fallback / .data$n_est,
        NA_real_
      ),
      bw_mtd_base = gsub("Norm$", "", .data$bw_mtd),
      bw_mtd_norm = ifelse(grepl("Norm$", .data$bw_mtd), "norm", "non-norm")
    )
  if (any(results$n_est < results$n_total)) {
    warning(
      "Some estimator/scenario combinations have incomplete stim/unstim ",
      "bandwidth pairs; estimate_rate records the complete-pair fraction."
    )
  }
  list(bw_list_raw_mtd = tbl, bw_tbl_results = results)
}

# ---------------------------------------------------------------------------
# Analysis 5: adaptive bandwidth estimates
# ---------------------------------------------------------------------------

#' Analysis 5 scenario: one adaptive bandwidth estimation simulation
#'
#' @param row data.frame One row of the analysis 5 `sim_grid`.
#' @param settings list Fixed `.simBandwidthEstBwDirectAdaptive()` arguments
#'   (e.g. `nSample`, `nIter`, `normAdaptiveNcell`, `bwFallback`).
#' @return tibble `.simBandwidthEstBwDirectAdaptive()` output, one row per
#'   sample and iteration.
.simBandwidthEstAdaptiveScenario <- function(row, settings) {
  do.call(.simBandwidthEstBwDirectAdaptive, c(settings, list(
    biasUns = row$bias_uns[[1]],
    bwMtd = row$bw_mtd[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]]
  )))
}

#' Analysis 5 validation: expected number of sample-level rows per sim_id
#'
#' @param tbl data.frame Collated analysis 5 outputs.
#' @param rows_per_sim integer Expected rows per `sim_id`
#'   (`nSample * nIter`).
#' @return character Problem strings (`character(0)` when valid).
.simBandwidthEstAdaptiveValidate <- function(tbl, rows_per_sim) {
  n_rows <- table(tbl$sim_id)
  bad_ids <- sort(as.integer(names(n_rows)[n_rows != rows_per_sim]))
  if (length(bad_ids) == 0L) {
    return(character())
  }
  paste0(
    "unexpected sample-row counts for sim_id: ",
    paste(bad_ids, collapse = ", ")
  )
}

#' Analysis 5 collated results
#'
#' Summarises the four adaptive bandwidth estimates (core and extra, for the
#' unstimulated and stimulated samples) per grid row: the number and
#' proportion of finite estimates and their mean. Means are conditional on
#' finite estimates.
#'
#' @param tbl data.frame Collated analysis 5 outputs.
#' @param grid_cols character Grid column names.
#' @return list `bw_list_raw` (the raw outputs) and `bw_tbl_results`.
.simBandwidthEstAdaptiveCollate <- function(tbl, grid_cols) {
  bw_measure_tbl <- tibble::tribble(
    ~bw_component, ~bw_condition, ~bw_col,
    "core", "unstim", "bw_uns_core",
    "core", "stim", "bw_stim_core",
    "extra", "unstim", "bw_uns_extra",
    "extra", "stim", "bw_stim_extra"
  )
  long <- purrr::map_dfr(seq_len(nrow(bw_measure_tbl)), function(i) {
    dplyr::mutate(
      tbl,
      bw_component = bw_measure_tbl$bw_component[[i]],
      bw_condition = bw_measure_tbl$bw_condition[[i]],
      bw_est = .data[[bw_measure_tbl$bw_col[[i]]]]
    )
  }) |>
    dplyr::mutate(
      bw_component = factor(.data$bw_component, levels = c("core", "extra")),
      bw_condition = factor(.data$bw_condition, levels = c("unstim", "stim"))
    )
  group_cols <- c(grid_cols, "bw_component", "bw_condition")
  results <- long |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(group_cols))) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_est = sum(is.finite(.data$bw_est)),
      prop_est = .data$n_est / .data$n_total,
      mean_bw = .simBandwidthFiniteMean(.data$bw_est),
      .groups = "drop"
    ) |>
    dplyr::mutate(bw_mtd_base = gsub("Norm$", "", .data$bw_mtd))
  list(bw_list_raw = tbl, bw_tbl_results = results)
}

# ---------------------------------------------------------------------------
# Analysis 6: manually adaptive background-subtracted frequency
# ---------------------------------------------------------------------------

# One scientific scenario; RNG and grid metadata belong to the shared runner.
.simBandwidthFreqBsAdaptiveScenario <- function(row, settings) {
  crossover <- row$bw_crossover[[1]]
  do.call(.simBandwidthBsFreq, c(settings, list(
    biasUns = row$bias_uns[[1]],
    bwAdaptiveCore = row$bw_core[[1]],
    bwAdaptiveExtra = row$bw_extra[[1]],
    bwAdaptiveCrossover = if (is.finite(crossover)) crossover else NULL,
    bwAdaptiveTransitionWidth = row$bw_transition_width[[1]],
    bwFallback = row$bw_fallback[[1]],
    nCellStim = row$n_cell[[1]],
    probResponse = row$prob_response[[1]],
    meanPos = row$mean_pos[[1]],
    transformation = row$transformation[[1]],
    samplePerturbationSd = row$sample_perturbation_sd[[1]],
    conditionPerturbationSd = row$condition_perturbation_sd[[1]],
    clusterPerturbationSd = row$cluster_perturbation_sd[[1]],
    backgroundRelativeToResponse = row$background_relative_to_response[[1]],
    ncellUnsRelativeToStim = row$n_cell_uns_relative_to_stim[[1]]
  )))
}

# Infrastructure failures are checked by the shared validator first.
.simBandwidthFreqBsAdaptiveValidate <- function(tbl) {
  required <- c(
    "method", "threshold", "propRespTruth", "propRespEst",
    "sim_id", "iter", "ind"
  )
  if (!all(required %in% names(tbl))) {
    return("Missing final sample result columns")
  }
  final <- tbl |>
    dplyr::filter(
      .data$method == "loc_sample",
      is.finite(.data$threshold),
      is.finite(.data$propRespTruth),
      is.finite(.data$propRespEst)
    )
  problems <- character()
  if (!setequal(final$sim_id, tbl$sim_id)) {
    problems <- "At least one sim_id had no finite final loc_sample result"
  }
  if (anyDuplicated(final[c("sim_id", "iter", "ind")]) > 0L) {
    problems <- c(problems, "Duplicate final sim_id/iter/ind result keys")
  }
  problems
}

# The frequency estimand and scenario summary match analysis 2. Persist the
# additional threshold summary during promotion, never during read-only plots.
.simBandwidthFreqBsAdaptiveCollate <- function(tbl, grid_cols) {
  results <- .simBandwidthFreqBsGlobalCollate(tbl, grid_cols)
  summary_tbl <- results$bw_tbl_results_raw |>
    dplyr::group_by(
      dplyr::pick(dplyr::any_of(grid_cols))
    ) |>
    dplyr::summarise(
      threshold_min = min(threshold, na.rm = TRUE),
      threshold_max = max(threshold, na.rm = TRUE),
      threshold_median = stats::median(threshold, na.rm = TRUE),
      threshold_iqr_lower = stats::quantile(threshold, 0.25, na.rm = TRUE, names = FALSE),
      threshold_iqr_upper = stats::quantile(threshold, 0.75, na.rm = TRUE, names = FALSE),
      threshold_sd = stats::sd(threshold, na.rm = TRUE),
      threshold_mad = stats::mad(threshold, na.rm = TRUE),
      threshold_iqr_length = .data$threshold_iqr_upper - .data$threshold_iqr_lower,
      propRespTruth_median = stats::median(propRespTruth, na.rm = TRUE),
      propRespEst_median = stats::median(propRespEst, na.rm = TRUE),
      propRespEst_median_diff =
        .data$propRespEst_median - .data$propRespTruth_median,
      propRespEst_median_diff_rel =
        (.data$propRespEst_median - .data$propRespTruth_median) /
          .data$propRespTruth_median,
      .groups = "drop"
    )
  c(results, list(summary_tbl = summary_tbl))
}
