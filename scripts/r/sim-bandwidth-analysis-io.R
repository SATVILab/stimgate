.path_sim_output <- function(
  sim_id,
  dir_output = NULL,
  sim_grid_chunk_index = NULL,
  sim_grid_n_chunks = NULL
) {
  if (is.null(dir_output) || !nzchar(dir_output)) {
    dir_output <- if (exists("dir_output", inherits = TRUE)) {
      get("dir_output", mode = "any", inherits = TRUE)
    } else {
      NULL
    }
  }
  if (is.null(dir_output) || !nzchar(dir_output)) {
    return(character(0))
  }

  if (is.null(sim_grid_chunk_index) || is.na(sim_grid_chunk_index)) {
    sim_grid_chunk_index <- if (exists("sim_grid_chunk_index", inherits = TRUE)) {
      get("sim_grid_chunk_index", mode = "any", inherits = TRUE)
    } else {
      1L
    }
  }
  if (is.null(sim_grid_n_chunks) || is.na(sim_grid_n_chunks)) {
    sim_grid_n_chunks <- if (exists("sim_grid_n_chunks", inherits = TRUE)) {
      get("sim_grid_n_chunks", mode = "any", inherits = TRUE)
    } else {
      1L
    }
  }

  file.path(
    dir_output,
    sprintf(
      "bw_list_raw-chunk_%03d-of_%03d-sim_id_%06d.rds",
      as.integer(sim_grid_chunk_index),
      as.integer(sim_grid_n_chunks),
      as.integer(sim_id)
    )
  )
}

path_sim_output <- .path_sim_output

.update_progress_summary <- function(
  path_progress_file,
  dir_jobs_chunk,
  total_sims,
  sim_grid_chunk_index = NULL,
  sim_grid_n_chunks = NULL,
  dir_output = NULL,
  heading = "ADAPTIVE SIMULATION PROGRESS DASHBOARD"
) {
  if (!dir.exists(dir_jobs_chunk)) {
    return(invisible(NULL))
  }

  all_job_files <- list.files(dir_jobs_chunk, full.names = FALSE)
  if (length(all_job_files) == 0L) {
    all_job_files <- character(0)
  }

  n_running <- sum(grepl("^running-", all_job_files))
  n_completed <- sum(grepl("^completed-", all_job_files))
  n_error <- sum(grepl("^error-", all_job_files))
  n_done <- n_completed + n_error

  running_files <- all_job_files[grepl("^running-", all_job_files)]
  running_ids <- sub("^running-", "", running_files)
  running_list_str <- if (length(running_ids) > 0L) {
    paste(running_ids, collapse = ", ")
  } else {
    "None"
  }

  has_chunk_info <- !is.null(sim_grid_chunk_index) && !is.null(sim_grid_n_chunks)
  has_output_dir <- !is.null(dir_output) && nzchar(dir_output)

  summary_text <- paste0(
    "========================================\n",
    " ", heading, " \n",
    "========================================\n",
    if (has_chunk_info) {
      paste0("Chunk              : ", sim_grid_chunk_index, " of ", sim_grid_n_chunks, "\n")
    } else {
      ""
    },
    "Total Simulations  : ", total_sims, "\n",
    "In Progress        : ", n_running, "\n",
    "Completed (Success): ", n_completed, "\n",
    "Completed (Errors) : ", n_error, "\n",
    "Total Done          : ", n_done, "\n",
    if (has_output_dir) {
      paste0("Output directory   : ", dir_output, "\n")
    } else {
      ""
    },
    "----------------------------------------\n",
    "Currently Running IDs:\n", running_list_str, "\n",
    "========================================\n"
  )

  # Tolerate concurrent-write conflicts across parallel workers writing the
  # same progress file; a transient dashboard write collision should not
  # fail an otherwise valid simulation.
  tryCatch(
    {
      cat(summary_text, file = path_progress_file, append = FALSE)
    },
    error = function(e) NULL
  )
  invisible(summary_text)
}

update_progress_summary <- .update_progress_summary

.find_bw_list_output_files <- function(output_dir = NULL, cache_dir = NULL, cache_path = NULL) {
  output_dir_vec <- unique(c(
    if (!is.null(output_dir) && nzchar(output_dir)) output_dir else character(),
    if (!is.null(cache_dir) && nzchar(cache_dir)) cache_dir else character()
  ))

  output_dir_vec <- output_dir_vec[dir.exists(output_dir_vec)]
  if (length(output_dir_vec) == 0L) {
    return(character(0))
  }

  unique(unlist(lapply(
    output_dir_vec,
    function(path) {
      list.files(
        path,
        pattern = "^bw_list_raw-chunk_[0-9]+-of_[0-9]+-sim_id_[0-9]+[.]rds$",
        full.names = TRUE,
        recursive = TRUE
      )
    }
  ), use.names = FALSE))
}

find_bw_list_output_files <- .find_bw_list_output_files
