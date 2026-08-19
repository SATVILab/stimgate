.path_sim_output <- function(sim_id, dir_output, sim_grid_chunk_index, sim_grid_n_chunks) {
  if (is.null(dir_output) || !nzchar(dir_output)) {
    return(character(0))
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
  sim_grid_chunk_index,
  sim_grid_n_chunks,
  dir_output
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

  summary_text <- paste0(
    "========================================\n",
    " ADAPTIVE SIMULATION PROGRESS DASHBOARD \n",
    "========================================\n",
    "Chunk              : ", sim_grid_chunk_index, " of ", sim_grid_n_chunks, "\n",
    "Total Simulations  : ", total_sims, "\n",
    "In Progress        : ", n_running, "\n",
    "Completed (Success): ", n_completed, "\n",
    "Completed (Errors) : ", n_error, "\n",
    "Total Done          : ", n_done, "\n",
    "Output directory   : ", dir_output, "\n",
    "----------------------------------------\n",
    "Currently Running IDs:\n", running_list_str, "\n",
    "========================================\n"
  )

  cat(summary_text, file = path_progress_file, append = FALSE)
  invisible(summary_text)
}

update_progress_summary <- .update_progress_summary

.find_bw_list_output_files <- function(output_dir = NULL, cache_dir = NULL) {
  output_dir_vec <- unique(c(
    if (!is.null(output_dir) && nzchar(output_dir)) output_dir else character(),
    if (!is.null(cache_dir) && nzchar(cache_dir)) file.path(cache_dir, "output") else character()
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
        full.names = TRUE
      )
    }
  ), use.names = FALSE))
}

find_bw_list_output_files <- .find_bw_list_output_files
