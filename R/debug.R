#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom  utils head

# Global variable bindings to avoid R CMD check notes
# These are primarily used in dplyr and ggplot2 contexts
globalVariables(c(
  # Variables used in dplyr operations
  "marker",
  "batch",
  "ind",
  "gate",
  "gate_cyt",
  "gate_single",
  "gate_name",
  "chnl",
  "gate_use",
  "gate_type",
  "gate_combn",
  "pop_gate",
  "gate_tbl",
  "chnl_cut",
  "tol",
  "ind_in_batch_gate",
  "tol_clust_single",
  "ind_batch_gate",
  "cp",
  "grp",
  "cp_join_lse_orig_mean_tg",
  "cp_orig_quant_min",
  "cp_join",
  "cp_join_lse",
  "cp_join_lse_orig",
  "cp_join_lse_orig_mean",
  "cp_join_tg_orig",
  "cp_join_tg_orig_mean",
  "prop_bs_orig",
  "prop_bs_cp_diff",
  "prop_bs_cp_diff_sd",
  "prop_bs_cp",
  "prop_l1se",
  "pred",
  "der",
  "cp_orig",
  "max_expr",
  "gate_05",
  "prop_bs_cp_diff_sd_max",
  "grp_level",
  "ind_vec",
  "x1",
  "x",
  "y",
  "x_ind",
  "count_stim",
  "n_cell_stim",
  "count_uns",
  "n_cell_uns",
  "prop_stim_pos",
  "prop_uns_pos",
  "prop_stim_sd",
  "prop_uns_sd",
  "x512",
  "x_vec",
  "freq_bs",
  "freq_stim",
  "pop_gate_curr",
  "cp_join_tg",
  "cp_orig_quant_min",
  "lse_orig",
  "cp_tg_ctrl",
  "cp_join_lse_orig_mean",
  "cp_join_tg_orig",
  "chnl_pos",
  "dir_save",
  "is.null_gate_tbl",
  "path_project",
  "exc_min",
  "prop_bs_diff",
  "prop_stim",
  "prop_uns",
  "prop_bs",
  "prob_smooth",
  "n_row",
  "y_stim",
  "y_uns",
  "stim",
  "x_stim",
  "prob",
  "x_uns",
  "prop_lab",
  "type",
  "line_id",
  "dens",
  "no",
  "yes",
  "prob_stim",
  "prob_stim_norm",
  "prop_pos",
  # Variables from .get_cp_uns_loc_prob_tbl_filter
  "minor_response_ind",
  "moderate_response_ind",
  "n_remaining",
  "prob_larger_count",
  "prob_larger_prop",
  # Variables from .get_prop_bs_by_cp_tbl_ind_calc
  "count_stim_cp",
  "count_uns_cp",
  "prop_stim_cp",
  "prop_uns_cp",
  "prop_bs_sd",
  "prop_stim_pos_cp",
  "prop_uns_pos_cp",
  "prop_stim_sd_cp",
  "prop_uns_sd_cp",
  "prop_bs_sd_cp",
  # Variables from other functions
  "cyt_combn",
  "freq_uns",
  "V1",
  "V2",
  "i",
  "tol_gate_single",
  # Variables used in plots and ggplot2 context
  ".debug",
  "prop_l1se",
  "rbeta",
  # Variables from fcs_write.R
  "concat",
  "gate_concat"
))

#' Create a debug file in tempdir()
#'
#' @return character Path to the created debug file (invisibly).
#' @keywords internal
.debug_file_create <- function() {
  dir_path <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  timestamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
  file_path <- file.path(dir_path, paste0("stimgate_", timestamp, ".txt"))
  file.create(file_path, showWarnings = FALSE)
  invisible(file_path)
}

#' Get the most recent debug file path
#'
#' @return character Path to the most recent debug file, or NULL if none.
#' @keywords internal
.debug_file_get_path <- function() {
  dir_path <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dir_path)) {
    return(NULL)
  }
  files <- list.files(
    dir_path,
    pattern = "^stimgate_\\d{4}-\\d{2}-\\d{2}-\\d{6}\\.txt$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(NULL)
  }
  files[which.max(file.info(files)$mtime)]
}

#' Print debug message conditionally
#'
#' Writes debug output to the most recent stimgate debug file when
#' STIMGATE_DEBUG is enabled.
#'
#' @param msg character Message to print
#' @param val object Optional value to append to message. Default is NULL.
#' @return logical invisibly TRUE if message was written, FALSE otherwise
#' @keywords internal
.debug <- function(msg, val = NULL) {
  must_debug <- tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")
  if (!must_debug) {
    return(invisible(FALSE))
  }
  if (!is.null(val)) {
    msg <- paste0(msg, ": ", val)
  }
  path_debug <- .debug_file_get_path()
  if (is.null(path_debug)) {
    path_debug <- .debug_file_create()
  }
  cat(msg, file = path_debug, sep = "\n", append = TRUE)
  invisible(TRUE)
}

#' Copy the latest stimgate debug file to the working directory
#'
#' @description Copies the most recent debug file created by stimgate
#' to the current working directory. The copied file uses the same
#' filename as the source.
#' @param path_dir character. Directory to copy the debug file to.
#' Default is `getwd()` (i.e. the working directory).
#'
#' @return character Path to the copied file (invisibly), or NULL if no
#'   debug file exists.
#' @export
stimgate_debug_copy <- function(path_dir = getwd()) {
  src <- .debug_file_get_path()
  if (is.null(src)) {
    message("No stimgate debug file found.")
    return(invisible(NULL))
  }
  if (!dir.exists(path_dir)) {
    dir.create(path_dir, recursive = TRUE)
  }
  dest <- file.path(path_dir, basename(src))
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) {
    message("Failed to copy debug file to working directory.")
    return(invisible(NULL))
  }
  invisible(dest)
}

#' Print the latest stimgate debug file to the console
#'
#' @description Prints the most recent debug file created by stimgate
#' to the console.
#'
#' @return character The debug text (invisibly), or NULL if no debug file exists.
#' @export
stimgate_debug_print <- function() {
  src <- .debug_file_get_path()
  if (is.null(src)) {
    message("No stimgate debug file found.")
    return(invisible(NULL))
  }
  txt <- readLines(src, warn = FALSE)
  cat(txt, sep = "\n")
  invisible(txt)
}

.int_save_nm <- function(name, obj, ind, stage, path_project) {
  if (!.int_save_check(ind)) {
    return(invisible(FALSE))
  }
  path_save <- .int_save_path_save(path_project, stage, ind, name)
  saveRDS(obj, path_save)
  invisible(TRUE)
}

.int_save <- function(ind, stage, path_project, ...) {
  if (!.int_save_check(ind)) {
    return(invisible(FALSE))
  }

  dots <- list(...)
  dot_names <- names(dots)

  call_names <- as.list(substitute(list(...)))[-1]
  call_names <- vapply(
    call_names,
    function(x) paste(deparse(x), collapse = ""),
    character(1)
  )

  if (is.null(dot_names)) {
    dot_names <- call_names
  } else {
    dot_names[dot_names == ""] <- call_names[dot_names == ""]
  }

  for (i in seq_along(dots)) {
    .int_save_nm(dot_names[[i]], dots[[i]], ind, stage, path_project)
  }

  invisible(TRUE)
}

#' @keywords internal
.is_invalid_ind <- function(ind) {
  is.null(ind) || length(ind) == 0 || all(is.na(ind))
}

.int_save_check <- function(ind) {
  # Return FALSE if ind is NULL or invalid
  if (.is_invalid_ind(ind)) {
    return(FALSE)
  }

  env_var <- Sys.getenv("STIMGATE_INTERMEDIATE") |>
    trimws() |>
    tolower()
  if (is.null(env_var) || length(env_var) == 0 || env_var == "") {
    return(FALSE)
  }
  if (env_var %in% c("y", "true", "yes", "all")) {
    return(TRUE)
  }
  env_var_split <- strsplit(env_var, ",|;") |>
    unlist() |>
    trimws()
  any(as.character(ind) %in% env_var_split)
}

.int_save_path_save <- function(path_project, stage, ind, name) {
  name <- paste0(name, ".rds")
  name <- gsub("\\.rds(\\.rds)*$", ".rds", name, ignore.case = TRUE)
  path_save <- file.path(
    path_project,
    "intermediate_data",
    stage,
    "ind",
    paste0(as.character(ind), collapse = "_"),
    name
  )
  if (!dir.exists(dirname(path_save))) {
    dir.create(dirname(path_save), recursive = TRUE, showWarnings = FALSE)
  }
  path_save
}

.browse <- function(ind) {
  if (!.browse_check(ind)) {
    return(invisible(FALSE))
  }
  eval(quote(browser()), envir = parent.frame())
  invisible(TRUE)
}


.browse_check <- function(ind) {
  # Return FALSE if ind is NULL or invalid
  if (.is_invalid_ind(ind)) {
    return(FALSE)
  }

  env_var <- Sys.getenv("STIMGATE_BROWSE") |>
    trimws() |>
    tolower()
  if (is.null(env_var) || length(env_var) == 0 || env_var == "") {
    return(FALSE)
  }
  if (env_var %in% c("y", "true", "yes", "all")) {
    return(TRUE)
  }
  env_var_split <- strsplit(env_var, ",|;") |>
    unlist() |>
    trimws()
  any(as.character(ind) %in% env_var_split)
}
