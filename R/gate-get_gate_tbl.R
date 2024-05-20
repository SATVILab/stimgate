#' @param data GatingSet.
#' Original GatingSet algorithm was applied to, with the same name.
#' @param pop_gate character.
#' Name of population to gate.
#' @param cut character vector.
#' Channels to gate on.
#' @param ind_in_batch_lab_vec character.
#' I-th element specifies name of i-th sample within a batch.
#' @param gate_name character. Name of gate.
#' @param ind_in_batch_gate numeric vector.
#' Specifies indices for each batch to gate.
#' @param ind_in_batch_uns numeric.
#' Specifies which index is unstim within each batch.
#'
#' @return
#' Gate table with gates for each sample for each marker.
#' @export
get_gate_tbl <- function(data,
                         pop_gate,
                         cut,
                         ind_in_batch_lab_vec,
                         gate_name,
                         ind_in_batch_gate,
                         ind_in_batch_uns,
                         path_project) {
  # ==============================
  # Preparation
  # ==============================

  # Create params list
  # -----------------------------

  # dataset_name
  data_name <- deparse(substitute(data))

  # chnl_lab
  chnl_lab_vec <- .get_labs(
    data = data[[1]],
    cut = cut
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name,
    cut = cut[1]
  )

  # get directory to save to
  # ------------------------------

  # get base directory
  dir_base <- stim_gate_dir_base_create(
    dir_base_init = path_project,
    params = params
  )

  dir_base <- stringr::str_remove(dir_base, paste0("/", chnl_lab_vec[cut[1]]))

  # get directory to save to
  marker_vec <- chnl_lab_vec[cut]
  markers <- paste0(marker_vec, collapse = "_")
  dir_save <- file.path(dir_base, markers, gate_name, "fcs")
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # ====================================
  # Get gates
  # ====================================

  purrr::map_df(cut, function(cut_curr) {
    params <- params[-which(names(params) == "cut")]
    # get base directory
    dir_base <- stim_gate_dir_base_create(
      dir_base_init = path_project,
      params = params |> append(list(cut = cut_curr))
    )
    # get stats tbl
    gate_stats <- readRDS(file.path(dir_base, "stats", "gate_stats_tbl"))

    gate_stats <- gate_stats |>
      dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())

    gate_stats
  })
}
