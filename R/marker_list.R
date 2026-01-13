# Complete marker parameter list with all required settings
# Ensures all parameters for each marker requiring a gate are properly defined

.complete_marker_list <- function(marker,
                                  bias_uns,
                                  bias_uns_factor,
                                  exc_min,
                                  .data,
                                  pop_gate,
                                  ind_batch_list,
                                  bw_min,
                                  cp_min,
                                  min_cell,
                                  tol_clust,
                                  max_pos_prob_x,
                                  gate_combn,
                                  marker_settings,
                                  path_project,
                                  .debug) {
  marker_settings_common <- list(
    bias_uns = bias_uns, bias_uns_factor = bias_uns_factor,
    exc_min = exc_min, cp_min = cp_min, bw_min = bw_min,
    min_cell = min_cell, tol_clust = tol_clust, gate_combn = gate_combn,
    max_pos_prob_x = max_pos_prob_x
  )
  chnl_lab <- stimgate_meta_settings_get_chnl_lab(path_project)
  marker_list <- purrr::map(marker, function(marker_curr) {
    .complete_marker_list_ind(
      marker = marker_curr,
      marker_settings_common = marker_settings_common,
      marker_settings_spec = list(
        marker = chnl_lab[[marker_curr]],
        chnl_cut = marker_curr
      ) |>
        append(marker_settings[[marker_curr]]),
      .data = .data,
      pop_gate = pop_gate,
      .debug = .debug,
      ind_batch_list = ind_batch_list,
      path_project = path_project
    )
  }) |>
    stats::setNames(chnl_lab[marker])

  .complete_marker_list_save(
    marker_list = marker_list,
    path_project = path_project
  )

  marker_list
}

.complete_marker_list_ind <- function(marker_settings_common,
                                      marker_settings_spec,
                                      marker,
                                      .data,
                                      pop_gate,
                                      .debug,
                                      ind_batch_list,
                                      path_project) {
  marker_settings <- .complete_marker_list_add_common(
    marker_settings_common = marker_settings_common,
    marker_settings = marker_settings_spec
  )

  marker_settings$bias_uns <- .complete_marker_list_bias_uns(
    bias_uns = marker_settings$bias_uns,
    bias_uns_factor = marker_settings$bias_uns_factor,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    .debug = .debug,
    ind_batch_list = ind_batch_list,
    path_project = path_project
  )

  marker_settings$bw_min <- .complete_marker_list_min_bw(
    bw_min = marker_settings$bw_min,
    bias_uns = max(marker_settings$bias_uns)
  )

  marker_settings$cp_min <- .complete_marker_list_cp_min(
    cp_min = marker_settings$cp_min,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    .debug = .debug,
    ind_batch_list = ind_batch_list,
    path_project = path_project
  )

  marker_settings
}

.complete_marker_list_add_common <- function(marker_settings_common,
                                             marker_settings) {
  marker_settings |>
    append(marker_settings_common[
      setdiff(names(marker_settings_common), names(marker_settings))
    ])
}

.complete_marker_list_bias_uns <- function(bias_uns,
                                           bias_uns_factor,
                                           .data,
                                           pop_gate,
                                           chnl_cut,
                                           .debug,
                                           ind_batch_list,
                                           path_project) {
  if (!is.null(bias_uns)) {
    return(bias_uns)
  }
  .debug_msg(.debug, "calculating bias_uns automatically") # nolint
  mean_range <- .complete_marker_list_bias_uns_get_mean_range(
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut,
    path_project = path_project
  )
  (mean_range / 12 * bias_uns_factor) |> signif(3)
}

.complete_marker_list_bias_uns_get_mean_range <- function(ind_batch_list,
                                                          .data,
                                                          pop_gate,
                                                          chnl_cut,
                                                          path_project) {
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i],
        path_project = path_project
      )
      purrr::map_dbl(ex_list, function(ex) {
        abs(
          diff(quantile(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], c(0.99, 0.01)),
            na.rm = TRUE
          )
        )[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


.complete_marker_list_min_bw <- function(bw_min, bias_uns) {
  if (!is.null(bw_min)) {
    return(bw_min)
  }
  bias_uns * 2.25
}

.complete_marker_list_cp_min <- function(cp_min,
                                         .data,
                                         pop_gate,
                                         chnl_cut,
                                         .debug,
                                         ind_batch_list,
                                         path_project) {
  if (!is.null(cp_min)) {
    return(cp_min)
  }
  .debug_msg(.debug, "calculating cp_min automatically") # nolint
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i],
        path_project = path_project
      )
      purrr::map_dbl(ex_list, function(ex) {
        median(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], na.rm = TRUE)[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}



# Get all cutpoint type names
# Returns character vector of all available cutpoint names
.get_full_cp_type_vec <- function(fdr) {
  # Get cutpoint names for unstim-based cuts
  # cp_name_vec_uns <- .get_cp_uns_name_vec(fdr)
  # cp_name_vec_uns_root <- paste0(cp_name_vec_uns, 'root')

  # output all cutpoint names
  c(
    "man", "tg", "dcp", "midp", "scp",
    "uns", "unsr", "loc"
  )
}

.complete_marker_list_save <- function(marker_list, path_project) {
  path_save <- file.path(path_project, "meta_data", "marker_list.rds")
  if (file.exists(path_save)) {
    invisible(file.remove(path_save))
  }
  if (!dir.exists(path_project)) {
    dir.create(path_project, recursive = TRUE)
  }
  saveRDS(
    marker_list,
    file = path_save
  )
}

.chnl_lab <- function(.data) {
  adf <- switch(class(.data)[1],
    "GatingSet" = {
      gh <- .data[[1]]
      fr <- flowWorkspace::gh_pop_get_data(gh)
      flowCore::parameters(fr)@data
    },
    "GatingHierarchy" = {
      fr <- flowWorkspace::gh_pop_get_data(.data)
      flowCore::parameters(fr)@data
    },
    "flowFrame" = flowCore::parameters(.data)@data,
    "flowSet" = flowCore::parameters(.data[[1]])@data,
    "cytoframe" = flowCore::parameters(.data)@data,
    "cytoset" = flowCore::parameters(.data[[1]])@data,
    stop(paste0("Unsupported data class: ", class(.data)[1]))
  )

  lab_vec <- setNames(adf$desc, adf$name)
  for (i in seq_along(lab_vec)) {
    if (is.na(lab_vec[i])) {
      lab_vec[i] <- names(lab_vec)[i]
    }
  }

  lab_vec
}

#' @title Read marker settings from project
#' @description Read the saved marker settings list from the project's meta_data folder.
#' @param path_project character Path to project.
#' @return A named list of marker settings (as saved by .complete_marker_list_save()).
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
#' stimgate_meta_settings_get_markers(tmp)
#' }
#' @export
stimgate_meta_settings_get_markers <- function(path_project) {
  path_marker_list <- file.path(path_project, "meta_data", "marker_list.rds")
  if (!file.exists(path_marker_list)) {
    stop("Marker list file not found in project meta_data folder")
  }
  readRDS(path_marker_list)
}
# ...existing code...
#' @title Read marker list with channel labels
#' @description Read the project's marker list and return it with element names
#'   replaced by channel labels (from chnl_lab).
#' @param path_project character Path to project.
#' @return A named list of marker settings where names are channel labels.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
#' stimgate_meta_settings_get_chnls(tmp)
#' }
#' @export
stimgate_meta_settings_get_chnls <- function(path_project) {
  marker_list <- stimgate_meta_settings_get_markers(path_project)
  chnl_lab <- stimgate_meta_settings_get_chnl_lab(path_project)
  # chnl_lab maps channel code -> label; rename marker_list elements using labels
  names(marker_list) <- chnl_lab[names(marker_list)]
  marker_list
}
# ...existing code...
#' @title Get marker settings for a single channel
#' @description Retrieve the marker settings for a single channel. The function
#'   accepts either a channel label (as returned by stimgate_meta_settings_get_chnl_lab)
#'   or the original channel name/key used in marker_list.
#' @param path_project character Path to project.
#' @param chnl character Channel label or channel name.
#' @return A list of settings for the requested channel.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
#' stimgate_meta_settings_get_chnl(tmp, "BC1 label")
#' stimgate_meta_settings_get_chnl(tmp, "BC1")
#' }
#' @export
stimgate_meta_settings_get_chnl <- function(path_project, chnl) {
  chnl_list <- stimgate_meta_settings_get_chnls(path_project)
  # If user supplied the label (names of chnl_list)
  if (chnl %in% names(chnl_list)) {
    return(chnl_list[[chnl]])
  }
  # If user supplied the original channel key/name, map via chnl_lab
  marker_list_orig <- stimgate_meta_settings_get_markers(path_project)
  if (chnl %in% names(marker_list_orig)) {
    chnl_lab <- stimgate_meta_settings_get_chnl_lab(path_project)
    label <- chnl_lab[[chnl]]
    if (!is.null(label) && label %in% names(chnl_list)) {
      return(chnl_list[[label]])
    }
  }
  stop(sprintf("Channel %s not found in marker list", chnl))
}

#' @title Get settings for a named marker
#' @description Retrieve the settings for a marker by its original name/key.
#' @param path_project character Path to project.
#' @param marker character Marker name/key as stored in marker_list.
#' @return A list of settings for the requested marker.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
#' stimgate_meta_settings_get_marker(tmp, "BC1")
#' }
#' @export
stimgate_meta_settings_get_marker <- function(path_project, marker) {
  marker_list <- stimgate_meta_settings_get_markers(path_project)
  if (!marker %in% names(marker_list)) {
    stop(paste0("Marker ", marker, " not found in marker list"))
  }
  marker_list[[marker]]
}

#' @title Read channel label mapping
#' @description Read the saved channel label mapping (chnl_lab.rds) from the project's meta_data folder.
#' @param path_project character Path to project.
#' @return Named character vector mapping channel names to labels.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
#' stimgate_meta_settings_get_chnl_lab(tmp)
#' }
#' @export
stimgate_meta_settings_get_chnl_lab <- function(path_project) {
  path_chnl_lab <- file.path(path_project, "meta_data", "chnl_lab.rds")
  if (!file.exists(path_chnl_lab)) {
    stop("Channel label file not found in project meta_data folder")
  }
  readRDS(path_chnl_lab)
}

.save_meta_data <- function(.data, path_project) {
  path_dir_meta_data <- file.path(path_project, "meta_data")
  if (!dir.exists(path_dir_meta_data)) {
    dir.create(path_dir_meta_data, recursive = TRUE)
  }
  chnl_lab <- .chnl_lab(.data)
  saveRDS(
    chnl_lab,
    file = file.path(path_dir_meta_data, "chnl_lab.rds")
  )
}