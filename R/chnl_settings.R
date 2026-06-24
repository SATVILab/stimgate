# Complete marker parameter list with all required settings
# Ensures all parameters for each marker requiring a gate are properly defined

#' @keywords internal
.complete_chnl_settings <- function(
  chnl,
  marker,
  chnl_settings,
  marker_settings,
  bias_uns,
  bias_uns_factor,
  exc_min,
  .data,
  pop_gate,
  ind_batch_list,
  bw,
  bw_min,
  bw_max,
  bw_mtd,
  bw_adj,
  bw_ncell_min,
  bw_ncell_max,
  cp_min,
  min_cell,
  tol_clust,
  max_pos_prob_x,
  gate_combn,
  gate_quant,
  path_project
) {
  chnl_settings <- .extract_chnl_settings(
    chnl = chnl,
    marker = marker,
    chnl_settings = chnl_settings,
    marker_settings = marker_settings,
    path_project = path_project
  )
  chnl_settings_common <- list(
    bias_uns = bias_uns,
    bias_uns_factor = bias_uns_factor,
    exc_min = exc_min,
    cp_min = cp_min,
    bw_min = bw_min,
    bw = bw,
    bw_max = bw_max,
    bw_mtd = bw_mtd,
    bw_adj = bw_adj,
    bw_ncell_min = bw_ncell_min,
    bw_ncell_max = bw_ncell_max,
    min_cell = min_cell,
    tol_clust = tol_clust,
    gate_combn = gate_combn,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant
  )
  chnl_lab <- stimgate_meta_read_chnl_lab(path_project)
  chnl_list <- purrr::map(chnl, function(chnl_curr) {
    chnl_settings_spec_curr <- list(
      marker = chnl_lab[[chnl_curr]],
      chnl_cut = chnl_curr
    ) |>
      append(chnl_settings[[chnl_curr]])

    .verify_chnl_settings_chnl(chnl_curr, chnl_settings_spec_curr)
    .complete_chnl_settings_ind(
      chnl = chnl_curr,
      chnl_settings_common = chnl_settings_common,
      chnl_settings_spec = chnl_settings_spec_curr,
      .data = .data,
      pop_gate = pop_gate,
      ind_batch_list = ind_batch_list,
      path_project = path_project
    )
  }) |>
    stats::setNames(chnl_lab[chnl])

  .complete_chnl_settings_save(
    chnl_list = chnl_list,
    path_project = path_project
  )

  chnl_list
}

#' @keywords internal
.complete_chnl_settings_ind <- function(
  chnl_settings_common,
  chnl_settings_spec,
  chnl,
  .data,
  pop_gate,
  ind_batch_list,
  path_project
) {
  chnl_settings <- .complete_chnl_settings_add_common(
    chnl_settings_common = chnl_settings_common,
    chnl_settings = chnl_settings_spec
  )

  chnl_settings$bias_uns <- .complete_chnl_settings_bias_uns(
    bias_uns = chnl_settings$bias_uns,
    bias_uns_factor = chnl_settings$bias_uns_factor,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl,
    ind_batch_list = ind_batch_list,
    path_project = path_project
  )

  chnl_settings$bw_min <- .complete_chnl_settings_min_bw(
    bw_min = chnl_settings$bw_min,
    bias_uns = max(chnl_settings$bias_uns)
  )

  chnl_settings$cp_min <- .complete_chnl_settings_cp_min(
    cp_min = chnl_settings$cp_min,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl,
    ind_batch_list = ind_batch_list,
    path_project = path_project
  )

  chnl_settings
}

#' @keywords internal
.complete_chnl_settings_add_common <- function(
  chnl_settings_common,
  chnl_settings
) {
  chnl_settings |>
    append(chnl_settings_common[
      setdiff(names(chnl_settings_common), names(chnl_settings))
    ])
}

#' @keywords internal
.complete_chnl_settings_bias_uns <- function(
  bias_uns,
  bias_uns_factor,
  .data,
  pop_gate,
  chnl_cut,
  ind_batch_list,
  path_project
) {
  if (!is.null(bias_uns)) {
    return(bias_uns)
  }
  .debug("calculating bias_uns automatically") # nolint
  mean_range <- .complete_chnl_settings_bias_uns_get_mean_range(
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut,
    path_project = path_project
  )
  (mean_range / 12 * bias_uns_factor) |> signif(3)
}

#' @keywords internal
.complete_chnl_settings_bias_uns_get_mean_range <- function(
  ind_batch_list,
  .data,
  pop_gate,
  chnl_cut,
  path_project
) {
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list(
        # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i],
        path_project = path_project
      )
      purrr::map_dbl(ex_list, function(ex) {
        abs(
          diff(
            quantile(
              .get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))],
              c(0.99, 0.01)
            ),
            na.rm = TRUE
          )
        )[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


#' @keywords internal
.complete_chnl_settings_min_bw <- function(bw_min, bias_uns) {
  if (!is.null(bw_min)) {
    return(bw_min)
  }
  bias_uns * 2.25
}

#' @keywords internal
.complete_chnl_settings_cp_min <- function(
  cp_min,
  .data,
  pop_gate,
  chnl_cut,
  ind_batch_list,
  path_project
) {
  if (!is.null(cp_min)) {
    return(cp_min)
  }
  .debug("calculating cp_min automatically") # nolint
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list(
        # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i],
        path_project = path_project
      )
      purrr::map_dbl(ex_list, function(ex) {
        median(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], na.rm = TRUE)[[
          1
        ]] # nolint
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


# Get all cutpoint type names
# Returns character vector of all available cutpoint names
#' @keywords internal
.get_full_cp_type_vec <- function(fdr) {
  # Get cutpoint names for unstim-based cuts
  # cp_name_vec_uns <- .get_cp_uns_name_vec(fdr)
  # cp_name_vec_uns_root <- paste0(cp_name_vec_uns, 'root')

  # output all cutpoint names
  c(
    "man",
    "tg",
    "dcp",
    "midp",
    "scp",
    "uns",
    "unsr",
    "loc"
  )
}

#' @keywords internal
.complete_chnl_settings_save <- function(chnl_list, path_project) {
  path_save <- file.path(path_project, "meta_data", "chnl_list.rds")
  if (file.exists(path_save)) {
    invisible(file.remove(path_save))
  }
  if (!dir.exists(path_project)) {
    dir.create(path_project, recursive = TRUE)
  }
  saveRDS(
    chnl_list,
    file = path_save
  )
}

#' @keywords internal
.chnl_lab <- function(.data) {
  adf <- switch(
    class(.data)[1],
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
#' @return A named list of marker settings (as saved by .complete_chnl_settings_save()).
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
#' stimgate_meta_read_settings_markers(tmp)
#' }
#' @export
stimgate_meta_read_settings_chnls <- function(path_project) {
  path_chnl_list <- file.path(path_project, "meta_data", "chnl_list.rds")
  if (!file.exists(path_chnl_list)) {
    stop("Channel list file not found in project meta_data folder")
  }
  readRDS(path_chnl_list)
}

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
#' stimgate_meta_read_settings_chnls(tmp)
#' }
#' @export
stimgate_meta_read_settings_markers <- function(path_project) {
  marker_list <- stimgate_meta_read_settings_chnls(path_project)
  chnl_lab <- stimgate_meta_read_chnl_lab(path_project)
  # chnl_lab maps channel code -> label; rename marker_list elements using labels
  names(marker_list) <- chnl_lab[names(marker_list)]
  marker_list
}

#' @title Get marker settings for a single channel
#' @description Retrieve the marker settings for a single channel. The function
#'   accepts either a channel label (as returned by stimgate_meta_read_chnl_lab)
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
#' stimgate_meta_read_settings_chnl(tmp, "BC1 label")
#' stimgate_meta_read_settings_chnl(tmp, "BC1")
#' }
#' @export
stimgate_meta_read_settings_chnl <- function(path_project, chnl) {
  chnl_list <- stimgate_meta_read_settings_chnls(path_project)
  # If user supplied the label (names of chnl_list)
  if (!chnl %in% names(chnl_list)) {
    stop(sprintf("Channel %s not found in marker list", chnl))
  }
  chnl_list[[chnl]]
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
#' stimgate_meta_read_settings_marker(tmp, "BC1")
#' }
#' @export
stimgate_meta_read_settings_marker <- function(path_project, marker) {
  marker_list <- stimgate_meta_read_settings_chnls(path_project)
  if (!marker %in% names(marker_list)) {
    stop(sprintf("Marker %s not found in marker list", marker))
  }
  marker_list[[marker]]
}

#' @rdname stimgate_meta_read_lab
#' @title Read channel or marker label mapping
#' @description Read the saved channel label mapping (chnl_lab.rds) from the project's meta_data folder.
#' @param path_project character Path to project.
#' @return Named character vector mapping channel names to labels.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
#' stimgate_meta_read_chnl_lab(tmp)
#' }
#' @export
stimgate_meta_read_chnl_lab <- function(path_project) {
  path_chnl_lab <- file.path(path_project, "meta_data", "chnl_lab.rds")
  if (!file.exists(path_chnl_lab)) {
    stop("Channel label file not found in project meta_data folder")
  }
  readRDS(path_chnl_lab)
}

#' @rdname stimgate_meta_read_lab
stimgate_meta_read_marker_lab <- function(path_project) {
  chnl_lab <- stimgate_meta_read_chnl_lab(path_project)
  stats::setNames(names(chnl_lab), chnl_lab)
}

#' @keywords internal
.save_meta_data <- function(.data, batch_list, path_project) {
  path_dir_meta_data <- file.path(path_project, "meta_data")
  if (!dir.exists(path_dir_meta_data)) {
    dir.create(path_dir_meta_data, recursive = TRUE)
  }
  .save_meta_data_chnl_lab(.data, path_dir_meta_data)
  .save_meta_data_batch_list(batch_list, path_dir_meta_data)
}

#' @keywords internal
.save_meta_data_chnl_lab <- function(.data, path_dir) {
  chnl_lab <- .chnl_lab(.data)
  saveRDS(
    chnl_lab,
    file = file.path(path_dir, "chnl_lab.rds")
  )
}

#' @keywords internal
.save_meta_data_batch_list <- function(batch_list, path_dir) {
  saveRDS(
    batch_list,
    file = file.path(path_dir, "batch_list.rds")
  )
}

#' @title Read batch list from project
#' @description Read the saved batch_list object from the project's meta_data folder.
#' @param path_project character Path to project.
#' @return A list describing sample grouping into batches (as saved by .save_meta_data_batch_list()).
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
#' saveRDS(list(batch1 = c(1, 2)), file.path(tmp, "meta_data", "batch_list.rds"))
#' stimgate_meta_read_batch_list(tmp)
#' }
#' @export
stimgate_meta_read_batch_list <- function(path_project) {
  path_batch_list <- file.path(path_project, "meta_data", "batch_list.rds")
  if (!file.exists(path_batch_list)) {
    stop("Batch list file not found in project meta_data folder")
  }
  readRDS(path_batch_list)
}

.extract_chnl <- function(chnl, marker, path_project) {
  if (!is.null(chnl)) {
    if (!length(chnl) == length(unique(chnl))) {
      stop(
        "Duplicate channel names found in `chnl`. Please ensure that each channel is unique."
      )
    }
    return(chnl)
  }
  marker_lab <- stimgate_meta_read_marker_lab(path_project)
  chnl_vec <- marker_lab[marker] |> stats::setNames(NULL)
  if (!length(chnl_vec) == length(unique(chnl_vec))) {
    stop(
      "Duplicate channel labels found for the specified markers. ",
      "Please ensure that each marker has a unique channel label. ",
      "Otherwise, simply specify `chnl` instead of `marker` ",
      "(and then also `chnl_settings` instead of `marker_settings` if applicable). "
    )
  }
  chnl_vec
}

.extract_chnl_settings <- function(
  chnl_settings,
  marker_settings,
  chnl,
  marker,
  path_project
) {
  .verify_chnl_settings(
    chnl_settings = chnl_settings,
    marker_settings = marker_settings,
    chnl = chnl,
    marker = marker
  )
  if (!is.null(chnl_settings)) {
    return(chnl_settings)
  }
  if (is.null(marker_settings)) {
    return(lapply(chnl, function(x) list()) |> stats::setNames(chnl))
  }
  marker_lab <- stimgate_meta_read_marker_lab(path_project)
  chnl_settings <- marker_settings
  names(chnl_settings) <- marker_lab[names(marker_settings)]
  for (i in seq_along(chnl_settings)) {
    chnl_settings[[i]]$chnl_cut <- marker_lab[[chnl_settings[[i]]$marker_cut]]
  }
  chnl_settings
}
