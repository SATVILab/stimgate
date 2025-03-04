#' @title Get list where each element is list of indices for a batch
#'
#' @inheritParams get_cp
#'
#' @return A list, where each element is list of indices for a batch.
.get_ind_batch_list <- function(data,
                                batch_size,
                                fcs = NULL,
                                ind_skip) {
  if (is.null(fcs)) {
    batch_vec <- vapply(
      seq_along(data), function(i) {
        keyword_vec <- flowWorkspace::gh_pop_get_data(data[[i]]) |>
          flowCore::keyword()
        keyword_vec[["GUID"]] |> basename()
      },
      character(1)
    )
  } else {
    batch_vec <- fcs
  }

  batch_vec <- gsub("_p1|_p4|_mtbaux|_ebv|_uns", "", batch_vec)
  split(seq_along(batch_vec), batch_vec) |>
    .get_ind_batch_list_skip(ind_skip = ind_skip)
}

.get_ind_batch_list_skip <- function(ind_batch_list_init,
                                     ind_skip) {
  if (is.null(ind_skip)) {
    return(ind_batch_list_init)
  }
  ind_batch_list <- list()
  k <- 0
  for (i in seq_along(ind_batch_list_init)) {
    any_in_ind_skip <- purrr::map_lgl(
      unlist(ind_skip), function(x) x %in% ind_batch_list_init[[i]]
    ) |>
      any()
    if (!any_in_ind_skip) {
      k <- k + 1
      ind_batch_list <- ind_batch_list |>
        append(list(ind_batch_list_init[[i]]))
    }
  }
  ind_batch_list
}

str_detect_any <- function(string, pattern) {
  purrr::map_lgl(
    pattern,
    function(pattern_curr) stringr::str_detect(string, pattern_curr)
  ) |>
    any()
}

#' @title Get expression matrix from multiple
#'
#' @description
#' Takes a GatingSet and returns the expression data for the
#' columns selected in \code{cut} and \code{high} as a tibble
#' for selected samples in the GatingHierarchy.
#'
#' @inheritParams get_cp
#' @param ind_data numeric vector.
#' Indices in GatingSet for which data should be drawn.
#' @param ind_uns numeric.
#' Index in GatingSet of sample containing unstim measurements for this batch.
#' @param data_name character. Name of GatingSet in workspace.
#'
#' @return A list, where each element is a tibble containing the cell
#' expression measurements from
#' a single sample for the channels named in \code{cut} and \code{high},
#' as well as the following columns:
#' - ind: index of sample in original GatingSet.
.get_ex_list <- function(data, ind_batch, ind_in_batch_gate,
                         ind_in_batch_uns, ind_in_batch_lab_vec,
                         pop, cut, high, data_name) {
  # Preparation
  # -------------------------

  # get indices of data
  ind_data <- ind_batch[c(
    ind_in_batch_gate[!ind_in_batch_gate == ind_in_batch_uns],
    ind_in_batch_uns
  )]
  ind_uns <- ind_batch[ind_in_batch_uns]

  # get the names of the associated fcs files
  purrr::map(ind_data, function(i) {
    ex <- .get_ex(
      data = data[[i]], pop = pop,
      cut = cut,
      high = high,
      ind = i,
      is_uns = ifelse(i == ind_uns, TRUE, FALSE),
      stim = ind_in_batch_lab_vec[which(ind_batch == i)],
      ind_in_batch = which(ind_batch == i),
      data_name = data_name
    ) |>
      tibble::as_tibble()
    ex |>
      dplyr::mutate(high = .check_if_high(ex = .env$ex, .env$high)) # nolint
  }) |>
    stats::setNames(as.character(ind_data))
}

#' @title Get expression matrix
#'
#' @description
#' Takes a GatingHierarchy and returns the expression data
#' for the columns selected in \code{cut} and \code{high} as a tibble.
#'
#' @inheritParams get_cp
#'
#' @return A tibble containing the cell expression measurements
#' from a single sample for the channels named in \code{cut} and \code{high}.
.get_ex <- function(data,
                    pop,
                    cut,
                    high = NULL,
                    ind = NULL,
                    is_uns = NULL,
                    stim = NULL,
                    ind_in_batch = NULL,
                    data_name) {
  # functions
  str_remove_all_v <- function(str, pattern) { # nolint
    for (pattern_curr in pattern) {
      str <- str |> stringr::str_remove_all(pattern_curr)
    }
    str
  }

  # get data
  fr <- flowWorkspace::gh_pop_get_data(data, y = pop)
  ex <- flowCore::exprs(fr)

  # filter ex data here

  # add in stim-related columns
  if (str_detect_any(data_name, c(
    "gs_cytof", "gs_proto", "gs_cd8_base",
    "gs_cytof_acs"
  ))) {
    out_tbl <- tibble::as_tibble(ex) |>
      dplyr::mutate(ind_cell = seq_len(dplyr::n())) |>
      dplyr::select(ind_cell, !!c( # nolint
        cut, names(high), "Ce140Di", "Lu175Di",
        "Eu151Di", "Eu153Di", "Ho165Di"
      ))
  } else {
    out_tbl <- tibble::as_tibble(ex) |>
      dplyr::mutate(ind_cell = seq_len(dplyr::n())) |>
      dplyr::select(ind_cell, !!c(cut, names(high))) # nolint
  }


  if (!is.null(ind)) {
    out_tbl <- out_tbl |> dplyr::mutate(ind = ind)
  }
  if (!is.null(is_uns)) {
    out_tbl <- out_tbl |> dplyr::mutate(is_uns = is_uns)
  }
  if (!is.null(stim)) {
    out_tbl <- out_tbl |> dplyr::mutate(stim = stim)
  }
  if (!is.null(ind_in_batch)) {
    out_tbl <- out_tbl |> dplyr::mutate(ind_in_batch = ind_in_batch)
  }

  out_tbl$cut <- out_tbl[[cut[1]]]

  # get the fcs file name and batch name
  if (str_detect_any(data_name, c(
    "gs_cytof", "gs_proto", "gs_cd8_base",
    "gs_cytof_acs"
  ))) {
    fcs <- flowCore::keyword(fr)[["GUID"]] |> basename()

    stim_vec <- c("ebv", "mtbaux", "p1", "p4", "uns")
    stim_ind_vec <- purrr::map_lgl(
      paste0("_", stim_vec),
      function(stim) {
        stringr::str_detect(fcs, stim)
      }
    )
    stim <- stim_vec[stim_ind_vec]

    underscore_loc_tbl <- stringr::str_locate_all(fcs, "_")[[1]]
    second_underscore_loc <- underscore_loc_tbl[2, "end"][[1]]
    batch_sh <- stringr::str_sub(fcs, end = second_underscore_loc - 1)

    # get batch and sample
    batch <- flowCore::keyword(fr)[["GUID"]] |>
      basename() |>
      stringr::str_remove("_ebv") |>
      stringr::str_remove("_mtbaux") |>
      stringr::str_remove("_p1") |>
      stringr::str_remove("_p4") |>
      stringr::str_remove("_uns") |>
      stringr::str_remove_all("_") |>
      stringr::str_remove_all("-") |>
      stringr::str_remove_all(" ") |>
      stringr::str_remove_all(".fcs") |>
      stringr::str_remove("2.1GB") |>
      stringr::str_remove("normalized") |>
      stringr::str_remove("timecut") |>
      stringr::str_remove("debeaded[2]?") |>
      stringr::str_remove("pid[12]")
  }


  # re-arrange
  out_tbl |>
    dplyr::mutate(batch = batch, batch_sh = batch_sh, fcs = fcs, stim = stim) |>
    dplyr::select(
      batch, batch_sh, fcs, stim, ind_cell, # nolint
      !!c(cut, names(high)), cut, dplyr::everything()
    )
}




#' @title Winsorises a numeric vector
#'
#' @description
#'
#' @param num numeric
#' vector. Contains data to winsorise.
#' @param min numeric.
#' Value to change values below this to.
#' @param quant numeric.
#' Quantile to select value to change values below this to.
.wins_vec <- function(num, min, quant) {
  if (missing(min)) min <- quantile(num, quant)
  vapply(num, function(x) ifelse(x > min, x, min), numeric(1))
}

#' @title Check each cell for being high on
#' at least one other marker
#'
#' @description
#' Indicates for each cell in a set whether
#' it has a high expression levelf or at least one other marker.
#'
#' @inheritParams get_cp
#' @param ex dataframe A dataframe containing the expression data.
#'
#' @return A logical vector with \code{TRUE} in a given
#' index if the corresponding cell was high for at least one of the markers
#' specified by the \code{high} parameter in the sense specified by \code{high}.
.check_if_high <- function(ex, high) {
  high_ind_vec <- rep(FALSE, nrow(ex))
  for (i in seq_along(high)) {
    chnl_high_vec <- ex[[names(high)[i]]]
    high_ind_vec <- high_ind_vec | chnl_high_vec > high[i]
  }
  high_ind_vec
}

#' @title Create and get the base directory
#'
#' @description
#'
#' @inheritParams save_plot_cp
#' @param dir_base_init character.
#' If provided, then this is the base directory
#' in which the plots are saved.
#' If \code{NULL}, then the base directory is
#' set to the package directory (i.e. DataPackageR::project_path()).
#' @param empty_dir logical.
#' If \code{TRUE}, then old contents from the directory the
#' files are saved to are deleted. Default is \code{FALSE}.
#'
#' @return String of the base directory
stimgate_dir_base_create <- function(params,
                                     dir_base_init = NULL,
                                     empty_dir = FALSE) {
  # create a base directory if not provided
  if (is.null(dir_base_init)) dir_base_init <- here::here("data-raw/gating")

  # create sub-folder to save results to
  cut <- params$cut
  chnl_lab_vec <- params$chnl_lab
  cut_marker <- chnl_lab_vec[[cut]]
  pop_gate <- stringr::str_replace_all(params$pop_gate, "[[//]]", "")
  pop_gate <- pop_gate |> stringr::str_remove_all(" ")
  while (stringr::str_sub(pop_gate) == "_") {
    pop_gate <- stringr::str_sub(pop_gate, 2)
  }
  data_name <- params$data_name
  batch <- params$ind_in_batch_lab_vec[params$ind_in_batch_gate] |>
    paste0(collapse = "")

  dir_base_save <- file.path(
    dir_base_init, data_name,
    pop_gate, batch, cut_marker
  )

  # empty the directory if it already exists and empty_dir == TRUE
  if (dir.exists(dir_base_save)) {
    if (empty_dir) unlink(dir_base_save, recursive = FALSE)
  }
  # create the folder if it does not exist
  if (!dir.exists(dir_base_save)) dir.create(dir_base_save, recursive = TRUE)

  # return the path to the folder
  dir_base_save
}

.spread_nearby_obs_tbl <- function(gate_tbl, width) {
  gate_tbl_adj <- purrr::map_df(
    split(gate_tbl, gate_tbl$gate_combn == "no"),
    function(gate_tbl_curr) {
      if (gate_tbl_curr$gate_combn[1] == "no") {
        return(gate_tbl_curr)
      }

      gate_tbl_curr |>
        dplyr::group_by(gate_name, batch) |> # nolint
        dplyr::slice(1) |>
        dplyr::ungroup()
    }
  )

  # calculate adjusted gates
  gate_tbl_adj <- gate_tbl_adj |>
    dplyr::group_by(batch) |> # nolint
    dplyr::arrange(batch, gate) |> # nolint
    dplyr::mutate(gate_adj = .spread_nearby_obs(
      num = gate, width = width, sort = FALSE, batch = .data$batch # nolint
    )) |>
    dplyr::ungroup()

  # get names for adjusted gates
  gate_tbl_adj <- gate_tbl_adj |>
    dplyr::mutate(name = purrr::map_chr(seq_len(dplyr::n()), function(i) {
      gate_combn <- gate_tbl_adj$gate_combn[i]
      gate_name <- gate_tbl_adj$gate_name[i]
      batch <- gate_tbl_adj$batch[i]
      ind <- gate_tbl_adj$ind[i]
      if (gate_combn == "no") {
        return(paste0(gate_name, batch, ind))
      }
      paste0(gate_name, batch)
    })) |>
    dplyr::group_by(gate_adj) |> # nolint
    dplyr::mutate(gate_adj_2 = purrr::map_dbl(dplyr::n(), function(N) {
      if (N == 1) {
        return(gate_adj) # nolint
      }
      gate_adj + runif(N, -width / 4, max = width / 4) # nolint
    })) |>
    dplyr::ungroup()

  gate_lab_vec_adj <- stats::setNames(
    gate_tbl_adj$gate_adj, gate_tbl_adj$name
  ) +
    runif(nrow(gate_tbl_adj), min = -width / 4, max = width / 4)

  gate_tbl |>
    dplyr::mutate(
      name = purrr::map_chr(seq_len(dplyr::n()), function(i) {
        gate_combn <- gate_tbl$gate_combn[i]
        gate_name <- gate_tbl$gate_name[i]
        batch <- gate_tbl$batch[i]
        ind <- gate_tbl$ind[i]
        if (gate_combn == "no") {
          return(paste0(gate_name, batch, ind))
        }
        paste0(gate_name, batch)
      }),
      gate_adj = gate_lab_vec_adj[.data$name] # nolint
    ) |>
    dplyr::select(-name) |>
    dplyr::rename(
      gate_orig = gate, # nolint
      gate = gate_adj # nolint
    )
}

#' @title Spread apart obs that are too close to each other
.spread_nearby_obs <- function(num, width = 0.25, sort = FALSE, batch = NULL) {
  if (sort) num <- sort(num)
  n <- length(num)
  for (i in seq_along(num)) {
    if (i == n) next
    val <- num[i]
    # get maximum index close to this one
    upper_dist_vec <- abs(val - num[(i + 1):n])
    if (any(upper_dist_vec < width)) {
      num[i] <- num[i] - width
    }
    if (i == 1) next
    lower_dist_vec <- abs(num[i] - num[1:(i - 1)])
    if (any(lower_dist_vec < width)) {
      num[which(lower_dist_vec < width)] <-
        num[which(lower_dist_vec < width)] - width
    }
  }
  dist_vec <- dist(num)[lower.tri(dist(num))]
  dist_vec <- dist_vec[!is.na(dist_vec)]
  if (any(dist_vec < width)) num <- .spread_nearby_obs(num)
  num

  purrr::map(split(num, num), function(x) {
    if (length(x) == 1) x
    x + runif(length(x), min = -width / 4, max = width / 4)
  }) |>
    unlist() |>
    stats::setNames(NULL)
}

.ensure_cytoutils <- function() {
  if (requireNamespace("cytoUtils", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    BiocManager::install("remotes")
  }
  remotes::install_github("RGLab/cytoUtils")
  invisible(TRUE)
}

.ensure_flowstats <- function() {
  if (requireNamespace("flowStats", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("flowStats")
}
.get_batch_from_ind <- function(ind,
                                ind_batch_list) {
  ind_batch_list[[
    which(
      purrr::map_lgl(ind_batch_list, function(x) ind %in% x)
    )
  ]]
}
