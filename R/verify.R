#' @keywords internal
.verify_gate_inputs <- function(
  path_project,
  .data,
  batch_list,
  pop_gate,
  chnl,
  marker,
  chnl_settings,
  marker_settings,
  calc_cyt_pos_gates,
  calc_single_pos_gates,
  bw_mtd,
  gate_combn,
  gate_quant,
  bw_ncell_min,
  bw_ncell_max,
  min_cell,
  tol_clust
) {
  # 1. Channel / Marker Mutual Exclusivity Checks
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Specify only one of 'chnl' or 'marker', not both.")
  }
  if (is.null(chnl) && is.null(marker)) {
    stop("Must specify one of 'chnl' or 'marker'.")
  }
  if (!is.null(chnl) && !is.null(marker_settings)) {
    stop("When 'chnl' is specified, 'marker_settings' must be NULL.")
  }
  if (!is.null(marker) && !is.null(chnl_settings)) {
    stop("When 'marker' is specified, 'chnl_settings' must be NULL.")
  }

  # 2. Structural & Type Checks
  if (
    !is.character(path_project) ||
      length(path_project) != 1 ||
      path_project == ""
  ) {
    stop(
      "`path_project` must be a single, non-empty character string specifying a directory."
    )
  }
  if (!is.character(pop_gate) || length(pop_gate) != 1) {
    stop("`pop_gate` must be a single character string (e.g., 'root').")
  }
  if (
    !inherits(
      .data,
      c(
        "GatingSet",
        "GatingHierarchy",
        "flowFrame",
        "flowSet",
        "cytoframe",
        "cytoset"
      )
    )
  ) {
    stop(
      "`.data` must be a valid flow core/workspace object (e.g., GatingSet, flowFrame)."
    )
  }
  if (!is.list(batch_list) || length(batch_list) == 0) {
    stop("`batch_list` must be a non-empty list of sample indices.")
  }

  # 3. Logical & Hyperparameter Checks
  if (!is.logical(calc_cyt_pos_gates) || length(calc_cyt_pos_gates) != 1) {
    stop("`calc_cyt_pos_gates` must be a single logical value (TRUE/FALSE).")
  }
  if (
    !is.logical(calc_single_pos_gates) || length(calc_single_pos_gates) != 1
  ) {
    stop("`calc_single_pos_gates` must be a single logical value (TRUE/FALSE).")
  }

  valid_bw_mtds <- c("nrd0", "sj", "hpi_0", "hpi_1", "hpi_2", "hpi_3")
  if (!bw_mtd %in% valid_bw_mtds) {
    stop(sprintf(
      "`bw_mtd` must be one of: %s",
      paste(valid_bw_mtds, collapse = ", ")
    ))
  }

  valid_gate_combns <- c("min", "median", "max")
  if (!gate_combn %in% valid_gate_combns) {
    stop(sprintf(
      "`gate_combn` must be one of: %s",
      paste(valid_gate_combns, collapse = ", ")
    ))
  }

  if (
    !is.numeric(gate_quant) ||
      length(gate_quant) != 2 ||
      any(gate_quant < 0 | gate_quant > 1)
  ) {
    stop(
      "`gate_quant` must be a numeric vector of two probabilities between 0 and 1."
    )
  }
  if (!is.numeric(bw_ncell_min) || bw_ncell_min <= 0) {
    stop("`bw_ncell_min` must be a positive number.")
  }
  if (!is.numeric(bw_ncell_max) || bw_ncell_max < bw_ncell_min) {
    stop(
      "`bw_ncell_max` must be a positive number greater than or equal to `bw_ncell_min`."
    )
  }
  if (!is.numeric(min_cell) || min_cell <= 0) {
    stop("`min_cell` must be a positive number.")
  }
  if (!is.numeric(tol_clust) || tol_clust <= 0) {
    stop("`tol_clust` must be a positive numeric threshold.")
  }
  chnl_lab <- .chnl_lab(.data)
  if (!is.null(chnl)) {
    all_chnl_valid <- all(chnl %in% names(chnl_lab))
    if (!all_chnl_valid) {
      stop(
        "The following channels are not found in the GatingSet: ",
        paste(setdiff(chnl, names(chnl_lab)), collapse = ", ")
      )
    }
  }
  if (!is.null(marker)) {
    all_marker_valid <- all(marker %in% chnl_lab)
    if (!all_marker_valid) {
      stop(
        "The following markers are not found in the GatingSet: ",
        paste(setdiff(marker, chnl_lab), collapse = ", ")
      )
    }
  }

  invisible(TRUE)
}

.verify_chnl_settings <- function(chnl_settings, chnl, marker_settings, marker) {
  if (!is.null(chnl_settings) && !is.null(marker_settings)) {
    stop("Specify only one of `chnl_settings` or `marker_settings`, not both.")
  }
  if (!is.null(chnl_settings) && !is.list(chnl_settings)) {
    stop("`chnl_settings` must be a list of channel-specific settings.")
  }
  if (!is.null(marker_settings) && !is.list(marker_settings)) {
    stop("`marker_settings` must be a list of marker-specific settings.")
  }
  if (!is.null(chnl_settings) && length(chnl_settings) > 0L) {
    chnl_vec_from_settings <- sapply(chnl_settings, function(x) x$chnl_cut)
    if (!all(chnl_vec_from_settings %in% chnl)) {
      stop(
        "All channels specified in `chnl_settings` must be included in `chnl`"
      )
    }
  } else if (!is.null(marker_settings) && length(marker_settings) > 0L) {
    marker_vec_from_settings <- sapply(marker_settings, function(x) x$marker_cut)
    if (!all(marker_vec_from_settings %in% marker)) {
      stop(
        "All markers specified in `marker_settings` must be included in `marker`"
      )
    }
  }
}

#' @keywords internal
.verify_chnl_settings_chnl <- function(chnl_curr, settings) {
  prefix <- sprintf("Channel '%s' setting error: ", chnl_curr)

  # Validate logic flags if present
  if (
    !is.null(settings$exc_min) &&
      (!is.logical(settings$exc_min) || length(settings$exc_min) != 1)
  ) {
    stop(paste0(prefix, "`exc_min` must be a single logical value."))
  }

  # Validate numeric limits and hyper-parameters
  if (
    !is.null(settings$bias_uns) &&
      (!is.numeric(settings$bias_uns) || length(settings$bias_uns) != 1)
  ) {
    stop(paste0(prefix, "`bias_uns` must be a single numeric value."))
  }
  if (
    !is.null(settings$bias_uns_factor) &&
      (!is.numeric(settings$bias_uns_factor) || settings$bias_uns_factor <= 0)
  ) {
    stop(paste0(prefix, "`bias_uns_factor` must be a positive number."))
  }
  if (!is.null(settings$bw) && (!is.numeric(settings$bw) || settings$bw <= 0)) {
    stop(paste0(prefix, "`bw` must be a positive numeric value."))
  }
  if (
    !is.null(settings$bw_adj) &&
      (!is.numeric(settings$bw_adj) || settings$bw_adj <= 0)
  ) {
    stop(paste0(prefix, "`bw_adj` must be a positive numeric multiplier."))
  }
  if (
    !is.null(settings$cp_min) &&
      (!is.numeric(settings$cp_min) || length(settings$cp_min) != 1)
  ) {
    stop(paste0(prefix, "`cp_min` must be a single numeric value."))
  }
  if (
    !is.null(settings$max_pos_prob_x) &&
      (!is.numeric(settings$max_pos_prob_x) ||
        length(settings$max_pos_prob_x) != 1)
  ) {
    stop(paste0(prefix, "`max_pos_prob_x` must be a single numeric value."))
  }

  invisible(TRUE)
}
