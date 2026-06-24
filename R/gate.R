#' @title Identify cytokine-positive cells through automated gating
#'
#' @description
#' Main function for identifying cytokine-positive cells using outlier-based gating
#' to compare stimulated versus unstimulated samples. This function implements a
#' comprehensive workflow that identifies cells responding to stimulation by
#' detecting outliers in cytokine expression distributions. The process includes
#' density estimation, threshold identification, clustering-based gate refinement,
#' and generation of comprehensive statistics and visualizations.
#'
#' The function operates by comparing cytokine expression in stimulated samples
#' against corresponding unstimulated controls from the same donor/batch to identify
#' cells that have likely responded to stimulation. It accounts for batch effects,
#' background cytokine production, and technical variability.
#'
#' @param path_project character. Path to project directory where all results will be saved.
#'   This directory will contain subdirectories for each marker with gate tables,
#'   statistics, and plots. The directory will be created if it doesn't exist.
#' @param .data GatingSet. A flowWorkspace GatingSet object containing the flow cytometry
#'   data with both stimulated and unstimulated samples. The GatingSet should have
#'   consistent channel names across all samples and include proper sample annotations.
#' @param batch_list list. Named list where each element contains indices of samples
#'   belonging to the same batch/donor. Names will be used for batch identification.
#'   Example: list(donor1 = 1:10, donor2 = 11:20). Proper batching is crucial for
#'   accurate background subtraction and gate identification.
#' @param chnl list. List where each element specifies parameters for gating a
#'   specific channel Each element should be a list containing at minimum the channel
#'   name (e.g., list(cut = "IL2")). Additional channel-specific parameters can
#'   override global defaults. The marker name should match channel names in the GatingSet.
#' @param marker character vector. Alternative way to specify markers to gate on.
#'   When provided, this is used instead of chnl to determine which markers to analyze.
#'   Default is NULL.
#' @param pop_gate character vector. Population(s) within which to perform gating.
#'   Default is "root" to gate on all cells. Can specify other populations like
#'   "CD3+" or "CD4+" if these gates already exist in the GatingSet.
#' @param bias_uns numeric. Bias adjustment for unstimulated samples to account for
#'   background cytokine production. When NULL (default), no bias correction is applied.
#'   Positive values shift the unstimulated distribution higher, making gates more
#'   conservative.
#' @param bias_uns_factor numeric. Multiplicative factor applied to bias_uns.
#'   Default is 1. Values > 1 increase the bias effect, values < 1 decrease it.
#'   This provides fine-tuning of the bias correction.
#' @param exc_min logical. Whether to exclude minimum expression values during
#'   analysis. Default is TRUE. Minimum values often represent technical artifacts
#'   or compensation spillover and should typically be excluded.
#' @param cp_min numeric. Minimum allowable cutpoint value. When NULL (default),
#'   no minimum is enforced. Useful for ensuring gates don't fall below known
#'   technical thresholds or background levels.
#' @param bw numeric. Specify the bandwith for density estimation. When NULL (default), bandwidth is estimated automatically. Default is `NULL`.
#' @param bw_min numeric. Minimum bandwidth for density estimation. Ignored if `bw` is set. Default is `NULL`.
#' @param bw_max numeric. Maximum bandwidth for density estimation. Ignored if `bw` is set. Default is `NULL`.
#' @param bw_mtd character. Method for automated bandwidth selection. Options include "nrd0", "sj", "hpi_0", "hpi_1", "hpi_2" and "hpi_3", which corresponds to the Silverman rule of thumbg (`"nrd0"`), the Sheather-Jones plug-in estimator (`"sj"`) and the Wand & Jones plugin-estimator for the 0-th, 1st, 2nd and 3rd derivatives of the density (`"hpi_1"`, `"hpi_1"`, `"hpi_2"` and `"hpi_3"`). Default is "nrd0". Ignored if `bw` is set. Default is `"hpi_1"`.
#' @param bw_adj numeric. Adjustment factor for bandwidth. Default is 1. Ignored if `bw` is set. Default is 1.
#' @param bw_ncell_min numeric. Minimum number of cells required for bandwidth estimation. If a sample has fewer cells than `bw_ncell_min`, cells are sampled with replacement to reach the minimum, with noise subsequently added. Ignored if `bw` is set. Default is 100.
#' @param bw_ncell_max numeric. Maximum number of cells used for bandwidth estimation. If a sample has more cells than `bw_ncell_max`, cells are sampled without replacement to reach the maximum. Ignored if `bw` is set. Default is 100 000.
#' @param min_cell numeric. Minimum number of cells required for reliable gating.
#'   Default is 100. Samples with fewer cells will be skipped as they don't provide
#'   sufficient statistical power for accurate gate identification.
#' @param max_pos_prob_x numeric. Maximum x-value (expression level) to consider
#'   when calculating positive probabilities. Default is Inf (no limit). Can be used
#'   to exclude extremely high expression values that may represent doublets or artifacts.
#' @param gate_quant numeric vector. Quantiles used for gate combination when multiple
#'   gates are identified. Default is c(0.25, 0.75). The method specified in gate_combn
#'   determines how these quantiles are used (e.g., minimum of 25th percentiles).
#' @param tol_clust numeric. Convergence tolerance for clustering algorithms used in
#'   gate refinement. Default is 1e-7. Smaller values require more precise convergence
#'   but may increase computation time.
#' @param gate_combn character. Method for combining multiple gate candidates.
#'   Default is "min" to use the most conservative (lowest) gate. Other options may
#'   include "median" or "max" depending on the desired stringency.
#' @param marker_settings list. Optional list of additional marker-specific settings
#'   that override global defaults. Each element should be named with the marker name
#'   and contain parameter overrides. Default is NULL.
#' @param chnl_settings list. Optional list of channel-specific settings that override
#'   global defaults. Similar to marker_settings but keyed by channel names.
#'   Default is NULL.
#' @param calc_cyt_pos_gates logical. Whether to calculate refined cytokine-positive
#'   gates using more sophisticated algorithms. Default is TRUE. When FALSE, only
#'   basic gates are calculated, which may be less accurate but faster.
#' @param calc_single_pos_gates logical. Whether to calculate single-positive gates
#'   for individual markers in addition to combination gates. Default is FALSE.
#'   Useful for detailed analysis of individual marker responses.
#' @return character. Returns the path to the project directory where all results
#'   have been saved. The directory structure created includes:
#'   \itemize{
#'     \item \code{path_project/[marker_name]/}: Directory for each marker containing:
#'     \item \code{gate_tbl_init.rds}: Initial gate table with preliminary gates
#'     \item \code{gate_tbl.rds}: Final refined gate table
#'     \item \code{stats/}: Directory containing statistics files
#'     \item \code{plots/}: Directory containing visualization plots (if generated)
#'   }
#'
#' @details
#' The function implements a multi-step workflow for identifying cytokine-positive cells:
#'
#' \strong{Step 1: Data Preparation}
#' \itemize{
#'   \item Validates input parameters and GatingSet structure
#'   \item Completes marker specifications with default values
#'   \item Organizes samples by batch for proper background subtraction
#' }
#'
#' \strong{Step 2: Initial Gate Identification}
#' \itemize{
#'   \item Extracts expression data for each marker within specified populations
#'   \item Estimates probability densities for stimulated and unstimulated samples
#'   \item Identifies candidate cutpoints using outlier detection algorithms
#'   \item Applies clustering to refine gate positions across batches
#' }
#'
#' \strong{Step 3: Cytokine-Positive Gate Refinement (if calc_cyt_pos_gates = TRUE)}
#' \itemize{
#'   \item Applies more sophisticated algorithms to refine initial gates
#'   \item Accounts for background cytokine production and technical variability
#'   \item Optimizes gates to minimize false positives while maintaining sensitivity
#' }
#'
#' \strong{Step 4: Single-Positive Gates (if calc_single_pos_gates = TRUE)}
#' \itemize{
#'   \item Calculates gates for individual markers independent of other markers
#'   \item Useful for understanding single-marker responses
#' }
#'
#' \strong{Step 5: Statistics Generation}
#' \itemize{
#'   \item Calculates comprehensive statistics including frequencies and combinations
#'   \item Generates cross-tabulations of cytokine-positive populations
#'   \item Saves results in structured format for downstream analysis
#' }
#'
#'
#' @seealso
#' \code{\link{stimgate_gate_get}} for extracting gate information,
#' \code{\link{get_stats}} for generating statistics from results,
#' \code{\link{stimgate_plot}} for visualizing identified gates,
#' \code{\link{stimgate_fcs_write}} for exporting cytokine-positive cells,
#' \code{\link[flowWorkspace]{GatingSet}} for GatingSet documentation
#'
#' @examples{
#' example_data <- get_example_data()
#' gs <- flowWorkspace::load_gs(example_data$path_gs)
#' path_project <- file.path(tempdir(), "demonstration")
#'
#' # Run gating
#' stimgate::stimgate_gate(
#'   .data = gs,
#'   path_project = path_project,
#'   pop_gate = "root",
#'   batch_list = example_data$batch_list,
#'   marker = example_data$marker
#' )
#'
#' # Create plots
#' plots <- stimgate_plot(
#'   ind = example_data$batch_list[[1]], # indices in `gs` to plot
#'   .data = gs, # GatingSet
#'   path_project = path_project,
#'   marker = example_data$marker,
#'   grid = TRUE
#' )
#'
#' # Advanced usage with parameter customization
#'
#' }
#' @importFrom flowCore exprs<- parameters<-
#' @importFrom stats approx as.formula binomial density glm kmeans median optim predict quantile rnorm sd dnorm fft
#' @importFrom utils read.csv write.csv
#' @importFrom cluster clusGap maxSE
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth geom_vline geom_hline
#' @importFrom graphics points
#' @importFrom dplyr everything
#' @export
stimgate_gate <- function(
  path_project,
  .data,
  pop_gate = "root",
  batch_list,
  chnl = NULL,
  marker = NULL,
  calc_cyt_pos_gates = TRUE,
  calc_single_pos_gates = FALSE,
  bias_uns = NULL,
  bias_uns_factor = 1,
  exc_min = TRUE,
  cp_min = NULL,
  bw = NULL,
  bw_min = NULL,
  bw_max = NULL,
  bw_mtd = "hpi_1",
  bw_adj = 1,
  bw_ncell_min = 1e2,
  bw_ncell_max = 1e5,
  min_cell = 1e2,
  max_pos_prob_x = Inf,
  gate_quant = c(0.25, 0.75),
  tol_clust = 1e-7,
  gate_combn = "min",
  marker_settings = NULL,
  chnl_settings = NULL
) {
  force(.data)
  if (Sys.getenv("STIMGATE_DEBUG") == "") {
    Sys.setenv("STIMGATE_DEBUG" = "FALSE")
    on.exit(Sys.unsetenv("STIMGATE_DEBUG"), add = TRUE)
  }
  if (Sys.getenv("STIMGATE_DEBUG") == "TRUE") {
    .debug_file_create()
    path_debug <- .debug_file_create()
    message(paste0("Saving debug output to ", path_debug))
    message(
      "Can copy it after the run to working directory with stimgate_debug_copy()"
    ) # nolint
    message(
      "Can print the output after the run to console with stimgate_debug_print()"
    ) # nolint
  }

  # Verify global function inputs
  .verify_gate_inputs(
    path_project = path_project,
    .data = .data,
    batch_list = batch_list,
    pop_gate = pop_gate,
    chnl = chnl,
    marker = marker,
    chnl_settings = chnl_settings,
    marker_settings = marker_settings,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    bw_mtd = bw_mtd,
    gate_combn = gate_combn,
    gate_quant = gate_quant,
    bw_ncell_min = bw_ncell_min,
    bw_ncell_max = bw_ncell_max,
    min_cell = min_cell,
    tol_clust = tol_clust
  )

  if (is.null(names(batch_list))) {
    batch_list <- batch_list |>
      stats::setNames(paste0("batch_", seq_along(batch_list)))
  }

  # get unspecified levels in marker elements
  .save_meta_data(.data, batch_list, path_project)
  chnl <- .extract_chnl(chnl, marker, path_project)

  chnl_settings <- .complete_chnl_settings(
    chnl = chnl,
    marker = marker,
    chnl_settings = chnl_settings,
    marker_settings = marker_settings,
    bias_uns = bias_uns,
    bias_uns_factor = bias_uns_factor,
    exc_min = exc_min,
    .data = .data,
    pop_gate = pop_gate,
    ind_batch_list = batch_list,
    bw_min = bw_min,
    bw_max = bw_max,
    bw = bw,
    bw_mtd = bw_mtd,
    bw_adj = bw_adj,
    bw_ncell_min = bw_ncell_min,
    bw_ncell_max = bw_ncell_max,
    cp_min = cp_min,
    min_cell = min_cell,
    tol_clust = tol_clust,
    max_pos_prob_x = max_pos_prob_x,
    gate_combn = gate_combn,
    gate_quant = gate_quant,
    path_project = path_project
  )

  # inital gates
  .gate_init(
    pop_gate = pop_gate,
    chnl_settings = chnl_settings,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project
  )

  # cytokine-positive gates
  gate_tbl <- .gate_cyt_pos(
    chnl_settings = chnl_settings,
    ind_batch_list = batch_list,
    pop_gate = pop_gate,
    .data = .data,
    calc_cyt_pos = calc_cyt_pos_gates,
    stage = "cyt_pos",
    path_project = path_project
  )

  # single-positive gates
  .gate_single(
    pop_gate = pop_gate,
    chnl_settings = chnl_settings,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    gate_tbl = gate_tbl
  )

  message("")
  message("")
  message("")
  message("getting cyt combn frequencies")

  .gate_stats(
    .data = .data,
    gate_tbl = gate_tbl,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    pop_gate = pop_gate,
    chnl_settings = chnl_settings,
    ind_batch_list = batch_list,
    path_project = path_project
  )

  path_project
}

#' @keywords internal
.gate_init <- function(
  pop_gate,
  chnl_settings,
  .data,
  ind_batch_list,
  path_project
) {
  message("----")
  message("getting base gates")
  message("----")
  message("")

  purrr::walk(chnl_settings, function(chnl_settings_curr) {
    txt <- paste0("chnl: ", chnl_settings_curr$chnl_cut)
    message(txt)

    gate_obj <- .gate_chnl(
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate,
      chnl_settings = chnl_settings_curr,
      path_project = path_project,
      stage = "init"
    )

    .gate_init_save(
      path_project = path_project,
      pop = pop_gate,
      chnl_cut = chnl_settings_curr$chnl_cut,
      gate_tbl = gate_obj$gate_tbl
    )
  })
}

#' @keywords internal
.gate_init_save <- function(path_project, pop, chnl_cut, gate_tbl) {
  path_save <- .gates_get_path_all(
    path_project = path_project,
    pop = pop,
    chnl_cut = chnl_cut,
    init = TRUE
  )
  if (!dir.exists(dirname(path_save))) {
    dir.create(dirname(path_save), recursive = TRUE, showWarnings = TRUE)
  }
  saveRDS(gate_tbl, path_save)
}

#' @keywords internal
.gate_single <- function(
  pop_gate,
  chnl_settings,
  .data,
  ind_batch_list,
  path_project,
  calc_cyt_pos_gates,
  calc_single_pos_gates,
  gate_tbl
) {
  message("")
  message("")
  message("----")
  message("getting single+ gates")
  message("----")
  message("")

  if (!calc_single_pos_gates) {
    .debug("Not gating single-pos gates")
    purrr::walk(chnl_settings, function(chnl_settings_curr) {
      path_save <- .gates_get_path_all(
        path_project = path_project,
        pop = pop_gate,
        chnl_cut = chnl_settings_curr$chnl_cut,
        init = FALSE
      )
      if (!dir.exists(dirname(path_save))) {
        dir.create(dirname(path_save), recursive = TRUE, showWarnings = TRUE)
      }
      saveRDS(
        gate_tbl |>
          dplyr::filter(chnl == chnl_settings_curr$chnl_cut) |>
          dplyr::mutate(gate_single = gate),
        path_save
      )
    })
    return(invisible(TRUE))
  } else {
    .debug("Gating single-pos gates")
  }

  purrr::walk(chnl_settings, function(chnl_settings_curr) {
    txt <- paste0("chnl: ", chnl_settings_curr$chnl_cut)
    message(txt)

    gate_obj <- .gate_chnl(
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate,
      chnl_cut = chnl_settings_curr$chnl_cut,
      chnl_settings = chnl_settings_curr,
      gate_tbl = gate_tbl,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      stage = "single"
    )

    .gate_single_save(
      path_project = path_project,
      pop = pop_gate,
      chnl = chnl_settings_curr$chnl_cut,
      gate_tbl = gate_obj$gate_tbl
    )
  })
}

#' @keywords internal
.gate_single_save <- function(path_project, pop, chnl, gate_tbl) {
  path_save <- .gates_get_path_all(
    path_project = path_project,
    pop = pop,
    chnl_cut = chnl,
    init = FALSE
  )
  if (!dir.exists(dirname(path_save))) {
    dir.create(dirname(path_save), recursive = TRUE, showWarnings = TRUE)
  }
  saveRDS(gate_tbl, path_save)
}

#' @keywords internal
.gate_stats <- function(
  .data,
  gate_tbl = NULL,
  calc_cyt_pos_gates,
  calc_single_pos_gates,
  pop_gate,
  chnl_settings,
  ind_batch_list,
  path_project
) {
  force(.data)
  .get_stats(
    gate_tbl = gate_tbl,
    filter_other_cyt_pos = FALSE,
    combn = TRUE,
    gate_type_cyt_pos_filter = if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_filter = if (calc_single_pos_gates) {
      "single"
    } else {
      "base"
    },
    gate_type_cyt_pos_calc = if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_calc = if (calc_single_pos_gates) "single" else "base",
    save = TRUE,
    pop_gate = pop_gate,
    chnl = purrr::map_chr(chnl_settings, function(x) x$chnl_cut),
    ind_batch_list = ind_batch_list,
    .data = .data,
    save_gate_tbl = TRUE,
    path_project = path_project
  )
}
