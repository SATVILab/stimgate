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
#' @param marker list. List where each element specifies parameters for gating a 
#'   specific marker. Each element should be a list containing at minimum the channel 
#'   name (e.g., list(cut = "IL2")). Additional marker-specific parameters can 
#'   override global defaults. The marker name should match channel names in the GatingSet.
#' @param pop_gate character vector. Population(s) within which to perform gating.
#'   Default is "root" to gate on all cells. Can specify other populations like 
#'   "CD3+" or "CD4+" if these gates already exist in the GatingSet.
#' @param bias_uns numeric. Bias adjustment for unstimulated samples to account for 
#'   background cytokine production. When NULL (default), no bias correction is applied.
#'   Positive values shift the unstimulated distribution higher, making gates more 
#'   conservative. Typically ranges from 0.1 to 1.0 when used.
#' @param bias_uns_factor numeric. Multiplicative factor applied to bias_uns. 
#'   Default is 1. Values > 1 increase the bias effect, values < 1 decrease it.
#'   This provides fine-tuning of the bias correction.
#' @param exc_min logical. Whether to exclude minimum expression values during 
#'   analysis. Default is TRUE. Minimum values often represent technical artifacts 
#'   or compensation spillover and should typically be excluded.
#' @param cp_min numeric. Minimum allowable cutpoint value. When NULL (default), 
#'   no minimum is enforced. Useful for ensuring gates don't fall below known 
#'   technical thresholds or background levels.
#' @param bw_min numeric. Minimum bandwidth for density estimation. When NULL (default),
#'   bandwidth is estimated automatically. Smaller values create more detailed density 
#'   estimates but may be noisier. Typical range is 0.01 to 0.1 on log-transformed data.
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
#' @param calc_cyt_pos_gates logical. Whether to calculate refined cytokine-positive 
#'   gates using more sophisticated algorithms. Default is TRUE. When FALSE, only 
#'   basic gates are calculated, which may be less accurate but faster.
#' @param calc_single_pos_gates logical. Whether to calculate single-positive gates 
#'   for individual markers in addition to combination gates. Default is FALSE. 
#'   Useful for detailed analysis of individual marker responses.
#' @param debug logical. Whether to enable detailed debug output and save intermediate 
#'   results. Default is FALSE. When TRUE, additional files and verbose output are 
#'   generated, useful for troubleshooting and method development.
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
#' \strong{Important Considerations:}
#' \itemize{
#'   \item Ensure stimulated and unstimulated samples are properly paired by batch
#'   \item Channel names in marker specifications must match GatingSet channels exactly
#'   \item Sufficient cell numbers (min_cell) are crucial for reliable gate identification
#'   \item Background bias correction (bias_uns) should be used cautiously and validated
#'   \item Debug mode generates extensive output useful for method validation
#' }
#'
#' @seealso 
#' \code{\link{get_gate_tbl}} for extracting gate information,
#' \code{\link{get_stats}} for generating statistics from results,
#' \code{\link{stimgate_plot}} for visualizing identified gates,
#' \code{\link{stimgate_fcs_write}} for exporting cytokine-positive cells,
#' \code{\link[flowWorkspace]{GatingSet}} for GatingSet documentation
#'
#' @examples
#' \donttest{
#'   # Load required libraries
#'   library(flowWorkspace)
#'   library(stimgate)
#'   
#'   # Basic usage with cytokine markers
#'   # Assuming 'gs' is a GatingSet with samples 1-20
#'   batch_list <- list(
#'     donor1 = c(1:3, 11:13),   # unstim=1:3, stim=11:13 for donor 1
#'     donor2 = c(4:6, 14:16),   # unstim=4:6, stim=14:16 for donor 2  
#'     donor3 = c(7:9, 17:19)    # unstim=7:9, stim=17:19 for donor 3
#'   )
#'   
#'   # Define markers to gate
#'   markers <- list(
#'     list(cut = "IL2"),        # Basic IL2 gating
#'     list(cut = "TNFa"),       # Basic TNF-alpha gating
#'     list(cut = "IFNg")        # Basic IFN-gamma gating
#'   )
#'   
#'   # Run basic gating
#'   result_path <- stimgate_gate(
#'     path_project = "/path/to/results",
#'     .data = gs,
#'     batch_list = batch_list,
#'     marker = markers
#'   )
#'   
#'   # Advanced usage with parameter customization
#'   markers_advanced <- list(
#'     list(cut = "IL2", tol = 0.5e-8),     # Custom tolerance for IL2
#'     list(cut = "TNFa", bias_uns = 0.2),  # Background bias for TNF-alpha
#'     list(cut = "IFNg", cp_min = 0.1)     # Minimum cutpoint for IFN-gamma
#'   )
#'   
#'   result_path <- stimgate_gate(
#'     path_project = "/path/to/results_advanced",
#'     .data = gs,
#'     batch_list = batch_list,
#'     marker = markers_advanced,
#'     pop_gate = "CD3+",                    # Gate within CD3+ population
#'     bias_uns = 0.1,                      # Global background bias
#'     min_cell = 200,                      # Require more cells for gating
#'     calc_cyt_pos_gates = TRUE,           # Use refined gating algorithms
#'     calc_single_pos_gates = TRUE,        # Calculate single-positive gates
#'     debug = TRUE                         # Enable debug output
#'   )
#'   
#'   # Example with quality control parameters
#'   result_path <- stimgate_gate(
#'     path_project = "/path/to/results_qc",
#'     .data = gs,
#'     batch_list = batch_list,
#'     marker = markers,
#'     exc_min = TRUE,                      # Exclude minimum values
#'     bw_min = 0.05,                       # Minimum bandwidth for density
#'     max_pos_prob_x = 6,                  # Exclude very high expression
#'     gate_quant = c(0.1, 0.9),           # More extreme quantiles
#'     gate_combn = "min"                   # Conservative gate combination
#'   )
#'   
#'   # Access results
#'   gates <- get_gate_tbl(result_path)     # Get gate table
#'   stats <- get_stats(result_path)        # Get statistics
#'   
#'   # Create visualizations
#'   plots <- stimgate_plot(
#'     ind = 1:3,
#'     .data = gs,
#'     path_project = result_path,
#'     marker = c("IL2", "TNFa")
#'   )
#' }
#' @importFrom flowCore exprs<- parameters<-
#' @importFrom stats approx as.formula binomial density glm kmeans median optim predict quantile rnorm sd
#' @importFrom utils read.csv write.csv
#' @importFrom cluster clusGap maxSE
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth geom_vline geom_hline
#' @importFrom dplyr everything
#' @export
stimgate_gate <- function(path_project,
                          .data,
                          batch_list,
                          marker,
                          pop_gate = "root",
                          bias_uns = NULL,
                          bias_uns_factor = 1,
                          exc_min = TRUE,
                          cp_min = NULL,
                          bw_min = NULL,
                          min_cell = 1e2,
                          max_pos_prob_x = Inf,
                          gate_quant = c(0.25, 0.75),
                          tol_clust = 1e-7,
                          gate_combn = "min",
                          marker_settings = NULL,
                          calc_cyt_pos_gates = TRUE,
                          calc_single_pos_gates = FALSE,
                          debug = FALSE) {
  force(.data)
  # capture and force-evaluate the .debug flag into a local .debug object
  .debug <- debug

  if (is.null(names(batch_list))) {
    batch_list <- batch_list |>
      stats::setNames(paste0("batch_", seq_along(batch_list)))
  }

  # get unspecified levels in marker elements
  marker <- .complete_marker_list( # nolint
    marker = marker,
    bias_uns = bias_uns,
    bias_uns_factor = bias_uns_factor,
    exc_min = exc_min,
    .data = .data,
    pop_gate = pop_gate,
    ind_batch_list = batch_list,
    bw_min = bw_min,
    cp_min = cp_min,
    min_cell = min_cell,
    tol_clust = tol_clust,
    max_pos_prob_x = max_pos_prob_x,
    gate_combn = gate_combn,
    marker_settings = marker_settings,
    path_project = path_project,
    .debug = .debug
  )

  # inital gates
  .gate_init(
    pop_gate = pop_gate,
    marker = marker,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project,
    noise_sd = NULL,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_clust = tol_clust,
    tol_gate_single = tol_clust * 1e-1,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    .debug = .debug
  )

  # cytokine-positive gates
  gate_tbl <- .gate_cyt_pos( # nolint
    marker_list = marker,
    ind_batch_list = batch_list,
    pop_gate = pop_gate,
    .data = .data,
    calc_cyt_pos = calc_cyt_pos_gates,
    .debug = .debug,
    path_project = path_project
  )

  # single-positive gates
  .gate_single(
    pop_gate = pop_gate,
    marker = marker,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project,
    noise_sd = NULL,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_clust = tol_clust,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    .debug = .debug,
    gate_tbl = gate_tbl
  )

  message("")
  message("")
  message("")
  message("getting cyt combn frequencies")

  path_dir_stats <- .gate_stats(
    .data = .data,
    params = NULL,
    gate_tbl = NULL,
    filter_other_cyt_pos = FALSE,
    combn = TRUE,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    .debug = .debug,
    save = TRUE,
    pop_gate = pop_gate,
    marker = marker,
    ind_batch_list = batch_list,
    path_project = path_project,
    tol_clust = tol_clust,
    save_gate_tbl = TRUE
  )

  path_project
}

.gate_init <- function(pop_gate,
                       marker,
                       .data,
                       ind_batch_list,
                       path_project,
                       noise_sd,
                       max_pos_prob_x,
                       gate_quant,
                       tol_clust,
                       tol_gate_single,
                       calc_cyt_pos_gates,
                       .debug) {
  # loop over populations
  message("----")
  message("getting base gates")
  message("----")
  message("")
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    txt <- paste0("chnl: ", marker_curr$chnl_cut)
    message(txt)
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate,
      chnl_cut = marker_curr$chnl_cut,
      gate_combn = marker_curr$gate_combn,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      exc_min = marker_curr$exc_min,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x = marker_curr$max_pos_prob_x,
      gate_quant = gate_quant,
      tol_clust = tol_clust,
      tol_gate_single = tol_gate_single,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      .debug = .debug
    )

    path_dir_save <- file.path(path_project, marker_curr$chnl_cut)
    dir.create(path_dir_save, recursive = TRUE, showWarnings = TRUE)

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(path_dir_save, "gate_tbl_init.rds")
    )
  })
}

.gate_single <- function(pop_gate,
                         marker,
                         .data,
                         ind_batch_list,
                         path_project,
                         noise_sd,
                         max_pos_prob_x,
                         gate_quant,
                         tol_clust,
                         tol_gate_single,
                         calc_cyt_pos_gates,
                         calc_single_pos_gates,
                         .debug,
                         gate_tbl) {
  # loop over populations
  message("")
  message("")
  message("----")
  message("getting single+ gates")
  message("----")
  message("")
  if (!calc_single_pos_gates) {
    .debug_msg(.debug, "Not gating single-pos gates") # nolint
    purrr::walk(marker, function(marker_curr) {
      saveRDS(
        gate_tbl |>
          dplyr::filter(chnl == marker_curr$chnl_cut) |>
          dplyr::mutate(gate_single = gate),
        file.path(path_project, marker_curr$chnl_cut, "gate_tbl.rds")
      )
    })
    return(invisible(TRUE))
  } else {
    .debug_msg(.debug, "Gating single-pos gates") # nolint
  }
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    txt <- paste0("chnl: ", marker_curr$chnl_cut)
    message(txt)
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate_curr,
      chnl_cut = marker_curr$chnl_cut,
      gate_combn = marker_curr$gate_combn,
      tol = marker_curr$tol,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      exc_min = marker_curr$exc_min,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x = marker_curr$max_pos_prob_x,
      gate_quant = gate_quant,
      tol_clust = tol_clust,
      tol_gate_single = tol_gate_single,
      gate_tbl = gate_tbl,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      .debug = .debug
    )

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(path_project, marker_curr$chnl_cut, "gate_tbl.rds")
    )
  })
}

.gate_stats <- function(.data,
                        params = NULL,
                        gate_tbl = NULL,
                        filter_other_cyt_pos = FALSE,
                        combn = TRUE,
                        calc_cyt_pos_gates,
                        calc_single_pos_gates,
                        .debug,
                        save = TRUE,
                        pop_gate,
                        marker,
                        ind_batch_list,
                        path_project,
                        tol_clust,
                        save_gate_tbl = TRUE) {
  force(.data)
  .get_stats( # nolint
    params = params,
    gate_tbl = gate_tbl,
    filter_other_cyt_pos = filter_other_cyt_pos,
    combn = combn,
    gate_type_cyt_pos_filter =
      if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_filter =
      if (calc_single_pos_gates) "single" else "base",
    gate_type_cyt_pos_calc =
      if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_calc =
      if (calc_single_pos_gates) "single" else "base",
    .debug = .debug,
    save = save,
    pop_gate = pop_gate,
    chnl = purrr::map_chr(marker, function(x) x$chnl_cut),
    ind_batch_list = ind_batch_list,
    .data = .data,
    save_gate_tbl = save_gate_tbl,
    path_project = path_project,
    tol_clust = tol_clust
  )
}
