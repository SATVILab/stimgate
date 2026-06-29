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
#' @param pathProject character. Path to project directory where all results will be saved.
#'   This directory will contain subdirectories for each marker with gate tables,
#'   statistics, and plots. The directory will be created if it doesn't exist.
#' @param .data GatingSet. A flowWorkspace GatingSet object containing the flow cytometry
#'   data with both stimulated and unstimulated samples. The GatingSet should have
#'   consistent channel names across all samples and include proper sample annotations.
#' @param batchList list. List where each element contains indices of samples belonging to the same batch/donor. Last index per element is the unstimulated vector, e.g. if `batchList = list(1:3, 4:6)`, then indices 3 and 6 correspond to the unstimulated samples for batches 1 and 2, respectively. If `batchList` is named, e.g. `list(pid1 = 1:3, pid2 = 4:6)`, then these names will be used for batch identification.
#' @param chnl list. List where each element specifies parameters for gating a
#'   specific channel Each element should be a list containing at minimum the channel
#'   name (e.g., list(cut = "IL2")). Additional channel-specific parameters can
#'   override global defaults. The marker name should match channel names in the GatingSet.
#' @param marker character vector. Alternative way to specify markers to gate on.
#'   When provided, this is used instead of chnl to determine which markers to analyze.
#'   Default is NULL.
#' @param popGate character vector. Population(s) within which to perform gating.
#'   Default is "root" to gate on all cells. Can specify other populations like
#'   "CD3+" or "CD4+" if these gates already exist in the GatingSet.
#' @param biasUns numeric. Bias adjustment for unstimulated samples to account for
#'   background cytokine production. When NULL (default), no bias correction is applied.
#'   Positive values shift the unstimulated distribution higher, making gates more
#'   conservative.
#' @param biasUnsFactor numeric. Multiplicative factor applied to biasUns.
#'   Default is 1. Values > 1 increase the bias effect, values < 1 decrease it.
#'   This provides fine-tuning of the bias correction.
#' @param excMin logical. Whether to exclude minimum expression values during
#'   analysis. Default is TRUE. Minimum values often represent technical artifacts
#'   or compensation spillover and should typically be excluded.
#' @param cpMin numeric. Minimum allowable cutpoint value. When NULL (default),
#'   no minimum is enforced. Useful for ensuring gates don't fall below known
#'   technical thresholds or background levels.
#' @param bw numeric. Specify the bandwith for density estimation. When NULL (default), bandwidth is estimated automatically. Default is `NULL`.
#' @param bwMin numeric or character. Minimum bandwidth for density estimation.
#'   Ignored if `bw` is set. Use `"auto"` to calculate automatically, `"none"`
#'   to apply no lower bound, or a numeric value to specify the lower bound.
#'   Default is `"auto"`.
#' @param bwMax numeric or character. Maximum bandwidth for density estimation.
#'   Ignored if `bw` is set. Use `"auto"` to calculate automatically, `"none"`
#'   to apply no upper bound, or a numeric value to specify the upper bound.
#'   Default is `"auto"`.
#' @param bwMtd character. Method for automated bandwidth selection. Options include "nrd0", "sj", "hpi0", "hpi1", "hpi2" and "hpi3", which corresponds to the Silverman rule of thumbg (`"nrd0"`), the Sheather-Jones plug-in estimator (`"sj"`) and the Wand & Jones plugin-estimator for the 0-th, 1st, 2nd and 3rd derivatives of the density (`"hpi1"`, `"hpi1"`, `"hpi2"` and `"hpi3"`). Default is "nrd0". Ignored if `bw` is set. Default is `"hpi1"`.
#' @param bwAdj numeric. Adjustment factor for bandwidth. Default is 1. Ignored if `bw` is set. Default is 1.
#' @param bwNcellMin numeric. Minimum number of cells required for bandwidth estimation. If a sample has fewer cells than `bwNcellMin`, cells are sampled with replacement to reach the minimum, with noise subsequently added. Ignored if `bw` is set. Default is 100.
#' @param bwNcellMax numeric. Maximum number of cells used for bandwidth estimation. If a sample has more cells than `bwNcellMax`, cells are sampled without replacement to reach the maximum. Ignored if `bw` is set. Default is 100 000.
#' @param bwCluster numeric. Optional fallback bandwidth for cluster-based local-FDR imputation. The current cluster step first tries to use a common bandwidth calculated as the median bandwidth across samples with generated local-FDR thresholds. `bwCluster` is used as a fallback when that common bandwidth cannot be estimated. Default is `NULL`.
#' @param minCell numeric. Minimum number of cells required for reliable gating.
#'   Default is 100. Samples with fewer cells will be skipped as they don't provide
#'   sufficient statistical power for accurate gate identification.
#' @param maxPosProbX numeric. Maximum x-value (expression level) to consider
#'   when calculating positive probabilities. Default is Inf (no limit). Can be used
#'   to exclude extremely high expression values that may represent doublets or artifacts.
#' @param gateQuant numeric vector. Quantiles used for gate combination when multiple
#'   gates are identified. Default is c(0.25, 0.75). The method specified in gateCombn
#'   determines how these quantiles are used (e.g., minimum of 25th percentiles).
#' @param tolClust numeric or NULL. Backwards-compatible switch for calculating
#'   cluster-adjusted local-FDR gates. When `NULL`, cluster adjustment is skipped.
#'   When non-NULL, samples without a generated local-FDR threshold may receive an
#'   imputed threshold from similar conditions using cluster-level equivalent
#'   tolerance values. The numeric value is retained for compatibility and is not
#'   used as the derivative tolerance in the current local-FDR cluster imputation.
#' @param locProbCol character. Probability column used by the local-FDR trimming
#'   step. Defaults to `"pred"`, the monotone smoothed response-probability
#'   estimate. Use `"probSmooth"` to force the raw/interpolated probability.
#' @param locMinPeakProb numeric. Minimum peak estimated response probability
#'   required before a local-FDR gate is considered credible. If the maximum
#'   probability is below this value, no true local-FDR threshold is marked as
#'   generated.
#' @param locDipAlpha numeric. Liberal dip-test p-value cutoff used to decide
#'   whether to inspect expression-density antimodes before thresholding.
#' @param locAntimodeHeightFrac numeric. Maximum allowed antimode height as a
#'   fraction of the highest density peak when identifying deep antimodes.
#' @param locAntimodeLowRel numeric. Antimode-separated regions to the left of the
#'   highest-response region are excluded when their mean response probability is
#'   below this fraction of the highest region's mean probability.
#' @param locAntimodeLowAbs numeric. Absolute response-probability cutoff for
#'   excluding low-response antimode-separated regions to the left of the
#'   highest-response region.
#' @param locFlatDerivFrac numeric. Fraction of the maximum positive derivative
#'   used to define the marginal-trimming anchor. The marginal scan then decides
#'   how far left of this anchor to extend, one bin at a time. Default is 0.5.
#' @param locFlatHardDerivFrac numeric. Lower derivative fraction used for a
#'   conservative hard exclusion of the very flat far-left region before the
#'   marginal bin scan. Default is 0.25.
#' @param locLeftLowRel numeric. Retained for backwards compatibility. The current
#'   derivative/marginal left-tail trim no longer uses this overall-region check.
#'   Candidate left-tail regions not separated by an
#'   antimode are considered low response when their mean response probability is
#'   below this fraction of the peak response probability.
#' @param locLeftLowAbs numeric. Retained for backwards compatibility. The
#'   current derivative/marginal left-tail trim no longer uses this overall-region
#'   check. Absolute response-probability cutoff for the non-antimode left-tail
#'   trimming rule.
#' @param locLeftCellFrac numeric. Retained for backwards compatibility. The
#'   current derivative/marginal left-tail trim no longer uses this overall-region
#'   check. Minimum size of the candidate low-response
#'   left-tail region, expressed as a fraction of the number of cells to the right
#'   of the start of the main probability rise.
#' @param locLeftLengthFrac numeric. Retained for backwards compatibility. The
#'   current derivative/marginal left-tail trim no longer uses this overall-region
#'   check. Minimum length of the candidate low-response left-tail region,
#'   expressed as a fraction of the expression interval over which the response
#'   probability rises from its minimum to its maximum.
#' @param locMarginalPurityRel numeric. Minimum allowed purity of each additional
#'   leftward bin, expressed as a fraction of the average response probability
#'   among cells to the right of the initial derivative-based local-FDR boundary.
#'   Default is 0.5.
#' @param locMarginalCellBinRatio numeric. Maximum number of cells allowed in each
#'   additional leftward bin, expressed as a multiple of the average number of
#'   cells per bin in the right-side reference interval. Empty reference bins are
#'   counted. Default is 2.
#' @param locMarginalRefQuantile numeric. Upper quantile of cells to the right of
#'   the initial derivative-based boundary used to define the reference interval
#'   for cells-per-bin calculations. Purity is still calculated using all cells to
#'   the right of the initial boundary. Default is 0.75.
#' @param locTolRefPeak character. Reference peak used by cluster imputation when
#'   calculating signed equivalent tolerance values from generated local-FDR
#'   thresholds. Default is `"highest"`.
#' @param gateCombn character vector. Method(s) for combining condition-level
#'   local-FDR gates within a batch. Supported values are `"no"`, `"min"`,
#'   `"median"`, `"max"`, and `"prejoin"`. Combination uses only thresholds that
#'   were actually generated by the local-FDR procedure, not fallback above-range
#'   cutpoints.
#' @param markerSettings list. Optional list of additional marker-specific settings
#'   that override global defaults. Each element should be named with the marker name
#'   and contain parameter overrides. Default is NULL.
#' @param chnlSettings list. Optional list of channel-specific settings that override
#'   global defaults. Similar to markerSettings but keyed by channel names.
#'   Default is NULL.
#' @param calcCytPosGates logical. Whether to calculate refined cytokine-positive
#'   gates using more sophisticated algorithms. Default is TRUE. When FALSE, only
#'   basic gates are calculated, which may be less accurate but faster.
#' @param calcSinglePosGates logical. Whether to calculate single-positive gates
#'   for individual markers in addition to combination gates. Default is FALSE.
#'   Useful for detailed analysis of individual marker responses.
#' @return character. Returns the path to the project directory where all results
#'   have been saved. The directory structure created includes:
#'   \itemize{
#'     \item \code{pathProject/[markerName]/}: Directory for each marker containing:
#'     \item \code{gateTblInit.rds}: Initial gate table with preliminary gates
#'     \item \code{gateTbl.rds}: Final refined gate table
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
#' \strong{Step 3: Cytokine-Positive Gate Refinement (if calcCytPosGates = TRUE)}
#' \itemize{
#'   \item Applies more sophisticated algorithms to refine initial gates
#'   \item Accounts for background cytokine production and technical variability
#'   \item Optimizes gates to minimize false positives while maintaining sensitivity
#' }
#'
#' \strong{Step 4: Single-Positive Gates (if calcSinglePosGates = TRUE)}
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
#' \code{\link{getStimGates}} for extracting gate information,
#' \code{\link{getStimStats}} for generating statistics from results,
#' \code{\link{plotStim}} for visualizing identified gates,
#' \code{\link{writeStimFCS}} for exporting cytokine-positive cells,
#' \code{\link[flowWorkspace]{GatingSet}} for GatingSet documentation
#'
#' @examples{
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$path_gs)
#' pathProject <- file.path(tempdir(), "demonstration")
#'
#' # Run gating
#' gateStim(
#'   .data = gs,
#'   pathProject = pathProject,
#'   popGate = "root",
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker
#' )
#'
#' # Create plots
#' plots <- plotStim(
#'   ind = exampleData$batchList[[1]], # indices in `gs` to plot
#'   .data = gs, # GatingSet
#'   pathProject = pathProject,
#'   marker = exampleData$marker,
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
gateStim <- function(
  pathProject,
  .data,
  popGate = "root",
  batchList,
  chnl = NULL,
  marker = NULL,
  calcCytPosGates = TRUE,
  calcSinglePosGates = FALSE,
  biasUns = NULL,
  biasUnsFactor = 1,
  excMin = TRUE,
  cpMin = NULL,
  bw = NULL,
  bwMin = NULL,
  bwMax = NULL,
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 1e2,
  bwNcellMax = 1e5,
  bwCluster = NULL,
  minCell = 1e2,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  tolClust = 1e-7,
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locLeftLowRel = 0.25,
  locLeftLowAbs = 0.15,
  locLeftCellFrac = 0.5,
  locLeftLengthFrac = 0.5,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  locTolRefPeak = "highest",
  gateCombn = "min",
  markerSettings = NULL,
  chnlSettings = NULL
) {
  force(.data)
  if (Sys.getenv("STIMGATE_DEBUG") == "") {
    Sys.setenv("STIMGATE_DEBUG" = "FALSE")
    on.exit(Sys.unsetenv("STIMGATE_DEBUG"), add = TRUE)
  }
  if (Sys.getenv("STIMGATE_DEBUG") == "TRUE") {
    .debugFileCreate()
    pathDebug <- .debugFileCreate()
    message(paste0("Saving debug output to ", pathDebug))
    message(
      "Can copy it after the run to working directory with stimgateDebugCopy()"
    ) # nolint
    message(
      "Can print the output after the run to console with stimgate_debug_print()"
    ) # nolint
  }

  # Verify global function inputs
  .verifyGateInputs(
    pathProject = pathProject,
    .data = .data,
    batchList = batchList,
    popGate = popGate,
    chnl = chnl,
    marker = marker,
    chnlSettings = chnlSettings,
    markerSettings = markerSettings,
    calcCytPosGates = calcCytPosGates,
    calcSinglePosGates = calcSinglePosGates,

    # NEWLY ADDED TO PASS TO VERIFIER:
    biasUns = biasUns,
    biasUnsFactor = biasUnsFactor,
    excMin = excMin,
    cpMin = cpMin,
    bw = bw,
    bwMin = "auto",
    bwMax = "auto",
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    bwCluster = bwCluster,
    minCell = minCell,
    maxPosProbX = maxPosProbX,
    gateQuant = gateQuant,
    tolClust = tolClust,
    locProbCol = locProbCol,
    locMinPeakProb = locMinPeakProb,
    locDipAlpha = locDipAlpha,
    locAntimodeHeightFrac = locAntimodeHeightFrac,
    locAntimodeLowRel = locAntimodeLowRel,
    locAntimodeLowAbs = locAntimodeLowAbs,
    locFlatDerivFrac = locFlatDerivFrac,
    locFlatHardDerivFrac = locFlatHardDerivFrac,
    locLeftLowRel = locLeftLowRel,
    locLeftLowAbs = locLeftLowAbs,
    locLeftCellFrac = locLeftCellFrac,
    locLeftLengthFrac = locLeftLengthFrac,
    locMarginalPurityRel = locMarginalPurityRel,
    locMarginalCellBinRatio = locMarginalCellBinRatio,
    locMarginalRefQuantile = locMarginalRefQuantile,
    locTolRefPeak = locTolRefPeak,
    gateCombn = gateCombn
  )

  if (is.null(names(batchList))) {
    batchList <- batchList |>
      stats::setNames(paste0("batch", seq_along(batchList)))
  }

  # get unspecified levels in marker elements
  .saveMetaData(.data, batchList, pathProject)
  chnl <- .extractChnl(chnl, marker, pathProject)

  chnlSettings <- .completeChnlSettings(
    chnl = chnl,
    marker = marker,
    chnlSettings = chnlSettings,
    markerSettings = markerSettings,
    biasUns = biasUns,
    biasUnsFactor = biasUnsFactor,
    excMin = excMin,
    .data = .data,
    popGate = popGate,
    indBatchList = batchList,
    bwMin = bwMin,
    bwMax = bwMax,
    bw = bw,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    bwCluster = bwCluster,
    cpMin = cpMin,
    minCell = minCell,
    tolClust = tolClust,
    locProbCol = locProbCol,
    locMinPeakProb = locMinPeakProb,
    locDipAlpha = locDipAlpha,
    locAntimodeHeightFrac = locAntimodeHeightFrac,
    locAntimodeLowRel = locAntimodeLowRel,
    locAntimodeLowAbs = locAntimodeLowAbs,
    locFlatDerivFrac = locFlatDerivFrac,
    locFlatHardDerivFrac = locFlatHardDerivFrac,
    locLeftLowRel = locLeftLowRel,
    locLeftLowAbs = locLeftLowAbs,
    locLeftCellFrac = locLeftCellFrac,
    locLeftLengthFrac = locLeftLengthFrac,
    locMarginalPurityRel = locMarginalPurityRel,
    locMarginalCellBinRatio = locMarginalCellBinRatio,
    locMarginalRefQuantile = locMarginalRefQuantile,
    locTolRefPeak = locTolRefPeak,
    maxPosProbX = maxPosProbX,
    gateCombn = gateCombn,
    gateQuant = gateQuant,
    pathProject = pathProject
  )

  # inital gates
  .gateInit(
    chnlSettings = chnlSettings,
    .data = .data,
    indBatchList = batchList,
    pathProject = pathProject
  )

  # cytokine-positive gates
  gateTbl <- .gateCytPos(
    chnlSettings = chnlSettings,
    indBatchList = batchList,
    .data = .data,
    calcCytPos = calcCytPosGates,
    stage = "cytPos",
    pathProject = pathProject
  )

  # single-positive gates
  .gateSingle(
    chnlSettings = chnlSettings,
    .data = .data,
    indBatchList = batchList,
    pathProject = pathProject,
    calcCytPosGates = calcCytPosGates,
    calcSinglePosGates = calcSinglePosGates,
    gateTbl = gateTbl
  )

  message("")
  message("")
  message("")
  message("getting cyt combn frequencies")

  .gateStats(
    .data = .data,
    gateTbl = gateTbl,
    calcCytPosGates = calcCytPosGates,
    calcSinglePosGates = calcSinglePosGates,
    chnlSettings = chnlSettings,
    indBatchList = batchList,
    pathProject = pathProject
  )

  pathProject
}

#' @keywords internal
.gateInit <- function(
  chnlSettings,
  .data,
  indBatchList,
  pathProject
) {
  message("----")
  message("getting base gates")
  message("----")
  message("")

  purrr::walk(chnlSettings, function(chnlSettingsCurr) {
    txt <- paste0("chnl: ", chnlSettingsCurr$chnlCut)
    message(txt)

    gateObj <- .gateChnl(
      .data = .data,
      indBatchList = indBatchList,
      chnlSettings = chnlSettingsCurr,
      pathProject = pathProject,
      stage = "init",
      calcCytPosGates = FALSE
    )

    .gateInitSave(
      pathProject = pathProject,
      chnlSettings = chnlSettingsCurr,
      gateTbl = gateObj$gateTbl
    )
  })
}

#' @keywords internal
.gateInitSave <- function(pathProject, chnlSettings, gateTbl) {
  pathSave <- .gatesGetPathAll(
    pathProject = pathProject,
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    init = TRUE
  )
  if (!dir.exists(dirname(pathSave))) {
    dir.create(dirname(pathSave), recursive = TRUE, showWarnings = TRUE)
  }
  saveRDS(gateTbl, pathSave)
}

#' @keywords internal
.gateSingle <- function(
  chnlSettings,
  .data,
  indBatchList,
  pathProject,
  calcCytPosGates,
  calcSinglePosGates,
  gateTbl
) {
  message("")
  message("")
  message("----")
  message("getting single+ gates")
  message("----")
  message("")

  if (!calcSinglePosGates) {
    .debug("Not gating single-pos gates")
    purrr::walk(chnlSettings, function(chnlSettingsCurr) {
      pathSave <- .gatesGetPathAll(
        pathProject = pathProject,
        pop = chnlSettings$popGate,
        chnlCut = chnlSettingsCurr$chnlCut,
        init = FALSE
      )
      if (!dir.exists(dirname(pathSave))) {
        dir.create(dirname(pathSave), recursive = TRUE, showWarnings = TRUE)
      }
      saveRDS(
        gateTbl |>
          dplyr::filter(chnl == chnlSettingsCurr$chnlCut) |>
          dplyr::mutate(gateSingle = gate),
        pathSave
      )
    })
    return(invisible(TRUE))
  } else {
    .debug("Gating single-pos gates")
  }

  purrr::walk(chnlSettings, function(chnlSettingsCurr) {
    txt <- paste0("chnl: ", chnlSettingsCurr$chnlCut)
    message(txt)

    gateObj <- .gateChnl(
      .data = .data,
      indBatchList = indBatchList,
      chnlSettings = chnlSettingsCurr,
      gateTbl = gateTbl,
      calcCytPosGates = calcCytPosGates,
      pathProject = pathProject,
      stage = "single"
    )

    .gateSingleSave(
      pathProject = pathProject,
      chnlSettings = chnlSettings,
      gateTbl = gateObj$gateTbl
    )
  })
}

#' @keywords internal
.gateSingleSave <- function(pathProject, chnlSettings, gateTbl) {
  pathSave <- .gatesGetPathAll(
    pathProject = pathProject,
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    init = FALSE
  )
  if (!dir.exists(dirname(pathSave))) {
    dir.create(dirname(pathSave), recursive = TRUE, showWarnings = TRUE)
  }
  saveRDS(gateTbl, pathSave)
}

#' @keywords internal
.gateStats <- function(
  .data,
  gateTbl = NULL,
  calcCytPosGates,
  calcSinglePosGates,
  chnlSettings,
  indBatchList,
  pathProject
) {
  force(.data)
  .getStats(
    gateTbl = gateTbl,
    filterOtherCytPos = FALSE,
    combn = TRUE,
    gateTypeCytPosFilter = if (calcCytPosGates) "cyt" else "base",
    gateTypeSinglePosFilter = if (calcSinglePosGates) {
      "single"
    } else {
      "base"
    },
    gateTypeCytPosCalc = if (calcCytPosGates) "cyt" else "base",
    gateTypeSinglePosCalc = if (calcSinglePosGates) "single" else "base",
    save = TRUE,
    popGate = chnlSettings[[1]]$popGate,
    chnl = purrr::map_chr(chnlSettings, function(x) x$chnlCut),
    indBatchList = indBatchList,
    .data = .data,
    saveGateTbl = TRUE,
    pathProject = pathProject
  )
}
