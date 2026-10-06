#' @title Tune stimulation gating
#' @description Set tuning options for [gateStim()]. Start with the defaults;
#'   set only the options you need to change.
#' @param calcCytPosGates logical Refine gates using cells positive for another
#'   cytokine. Lower a clustered gate to the leftmost internal antimode between
#'   the stimulated marginal peak plus one third of its left-window width and
#'   that gate; keep the gate if none exists. Default: TRUE.
#' @param minCell numeric Minimum cells required for gating; samples below this
#'   count are skipped. Default: 100.
#' @param biasUnsFactor numeric Automatic `biasUns` (when `biasUns` is NULL)
#'   equals this factor times the fallback bandwidth `bwFallback`. Default: 1.
#' @param excMin logical Exclude cells with minimum expression during gating.
#'   Default: TRUE.
#' @param cpMin numeric or NULL Minimum cutpoint. NULL estimates a 10% trimmed
#'   mean of tube medians after excluding minimum expression. Default: NULL.
#' @param bwMtd character Bandwidth selector: "nrd0", "sj", "hpi0", "hpi1",
#'   "hpi2", "hpi3", or any of these with a "Norm" suffix for background
#'   normalisation. Ignored with fixed `bw`. Default: "nrd0".
#' @param bwScope character Scalar bandwidth sharing: "cytokine" uses a 10%
#'   trimmed mean over about 100 tubes spread across batches; "cluster" groups
#'   tube densities up to the right shoulder of the background modal complex
#'   and uses cluster medians from about 100 tubes (at least one per cluster);
#'   "sample" estimates each stimulated/control pair separately. Smaller
#'   datasets use all tubes. Shared estimates exclude tubes below `minCell`.
#'   Each pair uses the smaller tube bandwidth. Inspect `bwShared` and, for
#'   clusters, `bwSharedTbl` with [stimgateMetaReadSettingsChnls()]. Ignored
#'   with fixed `bw` or adaptive bandwidths. Default: "cytokine".
#' @param bwAdj numeric Bandwidth multiplier; ignored with fixed `bw`. Default: 1.
#' @param bwMin numeric or character Lower bandwidth limit: "auto"
#'   estimates it, "none" disables it, or supply a number. Ignored with
#'   fixed `bw`; estimation failures use `bwFallback`. Default: "auto".
#' @param bwMax numeric or character Upper bandwidth limit: "auto"
#'   estimates it, "none" disables it, or supply a number. Ignored with
#'   fixed `bw`; estimation failures use `bwFallback`. Default: "auto".
#' @param bwFallback numeric or character Bandwidth when automatic estimation
#'   fails: "auto" estimates it from spread-out samples using `bwMtd`,
#'   or supply one finite positive number. NULL and "none" are invalid.
#'   Default: "auto".
#' @param bwNcellMin numeric Minimum selector sample size: a tube with fewer
#'   cells is upsampled (with jitter) to this many before its bandwidth is
#'   chosen. Scalar "Norm" methods use it when selecting background-core and
#'   right-excess cells. Ignored with fixed `bw`. Must not exceed `bwNcellMax`.
#'   Default: `bwNcellMax`, so every tube is resampled to exactly
#'   `bwNcellMax` cells.
#' @param bwNcellMax numeric Maximum selector sample size. Ordinary methods
#'   downsample; scalar "Norm" methods cap the constructed sample after
#'   inspecting the full distribution. Adaptive methods use
#'   `normAdaptiveNcell`. Shared bandwidths (`bwScope` "cytokine" or
#'   "cluster") prefer tubes with at least this many cells, then larger tubes
#'   down to half of it, and estimate each selected tube's bandwidth on exactly
#'   this many cells, upsampling smaller tubes. Ignored with fixed `bw`.
#'   Default: 10000.
#' @param bwCluster numeric or NULL Density bandwidth for threshold clustering.
#'   NULL uses the shared local-FDR bandwidth, or with `bwScope = "sample"`,
#'   the median bandwidth of samples with direct thresholds. Default: NULL.
#' @param clusterGates logical Share thresholds across clusters of paired
#'   stimulated/control densities on a common expression grid. With at least
#'   three direct thresholds, clip them to the cluster's 15th and 85th
#'   percentiles. Replace non-direct thresholds by the 60th percentile of
#'   available direct thresholds; retain high thresholds when none are
#'   available. Default: TRUE.
#' @param gateCombn character vector Combine direct condition-level thresholds
#'   within a batch: "no", "min", "median", "max", or "prejoin". Excludes
#'   fallback above-range cutpoints. Default: "min".
#' @param locProbCol character Probability used for trimming: "pred" (monotone
#'   smoothed response probability) or "probSmooth" (raw/interpolated
#'   probability). Default: "pred".
#' @param locMinPeakProb numeric Minimum peak response probability for a
#'   directly generated local-FDR threshold. Default: 0.25.
#' @param locEnforceShapeThreshold logical Refit densities and probabilities
#'   above the lower of the first stimulated-density antimode right of the main
#'   negative peak and the adjusted stimulated-density tailgate. Restrict later
#'   marginal filtering to this region. Default: FALSE.
#' @param locDipAlpha numeric Dip-test p-value cutoff for inspecting density
#'   antimodes before thresholding. Default: 0.2.
#' @param locAntimodeHeightFrac numeric Maximum antimode height as a fraction
#'   of the highest density peak. Default: 1/6.
#' @param locAntimodeLowRel numeric Exclude antimode-separated regions left of
#'   the highest-response region when their mean response probability is below
#'   this fraction of its mean. Default: 0.25.
#' @param locAntimodeLowAbs numeric Absolute mean response-probability cutoff
#'   for excluding those left-hand regions. Default: 0.15.
#' @param locFlatDerivFrac numeric Fraction of the maximum positive derivative
#'   defining the marginal-trimming anchor; scan leftward from there one bin at
#'   a time. Default: 0.5.
#' @param locFlatHardDerivFrac numeric Derivative fraction for hard exclusion
#'   of the flat far-left region before the marginal scan. Default: 0.25.
#' @param locMarginalPurityRel numeric Minimum purity of each added leftward
#'   bin, relative to average response probability right of the initial
#'   derivative boundary. Default: 0.5.
#' @param locMarginalCellBinRatio numeric Maximum cells per added leftward bin,
#'   as a multiple of average cells per reference bin, including empty bins.
#'   Default: 2.
#' @param locMarginalRefQuantile numeric Upper cell quantile right of the
#'   initial derivative boundary defining the cells-per-bin reference interval.
#'   Purity uses all cells right of that boundary. Default: 0.75.
#' @param bwAdaptive logical Use location-specific bandwidths when `bw` is
#'   NULL. Blend stimulated and control normalised bandwidth curves by
#'   preliminary density heights; use the shared curve for both final
#'   densities. Default: FALSE.
#' @param bwAdaptiveDensityN numeric or NULL Adaptive density grid points; NULL
#'   uses `normDensityN`. Default: NULL.
#' @param bwAdaptivePadFrac numeric Extend each end of the adaptive grid by
#'   this fraction of the combined expression range before area normalisation.
#'   Default: 0.15.
#' @param bwAdaptiveCore numeric or NULL Override the estimated
#'   background-component bandwidth in the adaptive curve. Default: NULL.
#' @param bwAdaptiveExtra numeric or NULL Override the estimated
#'   high-expression-component bandwidth in the adaptive curve. Default: NULL.
#' @param bwAdaptiveCrossover numeric or NULL Expression value at which the
#'   adaptive curve switches between components; NULL uses component-density
#'   weighting. Default: NULL.
#' @param bwAdaptiveTransitionWidth numeric Transition width in expression
#'   units around `bwAdaptiveCrossover`; 0 makes a hard switch. Default: 0.
#' @param normPeakMinRel numeric Relative peak/trough threshold identifying the
#'   background modal complex for "Norm" methods. Default: 0.75.
#' @param normExtraFrac numeric Target fraction of additional high-side values
#'   sampled for the normalised high component. Default: 0.2.
#' @param normExtraMax numeric Maximum additional high-side values for
#'   normalised methods; Inf allows no cap. Default: Inf.
#' @param normLambda numeric vector Box-Cox search grid; used only with
#'   `normMtd = "boxcox"`. Default: seq(-2, 2, length.out = 81).
#' @param normDensityN numeric Grid points for normalised bandwidth estimation.
#'   Default: 512.
#' @param normExcessBwMtd character Ordinary bandwidth selector for right-side
#'   excess density: "nrd0", "sj", "hpi0", "hpi1", "hpi2", or "hpi3".
#'   Default: "hpi3".
#' @param normExcessNcell numeric Maximum cells for estimating the
#'   excess-density bandwidth. Default: 10000.
#' @param normAdaptiveNcell numeric Fixed synthetic normal sample size per
#'   component for adaptive normalised bandwidths. Default: 2500.
#' @param normMtd character Normalisation for "Norm" selectors: "moments"
#'   replaces components with normals having matching moments; "boxcox" uses
#'   Box-Cox transformations for scalar bandwidths only. Default: "moments".
#' @details
#' Settings follow the gating workflow: cell and expression limits, bandwidth
#' selection and sharing, threshold sharing, local-FDR filtering, then adaptive
#' and normalised bandwidths. A fixed `bw` on [gateStim()] overrides automatic
#' bandwidth selection. Use its `markerControl` for per-marker settings.
#' @return A named list of class `stimControl`, with one element per setting.
#' @examples
#' stimControl()
#' stimControl(bwAdj = 1.5, clusterGates = FALSE)
#' @export
stimControl <- function(
  calcCytPosGates = TRUE,
  minCell = 1e2,
  biasUnsFactor = 1,
  excMin = TRUE,
  cpMin = NULL,
  bwMtd = "nrd0",
  bwScope = "cytokine",
  bwAdj = 1,
  bwMin = "auto",
  bwMax = "auto",
  bwFallback = "auto",
  bwNcellMin = bwNcellMax,
  bwNcellMax = 1e4,
  bwCluster = NULL,
  clusterGates = TRUE,
  gateCombn = "min",
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locEnforceShapeThreshold = FALSE,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  bwAdaptive = FALSE,
  bwAdaptiveDensityN = NULL,
  bwAdaptivePadFrac = 0.15,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  normMtd = "moments"
) {
  ctrl <- as.list(environment())

  if (
    !is.logical(calcCytPosGates) || length(calcCytPosGates) != 1L ||
      is.na(calcCytPosGates)
  ) {
    stop("`calcCytPosGates` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.numeric(minCell) || length(minCell) != 1L || minCell <= 0) {
    stop("`minCell` must be a positive number.")
  }
  if (
    !is.logical(locEnforceShapeThreshold) ||
      length(locEnforceShapeThreshold) != 1L ||
      is.na(locEnforceShapeThreshold)
  ) {
    stop("`locEnforceShapeThreshold` must be TRUE or FALSE")
  }

  if (
    !is.logical(clusterGates) || length(clusterGates) != 1L ||
      is.na(clusterGates)
  ) {
    stop("`clusterGates` must be TRUE or FALSE")
  }

  # Required global settings must be supplied explicitly; per-channel NULL/NA
  # overrides inherit these values.
  required <- c(
    "excMin", "biasUnsFactor", "bwAdj", "gateCombn", "bwFallback",
    "bwMtd", "bwScope"
  )
  isMissing <- vapply(ctrl[required], .verifyIsNullOrNa, logical(1))
  if (any(isMissing)) {
    stop(
      "Must be supplied (not NULL or NA): ",
      paste0("`", required[isMissing], "`", collapse = ", ")
    )
  }
  .verifyChnlSettingsChnl(settings = ctrl, prefix = "")

  structure(ctrl, class = "stimControl")
}
