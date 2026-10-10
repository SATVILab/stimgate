#' @title Tune stimulation gating
#' @description Set tuning options for [gateStim()]. Start with the defaults;
#'   set only the options you need to change.
#' @param cytPosMethod character Cytokine-positive rule, "refine" or
#'   "coexpression"; used only when `calcCytPosGates = TRUE`. Default: "refine".
#' @param coexNBin integer Number of bins used to lower and trim pairwise gates.
#'   Default: 20L.
#' @param coexResidualMin numeric Minimum Pearson residual beyond independence.
#'   Default: 3.5.
#' @param coexPurityFrac numeric Required fraction of double-positive purity,
#'   greater than zero and at most one. Default: 0.75.
#' @param coexZMin numeric Minimum net double-positive z score. Default: 2.
#' @param calcCytPosGates logical Apply the `cytPosMethod` rule. With "refine",
#'   refine gates using cells positive for another
#'   cytokine. Lower a clustered gate to the leftmost internal antimode between
#'   the stimulated marginal peak plus half its negative-population width
#'   (at least one marginal density bandwidth) and
#'   that gate; keep the gate if none exists. Default: TRUE.
#' @param minCell numeric Minimum cells required for gating; samples below this
#'   count are skipped. Default: 100.
#' @param biasUnsFactor numeric Automatic `biasUns` (when `biasUns` is NULL)
#'   equals this factor times the bandwidth every sample of the marker uses
#'   (a fixed `bw`, or the shared bandwidth with `bwScope = "cytokine"`), and
#'   otherwise times the fallback bandwidth `bwFallback`. Default: 1.
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
#' @param bwScaleNcell logical Widen a shared bandwidth (`bwScope`
#'   "cytokine" or "cluster") for samples whose smaller tube has fewer than
#'   `bwNcellMax` cells, by `(bwNcellMax / n)^(1/5)`, the rate at which a
#'   normal-reference bandwidth grows as cells decrease. The shared bandwidth
#'   is estimated on `bwNcellMax` cells; larger samples keep it. Default:
#'   TRUE.
#' @param bwCluster numeric or NULL Density bandwidth for threshold clustering.
#'   NULL uses the shared local-FDR bandwidth, or with `bwScope = "sample"`,
#'   the median bandwidth of samples with direct thresholds. Default: NULL.
#' @param clusterGates logical Share thresholds across clusters of paired
#'   stimulated/control densities on a common expression grid. Only
#'   responders (see Details) share their thresholds. With at least three
#'   responders in a cluster, clip their thresholds to the cluster's 15th and
#'   85th percentiles. Give other tubes the 60th percentile of the responders'
#'   thresholds; retain high thresholds when a cluster has no responder.
#'   Lowered thresholds are limited by `locShareCap` and `locShareCellCap`.
#'   Default: TRUE.
#' @param gateCombn character vector Combine the thresholds of responders (see
#'   Details) within a batch: "no", "min", "median", "max", or "prejoin".
#'   Lowered thresholds are limited by `locShareCap` and `locShareCellCap`.
#'   Default: "min".
#' @param locProbCol character Probability used for trimming: "pred" (monotone
#'   smoothed response probability) or "probSmooth" (raw/interpolated
#'   probability). Default: "pred".
#' @param locMinPeakProb numeric Minimum peak response probability for a
#'   directly generated local-FDR threshold. Default: 0.25.
#' @param locMinRiseProb numeric Minimum response probability that a rise in
#'   the fitted probability must reach, where it levels off, before the next
#'   rise, for its steepest point to start the response region. Rises that
#'   level off lower are skipped in favour of the next one to the right. 0
#'   turns the check off. Default: 1/3.
#' @param locThresholdMethod character How the condition-level gate is chosen
#'   from the local-FDR filtering: "region" uses the lower boundary of the
#'   region kept by filtering (`xSum` in diagnostics) as the gate; "match"
#'   moves the gate to where the background-subtracted frequency equals the
#'   sum of fitted response probabilities, which was the only method before
#'   stimgate 0.99.26; "cap" keeps the region boundary unless the
#'   background-subtracted frequency above it exceeds the sum of fitted
#'   response probabilities by more than the factor `locThresholdCap`, in
#'   which case the gate moves to the lowest value at or above the boundary
#'   where it no longer does. Cells count as positive when strictly above the
#'   gate. Default: "cap".
#' @param locThresholdCap numeric Largest allowed ratio of the
#'   background-subtracted frequency to the sum of fitted response
#'   probabilities under `locThresholdMethod = "cap"`; at least 1. Default:
#'   1.3.
#' @param locShareCap numeric How far a shared gate may lower a responder's
#'   gate: only while the background-subtracted frequency stays at most this
#'   multiple of the tube's sum of fitted response probabilities. At least 1;
#'   Inf accepts any lower shared gate. See Details. Default: 1.5.
#' @param locShareCellCap numeric How far a shared gate may lower the gate of
#'   a tube that is not a responder: only while its background-subtracted
#'   frequency stays at most this many cells divided by its number of
#'   stimulated cells, and at most the median frequency of the responders
#'   sharing their gates. At least 0; Inf switches this limit off. See
#'   Details. Default: 0.5.
#' @param locEnforceShapeThreshold logical Refit densities and probabilities
#'   above the lower of the first stimulated-density antimode right of the main
#'   negative peak and the adjusted stimulated-density tailgate. Restrict later
#'   marginal filtering to this region. Default: FALSE.
#' @param locShiftedPeakRef logical Handle a stimulated tube whose main peak
#'   has moved right because most of its cells respond. The search for
#'   responding cells normally starts half the larger negative-population width
#'   above the higher main peak, with a minimum of one local-FDR bandwidth.
#'   Width is measured left from each peak to its half-height point or first
#'   clear dip, whichever is closer; if neither exists, to the data minimum.
#'   When `TRUE` and the stimulated main peak lies more than `locShiftedPeakBwMult`
#'   bandwidths to the right of the unstimulated main peak, the search instead
#'   starts above the unstimulated peak, using only the unstimulated tube's
#'   negative width with the same one-density-bandwidth minimum. The bandwidth
#'   for triggering this rule is the shared bandwidth estimated on
#'   `bwNcellMax` cells, without the widening for small tubes from
#'   `bwScaleNcell`; a fixed `bw` when supplied; the sample's own bandwidth
#'   with `bwScope = "sample"`; and, for adaptive bandwidths, the bandwidth
#'   curve at the unstimulated peak. Tubes where the rule applied are flagged
#'   in the `locShiftedPeakRef` column of [getStimGates()]. Default: FALSE.
#' @param locShiftedPeakBwMult numeric Number of bandwidths by which the
#'   stimulated main peak must exceed the unstimulated main peak for
#'   `locShiftedPeakRef` to apply. Must be positive. Default: 2.
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
#'
#' Thresholds are shared within a batch (`gateCombn`) and then within
#' clusters (`clusterGates`). Only responders share their thresholds: tubes
#' whose own threshold was found by local FDR and that have more stimulated
#' than unstimulated cells above it (counting unshifted unstimulated
#' expression, as [getStimStats()] does). When a shared threshold is lower than
#' a responder's threshold, the responder accepts it only down to the lowest
#' value where its background-subtracted frequency is at most `locShareCap`
#' times its sum of fitted response probabilities. Any other tube accepts a
#' shared threshold only down to the lowest value where its frequency is at
#' most `locShareCellCap` cells divided by its number of stimulated cells, and
#' at most the median frequency of the responders sharing their thresholds.
#' Higher shared thresholds are accepted unchanged by responders. Setting both
#' limits to Inf shares responders' thresholds without limits.
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
  bwScaleNcell = TRUE,
  clusterGates = TRUE,
  gateCombn = "min",
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locMinRiseProb = 1 / 3,
  locThresholdMethod = "cap",
  locThresholdCap = 1.3,
  locShareCap = 1.5,
  locShareCellCap = 0.5,
  locEnforceShapeThreshold = FALSE,
  locShiftedPeakRef = FALSE,
  locShiftedPeakBwMult = 2,
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
  normMtd = "moments",
  cytPosMethod = "refine",
  coexNBin = 20L,
  coexResidualMin = 3.5,
  coexPurityFrac = 0.75,
  coexZMin = 2
) {
  ctrl <- as.list(environment())

  if (!is.character(cytPosMethod) || length(cytPosMethod) != 1L ||
      is.na(cytPosMethod) || !cytPosMethod %in% c("refine", "coexpression")) {
    stop("`cytPosMethod` must be 'refine' or 'coexpression'.")
  }
  for (nm in c("coexNBin", "coexResidualMin", "coexPurityFrac", "coexZMin")) {
    value <- ctrl[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || (nm == "coexNBin" && (value != floor(value) || value > .Machine$integer.max)) ||
        (nm == "coexPurityFrac" && value > 1)) {
      stop("`", nm, "` must be a finite positive number",
        if (nm == "coexNBin") " of whole bins." else if (nm == "coexPurityFrac") " at most 1." else ".")
    }
  }

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
    !is.logical(locShiftedPeakRef) || length(locShiftedPeakRef) != 1L ||
      is.na(locShiftedPeakRef)
  ) {
    stop("`locShiftedPeakRef` must be TRUE or FALSE")
  }

  if (
    !is.logical(clusterGates) || length(clusterGates) != 1L ||
      is.na(clusterGates)
  ) {
    stop("`clusterGates` must be TRUE or FALSE")
  }
  if (
    !is.logical(bwScaleNcell) || length(bwScaleNcell) != 1L ||
      is.na(bwScaleNcell)
  ) {
    stop("`bwScaleNcell` must be TRUE or FALSE")
  }

  # Required global settings must be supplied explicitly; per-channel NULL/NA
  # overrides inherit these values.
  required <- c(
    "excMin", "biasUnsFactor", "bwAdj", "gateCombn", "bwFallback",
    "bwMtd", "bwScope", "locThresholdMethod"
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
