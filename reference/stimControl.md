# Tune stimulation gating

Set tuning options for
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).
Start with the defaults; set only the options you need to change.

## Usage

``` r
stimControl(
  calcCytPosGates = TRUE,
  minCell = 100,
  biasUnsFactor = 1,
  excMin = TRUE,
  cpMin = NULL,
  bwMtd = "hpi1",
  bwScope = "cytokine",
  bwAdj = 1,
  bwMin = "auto",
  bwMax = "auto",
  bwFallback = "auto",
  bwNcellMin = 100,
  bwNcellMax = 1e+05,
  bwCluster = NULL,
  clusterGates = TRUE,
  gateCombn = "min",
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locEnforceShapeThreshold = FALSE,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1/6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1/2,
  locFlatHardDerivFrac = 1/4,
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
)
```

## Arguments

- calcCytPosGates:

  logical Refine gates using cells positive for another cytokine. Lower
  a clustered gate to the leftmost internal antimode between the
  stimulated marginal peak plus one third of its left-window width and
  that gate; keep the gate if none exists. Default: TRUE.

- minCell:

  numeric Minimum cells required for gating; samples below this count
  are skipped. Default: 100.

- biasUnsFactor:

  numeric Multiplier for automatically chosen `biasUns`. Default: 1.

- excMin:

  logical Exclude cells with minimum expression during gating. Default:
  TRUE.

- cpMin:

  numeric or NULL Minimum cutpoint. NULL estimates a 10% trimmed mean of
  tube medians after excluding minimum expression. Default: NULL.

- bwMtd:

  character Bandwidth selector: "nrd0", "sj", "hpi0", "hpi1", "hpi2",
  "hpi3", or any of these with a "Norm" suffix for background
  normalisation. Ignored with fixed `bw`. Default: "hpi1".

- bwScope:

  character Scalar bandwidth sharing: "cytokine" uses a 10% trimmed mean
  over about 100 tubes spread across batches; "cluster" groups tube
  densities up to the right shoulder of the background modal complex and
  uses cluster medians from about 100 tubes (at least one per cluster);
  "sample" estimates each stimulated/control pair separately. Smaller
  datasets use all tubes. Shared estimates exclude tubes below
  `minCell`. Each pair uses the smaller tube bandwidth. Inspect
  `bwShared` and, for clusters, `bwSharedTbl` with
  [`stimgateMetaReadSettingsChnls()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnls.md).
  Ignored with fixed `bw` or adaptive bandwidths. Default: "cytokine".

- bwAdj:

  numeric Bandwidth multiplier; ignored with fixed `bw`. Default: 1.

- bwMin:

  numeric or character Lower bandwidth limit: "auto" estimates it,
  "none" disables it, or supply a number. Ignored with fixed `bw`;
  estimation failures use `bwFallback`. Default: "auto".

- bwMax:

  numeric or character Upper bandwidth limit: "auto" estimates it,
  "none" disables it, or supply a number. Ignored with fixed `bw`;
  estimation failures use `bwFallback`. Default: "auto".

- bwFallback:

  numeric or character Bandwidth when automatic estimation fails: "auto"
  estimates it from spread-out samples using `bwMtd`, or supply one
  finite positive number. NULL and "none" are invalid. Default: "auto".

- bwNcellMin:

  numeric Minimum selector sample size. Ordinary methods upsample with
  jitter; scalar "Norm" methods use it when selecting background-core
  and right-excess cells. Ignored with fixed `bw`. Default: 100.

- bwNcellMax:

  numeric Maximum selector sample size. Ordinary methods downsample;
  scalar "Norm" methods cap the constructed sample after inspecting the
  full distribution. Adaptive methods use `normAdaptiveNcell`. Ignored
  with fixed `bw`. Default: 100000.

- bwCluster:

  numeric or NULL Density bandwidth for threshold clustering. NULL uses
  the shared local-FDR bandwidth, or with `bwScope = "sample"`, the
  median bandwidth of samples with direct thresholds. Default: NULL.

- clusterGates:

  logical Share thresholds across clusters of paired stimulated/control
  densities on a common expression grid. With at least three direct
  thresholds, clip them to the cluster's 15th and 85th percentiles.
  Replace non-direct thresholds by the 60th percentile of available
  direct thresholds; retain high thresholds when none are available.
  Default: TRUE.

- gateCombn:

  character vector Combine direct condition-level thresholds within a
  batch: "no", "min", "median", "max", or "prejoin". Excludes fallback
  above-range cutpoints. Default: "min".

- locProbCol:

  character Probability used for trimming: "pred" (monotone smoothed
  response probability) or "probSmooth" (raw/interpolated probability).
  Default: "pred".

- locMinPeakProb:

  numeric Minimum peak response probability for a directly generated
  local-FDR threshold. Default: 0.25.

- locEnforceShapeThreshold:

  logical Refit densities and probabilities above the lower of the first
  stimulated-density antimode right of the main negative peak and the
  adjusted stimulated-density tailgate. Restrict later marginal
  filtering to this region. Default: FALSE.

- locDipAlpha:

  numeric Dip-test p-value cutoff for inspecting density antimodes
  before thresholding. Default: 0.2.

- locAntimodeHeightFrac:

  numeric Maximum antimode height as a fraction of the highest density
  peak. Default: 1/6.

- locAntimodeLowRel:

  numeric Exclude antimode-separated regions left of the
  highest-response region when their mean response probability is below
  this fraction of its mean. Default: 0.25.

- locAntimodeLowAbs:

  numeric Absolute mean response-probability cutoff for excluding those
  left-hand regions. Default: 0.15.

- locFlatDerivFrac:

  numeric Fraction of the maximum positive derivative defining the
  marginal-trimming anchor; scan leftward from there one bin at a time.
  Default: 0.5.

- locFlatHardDerivFrac:

  numeric Derivative fraction for hard exclusion of the flat far-left
  region before the marginal scan. Default: 0.25.

- locMarginalPurityRel:

  numeric Minimum purity of each added leftward bin, relative to average
  response probability right of the initial derivative boundary.
  Default: 0.5.

- locMarginalCellBinRatio:

  numeric Maximum cells per added leftward bin, as a multiple of average
  cells per reference bin, including empty bins. Default: 2.

- locMarginalRefQuantile:

  numeric Upper cell quantile right of the initial derivative boundary
  defining the cells-per-bin reference interval. Purity uses all cells
  right of that boundary. Default: 0.75.

- bwAdaptive:

  logical Use location-specific bandwidths when `bw` is NULL. Blend
  stimulated and control normalised bandwidth curves by preliminary
  density heights; use the shared curve for both final densities.
  Default: FALSE.

- bwAdaptiveDensityN:

  numeric or NULL Adaptive density grid points; NULL uses
  `normDensityN`. Default: NULL.

- bwAdaptivePadFrac:

  numeric Extend each end of the adaptive grid by this fraction of the
  combined expression range before area normalisation. Default: 0.15.

- bwAdaptiveCore:

  numeric or NULL Override the estimated background-component bandwidth
  in the adaptive curve. Default: NULL.

- bwAdaptiveExtra:

  numeric or NULL Override the estimated high-expression-component
  bandwidth in the adaptive curve. Default: NULL.

- bwAdaptiveCrossover:

  numeric or NULL Expression value at which the adaptive curve switches
  between components; NULL uses component-density weighting. Default:
  NULL.

- bwAdaptiveTransitionWidth:

  numeric Transition width in expression units around
  `bwAdaptiveCrossover`; 0 makes a hard switch. Default: 0.

- normPeakMinRel:

  numeric Relative peak/trough threshold identifying the background
  modal complex for "Norm" methods. Default: 0.75.

- normExtraFrac:

  numeric Target fraction of additional high-side values sampled for the
  normalised high component. Default: 0.2.

- normExtraMax:

  numeric Maximum additional high-side values for normalised methods;
  Inf allows no cap. Default: Inf.

- normLambda:

  numeric vector Box-Cox search grid; used only with
  `normMtd = "boxcox"`. Default: seq(-2, 2, length.out = 81).

- normDensityN:

  numeric Grid points for normalised bandwidth estimation. Default: 512.

- normExcessBwMtd:

  character Ordinary bandwidth selector for right-side excess density:
  "nrd0", "sj", "hpi0", "hpi1", "hpi2", or "hpi3". Default: "hpi3".

- normExcessNcell:

  numeric Maximum cells for estimating the excess-density bandwidth.
  Default: 10000.

- normAdaptiveNcell:

  numeric Fixed synthetic normal sample size per component for adaptive
  normalised bandwidths. Default: 2500.

- normMtd:

  character Normalisation for "Norm" selectors: "moments" replaces
  components with normals having matching moments; "boxcox" uses Box-Cox
  transformations for scalar bandwidths only. Default: "moments".

## Value

A named list of class `stimControl`, with one element per setting.

## Details

Settings follow the gating workflow: cell and expression limits,
bandwidth selection and sharing, threshold sharing, local-FDR filtering,
then adaptive and normalised bandwidths. A fixed `bw` on
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
overrides automatic bandwidth selection. Use its `markerControl` for
per-marker settings.

## Examples

``` r
stimControl()
#> $calcCytPosGates
#> [1] TRUE
#> 
#> $minCell
#> [1] 100
#> 
#> $biasUnsFactor
#> [1] 1
#> 
#> $excMin
#> [1] TRUE
#> 
#> $cpMin
#> NULL
#> 
#> $bwMtd
#> [1] "hpi1"
#> 
#> $bwScope
#> [1] "cytokine"
#> 
#> $bwAdj
#> [1] 1
#> 
#> $bwMin
#> [1] "auto"
#> 
#> $bwMax
#> [1] "auto"
#> 
#> $bwFallback
#> [1] "auto"
#> 
#> $bwNcellMin
#> [1] 100
#> 
#> $bwNcellMax
#> [1] 1e+05
#> 
#> $bwCluster
#> NULL
#> 
#> $clusterGates
#> [1] TRUE
#> 
#> $gateCombn
#> [1] "min"
#> 
#> $locProbCol
#> [1] "pred"
#> 
#> $locMinPeakProb
#> [1] 0.25
#> 
#> $locEnforceShapeThreshold
#> [1] FALSE
#> 
#> $locDipAlpha
#> [1] 0.2
#> 
#> $locAntimodeHeightFrac
#> [1] 0.1666667
#> 
#> $locAntimodeLowRel
#> [1] 0.25
#> 
#> $locAntimodeLowAbs
#> [1] 0.15
#> 
#> $locFlatDerivFrac
#> [1] 0.5
#> 
#> $locFlatHardDerivFrac
#> [1] 0.25
#> 
#> $locMarginalPurityRel
#> [1] 0.5
#> 
#> $locMarginalCellBinRatio
#> [1] 2
#> 
#> $locMarginalRefQuantile
#> [1] 0.75
#> 
#> $bwAdaptive
#> [1] FALSE
#> 
#> $bwAdaptiveDensityN
#> NULL
#> 
#> $bwAdaptivePadFrac
#> [1] 0.15
#> 
#> $bwAdaptiveCore
#> NULL
#> 
#> $bwAdaptiveExtra
#> NULL
#> 
#> $bwAdaptiveCrossover
#> NULL
#> 
#> $bwAdaptiveTransitionWidth
#> [1] 0
#> 
#> $normPeakMinRel
#> [1] 0.75
#> 
#> $normExtraFrac
#> [1] 0.2
#> 
#> $normExtraMax
#> [1] Inf
#> 
#> $normLambda
#>  [1] -2.00 -1.95 -1.90 -1.85 -1.80 -1.75 -1.70 -1.65 -1.60 -1.55 -1.50 -1.45
#> [13] -1.40 -1.35 -1.30 -1.25 -1.20 -1.15 -1.10 -1.05 -1.00 -0.95 -0.90 -0.85
#> [25] -0.80 -0.75 -0.70 -0.65 -0.60 -0.55 -0.50 -0.45 -0.40 -0.35 -0.30 -0.25
#> [37] -0.20 -0.15 -0.10 -0.05  0.00  0.05  0.10  0.15  0.20  0.25  0.30  0.35
#> [49]  0.40  0.45  0.50  0.55  0.60  0.65  0.70  0.75  0.80  0.85  0.90  0.95
#> [61]  1.00  1.05  1.10  1.15  1.20  1.25  1.30  1.35  1.40  1.45  1.50  1.55
#> [73]  1.60  1.65  1.70  1.75  1.80  1.85  1.90  1.95  2.00
#> 
#> $normDensityN
#> [1] 512
#> 
#> $normExcessBwMtd
#> [1] "hpi3"
#> 
#> $normExcessNcell
#> [1] 10000
#> 
#> $normAdaptiveNcell
#> [1] 2500
#> 
#> $normMtd
#> [1] "moments"
#> 
#> attr(,"class")
#> [1] "stimControl"
stimControl(bwAdj = 1.5, clusterGates = FALSE)
#> $calcCytPosGates
#> [1] TRUE
#> 
#> $minCell
#> [1] 100
#> 
#> $biasUnsFactor
#> [1] 1
#> 
#> $excMin
#> [1] TRUE
#> 
#> $cpMin
#> NULL
#> 
#> $bwMtd
#> [1] "hpi1"
#> 
#> $bwScope
#> [1] "cytokine"
#> 
#> $bwAdj
#> [1] 1.5
#> 
#> $bwMin
#> [1] "auto"
#> 
#> $bwMax
#> [1] "auto"
#> 
#> $bwFallback
#> [1] "auto"
#> 
#> $bwNcellMin
#> [1] 100
#> 
#> $bwNcellMax
#> [1] 1e+05
#> 
#> $bwCluster
#> NULL
#> 
#> $clusterGates
#> [1] FALSE
#> 
#> $gateCombn
#> [1] "min"
#> 
#> $locProbCol
#> [1] "pred"
#> 
#> $locMinPeakProb
#> [1] 0.25
#> 
#> $locEnforceShapeThreshold
#> [1] FALSE
#> 
#> $locDipAlpha
#> [1] 0.2
#> 
#> $locAntimodeHeightFrac
#> [1] 0.1666667
#> 
#> $locAntimodeLowRel
#> [1] 0.25
#> 
#> $locAntimodeLowAbs
#> [1] 0.15
#> 
#> $locFlatDerivFrac
#> [1] 0.5
#> 
#> $locFlatHardDerivFrac
#> [1] 0.25
#> 
#> $locMarginalPurityRel
#> [1] 0.5
#> 
#> $locMarginalCellBinRatio
#> [1] 2
#> 
#> $locMarginalRefQuantile
#> [1] 0.75
#> 
#> $bwAdaptive
#> [1] FALSE
#> 
#> $bwAdaptiveDensityN
#> NULL
#> 
#> $bwAdaptivePadFrac
#> [1] 0.15
#> 
#> $bwAdaptiveCore
#> NULL
#> 
#> $bwAdaptiveExtra
#> NULL
#> 
#> $bwAdaptiveCrossover
#> NULL
#> 
#> $bwAdaptiveTransitionWidth
#> [1] 0
#> 
#> $normPeakMinRel
#> [1] 0.75
#> 
#> $normExtraFrac
#> [1] 0.2
#> 
#> $normExtraMax
#> [1] Inf
#> 
#> $normLambda
#>  [1] -2.00 -1.95 -1.90 -1.85 -1.80 -1.75 -1.70 -1.65 -1.60 -1.55 -1.50 -1.45
#> [13] -1.40 -1.35 -1.30 -1.25 -1.20 -1.15 -1.10 -1.05 -1.00 -0.95 -0.90 -0.85
#> [25] -0.80 -0.75 -0.70 -0.65 -0.60 -0.55 -0.50 -0.45 -0.40 -0.35 -0.30 -0.25
#> [37] -0.20 -0.15 -0.10 -0.05  0.00  0.05  0.10  0.15  0.20  0.25  0.30  0.35
#> [49]  0.40  0.45  0.50  0.55  0.60  0.65  0.70  0.75  0.80  0.85  0.90  0.95
#> [61]  1.00  1.05  1.10  1.15  1.20  1.25  1.30  1.35  1.40  1.45  1.50  1.55
#> [73]  1.60  1.65  1.70  1.75  1.80  1.85  1.90  1.95  2.00
#> 
#> $normDensityN
#> [1] 512
#> 
#> $normExcessBwMtd
#> [1] "hpi3"
#> 
#> $normExcessNcell
#> [1] 10000
#> 
#> $normAdaptiveNcell
#> [1] 2500
#> 
#> $normMtd
#> [1] "moments"
#> 
#> attr(,"class")
#> [1] "stimControl"
```
