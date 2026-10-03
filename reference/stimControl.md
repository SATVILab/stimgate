# Tuning settings for stimulation gating

Construct validated tuning settings for
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

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

  logical. Whether to refine each clustered one-marker gate using the
  target-marker distribution among cells positive for at least one other
  cytokine. A taut-string density is fitted to those cells. The
  clustered gate is lowered to the leftmost internal antimode strictly
  between the full stimulated marginal peak plus one third of its
  left-window width and the clustered gate. If no eligible antimode
  exists, the clustered gate is retained. Default is TRUE.

- minCell:

  numeric. Minimum number of cells required for reliable gating. Default
  is 100. Samples with fewer cells will be skipped as they don't provide
  sufficient statistical power for accurate gate identification.

- biasUnsFactor:

  numeric. Multiplicative factor applied to biasUns. Default is 1.
  Values \> 1 increase the bias effect, values \< 1 decrease it. This
  provides fine-tuning of the bias correction.

- excMin:

  logical. Whether to exclude minimum expression values during analysis.
  Default is TRUE. Minimum values often represent technical artifacts or
  compensation spillover and should typically be excluded.

- cpMin:

  numeric. Minimum allowable cutpoint value. When NULL (default), no
  minimum is enforced. Useful for ensuring gates don't fall below known
  technical thresholds or background levels.

- bwMtd:

  character. Method for automated bandwidth selection. Options include
  `"nrd0"`, `"sj"`, `"hpi0"`, `"hpi1"`, `"hpi2"` and `"hpi3"`, plus
  background-normalised variants `"nrd0Norm"`, `"sjNorm"`, `"hpi0Norm"`,
  `"hpi1Norm"`, `"hpi2Norm"` and `"hpi3Norm"`. The normalised variants
  first identify a background core and a high-side component, then
  estimate bandwidths after component normalisation. By default this
  uses moment-matched normal components (`normMtd = "moments"`); the
  older Box-Cox route remains available for scalar bandwidths with
  `normMtd = "boxcox"`. Ignored if `bw` is set. Default is `"hpi1"`.

- bwScope:

  "cytokine", "cluster" or "sample". Which samples share the scalar
  local-FDR bandwidth. `"cytokine"` estimates the bandwidths of about
  100 tubes spread across the batches (all tubes when there are fewer)
  and uses their 10% trimmed mean for every sample of the channel.
  `"cluster"` clusters all tubes on their densities up to the right
  shoulder of the left modal complex, estimates bandwidths for about 100
  tubes spread across the clusters (all tubes when there are fewer, and
  at least one per cluster) and gives each tube its cluster's median
  bandwidth. `"sample"` estimates the bandwidth separately for every
  stimulated sample. Tubes with fewer than `minCell` cells are excluded
  from shared bandwidths. In every case a sample uses the smaller of its
  stimulated and unstimulated tube bandwidths. The chosen values are
  reported while gating and saved as `bwShared` (and, for `"cluster"`,
  the per-tube table `bwSharedTbl`) in
  `stimgateMetaReadSettingsChnls(pathProject)`, so an automatic value
  can be inspected and then fixed for a marker through `bw` in
  `markerControl`. Ignored if `bw` is set or the adaptive bandwidth is
  used. Default is `"cytokine"`.

- bwAdj:

  numeric. Adjustment factor for bandwidth. Default is 1. Ignored if
  `bw` is set. Default is 1.

- bwMin:

  numeric or character. Minimum bandwidth for density estimation.
  Ignored if `bw` is set. Use `"auto"` to calculate automatically,
  `"none"` to apply no lower bound, or a numeric value to specify the
  lower bound. This is only a clipping limit; if automatic bandwidth
  estimation fails, `bwFallback` is used instead. Default is `"auto"`.

- bwMax:

  numeric or character. Maximum bandwidth for density estimation.
  Ignored if `bw` is set. Use `"auto"` to calculate automatically,
  `"none"` to apply no upper bound, or a numeric value to specify the
  upper bound. This is only a clipping limit; if automatic bandwidth
  estimation fails, `bwFallback` is used instead. Default is `"auto"`.

- bwFallback:

  numeric or character. Fallback bandwidth used whenever `bw` is `NULL`
  and automatic bandwidth estimation fails for a sample. Must be either
  `"auto"` or a single positive numeric value. Unlike `bwMin` and
  `bwMax`, `bwFallback` cannot be `NULL`, `"none"`, non-finite, zero, or
  negative, because a valid fallback bandwidth is always required. When
  `"auto"`, the fallback is calculated from randomly selected samples
  using the same bandwidth selector specified by `bwMtd`. Default is
  `"auto"`.

- bwNcellMin:

  numeric. Minimum number of cells requested by the bandwidth selector.
  For ordinary methods this controls internal up-sampling with jitter.
  For `*Norm` methods it is passed into the background-core/right-excess
  selector so rare right-tail cells are considered before any sampling
  is done. Ignored if `bw` is set. Default is 100.

- bwNcellMax:

  numeric. Maximum number of cells requested by the bandwidth selector.
  For ordinary methods this controls internal down-sampling. For `*Norm`
  methods it limits the constructed background-core/right-excess
  bandwidth sample after the full distribution has been inspected.
  Ignored if `bw` is set. Default is 100 000.

- bwCluster:

  numeric or NULL. Bandwidth for the densities clustered by the
  cluster-based threshold sharing (`clusterGates`). When `NULL`, the
  shared local-FDR bandwidth (see `bwScope`) is used, or, with
  `bwScope = "sample"`, the median bandwidth across samples with
  directly generated local-FDR thresholds. Default is `NULL`.

- clusterGates:

  logical. Whether to calculate cluster-adjusted local-FDR gates. When
  `FALSE`, cluster adjustment is skipped. When `TRUE`, paired stimulated
  and unstimulated densities are clustered on a common
  absolute-expression grid. Direct thresholds are winsorised within each
  cluster to its 15th and 85th percentiles when at least three direct
  thresholds are available. Every non-direct threshold is replaced by
  the cluster's 60th percentile when at least one direct threshold is
  available. A cluster without a direct threshold retains its original
  high thresholds. Default is `TRUE`.

- gateCombn:

  character vector. Method(s) for combining condition-level local-FDR
  gates within a batch. Supported values are `"no"`, `"min"`,
  `"median"`, `"max"`, and `"prejoin"`. Combination uses only thresholds
  that were actually generated by the local-FDR procedure, not fallback
  above-range cutpoints.

- locProbCol:

  character. Probability column used by the local-FDR trimming step.
  Defaults to `"pred"`, the monotone smoothed response-probability
  estimate. Use `"probSmooth"` to force the raw/interpolated
  probability.

- locMinPeakProb:

  numeric. Minimum peak estimated response probability required before a
  local-FDR gate is considered credible. If the maximum probability is
  below this value, no true local-FDR threshold is marked as generated.

- locEnforceShapeThreshold:

  logical. Whether density shape must first define the lowest expression
  value allowed to inform local-FDR thresholding. When `TRUE`, the lower
  of the first stimulated-density antimode to the right of the main
  negative peak and the adjusted stimulated-density tailgate is applied
  to both samples before densities, response probabilities, and the
  monotone probability curve are refitted. All subsequent marginal
  filtering is restricted to this refitted region. Default is `FALSE`.

- locDipAlpha:

  numeric. Liberal dip-test p-value cutoff used to decide whether to
  inspect expression-density antimodes before thresholding.

- locAntimodeHeightFrac:

  numeric. Maximum allowed antimode height as a fraction of the highest
  density peak when identifying deep antimodes.

- locAntimodeLowRel:

  numeric. Antimode-separated regions to the left of the
  highest-response region are excluded when their mean response
  probability is below this fraction of the highest region's mean
  probability.

- locAntimodeLowAbs:

  numeric. Absolute response-probability cutoff for excluding
  low-response antimode-separated regions to the left of the
  highest-response region.

- locFlatDerivFrac:

  numeric. Fraction of the maximum positive derivative used to define
  the marginal-trimming anchor. The marginal scan then decides how far
  left of this anchor to extend, one bin at a time. Default is 0.5.

- locFlatHardDerivFrac:

  numeric. Lower derivative fraction used for a conservative hard
  exclusion of the very flat far-left region before the marginal bin
  scan. Default is 0.25.

- locMarginalPurityRel:

  numeric. Minimum allowed purity of each additional leftward bin,
  expressed as a fraction of the average response probability among
  cells to the right of the initial derivative-based local-FDR boundary.
  Default is 0.5.

- locMarginalCellBinRatio:

  numeric. Maximum number of cells allowed in each additional leftward
  bin, expressed as a multiple of the average number of cells per bin in
  the right-side reference interval. Empty reference bins are counted.
  Default is 2.

- locMarginalRefQuantile:

  numeric. Upper quantile of cells to the right of the initial
  derivative-based boundary used to define the reference interval for
  cells-per-bin calculations. Purity is still calculated using all cells
  to the right of the initial boundary. Default is 0.75.

- bwAdaptive:

  logical. Whether local-FDR density estimation should use an adaptive
  location-specific bandwidth curve when `bw` is `NULL`. The adaptive
  path estimates separate normalised bandwidth curves for the stimulated
  and unstimulated samples, blends them by their preliminary density
  heights on a shared padded grid, and then evaluates both final
  densities with that shared bandwidth vector. Default is `FALSE`.

- bwAdaptiveDensityN:

  numeric. Number of grid points for the adaptive local-FDR density
  grid. When `NULL`, `normDensityN` is used. Default is `NULL`.

- bwAdaptivePadFrac:

  numeric. Fraction of the combined expression range by which the
  adaptive density grid is extended on both sides before area
  normalisation. Default is `0.15`.

- bwAdaptiveCore:

  numeric or NULL. Optional manually specified bandwidth for the
  background-core side of the adaptive local-FDR density curve. When
  supplied, it overrides the estimated core-component bandwidth for
  adaptive bandwidth construction. Default is `NULL`.

- bwAdaptiveExtra:

  numeric or NULL. Optional manually specified bandwidth for the
  high-expression/extra side of the adaptive local-FDR density curve.
  When supplied, it overrides the estimated extra-component bandwidth
  for adaptive bandwidth construction. Default is `NULL`.

- bwAdaptiveCrossover:

  numeric or NULL. Optional expression value at which the adaptive
  bandwidth curve crosses from the core bandwidth to the extra
  bandwidth. When `NULL`, component-density weighting is used. Default
  is `NULL`.

- bwAdaptiveTransitionWidth:

  numeric. Width, in expression units, of the optional smooth transition
  around `bwAdaptiveCrossover`. Use `0` for a hard switch at the
  crossover. Default is `0`.

- normPeakMinRel:

  numeric. Relative peak/trough threshold used to identify the main
  background modal complex for `*Norm` bandwidth methods. Default is
  `0.75`.

- normExtraFrac:

  numeric. Target fraction of additional high-side values sampled for
  the normalised high component. Default is `0.2`.

- normExtraMax:

  numeric. Maximum number of additional high-side values used by
  normalised bandwidth methods. May be `Inf`. Default is `Inf`.

- normLambda:

  numeric vector. Box-Cox lambda search grid used only when
  `normMtd = "boxcox"`. Default is `seq(-2, 2, length.out = 81)`.

- normDensityN:

  numeric. Number of grid points used inside normalised bandwidth
  helpers. Default is `512`.

- normExcessBwMtd:

  character. Ordinary bandwidth selector used for the right-side
  excess-density helper. Default is `"hpi3"`.

- normExcessNcell:

  numeric. Maximum number of cells used when estimating the
  excess-density helper bandwidth. Default is `10000`.

- normAdaptiveNcell:

  numeric. Fixed number of simulated normal-component values used per
  component when estimating adaptive normalised bandwidths. Default is
  `2500`.

- normMtd:

  character. Normalisation method for `*Norm` bandwidth selectors.
  `"moments"` replaces core and high components by normal components
  with matching moments; `"boxcox"` uses the older Box-Cox route for
  scalar bandwidths only. Default is `"moments"`.

## Value

A list of class `stimControl`.

## Details

Most users never need to change these settings. Arguments are grouped
into gating behaviour (`calcCytPosGates`), cell-count limits
(`minCell`), bias and expression limits (`biasUnsFactor`, `excMin`,
`cpMin`), bandwidth selection (`bwMtd` through `bwCluster`), threshold
sharing (`clusterGates`, `gateCombn`), local-FDR thresholding (`loc*`),
and advanced/experimental adaptive and normalised bandwidths
(`bwAdaptive*`, `norm*`). A fixed bandwidth is set on
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
through `bw`, or per marker through `markerControl`; when it is fixed
the bandwidth-selector settings are ignored. Use `markerControl` in
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
for per-marker overrides.

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
