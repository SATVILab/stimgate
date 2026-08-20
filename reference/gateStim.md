# Identify cytokine-positive cells through automated gating

Main function for identifying cytokine-positive cells using
outlier-based gating to compare stimulated versus unstimulated samples.
This function implements a comprehensive workflow that identifies cells
responding to stimulation by detecting outliers in cytokine expression
distributions. The process includes density estimation, threshold
identification, clustering-based gate refinement, and generation of
comprehensive statistics and visualizations.

The function operates by comparing cytokine expression in stimulated
samples against corresponding unstimulated controls from the same
donor/batch to identify cells that have likely responded to stimulation.
It accounts for batch effects, background cytokine production, and
technical variability.

## Usage

``` r
gateStim(
  pathProject,
  .data,
  popGate = "root",
  batchList,
  chnl = NULL,
  marker = NULL,
  calcCytPosGates = TRUE,
  biasUns = NULL,
  biasUnsFactor = 1,
  excMin = TRUE,
  cpMin = NULL,
  bw = NULL,
  bwMin = "auto",
  bwMax = "auto",
  bwFallback = "auto",
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 100,
  bwNcellMax = 1e+05,
  bwCluster = NULL,
  bwAdaptive = FALSE,
  bwAdaptiveDensityN = NULL,
  bwAdaptivePadFrac = 0.15,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normPeakFrac = 0.1,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  normMtd = "moments",
  minCell = 100,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  tolClust = 1e-07,
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locEnforceShapeThreshold = FALSE,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1/6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1/2,
  locFlatHardDerivFrac = 1/4,
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
)
```

## Arguments

- pathProject:

  character. Path to project directory where all results will be saved.
  This directory will contain subdirectories for each marker with gate
  tables, statistics, and plots. The directory will be created if it
  doesn't exist.

- .data:

  GatingSet. A flowWorkspace GatingSet object containing the flow
  cytometry data with both stimulated and unstimulated samples. The
  GatingSet should have consistent channel names across all samples and
  include proper sample annotations.

- popGate:

  character vector. Population(s) within which to perform gating.
  Default is "root" to gate on all cells. Can specify other populations
  like "CD3+" or "CD4+" if these gates already exist in the GatingSet.

- batchList:

  list. List where each element contains indices of samples belonging to
  the same batch/donor. The first index per element is the unstimulated
  control sample, e.g. if `batchList = list(c(3, 1, 2), c(6, 4, 5))`,
  then indices 3 and 6 correspond to the unstimulated samples for
  batches 1 and 2, respectively. If `batchList` is named, e.g.
  `list(pid1 = c(3, 1, 2), pid2 = c(6, 4, 5))`, then these names will be
  used for batch identification.

- chnl:

  list. List where each element specifies parameters for gating a
  specific channel Each element should be a list containing at minimum
  the channel name (e.g., list(cut = "IL2")). Additional
  channel-specific parameters can override global defaults. The marker
  name should match channel names in the GatingSet.

- marker:

  character vector. Alternative way to specify markers to gate on. When
  provided, this is used instead of chnl to determine which markers to
  analyze. Default is NULL.

- calcCytPosGates:

  logical. Whether to refine each clustered one-marker gate using the
  target-marker distribution among cells positive for at least one other
  cytokine. A taut-string density is fitted to those cells. The
  clustered gate is lowered to the leftmost internal antimode strictly
  between the full stimulated marginal peak plus one third of its
  left-window width and the clustered gate. If no eligible antimode
  exists, the clustered gate is retained. Default is TRUE.

- biasUns:

  numeric. Bias adjustment for unstimulated samples to account for
  background cytokine production. When NULL (default), no bias
  correction is applied. Positive values shift the unstimulated
  distribution higher, making gates more conservative.

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

- bw:

  numeric. Specify the bandwith for density estimation. When NULL
  (default), bandwidth is estimated automatically. Default is `NULL`.

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

- bwAdj:

  numeric. Adjustment factor for bandwidth. Default is 1. Ignored if
  `bw` is set. Default is 1.

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

  numeric. Optional fallback bandwidth for cluster-based local-FDR
  refinement. The cluster step first tries to use the median bandwidth
  across samples with directly generated local-FDR thresholds.
  `bwCluster` is used when that common bandwidth cannot be estimated.
  Default is `NULL`.

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

- normPeakFrac:

  numeric. Fraction of the selected background-core peak height used by
  normalised bandwidth helpers when identifying low-density tail
  regions. Default is `0.1`.

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

- normExtraJitterFrac:

  numeric. Jitter scale, as a fraction of a robust expression-scale
  standard deviation, applied to sampled high-side values. Default is
  `0.25`.

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

- minCell:

  numeric. Minimum number of cells required for reliable gating. Default
  is 100. Samples with fewer cells will be skipped as they don't provide
  sufficient statistical power for accurate gate identification.

- maxPosProbX:

  numeric. Maximum x-value (expression level) to consider when
  calculating positive probabilities. Default is Inf (no limit). Can be
  used to exclude extremely high expression values that may represent
  doublets or artifacts.

- gateQuant:

  numeric vector. Quantiles used for gate combination when multiple
  gates are identified. Default is c(0.25, 0.75). The method specified
  in gateCombn determines how these quantiles are used (e.g., minimum of
  25th percentiles).

- tolClust:

  numeric or NULL. Backwards-compatible switch for calculating
  cluster-adjusted local-FDR gates. When `NULL`, cluster adjustment is
  skipped. When non-NULL, paired stimulated and unstimulated densities
  are clustered on a common absolute-expression grid. Direct thresholds
  are winsorised within each cluster to its 15th and 85th percentiles
  when at least three direct thresholds are available. Every non-direct
  threshold is replaced by the cluster's 60th percentile when at least
  one direct threshold is available. A cluster without a direct
  threshold retains its original high thresholds. The numeric value is
  retained only as a backwards-compatible on/off switch.

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

- locLeftLowRel:

  numeric. Retained for backwards compatibility. The current
  derivative/marginal left-tail trim no longer uses this overall-region
  check. Candidate left-tail regions not separated by an antimode are
  considered low response when their mean response probability is below
  this fraction of the peak response probability.

- locLeftLowAbs:

  numeric. Retained for backwards compatibility. The current
  derivative/marginal left-tail trim no longer uses this overall-region
  check. Absolute response-probability cutoff for the non-antimode
  left-tail trimming rule.

- locLeftCellFrac:

  numeric. Retained for backwards compatibility. The current
  derivative/marginal left-tail trim no longer uses this overall-region
  check. Minimum size of the candidate low-response left-tail region,
  expressed as a fraction of the number of cells to the right of the
  start of the main probability rise.

- locLeftLengthFrac:

  numeric. Retained for backwards compatibility. The current
  derivative/marginal left-tail trim no longer uses this overall-region
  check. Minimum length of the candidate low-response left-tail region,
  expressed as a fraction of the expression interval over which the
  response probability rises from its minimum to its maximum.

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

- locTolRefPeak:

  character. Deprecated clustering setting retained for backwards
  compatibility. Joint-density quantile transfer does not use a
  derivative-tolerance reference peak. Default is `"highest"`.

- gateCombn:

  character vector. Method(s) for combining condition-level local-FDR
  gates within a batch. Supported values are `"no"`, `"min"`,
  `"median"`, `"max"`, and `"prejoin"`. Combination uses only thresholds
  that were actually generated by the local-FDR procedure, not fallback
  above-range cutpoints.

- markerSettings:

  list. Optional list of additional marker-specific settings that
  override global defaults. Each element should be named with the marker
  name and contain parameter overrides. Default is NULL.

- chnlSettings:

  list. Optional list of channel-specific settings that override global
  defaults. Similar to markerSettings but keyed by channel names.
  Default is NULL.

## Value

character. Returns the path to the project directory where all results
have been saved. The directory structure created includes:

- `pathProject/[markerName]/`: Directory for each marker containing:

- `gateTblInit.rds`: Initial gate table with preliminary gates

- `gateTbl.rds`: Final refined gate table

- `stats/`: Directory containing statistics files

- `plots/`: Directory containing visualization plots (if generated)

## Details

The function implements a multi-step workflow for identifying
cytokine-positive cells:

**Step 1: Data Preparation**

- Validates input parameters and GatingSet structure

- Completes marker specifications with default values

- Organizes samples by batch for proper background subtraction

**Step 2: Initial Gate Identification**

- Extracts expression data for each marker within specified populations

- Estimates probability densities for stimulated and unstimulated
  samples

- Identifies candidate cutpoints using outlier detection algorithms

- Applies clustering to refine gate positions across batches

**Step 3: Cytokine-Positive Gate Refinement (if calcCytPosGates =
TRUE)**

- Fits a taut-string density to target-marker expression among cells
  positive for another cytokine

- Finds internal antimodes between the protected negative-component
  region and the clustered one-marker gate

- Lowers the gate to the leftmost eligible antimode, or retains the
  clustered gate when none is available

**Step 4: Statistics Generation**

- Calculates comprehensive statistics including frequencies and
  combinations

- Generates cross-tabulations of cytokine-positive populations

- Saves results in structured format for downstream analysis

## See also

[`getStimGates`](https://satvilab.github.io/stimgate/reference/getStimGates.md)
for extracting gate information,
[`getStimStats`](https://satvilab.github.io/stimgate/reference/getStimStats.md)
for generating statistics from results,
[`plotStim`](https://satvilab.github.io/stimgate/reference/plotStim.md)
for visualizing identified gates,
[`writeStimFCS`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
for exporting cytokine-positive cells,
[`GatingSet`](https://rdrr.io/pkg/flowWorkspace/man/GatingSet-class.html)
for GatingSet documentation

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- file.path(tempdir(), "demonstration")

# Run gating
gateStim(
  .data = gs,
  pathProject = pathProject,
  popGate = "root",
  batchList = exampleData$batchList,
  marker = exampleData$marker
)
#> ----
#> getting base gates
#> ----
#> 
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> 
#> 
#> 
#> getting cyt combn frequencies
#> batch 2 of 2
#> [1] "/tmp/Rtmpl3WoLo/demonstration"

# Create plots
if (requireNamespace("hexbin", quietly = TRUE)) {
  plots <- plotStim(
    ind = exampleData$batchList[[1]], # indices in `gs` to plot
    .data = gs, # GatingSet
    pathProject = pathProject,
    marker = exampleData$marker,
    grid = TRUE
  )
}
```
