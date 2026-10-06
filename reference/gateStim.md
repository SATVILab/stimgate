# Gate cells responding to stimulation

Compare cytokine expression in stimulated samples with an unstimulated
control from the same donor or batch. Save gates and
background-subtracted statistics to a project directory.

## Usage

``` r
gateStim(
  pathProject,
  .data,
  batchList,
  marker = NULL,
  chnl = NULL,
  popGate = "root",
  biasUns = NULL,
  bw = NULL,
  control = stimControl(),
  markerControl = NULL,
  parallel = FALSE
)
```

## Arguments

- pathProject:

  character Directory for results; created if needed.

- .data:

  GatingSet, flowSet, cytoset, flowFrame, cytoframe, character, list or
  data.frame Cytometry samples: a GatingSet or other flowCore/
  flowWorkspace object, FCS file paths or one FCS directory, a list of
  numeric matrices or data frames (cells by channels), or one data frame
  with a `sample` column. See Details.

- batchList:

  list Samples grouped by donor or batch, as indices or sample names,
  with the unstimulated control first in each vector, e.g.
  `list(donor1 = c(3, 1, 2))`. List names identify batches. Each batch
  needs at least one stimulated sample; a control may be shared by
  batches (first in each), but a stimulated sample may belong to only
  one.

- marker:

  character vector or NULL Marker labels to gate; supply either `marker`
  or `chnl`. Default: NULL.

- chnl:

  character vector or NULL Channel names to gate. Default: NULL.

- popGate:

  character Population(s) already present in `.data`; only GatingSets
  have populations other than "root". Default: "root" (all cells).

- biasUns:

  numeric or NULL Upward shift of unstimulated expression. NULL uses one
  quarter of `bwFallback`, scaled by `biasUnsFactor`. Positive shifts
  make gating more conservative. Default: NULL.

- bw:

  numeric or NULL Fixed density bandwidth; NULL estimates it
  automatically. Per-marker values go in `markerControl`. Default: NULL.

- control:

  stimControl Tuning settings from
  [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md).
  Default: stimControl().

- markerControl:

  list or NULL Settings for single markers, named by marker label or
  channel, e.g. `list(IL2 = list(bw = 0.12, biasUns = 0))`. Accepts
  [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
  settings except `locEnforceShapeThreshold` and `calcCytPosGates`, plus
  `bw`, `biasUns` and `popGate`. Default: NULL.

- parallel:

  logical Use the active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  for initial channel gating. Default: FALSE.

## Value

A character string: `pathProject`. Gates are saved under `gates/`,
expression under `sampleData/`, settings under `metaData/`, and
statistics in `gateStats.rds` and `gateStats.csv`.

## Details

Thresholds can be shared across similar distributions, then refined
using cells positive for another cytokine. Read results with
[`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md),
[`getStimStats()`](https://satvilab.github.io/stimgate/reference/getStimStats.md)
and
[`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md);
inspect them with
[`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md).

**Input data.** Non-GatingSet inputs are converted to a GatingSet with
only the "root" population. StimGate does not transform data (e.g.
arcsinh), so supply values on the scale to gate. Sample order, which
`batchList` indices refer to, is:

- FCS directory: `.fcs` files (any case, not recursive) in sorted order;
  a vector of paths keeps its order. Sample names are file basenames.

- List of matrices/data frames: list order; names are sample names
  (default `sample1`, `sample2`, ...). Columns must have the same unique
  names in every sample and serve as both channels and markers.

- Data frame with `sample`: factor-level order, otherwise first
  appearance.

Pass the same data, in the same order, to
[`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md),
[`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
and
[`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md).

For parallel gating, set `parallel = TRUE` and choose a future plan,
e.g. `future::plan(future::multisession, workers = 4)`. All workers must
be able to access `pathProject`. Later stages run sequentially. Set a
seed for reproducible parallel subsampling; results may differ from
sequential runs.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- gateStim(
  tempfile("stimgate_"), gs, exampleData$batchList,
  marker = exampleData$marker
)
#> shared bandwidth for MarkerF1: 0.329
#> shared bandwidth for MarkerF2: 0.334
#> getting base gates
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> getting cyt combn frequencies
#> batch 2 of 2
getStimGates(pathProject)
#> # A tibble: 4 × 12
#>   pop   gateName  chnl  marker ind   batch  gate locGenerated locGeneratedDirect
#>   <chr> <chr>     <chr> <I<ch> <chr> <chr> <dbl> <lgl>        <lgl>             
#> 1 root  locminCl… BC1(… Marke… 2     batc…  4.42 TRUE         TRUE              
#> 2 root  locminCl… BC1(… Marke… 4     batc…  3.79 TRUE         TRUE              
#> 3 root  locminCl… BC2(… Marke… 2     batc…  3.40 TRUE         TRUE              
#> 4 root  locminCl… BC2(… Marke… 4     batc…  2.91 TRUE         TRUE              
#> # ℹ 3 more variables: locSource <chr>, locReason <chr>, gateCyt <dbl>

# Disable gate sharing and fix the first marker's bandwidth
gateStim(
  tempfile("custom_gating_"), gs, exampleData$batchList,
  marker = exampleData$marker, control = stimControl(clusterGates = FALSE),
  markerControl = stats::setNames(
    list(list(bw = 0.12, biasUns = 0)), exampleData$marker[1]
  )
)
#> shared bandwidth for MarkerF2: 0.334
#> getting base gates
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> getting cyt combn frequencies
#> batch 2 of 2
#> [1] "/tmp/Rtmp8ehnUQ/custom_gating_1a7e268f9d5"

# Gate in-memory matrices; column names act as channels and markers
matrices <- lapply(seq_along(gs), function(i) {
  flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
})
gateStim(
  tempfile("matrix_gating_"), matrices, exampleData$batchList,
  chnl = exampleData$chnl, control = stimControl(calcCytPosGates = FALSE)
)
#> shared bandwidth for BC1(La139)Dd: 0.329
#> shared bandwidth for BC2(Pr141)Dd: 0.334
#> getting base gates
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 2 of 2
#> getting clustered and/or controlled gates
#> getting cyt combn frequencies
#> batch 2 of 2
#> [1] "/tmp/Rtmp8ehnUQ/matrix_gating_1a7e18ca17c5"
```
