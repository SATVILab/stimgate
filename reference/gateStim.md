# Identify cytokine-positive cells through automated gating

Identify cells responding to stimulation by comparing cytokine
expression in stimulated samples with unstimulated controls from the
same donor/batch. Saves gates and background-subtracted statistics to a
project directory.

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

  character. Path to project directory where all results will be saved.
  This directory will contain subdirectories for each marker with gate
  tables, statistics, and plots. The directory will be created if it
  doesn't exist.

- .data:

  GatingSet. A flowWorkspace GatingSet object containing the flow
  cytometry data with both stimulated and unstimulated samples. The
  GatingSet should have consistent channel names across all samples and
  include proper sample annotations.

- batchList:

  list. List where each element contains indices of samples belonging to
  the same batch/donor. The first index per element is the unstimulated
  control sample, e.g. if `batchList = list(c(3, 1, 2), c(6, 4, 5))`,
  then indices 3 and 6 correspond to the unstimulated samples for
  batches 1 and 2, respectively. If `batchList` is named, e.g.
  `list(pid1 = c(3, 1, 2), pid2 = c(6, 4, 5))`, then these names will be
  used for batch identification.

- marker:

  character vector. Alternative way to specify markers to gate on. When
  provided, this is used instead of chnl to determine which markers to
  analyze. Default is NULL.

- chnl:

  character vector. Channel names to gate on. Specify either `chnl` or
  `marker`. Default is NULL.

- popGate:

  character vector. Population(s) within which to perform gating.
  Default is "root" to gate on all cells. Can specify other populations
  like "CD3+" or "CD4+" if these gates already exist in the GatingSet.

- biasUns:

  numeric. Bias adjustment for unstimulated samples to account for
  background cytokine production. When NULL (default), 1/4 of
  `bwFallback` is used (scaled by `biasUnsFactor`). Positive values
  shift the unstimulated distribution higher, making gates more
  conservative. Default is NULL.

- bw:

  numeric. Specify the bandwith for density estimation. When NULL
  (default), bandwidth is estimated automatically. A bandwidth may also
  be set per marker through `markerControl`. Default is `NULL`.

- control:

  stimControl Tuning settings from
  [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md).
  Most users do not need to change these. Default:
  [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md).

- markerControl:

  list or NULL. Named per-marker overrides, keyed by marker label or
  channel name, for example `list(IL2 = list(bw = 0.12, biasUns = 0))`.
  Settings from
  [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
  (except the global-only `locEnforceShapeThreshold` and
  `calcCytPosGates`), plus `bw`, `biasUns` and `popGate`, can be
  overridden. Default: NULL.

- parallel:

  logical If TRUE, gate channels in parallel during the initial gating
  stage using the active future::plan(). Default: FALSE.

## Value

character. Returns the path to the project directory where all results
have been saved. The directory structure created includes:

- `pathProject/[markerName]/`: Directory for each marker containing:

- `gateTblInit.rds`: Initial gate table with preliminary gates

- `gateTbl.rds`: Final refined gate table

- `stats/`: Directory containing statistics files

- `plots/`: Directory containing visualization plots (if generated)

## Details

Initial local-FDR gates compare each stimulated sample with its batch's
unstimulated control. Thresholds may be shared across similar
distributions, then refined using cells positive for another cytokine.
Results and statistics are saved to `pathProject` for
[`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md),
[`getStimStats()`](https://satvilab.github.io/stimgate/reference/getStimStats.md),
and
[`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md).
Use
[`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
for tuning and `markerControl` for per-marker overrides.

To gate channels in parallel, set `parallel = TRUE` and select a future
plan, for example `future::plan(future::multisession, workers = 4)`. The
default `parallel = FALSE` runs sequentially regardless of the active
plan. Only the initial per-channel gating stage is parallel; subsequent
cytokine-positive gating and statistics remain sequential. Workers read
expression data from the project disk cache rather than a GatingSet, so
the project directory must be accessible to all workers. With
`parallel = TRUE`, RNG-dependent subsampling uses parallel-safe L'Ecuyer
streams (`future.seed = TRUE`). Results are reproducible for a given
[`set.seed()`](https://rdrr.io/r/base/Random.html) and independent of
the chosen non-sequential plan, but may differ slightly from a
sequential run.

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
#> [1] "/tmp/Rtmp29tuQd/demonstration"

# Customise tuning and override the bandwidth for the first marker
gateStim(
  pathProject = file.path(tempdir(), "custom-gating"),
  .data = gs,
  batchList = exampleData$batchList,
  marker = exampleData$marker,
  bw = 0.1,
  control = stimControl(bwAdj = 1.5, clusterGates = FALSE),
  markerControl = stats::setNames(
    list(list(bw = 0.12, biasUns = 0)), exampleData$marker[1]
  )
)
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
#> [1] "/tmp/Rtmp29tuQd/custom-gating"

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
