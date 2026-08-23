# Plot stimulation gate

Plot bivariate hex and univariate density plots for batches of samples,
along with their gates.

## Usage

``` r
plotStim(
  ind,
  .data,
  pathProject,
  marker = NULL,
  chnl = NULL,
  pop = NULL,
  indLab = NULL,
  axisLab = NULL,
  excMin = TRUE,
  limitsExpand = NULL,
  limitsEqual = FALSE,
  grid = TRUE,
  gridNCol = 2,
  showGate = TRUE,
  minCell = 10,
  bias = FALSE,
  combnExc = NULL,
  chnlGate = NULL,
  markerGate = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE,
  gateUnsMethod = "min"
)
```

## Arguments

- ind:

  numeric vector. Specifies indices in `.data` to plot.

- .data:

  GatingSet. Same GatingSet passed to `gateStim`.

- pathProject:

  character. Path to the project directory used for `gateStim`.

- marker:

  character vector of length one or two. Specifies markers to be
  plotted. If only one is passed, then only univariate plots are
  created.

- chnl:

  character vector of length one or two. Specifies channels to be
  plotted. Ignored if `marker` is provided.

- pop:

  character. Specifies population within GatingSet that gates were
  calculated on. If `NULL`, defaults to population specified by folder
  name in `project_path/gates/pop_<pop>`, but throws an error if more
  than one population is detected (i.e. more than one directory in
  `gates/`). Default is `NULL`.

- indLab:

  named character vector. Labels for `ind` used in plot. Optional.

- axisLab:

  named character vector. Labels for axis titles, applied to `marker` or
  `chnl`. Optional.

- excMin:

  Logical. If `TRUE`, excludes the minimum expression values when
  processing the data. Default is `TRUE`.

- limitsExpand:

  list. Expand the limits of the plot axes. Default is `NULL`.

- limitsEqual:

  Logical. If TRUE, forces equal lengths of the limits.

- grid:

  Logical. If TRUE, arranges the resulting plots in a grid format using
  [`cowplot::plot_grid`](https://wilkelab.org/cowplot/reference/plot_grid.html).
  Default is `TRUE`.

- gridNCol:

  Integer. Number of columns in grid layout.

- showGate:

  Logical. If `TRUE`, overlays gate lines on the plots. Default is
  `TRUE`.

- minCell:

  integer. Minimum number of cells to be plotted. Will skip plots with
  fewer cells. Default is 10.

- bias:

  logical Whether to add bias to unstimulated sample used in the gating.
  Default is `FALSE`.

- combnExc:

  list or NULL Combinations of channels to exclude. Default is NULL.

- chnlGate:

  character or NULL Channel name(s) to use for gating. Cannot be
  specified with `marker_gate`. Default is NULL.

- markerGate:

  character or NULL Marker name(s) to use for gating. Cannot be
  specified with `chnl_gate`. Default is NULL.

- gateTypeCytPos:

  character Gate type to use for cytokine-positive cells. Default is
  "cyt".

- mult:

  logical Whether to return only multi-functional cells (positive for
  multiple markers). Default is FALSE.

- gateUnsMethod:

  character Method for gating unstimulated cells. Default is "min".

## Value

A grid of plots if `grid` is TRUE, otherwise a list of ggplot objects.

## Examples

``` r
# Create example data and run gating
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")

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
#> [1] "/tmp/RtmpIO1tZy/stimgate_example_data_4aa97266c797/stimgate"

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
