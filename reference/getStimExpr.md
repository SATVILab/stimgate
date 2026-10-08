# Read cell expression values

Read expression saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
optionally selecting stimulation-positive cells. Supply marker labels or
channel names.

## Usage

``` r
getStimExpr(
  pathProject,
  .data = NULL,
  pop = NULL,
  ind = NULL,
  chnl = NULL,
  marker = NULL,
  bias = FALSE,
  excMin = FALSE,
  combnExc = NULL,
  chnlGate = NULL,
  markerGate = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE,
  transFn = NULL,
  transChnl = NULL,
  transMarker = NULL
)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- .data:

  GatingSet, other input accepted by
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
  or NULL Data passed to
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
  in the same sample order; used only when the expression values were
  not saved. Default: NULL.

- pop:

  character or NULL Population names; NULL selects all saved
  populations. Default: NULL.

- ind:

  character or numeric vector or NULL Sample indices; NULL selects all
  saved samples, including controls. Default: NULL.

- chnl:

  character or NULL Channels to return; NULL selects all saved channels
  unless `marker` is supplied. Default: NULL.

- marker:

  character or NULL Marker labels to return; cannot be combined with
  `chnl`. Default: NULL.

- bias:

  logical Add the saved `biasUns` shift to controls. Default: FALSE.

- excMin:

  logical Exclude cells at the minimum of any requested channel.
  Default: FALSE.

- combnExc:

  list or NULL Channel combinations to exclude: each vector specifies
  positive channels, with other gating channels negative. Applies only
  when gating channels are supplied. Default: NULL.

- chnlGate:

  character or NULL Channels used to select positive cells; include
  these in the requested expression columns. Default: NULL.

- markerGate:

  character or NULL Marker labels used to select positive cells; cannot
  be combined with `chnlGate`. Default: NULL.

- gateTypeCytPos:

  character Positivity rule: "base" uses the main gate; "cyt" also
  admits cells above a refined gate when another marker clears its main
  gate. Default: "cyt".

- mult:

  logical Require positivity for at least two gating markers. Applies
  only when `chnlGate` or `markerGate` is supplied. Default: FALSE.

- transFn:

  function or NULL Transformation applied to the expression tibble
  before adding metadata columns. Default: NULL.

- transChnl:

  character or NULL Columns to transform when using channels; NULL
  transforms all expression columns. Default: NULL.

- transMarker:

  character or NULL Columns to transform when using markers; NULL
  transforms all expression columns. Default: NULL.

## Value

A tibble with one row per selected cell, `pop`, `ind`, and expression
columns named by channel (or marker when `marker` is supplied). Empty
selections have zero rows. The `nCellPos` attribute is a tibble with
`pop`, `ind`, `nCellPos` for every requested population/sample pair;
`probGMin` records the fraction of cells kept after removing
minimum-expression cells.

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
#> shared bandwidth for MarkerF1: 0.188
#> shared bandwidth for MarkerF2: 0.188
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
getStimExpr(pathProject, marker = exampleData$marker)
#> # A tibble: 40,000 × 4
#>    pop   ind   MarkerF1 MarkerF2
#>    <chr> <chr>    <dbl>    <dbl>
#>  1 root  1      -0.443     0.123
#>  2 root  1       0.635     1.12 
#>  3 root  1       1.06     -0.607
#>  4 root  1       0.0511   -1.92 
#>  5 root  1      -0.186     1.19 
#>  6 root  1       1.49     -0.438
#>  7 root  1      -0.479    -1.10 
#>  8 root  1       2.85      0.354
#>  9 root  1       0.721     1.69 
#> 10 root  1      -1.22      0.713
#> # ℹ 39,990 more rows
```
