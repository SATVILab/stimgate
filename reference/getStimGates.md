# Read stimulation gates

Read final gates saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
optionally selecting populations and markers. For threshold diagnostics,
use
[`getStimGatesDetailed()`](https://satvilab.github.io/stimgate/reference/getStimGatesDetailed.md).

## Usage

``` r
getStimGates(pathProject, pop = NULL, marker = NULL, chnl = NULL)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- pop:

  character or NULL Populations to retain; NULL selects all. Default:
  NULL.

- marker:

  character or NULL Marker labels to retain; takes precedence over
  `chnl`. Default: NULL (all markers).

- chnl:

  character or NULL Channels to retain. Default: NULL (all channels).

## Value

A tibble of stimulated-sample gates with identifiers `pop`, `marker`,
`chnl`, `batch`, `ind`, `gateName`, threshold `gate`, and refinement and
threshold-provenance columns when available.

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
#> # A tibble: 4 × 8
#>   pop   gateName    chnl         marker   ind   batch   gate gateCyt
#>   <chr> <chr>       <chr>        <I<chr>> <chr> <chr>  <dbl>   <dbl>
#> 1 root  locminClust BC1(La139)Dd MarkerF1 2     batch1  4.42    4.42
#> 2 root  locminClust BC1(La139)Dd MarkerF1 4     batch2  3.79    3.79
#> 3 root  locminClust BC2(Pr141)Dd MarkerF2 2     batch1  3.40    2.17
#> 4 root  locminClust BC2(Pr141)Dd MarkerF2 4     batch2  2.91    2.91
```
