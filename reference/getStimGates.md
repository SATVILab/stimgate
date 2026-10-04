# Read stimulation gates

Read final gates saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
optionally selecting populations and markers. For more detail on how
gates were chosen, use
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

  character or NULL Populations to keep; NULL keeps all. Default: NULL.

- marker:

  character or NULL Marker labels to keep; takes precedence over `chnl`.
  Default: NULL (all markers).

- chnl:

  character or NULL Channels to keep. Default: NULL (all channels).

## Value

A tibble of stimulated-sample gates with identifiers `pop`, `marker`,
`chnl`, `batch`, `ind`, `gateName`, the gate value `gate`, and, when
available, the cytokine-positive gate `gateCyt` and columns recording
how each gate was found.

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
#> 1 root  locminClust BC1(La139)Dd MarkerF1 2     batch1  4.40    4.40
#> 2 root  locminClust BC1(La139)Dd MarkerF1 4     batch2  3.87    3.87
#> 3 root  locminClust BC2(Pr141)Dd MarkerF2 2     batch1  3.40    2.17
#> 4 root  locminClust BC2(Pr141)Dd MarkerF2 4     batch2  2.90    2.90
```
