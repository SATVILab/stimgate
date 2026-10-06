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
#> 1 root  locminCl… BC1(… Marke… 2     batc…  4.41 TRUE         TRUE              
#> 2 root  locminCl… BC1(… Marke… 4     batc…  3.79 TRUE         TRUE              
#> 3 root  locminCl… BC2(… Marke… 2     batc…  3.40 TRUE         TRUE              
#> 4 root  locminCl… BC2(… Marke… 4     batc…  2.91 TRUE         TRUE              
#> # ℹ 3 more variables: locSource <chr>, locReason <chr>, gateCyt <dbl>
```
