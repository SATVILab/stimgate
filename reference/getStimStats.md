# Read gating statistics

Read cell counts and background-subtracted frequencies saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Usage

``` r
getStimStats(pathProject)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Value

A tibble (or data.frame when read from CSV) with sample and gate
identifiers, `countStim`, `countUns`, `nCellStim`, `nCellUns`,
proportions `propStim`, `propUns`, `propBs`, and percentages `freqStim`,
`freqUns`, `freqBs`. Background subtraction is stimulated minus
unstimulated.

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
#> shared bandwidth for MarkerF1: 0.185
#> shared bandwidth for MarkerF2: 0.189
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
getStimStats(pathProject)
#> # A tibble: 8 × 13
#>   gateName ind   cytCombn countStim nCellStim countUns nCellUns propStim propUns
#>   <chr>    <chr> <chr>        <int>     <int>    <int>    <int>    <dbl>   <dbl>
#> 1 loc_min… 2     BC1(La1…        45     10000        5    10000   0.0045  0.0005
#> 2 loc_min… 2     BC1(La1…         5     10000        0    10000   0.0005  0     
#> 3 loc_min… 2     BC1(La1…        45     10000        2    10000   0.0045  0.0002
#> 4 loc_min… 2     BC1(La1…      9905     10000     9993    10000   0.990   0.999 
#> 5 loc_min… 4     BC1(La1…       109     10000       33    10000   0.0109  0.0033
#> 6 loc_min… 4     BC1(La1…       560     10000      136    10000   0.056   0.0136
#> 7 loc_min… 4     BC1(La1…        59     10000        9    10000   0.0059  0.0009
#> 8 loc_min… 4     BC1(La1…      9272     10000     9822    10000   0.927   0.982 
#> # ℹ 4 more variables: propBs <dbl>, freqStim <dbl>, freqUns <dbl>, freqBs <dbl>
```
