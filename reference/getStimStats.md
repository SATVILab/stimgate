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
#> shared bandwidth for MarkerF1: 0.187
#> shared bandwidth for MarkerF2: 0.191
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
#> 1 loc_min… 2     BC1(La1…        38     10000        3    10000   0.0038  0.0003
#> 2 loc_min… 2     BC1(La1…       310     10000       80    10000   0.031   0.008 
#> 3 loc_min… 2     BC1(La1…        39     10000        1    10000   0.0039  0.0001
#> 4 loc_min… 2     BC1(La1…      9613     10000     9916    10000   0.961   0.992 
#> 5 loc_min… 4     BC1(La1…        61     10000       17    10000   0.0061  0.0017
#> 6 loc_min… 4     BC1(La1…       456     10000      100    10000   0.0456  0.01  
#> 7 loc_min… 4     BC1(La1…        35     10000        4    10000   0.0035  0.0004
#> 8 loc_min… 4     BC1(La1…      9448     10000     9879    10000   0.945   0.988 
#> # ℹ 4 more variables: propBs <dbl>, freqStim <dbl>, freqUns <dbl>, freqBs <dbl>
```
