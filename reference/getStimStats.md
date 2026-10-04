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
#> 1 locminC… 2     BC1(La1…        51     10000        6    10000   0.0051  0.0006
#> 2 locminC… 2     BC1(La1…       425     10000      154    10000   0.0425  0.0154
#> 3 locminC… 2     BC1(La1…        47     10000        2    10000   0.0047  0.0002
#> 4 locminC… 2     BC1(La1…      9477     10000     9838    10000   0.948   0.984 
#> 5 locminC… 4     BC1(La1…        94     10000       25    10000   0.0094  0.0025
#> 6 locminC… 4     BC1(La1…       516     10000      119    10000   0.0516  0.0119
#> 7 locminC… 4     BC1(La1…        50     10000        6    10000   0.005   0.0006
#> 8 locminC… 4     BC1(La1…      9340     10000     9850    10000   0.934   0.985 
#> # ℹ 4 more variables: propBs <dbl>, freqStim <dbl>, freqUns <dbl>, freqBs <dbl>
```
