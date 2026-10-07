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
#> 1 loc_min… 2     BC1(La1…        39     10000        3    10000   0.0039  0.0003
#> 2 loc_min… 2     BC1(La1…       262     10000       62    10000   0.0262  0.0062
#> 3 loc_min… 2     BC1(La1…        40     10000        1    10000   0.004   0.0001
#> 4 loc_min… 2     BC1(La1…      9659     10000     9934    10000   0.966   0.993 
#> 5 loc_min… 4     BC1(La1…        92     10000       23    10000   0.0092  0.0023
#> 6 loc_min… 4     BC1(La1…       484     10000      112    10000   0.0484  0.0112
#> 7 loc_min… 4     BC1(La1…        47     10000        6    10000   0.0047  0.0006
#> 8 loc_min… 4     BC1(La1…      9377     10000     9859    10000   0.938   0.986 
#> # ℹ 4 more variables: propBs <dbl>, freqStim <dbl>, freqUns <dbl>, freqBs <dbl>
```
