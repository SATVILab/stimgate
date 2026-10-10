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
#> 2 loc_min… 2     BC1(La1…       265     10000       67    10000   0.0265  0.0067
#> 3 loc_min… 2     BC1(La1…        41     10000        1    10000   0.0041  0.0001
#> 4 loc_min… 2     BC1(La1…      9655     10000     9929    10000   0.966   0.993 
#> 5 loc_min… 4     BC1(La1…        61     10000       14    10000   0.0061  0.0014
#> 6 loc_min… 4     BC1(La1…       455     10000       97    10000   0.0455  0.0097
#> 7 loc_min… 4     BC1(La1…        35     10000        4    10000   0.0035  0.0004
#> 8 loc_min… 4     BC1(La1…      9449     10000     9885    10000   0.945   0.988 
#> # ℹ 4 more variables: propBs <dbl>, freqStim <dbl>, freqUns <dbl>, freqBs <dbl>
```
