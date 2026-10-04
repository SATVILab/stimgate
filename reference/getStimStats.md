# Get gating statistics

Read and return gating statistics computed during gating.

## Usage

``` r
getStimStats(pathProject)
```

## Arguments

- pathProject:

  character. Path to the project directory.

## Value

A data frame with gating statistics.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- gateStim(
  pathProject = file.path(tempdir(), "getStimStatsExample"),
  .data = gs,
  batchList = exampleData$batchList,
  marker = exampleData$marker,
  popGate = "root"
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

# Get gating statistics
statTbl <- getStimStats(pathProject)
```
