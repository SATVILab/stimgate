# Get detailed gate diagnostics

Read detailed local-FDR threshold diagnostics saved during gating,
including condition-level, sample-level and final cluster-level
thresholds with the corresponding background-subtracted frequencies.

## Usage

``` r
getStimGatesDetailed(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL,
  save = FALSE,
  pathSave = NULL
)
```

## Arguments

- pathProject:

  character. Path to the project directory.

- pop:

  character. Optional population name(s) to filter gates by. Default is
  NULL (all populations).

- marker:

  character. Optional marker name(s) to retain.

- chnl:

  character. Optional channel name(s) to retain.

- save:

  logical. If TRUE, save the detailed table as an RDS file.

- pathSave:

  character. Optional path for the saved RDS file. Defaults to
  `file.path(pathProject, "gatesDetailed.rds")`.

## Value

A tibble with one row per saved threshold diagnostic.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- gateStim(
  pathProject = file.path(tempdir(), "getStimGatesDetailedExample"),
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

# Get threshold diagnostics
detailTbl <- getStimGatesDetailed(pathProject)
```
