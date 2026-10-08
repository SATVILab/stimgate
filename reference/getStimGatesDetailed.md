# Read details of how gates were chosen

Read the intermediate gates saved at each step (per condition, per
sample and per cluster) with background-subtracted frequencies. Use
[`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md)
for final gates.

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

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- pop:

  character or NULL Populations to keep. Default: NULL (all).

- marker:

  character or NULL Marker labels to keep. Default: NULL (all).

- chnl:

  character or NULL Channels to keep. Default: NULL (all).

- save:

  logical Save the returned table as RDS. Default: FALSE.

- pathSave:

  character or NULL Output path; NULL uses
  `file.path(pathProject, "gatesDetailed.rds")`. Default: NULL.

## Value

A tibble with one row per saved gate, including `pop`, `marker`, `chnl`,
gate and frequency columns, and the file each row was read from.

## Details

To save these details, run `Sys.setenv(STIMGATE_INTERMEDIATE = "all")`
before
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).
If nothing was saved, an empty tibble is returned.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
Sys.setenv(STIMGATE_INTERMEDIATE = "all")
pathProject <- gateStim(
  tempfile("stimgate_"), gs, exampleData$batchList,
  marker = exampleData$marker
)
#> shared bandwidth for MarkerF1: 0.188
#> shared bandwidth for MarkerF2: 0.19
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
Sys.unsetenv("STIMGATE_INTERMEDIATE")
getStimGatesDetailed(pathProject)
#> # A tibble: 12 × 60
#>    pop   marker   chnl         detailLevel stage ind   threshold thresholdOrigin
#>    <chr> <I<chr>> <chr>        <chr>       <chr> <chr>     <dbl> <chr>          
#>  1 root  MarkerF1 BC1(La139)Dd condition   init  2          4.24 condition_dete…
#>  2 root  MarkerF1 BC1(La139)Dd sample      init  2          4.24 condition_dete…
#>  3 root  MarkerF1 BC1(La139)Dd condition   init  4          3.68 condition_dete…
#>  4 root  MarkerF1 BC1(La139)Dd sample      init  4          3.68 condition_dete…
#>  5 root  MarkerF1 BC1(La139)Dd cluster_fi… NA    2          4.24 NA             
#>  6 root  MarkerF1 BC1(La139)Dd cluster_fi… NA    4          3.68 NA             
#>  7 root  MarkerF2 BC2(Pr141)Dd condition   init  2          3.99 condition_dete…
#>  8 root  MarkerF2 BC2(Pr141)Dd sample      init  2          3.99 condition_dete…
#>  9 root  MarkerF2 BC2(Pr141)Dd condition   init  4          2.78 condition_dete…
#> 10 root  MarkerF2 BC2(Pr141)Dd sample      init  4          2.78 condition_dete…
#> 11 root  MarkerF2 BC2(Pr141)Dd cluster_fi… NA    2          3.99 NA             
#> 12 root  MarkerF2 BC2(Pr141)Dd cluster_fi… NA    4          2.78 NA             
#> # ℹ 52 more variables: locGenerated <lgl>, locGeneratedDirect <lgl>,
#> #   locSource <chr>, locReason <chr>, bias <dbl>, propBsEst <dbl>,
#> #   propBsDiff <dbl>, nCellStim <int>, nCellUns <int>, propStim <dbl>,
#> #   propUns <dbl>, propBs <dbl>, locThresholdMethod <chr>, locRegionX <dbl>,
#> #   detailObject <chr>, detailPathStage <chr>, detailPathChnl <chr>,
#> #   detailPathInd <chr>, detailPathPop <chr>, detailSourceFile <chr>,
#> #   grp <chr>, grpUns <chr>, grpStim <chr>, cpOrigQuantMin <dbl>, …
```
