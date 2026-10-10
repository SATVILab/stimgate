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
#> # A tibble: 20 × 66
#>    pop   marker   chnl        ind   detailLevel gateCombn threshold locGenerated
#>    <chr> <I<chr>> <chr>       <chr> <chr>       <chr>         <dbl> <lgl>       
#>  1 root  MarkerF1 BC1(La139)… 2     batch_share min            4.67 TRUE        
#>  2 root  MarkerF1 BC1(La139)… 1     batch_share min            4.67 TRUE        
#>  3 root  MarkerF1 BC1(La139)… 2     condition   NA             4.67 TRUE        
#>  4 root  MarkerF1 BC1(La139)… 2     sample      NA             4.67 TRUE        
#>  5 root  MarkerF1 BC1(La139)… 4     batch_share min            4.29 TRUE        
#>  6 root  MarkerF1 BC1(La139)… 3     batch_share min            4.29 TRUE        
#>  7 root  MarkerF1 BC1(La139)… 4     condition   NA             4.29 TRUE        
#>  8 root  MarkerF1 BC1(La139)… 4     sample      NA             4.29 TRUE        
#>  9 root  MarkerF1 BC1(La139)… 2     cluster_fi… NA             4.67 TRUE        
#> 10 root  MarkerF1 BC1(La139)… 4     cluster_fi… NA             4.29 TRUE        
#> 11 root  MarkerF2 BC2(Pr141)… 2     batch_share min            4.09 TRUE        
#> 12 root  MarkerF2 BC2(Pr141)… 1     batch_share min            4.09 TRUE        
#> 13 root  MarkerF2 BC2(Pr141)… 2     condition   NA             4.09 TRUE        
#> 14 root  MarkerF2 BC2(Pr141)… 2     sample      NA             4.09 TRUE        
#> 15 root  MarkerF2 BC2(Pr141)… 4     batch_share min            3.40 TRUE        
#> 16 root  MarkerF2 BC2(Pr141)… 3     batch_share min            3.40 TRUE        
#> 17 root  MarkerF2 BC2(Pr141)… 4     condition   NA             3.40 TRUE        
#> 18 root  MarkerF2 BC2(Pr141)… 4     sample      NA             3.40 TRUE        
#> 19 root  MarkerF2 BC2(Pr141)… 2     cluster_fi… NA             4.09 TRUE        
#> 20 root  MarkerF2 BC2(Pr141)… 4     cluster_fi… NA             3.40 TRUE        
#> # ℹ 58 more variables: locGeneratedDirect <lgl>, locSource <chr>,
#> #   locReason <chr>, locResponder <lgl>, propBsEst <dbl>, locOwnFreq <dbl>,
#> #   locShareLimit <chr>, locShareProposed <dbl>, detailObject <chr>,
#> #   detailPathStage <chr>, detailPathChnl <chr>, detailPathInd <chr>,
#> #   detailPathPop <chr>, detailSourceFile <chr>, stage <chr>,
#> #   thresholdOrigin <chr>, bias <dbl>, propBsDiff <dbl>, nCellStim <int>,
#> #   nCellUns <int>, propStim <dbl>, propUns <dbl>, propBs <dbl>, …
```
