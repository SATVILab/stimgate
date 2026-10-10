# Read pairwise cytokine coexpression gates

Read the conditional thresholds and diagnostics saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
with `stimControl(cytPosMethod = "coexpression")` and cytokine-positive
gating enabled.

## Usage

``` r
getStimGatesCoexpression(pathProject, pop = NULL)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- pop:

  character or NULL Populations to keep; NULL keeps all. Default: NULL.

## Value

A tibble with one row per stimulated sample and ordered channel pair:
`pop`, `batch`, `ind`; `chnlCond`, `markerCond` identify the
conditioning cytokine; `chnl`, `marker` identify the target cytokine;
`gate`, `cut`, `condCut` are its ordinary gate, lowered gate and raised
conditioning cut; `floor`, `floorCond` are the target and conditioning
floors; `z` is the net double-positive z score; `purityDp` is
double-positive purity (NA when the initial tests fail); and `lowered`
indicates whether `cut < gate`. Projects using refinement or with
cytokine-positive gating disabled error.

## Details

A cell is positive for `chnl` if its expression exceeds `gate`, or if
any lowered pair for that channel has conditioning-channel expression
strictly above `condCut` and target-channel expression strictly above
`cut`. Conditioning uses expression, never recursively inferred
positivity. Control expression is raw (without `biasUns`). The floor is
the control's main nonzero negative peak plus 1.5 robust standard
deviations, estimated with an nrd0 Gaussian density on 2048 points.
Lowering requires a net double-positive response, excess coexpression
beyond independence, and adequate purity; the added band is trimmed
along the conditioning channel. The four `coex*` controls in
[`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
set these tests and bin counts.

## Examples

``` r
# \donttest{
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- gateStim(tempfile("stimgate_"), gs, exampleData$batchList,
  marker = exampleData$marker,
  control = stimControl(cytPosMethod = "coexpression"))
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
getStimGatesCoexpression(pathProject)
#> # A tibble: 4 × 15
#>   pop   batch  ind   chnlCond  markerCond chnl  marker  gate   cut condCut floor
#>   <chr> <chr>  <chr> <chr>     <chr>      <chr> <chr>  <dbl> <dbl>   <dbl> <dbl>
#> 1 root  batch1 2     BC2(Pr14… MarkerF2   BC1(… Marke…  4.67  4.27    4.43 0.697
#> 2 root  batch1 2     BC1(La13… MarkerF1   BC2(… Marke…  4.09  2.92    4.67 0.738
#> 3 root  batch2 4     BC2(Pr14… MarkerF2   BC1(… Marke…  4.29  3.74    3.54 0.584
#> 4 root  batch2 4     BC1(La13… MarkerF1   BC2(… Marke…  3.40  3.40    4.29 0.625
#> # ℹ 4 more variables: floorCond <dbl>, z <dbl>, purityDp <dbl>, lowered <lgl>
# }
```
