# Export stimulation-positive cells as FCS files

Select positive cells using saved or supplied gates and write one FCS
file per sample that has positive cells, plus `manifest.csv`.

## Usage

``` r
writeStimFCS(
  pathProject,
  .data,
  pop = NULL,
  indBatchList,
  pathDirSave,
  chnl = NULL,
  gateTbl = NULL,
  transFn = NULL,
  transChnl = NULL,
  combnExc = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE,
  gateUnsMethod = "min"
)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- .data:

  GatingSet, flowSet, cytoset, flowFrame, cytoframe, character, list or
  data.frame Data passed to
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
  in the same sample order.

- pop:

  character or NULL Population to export. NULL uses the single saved
  population, or "root" when `gateTbl` is supplied. Default: NULL.

- indBatchList:

  list Sample indices or names grouped by batch, unstimulated sample
  first, as for `batchList` in
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

- pathDirSave:

  character Output directory; existing contents are deleted.

- chnl:

  character vector or NULL Channels used to select positive cells; NULL
  uses all channels in the gate table. Default: NULL.

- gateTbl:

  data.frame or NULL Gates with `chnl`, `batch`, `ind`, `gate` and, for
  refined gates, `gateCyt`. NULL reads saved gates. Default: NULL.

- transFn:

  function or NULL Transformation applied to the selected cells before
  export. Default: NULL.

- transChnl:

  character vector or NULL Columns to transform; NULL transforms all
  expression columns. Default: NULL.

- combnExc:

  list or NULL Channel combinations to exclude: each vector specifies
  positive channels, with other selected channels negative. Default:
  NULL.

- gateTypeCytPos:

  character Positivity rule: "base" uses main gates; "cyt" uses the
  saved refinement or coexpression rule. Default: "cyt".

- mult:

  logical Require positivity for at least two markers. Default: FALSE.

- gateUnsMethod:

  character Summary of stimulated gates used for missing control gates:
  "min", "max", "mean", "tmean" (20% trimmed mean), or "med". Default:
  "min".

## Value

Invisibly, a tibble with one row per sample and columns `ind`, `batch`,
`fileName`, `nCellPos`, `written`, `reason`. The `pathDirSave` attribute
holds the output path. Samples with no positive cells get no FCS file.

## Details

With coexpression, `coexpression.csv` records pairwise thresholds.
Stimulated samples use their saved rules. For control exports, the
chosen `gateUnsMethod` summarises ordinary, lowered and conditioning
thresholds across the batch's stimulated samples.

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
#> shared bandwidth for MarkerF1: 0.188
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
manifest <- writeStimFCS(
  pathProject, gs, indBatchList = exampleData$batchList,
  pathDirSave = tempfile("positive_fcs_")
)
#> Using population 'root' from project directory.
#> Using channels 'BC1(La139)Dd, BC2(Pr141)Dd' from project directory.
#> Writing 1 of 4 files
#> Wrote sample001_unstim.fcs
#> Writing 2 of 4 files
#> Wrote sample001_stim1.fcs
#> Writing 3 of 4 files
#> Wrote sample002_unstim.fcs
#> Writing 4 of 4 files
#> Wrote sample002_stim1.fcs
```
