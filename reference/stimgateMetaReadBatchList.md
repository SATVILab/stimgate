# Read saved batches

Read the sample grouping saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Usage

``` r
stimgateMetaReadBatchList(pathProject)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Value

A list of sample indices by batch, with unstimulated controls first.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(list(batch1 = c(1, 2)),
  file.path(pathProject, "metaData", "batchList.rds")
)
stimgateMetaReadBatchList(pathProject)
#> $batch1
#> [1] 1 2
#> 
```
