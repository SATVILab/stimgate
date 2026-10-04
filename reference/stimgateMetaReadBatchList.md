# Read batch list from project

Read the saved batchList object from the project's metaData folder.

## Usage

``` r
stimgateMetaReadBatchList(pathProject)
```

## Arguments

- pathProject:

  character Path to project.

## Value

A list describing sample grouping into batches (as saved by
.saveMetaDataBatchList()).

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  list(batch1 = c(1, 2)),
  file.path(pathProject, "metaData", "batchList.rds")
)
stimgateMetaReadBatchList(pathProject)
#> $batch1
#> [1] 1 2
#> 
```
