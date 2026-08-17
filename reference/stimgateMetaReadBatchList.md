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
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(list(batch1 = c(1, 2)), file.path(tmp, "metaData", "batchList.rds"))
stimgateMetaReadBatchList(tmp)
} # }
```
