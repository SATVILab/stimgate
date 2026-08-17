# Read channel or marker label mapping

Read the saved channel label mapping (chnlLab.rds) from the project's
metaData folder.

## Usage

``` r
stimgateMetaReadChnlLab(pathProject)

stimgateMetaReadMarkerLab(pathProject)
```

## Arguments

- pathProject:

  character Path to project.

## Value

Named character vector mapping channel names to labels.

## Examples

``` r
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
stimgateMetaReadChnlLab(tmp)
} # }
```
