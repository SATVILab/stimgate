# Read marker list with channel labels

Read the project's marker list and return it with element names replaced
by channel labels (from chnlLab).

## Usage

``` r
stimgateMetaReadSettingsMarkers(pathProject)
```

## Arguments

- pathProject:

  character Path to project.

## Value

A named list of marker settings where names are channel labels.

## Examples

``` r
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
stimgateMetaReadSettingsChnls(tmp)
} # }
```
