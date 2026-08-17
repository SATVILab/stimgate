# Read marker settings from project

Read the saved marker settings list from the project's metaData folder.

## Usage

``` r
stimgateMetaReadSettingsChnls(pathProject)
```

## Arguments

- pathProject:

  character Path to project.

## Value

A named list of marker settings (as saved by
.completeChnlSettingsSave()).

## Examples

``` r
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
stimgateMetaReadSettingsChnls(tmp)
} # }
```
