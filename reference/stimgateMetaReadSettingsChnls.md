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
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  list(BC1 = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsChnls(pathProject)
#> $BC1
#> $BC1$bw
#> [1] 0.1
#> 
#> 
```
