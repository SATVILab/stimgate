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
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  list(BC1 = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
saveRDS(
  c(BC1 = "IFNg"),
  file.path(pathProject, "metaData", "chnlLab.rds")
)
stimgateMetaReadSettingsMarkers(pathProject)
#> $IFNg
#> $IFNg$bw
#> [1] 0.1
#> 
#> 
```
