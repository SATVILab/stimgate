# Relabel saved settings

Read
[`stimgateMetaReadSettingsChnls()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnls.md)
and map its names through
[`stimgateMetaReadChnlLab()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadLab.md).
Use only when the saved keys are channel names.

## Usage

``` r
stimgateMetaReadSettingsMarkers(pathProject)
```

## Arguments

- pathProject:

  character Project directory.

## Value

A list with names replaced by marker labels; unmatched keys become NA.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(list(BC1 = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
saveRDS(c(BC1 = "IFNg"), file.path(pathProject, "metaData", "chnlLab.rds"))
stimgateMetaReadSettingsMarkers(pathProject)
#> $IFNg
#> $IFNg$bw
#> [1] 0.1
#> 
#> 
```
