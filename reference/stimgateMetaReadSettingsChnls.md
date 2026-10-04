# Read saved gating settings

Read settings saved by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Usage

``` r
stimgateMetaReadSettingsChnls(pathProject)
```

## Arguments

- pathProject:

  character Project directory.

## Value

A list of settings per marker, named by the saved keys (marker labels
for projects created by
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)).

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(list(IFNg = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsChnls(pathProject)
#> $IFNg
#> $IFNg$bw
#> [1] 0.1
#> 
#> 
```
