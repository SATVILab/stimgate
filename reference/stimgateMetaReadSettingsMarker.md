# Read settings for a marker

Extract one entry from
[`stimgateMetaReadSettingsChnls()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnls.md).

## Usage

``` r
stimgateMetaReadSettingsMarker(pathProject, marker)
```

## Arguments

- pathProject:

  character Project directory.

- marker:

  character Exact saved marker key.

## Value

A list of marker settings; an unknown key raises an error.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(list(IFNg = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsMarker(pathProject, "IFNg")
#> $bw
#> [1] 0.1
#> 
```
