# Read settings by saved key

Extract one entry from
[`stimgateMetaReadSettingsChnls()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnls.md).
The key must match exactly; channel names are not converted to marker
labels.

## Usage

``` r
stimgateMetaReadSettingsChnl(pathProject, chnl)
```

## Arguments

- pathProject:

  character Project directory.

- chnl:

  character Exact saved key, usually a marker label.

## Value

A list of settings for the key; an unknown key raises an error.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(list(IFNg = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsChnl(pathProject, "IFNg")
#> $bw
#> [1] 0.1
#> 
```
