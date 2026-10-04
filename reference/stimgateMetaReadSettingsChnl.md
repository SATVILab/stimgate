# Get marker settings for a single channel

Retrieve the marker settings for a single channel. The function accepts
either a channel label (as returned by stimgateMetaReadChnlLab) or the
original channel name/key used in markerList.

## Usage

``` r
stimgateMetaReadSettingsChnl(pathProject, chnl)
```

## Arguments

- pathProject:

  character Path to project.

- chnl:

  character Channel label or channel name.

## Value

A list of settings for the requested channel.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  list(BC1 = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsChnl(pathProject, "BC1")
#> $bw
#> [1] 0.1
#> 
```
