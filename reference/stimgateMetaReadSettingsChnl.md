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
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
stimgateMetaReadSettingsChnl(tmp, "BC1 label")
stimgateMetaReadSettingsChnl(tmp, "BC1")
} # }
```
