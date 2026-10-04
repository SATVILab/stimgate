# Get settings for a named marker

Retrieve the settings for a marker by its original name/key.

## Usage

``` r
stimgateMetaReadSettingsMarker(pathProject, marker)
```

## Arguments

- pathProject:

  character Path to project.

- marker:

  character Marker name/key as stored in markerList.

## Value

A list of settings for the requested marker.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  list(BC1 = list(bw = 0.1)),
  file.path(pathProject, "metaData", "chnlSettings.rds")
)
stimgateMetaReadSettingsMarker(pathProject, "BC1")
#> $bw
#> [1] 0.1
#> 
```
