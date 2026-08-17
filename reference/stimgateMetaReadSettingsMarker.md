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
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
stimgateMetaReadSettingsMarker(tmp, "BC1")
} # }
```
