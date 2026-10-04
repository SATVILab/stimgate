# Read channel and marker mappings

Read channel-to-marker labels with `stimgateMetaReadChnlLab()`; read the
reverse mapping with `stimgateMetaReadMarkerLab()`.

## Usage

``` r
stimgateMetaReadChnlLab(pathProject)

stimgateMetaReadMarkerLab(pathProject)
```

## Arguments

- pathProject:

  character Project directory from
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Value

A named character vector: channel names to marker labels for
`stimgateMetaReadChnlLab()`, marker labels to channels for
`stimgateMetaReadMarkerLab()`.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(c(BC1 = "IFNg"), file.path(pathProject, "metaData", "chnlLab.rds"))
stimgateMetaReadChnlLab(pathProject)
#>    BC1 
#> "IFNg" 
stimgateMetaReadMarkerLab(pathProject)
#>  IFNg 
#> "BC1" 
```
