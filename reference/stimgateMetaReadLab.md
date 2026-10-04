# Read channel or marker label mapping

Read the saved channel label mapping (chnlLab.rds) from the project's
metaData folder.

## Usage

``` r
stimgateMetaReadChnlLab(pathProject)

stimgateMetaReadMarkerLab(pathProject)
```

## Arguments

- pathProject:

  character Path to project.

## Value

Named character vector mapping channel names to labels.

## Examples

``` r
pathProject <- tempfile("stimgate_meta_")
dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
saveRDS(
  c(BC1 = "IFNg"),
  file.path(pathProject, "metaData", "chnlLab.rds")
)
stimgateMetaReadChnlLab(pathProject)
#>    BC1 
#> "IFNg" 
stimgateMetaReadMarkerLab(pathProject)
#>  IFNg 
#> "BC1" 
```
