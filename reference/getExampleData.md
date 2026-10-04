# Load example cytometry data

Load the packaged dataset and save a GatingSet in a temporary directory.
Use the returned paths and labels to try
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## Usage

``` r
getExampleData()
```

## Value

A list with `pathGs` (saved GatingSet path), `batchList` (sample indices
by batch, control first), `chnl` (channels) and `marker` (labels).

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
exampleData$batchList
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3 4
#> 
```
