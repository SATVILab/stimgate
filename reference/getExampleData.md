# Get example GatingSet

Load the canonical packaged example dataset shipped under
`inst/extdata/stimgate_example_data/`. This keeps the regular package
examples and tests on a deterministic, realistic dataset without
requiring the simulation machinery to be installed with the package.

## Usage

``` r
getExampleData()
```

## Value

A list with the saved example-data path, channel labels, marker labels,
and sample-to-condition mapping.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
gs
#> A GatingSet with 4 samples
exampleData$batchList
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 3 4
#> 
exampleData$marker
#> [1] "MarkerF1" "MarkerF2"
```
