# Get example GatingSet

Create and save a complete synthetic GatingSet for testing and examples
using the simCytExperiment functions, avoiding external data
dependencies.

## Usage

``` r
getExampleData(
  scenario = "default",
  dirCache = NULL,
  clear = FALSE,
  nCell = 10000,
  nInd = 8
)

getTestData(
  scenario = "default",
  dirCache = NULL,
  clear = FALSE,
  nCell = 10000,
  nInd = 8
)
```

## Arguments

- scenario:

  Character specifying the simulation scenario ("default", "easy",
  "poorSeparation", "cytPos").

- dirCache:

  Directory to save the GatingSet. If NULL, uses a temporary directory.

- clear:

  Logical indicating whether to force cache clearing.

- nInd:

  Integer. Number of biological samples to simulate.

## Value

A list containing the path to the saved GatingSet, batchList, chnl, and
marker names.
