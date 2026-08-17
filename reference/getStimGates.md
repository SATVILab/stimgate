# Get gates

Get all the gates for each of the markers gated.

## Usage

``` r
getStimGates(pathProject)

getStimGates(pathProject)
```

## Arguments

- pathProject:

  character. Path to the project directory.

- pop:

  character. Optional population name(s) to filter gates by. Default is
  NULL (all populations).

- marker:

  character. Optional marker name(s) to filter gates by. Default is NULL
  (all markers).

- chnl:

  character. Optional channel name(s) to filter gates by. Default is
  NULL (all channels).

## Value

Gate table with gates for each sample for each marker.

A data frame with gating statistics.

## Examples

``` r
{
# Get example dataset
exampleData <- get_example_data()
gs <- flowWorkspace::load_gs(exampleData$pathGs)

# Run the stimgate pipeline
pathProject <- gateStim(
  pathProject = file.path(tempdir(), "getGateExample"),
  .data = gs,
  batchList = exampleData$batchList,
  marker = exampleData$marker,
  popGate = "root"
)

# Get statistics for the identified gates
gates <- getStimGates(pathProject)
}
#> Error in get_example_data(): could not find function "get_example_data"
```
