# Get gates

Get all the gates for each of the markers gated.

## Usage

``` r
getStimGates(pathProject, pop = NULL, marker = NULL, chnl = NULL)
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

## Examples

``` r
# Get example dataset
exampleData <- getExampleData()
#> Cache incomplete, regenerating synthetic test data...
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)

# Run the stimgate pipeline
pathProject <- gateStim(
  pathProject = file.path(tempdir(), "getGateExample"),
  .data = gs,
  batchList = exampleData$batchList,
  marker = exampleData$marker,
  popGate = "root"
)
#> ----
#> getting base gates
#> ----
#> 
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 8 of 8
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 8 of 8
#> getting clustered and/or controlled gates
#> 
#> 
#> ----
#> getting single+ gates
#> ----
#> 
#> 
#> 
#> 
#> getting cyt combn frequencies
#> batch 8 of 8

# Get identified gates
gates <- getStimGates(pathProject)
#> Warning: cannot open compressed file '/tmp/RtmpZlgykD/getGateExample/gates/poproot/chnlBC2(Pr141)Dd/all/gateTbl.rds', probable reason 'No such file or directory'
#> Error in map(.x, .f, ...): ℹ In index: 2.
#> Caused by error in `map()`:
#> ℹ In index: 2.
#> Caused by error in `gzfile()`:
#> ! cannot open the connection
```
