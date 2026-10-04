# Get channel-to-marker labels

Map channel names to marker labels in a cytometry object. Channels
without a marker label use their channel name.

## Usage

``` r
chnlLab(data)
```

## Arguments

- data:

  flowFrame, flowSet, GatingSet, GatingHierarchy, cytoframe or cytoset
  Cytometry data. For sets, labels come from the first sample.

## Value

A character vector of marker labels, named by channel.

## Examples

``` r
exampleData <- getExampleData()
#> Done
#> To reload it, use 'load_gs' function
gs <- flowWorkspace::load_gs(exampleData$pathGs)
chnlLab(gs)
#> BC1(La139)Dd BC2(Pr141)Dd 
#>   "MarkerF1"   "MarkerF2" 
```
