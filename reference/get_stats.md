# Get gating statistics

Get gating statistics

## Usage

``` r
get_stats(path_project)
```

## Arguments

- path_project:

  character. Path to the project directory.

## Value

A data frame with gating statistics.

## Examples

``` r
{
# Get example dataset
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)

# Run the stimgate pipeline
path_project <- stimgate_gate(
  path_project = file.path(tempdir(), "stimgate_example"),
  .data = gs,
  batch_list = example_data$batch_list,
  marker = example_data$marker,
  pop_gate = "root"
)

# Get statistics for the identified gates
stats <- get_stats(path_project)
}
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
#> 
#> see ?HDCytoData and browseVignettes('HDCytoData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> Done
#> To reload it, use 'load_gs' function
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
```
