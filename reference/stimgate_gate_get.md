# Get gates

Get all the gates for each of the markers gated.

## Usage

``` r
stimgate_gate_get(path_project)
```

## Arguments

- path_project:

  character. Path to the project directory.

## Value

Gate table with gates for each sample for each marker.

## Examples

``` r
{
# Get example dataset
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)

# Run the stimgate pipeline
path_project <- stimgate_gate(
  path_project = file.path(tempdir(), "get_gate_example"),
  .data = gs,
  batch_list = example_data$batch_list,
  marker = example_data$marker,
  pop_gate = "root"
)

# Get statistics for the identified gates
gates <- stimgate_gate_get(path_project)
}
#> see ?HDCytoData and browseVignettes('HDCytoData') for documentation
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
