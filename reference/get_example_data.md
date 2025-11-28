# Get example GatingSet

Create and save a complete example GatingSet for testing and examples.
This function internally creates a flowSet, samples channels, and saves
the GatingSet.

## Usage

``` r
get_example_data(dir_cache = NULL)
```

## Arguments

- dir_cache:

  Directory to save the GatingSet. If NULL, uses a temporary directory.

## Value

A list containing the path to the saved GatingSet, batch_list, and
marker names
