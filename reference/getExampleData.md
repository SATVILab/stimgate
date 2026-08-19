# Get example GatingSet

Load the canonical packaged example dataset shipped under
`inst/extdata/stimgate_example_data/`. This keeps the regular package
examples and tests on a deterministic, realistic dataset without
requiring the simulation machinery to be installed with the package.

## Usage

``` r
getExampleData()

getTestData()
```

## Value

A list with the saved example-data path, channel labels, marker labels,
and sample-to-condition mapping.
