# Generate a batch list of sample indices

Groups sample rows by batch/donor identifiers, screens out samples
falling below a minimum cell count threshold, and structures the output
so that the unstimulated control index is always positioned as the first
element of each batch.

## Usage

``` r
getBatchList(fnTblInfo, colGrp, colStim, unsChr, colNCell, minCell)
```

## Arguments

- fnTblInfo:

  data.frame. Sample metadata containing annotations.

- colGrp:

  character vector. One or more column names used to define
  batches/groups.

- colStim:

  character. Column name containing stimulation identifiers.

- unsChr:

  character. Name/string specifying the unstimulated control sample.

- colNCell:

  character. Column name containing the cell count for each sample.

- minCell:

  numeric. Minimum number of cells required to retain a sample.

## Value

A named list where each element contains a numeric vector of sample
indices representing a batch, with the unstimulated control index at the
beginning.

## Examples

``` r
fnTblInfo <- data.frame(
  donor = c("d1", "d1", "d2", "d2", "d3", "d3"),
  stim = c("stim", "uns", "uns", "stim", "uns", "stim"),
  nCell = c(5000, 4000, 6000, 5500, 3000, 50)
)
# Donor d3 is dropped because its stimulated sample has too few cells
getBatchList(
  fnTblInfo,
  colGrp = "donor",
  colStim = "stim",
  unsChr = "uns",
  colNCell = "nCell",
  minCell = 100
)
#> $d1
#> [1] 2 1
#> 
#> $d2
#> [1] 3 4
#> 
```
