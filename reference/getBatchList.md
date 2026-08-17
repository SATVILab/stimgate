# Generate a batch list of sample indices

Groups sample rows by batch/donor identifiers, screens out samples
falling below a minimum cell count threshold, and structures the output
so that the unstimulated control index is always positioned as the final
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
end.
