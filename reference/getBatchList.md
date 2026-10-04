# Group samples into batches

Group metadata rows by donor or batch, keeping the unstimulated control
first. Drop samples below `minCell` and groups with fewer than two
remaining samples or no control. A group with more than one control is
an error.

## Usage

``` r
getBatchList(fnTblInfo, colGrp, colStim, unsChr, colNCell, minCell)
```

## Arguments

- fnTblInfo:

  data.frame Sample metadata, one row per sample.

- colGrp:

  character vector Column names defining batches.

- colStim:

  character Column containing stimulation labels.

- unsChr:

  character Label identifying unstimulated controls.

- colNCell:

  character Column containing sample cell counts.

- minCell:

  numeric Minimum cell count to retain a sample.

## Value

A named list of integer row indices per batch. Names join group values
with underscores; control indices precede stimulated indices.

## Examples

``` r
samples <- data.frame(
  donor = c("d1", "d1", "d2", "d2"),
  stim = c("stim", "uns", "uns", "stim"),
  nCell = c(5000, 4000, 3000, 50)
)
# d2 is dropped: only its control meets the cell-count limit
getBatchList(samples, "donor", "stim", "uns", "nCell", minCell = 100)
#> $d1
#> [1] 2 1
#> 
```
