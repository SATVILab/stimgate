# Discount disagreement at the start of the responding region

Within half a density bandwidth of the lowest estimate value, the
initial contiguous run with `pred > probSmooth` receives weights
increasing linearly from zero at `probSmooth / pred = 0.75` to one at
ratio one. Only zero-weight rows in that window are removed as candidate
gates. The existing exclusion of minimum-expression rows (except a
single row) is retained.

## Usage

``` r
.getCpUnsLocGetCpDataThresholdCount(dataMod, densityBw = NULL)
```

## Arguments

- dataMod:

  data.frame Estimate rows after lower-margin exclusion.

- densityBw:

  numeric, list or NULL Density bandwidth. Default: NULL.

## Value

data.frame Candidate rows with a response-estimate `weight` column.
