# Restrict response search using the main negative-population widths

Start above the higher main peak plus half the larger tube width, with a
minimum offset of one local-FDR bandwidth at that peak. The shifted-peak
rule instead uses only the unstimulated peak, width and density
bandwidth.

## Usage

``` r
.getCpUnsLocProbTblFilter(
  probTbl,
  exVecStim,
  exVecUns,
  stage,
  peakStimX,
  peakUnsX,
  shiftedPeak = NULL,
  densityStim,
  densityUns,
  densityBw
)
```

## Arguments

- densityStim, densityUns:

  list Density grids used to identify the peaks.

- densityBw:

  numeric or list Local-FDR bandwidth or adaptive shared curve.
