# Measure the left half-width of a tube's main negative peak

Scan left for the first half-height crossing (linearly interpolated) and
the first local minimum at least one bandwidth from the peak with height
at most 75% of the peak. Use the closer candidate; ties use half-height.
A flat local minimum is represented by its rightmost grid point. If
neither candidate exists, use the minimum finite observed expression.

## Usage

``` r
.getCpUnsLocNegWidth(density, peakX, exVec, bw)
```

## Arguments

- density:

  list Density grid with increasing `x` and corresponding `y`.

- peakX:

  numeric Main peak location on that density grid.

- exVec:

  numeric Tube expression values, used for the fallback only.

- bw:

  numeric Local-FDR bandwidth at this peak (shared curve if adaptive).

## Value

list Width, source (`half_height`, `dip`, or `fallback`), selected
boundary, both candidate locations, peak height and bandwidth.
