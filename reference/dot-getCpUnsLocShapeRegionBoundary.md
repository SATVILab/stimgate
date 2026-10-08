# Lower boundary of the region kept by the shape-enforced filters

The shape route filters in turn at the pre-fit shape threshold, the
global derivative threshold and the final marginal cut, each keeping
`x >= cut`. The kept region therefore starts at the largest of the cuts
that were applied. If none was applied, it starts at the lowest kept
value.

## Usage

``` r
.getCpUnsLocShapeRegionBoundary(
  dataMod,
  shapeLowerBoundX,
  globalInfo,
  marginalInfo
)
```

## Value

list with `xSum` (the boundary), `source` and the candidate cuts.
