# Identify the current contiguous density-dominance boundary

The stimulated and unstimulated densities are first stabilised in the
same way as the earlier implementation: each density is collapsed to its
mean to the right of x_clear_init, then both are Gaussian-kernel
smoothed with half the original fixed density bandwidth. The 2:1
dominance region used here is specifically the contiguous region
containing x_clear_init. Its left onset is tempered by the strongest
dominance point between that onset and x_clear_init.

## Usage

``` r
.getCpUnsLocDominanceBoundaryCurrent(
  density,
  startX,
  densityBw = NULL,
  dominanceRatio = 2,
  onsetWeight = 2/3,
  lowerBoundX = NA_real_
)
```
