# Identify a smoothed left-to-right rise in density dominance

Stimulated and unstimulated densities are smoothed by Gaussian kernel
regression using half the fixed density bandwidth. When a derivative
threshold exists, densities to its right are first collapsed to their
respective means. A dominance threshold requires an observed transition
from non-dominance to dominance. Its location is two-thirds of the onset
location plus one-third of the subsequent dominance-score peak location.

## Usage

``` r
.getCpUnsLocMarginalDominanceStart(
  density,
  startX = NA_real_,
  densityBw = NULL,
  dominanceRatio = 2,
  onsetWeight = 2/3,
  lowerBoundX = NA_real_
)
```
