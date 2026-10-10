# Get the density grids used for both main peaks and negative widths

Fixed-bandwidth raw densities span the stimulated cells. If unstimulated
cells extend beyond that grid, recompute both Gaussian kernel densities
over the joint range at the same local-FDR bandwidth. Adaptive densities
already span both tubes and are reused without changing their bandwidth
curve.

## Usage

``` r
.getCpUnsLocGetPeakDensities(densTblRaw, exVecStim, exVecUns)
```
