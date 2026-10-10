# Get the safe lower boundary from the full stimulated marginal distribution

The left/main modal-complex peak is identified in the same way as in the
initial one-marker procedure. The negative width uses the closest left
half-height point or clear dip on that same marginal density (SJ
bandwidth with the `bwMin` floor). Refinement is not allowed at or below
`peakX + max(0.5 * windowWidth, densityBw)`.

## Usage

``` r
.getCytPosMarginalReference(ex, chnl, bwMin = NA_real_)
```
