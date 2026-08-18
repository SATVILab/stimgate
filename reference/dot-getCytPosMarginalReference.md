# Get the safe lower boundary from the full stimulated marginal distribution

The left/main modal-complex peak is identified in the same way as in the
initial one-marker procedure. `windowWidth` is the span from the 5th
percentile of values below that peak to the peak itself.
Cytokine-positive refinement is not allowed at or below
`peakX + windowWidth / 3`.

## Usage

``` r
.getCytPosMarginalReference(ex, chnl, bwMin = NA_real_)
```
