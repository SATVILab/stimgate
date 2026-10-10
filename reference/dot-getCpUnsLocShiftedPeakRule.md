# Apply the shifted-peak rule: it fires when the stimulated main peak lies more than `mult` reference bandwidths to the right of the unstimulated peak and the unstimulated negative width is available. NULL when it is off. Its width is measured by `.getCpUnsLocNegWidth()` on the peak density; the search offset uses half that width with a one-density-bandwidth minimum.

Apply the shifted-peak rule: it fires when the stimulated main peak lies
more than `mult` reference bandwidths to the right of the unstimulated
peak and the unstimulated negative width is available. NULL when it is
off. Its width is measured by
[`.getCpUnsLocNegWidth()`](https://satvilab.github.io/stimgate/reference/dot-getCpUnsLocNegWidth.md)
on the peak density; the search offset uses half that width with a
one-density-bandwidth minimum.

## Usage

``` r
.getCpUnsLocShiftedPeakRule(shiftedPeak, peakStimX, peakUnsX, windowWidthUns)
```
