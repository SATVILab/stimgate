# Return the right-most peak belonging to the left/main modal complex.

Peaks whose height is at least `peakMinRel * max(y)` are treated as
meaningful. If meaningful peaks are separated by a trough that is low
relative to both adjacent peaks and to the absolute peak, the first such
trough ends the left/main modal complex. Otherwise, shoulders and
unresolved peaks are allowed to belong to the same background complex.

## Usage

``` r
.getPeakMainLeftIdx(y, peakMinRel = 0.75, troughMaxRel = 0.75)
```
