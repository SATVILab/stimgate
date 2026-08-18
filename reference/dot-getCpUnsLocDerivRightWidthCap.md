# Cap a right-side derivative threshold using the selected peak's left width

The left width is measured at half height. The multiplier is chosen so
that a Gaussian-shaped peak reaches `rightFrac` at the same standardised
distance.

## Usage

``` r
.getCpUnsLocDerivRightWidthCap(peakData, iPeak, rightFrac, leftFrac = 0.5)
```
