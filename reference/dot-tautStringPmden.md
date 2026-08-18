# Piecewise-constant taut-string density via native C++ implementation

Returns a list with element `$y`: a numeric vector of length `n - 1`
(one value per interval between consecutive sorted data points)
representing the piecewise-constant taut-string density from the
FAUST-derived `cpPmden()` algorithm. Mode regions have higher values
than antimode regions, enabling downstream antimode detection.

## Usage

``` r
.tautStringPmden(x_sorted)
```
