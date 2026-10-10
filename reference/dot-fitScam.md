# Fit a monotone increasing SCAM to the modelled probabilities

Returns NULL when there are too few points and a try-error when the fit
fails. quiet suppresses warnings raised while fitting. The expression
values are fitted under the fixed column name `x`, because channel names
such as `PE-A` are not valid in a model formula.

## Usage

``` r
.fitScam(dataMod, bs, family, quiet)
```
