# Resolve the local-FDR threshold method

A missing (NULL or NA) per-marker value inherits the global value; a
control object saved before the option existed has no global value and
uses "region".

## Usage

``` r
.completeChnlSettingsLocThresholdMethod(method, methodCommon)
```
