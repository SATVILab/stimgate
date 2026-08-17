# Detect whether a transformation should receive upper-population ratio correction

This is currently applied to the raw gamma transformation and to the
skew transformation. The lower component is never rescaled by this
correction; only upper components are moved and rescaled so that their
observed-scale separation from the pre-perturbation lower reference
matches the corresponding raw-scale separation after sample, condition,
and cluster perturbations.

## Usage

``` r
.simCytUsesUpperRatioCorrection(transformationFunc)
```
