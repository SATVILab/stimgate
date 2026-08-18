# Native StimGate right-tail density shoulder gate

This helper defines the native StimGate tailgate rule: after the main
peak, locate the steepest negative derivative on the descending shoulder
and report the first point at which the density has flattened to the
requested fraction of that slope. This relative-derivative rule
preserves the local-FDR shape threshold behaviour and deliberately does
not use a fixed absolute tolerance like the
cytoUtils:::.cytokine_cutpoint() `tol` gate.

## Usage

``` r
.getStimGateTailgate(density, peakX = NULL, fraction = 1/200)
```
