# Place a gate in the gap below a selected cell value

Bandwidth resolution is shared with the leading-run response discount.

## Usage

``` r
.getCpUnsLocGateBelowCell(cp, x, densityBw = NULL)
```

## Arguments

- cp:

  numeric Selected cell value.

- x:

  numeric Stimulated and unstimulated expression; NULL keeps `cp`.

- densityBw:

  numeric, list or NULL Density bandwidth. An adaptive bandwidth object
  is interpolated at `cp`; an unavailable bandwidth leaves only the
  half-gap limit.

## Value

numeric Gate strictly below `cp` and above every lower value of `x`.
