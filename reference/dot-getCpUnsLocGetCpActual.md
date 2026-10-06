# Select the local-FDR threshold and place the gate below it

The selected threshold is the candidate cell value whose tail frequency
(cells at or above it) best matches the estimate. Gates are applied with
a strict `x > gate`, so the gate is placed below that cell, within the
gap to the next lower stimulated or unstimulated cell: it moves down by
the smaller of twice the density bandwidth and half that gap. The
applied gate then counts exactly the cells the selection counted.

## Usage

``` r
.getCpUnsLocGetCpActual(
  dataThreshold,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL,
  densityBw = NULL
)
```

## Arguments

- exTblStimOrig, exTblUnsOrig:

  data.frame or NULL Expression used to find the gap below the selected
  cell; NULL keeps the gate at the cell.

- densityBw:

  numeric, list or NULL Local-FDR density bandwidth; an adaptive
  bandwidth uses the shared bandwidth at the selected cell.

## Value

list from `.getCpUnsLocConditionOut()`, with the selected cell value as
attribute `cpSelected`.
