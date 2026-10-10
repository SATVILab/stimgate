# Package index

## Gating Functions

Main functions for identifying cytokine-positive cells

- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  : Gate cells responding to stimulation
- [`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
  : Tune stimulation gating
- [`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md)
  : Read stimulation gates
- [`getStimGatesCoexpression()`](https://satvilab.github.io/stimgate/reference/getStimGatesCoexpression.md)
  : Read pairwise cytokine coexpression gates
- [`getStimGatesDetailed()`](https://satvilab.github.io/stimgate/reference/getStimGatesDetailed.md)
  : Read details of how gates were chosen
- [`getStimStats()`](https://satvilab.github.io/stimgate/reference/getStimStats.md)
  : Read gating statistics

## Visualization

Functions for plotting gates and results

- [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md)
  : Plot stimulation gates

## Data Export

Functions for exporting gated data

- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  : Export stimulation-positive cells as FCS files

## Utilities

Helper and utility functions

- [`chnlLab()`](https://satvilab.github.io/stimgate/reference/chnlLab.md)
  : Get channel-to-marker labels
- [`getBatchList()`](https://satvilab.github.io/stimgate/reference/getBatchList.md)
  : Group samples into batches
- [`getExampleData()`](https://satvilab.github.io/stimgate/reference/getExampleData.md)
  : Load example cytometry data
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  : Read cell expression values
- [`stimgateMetaReadBatchList()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadBatchList.md)
  : Read saved batches
- [`stimgateMetaReadChnlLab()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadLab.md)
  [`stimgateMetaReadMarkerLab()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadLab.md)
  : Read channel and marker mappings
- [`stimgateMetaReadSettingsChnl()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnl.md)
  : Read settings by saved key
- [`stimgateMetaReadSettingsChnls()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsChnls.md)
  : Read saved gating settings
- [`stimgateMetaReadSettingsMarker()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsMarker.md)
  : Read settings for a marker
- [`stimgateMetaReadSettingsMarkers()`](https://satvilab.github.io/stimgate/reference/stimgateMetaReadSettingsMarkers.md)
  : Relabel saved settings
