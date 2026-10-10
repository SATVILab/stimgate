# Settings for the optional shifted-peak rule (`locShiftedPeakRef`), or NULL when it is off. The reference bandwidth is the per-channel shared bandwidth at the `bwNcellMax` reference size, before any widening for small tubes (smaller of the stimulated and unstimulated tubes' values with `bwScope = "cluster"`). Without a shared bandwidth it is the bandwidth of the local-FDR densities: a fixed `bw`, the per-sample estimate with `bwScope = "sample"`, or for adaptive bandwidths the blended bandwidth curve, read at the unstimulated peak. This trigger bandwidth is separate from the actual density bandwidth used for negative widths and search offsets.

Settings for the optional shifted-peak rule (`locShiftedPeakRef`), or
NULL when it is off. The reference bandwidth is the per-channel shared
bandwidth at the `bwNcellMax` reference size, before any widening for
small tubes (smaller of the stimulated and unstimulated tubes' values
with `bwScope = "cluster"`). Without a shared bandwidth it is the
bandwidth of the local-FDR densities: a fixed `bw`, the per-sample
estimate with `bwScope = "sample"`, or for adaptive bandwidths the
blended bandwidth curve, read at the unstimulated peak. This trigger
bandwidth is separate from the actual density bandwidth used for
negative widths and search offsets.

## Usage

``` r
.getCpUnsLocShiftedPeakSettings(
  chnlSettings,
  exTblStimThreshold,
  exTblUnsThreshold,
  densityBw
)
```
