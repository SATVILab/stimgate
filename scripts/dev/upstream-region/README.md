# Region-selection experiments (issue #612)

Exploratory scripts from October 2026 that tried to fix StimGate's region
selection before the cap at low separation. They were rebuilt from the session
transcript, since the original scratch copies were lost. They are not part of
any analysis. They ran inside the stimgate container from the main checkout
root, with the package at 0.99.36 (local-FDR code as in 0.99.34), and each
reads a scratch folder from the environment variable `P`.

## Background

At low separation the gate starts too high and rarely moves:

- The region starts at `dominanceX = 2/3 onset + 1/3 dominance-score peak`
  (`.getCpUnsLocDominanceBoundaryCurrent()`). The onset is the left end of the
  contiguous stretch where stim is at least twice unstim. At low separation it
  already lies right of the responders' density peak.
- The quality scan (`.getCpUnsLocQualityBoundaryCurrent()` /
  `.getCpUnsLocFilterMarginalBins()`) rarely extends below that start. Its
  crowding reference uses bins right of the start, on the responders' falling
  side, and its fitted purity is pulled down by the unstim bias.

See #612 for the results.

## The experiment harness

- `rerun-legacy.R`: re-runs one Analysis 7 scenario (all 20 tubes of dataset
  1) through `.simCompareRunScenario()`, wrapped in `.simDebugLoc()`. It saves
  `dbg_<tag>.rds` (per-tube summary, gate lines, simulated cells, the
  filter's `final` and `marginal` info, the antimode info and `minProbXPos`)
  and `res_<tag>.rds` (all methods' results).
  - Arguments: response probability, cell count, tag, an optional R list of
    call overrides (e.g. `"list(locMarginalPurityRel = 0.35,
    locThresholdMethod = 'region')"`), and an optional R list of row
    overrides (e.g. a batch shift: `"list(stim_sd_multiplier = 1.1,
    stim_sd_multiplier_clusters = 'gn', mismatch_type = 'sd_inflation')"`).
  - `TRANS` and `SEP` select the transformation and separation.
  - It sources `q7setup.R`: the setup chunks of
    `analysis/7-sim-compare-freq_bs.qmd` (commit 824202a8), purled with
    `knitr::purl(..., documentation = 0)` and cut after the
    `Total simulation scenarios` line.
- `mb.R`: a modified copy of `.getCpUnsLocFilterMarginalBins()`
  (`R/cp_uns_loc_filtering.R`, lines 1049–1333 at that time) defining
  `mbMod()`, with the variants below. `rerun-legacy.R` sources it when
  `QUAL_NEW=1`.

Environment switches read by `rerun-legacy.R`:

| Variable | Meaning |
|---|---|
| `QUAL_NEW=1` | Replace the quality boundary with the experimental one. The scan may go down to a floor: the lower of the dominance start and a Tailgate-style flattening point on one tube plus half a bandwidth. |
| `TG_TOL` | Flattening tolerance, as a share of the steepest slope (default 10^-2.5; 0.01 was also used). |
| `FLOOR_TUBE` | `stim` or `unstim`: the tube used for the flattening point. |
| `CROWD`, `CROWD_FRAC` | Crowding check: `ref` (three largest non-empty bins right of the start), `global` or `none`; cell limit as a multiple of their mean (default 2/3; 2 was also tried). |
| `PAIR=1` | Also test each bin together with the next one. |
| `PARTIAL=1` | Clip, not drop, the bin that straddles the lower bound (the partial-bin fix). |
| `MULTIRES=1` | Multi-scale scan: spans of 1/16, 1/8, 1/4 and 1/2 of the distance to the floor, plus the whole distance; take the smallest accepted span and repeat. |
| `PUR_SRC`, `PUR_EXP` | Purity from fitted probabilities (`fit`) or counts (`count`); the size-dependent bar is the purity bar × min(1, (span cells / reference cells)^`PUR_EXP`). |
| `DOM_ONSET_HALFBW=1` | Start at the dominance onset plus half a bandwidth instead of 2/3 onset + 1/3 peak. |
| `PLATEAU=1` (`PL_K`, `PL_BAR`, `PL_REJ`) | Prototype of the plateau step in the antimode slot. Superseded by `stimControl(locPlateauStep = )` (#611); see `../plateau/`. |

The best low-separation variant was:

```bash
P=<scratch> TRANS=gaussian SEP=low QUAL_NEW=1 MULTIRES=1 FLOOR_TUBE=stim \
  DOM_ONSET_HALFBW=1 CROWD_FRAC=2 TG_TOL=0.01 PUR_SRC=count PUR_EXP=0.5 CROWD=none \
  Rscript rerun-legacy.R 0.002 1e5 g002_nc_reg "list(locMarginalPurityRel = 0.35, locThresholdMethod = 'region')"
```

It gave median F1 0.50 → 0.65 at 0.2% response and 0.35 → 0.54 at 0.05%, but
at high separation and 0.01% it dropped from 1.000 to 0.960: the multi-scale
scan jumped straight to the floor on noisy small-count spans.

## Analysis scripts

All read the `dbg_*`/`res_*` files, taking tags as arguments unless noted:

- `cmp.R`: median gate, gate SD, F1 and frequency error by variant, with the
  region boundaries (`xClear`, `xDom`, Tailgate point, floor, `xQual`) and
  stop reasons.
- `scan.R`, `scan2.R`, `why.R`, `why2.R`, `stop.R`: the quality scan bin by
  bin, and why it stopped (crowding against purity).
- `dom.R`, `mkdom.R`: dominance onset, peak and start per tube.
- `fl.R`, `bw.R`: the floor against Tailgate's gate; flattening bandwidth
  against `ks::hpi()`.
- `pool.R`, `pool2.R`: own against shared gates.
- `over.R`: the spans accepted in over-calling tubes, against a baseline tag.
- `a.R`, `b.R`, `c.R`, `an.R`, `peek.R`: first looks at the cached
  Analysis 7 results.
- `tglow.R`, `tgdom.R`: where Tailgate cuts too low (Analyses 7, 7b, 8 and
  ACS), and how its flattening point compares with the dominance start.
