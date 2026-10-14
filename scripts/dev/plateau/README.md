# Plateau-step and responder-gate experiments

Developer scripts behind #611 (`locPlateauStep`), #612 (dominance onset) and
#613 (25th-percentile responder gate). They are exploratory, not part of any
analysis, and were run inside the stimgate container from the
`stimgate-plateau` worktree root. Each reads a scratch folder from the
environment variable `P`.

## Re-running scenarios

- `rerun.R`: re-runs one Analysis 7 scenario (all 20 tubes of dataset 1) with
  `.simCompareRunScenario()`, wrapped in `.simDebugLoc()`, and saves
  `rr_<tag>.rds`: per tube the own (pre-sharing) gate, direct flag, reason,
  preliminary lower bound (`minProbXPos`), bandwidth and simulated cells, plus
  the scenario's final results. It needs `q7setup_wt.R`, the setup chunks of
  `analysis/7-sim-compare-freq_bs.qmd` purled up to the
  `Total simulation scenarios` line. Arguments: response probability, cell
  count, tag and an optional R list of call overrides (for example
  `"list(locPlateauStep = FALSE)"`); `TRANS` and `SEP` select transformation
  and separation.
- `queue.sh`: runs the 14 Gaussian scenarios used for #613 with the step off,
  at most four R processes at once.

## Comparing rules offline

- `rules.R`: on the `rr_off_*` files, compares the final gates with the step
  off, the plateau step (from the final Analysis 7 results, iteration 1),
  moving to the higher of the cluster's 25th-percentile responder gate and the
  taut-string edge, the 25th-percentile rule when the edge lies below it (all
  tubes, responders only, with at least 5 responders), the step combined with
  that rule, and F-beta constrained to the interval between the edge and the
  gate (the binned comparator form and a cumulative-count form). It checks
  that the re-runs reproduce the saved pre-step final gates.
- `why.R`: why the step barely moved the lowest-frequency scenarios.

## Comparing full analysis results

`res7.R`, `res78b.R`, `resreal.R` and `resomip.R` compare Analyses 7, 7b/8,
ACS (9b) and OMIP-111 before and after the plateau step. They were written
when the plateau results were the current results. Since 2026-10-14 the
current results folders hold the pre-plateau (0.99.34-equivalent) results
again, and the plateau results are in
`stimgate-store/_backup/2026-10-14-plateau-step-final/`; swap the paths
before re-running them.

## ACS checks from saved gates

- `acsedge.R`: background-subtracted frequency at the taut-string edge plus a
  bandwidth for every ACS tube.
- `acsq25.R`: moving every tube to the higher of the 25th-percentile responder
  gate (population × cytokine as the cluster) and the edge. This is the loose
  version; the responder-only rule of #613 has not been run on ACS.

## History (`history/`)

The earlier steps that led to the plateau step and #613, in rough order. They
read `dbg_*`/`res_*` files written by `../upstream-region/rerun-legacy.R`,
whose `PLATEAU=1` switch prototyped the step in the antimode slot before it
went into the package (`am_*` tags = current rule, `pb_*` = step before the
cap; the after-cap version was applied offline).

- `am.R`, `am2.R`, `am3.R`: why the taut-string antimode step never fires
  (the fit starts just right of the negative peak and flattens the response).
- `tsplot.R`, `tsf1.R`: plots of the taut string fitted from the flattening
  point, with true-positive rugs; F1 at the antimode and the flattening point.
- `plateau.R`, `bwstep.R`: moving to the plateau edge, and bandwidth-sized
  steps towards it, with purity checks.
- `cmpab.R`: the step before and after the cap.
- `bias.R`, `ratio.R`, `base.R`, `shiftu.R`, `ls.R`: batch-shift handling
  (edge margins, the stim/unstim tail ratio, shifting or rescaling unstim).
- `gate_r0.R`, `grouprule.R`: switching the step off by the tail ratio, per
  tube and pooled per group with between-tube spread; `realratio.R` applies
  the ratio to ACS and OMIP-111 tubes.
- `guards.R`, `guards2.R`, `half.R`: small-count guards (cautious purity,
  minimum cells, empty-gap rule), graded unstim shifts and the automatic
  half-step.
- `diag.R`, `edgecmp.R`, `direct.R`, `pk2.R`, `smoke.R`, `smoke2.R`: checks of
  the package implementation (example data, saved tubes gated directly,
  package results against the offline rule).
- `share.R`: own against shared gates (needs records that are no longer kept).
- `batchq.sh`, `skewq.sh`: the queued batch-shift and skew re-runs.
