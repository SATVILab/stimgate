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
