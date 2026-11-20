# RQ4E results: intercept-addressable error diagnostic

Generated: 2026-11-17T16:29:29.503278+08:00  
Status: post-hoc derivation of frozen RQ3 held-out-target results; no model was refitted.

## Target-equal macro summary

| Analysis set | Targets | Brier U | Brier I | Brier F | I−U | F−I | Relative reduction I vs U | Log loss I−U | Log loss F−I |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary_six | 6 | 0.308627 | 0.169348 | 0.161723 | -0.139279 | -0.007625 | 45.1% | -0.353319 | -0.027427 |
| primary_five | 5 | 0.314002 | 0.179272 | 0.180676 | -0.134730 | 0.001404 | 42.9% | -0.307330 | 0.004047 |

## primary_six: target-level Brier diagnostic

| Target | n | Positives | U | I | F | I−U | F−I | Relative reduction I vs U |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AAH-I | 133 | 19 | 0.265536 | 0.111674 | 0.068714 | -0.153861 | -0.042961 | 57.9% |
| VEGF110 | 268 | 158 | 0.322558 | 0.236011 | 0.236044 | -0.086547 | 0.000033 | 26.8% |
| Nipah G | 1028 | 101 | 0.300146 | 0.087752 | 0.087732 | -0.212394 | -0.000020 | 70.8% |
| PD-L1 | 277 | 198 | 0.226346 | 0.211639 | 0.205524 | -0.014707 | -0.006115 | 6.5% |
| SARS-CoV-2 RBD | 268 | 226 | 0.489905 | 0.138527 | 0.131913 | -0.351378 | -0.006614 | 71.7% |
| TREM2 | 110 | 46 | 0.247269 | 0.230484 | 0.240408 | -0.016786 | 0.009924 | 6.8% |

## primary_five: target-level Brier diagnostic

| Target | n | Positives | U | I | F | I−U | F−I | Relative reduction I vs U |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| VEGF110 | 268 | 158 | 0.267236 | 0.236373 | 0.236371 | -0.030864 | -0.000002 | 11.5% |
| Nipah G | 1028 | 101 | 0.410135 | 0.087699 | 0.087740 | -0.322437 | 0.000042 | 78.6% |
| PD-L1 | 277 | 198 | 0.216353 | 0.206520 | 0.205469 | -0.009834 | -0.001051 | 4.5% |
| SARS-CoV-2 RBD | 268 | 226 | 0.402299 | 0.131453 | 0.132744 | -0.270845 | 0.001291 | 67.3% |
| TREM2 | 110 | 46 | 0.273987 | 0.234315 | 0.241054 | -0.039673 | 0.006739 | 14.5% |

## Interpretation boundary

In the six-target analysis, intercept updating reduced target-equal macro Brier by 0.139279 (45.1% relative to U), with improvement in 6/6 targets. Further slope updating changed macro Brier by -0.007625, improving 4 and worsening 2 targets.

In the separately refitted five-target sensitivity analysis, the corresponding intercept-updating change was -0.134730 (42.9% relative reduction), with improvement in 5/5 targets. Further slope updating changed macro Brier by 0.001404, improving 2 and worsening 3 targets.

These results quantify an intercept-addressable error component and the incremental benefit of slope updating. They are not a causal decomposition and do not establish that base rate is the largest mechanistic cause of transport failure. Post-I Brier and BSS describe remaining error after intercept adaptation; they do not identify its biological or experimental cause. No new confidence intervals, bootstrap tests, or pilot recommendation are introduced by this diagnostic.

## Quality control

Checks passed: 92/92. Source SHA-256 values were unchanged.
