# C01-C02 results for manuscript integration

Date: 2026-11-21

## C01: ipSAE extension

### Feasibility

- Common complete-case set: 2,074 pairs, 740 primary positives, six targets.
- The same 10 PD-L1 pairs lack Pinc and both released ipSAE fields; ipSAE caused no additional row loss.
- PAE dimensions, finite values, chain fields, and residue-token sums passed for all 2,074 pairs.
- Official `DunbrackLab/IPSAE` implementation was pinned to commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`, PAE cutoff 10 A, and distance cutoff 10 A.
- Official recomputation succeeded for 2,074/2,074 pairs.
- The released Proteinbase `boltz2_min_ipsae` field reproduced the official directional minimum in only 1,048/2,074 pairs. It was therefore retained as provenance but not used as the analysis score.
- No Boltz-2 affinity metric was present in the frozen Proteinbase release, so an affinity-score comparison could not be performed.

### Target-equal results

| Score/model | Layer | AUC | AP | Brier | Log loss |
|---|---|---:|---:|---:|---:|
| ipTM | raw | 0.6770 | 0.6132 | 0.4154 | 1.2224 |
| ipSAE | raw | 0.6658 | 0.5903 | 0.3341 | 1.6792 |
| ipTM | U | 0.6770 | 0.6132 | 0.3139 | 0.9204 |
| ipSAE | U | 0.6658 | 0.5903 | 0.3067 | 0.8269 |
| ipTM | I | 0.6761 | 0.6195 | 0.1686 | 0.5124 |
| ipSAE | I | 0.6271 | 0.5732 | 0.1714 | 0.5169 |
| ipTM | F | 0.6594 | 0.5914 | 0.1613 | 0.4822 |
| ipSAE | F | 0.6315 | 0.5597 | 0.1692 | 0.5095 |

### Interpretation

- ipSAE did not improve ranking over ipTM: target-equal AUC and AP were lower at every evaluated layer.
- ipSAE raw Brier was lower, but its log loss was worse; after local I/F adaptation, its Brier and log loss were also worse than ipTM.
- Under pair source weighting, the small ipSAE Brier advantage in U disappeared.
- The supported statement is that ipSAE_min in these data did not provide a stable incremental improvement over Boltz-2 ipTM for the evaluated ranking, probability-error, or local-adaptation outcomes.

## C02: paired cluster-bootstrap intervals

All intervals are 95% percentile intervals from 2,000 valid within-target sequence-cluster bootstrap replicates. U, I, and F used the same sampled rows in each replicate. Negative Brier/log-loss differences mean the left, later layer has lower error.

### Main primary-six

- `I-U` Brier difference: -0.1393 (95% interval -0.1572 to -0.1192).
- `I-U` log-loss difference: -0.3533 (-0.3987 to -0.3019).
- `F-I` Brier difference: -0.0076 (-0.0141 to -0.0015).
- `F-I` log-loss difference: -0.0274 (-0.0427 to -0.0131).
- Macro BSS-transport: U -0.0809 (-0.1133 to -0.0491); I 0.3859 (0.3332 to 0.4308); F 0.4121 (0.3679 to 0.4522).
- Relative macro-Brier reduction from U to I: 45.1% (39.5% to 49.8%).

### Main primary-five

- `I-U` Brier difference: -0.1347 (-0.1503 to -0.1173).
- `F-I` Brier difference: +0.0014 (-0.0026 to +0.0055).
- Macro BSS-transport: U -0.0921 (-0.1221 to -0.0600); I 0.3304 (0.2871 to 0.3679); F 0.3252 (0.2783 to 0.3661).
- Relative macro-Brier reduction from U to I: 42.9% (38.4% to 46.8%).

### External eight-target cohort

- `I-U` Brier difference: -0.0671 (-0.0837 to -0.0489).
- `F-I` Brier difference: -0.0057 (-0.0105 to -0.0008); the log-loss interval crossed zero.
- Macro BSS-transport: U -0.0687 (-0.1131 to -0.0296); I 0.2170 (0.1622 to 0.2680); F 0.2391 (0.1787 to 0.2960).
- Relative macro-Brier reduction from U to I: 25.9% (19.3% to 31.9%).

### Interpretation boundary

- Intercept adaptation has a consistently favorable macro interval in all three datasets.
- Additional slope adaptation is favorable at the primary-six macro level, uncertain in primary-five, and heterogeneous by target; it should not be described as a stable per-target gain.
- These intervals quantify conditional precision within the observed targets and do not represent uncertainty over a future population of targets.

## Files

- C01 report: `analysis/postreview_revision_2026-11-21/ipSAE/feasibility_report.md`
- C01 manifest: `analysis/postreview_revision_2026-11-21/ipSAE/run_manifest.json`
- C01 macros: `analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_macro_metrics.csv`
- C01 comparisons: `analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_paired_comparisons.csv`
- C02 method/QC: `analysis/postreview_revision_2026-11-21/bootstrap/bootstrap_method_and_qc.md`
- C02 macro intervals: `analysis/postreview_revision_2026-11-21/bootstrap/bootstrap_macro_intervals.csv`
- C02 target intervals: `analysis/postreview_revision_2026-11-21/bootstrap/bootstrap_target_intervals.csv`
