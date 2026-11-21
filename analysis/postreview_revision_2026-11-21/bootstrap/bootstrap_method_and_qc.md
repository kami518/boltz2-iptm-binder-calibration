# C02 paired cluster bootstrap

- Inputs were frozen pair-level U/I/F predictions; no model was refitted.
- Within each held-out target, 90%-identity sequence clusters were resampled with replacement.
- U, I, and F were evaluated on the same sampled clusters and rows in every replicate.
- Each target metric was recalculated before taking the target-equal macro average.
- Requested valid replicates: 2000; RNG seed: 20261121; percentile interval: 2.5% to 97.5%.
- Intervals describe conditional precision in the observed targets and do not represent uncertainty over future targets.

## Replicate QC

| Dataset | Pairs | Targets | Valid | Failed attempts |
|---|---:|---:|---:|---:|
| main_primary_six | 2084 | 6 | 2000 | 0 |
| main_primary_five | 1951 | 5 | 2000 | 0 |
| external_primary_weak_plus | 720 | 8 | 2000 | 0 |

The interval for `relative_macro_U_to_I_brier_reduction` is calculated by recomputing `(Brier_U-Brier_I)/Brier_U` in every replicate.
