# C01B: ipTM + ipSAE complementarity

## Design

- Cohort: the same 2,074 complete-case pairs, including 740 primary positives and six targets.
- Models compared: `ipTM_only`, `ipSAE_only`, and `ipTM_plus_ipSAE`.
- Evidence level: post hoc sensitivity extension added after the single-score C01 results had been examined. The combination execution rules were frozen before the combination result was inspected.
- U used one intercept per source target and two jointly estimated common slopes on `logit(ipTM)` and `logit(ipSAE_min)`.
- I updated only an overall intercept; F updated an overall intercept plus one scalar calibration slope.
- The relative ipTM/ipSAE weights were not re-estimated locally.
- Frozen primary-six LOTO and target-CV5 assignments, 90% sequence clusters, and both existing source-weighting schemes were reused.
- No feature selection, penalty tuning, Pinc combination, or random forest was run.

## Target-equal combination metrics

| Layer | AUC | AP | Brier | Log loss | BSS-transport |
|---|---:|---:|---:|---:|---:|
| U | 0.6708 | 0.6139 | 0.3128 | 0.8992 | -0.0974 |
| I | 0.6701 | 0.6174 | 0.1676 | 0.5077 | 0.3928 |
| F | 0.6551 | 0.5910 | 0.1606 | 0.4800 | 0.4163 |

## Target-equal paired comparisons

| Comparison/layer | AUC difference | AP difference | Brier difference | Log-loss difference | BSS difference |
|---|---:|---:|---:|---:|---:|
| ipTM_plus_ipSAE_vs_ipSAE_only/F | 0.0236 | 0.0313 | -0.0087 | -0.0295 | 0.0320 |
| ipTM_plus_ipSAE_vs_ipSAE_only/I | 0.0430 | 0.0442 | -0.0038 | -0.0092 | 0.0166 |
| ipTM_plus_ipSAE_vs_ipSAE_only/U | 0.0050 | 0.0236 | 0.0061 | 0.0723 | -0.0225 |
| ipTM_plus_ipSAE_vs_ipTM_only/F | -0.0043 | -0.0005 | -0.0007 | -0.0021 | 0.0031 |
| ipTM_plus_ipSAE_vs_ipTM_only/I | -0.0060 | -0.0021 | -0.0009 | -0.0047 | 0.0034 |
| ipTM_plus_ipSAE_vs_ipTM_only/U | -0.0062 | 0.0008 | -0.0011 | -0.0211 | 0.0017 |

Negative Brier/log-loss differences favor the combination; positive AUC/AP/BSS differences favor it.
These are point estimates from the existing folds, not a new model-selection exercise.

## Direction consistency versus ipTM-only

| Weighting | Layer | AUC improved | AP improved | Brier improved | Log loss improved | BSS improved |
|---|---|---:|---:|---:|---:|---:|
| target_equal | U | 1/6 | 1/6 | 4/6 | 4/6 | 4/6 |
| target_equal | I | 1/6 | 2/6 | 3/6 | 3/6 | 3/6 |
| target_equal | F | 3/6 | 3/6 | 3/6 | 3/6 | 3/6 |
| pair | U | 2/6 | 1/6 | 4/6 | 4/6 | 4/6 |
| pair | I | 1/6 | 1/6 | 4/6 | 4/6 | 4/6 |
| pair | F | 0/6 | 3/6 | 2/6 | 2/6 | 2/6 |

## Conditional source slopes

| Weighting | ipTM slope range | ipSAE slope range |
|---|---:|---:|
| pair | 0.5277 to 0.8958 | -0.1117 to -0.0339 |
| target_equal | 0.7180 to 1.3532 | -0.2018 to -0.0765 |

The conditional ipSAE slope was negative in every source fold and under both weighting schemes. Because ipTM and ipSAE are strongly correlated, this sign should not be interpreted as a causal or mechanistic effect. It does show that adding ipSAE did not provide a stable positive marginal signal after ipTM.

## Interpretation

The combination improved at least one probability-error metric in some configurations, but ranking was not consistently improved and the gains were small and target-dependent. The appropriate manuscript statement is that ipSAE, either alone or jointly with ipTM, did not show a stable incremental improvement in this evaluated cohort; this does not establish that ipSAE is generally uninformative or that no other calibration map could work.
