# C03 evidence and wording lock

Date: 2026-11-21

## Scope

This record fixes the scientific values and interpretation boundaries to be used in C04. It does not introduce a new analysis and does not replace any frozen result.

## Evidence hierarchy

- The original ipTM ranking, probability-error, target-held-out U/I/F, external-transport, and shared-budget analyses retain their existing evidence levels.
- ipSAE-only and ipTM-plus-ipSAE analyses are post hoc sensitivity extensions added after the original ipTM analyses had been examined.
- The bootstrap resampling in C02 describes conditional precision in the observed targets. It does not estimate variability across future targets.
- Local-label I/F results remain distinct from zero-label U evaluation.

## Canonical C01B wording

Use:

> On the same complete-case pairs, neither ipSAE alone nor a post hoc two-feature logistic model combining ipSAE with ipTM consistently improved discrimination or probability error across target-held-out transport and local-adaptation analyses. The small average reductions in probability error were inconsistent across targets and weighting schemes and were accompanied by slightly lower AUC.

Do not describe the two-feature model as prespecified. Its execution rules were frozen before examining the combination result, but the analysis was added after the single-score results.

Do not write that ipSAE and ipTM are not complementary, that ipSAE is generally ineffective, or that the released ipSAE field is erroneous.

## C01 numerical lock

- Cohort: 2,074 complete-case pairs, 740 primary positives, six targets.
- The same ten PD-L1 pairs lack Pinc and both released ipSAE fields. No additional loss occurred.
- Official implementation: DunbrackLab/IPSAE commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`, `ipsae.py`, PAE cutoff 10 A, distance cutoff 10 A.
- Official recomputation succeeded for 2,074/2,074 pairs.
- Released `boltz2_min_ipsae` reproduced the official directional minimum in 1,048/2,074 pairs and the other released ipSAE field reproduced the official directional maximum in 821/2,074 pairs.
- Main analysis therefore uses recomputed `ipSAE_min`. The released field remains valid provenance, but its generation protocol should not be assumed identical to the current official implementation.
- No Boltz-2 affinity metric was available in the frozen Proteinbase release.

### ipTM-only versus ipSAE-only

Target-equal results on the 2,074 common pairs:

| Layer/model | AUC | AP | Brier | Log loss |
|---|---:|---:|---:|---:|
| ipTM raw | 0.6770 | 0.6132 | 0.4154 | 1.2224 |
| ipSAE raw | 0.6658 | 0.5903 | 0.3341 | 1.6792 |
| ipTM U | 0.6770 | 0.6132 | 0.3139 | 0.9204 |
| ipSAE U | 0.6658 | 0.5903 | 0.3067 | 0.8269 |
| ipTM I | 0.6761 | 0.6195 | 0.1686 | 0.5124 |
| ipSAE I | 0.6271 | 0.5732 | 0.1714 | 0.5169 |
| ipTM F | 0.6594 | 0.5914 | 0.1613 | 0.4822 |
| ipSAE F | 0.6315 | 0.5597 | 0.1692 | 0.5095 |

- ipSAE did not improve AUC or AP at any layer.
- The lower ipSAE raw/U Brier was accompanied by worse log loss, did not survive source-pair weighting in U, and reversed after local I/F adaptation.
- Report the full trade-off; do not present an isolated Brier improvement.

## C01B numerical lock

Target-equal combined model versus ipTM-only:

| Layer | Delta AUC | Delta AP | Delta Brier | Delta log loss | Delta BSS-transport |
|---|---:|---:|---:|---:|---:|
| U | -0.0062 | +0.0008 | -0.0011 | -0.0211 | +0.0017 |
| I | -0.0060 | -0.0021 | -0.0009 | -0.0047 | +0.0034 |
| F | -0.0043 | -0.0005 | -0.0007 | -0.0021 | +0.0031 |

- Negative Brier/log-loss differences favor the combination. Positive AUC/AP/BSS differences favor it.
- Under target-equal weighting, Brier improved in 4/6 targets at U, 3/6 at I, and 3/6 at F.
- Under source-pair weighting, Brier improved in 4/6, 4/6, and 2/6 targets at U, I, and F, respectively.
- The conditional ipSAE slope was negative in all six source folds and under both weightings. This may be reported descriptively in the supplement because ipTM and ipSAE are strongly correlated; it must not be interpreted causally or mechanistically.

## C02 interval lock

All intervals are 95% percentile intervals from 2,000 valid within-target 90%-sequence-cluster bootstrap replicates.

| Dataset/quantity | Point | 95% interval |
|---|---:|---:|
| Primary-six relative macro-Brier reduction, U to I | 45.1% | 39.5% to 49.8% |
| Primary-six I-U Brier difference | -0.1393 | -0.1572 to -0.1192 |
| Primary-six I-U log-loss difference | -0.3533 | -0.3987 to -0.3019 |
| Primary-six F-I Brier difference | -0.0076 | -0.0141 to -0.0015 |
| Primary-six F-I log-loss difference | -0.0274 | -0.0427 to -0.0131 |
| Primary-six U BSS-transport | -0.0809 | -0.1133 to -0.0491 |
| Primary-six I BSS-transport | +0.3859 | +0.3332 to +0.4308 |
| Primary-six F BSS-transport | +0.4121 | +0.3679 to +0.4522 |
| Primary-five relative macro-Brier reduction, U to I | 42.9% | 38.4% to 46.8% |
| Primary-five F-I Brier difference | +0.0014 | -0.0026 to +0.0055 |
| External eight-target relative macro-Brier reduction, U to I | 25.9% | 19.3% to 31.9% |
| External I-U Brier difference | -0.0671 | -0.0837 to -0.0489 |
| External F-I Brier difference | -0.0057 | -0.0105 to -0.0008 |
| External F-I log-loss difference | -0.0118 | -0.0255 to +0.0025 |

- Intercept adaptation has a consistently favorable macro interval in all three datasets.
- Additional slope adaptation is not stable across the five-target analysis and is heterogeneous by target.
- The external F-I Brier interval is favorable, but its log-loss interval crosses zero.
- Do not call these intervals future-target prediction intervals.

## External report contribution boundary

Verified original report:

> Claude Science and Shanehsazzadeh, A. (2026) Autonomous de novo protein binder design with Claude. Technical report, Anthropic. https://assets.claude.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf

- The report described autonomous binder campaigns against 16 targets, two-vendor experimental validation, post hoc re-scoring under ten co-folding predictors, rank-based analyses, and relationships between co-folding confidence and binding or affinity.
- The present study asks whether a released score or target-transferred mapping supplies calibrated experimental-binding probabilities and whether local labels and pilot costs change experimental allocation.
- Cite both the technical report and the data release. Do not imply that the original report evaluated the present target-held-out probability-transport or shared-budget questions.

## Abstract decision

Do not add ipSAE-only or ipTM-plus-ipSAE results to the Abstract. The small and inconsistent differences would displace the main transport, local-adaptation, and decision-value findings.

## C04 wording constraints

- Prefer "did not show a stable incremental improvement in the evaluated cohort" over universal negative claims.
- Keep probability error and discrimination together whenever a small beneficial error difference is mentioned.
- State whether a result uses destination-target labels.
- State that external transport changes target/campaign and score-generation conditions together.
- Define the estimand as the probability of experimental binding among designs already selected and submitted for assay, not the success probability of an unselected candidate library.
- Replace "not a universal selection threshold" with "cannot be interpreted as a common cross-target probability threshold for experimental success".
