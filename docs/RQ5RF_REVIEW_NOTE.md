# RQ5-RF English results and supplementary-table note — for review only

Status: incorporated into the main manuscript, supplementary guide, and Supplementary Tables S1–S16 workbook on 18 November 2026.

## Proposed Results subsection

### Random forests alter probability error but do not consistently outperform pooled logistic regression on held-out targets

To examine whether the unstable Pinc increment could be explained by the linear functional form, we performed a post hoc two-by-two sensitivity analysis comparing pooled ridge logistic regression and random forests using ipTM alone or ipTM plus Pinc. All models used the same complete-case pairs, outer folds, source weights, and adaptation rules. The pooled logistic comparators omitted source-target intercepts so that matched comparisons varied model family without simultaneously changing the intercept structure.

Under six-target, target-equal leave-one-target-out evaluation without labels from the held-out target, the ipTM-only random forest lowered Brier from 0.2826 to 0.2624 and log loss from 0.7590 to 0.7254 relative to its logistic comparator, while BSS-transport increased from −0.0001 to 0.0690. Brier was lower in five of six targets, although the difference in SARS-CoV-2 RBD was negligible. These probability-error gains were accompanied by poorer discrimination: AUC decreased from 0.6515 to 0.5684 and AP from 0.6002 to 0.5536. BSS-target-oracle also remained negative (−0.7152). Thus, lower average probability error did not establish improved ranking or transportable calibration.

The two-feature random forest likewise did not yield a stable improvement. Relative to the two-feature logistic model, it lowered Brier and log loss in the primary six-target analysis, but this direction reversed after refitting on the five-target analysis set (Brier difference, +0.0466; log-loss difference, +0.1725; BSS-transport difference, −0.1294). Within the random-forest family, adding Pinc worsened macro-averaged Brier, log loss, and BSS-transport in all six analysis-set-by-weighting configurations. In the primary analysis, Brier increased from 0.2624 for ipTM alone to 0.2894 for ipTM plus Pinc, log loss increased from 0.7254 to 0.7870, and BSS-transport decreased from 0.0690 to −0.0240. In the primary analysis, the full random-forest-plus-Pinc model also failed to outperform the simplest pooled logistic ipTM model. With destination-target labels, the fully adapted two-feature random forest achieved Brier 0.1656 and BSS-target-oracle 0.0452 in the new pooled comparison; this local-label result does not alter the interpretation for prediction without destination-target labels and is not a pure model-family contrast with the previous target-intercept analysis. Overall, these exploratory findings show metric- and target-dependent changes under the evaluated random-forest grid, with no consistent improvement over pooled logistic regression when predicting held-out targets without their experimental labels, and no stable incremental contribution from Pinc (Supplementary Table S16).

## Proposed supplementary-table entry

### Table S16. Random-forest model-complexity sensitivity on common complete cases

**Purpose.** Determine whether a more flexible functional form changes the cross-target probability-error and incremental-Pinc conclusions when samples, folds, feature sets, source weights, and local-adaptation rules are held fixed.

**Design.** A post hoc exploratory two-by-two comparison of pooled ridge logistic regression and random forests, each using ipTM alone or ipTM plus Pinc. Random-forest hyperparameters were selected within source training data using leave-one-source-target-out log loss. The six configurations comprised primary-six, primary-five, and strict-six analysis sets under target-equal and source-pair weighting. The primary-six, primary-five, and strict-six sets contained 2,074 pairs and 740 positives, 1,941 pairs and 721 positives, and 2,074 pairs and 674 positives, respectively. U used no destination-target labels; I and F used the corresponding frozen local-adaptation procedures.

**Evidence boundary.** This was not an algorithm search or a confirmatory model ranking. Only the prespecified random-forest family was added. The pooled logistic models omitted source-target intercepts; consequently, the earlier target-intercept RQ5 models are descriptive anchors rather than pure model-family comparators. Differences between ipTM-plus-Pinc and ipTM-only random forests combine possible nonlinear main effects and interactions and do not identify a pure interaction. No independent external Pinc cohort was available, and local-label I/F results cannot establish transport to held-out targets without their experimental labels.

**Workbook panels.** A. Analysis-set, sample, and model-definition audit; B. Complete six-configuration U macro-averaged metrics; C. Target-level U metrics and prespecified pairwise differences; D. Local-adaptation I/F results on identical outer-test rows; E. Reliability summaries and non-evaluable calibration diagnostics; F. Random-forest candidate selection, implementation provenance, and QC summary.

**Interpretation.** Negative Brier or log-loss differences and positive BSS, AUC, or AP differences favour the left-hand model in each prespecified contrast. The prespecified directional-consistency screen was not met. Results therefore support a statement of local or inconsistent changes under the evaluated models, not a claim that all nonlinear methods are ineffective or that random forests resolve calibration transport.

## Numerical basis — internal review only

- `analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855/rq5rf_six_configuration_summary.csv`
- `analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855/rq5rf_macro_metrics.csv`
- `analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855/rq5rf_target_metrics.csv`
- `analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/rq5_adaptation_macro_summary.csv`
- `analysis/RQ5RF_model_complexity_frozen_protocol_2026-11-17.md`

## Integration notes — not manuscript prose

1. Insert the Results subsection after the existing complete-case Pinc paragraphs; do not replace the original target-intercept RQ5 results.
2. Add Table S16 rather than silently changing the existing S1–S15 workbook. If incorporated, update the supplement guide, workbook cover, source manifest, and the main-text supplementary-table inventory together.
3. Do not cite the old and new absolute F-layer values as a head-to-head random-forest-versus-logistic comparison; their intercept structures differ.
4. Keep the analysis labelled post hoc and exploratory. Do not move it into the Abstract or use it to revise the RQ1–RQ4 conclusions.
