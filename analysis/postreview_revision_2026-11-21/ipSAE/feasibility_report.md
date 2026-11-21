# C01 ipSAE feasibility and extension

**Evidence level:** post hoc sensitivity extension added after the original ipTM analyses had been examined.

## Feasibility

- Complete-case cohort: 2074 pairs, 740 primary positives, six targets.
- Missing set: the same 10 PD-L1 pairs lack both Pinc and Boltz-2 ipSAE inputs; no additional loss was introduced by the ipSAE comparison.
- Released fields audited: `boltz2_ipsae` and `boltz2_min_ipsae`.
- Official implementation: DunbrackLab/IPSAE commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`, `ipsae.py`, PAE cutoff 10 Å, distance cutoff 10 Å.
- Official script SHA-256: `10cf9b08c68c91e06cb28526cf2026f47a3980c9048fd3226d13e3304eaf1c27`.
- Stratified sample QC: 18/18 passed.
- Released minimum field exactly matched the official directional minimum in 8/18 sampled pairs; the comparison is recorded but the released field was not used as the algorithmic gold standard.
- Across the full batch, the released minimum field matched in 1048/2074 pairs and the released non-minimum field matched in 821/2074 pairs.
- The extension therefore uses freshly recomputed `ipSAE_min` from the stored Boltz-2 PAE/mmCIF pair with the pinned official implementation.
- PAE dimensions, finite values, chain labels, and residue-token sums were checked for all 2074 pairs.
- Boltz-2 affinity metrics found in the frozen Proteinbase file: none. No affinity output was available for this extension.

## Extension design

- The ipTM comparator and ipSAE model were both refitted on the same 2,074 rows, frozen LOTO and target-CV5 assignments, 90% sequence clusters, and source weighting.
- The source model used one intercept per source target and one common slope.
- I retained the frozen U prediction as an offset and estimated an intercept only.
- F estimated intercept plus a ridge-shrunk slope change; lambda was selected only within adaptation-training folds by pooled inner log loss.
- These results are an exploratory extension on complete cases and do not replace the frozen 2,084-pair ipTM analysis.
