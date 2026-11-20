# Reproduction guide

## 1. Environment

Use Python 3.12. The accepted numerical analyses recorded NumPy 2.4.4, pandas 3.0.2, SciPy 1.18.1, scikit-learn 1.9.0, Matplotlib 3.10.9, and DuckDB 1.5.3.

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

## 2. Available without upstream downloads

The repository includes accepted summary tables, QC records, figures, and frozen derived sources. Verify their integrity with:

```bash
python scripts/verify_publication.py
```

This checks published hashes and Python syntax. It does not refit models.

## 3. Recreate upstream table locations

```bash
python scripts/prepare_inputs.py --anthropic --proteinbase-parquet
```

Expected local locations:

- `data/proteinbase_latest.csv`
- `data/external_target_audit_2026-11-13/anthropic_design_summary.parquet`
- `data/external_target_audit_2026-11-13/anthropic_wetlab_summary.parquet`
- `data/external_target_audit_2026-11-13/anthropic_cofold_predictions.parquet`

These files are ignored by Git.

## 4. Execution order

The accepted analysis sequence is:

1. Build the pair-level master table and sequence-cluster assignments.
2. Run RQ1/RQ2 within-target discrimination and calibration analyses.
3. Run RQ3 target-held-out U/I/F and frozen external evaluation.
4. Run RQ4 selection and shared-budget decision analyses.
5. Run RQ5 complete-case Pinc analyses.
6. Run RQ6 predictor-pipeline comparisons.
7. Run the post-hoc RQ4E intercept diagnostic and RQ5RF random-forest sensitivity.
8. Rebuild the supplementary workbook, manuscript PDF, and QC summaries.

Each accepted run directory retains its scripts, commands, configuration snapshot, hash manifest, and QC result. Exact command strings are in the corresponding `COMMANDS.md`.

## 5. Important evidence boundaries

- Raw ipTM analysis tests the candidate probability interpretation; it does not imply that Boltz-2 developers defined ipTM as a wet-lab binding probability.
- External U is evaluated without destination-target labels. External I/F and shared-budget analyses are post hoc and label-dependent.
- RQ5RF is a post hoc sensitivity analysis over one prespecified random-forest family, not an algorithm search or confirmatory ranking.
- Full reproduction of Pinc requires the official implementation and source data under their own terms.
- Some frozen data include third-party identifiers and scores. Do not infer a new redistribution right from this repository's CC BY 4.0 data license.
