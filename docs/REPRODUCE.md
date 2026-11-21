# Reproduction guide

## 1. Environment

Use Python 3.12 with the pinned dependencies:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

The accepted runs recorded Python 3.12.16, NumPy 2.4.4, pandas 3.0.2, SciPy 1.18.1, scikit-learn 1.9.0, Matplotlib 3.10.9, and DuckDB 1.5.3. Exact environment records are retained in the accepted run manifests.

## 2. Repository integrity

```bash
python scripts/verify_publication.py
```

This checks the complete manifest, file sizes, SHA-256 hashes, reverse coverage, Python syntax, Markdown links, local paths, file URIs, credentials, and negative privacy test samples.

## 3. Upstream data

```bash
python scripts/prepare_inputs.py --proteinbase-parquet --anthropic
```

Pinned sources:

- Proteinbase `v2026-11-01`, revision `2d2825070da5ba77e0c7e1a0942146e8307e9282`, ODC-BY;
- Anthropic Claude binder-design release v1.0, revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`, CC BY 4.0 for data and documentation and MIT for included scripts.

Raw upstream datasets are ignored and are not redistributed.

## 4. Frozen derived analysis

The repository contains the input tables needed to rerun C01B, C02, figure generation, and PDF/QC generation. Run:

```bash
python scripts/reproduce_postreview.py --check
python scripts/reproduce_postreview.py --execute
```

The equivalent individual commands are:

```bash
python analysis/postreview_revision_2026-11-21/scripts/run_c01b_postreview.py
python analysis/postreview_revision_2026-11-21/scripts/run_c01_c02_postreview.py --bootstrap-only
python analysis/postreview_revision_2026-11-21/scripts/make_postreview_figures.py
python scripts/redraw_supplementary_figure_s7.py
python publication/manuscript/build_review_pdf.py
python publication/manuscript/qc_review_pdf.py
python analysis/postreview_revision_2026-11-21/scripts/qc_postreview_delivery.py
```

The exact script options are recorded in each run manifest. `run_c01_c02_postreview.py --bootstrap-only` reads the frozen U/I/F predictions and does not refit a model. `reproduce_postreview.py --check` verifies every required derived input before changing any output.

`--execute` writes regenerated outputs into the working tree. The result tables are expected to be byte-identical to the frozen tables. Run timestamps, Python patch versions, figure-library version strings, SVG object identifiers, and PDF creation dates can change without changing the plotted pixels, PDF page content, or manuscript text. The wrapper updates `RUN_STATUS.json` to identify its two supported tasks (`C01B` and `C02`) and should be run in a disposable copy when the released file hashes must remain unchanged.

The restored accepted-run records are also included:

- RQ3-F commands and stages: `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/COMMANDS.md`;
- RQ5 commands and stages: `analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/COMMANDS.md`;
- RQ6 command: `analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147/COMMANDS.md`;
- RQ5-RF commands: `analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855/COMMANDS.md`;
- RQ3-E permanent and external runs: `analysis/rq3e_outputs/run_20261114_164950/RUN_STATUS.json` and `analysis/rq3e_if_outputs/run_20261114_174403/RUN_STATUS.json`.

These records include the accepted commands, configuration snapshots, environment records, hash manifests, and QC outputs. Paths in copied text records must be resolved from the repository root; no external source tree is required for the frozen derived runs.

**Execution boundary:** these legacy RQ1/RQ2, RQ3, RQ4, RQ5, RQ6, and RQ5-RF files are accepted historical evidence, not a claim that every original upstream-to-result command can be rerun from this public snapshot without the documented external inputs. The supported frozen-derived entry point is `scripts/reproduce_postreview.py`; the detailed boundary is in `analysis/HISTORICAL_EVIDENCE_ONLY.md`.

## 5. C01 full upstream recalculation

The official IPSAE source is not redistributed. Clone the official repository at commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3` into a local directory and set:

```bash
export IPSAE_REPO=/path/to/IPSAE
export PROTEINBASE_CSV=/path/to/proteinbase_latest.csv
```

The complete recalculation then uses stored Boltz-2 PAE and mmCIF inputs under their source terms. Parameters are PAE cutoff 10 Å and CA-distance cutoff 10 Å. The released `boltz2_ipsae` and `boltz2_min_ipsae` fields are descriptive provenance fields, not the analysis endpoint.

## 6. Accepted analysis sequence

1. Build the pair-level master table and sequence-cluster assignments.
2. Run the within-target discrimination and hierarchical calibration analyses.
3. Run target-held-out U/I/F and frozen external evaluation.
4. Run fixed-budget and retrospective shared-budget analyses.
5. Run complete-case Pinc analyses.
6. Run predictor-pipeline descriptions and random-forest sensitivity.
7. Run C01 official ipSAE recomputation.
8. Run C01B ipTM-plus-ipSAE complementarity analysis.
9. Run C02 paired sequence-cluster bootstrap intervals.
10. Rebuild figures, the workbook, the PDF, and all QC records.

## 7. Evidence boundaries

- Raw ipTM is evaluated as a candidate experimental-binding probability only among designs already selected for assay.
- External U uses no destination-target labels; external I/F is post hoc and label-dependent.
- C01 and C01B are post hoc sensitivity extensions.
- C02 intervals describe conditional precision within the observed targets and are not future-target prediction intervals.
- The shared-budget analyses are retrospective simulations and do not validate a general pilot strategy or prospective assay savings.
- Pinc and ipSAE complete-case comparisons are restricted to the evaluated rows, folds, and weighting schemes.
