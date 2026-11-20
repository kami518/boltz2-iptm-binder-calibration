# Boltz-2 ipTM binder-calibration reproducibility release

Release-candidate repository for the manuscript:

> Evaluating Boltz-2 ipTM as an experimental-success probability across binder-design targets

This repository contains the analysis code, accepted numerical summaries, publication figures, and manuscript-facing outputs needed to audit the reported results. It is prepared for a public GitHub release and a later archival DOI, but no public remote or DOI has been created yet.

## What is included

- analyses for within-target discrimination and raw probability error;
- target-held-out transport and local intercept/slope adaptation;
- frozen external evaluation and post-hoc external adaptation;
- experimental-selection and shared-budget analyses;
- Pinc complete-case comparisons and random-forest sensitivity summaries;
- paired predictor-pipeline descriptions;
- accepted run scripts, commands, QC records, and hash manifests;
- main and supplementary figures and manuscript source text.

Three generated files are distributed as GitHub Release assets rather than Git blobs:

- the complete review PDF;
- the supplementary Excel workbook;
- the full pilot-parameter table.

Their filenames and SHA-256 hashes are recorded in `publication/RELEASE_ASSETS.md`.

## What is not included

- the raw Proteinbase CSV;
- the Anthropic release's raw per-design data, structure files, vendor reports, or PAE matrices;
- the Overath mirror or Overath-derived row-level tables;
- third-party model weights or bundled third-party software;
- the official Pinc implementation;
- local execution-environment paths or credentials.

See `docs/DATA_SOURCES.md` for the official pinned releases and `docs/REPRODUCE.md` for acquisition and verification steps.

## Quick integrity check

```bash
python scripts/verify_publication.py
```

The command checks every published file against `data/PUBLICATION_MANIFEST.csv`, compiles all included Python files, and scans text files for local home paths and credential-like content.

## Licenses

- Analysis code written for this study: MIT, see `LICENSE`.
- Derived data, documentation, figures, and manuscript outputs: CC BY 4.0, see `DATA_LICENSE.md`.
- Upstream data retain their own attribution and third-party terms; see `THIRD_PARTY_NOTICES.md`.

## Citation

The author list, repository URL, and DOI are intentionally pending. Do not cite this candidate release as a final archived record. Fill `CITATION.cff.template` and `ZENODO_METADATA.template.json` only after author approval and repository publication.
