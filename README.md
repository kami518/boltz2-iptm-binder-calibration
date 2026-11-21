# Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design

Reproducibility materials for the manuscript above, prepared as the `v1.1.0` submission candidate.

Author: KangMing Chen, College of Biological Science and Medical Engineering, Donghua University, Shanghai, China.

Repository: https://github.com/kami518/boltz2-iptm-binder-calibration

Reserved archival DOI: `10.5281/zenodo.24400262` (not activated in this candidate).

## Included evidence

- frozen Proteinbase and Claude binder-design dataset references;
- within-target discrimination and raw probability error;
- target-held-out U/I/F transport and local adaptation;
- external evaluation and nonconfirmatory low-event stress test;
- fixed-budget and retrospective shared-budget analyses;
- complete-case Pinc analyses;
- C01 official ipSAE recomputation;
- C01B post hoc ipTM-plus-ipSAE complementarity analysis;
- C02 paired sequence-cluster bootstrap intervals;
- random-forest model-complexity sensitivity;
- paired predictor-pipeline descriptions;
- publication figures, manuscript source, supplement guide, run manifests, and QC records.

## Large release assets

Three generated files are distributed through the GitHub Release rather than as ordinary Git blobs:

- `PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf`;
- `Supplementary_Tables_S1_S16_v1.1.0.xlsx`;
- `rq4c_pilot_parameters_v1.1.0.csv`.

Their pieces, hashes, and reassembly instructions are recorded in `data/release_assets_manifest.json` and `publication/RELEASE_ASSETS.md`.

## Exclusions

- the raw Proteinbase CSV;
- raw Anthropic per-design files, vendor reports, structure files, and PAE matrices;
- the Overath mirror and row-level Overath-derived tables whose redistribution basis is not established;
- model weights;
- the official IPSAE and Pinc source repositories;
- local execution paths and credentials.

See `docs/DATA_SOURCES.md`, `docs/REPRODUCE.md`, and `THIRD_PARTY_NOTICES.md`.

## Integrity check

```bash
python scripts/verify_publication.py --asset-dir /path/to/release-assets
```

The verifier checks the complete hash manifest, reverse file coverage, Python syntax, Markdown links, local paths, credentials, XLSX strings, PDF text, and PDF link annotations.

For the supported frozen-derived workflow:

```bash
python scripts/reproduce_postreview.py --check
python scripts/reproduce_postreview.py --execute
```

Legacy RQ1/RQ2, RQ3, RQ4, RQ5, RQ6, and RQ5-RF files are retained as historical evidence. They are not presented as an unrestricted upstream-to-result executable pipeline.

## Licenses

- analysis code: MIT;
- derived data, documentation, figures, and manuscript outputs: CC BY 4.0;
- upstream data and third-party outputs: their original terms.

## Release status

The source, figures, workbook, PDF, release assets, and Zenodo archive are prepared locally. No Git tag, GitHub Release, Zenodo record, or active DOI was created during candidate assembly. Remote publication remains subject to mainline approval.
