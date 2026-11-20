# Data-source verification record

Verification date: 18 November 2026

This record documents the verified external release metadata used in the manuscript. It distinguishes repository metadata from local download and analysis timestamps. It does not infer an author list, DOI, or license grant that is absent from the official source.

## Anthropic external dataset

- Official dataset repository: `Anthropic/claude-protein-binder-design`
- Dataset page: https://huggingface.co/datasets/Anthropic/claude-protein-binder-design
- Official release title: "Protein binder design data release: de novo miniprotein binders against 16 targets with two-vendor binding data and co-folded structure predictions"
- Official version: `v1.0`
- Verified repository revision: `9e1b81696da46835e9e9cde9a3da976e0abc92ab`
- Repository creation timestamp reported by Hugging Face: `2026-08-17T19:36:03Z`
- Last modification timestamp reported by Hugging Face: `2026-08-18T21:35:41Z`
- Citation author string in the official `data/CITATION.cff`: `The protein binder design data release authors`
- Data and documentation license: Creative Commons Attribution 4.0 International (`CC-BY-4.0`)
- Included scripts license: MIT
- Third-party material: retains its own terms, including conditions attached to some structure-prediction outputs
- Local source files were retrieved on 13 November 2026; the official repository and pinned files were rechecked on 18 November 2026.

Official citation evidence:

- `data/CITATION.cff` at revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`
- `data/LICENSE.md` at the same revision
- Dataset card at the same revision

Local integrity checks:

| File | SHA-256 |
|---|---|
| `data/external_target_audit_2026-11-13/anthropic_design_summary.parquet` | `3e796869b8dc0c90d7ab0e12daf25308c0cf272d04c6a5e27281f903ecd18744` |
| `data/external_target_audit_2026-11-13/anthropic_wetlab_summary.parquet` | `f80c3cfce67cb7da53576a020d476c5b74efed44e61c6fceaec9ee7de5722f88` |
| `data/external_target_audit_2026-11-13/anthropic_cofold_predictions.parquet` | `ba48010610daa9b9d184aa1ddb419beebdf728ccf745f5a8f3b48e4041e1a87d` |

## Proteinbase

- Official dataset repository: `adaptyvbio/proteinbase`
- Dataset page: https://huggingface.co/datasets/adaptyvbio/proteinbase
- Official snapshot tag and version: `v2026-11-01`
- Verified repository revision for the tag: `2d2825070da5ba77e0c7e1a0942146e8307e9282`
- Revision timestamp reported by Hugging Face: `2026-11-01T03:00:36Z`
- License tag and dataset-card license: Open Data Commons Attribution License (`ODC-BY`)
- The dataset card requires this attribution text: "This work used Proteinbase by Adaptyv Bio under ODC-BY license"
- The local CSV was saved on 5 November 2026. The study records analysis access on 7 November 2026. The manuscript states both dates separately to avoid conflating download and analysis access.
- Local CSV SHA-256: `8eac08bfc931409e0703eda904cdd1066c9d4a2ceb46af43763f94c74dabb9cc`
- The frozen CSV and the dataset-card parquet were compared across IDs, names, sequences, and normalized evaluation records; all 4,280 IDs matched with zero content discrepancies.

## Redistribution conclusion

- Proteinbase-derived material may be redistributed and adapted under ODC-BY with the required attribution.
- Anthropic release data and documentation may be redistributed and adapted under CC BY 4.0 with attribution; included scripts use MIT.
- Third-party structure-prediction outputs, reference sequences, vendor images, and other third-party materials must remain subject to their own terms.
- The public reproducibility repository should distribute analysis code, frozen derived tables permitted by these licenses, software versions, seeds, and hashes. It should provide acquisition instructions rather than re-host any third-party material whose terms do not permit redistribution.
