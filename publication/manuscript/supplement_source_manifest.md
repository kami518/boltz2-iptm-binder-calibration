# Supplement source and output manifest

Release: `v1.1.0`

This manifest records the portable dependency closure and generated supplement outputs. Raw upstream data and third-party source code are excluded.

## Pinned upstream sources

- Proteinbase: `v2026-11-01, revision 2d2825070da5ba77e0c7e1a0942146e8307e9282`
- Anthropic binder-design release: `v1.0, revision 9e1b81696da46835e9e9cde9a3da976e0abc92ab`
- IPSAE: `DunbrackLab/IPSAE commit 6174cf9e71cb1bd660cc805856a18c4871a6dec3`
- Pinc: `official implementation commit 57afb08e29fe97379db5b155da5fb78c4c8715d1`

## Dependency closure

The complete machine-readable dependency manifest is `data/postreview_dependency_manifest.json`.

| Path | Source bytes | Source SHA-256 | Published bytes | Published SHA-256 |
|---|---:|---|---:|---|
| `analysis/v3_outputs/pair_master_v3.csv` | 1,071,826 | `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91` | 1,071,826 | `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91` |
| `analysis/audit_outputs/pair_audit.csv` | 1,847,747 | `1ce35c8ee601a05dd390d43a76814800d8429aa64fff30f60c196ad30d74e571` | 1,847,747 | `1ce35c8ee601a05dd390d43a76814800d8429aa64fff30f60c196ad30d74e571` |
| `analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv` | 2,525,569 | `0a549cb890f42948ced456c4a128d480ea06514d4d431d069bcc0f7dc10acaab` | 2,598,159 | `adce8de974e9d3c4ea2a1a4bf8b753f2460d842a9cb637aaf22f38a234ec0b18` |
| `analysis/rq12_outputs/sequence_clusters_rq12.csv` | 422,726 | `45283741709016fb5db3a114d2a108a236c73b73750f01c959ec2753d3ff9efb` | 422,726 | `45283741709016fb5db3a114d2a108a236c73b73750f01c959ec2753d3ff9efb` |
| `analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/loto_pair_assignments.csv` | 18,762,175 | `6ab011035bd0314d6fbaffe9cbe0548c52668ec2fc1695ccf262cbda984194ea` | 18,762,175 | `6ab011035bd0314d6fbaffe9cbe0548c52668ec2fc1695ccf262cbda984194ea` |
| `analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/target_cv5_pair_assignments.csv` | 13,964,456 | `ac837a1b7e3710ab51ce20cf01c451675092fe06f9ebd86c830e7fd50bf61796` | 13,964,456 | `ac837a1b7e3710ab51ce20cf01c451675092fe06f9ebd86c830e7fd50bf61796` |
| `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_adaptation_pair_predictions.csv` | 16,344,953 | `d6cd2ddaa71f4424d41172663f051d44a340f03bfee084a6e562bbbfcda21456` | 16,344,953 | `d6cd2ddaa71f4424d41172663f051d44a340f03bfee084a6e562bbbfcda21456` |
| `analysis/rq3e_if_outputs/run_20261114_174403/oof_predictions.csv` | 231,295 | `030c5ced908bea93e7eb9cebd4d94f230891dad8597abb060fb659ce502b79bf` | 231,295 | `030c5ced908bea93e7eb9cebd4d94f230891dad8597abb060fb659ce502b79bf` |
| `analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_pair_extension.csv` | 657,429 | `8fe6dbe5af786c6e4a8275bf9f54fa8e285a5e84f9e29c7d67f54bda66afd972` | 657,429 | `8fe6dbe5af786c6e4a8275bf9f54fa8e285a5e84f9e29c7d67f54bda66afd972` |
| `analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_target_metrics.csv` | 17,857 | `84cf972f48fa34ea6f84c0e48c95db6158519aeccc9f262e1099c22bd67db7ff` | 17,857 | `84cf972f48fa34ea6f84c0e48c95db6158519aeccc9f262e1099c22bd67db7ff` |
| `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_common.py` | 39,171 | `282a039571a6a1498949f1affa92fe174727608f9ab758c9f64549b8552e6221` | 40,575 | `f448ef64af34c04bd4041b5571c6f99d9503b21ab6dfd20335d107313866b751` |

## Release workbook

- File: `Supplementary_Tables_S1_S16_v1.1.0.xlsx`
- Bytes: 12,769,500
- SHA-256: `59a2d4d260527e1c9b3a0c73aa41f18799afa9802a3cf5a8a55f7e28b433f301`

## Supplementary figures

| File | Format | Bytes | SHA-256 |
|---|---|---:|---|
| `publication/figures/supplementary/Supplementary_Figure_S1.png` | PNG | 383,173 | `8f90f2a12a4a0fcdcfbde68c673751e330c6f3b0ec4116876a02fc3b76f83923` |
| `publication/figures/supplementary/Supplementary_Figure_S2.png` | PNG | 301,130 | `91b78a6d90552cb3ee5bf7492c21b5f944cd57d14514c80a6c12c65fedeaca62` |
| `publication/figures/supplementary/Supplementary_Figure_S3.png` | PNG | 96,230 | `a474af3d64d9c279021deb11bf9330b5771d86b8ca104cad25b42f959a239923` |
| `publication/figures/supplementary/Supplementary_Figure_S4.png` | PNG | 275,014 | `cda72ac8051498b579e0f21645f3c1104448764c26b93f195fc8b04c90738b7e` |
| `publication/figures/supplementary/Supplementary_Figure_S5.png` | PNG | 244,159 | `0ae1ab8604374b9d56b86627d7e078ce2f84c7410d3b0c169b4857ef5dcc770e` |
| `publication/figures/supplementary/Supplementary_Figure_S6.png` | PNG | 84,424 | `1b50a043c78e2120189bab076fa3e1cac81a8176196f7c98dbedd7a2981aeae8` |
| `publication/figures/supplementary/Supplementary_Figure_S7.png` | PNG | 255,699 | `337771cf166ec1bc701b2c2a60c1ae27f6c3d0cf5ecf52a2d430758b51ff9885` |
| `publication/figures/supplementary/Supplementary_Figure_S8.png` | PNG | 290,070 | `9e7bd5357af6d1c1c522865f24a27fa08896008470469af99747dd8f6879102b` |
| `publication/figures/supplementary/Supplementary_Figure_S9.png` | PNG | 163,488 | `6d06dd0d19084f7bdb9fa6713d7cf7152ff8cbc0db176b3f65c03f2a36f2ba4b` |
| `publication/figures/supplementary/Supplementary_Figure_S10.png` | PNG | 357,238 | `a3d578e3c1848da2002f3ca9b0ead9a258700500dbb0daafc89cc3db0218349f` |
| `publication/figures/supplementary/Supplementary_Figure_S11_ipSAE_postreview.png` | PNG | 271,207 | `64249b76313d120bbc4d2ecad9efa66b17174a66f6da439d1f8db342c35dd640` |
| `publication/figures/supplementary/Supplementary_Figure_S11_ipSAE_postreview.svg` | SVG | 101,246 | `433f5b08907ac9912c02645ec3215388d4da8975eb5010fd703659da2da3895f` |
| `publication/figures/supplementary/Supplementary_Figure_S11_ipSAE_postreview.pdf` | PDF | 30,746 | `f0cb7c397b69f55f7ef7e24dce021b9a3462b450f373d87e652847e7a10e5b8a` |
