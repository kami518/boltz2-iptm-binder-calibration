# Data sources and pinned releases

## Proteinbase

- Official page: https://huggingface.co/datasets/adaptyvbio/proteinbase
- Version/tag: `v2026-11-01`
- Revision: `2d2825070da5ba77e0c7e1a0942146e8307e9282`
- License: ODC-BY
- Required attribution: "This work used Proteinbase by Adaptyv Bio under ODC-BY license"

The local source CSV used in the manuscript was saved on 5 November 2026 and recorded as accessed for analysis on 7 November 2026. Its SHA-256 was:

`8eac08bfc931409e0703eda904cdd1066c9d4a2ceb46af43763f94c74dabb9cc`

The raw 51.9 MB CSV is not redistributed here. Run:

```bash
python scripts/prepare_inputs.py --proteinbase-parquet
```

to download the pinned official parquet and create a semantically equivalent CSV under the ignored local-data path. For an exact byte-level re-creation of the working snapshot, obtain the public CSV from Proteinbase and compare the required identifiers and evaluation records.

## Anthropic Claude protein binder design data

- Official page: https://huggingface.co/datasets/Anthropic/claude-protein-binder-design
- Release/version: v1.0
- Revision: `9e1b81696da46835e9e9cde9a3da976e0abc92ab`
- Data and documentation: CC BY 4.0
- Included scripts: MIT
- Official citation author: The protein binder design data release authors

The following upstream tables can be downloaded at the pinned revision:

- `data/tables/design_summary.parquet`
- `data/tables/wetlab/summary.parquet`
- `data/tables/insilico/cofold_predictions.parquet`

Run:

```bash
python scripts/prepare_inputs.py --anthropic
```

The raw release, structure files, PAE matrices, vendor images, and local audit copies are not redistributed.

## Overath mirror

The Overath mirror and derived row-level overlap tables are excluded because a redistribution license for the mirror was not established. The manuscript's overlap audit must therefore be reproduced using a user-obtained copy that follows the source's original terms.

## Official Pinc implementation

Pinc is not bundled. The accepted analysis records the official repository:

https://git.mpi-cbg.de/tothpetroczylab/Pinc.git

and commit `57afb08e29fe97379db5b155da5fb78c4c8715d1`, run with R 4.6.2 and jsonlite 2.0.0. The accepted run used symmetrised Pinc only.

## Official IPSAE implementation

IPSAE source code is not bundled. The accepted C01 recomputation used:

https://github.com/DunbrackLab/IPSAE

at commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`, with a 10 Å PAE cutoff and a 10 Å CA-distance cutoff. The released Proteinbase `boltz2_ipsae` and `boltz2_min_ipsae` fields were retained for provenance but were not used as the analysis endpoint.
