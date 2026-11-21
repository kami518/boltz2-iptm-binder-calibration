# GitHub Release assets: v1.1.0

Three large generated files are distributed as deterministic ZIP archives split into 1 MB pieces:

- `PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf`;
- `Supplementary_Tables_S1_S16_v1.1.0.xlsx`;
- `rq4c_pilot_parameters_v1.1.0.csv`.

The machine-readable archive, piece, and final-file hashes are in:

`data/release_assets_manifest.json`

## Reassembly

Download every `v1.1.0` asset into one directory, then run:

```bash
python scripts/reassemble_release_assets.py /path/to/downloaded/assets
```

The script verifies every piece, rebuilds each archive, verifies archive hashes, safely extracts the archive, and verifies the three final outputs.

The Git tag identifies the exact source snapshot. The Release assets identify the PDF, supplementary workbook, and pilot-parameter table generated from that snapshot.
