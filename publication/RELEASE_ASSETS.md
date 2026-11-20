# GitHub Release assets

The complete review PDF and supplementary Excel workbook are distributed as ZIP archives split into 1 MB pieces to remain reliable over constrained connections. The compressed pilot-parameter CSV is distributed as one ZIP asset.

The machine-readable list of every archive, piece, byte count, and SHA-256 is:

`data/release_assets_manifest.json`

## Reassembly

Download all assets attached to release `v1.0.0` into one directory, then run:

```bash
python scripts/reassemble_release_assets.py /path/to/downloaded/assets
```

The script:

1. verifies every download piece;
2. concatenates pieces for archives that were split;
3. verifies each complete ZIP archive;
4. safely extracts the files; and
5. verifies the final PDF, XLSX, and CSV against their original SHA-256 values.

The Git tag identifies the exact code and summary-table state. The asset manifest identifies the accompanying large generated outputs.
