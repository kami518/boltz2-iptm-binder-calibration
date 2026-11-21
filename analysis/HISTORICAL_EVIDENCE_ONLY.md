# Historical legacy analysis records

The restored RQ1/RQ2, RQ3, RQ4, RQ5, RQ6, and RQ5-RF files under `analysis/` are **accepted historical evidence**, not the public reproduction entry point.

They are retained for auditability: scripts, command histories, configuration snapshots, hash manifests, status files, and QC records document how the earlier analyses were performed. Some paths and expected hashes were transformed during the public-release path cleanup. In particular, path-bearing CSV and JSON files can legitimately have different source and published hashes.

Do not present these legacy scripts as a complete runnable upstream-to-result pipeline. Use:

```bash
python scripts/reproduce_postreview.py --check
python scripts/reproduce_postreview.py --execute
```

for the verified frozen-derived chain covering C01B, C02, figures, PDF, and delivery QC. Full C01 recomputation additionally requires the pinned official IPSAE source and the user-acquired raw Proteinbase inputs described in `docs/REPRODUCE.md`.

The active wrapper intentionally refuses to start unless all of its frozen inputs are present.
