# RQ5-RF execution commands

Working directory: `<REPRODUCTION_ROOT>`

The executed sequence was:

```powershell
python scripts/stage0_preflight.py
python scripts/test_rq5rf.py
python scripts/stage1_shared_rows.py
python scripts/stage2_u_models.py
python scripts/stage3_adaptation.py
python scripts/stage4_metrics_figures.py
python scripts/independent_verification.py
python scripts/stage5_qc_finalize.py
```

RF parallel determinism was tested. `n_jobs=12` produced non-elementwise-identical
probabilities at floating-point roundoff, so the frozen recovery rule was applied:
the tested serial setting `n_jobs=1` was used. No network, download, dependency
installation, upgrade, or cross-version result splicing was used.
