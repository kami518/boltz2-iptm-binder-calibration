# RQ3 execution commands

All commands were run with working directory:

`<REPRODUCTION_ROOT>/analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311`

```powershell
python stage_0_input_and_fold_qc.py 2>&1 | Tee-Object -FilePath stage_0.log
python stage_1_u_layer.py 2>&1 | Tee-Object -FilePath stage_1.log
python stage_2_oracle_descriptive.py 2>&1 | Tee-Object -FilePath stage_2.log
python stage_3_predictive_I_F.py 2>&1 | Tee-Object -FilePath stage_3.log
python stage_4_finalize.py 2>&1 | Tee-Object -FilePath stage_4.log
python finalize_hash_manifest.py
```

`stage_4_finalize.py` invokes the complete F target-to-prediction regression audit and writes its JSON/CSV artifacts before final QC. The final hash manifest is run only after stage 4 exits and its log is closed.

No dependency installation or upgrade command was run.
