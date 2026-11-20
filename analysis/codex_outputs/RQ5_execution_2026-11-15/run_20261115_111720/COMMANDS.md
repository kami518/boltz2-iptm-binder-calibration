# RQ5 execution commands

Run directory:

```text
<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ5_execution_2026-11-15\run_20261115_111720
```

Input and fold preflight:

```powershell
python -m py_compile rq5_common.py stage_0_preflight.py stage_a_diagnostics.py
python stage_0_preflight.py 2>&1 | Tee-Object -FilePath stage_0.log
```

Stage A:

```powershell
python stage_a_diagnostics.py 2>&1 | Tee-Object -FilePath stage_a.log
```

Stage B:

```powershell
python stage_b_transport_adaptation.py 2>&1 | Tee-Object -FilePath stage_b.log
```

Required unit tests:

```powershell
python stage_c_unit_tests.py 2>&1 | Tee-Object -FilePath stage_c.log
```

Independent verification:

```powershell
python independent_verification.py 2>&1 | Tee-Object -FilePath independent_verification.log
```

Finalization:

```powershell
python stage_d_finalize.py 2>&1 | Tee-Object -FilePath stage_d.log
```

During Stage B development, two outer combo fits first stopped just above the
frozen gradient gate. The same objective was then refined by deterministic
damped Newton steps. The final run completed all 102 source-model cells and all
36,534 adaptation prediction rows. The failed trial entries remain in the
failure ledger with resolution labels.
