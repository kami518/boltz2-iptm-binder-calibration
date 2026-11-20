"""Diagnose conditional-mode convergence without repeating RQ1/RQ2 bootstraps."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
sys.path.insert(0, str(ROOT / "analysis"))
from rq12_adaptive import AdaptiveHierarchicalLogistic
from rq12_core import clipped_logit

config = json.loads((ROOT / "analysis" / "rq12_analysis_config.json").read_text(encoding="utf-8"))["hierarchical_model"]
df = pd.read_csv(ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv")
target_names = df.target.drop_duplicates().tolist()
groups = [(clipped_logit(part.boltz2_iptm, 1e-6), part.primary_outcome.to_numpy(dtype=int)) for _, part in df.groupby("target", sort=False)]
model = AdaptiveHierarchicalLogistic(groups, config)
original = model.objective
calls = 0

def wrapped(parameters):
    global calls
    calls += 1
    try:
        return original(parameters)
    except Exception:
        print(f"failure_at_call={calls} parameters={parameters.tolist()}", flush=True)
        for target, (x, y) in zip(target_names, groups):
            try:
                mode, information = model._mode_2d(x, y, *parameters)
                print(f"target={target} mode={mode.tolist()} information={information.tolist()}", flush=True)
            except Exception as error:
                print(f"target={target} error={error}", flush=True)
        raise

model.objective = wrapped
result = model.fit()
print(json.dumps(result, indent=2), flush=True)
