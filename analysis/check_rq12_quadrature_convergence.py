"""Check higher-order Gauss-Hermite convergence for the primary-six hierarchical model."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
sys.path.insert(0, str(ROOT / "analysis"))
from rq12_core import HierarchicalLogistic, clipped_logit

CONFIG = json.loads((ROOT / "analysis" / "rq12_analysis_config.json").read_text(encoding="utf-8"))
HCONFIG = CONFIG["hierarchical_model"]
OUT = ROOT / "analysis" / "rq12_outputs"
TARGETS = ["aah-i", "human-vegf110", "nipah-glycoprotein-g", "pd-l1", "sars-cov-2-rbd", "trem2"]
df = pd.read_csv(ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv")
groups = []
for target in TARGETS:
    selected = df[df.target == target]
    groups.append((clipped_logit(selected.boltz2_iptm, CONFIG["probability_clip"]), selected.primary_outcome.to_numpy(dtype=int)))
previous = json.loads((OUT / "rq2_hierarchical_full_results.json").read_text(encoding="utf-8"))
initial = np.asarray(previous["models"]["primary_six"]["parameters"], dtype=float)
results = {}
started = time.time()
for nodes in [31, 51, 71]:
    model = HierarchicalLogistic(groups, HCONFIG, nodes=nodes)
    bounds = [(-30.0, 30.0), (-30.0, 30.0), (np.log(0.01), np.log(10.0)), (np.log(0.01), np.log(10.0))]
    optimum = minimize(model.objective, initial, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6, "maxls": 50})
    covariance, positive_definite, eigenvalues = model.observed_information(model.objective, optimum.x)
    targets = model._posterior_targets(optimum.x)
    result = model._format_result(optimum, covariance, positive_definite, eigenvalues, targets)
    results[str(nodes)] = result
    initial = optimum.x
    print(f"nodes={nodes} success={optimum.success} pd={positive_definite} params={optimum.x.tolist()} elapsed={time.time()-started:.1f}s", flush=True)
    (OUT / "rq12_quadrature_convergence_diagnostic.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

reference = results["71"]
for nodes in ["31", "51", "71"]:
    current = results[nodes]
    differences = {
        key: abs(current["hyperparameters"][key] - reference["hyperparameters"][key])
        for key in ["mu_intercept", "mu_slope", "tau_intercept", "tau_slope"]
    }
    target_difference = max(
        abs(a[field] - b[field])
        for a, b in zip(current["targets"], reference["targets"])
        for field in ["alpha_partial", "beta_partial"]
    )
    print(f"comparison nodes={nodes} vs71 hyper={differences} max_target={target_difference}", flush=True)
