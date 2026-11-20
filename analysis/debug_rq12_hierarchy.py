"""Diagnose observed-information stability for the hierarchical quadrature model."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
sys.path.insert(0, str(ROOT / "analysis"))
from rq12_core import HierarchicalLogistic

CONFIG = json.loads((ROOT / "analysis" / "rq12_analysis_config.json").read_text(encoding="utf-8"))["hierarchical_model"]
rng = np.random.default_rng(29)
groups = []
for alpha, beta in [(-1.0, 0.8), (0.0, 1.2), (0.8, 0.6)]:
    x = rng.normal(size=160)
    outcomes = rng.binomial(1, 1 / (1 + np.exp(-(alpha + beta * x))))
    groups.append((x, outcomes))
model = HierarchicalLogistic(groups, CONFIG, nodes=11)
initial = model._initial_parameters()
bounds = [(-30.0, 30.0), (-30.0, 30.0), (np.log(0.01), np.log(10.0)), (np.log(0.01), np.log(10.0))]
result = minimize(model.objective, initial, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6, "maxls": 50})
print("optimum", result.x, result.fun, result.success, result.message)


def hessian_at(step_size):
    parameters = result.x
    n = len(parameters)
    hessian = np.zeros((n, n))
    steps = np.array([step_size * max(1.0, abs(v)) for v in parameters])
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = steps[i]
        hessian[i, i] = (model.objective(parameters + ei) - 2 * model.objective(parameters) + model.objective(parameters - ei)) / steps[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = steps[j]
            hessian[i, j] = hessian[j, i] = (
                model.objective(parameters + ei + ej)
                - model.objective(parameters + ei - ej)
                - model.objective(parameters - ei + ej)
                + model.objective(parameters - ei - ej)
            ) / (4 * steps[i] * steps[j])
    return hessian


for step in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5]:
    information = -hessian_at(step)
    print("step", step, "eigenvalues", np.linalg.eigvalsh(information))
print("gradient", np.array([(model.objective(result.x + np.eye(4)[i] * 1e-5) - model.objective(result.x - np.eye(4)[i] * 1e-5)) / 2e-5 for i in range(4)]))
