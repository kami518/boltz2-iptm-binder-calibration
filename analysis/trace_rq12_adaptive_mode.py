"""Trace Newton updates for the observed adaptive-mode numerical stall."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.special import expit

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
sys.path.insert(0, str(ROOT / "analysis"))
from rq12_core import clipped_logit

part = pd.read_csv(ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv").query("target == 'sars-cov-2-rbd'")
x = clipped_logit(part.boltz2_iptm, 1e-6)
y = part.primary_outcome.to_numpy(dtype=float)
mu_a, mu_b, log_tau_a, log_tau_b = [-3.9124716440567715, 5.391344178957667, 2.302585092994046, -3.1290876944826094]
mu = np.array([mu_a, mu_b])
theta = mu.copy()
precision = np.exp(-2 * np.array([log_tau_a, log_tau_b]))
history = []
for iteration in range(100):
    eta = theta[0] + theta[1] * x
    p = expit(eta)
    residual = y - p
    gradient = np.array([residual.sum(), (residual * x).sum()]) - precision * (theta - mu)
    w = p * (1 - p)
    info = np.array([[w.sum() + precision[0], (w * x).sum()], [(w * x).sum(), (w * x * x).sum() + precision[1]]])
    step = np.linalg.solve(info, gradient)
    logpost = float((eta * y).sum() - np.logaddexp(0, eta).sum() - 0.5 * np.sum(precision * (theta - mu) ** 2))
    history.append((iteration, logpost, np.linalg.norm(gradient), np.linalg.norm(step), theta.copy()))
    if np.linalg.norm(step) < 1e-9 * (1 + np.linalg.norm(theta)):
        break
    scale = 1.0
    improved = False
    for _ in range(60):
        candidate = theta + scale * step
        eta_new = candidate[0] + candidate[1] * x
        lp = float((eta_new * y).sum() - np.logaddexp(0, eta_new).sum() - 0.5 * np.sum(precision * (candidate - mu) ** 2))
        if lp >= logpost:
            theta = candidate
            improved = True
            break
        scale *= 0.5
    if not improved:
        print(f"line_search_failed iteration={iteration} last_scale={scale}")
        break
for row in history[-15:]:
    print(f"iteration={row[0]} logpost={row[1]:.15f} gradient={row[2]:.12g} step={row[3]:.12g} theta={row[4].tolist()}")
print(f"information={info.tolist()}")
