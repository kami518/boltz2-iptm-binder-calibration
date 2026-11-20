"""Post-hoc exploratory RQ3-E external predictive I/F validation."""
from __future__ import annotations

import hashlib
import json
import math
import platform
import re
import sys
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq
from scipy.special import expit
from sklearn.model_selection import StratifiedGroupKFold

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
PROTOCOL = ROOT / "analysis/RQ3E_IF_external_adaptation_protocol_2026-11-14.md"
PROTOCOL_SHA = "919a8010610e8b1abd9ac4a6c20f4516e60c1156c9a85e19ebcd4682d2b54df1"
UPSTREAM = ROOT / "analysis/rq3e_outputs/run_20261114_164950"
U_FILE = UPSTREAM / "rq3e_u_pair_predictions.csv"
PAIR_FILE = UPSTREAM / "rq3e_external_pairs.csv"
RAW_FILE = ROOT / "data/external_target_audit_2026-11-13/anthropic_audited_pair_table.csv"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
OUTPUT_ROOT = ROOT / "analysis/rq3e_if_outputs"
PRIMARY_CONFIG = "primary_target_equal_weak_plus_best_ipsae"
TARGETS = ["BHRF1", "IL-7Ra", "RBX1", "TrkA", "Nipah-G", "PD-L1", "TREM2", "VEGF-A"]
NEW_TARGETS = {"BHRF1", "IL-7Ra", "RBX1", "TrkA"}
LAMBDAS = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
CLIP = 1e-6
EDGES = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
CIGAR_RE = re.compile(r"(\d+)([=XMDI])")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False)
    temp.replace(path)


def make_run_dir() -> Path:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    candidate = OUTPUT_ROOT / stamp
    suffix = 1
    while candidate.exists():
        candidate = OUTPUT_ROOT / f"{stamp}_{suffix}"
        suffix += 1
    candidate.mkdir()
    return candidate


def alignment_values(result, qlen: int, tlen: int) -> tuple[float, float, float]:
    cigar = result.cigar.decode
    if isinstance(cigar, bytes):
        cigar = cigar.decode("ascii")
    matches = qc = tc = columns = 0
    for number, operation in CIGAR_RE.findall(cigar):
        n = int(number)
        columns += n
        if operation == "=":
            matches += n
            qc += n
            tc += n
        elif operation in ("X", "M"):
            if operation == "M":
                matches += n
            qc += n
            tc += n
        elif operation == "I":
            qc += n
        elif operation == "D":
            tc += n
    if not columns or not qc or not tc:
        return 0.0, 0.0, 0.0
    return matches / columns, qc / qlen, tc / tlen


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]


def fit_i(eta: np.ndarray, y: np.ndarray) -> dict:
    eta = np.asarray(eta, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(np.unique(y)) != 2:
        return {"status": "failed", "failure_type": "single_class", "a": np.nan}
    target = float(y.mean())
    def equation(a):
        return float(np.mean(expit(a + eta)) - target)
    if equation(-40.0) * equation(40.0) > 0:
        return {"status": "failed", "failure_type": "intercept_root_not_bracketed", "a": np.nan}
    a = float(brentq(equation, -40.0, 40.0, xtol=1e-13, rtol=1e-14))
    eta_i = a + eta
    return {
        "status": "evaluable", "failure_type": "", "a": a,
        "objective": float(np.mean(np.logaddexp(0.0, eta_i) - y * eta_i)),
        "gradient_max_abs": abs(float(np.mean(expit(eta_i) - y))),
    }


def f_objective(theta, eta, y, lam):
    a, delta = theta
    z = a + eta + delta * eta
    return float(np.mean(np.logaddexp(0.0, z) - y * z) + 0.5 * lam * delta * delta)


def fit_f(eta: np.ndarray, y: np.ndarray, lam: float) -> dict:
    eta = np.asarray(eta, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(np.unique(y)) != 2:
        return {"status": "failed", "failure_type": "single_class"}
    start_i = fit_i(eta, y)
    if start_i["status"] != "evaluable":
        return {"status": "failed", "failure_type": "I_start_failed"}
    theta = np.array([start_i["a"], 0.0], dtype=float)
    current = f_objective(theta, eta, y, lam)
    converged = False
    gradient = np.array([np.inf, np.inf])
    for _ in range(150):
        a, delta = theta
        z = a + eta + delta * eta
        p = expit(z)
        residual = p - y
        w = p * (1 - p)
        gradient = np.array([
            np.mean(residual),
            np.mean(residual * eta) + lam * delta,
        ])
        hessian = np.array([
            [np.mean(w), np.mean(w * eta)],
            [np.mean(w * eta), np.mean(w * eta * eta) + lam],
        ])
        try:
            step = np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError:
            return {"status": "failed", "failure_type": "singular_hessian"}
        if not np.isfinite(step).all():
            return {"status": "failed", "failure_type": "nonfinite_step"}
        scale = 1.0
        accepted = False
        for _ in range(60):
            candidate = theta - scale * step
            value = f_objective(candidate, eta, y, lam)
            if np.isfinite(value) and value <= current + 1e-13:
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            return {"status": "failed", "failure_type": "line_search_failed"}
        change = float(np.max(np.abs(candidate - theta)))
        theta = candidate
        current = value
        if change < 1e-10 and np.max(np.abs(gradient)) < 1e-8:
            converged = True
            break
    a, delta = theta
    z = a + eta + delta * eta
    p = expit(z)
    gradient = np.array([
        np.mean(p - y),
        np.mean((p - y) * eta) + lam * delta,
    ])
    gradient_max = float(np.max(np.abs(gradient)))
    bound_hit = bool(np.any(np.abs(theta) >= 39.999999))
    if not converged or gradient_max > 1e-6 or bound_hit or not np.isfinite(theta).all():
        return {"status": "failed", "failure_type": "nonconvergence_or_bound", "a": a, "delta": delta, "gradient_max_abs": gradient_max}
    return {
        "status": "evaluable", "failure_type": "", "a": float(a), "delta": float(delta),
        "gamma": float(1 + delta), "lambda": float(lam), "objective": current,
        "gradient_max_abs": gradient_max, "iterations": _ + 1,
    }


def fit_calibration(y, p):
    from scipy.optimize import minimize

    p = np.asarray(p, dtype=float)
    y = np.asarray(y, dtype=float)
    failure = {
        "calibration_intercept": np.nan, "calibration_slope": np.nan,
        "calibration_intercept_se": np.nan, "calibration_slope_se": np.nan,
        "calibration_gradient_max_abs": np.nan, "calibration_status": "failed",
    }
    if len(np.unique(y)) != 2 or np.std(p) <= 0:
        return {**failure, "calibration_failure_type": "single_class_or_constant_prediction"}
    q = np.log(np.clip(p, CLIP, 1 - CLIP) / (1 - np.clip(p, CLIP, 1 - CLIP)))
    X = np.column_stack([np.ones(len(q)), q])
    try:
        def objective_gradient(value):
            linear = X @ value
            prob = expit(linear)
            objective = float(np.sum(np.logaddexp(0.0, linear) - y * linear))
            return objective, X.T @ (prob - y)

        # A quasi-Newton solver with line search avoids confusing an oversized
        # unmodified Newton step with genuine separation. No penalty is added.
        fit = minimize(objective_gradient, np.array([0.0, 1.0]), jac=True,
                       method="BFGS", options={"maxiter": 1000, "gtol": 1e-8})
        theta = fit.x
        prob = expit(X @ theta)
        grad = X.T @ (prob - y)
        info = X.T @ ((prob * (1 - prob))[:, None] * X)
        if (not np.isfinite(theta).all() or np.max(np.abs(grad)) > 1e-5
                or not np.isfinite(info).all() or np.linalg.cond(info) > 1e14):
            return {**failure, "calibration_failure_type": "separation_or_nonconvergence"}
        covariance = np.linalg.inv(info)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return {**failure, "calibration_failure_type": "singular_or_separation"}
    return {
        "calibration_intercept": float(theta[0]), "calibration_slope": float(theta[1]),
        "calibration_intercept_se": float(np.sqrt(covariance[0, 0])),
        "calibration_slope_se": float(np.sqrt(covariance[1, 1])),
        "calibration_gradient_max_abs": float(np.max(np.abs(grad))),
        "calibration_status": "evaluable", "calibration_failure_type": "",
    }


def metric_row(g: pd.DataFrame, layer: str, calibration=True) -> dict:
    y = g.outcome_value.to_numpy(dtype=float)
    p = g[f"p_{layer}"].to_numpy(dtype=float)
    ref_train = g.adaptation_train_reference.to_numpy(dtype=float)
    ref_oracle = g.target_oracle_reference.to_numpy(dtype=float)
    ref_transport = g.bss_transport_reference.to_numpy(dtype=float)
    brier = float(np.mean((p - y) ** 2))
    ll_p = np.clip(p, CLIP, 1 - CLIP)
    result = {
        "layer": layer, "n": len(g), "positives": int(y.sum()), "prevalence": float(y.mean()),
        "mean_prediction": float(p.mean()), "brier": brier,
        "log_loss": float(-np.mean(y * np.log(ll_p) + (1 - y) * np.log1p(-ll_p))),
        "bss_adaptation_train": float(1 - np.sum((p - y) ** 2) / np.sum((ref_train - y) ** 2)),
        "bss_target_oracle": float(1 - np.sum((p - y) ** 2) / np.sum((ref_oracle - y) ** 2)),
        "bss_transport": float(1 - np.sum((p - y) ** 2) / np.sum((ref_transport - y) ** 2)),
    }
    if calibration:
        result.update(fit_calibration(y, p))
    return result


def synthetic_tests():
    rng = np.random.default_rng(20261114)
    eta = rng.normal(size=300)
    y = rng.binomial(1, expit(-0.4 + 1.3 * eta))
    fits = [fit_f(eta, y, lam) for lam in LAMBDAS]
    assert all(f["status"] == "evaluable" for f in fits)
    assert abs(fits[-1]["delta"]) < abs(fits[0]["delta"])
    fit = fits[2]
    theta = np.array([fit["a"], fit["delta"]])
    eps = 1e-5
    fd = []
    for j in range(2):
        plus = theta.copy(); plus[j] += eps
        minus = theta.copy(); minus[j] -= eps
        fd.append((f_objective(plus, eta, y, fit["lambda"]) - f_objective(minus, eta, y, fit["lambda"])) / (2 * eps))
    z = fit["a"] + eta + fit["delta"] * eta
    p = expit(z)
    analytic = np.array([np.mean(p - y), np.mean((p - y) * eta) + fit["lambda"] * fit["delta"]])
    assert np.max(np.abs(analytic - fd)) < 1e-6
    assert fit["gamma"] == 1 + fit["delta"]
    return {"status": "passed", "finite_difference_max_abs": float(np.max(np.abs(analytic - fd))), "large_lambda_delta_abs": abs(fits[-1]["delta"]), "small_lambda_delta_abs": abs(fits[0]["delta"])}


def main():
    run = make_run_dir()
    log_path = run / "execution.log"
    def log(message):
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")
    atomic_json(run / "RUN_STATUS.json", {"status": "running", "run_dir": str(run), "analysis_role": "post_hoc_exploratory_external_IF"})
    failures = []
    try:
        log("Checking frozen protocol and input hashes.")
        expected_hashes = {
            PROTOCOL: PROTOCOL_SHA,
            U_FILE: "3b256c3fd0f813885ce817bfa6df8f80cd23ac57473047ff976ee87e0b70b2c1",
            PAIR_FILE: "72f4127d410160224701729fa401b4076297641402aa0965a3f9aec250a8002e",
            MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
        }
        for path, expected in expected_hashes.items():
            if sha256(path) != expected:
                raise RuntimeError(f"Hash mismatch: {path}")
        tests = synthetic_tests()
        atomic_json(run / "synthetic_model_tests.json", tests)

        u = pd.read_csv(U_FILE)
        u = u.loc[u.config_id.eq(PRIMARY_CONFIG)].copy()
        raw = pd.read_csv(RAW_FILE, usecols=["uuid", "sequence"])
        data = u[["uuid", "target", "stratum", "eta_U", "predicted_probability", "bss_transport_reference", "bss_target_oracle_reference", "outcome_value"]].merge(raw, on="uuid", how="left", validate="one_to_one")
        data = data.rename(columns={"predicted_probability": "p_U"})
        if len(data) != 720 or data.uuid.nunique() != 720 or set(data.target) != set(TARGETS) or data.sequence.isna().any():
            raise RuntimeError("External primary cohort mismatch")
        if not np.allclose(data.p_U, expit(data.eta_U), rtol=0, atol=1e-14):
            raise RuntimeError("Frozen U identity failed")

        log("Building within-target 90% sequence clusters before any adaptation fit.")
        sys.path.insert(0, str(ROOT / "analysis/rq12_lib"))
        import parasail
        edge_rows = []
        cluster_frames = []
        cluster_summaries = []
        for target in TARGETS:
            block = data.loc[data.target.eq(target)].reset_index(drop=True)
            sequences = block.sequence.astype(str).str.upper().tolist()
            uf = UnionFind(len(block))
            for i in range(len(block)):
                for j in range(i + 1, len(block)):
                    result = parasail.sw_trace_striped_16(sequences[i], sequences[j], 11, 1, parasail.blosum62)
                    identity, qcov, tcov = alignment_values(result, len(sequences[i]), len(sequences[j]))
                    if identity >= 0.90 and qcov >= 0.80 and tcov >= 0.80:
                        uf.union(i, j)
                        edge_rows.append({"target": target, "uuid_1": block.uuid[i], "uuid_2": block.uuid[j], "identity": identity, "query_coverage": qcov, "target_coverage": tcov})
            roots = [uf.find(i) for i in range(len(block))]
            labels = {root: f"{target}|cluster_{k+1:03d}" for k, root in enumerate(sorted(set(roots), key=lambda r: roots.index(r)))}
            block["sequence_cluster90"] = [labels[r] for r in roots]
            cluster_frames.append(block[["uuid", "target", "sequence_cluster90"]])
            sizes = block.groupby("sequence_cluster90").size()
            cluster_summaries.append({"target": target, "n_pairs": len(block), "n_clusters": len(sizes), "largest_cluster": int(sizes.max()), "singleton_clusters": int((sizes == 1).sum())})
        clusters = pd.concat(cluster_frames, ignore_index=True)
        data = data.merge(clusters, on=["uuid", "target"], how="left", validate="one_to_one")
        atomic_csv(run / "sequence_cluster90_edges.csv", pd.DataFrame(edge_rows, columns=["target", "uuid_1", "uuid_2", "identity", "query_coverage", "target_coverage"]))
        atomic_csv(run / "sequence_cluster90_assignments.csv", clusters.sort_values(["target", "uuid"]))
        atomic_csv(run / "sequence_cluster90_summary.csv", pd.DataFrame(cluster_summaries))

        log("Freezing outer grouped folds before adaptation fitting.")
        data["outer_fold"] = -1
        outer_fold_rows = []
        for target in TARGETS:
            idx = np.flatnonzero(data.target.eq(target).to_numpy())
            block = data.iloc[idx]
            splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=20261112)
            for fold, (train_pos, test_pos) in enumerate(splitter.split(np.zeros(len(block)), block.outcome_value.to_numpy(int), block.sequence_cluster90), start=1):
                data.iloc[idx[test_pos], data.columns.get_loc("outer_fold")] = fold
                outer_fold_rows.append({"target": target, "outer_fold": fold, "train_n": len(train_pos), "train_positives": int(block.outcome_value.to_numpy(int)[train_pos].sum()), "test_n": len(test_pos), "test_positives": int(block.outcome_value.to_numpy(int)[test_pos].sum())})
        if (data.outer_fold < 1).any():
            raise RuntimeError("Incomplete outer fold assignment")
        outer_assignments = data[["uuid", "target", "stratum", "sequence_cluster90", "outer_fold", "outcome_value", "eta_U", "p_U"]].copy()
        atomic_csv(run / "outer_assignments_frozen.csv", outer_assignments.sort_values(["target", "outer_fold", "uuid"]))
        atomic_csv(run / "outer_fold_manifest.csv", pd.DataFrame(outer_fold_rows))
        for (target, fold), g in data.groupby(["target", "outer_fold"]):
            train_groups = set(data.loc[data.target.eq(target) & data.outer_fold.ne(fold), "sequence_cluster90"])
            assert set(g.sequence_cluster90).isdisjoint(train_groups)

        log("Running target-internal predictive I/F with nested lambda selection.")
        inner_assignment_rows = []
        inner_prediction_rows = []
        lambda_rows = []
        fold_parameter_rows = []
        oof_frames = []
        fold_metric_rows = []
        for target in TARGETS:
            target_data = data.loc[data.target.eq(target)].copy()
            for fold in range(1, 6):
                train = target_data.loc[target_data.outer_fold.ne(fold)].copy()
                test = target_data.loc[target_data.outer_fold.eq(fold)].copy()
                y_train = train.outcome_value.to_numpy(dtype=float)
                eta_train = train.eta_U.to_numpy(dtype=float)
                fit_i_outer = fit_i(eta_train, y_train)
                fold_parameter_rows.append({"target": target, "outer_fold": fold, "layer": "I", "status": fit_i_outer["status"], "failure_type": fit_i_outer["failure_type"], "a": fit_i_outer["a"], "delta": 0.0, "gamma": 1.0, "lambda": np.nan, "objective": fit_i_outer.get("objective", np.nan), "gradient_max_abs": fit_i_outer.get("gradient_max_abs", np.nan), "train_n": len(train), "train_positives": int(y_train.sum()), "test_n": len(test), "test_positives": int(test.outcome_value.sum())})
                if fit_i_outer["status"] != "evaluable":
                    raise RuntimeError(f"Unexpected I failure under frozen complete-case design: {target} fold {fold}")
                test["a_I"] = fit_i_outer["a"]
                test["p_I"] = expit(fit_i_outer["a"] + test.eta_U.to_numpy(float))

                inner = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=20261113)
                inner_splits = list(inner.split(np.zeros(len(train)), y_train.astype(int), train.sequence_cluster90))
                for inner_fold, (itr, iva) in enumerate(inner_splits, start=1):
                    for row_index in iva:
                        inner_assignment_rows.append({"target": target, "outer_fold": fold, "inner_fold": inner_fold, "uuid": train.uuid.iloc[row_index], "sequence_cluster90": train.sequence_cluster90.iloc[row_index]})
                candidate_results = {}
                for lam in LAMBDAS:
                    candidate_ok = True
                    candidate_message = ""
                    val_frames = []
                    for inner_fold, (itr, iva) in enumerate(inner_splits, start=1):
                        fit = fit_f(eta_train[itr], y_train[itr], lam)
                        if fit["status"] != "evaluable":
                            candidate_ok = False
                            candidate_message = fit["failure_type"]
                            failures.append({"stage": "inner_F", "item": f"{target}|{fold}|{inner_fold}|{lam}", "failure_type": fit["failure_type"]})
                            break
                        p_val = expit(fit["a"] + fit["gamma"] * eta_train[iva])
                        val_frames.append(pd.DataFrame({"target": target, "outer_fold": fold, "inner_fold": inner_fold, "uuid": train.uuid.iloc[iva].to_numpy(), "lambda": lam, "outcome_value": y_train[iva], "predicted_probability": p_val, "a": fit["a"], "delta": fit["delta"], "gamma": fit["gamma"], "gradient_max_abs": fit["gradient_max_abs"]}))
                    if candidate_ok:
                        vals = pd.concat(val_frames, ignore_index=True)
                        inner_prediction_rows.append(vals)
                        p = np.clip(vals.predicted_probability.to_numpy(float), CLIP, 1 - CLIP)
                        yy = vals.outcome_value.to_numpy(float)
                        pooled_ll = float(-np.mean(yy * np.log(p) + (1 - yy) * np.log1p(-p)))
                        lambda_rows.append({"target": target, "outer_fold": fold, "lambda": lam, "candidate_status": "evaluable", "failure_type": "", "pooled_inner_log_loss": pooled_ll, "n_inner_validation": len(vals), "selected": False})
                        candidate_results[lam] = pooled_ll
                    else:
                        lambda_rows.append({"target": target, "outer_fold": fold, "lambda": lam, "candidate_status": "failed", "failure_type": candidate_message, "pooled_inner_log_loss": np.nan, "n_inner_validation": 0, "selected": False})
                if not candidate_results:
                    raise RuntimeError(f"All F candidates failed: {target} fold {fold}")
                minimum = min(candidate_results.values())
                selected_lambda = max(lam for lam, value in candidate_results.items() if abs(value - minimum) <= 1e-12)
                for row in lambda_rows:
                    if row["target"] == target and row["outer_fold"] == fold and row["lambda"] == selected_lambda:
                        row["selected"] = True
                final_f = fit_f(eta_train, y_train, selected_lambda)
                fold_parameter_rows.append({"target": target, "outer_fold": fold, "layer": "F", "status": final_f["status"], "failure_type": final_f.get("failure_type", ""), "a": final_f.get("a", np.nan), "delta": final_f.get("delta", np.nan), "gamma": final_f.get("gamma", np.nan), "lambda": selected_lambda, "objective": final_f.get("objective", np.nan), "gradient_max_abs": final_f.get("gradient_max_abs", np.nan), "train_n": len(train), "train_positives": int(y_train.sum()), "test_n": len(test), "test_positives": int(test.outcome_value.sum())})
                if final_f["status"] != "evaluable":
                    raise RuntimeError(f"Unexpected final F failure under frozen complete-case design: {target} fold {fold}")
                test["a_F"] = final_f["a"]
                test["delta_F"] = final_f["delta"]
                test["gamma_F"] = final_f["gamma"]
                test["selected_lambda"] = selected_lambda
                test["p_F"] = expit(final_f["a"] + final_f["gamma"] * test.eta_U.to_numpy(float))
                test["adaptation_train_reference"] = float(y_train.mean())
                test["target_oracle_reference"] = float(target_data.outcome_value.mean())
                oof_frames.append(test)
                for layer in ["U", "I", "F"]:
                    fold_metric_rows.append({"target": target, "outer_fold": fold, **metric_row(test, layer, calibration=False)})

        atomic_csv(run / "inner_assignments.csv", pd.DataFrame(inner_assignment_rows))
        atomic_csv(run / "inner_validation_predictions.csv", pd.concat(inner_prediction_rows, ignore_index=True))
        atomic_csv(run / "lambda_selection.csv", pd.DataFrame(lambda_rows))
        atomic_csv(run / "fold_parameters.csv", pd.DataFrame(fold_parameter_rows))
        oof = pd.concat(oof_frames, ignore_index=True).sort_values(["target", "uuid"])
        if len(oof) != 720 or oof.uuid.nunique() != 720:
            raise RuntimeError("OOF output is not one complete row per external pair")
        oof_columns = ["uuid", "target", "stratum", "sequence_cluster90", "outer_fold", "outcome_value", "eta_U", "p_U", "a_I", "p_I", "a_F", "delta_F", "gamma_F", "selected_lambda", "p_F", "adaptation_train_reference", "target_oracle_reference", "bss_transport_reference", "bss_target_oracle_reference"]
        atomic_csv(run / "oof_predictions.csv", oof[oof_columns])
        atomic_csv(run / "fold_metrics.csv", pd.DataFrame(fold_metric_rows))

        log("Evaluating frozen out-of-fold predictions.")
        target_rows = []
        reliability_rows = []
        comparison_rows = []
        for target in TARGETS:
            g = oof.loc[oof.target.eq(target)]
            target_result = {"target": target, "stratum": g.stratum.iloc[0]}
            for layer in ["U", "I", "F"]:
                target_result.update(metric_row(g, layer))
                target_rows.append({"target": target, "stratum": g.stratum.iloc[0], **metric_row(g, layer)})
                p = g[f"p_{layer}"].to_numpy(float)
                y = g.outcome_value.to_numpy(int)
                bins = np.digitize(p, EDGES[1:-1], right=False)
                for b, (lower, upper) in enumerate(zip(EDGES[:-1], EDGES[1:])):
                    take = bins == b
                    reliability_rows.append({"target": target, "layer": layer, "bin_index": b, "bin_lower": lower, "bin_upper": upper, "n": int(take.sum()), "positives": int(y[take].sum()), "mean_prediction": float(p[take].mean()) if take.any() else np.nan, "observed_rate": float(y[take].mean()) if take.any() else np.nan})
            layer_metrics = {layer: metric_row(g, layer, calibration=False) for layer in ["U", "I", "F"]}
            for comparison, left, right in [("I_minus_U", "I", "U"), ("F_minus_U", "F", "U"), ("F_minus_I", "F", "I")]:
                comparison_rows.append({"target": target, "stratum": g.stratum.iloc[0], "comparison": comparison, "brier_difference": layer_metrics[left]["brier"] - layer_metrics[right]["brier"], "log_loss_difference": layer_metrics[left]["log_loss"] - layer_metrics[right]["log_loss"]})
        target_metrics = pd.DataFrame(target_rows)
        for row in target_metrics.loc[target_metrics.calibration_status.eq("failed")].itertuples(index=False):
            failures.append({"stage": "OOF_calibration_diagnostic", "item": f"{row.target}|{row.layer}", "failure_type": row.calibration_failure_type})
        atomic_csv(run / "target_metrics.csv", target_metrics)
        atomic_csv(run / "reliability.csv", pd.DataFrame(reliability_rows))
        comparisons = pd.DataFrame(comparison_rows)
        atomic_csv(run / "target_comparisons.csv", comparisons)

        group_rows = []
        for stratum, targets in [("all_8", set(TARGETS)), ("new_target_4", NEW_TARGETS), ("same_target_campaign_4", set(TARGETS) - NEW_TARGETS)]:
            for layer in ["U", "I", "F"]:
                subset = target_metrics.loc[target_metrics.layer.eq(layer) & target_metrics.target.isin(targets)]
                row = {"aggregation": "target_equal_macro", "stratum": stratum, "layer": layer, "n_targets": len(subset), "n_pairs": int(subset.n.sum()), "positives": int(subset.positives.sum())}
                for name in ["prevalence", "mean_prediction", "brier", "log_loss", "bss_adaptation_train", "bss_target_oracle", "bss_transport"]:
                    row[name] = float(subset[name].mean())
                    row[f"{name}_min"] = float(subset[name].min())
                    row[f"{name}_max"] = float(subset[name].max())
                group_rows.append(row)
            for comparison, left, right in [("I_minus_U", "I", "U"), ("F_minus_U", "F", "U"), ("F_minus_I", "F", "I")]:
                subset = comparisons.loc[comparisons.comparison.eq(comparison) & comparisons.target.isin(targets)]
                group_rows.append({"aggregation": "target_equal_macro_comparison", "stratum": stratum, "layer": comparison, "n_targets": len(subset), "n_pairs": int(target_metrics.loc[target_metrics.layer.eq("U") & target_metrics.target.isin(targets), "n"].sum()), "positives": int(target_metrics.loc[target_metrics.layer.eq("U") & target_metrics.target.isin(targets), "positives"].sum()), "brier": float(subset.brier_difference.mean()), "log_loss": float(subset.log_loss_difference.mean()), "brier_min": float(subset.brier_difference.min()), "brier_max": float(subset.brier_difference.max()), "log_loss_min": float(subset.log_loss_difference.min()), "log_loss_max": float(subset.log_loss_difference.max()), "brier_improved_targets": int(subset.brier_difference.lt(0).sum()), "log_loss_improved_targets": int(subset.log_loss_difference.lt(0).sum())})
            pooled = oof.loc[oof.target.isin(targets)]
            for layer in ["U", "I", "F"]:
                group_rows.append({"aggregation": "pair_pooled", "stratum": stratum, **metric_row(pooled, layer)})
        group_metrics = pd.DataFrame(group_rows)
        atomic_csv(run / "group_metrics.csv", group_metrics)

        log("Drawing OOF adaptation figures.")
        order = TARGETS
        fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
        x = np.arange(8); width = 0.25
        for k, layer in enumerate(["U", "I", "F"]):
            values = target_metrics.loc[target_metrics.layer.eq(layer)].set_index("target").loc[order, "brier"]
            axes[0].bar(x + (k - 1) * width, values, width, label=layer)
        axes[0].set_xticks(x, order, rotation=35, ha="right"); axes[0].set_ylabel("OOF Brier"); axes[0].legend(title="Layer")
        for comparison, color in [("I_minus_U", "#1b9e77"), ("F_minus_I", "#7570b3")]:
            values = comparisons.loc[comparisons.comparison.eq(comparison)].set_index("target").loc[order, "brier_difference"]
            axes[1].bar(x + (-0.18 if comparison == "I_minus_U" else 0.18), values, 0.36, label=comparison, color=color)
        axes[1].axhline(0, color="black", linewidth=0.8); axes[1].set_xticks(x, order, rotation=35, ha="right"); axes[1].set_ylabel("OOF Brier difference"); axes[1].legend()
        fig.savefig(run / "RQ3E_IF_OOF_metrics.png", dpi=220)
        plt.close(fig)
        fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
        for ax, target in zip(axes.ravel(), TARGETS):
            reliability = pd.DataFrame(reliability_rows)
            for layer, color in [("U", "#4c78a8"), ("I", "#1b9e77"), ("F", "#7570b3")]:
                r = reliability.loc[reliability.target.eq(target) & reliability.layer.eq(layer) & reliability.n.gt(0)]
                ax.plot(r.mean_prediction, r.observed_rate, marker="o", label=layer, color=color)
            ax.plot([0, 1], [0, 1], "--", color="0.7"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title(target)
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.savefig(run / "RQ3E_IF_OOF_reliability.png", dpi=220)
        plt.close(fig)

        log("Writing result summary, QC and final hashes.")
        primary_macro = group_metrics.loc[group_metrics.aggregation.eq("target_equal_macro") & group_metrics.stratum.eq("all_8")].set_index("layer")
        comparison_macro = group_metrics.loc[group_metrics.aggregation.eq("target_equal_macro_comparison") & group_metrics.stratum.eq("all_8")].set_index("layer")
        lines = [
            "# RQ3‑E‑IF 外部有标签适配结果", "",
            "**定位：看到零标签外部 U 结果后新增的探索性分析；不改变确认性 U 结果。**", "",
            "主分析：八靶点／720 行，Adaptyv Weak+，冻结 target-equal U；靶点内序列簇 5 折，F 使用训练折内 3 折选择 lambda。", "",
            "## 八靶点等权结果", "",
            "| 层级 | Brier | Log loss | BSS-adaptation-train | BSS-target-oracle | BSS-transport |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for layer in ["U", "I", "F"]:
            r = primary_macro.loc[layer]
            lines.append(f"| {layer} | {r.brier:.4f} | {r.log_loss:.4f} | {r.bss_adaptation_train:.4f} | {r.bss_target_oracle:.4f} | {r.bss_transport:.4f} |")
        lines += ["", "## 跨层比较", "", "| 比较 | Brier 宏平均差 | Log loss 宏平均差 | Brier 改善靶点 | Log loss 改善靶点 |", "|---|---:|---:|---:|---:|"]
        for comparison in ["I_minus_U", "F_minus_U", "F_minus_I"]:
            r = comparison_macro.loc[comparison]
            lines.append(f"| {comparison} | {r.brier:.4f} | {r.log_loss:.4f} | {int(r.brier_improved_targets)}／8 | {int(r.log_loss_improved_targets)}／8 |")
        lines += ["", "## 逐靶点结果", "", "| 靶点 | U Brier | I Brier | F Brier | I−U | F−I | I BSS-train | F BSS-train |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
        for target in TARGETS:
            m = target_metrics.loc[target_metrics.target.eq(target)].set_index("layer")
            c = comparisons.loc[comparisons.target.eq(target)].set_index("comparison")
            lines.append(f"| {target} | {m.loc['U','brier']:.4f} | {m.loc['I','brier']:.4f} | {m.loc['F','brier']:.4f} | {c.loc['I_minus_U','brier_difference']:.4f} | {c.loc['F_minus_I','brier_difference']:.4f} | {m.loc['I','bss_adaptation_train']:.4f} | {m.loc['F','bss_adaptation_train']:.4f} |")
        lambda_counts = oof.groupby("selected_lambda").uuid.nunique().to_dict()
        lines += ["", "## 解释边界", "", "- 所有 I／F 数值均为靶点内折外预测；外部测试折标签没有参与拟合或 lambda 选择。", "- I／F 使用本地标签，不能挽救或改写零标签 U 的确认性结论。", "- 每靶点仅 90 行，F 的斜率适配仍可能高方差；八靶点方向计数不是未来靶点总体推断。", f"- 外测行选定 lambda 的配对计数：{lambda_counts}。", "- 未运行 Pinc、其他结局、其他 score 或其他 ML。"]
        (run / "RQ3E_IF_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        atomic_csv(run / "failed_item_ledger.csv", pd.DataFrame(failures, columns=["stage", "item", "failure_type"]))
        qc = {
            "status": "passed", "analysis_role": "post_hoc_exploratory_external_IF", "protocol_sha256": sha256(PROTOCOL),
            "external_pairs": len(oof), "unique_pairs": int(oof.uuid.nunique()), "targets": int(oof.target.nunique()),
            "outer_folds": 40, "lambda_candidate_rows": len(lambda_rows), "inner_prediction_rows": sum(len(x) for x in inner_prediction_rows),
            "all_U_I_F_probabilities_finite": bool(np.isfinite(oof[["p_U", "p_I", "p_F"]].to_numpy()).all()),
            "outer_group_leakage_detected": False, "outer_test_labels_used_for_fitting_or_selection": False,
            "U_refit": False, "Pinc_or_ML_run": False, "failed_items": len(failures),
            "master_sha256_after": sha256(MASTER), "result_driven_model_changes": False,
            "calibration_SE_definition": "inverse unnormalised Hessian X'WX",
        }
        atomic_json(run / "qc_checks.json", qc)
        config = {
            "protocol": str(PROTOCOL), "protocol_sha256": sha256(PROTOCOL), "analysis_role": "post_hoc_exploratory_external_IF",
            "input_hashes": {str(p): sha256(p) for p in expected_hashes}, "primary_config": PRIMARY_CONFIG,
            "targets": TARGETS, "new_targets": sorted(NEW_TARGETS), "lambda_grid": LAMBDAS,
            "outer_cv": {"algorithm": "StratifiedGroupKFold", "n_splits": 5, "shuffle": True, "random_state": 20261112, "groups": "within-target sequence_cluster90"},
            "inner_cv": {"algorithm": "StratifiedGroupKFold", "n_splits": 3, "shuffle": True, "random_state": 20261113, "groups": "within-target sequence_cluster90"},
            "I_model": "eta_I=a+eta_U", "F_model": "eta_F=a+(1+delta)*eta_U; mean NLL + 0.5*lambda*delta^2",
            "lambda_selection": "minimum pooled inner log loss; abs difference <=1e-12 ties choose larger lambda",
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "parasail": parasail.__version__},
            "execution_script_sha256": sha256(Path(__file__)),
        }
        atomic_json(run / "run_config.snapshot.json", config)
        atomic_json(run / "RUN_STATUS.json", {"status": "complete", "qc": "passed", "run_dir": str(run), "analysis_role": "post_hoc_exploratory_external_IF", "external_pairs": 720, "outer_folds": 40})
        log(f"Complete with QC passed: {run}")
        files = sorted(p for p in run.iterdir() if p.is_file() and p.name != "final_hash_manifest.csv")
        atomic_csv(run / "final_hash_manifest.csv", pd.DataFrame([{"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files]))
    except Exception as error:
        atomic_json(run / "RUN_STATUS.json", {"status": "failed", "run_dir": str(run), "error": str(error), "traceback": traceback.format_exc()})
        log(f"FAILED: {error}")
        raise


if __name__ == "__main__":
    main()
