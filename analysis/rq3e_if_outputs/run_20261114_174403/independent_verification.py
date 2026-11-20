"""Independent reconstruction and numerical audit of the frozen external I/F run.

This file deliberately does not import the production analysis script.  Sequence
edges, connected components, folds, model fits, predictions and summaries are
reconstructed from the frozen upstream inputs.  Existing outputs are read only
for comparison, never used as fitting targets.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import expit
from sklearn.model_selection import StratifiedGroupKFold

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = ROOT / "analysis/rq3e_if_outputs/run_20261114_174403"
UPSTREAM = ROOT / "analysis/rq3e_outputs/run_20261114_164950"
U_FILE = UPSTREAM / "rq3e_u_pair_predictions.csv"
RAW_FILE = ROOT / "data/external_target_audit_2026-11-13/anthropic_audited_pair_table.csv"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
PROTOCOL = ROOT / "analysis/RQ3E_IF_external_adaptation_protocol_2026-11-14.md"
CONFIG = "primary_target_equal_weak_plus_best_ipsae"
TARGETS = ["BHRF1", "IL-7Ra", "RBX1", "TrkA", "Nipah-G", "PD-L1", "TREM2", "VEGF-A"]
NEW_TARGETS = {"BHRF1", "IL-7Ra", "RBX1", "TrkA"}
LAMBDAS = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
CLIP = 1e-6
EDGES = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
REPORT = {"status": "running", "run": str(RUN), "checks": {}, "maxima": {}}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def record(name: str, **details) -> None:
    REPORT["checks"][name] = {"status": "passed", **details}
    print(f"PASS {name}: {details}", flush=True)


def maximum(name: str, value: float) -> None:
    value = float(value)
    REPORT["maxima"][name] = value
    if not np.isfinite(value):
        raise AssertionError(f"Nonfinite audit maximum: {name}")


def close(actual, expected, label: str, atol=1e-9, rtol=1e-10) -> None:
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    if a.shape != e.shape:
        raise AssertionError(f"{label}: shape {a.shape} != {e.shape}")
    finite = np.isfinite(a) & np.isfinite(e)
    both_nan = np.isnan(a) & np.isnan(e)
    if not np.all(finite | both_nan):
        raise AssertionError(f"{label}: nonfinite-pattern mismatch")
    difference = np.abs(a[finite] - e[finite])
    err = float(difference.max()) if difference.size else 0.0
    maximum(label, err)
    if not np.allclose(a, e, atol=atol, rtol=rtol, equal_nan=True):
        raise AssertionError(f"{label}: numerical mismatch, maximum absolute error {err}")


def compare_frames(name: str, expected: pd.DataFrame, saved: pd.DataFrame,
                   keys: list[str], atol=1e-9) -> None:
    if list(expected.columns) != list(saved.columns):
        raise AssertionError(f"{name}: columns differ")
    if len(expected) != len(saved):
        raise AssertionError(f"{name}: row counts {len(expected)} != {len(saved)}")
    left = expected.sort_values(keys).reset_index(drop=True)
    right = saved.sort_values(keys).reset_index(drop=True)
    for key in keys:
        if not left[key].equals(right[key]):
            raise AssertionError(f"{name}: key mismatch in {key}")
    max_error = 0.0
    for column in expected.columns:
        if column in keys:
            continue
        a = left[column]
        b = right[column]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            av = a.to_numpy(dtype=float)
            bv = b.to_numpy(dtype=float)
            finite = np.isfinite(av) & np.isfinite(bv)
            if finite.any():
                max_error = max(max_error, float(np.max(np.abs(av[finite] - bv[finite]))))
            if not np.allclose(av, bv, atol=atol, rtol=1e-10, equal_nan=True):
                bad = np.flatnonzero(~np.isclose(av, bv, atol=atol, rtol=1e-10, equal_nan=True))
                raise AssertionError(f"{name}: mismatch in {column}, first row {bad[0]}")
        else:
            # CSV reads empty failure_type cells back as NA; both encode no failure.
            av = a.fillna("").astype(str)
            bv = b.fillna("").astype(str)
            if not av.equals(bv):
                raise AssertionError(f"{name}: categorical mismatch in {column}")
    maximum(f"{name}_max_abs", max_error)
    record(name, rows=len(saved), columns=len(saved.columns), max_abs_error=max_error)


def alignment_metrics(result, query_length: int, target_length: int) -> tuple[float, float, float]:
    """Parse local-alignment CIGAR independently of the production function."""
    cigar = result.cigar.decode
    if isinstance(cigar, bytes):
        cigar = cigar.decode("ascii")
    exact = query_aligned = target_aligned = total = 0
    for token in re.finditer(r"([0-9]+)([=XMDI])", cigar):
        count = int(token.group(1))
        operation = token.group(2)
        total += count
        if operation in "=XM":
            query_aligned += count
            target_aligned += count
            if operation in "=M":
                exact += count
        elif operation == "I":
            query_aligned += count
        elif operation == "D":
            target_aligned += count
    if min(total, query_aligned, target_aligned) == 0:
        return 0.0, 0.0, 0.0
    return exact / total, query_aligned / query_length, target_aligned / target_length


def components_from_edges(n: int, edges: list[tuple[int, int]]) -> list[int]:
    adjacency = [set() for _ in range(n)]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    labels = [-1] * n
    component = 0
    for start in range(n):
        if labels[start] >= 0:
            continue
        stack = [start]
        labels[start] = component
        while stack:
            node = stack.pop()
            for neighbor in adjacency[node]:
                if labels[neighbor] < 0:
                    labels[neighbor] = component
                    stack.append(neighbor)
        component += 1
    return labels


def objective_gradient(theta: np.ndarray, eta: np.ndarray, y: np.ndarray,
                       lam: float) -> tuple[float, np.ndarray]:
    a, delta = theta
    z = a + (1.0 + delta) * eta
    p = expit(z)
    value = float(np.mean(np.logaddexp(0.0, z) - y * z) + 0.5 * lam * delta ** 2)
    gradient = np.array([
        np.mean(p - y),
        np.mean((p - y) * eta) + lam * delta,
    ])
    return value, gradient


def fit_f_independent(eta: np.ndarray, y: np.ndarray, lam: float) -> dict:
    if np.unique(y).size != 2:
        raise AssertionError("Unexpected single-class F training set")
    eta = np.asarray(eta, dtype=float)
    y = np.asarray(y, dtype=float)
    prevalence = float(y.mean())
    a0 = brentq(lambda a: float(np.mean(expit(a + eta)) - prevalence), -40, 40,
                xtol=1e-13, rtol=1e-14)
    result = minimize(
        objective_gradient, np.array([a0, 0.0]), args=(eta, y, lam), jac=True,
        method="L-BFGS-B", bounds=[(-40.0, 40.0), (-40.0, 40.0)],
        options={"maxiter": 2000, "maxls": 100, "ftol": 1e-15, "gtol": 1e-11},
    )
    a, delta = result.x
    value, gradient = objective_gradient(result.x, eta, y, lam)
    if not np.isfinite(value) or np.max(np.abs(gradient)) > 2e-7:
        raise AssertionError(f"Independent F fit failed: {result.message}; gradient={gradient}")
    return {"a": float(a), "delta": float(delta), "gamma": float(1 + delta),
            "objective": value, "gradient": gradient, "optimizer_success": bool(result.success)}


def calibration_independent(y: np.ndarray, p: np.ndarray) -> dict:
    failure = {
        "calibration_intercept": np.nan, "calibration_slope": np.nan,
        "calibration_intercept_se": np.nan, "calibration_slope_se": np.nan,
        "calibration_gradient_max_abs": np.nan, "calibration_status": "failed",
    }
    if np.unique(y).size != 2 or np.std(p) <= 0:
        return {**failure, "calibration_failure_type": "single_class_or_constant_prediction"}
    clipped = np.clip(p, CLIP, 1 - CLIP)
    q = np.log(clipped / (1 - clipped))
    X = np.column_stack([np.ones(len(q)), q])

    def fun(theta):
        z = X @ theta
        prob = expit(z)
        return float(np.sum(np.logaddexp(0.0, z) - y * z)), X.T @ (prob - y)

    result = minimize(fun, np.array([0.0, 1.0]), jac=True, method="BFGS",
                      options={"maxiter": 1000, "gtol": 1e-8})
    theta = result.x
    z = X @ theta
    prob = expit(z)
    gradient = X.T @ (prob - y)
    hessian = X.T @ ((prob * (1 - prob))[:, None] * X)
    # Complete separation must remain a failed diagnostic; never use a pseudo-inverse.
    if (not np.isfinite(theta).all() or np.max(np.abs(gradient)) > 1e-5
            or np.linalg.cond(hessian) > 1e14):
        return {**failure, "calibration_failure_type": "singular_or_separation"}
    try:
        covariance = np.linalg.inv(hessian)
    except np.linalg.LinAlgError:
        return {**failure, "calibration_failure_type": "singular_or_separation"}
    return {
        "calibration_intercept": float(theta[0]), "calibration_slope": float(theta[1]),
        "calibration_intercept_se": float(np.sqrt(covariance[0, 0])),
        "calibration_slope_se": float(np.sqrt(covariance[1, 1])),
        "calibration_gradient_max_abs": float(np.max(np.abs(gradient))),
        "calibration_status": "evaluable", "calibration_failure_type": "",
    }


def metrics(frame: pd.DataFrame, layer: str, calibration=True) -> dict:
    y = frame.outcome_value.to_numpy(dtype=float)
    p = frame[f"p_{layer}"].to_numpy(dtype=float)
    clipped = np.clip(p, CLIP, 1 - CLIP)
    result = {
        "layer": layer, "n": len(frame), "positives": int(y.sum()),
        "prevalence": float(y.mean()), "mean_prediction": float(p.mean()),
        "brier": float(np.mean((p - y) ** 2)),
        "log_loss": float(-np.mean(y * np.log(clipped) + (1 - y) * np.log1p(-clipped))),
    }
    for name, reference in [
        ("bss_adaptation_train", "adaptation_train_reference"),
        ("bss_target_oracle", "target_oracle_reference"),
        ("bss_transport", "bss_transport_reference"),
    ]:
        r = frame[reference].to_numpy(dtype=float)
        result[name] = float(1 - np.sum((p - y) ** 2) / np.sum((r - y) ** 2))
    if calibration:
        result.update(calibration_independent(y, p))
    return result


def verify_hashes_at_start() -> None:
    fixed = {
        PROTOCOL: "919a8010610e8b1abd9ac4a6c20f4516e60c1156c9a85e19ebcd4682d2b54df1",
        U_FILE: "3b256c3fd0f813885ce817bfa6df8f80cd23ac57473047ff976ee87e0b70b2c1",
        UPSTREAM / "rq3e_external_pairs.csv": "72f4127d410160224701729fa401b4076297641402aa0965a3f9aec250a8002e",
        MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
    }
    for path, digest in fixed.items():
        if sha256(path) != digest:
            raise AssertionError(f"Frozen input hash mismatch: {path}")
    record("frozen_input_and_master_hashes", files=len(fixed), master_unchanged=True)


def main() -> None:
    verify_hashes_at_start()
    sys.path.insert(0, str(ROOT / "analysis/rq12_lib"))
    import parasail

    u_all = pd.read_csv(U_FILE)
    u = u_all.loc[u_all.config_id.eq(CONFIG)].copy()
    raw = pd.read_csv(RAW_FILE, usecols=["uuid", "sequence"])
    data = u[["uuid", "target", "stratum", "eta_U", "predicted_probability",
              "bss_transport_reference", "bss_target_oracle_reference", "outcome_value"]].merge(
                  raw, on="uuid", how="left", validate="one_to_one")
    data = data.rename(columns={"predicted_probability": "p_U"})
    assert len(data) == 720 and data.uuid.nunique() == 720
    assert data.target.value_counts().to_dict() == {target: 90 for target in TARGETS}
    assert int(data.outcome_value.sum()) == 287
    assert data.sequence.notna().all()
    close(data.p_U, expit(data.eta_U), "frozen_U_expit_identity", atol=1e-14, rtol=0)

    saved_oof = pd.read_csv(RUN / "oof_predictions.csv")
    saved_outer = pd.read_csv(RUN / "outer_assignments_frozen.csv")
    joined = saved_oof[["uuid", "eta_U", "p_U", "outcome_value", "bss_transport_reference",
                        "bss_target_oracle_reference"]].merge(
        u[["uuid", "eta_U", "predicted_probability", "outcome_value", "bss_transport_reference",
           "bss_target_oracle_reference"]], on="uuid", suffixes=("_saved", "_upstream"),
        validate="one_to_one")
    for saved, upstream in [("eta_U_saved", "eta_U_upstream"), ("p_U", "predicted_probability"),
                            ("outcome_value_saved", "outcome_value_upstream"),
                            ("bss_transport_reference_saved", "bss_transport_reference_upstream"),
                            ("bss_target_oracle_reference_saved", "bss_target_oracle_reference_upstream")]:
        close(joined[saved], joined[upstream], f"unchanged_{saved}", atol=1e-14, rtol=0)
    record("U_and_labels_frozen", pairs=len(joined), U_refit_detected=False)

    # Recompute every within-target pairwise local alignment and every qualifying edge.
    edge_rows = []
    cluster_parts = []
    cluster_summary_rows = []
    for target in TARGETS:
        block = data.loc[data.target.eq(target)].reset_index(drop=True)
        sequences = block.sequence.astype(str).str.upper().tolist()
        pairs = []
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                result = parasail.sw_trace_striped_16(sequences[i], sequences[j], 11, 1, parasail.blosum62)
                identity, query_coverage, target_coverage = alignment_metrics(
                    result, len(sequences[i]), len(sequences[j]))
                if identity >= 0.90 and query_coverage >= 0.80 and target_coverage >= 0.80:
                    pairs.append((i, j))
                    edge_rows.append({"target": target, "uuid_1": block.uuid[i], "uuid_2": block.uuid[j],
                                      "identity": identity, "query_coverage": query_coverage,
                                      "target_coverage": target_coverage})
        component = components_from_edges(len(block), pairs)
        labels = [f"{target}|cluster_{value + 1:03d}" for value in component]
        part = block[["uuid", "target"]].copy()
        part["sequence_cluster90"] = labels
        cluster_parts.append(part)
        sizes = part.groupby("sequence_cluster90").size()
        cluster_summary_rows.append({"target": target, "n_pairs": len(part), "n_clusters": len(sizes),
                                     "largest_cluster": int(sizes.max()),
                                     "singleton_clusters": int((sizes == 1).sum())})
    clusters = pd.concat(cluster_parts, ignore_index=True)
    compare_frames("sequence_edges_rebuilt", pd.DataFrame(edge_rows),
                   pd.read_csv(RUN / "sequence_cluster90_edges.csv"), ["target", "uuid_1", "uuid_2"])
    compare_frames("sequence_clusters_rebuilt", clusters.sort_values(["target", "uuid"]).reset_index(drop=True),
                   pd.read_csv(RUN / "sequence_cluster90_assignments.csv"), ["uuid", "target"])
    compare_frames("sequence_cluster_summaries", pd.DataFrame(cluster_summary_rows),
                   pd.read_csv(RUN / "sequence_cluster90_summary.csv"), ["target"])
    data = data.merge(clusters, on=["uuid", "target"], validate="one_to_one")

    # Rebuild outer folds from upstream data and independently rebuilt sequence clusters.
    data["outer_fold"] = -1
    outer_manifest_rows = []
    for target in TARGETS:
        positions = np.flatnonzero(data.target.eq(target).to_numpy())
        block = data.iloc[positions]
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=20261112)
        for fold, (train_pos, test_pos) in enumerate(splitter.split(
                np.zeros(len(block)), block.outcome_value.to_numpy(int), block.sequence_cluster90), start=1):
            data.iloc[positions[test_pos], data.columns.get_loc("outer_fold")] = fold
            outer_manifest_rows.append({"target": target, "outer_fold": fold,
                "train_n": len(train_pos), "train_positives": int(block.outcome_value.to_numpy(int)[train_pos].sum()),
                "test_n": len(test_pos), "test_positives": int(block.outcome_value.to_numpy(int)[test_pos].sum())})
    assert (data.outer_fold > 0).all()
    outer = data[["uuid", "target", "stratum", "sequence_cluster90", "outer_fold", "outcome_value", "eta_U", "p_U"]]
    compare_frames("outer_assignments_rebuilt", outer.sort_values(["target", "outer_fold", "uuid"]).reset_index(drop=True),
                   saved_outer, ["uuid", "target"])
    compare_frames("outer_fold_manifest", pd.DataFrame(outer_manifest_rows),
                   pd.read_csv(RUN / "outer_fold_manifest.csv"), ["target", "outer_fold"])
    for (target, fold), test in data.groupby(["target", "outer_fold"]):
        train = data.loc[data.target.eq(target) & data.outer_fold.ne(fold)]
        assert set(test.sequence_cluster90).isdisjoint(set(train.sequence_cluster90))
    record("outer_sequence_cluster_isolation", outer_folds=40, leakage_pairs=0)

    # Independently fit every inner candidate and final outer model.
    saved_inner_predictions = pd.read_csv(RUN / "inner_validation_predictions.csv")
    saved_lambda = pd.read_csv(RUN / "lambda_selection.csv")
    saved_parameters = pd.read_csv(RUN / "fold_parameters.csv")
    inner_assignment_rows = []
    inner_prediction_rows = []
    lambda_rows = []
    parameter_rows = []
    oof_parts = []
    fold_metric_rows = []
    intercept_score_errors = []
    final_parameter_errors = []
    inner_parameter_errors = []
    for target in TARGETS:
        target_data = data.loc[data.target.eq(target)].copy()
        for fold in range(1, 6):
            train = target_data.loc[target_data.outer_fold.ne(fold)].copy()
            test = target_data.loc[target_data.outer_fold.eq(fold)].copy()
            eta_train = train.eta_U.to_numpy(dtype=float)
            y_train = train.outcome_value.to_numpy(dtype=float)
            a_i = brentq(lambda a: float(np.mean(expit(a + eta_train)) - y_train.mean()),
                         -40, 40, xtol=1e-13, rtol=1e-14)
            eta_i = a_i + eta_train
            p_i_train = expit(eta_i)
            intercept_score_errors.append(abs(float(np.mean(p_i_train) - y_train.mean())))
            i_objective = float(np.mean(np.logaddexp(0.0, eta_i) - y_train * eta_i))
            parameter_rows.append({"target": target, "outer_fold": fold, "layer": "I", "status": "evaluable",
                "failure_type": "", "a": a_i, "delta": 0.0, "gamma": 1.0, "lambda": np.nan,
                "objective": i_objective, "gradient_max_abs": abs(float(np.mean(p_i_train - y_train))),
                "train_n": len(train), "train_positives": int(y_train.sum()), "test_n": len(test),
                "test_positives": int(test.outcome_value.sum())})
            test["a_I"] = a_i
            test["p_I"] = expit(a_i + test.eta_U.to_numpy(dtype=float))

            inner = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=20261113)
            splits = list(inner.split(np.zeros(len(train)), y_train.astype(int), train.sequence_cluster90))
            for inner_fold, (itr, iva) in enumerate(splits, start=1):
                assert set(train.sequence_cluster90.iloc[itr]).isdisjoint(set(train.sequence_cluster90.iloc[iva]))
                for position in iva:
                    inner_assignment_rows.append({"target": target, "outer_fold": fold, "inner_fold": inner_fold,
                        "uuid": train.uuid.iloc[position], "sequence_cluster90": train.sequence_cluster90.iloc[position]})
            candidates = {}
            for lam in LAMBDAS:
                validation_parts = []
                for inner_fold, (itr, iva) in enumerate(splits, start=1):
                    fit = fit_f_independent(eta_train[itr], y_train[itr], lam)
                    p_validation = expit(fit["a"] + fit["gamma"] * eta_train[iva])
                    validation_parts.append(pd.DataFrame({"target": target, "outer_fold": fold,
                        "inner_fold": inner_fold, "uuid": train.uuid.iloc[iva].to_numpy(), "lambda": lam,
                        "outcome_value": y_train[iva], "predicted_probability": p_validation,
                        "a": fit["a"], "delta": fit["delta"], "gamma": fit["gamma"],
                        "gradient_max_abs": float(np.max(np.abs(fit["gradient"])))}))
                    saved_fit = saved_inner_predictions.loc[
                        saved_inner_predictions.target.eq(target) & saved_inner_predictions.outer_fold.eq(fold)
                        & saved_inner_predictions.inner_fold.eq(inner_fold) & saved_inner_predictions["lambda"].eq(lam)]
                    assert len(saved_fit) == len(iva)
                    inner_parameter_errors.append(max(
                        abs(fit["a"] - float(saved_fit.a.iloc[0])),
                        abs(fit["delta"] - float(saved_fit.delta.iloc[0]))))
                values = pd.concat(validation_parts, ignore_index=True)
                inner_prediction_rows.append(values)
                p = np.clip(values.predicted_probability.to_numpy(dtype=float), CLIP, 1 - CLIP)
                yy = values.outcome_value.to_numpy(dtype=float)
                loss = float(-np.mean(yy * np.log(p) + (1 - yy) * np.log1p(-p)))
                candidates[lam] = loss
                lambda_rows.append({"target": target, "outer_fold": fold, "lambda": lam,
                    "candidate_status": "evaluable", "failure_type": "", "pooled_inner_log_loss": loss,
                    "n_inner_validation": len(values), "selected": False})
            minimum_loss = min(candidates.values())
            selected = max(lam for lam, loss in candidates.items() if abs(loss - minimum_loss) <= 1e-12)
            for row in lambda_rows:
                if row["target"] == target and row["outer_fold"] == fold and row["lambda"] == selected:
                    row["selected"] = True
            fit_f = fit_f_independent(eta_train, y_train, selected)
            parameter_rows.append({"target": target, "outer_fold": fold, "layer": "F", "status": "evaluable",
                "failure_type": "", "a": fit_f["a"], "delta": fit_f["delta"], "gamma": fit_f["gamma"],
                "lambda": selected, "objective": fit_f["objective"],
                "gradient_max_abs": float(np.max(np.abs(fit_f["gradient"]))),
                "train_n": len(train), "train_positives": int(y_train.sum()), "test_n": len(test),
                "test_positives": int(test.outcome_value.sum())})
            saved_f = saved_parameters.loc[saved_parameters.target.eq(target) & saved_parameters.outer_fold.eq(fold)
                                           & saved_parameters.layer.eq("F")].iloc[0]
            final_parameter_errors.append(max(abs(fit_f["a"] - saved_f.a), abs(fit_f["delta"] - saved_f.delta)))
            test["a_F"] = fit_f["a"]
            test["delta_F"] = fit_f["delta"]
            test["gamma_F"] = fit_f["gamma"]
            test["selected_lambda"] = selected
            test["p_F"] = expit(fit_f["a"] + fit_f["gamma"] * test.eta_U.to_numpy(dtype=float))
            test["adaptation_train_reference"] = float(y_train.mean())
            test["target_oracle_reference"] = float(target_data.outcome_value.mean())
            oof_parts.append(test)
            for layer in ["U", "I", "F"]:
                fold_metric_rows.append({"target": target, "outer_fold": fold, **metrics(test, layer, calibration=False)})

    maximum("I_intercept_training_score_error", max(intercept_score_errors))
    maximum("final_F_independent_optimizer_max_parameter_error", max(final_parameter_errors))
    maximum("inner_F_independent_optimizer_max_parameter_error", max(inner_parameter_errors))
    record("I_score_equations", outer_fits=40, max_abs_error=max(intercept_score_errors))
    record("F_independent_optimizer", final_fits=40, inner_fits=720,
           max_final_parameter_error=max(final_parameter_errors),
           max_inner_parameter_error=max(inner_parameter_errors))
    compare_frames("inner_assignments_rebuilt", pd.DataFrame(inner_assignment_rows),
                   pd.read_csv(RUN / "inner_assignments.csv"), ["target", "outer_fold", "uuid"])
    inner_predictions = pd.concat(inner_prediction_rows, ignore_index=True)
    compare_frames("inner_predictions_and_loss_inputs", inner_predictions, saved_inner_predictions,
                   ["target", "outer_fold", "inner_fold", "uuid", "lambda"], atol=2e-6)
    compare_frames("lambda_losses_and_selection", pd.DataFrame(lambda_rows), saved_lambda,
                   ["target", "outer_fold", "lambda"], atol=2e-6)
    compare_frames("outer_I_F_parameters", pd.DataFrame(parameter_rows), saved_parameters,
                   ["target", "outer_fold", "layer"], atol=2e-6)

    oof = pd.concat(oof_parts, ignore_index=True).sort_values(["target", "uuid"])
    oof_columns = ["uuid", "target", "stratum", "sequence_cluster90", "outer_fold", "outcome_value",
        "eta_U", "p_U", "a_I", "p_I", "a_F", "delta_F", "gamma_F", "selected_lambda", "p_F",
        "adaptation_train_reference", "target_oracle_reference", "bss_transport_reference", "bss_target_oracle_reference"]
    compare_frames("complete_OOF_predictions", oof[oof_columns], saved_oof, ["uuid", "target"], atol=2e-6)
    compare_frames("fold_metrics", pd.DataFrame(fold_metric_rows), pd.read_csv(RUN / "fold_metrics.csv"),
                   ["target", "outer_fold", "layer"], atol=2e-6)

    # Recompute all target, comparison, reliability and grouped summaries from OOF predictions.
    target_rows = []
    reliability_rows = []
    comparison_rows = []
    for target in TARGETS:
        group = oof.loc[oof.target.eq(target)]
        layer_metrics = {}
        for layer in ["U", "I", "F"]:
            layer_metrics[layer] = metrics(group, layer)
            target_rows.append({"target": target, "stratum": group.stratum.iloc[0], **layer_metrics[layer]})
            p = group[f"p_{layer}"].to_numpy(dtype=float)
            y = group.outcome_value.to_numpy(dtype=int)
            bins = np.digitize(p, EDGES[1:-1], right=False)
            for bin_index, (lower, upper) in enumerate(zip(EDGES[:-1], EDGES[1:])):
                take = bins == bin_index
                reliability_rows.append({"target": target, "layer": layer, "bin_index": bin_index,
                    "bin_lower": lower, "bin_upper": upper, "n": int(take.sum()), "positives": int(y[take].sum()),
                    "mean_prediction": float(p[take].mean()) if take.any() else np.nan,
                    "observed_rate": float(y[take].mean()) if take.any() else np.nan})
        for comparison, left, right in [("I_minus_U", "I", "U"), ("F_minus_U", "F", "U"), ("F_minus_I", "F", "I")]:
            comparison_rows.append({"target": target, "stratum": group.stratum.iloc[0], "comparison": comparison,
                "brier_difference": layer_metrics[left]["brier"] - layer_metrics[right]["brier"],
                "log_loss_difference": layer_metrics[left]["log_loss"] - layer_metrics[right]["log_loss"]})
    target_metrics = pd.DataFrame(target_rows)
    compare_frames("target_metrics_and_calibration", target_metrics, pd.read_csv(RUN / "target_metrics.csv"),
                   ["target", "layer"], atol=2e-6)
    compare_frames("target_comparisons", pd.DataFrame(comparison_rows), pd.read_csv(RUN / "target_comparisons.csv"),
                   ["target", "comparison"], atol=2e-6)
    compare_frames("reliability_records", pd.DataFrame(reliability_rows), pd.read_csv(RUN / "reliability.csv"),
                   ["target", "layer", "bin_index"], atol=2e-6)
    failed = target_metrics.loc[target_metrics.calibration_status.eq("failed"), ["target", "layer"]]
    record("calibration_diagnostics_rebuilt", failed=failed.to_dict("records"),
           standard_error="inverse unnormalised Hessian X'WX", pseudo_inverse_used=False)

    group_rows = []
    for stratum, targets in [("all_8", set(TARGETS)), ("new_target_4", NEW_TARGETS),
                             ("same_target_campaign_4", set(TARGETS) - NEW_TARGETS)]:
        for layer in ["U", "I", "F"]:
            subset = target_metrics.loc[target_metrics.layer.eq(layer) & target_metrics.target.isin(targets)]
            row = {"aggregation": "target_equal_macro", "stratum": stratum, "layer": layer,
                   "n_targets": len(subset), "n_pairs": int(subset.n.sum()), "positives": int(subset.positives.sum())}
            for name in ["prevalence", "mean_prediction", "brier", "log_loss", "bss_adaptation_train",
                         "bss_target_oracle", "bss_transport"]:
                row[name] = float(subset[name].mean())
                row[f"{name}_min"] = float(subset[name].min())
                row[f"{name}_max"] = float(subset[name].max())
            group_rows.append(row)
        comparisons = pd.DataFrame(comparison_rows)
        for comparison, left, right in [("I_minus_U", "I", "U"), ("F_minus_U", "F", "U"), ("F_minus_I", "F", "I")]:
            subset = comparisons.loc[comparisons.comparison.eq(comparison) & comparisons.target.isin(targets)]
            reference = target_metrics.loc[target_metrics.layer.eq("U") & target_metrics.target.isin(targets)]
            group_rows.append({"aggregation": "target_equal_macro_comparison", "stratum": stratum,
                "layer": comparison, "n_targets": len(subset), "n_pairs": int(reference.n.sum()),
                "positives": int(reference.positives.sum()), "brier": float(subset.brier_difference.mean()),
                "log_loss": float(subset.log_loss_difference.mean()),
                "brier_min": float(subset.brier_difference.min()), "brier_max": float(subset.brier_difference.max()),
                "log_loss_min": float(subset.log_loss_difference.min()), "log_loss_max": float(subset.log_loss_difference.max()),
                "brier_improved_targets": int(subset.brier_difference.lt(0).sum()),
                "log_loss_improved_targets": int(subset.log_loss_difference.lt(0).sum())})
        pooled = oof.loc[oof.target.isin(targets)]
        for layer in ["U", "I", "F"]:
            group_rows.append({"aggregation": "pair_pooled", "stratum": stratum, **metrics(pooled, layer)})
    group_metrics = pd.DataFrame(group_rows)
    compare_frames("group_metrics", group_metrics, pd.read_csv(RUN / "group_metrics.csv"),
                   ["aggregation", "stratum", "layer"], atol=2e-6)

    # Independently check the exact numerical claims intended for reporting.
    macro = group_metrics.loc[group_metrics.aggregation.eq("target_equal_macro") & group_metrics.stratum.eq("all_8")].set_index("layer")
    macro_c = group_metrics.loc[group_metrics.aggregation.eq("target_equal_macro_comparison") & group_metrics.stratum.eq("all_8")].set_index("layer")
    expected = {
        ("U", "brier"): 0.2590, ("U", "log_loss"): 0.7120, ("U", "bss_transport"): -0.0687,
        ("I", "brier"): 0.1918, ("I", "bss_adaptation_train"): 0.0561,
        ("F", "brier"): 0.1861, ("F", "bss_adaptation_train"): 0.0807,
        ("F_minus_I", "brier"): -0.0057,
    }
    for (layer, metric), value in expected.items():
        table = macro_c if layer.endswith(("_minus_U", "_minus_I")) else macro
        assert round(float(table.loc[layer, metric]), 4) == value
    assert int(macro_c.loc["I_minus_U", "brier_improved_targets"]) == 7
    assert int(macro_c.loc["F_minus_I", "brier_improved_targets"]) == 2
    lambda_counts = oof.groupby("selected_lambda").uuid.nunique().to_dict()
    assert lambda_counts == {0.001: 214, 0.01: 144, 0.1: 71, 1.0: 18, 100.0: 273}
    record("headline_numbers_reconstructed", macro_rows=9, comparison_rows=9, pooled_rows=9,
           I_brier_improved_targets=7, F_minus_I_brier_improved_targets=2, lambda_counts=lambda_counts)

    # Verify the frozen manifest before adding this audit's own two files.
    manifest = pd.read_csv(RUN / "final_hash_manifest.csv")
    assert len(manifest) == 24
    for row in manifest.itertuples(index=False):
        path = RUN / row.file
        assert path.stat().st_size == row.bytes
        assert sha256(path) == row.sha256
    assert sha256(MASTER) == "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91"
    record("final_hash_manifest", files=len(manifest), mismatches=0, master_unchanged=True)
    REPORT["output_row_counts"] = {"oof": 720, "fold_metrics": 120, "target_metrics": 24,
        "target_comparisons": 24, "group_metrics": 27, "reliability": 168,
        "inner_validation_predictions": 17280, "lambda_candidates": 240}
    REPORT["status"] = "passed"


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        REPORT["status"] = "failed"
        REPORT["error"] = str(error)
        REPORT["traceback"] = traceback.format_exc()
        raise
    finally:
        (RUN / "independent_verification.json").write_text(
            json.dumps(REPORT, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
