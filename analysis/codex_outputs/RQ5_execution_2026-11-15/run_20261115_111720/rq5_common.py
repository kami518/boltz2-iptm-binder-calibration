"""Shared, deterministic utilities for the frozen RQ5 execution.

RQ5 is an exploratory incremental analysis.  This module reuses the frozen
RQ3 numerical formulas where explicitly requested, while keeping all RQ5
writes inside this run directory.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

sys.dont_write_bytecode = True

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score

RUN_DIR = Path(__file__).resolve().parent
ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RQ5_HANDOFF_DIR = (
    ROOT / "analysis" / "codex_tasks" / "RQ5_EXECUTION_CODEX_HANDOFF_2026-11-15"
)
RQ5_HANDOFF_PATH = RQ5_HANDOFF_DIR / "START_HERE_CODEX.md"
RQ3_HANDOFF_DIR = (
    ROOT / "analysis" / "codex_tasks" / "RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13"
)
RQ3_CONFIG_PATH = RQ3_HANDOFF_DIR / "rq3_execution_config.json"
RQ3_IMPLEMENTATION_DIR = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_F_correction_2026-11-13"
    / "run_20261113_182311"
)
CONFIG_PATH = RUN_DIR / "run_config.snapshot.json"
STATUS_PATH = RUN_DIR / "RUN_STATUS.json"
FAILURE_PATH = RUN_DIR / "failed_item_ledger.csv"

# Import the accepted frozen implementation only for formula reuse.  RQ5 does
# not modify, regenerate, or reuse its 2,084-row result predictions.
if str(RQ3_IMPLEMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(RQ3_IMPLEMENTATION_DIR))
import rq3_common as rq3  # noqa: E402

ANALYSIS_SETS = ("primary_six", "primary_five", "strict_six")
OUTCOMES = {
    "primary_six": "primary_outcome",
    "primary_five": "primary_outcome",
    "strict_six": "strict_outcome",
}
CLUSTER_RULE = "sequence_cluster90"
SOURCE_WEIGHTS = ("target_equal", "pair")
SCORES = {
    "ipTM": ("boltz2_iptm", "x_ipTM"),
    "Pinc": ("pinc", "x_Pinc"),
}
EXPOSURES = ("ipTM_only", "Pinc_only", "ipTM_plus_Pinc")
ADAPTATION_LAYERS = ("U_transport_zero_label", "predictive_I", "predictive_F")
LAMBDA_GRID = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
RELIABILITY_EDGES = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0)
BOOTSTRAP_SEED = 20261115
BOOTSTRAP_REPLICATES = 2000


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, set):
        return sorted(_jsonable(item) for item in value)
    return value


def atomic_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(
            _jsonable(payload),
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        handle.write("\n")
    os.replace(temporary, path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8")
    os.replace(temporary, path)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_config() -> dict:
    return read_json(CONFIG_PATH)


def update_status(**updates: Any) -> dict:
    if STATUS_PATH.exists():
        status = read_json(STATUS_PATH)
    else:
        status = {
            "handoff_id": "RQ5_EXECUTION_2026-11-15",
            "run_directory": str(RUN_DIR),
            "overall_status": "initializing",
            "created_at": now_iso(),
            "analysis_scope": "exploratory_incremental_only",
            "stages": {},
        }
    status.update(updates)
    status["updated_at"] = now_iso()
    atomic_json(STATUS_PATH, status)
    return status


def set_stage_status(stage: str, status: str, **details: Any) -> None:
    current = update_status()
    stages = current.setdefault("stages", {})
    stages[stage] = {"status": status, "updated_at": now_iso(), **details}
    values = [item.get("status") for item in stages.values()]
    if status == "running":
        overall = "running"
    elif "blocked" in values:
        overall = "blocked"
    elif "failed" in values:
        overall = "failed"
    elif any(value in {"failed", "blocked"} for value in values):
        overall = "partial"
    elif all(value == "complete" for value in values):
        overall = "complete"
    else:
        overall = "partial"
    update_status(current_stage=stage, stages=stages, overall_status=overall)


def append_failure(
    stage: str,
    item_id: str,
    failure_type: str,
    detail: str,
    **context: Any,
) -> None:
    row = {
        "timestamp": now_iso(),
        "stage": stage,
        "item_id": item_id,
        "failure_type": failure_type,
        "detail": detail,
        "context_json": json.dumps(context, ensure_ascii=False, sort_keys=True),
    }
    if FAILURE_PATH.exists():
        existing = pd.read_csv(FAILURE_PATH)
        duplicate = (
            existing["stage"].astype(str).eq(stage)
            & existing["item_id"].astype(str).eq(item_id)
            & existing["failure_type"].astype(str).eq(failure_type)
            & existing["detail"].astype(str).eq(detail)
        )
        if duplicate.any():
            return
        frame = pd.concat([existing, pd.DataFrame([row])], ignore_index=True)
    else:
        frame = pd.DataFrame([row])
    atomic_csv(FAILURE_PATH, frame)


def ensure_failure_ledger() -> None:
    if not FAILURE_PATH.exists():
        atomic_csv(
            FAILURE_PATH,
            pd.DataFrame(
                columns=[
                    "timestamp",
                    "stage",
                    "item_id",
                    "failure_type",
                    "detail",
                    "context_json",
                ]
            ),
        )


def environment_snapshot() -> dict:
    return {
        "captured_at": now_iso(),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "network_access_used": False,
        "dependency_install_or_upgrade_performed": False,
    }


def logit_clip(probability: Iterable[float], clip: float = 1e-6) -> np.ndarray:
    p = np.asarray(list(probability) if not isinstance(probability, np.ndarray) else probability, dtype=float)
    clipped = np.clip(p, clip, 1.0 - clip)
    return np.log(clipped) - np.log1p(-clipped)


def load_core_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    specs = config["inputs"]
    master = pd.read_csv(specs["pair_master"]["path"], low_memory=False)
    pinc = pd.read_csv(specs["pinc_extension"]["path"], low_memory=False)
    clusters = pd.read_csv(specs["sequence_clusters"]["path"], low_memory=False)
    if len(master) != int(specs["pair_master"]["expected_rows"]):
        raise RuntimeError("Frozen pair master row count mismatch")
    if len(pinc) != int(specs["pinc_extension"]["expected_rows"]):
        raise RuntimeError("Pinc extension row count mismatch")
    if len(clusters) != int(specs["sequence_clusters"]["expected_rows"]):
        raise RuntimeError("Sequence cluster row count mismatch")
    for name, frame in (("master", master), ("pinc", pinc), ("clusters", clusters)):
        if frame["pair_id"].duplicated().any():
            raise RuntimeError(f"pair_id is not unique in {name}")
        if frame["pair_id"].isna().any():
            raise RuntimeError(f"pair_id is missing in {name}")
    if set(master["pair_id"]) != set(pinc["pair_id"]):
        raise RuntimeError("Pinc pair_id set differs from master")
    if set(master["pair_id"]) != set(clusters["pair_id"]):
        raise RuntimeError("Cluster pair_id set differs from master")

    joined = master.merge(
        pinc[
            [
                "pair_id",
                "design_id",
                "target",
                "pinc",
                "pinc_status",
                "pinc_reason",
            ]
        ],
        on="pair_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_pinc"),
    )
    if len(joined) != len(master):
        raise RuntimeError("Pinc one-to-one join lost rows")
    for column in ("design_id", "target"):
        if not joined[column].astype(str).eq(
            joined[f"{column}_pinc"].astype(str)
        ).all():
            raise RuntimeError(f"Master and Pinc disagree on {column}")
    joined = joined.drop(columns=["design_id_pinc", "target_pinc"])

    joined = joined.merge(
        clusters[
            [
                "pair_id",
                "target",
                "design_id",
                "sequence",
                "sequence_cluster90",
                "sequence_cluster80",
            ]
        ],
        on="pair_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_cluster"),
    )
    if len(joined) != len(master):
        raise RuntimeError("Cluster one-to-one join lost rows")
    for column in ("target", "design_id", "sequence"):
        if not joined[column].astype(str).eq(
            joined[f"{column}_cluster"].astype(str)
        ).all():
            raise RuntimeError(f"Master and cluster inputs disagree on {column}")
    joined = joined.drop(
        columns=["target_cluster", "design_id_cluster", "sequence_cluster"]
    )
    joined["x_ipTM"] = np.nan
    joined["x_Pinc"] = np.nan
    success = joined["pinc_status"].eq("success")
    joined.loc[success, "x_ipTM"] = logit_clip(
        joined.loc[success, "boltz2_iptm"].to_numpy(dtype=float),
        float(config["probability_clip"]),
    )
    joined.loc[success, "x_Pinc"] = logit_clip(
        joined.loc[success, "pinc"].to_numpy(dtype=float),
        float(config["probability_clip"]),
    )
    return master, pinc, clusters, joined


def analysis_frame(config: dict, joined: pd.DataFrame, analysis_set: str) -> pd.DataFrame:
    if analysis_set not in ANALYSIS_SETS:
        raise ValueError(f"Unknown RQ5 analysis set: {analysis_set}")
    frame = joined[joined["pinc_status"].eq("success")].copy()
    if analysis_set == "primary_five":
        frame = frame[frame["target"].astype(str).ne("aah-i")].copy()
    outcome = OUTCOMES[analysis_set]
    frame = frame[frame[outcome].notna()].copy()
    frame["outcome_value"] = frame[outcome].astype(int)
    if not np.isfinite(
        frame[["boltz2_iptm", "pinc", "x_ipTM", "x_Pinc"]].to_numpy(dtype=float)
    ).all():
        raise RuntimeError(f"Nonfinite RQ5 score or transform in {analysis_set}")
    return frame.sort_values("pair_id").reset_index(drop=True)


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    yy = np.asarray(y, dtype=float)
    ss = np.asarray(score, dtype=float)
    if np.unique(yy).size < 2 or not np.isfinite(ss).all():
        return np.nan
    return float(roc_auc_score(yy, ss))


def safe_ap(y: np.ndarray, score: np.ndarray) -> float:
    yy = np.asarray(y, dtype=float)
    ss = np.asarray(score, dtype=float)
    if np.sum(yy) <= 0 or not np.isfinite(ss).all():
        return np.nan
    return float(average_precision_score(yy, ss))


def brier(y: np.ndarray, probability: np.ndarray) -> float:
    return float(np.mean((np.asarray(probability, dtype=float) - np.asarray(y, dtype=float)) ** 2))


def log_loss(y: np.ndarray, probability: np.ndarray, clip: float) -> float:
    yy = np.asarray(y, dtype=float)
    pp = np.clip(np.asarray(probability, dtype=float), clip, 1.0 - clip)
    return float(-np.mean(yy * np.log(pp) + (1.0 - yy) * np.log1p(-pp)))


def bss(y: np.ndarray, probability: np.ndarray, reference: np.ndarray) -> float:
    return rq3.bss_from_reference(np.asarray(y), np.asarray(probability), np.asarray(reference))


def score_point_metrics(
    y: np.ndarray,
    score: np.ndarray,
    *,
    clip: float,
    include_calibration: bool = True,
) -> dict:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    baseline = float(np.mean(y))
    result = {
        "n": int(len(y)),
        "positives": int(np.sum(y)),
        "prevalence": baseline,
        "score_mean": float(np.mean(score)),
        "score_sd": float(np.std(score, ddof=1)) if len(score) > 1 else 0.0,
        "score_median": float(np.median(score)),
        "score_min": float(np.min(score)),
        "score_max": float(np.max(score)),
        "roc_auc": safe_auc(y, score),
        "ap": safe_ap(y, score),
        "brier": brier(y, score),
        "log_loss": log_loss(y, score, clip),
        "bss_target_constant": bss(y, score, np.full(len(y), baseline)),
    }
    if include_calibration:
        result.update(rq3.fit_calibration(y, score, clip))
    return result


def performance_point_metrics(
    y: np.ndarray,
    probability: np.ndarray,
    *,
    clip: float,
    transport_reference: np.ndarray,
    target_oracle_reference: np.ndarray,
) -> dict:
    y = np.asarray(y, dtype=float)
    probability = np.asarray(probability, dtype=float)
    result = {
        "n": int(len(y)),
        "positives": int(np.sum(y)),
        "prevalence": float(np.mean(y)),
        "roc_auc": safe_auc(y, probability),
        "ap": safe_ap(y, probability),
        "brier": brier(y, probability),
        "log_loss": log_loss(y, probability, clip),
        "bss_transport": bss(y, probability, transport_reference),
        "bss_target_oracle": bss(y, probability, target_oracle_reference),
    }
    if np.isfinite(probability).all() and np.unique(y).size == 2:
        result.update(rq3.fit_calibration(y, probability, clip))
    else:
        result.update(
            {
                "calibration_intercept": np.nan,
                "calibration_slope": np.nan,
                "calibration_intercept_se": np.nan,
                "calibration_slope_se": np.nan,
                "calibration_status": "failed",
                "calibration_failure_type": "single_class_or_nonfinite",
                "calibration_message": "Calibration diagnostic not evaluable",
            }
        )
    return result


def reliability_frame(
    y: np.ndarray,
    probability: np.ndarray,
    edges: Sequence[float],
) -> list[dict]:
    frame = pd.DataFrame(
        {"outcome_value": np.asarray(y, dtype=int), "predicted_probability": np.asarray(probability, dtype=float)}
    )
    return rq3.reliability_rows(frame, "predicted_probability", edges)


def weighted_center_scale(
    values: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray, bool]:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.ndim != 2:
        raise ValueError("values must be two-dimensional")
    if not (np.isfinite(x).all() and np.isfinite(w).all()):
        return np.full(x.shape[1], np.nan), np.full(x.shape[1], np.nan), False
    if len(x) == 0 or (w <= 0).any() or float(w.sum()) <= 0:
        return np.full(x.shape[1], np.nan), np.full(x.shape[1], np.nan), False
    mean = np.sum(x * w[:, None], axis=0) / np.sum(w)
    variance = np.sum(w[:, None] * (x - mean) ** 2, axis=0) / np.sum(w)
    sd = np.sqrt(np.maximum(variance, 0.0))
    return mean, sd, bool(np.isfinite(sd).all() and (sd > 0).all())


def fit_combo_ridge(
    x: np.ndarray,
    y: np.ndarray,
    source_targets: np.ndarray,
    weights: np.ndarray,
    penalty_lambda: float,
    *,
    label: str,
    bound: float = 40.0,
) -> dict:
    """Fit source intercepts plus two jointly ridge-penalized slopes."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    targets = np.asarray(source_targets).astype(str)
    w = np.asarray(weights, dtype=float)
    levels = sorted(pd.unique(targets))
    n_slope = x.shape[1]
    design = np.column_stack([rq3.target_one_hot(targets, levels), x])
    if not (np.isfinite(design).all() and np.isfinite(y).all() and np.isfinite(w).all()):
        return {"status": "failed", "estimable": False, "failure_type": "nonfinite_input"}
    if np.unique(y).size < 2:
        return {"status": "failed", "estimable": False, "failure_type": "single_class"}
    if len(y) == 0 or (w <= 0).any() or float(w.sum()) <= 0:
        return {"status": "failed", "estimable": False, "failure_type": "invalid_weight"}

    penalty_mask = np.zeros(design.shape[1], dtype=float)
    penalty_mask[len(levels) :] = 1.0
    weight_sum = float(w.sum())

    def objective(theta: np.ndarray) -> float:
        eta = design @ theta
        nll = float(np.sum(w * rq3.stable_nll(eta, y)) / weight_sum)
        return nll + 0.5 * float(penalty_lambda) * float(np.sum(theta * theta * penalty_mask))

    def gradient(theta: np.ndarray) -> np.ndarray:
        eta = design @ theta
        return (
            design.T @ (w * (expit(eta) - y)) / weight_sum
            + float(penalty_lambda) * theta * penalty_mask
        )

    initial = np.zeros(design.shape[1], dtype=float)
    for index, level in enumerate(levels):
        mask = targets == level
        prevalence = np.clip(float(np.average(y[mask], weights=w[mask])), 1e-6, 1.0 - 1e-6)
        initial[index] = math.log(prevalence / (1.0 - prevalence))
    starts = [
        initial,
        np.r_[initial[: len(levels)], np.full(n_slope, 0.25)],
        np.r_[initial[: len(levels)], np.full(n_slope, -0.25)],
    ]
    records = []
    best = None
    for index, start in enumerate(starts):
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            jac=gradient,
            bounds=[(-bound, bound)] * design.shape[1],
            options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8, "maxls": 100},
        )
        value = float(result.fun)
        records.append(
            {
                "start_index": index,
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "objective": value,
            }
        )
        if best is None or value < best[0]:
            best = (value, result)
    if best is None:
        return {"status": "failed", "estimable": False, "failure_type": "no_valid_start"}
    objective_value, result = best
    theta = np.asarray(result.x, dtype=float)
    # L-BFGS-B can stop a few ulps above the frozen gradient gate when the
    # lambda=100 objective is strongly scaled.  Refine the same convex
    # objective with deterministic damped Newton steps; this changes only the
    # numerical solver, not the model, penalty, features, or selection rule.
    newton_iterations = 0
    newton_converged = False
    for newton_iterations in range(1, 101):
        eta = design @ theta
        gradient_value = (
            design.T @ (w * (expit(eta) - y)) / weight_sum
            + float(penalty_lambda) * theta * penalty_mask
        )
        if float(np.max(np.abs(gradient_value))) <= 1e-10:
            newton_converged = True
            break
        probability_value = expit(eta)
        hessian_value = (
            design.T
            * (w * probability_value * (1.0 - probability_value))
        ) @ design / weight_sum
        hessian_value[np.diag_indices_from(hessian_value)] += (
            penalty_mask * float(penalty_lambda)
        )
        try:
            direction = np.linalg.solve(hessian_value, -gradient_value)
        except np.linalg.LinAlgError:
            break
        current_objective = objective(theta)
        step_size = 1.0
        accepted = False
        while step_size >= 1e-12:
            candidate = theta + step_size * direction
            if np.all(np.abs(candidate) < bound - 1e-9):
                candidate_objective = objective(candidate)
                if (
                    np.isfinite(candidate_objective)
                    and candidate_objective <= current_objective + 1e-15
                ):
                    theta = candidate
                    objective_value = float(candidate_objective)
                    accepted = True
                    break
            step_size *= 0.5
        if not accepted:
            break
    if not newton_converged:
        eta = design @ theta
        gradient_value = (
            design.T @ (w * (expit(eta) - y)) / weight_sum
            + float(penalty_lambda) * theta * penalty_mask
        )
        newton_converged = bool(float(np.max(np.abs(gradient_value))) <= 1e-10)
    objective_value = float(objective(theta))
    probability = expit(design @ theta)
    hessian = (design.T * (w * probability * (1.0 - probability))) @ design / weight_sum
    hessian[np.diag_indices_from(hessian)] += penalty_mask * float(penalty_lambda)
    eigenvalues = np.linalg.eigvalsh((hessian + hessian.T) / 2.0)
    information_ok = bool(
        np.isfinite(eigenvalues).all()
        and float(eigenvalues.min()) > np.finfo(float).eps * 100.0 * max(1.0, float(eigenvalues.max()))
    )
    final_gradient = gradient(theta)
    bound_hit = bool(np.any(np.isclose(np.abs(theta), bound, atol=1e-6, rtol=0.0)))
    finite = bool(np.isfinite(theta).all() and np.isfinite(objective_value) and np.isfinite(final_gradient).all())
    estimable = bool(
        (result.success or newton_converged)
        and finite
        and information_ok
        and not bound_hit
        and float(np.max(np.abs(final_gradient))) <= 1e-5
    )
    if not result.success:
        failure_type = "optimizer_nonconvergence"
    elif not finite:
        failure_type = "nonfinite_solution"
    elif bound_hit:
        failure_type = "parameter_bound_hit"
    elif not information_ok:
        failure_type = "singular_or_ill_conditioned_information"
    elif float(np.max(np.abs(final_gradient))) > 1e-5:
        failure_type = "gradient_not_small"
    else:
        failure_type = ""
    return {
        "status": "evaluable" if estimable else "failed",
        "estimable": estimable,
        "failure_type": failure_type,
        "message": f"{label}: {'passed' if estimable else failure_type}",
        "theta": theta,
        "levels": levels,
        "objective": objective_value,
        "gradient_max_abs": float(np.max(np.abs(final_gradient))),
        "iterations": int(result.nit),
        "optimizer_success": bool(result.success),
        "newton_refinement": {
            "iterations": int(newton_iterations),
            "converged": bool(newton_converged),
        },
        "information_ok": information_ok,
        "starts": records,
    }


def source_weights(frame: pd.DataFrame, weighting: str) -> np.ndarray:
    if weighting == "target_equal":
        values = frame["weight_target_equal"].to_numpy(dtype=float)
    elif weighting == "pair":
        values = frame["weight_pair"].to_numpy(dtype=float)
    else:
        raise ValueError(f"Unknown source weighting: {weighting}")
    if not np.isfinite(values).all() or (values <= 0).any():
        raise RuntimeError(f"Invalid {weighting} source weights")
    return values


def transport_intercept(
    levels: Sequence[str],
    intercepts: np.ndarray,
    targets: np.ndarray,
    weighting: str,
) -> float:
    intercepts = np.asarray(intercepts, dtype=float)
    if weighting == "target_equal":
        return float(np.mean(intercepts))
    counts = pd.Series(np.asarray(targets).astype(str)).value_counts()
    source_counts = np.asarray([float(counts.loc[level]) for level in levels])
    return float(np.sum(source_counts * intercepts) / np.sum(source_counts))


def fit_linear_exposure(
    x: np.ndarray,
    y: np.ndarray,
    targets: np.ndarray,
    weights: np.ndarray,
    weighting: str,
    label: str,
) -> dict:
    fit = rq3.fit_source_linear(np.asarray(x), np.asarray(y), np.asarray(targets), np.asarray(weights), label)
    if fit["estimable"]:
        levels = list(fit["levels"])
        source_intercepts = np.asarray(fit["theta"][: len(levels)], dtype=float)
        slope = float(fit["theta"][len(levels)])
        fit["transport_intercept"] = transport_intercept(levels, source_intercepts, targets, weighting)
        fit["common_slope"] = slope
    else:
        fit["transport_intercept"] = np.nan
        fit["common_slope"] = np.nan
    return fit


def fit_combo_model(
    raw_x: np.ndarray,
    y: np.ndarray,
    targets: np.ndarray,
    weights: np.ndarray,
    weighting: str,
    penalty_lambda: float,
    label: str,
) -> dict:
    raw_x = np.asarray(raw_x, dtype=float)
    mean, sd, scale_ok = weighted_center_scale(raw_x, weights)
    if not scale_ok:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "zero_or_invalid_feature_sd",
            "feature_mean": mean,
            "feature_sd": sd,
        }
    transformed = (raw_x - mean) / sd
    fit = fit_combo_ridge(
        transformed,
        y,
        targets,
        weights,
        penalty_lambda,
        label=label,
    )
    fit.update({"feature_mean": mean, "feature_sd": sd})
    if fit["estimable"]:
        levels = list(fit["levels"])
        source_intercepts = np.asarray(fit["theta"][: len(levels)], dtype=float)
        slopes = np.asarray(fit["theta"][len(levels) :], dtype=float)
        fit["transport_intercept"] = transport_intercept(levels, source_intercepts, targets, weighting)
        fit["common_slopes_standardized"] = slopes
    else:
        fit["transport_intercept"] = np.nan
        fit["common_slopes_standardized"] = np.full(raw_x.shape[1], np.nan)
    return fit


def combo_eta(
    raw_x: np.ndarray,
    transport_intercept_value: float,
    slopes_standardized: np.ndarray,
    feature_mean: np.ndarray,
    feature_sd: np.ndarray,
) -> np.ndarray:
    transformed = (np.asarray(raw_x, dtype=float) - np.asarray(feature_mean, dtype=float)) / np.asarray(feature_sd, dtype=float)
    return float(transport_intercept_value) + transformed @ np.asarray(slopes_standardized, dtype=float)


def summarize_target_metrics(
    frame: pd.DataFrame,
    metrics: Sequence[str],
    group_columns: Sequence[str],
) -> pd.DataFrame:
    rows = []
    for keys, group in frame.groupby(list(group_columns), dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_columns, keys))
        for metric in metrics:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            rows.append(
                {
                    **base,
                    "summary_type": "target_equal_macro",
                    "metric": metric,
                    "mean": float(values.mean()) if len(values) else np.nan,
                    "target_count": int(len(values)),
                    "min": float(values.min()) if len(values) else np.nan,
                    "max": float(values.max()) if len(values) else np.nan,
                    "positive_sign_targets": int((values > 0).sum()),
                    "negative_sign_targets": int((values < 0).sum()),
                    "zero_sign_targets": int((values == 0).sum()),
                }
            )
    return pd.DataFrame(rows)


def save_figure(path: Path, figure: plt.Figure, dpi: int = 220) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    figure.savefig(
        temporary,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        format=path.suffix.lstrip("."),
    )
    plt.close(figure)
    os.replace(temporary, path)

