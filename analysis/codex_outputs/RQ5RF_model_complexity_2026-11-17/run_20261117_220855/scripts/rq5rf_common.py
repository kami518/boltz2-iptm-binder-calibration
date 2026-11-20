"""Shared deterministic implementation for RQ5-RF v1 frozen protocol."""
from __future__ import annotations

import gzip
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

RUN_DIR = Path(__file__).resolve().parent.parent
ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
CONFIG_PATH = RUN_DIR / "run_config.snapshot.json"
STATUS_PATH = RUN_DIR / "RUN_STATUS.json"
FAILURE_PATH = RUN_DIR / "failed_item_ledger.csv"
RQ3_IMPLEMENTATION_DIR = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_F_correction_2026-11-13"
    / "run_20261113_182311"
)

if str(RQ3_IMPLEMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(RQ3_IMPLEMENTATION_DIR))
import rq3_common as rq3  # noqa: E402

ANALYSIS_SETS = ("primary_six", "primary_five", "strict_six")
SOURCE_WEIGHTS = ("target_equal", "pair")
EXPOSURES = (
    "pooled_logistic_ipTM",
    "pooled_logistic_combo",
    "rf_ipTM",
    "rf_combo",
)
LOGISTIC_EXPOSURES = ("pooled_logistic_ipTM", "pooled_logistic_combo")
RF_EXPOSURES = ("rf_ipTM", "rf_combo")
FEATURE_MAP = {
    "pooled_logistic_ipTM": ("x_ipTM",),
    "pooled_logistic_combo": ("x_ipTM", "x_Pinc"),
    "rf_ipTM": ("x_ipTM",),
    "rf_combo": ("x_ipTM", "x_Pinc"),
}
LAYERS = ("U", "I", "F")
PRESET_COMPARISONS = (
    ("pooled_logistic_combo", "pooled_logistic_ipTM"),
    ("rf_ipTM", "pooled_logistic_ipTM"),
    ("rf_combo", "pooled_logistic_combo"),
    ("rf_combo", "rf_ipTM"),
    ("rf_combo", "pooled_logistic_ipTM"),
)
METRICS = (
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "roc_auc",
    "ap",
    "calibration_intercept",
    "calibration_slope",
)


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_config() -> dict:
    return read_json(CONFIG_PATH)


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


def atomic_gzip_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", newline="\n") as handle:
        json.dump(
            _jsonable(payload),
            handle,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=True,
        )
    os.replace(temporary, path)


def read_gzip_json(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8")
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        _jsonable(payload),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def pair_set_sha256(pair_ids: Iterable[str]) -> str:
    values = sorted(str(value) for value in pair_ids)
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def weight_vector_sha256(pair_ids: Iterable[str], weights: Iterable[float]) -> str:
    rows = [(str(pair_id), float(weight)) for pair_id, weight in zip(pair_ids, weights)]
    rows.sort(key=lambda item: item[0])
    payload = "".join(f"{pair_id}\t{weight:.17g}\n" for pair_id, weight in rows)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def task_seed(
    analysis_set: str,
    weighting: str,
    outer_fold: str,
    inner_fold: str,
    stage: str,
) -> int:
    if stage not in {"inner", "final"}:
        raise ValueError("RF stage must be inner or final")
    canonical = (
        f"RQ5RF|v1|{analysis_set}|{weighting}|{outer_fold}|"
        f"{inner_fold}|{stage}"
    )
    return int.from_bytes(hashlib.sha256(canonical.encode("utf-8")).digest()[:4], "big")


def update_status(**updates: Any) -> dict:
    status = read_json(STATUS_PATH) if STATUS_PATH.exists() else {}
    status.update(updates)
    status.setdefault(
        "handoff_id", "RQ5RF_EXECUTION_CODEX_HANDOFF_2026-11-17"
    )
    status.setdefault("run_directory", str(RUN_DIR))
    status.setdefault("analysis_scope", "post_hoc_exploratory_model_complexity_sensitivity")
    status.setdefault("created_at", now_iso())
    status["updated_at"] = now_iso()
    atomic_json(STATUS_PATH, status)
    return status


def set_stage_status(
    stage: str,
    status: str,
    *,
    overall_status: str | None = None,
    **details: Any,
) -> None:
    current = read_json(STATUS_PATH) if STATUS_PATH.exists() else {}
    stages = current.get("stages", {})
    stages[stage] = {"status": status, "updated_at": now_iso(), **details}
    current["stages"] = stages
    current["current_stage"] = stage
    if overall_status is not None:
        current["overall_status"] = overall_status
    current["updated_at"] = now_iso()
    atomic_json(STATUS_PATH, current)


def ensure_failure_ledger() -> None:
    if FAILURE_PATH.exists():
        return
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
                "resolution_status",
            ]
        ),
    )


def append_failure(
    stage: str,
    item_id: str,
    failure_type: str,
    detail: str,
    *,
    resolution_status: str = "unresolved",
    **context: Any,
) -> None:
    ensure_failure_ledger()
    existing = pd.read_csv(FAILURE_PATH)
    key = (
        existing["stage"].astype(str).eq(str(stage))
        & existing["item_id"].astype(str).eq(str(item_id))
        & existing["failure_type"].astype(str).eq(str(failure_type))
        & existing["detail"].astype(str).eq(str(detail))
    )
    if key.any():
        return
    row = pd.DataFrame(
        [
            {
                "timestamp": now_iso(),
                "stage": stage,
                "item_id": item_id,
                "failure_type": failure_type,
                "detail": detail,
                "context_json": json.dumps(
                    _jsonable(context), ensure_ascii=False, sort_keys=True
                ),
                "resolution_status": resolution_status,
            }
        ]
    )
    atomic_csv(FAILURE_PATH, pd.concat([existing, row], ignore_index=True))


def resolve_failure(stage: str, item_id: str, failure_type: str) -> None:
    if not FAILURE_PATH.exists():
        return
    frame = pd.read_csv(FAILURE_PATH)
    mask = (
        frame["stage"].astype(str).eq(str(stage))
        & frame["item_id"].astype(str).eq(str(item_id))
        & frame["failure_type"].astype(str).eq(str(failure_type))
    )
    if mask.any():
        frame.loc[mask, "resolution_status"] = "resolved"
        atomic_csv(FAILURE_PATH, frame)


def environment_snapshot() -> dict:
    return {
        "captured_at": now_iso(),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "network_access_used": False,
        "dependency_install_or_upgrade_performed": False,
    }


def stable_sort(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.assign(
        _pair_id_sort=frame["pair_id"].astype(str)
    ).sort_values("_pair_id_sort", kind="mergesort").drop(
        columns="_pair_id_sort"
    ).reset_index(drop=True)


def logit_clip(probability: Iterable[float], clip: float = 1e-6) -> np.ndarray:
    values = np.asarray(probability, dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError("Nonfinite probability before logit transform")
    if (values < 0.0).any() or (values > 1.0).any():
        raise RuntimeError("Probability outside [0,1] before logit transform")
    clipped = np.clip(values, float(clip), 1.0 - float(clip))
    return np.log(clipped) - np.log1p(-clipped)


def load_joined_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    specs = config["inputs"]
    master = pd.read_csv(specs["pair_master"]["path"], low_memory=False)
    pinc = pd.read_csv(specs["pinc_extension"]["path"], low_memory=False)
    clusters = pd.read_csv(specs["sequence_clusters"]["path"], low_memory=False)
    expected_lengths = (
        ("master", master, int(specs["pair_master"]["expected_rows"])),
        ("pinc", pinc, int(specs["pinc_extension"]["expected_rows"])),
        ("clusters", clusters, int(specs["sequence_clusters"]["expected_rows"])),
    )
    for label, frame, expected in expected_lengths:
        if len(frame) != expected:
            raise RuntimeError(f"{label} row count mismatch: {len(frame)} != {expected}")
        if frame["pair_id"].isna().any() or frame["pair_id"].duplicated().any():
            raise RuntimeError(f"{label} pair_id is missing or nonunique")
    if not (set(master["pair_id"]) == set(pinc["pair_id"]) == set(clusters["pair_id"])):
        raise RuntimeError("Core input pair_id sets differ")

    joined = master.merge(
        pinc[["pair_id", "design_id", "target", "pinc", "pinc_status", "pinc_reason"]],
        on="pair_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_pinc"),
    )
    if len(joined) != len(master):
        raise RuntimeError("Pinc join changed row count")
    for column in ("design_id", "target"):
        if not joined[column].astype(str).eq(joined[f"{column}_pinc"].astype(str)).all():
            raise RuntimeError(f"Master and Pinc disagree on {column}")
    joined = joined.drop(columns=["design_id_pinc", "target_pinc"])
    joined = joined.merge(
        clusters[
            [
                "pair_id",
                "target",
                "design_id",
                "sequence_cluster90",
                "sequence_cluster80",
            ]
        ],
        on="pair_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_cluster"),
    )
    if len(joined) != len(master):
        raise RuntimeError("Cluster join changed row count")
    for column in ("target", "design_id"):
        if not joined[column].astype(str).eq(joined[f"{column}_cluster"].astype(str)).all():
            raise RuntimeError(f"Master and cluster inputs disagree on {column}")
    joined = joined.drop(columns=["target_cluster", "design_id_cluster"])
    joined["sequence_cluster_value"] = joined["sequence_cluster90"].astype(str)

    success = joined["pinc_status"].eq("success")
    for column in ("boltz2_iptm", "pinc"):
        values = joined.loc[success, column].astype(float)
        if not np.isfinite(values).all():
            raise RuntimeError(f"Nonfinite success value in {column}")
        if (values < 0.0).any() or (values > 1.0).any():
            raise RuntimeError(f"Out-of-range success value in {column}")
    joined["x_ipTM"] = np.nan
    joined["x_Pinc"] = np.nan
    joined.loc[success, "x_ipTM"] = logit_clip(
        joined.loc[success, "boltz2_iptm"].to_numpy(dtype=float),
        float(config["probability_clip"]),
    )
    joined.loc[success, "x_Pinc"] = logit_clip(
        joined.loc[success, "pinc"].to_numpy(dtype=float),
        float(config["probability_clip"]),
    )
    return master, pinc, clusters, stable_sort(joined)


def analysis_frame(config: dict, joined: pd.DataFrame, analysis_set: str) -> pd.DataFrame:
    if analysis_set not in ANALYSIS_SETS:
        raise ValueError(f"Unknown analysis set: {analysis_set}")
    spec = config["analysis_sets"][analysis_set]
    frame = joined[joined["pinc_status"].eq("success")].copy()
    if spec["exclude_target"] is not None:
        frame = frame[frame["target"].astype(str).ne(str(spec["exclude_target"]))].copy()
    outcome_column = spec["outcome"]
    frame = frame[frame[outcome_column].notna()].copy()
    frame["outcome_value"] = frame[outcome_column].astype(int)
    if not np.isin(frame["outcome_value"].unique(), [0, 1]).all():
        raise RuntimeError(f"Nonbinary outcome in {analysis_set}")
    frame = stable_sort(frame)
    expected_rows = int(spec["expected_rows"])
    expected_positives = int(spec["expected_positives"])
    if len(frame) != expected_rows or int(frame["outcome_value"].sum()) != expected_positives:
        raise RuntimeError(
            f"{analysis_set} sample mismatch: {len(frame)}/{int(frame['outcome_value'].sum())}"
        )
    return frame


def current_source_weights(frame: pd.DataFrame, weighting: str) -> tuple[np.ndarray, dict]:
    pair_ids = frame["pair_id"].astype(str).to_numpy()
    targets = frame["target"].astype(str)
    levels = sorted(targets.unique())
    counts = targets.value_counts().to_dict()
    if weighting == "target_equal":
        if not levels:
            raise RuntimeError("Cannot form target_equal weights without source targets")
        weights = np.asarray(
            [1.0 / (len(levels) * int(counts[target])) for target in targets],
            dtype=float,
        )
        target_sums = {
            target: float(weights[np.asarray(targets == target)].sum())
            for target in levels
        }
        max_target_error = max(
            abs(value - 1.0 / len(levels)) for value in target_sums.values()
        )
    elif weighting == "pair":
        weights = np.full(len(frame), 1.0 / len(frame), dtype=float)
        target_sums = {
            target: float(weights[np.asarray(targets == target)].sum())
            for target in levels
        }
        max_target_error = np.nan
    else:
        raise ValueError(f"Unknown weighting: {weighting}")
    if not np.isfinite(weights).all() or (weights <= 0).any():
        raise RuntimeError("Source weights are nonfinite or nonpositive")
    normalization_error = abs(float(weights.sum()) - 1.0)
    if normalization_error > 1e-12:
        raise RuntimeError("Source weights do not sum to one")
    if weighting == "target_equal" and max_target_error > 1e-12:
        raise RuntimeError("Target-total source weight tolerance exceeded")
    audit = {
        "pair_id_set_sha256": pair_set_sha256(pair_ids),
        "weight_vector_sha256": weight_vector_sha256(pair_ids, weights),
        "weight_normalization_error": normalization_error,
        "source_n_targets": len(levels),
        "source_target_counts_json": json.dumps(counts, sort_keys=True),
        "source_target_weight_sums_json": json.dumps(target_sums, sort_keys=True),
        "max_target_weight_sum_error": max_target_error,
    }
    return weights, audit


def weighted_center_scale(
    values: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.ndim != 2 or len(x) == 0 or x.shape[0] != len(w):
        raise RuntimeError("Invalid matrix or weights for standardization")
    if not np.isfinite(x).all() or not np.isfinite(w).all() or (w <= 0).any():
        raise RuntimeError("Invalid values or weights for standardization")
    mean = np.sum(x * w[:, None], axis=0) / w.sum()
    variance = np.sum(w[:, None] * (x - mean) ** 2, axis=0) / w.sum()
    sd = np.sqrt(np.maximum(variance, 0.0))
    if not np.isfinite(sd).all() or not (sd > 0.0).all():
        raise RuntimeError("Zero or nonfinite training-only feature SD")
    return mean, sd


def logistic_design(
    values: np.ndarray, mean: np.ndarray, sd: np.ndarray
) -> np.ndarray:
    transformed = (np.asarray(values, dtype=float) - mean) / sd
    return np.column_stack([np.ones(len(transformed)), transformed])


def fit_pooled_logistic(
    values: np.ndarray,
    outcome: np.ndarray,
    weights: np.ndarray,
    penalty_lambda: float,
    *,
    label: str,
) -> dict:
    x = np.asarray(values, dtype=float)
    y = np.asarray(outcome, dtype=float)
    w = np.asarray(weights, dtype=float)
    if len(y) == 0 or len(y) != len(x) or len(y) != len(w):
        return {"status": "failed", "estimable": False, "failure_type": "shape_error"}
    if np.unique(y).size != 2:
        return {"status": "failed", "estimable": False, "failure_type": "single_class"}
    if not (
        np.isfinite(x).all()
        and np.isfinite(y).all()
        and np.isfinite(w).all()
        and (w > 0).all()
    ):
        return {"status": "failed", "estimable": False, "failure_type": "nonfinite_or_invalid_weight"}
    try:
        mean, sd = weighted_center_scale(x, w)
    except RuntimeError as error:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "zero_or_invalid_feature_sd",
            "message": str(error),
        }
    design = logistic_design(x, mean, sd)
    penalty_mask = np.zeros(design.shape[1], dtype=float)
    penalty_mask[1:] = 1.0
    prevalence = float(np.average(y, weights=w))
    prevalence = min(max(prevalence, 1e-8), 1.0 - 1e-8)
    intercept_start = math.log(prevalence / (1.0 - prevalence))
    starts = [
        np.r_[intercept_start, np.zeros(design.shape[1] - 1)],
        np.r_[intercept_start, np.full(design.shape[1] - 1, 0.25)],
        np.r_[intercept_start, np.full(design.shape[1] - 1, -0.25)],
    ]

    def objective(theta: np.ndarray) -> float:
        eta = design @ theta
        nll = float(np.sum(w * rq3.stable_nll(eta, y)) / w.sum())
        return nll + 0.5 * float(penalty_lambda) * float(
            np.sum(theta * theta * penalty_mask)
        )

    def gradient(theta: np.ndarray) -> np.ndarray:
        eta = design @ theta
        return (
            design.T @ (w * (expit(eta) - y)) / w.sum()
            + float(penalty_lambda) * theta * penalty_mask
        )

    candidates = []
    for index, start in enumerate(starts):
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            jac=gradient,
            bounds=[(-40.0, 40.0)] * design.shape[1],
            options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-9, "maxls": 100},
        )
        candidates.append((float(result.fun), index, result))
    objective_value, _, result = min(candidates, key=lambda item: (item[0], item[1]))
    theta = np.asarray(result.x, dtype=float)
    newton_converged = False
    newton_iterations = 0
    for newton_iterations in range(1, 101):
        probabilities = expit(design @ theta)
        grad = (
            design.T @ (w * (probabilities - y)) / w.sum()
            + float(penalty_lambda) * theta * penalty_mask
        )
        if float(np.max(np.abs(grad))) <= 1e-10:
            newton_converged = True
            break
        information = (
            design.T * (w * probabilities * (1.0 - probabilities))
        ) @ design / w.sum()
        information[np.diag_indices_from(information)] += (
            float(penalty_lambda) * penalty_mask
        )
        try:
            direction = np.linalg.solve(information, -grad)
        except np.linalg.LinAlgError:
            break
        current = objective(theta)
        step = 1.0
        accepted = False
        while step >= 1e-12:
            candidate = theta + step * direction
            if np.all(np.abs(candidate) < 39.999999):
                candidate_value = objective(candidate)
                if np.isfinite(candidate_value) and candidate_value <= current + 1e-15:
                    theta = candidate
                    objective_value = float(candidate_value)
                    accepted = True
                    break
            step *= 0.5
        if not accepted:
            break
    probabilities = expit(design @ theta)
    information = (design.T * (w * probabilities * (1.0 - probabilities))) @ design / w.sum()
    information[np.diag_indices_from(information)] += float(penalty_lambda) * penalty_mask
    eigenvalues = np.linalg.eigvalsh((information + information.T) / 2.0)
    information_ok = bool(
        np.isfinite(eigenvalues).all()
        and float(eigenvalues.min())
        > np.finfo(float).eps * 100.0 * max(1.0, float(eigenvalues.max()))
    )
    final_gradient = gradient(theta)
    bound_hit = bool(np.any(np.isclose(np.abs(theta), 40.0, atol=1e-6, rtol=0.0)))
    estimable = bool(
        (bool(result.success) or newton_converged)
        and np.isfinite(theta).all()
        and np.isfinite(objective_value)
        and np.isfinite(final_gradient).all()
        and information_ok
        and not bound_hit
        and float(np.max(np.abs(final_gradient))) <= 1e-5
    )
    failure_type = ""
    if not bool(result.success) and not newton_converged:
        failure_type = "optimizer_nonconvergence"
    elif not information_ok:
        failure_type = "singular_or_ill_conditioned_information"
    elif bound_hit:
        failure_type = "parameter_bound_hit"
    elif float(np.max(np.abs(final_gradient))) > 1e-5:
        failure_type = "gradient_not_small"
    return {
        "status": "evaluable" if estimable else "failed",
        "estimable": estimable,
        "failure_type": failure_type,
        "message": label,
        "theta": theta,
        "intercept": float(theta[0]),
        "slopes_standardized": theta[1:].copy(),
        "feature_mean": mean,
        "feature_sd": sd,
        "objective": objective_value,
        "gradient_max_abs": float(np.max(np.abs(final_gradient))),
        "optimizer_success": bool(result.success),
        "newton_converged": newton_converged,
        "newton_iterations": int(newton_iterations),
        "information_ok": information_ok,
    }


def logistic_eta(values: np.ndarray, fit: dict) -> np.ndarray:
    mean = np.asarray(fit["feature_mean"], dtype=float)
    sd = np.asarray(fit["feature_sd"], dtype=float)
    design = logistic_design(np.asarray(values, dtype=float), mean, sd)
    return design @ np.r_[float(fit["intercept"]), np.asarray(fit["slopes_standardized"], dtype=float)]


def effective_rf_grid(n_features: int) -> list[dict]:
    grid = []
    seen = set()
    for max_depth in (2, 4, None):
        for min_samples_leaf in (10, 30):
            for requested_max_features in (1, 2):
                effective = min(requested_max_features, int(n_features))
                key = (max_depth, min_samples_leaf, effective)
                if key in seen:
                    continue
                seen.add(key)
                grid.append(
                    {
                        "max_depth": max_depth,
                        "min_samples_leaf": min_samples_leaf,
                        "max_features": effective,
                        "requested_max_features": requested_max_features,
                    }
                )
    return grid


def fit_random_forest(
    values: np.ndarray,
    outcome: np.ndarray,
    weights: np.ndarray,
    parameters: dict,
    *,
    seed: int,
    n_jobs: int,
    label: str,
) -> dict:
    x = np.asarray(values, dtype=float)
    y = np.asarray(outcome, dtype=int)
    w = np.asarray(weights, dtype=float)
    if len(y) == 0 or len(y) != len(x) or len(y) != len(w):
        return {"status": "failed", "estimable": False, "failure_type": "shape_error"}
    if not np.isin(y, [0, 1]).all() or np.unique(y).size != 2:
        return {"status": "failed", "estimable": False, "failure_type": "single_class"}
    if not (
        np.isfinite(x).all()
        and np.isfinite(w).all()
        and (w > 0).all()
    ):
        return {"status": "failed", "estimable": False, "failure_type": "nonfinite_or_invalid_weight"}
    estimator = RandomForestClassifier(
        n_estimators=500,
        criterion="gini",
        max_depth=parameters["max_depth"],
        min_samples_split=2,
        min_samples_leaf=int(parameters["min_samples_leaf"]),
        min_weight_fraction_leaf=0.0,
        max_features=int(parameters["max_features"]),
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=True,
        oob_score=False,
        n_jobs=int(n_jobs),
        random_state=int(seed),
        verbose=0,
        warm_start=False,
        class_weight=None,
        max_samples=None,
        ccp_alpha=0.0,
    )
    try:
        estimator.fit(x, y, sample_weight=w)
    except Exception as error:  # Preserve every estimator failure in the ledger.
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "fit_exception",
            "message": f"{label}: {type(error).__name__}: {error}",
        }
    classes = np.asarray(estimator.classes_)
    if not (classes.shape == (2,) and np.array_equal(classes, np.asarray([0, 1]))):
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "invalid_classes_mapping",
            "classes": classes.tolist(),
        }
    sample_lists = getattr(estimator, "estimators_samples_", None)
    if sample_lists is None or len(sample_lists) != 500:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "missing_estimators_samples",
        }
    in_bag_counts = np.asarray(
        [len(np.unique(np.asarray(sample, dtype=int))) for sample in sample_lists],
        dtype=float,
    )
    return {
        "status": "evaluable",
        "estimable": True,
        "failure_type": "",
        "message": label,
        "estimator": estimator,
        "classes": [0, 1],
        "n_in_bag": float(in_bag_counts.mean()),
        "n_in_bag_min": float(in_bag_counts.min()),
        "n_in_bag_max": float(in_bag_counts.max()),
        "seed": int(seed),
    }


def rf_positive_probability(estimator: Any, values: np.ndarray) -> np.ndarray:
    classes = np.asarray(estimator.classes_)
    if classes.shape != (2,) or not np.array_equal(classes, np.asarray([0, 1])):
        raise RuntimeError("RF classes_ must be exactly [0,1]")
    positive_index = int(np.flatnonzero(classes == 1)[0])
    probability = np.asarray(estimator.predict_proba(np.asarray(values, dtype=float))[:, positive_index], dtype=float)
    if not np.isfinite(probability).all():
        raise RuntimeError("RF returned nonfinite probability")
    return probability


def selection_tie_value(model_family: str, candidate: dict) -> tuple:
    max_depth_rank = int(candidate["depth_rank"])
    return (
        max_depth_rank,
        -int(candidate["min_samples_leaf"]),
        int(candidate["max_features"]),
    )


def minimum_selection(
    records: Sequence[dict], model_family: str
) -> float:
    evaluable = [
        record
        for record in records
        if record["candidate_status"] == "evaluable"
        and np.isfinite(float(record["macro_mean_inner_log_loss"]))
    ]
    if not evaluable:
        return float("nan")
    best = min(float(record["macro_mean_inner_log_loss"]) for record in evaluable)
    ties = [
        record
        for record in evaluable
        if float(record["macro_mean_inner_log_loss"]) - best <= 1e-12
    ]
    if model_family == "pooled_logistic":
        return float(max(float(record["lambda"]) for record in ties))
    selected = min(ties, key=lambda record: selection_tie_value(model_family, record))
    return float(selected["candidate_id"])


def log_loss(y: np.ndarray, probability: np.ndarray, clip: float) -> float:
    yy = np.asarray(y, dtype=float)
    pp = np.clip(np.asarray(probability, dtype=float), float(clip), 1.0 - float(clip))
    if not np.isfinite(pp).all():
        return float("nan")
    return float(-np.mean(yy * np.log(pp) + (1.0 - yy) * np.log1p(-pp)))


def brier(y: np.ndarray, probability: np.ndarray) -> float:
    return float(np.mean((np.asarray(probability, dtype=float) - np.asarray(y, dtype=float)) ** 2))


def bss(y: np.ndarray, probability: np.ndarray, reference: np.ndarray) -> float:
    return rq3.bss_from_reference(np.asarray(y), np.asarray(probability), np.asarray(reference))


def safe_auc(y: np.ndarray, probability: np.ndarray) -> float:
    if np.unique(np.asarray(y)).size < 2 or not np.isfinite(np.asarray(probability)).all():
        return float("nan")
    return float(roc_auc_score(np.asarray(y), np.asarray(probability)))


def safe_ap(y: np.ndarray, probability: np.ndarray) -> float:
    if np.sum(np.asarray(y)) <= 0 or not np.isfinite(np.asarray(probability)).all():
        return float("nan")
    return float(average_precision_score(np.asarray(y), np.asarray(probability)))


def performance_metrics(
    y: np.ndarray,
    probability: np.ndarray,
    *,
    clip: float,
    transport_reference: np.ndarray,
    target_oracle_reference: np.ndarray,
    adaptation_reference: np.ndarray | None = None,
) -> dict:
    yy = np.asarray(y, dtype=float)
    pp = np.asarray(probability, dtype=float)
    result = {
        "n": int(len(yy)),
        "positives": int(np.sum(yy)),
        "prevalence": float(np.mean(yy)),
        "roc_auc": safe_auc(yy, pp),
        "ap": safe_ap(yy, pp),
        "brier": brier(yy, pp) if np.isfinite(pp).all() else float("nan"),
        "log_loss": log_loss(yy, pp, clip),
        "bss_transport": bss(yy, pp, transport_reference) if np.isfinite(pp).all() else float("nan"),
        "bss_target_oracle": bss(yy, pp, target_oracle_reference) if np.isfinite(pp).all() else float("nan"),
        "bss_adaptation_train": (
            bss(yy, pp, adaptation_reference)
            if adaptation_reference is not None and np.isfinite(pp).all()
            else float("nan")
        ),
    }
    if np.isfinite(pp).all() and np.unique(yy).size == 2:
        result.update(rq3.fit_calibration(yy, pp, clip))
    else:
        result.update(
            {
                "calibration_intercept": float("nan"),
                "calibration_slope": float("nan"),
                "calibration_intercept_se": float("nan"),
                "calibration_slope_se": float("nan"),
                "calibration_status": "failed",
                "calibration_failure_type": "single_class_or_nonfinite",
                "calibration_message": "Calibration diagnostic not evaluable",
            }
        )
    return result


def reliability_rows(
    y: np.ndarray, probability: np.ndarray, edges: Sequence[float]
) -> list[dict]:
    frame = pd.DataFrame(
        {
            "outcome_value": np.asarray(y, dtype=int),
            "predicted_probability": np.asarray(probability, dtype=float),
        }
    )
    return rq3.reliability_rows(frame, "predicted_probability", edges)


def solve_offset_intercept(offset: np.ndarray, outcome: np.ndarray, label: str) -> dict:
    return rq3.solve_offset_intercept(
        np.asarray(offset, dtype=float), np.asarray(outcome, dtype=float), label
    )


def fit_offset_slope(
    offset: np.ndarray,
    outcome: np.ndarray,
    *,
    penalty_lambda: float,
    label: str,
) -> dict:
    return rq3.fit_offset_slope(
        np.asarray(offset, dtype=float),
        np.asarray(outcome, dtype=float),
        penalty_lambda=float(penalty_lambda),
        label=label,
    )


def adaptation_inner_splits(
    train: pd.DataFrame, *, outcome_column: str = "outcome_value"
) -> tuple[list[tuple[np.ndarray, np.ndarray]], str]:
    ordered = stable_sort(train)
    splitter = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=20261113)
    splits = list(
        splitter.split(
            np.zeros((len(ordered), 1)),
            ordered[outcome_column].astype(int).to_numpy(),
            groups=ordered["sequence_cluster_value"].astype(str).to_numpy(),
        )
    )
    if len(splits) != 3:
        raise RuntimeError("StratifiedGroupKFold did not produce three folds")
    assignment_payload = []
    for number, (train_index, test_index) in enumerate(splits, start=1):
        train_groups = set(ordered.iloc[train_index]["sequence_cluster_value"].astype(str))
        test_groups = set(ordered.iloc[test_index]["sequence_cluster_value"].astype(str))
        if train_groups & test_groups:
            raise RuntimeError("Adaptation inner split has cluster leakage")
        if ordered.iloc[train_index]["outcome_value"].nunique() < 2:
            raise RuntimeError("Adaptation inner split has single-class training")
        if ordered.iloc[test_index]["outcome_value"].nunique() < 1:
            raise RuntimeError("Adaptation inner split has no validation rows")
        assignment_payload.append(
            {
                "inner_fold_number": number,
                "train_pair_ids": sorted(ordered.iloc[train_index]["pair_id"].astype(str)),
                "test_pair_ids": sorted(ordered.iloc[test_index]["pair_id"].astype(str)),
            }
        )
    return splits, canonical_json_sha256(assignment_payload)


def source_reference(
    train: pd.DataFrame, weighting: str
) -> tuple[float, dict[str, float]]:
    prevalence = train.groupby("target", sort=True)["outcome_value"].mean().astype(float)
    if weighting == "target_equal":
        reference = float(prevalence.mean())
    elif weighting == "pair":
        reference = float(train["outcome_value"].mean())
    else:
        raise ValueError(f"Unknown weighting: {weighting}")
    return reference, prevalence.to_dict()


def source_parameter_row(
    *,
    base: dict,
    exposure: str,
    family: str,
    fit: dict,
    train: pd.DataFrame,
    weight_audit: dict,
    selected_lambda: float,
    seed: int | None,
    parameters: dict | None,
) -> dict:
    features = FEATURE_MAP[exposure]
    row = {
        **base,
        **weight_audit,
        "model_exposure": exposure,
        "model_family": family,
        "features_json": json.dumps(features),
        "fit_status": fit.get("status", "failed"),
        "failure_type": fit.get("failure_type", ""),
        "selected_lambda": selected_lambda,
        "seed": seed if seed is not None else np.nan,
        "training_count": int(len(train)),
        "training_positives": int(train["outcome_value"].sum()),
        "classes_json": json.dumps(fit.get("classes", [])),
        "n_in_bag": fit.get("n_in_bag", np.nan),
        "n_in_bag_min": fit.get("n_in_bag_min", np.nan),
        "n_in_bag_max": fit.get("n_in_bag_max", np.nan),
        "objective": fit.get("objective", np.nan),
        "gradient_max_abs": fit.get("gradient_max_abs", np.nan),
        "parameters_json": json.dumps(parameters or {}, sort_keys=True),
    }
    if family == "pooled_logistic":
        row.update(
            {
                "intercept": fit.get("intercept", np.nan),
                "slopes_standardized_json": json.dumps(
                    np.asarray(
                        fit.get("slopes_standardized", []), dtype=float
                    ).tolist(),
                    separators=(",", ":"),
                ),
                "feature_mean_json": json.dumps(
                    np.asarray(fit.get("feature_mean", []), dtype=float).tolist(),
                    separators=(",", ":"),
                ),
                "feature_sd_json": json.dumps(
                    np.asarray(fit.get("feature_sd", []), dtype=float).tolist(),
                    separators=(",", ":"),
                ),
            }
        )
    return row


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
