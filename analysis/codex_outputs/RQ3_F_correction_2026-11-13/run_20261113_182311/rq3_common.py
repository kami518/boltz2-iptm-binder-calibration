"""Shared, deterministic utilities for the frozen RQ3 execution.

All writes are confined to this execution run directory. The module implements
the linear and restricted cubic spline source models with NumPy/SciPy only.
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

RUN_DIR = Path(__file__).resolve().parent
ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
HANDOFF_DIR = ROOT / "analysis" / "codex_tasks" / "RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13"
CONFIG_PATH = HANDOFF_DIR / "rq3_execution_config.json"
STATUS_PATH = RUN_DIR / "RUN_STATUS.json"
FAILURE_PATH = RUN_DIR / "failed_item_ledger.csv"

ANALYSIS_SETS = ("primary_six", "primary_five", "strict_six", "strict_five")
CLUSTER_RULES = ("sequence_cluster90", "sequence_cluster80")
SOURCE_WEIGHTS = ("target_equal", "pair")
LINEAR = "linear"
RCS = "restricted_cubic_spline"
U_LAYER = "U_transport_zero_label"
ORACLE_LAYER = "oracle_base_rate_scenario"
DESCRIPTIVE_LAYER = "descriptive_in_sample_adaptation"
PRED_I_LAYER = "predictive_I"
PRED_F_LAYER = "predictive_F"


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


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


def update_status(**updates: Any) -> dict:
    if STATUS_PATH.exists():
        with STATUS_PATH.open("r", encoding="utf-8") as handle:
            status = json.load(handle)
    else:
        config = load_config()
        status = {
            "handoff_id": config["handoff_id"],
            "run_directory": str(RUN_DIR),
            "overall_status": "initializing",
            "created_at": now_iso(),
            "config_frozen_before_any_rq3_result": True,
            "stages": {
                name: {"status": "not_started"}
                for name in config["execution_stages"]
            },
        }
    status.update(updates)
    status["updated_at"] = now_iso()
    atomic_json(STATUS_PATH, status)
    return status


def set_stage_status(stage: str, status: str, **details: Any) -> None:
    current = update_status()
    stages = current.setdefault("stages", {})
    stages[stage] = {"status": status, "updated_at": now_iso(), **details}
    if status == "running":
        overall = "running"
    else:
        values = [item.get("status") for item in stages.values()]
        if any(value in {"failed", "blocked"} for value in values):
            overall = "partial"
        elif all(value == "complete" for value in values):
            overall = "complete"
        else:
            overall = "partial"
    update_status(
        current_stage=stage,
        stages=stages,
        overall_status=overall,
    )


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
    frame = pd.DataFrame([row])
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
        frame = pd.concat([existing, frame], ignore_index=True)
    atomic_csv(FAILURE_PATH, frame)


def environment_snapshot() -> dict:
    return {
        "captured_at": now_iso(),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "dependency_install_or_upgrade_performed": False,
    }


def verify_input_hashes(config: dict) -> tuple[pd.DataFrame, list[str]]:
    rows = []
    mismatches = []
    specs = {"protocol": config["protocol"], **config["inputs"]}
    for key, spec in specs.items():
        path = Path(spec["path"])
        if not path.is_file():
            actual = ""
            rows.append(
                {
                    "input_key": key,
                    "path": str(path),
                    "expected_sha256": spec["sha256"],
                    "actual_sha256": actual,
                    "sha256_matches": False,
                    "bytes": 0,
                }
            )
            mismatches.append(key)
            continue
        actual = sha256_file(path)
        matches = actual.lower() == spec["sha256"].lower()
        rows.append(
            {
                "input_key": key,
                "path": str(path),
                "expected_sha256": spec["sha256"],
                "actual_sha256": actual,
                "sha256_matches": matches,
                "bytes": path.stat().st_size,
            }
        )
        if not matches:
            mismatches.append(key)
    return pd.DataFrame(rows), mismatches


def require_hash_gate(config: dict) -> pd.DataFrame:
    hashes, mismatches = verify_input_hashes(config)
    if mismatches:
        append_failure(
            "preflight",
            "frozen_input_hash_gate",
            "sha256_mismatch",
            "Frozen input hash gate failed; no analysis should proceed.",
            mismatches=mismatches,
        )
        raise RuntimeError(f"Frozen input SHA-256 mismatch: {mismatches}")
    return hashes


def load_joined_data(config: dict) -> pd.DataFrame:
    pair_path = Path(config["inputs"]["pair_master"]["path"])
    cluster_path = Path(config["inputs"]["sequence_clusters"]["path"])
    pairs = pd.read_csv(pair_path, low_memory=False)
    clusters = pd.read_csv(cluster_path, low_memory=False)
    if len(pairs) != int(config["inputs"]["pair_master"]["expected_rows"]):
        raise RuntimeError("pair_master_v3.csv row count mismatch")
    if len(clusters) != int(config["inputs"]["sequence_clusters"]["expected_rows"]):
        raise RuntimeError("sequence_clusters_rq12.csv row count mismatch")
    if pairs["pair_id"].duplicated().any() or clusters["pair_id"].duplicated().any():
        raise RuntimeError("pair_id is not unique in a frozen input")
    joined = pairs.merge(
        clusters[
            ["pair_id", "target", "design_id", "sequence", *CLUSTER_RULES]
        ],
        on="pair_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_cluster_source"),
    )
    if len(joined) != len(pairs):
        raise RuntimeError("Pair/cluster inner join lost rows")
    for column in ("target", "design_id", "sequence"):
        if not joined[column].astype(str).eq(
            joined[f"{column}_cluster_source"].astype(str)
        ).all():
            raise RuntimeError(f"Frozen inputs disagree on {column}")
    joined = joined.drop(
        columns=[
            "target_cluster_source",
            "design_id_cluster_source",
            "sequence_cluster_source",
        ]
    )
    joined["logit_iptm"] = logit_clip(
        joined["boltz2_iptm"].to_numpy(dtype=float),
        float(config["probability_clip"]),
    )
    return joined


def logit_clip(probability: np.ndarray, clip: float) -> np.ndarray:
    p = np.asarray(probability, dtype=float)
    clipped = np.clip(p, clip, 1.0 - clip)
    return np.log(clipped) - np.log1p(-clipped)


def validate_fold_structure(
    config: dict, joined: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    input_paths = config["inputs"]
    loto = pd.read_csv(input_paths["loto_assignments"]["path"], low_memory=False)
    loto_manifest = pd.read_csv(input_paths["loto_manifest"]["path"], low_memory=False)
    cv5 = pd.read_csv(input_paths["target_cv5_assignments"]["path"], low_memory=False)
    cv5_manifest = pd.read_csv(
        input_paths["target_cv5_manifest"]["path"], low_memory=False
    )
    exclusions = pd.read_csv(
        input_paths["cross_target_exclusions"]["path"], low_memory=False
    )

    pair_ids = set(joined["pair_id"].astype(str))
    checks: list[dict] = []

    def check(name: str, passed: bool, observed: Any, expected: Any) -> None:
        checks.append(
            {
                "check": name,
                "passed": bool(passed),
                "observed": observed,
                "expected": expected,
            }
        )
        if not passed:
            raise RuntimeError(f"Stage 0 check failed: {name}; observed={observed}")

    check("master_pair_count", len(joined) == 2084, len(joined), 2084)
    check("unique_pair_ids", joined["pair_id"].nunique() == 2084, joined["pair_id"].nunique(), 2084)
    check("cluster_pair_count", len(joined) == 2084, len(joined), 2084)
    check("loto_manifest_rows", len(loto_manifest) == 44, len(loto_manifest), 44)
    check("cv5_manifest_rows", len(cv5_manifest) == 220, len(cv5_manifest), 220)
    check("loto_manifest_qc", loto_manifest["qc_passed"].all(), int(loto_manifest["qc_passed"].sum()), 44)
    check("cv5_manifest_qc", cv5_manifest["qc_passed"].all(), int(cv5_manifest["qc_passed"].sum()), 220)
    check(
        "analysis_set_names",
        set(loto_manifest["analysis_set"]) == set(ANALYSIS_SETS),
        sorted(set(loto_manifest["analysis_set"])),
        sorted(ANALYSIS_SETS),
    )
    check(
        "cluster_rule_names",
        set(loto_manifest["cluster_rule"]) == set(CLUSTER_RULES),
        sorted(set(loto_manifest["cluster_rule"])),
        sorted(CLUSTER_RULES),
    )
    check(
        "all_assignment_pair_ids",
        set(loto["pair_id"].astype(str)) == pair_ids
        and set(cv5["pair_id"].astype(str)) == pair_ids,
        {
            "loto": len(set(loto["pair_id"].astype(str))),
            "cv5": len(set(cv5["pair_id"].astype(str))),
        },
        len(pair_ids),
    )
    check(
        "loto_expected_rows_per_fold",
        loto.groupby("fold_id").size().to_dict()
        == loto_manifest.set_index("fold_id")["n_pairs_total"].to_dict(),
        "manifest-consistent",
        "manifest-consistent",
    )
    check(
        "cv5_expected_rows_per_fold",
        cv5.groupby("fold_id").size().to_dict()
        == cv5_manifest.set_index("fold_id")["n_pairs_total"].to_dict(),
        "manifest-consistent",
        "manifest-consistent",
    )

    loto_leaks = []
    loto_test_incomplete = []
    for row in loto_manifest.itertuples(index=False):
        fold = loto[loto["fold_id"].eq(row.fold_id)]
        target = str(row.held_out_target)
        train = fold[fold["split_status"].eq("train")]
        test = fold[fold["split_status"].eq("test")]
        removed = fold[fold["split_status"].eq("removed_cross_target_cluster")]
        expected_test = joined[joined["target"].astype(str).eq(target)]
        expected_test = expected_test[
            expected_test["pair_id"].isin(
                joined.loc[
                    joined["target"].astype(str).ne(target), "pair_id"
                ].pipe(lambda s: joined["pair_id"])
            )
        ]
        if train["target"].astype(str).eq(target).any():
            loto_leaks.append(row.fold_id)
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            loto_leaks.append(row.fold_id)
        if len(test) != int(row.n_test_pairs):
            loto_test_incomplete.append(row.fold_id)
        if len(train) != int(row.n_train_pairs_final):
            loto_test_incomplete.append(row.fold_id)
        if len(removed) != int(row.n_removed_cross_target_cluster):
            loto_test_incomplete.append(row.fold_id)
        if set(test["target"].astype(str)) != {target}:
            loto_test_incomplete.append(row.fold_id)
    check("loto_training_target_absent", not loto_leaks, len(loto_leaks), 0)
    check("loto_test_complete", not loto_test_incomplete, sorted(set(loto_test_incomplete)), [])

    cv_leaks = []
    cv_test_once = []
    for row in cv5_manifest.itertuples(index=False):
        fold = cv5[cv5["fold_id"].eq(row.fold_id)]
        train = fold[fold["split"].eq("train")]
        test = fold[fold["split"].eq("test")]
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            cv_leaks.append(row.fold_id)
        if len(train) != int(row.n_train_pairs) or len(test) != int(row.n_test_pairs):
            cv_test_once.append(row.fold_id)
    test_counts = cv5[cv5["split"].eq("test")].groupby(
        ["analysis_set", "outcome_name", "cluster_rule", "target", "pair_id"]
    ).size()
    check("cv5_group_isolation", not cv_leaks, len(cv_leaks), 0)
    check("cv5_fold_counts", not cv_test_once, sorted(set(cv_test_once)), [])
    check("cv5_each_pair_tested_once", test_counts.eq(1).all(), int(test_counts.eq(1).sum()), len(test_counts))
    check(
        "exclusion_ledger_nonempty_and_accounted",
        len(exclusions) == int(loto["split_status"].eq("removed_cross_target_cluster").sum()),
        len(exclusions),
        int(loto["split_status"].eq("removed_cross_target_cluster").sum()),
    )

    # One row per expected analysis-set x cluster-rule x target combination.
    expected_loto = set()
    for set_name in ANALYSIS_SETS:
        for rule in CLUSTER_RULES:
            manifest_subset = loto_manifest[
                loto_manifest["analysis_set"].eq(set_name)
                & loto_manifest["cluster_rule"].eq(rule)
            ]
            for target in manifest_subset["held_out_target"]:
                expected_loto.add((set_name, rule, target))
    actual_loto = set(
        loto_manifest[["analysis_set", "cluster_rule", "held_out_target"]]
        .itertuples(index=False, name=None)
    )
    check("loto_combination_complete", actual_loto == expected_loto, len(actual_loto), len(expected_loto))

    summary = {
        "checks": checks,
        "all_passed": all(item["passed"] for item in checks),
        "n_checks": len(checks),
        "failure_count": 0,
        "note": "Tests above are structural and do not fit or evaluate any model.",
    }
    return loto, loto_manifest, cv5, cv5_manifest, exclusions, summary


def model_configurations() -> Iterable[dict]:
    for analysis_set in ANALYSIS_SETS:
        for cluster_rule in CLUSTER_RULES:
            for source_weighting in SOURCE_WEIGHTS:
                yield {
                    "analysis_set": analysis_set,
                    "cluster_rule": cluster_rule,
                    "source_weighting": source_weighting,
                }


def stable_nll(eta: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.logaddexp(0.0, eta) - np.asarray(y, dtype=float) * eta


def _information_qc(hessian: np.ndarray, label: str) -> tuple[bool, str]:
    if not np.isfinite(hessian).all():
        return False, f"{label}: nonfinite Hessian"
    eigenvalues = np.linalg.eigvalsh((hessian + hessian.T) / 2.0)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    if eigenvalues.min() <= np.finfo(float).eps * 100.0 * scale:
        return False, f"{label}: singular/unidentifiable information; eig_min={eigenvalues.min():.6g}"
    condition = float(np.linalg.cond(hessian))
    if not math.isfinite(condition) or condition > 1e14:
        return False, f"{label}: ill-conditioned information; cond={condition:.6g}"
    return True, f"{label}: Hessian finite, eig_min={eigenvalues.min():.6g}, cond={condition:.6g}"


def fit_weighted_logistic(
    design: np.ndarray,
    outcome: np.ndarray,
    weights: np.ndarray,
    *,
    label: str,
    starts: Sequence[np.ndarray],
    penalty_lambda: float | None = None,
    delta_index: int | None = None,
    bound: float = 40.0,
    offset: np.ndarray | None = None,
) -> dict:
    X = np.asarray(design, dtype=float)
    y = np.asarray(outcome, dtype=float)
    w = np.asarray(weights, dtype=float)
    offset_values = (
        np.zeros(X.shape[0], dtype=float)
        if offset is None
        else np.asarray(offset, dtype=float)
    )
    if X.ndim != 2 or X.shape[0] != y.size or y.size != w.size:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "shape_error",
            "message": f"{label}: design/outcome/weight shape mismatch",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
        }
    if (
        offset_values.ndim != 1
        or offset_values.size != X.shape[0]
        or not (
            np.isfinite(X).all()
            and np.isfinite(y).all()
            and np.isfinite(w).all()
            and np.isfinite(offset_values).all()
        )
    ):
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "nonfinite_input",
            "message": f"{label}: invalid or nonfinite design/outcome/weight/offset",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
        }
    if not np.isin(y, [0.0, 1.0]).all() or np.unique(y).size < 2:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "single_class",
            "message": f"{label}: outcome must contain both classes",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
        }
    if (w <= 0).any():
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "invalid_weight",
            "message": f"{label}: all fitting weights must be positive",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
        }
    weight_sum = float(w.sum())
    if not math.isfinite(weight_sum) or weight_sum <= 0:
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "invalid_weight_sum",
            "message": f"{label}: invalid weight sum",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
        }

    def objective(theta: np.ndarray) -> float:
        eta = X @ theta + offset_values
        value = float(np.sum(w * stable_nll(eta, y)) / weight_sum)
        if penalty_lambda is not None and delta_index is not None:
            value += 0.5 * float(penalty_lambda) * float(theta[delta_index] ** 2)
        return value

    def gradient(theta: np.ndarray) -> np.ndarray:
        eta = X @ theta + offset_values
        grad = X.T @ (w * (expit(eta) - y)) / weight_sum
        if penalty_lambda is not None and delta_index is not None:
            grad = grad.copy()
            grad[delta_index] += float(penalty_lambda) * theta[delta_index]
        return grad

    bounds = [(-bound, bound)] * X.shape[1]
    best = None
    start_records = []
    for index, start in enumerate(starts):
        initial = np.asarray(start, dtype=float)
        if initial.shape != (X.shape[1],) or not np.isfinite(initial).all():
            continue
        initial = np.clip(initial, -bound + 1e-9, bound - 1e-9)
        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=gradient,
            bounds=bounds,
            options={
                "maxiter": 5000,
                "ftol": 1e-12,
                "gtol": 1e-8,
                "maxls": 100,
            },
        )
        value = float(result.fun)
        start_records.append(
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
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "no_valid_start",
            "message": f"{label}: no finite deterministic start",
            "theta": np.full(X.shape[1], np.nan),
            "covariance": np.full((X.shape[1], X.shape[1]), np.nan),
            "starts": start_records,
        }
    objective_value, result = best
    theta = np.asarray(result.x, dtype=float)
    probabilities = expit(X @ theta + offset_values)
    hessian = (X.T * (w * probabilities * (1.0 - probabilities))) @ X / weight_sum
    if penalty_lambda is not None and delta_index is not None:
        hessian = hessian.copy()
        hessian[delta_index, delta_index] += float(penalty_lambda)
    information_ok, information_message = _information_qc(hessian, label)
    bound_hit = bool(np.any(np.isclose(np.abs(theta), bound, atol=1e-6, rtol=0.0)))
    final_gradient = gradient(theta)
    finite = bool(
        np.isfinite(theta).all()
        and np.isfinite(objective_value)
        and np.isfinite(final_gradient).all()
    )
    estimable = bool(
        result.success
        and finite
        and information_ok
        and not bound_hit
        and np.max(np.abs(final_gradient)) <= 1e-5
    )
    failure_type = None
    if not result.success:
        failure_type = "optimizer_nonconvergence"
    elif not finite:
        failure_type = "nonfinite_solution"
    elif bound_hit:
        failure_type = "parameter_bound_hit"
    elif not information_ok:
        failure_type = "singular_or_ill_conditioned_information"
    elif np.max(np.abs(final_gradient)) > 1e-5:
        failure_type = "gradient_not_small"
    covariance = (
        np.linalg.pinv(hessian)
        if np.isfinite(hessian).all()
        else np.full((X.shape[1], X.shape[1]), np.nan)
    )
    return {
        "status": "evaluable" if estimable else "failed",
        "estimable": estimable,
        "failure_type": failure_type,
        "message": (
            f"{label}: converged and passed numerical QC"
            if estimable
            else f"{label}: {failure_type}"
        ),
        "theta": theta,
        "covariance": covariance,
        "objective": objective_value,
        "gradient": final_gradient,
        "gradient_max_abs": float(np.max(np.abs(final_gradient))),
        "iterations": int(result.nit),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "bound_hit": bound_hit,
        "information_qc_message": information_message,
        "starts": start_records,
    }


def target_one_hot(targets: np.ndarray, levels: Sequence[str]) -> np.ndarray:
    values = np.asarray(targets).astype(str)
    return np.column_stack([values == level for level in levels]).astype(float)


def fit_source_linear(
    x: np.ndarray,
    y: np.ndarray,
    training_targets: np.ndarray,
    weights: np.ndarray,
    label: str,
) -> dict:
    levels = sorted(pd.unique(np.asarray(training_targets).astype(str)))
    design = np.column_stack([target_one_hot(training_targets, levels), x])
    start = np.zeros(design.shape[1])
    for index, level in enumerate(levels):
        mask = np.asarray(training_targets).astype(str) == level
        prevalence = np.clip(
            float(np.average(y[mask], weights=weights[mask])),
            1e-6,
            1.0 - 1e-6,
        )
        start[index] = math.log(prevalence / (1.0 - prevalence))
    starts = [
        start,
        np.r_[start[:-1], 0.25],
        np.r_[start[:-1], -0.25],
    ]
    fit = fit_weighted_logistic(
        design,
        y,
        weights,
        label=label,
        starts=starts,
    )
    fit.update({"design": design, "levels": levels})
    return fit


def weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantiles: Sequence[float]
) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    order = np.argsort(x, kind="mergesort")
    x_sorted = x[order]
    w_sorted = w[order]
    cumulative = (np.cumsum(w_sorted) - 0.5 * w_sorted) / np.sum(w_sorted)
    return np.interp(
        np.asarray(quantiles, dtype=float),
        cumulative,
        x_sorted,
        left=float(x_sorted.min()),
        right=float(x_sorted.max()),
    )


def rcs_basis(x: np.ndarray, knots: Sequence[float]) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    k1, k2, k3 = (float(value) for value in knots)
    if not k1 < k2 < k3:
        raise ValueError("RCS knots must be strictly increasing")
    c = (k3 - k1) / (k3 - k2)
    d = (k2 - k1) / (k3 - k2)
    denominator = (k3 - k1) ** 2
    b1 = values
    b2 = (
        np.maximum(values - k1, 0.0) ** 3
        - c * np.maximum(values - k2, 0.0) ** 3
        + d * np.maximum(values - k3, 0.0) ** 3
    ) / denominator
    return np.column_stack([b1, b2])


def rcs_second_derivative(x: np.ndarray, knots: Sequence[float]) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    k1, k2, k3 = (float(value) for value in knots)
    c = (k3 - k1) / (k3 - k2)
    d = (k2 - k1) / (k3 - k2)
    denominator = (k3 - k1) ** 2
    return (
        6.0
        * (
            np.maximum(values - k1, 0.0)
            - c * np.maximum(values - k2, 0.0)
            + d * np.maximum(values - k3, 0.0)
        )
        / denominator
    )


def fit_source_rcs(
    x: np.ndarray,
    y: np.ndarray,
    training_targets: np.ndarray,
    weights: np.ndarray,
    quantiles: Sequence[float],
    label: str,
) -> dict:
    knots = weighted_quantile(x, weights, quantiles)
    levels = sorted(pd.unique(np.asarray(training_targets).astype(str)))
    if not np.all(np.diff(knots) > 0):
        return {
            "status": "failed",
            "estimable": False,
            "failure_type": "non_strict_rcs_knots",
            "message": f"{label}: RCS knots are not strictly increasing",
            "theta": np.full(len(levels) + 2, np.nan),
            "covariance": np.full(
                (len(levels) + 2,) * 2, np.nan
            ),
            "knots": knots,
            "levels": levels,
        }
    basis = rcs_basis(x, knots)
    design = np.column_stack([target_one_hot(training_targets, levels), basis])
    start = np.zeros(design.shape[1])
    for index, level in enumerate(levels):
        mask = np.asarray(training_targets).astype(str) == level
        prevalence = np.clip(
            float(np.average(y[mask], weights=weights[mask])),
            1e-6,
            1.0 - 1e-6,
        )
        start[index] = math.log(prevalence / (1.0 - prevalence))
    starts = [
        start,
        np.r_[start[:-2], 0.25, 0.0],
        np.r_[start[:-2], -0.25, 0.0],
    ]
    fit = fit_weighted_logistic(
        design,
        y,
        weights,
        label=label,
        starts=starts,
    )
    fit.update({"design": design, "levels": levels, "knots": knots})
    return fit


def fit_calibration(
    outcome: np.ndarray,
    prediction: np.ndarray,
    clip: float,
) -> dict:
    q = logit_clip(np.asarray(prediction, dtype=float), clip)
    design = np.column_stack([np.ones(len(q)), q])
    starts = [np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.5])]
    fit = fit_weighted_logistic(
        design,
        np.asarray(outcome, dtype=float),
        np.ones(len(q)),
        label="calibration",
        starts=starts,
    )
    if not fit["estimable"]:
        return {
            "calibration_intercept": np.nan,
            "calibration_slope": np.nan,
            "calibration_intercept_se": np.nan,
            "calibration_slope_se": np.nan,
            "calibration_status": "failed",
            "calibration_failure_type": fit["failure_type"],
            "calibration_message": fit["message"],
        }
    covariance = fit["covariance"]
    return {
        "calibration_intercept": float(fit["theta"][0]),
        "calibration_slope": float(fit["theta"][1]),
        "calibration_intercept_se": float(math.sqrt(max(covariance[0, 0], 0.0))),
        "calibration_slope_se": float(math.sqrt(max(covariance[1, 1], 0.0))),
        "calibration_status": "evaluable",
        "calibration_failure_type": "",
        "calibration_message": fit["message"],
    }


def brier_score(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean((np.asarray(p, dtype=float) - np.asarray(y, dtype=float)) ** 2))


def log_loss_score(y: np.ndarray, p: np.ndarray, clip: float) -> float:
    yy = np.asarray(y, dtype=float)
    pp = np.clip(np.asarray(p, dtype=float), clip, 1.0 - clip)
    return float(-np.mean(yy * np.log(pp) + (1.0 - yy) * np.log1p(-pp)))


def bss_from_reference(y: np.ndarray, p: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.sum((np.asarray(reference) - np.asarray(y)) ** 2))
    if denominator <= 0 or not math.isfinite(denominator):
        return np.nan
    return float(1.0 - np.sum((np.asarray(p) - np.asarray(y)) ** 2) / denominator)


def compute_performance_metrics(
    y: np.ndarray,
    p: np.ndarray,
    *,
    clip: float,
    transport_reference: np.ndarray,
    target_oracle_reference: np.ndarray,
) -> dict:
    result = {
        "n": int(len(y)),
        "positives": int(np.sum(y)),
        "prevalence": float(np.mean(y)),
        "brier": brier_score(y, p),
        "log_loss": log_loss_score(y, p, clip),
        "bss_transport": bss_from_reference(y, p, transport_reference),
        "bss_target_oracle": bss_from_reference(y, p, target_oracle_reference),
    }
    if np.isfinite(p).all() and np.unique(y).size == 2:
        result.update(fit_calibration(y, p, clip))
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


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return np.nan, np.nan
    proportion = k / n
    denominator = 1.0 + z * z / n
    center = (proportion + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        proportion * (1.0 - proportion) / n + z * z / (4.0 * n * n)
    ) / denominator
    return center - half, center + half


def reliability_rows(
    frame: pd.DataFrame,
    probability_column: str,
    edges: Sequence[float],
) -> list[dict]:
    p = frame[probability_column].to_numpy(dtype=float)
    y = frame["outcome_value"].to_numpy(dtype=float)
    bin_indices = np.digitize(p, np.asarray(edges[1:-1], dtype=float), right=False)
    rows = []
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
        mask = bin_indices == index
        count = int(mask.sum())
        positives = int(np.sum(y[mask])) if count else 0
        low, high = wilson_interval(positives, count)
        rows.append(
            {
                "bin_index": index,
                "bin_lower": float(lower),
                "bin_upper": float(upper),
                "n": count,
                "positives": positives,
                "mean_prediction": float(np.mean(p[mask])) if count else np.nan,
                "observed_rate": float(np.mean(y[mask])) if count else np.nan,
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )
    return rows


def solve_offset_intercept(
    offset: np.ndarray, outcome: np.ndarray, label: str
) -> dict:
    y = np.asarray(outcome, dtype=float)
    eta = np.asarray(offset, dtype=float)
    if not (np.isfinite(eta).all() and np.isfinite(y).all()):
        return {
            "status": "failed",
            "failure_type": "nonfinite_input",
            "intercept": np.nan,
            "message": f"{label}: nonfinite input",
        }
    if np.unique(y).size < 2:
        return {
            "status": "failed",
            "failure_type": "single_class",
            "intercept": np.nan,
            "message": f"{label}: single-class outcome",
        }
    prevalence = float(np.mean(y))

    def root_function(intercept: float) -> float:
        return float(np.mean(expit(intercept + eta)) - prevalence)

    grid = np.linspace(-80.0, 80.0, 321)
    values = np.asarray([root_function(float(value)) for value in grid])
    brackets = np.where(np.sign(values[:-1]) * np.sign(values[1:]) <= 0)[0]
    if brackets.size == 0:
        return {
            "status": "failed",
            "failure_type": "root_not_bracketed",
            "intercept": np.nan,
            "message": f"{label}: intercept root not bracketed",
        }
    index = int(brackets[0])
    lower = float(grid[index])
    upper = float(grid[index + 1])
    if abs(values[index]) <= 1e-13:
        intercept = lower
    elif abs(values[index + 1]) <= 1e-13:
        intercept = upper
    else:
        for _ in range(200):
            middle = (lower + upper) / 2.0
            value = root_function(middle)
            if abs(value) <= 1e-13 or (upper - lower) <= 1e-12:
                intercept = middle
                break
            if np.sign(value) == np.sign(values[index]):
                lower = middle
            else:
                upper = middle
        else:
            intercept = (lower + upper) / 2.0
    residual = float(np.mean(expit(intercept + eta)) - prevalence)
    return {
        "status": "evaluable",
        "failure_type": "",
        "intercept": float(intercept),
        "residual": residual,
        "prevalence": prevalence,
        "message": f"{label}: base-rate intercept solved",
    }


def fit_offset_slope(
    offset: np.ndarray,
    outcome: np.ndarray,
    *,
    penalty_lambda: float,
    label: str,
) -> dict:
    eta_u = np.asarray(offset, dtype=float)
    y = np.asarray(outcome, dtype=float)
    # The frozen F parameterization is eta_F = a + (1 + delta) * eta_U.
    # eta_U is an explicit model offset: the design supplies [a, delta].
    design = np.column_stack([np.ones(len(eta_u)), eta_u])
    intercept_probe = solve_offset_intercept(eta_u, y, f"{label}|initial_intercept")
    initial_intercept = (
        float(intercept_probe["intercept"])
        if intercept_probe["status"] == "evaluable"
        else 0.0
    )
    starts = [
        np.array([initial_intercept, 0.0]),
        np.array([initial_intercept, 0.25]),
        np.array([initial_intercept, -0.25]),
    ]
    fit = fit_weighted_logistic(
        design,
        y,
        np.ones(len(y)),
        label=label,
        starts=starts,
        penalty_lambda=penalty_lambda,
        delta_index=1,
        offset=eta_u,
    )
    if fit["estimable"]:
        fit["intercept"] = float(fit["theta"][0])
        fit["delta"] = float(fit["theta"][1])
        fit["gamma"] = float(1.0 + fit["delta"])
    else:
        fit["intercept"] = np.nan
        fit["gamma"] = np.nan
        fit["delta"] = np.nan
    return fit


def macro_summary(
    frame: pd.DataFrame,
    metric_columns: Sequence[str],
    group_columns: Sequence[str],
) -> pd.DataFrame:
    rows = []
    for keys, group in frame.groupby(list(group_columns), dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_columns, keys))
        for metric in metric_columns:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            row = {
                **base,
                "summary_type": "target_equal_macro",
                "metric": metric,
                "mean": float(values.mean()) if len(values) else np.nan,
                "target_count": int(len(values)),
                "min": float(values.min()) if len(values) else np.nan,
                "max": float(values.max()) if len(values) else np.nan,
                "n_positive_sign": int((values > 0).sum()),
                "n_negative_sign": int((values < 0).sum()),
                "n_zero_sign": int((values == 0).sum()),
            }
            rows.append(row)
    return pd.DataFrame(rows)


def pair_summary_metrics(
    predictions: pd.DataFrame,
    *,
    probability_column: str,
    transport_column: str,
    oracle_column: str,
    clip: float,
) -> dict:
    return compute_performance_metrics(
        predictions["outcome_value"].to_numpy(dtype=float),
        predictions[probability_column].to_numpy(dtype=float),
        clip=clip,
        transport_reference=predictions[transport_column].to_numpy(dtype=float),
        target_oracle_reference=predictions[oracle_column].to_numpy(dtype=float),
    )


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


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def safe_json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
