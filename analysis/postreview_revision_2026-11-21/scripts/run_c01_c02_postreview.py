"""C01 ipSAE extension and C02 frozen-prediction bootstrap.

The script writes only under analysis/postreview_revision_2026-11-21 and never
rewrites a frozen manuscript, data, fold, Pinc, or model output.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis" / "postreview_revision_2026-11-21"
IPSAE_OUT = OUT / "ipSAE"
BOOT_OUT = OUT / "bootstrap"
SCRIPT_DIR = OUT / "scripts"
PYTHON = Path(sys.executable)
IPSAE_REPO = Path(os.environ.get("IPSAE_REPO", str(ROOT / "external" / "ipsae" / "IPSAE")))
IPSAE_SCRIPT = IPSAE_REPO / "ipsae.py"
IPSAE_COMMIT = "6174cf9e71cb1bd660cc805856a18c4871a6dec3"

PAIR_MASTER = ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv"
PAIR_AUDIT = ROOT / "analysis" / "audit_outputs" / "pair_audit.csv"
PINC = ROOT / "analysis" / "pinc_batch" / "run_20261110_181500" / "pinc_pair_extension.csv"
CLUSTERS = ROOT / "analysis" / "rq12_outputs" / "sequence_clusters_rq12.csv"
LOTO = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_prep_2026-11-13"
    / "run_20261113_105301"
    / "loto_pair_assignments.csv"
)
CV5 = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_prep_2026-11-13"
    / "run_20261113_105301"
    / "target_cv5_pair_assignments.csv"
)
PROTEINBASE = Path(os.environ.get("PROTEINBASE_CSV", str(ROOT / "data" / "proteinbase_latest.csv")))
RQ3_COMMON = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_F_correction_2026-11-13"
    / "run_20261113_182311"
    / "rq3_common.py"
)
MAIN_FROZEN_PREDICTIONS = (
    ROOT
    / "analysis"
    / "codex_outputs"
    / "RQ3_F_correction_2026-11-13"
    / "run_20261113_182311"
    / "rq3_adaptation_pair_predictions.csv"
)
EXTERNAL_PREDICTIONS = (
    ROOT
    / "analysis"
    / "rq3e_if_outputs"
    / "run_20261114_174403"
    / "oof_predictions.csv"
)

CLIP = 1e-6
BOOT_REPLICATES = 2000
BOOT_SEED = 20261121
LAMBDAS = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
LAYERS = ("U", "I", "F")
SCORE_MODELS = ("ipTM_only", "ipSAE_only")
WEIGHTINGS = ("target_equal", "pair")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def write_json(path: Path, payload) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(jsonable(payload), ensure_ascii=False, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    os.replace(tmp, path)


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False, encoding="utf-8")
    os.replace(tmp, path)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rq3 = load_module("postreview_rq3_common", RQ3_COMMON)


def logit_clip(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    clipped = np.clip(p, CLIP, 1.0 - CLIP)
    return np.log(clipped) - np.log1p(-clipped)


def safe_auc(y, score) -> float:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    if np.unique(y).size < 2 or not np.isfinite(score).all():
        return np.nan
    return float(roc_auc_score(y, score))


def safe_ap(y, score) -> float:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    if y.sum() <= 0 or not np.isfinite(score).all():
        return np.nan
    return float(average_precision_score(y, score))


def brier(y, p) -> float:
    return float(np.mean((np.asarray(p, float) - np.asarray(y, float)) ** 2))


def log_loss(y, p) -> float:
    y = np.asarray(y, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), CLIP, 1.0 - CLIP)
    return float(-np.mean(y * np.log(p) + (1.0 - y) * np.log1p(-p)))


def bss(y, p, reference) -> float:
    return float(1.0 - brier(y, p) / brier(y, reference))


def metric_set(y, p, transport_reference=None, target_reference=None) -> dict:
    result = {
        "n": int(len(y)),
        "positives": int(np.sum(y)),
        "prevalence": float(np.mean(y)),
        "roc_auc": safe_auc(y, p),
        "ap": safe_ap(y, p),
        "brier": brier(y, p),
        "log_loss": log_loss(y, p),
    }
    if transport_reference is not None:
        result["bss_transport"] = bss(y, p, transport_reference)
    if target_reference is not None:
        result["bss_target_oracle"] = bss(y, p, target_reference)
    return result


def source_weights(frame: pd.DataFrame, weighting: str) -> np.ndarray:
    column = "weight_target_equal" if weighting == "target_equal" else "weight_pair"
    weights = frame[column].to_numpy(dtype=float)
    if not np.isfinite(weights).all() or (weights <= 0).any():
        raise RuntimeError(f"Invalid {weighting} weights")
    return weights


def transport_intercept(
    levels: list[str],
    intercepts: np.ndarray,
    targets: np.ndarray,
    weighting: str,
) -> float:
    intercepts = np.asarray(intercepts, dtype=float)
    if weighting == "target_equal":
        return float(np.mean(intercepts))
    counts = pd.Series(np.asarray(targets, dtype=str)).value_counts()
    source_counts = np.asarray([float(counts.loc[level]) for level in levels])
    return float(np.sum(source_counts * intercepts) / source_counts.sum())


def source_reference(frame: pd.DataFrame, weighting: str) -> tuple[float, dict]:
    prevalences = frame.groupby("target", sort=True)["outcome_value"].mean().astype(float)
    reference = (
        float(prevalences.mean())
        if weighting == "target_equal"
        else float(frame["outcome_value"].mean())
    )
    return reference, prevalences.to_dict()


def fit_source(
    train: pd.DataFrame,
    x_column: str,
    weighting: str,
    label: str,
) -> dict:
    fit = rq3.fit_source_linear(
        train[x_column].to_numpy(dtype=float),
        train["outcome_value"].to_numpy(dtype=float),
        train["target"].astype(str).to_numpy(),
        source_weights(train, weighting),
        label,
    )
    if fit.get("estimable", False):
        levels = list(fit["levels"])
        intercepts = np.asarray(fit["theta"][: len(levels)], dtype=float)
        fit["common_slope"] = float(fit["theta"][len(levels)])
        fit["transport_intercept"] = transport_intercept(
            levels,
            intercepts,
            train["target"].astype(str).to_numpy(),
            weighting,
        )
    else:
        fit["common_slope"] = np.nan
        fit["transport_intercept"] = np.nan
    return fit


def choose_f_lambda(train: pd.DataFrame) -> tuple[float, pd.DataFrame]:
    y = train["outcome_value"].to_numpy(dtype=int)
    eta = train["eta_U"].to_numpy(dtype=float)
    groups = train["sequence_cluster_value"].astype(str).to_numpy()
    splitter = StratifiedGroupKFold(
        n_splits=3,
        shuffle=True,
        random_state=20261113,
    )
    splits = list(splitter.split(np.zeros(len(train)), y, groups))
    rows = []
    for penalty_lambda in LAMBDAS:
        preds = []
        outcomes = []
        status = "evaluable"
        failure_type = ""
        for inner_fold, (itr, iva) in enumerate(splits, start=1):
            fit = rq3.fit_offset_slope(
                eta[itr],
                y[itr],
                penalty_lambda=penalty_lambda,
                label=f"inner|lambda_{penalty_lambda:g}|fold_{inner_fold}",
            )
            if not fit.get("estimable", False):
                status = "failed"
                failure_type = fit.get("failure_type", "inner_fit_failure")
                break
            preds.append(expit(fit["intercept"] + fit["gamma"] * eta[iva]))
            outcomes.append(y[iva])
        rows.append(
            {
                "lambda": penalty_lambda,
                "status": status,
                "failure_type": failure_type,
                "pooled_inner_log_loss": (
                    log_loss(np.concatenate(outcomes), np.concatenate(preds))
                    if status == "evaluable"
                    else np.nan
                ),
            }
        )
    candidates = pd.DataFrame(rows)
    valid = candidates[candidates["status"].eq("evaluable")]
    if valid.empty:
        return np.nan, candidates
    minimum = float(valid["pooled_inner_log_loss"].min())
    tied = valid[
        (valid["pooled_inner_log_loss"] - minimum).abs().le(1e-12)
    ]["lambda"].astype(float)
    selected = float(tied.max())
    candidates["selected"] = candidates["lambda"].eq(selected) & candidates[
        "status"
    ].eq("evaluable")
    return selected, candidates


def fit_adaptation(
    target_rows: pd.DataFrame,
    cv_target: pd.DataFrame,
    *,
    score_model: str,
    weighting: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    predictions = []
    parameters = []
    lambda_rows = []
    fold_ids = sorted(cv_target["fold_id"].unique())
    if len(fold_ids) != 5:
        raise RuntimeError(f"Expected five CV folds; observed {len(fold_ids)}")
    for fold_id in fold_ids:
        assignments = cv_target[cv_target["fold_id"].eq(fold_id)][
            ["pair_id", "split"]
        ].drop_duplicates()
        frame = target_rows.merge(
            assignments,
            on="pair_id",
            how="inner",
            validate="one_to_one",
        )
        train = frame[frame["split"].eq("train")].copy().reset_index(drop=True)
        test = frame[frame["split"].eq("test")].copy().reset_index(drop=True)
        if set(train["sequence_cluster_value"]) & set(
            test["sequence_cluster_value"]
        ):
            raise RuntimeError(f"Cluster leakage in {score_model} {weighting} {fold_id}")
        y_train = train["outcome_value"].to_numpy(dtype=float)
        y_test = test["outcome_value"].to_numpy(dtype=float)
        eta_train = train["eta_U"].to_numpy(dtype=float)
        eta_test = test["eta_U"].to_numpy(dtype=float)

        i_fit = rq3.solve_offset_intercept(
            eta_train,
            y_train,
            f"{score_model}|{weighting}|{fold_id}|I",
        )
        if i_fit["status"] != "evaluable":
            raise RuntimeError(f"I fit failed: {i_fit}")
        p_i = expit(i_fit["intercept"] + eta_test)

        selected_lambda, candidates = choose_f_lambda(train)
        candidates.insert(0, "fold_id", fold_id)
        candidates.insert(0, "source_weighting", weighting)
        candidates.insert(0, "model", score_model)
        lambda_rows.append(candidates)
        if not np.isfinite(selected_lambda):
            raise RuntimeError(f"All F lambdas failed for {score_model} {weighting} {fold_id}")
        f_fit = rq3.fit_offset_slope(
            eta_train,
            y_train,
            penalty_lambda=selected_lambda,
            label=f"{score_model}|{weighting}|{fold_id}|F",
        )
        if not f_fit.get("estimable", False):
            raise RuntimeError(f"Final F fit failed: {f_fit}")
        p_f = expit(f_fit["intercept"] + f_fit["gamma"] * eta_test)

        predictions.append(
            pd.DataFrame(
                {
                    "pair_id": test["pair_id"].astype(str),
                    "target": test["target"].astype(str),
                    "sequence_cluster_value": test["sequence_cluster_value"].astype(str),
                    "outer_fold_id": fold_id,
                    "model": score_model,
                    "source_weighting": weighting,
                    "outcome_value": y_test.astype(int),
                    "raw_score": test["raw_score"].to_numpy(dtype=float),
                    "eta_U": eta_test,
                    "p_U": test["p_U"].to_numpy(dtype=float),
                    "p_I": p_i,
                    "p_F": p_f,
                    "adaptation_train_reference": float(y_train.mean()),
                    "bss_transport_reference": test[
                        "bss_transport_reference"
                    ].to_numpy(dtype=float),
                    "bss_target_oracle_reference": test[
                        "bss_target_oracle_reference"
                    ].to_numpy(dtype=float),
                    "selected_lambda": selected_lambda,
                }
            )
        )
        parameters.append(
            {
                "model": score_model,
                "source_weighting": weighting,
                "fold_id": fold_id,
                "I_intercept": i_fit["intercept"],
                "F_intercept": f_fit["intercept"],
                "F_delta": f_fit["delta"],
                "F_gamma": f_fit["gamma"],
                "selected_lambda": selected_lambda,
                "train_n": len(train),
                "test_n": len(test),
            }
        )
    result = pd.concat(predictions, ignore_index=True)
    if result["pair_id"].duplicated().any():
        raise RuntimeError("Each adaptation pair must be tested exactly once")
    return result, pd.DataFrame(parameters), pd.concat(lambda_rows, ignore_index=True)


def parse_ipsae_output(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    start = next(index for index, line in enumerate(lines) if line.lstrip().startswith("Chn1"))
    from io import StringIO

    return pd.read_csv(StringIO("\n".join(lines[start:])), sep=r"\s+")


def run_ipsae_sample_qc(data: pd.DataFrame) -> pd.DataFrame:
    sample_parts = []
    for _, group in data.sort_values(["target", "boltz2_min_ipsae"]).groupby(
        "target", sort=True
    ):
        sample_parts.append(
            group.iloc[sorted({0, len(group) // 2, len(group) - 1})]
        )
    sample = pd.concat(sample_parts, ignore_index=True)
    sample_dir = IPSAE_OUT / "sample_qc"
    sample_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in sample.itertuples(index=False):
        safe = str(row.pair_id).replace("|", "__")
        source = ROOT / "data" / "pinc_batch" / "inputs" / safe
        case_dir = sample_dir / safe
        case_dir.mkdir(parents=True, exist_ok=True)
        pae_source = source / "pae.json"
        cif_source = source / "structure.cif"
        if not pae_source.is_file() or not cif_source.is_file():
            rows.append(
                {
                    "pair_id": row.pair_id,
                    "target": row.target,
                    "status": "failed",
                    "failure_type": "missing_local_input",
                }
            )
            continue
        payload = json.loads(pae_source.read_text(encoding="utf-8"))
        if isinstance(payload, list):
            payload = payload[0]
        wrapper = case_dir / "pae.af3.json"
        wrapper.write_text(
            json.dumps({"pae": payload["predicted_aligned_error"]}),
            encoding="utf-8",
        )
        cif_copy = case_dir / "structure.cif"
        shutil.copy2(cif_source, cif_copy)
        completed = subprocess.run(
            [
                str(PYTHON),
                str(IPSAE_SCRIPT),
                str(wrapper),
                str(cif_copy),
                "10",
                "10",
            ],
            cwd=case_dir,
            capture_output=True,
            text=True,
            timeout=180,
        )
        (case_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (case_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
        output = case_dir / "structure_10_10.txt"
        if completed.returncode != 0 or not output.is_file():
            rows.append(
                {
                    "pair_id": row.pair_id,
                    "target": row.target,
                    "status": "failed",
                    "failure_type": "ipsae_execution_failed",
                    "return_code": completed.returncode,
                    "stderr": completed.stderr[-1000:],
                }
            )
            continue
        scores = parse_ipsae_output(output)
        scores = scores[scores["Type"].astype(str).eq("asym")].copy()
        d2t = scores[
            scores["Chn1"].astype(str).eq(str(row.design_chain))
            & scores["Chn2"].astype(str).eq(str(row.target_chain))
        ]["ipSAE"]
        t2d = scores[
            scores["Chn1"].astype(str).eq(str(row.target_chain))
            & scores["Chn2"].astype(str).eq(str(row.design_chain))
        ]["ipSAE"]
        if len(d2t) != 1 or len(t2d) != 1:
            rows.append(
                {
                    "pair_id": row.pair_id,
                    "target": row.target,
                    "status": "failed",
                    "failure_type": "directional_rows_not_unique",
                }
            )
            continue
        design_to_target = float(d2t.iloc[0])
        target_to_design = float(t2d.iloc[0])
        official_max = max(design_to_target, target_to_design)
        official_min = min(design_to_target, target_to_design)
        match_max = abs(official_max - float(row.boltz2_ipsae)) <= 1e-6
        match_min = abs(official_min - float(row.boltz2_min_ipsae)) <= 1e-6
        rows.append(
            {
                "pair_id": row.pair_id,
                "target": row.target,
                "design_chain": row.design_chain,
                "target_chain": row.target_chain,
                "official_design_to_target": design_to_target,
                "official_target_to_design": target_to_design,
                "official_max": official_max,
                "official_min": official_min,
                "released_ipsae": row.boltz2_ipsae,
                "released_min_ipsae": row.boltz2_min_ipsae,
                "max_matches_release": match_max,
                "min_matches_release": match_min,
                "status": "passed",
                "failure_type": "",
            }
        )
    result = pd.DataFrame(rows)
    write_csv(IPSAE_OUT / "ipsae_sample_qc.csv", result)
    return result


def run_ipsae_one(row, batch_root: Path) -> dict:
    safe = str(row.pair_id).replace("|", "__")
    source = ROOT / "data" / "pinc_batch" / "inputs" / safe
    pae_source = source / "pae.json"
    cif_source = source / "structure.cif"
    if not pae_source.is_file() or not cif_source.is_file():
        return {
            "pair_id": row.pair_id,
            "target": row.target,
            "status": "failed",
            "failure_type": "missing_local_input",
        }
    try:
        with TemporaryDirectory(prefix=f"{safe}__", dir=batch_root) as tmp:
            case_dir = Path(tmp)
            payload = json.loads(pae_source.read_text(encoding="utf-8"))
            if isinstance(payload, list):
                payload = payload[0]
            wrapper = case_dir / "pae.af3.json"
            wrapper.write_text(
                json.dumps({"pae": payload["predicted_aligned_error"]}),
                encoding="utf-8",
            )
            cif_copy = case_dir / "structure.cif"
            shutil.copy2(cif_source, cif_copy)
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(IPSAE_SCRIPT),
                    str(wrapper),
                    str(cif_copy),
                    "10",
                    "10",
                ],
                cwd=case_dir,
                capture_output=True,
                text=True,
                timeout=180,
            )
            output = case_dir / "structure_10_10.txt"
            if completed.returncode != 0 or not output.is_file():
                return {
                    "pair_id": row.pair_id,
                    "target": row.target,
                    "status": "failed",
                    "failure_type": "ipsae_execution_failed",
                    "return_code": completed.returncode,
                    "stderr": completed.stderr[-1000:],
                }
            scores = parse_ipsae_output(output)
            scores = scores[scores["Type"].astype(str).eq("asym")].copy()
            d2t = scores[
                scores["Chn1"].astype(str).eq(str(row.design_chain))
                & scores["Chn2"].astype(str).eq(str(row.target_chain))
            ]["ipSAE"]
            t2d = scores[
                scores["Chn1"].astype(str).eq(str(row.target_chain))
                & scores["Chn2"].astype(str).eq(str(row.design_chain))
            ]["ipSAE"]
            if len(d2t) != 1 or len(t2d) != 1:
                return {
                    "pair_id": row.pair_id,
                    "target": row.target,
                    "status": "failed",
                    "failure_type": "directional_rows_not_unique",
                }
            design_to_target = float(d2t.iloc[0])
            target_to_design = float(t2d.iloc[0])
            return {
                "pair_id": row.pair_id,
                "target": row.target,
                "status": "passed",
                "failure_type": "",
                "official_design_to_target": design_to_target,
                "official_target_to_design": target_to_design,
                "official_ipsae_max": max(design_to_target, target_to_design),
                "official_ipsae_min": min(design_to_target, target_to_design),
                "released_ipsae": row.boltz2_ipsae,
                "released_min_ipsae": row.boltz2_min_ipsae,
                "max_matches_release": abs(
                    max(design_to_target, target_to_design)
                    - float(row.boltz2_ipsae)
                )
                <= 1e-6,
                "min_matches_release": abs(
                    min(design_to_target, target_to_design)
                    - float(row.boltz2_min_ipsae)
                )
                <= 1e-6,
            }
    except Exception as error:
        return {
            "pair_id": row.pair_id,
            "target": row.target,
            "status": "failed",
            "failure_type": type(error).__name__,
            "error": str(error)[:1000],
        }


def run_ipsae_batch(data: pd.DataFrame) -> pd.DataFrame:
    batch_root = IPSAE_OUT / "batch_work"
    batch_root.mkdir(parents=True, exist_ok=True)
    cached_path = IPSAE_OUT / "ipsae_official_pair_metrics.csv"
    if cached_path.is_file():
        cached = pd.read_csv(cached_path, low_memory=False)
        if (
            len(cached) == len(data)
            and cached["status"].eq("passed").all()
            and set(cached["pair_id"].astype(str))
            == set(data["pair_id"].astype(str))
        ):
            return cached.sort_values("pair_id").reset_index(drop=True)
    rows = []
    workers = min(8, max(1, os.cpu_count() or 1))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(run_ipsae_one, row, batch_root)
            for row in data.itertuples(index=False)
        ]
        for future in as_completed(futures):
            rows.append(future.result())
    result = pd.DataFrame(rows).sort_values("pair_id").reset_index(drop=True)
    write_csv(cached_path, result)
    if len(result) != len(data) or not result["status"].eq("passed").all():
        failures = result[result["status"].ne("passed")]
        raise RuntimeError(
            f"Official ipSAE batch failures: {len(failures)}"
        )
    return result


def build_c01_extension() -> None:
    IPSAE_OUT.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(PAIR_MASTER, low_memory=False)
    pinc = pd.read_csv(PINC, low_memory=False)
    audit = pd.read_csv(PAIR_AUDIT, low_memory=False)
    clusters = pd.read_csv(CLUSTERS, low_memory=False)
    loto = pd.read_csv(LOTO, low_memory=False)
    cv5 = pd.read_csv(CV5, low_memory=False)

    audit["pair_id"] = (
        audit["design_id"].astype(str) + "|" + audit["target"].astype(str)
    )
    if audit["pair_id"].duplicated().any():
        raise RuntimeError("Pair-audit pair_id is not unique")
    data = (
        master[
            [
                "pair_id",
                "design_id",
                "target",
                "design_class",
                "boltz2_iptm",
                "primary_outcome",
                "strict_outcome",
            ]
        ]
        .merge(
            pinc[
                [
                    "pair_id",
                    "pinc",
                    "pinc_status",
                    "cif_url",
                    "pae_url",
                    "pae_n_rows",
                    "pae_n_cols",
                    "pae_all_finite",
                    "chain_order",
                    "design_chain",
                    "target_chain",
                    "design_chain_residue_n",
                    "target_chain_residue_n",
                    "total_residue_tokens",
                    "cif_sha256",
                    "pae_sha256",
                ]
            ],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            audit[["pair_id", "boltz2_ipsae", "boltz2_min_ipsae"]],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            clusters[
                [
                    "pair_id",
                    "sequence_cluster90",
                    "sequence_cluster80",
                ]
            ],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
    )
    if len(data) != 2084:
        raise RuntimeError("C01 master row count is not 2084")
    complete = data[
        data[["boltz2_iptm", "boltz2_min_ipsae", "pinc"]].notna().all(axis=1)
    ].copy()
    if len(complete) != 2074 or int(complete["primary_outcome"].sum()) != 740:
        raise RuntimeError("C01 complete-case count or positive count changed")
    complete["outcome_value"] = complete["primary_outcome"].astype(int)

    sample_qc = run_ipsae_sample_qc(complete)
    if len(sample_qc) == 0 or not sample_qc["status"].eq("passed").all():
        raise RuntimeError("Official ipSAE sample QC failed")
    official = run_ipsae_batch(complete)
    complete = complete.merge(
        official[
            [
                "pair_id",
                "official_design_to_target",
                "official_target_to_design",
                "official_ipsae_max",
                "official_ipsae_min",
                "max_matches_release",
                "min_matches_release",
            ]
        ],
        on="pair_id",
        how="left",
        validate="one_to_one",
    )
    if complete[["official_ipsae_min", "official_ipsae_max"]].isna().any().any():
        raise RuntimeError("Official ipSAE batch merge lost rows")
    complete["x_ipTM"] = logit_clip(complete["boltz2_iptm"])
    complete["x_ipSAE"] = logit_clip(complete["official_ipsae_min"])
    complete["x_Pinc"] = logit_clip(complete["pinc"])
    write_csv(
        IPSAE_OUT / "ipSAE_pair_extension.csv",
        complete[
            [
                "pair_id",
                "design_id",
                "target",
                "design_class",
                "sequence_cluster90",
                "primary_outcome",
                "strict_outcome",
                "pinc",
                "boltz2_iptm",
                "official_ipsae_min",
                "official_ipsae_max",
                "boltz2_ipsae",
                "boltz2_min_ipsae",
                "max_matches_release",
                "min_matches_release",
                "chain_order",
                "design_chain",
                "target_chain",
                "pae_n_rows",
                "pae_n_cols",
                "total_residue_tokens",
                "cif_sha256",
                "pae_sha256",
            ]
        ],
    )

    mapping_checks = {
        "complete_cases": int(len(complete)),
        "primary_positives": int(complete["outcome_value"].sum()),
        "targets": complete.groupby("target").size().to_dict(),
        "pae_square_and_matches_tokens": bool(
            complete["pae_n_rows"].eq(complete["pae_n_cols"]).all()
            and complete["pae_n_rows"].eq(complete["total_residue_tokens"]).all()
        ),
        "pae_all_finite": bool(complete["pae_all_finite"].eq(True).all()),
        "chain_fields_complete": bool(
            complete[["chain_order", "design_chain", "target_chain"]]
            .notna()
            .all()
            .all()
        ),
        "residue_counts_sum_to_tokens": bool(
            complete["design_chain_residue_n"]
            .add(complete["target_chain_residue_n"])
            .eq(complete["total_residue_tokens"])
            .all()
        ),
        "sample_qc_n": int(len(sample_qc)),
        "sample_qc_passed": int(sample_qc["status"].eq("passed").sum()),
    }
    if not all(
        value
        for key, value in mapping_checks.items()
        if key
        in {
            "pae_square_and_matches_tokens",
            "pae_all_finite",
            "chain_fields_complete",
            "residue_counts_sum_to_tokens",
        }
    ):
        raise RuntimeError("C01 PAE/CIF/chain mapping audit failed")

    affinity_metrics = set()
    for value in pd.read_csv(PROTEINBASE, usecols=["evaluations"])["evaluations"].dropna():
        try:
            evaluations = json.loads(value)
        except Exception:
            continue
        for evaluation in evaluations:
            if "affin" in str(evaluation.get("metric", "")).lower():
                affinity_metrics.add(str(evaluation.get("metric")))

    fold_ids = (
        loto[
            loto["analysis_set"].eq("primary_six")
            & loto["cluster_rule"].eq("sequence_cluster90")
        ]["fold_id"]
        .drop_duplicates()
        .tolist()
    )
    if len(fold_ids) != 6:
        raise RuntimeError("Expected six primary-six LOTO folds")

    u_frames = []
    source_rows = []
    complete_ids = set(complete["pair_id"].astype(str))
    for fold_id in fold_ids:
        fold = loto[loto["fold_id"].eq(fold_id)].copy()
        fold = fold[fold["pair_id"].astype(str).isin(complete_ids)]
        train = fold[fold["split_status"].eq("train")].merge(
            complete[
                [
                    "pair_id",
                    "boltz2_iptm",
                    "boltz2_min_ipsae",
                    "official_ipsae_min",
                    "x_ipTM",
                    "x_ipSAE",
                ]
            ],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        test = fold[fold["split_status"].eq("test")].merge(
            complete[
                [
                    "pair_id",
                    "sequence_cluster90",
                    "boltz2_iptm",
                    "boltz2_min_ipsae",
                    "official_ipsae_min",
                    "x_ipTM",
                    "x_ipSAE",
                ]
            ],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        if train[["x_ipTM", "x_ipSAE"]].isna().any().any():
            raise RuntimeError("Missing source score after complete-case filter")
        if test[["x_ipTM", "x_ipSAE"]].isna().any().any():
            raise RuntimeError("Missing held-out score after complete-case filter")
        held_target = str(test["target"].unique()[0])
        for weighting in WEIGHTINGS:
            ref, _ = source_reference(train, weighting)
            for model, column in (("ipTM_only", "x_ipTM"), ("ipSAE_only", "x_ipSAE")):
                fit = fit_source(
                    train,
                    column,
                    weighting,
                    f"{fold_id}|{weighting}|{model}",
                )
                source_rows.append(
                    {
                        "fold_id": fold_id,
                        "held_out_target": held_target,
                        "source_weighting": weighting,
                        "model": model,
                        "fit_status": fit.get("status"),
                        "failure_type": fit.get("failure_type", ""),
                        "transport_intercept": fit.get("transport_intercept"),
                        "common_slope": fit.get("common_slope"),
                        "source_weights_sum": float(
                            source_weights(train, weighting).sum()
                        ),
                        "source_target_count": int(train["target"].nunique()),
                        "source_n": int(len(train)),
                    }
                )
                if not fit.get("estimable", False):
                    raise RuntimeError(f"C01 source fit failed: {fold_id} {weighting} {model}")
                eta = (
                    fit["transport_intercept"]
                    + fit["common_slope"] * test[column].to_numpy(dtype=float)
                )
                u_frames.append(
                    pd.DataFrame(
                        {
                            "pair_id": test["pair_id"].astype(str),
                            "target": test["target"].astype(str),
                            "sequence_cluster_value": test[
                                "sequence_cluster90"
                            ].astype(str),
                            "outer_fold_id": fold_id,
                            "model": model,
                            "source_weighting": weighting,
                            "outcome_value": test["outcome_value"].astype(int),
                            "raw_score": test[
                                "boltz2_iptm"
                                if model == "ipTM_only"
                                else "official_ipsae_min"
                            ].to_numpy(dtype=float),
                            "eta_U": eta,
                            "p_U": expit(eta),
                            "bss_transport_reference": ref,
                            "bss_target_oracle_reference": float(
                                test["outcome_value"].mean()
                            ),
                        }
                    )
                )
    u_predictions = pd.concat(u_frames, ignore_index=True)
    if len(u_predictions) != 2074 * 2 * 2:
        raise RuntimeError("U prediction row count mismatch")

    adaptation_frames = []
    parameter_frames = []
    lambda_frames = []
    for (model, weighting, target), target_rows in u_predictions.groupby(
        ["model", "source_weighting", "target"], sort=True
    ):
        cv_target = cv5[
            cv5["analysis_set"].eq("primary_six")
            & cv5["outcome_name"].eq("primary_outcome")
            & cv5["cluster_rule"].eq("sequence_cluster90")
            & cv5["target"].astype(str).eq(str(target))
        ].copy()
        adapted, parameters, lambdas = fit_adaptation(
            target_rows,
            cv_target,
            score_model=model,
            weighting=weighting,
        )
        adaptation_frames.append(adapted)
        parameter_frames.append(parameters)
        lambda_frames.append(lambdas)
    adapted = pd.concat(adaptation_frames, ignore_index=True)
    if len(adapted) != 2074 * 2 * 2:
        raise RuntimeError("Adaptation prediction row count mismatch")
    write_csv(IPSAE_OUT / "ipSAE_UIF_pair_predictions.csv", adapted)
    write_csv(
        IPSAE_OUT / "ipSAE_source_parameters.csv",
        pd.DataFrame(source_rows),
    )
    write_csv(
        IPSAE_OUT / "ipSAE_adaptation_parameters.csv",
        pd.concat(parameter_frames, ignore_index=True),
    )
    write_csv(
        IPSAE_OUT / "ipSAE_lambda_selection.csv",
        pd.concat(lambda_frames, ignore_index=True),
    )

    target_rows = []
    raw_rows = []
    for (model, weighting, target), group in adapted.groupby(
        ["model", "source_weighting", "target"], sort=True
    ):
        y = group["outcome_value"].to_numpy(dtype=float)
        raw_rows.append(
            {
                "model": model,
                "source_weighting": weighting,
                "target": target,
                "layer": "raw",
                **metric_set(
                    y,
                    group["raw_score"].to_numpy(dtype=float),
                    target_reference=group[
                        "bss_target_oracle_reference"
                    ].to_numpy(dtype=float),
                ),
            }
        )
        for layer in LAYERS:
            target_rows.append(
                {
                    "model": model,
                    "source_weighting": weighting,
                    "target": target,
                    "layer": layer,
                    **metric_set(
                        y,
                        group[f"p_{layer}"].to_numpy(dtype=float),
                        transport_reference=group[
                            "bss_transport_reference"
                        ].to_numpy(dtype=float),
                        target_reference=group[
                            "bss_target_oracle_reference"
                        ].to_numpy(dtype=float),
                    ),
                    "bss_adaptation_train": (
                        bss(
                            y,
                            group[f"p_{layer}"].to_numpy(dtype=float),
                            group["adaptation_train_reference"].to_numpy(dtype=float),
                        )
                        if layer in {"I", "F"}
                        else np.nan
                    ),
                }
            )
    target_metrics = pd.DataFrame(raw_rows + target_rows)
    write_csv(IPSAE_OUT / "ipSAE_target_metrics.csv", target_metrics)

    macro_rows = []
    for (model, weighting, layer), group in target_metrics.groupby(
        ["model", "source_weighting", "layer"], sort=True
    ):
        row = {
            "model": model,
            "source_weighting": weighting,
            "layer": layer,
            "summary_type": "target_equal_macro",
            "target_count": len(group),
            "n_pairs": int(group["n"].sum()),
            "positives": int(group["positives"].sum()),
        }
        for metric in (
            "prevalence",
            "roc_auc",
            "ap",
            "brier",
            "log_loss",
            "bss_transport",
            "bss_target_oracle",
            "bss_adaptation_train",
        ):
            if metric in group.columns:
                row[metric] = float(group[metric].mean())
        macro_rows.append(row)
    macro_metrics = pd.DataFrame(macro_rows)
    write_csv(IPSAE_OUT / "ipSAE_macro_metrics.csv", macro_metrics)

    comparison_rows = []
    for (weighting, layer, target), group in target_metrics.groupby(
        ["source_weighting", "layer", "target"], sort=True
    ):
        by_model = group.set_index("model")
        if set(SCORE_MODELS) - set(by_model.index):
            raise RuntimeError("Missing ipTM/ipSAE target metric for comparison")
        comparison_rows.append(
            {
                "source_weighting": weighting,
                "layer": layer,
                "target": target,
                "auc_ipSAE_minus_ipTM": by_model.loc["ipSAE_only", "roc_auc"]
                - by_model.loc["ipTM_only", "roc_auc"],
                "ap_ipSAE_minus_ipTM": by_model.loc["ipSAE_only", "ap"]
                - by_model.loc["ipTM_only", "ap"],
                "brier_ipSAE_minus_ipTM": by_model.loc["ipSAE_only", "brier"]
                - by_model.loc["ipTM_only", "brier"],
                "log_loss_ipSAE_minus_ipTM": by_model.loc[
                    "ipSAE_only", "log_loss"
                ]
                - by_model.loc["ipTM_only", "log_loss"],
                "bss_transport_ipSAE_minus_ipTM": by_model.loc[
                    "ipSAE_only", "bss_transport"
                ]
                - by_model.loc["ipTM_only", "bss_transport"],
            }
        )
    target_comparisons = pd.DataFrame(comparison_rows)
    macro_comparisons = (
        target_comparisons.groupby(["source_weighting", "layer"], as_index=False)
        .agg(
            target_count=("target", "nunique"),
            auc_ipSAE_minus_ipTM=("auc_ipSAE_minus_ipTM", "mean"),
            ap_ipSAE_minus_ipTM=("ap_ipSAE_minus_ipTM", "mean"),
            brier_ipSAE_minus_ipTM=("brier_ipSAE_minus_ipTM", "mean"),
            log_loss_ipSAE_minus_ipTM=("log_loss_ipSAE_minus_ipTM", "mean"),
            bss_transport_ipSAE_minus_ipTM=(
                "bss_transport_ipSAE_minus_ipTM",
                "mean",
            ),
        )
    )
    macro_comparisons["summary_type"] = "target_equal_macro"
    paired = pd.concat(
        [
            target_comparisons.assign(summary_type="per_target"),
            macro_comparisons,
        ],
        ignore_index=True,
        sort=False,
    )
    write_csv(IPSAE_OUT / "ipSAE_paired_comparisons.csv", paired)

    report = [
        "# C01 ipSAE feasibility and extension",
        "",
        "**Evidence level:** post hoc sensitivity extension added after the original ipTM analyses had been examined.",
        "",
        "## Feasibility",
        "",
        f"- Complete-case cohort: {len(complete)} pairs, {int(complete['outcome_value'].sum())} primary positives, six targets.",
        "- Missing set: the same 10 PD-L1 pairs lack both Pinc and Boltz-2 ipSAE inputs; no additional loss was introduced by the ipSAE comparison.",
        "- Released fields audited: `boltz2_ipsae` and `boltz2_min_ipsae`.",
        f"- Official implementation: DunbrackLab/IPSAE commit `{IPSAE_COMMIT}`, `ipsae.py`, PAE cutoff 10 Å, distance cutoff 10 Å.",
        f"- Official script SHA-256: `{sha256(IPSAE_SCRIPT)}`.",
        f"- Stratified sample QC: {int(sample_qc['status'].eq('passed').sum())}/{len(sample_qc)} passed.",
        f"- Released minimum field exactly matched the official directional minimum in {int(sample_qc['min_matches_release'].sum())}/{len(sample_qc)} sampled pairs; the comparison is recorded but the released field was not used as the algorithmic gold standard.",
        f"- Across the full batch, the released minimum field matched in {int(complete['min_matches_release'].sum())}/{len(complete)} pairs and the released non-minimum field matched in {int(complete['max_matches_release'].sum())}/{len(complete)} pairs.",
        "- The extension therefore uses freshly recomputed `ipSAE_min` from the stored Boltz-2 PAE/mmCIF pair with the pinned official implementation.",
        f"- PAE dimensions, finite values, chain labels, and residue-token sums were checked for all {len(complete)} pairs.",
        f"- Boltz-2 affinity metrics found in the frozen Proteinbase file: {sorted(affinity_metrics) if affinity_metrics else 'none'}. No affinity output was available for this extension.",
        "",
        "## Extension design",
        "",
        "- The ipTM comparator and ipSAE model were both refitted on the same 2,074 rows, frozen LOTO and target-CV5 assignments, 90% sequence clusters, and source weighting.",
        "- The source model used one intercept per source target and one common slope.",
        "- I retained the frozen U prediction as an offset and estimated an intercept only.",
        "- F estimated intercept plus a ridge-shrunk slope change; lambda was selected only within adaptation-training folds by pooled inner log loss.",
        "- These results are an exploratory extension on complete cases and do not replace the frozen 2,084-pair ipTM analysis.",
        "",
    ]
    (IPSAE_OUT / "feasibility_report.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )

    manifest = {
        "status": "complete",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "analysis": "C01",
        "analysis_role": "post_hoc_sensitivity_extension",
        "complete_cases": len(complete),
        "positives": int(complete["outcome_value"].sum()),
        "ipSAE_definition": "minimum directional ipSAE at PAE cutoff 10 A and distance cutoff 10 A",
        "released_metric_fields": ["boltz2_ipsae", "boltz2_min_ipsae"],
        "official_implementation": {
            "repository": "https://github.com/DunbrackLab/IPSAE",
            "commit": IPSAE_COMMIT,
            "script": str(IPSAE_SCRIPT.relative_to(ROOT)),
            "script_sha256": sha256(IPSAE_SCRIPT),
            "parameters": {"pae_cutoff": 10, "distance_cutoff": 10},
        },
        "affinity_metrics_found": sorted(affinity_metrics),
        "sample_qc": {
            "n": int(len(sample_qc)),
            "passed": int(sample_qc["status"].eq("passed").sum()),
            "failed": int(sample_qc["status"].eq("failed").sum()),
        },
        "mapping_checks": mapping_checks,
        "released_field_agreement": {
            "minimum_field_exact_matches": int(complete["min_matches_release"].sum()),
            "maximum_field_exact_matches": int(complete["max_matches_release"].sum()),
            "complete_cases": len(complete),
        },
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                PAIR_MASTER,
                PAIR_AUDIT,
                PINC,
                CLUSTERS,
                LOTO,
                CV5,
                PROTEINBASE,
                IPSAE_SCRIPT,
            )
        },
        "outputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                IPSAE_OUT / "ipSAE_pair_extension.csv",
                IPSAE_OUT / "ipsae_official_pair_metrics.csv",
                IPSAE_OUT / "ipSAE_UIF_pair_predictions.csv",
                IPSAE_OUT / "ipSAE_source_parameters.csv",
                IPSAE_OUT / "ipSAE_adaptation_parameters.csv",
                IPSAE_OUT / "ipSAE_lambda_selection.csv",
                IPSAE_OUT / "ipSAE_target_metrics.csv",
                IPSAE_OUT / "ipSAE_macro_metrics.csv",
                IPSAE_OUT / "ipSAE_paired_comparisons.csv",
                IPSAE_OUT / "ipsae_sample_qc.csv",
                IPSAE_OUT / "feasibility_report.md",
            )
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "script_sha256": sha256(Path(__file__)),
        "frozen_inputs_modified": False,
    }
    write_json(IPSAE_OUT / "run_manifest.json", manifest)


def build_dataset_for_bootstrap() -> dict[str, pd.DataFrame]:
    main = pd.read_csv(MAIN_FROZEN_PREDICTIONS, low_memory=False)
    external = pd.read_csv(EXTERNAL_PREDICTIONS, low_memory=False)
    datasets = {}
    for dataset, analysis_set, expected_n, expected_pos in (
        ("main_primary_six", "primary_six", 2084, 748),
        ("main_primary_five", "primary_five", 1951, 729),
    ):
        frame = main[
            main["analysis_set"].eq(analysis_set)
            & main["cluster_rule"].eq("sequence_cluster90")
            & main["source_weighting"].eq("target_equal")
        ].copy()
        if len(frame) != expected_n or int(frame["outcome_value"].sum()) != expected_pos:
            raise RuntimeError(f"Frozen {dataset} row/positive count mismatch")
        frame = frame[
            [
                "pair_id",
                "target",
                "sequence_cluster_value",
                "outcome_value",
                "p_U",
                "p_I",
                "p_F",
                "bss_transport_reference",
            ]
        ].rename(columns={"pair_id": "id", "sequence_cluster_value": "cluster"})
        datasets[dataset] = frame
    external = external[
        [
            "uuid",
            "target",
            "sequence_cluster90",
            "outcome_value",
            "p_U",
            "p_I",
            "p_F",
            "bss_transport_reference",
        ]
    ].rename(columns={"uuid": "id", "sequence_cluster90": "cluster"})
    if len(external) != 720 or external["id"].nunique() != 720:
        raise RuntimeError("Frozen external cohort row count mismatch")
    datasets["external_primary_weak_plus"] = external
    return datasets


def bootstrap_metric_bundle(frame: pd.DataFrame) -> dict[str, float]:
    y = frame["outcome_value"].to_numpy(dtype=float)
    ref = frame["bss_transport_reference"].to_numpy(dtype=float)
    values = {}
    for layer in LAYERS:
        p = frame[f"p_{layer}"].to_numpy(dtype=float)
        values[f"{layer}_brier"] = brier(y, p)
        values[f"{layer}_log_loss"] = log_loss(y, p)
        values[f"{layer}_bss_transport"] = bss(y, p, ref)
    values["I_minus_U_brier"] = values["I_brier"] - values["U_brier"]
    values["I_minus_U_log_loss"] = values["I_log_loss"] - values["U_log_loss"]
    values["F_minus_I_brier"] = values["F_brier"] - values["I_brier"]
    values["F_minus_I_log_loss"] = values["F_log_loss"] - values["I_log_loss"]
    values["I_minus_U_bss_transport"] = (
        values["I_bss_transport"] - values["U_bss_transport"]
    )
    values["F_minus_I_bss_transport"] = (
        values["F_bss_transport"] - values["I_bss_transport"]
    )
    values["relative_target_U_to_I_brier_reduction"] = (
        (values["U_brier"] - values["I_brier"]) / values["U_brier"]
    )
    return values


def bootstrap_dataset(
    dataset: str,
    frame: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    targets = sorted(frame["target"].unique())
    by_target = {target: frame[frame["target"].eq(target)].copy() for target in targets}
    cluster_indices = {}
    for target, block in by_target.items():
        groups = {}
        for cluster, group in block.groupby("cluster", sort=True):
            groups[str(cluster)] = group.index.to_numpy()
        cluster_indices[target] = list(groups.values())
    metric_names = list(bootstrap_metric_bundle(next(iter(by_target.values()))).keys())
    macro_metric_names = metric_names + [
        "relative_macro_U_to_I_brier_reduction"
    ]

    point_target = {
        target: bootstrap_metric_bundle(by_target[target]) for target in targets
    }
    point_macro = {
        metric: float(np.mean([point_target[target][metric] for target in targets]))
        for metric in metric_names
    }
    point_macro["relative_macro_U_to_I_brier_reduction"] = (
        (point_macro["U_brier"] - point_macro["I_brier"])
        / point_macro["U_brier"]
    )

    target_replicates = {
        target: {metric: [] for metric in metric_names} for target in targets
    }
    macro_replicates = {metric: [] for metric in macro_metric_names}
    valid = 0
    attempts = 0
    failed = 0
    index_lookup = {target: frame.index for target in targets}
    while valid < BOOT_REPLICATES and attempts < BOOT_REPLICATES * 5:
        attempts += 1
        replicate_target = {}
        try:
            for target in targets:
                blocks = cluster_indices[target]
                sampled = rng.integers(0, len(blocks), size=len(blocks))
                indices = np.concatenate([blocks[position] for position in sampled])
                values = bootstrap_metric_bundle(frame.loc[indices])
                if not all(np.isfinite(value) for value in values.values()):
                    raise ValueError("nonfinite replicate")
                replicate_target[target] = values
            macro = {
                metric: float(
                    np.mean([replicate_target[target][metric] for target in targets])
                )
                for metric in metric_names
            }
            macro["relative_macro_U_to_I_brier_reduction"] = (
                (macro["U_brier"] - macro["I_brier"]) / macro["U_brier"]
            )
            if not all(np.isfinite(value) for value in macro.values()):
                raise ValueError("nonfinite macro replicate")
        except Exception:
            failed += 1
            continue
        valid += 1
        for target in targets:
            for metric in metric_names:
                target_replicates[target][metric].append(
                    replicate_target[target][metric]
                )
        for metric in macro_metric_names:
            macro_replicates[metric].append(macro[metric])

    if valid < BOOT_REPLICATES:
        raise RuntimeError(
            f"Only {valid} valid bootstrap replicates for {dataset}; required {BOOT_REPLICATES}"
        )

    target_rows = []
    for target in targets:
        for metric in metric_names:
            values = np.asarray(target_replicates[target][metric], dtype=float)
            target_rows.append(
                {
                    "dataset": dataset,
                    "target": target,
                    "metric": metric,
                    "point_estimate": point_target[target][metric],
                    "bootstrap_95_low": float(np.quantile(values, 0.025)),
                    "bootstrap_95_high": float(np.quantile(values, 0.975)),
                    "n_valid_boot": valid,
                    "n_failed_boot": failed,
                }
            )
    macro_rows = []
    for metric in macro_metric_names:
        values = np.asarray(macro_replicates[metric], dtype=float)
        macro_rows.append(
            {
                "dataset": dataset,
                "target": "__macro_target_equal__",
                "metric": metric,
                "point_estimate": point_macro[metric],
                "bootstrap_95_low": float(np.quantile(values, 0.025)),
                "bootstrap_95_high": float(np.quantile(values, 0.975)),
                "n_valid_boot": valid,
                "n_failed_boot": failed,
            }
        )
    qc = {
        "dataset": dataset,
        "n_pairs": len(frame),
        "n_targets": len(targets),
        "targets": targets,
        "n_clusters_by_target": {
            target: len(cluster_indices[target]) for target in targets
        },
        "requested_replicates": BOOT_REPLICATES,
        "valid_replicates": valid,
        "failed_attempts": failed,
        "all_core_estimates_finite": True,
    }
    return pd.DataFrame(target_rows), pd.DataFrame(macro_rows), qc


def build_c02_bootstrap() -> None:
    BOOT_OUT.mkdir(parents=True, exist_ok=True)
    datasets = build_dataset_for_bootstrap()
    rng = np.random.default_rng(BOOT_SEED)
    target_frames = []
    macro_frames = []
    qc_rows = []
    for dataset, frame in datasets.items():
        target, macro, qc = bootstrap_dataset(dataset, frame, rng)
        target_frames.append(target)
        macro_frames.append(macro)
        qc_rows.append(qc)
    target_intervals = pd.concat(target_frames, ignore_index=True)
    macro_intervals = pd.concat(macro_frames, ignore_index=True)
    write_csv(BOOT_OUT / "bootstrap_target_intervals.csv", target_intervals)
    write_csv(BOOT_OUT / "bootstrap_macro_intervals.csv", macro_intervals)

    method = [
        "# C02 paired cluster bootstrap",
        "",
        "- Inputs were frozen pair-level U/I/F predictions; no model was refitted.",
        "- Within each held-out target, 90%-identity sequence clusters were resampled with replacement.",
        "- U, I, and F were evaluated on the same sampled clusters and rows in every replicate.",
        "- Each target metric was recalculated before taking the target-equal macro average.",
        f"- Requested valid replicates: {BOOT_REPLICATES}; RNG seed: {BOOT_SEED}; percentile interval: 2.5% to 97.5%.",
        "- Intervals describe conditional precision in the observed targets and do not represent uncertainty over future targets.",
        "",
        "## Replicate QC",
        "",
        "| Dataset | Pairs | Targets | Valid | Failed attempts |",
        "|---|---:|---:|---:|---:|",
    ]
    for qc in qc_rows:
        method.append(
            f"| {qc['dataset']} | {qc['n_pairs']} | {qc['n_targets']} | "
            f"{qc['valid_replicates']} | {qc['failed_attempts']} |"
        )
    method += [
        "",
        "The interval for `relative_macro_U_to_I_brier_reduction` is calculated by "
        "recomputing `(Brier_U-Brier_I)/Brier_U` in every replicate.",
    ]
    (BOOT_OUT / "bootstrap_method_and_qc.md").write_text(
        "\n".join(method) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "status": "complete",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "analysis": "C02",
        "replicates": BOOT_REPLICATES,
        "seed": BOOT_SEED,
        "interval": "percentile 2.5-97.5",
        "cluster_rule": "within-target sequence_cluster90",
        "same_rows_for_U_I_F": True,
        "model_refit": False,
        "caveat": "Conditional precision within observed targets; not a future-target prediction interval.",
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (MAIN_FROZEN_PREDICTIONS, EXTERNAL_PREDICTIONS)
        },
        "qc": qc_rows,
        "outputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                BOOT_OUT / "bootstrap_target_intervals.csv",
                BOOT_OUT / "bootstrap_macro_intervals.csv",
                BOOT_OUT / "bootstrap_method_and_qc.md",
            )
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "script_sha256": sha256(Path(__file__)),
    }
    write_json(BOOT_OUT / "run_manifest.json", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bootstrap-only",
        action="store_true",
        help="Recompute C02 from frozen U/I/F prediction tables without rerunning C01.",
    )
    args = parser.parse_args()
    if not args.bootstrap_only:
        IPSAE_OUT.mkdir(parents=True, exist_ok=True)
    BOOT_OUT.mkdir(parents=True, exist_ok=True)
    write_json(
        OUT / "RUN_STATUS.json",
        {
            "status": "running",
            "started_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "tasks": ["C02"] if args.bootstrap_only else ["C01", "C02"],
        },
    )
    if not args.bootstrap_only:
        build_c01_extension()
    build_c02_bootstrap()
    write_json(
        OUT / "RUN_STATUS.json",
        {
            "status": "complete",
            "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "tasks": ["C02"] if args.bootstrap_only else ["C01", "C02"],
        },
    )


if __name__ == "__main__":
    main()
