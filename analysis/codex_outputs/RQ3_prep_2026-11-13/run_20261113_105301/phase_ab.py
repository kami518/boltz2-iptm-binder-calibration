"""Deterministic RQ3 preparation phases A and B.

This script only creates fold assignments and source annotations. It does not
fit, tune, evaluate, or select any model.
"""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True

import duckdb
import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import StratifiedGroupKFold

RUN_DIR = Path(__file__).resolve().parent
ROOT = Path(r"<REPRODUCTION_ROOT>")
CONFIG_PATH = (
    ROOT
    / "analysis"
    / "codex_tasks"
    / "RQ3_PREP_CODEX_HANDOFF_2026-11-13"
    / "rq3_prep_config.json"
)
PAIR_PATH = ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv"
CLUSTER_PATH = ROOT / "analysis" / "rq12_outputs" / "sequence_clusters_rq12.csv"
PROTEINBASE_PATH = ROOT / "data" / "proteinbase_latest.csv"


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def collect_environment(config: dict) -> dict:
    disk = shutil.disk_usage(RUN_DIR)
    software = {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scikit-learn": sklearn.__version__,
        "duckdb": duckdb.__version__,
        "pyarrow": package_version("pyarrow"),
        "requests": package_version("requests"),
        "parasail": "1.3.4 (vendored; imported in phase C from analysis/rq12_lib)",
    }
    current = {
        "run_directory": str(RUN_DIR),
        "captured_at": now_iso(),
        "software": software,
        "input_paths_and_hashes": {},
        "disk": {
            "root": str(RUN_DIR.anchor),
            "total_bytes": disk.total,
            "used_bytes": disk.used,
            "free_bytes": disk.free,
            "free_gib": round(disk.free / (1024**3), 3),
        },
    }
    mismatches = []
    for key, spec in config["inputs"].items():
        path = Path(spec["path"])
        actual_sha = sha256_file(path)
        actual_bytes = path.stat().st_size
        record = {
            "path": str(path),
            "expected_sha256": spec["sha256"],
            "actual_sha256": actual_sha,
            "sha256_matches": actual_sha.lower() == spec["sha256"].lower(),
            "bytes": actual_bytes,
        }
        if "expected_rows" in spec:
            record["expected_rows_in_config"] = spec["expected_rows"]
        if "expected_rows_from_prior_audit" in spec:
            record["expected_rows_from_prior_audit"] = spec[
                "expected_rows_from_prior_audit"
            ]
        current["input_paths_and_hashes"][key] = record
        if not record["sha256_matches"]:
            mismatches.append(key)
    current["hash_gate_passed"] = not mismatches
    current["hash_mismatches"] = mismatches
    atomic_json(RUN_DIR / "preflight_environment.json", current)
    if mismatches:
        raise RuntimeError(f"Frozen input SHA-256 mismatch: {mismatches}")
    return current


def prepare_pair_data(config: dict) -> pd.DataFrame:
    pairs = pd.read_csv(PAIR_PATH, low_memory=False)
    clusters = pd.read_csv(CLUSTER_PATH, low_memory=False)
    if len(pairs) != int(config["inputs"]["pair_master"]["expected_rows"]):
        raise RuntimeError(f"Pair master row count is {len(pairs)}")
    if len(clusters) != int(config["inputs"]["sequence_clusters"]["expected_rows"]):
        raise RuntimeError(f"Cluster row count is {len(clusters)}")
    if pairs["pair_id"].nunique() != len(pairs):
        raise RuntimeError("pair_master_v3.csv pair_id is not unique")
    if clusters["pair_id"].nunique() != len(clusters):
        raise RuntimeError("sequence_clusters_rq12.csv pair_id is not unique")
    cluster_columns = [
        "pair_id",
        "design_id",
        "target",
        "sequence",
        *config["sequence_cluster_rules"],
    ]
    joined = pairs.merge(
        clusters[cluster_columns],
        on="pair_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_cluster_source"),
    )
    if len(joined) != len(pairs):
        raise RuntimeError("One-to-one pair/cluster join lost rows")
    for column in ["design_id", "target", "sequence"]:
        mismatch = joined[column].astype(str) != joined[
            f"{column}_cluster_source"
        ].astype(str)
        if mismatch.any():
            raise RuntimeError(
                f"{mismatch.sum()} rows disagree on {column} across frozen inputs"
            )
    joined = joined.drop(
        columns=[
            "design_id_cluster_source",
            "target_cluster_source",
            "sequence_cluster_source",
        ]
    )
    for rule in config["sequence_cluster_rules"]:
        if joined[rule].isna().any():
            raise RuntimeError(f"Missing sequence cluster values in {rule}")
    return joined


def analysis_subset(
    joined: pd.DataFrame, set_name: str, spec: dict
) -> pd.DataFrame:
    outcome_column = spec["outcome"]
    frame = joined.copy()
    excluded = spec.get("exclude_target")
    if excluded is not None:
        frame = frame[
            frame["target"].astype(str).str.casefold()
            != str(excluded).casefold()
        ].copy()
    if frame.empty:
        raise RuntimeError(f"Analysis set {set_name} is empty")
    outcome = pd.to_numeric(frame[outcome_column], errors="raise")
    if not outcome.isin([0, 1]).all():
        raise RuntimeError(f"Non-binary outcome in analysis set {set_name}")
    frame["analysis_set"] = set_name
    frame["outcome_name"] = outcome_column
    frame["outcome_value"] = outcome.astype(int)
    return frame


def build_loto(
    joined: pd.DataFrame, config: dict
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    manifest_rows: list[dict] = []
    assignment_frames: list[pd.DataFrame] = []
    exclusion_rows: list[dict] = []
    qc_folds: list[dict] = []

    for set_name, set_spec in config["analysis_sets"].items():
        subset = analysis_subset(joined, set_name, set_spec)
        targets = sorted(subset["target"].astype(str).unique())
        for rule in config["sequence_cluster_rules"]:
            for held_target in targets:
                fold_id = f"loto|{set_name}|{rule}|{held_target}"
                frame = subset.copy()
                frame["cluster_rule"] = rule
                frame["sequence_cluster_value"] = frame[rule].astype(str)
                is_test = frame["target"].astype(str).eq(held_target)
                test_clusters = set(
                    frame.loc[is_test, "sequence_cluster_value"].astype(str)
                )
                is_candidate = ~is_test
                is_removed = is_candidate & frame[
                    "sequence_cluster_value"
                ].astype(str).isin(test_clusters)
                is_train = is_candidate & ~is_removed
                if not is_test.any():
                    raise RuntimeError(f"No test rows for LOTO fold {fold_id}")
                if not is_train.any():
                    raise RuntimeError(f"No training rows for LOTO fold {fold_id}")

                test_pairs_by_cluster: dict[str, list[str]] = defaultdict(list)
                for row in frame.loc[is_test].itertuples(index=False):
                    test_pairs_by_cluster[
                        str(getattr(row, rule))
                    ].append(str(row.pair_id))

                counts = (
                    frame.loc[is_train]
                    .groupby("target", dropna=False)
                    .size()
                    .astype(float)
                )
                raw_map = {
                    str(target): 1.0 / float(count)
                    for target, count in counts.items()
                }
                raw_values = frame["target"].astype(str).map(raw_map)
                train_mean = float(raw_values.loc[is_train].mean())
                if not math.isfinite(train_mean) or train_mean <= 0:
                    raise RuntimeError(f"Invalid target-equal mean for {fold_id}")

                assignments = frame[
                    [
                        "pair_id",
                        "analysis_set",
                        "outcome_name",
                        "target",
                        "outcome_value",
                        "cluster_rule",
                        "sequence_cluster_value",
                    ]
                ].copy()
                assignments.insert(3, "fold_id", fold_id)
                assignments.insert(4, "held_out_target", held_target)
                assignments["split_status"] = np.select(
                    [is_train.to_numpy(), is_test.to_numpy()],
                    ["train", "test"],
                    default="removed_cross_target_cluster",
                )
                assignments["weight_pair"] = np.where(is_train, 1.0, np.nan)
                assignments["weight_target_equal_raw"] = np.where(
                    is_train, raw_values, np.nan
                )
                assignments["weight_target_equal"] = np.where(
                    is_train, raw_values / train_mean, np.nan
                )
                assignment_frames.append(assignments)

                removed_rows = frame.loc[is_removed]
                for row in removed_rows.itertuples(index=False):
                    cluster = str(getattr(row, rule))
                    shared_pairs = test_pairs_by_cluster.get(cluster, [])
                    exclusion_rows.append(
                        {
                            "fold_id": fold_id,
                            "analysis_set": set_name,
                            "cluster_rule": rule,
                            "held_out_target": held_target,
                            "pair_id": str(row.pair_id),
                            "source_target": str(row.target),
                            "sequence_cluster_value": cluster,
                            "n_test_pairs_in_shared_cluster": len(shared_pairs),
                            "shared_test_pair_ids_sample": ";".join(
                                sorted(shared_pairs)[:10]
                            ),
                            "removal_reason": (
                                "training_pair_cluster_occurs_in_held_out_target"
                            ),
                        }
                    )

                test_clusters_after = set(
                    assignments.loc[
                        assignments["split_status"].eq("test"),
                        "sequence_cluster_value",
                    ].astype(str)
                )
                train_clusters_after = set(
                    assignments.loc[
                        assignments["split_status"].eq("train"),
                        "sequence_cluster_value",
                    ].astype(str)
                )
                test_targets = set(
                    assignments.loc[
                        assignments["split_status"].eq("test"), "target"
                    ].astype(str)
                )
                per_target_raw_sums = (
                    assignments.loc[
                        assignments["split_status"].eq("train")
                    ]
                    .groupby("target")["weight_target_equal_raw"]
                    .sum()
                )
                raw_sum_error = float(
                    np.max(
                        np.abs(
                            per_target_raw_sums.to_numpy(dtype=float) - 1.0
                        )
                    )
                )
                normalized_mean_error = float(
                    abs(
                        assignments.loc[
                            assignments["split_status"].eq("train"),
                            "weight_target_equal",
                        ].mean()
                        - 1.0
                    )
                )
                tests = {
                    "test_pairs_not_removed": bool(
                        (assignments.loc[assignments["split_status"].eq("test"), "weight_pair"].isna()).all()
                    ),
                    "train_test_cluster_intersection_empty": not bool(
                        train_clusters_after & test_clusters_after
                    ),
                    "test_target_pure": test_targets == {held_target},
                    "target_equal_raw_sums_to_one_within_source_target": raw_sum_error
                    < 1e-12,
                    "target_equal_mean_is_one": normalized_mean_error < 1e-12,
                    "all_analysis_pairs_accounted_once": int(
                        len(assignments)
                    )
                    == len(frame),
                }
                passed = all(tests.values())
                qc_folds.append(
                    {
                        "fold_kind": "loto",
                        "fold_id": fold_id,
                        "analysis_set": set_name,
                        "cluster_rule": rule,
                        "held_out_target": held_target,
                        "passed": passed,
                        "tests": tests,
                        "raw_sum_max_abs_error": raw_sum_error,
                        "normalized_mean_abs_error": normalized_mean_error,
                    }
                )
                manifest_rows.append(
                    {
                        "fold_id": fold_id,
                        "analysis_set": set_name,
                        "outcome": set_spec["outcome"],
                        "cluster_rule": rule,
                        "held_out_target": held_target,
                        "n_pairs_total": int(len(frame)),
                        "n_test_pairs": int(is_test.sum()),
                        "n_train_pairs_candidate": int(is_candidate.sum()),
                        "n_removed_cross_target_cluster": int(
                            is_removed.sum()
                        ),
                        "n_train_pairs_final": int(is_train.sum()),
                        "n_test_positive": int(
                            frame.loc[is_test, "outcome_value"].sum()
                        ),
                        "n_train_positive": int(
                            frame.loc[is_train, "outcome_value"].sum()
                        ),
                        "n_source_targets_in_train": int(
                            frame.loc[is_train, "target"].nunique()
                        ),
                        "n_test_sequence_clusters": len(test_clusters_after),
                        "n_train_sequence_clusters": len(train_clusters_after),
                        "n_removed_sequence_clusters": len(
                            test_clusters_after - train_clusters_after
                        ),
                        "target_equal_mean_before_normalization": train_mean,
                        "qc_passed": passed,
                    }
                )

    assignments = pd.concat(assignment_frames, ignore_index=True)
    manifest = pd.DataFrame(manifest_rows)
    exclusions = pd.DataFrame(
        exclusion_rows,
        columns=[
            "fold_id",
            "analysis_set",
            "cluster_rule",
            "held_out_target",
            "pair_id",
            "source_target",
            "sequence_cluster_value",
            "n_test_pairs_in_shared_cluster",
            "shared_test_pair_ids_sample",
            "removal_reason",
        ],
    )
    qc = {
        "fold_kind": "loto",
        "n_folds": len(manifest),
        "n_passed": int(manifest["qc_passed"].sum()),
        "n_failed": int((~manifest["qc_passed"]).sum()),
        "n_exclusion_rows": len(exclusions),
        "all_fold_ids_unique": manifest["fold_id"].is_unique,
        "folds": qc_folds,
    }
    return manifest, assignments, exclusions, qc


def build_target_cv5(
    joined: pd.DataFrame, config: dict
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    manifest_rows: list[dict] = []
    assignment_frames: list[pd.DataFrame] = []
    qc_folds: list[dict] = []
    n_splits = int(config["target_internal_cv"]["n_splits"])
    seed = int(config["target_internal_cv"]["random_state"])

    for set_name, set_spec in config["analysis_sets"].items():
        subset = analysis_subset(joined, set_name, set_spec)
        targets = sorted(subset["target"].astype(str).unique())
        for rule in config["sequence_cluster_rules"]:
            for target in targets:
                target_frame = subset[
                    subset["target"].astype(str).eq(target)
                ].copy()
                fold_errors: list[str] = []
                try:
                    splitter = StratifiedGroupKFold(
                        n_splits=n_splits,
                        shuffle=True,
                        random_state=seed,
                    )
                    splits = list(
                        splitter.split(
                            np.zeros(len(target_frame), dtype=np.int8),
                            target_frame["outcome_value"].to_numpy(dtype=int),
                            target_frame[rule].astype(str).to_numpy(),
                        )
                    )
                except Exception as error:  # Explicit QC failure; no fallback.
                    splits = []
                    fold_errors.append(f"{type(error).__name__}: {error}")
                if len(splits) != n_splits and not fold_errors:
                    fold_errors.append(
                        f"Expected {n_splits} splits, received {len(splits)}"
                    )

                fold_assignments: list[pd.DataFrame] = []
                for fold_number, (train_index, test_index) in enumerate(
                    splits, start=1
                ):
                    fold_id = (
                        f"cv5|{set_name}|{rule}|{target}|fold_{fold_number:02d}"
                    )
                    status = np.full(len(target_frame), "train", dtype=object)
                    status[test_index] = "test"
                    train_frame = target_frame.iloc[train_index]
                    test_frame = target_frame.iloc[test_index]
                    train_groups = set(
                        train_frame[rule].astype(str)
                    )
                    test_groups = set(test_frame[rule].astype(str))
                    accounted = np.zeros(len(target_frame), dtype=int)
                    accounted[train_index] += 1
                    accounted[test_index] += 1
                    checks = {
                        "each_pair_assigned_once": bool((accounted == 1).all()),
                        "group_intersection_empty": not bool(
                            train_groups & test_groups
                        ),
                        "test_target_pure": set(
                            test_frame["target"].astype(str)
                        )
                        == {target},
                        "train_has_two_outcome_classes": train_frame[
                            "outcome_value"
                        ].nunique()
                        == 2,
                        "test_has_two_outcome_classes": test_frame[
                            "outcome_value"
                        ].nunique()
                        == 2,
                    }
                    hard_checks = [
                        checks["each_pair_assigned_once"],
                        checks["group_intersection_empty"],
                        checks["test_target_pure"],
                    ]
                    soft_checks = [
                        checks["train_has_two_outcome_classes"],
                        checks["test_has_two_outcome_classes"],
                    ]
                    passed = all(hard_checks) and all(soft_checks)
                    qc_folds.append(
                        {
                            "fold_kind": "target_cv5",
                            "fold_id": fold_id,
                            "analysis_set": set_name,
                            "cluster_rule": rule,
                            "target": target,
                            "passed": passed,
                            "tests": checks,
                            "failure_reason": None
                            if passed
                            else (
                                "class_or_structure_qc_failed_without_fallback"
                            ),
                        }
                    )
                    manifest_rows.append(
                        {
                            "fold_id": fold_id,
                            "analysis_set": set_name,
                            "outcome": set_spec["outcome"],
                            "cluster_rule": rule,
                            "target": target,
                            "n_pairs_total": len(target_frame),
                            "n_train_pairs": len(train_index),
                            "n_test_pairs": len(test_index),
                            "n_train_positive": int(
                                train_frame["outcome_value"].sum()
                            ),
                            "n_test_positive": int(
                                test_frame["outcome_value"].sum()
                            ),
                            "n_train_sequence_clusters": len(train_groups),
                            "n_test_sequence_clusters": len(test_groups),
                            "qc_passed": passed,
                        }
                    )
                    fold_frame = target_frame[
                        [
                            "pair_id",
                            "analysis_set",
                            "outcome_name",
                            "target",
                            "outcome_value",
                            rule,
                        ]
                    ].copy()
                    fold_frame.insert(3, "fold_id", fold_id)
                    fold_frame.insert(4, "cluster_rule", rule)
                    fold_frame["split"] = status
                    fold_frame = fold_frame.rename(
                        columns={rule: "sequence_cluster_value"}
                    )
                    fold_assignments.append(fold_frame)

                if fold_assignments:
                    combined = pd.concat(fold_assignments, ignore_index=True)
                    assignments_for_target = combined
                    assignment_frames.append(assignments_for_target)
                if fold_errors:
                    qc_folds.append(
                        {
                            "fold_kind": "target_cv5",
                            "fold_id": (
                                f"cv5|{set_name}|{rule}|{target}|ERROR"
                            ),
                            "analysis_set": set_name,
                            "cluster_rule": rule,
                            "target": target,
                            "passed": False,
                            "tests": {},
                            "failure_reason": "; ".join(fold_errors),
                        }
                    )
                    manifest_rows.append(
                        {
                            "fold_id": f"cv5|{set_name}|{rule}|{target}|ERROR",
                            "analysis_set": set_name,
                            "outcome": set_spec["outcome"],
                            "cluster_rule": rule,
                            "target": target,
                            "n_pairs_total": len(target_frame),
                            "n_train_pairs": None,
                            "n_test_pairs": None,
                            "n_train_positive": None,
                            "n_test_positive": None,
                            "n_train_sequence_clusters": None,
                            "n_test_sequence_clusters": None,
                            "qc_passed": False,
                        }
                    )

    assignments = (
        pd.concat(assignment_frames, ignore_index=True)
        if assignment_frames
        else pd.DataFrame()
    )
    manifest = pd.DataFrame(manifest_rows)
    qc = {
        "fold_kind": "target_cv5",
        "n_splitter_runs": int(
            manifest["fold_id"].astype(str).str.endswith("ERROR").eq(False).sum()
        ),
        "n_fold_records": len(
            manifest[
                ~manifest["fold_id"].astype(str).str.endswith("ERROR")
            ]
        ),
        "n_passed": int(
            (
                manifest["qc_passed"]
                & ~manifest["fold_id"].astype(str).str.endswith("ERROR")
            ).sum()
        ),
        "n_failed": int(
            (
                ~manifest["qc_passed"]
                | manifest["fold_id"].astype(str).str.endswith("ERROR")
            ).sum()
        ),
        "n_splits": n_splits,
        "shuffle": True,
        "random_state": seed,
        "folds": qc_folds,
    }
    return manifest, assignments, qc


def extract_evaluation_metadata(raw_text: object) -> dict:
    if raw_text is None or (isinstance(raw_text, float) and math.isnan(raw_text)):
        return {
            "parse_ok": False,
            "parse_error": "missing evaluations JSON",
            "experimental_methods": [],
            "experimental_metrics": [],
            "experimental_targets": [],
            "computational_methods": [],
            "experimental_date_min": None,
            "experimental_date_max": None,
            "iptm_metric_present": False,
        }
    try:
        values = json.loads(str(raw_text))
        if not isinstance(values, list):
            raise ValueError("evaluations is not a JSON list")
    except Exception as error:
        return {
            "parse_ok": False,
            "parse_error": f"{type(error).__name__}: {error}",
            "experimental_methods": [],
            "experimental_metrics": [],
            "experimental_targets": [],
            "computational_methods": [],
            "experimental_date_min": None,
            "experimental_date_max": None,
            "iptm_metric_present": False,
        }
    experimental_methods = set()
    experimental_metrics = set()
    experimental_targets = set()
    computational_methods = set()
    dates = []
    iptm_metric_present = False
    for item in values:
        if not isinstance(item, dict):
            continue
        record_type = item.get("type")
        method = item.get("method")
        metric = item.get("metric")
        target = item.get("target")
        measured_at = item.get("measuredAt")
        if method not in (None, ""):
            if record_type == "experimental":
                experimental_methods.add(str(method))
            elif record_type == "computational":
                computational_methods.add(str(method))
        if record_type == "experimental" and metric not in (None, ""):
            experimental_metrics.add(str(metric))
        if record_type == "experimental" and target not in (None, ""):
            experimental_targets.add(str(target))
        if metric is not None and "iptm" in str(metric).casefold():
            iptm_metric_present = True
        if record_type == "experimental" and measured_at not in (None, ""):
            dates.append(str(measured_at))
    return {
        "parse_ok": True,
        "parse_error": None,
        "experimental_methods": sorted(experimental_methods),
        "experimental_metrics": sorted(experimental_metrics),
        "experimental_targets": sorted(experimental_targets),
        "computational_methods": sorted(computational_methods),
        "experimental_date_min": min(dates) if dates else None,
        "experimental_date_max": max(dates) if dates else None,
        "iptm_metric_present": iptm_metric_present,
    }


def normalized_text(value: object) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    text = str(value).strip()
    return text if text else "unknown"


def build_source_annotations(
    joined: pd.DataFrame, config: dict
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    raw = pd.read_csv(PROTEINBASE_PATH, low_memory=False)
    if raw["id"].nunique() != len(raw):
        raise RuntimeError("proteinbase_latest.csv id is not unique")
    if len(joined) != int(config["inputs"]["pair_master"]["expected_rows"]):
        raise RuntimeError("Unexpected main-table row count for source audit")

    id_to_index = {str(value): index for index, value in enumerate(raw["id"])}
    missing_ids = [
        value
        for value in joined["design_id"].astype(str)
        if value not in id_to_index
    ]
    if missing_ids:
        raise RuntimeError(
            f"{len(missing_ids)} design_id values are absent from Proteinbase"
        )
    raw_index = joined["design_id"].astype(str).map(id_to_index).astype(int)
    raw_rows = raw.iloc[raw_index.to_numpy()].reset_index(drop=True)
    physical_line = raw_index.to_numpy() + 2
    sequence_match = joined["sequence"].astype(str).str.strip().str.upper().eq(
        raw_rows["sequence"].astype(str).str.strip().str.upper()
    )
    source_row_match = physical_line == joined["source_row"].to_numpy()
    design_id_match = joined["design_id"].astype(str).eq(
        raw_rows["id"].astype(str)
    )
    if not sequence_match.all():
        raise RuntimeError(
            f"{int((~sequence_match).sum())} design_id/sequence links disagree"
        )
    if not source_row_match.all():
        examples = joined.loc[
            ~source_row_match, ["pair_id", "design_id", "source_row"]
        ].head(5)
        raise RuntimeError(
            "source_row does not equal the validated CSV physical data line "
            f"for {int((~source_row_match).sum())} rows; examples:\n{examples}"
        )
    if not design_id_match.all():
        raise RuntimeError("Validated design_id linkage changed unexpectedly")

    author_match = raw_rows["author"].map(normalized_text) == joined[
        "author"
    ].map(normalized_text)
    design_method_match = raw_rows["designMethod"].map(
        normalized_text
    ) == joined["designMethod"].map(normalized_text)

    annotations = pd.DataFrame(
        {
            "pair_id": joined["pair_id"].astype(str),
            "design_id": joined["design_id"].astype(str),
            "target": joined["target"].astype(str),
            "source_row": joined["source_row"].astype(int),
            "verified_csv_physical_line": physical_line,
            "source_row_verified": source_row_match,
            "design_id_verified": design_id_match,
            "sequence_verified": sequence_match,
            "author_master": joined["author"],
            "author_raw": raw_rows["author"],
            "author_exact_match": author_match,
            "designMethod_master": joined["designMethod"],
            "designMethod_raw": raw_rows["designMethod"],
            "designMethod_exact_match": design_method_match,
            "design_class": joined["design_class"],
            "primary_assay": joined["primary_assay"],
            "assays_observed": joined["assays_observed"],
            "n_assays_observed": joined["n_assays_observed"],
            "primary_outcome": joined["primary_outcome"],
            "strict_outcome": joined["strict_outcome"],
            "raw_name": raw_rows["name"],
            "raw_sequence_sha256": raw_rows["sequence"].map(
                lambda value: sha256_text(str(value))
            ),
            "raw_evaluations_sha256": raw_rows["evaluations"].map(
                lambda value: sha256_text(str(value))
            ),
            "raw_evaluations_bytes": raw_rows["evaluations"].map(
                lambda value: len(str(value).encode("utf-8"))
            ),
        }
    )

    source_proxy = annotations["author_raw"].map(normalized_text)
    campaign_proxy = [
        f"{normalized_text(author)}|{normalized_text(target)}"
        for author, target in zip(
            annotations["author_raw"], annotations["target"], strict=True
        )
    ]
    annotations.insert(22, "source_proxy", source_proxy)
    annotations.insert(23, "campaign_proxy", campaign_proxy)

    evaluation_records = []
    for value in raw_rows["evaluations"]:
        evaluation_records.append(extract_evaluation_metadata(value))
    annotations["experimental_methods_json"] = [
        json.dumps(record["experimental_methods"], ensure_ascii=False)
        for record in evaluation_records
    ]
    annotations["experimental_metrics_json"] = [
        json.dumps(record["experimental_metrics"], ensure_ascii=False)
        for record in evaluation_records
    ]
    annotations["experimental_targets_json"] = [
        json.dumps(record["experimental_targets"], ensure_ascii=False)
        for record in evaluation_records
    ]
    annotations["computational_methods_json"] = [
        json.dumps(record["computational_methods"], ensure_ascii=False)
        for record in evaluation_records
    ]
    annotations["experimental_date_min"] = [
        record["experimental_date_min"] for record in evaluation_records
    ]
    annotations["experimental_date_max"] = [
        record["experimental_date_max"] for record in evaluation_records
    ]
    annotations["any_iptm_metric_present"] = [
        record["iptm_metric_present"] for record in evaluation_records
    ]

    dictionary_mappings = []
    grouped = annotations.groupby(
        ["author_raw", "target", "source_proxy", "campaign_proxy"],
        dropna=False,
        sort=True,
    )
    for key, group in grouped:
        raw_author, target, source_value, campaign_value = key
        unknown_fields = []
        if normalized_text(raw_author) == "unknown":
            unknown_fields.append("author_raw")
        if normalized_text(target) == "unknown":
            unknown_fields.append("target")
        dictionary_mappings.append(
            {
                "raw_author": raw_author,
                "target": target,
                "source_proxy": source_value,
                "campaign_proxy": campaign_value,
                "pair_count": int(len(group)),
                "unknown_fields": unknown_fields,
            }
        )
    source_dictionary = {
        "version_date": config["handoff_date"],
        "factual_status": (
            "Proxy labels only. No campaign membership, de novo status, or "
            "selection role is inferred beyond traceable raw fields."
        ),
        "rules": {
            "source_proxy": (
                "Use the raw Proteinbase author value exactly, trimmed; "
                "empty or missing becomes 'unknown'."
            ),
            "campaign_proxy": (
                "Concatenate raw Proteinbase author and frozen main-table "
                "target as '<author>|<target>'; either missing component "
                "remains literally 'unknown'. This is a deterministic "
                "author-target proxy, not a verified campaign identifier."
            ),
        },
        "raw_fields_used": [
            "id",
            "name",
            "sequence",
            "author",
            "designMethod",
            "evaluations",
        ],
        "mapping_count": len(dictionary_mappings),
        "unknown_pair_counts": {
            "source_proxy_unknown": int(
                annotations["source_proxy"].eq("unknown").sum()
            ),
            "campaign_proxy_component_unknown": int(
                annotations["campaign_proxy"].str.contains(
                    "unknown", regex=False
                ).sum()
            ),
            "author_raw_missing": int(annotations["author_raw"].isna().sum()),
            "designMethod_raw_missing": int(
                annotations["designMethod_raw"].isna().sum()
            ),
        },
        "mappings": dictionary_mappings,
    }

    summary_rows = []
    for target, group in annotations.groupby("target", sort=True):
        source_counts = group["source_proxy"].value_counts(dropna=False)
        campaign_counts = group["campaign_proxy"].value_counts(dropna=False)
        assay_counts = group["primary_assay"].map(normalized_text).value_counts()
        design_class_counts = group["design_class"].map(
            normalized_text
        ).value_counts()
        summary_rows.append(
            {
                "target": target,
                "n_pairs": len(group),
                "n_primary_positive": int(group["primary_outcome"].sum()),
                "n_strict_positive": int(group["strict_outcome"].sum()),
                "n_authors": int(
                    group["author_raw"].map(normalized_text).nunique()
                ),
                "n_designMethods": int(
                    group["designMethod_raw"]
                    .map(normalized_text)
                    .replace("unknown", np.nan)
                    .nunique()
                ),
                "n_source_proxies": int(source_counts.size),
                "n_campaign_proxies": int(campaign_counts.size),
                "max_source_proxy": str(source_counts.index[0]),
                "max_source_share": float(source_counts.iloc[0] / len(group)),
                "max_campaign_proxy": str(campaign_counts.index[0]),
                "max_campaign_share": float(campaign_counts.iloc[0] / len(group)),
                "assay_composition_json": json.dumps(
                    assay_counts.to_dict(), ensure_ascii=False
                ),
                "design_class_composition_json": json.dumps(
                    design_class_counts.to_dict(), ensure_ascii=False
                ),
                "author_composition_json": json.dumps(
                    group["author_raw"].map(normalized_text).value_counts().to_dict(),
                    ensure_ascii=False,
                ),
            }
        )
    target_summary = pd.DataFrame(summary_rows)

    evaluation_parse_failures = sum(
        not record["parse_ok"] for record in evaluation_records
    )
    qc = {
        "created_at": now_iso(),
        "pair_master_rows": len(joined),
        "pair_master_unique_pair_ids": int(joined["pair_id"].nunique()),
        "proteinbase_rows_parsed": len(raw),
        "proteinbase_unique_ids": int(raw["id"].nunique()),
        "proteinbase_unique_sequences": int(raw["sequence"].nunique()),
        "pair_design_ids_linked_to_proteinbase": int(design_id_match.sum()),
        "link_rate": float(design_id_match.mean()),
        "design_id_sequence_validation_passed": bool(sequence_match.all()),
        "source_row_validation_passed": bool(source_row_match.all()),
        "source_row_rule": (
            "source_row equals the CSV physical line number (header=1, "
            "first data line=2), validated only after unique design_id and "
            "exact sequence matching."
        ),
        "source_row_max_abs_difference": int(
            np.max(np.abs(physical_line - joined["source_row"].to_numpy()))
        ),
        "author_exact_match_count": int(author_match.sum()),
        "author_mismatch_count": int((~author_match).sum()),
        "designMethod_exact_match_count": int(design_method_match.sum()),
        "designMethod_mismatch_count": int((~design_method_match).sum()),
        "missing_master_author_count": int(joined["author"].isna().sum()),
        "missing_master_designMethod_count": int(
            joined["designMethod"].isna().sum()
        ),
        "missing_raw_author_count": int(raw_rows["author"].isna().sum()),
        "missing_raw_designMethod_count": int(
            raw_rows["designMethod"].isna().sum()
        ),
        "duplicate_pair_design_ids": int(
            joined["design_id"].duplicated().sum()
        ),
        "duplicate_raw_ids": int(raw["id"].duplicated().sum()),
        "duplicate_raw_sequences": int(raw["sequence"].duplicated().sum()),
        "evaluations_parse_failures": int(evaluation_parse_failures),
        "unknown_source_proxy_count": int(
            annotations["source_proxy"].eq("unknown").sum()
        ),
        "unknown_campaign_proxy_component_count": int(
            annotations["campaign_proxy"].str.contains(
                "unknown", regex=False
            ).sum()
        ),
        "iptm_metric_present_count": int(
            annotations["any_iptm_metric_present"].sum()
        ),
        "iptm_interpretation": (
            "Presence of an ipTM-named Proteinbase metric is recorded only. "
            "No inference is made that ipTM participated in candidate selection."
        ),
        "unable_to_confirm": [
            "Verified campaign membership is unavailable from the frozen raw fields.",
            "De novo status is not inferred.",
            "Candidate-selection use of ipTM is not inferred.",
        ],
    }
    return annotations, target_summary, qc, source_dictionary


def main() -> None:
    started = time.time()
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    preflight = collect_environment(config)
    joined = prepare_pair_data(config)

    loto_manifest, loto_assignments, exclusions, loto_qc = build_loto(
        joined, config
    )
    cv_manifest, cv_assignments, cv_qc = build_target_cv5(joined, config)
    atomic_csv(RUN_DIR / "loto_fold_manifest.csv", loto_manifest)
    atomic_csv(RUN_DIR / "loto_pair_assignments.csv", loto_assignments)
    atomic_csv(
        RUN_DIR / "cross_target_cluster_exclusions.csv", exclusions
    )
    atomic_csv(RUN_DIR / "target_cv5_fold_manifest.csv", cv_manifest)
    atomic_csv(
        RUN_DIR / "target_cv5_pair_assignments.csv", cv_assignments
    )

    expected_targets_by_set = {
        name: int(
            joined.loc[
                ~joined["target"].astype(str).str.casefold().eq(
                    str(spec.get("exclude_target")).casefold()
                )
                if spec.get("exclude_target")
                else np.ones(len(joined), dtype=bool),
                "target",
            ].nunique()
        )
        for name, spec in config["analysis_sets"].items()
    }
    expected_loto = (
        sum(expected_targets_by_set.values())
        * len(config["sequence_cluster_rules"])
    )
    expected_cv_splits = (
        sum(expected_targets_by_set.values())
        * len(config["sequence_cluster_rules"])
    )
    expected_cv_fold_records = (
        expected_cv_splits * int(config["target_internal_cv"]["n_splits"])
    )
    fold_qc = {
        "created_at": now_iso(),
        "inputs": {
            "pair_master": str(PAIR_PATH),
            "sequence_clusters": str(CLUSTER_PATH),
            "pair_master_sha256": preflight["input_paths_and_hashes"][
                "pair_master"
            ]["actual_sha256"],
            "sequence_clusters_sha256": preflight["input_paths_and_hashes"][
                "sequence_clusters"
            ]["actual_sha256"],
        },
        "joined_rows": len(joined),
        "analysis_set_targets": expected_targets_by_set,
        "loto": {
            **loto_qc,
            "expected_folds": expected_loto,
            "all_expected_folds_present": len(loto_manifest)
            == expected_loto,
            "each_analysis_set_rule_target_has_one_fold": not loto_manifest.duplicated(
                ["analysis_set", "cluster_rule", "held_out_target"]
            ).any(),
        },
        "target_cv5": {
            **cv_qc,
            "expected_splitter_runs": expected_cv_splits,
            "expected_complete_fold_records": expected_cv_fold_records,
            "all_complete_fold_records_present": int(
                cv_manifest["fold_id"]
                .astype(str)
                .str.endswith("ERROR")
                .eq(False)
                .sum()
            )
            == expected_cv_fold_records,
            "nested_inner_folds_generated": False,
            "fallback_used": False,
        },
        "outcome_columns_modified": False,
        "frozen_fields_modified": False,
    }
    atomic_json(RUN_DIR / "fold_qc.json", fold_qc)

    annotations, source_summary, source_qc, source_dictionary = (
        build_source_annotations(joined, config)
    )
    atomic_csv(RUN_DIR / "pair_source_annotations.csv", annotations)
    atomic_csv(RUN_DIR / "target_source_summary.csv", source_summary)
    atomic_json(
        RUN_DIR / "source_campaign_dictionary.json", source_dictionary
    )
    atomic_json(RUN_DIR / "source_metadata_qc.json", source_qc)

    a_passed = bool(
        fold_qc["loto"]["all_expected_folds_present"]
        and fold_qc["loto"]["n_failed"] == 0
        and fold_qc["target_cv5"]["all_complete_fold_records_present"]
        and fold_qc["target_cv5"]["n_failed"] == 0
    )
    b_passed = bool(
        source_qc["source_row_validation_passed"]
        and source_qc["design_id_sequence_validation_passed"]
        and source_qc["pair_design_ids_linked_to_proteinbase"] == len(joined)
        and source_qc["evaluations_parse_failures"] == 0
    )
    phase_summary = {
        "created_at": now_iso(),
        "elapsed_seconds": round(time.time() - started, 3),
        "A": {
            "status": "complete" if a_passed else "partial",
            "loto_folds": len(loto_manifest),
            "loto_assignment_rows": len(loto_assignments),
            "cross_target_exclusions": len(exclusions),
            "target_cv5_fold_records": int(
                cv_manifest["fold_id"]
                .astype(str)
                .str.endswith("ERROR")
                .eq(False)
                .sum()
            ),
            "target_cv5_assignment_rows": len(cv_assignments),
            "qc_passed": a_passed,
            "forbidden_model_fitting_or_metrics_performed": False,
        },
        "B": {
            "status": "complete" if b_passed else "partial",
            "pair_annotations": len(annotations),
            "targets_summarized": len(source_summary),
            "source_row_link_rate": source_qc["link_rate"],
            "qc_passed": b_passed,
        },
    }
    atomic_json(RUN_DIR / "phase_AB_summary.json", phase_summary)
    print(json.dumps(phase_summary, ensure_ascii=False, indent=2))
    if not a_passed or not b_passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
