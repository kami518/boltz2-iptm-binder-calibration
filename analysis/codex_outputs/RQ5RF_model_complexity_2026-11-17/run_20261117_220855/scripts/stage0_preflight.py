"""RQ5-RF Stage 0: validate frozen inputs and save the execution snapshot."""
from __future__ import annotations

import json
import os
import time
import ctypes
from ctypes import wintypes

import numpy as np
import pandas as pd

from rq5rf_common import (
    ANALYSIS_SETS,
    CONFIG_PATH,
    RUN_DIR,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_json,
    current_source_weights,
    effective_rf_grid,
    ensure_failure_ledger,
    environment_snapshot,
    load_config,
    load_joined_inputs,
    now_iso,
    resolve_failure,
    set_stage_status,
    sha256_file,
    stable_sort,
    update_status,
    weighted_center_scale,
)

STAGE = "stage_0_preflight"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class _ProcessMemoryCounters(ctypes.Structure):
    _fields_ = [
        ("cb", wintypes.DWORD),
        ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


def windows_process_memory() -> dict[str, int]:
    access = 0x0400 | 0x0010
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    psapi.GetProcessMemoryInfo.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(_ProcessMemoryCounters),
        wintypes.DWORD,
    ]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    handle = kernel32.OpenProcess(access, False, os.getpid())
    if not handle:
        raise OSError(ctypes.get_last_error(), "OpenProcess failed")
    try:
        counters = _ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        if not psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.cb):
            raise OSError(ctypes.get_last_error(), "GetProcessMemoryInfo failed")
        return {
            "rss_bytes": int(counters.WorkingSetSize),
            "peak_working_set_bytes": int(counters.PeakWorkingSetSize),
            "peak_pagefile_bytes": int(counters.PeakPagefileUsage),
        }
    finally:
        kernel32.CloseHandle(handle)


def main() -> None:
    ensure_failure_ledger()
    update_status(
        overall_status="running",
        status_semantics={
            "complete": "all required stages and QC completed",
            "blocked": "frozen-rule or environment blocker requiring renewal",
            "failed": "code error remains after permitted repair attempts",
        },
    )
    set_stage_status(STAGE, "running")
    config = load_config()
    environment = environment_snapshot()
    atomic_json(RUN_DIR / "environment.json", environment)
    require(
        environment["scikit_learn"] == "1.9.0",
        f"Blocked: scikit-learn must be 1.9.0, found {environment['scikit_learn']}",
    )

    hash_rows = []
    protocol_path = RUN_DIR / "RQ5RF_model_complexity_frozen_protocol_2026-11-17.md"
    handoff_path = RUN_DIR / "START_HERE_CODEX.md"
    source_protocol_path = config["protocol"]["path"]
    require(sha256_file(protocol_path) == config["protocol"]["sha256"].lower(), "Protocol copy hash mismatch")
    require(sha256_file(source_protocol_path) == config["protocol"]["sha256"].lower(), "Source protocol hash mismatch")
    hash_rows.extend(
        [
            {
                "role": "frozen_protocol",
                "path": source_protocol_path,
                "expected_sha256": config["protocol"]["sha256"],
                "actual_sha256": sha256_file(source_protocol_path),
                "status": "pass",
            },
            {
                "role": "execution_handoff_copy",
                "path": str(handoff_path),
                "expected_sha256": "",
                "actual_sha256": sha256_file(handoff_path),
                "status": "pass",
            },
        ]
    )

    for name, spec in config["inputs"].items():
        path = spec["path"]
        actual = sha256_file(path)
        expected = spec["sha256"].lower()
        require(actual == expected, f"Input hash mismatch for {name}: {actual} != {expected}")
        hash_rows.append(
            {
                "role": name,
                "path": path,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "status": "pass",
            }
        )
    atomic_csv(RUN_DIR / "final_input_hashes.csv", pd.DataFrame(hash_rows))

    master, pinc, clusters, joined = load_joined_inputs(config)
    success = int(pinc["pinc_status"].eq("success").sum())
    missing = int(pinc["pinc_status"].eq("missing_input").sum())
    require(success == 2074, f"Pinc success count mismatch: {success}")
    require(missing == 10, f"Pinc missing count mismatch: {missing}")
    require(
        set(pinc.loc[pinc["pinc_status"].ne("success"), "target"].astype(str)) == {"pd-l1"},
        "Pinc missing rows are not exclusively pd-l1",
    )

    set_rows = []
    frames = {}
    for analysis_set in ANALYSIS_SETS:
        frame = analysis_frame(config, joined, analysis_set)
        frames[analysis_set] = frame
        expected = config["analysis_sets"][analysis_set]
        set_rows.append(
            {
                "analysis_set": analysis_set,
                "outcome": expected["outcome"],
                "excluded_target": expected["exclude_target"],
                "rows": len(frame),
                "positives": int(frame["outcome_value"].sum()),
                "targets": int(frame["target"].nunique()),
                "unique_pair_ids": int(frame["pair_id"].nunique()),
                "iptm_min": float(frame["boltz2_iptm"].min()),
                "iptm_max": float(frame["boltz2_iptm"].max()),
                "pinc_min": float(frame["pinc"].min()),
                "pinc_max": float(frame["pinc"].max()),
                "status": "pass",
            }
        )
    atomic_csv(RUN_DIR / "analysis_set_audit.csv", pd.DataFrame(set_rows))

    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    loto_manifest = pd.read_csv(config["inputs"]["loto_manifest"]["path"], low_memory=False)
    cv5 = pd.read_csv(config["inputs"]["target_cv5_assignments"]["path"], low_memory=False)
    cv5_manifest = pd.read_csv(config["inputs"]["target_cv5_manifest"]["path"], low_memory=False)
    exclusions = pd.read_csv(config["inputs"]["cross_target_exclusions"]["path"], low_memory=False)
    require(len(loto_manifest) == 44, "LOTO manifest row count mismatch")
    require(len(cv5_manifest) == 220, "Target-CV5 manifest row count mismatch")
    cluster_rule = config["cluster_rule"]
    loto90 = loto[loto["cluster_rule"].eq(cluster_rule)].copy()
    loto_manifest90 = loto_manifest[loto_manifest["cluster_rule"].eq(cluster_rule)].copy()
    cv590 = cv5[cv5["cluster_rule"].eq(cluster_rule)].copy()
    cv5_manifest90 = cv5_manifest[cv5_manifest["cluster_rule"].eq(cluster_rule)].copy()
    exclusions90 = exclusions[exclusions["cluster_rule"].eq(cluster_rule)].copy()
    require(
        len(loto_manifest90[loto_manifest90["analysis_set"].isin(ANALYSIS_SETS)]) == 17,
        "Expected 17 selected-analysis-set frozen LOTO90 source tasks",
    )
    require(
        len(cv5_manifest90[cv5_manifest90["analysis_set"].isin(ANALYSIS_SETS)]) == 85,
        "Expected 85 selected-analysis-set frozen target-CV5 folds",
    )

    fold_rows = []
    for analysis_set in ANALYSIS_SETS:
        valid_ids = set(frames[analysis_set]["pair_id"].astype(str))
        manifest_set = loto_manifest90[loto_manifest90["analysis_set"].eq(analysis_set)]
        require(len(manifest_set) == int(frames[analysis_set]["target"].nunique()), f"Unexpected LOTO target count for {analysis_set}")
        for held_target in sorted(manifest_set["held_out_target"].astype(str)):
            fold_id = f"loto|{analysis_set}|{cluster_rule}|{held_target}"
            fold = loto90[loto90["fold_id"].eq(fold_id)].copy()
            require(len(fold) > 0, f"Missing LOTO fold {fold_id}")
            require(fold["pair_id"].is_unique, f"Duplicate pair rows in {fold_id}")
            fold = fold[fold["pair_id"].astype(str).isin(valid_ids)].copy()
            train = fold[fold["split_status"].eq("train")]
            test = fold[fold["split_status"].eq("test")]
            require(len(test) == int((frames[analysis_set]["target"].astype(str) == held_target).sum()), f"LOTO test does not cover held target {fold_id}")
            require(set(test["target"].astype(str)) == {held_target}, f"Wrong LOTO test target {fold_id}")
            require(held_target not in set(train["target"].astype(str)), f"Source target leak {fold_id}")
            require(
                not (
                    set(train["sequence_cluster_value"].astype(str))
                    & set(test["sequence_cluster_value"].astype(str))
                ),
                f"LOTO cluster leak {fold_id}",
            )
            explicit_exclusions = set(
                exclusions90.loc[exclusions90["fold_id"].eq(fold_id), "pair_id"].astype(str)
            )
            require(not explicit_exclusions.intersection(set(train["pair_id"].astype(str))), f"Frozen exclusion remains in train {fold_id}")
            fold_rows.append(
                {
                    "analysis_set": analysis_set,
                    "fold_id": fold_id,
                    "held_out_target": held_target,
                    "train_rows_after_pinc": len(train),
                    "test_rows_after_pinc": len(test),
                    "source_targets": train["target"].nunique(),
                    "train_positives": int(train["outcome_value"].sum()),
                    "test_positives": int(test["outcome_value"].sum()),
                    "frozen_exclusion_rows_total": len(explicit_exclusions),
                    "status": "pass",
                }
            )
    atomic_csv(RUN_DIR / "loto_fold_audit.csv", pd.DataFrame(fold_rows))

    cv_rows = []
    for analysis_set in ANALYSIS_SETS:
        valid_ids = set(frames[analysis_set]["pair_id"].astype(str))
        subset = cv590[cv590["analysis_set"].eq(analysis_set)].copy()
        subset = subset[subset["pair_id"].astype(str).isin(valid_ids)]
        expected_manifest = cv5_manifest90[cv5_manifest90["analysis_set"].eq(analysis_set)]
        require(
            subset["fold_id"].nunique() == len(expected_manifest),
            f"Target-CV5 fold count mismatch for {analysis_set}",
        )
        for fold_id, fold in subset.groupby("fold_id", sort=True):
            train = fold[fold["split"].eq("train")]
            test = fold[fold["split"].eq("test")]
            require(len(train) > 0 and len(test) > 0, f"Empty target-CV5 fold {fold_id}")
            require(fold["pair_id"].is_unique, f"Duplicate target-CV5 row in {fold_id}")
            require(
                not (
                    set(train["sequence_cluster_value"].astype(str))
                    & set(test["sequence_cluster_value"].astype(str))
                ),
                f"Target-CV5 cluster leak {fold_id}",
            )
            require(train["outcome_value"].nunique() == 2, f"Single-class target-CV5 train {fold_id}")
            cv_rows.append(
                {
                    "analysis_set": analysis_set,
                    "fold_id": fold_id,
                    "target": str(train["target"].iloc[0]),
                    "train_rows": len(train),
                    "test_rows": len(test),
                    "train_positives": int(train["outcome_value"].sum()),
                    "test_positives": int(test["outcome_value"].sum()),
                    "status": "pass",
                }
            )
    atomic_csv(RUN_DIR / "target_cv5_fold_audit.csv", pd.DataFrame(cv_rows))

    weight_rows = []
    for analysis_set in ANALYSIS_SETS:
        frame = frames[analysis_set]
        for weighting in config["source_weightings"]:
            _, audit = current_source_weights(frame, weighting)
            weight_rows.append(
                {
                    "analysis_set": analysis_set,
                    "source_weighting": weighting,
                    "n_rows": len(frame),
                    **audit,
                    "status": "pass",
                }
            )
    atomic_csv(RUN_DIR / "global_weight_audit.csv", pd.DataFrame(weight_rows))

    # Representative full-size 500-tree fit used only for resource planning.
    representative = frames["primary_six"]
    representative = representative.sample(n=min(1600, len(representative)), random_state=20261117)
    weights, _ = current_source_weights(representative, "target_equal")
    mean, sd = weighted_center_scale(
        representative[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float), weights
    )
    del mean, sd
    parameters = effective_rf_grid(2)[0]
    try:
        memory_before_info = windows_process_memory()
        memory_before = int(memory_before_info["rss_bytes"])
        started = time.perf_counter()
        from sklearn.ensemble import RandomForestClassifier

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
            n_jobs=int(config["random_forest"]["n_jobs"]),
            random_state=20261117,
            verbose=0,
            warm_start=False,
            class_weight=None,
            max_samples=None,
            ccp_alpha=0.0,
        ).fit(
            representative[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float),
            representative["outcome_value"].to_numpy(dtype=int),
            sample_weight=weights,
        )
        elapsed = float(time.perf_counter() - started)
        memory_after_info = windows_process_memory()
        memory_after = int(memory_after_info["rss_bytes"])
        peak = int(memory_after_info["peak_working_set_bytes"])
        del estimator
        representative_resource = {
            "status": "complete",
            "n_rows": int(len(representative)),
            "n_trees": 500,
            "n_jobs": int(config["random_forest"]["n_jobs"]),
            "elapsed_seconds": elapsed,
            "rss_before_bytes": memory_before,
            "rss_after_bytes": memory_after,
            "process_peak_working_set_bytes": peak,
            "process_peak_pagefile_bytes": int(
                memory_after_info["peak_pagefile_bytes"]
            ),
            "rss_delta_bytes": memory_after - memory_before,
            "note": "Representative resource probe only; not a result fit.",
        }
        resolve_failure(
            STAGE,
            "representative_rf_resource_probe",
            "resource_probe_failure",
        )
    except Exception as error:
        representative_resource = {
            "status": "failed",
            "failure": f"{type(error).__name__}: {error}",
        }
        append_failure(
            STAGE,
            "representative_rf_resource_probe",
            "resource_probe_failure",
            str(error),
        )
    atomic_json(RUN_DIR / "logs" / "representative_fit_resource.json", representative_resource)

    summary = {
        "stage": STAGE,
        "created_at": now_iso(),
        "status": "complete",
        "scikit_learn": environment["scikit_learn"],
        "protocol_sha256": sha256_file(source_protocol_path),
        "analysis_sets": set_rows,
        "frozen_loto90_tasks": int(
            loto_manifest90["analysis_set"].isin(ANALYSIS_SETS).sum()
        ),
        "frozen_target_cv5_90_folds": int(
            cv5_manifest90["analysis_set"].isin(ANALYSIS_SETS).sum()
        ),
        "core_input_rows": {
            "pair_master": len(master),
            "pinc_extension": len(pinc),
            "sequence_clusters": len(clusters),
        },
        "representative_resource_probe": representative_resource,
        "all_input_hash_checks_passed": True,
        "all_sample_and_fold_checks_passed": True,
        "network_access_used": False,
        "dependency_install_or_upgrade_performed": False,
    }
    atomic_json(RUN_DIR / "stage0_preflight.json", summary)
    set_stage_status(STAGE, "complete", details=summary)
    print(json.dumps(_summary_print(summary), ensure_ascii=False))


def _summary_print(summary: dict) -> dict:
    return {
        "stage": summary["stage"],
        "status": summary["status"],
        "analysis_sets": [
            {
                "set": row["analysis_set"],
                "rows": row["rows"],
                "positives": row["positives"],
                "targets": row["targets"],
            }
            for row in summary["analysis_sets"]
        ],
        "representative_fit_seconds": summary["representative_resource_probe"].get(
            "elapsed_seconds"
        ),
    }


if __name__ == "__main__":
    main()
