"""RQ5 stage 0: freeze the run, verify accepted inputs, and audit folds."""
from __future__ import annotations

import shutil
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

from rq5_common import (
    ANALYSIS_SETS,
    CLUSTER_RULE,
    OUTCOMES,
    RQ3_CONFIG_PATH,
    RQ3_IMPLEMENTATION_DIR,
    RQ5_HANDOFF_PATH,
    ROOT,
    RUN_DIR,
    atomic_csv,
    atomic_json,
    append_failure,
    ensure_failure_ledger,
    environment_snapshot,
    load_core_inputs,
    now_iso,
    read_json,
    set_stage_status,
    sha256_file,
)

STAGE = "stage_0_input_and_fold_qc"


def fail(message: str) -> None:
    append_failure(
        STAGE,
        "preflight_acceptance_gate",
        "acceptance_count_mismatch",
        message,
    )
    set_stage_status(STAGE, "blocked", blocker=message)
    raise RuntimeError(message)


def main() -> None:
    set_stage_status(STAGE, "running")
    rq3_config = read_json(RQ3_CONFIG_PATH)
    # Preserve the accepted JSON exactly rather than reserializing it.
    shutil.copyfile(RQ3_CONFIG_PATH, RUN_DIR / "rq3_execution_config.accepted.json")
    shutil.copyfile(RQ3_CONFIG_PATH, RUN_DIR / "run_config.snapshot.json")
    shutil.copyfile(RQ5_HANDOFF_PATH, RUN_DIR / "START_HERE_CODEX.copy.md")
    config = rq3_config
    rq5_config = {
        "handoff_id": "RQ5_EXECUTION_2026-11-15",
        "analysis_scope": "exploratory_incremental_only",
        "created_at": now_iso(),
        "workspace": r"<REPRODUCTION_ROOT>",
        "run_directory": str(RUN_DIR),
        "handoff_path": str(RQ5_HANDOFF_PATH),
        "handoff_sha256": sha256_file(RQ5_HANDOFF_PATH),
        "rq3_config_path": str(RQ3_CONFIG_PATH),
        "rq3_config_sha256": sha256_file(RQ3_CONFIG_PATH),
        "inputs": {
            "pair_master": {
                "path": r"<REPRODUCTION_ROOT>\analysis\v3_outputs\pair_master_v3.csv",
                "sha256": "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
                "expected_rows": 2084,
            },
            "pinc_extension": {
                "path": r"<REPRODUCTION_ROOT>\analysis\pinc_batch\run_20261110_181500\pinc_pair_extension.csv",
                "sha256": "0a549cb890f42948ced456c4a128d480ea06514d4d431d069bcc0f7dc10acaab",
                "expected_rows": 2084,
                "expected_success": 2074,
                "expected_missing_input": 10,
                "missing_target": "pd-l1",
            },
            "sequence_clusters": rq3_config["inputs"]["sequence_clusters"],
            "loto_assignments": rq3_config["inputs"]["loto_assignments"],
            "loto_manifest": rq3_config["inputs"]["loto_manifest"],
            "target_cv5_assignments": rq3_config["inputs"]["target_cv5_assignments"],
            "target_cv5_manifest": rq3_config["inputs"]["target_cv5_manifest"],
            "cross_target_exclusions": rq3_config["inputs"]["cross_target_exclusions"],
            "rq3_implementation": {
                "path": str(RQ3_IMPLEMENTATION_DIR / "rq3_common.py"),
                "sha256": sha256_file(RQ3_IMPLEMENTATION_DIR / "rq3_common.py"),
            },
        },
        "analysis_sets": {
            "primary_six": {"outcome": "primary_outcome", "exclude_target": None},
            "primary_five": {"outcome": "primary_outcome", "exclude_target": "aah-i"},
            "strict_six": {"outcome": "strict_outcome", "exclude_target": None},
        },
        "cluster_rule": CLUSTER_RULE,
        "source_weighting": {"main": "target_equal", "sensitivity": "pair"},
        "probability_clip": 1e-6,
        "reliability_bin_edges": [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0],
        "bootstrap": {"replicates": 2000, "seed": 20261115, "unit": "sequence_cluster90 whole cluster", "paired_scores": True},
        "combo_u": {
            "estimator": "source-target intercepts plus two common ridge-penalized slopes",
            "penalty": "0.5*lambda*(beta_iptm^2+beta_pinc^2); intercepts unpenalized",
            "lambda_grid": [0.001, 0.01, 0.1, 1.0, 10.0, 100.0],
            "inner_selection": "leave-one-source-target-out; equal-target macro mean log loss; 1e-12 ties choose larger lambda",
            "standardization": "weighted mean and population SD from each current source/inner training set only",
        },
        "adaptation": rq3_config["predictive_adaptation"],
        "forbidden": [
            "modify v1/v2/v3, frozen master, clusters, folds, Pinc extension, or prior results",
            "network, install, upgrade, RQ6, random forest, other ML, ipSAE, LIS, extra feature search",
            "change outcomes, features, lambda, weights, folds, or analysis sets based on results",
            "claim population uncertainty over future targets or rewrite RQ3 conclusions",
        ],
    }
    atomic_json(RUN_DIR / "run_config.snapshot.json", rq5_config)
    ensure_failure_ledger()

    expected_hashes = {}
    for key, spec in rq5_config["inputs"].items():
        if isinstance(spec, dict) and "sha256" in spec and "path" in spec:
            expected_hashes[key] = spec
    hash_rows = []
    mismatches = []
    for key, spec in expected_hashes.items():
        path = __import__("pathlib").Path(spec["path"])
        actual = sha256_file(path) if path.is_file() else ""
        matches = actual.lower() == str(spec["sha256"]).lower()
        hash_rows.append(
            {
                "input_key": key,
                "path": str(path),
                "expected_sha256": spec["sha256"],
                "actual_sha256": actual,
                "sha256_matches": matches,
                "bytes": path.stat().st_size if path.is_file() else 0,
                "rows": int(len(pd.read_csv(path, low_memory=False))) if path.suffix.lower() == ".csv" else "",
            }
        )
        if not matches:
            mismatches.append(key)
    atomic_csv(RUN_DIR / "final_input_hashes.csv", pd.DataFrame(hash_rows))
    if mismatches:
        fail(f"Input SHA-256 mismatch: {mismatches}")

    master, pinc, clusters, joined = load_core_inputs(rq5_config)
    statuses = pinc["pinc_status"].value_counts(dropna=False).to_dict()
    if statuses != {"success": 2074, "missing_input": 10}:
        fail(f"Unexpected Pinc status counts: {statuses}")
    missing = joined[joined["pinc_status"].ne("success")]
    if len(missing) != 10 or set(missing["target"].astype(str)) != {"pd-l1"}:
        fail("Pinc missing rows are not exactly 10 PD-L1 rows")
    if int(missing["primary_outcome"].sum()) != 8:
        fail("Pinc missing rows do not contain 8 primary positives")
    if int(missing["strict_outcome"].sum()) != 4:
        fail("Pinc missing rows do not contain 4 strict positives")
    pd_success = joined[
        joined["pinc_status"].eq("success") & joined["target"].astype(str).eq("pd-l1")
    ]
    if len(pd_success) != 267 or int(pd_success["primary_outcome"].sum()) != 190:
        fail("PD-L1 complete-case acceptance count failed")

    set_counts = {}
    for analysis_set in ANALYSIS_SETS:
        outcome = OUTCOMES[analysis_set]
        frame = joined[joined["pinc_status"].eq("success")].copy()
        if analysis_set == "primary_five":
            frame = frame[frame["target"].astype(str).ne("aah-i")]
        n = len(frame)
        positives = int(frame[outcome].sum())
        set_counts[analysis_set] = {
            "n": n,
            "positives": positives,
            "targets": frame.groupby("target").size().to_dict(),
        }
    if set_counts["primary_six"]["n"] != 2074 or set_counts["primary_six"]["positives"] != 740:
        fail("primary_six acceptance count failed")
    if set_counts["primary_five"]["n"] != 1941:
        fail("primary_five acceptance count failed")
    if set_counts["strict_six"]["positives"] != 674:
        fail(
            "strict_six acceptance positive count failed: "
            f"observed {set_counts['strict_six']['positives']} versus corrected handoff-expected 674"
        )

    loto = pd.read_csv(rq5_config["inputs"]["loto_assignments"]["path"], low_memory=False)
    loto_manifest = pd.read_csv(rq5_config["inputs"]["loto_manifest"]["path"], low_memory=False)
    cv5 = pd.read_csv(rq5_config["inputs"]["target_cv5_assignments"]["path"], low_memory=False)
    cv5_manifest = pd.read_csv(rq5_config["inputs"]["target_cv5_manifest"]["path"], low_memory=False)
    exclusions = pd.read_csv(rq5_config["inputs"]["cross_target_exclusions"]["path"], low_memory=False)

    checks = []

    def check(name: str, passed: bool, observed, expected) -> None:
        checks.append({"check": name, "passed": bool(passed), "observed": observed, "expected": expected})
        if not passed:
            fail(f"Stage 0 check failed: {name}; observed={observed}; expected={expected}")

    check("master_rows", len(master) == 2084, len(master), 2084)
    check("master_unique_pair_ids", master["pair_id"].nunique() == 2084, master["pair_id"].nunique(), 2084)
    check("joined_rows", len(joined) == 2084, len(joined), 2084)
    check("loto_manifest_rows", len(loto_manifest) == 44, len(loto_manifest), 44)
    check("target_cv5_manifest_rows", len(cv5_manifest) == 220, len(cv5_manifest), 220)
    check("loto_manifest_qc", loto_manifest["qc_passed"].all(), int(loto_manifest["qc_passed"].sum()), 44)
    check("cv5_manifest_qc", cv5_manifest["qc_passed"].all(), int(cv5_manifest["qc_passed"].sum()), 220)
    check(
        "analysis_rows_match_manifest_primary_six",
        loto[loto["fold_id"].eq("loto|primary_six|sequence_cluster90|aah-i")].shape[0] == 2084,
        loto[loto["fold_id"].eq("loto|primary_six|sequence_cluster90|aah-i")].shape[0],
        2084,
    )
    check(
        "analysis_rows_match_manifest_primary_five",
        loto[loto["fold_id"].eq("loto|primary_five|sequence_cluster90|human-vegf110")].shape[0] == 1951,
        loto[loto["fold_id"].eq("loto|primary_five|sequence_cluster90|human-vegf110")].shape[0],
        1951,
    )
    check(
        "exclusions_accounted",
        len(exclusions) == int(loto["split_status"].eq("removed_cross_target_cluster").sum()),
        len(exclusions),
        int(loto["split_status"].eq("removed_cross_target_cluster").sum()),
    )

    leaks = []
    for row in loto_manifest[
        loto_manifest["analysis_set"].isin(ANALYSIS_SETS)
        & loto_manifest["cluster_rule"].eq(CLUSTER_RULE)
    ].itertuples(index=False):
        fold = loto[loto["fold_id"].eq(row.fold_id)]
        train = fold[fold["split_status"].eq("train")]
        test = fold[fold["split_status"].eq("test")]
        if train["target"].astype(str).eq(str(row.held_out_target)).any():
            leaks.append((row.fold_id, "held_target_in_train"))
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            leaks.append((row.fold_id, "cluster_overlap"))
        if len(test) != int(row.n_test_pairs) or len(train) != int(row.n_train_pairs_final):
            leaks.append((row.fold_id, "fold_count"))
    check("loto_isolation", not leaks, leaks[:5], [])

    cv_leaks = []
    for fold_id, fold in cv5[
        cv5["analysis_set"].isin(ANALYSIS_SETS)
        & cv5["cluster_rule"].eq(CLUSTER_RULE)
    ].groupby("fold_id"):
        train = fold[fold["split"].eq("train")]
        test = fold[fold["split"].eq("test")]
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            cv_leaks.append(fold_id)
    check("cv5_isolation", not cv_leaks, cv_leaks[:5], [])

    qc = {
        "overall_status": "passed",
        "checks": checks,
        "set_counts": set_counts,
        "pinc_status_counts": statuses,
        "pinc_missing": {
            "n": int(len(missing)),
            "targets": missing["target"].astype(str).tolist(),
            "primary_positive": int(missing["primary_outcome"].sum()),
            "strict_positive": int(missing["strict_outcome"].sum()),
            "imputation_performed": False,
        },
        "pd_l1_complete_case": {"n": int(len(pd_success)), "primary_positive": int(pd_success["primary_outcome"].sum())},
        "fold_rows": {
            "loto": int(len(loto)),
            "target_cv5": int(len(cv5)),
            "cross_target_exclusions": int(len(exclusions)),
        },
    }
    atomic_json(RUN_DIR / "stage0_qc.json", qc)
    atomic_json(RUN_DIR / "environment.json", environment_snapshot())
    set_stage_status(
        STAGE,
        "complete",
        checks_passed=len(checks),
        failures=0,
        hash_gate_passed=True,
        config_snapshot_sha256=sha256_file(RUN_DIR / "run_config.snapshot.json"),
    )
    print(f"Stage 0 complete: {len(checks)} checks passed; acceptance counts and frozen folds verified.")


if __name__ == "__main__":
    main()

