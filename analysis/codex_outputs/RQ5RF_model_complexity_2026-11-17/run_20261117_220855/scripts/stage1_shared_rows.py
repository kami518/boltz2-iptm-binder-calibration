"""RQ5-RF Stage 1: freeze shared outer-source and inner-task manifests."""
from __future__ import annotations

import json

import pandas as pd
from scipy.special import expit

from rq5rf_common import (
    ANALYSIS_SETS,
    RUN_DIR,
    SOURCE_WEIGHTS,
    analysis_frame,
    atomic_csv,
    atomic_json,
    current_source_weights,
    load_config,
    load_joined_inputs,
    now_iso,
    pair_set_sha256,
    set_stage_status,
    stable_sort,
)

STAGE = "stage_1_shared_rows_and_unit_tests"


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    unit_tests_path = RUN_DIR / "rq5rf_unit_tests.json"
    if not unit_tests_path.exists():
        raise RuntimeError("Unit tests have not produced a result")
    unit_tests = json.loads(unit_tests_path.read_text(encoding="utf-8"))
    if not unit_tests.get("successful", False):
        raise RuntimeError("Unit tests did not pass; formal fitting is prohibited")

    _, _, _, joined = load_joined_inputs(config)
    joined_index = joined.set_index("pair_id")
    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    manifest = pd.read_csv(config["inputs"]["loto_manifest"]["path"], low_memory=False)
    exclusions = pd.read_csv(
        config["inputs"]["cross_target_exclusions"]["path"], low_memory=False
    )
    rule = config["cluster_rule"]
    tasks = []
    inner_rows = []
    for analysis_set in ANALYSIS_SETS:
        outcome_name = config["analysis_sets"][analysis_set]["outcome"]
        valid = set(analysis_frame(config, joined, analysis_set)["pair_id"].astype(str))
        subset_manifest = manifest[
            manifest["analysis_set"].eq(analysis_set)
            & manifest["cluster_rule"].eq(rule)
        ].sort_values("held_out_target")
        for held_target in subset_manifest["held_out_target"].astype(str):
            fold_id = f"loto|{analysis_set}|{rule}|{held_target}"
            fold = loto[loto["fold_id"].eq(fold_id)].copy()
            fold = fold[fold["pair_id"].astype(str).isin(valid)].copy()
            train = fold[fold["split_status"].eq("train")].copy()
            test = fold[fold["split_status"].eq("test")].copy()
            excluded_ids = set(
                exclusions.loc[
                    exclusions["fold_id"].eq(fold_id), "pair_id"
                ].astype(str)
            )
            train = train[~train["pair_id"].astype(str).isin(excluded_ids)].copy()
            train = stable_sort(
                train.merge(
                    joined_index[["x_ipTM", "x_Pinc"]].reset_index(),
                    on="pair_id",
                    how="left",
                    validate="one_to_one",
                )
            )
            test = stable_sort(
                test.merge(
                    joined_index[["x_ipTM", "x_Pinc"]].reset_index(),
                    on="pair_id",
                    how="left",
                    validate="one_to_one",
                )
            )
            if train.empty or test.empty:
                raise RuntimeError(f"Empty outer fold after filtering: {fold_id}")
            if test["target"].astype(str).nunique() != 1:
                raise RuntimeError(f"Outer test must contain one target: {fold_id}")
            source_targets = sorted(train["target"].astype(str).unique())
            source_counts = train["target"].astype(str).value_counts().to_dict()
            if any(count <= 0 for count in source_counts.values()):
                raise RuntimeError(f"Empty outer source target: {fold_id}")
            for weighting in SOURCE_WEIGHTS:
                _, weight_audit = current_source_weights(train, weighting)
                task_id = f"{analysis_set}|{weighting}|{held_target}"
                tasks.append(
                    {
                        "task_id": task_id,
                        "outer_fold_id": fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": rule,
                        "source_weighting": weighting,
                        "held_out_target": held_target,
                        "outer_train_rows": len(train),
                        "outer_train_positives": int(train["outcome_value"].sum()),
                        "outer_test_rows": len(test),
                        "outer_test_positives": int(test["outcome_value"].sum()),
                        "source_n_targets": len(source_targets),
                        "outer_train_pair_set_sha256": pair_set_sha256(
                            train["pair_id"].astype(str)
                        ),
                        "outer_test_pair_set_sha256": pair_set_sha256(
                            test["pair_id"].astype(str)
                        ),
                        "source_target_counts_json": json.dumps(
                            source_counts, sort_keys=True
                        ),
                        "weight_vector_sha256": weight_audit[
                            "weight_vector_sha256"
                        ],
                        "frozen_cross_target_exclusions_present_after_pinc": len(
                            excluded_ids & valid
                        ),
                        "checkpoint_file": f"checkpoints/{analysis_set}__{weighting}__{held_target}.json.gz",
                    }
                )
                for inner_held_target in source_targets:
                    validation = train[
                        train["target"].astype(str).eq(inner_held_target)
                    ]
                    validation_clusters = set(
                        validation["sequence_cluster_value"].astype(str)
                    )
                    candidate = train[
                        train["target"].astype(str).ne(inner_held_target)
                    ]
                    inner_train = candidate[
                        ~candidate["sequence_cluster_value"]
                        .astype(str)
                        .isin(validation_clusters)
                    ]
                    inner_rows.append(
                        {
                            "task_id": task_id,
                            "outer_fold_id": fold_id,
                            "analysis_set": analysis_set,
                            "source_weighting": weighting,
                            "held_out_target": held_target,
                            "inner_held_out_source_target": inner_held_target,
                            "inner_train_candidate_rows": len(candidate),
                            "inner_rows_removed_by_cluster_exclusion": len(candidate)
                            - len(inner_train),
                            "inner_train_rows": len(inner_train),
                            "inner_validation_rows": len(validation),
                            "inner_train_pair_set_sha256": pair_set_sha256(
                                inner_train["pair_id"].astype(str)
                            ),
                            "inner_validation_pair_set_sha256": pair_set_sha256(
                                validation["pair_id"].astype(str)
                            ),
                            "inner_train_target_count": inner_train["target"].nunique(),
                            "expected_inner_train_target_count": len(source_targets)
                            - 1,
                            "cluster_disjoint": not (
                                set(
                                    inner_train["sequence_cluster_value"].astype(str)
                                )
                                & set(
                                    validation["sequence_cluster_value"].astype(str)
                                )
                            ),
                        }
                    )
    task_frame = pd.DataFrame(tasks).sort_values(
        ["analysis_set", "source_weighting", "held_out_target"]
    )
    inner_frame = pd.DataFrame(inner_rows).sort_values(
        [
            "analysis_set",
            "source_weighting",
            "held_out_target",
            "inner_held_out_source_target",
        ]
    )
    if len(task_frame) != 34:
        raise RuntimeError(f"Expected 34 outer source tasks, found {len(task_frame)}")
    if len(inner_frame) != 160:
        raise RuntimeError(f"Expected 160 inner source-target tasks, found {len(inner_frame)}")
    if not inner_frame["cluster_disjoint"].all():
        raise RuntimeError("Inner task manifest contains cluster leakage")
    if not (
        inner_frame["inner_train_target_count"]
        == inner_frame["expected_inner_train_target_count"]
    ).all():
        raise RuntimeError("An inner task excluded an expected source target")
    atomic_csv(RUN_DIR / "shared_outer_task_manifest.csv", task_frame)
    atomic_csv(RUN_DIR / "shared_inner_task_manifest.csv", inner_frame)
    summary = {
        "status": "complete",
        "created_at": now_iso(),
        "unit_tests": unit_tests,
        "outer_source_tasks": len(task_frame),
        "inner_source_target_tasks": len(inner_frame),
        "shared_split_contract": "all four exposures, both weightings, and all lambdas reuse frozen assignments",
        "recoverable_checkpoint_root": str(RUN_DIR / "checkpoints"),
        "all_inner_splits_cluster_disjoint": bool(
            inner_frame["cluster_disjoint"].all()
        ),
    }
    atomic_json(RUN_DIR / "stage1_summary.json", summary)
    set_stage_status(STAGE, "complete", details=summary)
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
