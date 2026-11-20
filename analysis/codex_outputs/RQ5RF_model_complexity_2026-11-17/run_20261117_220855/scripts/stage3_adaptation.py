"""RQ5-RF Stage 3: target-CV5 I/F adaptation with shared inner assignments."""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from scipy.special import expit

from rq5rf_common import (
    EXPOSURES,
    RUN_DIR,
    SOURCE_WEIGHTS,
    adaptation_inner_splits,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_json,
    canonical_json_sha256,
    fit_offset_slope,
    load_config,
    load_joined_inputs,
    log_loss,
    pair_set_sha256,
    performance_metrics,
    set_stage_status,
    solve_offset_intercept,
    stable_sort,
)

STAGE = "stage_3_adaptation"


def evaluate_layers(
    prediction: dict,
    y_test: np.ndarray,
    *,
    clip: float,
    transport_reference: np.ndarray,
    target_oracle_reference: np.ndarray,
    adaptation_reference: np.ndarray,
) -> list[dict]:
    rows = []
    for layer, probability, status, failure in (
        ("U_transport_zero_label", prediction["p_U"], "evaluable", ""),
        ("predictive_I", prediction["p_I"], prediction["I_status"], prediction["I_failure_type"]),
        ("predictive_F", prediction["p_F"], prediction["F_status"], prediction["F_failure_type"]),
    ):
        metrics = performance_metrics(
            y_test,
            np.asarray(probability, dtype=float),
            clip=clip,
            transport_reference=transport_reference,
            target_oracle_reference=target_oracle_reference,
            adaptation_reference=adaptation_reference,
        )
        rows.append(
            {
                "outer_fold_id": prediction["outer_fold_id"],
                "analysis_set": prediction["analysis_set"],
                "outcome_name": prediction["outcome_name"],
                "target": prediction["target"],
                "source_weighting": prediction["source_weighting"],
                "model_exposure": prediction["model_exposure"],
                "model_layer_label": layer,
                "fit_status": status,
                "failure_type": failure,
                "selected_lambda": prediction["selected_lambda"] if layer == "predictive_F" else np.nan,
                **metrics,
            }
        )
    return rows


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    clip = float(config["probability_clip"])
    u_path = RUN_DIR / "rf_u_pair_predictions.csv"
    if not u_path.exists():
        raise RuntimeError("Stage 2 U predictions are required")
    u_predictions = pd.read_csv(u_path, low_memory=False)
    _, _, _, joined = load_joined_inputs(config)
    cv5 = pd.read_csv(config["inputs"]["target_cv5_assignments"]["path"], low_memory=False)
    cv5 = cv5[cv5["cluster_rule"].eq(config["cluster_rule"])].copy()
    lambda_grid = [float(value) for value in config["adaptation"]["lambda_grid"]]

    pair_prediction_rows = []
    fold_metric_rows = []
    parameter_rows = []
    lambda_rows = []
    inner_assignment_rows = []
    assignment_manifest = []

    for analysis_set in ("primary_six", "primary_five", "strict_six"):
        valid = analysis_frame(config, joined, analysis_set)
        valid_ids = set(valid["pair_id"].astype(str))
        target_u_all = u_predictions[
            u_predictions["analysis_set"].eq(analysis_set)
        ].copy()
        target_u_all = target_u_all[
            target_u_all["pair_id"].astype(str).isin(valid_ids)
        ].copy()
        cv_set = cv5[cv5["analysis_set"].eq(analysis_set)].copy()
        cv_set = cv_set[cv_set["pair_id"].astype(str).isin(valid_ids)].copy()
        for held_target in sorted(cv_set["target"].astype(str).unique()):
            target_cv = cv_set[cv_set["target"].astype(str).eq(held_target)].copy()
            source_outer_fold_id = f"loto|{analysis_set}|{config['cluster_rule']}|{held_target}"
            target_u = target_u_all[
                target_u_all["held_out_target"].astype(str).eq(held_target)
            ].copy()
            if target_u.empty:
                append_failure(
                    STAGE,
                    f"{analysis_set}|{held_target}",
                    "missing_target_u_predictions",
                    "No U predictions for target-CV5 adaptation",
                )
                continue
            expected_target_rows = set(
                valid.loc[valid["target"].astype(str).eq(held_target), "pair_id"].astype(str)
            )
            for cv_outer_fold_id in sorted(target_cv["fold_id"].unique()):
                fold_rows = target_cv[target_cv["fold_id"].eq(cv_outer_fold_id)].copy()
                train_base = fold_rows[fold_rows["split"].eq("train")].copy()
                test_base = fold_rows[fold_rows["split"].eq("test")].copy()
                if train_base.empty or test_base.empty:
                    append_failure(
                        STAGE,
                        cv_outer_fold_id,
                        "empty_frozen_adaptation_fold",
                        "Frozen target-CV5 train or test is empty after Pinc filtering",
                    )
                    continue
                train_base = stable_sort(train_base)
                test_base = stable_sort(test_base)
                if set(train_base["sequence_cluster_value"].astype(str)) & set(
                    test_base["sequence_cluster_value"].astype(str)
                ):
                    append_failure(
                        STAGE,
                        cv_outer_fold_id,
                        "adaptation_cluster_leak",
                        "Frozen target-CV5 train/test clusters overlap",
                    )
                    continue
                if train_base["outcome_value"].nunique() != 2:
                    append_failure(
                        STAGE,
                        cv_outer_fold_id,
                        "adaptation_train_single_class",
                        "Frozen adaptation training has one class",
                    )
                    continue
                try:
                    splits, assignment_hash = adaptation_inner_splits(train_base)
                except Exception as error:
                    append_failure(
                        STAGE,
                        cv_outer_fold_id,
                        "adaptation_inner_split_failure",
                        str(error),
                    )
                    continue
                assignment_manifest.append(
                    {
                        "analysis_set": analysis_set,
                        "held_out_target": held_target,
                        "outer_fold_id": cv_outer_fold_id,
                        "n_adaptation_train": len(train_base),
                        "n_adaptation_test": len(test_base),
                        "inner_assignment_sha256": assignment_hash,
                        "reused_exposures": len(EXPOSURES),
                        "reused_weightings": len(SOURCE_WEIGHTS),
                    }
                )
                for number, (train_index, test_index) in enumerate(splits, start=1):
                    assignment = np.full(len(train_base), "unassigned", dtype=object)
                    assignment[train_index] = "train"
                    assignment[test_index] = "test"
                    for row_number, row in train_base.iterrows():
                        inner_assignment_rows.append(
                            {
                                "analysis_set": analysis_set,
                                "held_out_target": held_target,
                                "outer_fold_id": cv_outer_fold_id,
                                "inner_fold_number": number,
                                "pair_id": str(row["pair_id"]),
                                "outcome_value": int(row["outcome_value"]),
                                "sequence_cluster_value": str(row["sequence_cluster_value"]),
                                "split": str(assignment[row_number]),
                            }
                        )
                ordered_train = stable_sort(train_base)
                ordered_test = stable_sort(test_base)
                y_train = ordered_train["outcome_value"].to_numpy(dtype=float)
                y_test = ordered_test["outcome_value"].to_numpy(dtype=float)
                adaptation_prevalence = float(np.mean(y_train))
                adaptation_reference = np.full(len(ordered_test), adaptation_prevalence)

                for weighting in SOURCE_WEIGHTS:
                    weight_u = target_u[
                        target_u["source_weighting"].eq(weighting)
                    ].copy()
                    for exposure in EXPOSURES:
                        exposure_u = weight_u[
                            weight_u["model_exposure"].eq(exposure)
                        ].copy()
                        train = ordered_train.merge(
                            exposure_u[
                                [
                                    "pair_id",
                                    "outer_fold_id",
                                    "eta_U",
                                    "p_U",
                                    "bss_transport_reference",
                                    "bss_target_oracle_reference",
                                ]
                            ].rename(columns={"outer_fold_id": "source_outer_fold_id"}),
                            on="pair_id",
                            how="left",
                            validate="one_to_one",
                        )
                        test = ordered_test.merge(
                            exposure_u[
                                [
                                    "pair_id",
                                    "outer_fold_id",
                                    "eta_U",
                                    "p_U",
                                    "bss_transport_reference",
                                    "bss_target_oracle_reference",
                                ]
                            ].rename(columns={"outer_fold_id": "source_outer_fold_id"}),
                            on="pair_id",
                            how="left",
                            validate="one_to_one",
                        )
                        needed = ["eta_U", "p_U", "bss_transport_reference", "bss_target_oracle_reference"]
                        if train[needed].isna().any().any() or test[needed].isna().any().any():
                            append_failure(
                                STAGE,
                                f"{cv_outer_fold_id}|{weighting}|{exposure}",
                                "missing_u_after_merge",
                                "U prediction or reference missing after exact pair_id merge",
                            )
                            continue
                        if train["source_outer_fold_id"].astype(str).nunique() != 1 or test["source_outer_fold_id"].astype(str).nunique() != 1:
                            append_failure(
                                STAGE,
                                f"{cv_outer_fold_id}|{weighting}|{exposure}",
                                "source_fold_mapping_failure",
                                "Adaptation rows map to multiple LOTO source folds",
                            )
                            continue
                        eta_train = train["eta_U"].to_numpy(dtype=float)
                        eta_test = test["eta_U"].to_numpy(dtype=float)
                        i_fit = solve_offset_intercept(
                            eta_train,
                            y_train,
                            f"{cv_outer_fold_id}|{weighting}|{exposure}|I",
                        )
                        if i_fit["status"] == "evaluable":
                            eta_i = float(i_fit["intercept"]) + eta_test
                            p_i = expit(eta_i)
                            i_status = "evaluable"
                            i_failure = ""
                        else:
                            eta_i = np.full(len(test), np.nan)
                            p_i = np.full(len(test), np.nan)
                            i_status = "failed"
                            i_failure = i_fit.get("failure_type", "I_fit_failure")
                            append_failure(
                                STAGE,
                                f"{cv_outer_fold_id}|{weighting}|{exposure}|I",
                                i_failure,
                                i_fit.get("message", "I fit failed"),
                            )
                        parameter_rows.append(
                            {
                                "outer_fold_id": cv_outer_fold_id,
                                "source_outer_fold_id": source_outer_fold_id,
                                "analysis_set": analysis_set,
                                "outcome_name": config["analysis_sets"][analysis_set]["outcome"],
                                "target": held_target,
                                "source_weighting": weighting,
                                "model_exposure": exposure,
                                "model_layer_label": "predictive_I",
                                "fit_status": i_status,
                                "failure_type": i_failure,
                                "intercept": i_fit.get("intercept", np.nan),
                                "gamma": 1.0,
                                "delta": 0.0,
                                "selected_lambda": np.nan,
                                "n_adaptation_train": len(train),
                                "n_adaptation_test": len(test),
                                "adaptation_train_positives": int(y_train.sum()),
                                "adaptation_prevalence": adaptation_prevalence,
                                "inner_assignment_sha256": assignment_hash,
                            }
                        )

                        candidates_for_exposure = []
                        for penalty_lambda in lambda_grid:
                            pooled_y = []
                            pooled_p = []
                            inner_objectives = []
                            candidate_failure = ""
                            for number, (inner_train_index, inner_test_index) in enumerate(splits, start=1):
                                inner_train = train.iloc[inner_train_index]
                                inner_test = train.iloc[inner_test_index]
                                if inner_train["outcome_value"].nunique() != 2:
                                    candidate_failure = f"inner_{number}_single_class_train"
                                    break
                                fit = fit_offset_slope(
                                    inner_train["eta_U"].to_numpy(dtype=float),
                                    inner_train["outcome_value"].to_numpy(dtype=float),
                                    penalty_lambda=penalty_lambda,
                                    label=(
                                        f"{cv_outer_fold_id}|{weighting}|{exposure}|"
                                        f"inner_{number}|lambda_{penalty_lambda:g}"
                                    ),
                                )
                                if not fit.get("estimable", False):
                                    candidate_failure = fit.get(
                                        "failure_type", "inner_fit_failure"
                                    )
                                    break
                                eta_inner = (
                                    float(fit["intercept"])
                                    + float(fit["gamma"])
                                    * inner_test["eta_U"].to_numpy(dtype=float)
                                )
                                pooled_y.append(
                                    inner_test["outcome_value"].to_numpy(dtype=float)
                                )
                                pooled_p.append(expit(eta_inner))
                                inner_objectives.append(float(fit["objective"]))
                            if candidate_failure:
                                value = np.nan
                                status = "failed"
                            else:
                                value = log_loss(
                                    np.concatenate(pooled_y),
                                    np.concatenate(pooled_p),
                                    clip,
                                )
                                status = "evaluable"
                            candidates_for_exposure.append(
                                {
                                    "outer_fold_id": cv_outer_fold_id,
                                    "source_outer_fold_id": source_outer_fold_id,
                                    "analysis_set": analysis_set,
                                    "outcome_name": config["analysis_sets"][analysis_set]["outcome"],
                                    "target": held_target,
                                    "source_weighting": weighting,
                                    "model_exposure": exposure,
                                    "lambda": penalty_lambda,
                                    "candidate_status": status,
                                    "failure_type": candidate_failure,
                                    "pooled_inner_log_loss": value,
                                    "n_pooled_inner_predictions": (
                                        int(sum(len(values) for values in pooled_y))
                                        if status == "evaluable"
                                        else 0
                                    ),
                                    "mean_inner_objective": (
                                        float(np.mean(inner_objectives))
                                        if status == "evaluable"
                                        else np.nan
                                    ),
                                    "inner_assignment_sha256": assignment_hash,
                                    "selected": False,
                                }
                            )
                        evaluable = [
                            row
                            for row in candidates_for_exposure
                            if row["candidate_status"] == "evaluable"
                        ]
                        if evaluable:
                            best = min(row["pooled_inner_log_loss"] for row in evaluable)
                            ties = [
                                row
                                for row in evaluable
                                if row["pooled_inner_log_loss"] - best <= 1e-12
                            ]
                            selected_lambda = float(
                                max(row["lambda"] for row in ties)
                            )
                            f_fit = fit_offset_slope(
                                eta_train,
                                y_train,
                                penalty_lambda=selected_lambda,
                                label=(
                                    f"{cv_outer_fold_id}|{weighting}|{exposure}|"
                                    f"final_F_lambda_{selected_lambda:g}"
                                ),
                            )
                            for row in candidates_for_exposure:
                                row["selected"] = bool(
                                    row["candidate_status"] == "evaluable"
                                    and float(row["lambda"]) == selected_lambda
                                )
                        else:
                            selected_lambda = np.nan
                            f_fit = {
                                "status": "failed",
                                "estimable": False,
                                "failure_type": "all_F_lambdas_failed",
                                "intercept": np.nan,
                                "gamma": np.nan,
                                "delta": np.nan,
                            }
                        lambda_rows.extend(candidates_for_exposure)
                        if f_fit.get("estimable", False):
                            eta_f = (
                                float(f_fit["intercept"])
                                + float(f_fit["gamma"]) * eta_test
                            )
                            p_f = expit(eta_f)
                            f_status = "evaluable"
                            f_failure = ""
                        else:
                            eta_f = np.full(len(test), np.nan)
                            p_f = np.full(len(test), np.nan)
                            f_status = "failed"
                            f_failure = f_fit.get(
                                "failure_type", "F_fit_failure"
                            )
                            append_failure(
                                STAGE,
                                f"{cv_outer_fold_id}|{weighting}|{exposure}|F",
                                f_failure,
                                f_fit.get("message", "F fit failed"),
                            )
                        parameter_rows.append(
                            {
                                "outer_fold_id": cv_outer_fold_id,
                                "source_outer_fold_id": source_outer_fold_id,
                                "analysis_set": analysis_set,
                                "outcome_name": config["analysis_sets"][analysis_set]["outcome"],
                                "target": held_target,
                                "source_weighting": weighting,
                                "model_exposure": exposure,
                                "model_layer_label": "predictive_F",
                                "fit_status": f_status,
                                "failure_type": f_failure,
                                "intercept": f_fit.get("intercept", np.nan),
                                "gamma": f_fit.get("gamma", np.nan),
                                "delta": f_fit.get("delta", np.nan),
                                "selected_lambda": selected_lambda,
                                "objective": f_fit.get("objective", np.nan),
                                "gradient_max_abs": f_fit.get("gradient_max_abs", np.nan),
                                "n_adaptation_train": len(train),
                                "n_adaptation_test": len(test),
                                "adaptation_train_positives": int(y_train.sum()),
                                "adaptation_prevalence": adaptation_prevalence,
                                "inner_assignment_sha256": assignment_hash,
                            }
                        )
                        prediction = {
                            "pair_id": test["pair_id"].astype(str).to_numpy(),
                            "analysis_set": analysis_set,
                            "outcome_name": config["analysis_sets"][analysis_set]["outcome"],
                            "outer_fold_id": cv_outer_fold_id,
                            "source_outer_fold_id": source_outer_fold_id,
                            "cluster_rule": config["cluster_rule"],
                            "target": held_target,
                            "sequence_cluster_value": test["sequence_cluster_value"].astype(str).to_numpy(),
                            "outcome_value": y_test.astype(int),
                            "source_weighting": weighting,
                            "model_exposure": exposure,
                            "eta_U": eta_test,
                            "p_U": test["p_U"].to_numpy(dtype=float),
                            "eta_I": eta_i,
                            "p_I": p_i,
                            "eta_F": eta_f,
                            "p_F": p_f,
                            "selected_lambda": selected_lambda,
                            "I_status": i_status,
                            "F_status": f_status,
                            "I_failure_type": i_failure,
                            "F_failure_type": f_failure,
                            "adaptation_train_prevalence": adaptation_prevalence,
                            "bss_adaptation_train_reference": np.full(len(test), adaptation_prevalence),
                            "bss_transport_reference": test["bss_transport_reference"].to_numpy(dtype=float),
                            "bss_target_oracle_reference": test["bss_target_oracle_reference"].to_numpy(dtype=float),
                            "source_fold_mapping_sha256": canonical_json_sha256(
                                {
                                    "source_fold": source_outer_fold_id,
                                    "target": held_target,
                                    "weighting": weighting,
                                }
                            ),
                        }
                        pair_prediction_rows.append(pd.DataFrame(prediction))
                        fold_metric_rows.extend(
                            evaluate_layers(
                                prediction,
                                y_test,
                                clip=clip,
                                transport_reference=test[
                                    "bss_transport_reference"
                                ].to_numpy(dtype=float),
                                target_oracle_reference=test[
                                    "bss_target_oracle_reference"
                                ].to_numpy(dtype=float),
                                adaptation_reference=prediction[
                                    "bss_adaptation_train_reference"
                                ],
                            )
                        )

    predictions = pd.concat(pair_prediction_rows, ignore_index=True)
    fold_metrics = pd.DataFrame(fold_metric_rows)
    parameters = pd.DataFrame(parameter_rows)
    lambdas = pd.DataFrame(lambda_rows)
    assignments = pd.DataFrame(inner_assignment_rows)
    manifests = pd.DataFrame(assignment_manifest)
    for column_set in (
        ["analysis_set", "source_weighting", "target", "pair_id"],
        ["outer_fold_id", "source_weighting", "model_exposure", "pair_id"],
    ):
        predictions.sort_values(column_set, kind="mergesort", inplace=True)
    atomic_csv(RUN_DIR / "rq5rf_pair_predictions.csv", predictions)
    atomic_csv(RUN_DIR / "rq5rf_adaptation_fold_metrics.csv", fold_metrics)
    atomic_csv(RUN_DIR / "rq5rf_adaptation_parameters.csv", parameters)
    atomic_csv(RUN_DIR / "rq5rf_adaptation_lambda_selection.csv", lambdas)
    atomic_csv(RUN_DIR / "rq5rf_adaptation_inner_assignments.csv", assignments)
    atomic_csv(RUN_DIR / "rq5rf_adaptation_inner_manifest.csv", manifests)

    # Hard checks: each target row appears exactly once in concatenated CV tests
    # for each analysis_set × weighting × exposure.
    coverage_failures = []
    for keys, group in predictions.groupby(
        ["analysis_set", "source_weighting", "model_exposure", "target"]
    ):
        duplicates = int(group["pair_id"].astype(str).duplicated().sum())
        if duplicates:
            coverage_failures.append((keys, duplicates, len(group)))
            append_failure(
                STAGE,
                "|".join(str(value) for value in keys),
                "duplicate_target_oof_rows",
                f"Duplicated pair rows in target OOF predictions: {duplicates}",
            )
    same_assignment_hashes = {
        key: group["inner_assignment_sha256"].nunique()
        for key, group in manifests.groupby(["analysis_set", "held_out_target", "outer_fold_id"])
    }
    summary = {
        "status": "complete" if not coverage_failures else "failed",
        "pair_prediction_rows": len(predictions),
        "fold_metric_rows": len(fold_metrics),
        "parameter_rows": len(parameters),
        "lambda_candidate_rows": len(lambdas),
        "inner_assignment_rows": len(assignments),
        "shared_adaptation_outer_folds": len(manifests),
        "coverage_failures": len(coverage_failures),
        "all_each_target_model_oof_unique": not bool(coverage_failures),
        "inner_assignment_hashes_shared": bool(
            len(same_assignment_hashes)
            and all(value == 1 for value in same_assignment_hashes.values())
        ),
        "F_selected_all_evaluable": bool(
            lambdas[lambdas["selected"]]["candidate_status"].eq("evaluable").all()
        ),
        "I_all_evaluable": bool(parameters[parameters["model_layer_label"].eq("predictive_I")]["fit_status"].eq("evaluable").all()),
        "F_all_evaluable": bool(parameters[parameters["model_layer_label"].eq("predictive_F")]["fit_status"].eq("evaluable").all()),
    }
    atomic_json(RUN_DIR / "stage3_adaptation_summary.json", summary)
    set_stage_status(STAGE, summary["status"], details=summary)
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
