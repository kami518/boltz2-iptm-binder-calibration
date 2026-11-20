"""RQ5-RF Stage 2: four pooled U models and complete inner model selection."""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from rq5rf_common import (
    EXPOSURES,
    FEATURE_MAP,
    LOGISTIC_EXPOSURES,
    RF_EXPOSURES,
    RUN_DIR,
    SOURCE_WEIGHTS,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_gzip_json,
    atomic_json,
    canonical_json_sha256,
    current_source_weights,
    effective_rf_grid,
    fit_pooled_logistic,
    fit_random_forest,
    load_config,
    load_joined_inputs,
    log_loss,
    logit_clip,
    logistic_eta,
    minimum_selection,
    pair_set_sha256,
    performance_metrics,
    read_gzip_json,
    rf_positive_probability,
    selection_tie_value,
    set_stage_status,
    source_parameter_row,
    source_reference,
    stable_sort,
    task_seed,
    weight_vector_sha256,
)

STAGE = "stage_2_u_models"
RF_JOBS = 1


def prepare_outer_tasks(config: dict, joined: pd.DataFrame) -> list[dict]:
    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    manifest = pd.read_csv(config["inputs"]["loto_manifest"]["path"], low_memory=False)
    exclusions = pd.read_csv(
        config["inputs"]["cross_target_exclusions"]["path"], low_memory=False
    )
    rule = config["cluster_rule"]
    tasks = []
    for analysis_set in ("primary_six", "primary_five", "strict_six"):
        outcome_name = config["analysis_sets"][analysis_set]["outcome"]
        valid = analysis_frame(config, joined, analysis_set).copy()
        valid_ids = set(valid["pair_id"].astype(str))
        useful = valid[
            [
                "pair_id",
                "target",
                "outcome_value",
                "sequence_cluster_value",
                "boltz2_iptm",
                "pinc",
                "x_ipTM",
                "x_Pinc",
            ]
        ]
        selected_manifest = manifest[
            manifest["analysis_set"].eq(analysis_set)
            & manifest["cluster_rule"].eq(rule)
        ].sort_values("held_out_target")
        for held_target in selected_manifest["held_out_target"].astype(str):
            fold_id = f"loto|{analysis_set}|{rule}|{held_target}"
            fold = loto[loto["fold_id"].eq(fold_id)].copy()
            fold = fold[fold["pair_id"].astype(str).isin(valid_ids)].copy()
            train = fold[fold["split_status"].eq("train")][["pair_id", "split_status"]]
            test = fold[fold["split_status"].eq("test")][["pair_id", "split_status"]]
            train = train.merge(useful, on="pair_id", how="left", validate="one_to_one")
            test = test.merge(useful, on="pair_id", how="left", validate="one_to_one")
            frozen_exclusions = set(
                exclusions.loc[
                    exclusions["fold_id"].eq(fold_id), "pair_id"
                ].astype(str)
            )
            train = train[~train["pair_id"].astype(str).isin(frozen_exclusions)].copy()
            train = stable_sort(train)
            test = stable_sort(test)
            if train.empty or test.empty:
                raise RuntimeError(f"Empty outer fold after Pinc filtering: {fold_id}")
            if test["target"].astype(str).nunique() != 1:
                raise RuntimeError(f"Wrong outer test target set: {fold_id}")
            if set(train["target"].astype(str)).intersection({held_target}):
                raise RuntimeError(f"Held target leaked into source: {fold_id}")
            if set(train["sequence_cluster_value"].astype(str)) & set(
                test["sequence_cluster_value"].astype(str)
            ):
                raise RuntimeError(f"Outer cluster leakage: {fold_id}")
            expected_test_count = int(
                valid["target"].astype(str).eq(held_target).sum()
            )
            if len(test) != expected_test_count:
                raise RuntimeError(
                    f"Outer test does not contain every complete-case target row: {fold_id}"
                )
            target_counts = train["target"].astype(str).value_counts().to_dict()
            if not target_counts or min(target_counts.values()) <= 0:
                raise RuntimeError(f"An expected source target is empty: {fold_id}")
            for weighting in SOURCE_WEIGHTS:
                tasks.append(
                    {
                        "task_id": f"{analysis_set}|{weighting}|{held_target}",
                        "outer_fold_id": fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "source_weighting": weighting,
                        "held_out_target": held_target,
                        "train": train,
                        "test": test,
                        "target_counts": target_counts,
                    }
                )
    if len(tasks) != 34:
        raise RuntimeError(f"Expected 34 outer tasks, built {len(tasks)}")
    return tasks


def task_signature(task: dict, config: dict) -> str:
    return canonical_json_sha256(
        {
            "protocol_version": config["protocol_version"],
            "task_id": task["task_id"],
            "outer_fold_id": task["outer_fold_id"],
            "train_pair_ids": task["train"]["pair_id"].astype(str).tolist(),
            "test_pair_ids": task["test"]["pair_id"].astype(str).tolist(),
            "source_weighting": task["source_weighting"],
            "n_jobs": config["random_forest"]["n_jobs"],
            "n_estimators": config["random_forest"]["n_estimators"],
            "models": EXPOSURES,
        }
    )


def candidate_base(task: dict, inner_target: str, phase: str) -> dict:
    return {
        "outer_fold_id": task["outer_fold_id"],
        "analysis_set": task["analysis_set"],
        "outcome_name": task["outcome_name"],
        "cluster_rule": config_rule(),
        "source_weighting": task["source_weighting"],
        "held_out_target": task["held_out_target"],
        "inner_held_out_source_target": inner_target,
        "fit_phase": phase,
    }


_CONFIG_CACHE: dict | None = None


def config_rule() -> str:
    global _CONFIG_CACHE
    if _CONFIG_CACHE is None:
        _CONFIG_CACHE = load_config()
    return str(_CONFIG_CACHE["cluster_rule"])


def fit_audit(
    *,
    base: dict,
    exposure: str,
    family: str,
    candidate_id: str,
    parameters: dict | None,
    fit: dict,
    train: pd.DataFrame,
    weights: np.ndarray,
    weight_audit: dict,
    selected_lambda: float,
    seed: int | None,
    validation_log_loss: float,
) -> dict:
    return {
        **base,
        **weight_audit,
        "model_exposure": exposure,
        "model_family": family,
        "candidate_id": candidate_id,
        "parameters_json": json.dumps(parameters or {}, sort_keys=True),
        "lambda": selected_lambda,
        "seed": seed if seed is not None else np.nan,
        "training_count": int(len(train)),
        "training_positives": int(train["outcome_value"].sum()),
        "classes_json": json.dumps(fit.get("classes", [])),
        "n_in_bag": fit.get("n_in_bag", np.nan),
        "n_in_bag_min": fit.get("n_in_bag_min", np.nan),
        "n_in_bag_max": fit.get("n_in_bag_max", np.nan),
        "fit_status": fit.get("status", "failed"),
        "failure_type": fit.get("failure_type", ""),
        "validation_log_loss": validation_log_loss,
    }


def aggregate_candidates(
    fit_records: list[dict],
    *,
    expected_inner_targets: int,
) -> tuple[list[dict], dict]:
    candidate_rows = []
    selected_by_exposure = {}
    for (exposure, family, candidate_id), group in pd.DataFrame(fit_records).groupby(
        ["model_exposure", "model_family", "candidate_id"], sort=True
    ):
        losses = pd.to_numeric(group["validation_log_loss"], errors="coerce")
        all_evaluable = bool(
            len(group) == expected_inner_targets
            and group["fit_status"].eq("evaluable").all()
            and losses.notna().all()
            and np.isfinite(losses.to_numpy(dtype=float)).all()
        )
        first = group.iloc[0].to_dict()
        parameter_values = json.loads(first["parameters_json"])
        candidate_rows.append(
            {
                "outer_fold_id": first["outer_fold_id"],
                "analysis_set": first["analysis_set"],
                "outcome_name": first["outcome_name"],
                "cluster_rule": first["cluster_rule"],
                "source_weighting": first["source_weighting"],
                "held_out_target": first["held_out_target"],
                "model_exposure": exposure,
                "model_family": family,
                "candidate_id": candidate_id,
                "parameters_json": first["parameters_json"],
                "max_depth": parameter_values.get("max_depth", np.nan),
                "min_samples_leaf": parameter_values.get("min_samples_leaf", np.nan),
                "max_features": parameter_values.get("max_features", np.nan),
                "depth_rank": {2: 0, 4: 1, None: 2}.get(
                    parameter_values.get("max_depth"), np.nan
                ),
                "lambda": first["lambda"],
                "seed": first["seed"],
                "candidate_status": "evaluable" if all_evaluable else "failed",
                "failure_type": (
                    ""
                    if all_evaluable
                    else "one_or_more_inner_fits_or_validation_not_evaluable"
                ),
                "macro_mean_inner_log_loss": (
                    float(losses.mean()) if all_evaluable else np.nan
                ),
                "inner_targets_evaluable": int(group["fit_status"].eq("evaluable").sum()),
                "inner_targets_expected": int(expected_inner_targets),
                "selected": False,
            }
        )
    candidate_frame = pd.DataFrame(candidate_rows)
    for exposure in EXPOSURES:
        subset = candidate_frame[
            candidate_frame["model_exposure"].eq(exposure)
        ].copy()
        family = str(subset["model_family"].iloc[0])
        records = subset.to_dict("records")
        selected = minimum_selection(records, family)
        selected_by_exposure[exposure] = selected
        if np.isfinite(selected):
            if family == "pooled_logistic":
                mask = subset["lambda"].astype(float).eq(float(selected))
            else:
                mask = subset["candidate_id"].astype(float).eq(float(selected))
            candidate_frame.loc[
                candidate_frame["model_exposure"].eq(exposure) & mask, "selected"
            ] = True
    return candidate_frame.to_dict("records"), selected_by_exposure


def select_record(
    records: list[dict], exposure: str, family: str, selected: float
) -> dict | None:
    if not np.isfinite(selected):
        return None
    for record in records:
        if record["model_exposure"] != exposure:
            continue
        if family == "pooled_logistic" and float(record["lambda"]) == float(selected):
            return record
        if family == "random_forest" and float(record["candidate_id"]) == float(selected):
            return record
    return None


def run_inner_selection(task: dict, config: dict) -> tuple[list[dict], list[dict], dict]:
    train = task["train"]
    source_targets = sorted(train["target"].astype(str).unique())
    fit_records = []
    candidate_aggregate = []
    selections = {}
    lambda_grid = [float(value) for value in config["logistic"]["lambda_grid"]]
    for inner_target in source_targets:
        validation = stable_sort(train[train["target"].astype(str).eq(inner_target)])
        candidate_pool = train[train["target"].astype(str).ne(inner_target)].copy()
        validation_clusters = set(validation["sequence_cluster_value"].astype(str))
        inner_train = candidate_pool[
            ~candidate_pool["sequence_cluster_value"]
            .astype(str)
            .isin(validation_clusters)
        ].copy()
        inner_train = stable_sort(inner_train)
        expected_targets = set(source_targets) - {inner_target}
        present_targets = set(inner_train["target"].astype(str))
        base = candidate_base(task, inner_target, "inner")
        if (
            inner_train.empty
            or len(validation) == 0
            or present_targets != expected_targets
            or validation_clusters & set(inner_train["sequence_cluster_value"].astype(str))
        ):
            raise RuntimeError(
                f"Invalid inner task after exclusions: {task['task_id']} | {inner_target}"
            )
        weights, weight_audit = current_source_weights(
            inner_train, task["source_weighting"]
        )
        outcome = inner_train["outcome_value"].to_numpy(dtype=int)
        validation_outcome = validation["outcome_value"].to_numpy(dtype=float)
        seed = task_seed(
            task["analysis_set"],
            task["source_weighting"],
            task["outer_fold_id"],
            inner_target,
            "inner",
        )

        for exposure in LOGISTIC_EXPOSURES:
            features = FEATURE_MAP[exposure]
            for penalty_lambda in lambda_grid:
                fit = fit_pooled_logistic(
                    inner_train[list(features)].to_numpy(dtype=float),
                    outcome,
                    weights,
                    penalty_lambda,
                    label=f"{task['task_id']}|inner|{inner_target}|{exposure}|{penalty_lambda:g}",
                )
                if fit["estimable"]:
                    eta = logistic_eta(
                        validation[list(features)].to_numpy(dtype=float), fit
                    )
                    value = log_loss(
                        validation_outcome,
                        expit(eta),
                        float(config["probability_clip"]),
                    )
                else:
                    value = float("nan")
                fit_records.append(
                    fit_audit(
                        base=base,
                        exposure=exposure,
                        family="pooled_logistic",
                        candidate_id=f"lambda_{penalty_lambda:g}",
                        parameters=None,
                        fit=fit,
                        train=inner_train,
                        weights=weights,
                        weight_audit=weight_audit,
                        selected_lambda=penalty_lambda,
                        seed=None,
                        validation_log_loss=value,
                    )
                )

        for exposure in RF_EXPOSURES:
            features = FEATURE_MAP[exposure]
            for candidate_id, parameters in enumerate(effective_rf_grid(len(features))):
                fit = fit_random_forest(
                    inner_train[list(features)].to_numpy(dtype=float),
                    outcome,
                    weights,
                    parameters,
                    seed=seed,
                    n_jobs=RF_JOBS,
                    label=f"{task['task_id']}|inner|{inner_target}|{exposure}|candidate_{candidate_id}",
                )
                if fit["estimable"]:
                    probability = rf_positive_probability(
                        fit["estimator"], validation[list(features)].to_numpy(dtype=float)
                    )
                    value = log_loss(
                        validation_outcome,
                        probability,
                        float(config["probability_clip"]),
                    )
                    del fit["estimator"]
                else:
                    value = float("nan")
                fit_records.append(
                    fit_audit(
                        base=base,
                        exposure=exposure,
                        family="random_forest",
                        candidate_id=str(candidate_id),
                        parameters=parameters,
                        fit=fit,
                        train=inner_train,
                        weights=weights,
                        weight_audit=weight_audit,
                        selected_lambda=np.nan,
                        seed=seed,
                        validation_log_loss=value,
                    )
                )
    aggregate, selections = aggregate_candidates(
        fit_records, expected_inner_targets=len(source_targets)
    )
    selected_keys = {
        (row["model_exposure"], row["candidate_id"])
        for row in aggregate
        if bool(row["selected"])
    }
    for record in fit_records:
        record["candidate_selected"] = (
            record["model_exposure"],
            record["candidate_id"],
        ) in selected_keys
    return fit_records, aggregate, selections


def fit_final_models(
    task: dict, config: dict, inner_aggregate: list[dict], selections: dict
) -> tuple[dict, list[dict], dict]:
    train = task["train"]
    weights, weight_audit = current_source_weights(train, task["source_weighting"])
    outcome = train["outcome_value"].to_numpy(dtype=int)
    final_seed = task_seed(
        task["analysis_set"],
        task["source_weighting"],
        task["outer_fold_id"],
        "none",
        "final",
    )
    base = candidate_base(task, "none", "final")
    selected_fits = {}
    source_rows = []
    final_audits = []
    for exposure in LOGISTIC_EXPOSURES:
        family = "pooled_logistic"
        selected = float(selections[exposure])
        selected_record = select_record(inner_aggregate, exposure, family, selected)
        features = FEATURE_MAP[exposure]
        if selected_record is None:
            fit = {
                "status": "failed",
                "estimable": False,
                "failure_type": "no_evaluable_inner_candidate",
            }
            selected = np.nan
        else:
            fit = fit_pooled_logistic(
                train[list(features)].to_numpy(dtype=float),
                outcome,
                weights,
                selected,
                label=f"{task['task_id']}|final|{exposure}|lambda_{selected:g}",
            )
        selected_fits[exposure] = fit
        source_rows.append(
            source_parameter_row(
                base=base,
                exposure=exposure,
                family=family,
                fit=fit,
                train=train,
                weight_audit=weight_audit,
                selected_lambda=selected,
                seed=None,
                parameters=None,
            )
        )
        final_audits.append(
            fit_audit(
                base=base,
                exposure=exposure,
                family=family,
                candidate_id=f"lambda_{selected:g}" if np.isfinite(selected) else "none",
                parameters=None,
                fit=fit,
                train=train,
                weights=weights,
                weight_audit=weight_audit,
                selected_lambda=selected,
                seed=None,
                validation_log_loss=np.nan,
            )
        )
        if not fit.get("estimable", False):
            append_failure(
                STAGE,
                f"{task['task_id']}|{exposure}",
                fit.get("failure_type", "final_fit_failure"),
                "Final pooled logistic fit failed",
            )

    # Fit all final RF candidates for the frozen 3,492-fit audit total. Only
    # the inner-selected candidate is used for prediction.
    for exposure in RF_EXPOSURES:
        features = FEATURE_MAP[exposure]
        selected = float(selections[exposure])
        selected_estimator = None
        selected_fit_record = None
        for candidate_id, parameters in enumerate(effective_rf_grid(len(features))):
            fit = fit_random_forest(
                train[list(features)].to_numpy(dtype=float),
                outcome,
                weights,
                parameters,
                seed=final_seed,
                n_jobs=RF_JOBS,
                label=f"{task['task_id']}|final|{exposure}|candidate_{candidate_id}",
            )
            final_audits.append(
                fit_audit(
                    base=base,
                    exposure=exposure,
                    family="random_forest",
                    candidate_id=str(candidate_id),
                    parameters=parameters,
                    fit=fit,
                    train=train,
                    weights=weights,
                    weight_audit=weight_audit,
                    selected_lambda=np.nan,
                    seed=final_seed,
                    validation_log_loss=np.nan,
                )
            )
            if np.isfinite(selected) and float(candidate_id) == selected:
                selected_estimator = fit.get("estimator")
                selected_fit_record = {
                    key: value for key, value in fit.items() if key != "estimator"
                }
                selected_fit_record["parameters"] = parameters
            elif "estimator" in fit:
                del fit["estimator"]
        if selected_estimator is None or selected_fit_record is None or not selected_fit_record.get("estimable", False):
            selected_fit_record = {
                "status": "failed",
                "estimable": False,
                "failure_type": "no_evaluable_inner_candidate_or_final_refit",
            }
            append_failure(
                STAGE,
                f"{task['task_id']}|{exposure}",
                "final_rf_fit_failure",
                "Selected final RF fit was not available or failed",
            )
        selected_fits[exposure] = selected_fit_record
        selected_fit_record["estimator"] = selected_estimator
        source_rows.append(
            source_parameter_row(
                base=base,
                exposure=exposure,
                family="random_forest",
                fit=selected_fit_record,
                train=train,
                weight_audit=weight_audit,
                selected_lambda=np.nan,
                seed=final_seed,
                parameters=selected_fit_record.get("parameters"),
            )
        )
    return selected_fits, source_rows, final_audits


def run_task(task: dict, config: dict) -> dict:
    outer_id = task["outer_fold_id"]
    task_id = task["task_id"]
    started = time.perf_counter()
    inner_fits, inner_aggregate, selections = run_inner_selection(task, config)
    selected_fits, source_rows, final_audits = fit_final_models(
        task, config, inner_aggregate, selections
    )
    train = task["train"]
    test = task["test"]
    source_ref, source_prevalences = source_reference(
        train, task["source_weighting"]
    )
    target_oracle_reference = float(test["outcome_value"].mean())
    prediction_frames = []
    fold_metrics = []
    for exposure in EXPOSURES:
        fit = selected_fits[exposure]
        features = FEATURE_MAP[exposure]
        if fit.get("estimable", False):
            if exposure in LOGISTIC_EXPOSURES:
                eta = logistic_eta(test[list(features)].to_numpy(dtype=float), fit)
                probability = expit(eta)
            else:
                probability = rf_positive_probability(
                    fit["estimator"], test[list(features)].to_numpy(dtype=float)
                )
                eta = logit_clip(probability, float(config["probability_clip"]))
            status = "evaluable"
            failure_type = ""
        else:
            eta = np.full(len(test), np.nan)
            probability = np.full(len(test), np.nan)
            status = "failed"
            failure_type = fit.get("failure_type", "final_model_failure")
        if not (np.isfinite(eta).all() and np.isfinite(probability).all()):
            status = "failed"
            failure_type = failure_type or "nonfinite_prediction"
            eta = np.full(len(test), np.nan)
            probability = np.full(len(test), np.nan)
        frame = pd.DataFrame(
            {
                "pair_id": test["pair_id"].astype(str).to_numpy(),
                "analysis_set": task["analysis_set"],
                "outcome_name": task["outcome_name"],
                "outer_fold_id": outer_id,
                "cluster_rule": config["cluster_rule"],
                "source_weighting": task["source_weighting"],
                "model_exposure": exposure,
                "held_out_target": task["held_out_target"],
                "target": test["target"].astype(str).to_numpy(),
                "sequence_cluster_value": test["sequence_cluster_value"].astype(str).to_numpy(),
                "outcome_value": test["outcome_value"].astype(int).to_numpy(),
                "boltz2_iptm": test["boltz2_iptm"].to_numpy(dtype=float),
                "pinc": test["pinc"].to_numpy(dtype=float),
                "x_ipTM": test["x_ipTM"].to_numpy(dtype=float),
                "x_Pinc": test["x_Pinc"].to_numpy(dtype=float),
                "eta_U": eta,
                "p_U": probability,
                "fit_status": status,
                "failure_type": failure_type,
                "bss_transport_reference": source_ref,
                "bss_target_oracle_reference": target_oracle_reference,
                "source_target_prevalences_json": json.dumps(
                    source_prevalences, sort_keys=True
                ),
            }
        )
        prediction_frames.append(frame)
        metrics = performance_metrics(
            test["outcome_value"].to_numpy(dtype=float),
            probability,
            clip=float(config["probability_clip"]),
            transport_reference=np.full(len(test), source_ref),
            target_oracle_reference=np.full(len(test), target_oracle_reference),
        )
        fold_metrics.append(
            {
                "outer_fold_id": outer_id,
                "analysis_set": task["analysis_set"],
                "outcome_name": task["outcome_name"],
                "cluster_rule": config["cluster_rule"],
                "source_weighting": task["source_weighting"],
                "model_exposure": exposure,
                "held_out_target": task["held_out_target"],
                "model_layer_label": "U_transport_zero_label",
                "fit_status": status,
                "failure_type": failure_type,
                **metrics,
            }
        )
    predictions = pd.concat(prediction_frames, ignore_index=True)
    expected_ids = set(test["pair_id"].astype(str))
    sets_equal = all(
        set(group["pair_id"].astype(str)) == expected_ids
        for _, group in predictions.groupby("model_exposure")
    )
    if not sets_equal:
        append_failure(
            STAGE,
            f"{task_id}|outer_rows",
            "outer_test_row_mismatch",
            "Four U exposures do not share the complete outer-test row set",
        )
    rf_inner_count = sum(
        row["model_family"] == "random_forest" for row in inner_fits
    )
    payload = {
        "checkpoint_status": "complete",
        "task_id": task_id,
        "task_signature": task_signature(task, config),
        "completed_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "elapsed_seconds": float(time.perf_counter() - started),
        "u_predictions": predictions.to_dict("records"),
        "u_fold_metrics": fold_metrics,
        "source_parameters": source_rows,
        "inner_fit_records": inner_fits,
        "inner_candidate_selection": inner_aggregate,
        "final_candidate_fit_audit": final_audits,
        "selections": selections,
        "same_outer_test_rows": sets_equal,
        "outer_test_pair_set_sha256": pair_set_sha256(expected_ids),
        "rf_inner_fits": int(rf_inner_count),
        "rf_final_candidate_fits": int(
            sum(row["model_family"] == "random_forest" for row in final_audits)
        ),
    }
    return payload


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    if float(config["random_forest"]["n_jobs"]) != float(RF_JOBS):
        raise RuntimeError("Configured RF n_jobs disagrees with serial tested setting")
    _, _, _, joined = load_joined_inputs(config)
    tasks = prepare_outer_tasks(config, joined)
    checkpoint_dir = RUN_DIR / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    payloads = []
    for index, task in enumerate(tasks, start=1):
        safe_name = (
            f"{task['analysis_set']}__{task['source_weighting']}__"
            f"{task['held_out_target']}.json.gz"
        )
        checkpoint_path = checkpoint_dir / safe_name
        signature = task_signature(task, config)
        if checkpoint_path.exists():
            payload = read_gzip_json(checkpoint_path)
            if (
                payload.get("checkpoint_status") == "complete"
                and payload.get("task_signature") == signature
            ):
                payloads.append(payload)
                print(
                    f"[{index:02d}/{len(tasks)}] resume {task['task_id']} "
                    f"({payload['rf_inner_fits'] + payload['rf_final_candidate_fits']} RF fits)",
                    flush=True,
                )
                continue
        print(f"[{index:02d}/{len(tasks)}] start {task['task_id']}", flush=True)
        payload = run_task(task, config)
        atomic_gzip_json(checkpoint_path, payload)
        payloads.append(payload)
        print(
            f"[{index:02d}/{len(tasks)}] complete {task['task_id']} "
            f"in {payload['elapsed_seconds']:.2f}s; "
            f"RF fits={payload['rf_inner_fits'] + payload['rf_final_candidate_fits']}",
            flush=True,
        )

    u_predictions = pd.DataFrame(
        [row for payload in payloads for row in payload["u_predictions"]]
    )
    u_fold_metrics = pd.DataFrame(
        [row for payload in payloads for row in payload["u_fold_metrics"]]
    )
    source_parameters = pd.DataFrame(
        [row for payload in payloads for row in payload["source_parameters"]]
    )
    inner_fit_records = pd.DataFrame(
        [row for payload in payloads for row in payload["inner_fit_records"]]
    )
    inner_candidate_selection = pd.DataFrame(
        [row for payload in payloads for row in payload["inner_candidate_selection"]]
    )
    final_candidate_audit = pd.DataFrame(
        [row for payload in payloads for row in payload["final_candidate_fit_audit"]]
    )
    same_rows = pd.DataFrame(
        [
            {
                "task_id": payload["task_id"],
                "same_outer_test_rows": payload["same_outer_test_rows"],
                "outer_test_pair_set_sha256": payload["outer_test_pair_set_sha256"],
            }
            for payload in payloads
        ]
    )
    for frame in (
        u_predictions,
        u_fold_metrics,
        source_parameters,
        inner_fit_records,
        inner_candidate_selection,
        final_candidate_audit,
    ):
        frame.sort_values(
            [
                column
                for column in (
                    "analysis_set",
                    "source_weighting",
                    "held_out_target",
                    "model_exposure",
                    "pair_id",
                    "inner_held_out_source_target",
                    "candidate_id",
                )
                if column in frame.columns
            ],
            inplace=True,
            kind="mergesort",
        )
    atomic_csv(RUN_DIR / "rf_u_pair_predictions.csv", u_predictions)
    atomic_csv(RUN_DIR / "rf_u_fold_metrics.csv", u_fold_metrics)
    atomic_csv(RUN_DIR / "rf_source_model_parameters.csv", source_parameters)
    atomic_csv(RUN_DIR / "rf_inner_candidate_fits.csv", inner_fit_records)
    atomic_csv(RUN_DIR / "rf_inner_candidate_selection.csv", inner_candidate_selection)
    atomic_csv(RUN_DIR / "rf_final_candidate_fit_audit.csv", final_candidate_audit)
    atomic_csv(RUN_DIR / "same_outer_test_rows.csv", same_rows)

    rf_inner_total = int(
        inner_fit_records["model_family"].eq("random_forest").sum()
    )
    rf_final_total = int(
        final_candidate_audit["model_family"].eq("random_forest").sum()
    )
    summary = {
        "status": "complete",
        "outer_tasks": len(tasks),
        "inner_source_target_tasks": int(
            inner_fit_records[
                [
                    "analysis_set",
                    "source_weighting",
                    "held_out_target",
                    "inner_held_out_source_target",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        ),
        "u_prediction_rows": len(u_predictions),
        "source_parameter_rows": len(source_parameters),
        "inner_candidate_fit_rows": len(inner_fit_records),
        "inner_candidate_selection_rows": len(inner_candidate_selection),
        "final_candidate_fit_rows": len(final_candidate_audit),
        "rf_inner_fits": rf_inner_total,
        "rf_final_candidate_fits": rf_final_total,
        "rf_total_fits": rf_inner_total + rf_final_total,
        "all_same_outer_test_rows": bool(same_rows["same_outer_test_rows"].all()),
        "all_source_fits_evaluable": bool(source_parameters["fit_status"].eq("evaluable").all()),
        "all_inner_selections_evaluable": int(
            inner_candidate_selection["selected"].fillna(False).sum()
        )
        == len(tasks) * len(EXPOSURES),
        "checkpoints_complete": len(payloads),
    }
    if summary["rf_total_fits"] != 3492:
        append_failure(
            STAGE,
            "rf_fit_total",
            "rf_fit_count_mismatch",
            f"Expected 3492 RF fits, observed {summary['rf_total_fits']}",
        )
    atomic_json(RUN_DIR / "stage2_u_models_summary.json", summary)
    set_stage_status(
        STAGE,
        "complete" if summary["all_same_outer_test_rows"] else "failed",
        details=summary,
    )
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
