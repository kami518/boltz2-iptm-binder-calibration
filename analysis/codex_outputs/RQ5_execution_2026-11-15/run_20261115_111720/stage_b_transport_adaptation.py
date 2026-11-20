"""RQ5 Stage B: U transport, fixed small-feature combo, and I/F adaptation."""
from __future__ import annotations

import json
import math
import sys

sys.dont_write_bytecode = True

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.model_selection import StratifiedGroupKFold

from rq5_common import (
    ADAPTATION_LAYERS,
    ANALYSIS_SETS,
    CLUSTER_RULE,
    EXPOSURES,
    LAMBDA_GRID,
    OUTCOMES,
    RELIABILITY_EDGES,
    RUN_DIR,
    SOURCE_WEIGHTS,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_json,
    combo_eta,
    fit_combo_model,
    fit_linear_exposure,
    load_config,
    load_core_inputs,
    log_loss,
    performance_point_metrics,
    reliability_frame,
    rq3,
    save_figure,
    set_stage_status,
    source_weights,
    summarize_target_metrics,
)

STAGE = "stage_b_transport_and_adaptation"
U_LAYER = "U_transport_zero_label"
I_LAYER = "predictive_I"
F_LAYER = "predictive_F"
METRICS = (
    "roc_auc",
    "ap",
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "calibration_intercept",
    "calibration_slope",
)


def source_reference(
    train: pd.DataFrame, weighting: str
) -> tuple[float, dict[str, float]]:
    prevalences = (
        train.groupby("target", sort=True)["outcome_value"].mean().astype(float)
    )
    if weighting == "target_equal":
        reference = float(prevalences.mean())
    else:
        reference = float(train["outcome_value"].mean())
    return reference, prevalences.to_dict()


def source_parameter_row(
    *,
    base: dict,
    exposure: str,
    fit: dict,
    train: pd.DataFrame,
    weighting: str,
    selected_lambda: float = np.nan,
) -> dict:
    levels = list(fit.get("levels", sorted(train["target"].unique())))
    if fit.get("estimable", False):
        if exposure == "ipTM_plus_Pinc":
            source_intercepts = np.asarray(fit["theta"][: len(levels)], dtype=float)
            slopes = np.asarray(fit["common_slopes_standardized"], dtype=float)
            intercept_json = json.dumps(
                {level: float(value) for level, value in zip(levels, source_intercepts)},
                sort_keys=True,
            )
            beta_iptm, beta_pinc = (float(slopes[0]), float(slopes[1]))
            common_slope = np.nan
        else:
            source_intercepts = np.asarray(fit["theta"][: len(levels)], dtype=float)
            common_slope = float(fit["theta"][len(levels)])
            beta_iptm = common_slope if exposure == "ipTM_only" else np.nan
            beta_pinc = common_slope if exposure == "Pinc_only" else np.nan
            intercept_json = json.dumps(
                {level: float(value) for level, value in zip(levels, source_intercepts)},
                sort_keys=True,
            )
        means = np.asarray(fit.get("feature_mean", [np.nan, np.nan]), dtype=float)
        sds = np.asarray(fit.get("feature_sd", [np.nan, np.nan]), dtype=float)
        # Linear models do not standardize their raw score.
        if exposure != "ipTM_plus_Pinc":
            means = np.full(2, np.nan)
            sds = np.full(2, np.nan)
        transport_value = float(fit["transport_intercept"])
    else:
        beta_iptm = np.nan
        beta_pinc = np.nan
        common_slope = np.nan
        means = np.full(2, np.nan)
        sds = np.full(2, np.nan)
        transport_value = np.nan
        intercept_json = "{}"
    return {
        **base,
        "exposure": exposure,
        "fit_status": fit.get("status", "failed"),
        "failure_type": fit.get("failure_type", ""),
        "selected_lambda": selected_lambda,
        "n_source_training_rows": int(len(train)),
        "n_source_targets": int(train["target"].nunique()),
        "transport_intercept": transport_value,
        "common_slope": common_slope,
        "beta_ipTM": beta_iptm,
        "beta_Pinc": beta_pinc,
        "feature_mean_ipTM": float(means[0]),
        "feature_sd_ipTM": float(sds[0]),
        "feature_mean_Pinc": float(means[1]),
        "feature_sd_Pinc": float(sds[1]),
        "source_intercepts_json": intercept_json,
        "objective": fit.get("objective", np.nan),
        "gradient_max_abs": fit.get("gradient_max_abs", np.nan),
        "optimizer_success": fit.get("optimizer_success", fit.get("estimable", False)),
    }


def choose_combo_lambda(
    train: pd.DataFrame,
    *,
    analysis_set: str,
    weighting: str,
    fold_id: str,
    held_target: str,
    outcome_name: str,
) -> tuple[dict, list[dict], list[dict]]:
    source_targets = sorted(train["target"].astype(str).unique())
    detailed_rows = []
    candidate_rows = []
    candidate_cache = {}
    for penalty_lambda in LAMBDA_GRID:
        heldout_losses = []
        candidate_failed = False
        for inner_target in source_targets:
            validation = train[train["target"].astype(str).eq(inner_target)].copy()
            inner_train = train[train["target"].astype(str).ne(inner_target)].copy()
            validation_clusters = set(validation["sequence_cluster_value"].astype(str))
            inner_train = inner_train[
                ~inner_train["sequence_cluster_value"].astype(str).isin(validation_clusters)
            ].copy()
            base = {
                "outer_fold_id": fold_id,
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "cluster_rule": CLUSTER_RULE,
                "source_weighting": weighting,
                "held_out_target": held_target,
                "inner_held_out_source_target": inner_target,
                "lambda": penalty_lambda,
            }
            if inner_train.empty or validation.empty:
                detailed_rows.append(
                    {
                        **base,
                        "candidate_status": "failed",
                        "failure_type": "empty_inner_fold",
                        "validation_log_loss": np.nan,
                        "n_inner_train": int(len(inner_train)),
                        "n_inner_validation": int(len(validation)),
                    }
                )
                candidate_failed = True
                continue
            inner_weights = source_weights(inner_train, weighting)
            fit = fit_combo_model(
                inner_train[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float),
                inner_train["outcome_value"].to_numpy(dtype=float),
                inner_train["target"].astype(str).to_numpy(),
                inner_weights,
                weighting,
                penalty_lambda,
                f"{fold_id}|{weighting}|inner_{inner_target}|lambda_{penalty_lambda:g}",
            )
            if fit["estimable"]:
                eta = combo_eta(
                    validation[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float),
                    fit["transport_intercept"],
                    fit["common_slopes_standardized"],
                    fit["feature_mean"],
                    fit["feature_sd"],
                )
                probability = expit(eta)
                validation_loss = log_loss(
                    validation["outcome_value"].to_numpy(dtype=float),
                    probability,
                    float(0.000001),
                )
                heldout_losses.append(validation_loss)
                detailed_rows.append(
                    {
                        **base,
                        "candidate_status": "evaluable",
                        "failure_type": "",
                        "validation_log_loss": validation_loss,
                        "n_inner_train": int(len(inner_train)),
                        "n_inner_validation": int(len(validation)),
                        "inner_feature_mean_ipTM": float(fit["feature_mean"][0]),
                        "inner_feature_sd_ipTM": float(fit["feature_sd"][0]),
                        "inner_feature_mean_Pinc": float(fit["feature_mean"][1]),
                        "inner_feature_sd_Pinc": float(fit["feature_sd"][1]),
                        "inner_transport_intercept": float(fit["transport_intercept"]),
                        "inner_beta_ipTM": float(fit["common_slopes_standardized"][0]),
                        "inner_beta_Pinc": float(fit["common_slopes_standardized"][1]),
                        "inner_objective": float(fit["objective"]),
                    }
                )
            else:
                detailed_rows.append(
                    {
                        **base,
                        "candidate_status": "failed",
                        "failure_type": fit.get("failure_type", "inner_fit_failure"),
                        "validation_log_loss": np.nan,
                        "n_inner_train": int(len(inner_train)),
                        "n_inner_validation": int(len(validation)),
                    }
                )
                candidate_failed = True
        if candidate_failed:
            candidate_rows.append(
                {
                    "outer_fold_id": fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "source_weighting": weighting,
                    "held_out_target": held_target,
                    "lambda": penalty_lambda,
                    "candidate_status": "failed",
                    "failure_type": "one_or_more_inner_folds_failed",
                    "macro_mean_inner_log_loss": np.nan,
                    "n_inner_targets_evaluable": len(heldout_losses),
                    "n_inner_targets_expected": len(source_targets),
                    "selected": False,
                }
            )
        else:
            value = float(np.mean(heldout_losses))
            candidate_cache[penalty_lambda] = value
            candidate_rows.append(
                {
                    "outer_fold_id": fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "source_weighting": weighting,
                    "held_out_target": held_target,
                    "lambda": penalty_lambda,
                    "candidate_status": "evaluable",
                    "failure_type": "",
                    "macro_mean_inner_log_loss": value,
                    "n_inner_targets_evaluable": len(heldout_losses),
                    "n_inner_targets_expected": len(source_targets),
                    "selected": False,
                }
            )
    if not candidate_cache:
        return (
            {
                "status": "failed",
                "estimable": False,
                "failure_type": "all_combo_lambdas_failed",
                "selected_lambda": np.nan,
            },
            candidate_rows,
            detailed_rows,
        )
    best = min(candidate_cache.values())
    tied = [
        penalty_lambda
        for penalty_lambda, value in candidate_cache.items()
        if value - best <= 1e-12
    ]
    selected = float(max(tied))
    for row in candidate_rows:
        row["selected"] = bool(
            row["candidate_status"] == "evaluable" and row["lambda"] == selected
        )
    return (
        {
            "selected_lambda": selected,
            "macro_mean_inner_log_loss": candidate_cache[selected],
            "selection_tie_rule": "1e-12 ties choose larger lambda",
        },
        candidate_rows,
        detailed_rows,
    )


def align_target_metrics(
    left: pd.DataFrame,
    right: pd.DataFrame,
    group_columns: list[str],
) -> pd.DataFrame:
    keys = list(group_columns) + ["held_out_target"]
    columns = keys + [metric for metric in METRICS if metric in left.columns or metric in right.columns]
    return left[columns].merge(
        right[columns],
        on=keys,
        how="inner",
        validate="one_to_one",
        suffixes=("__left", "__right"),
    )


def target_difference_rows(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    grouping: list[str],
    comparison_family: str,
    comparison: str,
    left_model: str,
    right_model: str,
    layer: str,
) -> list[dict]:
    aligned = align_target_metrics(left, right, grouping)
    rows = []
    for _, row in aligned.iterrows():
        base = {column: row[column] for column in grouping + ["held_out_target"]}
        for metric in METRICS:
            left_value = row.get(f"{metric}__left", np.nan)
            right_value = row.get(f"{metric}__right", np.nan)
            rows.append(
                {
                    **base,
                    "comparison_family": comparison_family,
                    "comparison": comparison,
                    "layer": layer,
                    "left_model": left_model,
                    "right_model": right_model,
                    "metric": metric,
                    "left_value": left_value,
                    "right_value": right_value,
                    "difference_left_minus_right": pd.to_numeric(
                        left_value, errors="coerce"
                    )
                    - pd.to_numeric(right_value, errors="coerce"),
                    "brier_logloss_convention": "negative difference means left is better",
                    "auc_ap_bss_convention": "positive difference means left is better",
                }
            )
    return rows


def macro_comparison_rows(
    target_differences: pd.DataFrame,
    grouping: list[str],
) -> list[dict]:
    rows = []
    keys = grouping + [
        "comparison_family",
        "comparison",
        "layer",
        "left_model",
        "right_model",
        "metric",
    ]
    for key_values, group in target_differences.groupby(keys, sort=True):
        if not isinstance(key_values, tuple):
            key_values = (key_values,)
        base = dict(zip(keys, key_values))
        values = pd.to_numeric(
            group["difference_left_minus_right"], errors="coerce"
        ).dropna()
        rows.append(
            {
                **base,
                "summary_type": "target_equal_macro",
                "mean_difference": float(values.mean()) if len(values) else np.nan,
                "target_count": int(len(values)),
                "min_difference": float(values.min()) if len(values) else np.nan,
                "max_difference": float(values.max()) if len(values) else np.nan,
                "positive_difference_targets": int((values > 0).sum()),
                "negative_difference_targets": int((values < 0).sum()),
                "zero_difference_targets": int((values == 0).sum()),
                "direction_convention": (
                    "Brier/log loss: negative means left better; "
                    "AUC/AP/BSS: positive means left better"
                ),
            }
        )
    return rows


def pair_prediction_difference(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    left_probability: str,
    right_probability: str,
    clip: float,
) -> dict:
    common = left[["pair_id"]].merge(
        right[["pair_id"]], on="pair_id", how="inner", validate="one_to_one"
    )["pair_id"]
    lf = left.set_index("pair_id").loc[common]
    rf = right.set_index("pair_id").loc[common]
    y = lf["outcome_value"].to_numpy(dtype=float)
    left_metrics = performance_point_metrics(
        y,
        lf[left_probability].to_numpy(dtype=float),
        clip=clip,
        transport_reference=lf["bss_transport_reference"].to_numpy(dtype=float),
        target_oracle_reference=lf["bss_target_oracle_reference"].to_numpy(dtype=float),
    )
    right_metrics = performance_point_metrics(
        y,
        rf[right_probability].to_numpy(dtype=float),
        clip=clip,
        transport_reference=rf["bss_transport_reference"].to_numpy(dtype=float),
        target_oracle_reference=rf["bss_target_oracle_reference"].to_numpy(dtype=float),
    )
    if "bss_adaptation_train_reference" in lf.columns:
        left_metrics["bss_adaptation_train"] = rq3.bss_from_reference(
            y,
            lf[left_probability].to_numpy(dtype=float),
            lf["bss_adaptation_train_reference"].to_numpy(dtype=float),
        )
        right_metrics["bss_adaptation_train"] = rq3.bss_from_reference(
            y,
            rf[right_probability].to_numpy(dtype=float),
            rf["bss_adaptation_train_reference"].to_numpy(dtype=float),
        )
    result = {"n_pairs": int(len(common))}
    for metric in METRICS:
        left_value = left_metrics.get(metric, np.nan)
        right_value = right_metrics.get(metric, np.nan)
        result[f"left_{metric}"] = left_value
        result[f"right_{metric}"] = right_value
        result[f"difference_{metric}"] = left_value - right_value
    return result


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    clip = float(config["probability_clip"])
    _, _, _, joined = load_core_inputs(config)
    joined_index = joined.set_index("pair_id")
    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    loto_manifest = pd.read_csv(config["inputs"]["loto_manifest"]["path"], low_memory=False)
    cv5 = pd.read_csv(config["inputs"]["target_cv5_assignments"]["path"], low_memory=False)

    u_prediction_frames = []
    u_fold_metric_rows = []
    source_parameter_rows = []
    combo_selection_rows = []
    combo_inner_rows = []
    same_test_row_checks = []

    for analysis_set in ANALYSIS_SETS:
        outcome_name = OUTCOMES[analysis_set]
        valid_pairs = set(
            analysis_frame(config, joined, analysis_set)["pair_id"].astype(str)
        )
        manifest_set = loto_manifest[
            loto_manifest["analysis_set"].eq(analysis_set)
            & loto_manifest["cluster_rule"].eq(CLUSTER_RULE)
        ].sort_values("held_out_target")
        for held_target in manifest_set["held_out_target"].astype(str):
            fold_id = f"loto|{analysis_set}|{CLUSTER_RULE}|{held_target}"
            fold = loto[loto["fold_id"].eq(fold_id)].copy()
            fold = fold[fold["pair_id"].astype(str).isin(valid_pairs)]
            train = fold[fold["split_status"].eq("train")].merge(
                joined_index[
                    ["boltz2_iptm", "pinc", "x_ipTM", "x_Pinc"]
                ].reset_index(),
                on="pair_id",
                how="left",
                validate="one_to_one",
            )
            test = fold[fold["split_status"].eq("test")].merge(
                joined_index[
                    ["boltz2_iptm", "pinc", "x_ipTM", "x_Pinc"]
                ].reset_index(),
                on="pair_id",
                how="left",
                validate="one_to_one",
            )
            train = train.sort_values("pair_id").reset_index(drop=True)
            test = test.sort_values("pair_id").reset_index(drop=True)
            if train.empty or test.empty:
                append_failure(STAGE, fold_id, "empty_fold_after_pinc_filter", "Train or test rows are empty")
                continue
            if train["target"].astype(str).eq(held_target).any():
                append_failure(STAGE, fold_id, "source_target_leak", "Held-out target appears in source training")
                continue
            if set(train["sequence_cluster_value"].astype(str)) & set(
                test["sequence_cluster_value"].astype(str)
            ):
                append_failure(STAGE, fold_id, "cluster_leak", "Source train/test clusters overlap")
                continue
            if train[["x_ipTM", "x_Pinc"]].isna().any().any() or test[
                ["x_ipTM", "x_Pinc"]
            ].isna().any().any():
                append_failure(STAGE, fold_id, "missing_score_after_filter", "A score is missing")
                continue
            source_ref, source_prevalences = source_reference(train, "target_equal")
            pair_source_ref, _ = source_reference(train, "pair")
            target_oracle_ref = float(test["outcome_value"].mean())
            test_ids = set(test["pair_id"].astype(str))

            for weighting in SOURCE_WEIGHTS:
                weights = source_weights(train, weighting)
                base = {
                    "outer_fold_id": fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "held_out_target": held_target,
                    "source_weighting": weighting,
                }
                fits = {}
                for exposure, column in (
                    ("ipTM_only", "x_ipTM"),
                    ("Pinc_only", "x_Pinc"),
                ):
                    fit = fit_linear_exposure(
                        train[column].to_numpy(dtype=float),
                        train["outcome_value"].to_numpy(dtype=float),
                        train["target"].astype(str).to_numpy(),
                        weights,
                        weighting,
                        f"{fold_id}|{weighting}|{exposure}",
                    )
                    fits[exposure] = fit
                    source_parameter_rows.append(
                        source_parameter_row(
                            base=base,
                            exposure=exposure,
                            fit=fit,
                            train=train,
                            weighting=weighting,
                        )
                    )
                    if not fit["estimable"]:
                        append_failure(
                            STAGE,
                            f"{fold_id}|{weighting}|{exposure}",
                            fit.get("failure_type", "source_model_failure"),
                            fit.get("message", "Single-score U fit failed"),
                        )

                selection, candidate_rows, detailed_rows = choose_combo_lambda(
                    train,
                    analysis_set=analysis_set,
                    weighting=weighting,
                    fold_id=fold_id,
                    held_target=held_target,
                    outcome_name=outcome_name,
                )
                combo_selection_rows.extend(candidate_rows)
                combo_inner_rows.extend(detailed_rows)
                selected_lambda = float(selection.get("selected_lambda", np.nan))
                if np.isfinite(selected_lambda):
                    combo_fit = fit_combo_model(
                        train[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float),
                        train["outcome_value"].to_numpy(dtype=float),
                        train["target"].astype(str).to_numpy(),
                        weights,
                        weighting,
                        selected_lambda,
                        f"{fold_id}|{weighting}|combo_final_lambda_{selected_lambda:g}",
                    )
                else:
                    combo_fit = {
                        "status": "failed",
                        "estimable": False,
                        "failure_type": "all_combo_lambdas_failed",
                    }
                fits["ipTM_plus_Pinc"] = combo_fit
                source_parameter_rows.append(
                    source_parameter_row(
                        base=base,
                        exposure="ipTM_plus_Pinc",
                        fit=combo_fit,
                        train=train,
                        weighting=weighting,
                        selected_lambda=selected_lambda,
                    )
                )
                if not combo_fit["estimable"]:
                    append_failure(
                        STAGE,
                        f"{fold_id}|{weighting}|ipTM_plus_Pinc",
                        combo_fit.get("failure_type", "combo_failure"),
                        combo_fit.get("message", "Combination U fit failed"),
                    )

                fold_prediction_frames = []
                for exposure, fit in fits.items():
                    if fit.get("estimable", False):
                        if exposure == "ipTM_plus_Pinc":
                            eta = combo_eta(
                                test[["x_ipTM", "x_Pinc"]].to_numpy(dtype=float),
                                fit["transport_intercept"],
                                fit["common_slopes_standardized"],
                                fit["feature_mean"],
                                fit["feature_sd"],
                            )
                        else:
                            column = "x_ipTM" if exposure == "ipTM_only" else "x_Pinc"
                            eta = (
                                float(fit["transport_intercept"])
                                + float(fit["common_slope"])
                                * test[column].to_numpy(dtype=float)
                            )
                        probability = expit(eta)
                        status = "evaluable"
                        failure_type = ""
                    else:
                        eta = np.full(len(test), np.nan)
                        probability = np.full(len(test), np.nan)
                        status = "failed"
                        failure_type = fit.get("failure_type", "fit_failure")
                    transport_reference_value = (
                        source_ref if weighting == "target_equal" else pair_source_ref
                    )
                    frame = pd.DataFrame(
                        {
                            "pair_id": test["pair_id"].astype(str).to_numpy(),
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "outer_fold_id": fold_id,
                            "cluster_rule": CLUSTER_RULE,
                            "held_out_target": held_target,
                            "target": test["target"].astype(str).to_numpy(),
                            "sequence_cluster_value": test["sequence_cluster_value"].astype(str).to_numpy(),
                            "outcome_value": test["outcome_value"].astype(int).to_numpy(),
                            "source_weighting": weighting,
                            "model_exposure": exposure,
                            "model_layer_label": U_LAYER,
                            "boltz2_iptm": test["boltz2_iptm"].to_numpy(dtype=float),
                            "pinc": test["pinc"].to_numpy(dtype=float),
                            "x_ipTM": test["x_ipTM"].to_numpy(dtype=float),
                            "x_Pinc": test["x_Pinc"].to_numpy(dtype=float),
                            "eta_U": eta,
                            "p_U": probability,
                            "selected_lambda": (
                                selected_lambda
                                if exposure == "ipTM_plus_Pinc"
                                else np.nan
                            ),
                            "fit_status": status,
                            "failure_type": failure_type,
                            "bss_transport_reference": transport_reference_value,
                            "bss_target_oracle_reference": target_oracle_ref,
                            "source_target_prevalences_json": json.dumps(
                                source_prevalences, sort_keys=True
                            ),
                        }
                    )
                    fold_prediction_frames.append(frame)
                    u_prediction_frames.append(frame)
                    metrics = performance_point_metrics(
                        frame["outcome_value"].to_numpy(dtype=float),
                        frame["p_U"].to_numpy(dtype=float),
                        clip=clip,
                        transport_reference=frame["bss_transport_reference"].to_numpy(dtype=float),
                        target_oracle_reference=frame["bss_target_oracle_reference"].to_numpy(dtype=float),
                    )
                    u_fold_metric_rows.append(
                        {
                            **base,
                            "model_exposure": exposure,
                            "model_layer_label": U_LAYER,
                            "selected_lambda": (
                                selected_lambda
                                if exposure == "ipTM_plus_Pinc"
                                else np.nan
                            ),
                            "fit_status": status,
                            "failure_type": failure_type,
                            **metrics,
                        }
                    )
                sets = [set(frame["pair_id"].astype(str)) for frame in fold_prediction_frames]
                same_rows = bool(sets and all(value == test_ids for value in sets))
                same_test_row_checks.append(
                    {
                        "analysis_set": analysis_set,
                        "source_weighting": weighting,
                        "fold_id": fold_id,
                        "n_outer_test_rows": len(test_ids),
                        "same_rows_across_exposures": same_rows,
                    }
                )
                if not same_rows:
                    append_failure(
                        STAGE,
                        f"{fold_id}|{weighting}",
                        "outer_test_row_mismatch",
                        "Exposure predictions do not use identical outer-test rows",
                    )

    u_predictions = pd.concat(u_prediction_frames, ignore_index=True)
    u_fold_metrics = pd.DataFrame(u_fold_metric_rows)
    source_parameters = pd.DataFrame(source_parameter_rows)
    combo_selection = pd.DataFrame(combo_selection_rows)
    combo_inner = pd.DataFrame(combo_inner_rows)
    atomic_csv(RUN_DIR / "rq5_u_pair_predictions.csv", u_predictions)
    atomic_csv(RUN_DIR / "rq5_u_fold_metrics.csv", u_fold_metrics)
    atomic_csv(RUN_DIR / "rq5_source_model_parameters.csv", source_parameters)
    atomic_csv(RUN_DIR / "rq5_combo_lambda_selection.csv", combo_selection)
    atomic_csv(RUN_DIR / "rq5_combo_lambda_inner_metrics.csv", combo_inner)
    atomic_csv(RUN_DIR / "rq5_same_outer_test_rows.csv", pd.DataFrame(same_test_row_checks))

    u_target_rows = []
    u_pair_summary_rows = []
    u_reliability_rows = []
    target_groups = [
        "analysis_set",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "held_out_target",
    ]
    for keys, group in u_predictions.groupby(target_groups, sort=True):
        analysis_set, cluster_rule, weighting, exposure, target = keys
        group = group.sort_values("pair_id").reset_index(drop=True)
        metrics = performance_point_metrics(
            group["outcome_value"].to_numpy(dtype=float),
            group["p_U"].to_numpy(dtype=float),
            clip=clip,
            transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
            target_oracle_reference=group["bss_target_oracle_reference"].to_numpy(dtype=float),
        )
        u_target_rows.append(
            {
                "analysis_set": analysis_set,
                "outcome_name": str(group["outcome_name"].iloc[0]),
                "cluster_rule": cluster_rule,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "held_out_target": target,
                "model_layer_label": U_LAYER,
                "fit_status": str(group["fit_status"].iloc[0]),
                "failure_type": str(group["failure_type"].iloc[0]),
                "bss_adaptation_train": np.nan,
                **metrics,
            }
        )
        for reliability in reliability_frame(
            group["outcome_value"].to_numpy(dtype=float),
            group["p_U"].to_numpy(dtype=float),
            RELIABILITY_EDGES,
        ):
            u_reliability_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": str(group["outcome_name"].iloc[0]),
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "held_out_target": target,
                    "model_layer_label": U_LAYER,
                    **reliability,
                }
            )
    for keys, group in u_predictions.groupby(
        ["analysis_set", "outcome_name", "cluster_rule", "source_weighting", "model_exposure"],
        sort=True,
    ):
        analysis_set, outcome_name, cluster_rule, weighting, exposure = keys
        pooled = performance_point_metrics(
            group["outcome_value"].to_numpy(dtype=float),
            group["p_U"].to_numpy(dtype=float),
            clip=clip,
            transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
            target_oracle_reference=group["bss_target_oracle_reference"].to_numpy(dtype=float),
        )
        u_pair_summary_rows.append(
            {
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "cluster_rule": cluster_rule,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "model_layer_label": U_LAYER,
                "summary_type": "pair_pooled",
                **pooled,
            }
        )
    u_target_metrics = pd.DataFrame(u_target_rows)
    u_macro = summarize_target_metrics(
        u_target_metrics,
        METRICS,
        ["analysis_set", "outcome_name", "cluster_rule", "source_weighting", "model_exposure"],
    )
    atomic_csv(RUN_DIR / "rq5_u_target_metrics.csv", u_target_metrics)
    atomic_csv(RUN_DIR / "rq5_u_macro_summary.csv", u_macro)
    atomic_csv(RUN_DIR / "rq5_u_pair_summary.csv", pd.DataFrame(u_pair_summary_rows))
    atomic_csv(RUN_DIR / "rq5_u_reliability.csv", pd.DataFrame(u_reliability_rows))

    adaptation_predictions = []
    adaptation_fold_metrics = []
    adaptation_parameters = []
    adaptation_inner_assignments = []
    adaptation_lambda_selection = []
    cv_outer_row_checks = []

    for keys, target_u in u_predictions.groupby(
        ["analysis_set", "source_weighting", "model_exposure", "held_out_target"],
        sort=True,
    ):
        analysis_set, weighting, exposure, held_target = keys
        outcome_name = OUTCOMES[analysis_set]
        target_u = target_u.sort_values("pair_id").reset_index(drop=True)
        target_cv = cv5[
            cv5["analysis_set"].eq(analysis_set)
            & cv5["outcome_name"].eq(outcome_name)
            & cv5["cluster_rule"].eq(CLUSTER_RULE)
            & cv5["target"].astype(str).eq(str(held_target))
        ].copy()
        present_ids = set(target_u["pair_id"].astype(str))
        target_cv = target_cv[target_cv["pair_id"].astype(str).isin(present_ids)].copy()
        outer_fold_ids = sorted(target_cv["fold_id"].unique())
        if len(outer_fold_ids) != 5:
            append_failure(
                STAGE,
                f"{analysis_set}|{weighting}|{exposure}|{held_target}",
                "wrong_frozen_outer_cv_count",
                "Expected five frozen target-CV5 folds",
                observed=len(outer_fold_ids),
            )
            continue
        target_index = target_u.set_index("pair_id")
        for outer_fold_id in outer_fold_ids:
            fold_rows = target_cv[target_cv["fold_id"].eq(outer_fold_id)].copy()
            cv_train = fold_rows[fold_rows["split"].eq("train")].sort_values("pair_id")
            cv_test = fold_rows[fold_rows["split"].eq("test")].sort_values("pair_id")
            if cv_train.empty or cv_test.empty:
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{weighting}|{exposure}",
                    "empty_frozen_outer_fold",
                    "Frozen target CV5 train or test is empty after Pinc filtering",
                )
                continue
            if set(cv_train["sequence_cluster_value"].astype(str)) & set(
                cv_test["sequence_cluster_value"].astype(str)
            ):
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{weighting}|{exposure}",
                    "outer_cv_cluster_leak",
                    "Frozen target CV5 train/test clusters overlap",
                )
                continue
            train = cv_train.merge(
                target_u,
                on="pair_id",
                how="left",
                validate="one_to_one",
                suffixes=("", "_u"),
            )
            test = cv_test.merge(
                target_u,
                on="pair_id",
                how="left",
                validate="one_to_one",
                suffixes=("", "_u"),
            )
            if train[["eta_U", "p_U"]].isna().any().any() or test[
                ["eta_U", "p_U"]
            ].isna().any().any():
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{weighting}|{exposure}",
                    "missing_fixed_u",
                    "Fixed U prediction missing from frozen outer CV row",
                )
                continue
            y_train = train["outcome_value"].to_numpy(dtype=float)
            y_test = test["outcome_value"].to_numpy(dtype=float)
            eta_train = train["eta_U"].to_numpy(dtype=float)
            eta_test = test["eta_U"].to_numpy(dtype=float)
            adaptation_prevalence = float(np.mean(y_train))
            adaptation_reference = np.full(len(test), adaptation_prevalence)
            i_fit = rq3.solve_offset_intercept(
                eta_train, y_train, f"{outer_fold_id}|{exposure}|predictive_I"
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
                    f"{outer_fold_id}|{weighting}|{exposure}|predictive_I",
                    i_failure,
                    i_fit.get("message", "I intercept fit failed"),
                )
            adaptation_parameters.append(
                {
                    "outer_fold_id": outer_fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "held_out_target": held_target,
                    "model_layer_label": I_LAYER,
                    "fit_status": i_status,
                    "failure_type": i_failure,
                    "selected_lambda": np.nan,
                    "intercept": i_fit.get("intercept", np.nan),
                    "gamma": 1.0,
                    "delta": 0.0,
                    "n_adaptation_train": int(len(train)),
                    "n_outer_test": int(len(test)),
                    "adaptation_train_prevalence": adaptation_prevalence,
                }
            )

            splitter = StratifiedGroupKFold(
                n_splits=3,
                shuffle=True,
                random_state=int(
                    config["adaptation"]["inner_cv"]["random_state"]
                ),
            )
            inner_splits = list(
                splitter.split(
                    np.zeros((len(train), 1)),
                    y_train.astype(int),
                    groups=train["sequence_cluster_value"].astype(str).to_numpy(),
                )
            )
            if len(inner_splits) != 3:
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{weighting}|{exposure}|inner_cv",
                    "inner_fold_formation_failure",
                    "StratifiedGroupKFold did not form three inner folds",
                )
            for inner_number, (inner_train_index, inner_test_index) in enumerate(
                inner_splits, start=1
            ):
                train_groups = set(
                    train.iloc[inner_train_index]["sequence_cluster_value"].astype(str)
                )
                test_groups = set(
                    train.iloc[inner_test_index]["sequence_cluster_value"].astype(str)
                )
                if train_groups & test_groups:
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{weighting}|{exposure}|inner_{inner_number}",
                        "inner_cluster_leak",
                        "F adaptation inner train/test clusters overlap",
                    )
                inner_train_set = set(inner_train_index.tolist())
                for row_number, row in train.reset_index(drop=True).iterrows():
                    adaptation_inner_assignments.append(
                        {
                            "outer_fold_id": outer_fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": CLUSTER_RULE,
                            "source_weighting": weighting,
                            "model_exposure": exposure,
                            "held_out_target": held_target,
                            "inner_fold_number": inner_number,
                            "pair_id": row["pair_id"],
                            "outcome_value": int(row["outcome_value"]),
                            "sequence_cluster_value": row["sequence_cluster_value"],
                            "split": (
                                "train"
                                if row_number in inner_train_set
                                else "test"
                            ),
                        }
                    )

            lambda_records = []
            for penalty_lambda in LAMBDA_GRID:
                pooled_y = []
                pooled_p = []
                candidate_failure = ""
                objectives = []
                for inner_number, (inner_train_index, inner_test_index) in enumerate(
                    inner_splits, start=1
                ):
                    inner_train = train.iloc[inner_train_index]
                    inner_test = train.iloc[inner_test_index]
                    if inner_train["outcome_value"].nunique() < 2:
                        candidate_failure = f"inner_{inner_number}_single_class_train"
                        break
                    fit = rq3.fit_offset_slope(
                        inner_train["eta_U"].to_numpy(dtype=float),
                        inner_train["outcome_value"].to_numpy(dtype=float),
                        penalty_lambda=penalty_lambda,
                        label=(
                            f"{outer_fold_id}|{weighting}|{exposure}|"
                            f"inner_{inner_number}|lambda_{penalty_lambda:g}"
                        ),
                    )
                    if not fit["estimable"]:
                        candidate_failure = (
                            f"inner_{inner_number}_{fit['failure_type']}"
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
                    objectives.append(float(fit["objective"]))
                if candidate_failure:
                    lambda_records.append(
                        {
                            "outer_fold_id": outer_fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": CLUSTER_RULE,
                            "source_weighting": weighting,
                            "model_exposure": exposure,
                            "held_out_target": held_target,
                            "lambda": penalty_lambda,
                            "candidate_status": "failed",
                            "failure_type": candidate_failure,
                            "pooled_inner_log_loss": np.nan,
                            "n_pooled_inner_predictions": 0,
                            "selected": False,
                        }
                    )
                    continue
                y_pooled = np.concatenate(pooled_y)
                p_pooled = np.concatenate(pooled_p)
                lambda_records.append(
                    {
                        "outer_fold_id": outer_fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": CLUSTER_RULE,
                        "source_weighting": weighting,
                        "model_exposure": exposure,
                        "held_out_target": held_target,
                        "lambda": penalty_lambda,
                        "candidate_status": "evaluable",
                        "failure_type": "",
                        "pooled_inner_log_loss": log_loss(
                            y_pooled, p_pooled, clip
                        ),
                        "n_pooled_inner_predictions": int(len(y_pooled)),
                        "mean_inner_objective": float(np.mean(objectives)),
                        "selected": False,
                    }
                )
            evaluable_lambdas = [
                row
                for row in lambda_records
                if row["candidate_status"] == "evaluable"
            ]
            selected_lambda = np.nan
            if evaluable_lambdas:
                best_loss = min(row["pooled_inner_log_loss"] for row in evaluable_lambdas)
                selected_row = max(
                    [
                        row
                        for row in evaluable_lambdas
                        if row["pooled_inner_log_loss"] - best_loss <= 1e-12
                    ],
                    key=lambda row: row["lambda"],
                )
                selected_lambda = float(selected_row["lambda"])
                for row in lambda_records:
                    row["selected"] = bool(
                        row["candidate_status"] == "evaluable"
                        and row["lambda"] == selected_lambda
                    )
                f_fit = rq3.fit_offset_slope(
                    eta_train,
                    y_train,
                    penalty_lambda=selected_lambda,
                    label=(
                        f"{outer_fold_id}|{weighting}|{exposure}|"
                        f"final_F_lambda_{selected_lambda:g}"
                    ),
                )
            else:
                f_fit = {
                    "status": "failed",
                    "estimable": False,
                    "failure_type": "all_lambdas_failed",
                    "message": "No F lambda candidate was evaluable",
                    "intercept": np.nan,
                    "gamma": np.nan,
                    "delta": np.nan,
                }
            adaptation_lambda_selection.extend(lambda_records)
            if f_fit.get("estimable", False):
                eta_f = float(f_fit["intercept"]) + float(f_fit["gamma"]) * eta_test
                p_f = expit(eta_f)
                f_status = "evaluable"
                f_failure = ""
            else:
                eta_f = np.full(len(test), np.nan)
                p_f = np.full(len(test), np.nan)
                f_status = "failed"
                f_failure = f_fit.get("failure_type", "F_fit_failure")
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{weighting}|{exposure}|predictive_F",
                    f_failure,
                    f_fit.get("message", "F fit failed"),
                )
            adaptation_parameters.append(
                {
                    "outer_fold_id": outer_fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "held_out_target": held_target,
                    "model_layer_label": F_LAYER,
                    "fit_status": f_status,
                    "failure_type": f_failure,
                    "selected_lambda": selected_lambda,
                    "intercept": f_fit.get("intercept", np.nan),
                    "gamma": f_fit.get("gamma", np.nan),
                    "delta": f_fit.get("delta", np.nan),
                    "n_adaptation_train": int(len(train)),
                    "n_outer_test": int(len(test)),
                    "adaptation_train_prevalence": adaptation_prevalence,
                    "objective": f_fit.get("objective", np.nan),
                    "gradient_max_abs": f_fit.get("gradient_max_abs", np.nan),
                }
            )
            base = {
                "pair_id": test["pair_id"].astype(str).to_numpy(),
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "outer_fold_id": outer_fold_id,
                "cluster_rule": CLUSTER_RULE,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "held_out_target": held_target,
                "target": test["target"].astype(str).to_numpy(),
                "sequence_cluster_value": test["sequence_cluster_value"].astype(str).to_numpy(),
                "outcome_value": y_test.astype(int),
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
                "bss_adaptation_train_reference": adaptation_reference,
                "bss_transport_reference": test[
                    "bss_transport_reference"
                ].to_numpy(dtype=float),
                "bss_target_oracle_reference": test[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
                "same_outer_test_rows": True,
            }
            adaptation_predictions.append(pd.DataFrame(base))
            for layer, probability, status in (
                (U_LAYER, test["p_U"].to_numpy(dtype=float), "evaluable"),
                (I_LAYER, p_i, i_status),
                (F_LAYER, p_f, f_status),
            ):
                metrics = performance_point_metrics(
                    y_test,
                    probability,
                    clip=clip,
                    transport_reference=test[
                        "bss_transport_reference"
                    ].to_numpy(dtype=float),
                    target_oracle_reference=test[
                        "bss_target_oracle_reference"
                    ].to_numpy(dtype=float),
                )
                metrics["bss_adaptation_train"] = rq3.bss_from_reference(
                    y_test, probability, adaptation_reference
                )
                adaptation_fold_metrics.append(
                    {
                        "outer_fold_id": outer_fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": CLUSTER_RULE,
                        "source_weighting": weighting,
                        "model_exposure": exposure,
                        "held_out_target": held_target,
                        "model_layer_label": layer,
                        "selected_lambda": (
                            selected_lambda if layer == F_LAYER else np.nan
                        ),
                        "fit_status": status,
                        **metrics,
                    }
                )

    adaptation_prediction_frame = pd.concat(adaptation_predictions, ignore_index=True)
    adaptation_fold_metric_frame = pd.DataFrame(adaptation_fold_metrics)
    adaptation_parameter_frame = pd.DataFrame(adaptation_parameters)
    adaptation_lambda_frame = pd.DataFrame(adaptation_lambda_selection)
    inner_assignment_frame = pd.DataFrame(adaptation_inner_assignments)
    atomic_csv(
        RUN_DIR / "rq5_adaptation_pair_predictions.csv",
        adaptation_prediction_frame,
    )
    atomic_csv(
        RUN_DIR / "rq5_adaptation_fold_metrics.csv",
        adaptation_fold_metric_frame,
    )
    atomic_csv(
        RUN_DIR / "rq5_adaptation_parameters.csv",
        adaptation_parameter_frame,
    )
    atomic_csv(
        RUN_DIR / "rq5_adaptation_lambda_selection.csv",
        adaptation_lambda_frame,
    )
    atomic_csv(
        RUN_DIR / "rq5_adaptation_inner_assignments.csv",
        inner_assignment_frame,
    )

    target_metric_rows = []
    reliability_rows = []
    pair_summary_rows = []
    target_groups = [
        "analysis_set",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "held_out_target",
    ]
    for keys, group in adaptation_prediction_frame.groupby(target_groups, sort=True):
        analysis_set, cluster_rule, weighting, exposure, target = keys
        group = group.sort_values("pair_id").reset_index(drop=True)
        y = group["outcome_value"].to_numpy(dtype=float)
        for layer, probability_column, status_column in (
            (U_LAYER, "p_U", None),
            (I_LAYER, "p_I", "I_status"),
            (F_LAYER, "p_F", "F_status"),
        ):
            metrics = performance_point_metrics(
                y,
                group[probability_column].to_numpy(dtype=float),
                clip=clip,
                transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
                target_oracle_reference=group[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
            )
            metrics["bss_adaptation_train"] = rq3.bss_from_reference(
                y,
                group[probability_column].to_numpy(dtype=float),
                group["bss_adaptation_train_reference"].to_numpy(dtype=float),
            )
            target_metric_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": str(group["outcome_name"].iloc[0]),
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "held_out_target": target,
                    "model_layer_label": layer,
                    "target_oof_fold_count": int(group["outer_fold_id"].nunique()),
                    "fit_status": (
                        "evaluable"
                        if status_column is None
                        else str(group[status_column].iloc[0])
                    ),
                    **metrics,
                }
            )
            for reliability in reliability_frame(
                y,
                group[probability_column].to_numpy(dtype=float),
                RELIABILITY_EDGES,
            ):
                reliability_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "outcome_name": str(group["outcome_name"].iloc[0]),
                        "cluster_rule": cluster_rule,
                        "source_weighting": weighting,
                        "model_exposure": exposure,
                        "held_out_target": target,
                        "model_layer_label": layer,
                        **reliability,
                    }
                )
    target_metrics = pd.DataFrame(target_metric_rows)
    adaptation_reliability = pd.DataFrame(reliability_rows)
    atomic_csv(RUN_DIR / "rq5_adaptation_target_metrics.csv", target_metrics)
    atomic_csv(RUN_DIR / "rq5_adaptation_reliability.csv", adaptation_reliability)
    first_adaptation = adaptation_prediction_frame.drop_duplicates(
        ["analysis_set", "source_weighting", "model_exposure", "held_out_target", "pair_id"]
    )
    for keys, group in first_adaptation.groupby(
        ["analysis_set", "outcome_name", "cluster_rule", "source_weighting", "model_exposure"],
        sort=True,
    ):
        analysis_set, outcome_name, cluster_rule, weighting, exposure = keys
        group = group.set_index("pair_id")
        for layer, probability_column in (
            (U_LAYER, "p_U"),
            (I_LAYER, "p_I"),
            (F_LAYER, "p_F"),
        ):
            metrics = performance_point_metrics(
                group["outcome_value"].to_numpy(dtype=float),
                group[probability_column].to_numpy(dtype=float),
                clip=clip,
                transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
                target_oracle_reference=group[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
            )
            metrics["bss_adaptation_train"] = rq3.bss_from_reference(
                group["outcome_value"].to_numpy(dtype=float),
                group[probability_column].to_numpy(dtype=float),
                group["bss_adaptation_train_reference"].to_numpy(dtype=float),
            )
            pair_summary_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "model_layer_label": layer,
                    "summary_type": "pair_pooled",
                    **metrics,
                }
            )
    adaptation_macro = summarize_target_metrics(
        target_metrics,
        METRICS,
        ["analysis_set", "outcome_name", "cluster_rule", "source_weighting", "model_exposure", "model_layer_label"],
    )
    pair_summary = pd.DataFrame(pair_summary_rows)
    atomic_csv(RUN_DIR / "rq5_adaptation_macro_summary.csv", adaptation_macro)
    atomic_csv(RUN_DIR / "rq5_adaptation_pair_summary.csv", pair_summary)

    # Verify that exposure and layer predictions use one identical row set per target/fold.
    row_set_failures = []
    for keys, group in adaptation_prediction_frame.groupby(
        ["analysis_set", "source_weighting", "held_out_target", "outer_fold_id"]
    ):
        sets = group.groupby(["model_exposure"])["pair_id"].apply(lambda s: frozenset(s.astype(str)))
        if sets.nunique() != 1:
            row_set_failures.append(keys)
    if row_set_failures:
        append_failure(
            STAGE,
            "adaptation_exposure_row_sets",
            "outer_test_row_mismatch",
            "Exposure outer-test row sets differ",
            failures=row_set_failures[:10],
        )

    target_difference_data = []
    pair_comparison_rows = []
    grouping_target = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
    ]
    # Prespecified single-exposure U comparisons.
    exposures = {
        "ipTM_only": u_target_metrics[u_target_metrics["model_exposure"].eq("ipTM_only")],
        "Pinc_only": u_target_metrics[u_target_metrics["model_exposure"].eq("Pinc_only")],
        "ipTM_plus_Pinc": u_target_metrics[u_target_metrics["model_exposure"].eq("ipTM_plus_Pinc")],
    }
    for left_name, right_name in (
        ("Pinc_only", "ipTM_only"),
        ("ipTM_plus_Pinc", "ipTM_only"),
        ("ipTM_plus_Pinc", "Pinc_only"),
    ):
        target_difference_data.extend(
            target_difference_rows(
                exposures[left_name],
                exposures[right_name],
                grouping=grouping_target,
                comparison_family="preset_U_exposure",
                comparison=f"{left_name}_minus_{right_name}",
                left_model=left_name,
                right_model=right_name,
                layer=U_LAYER,
            )
        )
        left_predictions = u_predictions[u_predictions["model_exposure"].eq(left_name)]
        right_predictions = u_predictions[u_predictions["model_exposure"].eq(right_name)]
        for keys, left_group in left_predictions.groupby(grouping_target, sort=True):
            analysis_set, outcome_name, cluster_rule, weighting = keys
            right_group = right_predictions[
                right_predictions["analysis_set"].eq(analysis_set)
                & right_predictions["outcome_name"].eq(outcome_name)
                & right_predictions["cluster_rule"].eq(cluster_rule)
                & right_predictions["source_weighting"].eq(weighting)
            ]
            result = pair_prediction_difference(
                left_group,
                right_group,
                left_probability="p_U",
                right_probability="p_U",
                clip=clip,
            )
            pair_comparison_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "comparison_family": "preset_U_exposure",
                    "comparison": f"{left_name}_minus_{right_name}",
                    "layer": U_LAYER,
                    "left_model": left_name,
                    "right_model": right_name,
                    "summary_type": "pair_pooled",
                    **result,
                }
            )

    # Prespecified I-U, F-U, F-I increments for each exposure.
    for exposure in EXPOSURES:
        subset = target_metrics[target_metrics["model_exposure"].eq(exposure)]
        layers = {
            layer: subset[subset["model_layer_label"].eq(layer)]
            for layer in ADAPTATION_LAYERS
        }
        for left_layer, right_layer in (
            (I_LAYER, U_LAYER),
            (F_LAYER, U_LAYER),
            (F_LAYER, I_LAYER),
        ):
            comparison = f"{left_layer}_minus_{right_layer}"
            target_difference_data.extend(
                target_difference_rows(
                    layers[left_layer],
                    layers[right_layer],
                    grouping=grouping_target,
                    comparison_family="preset_adaptation_increment",
                    comparison=comparison,
                    left_model=exposure,
                    right_model=exposure,
                    layer=f"{left_layer}|{right_layer}",
                )
            )
            for keys, left_group in adaptation_prediction_frame[
                adaptation_prediction_frame["model_exposure"].eq(exposure)
            ].groupby(grouping_target, sort=True):
                analysis_set, outcome_name, cluster_rule, weighting = keys
                right_group = left_group
                probability_left = {
                    I_LAYER: "p_I",
                    F_LAYER: "p_F",
                }[left_layer]
                probability_right = {
                    U_LAYER: "p_U",
                    I_LAYER: "p_I",
                }[right_layer]
                result = pair_prediction_difference(
                    left_group,
                    right_group,
                    left_probability=probability_left,
                    right_probability=probability_right,
                    clip=clip,
                )
                pair_comparison_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": cluster_rule,
                        "source_weighting": weighting,
                        "comparison_family": "preset_adaptation_increment",
                        "comparison": comparison,
                        "layer": f"{left_layer}|{right_layer}",
                        "left_model": exposure,
                        "right_model": exposure,
                        "summary_type": "pair_pooled",
                        **result,
                    }
                )

    target_differences = pd.DataFrame(target_difference_data)
    macro_pairwise = macro_comparison_rows(target_differences, grouping_target)
    pairwise = pd.concat(
        [pd.DataFrame(macro_pairwise), pd.DataFrame(pair_comparison_rows)],
        ignore_index=True,
        sort=False,
    )
    atomic_csv(RUN_DIR / "rq5_target_pairwise_differences.csv", target_differences)
    atomic_csv(RUN_DIR / "rq5_pairwise_comparisons.csv", pairwise)
    sensitivity = pairwise[
        pairwise["comparison_family"].eq("preset_U_exposure")
        | pairwise["comparison_family"].eq("preset_adaptation_increment")
    ].copy()
    sensitivity["analysis_label"] = sensitivity["analysis_set"].map(
        {
            "primary_six": "main_six",
            "primary_five": "five_target_sensitivity",
            "strict_six": "strict_outcome_sensitivity",
        }
    )
    sensitivity.loc[
        sensitivity["source_weighting"].eq("pair"),
        "analysis_label",
    ] = sensitivity.loc[
        sensitivity["source_weighting"].eq("pair"), "analysis_label"
    ] + "__pair_weighting_sensitivity"
    atomic_csv(RUN_DIR / "rq5_sensitivity_summary.csv", sensitivity)

    # Figure 1: U transport diagnostics across exposures, primary_six/target_equal.
    primary_u = u_target_metrics[
        u_target_metrics["analysis_set"].eq("primary_six")
        & u_target_metrics["source_weighting"].eq("target_equal")
    ]
    targets = sorted(primary_u["held_out_target"].unique())
    positions = np.arange(len(targets))
    width = 0.25
    colors = {
        "ipTM_only": "#2468A2",
        "Pinc_only": "#D97941",
        "ipTM_plus_Pinc": "#4D8B63",
    }
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4), constrained_layout=True)
    for axis, metric in zip(axes, ("brier", "bss_transport", "roc_auc")):
        for offset, exposure in zip((-width, 0, width), EXPOSURES):
            values = []
            for target in targets:
                value = primary_u[
                    primary_u["held_out_target"].eq(target)
                    & primary_u["model_exposure"].eq(exposure)
                ][metric]
                values.append(float(value.iloc[0]) if len(value) else np.nan)
            axis.bar(
                positions + offset,
                values,
                width=width,
                color=colors[exposure],
                alpha=0.86,
                label=exposure,
            )
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set_xticks(positions)
        axis.set_xticklabels(targets, rotation=30, ha="right")
        axis.set_title(metric.replace("_", " ").upper())
        axis.grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "RQ5 exploratory U transport by held-out target (primary_six, target_equal; LOTO)",
        fontsize=14,
    )
    save_figure(RUN_DIR / "RQ5_U_transport.png", fig)

    # Figure 2: adaptation increments on identical outer-test rows.
    primary_adapt = target_metrics[
        target_metrics["analysis_set"].eq("primary_six")
        & target_metrics["source_weighting"].eq("target_equal")
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), constrained_layout=True, sharey=True)
    for axis, exposure in zip(axes, EXPOSURES):
        subset = primary_adapt[primary_adapt["model_exposure"].eq(exposure)]
        pivot = subset.pivot(index="held_out_target", columns="model_layer_label", values="brier")
        delta_i = pivot[I_LAYER] - pivot[U_LAYER]
        delta_f = pivot[F_LAYER] - pivot[U_LAYER]
        order = list(pivot.index)
        y = np.arange(len(order))
        axis.scatter(delta_i, y, color="#6B7FD7", label="I - U", s=58)
        axis.scatter(delta_f, y, color="#C8553D", label="F - U", s=58)
        axis.axvline(0, color="black", linewidth=0.8)
        axis.set_yticks(y)
        axis.set_yticklabels(order)
        axis.set_xlabel("Brier difference (negative is better)")
        axis.set_title(exposure)
        axis.grid(axis="x", alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "RQ5 exploratory adaptation increment on identical outer-test rows "
        "(primary_six, target_equal)",
        fontsize=14,
    )
    save_figure(RUN_DIR / "RQ5_adaptation_increment.png", fig)

    summary = {
        "analysis_sets": list(ANALYSIS_SETS),
        "source_weightings": list(SOURCE_WEIGHTS),
        "exposures": list(EXPOSURES),
        "u_prediction_rows": int(len(u_predictions)),
        "adaptation_prediction_rows": int(len(adaptation_prediction_frame)),
        "source_parameter_rows": int(len(source_parameters)),
        "combo_lambda_selection_rows": int(len(combo_selection)),
        "combo_lambda_inner_rows": int(len(combo_inner)),
        "adaptation_parameter_rows": int(len(adaptation_parameter_frame)),
        "adaptation_lambda_selection_rows": int(len(adaptation_lambda_frame)),
        "adaptation_inner_assignment_rows": int(len(inner_assignment_frame)),
        "u_target_metric_rows": int(len(u_target_metrics)),
        "adaptation_target_metric_rows": int(len(target_metrics)),
        "pairwise_summary_rows": int(len(pairwise)),
        "same_outer_test_row_checks": int(len(same_test_row_checks)),
        "same_outer_test_rows_all_passed": all(
            row["same_rows_across_exposures"] for row in same_test_row_checks
        ),
        "adaptation_exposure_row_set_failures": len(row_set_failures),
        "figures": ["RQ5_U_transport.png", "RQ5_adaptation_increment.png"],
    }
    atomic_json(RUN_DIR / "stage_b_summary.json", summary)
    if not summary["same_outer_test_rows_all_passed"] or row_set_failures:
        set_stage_status(STAGE, "failed", **summary)
        raise RuntimeError("Outer-test row equality check failed")
    set_stage_status(STAGE, "complete", **summary)
    print(
        "Stage B complete: "
        f"{len(u_predictions)} U rows, {len(adaptation_prediction_frame)} adaptation rows, "
        f"{len(pairwise)} pairwise summary rows."
    )


if __name__ == "__main__":
    main()
