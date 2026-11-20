"""RQ3 stage 1: fixed linear and RCS U-layer LOTO transport analyses."""
from __future__ import annotations

import math
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from scipy.special import expit

from rq3_common import (
    ANALYSIS_SETS,
    CLUSTER_RULES,
    LINEAR,
    RCS,
    RUN_DIR,
    SOURCE_WEIGHTS,
    U_LAYER,
    append_failure,
    atomic_csv,
    atomic_json,
    compute_performance_metrics,
    fit_source_linear,
    fit_source_rcs,
    load_config,
    load_joined_data,
    model_configurations,
    now_iso,
    rcs_basis,
    rcs_second_derivative,
    read_csv,
    reliability_rows,
    require_hash_gate,
    set_stage_status,
    validate_fold_structure,
)

STAGE = "stage_1_U_linear_and_RCS_with_predictions_metrics_reliability"
METRIC_COLUMNS = [
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "calibration_intercept",
    "calibration_slope",
    "calibration_intercept_se",
    "calibration_slope_se",
]


def rcs_unit_test() -> dict:
    knots = [-1.0, 0.0, 1.0]
    below = np.linspace(-5.0, -1.0, 101)
    above = np.linspace(1.0, 5.0, 101)
    below_basis = rcs_basis(below, knots)
    above_basis = rcs_basis(above, knots)
    second_below = rcs_second_derivative(below, knots)
    second_above = rcs_second_derivative(above, knots)
    slope_above = np.diff(above_basis[:, 1]) / np.diff(above)
    checks = {
        "b2_zero_below_k1": bool(np.max(np.abs(below_basis[:, 1])) < 1e-14),
        "second_derivative_zero_below_k1": bool(
            np.max(np.abs(second_below)) < 1e-14
        ),
        "second_derivative_zero_above_k3": bool(
            np.max(np.abs(second_above)) < 1e-14
        ),
        "linear_extrapolation_above_k3": bool(
            np.ptp(slope_above) < 1e-12
        ),
    }
    return {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "below_k1_max_abs_B2": float(np.max(np.abs(below_basis[:, 1]))),
        "below_k1_max_abs_second_derivative": float(
            np.max(np.abs(second_below))
        ),
        "above_k3_max_abs_second_derivative": float(
            np.max(np.abs(second_above))
        ),
        "above_k3_slope_range": float(np.ptp(slope_above)),
    }


def parameter_rows(
    *,
    base: dict,
    fit: dict,
    model_form: str,
    source_targets: np.ndarray,
    source_counts: pd.Series,
    source_positive: pd.Series,
    weight_kind: str,
) -> tuple[list[dict], float, float]:
    rows: list[dict] = []
    source_target_count = int(len(source_counts))
    covariance = fit["covariance"]
    for index, level in enumerate(fit["levels"]):
        rows.append(
            {
                **base,
                "model_form": model_form,
                "parameter_name": "source_intercept",
                "source_target": level,
                "estimate": (
                    float(fit["theta"][index])
                    if fit["estimable"]
                    else np.nan
                ),
                "standard_error": (
                    float(math.sqrt(max(covariance[index, index], 0.0)))
                    if fit["estimable"] and np.isfinite(covariance[index, index])
                    else np.nan
                ),
                "n_pairs_source_target": int(source_counts.loc[level]),
                "n_positive_source_target": int(source_positive.loc[level]),
                "fit_status": fit["status"],
                "failure_type": fit["failure_type"] or "",
                "objective": fit.get("objective", np.nan),
                "iterations": fit.get("iterations", np.nan),
                "gradient_max_abs": fit.get("gradient_max_abs", np.nan),
                "information_qc_message": fit.get("information_qc_message", ""),
            }
        )
    coefficient_names = (
        ["common_slope_x"]
        if model_form == LINEAR
        else ["common_B1", "common_B2"]
    )
    coefficient_indices = list(
        range(source_target_count, source_target_count + len(coefficient_names))
    )
    for name, index in zip(coefficient_names, coefficient_indices):
        rows.append(
            {
                **base,
                "model_form": model_form,
                "parameter_name": name,
                "source_target": "",
                "estimate": (
                    float(fit["theta"][index])
                    if fit["estimable"]
                    else np.nan
                ),
                "standard_error": (
                    float(math.sqrt(max(covariance[index, index], 0.0)))
                    if fit["estimable"] and np.isfinite(covariance[index, index])
                    else np.nan
                ),
                "n_pairs_source_target": source_target_count,
                "n_positive_source_target": int(sum(source_positive)),
                "fit_status": fit["status"],
                "failure_type": fit["failure_type"] or "",
                "objective": fit.get("objective", np.nan),
                "iterations": fit.get("iterations", np.nan),
                "gradient_max_abs": fit.get("gradient_max_abs", np.nan),
                "information_qc_message": fit.get("information_qc_message", ""),
            }
        )
    if weight_kind == "target_equal":
        transport_weights = np.full(source_target_count, 1.0 / source_target_count)
    else:
        ordered_counts = source_counts.loc[fit["levels"]].to_numpy(dtype=float)
        transport_weights = ordered_counts / ordered_counts.sum()
    transport_intercept = (
        float(np.dot(transport_weights, fit["theta"][:source_target_count]))
        if fit["estimable"]
        else np.nan
    )
    if fit["estimable"]:
        coefficient = np.r_[transport_weights, np.zeros(len(coefficient_names))]
        transport_se = float(
            math.sqrt(max(float(coefficient @ covariance @ coefficient), 0.0))
        )
    else:
        transport_se = np.nan
    rows.append(
        {
            **base,
            "model_form": model_form,
            "parameter_name": "transport_intercept",
            "source_target": "",
            "estimate": transport_intercept,
            "standard_error": transport_se,
            "n_pairs_source_target": int(source_counts.sum()),
            "n_positive_source_target": int(source_positive.sum()),
            "fit_status": fit["status"],
            "failure_type": fit["failure_type"] or "",
            "objective": fit.get("objective", np.nan),
            "iterations": fit.get("iterations", np.nan),
            "gradient_max_abs": fit.get("gradient_max_abs", np.nan),
            "information_qc_message": fit.get("information_qc_message", ""),
        }
    )
    return rows, transport_intercept, transport_se


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    require_hash_gate(config)
    joined = load_joined_data(config)
    loto = read_csv(
        __import__("pathlib").Path(config["inputs"]["loto_assignments"]["path"])
    )
    manifest = read_csv(
        __import__("pathlib").Path(config["inputs"]["loto_manifest"]["path"])
    )
    validate_fold_structure(config, joined)[0:0]  # preserve the stage-0 gate

    rcs_test = rcs_unit_test()
    atomic_json(RUN_DIR / "rcs_unit_test.json", rcs_test)
    if rcs_test["status"] != "passed":
        append_failure(
            STAGE,
            "rcs_unit_test",
            "unit_test_failure",
            "Restricted cubic spline formula unit tests failed.",
            checks=rcs_test["checks"],
        )
        raise RuntimeError("RCS unit test failed")

    joined_index = joined.set_index("pair_id")
    prediction_frames: list[pd.DataFrame] = []
    parameter_frames: list[pd.DataFrame] = []
    metric_frames: list[pd.DataFrame] = []
    reliability_frames: list[pd.DataFrame] = []
    knot_frames: list[pd.DataFrame] = []
    expected_evaluations = 0

    for analysis_set in ANALYSIS_SETS:
        for cluster_rule in CLUSTER_RULES:
            manifest_set = manifest[
                manifest["analysis_set"].eq(analysis_set)
                & manifest["cluster_rule"].eq(cluster_rule)
            ].sort_values("held_out_target")
            for held_target in manifest_set["held_out_target"].astype(str):
                fold_id = (
                    f"loto|{analysis_set}|{cluster_rule}|{held_target}"
                )
                fold = loto[loto["fold_id"].eq(fold_id)]
                train = fold[fold["split_status"].eq("train")].copy()
                test = fold[fold["split_status"].eq("test")].copy()
                if train.empty or test.empty:
                    append_failure(
                        STAGE,
                        fold_id,
                        "missing_fold_rows",
                        "Frozen LOTO fold has no training or test rows.",
                    )
                    continue
                if train["target"].astype(str).eq(held_target).any():
                    append_failure(
                        STAGE,
                        fold_id,
                        "held_out_target_leak",
                        "Held-out target appears in U training rows.",
                    )
                    continue
                if set(train["sequence_cluster_value"].astype(str)) & set(
                    test["sequence_cluster_value"].astype(str)
                ):
                    append_failure(
                        STAGE,
                        fold_id,
                        "cluster_leak",
                        "Training and test sequence clusters overlap.",
                    )
                    continue
                x_train = train["sequence_cluster_value"].notna()
                if not x_train.all():
                    append_failure(
                        STAGE,
                        fold_id,
                        "missing_cluster",
                        "A frozen cluster value is missing.",
                    )
                    continue
                train_joined = train.join(
                    joined_index[["boltz2_iptm", "logit_iptm"]],
                    on="pair_id",
                    how="left",
                    validate="many_to_one",
                )
                test_joined = test.join(
                    joined_index[["boltz2_iptm", "logit_iptm"]],
                    on="pair_id",
                    how="left",
                    validate="many_to_one",
                )
                if (
                    train_joined[["boltz2_iptm", "logit_iptm"]].isna().any().any()
                    or test_joined[["boltz2_iptm", "logit_iptm"]]
                    .isna()
                    .any()
                    .any()
                ):
                    append_failure(
                        STAGE,
                        fold_id,
                        "missing_iptm",
                        "A LOTO fold contains missing ipTM values.",
                    )
                    continue

                x = train_joined["logit_iptm"].to_numpy(dtype=float)
                y = train_joined["outcome_value"].to_numpy(dtype=float)
                training_targets = train_joined["target"].astype(str).to_numpy()
                test_x = test_joined["logit_iptm"].to_numpy(dtype=float)
                source_counts = (
                    train_joined.groupby("target", sort=True).size().astype(int)
                )
                source_positive = (
                    train_joined.groupby("target", sort=True)["outcome_value"]
                    .sum()
                    .astype(int)
                )

                for weight_kind in SOURCE_WEIGHTS:
                    if weight_kind == "target_equal":
                        weights = train_joined[
                            "weight_target_equal"
                        ].to_numpy(dtype=float)
                        knot_weights = train_joined[
                            "weight_target_equal_raw"
                        ].to_numpy(dtype=float)
                    else:
                        weights = train_joined["weight_pair"].to_numpy(dtype=float)
                        knot_weights = weights.copy()
                    if (
                        not np.isfinite(weights).all()
                        or not np.isfinite(knot_weights).all()
                        or (weights <= 0).any()
                        or (knot_weights <= 0).any()
                    ):
                        append_failure(
                            STAGE,
                            f"{fold_id}|{weight_kind}",
                            "invalid_training_weight",
                            "Training weights are nonfinite, missing, or nonpositive.",
                        )
                        continue
                    base = {
                        "fold_id": fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": test_joined["outcome_name"].iloc[0],
                        "cluster_rule": cluster_rule,
                        "held_out_target": held_target,
                        "source_weighting": weight_kind,
                        "n_training_pairs": int(len(train_joined)),
                        "n_source_targets": int(len(source_counts)),
                    }
                    for model_form in (LINEAR, RCS):
                        expected_evaluations += 1
                        if model_form == LINEAR:
                            fit = fit_source_linear(
                                x,
                                y,
                                training_targets,
                                weights,
                                f"{fold_id}|{weight_kind}|{model_form}",
                            )
                        else:
                            fit = fit_source_rcs(
                                x,
                                y,
                                training_targets,
                                knot_weights,
                                config["restricted_cubic_spline"]["knot_quantiles"],
                                f"{fold_id}|{weight_kind}|{model_form}",
                            )
                        if not fit["estimable"]:
                            append_failure(
                                STAGE,
                                f"{fold_id}|{weight_kind}|{model_form}",
                                fit.get("failure_type") or "source_model_failure",
                                fit.get("message", "Source model failed."),
                                objective=fit.get("objective"),
                            )
                        rows, transport_intercept, transport_se = parameter_rows(
                            base=base,
                            fit=fit,
                            model_form=model_form,
                            source_targets=training_targets,
                            source_counts=source_counts,
                            source_positive=source_positive,
                            weight_kind=weight_kind,
                        )
                        parameter_frames.append(pd.DataFrame(rows))

                        if model_form == RCS:
                            knots = fit.get("knots", np.full(3, np.nan))
                            knot_frames.append(
                                pd.DataFrame(
                                    [
                                        {
                                            **base,
                                            "quantile": quantile,
                                            "knot_x": float(knot),
                                            "knot_status": (
                                                "evaluable"
                                                if np.isfinite(knot)
                                                else "failed"
                                            ),
                                            "knot_source": "current LOTO source training rows only",
                                        }
                                        for quantile, knot in zip(
                                            config["restricted_cubic_spline"][
                                                "knot_quantiles"
                                            ],
                                            knots,
                                        )
                                    ]
                                )
                            )

                        if fit["estimable"]:
                            source_target_count = len(fit["levels"])
                            if model_form == LINEAR:
                                eta = transport_intercept + fit["theta"][
                                    source_target_count
                                ] * test_x
                            else:
                                basis = rcs_basis(test_x, fit["knots"])
                                eta = (
                                    transport_intercept
                                    + basis
                                    @ fit["theta"][
                                        source_target_count : source_target_count + 2
                                    ]
                                )
                            probability = expit(eta)
                        else:
                            eta = np.full(len(test_joined), np.nan)
                            probability = np.full(len(test_joined), np.nan)

                        if weight_kind == "target_equal":
                            source_target_prevalences = (
                                source_positive.astype(float)
                                / source_counts.astype(float)
                            )
                            transport_reference_value = float(
                                source_target_prevalences.mean()
                            )
                        else:
                            transport_reference_value = float(
                                source_positive.sum() / source_counts.sum()
                            )
                        target_oracle_reference_value = float(
                            test_joined["outcome_value"].mean()
                        )
                        predictions = pd.DataFrame(
                            {
                                "pair_id": test_joined["pair_id"].astype(str).to_numpy(),
                                "analysis_set": analysis_set,
                                "outcome_name": test_joined[
                                    "outcome_name"
                                ].astype(str).to_numpy(),
                                "fold_id": fold_id,
                                "cluster_rule": cluster_rule,
                                "held_out_target": held_target,
                                "target": test_joined["target"].astype(str).to_numpy(),
                                "sequence_cluster_value": test_joined[
                                    "sequence_cluster_value"
                                ].astype(str).to_numpy(),
                                "outcome_value": test_joined[
                                    "outcome_value"
                                ].astype(int).to_numpy(),
                                "boltz2_iptm": test_joined[
                                    "boltz2_iptm"
                                ].to_numpy(dtype=float),
                                "x": test_x,
                                "source_weighting": weight_kind,
                                "model_form": model_form,
                                "model_layer_label": U_LAYER,
                                "transport_intercept": transport_intercept,
                                "transport_intercept_se": transport_se,
                                "eta_U": eta,
                                "predicted_probability": probability,
                                "bss_transport_reference": transport_reference_value,
                                "bss_target_oracle_reference": (
                                    target_oracle_reference_value
                                ),
                                "fit_status": fit["status"],
                                "failure_type": fit["failure_type"] or "",
                            }
                        )
                        prediction_frames.append(predictions)

                        metric_base = {
                            "fold_id": fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": str(test_joined["outcome_name"].iloc[0]),
                            "cluster_rule": cluster_rule,
                            "held_out_target": held_target,
                            "source_weighting": weight_kind,
                            "model_form": model_form,
                            "model_layer_label": U_LAYER,
                            "n": int(len(predictions)),
                            "positives": int(test_joined["outcome_value"].sum()),
                            "convergence_flags": fit["status"],
                            "fit_message": fit.get("message", ""),
                            "transport_intercept": transport_intercept,
                            "transport_intercept_se": transport_se,
                            "bss_transport_reference_value": (
                                transport_reference_value
                            ),
                            "bss_target_oracle_reference_value": (
                                target_oracle_reference_value
                            ),
                        }
                        if fit["estimable"]:
                            metrics = compute_performance_metrics(
                                test_joined["outcome_value"].to_numpy(dtype=float),
                                probability,
                                clip=float(config["probability_clip"]),
                                transport_reference=predictions[
                                    "bss_transport_reference"
                                ].to_numpy(dtype=float),
                                target_oracle_reference=predictions[
                                    "bss_target_oracle_reference"
                                ].to_numpy(dtype=float),
                            )
                            if metrics["calibration_status"] != "evaluable":
                                append_failure(
                                    STAGE,
                                    f"{fold_id}|{weight_kind}|{model_form}|calibration",
                                    "calibration_failure",
                                    metrics["calibration_message"],
                                )
                        else:
                            metrics = {
                                "n": int(len(predictions)),
                                "positives": int(
                                    test_joined["outcome_value"].sum()
                                ),
                                "prevalence": target_oracle_reference_value,
                                "brier": np.nan,
                                "log_loss": np.nan,
                                "bss_transport": np.nan,
                                "bss_target_oracle": np.nan,
                                "calibration_intercept": np.nan,
                                "calibration_slope": np.nan,
                                "calibration_intercept_se": np.nan,
                                "calibration_slope_se": np.nan,
                                "calibration_status": "failed",
                                "calibration_failure_type": "source_model_failure",
                                "calibration_message": fit.get("message", ""),
                            }
                        metric_frames.append(
                            pd.DataFrame([{**metric_base, **metrics}])
                        )

                        if fit["estimable"]:
                            for reliability in reliability_rows(
                                predictions,
                                "predicted_probability",
                                config["reliability_bin_edges"],
                            ):
                                reliability_frames.append(
                                    pd.DataFrame(
                                        [
                                            {
                                                **base,
                                                "model_form": model_form,
                                                "model_layer_label": U_LAYER,
                                                **reliability,
                                            }
                                        ]
                                    )
                                )
                        else:
                            for bin_index in range(
                                len(config["reliability_bin_edges"]) - 1
                            ):
                                reliability_frames.append(
                                    pd.DataFrame(
                                        [
                                            {
                                                **base,
                                                "model_form": model_form,
                                                "model_layer_label": U_LAYER,
                                                "bin_index": bin_index,
                                                "bin_lower": config[
                                                    "reliability_bin_edges"
                                                ][bin_index],
                                                "bin_upper": config[
                                                    "reliability_bin_edges"
                                                ][bin_index + 1],
                                                "n": 0,
                                                "positives": 0,
                                                "mean_prediction": np.nan,
                                                "observed_rate": np.nan,
                                                "wilson_95_low": np.nan,
                                                "wilson_95_high": np.nan,
                                            }
                                        ]
                                    )
                                )

    predictions = pd.concat(prediction_frames, ignore_index=True)
    parameters = pd.concat(parameter_frames, ignore_index=True)
    metrics = pd.concat(metric_frames, ignore_index=True)
    reliability = pd.concat(reliability_frames, ignore_index=True)
    knots = pd.concat(knot_frames, ignore_index=True)

    expected_rows = 16140 * 2 * 2
    if len(metrics) != expected_evaluations or expected_evaluations != 176:
        append_failure(
            STAGE,
            "stage_1_model_count",
            "incomplete_model_grid",
            "The expected U model grid was not fully enumerated.",
            expected=176,
            observed=int(expected_evaluations),
            metric_rows=int(len(metrics)),
        )
    if len(predictions) != expected_rows:
        append_failure(
            STAGE,
            "stage_1_prediction_count",
            "incomplete_pair_predictions",
            "The expected number of frozen held-out-target prediction rows is incorrect.",
            expected=expected_rows,
            observed=int(len(predictions)),
        )

    atomic_csv(RUN_DIR / "rq3_u_pair_predictions.csv", predictions)
    atomic_csv(RUN_DIR / "rq3_source_model_parameters.csv", parameters)
    atomic_csv(RUN_DIR / "rq3_loto_metrics.csv", metrics)
    atomic_csv(RUN_DIR / "rq3_loto_reliability.csv", reliability)
    atomic_csv(RUN_DIR / "rq3_rcs_knots.csv", knots)

    macro = []
    for keys, group in metrics.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "model_form",
            "model_layer_label",
        ],
        sort=True,
    ):
        base = dict(
            zip(
                [
                    "analysis_set",
                    "cluster_rule",
                    "source_weighting",
                    "model_form",
                    "model_layer_label",
                ],
                keys,
            )
        )
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            macro.append(
                {
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
            )
    macro_frame = pd.DataFrame(macro)
    atomic_csv(RUN_DIR / "rq3_loto_macro_summary.csv", macro_frame)

    pair_rows = []
    for keys, group in predictions.groupby(
        ["analysis_set", "cluster_rule", "source_weighting", "model_form"],
        sort=True,
    ):
        base = dict(
            zip(
                ["analysis_set", "cluster_rule", "source_weighting", "model_form"],
                keys,
            )
        )
        if group["predicted_probability"].isna().any():
            metrics_pair = {
                "n": int(len(group)),
                "positives": int(group["outcome_value"].sum()),
                "brier": np.nan,
                "log_loss": np.nan,
                "bss_transport": np.nan,
                "bss_target_oracle": np.nan,
            }
        else:
            metrics_pair = compute_performance_metrics(
                group["outcome_value"].to_numpy(dtype=float),
                group["predicted_probability"].to_numpy(dtype=float),
                clip=float(config["probability_clip"]),
                transport_reference=group["bss_transport_reference"].to_numpy(
                    dtype=float
                ),
                target_oracle_reference=group[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
            )
        pair_rows.append(
            {
                **base,
                "model_layer_label": U_LAYER,
                "summary_type": "pair_pooled_across_held_out_targets",
                **metrics_pair,
            }
        )
    pair_frame = pd.DataFrame(pair_rows)
    atomic_csv(RUN_DIR / "rq3_loto_pair_summary.csv", pair_frame)

    failures = int(metrics["convergence_flags"].ne("evaluable").sum())
    stage_summary = {
        "stage": STAGE,
        "completed_at": now_iso(),
        "expected_model_evaluations": 176,
        "observed_model_evaluations": int(len(metrics)),
        "observed_pair_prediction_rows": int(len(predictions)),
        "expected_pair_prediction_rows": expected_rows,
        "source_models_evaluable": int(metrics["convergence_flags"].eq("evaluable").sum()),
        "source_models_failed": failures,
        "all_predictions_finite_for_evaluable_models": bool(
            predictions.loc[
                predictions["fit_status"].eq("evaluable"), "predicted_probability"
            ].map(np.isfinite).all()
        ),
        "rcs_unit_test": rcs_test,
    }
    atomic_json(RUN_DIR / "stage1_summary.json", stage_summary)
    set_stage_status(
        STAGE,
        "complete",
        expected_models=176,
        evaluable_models=stage_summary["source_models_evaluable"],
        failed_models=failures,
        pair_prediction_rows=int(len(predictions)),
    )
    print(
        f"Stage 1 complete: {len(metrics)} model-target evaluations, "
        f"{len(predictions)} pair predictions, {failures} failed models."
    )


if __name__ == "__main__":
    main()
