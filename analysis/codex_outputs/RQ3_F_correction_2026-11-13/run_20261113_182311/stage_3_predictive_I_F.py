"""RQ3 stage 3: predictive I and nested-CV penalized F adaptation."""
from __future__ import annotations

import math
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.model_selection import StratifiedGroupKFold

from rq3_common import (
    LINEAR,
    PRED_F_LAYER,
    PRED_I_LAYER,
    RUN_DIR,
    append_failure,
    atomic_csv,
    atomic_json,
    bss_from_reference,
    compute_performance_metrics,
    fit_offset_slope,
    load_config,
    read_csv,
    reliability_rows,
    require_hash_gate,
    set_stage_status,
    solve_offset_intercept,
)

STAGE = "stage_3_predictive_I_F_with_nested_inner_cv"
METRIC_COLUMNS = [
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "calibration_intercept",
    "calibration_slope",
    "calibration_intercept_se",
    "calibration_slope_se",
]


def adaptation_unit_tests(config: dict) -> dict:
    eta = np.linspace(-2.5, 2.5, 500)
    delta_zero_p = expit(0.0 + eta + 0.0 * eta)
    exact_delta_zero = bool(
        np.array_equal(delta_zero_p, expit(eta))
    )
    rng = np.random.default_rng(20261113)
    true_probability = expit(0.15 + eta + 0.65 * eta)
    y = rng.binomial(1, true_probability)
    delta_abs = []
    statuses = []
    objective_errors = []
    identity_errors = []
    gradient_maxima = []
    fitted_records = []
    for penalty_lambda in config["predictive_adaptation"]["lambda_grid"]:
        fit = fit_offset_slope(
            eta,
            y,
            penalty_lambda=float(penalty_lambda),
            label=f"F_lambda_unit_test|{penalty_lambda}",
        )
        statuses.append(fit["status"])
        delta_abs.append(abs(float(fit["delta"])) if fit["estimable"] else np.nan)
        if fit["estimable"]:
            fitted_records.append((float(penalty_lambda), fit))
            eta_f = (
                float(fit["intercept"])
                + eta
                + float(fit["delta"]) * eta
            )
            reconstructed_objective = float(
                np.mean(np.logaddexp(0.0, eta_f) - y * eta_f)
                + 0.5 * float(penalty_lambda) * float(fit["delta"]) ** 2
            )
            objective_errors.append(
                abs(reconstructed_objective - float(fit["objective"]))
            )
            identity_errors.append(
                abs(
                    float(fit["gamma"])
                    - (1.0 + float(fit["delta"]))
                )
            )
            gradient_maxima.append(float(fit["gradient_max_abs"]))
    finite = np.asarray(delta_abs, dtype=float)
    nonincreasing = bool(
        np.isfinite(finite).all()
        and np.all(np.diff(finite) <= 1e-10)
    )

    def independent_objective(theta: np.ndarray, penalty_lambda: float) -> float:
        linear_predictor = theta[0] + eta + theta[1] * eta
        return float(
            np.mean(np.logaddexp(0.0, linear_predictor) - y * linear_predictor)
            + 0.5 * penalty_lambda * theta[1] ** 2
        )

    def independent_gradient(theta: np.ndarray, penalty_lambda: float) -> np.ndarray:
        linear_predictor = theta[0] + eta + theta[1] * eta
        residual = expit(linear_predictor) - y
        return np.asarray(
            [
                np.mean(residual),
                np.mean(residual * eta) + penalty_lambda * theta[1],
            ],
            dtype=float,
        )

    finite_difference_error = np.nan
    if fitted_records:
        _, central_fit = fitted_records[len(fitted_records) // 2]
        theta_fit = np.asarray(
            [central_fit["intercept"], central_fit["delta"]], dtype=float
        )
        step = 1e-6
        finite_difference = np.empty(2, dtype=float)
        for index in range(2):
            plus = theta_fit.copy()
            minus = theta_fit.copy()
            plus[index] += step
            minus[index] -= step
            finite_difference[index] = (
                independent_objective(plus, 1.0)
                - independent_objective(minus, 1.0)
            ) / (2.0 * step)
        analytic_gradient = independent_gradient(theta_fit, 1.0)
        finite_difference_error = float(
            np.max(np.abs(finite_difference - analytic_gradient))
        )

    lambda_zero_fit = fit_offset_slope(
        eta,
        y,
        penalty_lambda=0.0,
        label="F_lambda_zero_independent_check",
    )

    def gamma_objective(values: np.ndarray) -> float:
        linear_predictor = values[0] + values[1] * eta
        return float(
            np.mean(np.logaddexp(0.0, linear_predictor) - y * linear_predictor)
        )

    def gamma_gradient(values: np.ndarray) -> np.ndarray:
        residual = expit(values[0] + values[1] * eta) - y
        return np.asarray([np.mean(residual), np.mean(residual * eta)])

    independent_zero = minimize(
        gamma_objective,
        np.asarray(
            [
                lambda_zero_fit["intercept"]
                if lambda_zero_fit["estimable"]
                else 0.0,
                1.0,
            ]
        ),
        method="BFGS",
        jac=gamma_gradient,
        options={"gtol": 1e-10, "maxiter": 5000},
    )
    lambda_zero_prediction_error = np.nan
    if lambda_zero_fit["estimable"] and independent_zero.success:
        lambda_zero_p = expit(
            lambda_zero_fit["intercept"]
            + lambda_zero_fit["gamma"] * eta
        )
        independent_zero_p = expit(
            independent_zero.x[0] + independent_zero.x[1] * eta
        )
        lambda_zero_prediction_error = float(
            np.max(np.abs(lambda_zero_p - independent_zero_p))
        )

    large_lambda = 1e8
    large_lambda_fit = fit_offset_slope(
        eta,
        y,
        penalty_lambda=large_lambda,
        label=f"F_large_lambda_limit|{large_lambda:g}",
    )
    large_lambda_abs_delta = (
        abs(float(large_lambda_fit["delta"]))
        if large_lambda_fit["estimable"]
        else np.nan
    )

    checks = {
        "delta_zero_exactly_reproduces_U": exact_delta_zero,
        "abs_delta_non_increasing_with_lambda": nonincreasing,
        "all_synthetic_fits_evaluable": all(
            item == "evaluable" for item in statuses
        ),
        "saved_objective_reconstructs_from_parameters": bool(
            objective_errors
            and max(objective_errors) <= 1e-10
        ),
        "saved_gamma_equals_one_plus_delta": bool(
            identity_errors and max(identity_errors) == 0.0
        ),
        "saved_solution_gradient_near_zero": bool(
            gradient_maxima and max(gradient_maxima) <= 1e-5
        ),
        "analytic_gradient_matches_central_finite_difference": bool(
            np.isfinite(finite_difference_error)
            and finite_difference_error <= 1e-6
        ),
        "lambda_zero_matches_independent_gamma_fit": bool(
            lambda_zero_fit["estimable"]
            and independent_zero.success
            and np.isfinite(lambda_zero_prediction_error)
            and lambda_zero_prediction_error <= 1e-9
        ),
        "large_lambda_compresses_delta_toward_zero": bool(
            large_lambda_fit["estimable"]
            and np.isfinite(large_lambda_abs_delta)
            and large_lambda_abs_delta <= 1e-5
        ),
    }
    return {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "lambda_grid": config["predictive_adaptation"]["lambda_grid"],
        "abs_delta_by_lambda": finite.tolist(),
        "fit_statuses": statuses,
        "numerical_diagnostics": {
            "max_objective_reconstruction_error": (
                float(max(objective_errors)) if objective_errors else np.nan
            ),
            "max_gamma_identity_error": (
                float(max(identity_errors)) if identity_errors else np.nan
            ),
            "max_saved_gradient": (
                float(max(gradient_maxima)) if gradient_maxima else np.nan
            ),
            "max_finite_difference_gradient_error": finite_difference_error,
            "lambda_zero_max_probability_error": lambda_zero_prediction_error,
            "large_lambda": large_lambda,
            "large_lambda_abs_delta": large_lambda_abs_delta,
            "large_lambda_note": (
                "finite-lambda compression check only; exact equality with I "
                "is not required"
            ),
        },
        "tolerances": {
            "objective_reconstruction_abs": 1e-10,
            "gamma_identity_abs": 0.0,
            "saved_gradient_max_abs": 1e-5,
            "finite_difference_gradient_max_abs": 1e-6,
            "lambda_zero_probability_max_abs": 1e-9,
            "large_lambda_abs_delta": 1e-5,
        },
        "synthetic_seed": 20261113,
    }


def metric_row(
    *,
    base: dict,
    y: np.ndarray,
    prediction: np.ndarray,
    clip: float,
    transport_reference: np.ndarray,
    target_oracle_reference: np.ndarray,
    adaptation_reference: np.ndarray,
) -> dict:
    finite = bool(np.isfinite(prediction).all())
    if not finite:
        return {
            **base,
            "convergence_flags": "failed",
            "n": int(len(y)),
            "positives": int(np.sum(y)),
            "brier": np.nan,
            "log_loss": np.nan,
            "bss_transport": np.nan,
            "bss_target_oracle": np.nan,
            "bss_adaptation_train": np.nan,
            "calibration_intercept": np.nan,
            "calibration_slope": np.nan,
            "calibration_intercept_se": np.nan,
            "calibration_slope_se": np.nan,
            "calibration_status": "failed",
            "calibration_failure_type": "prediction_failure",
            "calibration_message": "Prediction is nonfinite",
        }
    performance = compute_performance_metrics(
        y,
        prediction,
        clip=clip,
        transport_reference=transport_reference,
        target_oracle_reference=target_oracle_reference,
    )
    performance["bss_adaptation_train"] = bss_from_reference(
        y, prediction, adaptation_reference
    )
    return {
        **base,
        "convergence_flags": "evaluable",
        **performance,
    }


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    require_hash_gate(config)
    if float(config["probability_clip"]) != 1e-6:
        raise RuntimeError("Unexpected frozen probability clip")
    clip = float(config["probability_clip"])
    lambda_grid = [float(value) for value in config["predictive_adaptation"]["lambda_grid"]]

    unit_tests = adaptation_unit_tests(config)
    atomic_json(RUN_DIR / "adaptation_unit_tests.json", unit_tests)
    if unit_tests["status"] != "passed":
        append_failure(
            STAGE,
            "adaptation_unit_tests",
            "unit_test_failure",
            "Predictive adaptation formula/penalty unit tests failed.",
            checks=unit_tests["checks"],
        )
        raise RuntimeError("Predictive adaptation unit tests failed")

    u = read_csv(RUN_DIR / "rq3_u_pair_predictions.csv")
    u = u[
        u["model_form"].eq(LINEAR) & u["fit_status"].eq("evaluable")
    ].copy()
    cv5 = read_csv(
        __import__("pathlib").Path(
            config["inputs"]["target_cv5_assignments"]["path"]
        )
    )

    prediction_frames: list[pd.DataFrame] = []
    inner_frames: list[pd.DataFrame] = []
    lambda_frames: list[pd.DataFrame] = []
    fold_metric_frames: list[dict] = []
    parameter_frames: list[dict] = []
    expected_outer_units = 44 * 2
    expected_cv_folds = expected_outer_units * 5
    observed_outer_units = 0

    group_columns = [
        "analysis_set",
        "cluster_rule",
        "source_weighting",
        "held_out_target",
    ]
    for keys, target_u in u.groupby(group_columns, sort=True):
        analysis_set, cluster_rule, source_weighting, held_target = keys
        observed_outer_units += 1
        target_u = target_u.sort_values("pair_id").reset_index(drop=True)
        if target_u["pair_id"].duplicated().any():
            append_failure(
                STAGE,
                f"{analysis_set}|{cluster_rule}|{source_weighting}|{held_target}",
                "duplicate_fixed_U_rows",
                "Fixed U predictions contain duplicate pair rows.",
            )
            continue
        outcome_name = str(target_u["outcome_name"].iloc[0])
        target_cv = cv5[
            cv5["analysis_set"].eq(analysis_set)
            & cv5["outcome_name"].eq(outcome_name)
            & cv5["cluster_rule"].eq(cluster_rule)
            & cv5["target"].astype(str).eq(str(held_target))
        ].copy()
        if target_cv.empty:
            append_failure(
                STAGE,
                f"{analysis_set}|{cluster_rule}|{source_weighting}|{held_target}",
                "missing_frozen_cv5",
                "No frozen CV5 assignments match the LOTO target.",
            )
            continue

        target_index = target_u.set_index("pair_id")
        outer_fold_ids = sorted(target_cv["fold_id"].unique())
        if len(outer_fold_ids) != 5:
            append_failure(
                STAGE,
                f"{analysis_set}|{cluster_rule}|{source_weighting}|{held_target}",
                "wrong_outer_fold_count",
                "Frozen CV5 target does not have five outer folds.",
                observed=len(outer_fold_ids),
            )
            continue

        for outer_fold_id in outer_fold_ids:
            fold_rows = target_cv[target_cv["fold_id"].eq(outer_fold_id)].copy()
            train_cv = fold_rows[fold_rows["split"].eq("train")].copy()
            test_cv = fold_rows[fold_rows["split"].eq("test")].copy()
            if train_cv.empty or test_cv.empty:
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}",
                    "empty_outer_fold",
                    "Outer adaptation train or test set is empty.",
                )
                continue
            if set(train_cv["sequence_cluster_value"].astype(str)) & set(
                test_cv["sequence_cluster_value"].astype(str)
            ):
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}",
                    "outer_cluster_leak",
                    "Outer train/test sequence clusters overlap.",
                )
            train = target_index.loc[train_cv["pair_id"].astype(str)].reset_index()
            test = target_index.loc[test_cv["pair_id"].astype(str)].reset_index()
            train = train.merge(
                train_cv[
                    [
                        "pair_id",
                        "outcome_value",
                        "sequence_cluster_value",
                    ]
                ].rename(
                    columns={
                        "outcome_value": "cv_outcome_value",
                        "sequence_cluster_value": "cv_group",
                    }
                ),
                on="pair_id",
                how="left",
                validate="one_to_one",
            )
            test = test.merge(
                test_cv[
                    [
                        "pair_id",
                        "outcome_value",
                        "sequence_cluster_value",
                    ]
                ].rename(
                    columns={
                        "outcome_value": "cv_outcome_value",
                        "sequence_cluster_value": "cv_group",
                    }
                ),
                on="pair_id",
                how="left",
                validate="one_to_one",
            )
            for role, frame in (("train", train), ("test", test)):
                if frame[["cv_outcome_value", "cv_group"]].isna().any().any():
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{source_weighting}|{role}",
                        "cv_join_missing",
                        "Frozen CV assignment did not join to fixed U row.",
                    )
                if not frame["outcome_value"].astype(int).eq(
                    frame["cv_outcome_value"].astype(int)
                ).all():
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{source_weighting}|{role}",
                        "outcome_mismatch",
                        "Frozen CV outcome disagrees with master outcome.",
                    )
            if (
                train[["cv_outcome_value", "cv_group"]].isna().any().any()
                or test[["cv_outcome_value", "cv_group"]].isna().any().any()
            ):
                continue

            eta_train = train["eta_U"].to_numpy(dtype=float)
            y_train = train["cv_outcome_value"].to_numpy(dtype=int)
            eta_test = test["eta_U"].to_numpy(dtype=float)
            y_test = test["cv_outcome_value"].to_numpy(dtype=int)
            adaptation_train_prevalence = float(np.mean(y_train))
            adaptation_reference_test = np.full(
                len(test), adaptation_train_prevalence, dtype=float
            )
            transport_reference_test = test[
                "bss_transport_reference"
            ].to_numpy(dtype=float)
            target_oracle_reference_test = test[
                "bss_target_oracle_reference"
            ].to_numpy(dtype=float)

            i_fit = solve_offset_intercept(
                eta_train,
                y_train,
                f"{outer_fold_id}|{source_weighting}|predictive_I",
            )
            if i_fit["status"] == "evaluable":
                eta_i_test = float(i_fit["intercept"]) + eta_test
                p_i_test = expit(eta_i_test)
            else:
                eta_i_test = np.full(len(test), np.nan)
                p_i_test = np.full(len(test), np.nan)
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}|predictive_I",
                    i_fit["failure_type"],
                    i_fit["message"],
                )
            parameter_frames.append(
                {
                    "outer_fold_id": outer_fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "held_out_target": held_target,
                    "source_weighting": source_weighting,
                    "model_layer_label": PRED_I_LAYER,
                    "fit_status": i_fit["status"],
                    "failure_type": i_fit.get("failure_type", ""),
                    "selected_lambda": np.nan,
                    "intercept": i_fit.get("intercept", np.nan),
                    "gamma": 1.0,
                    "delta": 0.0,
                    "objective": np.nan,
                    "iterations": np.nan,
                    "gradient_max_abs": np.nan,
                    "n_adaptation_train": int(len(train)),
                    "n_outer_test": int(len(test)),
                    "adaptation_train_prevalence": adaptation_train_prevalence,
                }
            )

            splitter = StratifiedGroupKFold(
                n_splits=int(config["predictive_adaptation"]["inner_cv"]["n_splits"]),
                shuffle=True,
                random_state=int(
                    config["predictive_adaptation"]["inner_cv"]["random_state"]
                ),
            )
            inner_splits = list(
                splitter.split(
                    np.zeros((len(train), 1)),
                    y_train,
                    groups=train["cv_group"].astype(str).to_numpy(),
                )
            )
            if len(inner_splits) != 3:
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}|inner_cv",
                    "inner_fold_formation_failure",
                    "StratifiedGroupKFold did not form three inner folds.",
                )
            for inner_number, (inner_train_indices, inner_test_indices) in enumerate(
                inner_splits, start=1
            ):
                inner_train_index_set = set(inner_train_indices.tolist())
                inner_train_groups = set(
                    train.iloc[inner_train_indices]["cv_group"].astype(str)
                )
                inner_test_groups = set(
                    train.iloc[inner_test_indices]["cv_group"].astype(str)
                )
                if inner_train_groups & inner_test_groups:
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{source_weighting}|inner_{inner_number}",
                        "inner_cluster_leak",
                        "Inner train/test sequence clusters overlap.",
                    )
                for row_number, (_, row) in enumerate(train.iterrows()):
                    inner_frames.append(
                        {
                            "outer_fold_id": outer_fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": cluster_rule,
                            "held_out_target": held_target,
                            "source_weighting": source_weighting,
                            "inner_fold_number": inner_number,
                            "pair_id": row["pair_id"],
                            "outcome_value": int(row["cv_outcome_value"]),
                            "sequence_cluster_value": row["cv_group"],
                            "split": (
                                "train"
                                if row_number in inner_train_index_set
                                else "test"
                            ),
                        }
                    )

            lambda_records = []
            for penalty_lambda in lambda_grid:
                pooled_y = []
                pooled_p = []
                candidate_failure = ""
                candidate_objectives = []
                for inner_number, (inner_train_indices, inner_test_indices) in enumerate(
                    inner_splits, start=1
                ):
                    inner_train = train.iloc[inner_train_indices]
                    inner_test = train.iloc[inner_test_indices]
                    if inner_train["cv_outcome_value"].nunique() < 2:
                        candidate_failure = f"inner_{inner_number}_single_class"
                        break
                    fit = fit_offset_slope(
                        inner_train["eta_U"].to_numpy(dtype=float),
                        inner_train["cv_outcome_value"].to_numpy(dtype=float),
                        penalty_lambda=penalty_lambda,
                        label=(
                            f"{outer_fold_id}|{source_weighting}|"
                            f"inner_{inner_number}|lambda_{penalty_lambda:g}"
                        ),
                    )
                    if not fit["estimable"]:
                        candidate_failure = (
                            f"inner_{inner_number}_{fit['failure_type']}"
                        )
                        break
                    eta_inner_test = (
                        fit["intercept"]
                        + fit["gamma"]
                        * inner_test["eta_U"].to_numpy(dtype=float)
                    )
                    pooled_y.append(
                        inner_test["cv_outcome_value"].to_numpy(dtype=float)
                    )
                    pooled_p.append(expit(eta_inner_test))
                    candidate_objectives.append(float(fit["objective"]))
                if candidate_failure:
                    lambda_records.append(
                        {
                            "outer_fold_id": outer_fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": cluster_rule,
                            "held_out_target": held_target,
                            "source_weighting": source_weighting,
                            "lambda": penalty_lambda,
                            "candidate_status": "unevaluable",
                            "failure_type": candidate_failure,
                            "pooled_inner_log_loss": np.nan,
                            "n_pooled_inner_predictions": 0,
                            "mean_inner_objective": np.nan,
                            "selected": False,
                            "tie_rule": "minimum pooled inner log loss; 1e-12 ties choose larger lambda",
                        }
                    )
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{source_weighting}|lambda_{penalty_lambda:g}",
                        candidate_failure,
                        "Candidate lambda failed nested inner validation QC.",
                    )
                    continue
                y_pooled = np.concatenate(pooled_y)
                p_pooled = np.concatenate(pooled_p)
                pooled_log_loss = float(
                    -np.mean(
                        y_pooled * np.log(np.clip(p_pooled, clip, 1.0 - clip))
                        + (1.0 - y_pooled)
                        * np.log1p(-np.clip(p_pooled, clip, 1.0 - clip))
                    )
                )
                lambda_records.append(
                    {
                        "outer_fold_id": outer_fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": cluster_rule,
                        "held_out_target": held_target,
                        "source_weighting": source_weighting,
                        "lambda": penalty_lambda,
                        "candidate_status": "evaluable",
                        "failure_type": "",
                        "pooled_inner_log_loss": pooled_log_loss,
                        "n_pooled_inner_predictions": int(len(y_pooled)),
                        "mean_inner_objective": float(
                            np.mean(candidate_objectives)
                        ),
                        "selected": False,
                        "tie_rule": "minimum pooled inner log loss; 1e-12 ties choose larger lambda",
                    }
                )

            evaluable_records = [
                row
                for row in lambda_records
                if row["candidate_status"] == "evaluable"
            ]
            selected_lambda = np.nan
            if evaluable_records:
                best_loss = min(
                    row["pooled_inner_log_loss"] for row in evaluable_records
                )
                tied = [
                    row
                    for row in evaluable_records
                    if row["pooled_inner_log_loss"] - best_loss <= 1e-12
                ]
                selected_record = max(tied, key=lambda row: row["lambda"])
                selected_lambda = float(selected_record["lambda"])
                for row in lambda_records:
                    row["selected"] = bool(
                        row["candidate_status"] == "evaluable"
                        and row["lambda"] == selected_lambda
                    )
                final_fit = fit_offset_slope(
                    eta_train,
                    y_train,
                    penalty_lambda=selected_lambda,
                    label=(
                        f"{outer_fold_id}|{source_weighting}|"
                        f"final_F_lambda_{selected_lambda:g}"
                    ),
                )
            else:
                final_fit = {
                    "status": "failed",
                    "estimable": False,
                    "failure_type": "all_lambdas_unevaluable",
                    "message": "No lambda candidate was evaluable.",
                    "intercept": np.nan,
                    "gamma": np.nan,
                    "delta": np.nan,
                    "objective": np.nan,
                    "iterations": np.nan,
                    "gradient_max_abs": np.nan,
                }
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}|predictive_F",
                    "all_lambdas_unevaluable",
                    "All prespecified lambda candidates failed inner validation.",
                )
            lambda_frames.extend(lambda_records)

            if final_fit.get("estimable", False):
                eta_f_test = (
                    float(final_fit["intercept"])
                    + float(final_fit["gamma"]) * eta_test
                )
                p_f_test = expit(eta_f_test)
                f_status = "evaluable"
                f_failure = ""
            else:
                eta_f_test = np.full(len(test), np.nan)
                p_f_test = np.full(len(test), np.nan)
                f_status = "failed"
                f_failure = final_fit.get("failure_type", "final_fit_failure")
                append_failure(
                    STAGE,
                    f"{outer_fold_id}|{source_weighting}|predictive_F",
                    f_failure,
                    final_fit.get("message", "Final F fit failed."),
                )
            parameter_frames.append(
                {
                    "outer_fold_id": outer_fold_id,
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "held_out_target": held_target,
                    "source_weighting": source_weighting,
                    "model_layer_label": PRED_F_LAYER,
                    "fit_status": f_status,
                    "failure_type": f_failure,
                    "selected_lambda": selected_lambda,
                    "intercept": final_fit.get("intercept", np.nan),
                    "gamma": final_fit.get("gamma", np.nan),
                    "delta": final_fit.get("delta", np.nan),
                    "objective": final_fit.get("objective", np.nan),
                    "iterations": final_fit.get("iterations", np.nan),
                    "gradient_max_abs": final_fit.get("gradient_max_abs", np.nan),
                    "n_adaptation_train": int(len(train)),
                    "n_outer_test": int(len(test)),
                    "adaptation_train_prevalence": adaptation_train_prevalence,
                }
            )

            prediction_frames.append(
                pd.DataFrame(
                    {
                        "pair_id": test["pair_id"].astype(str).to_numpy(),
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "outer_fold_id": outer_fold_id,
                        "cluster_rule": cluster_rule,
                        "held_out_target": held_target,
                        "target": test["target"].astype(str).to_numpy(),
                        "sequence_cluster_value": test["cv_group"].astype(str),
                        "outcome_value": y_test,
                        "boltz2_iptm": test["boltz2_iptm"].to_numpy(dtype=float),
                        "x": test["x"].to_numpy(dtype=float),
                        "source_weighting": source_weighting,
                        "source_U_fold_id": test["fold_id"].astype(str).to_numpy(),
                        "eta_U": eta_test,
                        "p_U": test["predicted_probability"].to_numpy(dtype=float),
                        "eta_I": eta_i_test,
                        "p_I": p_i_test,
                        "eta_F": eta_f_test,
                        "p_F": p_f_test,
                        "adaptation_train_prevalence": adaptation_train_prevalence,
                        "bss_adaptation_train_reference": adaptation_reference_test,
                        "bss_transport_reference": transport_reference_test,
                        "bss_target_oracle_reference": target_oracle_reference_test,
                        "selected_lambda": selected_lambda,
                        "I_status": i_fit["status"],
                        "F_status": f_status,
                        "I_failure_type": i_fit.get("failure_type", ""),
                        "F_failure_type": f_failure,
                        "same_outer_test_rows": True,
                    }
                )
            )

            for layer, probability, status in [
                ("U_transport_zero_label", test["predicted_probability"].to_numpy(dtype=float), "evaluable"),
                (PRED_I_LAYER, p_i_test, i_fit["status"]),
                (PRED_F_LAYER, p_f_test, f_status),
            ]:
                row = metric_row(
                    base={
                        "outer_fold_id": outer_fold_id,
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": cluster_rule,
                        "held_out_target": held_target,
                        "source_weighting": source_weighting,
                        "model_layer_label": layer,
                        "selected_lambda": (
                            selected_lambda if layer == PRED_F_LAYER else np.nan
                        ),
                        "adaptation_train_prevalence": adaptation_train_prevalence,
                    },
                    y=y_test.astype(float),
                    prediction=probability,
                    clip=clip,
                    transport_reference=transport_reference_test,
                    target_oracle_reference=target_oracle_reference_test,
                    adaptation_reference=adaptation_reference_test,
                )
                row["convergence_flags"] = status
                fold_metric_frames.append(row)
                if status == "evaluable" and row.get("calibration_status") != "evaluable":
                    append_failure(
                        STAGE,
                        f"{outer_fold_id}|{source_weighting}|{layer}|calibration",
                        "calibration_failure",
                        row.get("calibration_message", "Calibration failed."),
                    )

    predictions = pd.concat(prediction_frames, ignore_index=True)
    inner_assignments = pd.DataFrame(inner_frames)
    lambdas = pd.DataFrame(lambda_frames)
    fold_metrics = pd.DataFrame(fold_metric_frames)
    parameters = pd.DataFrame(parameter_frames)
    if observed_outer_units != expected_outer_units:
        append_failure(
            STAGE,
            "outer_unit_count",
            "incomplete_outer_grid",
            "Predictive adaptation did not enumerate all expected target/config units.",
            expected=expected_outer_units,
            observed=observed_outer_units,
        )
    if len(predictions) != 32280:
        append_failure(
            STAGE,
            "outer_test_prediction_count",
            "incomplete_pair_predictions",
            "Predictive outer-test pair count is incorrect.",
            expected=32280,
            observed=int(len(predictions)),
        )
    if len(fold_metrics) != expected_cv_folds * 3:
        append_failure(
            STAGE,
            "outer_fold_metric_count",
            "incomplete_fold_metrics",
            "Predictive outer-fold metric count is incorrect.",
            expected=expected_cv_folds * 3,
            observed=int(len(fold_metrics)),
        )
    if len(lambdas) != expected_outer_units * 5 * len(lambda_grid):
        append_failure(
            STAGE,
            "lambda_selection_count",
            "incomplete_lambda_grid",
            "Lambda selection grid count is incorrect.",
            expected=expected_outer_units * 5 * len(lambda_grid),
            observed=int(len(lambdas)),
        )

    atomic_csv(RUN_DIR / "rq3_adaptation_pair_predictions.csv", predictions)
    atomic_csv(RUN_DIR / "rq3_adaptation_inner_assignments.csv", inner_assignments)
    atomic_csv(RUN_DIR / "rq3_lambda_selection.csv", lambdas)
    atomic_csv(RUN_DIR / "rq3_adaptation_fold_metrics.csv", fold_metrics)
    atomic_csv(RUN_DIR / "rq3_adaptation_parameters.csv", parameters)

    target_metric_rows = []
    reliability_frames = []
    target_group_columns = [
        "analysis_set",
        "cluster_rule",
        "source_weighting",
        "held_out_target",
    ]
    for keys, group in predictions.groupby(target_group_columns, sort=True):
        analysis_set, cluster_rule, source_weighting, held_target = keys
        group = group.sort_values("pair_id").reset_index(drop=True)
        y_target = group["outcome_value"].to_numpy(dtype=float)
        for layer, probability_column in [
            ("U_transport_zero_label", "p_U"),
            (PRED_I_LAYER, "p_I"),
            (PRED_F_LAYER, "p_F"),
        ]:
            row = metric_row(
                base={
                    "analysis_set": analysis_set,
                    "outcome_name": str(group["outcome_name"].iloc[0]),
                    "cluster_rule": cluster_rule,
                    "held_out_target": held_target,
                    "source_weighting": source_weighting,
                    "model_layer_label": layer,
                    "target_oof_fold_count": int(
                        group["outer_fold_id"].nunique()
                    ),
                },
                y=y_target,
                prediction=group[probability_column].to_numpy(dtype=float),
                clip=clip,
                transport_reference=group["bss_transport_reference"].to_numpy(
                    dtype=float
                ),
                target_oracle_reference=group[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
                adaptation_reference=group[
                    "bss_adaptation_train_reference"
                ].to_numpy(dtype=float),
            )
            target_metric_rows.append(row)
            if row["convergence_flags"] == "evaluable":
                reliability_frame = pd.DataFrame(
                    {
                        "outcome_value": y_target,
                        "predicted_probability": group[
                            probability_column
                        ].to_numpy(dtype=float),
                    }
                )
                for reliability in reliability_rows(
                    reliability_frame,
                    "predicted_probability",
                    config["reliability_bin_edges"],
                ):
                    reliability_frames.append(
                        {
                            "analysis_set": analysis_set,
                            "outcome_name": str(group["outcome_name"].iloc[0]),
                            "cluster_rule": cluster_rule,
                            "held_out_target": held_target,
                            "source_weighting": source_weighting,
                            "model_layer_label": layer,
                            **reliability,
                        }
                    )
    target_metrics = pd.DataFrame(target_metric_rows)
    atomic_csv(RUN_DIR / "rq3_target_adaptation_metrics.csv", target_metrics)
    atomic_csv(
        RUN_DIR / "rq3_adaptation_reliability.csv",
        pd.DataFrame(reliability_frames),
    )

    macro_rows = []
    for keys, group in target_metrics.groupby(
        ["analysis_set", "cluster_rule", "source_weighting", "model_layer_label"],
        sort=True,
    ):
        base = dict(
            zip(
                [
                    "analysis_set",
                    "cluster_rule",
                    "source_weighting",
                    "model_layer_label",
                ],
                keys,
            )
        )
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            macro_rows.append(
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
    macro = pd.DataFrame(macro_rows)
    atomic_csv(RUN_DIR / "rq3_adaptation_macro_summary.csv", macro)

    pair_rows = []
    for keys, group in predictions.groupby(
        ["analysis_set", "cluster_rule", "source_weighting"], sort=True
    ):
        base = dict(
            zip(["analysis_set", "cluster_rule", "source_weighting"], keys)
        )
        for layer, probability_column in [
            ("U_transport_zero_label", "p_U"),
            (PRED_I_LAYER, "p_I"),
            (PRED_F_LAYER, "p_F"),
        ]:
            pair_rows.append(
                {
                    **base,
                    "model_layer_label": layer,
                    "summary_type": "pair_pooled_across_held_out_targets",
                    **metric_row(
                        base={},
                        y=group["outcome_value"].to_numpy(dtype=float),
                        prediction=group[probability_column].to_numpy(dtype=float),
                        clip=clip,
                        transport_reference=group[
                            "bss_transport_reference"
                        ].to_numpy(dtype=float),
                        target_oracle_reference=group[
                            "bss_target_oracle_reference"
                        ].to_numpy(dtype=float),
                        adaptation_reference=group[
                            "bss_adaptation_train_reference"
                        ].to_numpy(dtype=float),
                    ),
                }
            )
    pair_summary = pd.DataFrame(pair_rows)
    atomic_csv(RUN_DIR / "rq3_adaptation_pair_summary.csv", pair_summary)

    comparison_rows = []
    for keys, group in target_metrics.groupby(
        ["analysis_set", "cluster_rule", "source_weighting"], sort=True
    ):
        for comparison_name, left, right in [
            ("predictive_I_minus_U", PRED_I_LAYER, "U_transport_zero_label"),
            ("predictive_F_minus_U", PRED_F_LAYER, "U_transport_zero_label"),
            ("predictive_F_minus_I", PRED_F_LAYER, PRED_I_LAYER),
        ]:
            if not {left, right}.issubset(
                set(group["model_layer_label"])
            ):
                continue
            for metric in ("brier", "log_loss", "bss_target_oracle"):
                pivot = group.pivot(
                    index="held_out_target",
                    columns="model_layer_label",
                    values=metric,
                )
                difference = (
                    pd.to_numeric(pivot[left], errors="coerce")
                    - pd.to_numeric(pivot[right], errors="coerce")
                ).dropna()
                comparison_rows.append(
                    {
                        "analysis_set": keys[0],
                        "cluster_rule": keys[1],
                        "source_weighting": keys[2],
                        "comparison": comparison_name,
                        "metric": metric,
                        "difference_direction": (
                            "negative_is_better"
                            if metric in {"brier", "log_loss"}
                            else "positive_is_better"
                        ),
                        "target_count": int(len(difference)),
                        "macro_mean_difference": float(difference.mean()),
                        "n_improved": int(
                            (difference < 0).sum()
                            if metric in {"brier", "log_loss"}
                            else (difference > 0).sum()
                        ),
                        "n_worsened": int(
                            (difference > 0).sum()
                            if metric in {"brier", "log_loss"}
                            else (difference < 0).sum()
                        ),
                        "min_difference": float(difference.min()),
                        "max_difference": float(difference.max()),
                    }
                )
    comparisons = pd.DataFrame(comparison_rows)
    atomic_csv(RUN_DIR / "rq3_adaptation_comparisons.csv", comparisons)

    selected = lambdas[lambdas["selected"]].copy()
    failed_fold_layers = int(
        fold_metrics["convergence_flags"].ne("evaluable").sum()
    )
    stage_summary = {
        "stage": STAGE,
        "outer_target_config_units": int(observed_outer_units),
        "expected_outer_target_config_units": int(expected_outer_units),
        "outer_cv_fold_count": int(predictions["outer_fold_id"].nunique()),
        "outer_test_pair_predictions": int(len(predictions)),
        "inner_assignment_rows": int(len(inner_assignments)),
        "lambda_candidates": int(len(lambdas)),
        "lambda_candidates_evaluable": int(
            lambdas["candidate_status"].eq("evaluable").sum()
        ),
        "selected_lambda_count": int(len(selected)),
        "selected_lambda_frequencies": {
            str(key): int(value)
            for key, value in selected["lambda"].value_counts().sort_index().items()
        },
        "fold_metric_rows": int(len(fold_metrics)),
        "failed_fold_layer_evaluations": failed_fold_layers,
        "target_metric_rows": int(len(target_metrics)),
        "same_outer_test_rows_for_U_I_F": bool(
            predictions[["p_U", "p_I", "p_F"]].shape[0]
            == len(predictions)
        ),
        "adaptation_unit_tests": unit_tests,
        "labels": {
            "predictive_I": PRED_I_LAYER,
            "predictive_F": PRED_F_LAYER,
        },
    }
    atomic_json(RUN_DIR / "stage3_summary.json", stage_summary)
    set_stage_status(
        STAGE,
        "complete",
        outer_target_config_units=int(observed_outer_units),
        outer_cv_folds=int(predictions["outer_fold_id"].nunique()),
        lambda_candidates=int(len(lambdas)),
        failed_fold_layer_evaluations=failed_fold_layers,
    )
    print(
        f"Stage 3 complete: {len(predictions)} outer-test pair predictions, "
        f"{len(lambdas)} lambda candidates, {failed_fold_layers} failed fold-layer evaluations."
    )


if __name__ == "__main__":
    main()
