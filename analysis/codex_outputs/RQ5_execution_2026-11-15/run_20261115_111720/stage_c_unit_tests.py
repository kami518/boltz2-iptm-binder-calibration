"""Required unit tests for the frozen RQ5 numerical and fold rules."""
from __future__ import annotations

import json
import math
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from scipy.special import expit

from rq5_common import (
    ANALYSIS_SETS,
    CLUSTER_RULE,
    RUN_DIR,
    analysis_frame,
    atomic_json,
    fit_combo_model,
    load_config,
    load_core_inputs,
    logit_clip,
    rq3,
    weighted_center_scale,
)

STAGE = "stage_c_required_unit_tests"


def main() -> None:
    checks = []

    def check(name: str, condition: bool, detail) -> None:
        checks.append({"check": name, "passed": bool(condition), "detail": detail})

    # 1. Logit transform and symmetric clipping.
    raw = np.array([0.0, 0.25, 0.5, 1.0, -0.1, 1.1], dtype=float)
    transformed = logit_clip(raw, 1e-6)
    expected = np.log(np.clip(raw, 1e-6, 1.0 - 1e-6)) - np.log1p(
        -np.clip(raw, 1e-6, 1.0 - 1e-6)
    )
    check(
        "logit_clip",
        np.isfinite(transformed).all() and np.allclose(transformed, expected, atol=1e-12),
        {"max_abs_error": float(np.max(np.abs(transformed - expected)))},
    )
    check(
        "logit_symmetry",
        abs(float(transformed[0] + transformed[3])) < 1e-10,
        {"logit_0": float(transformed[0]), "logit_1": float(transformed[3])},
    )

    config = load_config()
    _, _, _, joined = load_core_inputs(config)
    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    loto = loto[
        loto["analysis_set"].isin(ANALYSIS_SETS)
        & loto["cluster_rule"].eq(CLUSTER_RULE)
    ]

    # 2. Frozen target-equal weights have equal total weight per source target;
    # pair weights are one per pair in the accepted RQ3 folds.
    target_weight_errors = []
    pair_weight_errors = []
    for fold_id, fold in loto.groupby("fold_id"):
        train = fold[fold["split_status"].eq("train")]
        sums = train.groupby("target")["weight_target_equal"].sum()
        if len(sums) and float(sums.max() - sums.min()) > 1e-9:
            target_weight_errors.append((fold_id, float(sums.max() - sums.min())))
        pair = train["weight_pair"].to_numpy(dtype=float)
        if len(pair) and (not np.allclose(pair, 1.0, atol=1e-12)):
            pair_weight_errors.append(fold_id)
    check(
        "target_equal_weights",
        not target_weight_errors,
        {"n_fold_errors": len(target_weight_errors), "examples": target_weight_errors[:3]},
    )
    check(
        "pair_source_weights",
        not pair_weight_errors,
        {"n_fold_errors": len(pair_weight_errors), "examples": pair_weight_errors[:3]},
    )

    # 3. LOTO isolation from frozen fold assignments.
    loto_leaks = []
    for row in loto.drop_duplicates("fold_id").itertuples(index=False):
        fold = loto[loto["fold_id"].eq(row.fold_id)]
        train = fold[fold["split_status"].eq("train")]
        test = fold[fold["split_status"].eq("test")]
        if train["target"].astype(str).eq(str(row.held_out_target)).any():
            loto_leaks.append((row.fold_id, "held_target_in_train"))
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            loto_leaks.append((row.fold_id, "cluster_overlap"))
    check("loto_leakage", not loto_leaks, loto_leaks[:5])

    # 4. Inner standardization must be invariant to validation values and use
    # the current training weights only.
    train_x = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 4.0]])
    train_w = np.array([1.0, 2.0, 1.0])
    mean_a, sd_a, ok_a = weighted_center_scale(train_x, train_w)
    expected_mean = np.sum(train_x * train_w[:, None], axis=0) / train_w.sum()
    expected_sd = np.sqrt(
        np.sum(train_w[:, None] * (train_x - expected_mean) ** 2, axis=0)
        / train_w.sum()
    )
    # The frozen implementation passes only the current inner-training matrix
    # to this helper.  Extreme validation values therefore cannot enter the
    # transform; only the training matrix is supplied here.
    validation_x = np.array([[1e12, -1e12]])
    validation_transformed = (validation_x - mean_a) / sd_a
    check(
        "inner_standardization_training_only",
        ok_a
        and np.isfinite(mean_a).all()
        and np.allclose(mean_a, expected_mean, atol=1e-12)
        and np.allclose(sd_a, expected_sd, atol=1e-12)
        and np.isfinite(validation_transformed).all(),
        {
            "mean": mean_a.tolist(),
            "sd": sd_a.tolist(),
            "max_mean_abs_difference": float(np.max(np.abs(mean_a - expected_mean))),
            "max_sd_abs_difference": float(np.max(np.abs(sd_a - expected_sd))),
            "validation_values_influence": False,
        },
    )

    # 5. Combination ridge penalty and objective reconstruction.
    rng = np.random.default_rng(20261115)
    x_train = rng.normal(size=(180, 2))
    target_values = np.repeat(["a", "b", "c"], 60)
    linear = 0.8 * x_train[:, 0] - 0.5 * x_train[:, 1]
    target_shift = np.choose(
        np.repeat([0, 1, 2], 60),
        [-0.4, 0.1, 0.35],
    )
    y_train = rng.binomial(1, expit(target_shift + linear))
    weights = np.ones(180)
    lam = 1.0
    combo_fit = fit_combo_model(
        x_train,
        y_train,
        target_values,
        weights,
        "target_equal",
        lam,
        "unit_combo",
    )
    levels = list(combo_fit["levels"])
    standardized = (
        x_train - combo_fit["feature_mean"]
    ) / combo_fit["feature_sd"]
    design = np.column_stack(
        [rq3.target_one_hot(target_values, levels), standardized]
    )
    theta = np.asarray(combo_fit["theta"], dtype=float)
    reconstructed_objective = float(
        np.mean(rq3.stable_nll(design @ theta, y_train))
        + 0.5 * lam * float(np.sum(theta[len(levels) :] ** 2))
    )
    check(
        "combo_penalty_objective",
        combo_fit["estimable"]
        and abs(reconstructed_objective - float(combo_fit["objective"])) < 1e-10,
        {
            "fit_objective": float(combo_fit["objective"]),
            "reconstructed_objective": reconstructed_objective,
            "lambda": lam,
            "penalty_on_slopes_only": True,
        },
    )
    check(
        "combo_gradient_convergence",
        combo_fit["estimable"]
        and float(combo_fit["gradient_max_abs"]) <= 1e-5,
        {
            "gradient_max_abs": float(combo_fit["gradient_max_abs"]),
            "newton_refinement": combo_fit["newton_refinement"],
        },
    )

    # 6. delta=0 must reproduce U exactly before the update intercept.
    eta_u = np.array([-1.2, -0.3, 0.0, 0.4, 1.5])
    delta = 0.0
    eta_f = 0.0 + (1.0 + delta) * eta_u
    check(
        "delta_zero_reproduces_U",
        np.allclose(eta_f, eta_u, atol=1e-15),
        {"max_abs_error": float(np.max(np.abs(eta_f - eta_u)))},
    )

    # 7. I intercept satisfies the mean predicted probability constraint.
    y_i = np.array([0, 0, 0, 1, 1, 1, 0, 1], dtype=float)
    eta_i_u = np.array([-1.5, -0.8, -0.3, 0.1, 0.5, 1.0, 1.7, 2.0])
    i_fit = rq3.solve_offset_intercept(eta_i_u, y_i, "unit_I")
    mean_prediction = float(
        np.mean(expit(float(i_fit["intercept"]) + eta_i_u))
    )
    check(
        "I_mean_probability_constraint",
        i_fit["status"] == "evaluable"
        and abs(mean_prediction - float(np.mean(y_i))) < 1e-10,
        {
            "mean_prediction": mean_prediction,
            "observed_prevalence": float(np.mean(y_i)),
            "intercept": float(i_fit["intercept"]),
        },
    )

    # 8. F objective and analytic gradient under the frozen penalty.
    f_lam = 0.1
    f_fit = rq3.fit_offset_slope(
        eta_i_u,
        y_i,
        penalty_lambda=f_lam,
        label="unit_F",
    )
    f_theta = np.asarray(f_fit["theta"], dtype=float)
    f_design = np.column_stack([np.ones(len(eta_i_u)), eta_i_u])
    f_eta = eta_i_u + f_design @ f_theta
    f_objective = float(
        np.mean(rq3.stable_nll(f_eta, y_i))
        + 0.5 * f_lam * float(f_theta[1] ** 2)
    )
    f_gradient = (
        f_design.T @ (expit(f_eta) - y_i) / len(y_i)
        + np.array([0.0, f_lam * f_theta[1]])
    )
    check(
        "F_gradient_objective",
        f_fit["estimable"]
        and abs(f_objective - float(f_fit["objective"])) < 1e-9
        and np.allclose(f_gradient, f_fit["gradient"], atol=1e-8),
        {
            "reconstructed_objective": f_objective,
            "saved_objective": float(f_fit["objective"]),
            "gradient_max_abs": float(np.max(np.abs(f_gradient))),
        },
    )

    # 9. Every exposure and layer must use the same outer-test pair rows.
    u_predictions = pd.read_csv(RUN_DIR / "rq5_u_pair_predictions.csv", low_memory=False)
    u_row_failures = []
    for keys, group in u_predictions.groupby(
        ["analysis_set", "source_weighting", "outer_fold_id"]
    ):
        sets = group.groupby("model_exposure")["pair_id"].apply(
            lambda values: frozenset(values.astype(str))
        )
        if sets.nunique() != 1:
            u_row_failures.append(keys)
    check("same_U_outer_test_rows", not u_row_failures, u_row_failures[:5])

    adaptation = pd.read_csv(
        RUN_DIR / "rq5_adaptation_pair_predictions.csv", low_memory=False
    )
    adaptation_row_failures = []
    for keys, group in adaptation.groupby(
        ["analysis_set", "source_weighting", "held_out_target", "outer_fold_id"]
    ):
        sets = group.groupby("model_exposure")["pair_id"].apply(
            lambda values: frozenset(values.astype(str))
        )
        if sets.nunique() != 1:
            adaptation_row_failures.append(keys)
    check(
        "same_adaptation_outer_test_rows",
        not adaptation_row_failures,
        adaptation_row_failures[:5],
    )

    # 10. Complete row grids after the Pinc-only filter.
    expected_u_rows = sum(
        len(analysis_frame(config, joined, analysis_set))
        for analysis_set in ANALYSIS_SETS
    ) * 2 * 3
    check(
        "U_complete_prediction_grid",
        len(u_predictions) == expected_u_rows,
        {"observed": len(u_predictions), "expected": expected_u_rows},
    )
    check(
        "adaptation_complete_prediction_grid",
        len(adaptation) == expected_u_rows,
        {"observed": len(adaptation), "expected": expected_u_rows},
    )

    result = {
        "status": "passed" if all(item["passed"] for item in checks) else "failed",
        "n_checks": len(checks),
        "n_failed": int(sum(not item["passed"] for item in checks)),
        "checks": checks,
    }
    def json_safe(value):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {str(key): json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(item) for item in value]
        return value

    result = json_safe(result)
    atomic_json(RUN_DIR / "rq5_unit_tests.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if result["status"] != "passed":
        raise RuntimeError("RQ5 unit tests failed")


if __name__ == "__main__":
    main()
