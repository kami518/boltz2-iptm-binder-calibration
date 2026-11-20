"""Independent regression and correction audit for the RQ3 F parameterization."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit

from rq3_common import (
    DESCRIPTIVE_LAYER,
    LINEAR,
    PRED_F_LAYER,
    PRED_I_LAYER,
    RUN_DIR,
    atomic_csv,
    atomic_json,
    read_csv,
)

OLD_RUN = Path(
    r"<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ3_execution_2026-11-13"
    r"\run_20261113_175022"
)
IDENTITY_TOLERANCE = 1e-14
OBJECTIVE_TOLERANCE = 1e-10
ETA_RECONSTRUCTION_TOLERANCE = 1e-12
PROBABILITY_RECONSTRUCTION_TOLERANCE = 1e-12
GRADIENT_TOLERANCE = 1e-5
OLD_NEW_ETA_TOLERANCE = 1e-10
OLD_NEW_PROBABILITY_TOLERANCE = 1e-12

FIT_KEYS = [
    "analysis_set",
    "outcome_name",
    "cluster_rule",
    "held_out_target",
    "source_weighting",
]
PAIR_KEYS = [
    "pair_id",
    "analysis_set",
    "outcome_name",
    "cluster_rule",
    "held_out_target",
    "source_weighting",
]


def _max_abs(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = numeric[np.isfinite(numeric)]
    return float(np.max(np.abs(finite))) if finite.size else 0.0


def _objective(eta: np.ndarray, y: np.ndarray, penalty_lambda: float, delta: float) -> float:
    return float(
        np.mean(np.logaddexp(0.0, eta) - y * eta)
        + 0.5 * penalty_lambda * delta**2
    )


def _descriptive_consistency(
    predictions: pd.DataFrame,
    parameters: pd.DataFrame,
) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    selected_parameters = parameters[
        parameters["model_layer_label"].eq(DESCRIPTIVE_LAYER)
    ].copy()
    selected_predictions = predictions[
        predictions["descriptive_F_status"].eq("evaluable")
    ].drop(
        columns=[
            "descriptive_intercept",
            "descriptive_gamma",
            "descriptive_delta",
        ]
    ).copy()
    merged = selected_predictions.merge(
        selected_parameters,
        on=FIT_KEYS,
        how="inner",
        validate="many_to_one",
        suffixes=("", "_parameter"),
    )
    eta_reconstructed = (
        merged["descriptive_intercept"].to_numpy(dtype=float)
        + merged["descriptive_gamma"].to_numpy(dtype=float)
        * merged["eta_U"].to_numpy(dtype=float)
    )
    merged["eta_reconstructed"] = eta_reconstructed
    merged["p_reconstructed"] = expit(eta_reconstructed)
    merged["eta_reconstruction_abs_error"] = np.abs(
        merged["eta_descriptive_F"].to_numpy(dtype=float) - eta_reconstructed
    )
    merged["probability_reconstruction_abs_error"] = np.abs(
        merged["p_descriptive_F"].to_numpy(dtype=float)
        - merged["p_reconstructed"].to_numpy(dtype=float)
    )
    merged["gamma_identity_abs_error"] = np.abs(
        merged["descriptive_gamma"].to_numpy(dtype=float)
        - 1.0
        - merged["descriptive_delta"].to_numpy(dtype=float)
    )

    rows = []
    independent_frames = []
    for keys, group in merged.groupby(FIT_KEYS, sort=True):
        first = group.iloc[0]
        y = group["outcome_value"].to_numpy(dtype=float)
        eta_u = group["eta_U"].to_numpy(dtype=float)
        reconstructed = _objective(
            group["eta_reconstructed"].to_numpy(dtype=float),
            y,
            0.0,
            float(first["descriptive_delta"]),
        )

        def direct_objective(values: np.ndarray) -> float:
            linear_predictor = values[0] + values[1] * eta_u
            return float(
                np.mean(
                    np.logaddexp(0.0, linear_predictor)
                    - y * linear_predictor
                )
            )

        def direct_gradient(values: np.ndarray) -> np.ndarray:
            residual = expit(values[0] + values[1] * eta_u) - y
            return np.asarray(
                [np.mean(residual), np.mean(residual * eta_u)], dtype=float
            )

        prevalence = float(np.clip(np.mean(y), 1e-9, 1.0 - 1e-9))
        direct_fit = minimize(
            direct_objective,
            np.asarray([np.log(prevalence / (1.0 - prevalence)), 1.0]),
            method="L-BFGS-B",
            jac=direct_gradient,
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-10},
        )
        direct_gradient_max_abs = float(
            np.max(np.abs(direct_gradient(direct_fit.x)))
        )
        direct_predictions = expit(
            direct_fit.x[0] + direct_fit.x[1] * eta_u
        )
        stored_predictions = group["p_descriptive_F"].to_numpy(dtype=float)
        direct_error = (
            float(np.max(np.abs(stored_predictions - direct_predictions)))
            if np.isfinite(direct_fit.x).all()
            else np.nan
        )
        direct_success = bool(
            np.isfinite(direct_fit.x).all()
            and direct_gradient_max_abs <= 1e-6
            and float(direct_fit.fun)
            <= float(first["descriptive_objective"]) + 1e-10
        )
        independent_frames.append(
            pd.DataFrame(
                [
                    {
                        **dict(zip(FIT_KEYS, keys)),
                        "n": len(group),
                        "stored_intercept": float(
                            first["descriptive_intercept"]
                        ),
                        "independent_intercept": float(direct_fit.x[0]),
                        "stored_gamma": float(first["descriptive_gamma"]),
                        "independent_gamma": float(direct_fit.x[1]),
                        "stored_objective": float(
                            first["descriptive_objective"]
                        ),
                        "independent_objective": float(direct_fit.fun),
                        "independent_gradient_max_abs": (
                            direct_gradient_max_abs
                        ),
                        "max_probability_error": direct_error,
                        "independent_optimizer_success": bool(
                            direct_fit.success
                        ),
                        "independent_success": direct_success,
                    }
                ]
            )
        )
        rows.append(
            {
                "record_type": "descriptive_F",
                "fit_id": "|".join(str(value) for value in keys),
                **dict(zip(FIT_KEYS, keys)),
                "n_training_rows": int(len(group)),
                "saved_objective": float(first["descriptive_objective"]),
                "reconstructed_objective": reconstructed,
                "objective_abs_error": abs(
                    reconstructed - float(first["descriptive_objective"])
                ),
                "saved_gradient_max_abs": float(
                    first["descriptive_gradient_max_abs"]
                ),
                "gamma_identity_abs_error": float(
                    group["gamma_identity_abs_error"].max()
                ),
                "max_eta_reconstruction_abs_error": float(
                    group["eta_reconstruction_abs_error"].max()
                ),
                "max_probability_reconstruction_abs_error": float(
                    group["probability_reconstruction_abs_error"].max()
                ),
                "lambda_zero_independent_max_probability_error": direct_error,
            }
        )
    return pd.DataFrame(rows), independent_frames


def _predictive_consistency(
    u_predictions: pd.DataFrame,
    adaptation_predictions: pd.DataFrame,
    inner_assignments: pd.DataFrame,
    parameters: pd.DataFrame,
) -> pd.DataFrame:
    selected_parameters = parameters[
        parameters["model_layer_label"].eq(PRED_F_LAYER)
    ].copy()
    training_pairs = inner_assignments[
        [
            "outer_fold_id",
            *FIT_KEYS,
            "pair_id",
        ]
    ].drop_duplicates()
    u_linear = u_predictions[
        u_predictions["model_form"].eq(LINEAR)
    ][
        [
            "pair_id",
            *FIT_KEYS,
            "eta_U",
            "outcome_value",
        ]
    ].copy()
    training = training_pairs.merge(
        u_linear,
        on=["pair_id", *FIT_KEYS],
        how="inner",
        validate="many_to_one",
    )
    training = training.merge(
        selected_parameters[
            [
                "outer_fold_id",
                *FIT_KEYS,
                "intercept",
                "gamma",
                "delta",
                "selected_lambda",
                "objective",
                "gradient_max_abs",
                "n_adaptation_train",
            ]
        ],
        on=["outer_fold_id", *FIT_KEYS],
        how="inner",
        validate="many_to_one",
        suffixes=("", "_parameter"),
    )
    training["eta_reconstructed"] = (
        training["intercept"].to_numpy(dtype=float)
        + training["gamma"].to_numpy(dtype=float)
        * training["eta_U"].to_numpy(dtype=float)
    )
    training["gamma_identity_abs_error"] = np.abs(
        training["gamma"].to_numpy(dtype=float)
        - 1.0
        - training["delta"].to_numpy(dtype=float)
    )

    test_consistency = adaptation_predictions[
        [
            "outer_fold_id",
            *FIT_KEYS,
            "eta_U",
            "eta_F",
            "p_F",
        ]
    ].merge(
        selected_parameters[
            [
                "outer_fold_id",
                *FIT_KEYS,
                "intercept",
                "gamma",
                "delta",
            ]
        ],
        on=["outer_fold_id", *FIT_KEYS],
        how="inner",
        validate="many_to_one",
    )
    test_consistency["eta_reconstructed"] = (
        test_consistency["intercept"].to_numpy(dtype=float)
        + test_consistency["gamma"].to_numpy(dtype=float)
        * test_consistency["eta_U"].to_numpy(dtype=float)
    )
    test_consistency["eta_reconstruction_abs_error"] = np.abs(
        test_consistency["eta_F"].to_numpy(dtype=float)
        - test_consistency["eta_reconstructed"].to_numpy(dtype=float)
    )
    test_consistency["probability_reconstruction_abs_error"] = np.abs(
        test_consistency["p_F"].to_numpy(dtype=float)
        - expit(test_consistency["eta_reconstructed"].to_numpy(dtype=float))
    )
    test_errors = test_consistency.groupby(
        ["outer_fold_id", *FIT_KEYS]
    ).agg(
        max_eta_reconstruction_abs_error=(
            "eta_reconstruction_abs_error",
            "max",
        ),
        max_probability_reconstruction_abs_error=(
            "probability_reconstruction_abs_error",
            "max",
        ),
    )

    rows = []
    for keys, group in training.groupby(
        ["outer_fold_id", *FIT_KEYS], sort=True
    ):
        first = group.iloc[0]
        eta_f = group["eta_reconstructed"].to_numpy(dtype=float)
        y = group["outcome_value"].to_numpy(dtype=float)
        reconstructed = _objective(
            eta_f,
            y,
            float(first["selected_lambda"]),
            float(first["delta"]),
        )
        fit_errors = test_errors.loc[keys]
        rows.append(
            {
                "record_type": "predictive_F",
                "fit_id": "|".join(str(value) for value in keys),
                **dict(zip(["outer_fold_id", *FIT_KEYS], keys)),
                "n_training_rows": int(len(group)),
                "saved_objective": float(first["objective"]),
                "reconstructed_objective": reconstructed,
                "objective_abs_error": abs(
                    reconstructed - float(first["objective"])
                ),
                "saved_gradient_max_abs": float(first["gradient_max_abs"]),
                "gamma_identity_abs_error": float(
                    group["gamma_identity_abs_error"].max()
                ),
                "max_eta_reconstruction_abs_error": float(
                    fit_errors["max_eta_reconstruction_abs_error"]
                ),
                "max_probability_reconstruction_abs_error": float(
                    fit_errors[
                        "max_probability_reconstruction_abs_error"
                    ]
                ),
                "lambda_zero_independent_max_probability_error": np.nan,
            }
        )
    result = pd.DataFrame(rows)
    expected_counts = selected_parameters.set_index(
        ["outer_fold_id", *FIT_KEYS]
    )["n_adaptation_train"]
    observed_counts = training.groupby(
        ["outer_fold_id", *FIT_KEYS]
    ).size()
    if not expected_counts.reindex(observed_counts.index).eq(observed_counts).all():
        raise RuntimeError("Outer F training-row reconstruction does not match saved counts")
    return result


def _compare_old_new(
    new: pd.DataFrame,
    old: pd.DataFrame,
    keys: list[str],
    pairs: list[tuple[str, str]],
    layer: str,
) -> tuple[dict, pd.DataFrame]:
    new_frame = new[keys + [column for column, _ in pairs]].copy()
    old_frame = old[keys + [column for _, column in pairs]].copy()
    merged = new_frame.merge(
        old_frame,
        on=keys,
        how="outer",
        validate="one_to_one",
        suffixes=("_new", "_old"),
        indicator=True,
    )
    unmatched = int(merged["_merge"].ne("both").sum())
    change_rows = []
    for new_column, old_column in pairs:
        merged[f"{new_column}_abs_difference"] = np.abs(
            pd.to_numeric(merged[f"{new_column}_new"], errors="coerce")
            - pd.to_numeric(merged[f"{old_column}_old"], errors="coerce")
        )
    eta_pairs = [
        (new_column, old_column)
        for new_column, old_column in pairs
        if new_column.startswith("eta_")
    ]
    probability_pairs = [
        (new_column, old_column)
        for new_column, old_column in pairs
        if new_column.startswith("p_")
        or new_column == "predicted_probability"
    ]
    max_eta = max(
        (
            _max_abs(merged[f"{new_column}_abs_difference"])
            for new_column, _ in eta_pairs
        ),
        default=0.0,
    )
    max_probability = max(
        (
            _max_abs(merged[f"{new_column}_abs_difference"])
            for new_column, _ in probability_pairs
        ),
        default=0.0,
    )
    material_mask = pd.Series(False, index=merged.index)
    for new_column, _ in eta_pairs:
        material_mask |= (
            merged[f"{new_column}_abs_difference"]
            > OLD_NEW_ETA_TOLERANCE
        ).fillna(False)
    for new_column, _ in probability_pairs:
        material_mask |= (
            merged[f"{new_column}_abs_difference"]
            > OLD_NEW_PROBABILITY_TOLERANCE
        ).fillna(False)
    material = int(material_mask.sum() + unmatched)
    summary = {
        "layer": layer,
        "new_rows": int(len(new)),
        "old_rows": int(len(old)),
        "joined_rows": int(merged["_merge"].eq("both").sum()),
        "unmatched_rows": unmatched,
        "max_absolute_eta_difference": max_eta,
        "max_absolute_probability_difference": max_probability,
        "material_change_rows": material,
        "eta_materiality_tolerance": OLD_NEW_ETA_TOLERANCE,
        "probability_materiality_tolerance": OLD_NEW_PROBABILITY_TOLERANCE,
    }
    if material:
        changed = merged.loc[material_mask, keys + [
            column for column in merged.columns
            if column.endswith("_abs_difference")
        ]].copy()
        changed.insert(0, "layer", layer)
        change_rows.append(changed)
    return summary, (
        pd.concat(change_rows, ignore_index=True)
        if change_rows
        else pd.DataFrame()
    )


def run_regression(
    *,
    u_predictions: pd.DataFrame,
    descriptive_predictions: pd.DataFrame,
    descriptive_parameters: pd.DataFrame,
    adaptation_predictions: pd.DataFrame,
    inner_assignments: pd.DataFrame,
    adaptation_parameters: pd.DataFrame,
) -> dict:
    descriptive, independent_frames = _descriptive_consistency(
        descriptive_predictions,
        descriptive_parameters,
    )
    predictive = _predictive_consistency(
        u_predictions,
        adaptation_predictions,
        inner_assignments,
        adaptation_parameters,
    )
    consistency = pd.concat([descriptive, predictive], ignore_index=True)
    atomic_csv(RUN_DIR / "rq3_f_objective_consistency.csv", consistency)
    independent = pd.concat(independent_frames, ignore_index=True)
    atomic_csv(RUN_DIR / "rq3_descriptive_lambda_zero_independent.csv", independent)

    old_u = read_csv(OLD_RUN / "rq3_u_pair_predictions.csv")
    old_descriptive = read_csv(
        OLD_RUN / "rq3_descriptive_pair_predictions.csv"
    )
    old_adaptation = read_csv(
        OLD_RUN / "rq3_adaptation_pair_predictions.csv"
    )
    u_summary, u_changes = _compare_old_new(
        u_predictions,
        old_u,
        [
            "pair_id",
            "analysis_set",
            "outcome_name",
            "fold_id",
            "cluster_rule",
            "held_out_target",
            "source_weighting",
            "model_form",
        ],
        [("eta_U", "eta_U"), ("predicted_probability", "predicted_probability")],
        "U",
    )
    oracle_summary, oracle_changes = _compare_old_new(
        descriptive_predictions,
        old_descriptive,
        PAIR_KEYS,
        [("eta_oracle_I", "eta_oracle_I"), ("p_oracle_I", "p_oracle_I")],
        "oracle_I",
    )
    predictive_i_summary, predictive_i_changes = _compare_old_new(
        adaptation_predictions,
        old_adaptation,
        [
            "pair_id",
            "analysis_set",
            "outcome_name",
            "outer_fold_id",
            "cluster_rule",
            "held_out_target",
            "source_weighting",
        ],
        [("eta_I", "eta_I"), ("p_I", "p_I")],
        "predictive_I",
    )
    descriptive_f_summary, descriptive_f_changes = _compare_old_new(
        descriptive_predictions,
        old_descriptive,
        PAIR_KEYS,
        [
            ("eta_descriptive_F", "eta_descriptive_F"),
            ("p_descriptive_F", "p_descriptive_F"),
        ],
        "descriptive_F",
    )
    predictive_f_summary, predictive_f_changes = _compare_old_new(
        adaptation_predictions,
        old_adaptation,
        [
            "pair_id",
            "analysis_set",
            "outcome_name",
            "outer_fold_id",
            "cluster_rule",
            "held_out_target",
            "source_weighting",
        ],
        [("eta_F", "eta_F"), ("p_F", "p_F")],
        "predictive_F",
    )
    old_new_summary = pd.DataFrame(
        [
            u_summary,
            oracle_summary,
            predictive_i_summary,
            descriptive_f_summary,
            predictive_f_summary,
        ]
    )
    atomic_csv(RUN_DIR / "rq3_old_new_consistency.csv", old_new_summary)
    changes = [
        frame
        for frame in [
            u_changes,
            oracle_changes,
            predictive_i_changes,
            descriptive_f_changes,
            predictive_f_changes,
        ]
        if not frame.empty
    ]
    atomic_csv(
        RUN_DIR / "rq3_old_new_material_changes.csv",
        pd.concat(changes, ignore_index=True) if changes else pd.DataFrame(),
    )

    checks = {
        "objective_prediction_rows_complete": bool(
            len(consistency) == 88 + 440
        ),
        "objective_reconstruction_within_tolerance": bool(
            consistency["objective_abs_error"].max() <= OBJECTIVE_TOLERANCE
        ),
        "gamma_identity_within_tolerance": bool(
            consistency["gamma_identity_abs_error"].max()
            <= IDENTITY_TOLERANCE
        ),
        "eta_reconstruction_within_tolerance": bool(
            consistency["max_eta_reconstruction_abs_error"].max()
            <= ETA_RECONSTRUCTION_TOLERANCE
        ),
        "probability_reconstruction_within_tolerance": bool(
            consistency[
                "max_probability_reconstruction_abs_error"
            ].max()
            <= PROBABILITY_RECONSTRUCTION_TOLERANCE
        ),
        "saved_gradient_within_tolerance": bool(
            consistency["saved_gradient_max_abs"].max()
            <= GRADIENT_TOLERANCE
        ),
        "all_lambda_zero_independent_fits_succeeded": bool(
            independent["independent_success"].all()
        ),
        "all_lambda_zero_independent_predictions_match": bool(
            independent["max_probability_error"].max() <= 1e-7
        ),
        "U_unchanged_within_materiality": bool(
            u_summary["material_change_rows"] == 0
        ),
        "oracle_I_unchanged_within_materiality": bool(
            oracle_summary["material_change_rows"] == 0
        ),
        "predictive_I_unchanged_within_materiality": bool(
            predictive_i_summary["material_change_rows"] == 0
        ),
    }
    result = {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "tolerances": {
            "objective_abs": OBJECTIVE_TOLERANCE,
            "gamma_identity_abs": IDENTITY_TOLERANCE,
            "eta_reconstruction_abs": ETA_RECONSTRUCTION_TOLERANCE,
            "probability_reconstruction_abs": (
                PROBABILITY_RECONSTRUCTION_TOLERANCE
            ),
            "saved_gradient_max_abs": GRADIENT_TOLERANCE,
            "lambda_zero_independent_probability_abs": 1e-7,
            "old_new_eta_abs": OLD_NEW_ETA_TOLERANCE,
            "old_new_probability_abs": OLD_NEW_PROBABILITY_TOLERANCE,
        },
        "record_counts": {
            "descriptive_F_fits": int(
                descriptive["record_type"].eq("descriptive_F").sum()
            ),
            "predictive_F_fits": int(
                predictive["record_type"].eq("predictive_F").sum()
            ),
            "total_F_fits": int(len(consistency)),
        },
        "numerical_maxima": {
            "objective_reconstruction_abs_error": float(
                consistency["objective_abs_error"].max()
            ),
            "gamma_identity_abs_error": float(
                consistency["gamma_identity_abs_error"].max()
            ),
            "eta_reconstruction_abs_error": float(
                consistency["max_eta_reconstruction_abs_error"].max()
            ),
            "probability_reconstruction_abs_error": float(
                consistency[
                    "max_probability_reconstruction_abs_error"
                ].max()
            ),
            "saved_gradient_max_abs": float(
                consistency["saved_gradient_max_abs"].max()
            ),
            "lambda_zero_independent_probability_abs_error": float(
                independent["max_probability_error"].max()
            ),
        },
        "old_new_consistency": old_new_summary.to_dict(orient="records"),
        "invalidated_old_F_scope": {
            "descriptive_F_fits": int(
                (
                    read_csv(OLD_RUN / "rq3_descriptive_parameters.csv")[
                        "model_layer_label"
                    ].eq(DESCRIPTIVE_LAYER)
                ).sum()
            ),
            "predictive_F_fits": int(
                (
                    read_csv(OLD_RUN / "rq3_adaptation_parameters.csv")[
                        "model_layer_label"
                    ].eq(PRED_F_LAYER)
                ).sum()
            ),
            "lambda_candidates": int(
                len(read_csv(OLD_RUN / "rq3_lambda_selection.csv"))
            ),
            "descriptive_F_pair_predictions": int(len(old_descriptive)),
            "predictive_F_pair_predictions": int(len(old_adaptation)),
            "reason": (
                "The old fit optimized a + theta_1*eta_U while saving "
                "delta=theta_1 and predicting with a+(1+delta)*eta_U."
            ),
        },
    }
    atomic_json(RUN_DIR / "rq3_f_regression_tests.json", result)
    return result
