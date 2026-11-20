"""RQ3 stage 2: oracle base-rate I and in-sample descriptive F diagnostics."""
from __future__ import annotations

import math
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from scipy.special import expit

from rq3_common import (
    DESCRIPTIVE_LAYER,
    LINEAR,
    ORACLE_LAYER,
    RUN_DIR,
    append_failure,
    atomic_csv,
    atomic_json,
    compute_performance_metrics,
    fit_offset_slope,
    load_config,
    read_csv,
    reliability_rows,
    require_hash_gate,
    set_stage_status,
    solve_offset_intercept,
)

STAGE = "stage_2_oracle_I_and_descriptive_F"
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


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    require_hash_gate(config)
    u_predictions = read_csv(RUN_DIR / "rq3_u_pair_predictions.csv")
    u_predictions = u_predictions[
        u_predictions["model_form"].eq(LINEAR)
        & u_predictions["fit_status"].eq("evaluable")
    ].copy()

    wide_frames: list[pd.DataFrame] = []
    metric_frames: list[pd.DataFrame] = []
    parameter_frames: list[pd.DataFrame] = []
    reliability_frames: list[pd.DataFrame] = []
    expected_units = 0

    group_columns = [
        "analysis_set",
        "cluster_rule",
        "source_weighting",
        "held_out_target",
        "fold_id",
    ]
    for keys, frame in u_predictions.groupby(group_columns, sort=True):
        analysis_set, cluster_rule, source_weighting, held_target, fold_id = keys
        expected_units += 1
        frame = frame.sort_values("pair_id").reset_index(drop=True)
        y = frame["outcome_value"].to_numpy(dtype=float)
        eta_u = frame["eta_U"].to_numpy(dtype=float)
        p_u = frame["predicted_probability"].to_numpy(dtype=float)
        clip = float(config["probability_clip"])
        transport_ref = frame["bss_transport_reference"].to_numpy(dtype=float)
        oracle_ref = frame["bss_target_oracle_reference"].to_numpy(dtype=float)

        oracle = solve_offset_intercept(
            eta_u,
            y,
            f"{fold_id}|{source_weighting}|oracle_I",
        )
        if oracle["status"] == "evaluable":
            oracle_intercept = float(oracle["intercept"])
            eta_oracle = oracle_intercept + eta_u
            p_oracle = expit(eta_oracle)
            oracle_status = "evaluable"
            oracle_failure = ""
        else:
            oracle_intercept = np.nan
            eta_oracle = np.full(len(frame), np.nan)
            p_oracle = np.full(len(frame), np.nan)
            oracle_status = "failed"
            oracle_failure = oracle["failure_type"]
            append_failure(
                STAGE,
                f"{fold_id}|{source_weighting}|oracle_I",
                oracle["failure_type"],
                oracle["message"],
            )

        descriptive = fit_offset_slope(
            eta_u,
            y,
            penalty_lambda=0.0,
            label=f"{fold_id}|{source_weighting}|descriptive_F",
        )
        if descriptive["estimable"]:
            eta_descriptive = (
                descriptive["intercept"] + descriptive["gamma"] * eta_u
            )
            p_descriptive = expit(eta_descriptive)
            descriptive_status = "evaluable"
            descriptive_failure = ""
        else:
            eta_descriptive = np.full(len(frame), np.nan)
            p_descriptive = np.full(len(frame), np.nan)
            descriptive_status = "failed"
            descriptive_failure = descriptive["failure_type"]
            append_failure(
                STAGE,
                f"{fold_id}|{source_weighting}|descriptive_F",
                descriptive["failure_type"],
                descriptive["message"],
            )

        wide = frame[
            [
                "pair_id",
                "analysis_set",
                "outcome_name",
                "fold_id",
                "cluster_rule",
                "held_out_target",
                "target",
                "sequence_cluster_value",
                "outcome_value",
                "boltz2_iptm",
                "x",
                "source_weighting",
                "eta_U",
                "predicted_probability",
                "bss_transport_reference",
                "bss_target_oracle_reference",
                "transport_intercept",
            ]
        ].rename(columns={"predicted_probability": "p_U"})
        wide["model_form"] = LINEAR
        wide["eta_oracle_I"] = eta_oracle
        wide["p_oracle_I"] = p_oracle
        wide["oracle_intercept"] = oracle_intercept
        wide["oracle_I_status"] = oracle_status
        wide["eta_descriptive_F"] = eta_descriptive
        wide["p_descriptive_F"] = p_descriptive
        wide["descriptive_intercept"] = descriptive.get("intercept", np.nan)
        wide["descriptive_gamma"] = descriptive.get("gamma", np.nan)
        wide["descriptive_delta"] = descriptive.get("delta", np.nan)
        wide["descriptive_F_status"] = descriptive_status
        wide_frames.append(wide)

        base = {
            "fold_id": fold_id,
            "analysis_set": analysis_set,
            "outcome_name": str(frame["outcome_name"].iloc[0]),
            "cluster_rule": cluster_rule,
            "held_out_target": held_target,
            "source_weighting": source_weighting,
            "model_form": LINEAR,
            "n": int(len(frame)),
            "positives": int(np.sum(y)),
            "bss_transport_reference_value": float(transport_ref[0]),
            "bss_target_oracle_reference_value": float(oracle_ref[0]),
        }
        for layer, prediction, eta, status, failure, parameter_summary in [
            (
                ORACLE_LAYER,
                p_oracle,
                eta_oracle,
                oracle_status,
                oracle_failure,
                {
                    "oracle_intercept": oracle_intercept,
                    "oracle_residual": oracle.get("residual", np.nan),
                },
            ),
            (
                DESCRIPTIVE_LAYER,
                p_descriptive,
                eta_descriptive,
                descriptive_status,
                descriptive_failure,
                {
                    "descriptive_intercept": descriptive.get("intercept", np.nan),
                    "descriptive_gamma": descriptive.get("gamma", np.nan),
                    "descriptive_delta": descriptive.get("delta", np.nan),
                    "descriptive_objective": descriptive.get(
                        "objective", np.nan
                    ),
                    "descriptive_gradient_max_abs": descriptive.get(
                        "gradient_max_abs", np.nan
                    ),
                    "descriptive_iterations": descriptive.get(
                        "iterations", np.nan
                    ),
                },
            ),
        ]:
            parameter_frames.append(
                pd.DataFrame(
                    [
                        {
                            "fold_id": fold_id,
                            "analysis_set": analysis_set,
                            "outcome_name": str(frame["outcome_name"].iloc[0]),
                            "cluster_rule": cluster_rule,
                            "held_out_target": held_target,
                            "source_weighting": source_weighting,
                            "model_form": LINEAR,
                            "model_layer_label": layer,
                            "fit_status": status,
                            "failure_type": failure,
                            "n_target_pairs": int(len(frame)),
                            "n_target_positives": int(np.sum(y)),
                            **parameter_summary,
                        }
                    ]
                )
            )
            if status == "evaluable":
                performance = compute_performance_metrics(
                    y,
                    prediction,
                    clip=clip,
                    transport_reference=transport_ref,
                    target_oracle_reference=oracle_ref,
                )
                if performance["calibration_status"] != "evaluable":
                    append_failure(
                        STAGE,
                        f"{fold_id}|{source_weighting}|{layer}|calibration",
                        "calibration_failure",
                        performance["calibration_message"],
                    )
            else:
                performance = {
                    "n": int(len(frame)),
                    "positives": int(np.sum(y)),
                    "prevalence": float(np.mean(y)),
                    "brier": np.nan,
                    "log_loss": np.nan,
                    "bss_transport": np.nan,
                    "bss_target_oracle": np.nan,
                    "calibration_intercept": np.nan,
                    "calibration_slope": np.nan,
                    "calibration_intercept_se": np.nan,
                    "calibration_slope_se": np.nan,
                    "calibration_status": "failed",
                    "calibration_failure_type": "adaptation_failure",
                    "calibration_message": failure,
                }
            metric_frames.append(
                pd.DataFrame(
                    [
                        {
                            **base,
                            "model_layer_label": layer,
                            "convergence_flags": status,
                            **performance,
                        }
                    ]
                )
            )

            reliability_frame = pd.DataFrame(
                {
                    "outcome_value": y,
                    "predicted_probability": prediction,
                }
            )
            for reliability in reliability_rows(
                reliability_frame,
                "predicted_probability",
                config["reliability_bin_edges"],
            ):
                reliability_frames.append(
                    pd.DataFrame(
                        [
                            {
                                "fold_id": fold_id,
                                "analysis_set": analysis_set,
                                "outcome_name": str(frame["outcome_name"].iloc[0]),
                                "cluster_rule": cluster_rule,
                                "held_out_target": held_target,
                                "source_weighting": source_weighting,
                                "model_form": LINEAR,
                                "model_layer_label": layer,
                                "n_training_pairs": np.nan,
                                "n_source_targets": np.nan,
                                **reliability,
                            }
                        ]
                    )
                )

    wide_predictions = pd.concat(wide_frames, ignore_index=True)
    metrics = pd.concat(metric_frames, ignore_index=True)
    parameters = pd.concat(parameter_frames, ignore_index=True)
    new_reliability = pd.concat(reliability_frames, ignore_index=True)
    expected_units = 44 * 2
    if len(metrics) != expected_units * 2:
        append_failure(
            STAGE,
            "stage_2_grid",
            "incomplete_model_grid",
            "Oracle/descriptive model grid is incomplete.",
            expected=expected_units * 2,
            observed=int(len(metrics)),
        )
    if len(wide_predictions) != 32280:
        append_failure(
            STAGE,
            "stage_2_prediction_count",
            "incomplete_pair_predictions",
            "Oracle/descriptive pair prediction count is incorrect.",
            expected=32280,
            observed=int(len(wide_predictions)),
        )

    atomic_csv(RUN_DIR / "rq3_descriptive_pair_predictions.csv", wide_predictions)
    atomic_csv(RUN_DIR / "rq3_descriptive_metrics.csv", metrics)
    atomic_csv(RUN_DIR / "rq3_descriptive_parameters.csv", parameters)

    existing_reliability = read_csv(RUN_DIR / "rq3_loto_reliability.csv")
    combined_reliability = pd.concat(
        [existing_reliability, new_reliability], ignore_index=True
    )
    atomic_csv(RUN_DIR / "rq3_loto_reliability.csv", combined_reliability)

    macro_rows = []
    for keys, group in metrics.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
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
                    "model_form": LINEAR,
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
    atomic_csv(RUN_DIR / "rq3_descriptive_macro_summary.csv", macro)

    pair_rows = []
    for layer, probability_column in [
        (ORACLE_LAYER, "p_oracle_I"),
        (DESCRIPTIVE_LAYER, "p_descriptive_F"),
    ]:
        for keys, group in wide_predictions.groupby(
            ["analysis_set", "cluster_rule", "source_weighting"], sort=True
        ):
            base = dict(
                zip(
                    ["analysis_set", "cluster_rule", "source_weighting"], keys
                )
            )
            if group[probability_column].isna().any():
                performance = {
                    "n": int(len(group)),
                    "positives": int(group["outcome_value"].sum()),
                    "brier": np.nan,
                    "log_loss": np.nan,
                    "bss_transport": np.nan,
                    "bss_target_oracle": np.nan,
                }
            else:
                performance = compute_performance_metrics(
                    group["outcome_value"].to_numpy(dtype=float),
                    group[probability_column].to_numpy(dtype=float),
                    clip=float(config["probability_clip"]),
                    transport_reference=group[
                        "bss_transport_reference"
                    ].to_numpy(dtype=float),
                    target_oracle_reference=group[
                        "bss_target_oracle_reference"
                    ].to_numpy(dtype=float),
                )
            pair_rows.append(
                {
                    **base,
                    "model_form": LINEAR,
                    "model_layer_label": layer,
                    "summary_type": "pair_pooled_across_held_out_targets",
                    **performance,
                }
            )
    pair_summary = pd.DataFrame(pair_rows)
    atomic_csv(RUN_DIR / "rq3_descriptive_pair_summary.csv", pair_summary)

    evaluable = metrics["convergence_flags"].eq("evaluable")
    stage_summary = {
        "stage": STAGE,
        "expected_model_evaluations": expected_units * 2,
        "observed_model_evaluations": int(len(metrics)),
        "evaluable_model_evaluations": int(evaluable.sum()),
        "failed_model_evaluations": int((~evaluable).sum()),
        "pair_rows": int(len(wide_predictions)),
        "oracle_mean_intercept": float(
            parameters.loc[
                parameters["model_layer_label"].eq(ORACLE_LAYER),
                "oracle_intercept",
            ].mean()
        ),
        "descriptive_mean_gamma": float(
            parameters.loc[
                parameters["model_layer_label"].eq(DESCRIPTIVE_LAYER),
                "descriptive_gamma",
            ].mean()
        ),
        "labels": {
            "oracle_I": ORACLE_LAYER,
            "descriptive_F": DESCRIPTIVE_LAYER,
        },
    }
    atomic_json(RUN_DIR / "stage2_summary.json", stage_summary)
    set_stage_status(
        STAGE,
        "complete",
        expected_models=expected_units * 2,
        evaluable_models=int(evaluable.sum()),
        failed_models=int((~evaluable).sum()),
    )
    print(
        f"Stage 2 complete: {len(metrics)} diagnostics, "
        f"{int(evaluable.sum())} evaluable."
    )


if __name__ == "__main__":
    main()
