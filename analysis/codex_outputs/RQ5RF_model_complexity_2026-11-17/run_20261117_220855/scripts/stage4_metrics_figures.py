"""RQ5-RF Stage 4: metrics, complete-sample gates, pairwise results, figures."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rq5rf_common import (
    ANALYSIS_SETS,
    EXPOSURES,
    PRESET_COMPARISONS,
    RUN_DIR,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_json,
    bss,
    brier,
    load_config,
    load_joined_inputs,
    log_loss,
    pair_set_sha256,
    performance_metrics,
    reliability_rows,
    save_figure,
    set_stage_status,
)

STAGE = "stage_4_metrics_and_figures"
LAYERS = (
    "U_transport_zero_label",
    "predictive_I",
    "predictive_F",
)
METRIC_COLUMNS = (
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "roc_auc",
    "ap",
    "calibration_intercept",
    "calibration_slope",
)
PAIR_METRICS = (
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "roc_auc",
    "ap",
    "calibration_intercept",
    "calibration_slope",
)


def layer_probability_column(layer: str) -> str:
    return {
        "U_transport_zero_label": "p_U",
        "predictive_I": "p_I",
        "predictive_F": "p_F",
    }[layer]


def metric_direction(metric: str) -> str:
    if metric in {"brier", "log_loss"}:
        return "negative_is_left_improvement"
    if metric in {"bss_transport", "bss_target_oracle", "bss_adaptation_train"}:
        return "positive_is_left_improvement"
    if metric in {"roc_auc", "ap"}:
        return "positive_is_left_improvement"
    return "descriptive_no_better_direction"


def is_improvement(metric: str, difference: float) -> bool:
    if not np.isfinite(difference):
        return False
    if metric in {"brier", "log_loss"}:
        return difference < -1e-12
    if metric in {
        "bss_transport",
        "bss_target_oracle",
        "bss_adaptation_train",
        "roc_auc",
        "ap",
    }:
        return difference > 1e-12
    return False


def is_reversal(metric: str, difference: float) -> bool:
    if not np.isfinite(difference):
        return False
    if metric in {"brier", "log_loss"}:
        return difference > 1e-12
    if metric in {
        "bss_transport",
        "bss_target_oracle",
        "bss_adaptation_train",
        "roc_auc",
        "ap",
    }:
        return difference < -1e-12
    return False


def target_metrics(predictions: pd.DataFrame, clip: float) -> pd.DataFrame:
    rows = []
    groups = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "target",
    ]
    for keys, group in predictions.groupby(groups, sort=True):
        analysis_set, outcome_name, cluster_rule, weighting, exposure, target = keys
        group = group.sort_values("pair_id", kind="mergesort").reset_index(drop=True)
        for layer in LAYERS:
            probability_column = layer_probability_column(layer)
            probability = group[probability_column].to_numpy(dtype=float)
            metrics = performance_metrics(
                group["outcome_value"].to_numpy(dtype=float),
                probability,
                clip=clip,
                transport_reference=group[
                    "bss_transport_reference"
                ].to_numpy(dtype=float),
                target_oracle_reference=group[
                    "bss_target_oracle_reference"
                ].to_numpy(dtype=float),
                adaptation_reference=(
                    group["bss_adaptation_train_reference"].to_numpy(dtype=float)
                    if layer != "U_transport_zero_label"
                    else None
                ),
            )
            status_column = {
                "U_transport_zero_label": None,
                "predictive_I": "I_status",
                "predictive_F": "F_status",
            }[layer]
            fit_status = (
                "evaluable"
                if status_column is None
                else str(group[status_column].iloc[0])
            )
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "model_layer_label": layer,
                    "target": target,
                    "target_outer_fold_count": int(
                        group["outer_fold_id"].nunique()
                    ),
                    "pair_id_set_sha256": pair_set_sha256(
                        group["pair_id"].astype(str)
                    ),
                    "fit_status": fit_status,
                    "failure_type": (
                        ""
                        if layer == "U_transport_zero_label"
                        else str(
                            group[
                                "I_failure_type"
                                if layer == "predictive_I"
                                else "F_failure_type"
                            ].iloc[0]
                        )
                    ),
                    "selected_lambda": (
                        float(group["selected_lambda"].iloc[0])
                        if layer == "predictive_F"
                        else np.nan
                    ),
                    **metrics,
                }
            )
    return pd.DataFrame(rows).sort_values(
        [
            "analysis_set",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
            "target",
        ],
        kind="mergesort",
    )


def expected_targets_for_set(analysis_set: str, config: dict) -> set[str]:
    # Frozen complete-case target membership, independent of observed rows.
    if analysis_set == "primary_five":
        return {
            "human-vegf110",
            "nipah-glycoprotein-g",
            "pd-l1",
            "sars-cov-2-rbd",
            "trem2",
        }
    return {
        "aah-i",
        "human-vegf110",
        "nipah-glycoprotein-g",
        "pd-l1",
        "sars-cov-2-rbd",
        "trem2",
    }


def macro_metrics(target: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows = []
    grouping = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "model_layer_label",
    ]
    for keys, group in target.groupby(grouping, sort=True):
        analysis_set, outcome_name, cluster_rule, weighting, exposure, layer = keys
        expected = expected_targets_for_set(analysis_set, config)
        observed = set(group["target"].astype(str))
        all_expected = observed == expected
        row = {
            "analysis_set": analysis_set,
            "outcome_name": outcome_name,
            "cluster_rule": cluster_rule,
            "source_weighting": weighting,
            "model_exposure": exposure,
            "model_layer_label": layer,
            "summary_type": "target_equal_macro",
            "expected_target_count": len(expected),
            "observed_target_count": len(observed),
            "targets_observed_json": json.dumps(sorted(observed)),
            "missing_targets_json": json.dumps(sorted(expected - observed)),
        }
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group[metric], errors="coerce")
            evaluable = bool(
                all_expected
                and values.notna().all()
                and np.isfinite(values.to_numpy(dtype=float)).all()
            )
            if metric == "bss_adaptation_train" and layer == "U_transport_zero_label":
                evaluable = False
            row[metric] = float(values.mean()) if evaluable else np.nan
            row[f"{metric}_status"] = (
                "evaluable"
                if evaluable
                else (
                    "not_applicable_to_U"
                    if metric == "bss_adaptation_train"
                    and layer == "U_transport_zero_label"
                    else "not_evaluable_missing_expected_target_or_metric"
                )
            )
        row["macro_status"] = (
            "evaluable"
            if all_expected
            and group["fit_status"].eq("evaluable").all()
            else "not_evaluable"
        )
        rows.append(row)
    return pd.DataFrame(rows)


def pair_pooled_metrics(
    predictions: pd.DataFrame, config: dict, joined: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    grouping = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
    ]
    for keys, group in predictions.groupby(grouping, sort=True):
        analysis_set, outcome_name, cluster_rule, weighting, exposure = keys
        valid = analysis_frame(config, joined, analysis_set)
        expected_ids = set(valid["pair_id"].astype(str))
        duplicated = int(group["pair_id"].astype(str).duplicated().sum())
        for layer in LAYERS:
            probability_column = layer_probability_column(layer)
            probability = group[probability_column].to_numpy(dtype=float)
            complete = bool(
                set(group["pair_id"].astype(str)) == expected_ids
                and not duplicated
                and np.isfinite(probability).all()
            )
            row = {
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "cluster_rule": cluster_rule,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "model_layer_label": layer,
                "summary_type": "pair_pooled",
                "n_expected_pairs": len(expected_ids),
                "n_observed_pairs": len(group),
                "n_unique_observed_pairs": group["pair_id"].astype(str).nunique(),
                "duplicate_pair_rows": duplicated,
                "pair_set_sha256": pair_set_sha256(group["pair_id"].astype(str)),
                "complete_sample_gate": "pass" if complete else "fail",
            }
            if complete:
                metrics = performance_metrics(
                    group["outcome_value"].to_numpy(dtype=float),
                    probability,
                    clip=float(config["probability_clip"]),
                    transport_reference=group[
                        "bss_transport_reference"
                    ].to_numpy(dtype=float),
                    target_oracle_reference=group[
                        "bss_target_oracle_reference"
                    ].to_numpy(dtype=float),
                    adaptation_reference=(
                        group["bss_adaptation_train_reference"].to_numpy(dtype=float)
                        if layer != "U_transport_zero_label"
                        else None
                    ),
                )
                row.update(metrics)
            else:
                row.update({metric: np.nan for metric in METRIC_COLUMNS})
            rows.append(row)
    return pd.DataFrame(rows)


def pairwise_differences(target: pd.DataFrame, macro: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_rows = []
    macro_rows = []
    config_keys = ["analysis_set", "outcome_name", "cluster_rule", "source_weighting"]
    for left_model, right_model in PRESET_COMPARISONS:
        comparison = f"{left_model}_minus_{right_model}"
        for layer in LAYERS:
            left = target[
                target["model_exposure"].eq(left_model)
                & target["model_layer_label"].eq(layer)
            ]
            right = target[
                target["model_exposure"].eq(right_model)
                & target["model_layer_label"].eq(layer)
            ]
            aligned = left.merge(
                right,
                on=config_keys + ["target"],
                how="outer",
                validate="one_to_one",
                suffixes=("__left", "__right"),
                indicator=True,
            )
            for _, row in aligned.iterrows():
                base = {
                    "analysis_set": row["analysis_set"],
                    "outcome_name": row["outcome_name"],
                    "cluster_rule": row["cluster_rule"],
                    "source_weighting": row["source_weighting"],
                    "comparison": comparison,
                    "left_model": left_model,
                    "right_model": right_model,
                    "model_layer_label": layer,
                    "target": row["target"],
                    "alignment_status": "aligned" if row["_merge"] == "both" else "missing_side",
                }
                for metric in PAIR_METRICS:
                    left_value = row[f"{metric}__left"]
                    right_value = row[f"{metric}__right"]
                    difference = (
                        float(left_value) - float(right_value)
                        if pd.notna(left_value) and pd.notna(right_value)
                        else np.nan
                    )
                    target_rows.append(
                        {
                            **base,
                            "metric": metric,
                            "left_value": left_value,
                            "right_value": right_value,
                            "difference_left_minus_right": difference,
                            "absolute_difference": abs(difference)
                            if np.isfinite(difference)
                            else np.nan,
                            "direction_convention": metric_direction(metric),
                            "left_improved": is_improvement(metric, difference),
                            "left_reversed": is_reversal(metric, difference),
                        }
                    )
            for keys, group in pd.DataFrame(target_rows).query(
                "comparison == @comparison and model_layer_label == @layer"
            ).groupby(config_keys, sort=True):
                analysis_set, outcome_name, cluster_rule, weighting = keys
                expected = expected_targets_for_set(analysis_set, pd.DataFrame())
                observed = set(group["target"].astype(str))
                for metric in PAIR_METRICS:
                    metric_group = group[group["metric"].eq(metric)]
                    values = pd.to_numeric(
                        metric_group["difference_left_minus_right"], errors="coerce"
                    )
                    complete = (
                        observed == expected
                        and len(metric_group) == len(expected)
                        and values.notna().all()
                        and np.isfinite(values.to_numpy(dtype=float)).all()
                        and metric_group["alignment_status"].eq("aligned").all()
                    )
                    macro_rows.append(
                        {
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": cluster_rule,
                            "source_weighting": weighting,
                            "comparison": comparison,
                            "left_model": left_model,
                            "right_model": right_model,
                            "model_layer_label": layer,
                            "metric": metric,
                            "summary_type": "target_equal_macro",
                            "mean_difference_left_minus_right": (
                                float(values.mean()) if complete else np.nan
                            ),
                            "absolute_mean_difference": (
                                abs(float(values.mean())) if complete else np.nan
                            ),
                            "min_difference": (
                                float(values.min()) if complete else np.nan
                            ),
                            "max_difference": (
                                float(values.max()) if complete else np.nan
                            ),
                            "improvement_target_count": int(
                                sum(is_improvement(metric, value) for value in values)
                            )
                            if complete
                            else 0,
                            "reversal_target_count": int(
                                sum(is_reversal(metric, value) for value in values)
                            )
                            if complete
                            else 0,
                            "target_count_expected": len(expected),
                            "target_count_observed": len(metric_group),
                            "complete_sample_gate": "pass" if complete else "fail",
                            "direction_convention": metric_direction(metric),
                        }
                    )
    return pd.DataFrame(target_rows), pd.DataFrame(macro_rows)


def pairwise_pair_pooled(pair_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for left_model, right_model in PRESET_COMPARISONS:
        comparison = f"{left_model}_minus_{right_model}"
        for layer in LAYERS:
            for keys, all_group in pair_metrics[
                pair_metrics["model_layer_label"].eq(layer)
            ].groupby(
                ["analysis_set", "outcome_name", "cluster_rule", "source_weighting"],
                sort=True,
            ):
                analysis_set, outcome_name, cluster_rule, weighting = keys
                left = all_group[all_group["model_exposure"].eq(left_model)]
                right = all_group[all_group["model_exposure"].eq(right_model)]
                if len(left) != 1 or len(right) != 1:
                    continue
                left = left.iloc[0]
                right = right.iloc[0]
                complete = (
                    left["complete_sample_gate"] == "pass"
                    and right["complete_sample_gate"] == "pass"
                    and left["pair_set_sha256"] == right["pair_set_sha256"]
                )
                row = {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "comparison": comparison,
                    "left_model": left_model,
                    "right_model": right_model,
                    "model_layer_label": layer,
                    "summary_type": "pair_pooled",
                    "complete_sample_gate": "pass" if complete else "fail",
                    "left_pair_set_sha256": left["pair_set_sha256"],
                    "right_pair_set_sha256": right["pair_set_sha256"],
                }
                for metric in PAIR_METRICS:
                    difference = (
                        float(left[metric]) - float(right[metric])
                        if complete and pd.notna(left[metric]) and pd.notna(right[metric])
                        else np.nan
                    )
                    row[f"left_{metric}"] = left[metric]
                    row[f"right_{metric}"] = right[metric]
                    row[f"difference_{metric}"] = difference
                    row[f"absolute_difference_{metric}"] = (
                        abs(difference) if np.isfinite(difference) else np.nan
                    )
                    row[f"direction_{metric}"] = metric_direction(metric)
                rows.append(row)
    return pd.DataFrame(rows)


def reliability(predictions: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows = []
    grouping = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "target",
    ]
    for keys, group in predictions.groupby(grouping, sort=True):
        analysis_set, outcome_name, cluster_rule, weighting, exposure, target = keys
        for layer in LAYERS:
            probability_column = layer_probability_column(layer)
            for record in reliability_rows(
                group["outcome_value"].to_numpy(dtype=int),
                group[probability_column].to_numpy(dtype=float),
                config["reliability_bin_edges"],
            ):
                rows.append(
                    {
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": cluster_rule,
                        "source_weighting": weighting,
                        "model_exposure": exposure,
                        "model_layer_label": layer,
                        "target": target,
                        **record,
                    }
                )
    return pd.DataFrame(rows)


def screening_10_2(
    pairwise_macro: pd.DataFrame, pairwise_target: pd.DataFrame
) -> dict:
    comparisons = (
        "rf_combo_minus_pooled_logistic_combo",
        "rf_combo_minus_rf_ipTM",
    )
    metrics = ("brier", "log_loss", "bss_transport")

    def macro_value(analysis_set: str, weighting: str, comparison: str, metric: str) -> float:
        rows = pairwise_macro[
            pairwise_macro["analysis_set"].eq(analysis_set)
            & pairwise_macro["source_weighting"].eq(weighting)
            & pairwise_macro["comparison"].eq(comparison)
            & pairwise_macro["model_layer_label"].eq("U_transport_zero_label")
            & pairwise_macro["metric"].eq(metric)
        ]
        return float(rows.iloc[0]["mean_difference_left_minus_right"]) if len(rows) == 1 else np.nan

    condition_1_values = {}
    condition_1 = True
    for comparison in comparisons:
        for metric in metrics:
            value = macro_value("primary_six", "target_equal", comparison, metric)
            condition_1_values[f"{comparison}|{metric}"] = value
            condition_1 &= is_improvement(metric, value)

    condition_2_values = {}
    condition_2 = True
    for analysis_set, required in (("primary_six", 5), ("primary_five", 4)):
        for weighting in ("target_equal", "pair"):
            for comparison in comparisons:
                rows = pairwise_target[
                    pairwise_target["analysis_set"].eq(analysis_set)
                    & pairwise_target["source_weighting"].eq(weighting)
                    & pairwise_target["comparison"].eq(comparison)
                    & pairwise_target["model_layer_label"].eq("U_transport_zero_label")
                    & pairwise_target["metric"].eq("brier")
                ]
                improved = int(rows["left_improved"].fillna(False).sum())
                expected = len(expected_targets_for_set(analysis_set, {}))
                key = f"{analysis_set}|{weighting}|{comparison}|brier_improved_targets"
                condition_2_values[key] = improved
                condition_2 &= bool(len(rows) == expected and improved >= required)

    condition_3_values = {}
    condition_3 = True
    for analysis_set in ANALYSIS_SETS:
        for weighting in ("target_equal", "pair"):
            for comparison in comparisons:
                for metric in metrics:
                    value = macro_value(analysis_set, weighting, comparison, metric)
                    key = f"{analysis_set}|{weighting}|{comparison}|{metric}"
                    condition_3_values[key] = value
                    condition_3 &= bool(np.isfinite(value) and not is_reversal(metric, value))
    return {
        "screening_name": "10.2 exploratory direction-consistency screen",
        "is_statistical_test": False,
        "is_effect_importance_test": False,
        "neutral_tolerance": 1e-12,
        "condition_1_primary_six_target_equal_strict_improvement": bool(condition_1),
        "condition_1_values": condition_1_values,
        "condition_2_target_direction_counts": bool(condition_2),
        "condition_2_values": condition_2_values,
        "condition_3_six_configurations_no_reversal": bool(condition_3),
        "condition_3_values": condition_3_values,
        "condition_4_qc_provisional": False,
        "all_conditions_provisional_before_qc": bool(
            condition_1 and condition_2 and condition_3
        ),
        "allowed_conclusion_if_condition_4_also_passes": (
            "Direction-consistent exploratory increment; repeat in independent complete-case data."
        ),
    }


def make_figures(
    target: pd.DataFrame,
    pairwise_macro: pd.DataFrame,
    reliability_frame: pd.DataFrame,
) -> list[str]:
    figures = []
    primary_target = target[
        target["analysis_set"].eq("primary_six")
        & target["source_weighting"].eq("target_equal")
    ].copy()
    targets = sorted(
        primary_target.loc[
            primary_target["model_layer_label"].eq("U_transport_zero_label"),
            "target",
        ].unique()
    )
    x = np.arange(len(targets))
    width = 0.19
    colors = {
        "pooled_logistic_ipTM": "#2C6E9B",
        "pooled_logistic_combo": "#57A773",
        "rf_ipTM": "#D66B3D",
        "rf_combo": "#8C5AA8",
    }
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.3), constrained_layout=True)
    for axis, metric in zip(axes, ("brier", "bss_transport")):
        for offset, exposure in zip(np.linspace(-1.5 * width, 1.5 * width, 4), EXPOSURES):
            subset = primary_target[
                primary_target["model_exposure"].eq(exposure)
                & primary_target["model_layer_label"].eq("U_transport_zero_label")
            ].set_index("target")
            values = [float(subset.loc[target_name, metric]) if target_name in subset.index else np.nan for target_name in targets]
            axis.bar(x + offset, values, width=width, label=exposure, color=colors[exposure])
        axis.set_xticks(x)
        axis.set_xticklabels(targets, rotation=28, ha="right")
        axis.axhline(0, color="black", linewidth=0.8)
        axis.grid(axis="y", alpha=0.2)
        axis.set_title(metric.replace("_", " ").upper())
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("U-layer performance by target: primary_six, target_equal", fontsize=14)
    path = RUN_DIR / "RQ5RF_U_four_model_transport.png"
    save_figure(path, fig)
    figures.append(path.name)

    macro = pairwise_macro[
        pairwise_macro["model_layer_label"].eq("U_transport_zero_label")
        & pairwise_macro["metric"].isin(["brier", "log_loss", "bss_transport"])
    ].copy()
    labels = []
    values = []
    for comparison in [f"{left}_minus_{right}" for left, right in PRESET_COMPARISONS]:
        subset = macro[
            macro["comparison"].eq(comparison)
            & macro["analysis_set"].eq("primary_six")
            & macro["source_weighting"].eq("target_equal")
            & macro["metric"].eq("brier")
        ]
        labels.append(comparison.replace("_minus_", "\nminus\n"))
        values.append(float(subset.iloc[0]["mean_difference_left_minus_right"]) if len(subset) else np.nan)
    fig, axis = plt.subplots(figsize=(11, 5.2), constrained_layout=True)
    axis.bar(np.arange(len(labels)), values, color=["#386CB0", "#4D8B63", "#C8553D", "#8C5AA8", "#7A6F5B"])
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xticks(np.arange(len(labels)))
    axis.set_xticklabels(labels, fontsize=8)
    axis.set_ylabel("Brier difference (left minus right; negative favors left)")
    axis.set_title("Primary U-layer pairwise Brier differences: primary_six, target_equal")
    axis.grid(axis="y", alpha=0.2)
    path = RUN_DIR / "RQ5RF_primary_pairwise_differences.png"
    save_figure(path, fig)
    figures.append(path.name)

    fig, axes = plt.subplots(1, 4, figsize=(18, 5.3), sharey=True, constrained_layout=True)
    for axis, exposure in zip(axes, EXPOSURES):
        subset = primary_target[primary_target["model_exposure"].eq(exposure)]
        pivot = subset.pivot(index="target", columns="model_layer_label", values="brier")
        order = list(targets)
        y = np.arange(len(order))
        axis.scatter(
            pivot.loc[order, "predictive_I"] - pivot.loc[order, "U_transport_zero_label"],
            y,
            color="#6B7FD7",
            label="I - U",
            s=55,
        )
        axis.scatter(
            pivot.loc[order, "predictive_F"] - pivot.loc[order, "U_transport_zero_label"],
            y,
            color="#C8553D",
            label="F - U",
            s=55,
        )
        axis.axvline(0, color="black", linewidth=0.8)
        axis.set_yticks(y)
        axis.set_yticklabels(order)
        axis.set_xlabel("Brier difference (negative is better)")
        axis.set_title(exposure, fontsize=10)
        axis.grid(axis="x", alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle("I/F same-layer Brier comparison: primary_six, target_equal", fontsize=14)
    path = RUN_DIR / "RQ5RF_IF_same_layer_comparison.png"
    save_figure(path, fig)
    figures.append(path.name)

    rel = reliability_frame[
        reliability_frame["analysis_set"].eq("primary_six")
        & reliability_frame["source_weighting"].eq("target_equal")
    ].copy()
    fig, axes = plt.subplots(4, 3, figsize=(13, 14), constrained_layout=True)
    layer_order = LAYERS
    for row_index, exposure in enumerate(EXPOSURES):
        for column_index, layer in enumerate(layer_order):
            axis = axes[row_index, column_index]
            subset = rel[
                rel["model_exposure"].eq(exposure)
                & rel["model_layer_label"].eq(layer)
                & rel["n"].gt(0)
            ].copy()
            if subset.empty:
                axis.text(0.5, 0.5, "No evaluable bins", ha="center", va="center")
            else:
                subset["weighted_mean"] = (
                    subset["mean_prediction"] * subset["n"]
                )
                pooled_bins = (
                    subset.groupby(["bin_index", "bin_lower", "bin_upper"], sort=True)
                    .agg(
                        n=("n", "sum"),
                        positives=("positives", "sum"),
                        weighted_mean=("weighted_mean", "sum"),
                    )
                    .reset_index()
                )
                pooled_bins["mean_prediction"] = (
                    pooled_bins["weighted_mean"] / pooled_bins["n"]
                )
                pooled_bins["observed_rate"] = (
                    pooled_bins["positives"] / pooled_bins["n"]
                )
                axis.plot(
                    pooled_bins["mean_prediction"],
                    pooled_bins["observed_rate"],
                    marker="o",
                    color=colors[exposure],
                    linewidth=1.2,
                    markersize=4,
                )
            axis.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=0.8)
            axis.set_xlim(0, 1)
            axis.set_ylim(0, 1)
            axis.grid(alpha=0.2)
            axis.set_title(f"{exposure} | {layer}", fontsize=8)
            axis.set_xlabel("Mean prediction")
            axis.set_ylabel("Observed rate")
    fig.suptitle(
        "Reliability, six targets pooled (all evaluable bins retained): primary_six, target_equal",
        fontsize=14,
    )
    path = RUN_DIR / "RQ5RF_reliability_primary_six.png"
    save_figure(path, fig)
    figures.append(path.name)
    return figures


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    clip = float(config["probability_clip"])
    predictions = pd.read_csv(RUN_DIR / "rq5rf_pair_predictions.csv", low_memory=False)
    _, _, _, joined = load_joined_inputs(config)
    target = target_metrics(predictions, clip)
    for _, row in target[
        target["calibration_status"].ne("evaluable")
    ].iterrows():
        append_failure(
            STAGE,
            "|".join(
                [
                    str(row["analysis_set"]),
                    str(row["source_weighting"]),
                    str(row["model_exposure"]),
                    str(row["model_layer_label"]),
                    str(row["target"]),
                ]
            ),
            "calibration_non_evaluable",
            str(row.get("calibration_message", "Calibration diagnostic not evaluable")),
            resolution_status="recorded_non_evaluable_diagnostic",
            failure_type_detail=str(row.get("calibration_failure_type", "")),
        )
    macro = macro_metrics(target, config)
    pair_metrics = pair_pooled_metrics(predictions, config, joined)
    pairwise_target, pairwise_macro = pairwise_differences(target, macro)
    pairwise_pair = pairwise_pair_pooled(pair_metrics)
    reliability_frame = reliability(predictions, config)
    screening = screening_10_2(pairwise_macro, pairwise_target)
    figures = make_figures(target, pairwise_macro, reliability_frame)

    atomic_csv(RUN_DIR / "rq5rf_target_metrics.csv", target)
    atomic_csv(RUN_DIR / "rq5rf_macro_metrics.csv", macro)
    atomic_csv(RUN_DIR / "rq5rf_pair_metrics.csv", pair_metrics)
    atomic_csv(RUN_DIR / "rq5rf_pairwise_target_differences.csv", pairwise_target)
    atomic_csv(RUN_DIR / "rq5rf_pairwise_macro_comparisons.csv", pairwise_macro)
    atomic_csv(RUN_DIR / "rq5rf_pairwise_pair_pooled.csv", pairwise_pair)
    atomic_csv(RUN_DIR / "rq5rf_reliability.csv", reliability_frame)
    atomic_json(RUN_DIR / "direction_consistency_10_2.json", screening)

    summary = {
        "status": "complete",
        "target_metric_rows": len(target),
        "macro_metric_rows": len(macro),
        "pair_metric_rows": len(pair_metrics),
        "pairwise_target_difference_rows": len(pairwise_target),
        "pairwise_macro_rows": len(pairwise_macro),
        "pairwise_pair_pooled_rows": len(pairwise_pair),
        "reliability_rows": len(reliability_frame),
        "macro_all_expected_targets": bool(
            macro["macro_status"].eq("evaluable").all()
        ),
        "pair_pooled_all_complete": bool(
            pair_metrics["complete_sample_gate"].eq("pass").all()
        ),
        "pairwise_pair_sets_all_match": bool(
            pairwise_pair["complete_sample_gate"].eq("pass").all()
        ),
        "figures": figures,
        "screening_10_2_conditions_1_to_3": screening,
    }
    atomic_json(RUN_DIR / "stage4_metrics_summary.json", summary)
    set_stage_status(STAGE, "complete", details=summary)
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
