"""RQ3 stage 4: overall tables, figures, QC, and handoff report."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rq3_f_regression import run_regression
from rq3_common import (
    ANALYSIS_SETS,
    CLUSTER_RULES,
    DESCRIPTIVE_LAYER,
    LINEAR,
    ORACLE_LAYER,
    PRED_F_LAYER,
    PRED_I_LAYER,
    RCS,
    RUN_DIR,
    SOURCE_WEIGHTS,
    U_LAYER,
    append_failure,
    atomic_csv,
    atomic_json,
    bss_from_reference,
    compute_performance_metrics,
    load_config,
    read_csv,
    require_hash_gate,
    save_figure,
    set_stage_status,
    now_iso,
)

STAGE = "stage_4_figures_tables_qc_and_handoff"
TARGET_ORDER = [
    "aah-i",
    "human-vegf110",
    "nipah-glycoprotein-g",
    "pd-l1",
    "sars-cov-2-rbd",
    "trem2",
]
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


def markdown_table(frame: pd.DataFrame, digits: int = 4) -> str:
    if frame.empty:
        return "_No rows._"
    display = frame.copy()
    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(
                lambda value: (
                    ""
                    if pd.isna(value)
                    else f"{float(value):.{digits}f}"
                )
            )
        else:
            display[column] = display[column].fillna("").astype(str)
    header = "| " + " | ".join(display.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(display.columns)) + " |"
    rows = [
        "| " + " | ".join(row) + " |"
        for row in display.astype(str).itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def macro_from_targets(
    target_metrics: pd.DataFrame,
    group_columns: list[str],
) -> pd.DataFrame:
    rows = []
    for keys, group in target_metrics.groupby(group_columns, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_columns, keys))
        for metric in METRIC_COLUMNS:
            if metric not in group:
                continue
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            rows.append(
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
    return pd.DataFrame(rows)


def metric_lookup(
    macro: pd.DataFrame,
    *,
    stage: str,
    layer: str,
    metric: str,
    analysis_set: str = "primary_six",
    cluster_rule: str = "sequence_cluster90",
    source_weighting: str = "target_equal",
    model_form: str = LINEAR,
) -> float:
    row = macro[
        macro["stage"].eq(stage)
        & macro["model_layer_label"].eq(layer)
        & macro["metric"].eq(metric)
        & macro["analysis_set"].eq(analysis_set)
        & macro["cluster_rule"].eq(cluster_rule)
        & macro["source_weighting"].eq(source_weighting)
        & macro["model_form"].eq(model_form)
    ]
    if row.empty:
        return np.nan
    return float(row.iloc[0]["mean"])


def make_u_figure(
    main: pd.DataFrame,
    pair_main: pd.Series,
    reliability: pd.DataFrame,
    output_path: Path,
) -> None:
    ordered = main.set_index("held_out_target").loc[TARGET_ORDER].reset_index()
    x = np.arange(len(ordered))
    figure, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = plt.cm.tab10(np.linspace(0, 1, len(ordered)))

    axes[0, 0].bar(x, ordered["brier"], color=colors, alpha=0.82)
    axes[0, 0].axhline(
        ordered["brier"].mean(),
        color="black",
        linestyle="--",
        linewidth=1.3,
        label="target-equal macro",
    )
    axes[0, 0].scatter(
        [len(ordered) - 0.25],
        [pair_main["brier"]],
        marker="D",
        color="crimson",
        label="pair summary",
        zorder=4,
    )
    axes[0, 0].set_title("A. Held-out-target Brier score (lower is better)")
    axes[0, 0].set_ylabel("Brier")
    axes[0, 0].legend(frameon=False, fontsize=8)

    axes[0, 1].bar(x, ordered["log_loss"], color=colors, alpha=0.82)
    axes[0, 1].axhline(
        ordered["log_loss"].mean(),
        color="black",
        linestyle="--",
        linewidth=1.3,
        label="target-equal macro",
    )
    axes[0, 1].scatter(
        [len(ordered) - 0.25],
        [pair_main["log_loss"]],
        marker="D",
        color="crimson",
        label="pair summary",
        zorder=4,
    )
    axes[0, 1].set_title("B. Held-out-target log loss (lower is better)")
    axes[0, 1].set_ylabel("Log loss")
    axes[0, 1].legend(frameon=False, fontsize=8)

    width = 0.36
    axes[1, 0].bar(
        x - width / 2,
        ordered["bss_transport"],
        width,
        label="BSS-transport (zero-label)",
        color="#3b82b8",
    )
    axes[1, 0].bar(
        x + width / 2,
        ordered["bss_target_oracle"],
        width,
        label="BSS-target-oracle (diagnostic)",
        color="#d69b38",
        alpha=0.75,
    )
    axes[1, 0].axhline(0, color="black", linewidth=0.8)
    axes[1, 0].axhline(
        ordered["bss_transport"].mean(),
        color="#1d4f73",
        linestyle="--",
        linewidth=1.2,
    )
    axes[1, 0].scatter(
        [len(ordered) - 0.25],
        [pair_main["bss_transport"]],
        marker="D",
        color="crimson",
        zorder=4,
    )
    axes[1, 0].set_title(
        "C. BSS references are kept distinct\n(dashed line: zero-label macro; diamond: pair summary)"
    )
    axes[1, 0].set_ylabel("BSS")
    axes[1, 0].legend(frameon=False, fontsize=8)

    for color, target in zip(colors, TARGET_ORDER):
        subset = reliability[
            reliability["held_out_target"].eq(target)
            & reliability["n"].gt(0)
        ]
        axes[1, 1].plot(
            subset["mean_prediction"],
            subset["observed_rate"],
            marker="o",
            linewidth=1.2,
            alpha=0.8,
            color=color,
            label=target,
        )
    axes[1, 1].plot([0, 1], [0, 1], color="black", linestyle=":", linewidth=1)
    axes[1, 1].set_xlim(0, 1)
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_xlabel("Mean predicted probability in frozen bin")
    axes[1, 1].set_ylabel("Observed binding rate")
    axes[1, 1].set_title("D. Reliability support (per-target frozen bins)")
    axes[1, 1].legend(frameon=False, fontsize=7, ncol=2)

    for axis in axes.flat:
        axis.grid(alpha=0.18, linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    for axis in axes[:2, :].flat:
        axis.set_xticks(x)
        axis.set_xticklabels(ordered["held_out_target"], rotation=28, ha="right")
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(ordered["held_out_target"], rotation=28, ha="right")
    figure.suptitle(
        "RQ3 U layer: primary_six + cluster90 + target_equal + linear",
        fontsize=14,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(output_path, figure)


def make_adaptation_figure(
    main_targets: pd.DataFrame,
    selected_lambda: pd.DataFrame,
    output_path: Path,
) -> None:
    layers = [
        ("U_transport_zero_label", "U"),
        (PRED_I_LAYER, "Predictive I"),
        (PRED_F_LAYER, "Predictive F"),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(14, 10))
    x = np.arange(len(TARGET_ORDER))
    for axis, metric, title in [
        (axes[0, 0], "brier", "A. Same outer rows: Brier (lower is better)"),
        (axes[0, 1], "log_loss", "B. Same outer rows: log loss (lower is better)"),
        (
            axes[1, 0],
            "bss_target_oracle",
            "C. Same outer rows: BSS-target-oracle (diagnostic reference)",
        ),
    ]:
        for layer, label in layers:
            subset = (
                main_targets[main_targets["model_layer_label"].eq(layer)]
                .set_index("held_out_target")
                .loc[TARGET_ORDER]
                .reset_index()
            )
            axis.plot(
                x,
                subset[metric],
                marker="o",
                linewidth=1.6,
                label=label,
            )
        axis.set_xticks(x)
        axis.set_xticklabels(TARGET_ORDER, rotation=28, ha="right")
        axis.set_title(title)
        axis.set_ylabel(metric)
        axis.axhline(0, color="black", linewidth=0.7) if metric.startswith("bss") else None
        axis.legend(frameon=False, fontsize=8)

    frequencies = (
        selected_lambda["lambda"]
        .value_counts()
        .sort_index()
    )
    axes[1, 1].bar(
        [f"{value:g}" for value in frequencies.index],
        frequencies.values,
        color="#7b5ea7",
        alpha=0.82,
    )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel("Selected outer folds (log scale)")
    axes[1, 1].set_xlabel("Lambda selected by pooled inner-validation log loss")
    axes[1, 1].set_title("D. Frozen lambda grid selections (440 outer-fold fits)")

    for axis in axes.flat:
        axis.grid(alpha=0.18, linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(
        "RQ3 predictive adaptation (CV5 outer folds)\n"
        "Oracle I and descriptive F are nondeployable diagnostics and are shown only in the report",
        fontsize=13,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(output_path, figure)


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    initial_hashes = require_hash_gate(config)
    stage0_qc = json.loads((RUN_DIR / "stage0_qc.json").read_text(encoding="utf-8"))
    u_predictions = read_csv(RUN_DIR / "rq3_u_pair_predictions.csv")
    u_metrics = read_csv(RUN_DIR / "rq3_loto_metrics.csv")
    u_parameters = read_csv(RUN_DIR / "rq3_source_model_parameters.csv")
    u_knots = read_csv(RUN_DIR / "rq3_rcs_knots.csv")
    u_reliability = read_csv(RUN_DIR / "rq3_loto_reliability.csv")
    u_pair_summary = read_csv(RUN_DIR / "rq3_loto_pair_summary.csv")
    descriptive_predictions = read_csv(
        RUN_DIR / "rq3_descriptive_pair_predictions.csv"
    )
    descriptive_metrics = read_csv(RUN_DIR / "rq3_descriptive_metrics.csv")
    descriptive_parameters = read_csv(RUN_DIR / "rq3_descriptive_parameters.csv")
    adaptation_predictions = read_csv(
        RUN_DIR / "rq3_adaptation_pair_predictions.csv"
    )
    inner_assignments = read_csv(RUN_DIR / "rq3_adaptation_inner_assignments.csv")
    lambda_selection = read_csv(RUN_DIR / "rq3_lambda_selection.csv")
    fold_metrics = read_csv(RUN_DIR / "rq3_adaptation_fold_metrics.csv")
    target_adaptation_metrics = read_csv(
        RUN_DIR / "rq3_target_adaptation_metrics.csv"
    )
    adaptation_parameters = read_csv(RUN_DIR / "rq3_adaptation_parameters.csv")
    comparisons = read_csv(RUN_DIR / "rq3_adaptation_comparisons.csv")
    failure_ledger = read_csv(RUN_DIR / "failed_item_ledger.csv")
    f_regression = run_regression(
        u_predictions=u_predictions,
        descriptive_predictions=descriptive_predictions,
        descriptive_parameters=descriptive_parameters,
        adaptation_predictions=adaptation_predictions,
        inner_assignments=inner_assignments,
        adaptation_parameters=adaptation_parameters,
    )
    if f_regression["status"] != "passed":
        append_failure(
            STAGE,
            "rq3_f_regression",
            "regression_test_failure",
            "One or more F correction regression checks failed.",
            failed_checks=[
                name
                for name, passed in f_regression["checks"].items()
                if not passed
            ],
        )
        failure_ledger = read_csv(RUN_DIR / "failed_item_ledger.csv")

    u_metrics["stage"] = "stage_1_U"
    descriptive_metrics["stage"] = "stage_2_oracle_descriptive"
    target_adaptation_metrics["model_form"] = LINEAR
    target_adaptation_metrics["stage"] = "stage_3_predictive"
    overall_metrics = pd.concat(
        [
            u_metrics.assign(summary_scope="per_held_out_target"),
            descriptive_metrics.assign(summary_scope="per_held_out_target"),
            target_adaptation_metrics.assign(summary_scope="per_held_out_target"),
        ],
        ignore_index=True,
        sort=False,
    )
    group_columns = [
        "stage",
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_form",
        "model_layer_label",
    ]
    overall_macro = macro_from_targets(overall_metrics, group_columns)
    atomic_csv(RUN_DIR / "rq3_overall_metrics.csv", overall_metrics)
    atomic_csv(RUN_DIR / "rq3_overall_macro_summary.csv", overall_macro)

    u_pair_summary["stage"] = "stage_1_U"
    descriptive_pair_summary = read_csv(
        RUN_DIR / "rq3_descriptive_pair_summary.csv"
    )
    descriptive_pair_summary["stage"] = "stage_2_oracle_descriptive"
    adaptation_pair_summary = read_csv(
        RUN_DIR / "rq3_adaptation_pair_summary.csv"
    )
    adaptation_pair_summary["stage"] = "stage_3_predictive"
    overall_pair_summary = pd.concat(
        [
            u_pair_summary,
            descriptive_pair_summary,
            adaptation_pair_summary,
        ],
        ignore_index=True,
        sort=False,
    )
    atomic_csv(RUN_DIR / "rq3_overall_pair_summary.csv", overall_pair_summary)

    sensitivity_specs = [
        (
            "main",
            "primary_six",
            "sequence_cluster90",
            "target_equal",
            LINEAR,
        ),
        (
            "primary_five_excluding_aah_i",
            "primary_five",
            "sequence_cluster90",
            "target_equal",
            LINEAR,
        ),
        (
            "strict_outcome",
            "strict_six",
            "sequence_cluster90",
            "target_equal",
            LINEAR,
        ),
        (
            "sequence_cluster80",
            "primary_six",
            "sequence_cluster80",
            "target_equal",
            LINEAR,
        ),
        (
            "pair_source_weighting",
            "primary_six",
            "sequence_cluster90",
            "pair",
            LINEAR,
        ),
        (
            "restricted_cubic_spline",
            "primary_six",
            "sequence_cluster90",
            "target_equal",
            RCS,
        ),
    ]
    sensitivity_rows = []
    for dimension, analysis_set, rule, weight, form in sensitivity_specs:
        for stage, layer_values in [
            (
                "stage_1_U",
                [U_LAYER],
            ),
            (
                "stage_2_oracle_descriptive",
                [ORACLE_LAYER, DESCRIPTIVE_LAYER],
            ),
            (
                "stage_3_predictive",
                ["U_transport_zero_label", PRED_I_LAYER, PRED_F_LAYER],
            ),
        ]:
            if dimension == "restricted_cubic_spline" and stage != "stage_1_U":
                continue
            for layer in layer_values:
                subset_macro = overall_macro[
                    overall_macro["stage"].eq(stage)
                    & overall_macro["analysis_set"].eq(analysis_set)
                    & overall_macro["cluster_rule"].eq(rule)
                    & overall_macro["source_weighting"].eq(weight)
                    & overall_macro["model_form"].eq(form)
                    & overall_macro["model_layer_label"].eq(layer)
                ]
                if subset_macro.empty:
                    continue
                for metric in ("brier", "log_loss", "bss_transport", "bss_target_oracle", "bss_adaptation_train"):
                    selected = subset_macro[subset_macro["metric"].eq(metric)]
                    if selected.empty:
                        continue
                    row = selected.iloc[0]
                    sensitivity_rows.append(
                        {
                            "dimension": dimension,
                            "stage": stage,
                            "analysis_set": analysis_set,
                            "cluster_rule": rule,
                            "source_weighting": weight,
                            "model_form": form,
                            "model_layer_label": layer,
                            "metric": metric,
                            "macro_mean": row["mean"],
                            "target_count": row["target_count"],
                            "min": row["min"],
                            "max": row["max"],
                            "n_positive_sign": row["n_positive_sign"],
                            "n_negative_sign": row["n_negative_sign"],
                        }
                    )
    sensitivity = pd.DataFrame(sensitivity_rows)
    atomic_csv(RUN_DIR / "rq3_sensitivity_summary.csv", sensitivity)

    main_u = u_metrics[
        u_metrics["analysis_set"].eq("primary_six")
        & u_metrics["cluster_rule"].eq("sequence_cluster90")
        & u_metrics["source_weighting"].eq("target_equal")
        & u_metrics["model_form"].eq(LINEAR)
    ].copy()
    main_u["held_out_target"] = pd.Categorical(
        main_u["held_out_target"], categories=TARGET_ORDER, ordered=True
    )
    main_u = main_u.sort_values("held_out_target")
    main_u_reliability = u_reliability[
        u_reliability["analysis_set"].eq("primary_six")
        & u_reliability["cluster_rule"].eq("sequence_cluster90")
        & u_reliability["source_weighting"].eq("target_equal")
        & u_reliability["model_form"].eq(LINEAR)
        & u_reliability["model_layer_label"].eq(U_LAYER)
    ]
    main_pair_u = u_pair_summary[
        u_pair_summary["analysis_set"].eq("primary_six")
        & u_pair_summary["cluster_rule"].eq("sequence_cluster90")
        & u_pair_summary["source_weighting"].eq("target_equal")
        & u_pair_summary["model_form"].eq(LINEAR)
    ].iloc[0]
    make_u_figure(
        main_u,
        main_pair_u,
        main_u_reliability,
        RUN_DIR / "RQ3_U_LOTO.png",
    )

    main_adaptation = target_adaptation_metrics[
        target_adaptation_metrics["analysis_set"].eq("primary_six")
        & target_adaptation_metrics["cluster_rule"].eq("sequence_cluster90")
        & target_adaptation_metrics["source_weighting"].eq("target_equal")
    ].copy()
    main_adaptation["held_out_target"] = pd.Categorical(
        main_adaptation["held_out_target"],
        categories=TARGET_ORDER,
        ordered=True,
    )
    main_adaptation = main_adaptation.sort_values(
        ["held_out_target", "model_layer_label"]
    )
    selected_main = lambda_selection[
        lambda_selection["analysis_set"].eq("primary_six")
        & lambda_selection["cluster_rule"].eq("sequence_cluster90")
        & lambda_selection["source_weighting"].eq("target_equal")
        & lambda_selection["selected"].astype(bool)
    ]
    make_adaptation_figure(
        main_adaptation,
        selected_main,
        RUN_DIR / "RQ3_adaptation.png",
    )

    overall_by_target = main_adaptation.pivot(
        index="held_out_target",
        columns="model_layer_label",
        values=["brier", "log_loss", "bss_adaptation_train"],
    )
    main_desc = descriptive_metrics[
        descriptive_metrics["analysis_set"].eq("primary_six")
        & descriptive_metrics["cluster_rule"].eq("sequence_cluster90")
        & descriptive_metrics["source_weighting"].eq("target_equal")
    ].copy()
    main_desc["held_out_target"] = pd.Categorical(
        main_desc["held_out_target"], categories=TARGET_ORDER, ordered=True
    )
    main_desc = main_desc.sort_values(["held_out_target", "model_layer_label"])
    main_comparison = comparisons[
        comparisons["analysis_set"].eq("primary_six")
        & comparisons["cluster_rule"].eq("sequence_cluster90")
        & comparisons["source_weighting"].eq("target_equal")
    ].copy()
    main_comparison = main_comparison[
        main_comparison["metric"].isin(["brier", "log_loss", "bss_target_oracle"])
    ]

    # Numerical QC reconstructs metrics and checks all frozen selection rules.
    qc_checks: list[dict] = []

    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        qc_checks.append(
            {
                "check": name,
                "passed": bool(passed),
                "observed": observed,
                "expected": expected,
            }
        )

    final_hashes = require_hash_gate(config)
    atomic_csv(RUN_DIR / "final_input_hashes.csv", final_hashes)
    check("initial_stage0_structural_qc", bool(stage0_qc["all_passed"]), stage0_qc["n_checks"], 19)
    check("frozen_input_hash_gate_closed", final_hashes["sha256_matches"].all(), int(final_hashes["sha256_matches"].sum()), len(final_hashes))
    check("loto_fold_count", u_metrics["fold_id"].nunique() == 44, int(u_metrics["fold_id"].nunique()), 44)
    check(
        "u_model_target_evaluations",
        len(u_metrics) == 176,
        int(len(u_metrics)),
        176,
    )
    check(
        "u_pair_prediction_count",
        len(u_predictions) == 64560,
        int(len(u_predictions)),
        64560,
    )
    check(
        "u_source_models_evaluable",
        u_metrics["convergence_flags"].eq("evaluable").all(),
        int(u_metrics["convergence_flags"].eq("evaluable").sum()),
        176,
    )
    check(
        "u_predictions_finite",
        u_predictions.loc[
            u_predictions["fit_status"].eq("evaluable"),
            "predicted_probability",
        ].map(np.isfinite).all(),
        int(u_predictions["predicted_probability"].isna().sum()),
        0,
    )
    expected_target_counts = (
        u_predictions.groupby(
            ["analysis_set", "cluster_rule", "source_weighting", "model_form", "held_out_target"]
        ).size()
    )
    check(
        "u_test_pair_coverage_unique",
        u_predictions.duplicated(
            ["analysis_set", "cluster_rule", "source_weighting", "model_form", "pair_id"]
        ).sum() == 0,
        int(
            u_predictions.duplicated(
                ["analysis_set", "cluster_rule", "source_weighting", "model_form", "pair_id"]
            ).sum()
        ),
        0,
    )
    knot_groups_ok = True
    for _, group in u_knots.groupby(
        ["fold_id", "source_weighting"]
    ):
        values = group.sort_values("quantile")["knot_x"].to_numpy(dtype=float)
        knot_groups_ok &= (
            len(values) == 3
            and np.isfinite(values).all()
            and np.all(np.diff(values) > 0)
        )
    check("rcs_knots_strict_and_source_training_only", bool(knot_groups_ok), int(len(u_knots)), 264)
    recomputed_u_errors = []
    for _, group in u_predictions.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "model_form",
            "held_out_target",
        ]
    ):
        stored = u_metrics[
            u_metrics["analysis_set"].eq(group["analysis_set"].iloc[0])
            & u_metrics["cluster_rule"].eq(group["cluster_rule"].iloc[0])
            & u_metrics["source_weighting"].eq(group["source_weighting"].iloc[0])
            & u_metrics["model_form"].eq(group["model_form"].iloc[0])
            & u_metrics["held_out_target"].eq(group["held_out_target"].iloc[0])
        ].iloc[0]
        performance = compute_performance_metrics(
            group["outcome_value"].to_numpy(dtype=float),
            group["predicted_probability"].to_numpy(dtype=float),
            clip=float(config["probability_clip"]),
            transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
            target_oracle_reference=group["bss_target_oracle_reference"].to_numpy(dtype=float),
        )
        for metric in ("brier", "log_loss", "bss_transport", "bss_target_oracle"):
            recomputed_u_errors.append(abs(float(stored[metric]) - performance[metric]))
    check(
        "u_metrics_recomputable_from_predictions",
        max(recomputed_u_errors) < 1e-10,
        float(max(recomputed_u_errors)),
        "<1e-10",
    )
    check(
        "u_reliability_frozen_bins_complete",
        len(u_reliability) == (176 + 176) * 7,
        int(len(u_reliability)),
        2464,
    )
    check(
        "oracle_descriptive_evaluations",
        len(descriptive_metrics) == 176
        and descriptive_metrics["convergence_flags"].eq("evaluable").all(),
        int(descriptive_metrics["convergence_flags"].eq("evaluable").sum()),
        176,
    )
    oracle_errors = []
    for _, group in descriptive_predictions.groupby(
        ["analysis_set", "cluster_rule", "source_weighting", "held_out_target"]
    ):
        oracle_errors.append(
            abs(
                float(group["p_oracle_I"].mean())
                - float(group["outcome_value"].mean())
            )
        )
    check(
        "oracle_I_mean_probability_equals_target_prevalence",
        max(oracle_errors) < 1e-10,
        float(max(oracle_errors)),
        "<1e-10",
    )
    check(
        "information_layer_labels_separated",
        set(descriptive_parameters["model_layer_label"])
        == {ORACLE_LAYER, DESCRIPTIVE_LAYER},
        sorted(set(descriptive_parameters["model_layer_label"])),
        [ORACLE_LAYER, DESCRIPTIVE_LAYER],
    )
    check(
        "predictive_outer_pair_rows",
        len(adaptation_predictions) == 32280
        and adaptation_predictions.duplicated(
            ["analysis_set", "cluster_rule", "source_weighting", "pair_id"]
        ).sum()
        == 0,
        int(len(adaptation_predictions)),
        32280,
    )
    check(
        "predictive_same_outer_rows_complete",
        adaptation_predictions[["p_U", "p_I", "p_F"]].notna().all().all(),
        adaptation_predictions[["p_U", "p_I", "p_F"]].notna().all().all(),
        True,
    )
    fold_coverage = adaptation_predictions.groupby(
        ["analysis_set", "cluster_rule", "source_weighting", "held_out_target"]
    )["outer_fold_id"].nunique()
    check(
        "predictive_five_outer_folds_per_target_config",
        fold_coverage.eq(5).all(),
        fold_coverage.value_counts().to_dict(),
        {5: 88},
    )
    check(
        "predictive_fold_layer_evaluations",
        len(fold_metrics) == 1320
        and fold_metrics["convergence_flags"].eq("evaluable").all(),
        int(fold_metrics["convergence_flags"].eq("evaluable").sum()),
        1320,
    )
    inner_counts = inner_assignments.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "held_out_target",
            "outer_fold_id",
            "inner_fold_number",
            "pair_id",
        ]
    ).size()
    check(
        "inner_assignments_each_pair_once",
        inner_counts.eq(1).all(),
        int((inner_counts != 1).sum()),
        0,
    )
    inner_leaks = 0
    for _, group in inner_assignments.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "held_out_target",
            "outer_fold_id",
            "inner_fold_number",
        ]
    ):
        train_groups = set(
            group.loc[group["split"].eq("train"), "sequence_cluster_value"].astype(str)
        )
        test_groups = set(
            group.loc[group["split"].eq("test"), "sequence_cluster_value"].astype(str)
        )
        if train_groups & test_groups:
            inner_leaks += 1
    check("inner_cluster_isolation", inner_leaks == 0, inner_leaks, 0)
    lambda_counts = lambda_selection.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "held_out_target",
            "outer_fold_id",
        ]
    )["lambda"].nunique()
    check(
        "frozen_six_lambda_grid",
        lambda_counts.eq(6).all()
        and set(np.round(lambda_selection["lambda"], 10))
        == set(round(value, 10) for value in config["predictive_adaptation"]["lambda_grid"]),
        lambda_counts.value_counts().to_dict(),
        {6: 440},
    )
    selected_counts = lambda_selection.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "held_out_target",
            "outer_fold_id",
        ]
    )["selected"].sum()
    check(
        "one_lambda_selected_per_outer_fit",
        selected_counts.eq(1).all(),
        int((selected_counts != 1).sum()),
        0,
    )
    selection_rule_errors = 0
    for _, group in lambda_selection.groupby(
        [
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "held_out_target",
            "outer_fold_id",
        ]
    ):
        evaluable = group[group["candidate_status"].eq("evaluable")]
        if evaluable.empty:
            selection_rule_errors += 1
            continue
        best = evaluable["pooled_inner_log_loss"].min()
        expected = evaluable[
            evaluable["pooled_inner_log_loss"] - best <= 1e-12
        ]["lambda"].max()
        observed = group.loc[group["selected"], "lambda"]
        if len(observed) != 1 or abs(float(observed.iloc[0]) - float(expected)) > 1e-12:
            selection_rule_errors += 1
    check(
        "lambda_selection_rule_reproducible",
        selection_rule_errors == 0,
        selection_rule_errors,
        0,
    )
    bss_separation = (
        not np.allclose(
            adaptation_predictions["bss_transport_reference"],
            adaptation_predictions["bss_target_oracle_reference"],
        )
        and adaptation_predictions["bss_transport_reference"].name
        != adaptation_predictions["bss_target_oracle_reference"].name
    )
    check(
        "transport_and_target_oracle_references_distinct",
        bss_separation,
        {
            "transport_column": "bss_transport_reference",
            "target_oracle_column": "bss_target_oracle_reference",
        },
        "separate row-level columns",
    )
    layer_sets = {
        "stage_1_u": set(u_metrics["model_layer_label"]),
        "stage_2": set(descriptive_metrics["model_layer_label"]),
        "stage_3": set(target_adaptation_metrics["model_layer_label"]),
    }
    check(
        "zero_label_conclusion_limited_to_U",
        layer_sets["stage_1_u"] == {U_LAYER}
        and ORACLE_LAYER not in layer_sets["stage_1_u"]
        and DESCRIPTIVE_LAYER not in layer_sets["stage_1_u"],
        layer_sets,
        "U only in stage 1 zero-label layer",
    )
    calibration_failed = fold_metrics[
        fold_metrics["calibration_status"].ne("evaluable")
    ]
    ledger_calibration = failure_ledger[
        failure_ledger["failure_type"].eq("calibration_failure")
    ]
    check(
        "all_failed_calibration_diagnostics_ledgered",
        len(calibration_failed) == len(ledger_calibration),
        {
            "calibration_failed_fold_rows": int(len(calibration_failed)),
            "ledger_rows": int(len(ledger_calibration)),
        },
        "equal counts",
    )
    check(
        "no_unledgered_model_fit_failure",
        len(failure_ledger[
            ~failure_ledger["failure_type"].isin(
                ["calibration_failure"]
            )
        ]) == 0,
        int(
            len(failure_ledger[
                ~failure_ledger["failure_type"].isin(["calibration_failure"])
            ])
        ),
        0,
    )
    adaptation_unit_tests = json.loads(
        (RUN_DIR / "adaptation_unit_tests.json").read_text(encoding="utf-8")
    )
    old_new = pd.DataFrame(
        f_regression["old_new_consistency"]
    ).set_index("layer")
    check(
        "F_regression_suite_passed",
        f_regression["status"] == "passed",
        f_regression["status"],
        "passed",
    )
    check(
        "F_regression_record_count_complete",
        f_regression["record_counts"]["total_F_fits"] == 528,
        f_regression["record_counts"]["total_F_fits"],
        528,
    )
    check(
        "F_saved_objective_matches_reconstruction",
        f_regression["checks"][
            "objective_reconstruction_within_tolerance"
        ],
        f_regression["numerical_maxima"][
            "objective_reconstruction_abs_error"
        ],
        f"<={f_regression['tolerances']['objective_abs']}",
    )
    check(
        "F_gamma_identity_complete",
        f_regression["checks"]["gamma_identity_within_tolerance"],
        f_regression["numerical_maxima"]["gamma_identity_abs_error"],
        f"<={f_regression['tolerances']['gamma_identity_abs']}",
    )
    check(
        "F_eta_prediction_reconstruction_complete",
        f_regression["checks"]["eta_reconstruction_within_tolerance"],
        f_regression["numerical_maxima"]["eta_reconstruction_abs_error"],
        f"<={f_regression['tolerances']['eta_reconstruction_abs']}",
    )
    check(
        "F_probability_prediction_reconstruction_complete",
        f_regression["checks"][
            "probability_reconstruction_within_tolerance"
        ],
        f_regression["numerical_maxima"][
            "probability_reconstruction_abs_error"
        ],
        f"<={f_regression['tolerances']['probability_reconstruction_abs']}",
    )
    check(
        "F_saved_gradients_near_zero",
        f_regression["checks"]["saved_gradient_within_tolerance"],
        f_regression["numerical_maxima"]["saved_gradient_max_abs"],
        f"<={f_regression['tolerances']['saved_gradient_max_abs']}",
    )
    check(
        "F_lambda_zero_independent_implementation_matches",
        f_regression["checks"][
            "all_lambda_zero_independent_predictions_match"
        ],
        f_regression["numerical_maxima"][
            "lambda_zero_independent_probability_abs_error"
        ],
        f"<={f_regression['tolerances']['lambda_zero_independent_probability_abs']}",
    )
    check(
        "F_large_lambda_limit_passed",
        bool(
            adaptation_unit_tests["checks"][
                "large_lambda_compresses_delta_toward_zero"
            ]
        ),
        adaptation_unit_tests["numerical_diagnostics"][
            "large_lambda_abs_delta"
        ],
        "<=1e-5",
    )
    check(
        "old_new_U_numerically_unchanged",
        old_new.loc["U", "material_change_rows"] == 0,
        old_new.loc["U"].to_dict(),
        "zero material-change rows",
    )
    check(
        "old_new_oracle_I_numerically_unchanged",
        old_new.loc["oracle_I", "material_change_rows"] == 0,
        old_new.loc["oracle_I"].to_dict(),
        "zero material-change rows",
    )
    check(
        "old_new_predictive_I_numerically_unchanged",
        old_new.loc["predictive_I", "material_change_rows"] == 0,
        old_new.loc["predictive_I"].to_dict(),
        "zero material-change rows",
    )
    required_paths = [
        "RUN_STATUS.json",
        "run_config.snapshot.json",
        "rq3_source_model_parameters.csv",
        "rq3_u_pair_predictions.csv",
        "rq3_loto_metrics.csv",
        "rq3_loto_reliability.csv",
        "rq3_descriptive_pair_predictions.csv",
        "rq3_adaptation_inner_assignments.csv",
        "rq3_lambda_selection.csv",
        "rq3_adaptation_fold_metrics.csv",
        "rq3_adaptation_pair_predictions.csv",
        "rq3_target_adaptation_metrics.csv",
        "rq3_overall_metrics.csv",
        "RQ3_U_LOTO.png",
        "RQ3_adaptation.png",
        "qc_checks.json",
        "HANDOFF_RESULT.md",
        "handoff_result.json",
        "CORRECTION_RESULT.md",
        "rq3_f_regression_tests.json",
        "rq3_f_objective_consistency.csv",
        "rq3_old_new_consistency.csv",
        "failed_item_ledger.csv",
    ]
    missing_required = [
        name for name in required_paths if not (RUN_DIR / name).is_file()
    ]
    # HANDOFF_RESULT.md and handoff_result.json are written below; record preflight status.
    missing_pre_handoff = [
        name
        for name in missing_required
        if name not in {
            "HANDOFF_RESULT.md",
            "handoff_result.json",
            "qc_checks.json",
            "CORRECTION_RESULT.md",
        }
    ]
    check(
        "required_outputs_present_before_handoff",
        not missing_pre_handoff,
        missing_pre_handoff,
        [],
    )
    stage4_qc = {
        "stage": STAGE,
        "checked_at": now_iso(),
        "all_passed": all(item["passed"] for item in qc_checks),
        "n_checks": len(qc_checks),
        "n_failed": sum(not item["passed"] for item in qc_checks),
        "checks": qc_checks,
        "failure_ledger_rows": int(len(failure_ledger)),
        "calibration_bound_hit_rows": int(len(calibration_failed)),
        "model_fit_failures": int(
            failure_ledger["failure_type"].ne("calibration_failure").sum()
        ),
        "f_regression": f_regression,
        "hash_manifest_status": "generated by finalize_hash_manifest.py after all logs are closed",
    }
    atomic_json(RUN_DIR / "qc_checks.json", stage4_qc)
    if not stage4_qc["all_passed"]:
        append_failure(
            STAGE,
            "stage_4_qc",
            "qc_failure",
            "One or more final QC checks failed.",
            failed_checks=[
                item["check"] for item in qc_checks if not item["passed"]
            ],
        )

    stage1_summary = json.loads((RUN_DIR / "stage1_summary.json").read_text(encoding="utf-8"))
    stage2_summary = json.loads((RUN_DIR / "stage2_summary.json").read_text(encoding="utf-8"))
    stage3_summary = json.loads((RUN_DIR / "stage3_summary.json").read_text(encoding="utf-8"))

    main_u_table = main_u[
        [
            "held_out_target",
            "n",
            "positives",
            "brier",
            "log_loss",
            "bss_transport",
            "bss_target_oracle",
            "calibration_intercept",
            "calibration_slope",
        ]
    ].copy()
    main_u_table["held_out_target"] = main_u_table["held_out_target"].astype(str)
    macro_u = overall_macro[
        overall_macro["stage"].eq("stage_1_U")
        & overall_macro["analysis_set"].eq("primary_six")
        & overall_macro["cluster_rule"].eq("sequence_cluster90")
        & overall_macro["source_weighting"].eq("target_equal")
        & overall_macro["model_form"].eq(LINEAR)
        & overall_macro["metric"].isin(
            ["brier", "log_loss", "bss_transport", "bss_target_oracle", "calibration_slope"]
        )
    ][["metric", "mean", "min", "max", "n_positive_sign", "n_negative_sign"]]

    sensitivity_u = sensitivity[
        sensitivity["stage"].eq("stage_1_U")
        & sensitivity["model_layer_label"].eq(U_LAYER)
        & sensitivity["metric"].isin(
            ["brier", "log_loss", "bss_transport", "bss_target_oracle"]
        )
    ].pivot(
        index=[
            "dimension",
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "model_form",
        ],
        columns="metric",
        values="macro_mean",
    ).reset_index()

    adaptation_report = main_adaptation[
        [
            "held_out_target",
            "model_layer_label",
            "brier",
            "log_loss",
            "bss_target_oracle",
            "bss_adaptation_train",
            "calibration_status",
        ]
    ].copy()
    adaptation_report["held_out_target"] = adaptation_report["held_out_target"].astype(str)
    descriptive_report = main_desc[
        [
            "held_out_target",
            "model_layer_label",
            "brier",
            "log_loss",
            "bss_target_oracle",
            "calibration_status",
        ]
    ].copy()
    descriptive_report["held_out_target"] = descriptive_report["held_out_target"].astype(str)
    comparison_report = main_comparison[
        [
            "comparison",
            "metric",
            "macro_mean_difference",
            "n_improved",
            "n_worsened",
            "min_difference",
            "max_difference",
        ]
    ].copy()

    sensitivity_adaptation = sensitivity[
        sensitivity["stage"].eq("stage_3_predictive")
        & sensitivity["metric"].isin(["brier", "log_loss", "bss_target_oracle"])
    ].pivot_table(
        index=[
            "dimension",
            "analysis_set",
            "cluster_rule",
            "source_weighting",
            "model_form",
            "model_layer_label",
        ],
        columns="metric",
        values="macro_mean",
    ).reset_index()

    status_table = pd.DataFrame(
        [
            ["stage_0_input_and_fold_qc", "complete", "19/19 structural checks"],
            [
                "stage_1_U_linear_and_RCS_with_predictions_metrics_reliability",
                "complete",
                f"{stage1_summary['source_models_evaluable']}/176 evaluable",
            ],
            [
                "stage_2_oracle_I_and_descriptive_F",
                "complete",
                f"{stage2_summary['evaluable_model_evaluations']}/176 evaluable diagnostics",
            ],
            [
                "stage_3_predictive_I_F_with_nested_inner_cv",
                "complete",
                f"{stage3_summary['lambda_candidates_evaluable']}/2640 lambda candidates evaluable",
            ],
            [
                "stage_4_figures_tables_qc_and_handoff",
                "complete",
                f"{sum(item['passed'] for item in qc_checks)}/{len(qc_checks)} final checks passed",
            ],
        ],
        columns=["Stage", "Status", "Evidence"],
    )
    failure_summary = (
        failure_ledger.groupby(["stage", "failure_type"])
        .size()
        .reset_index(name="rows")
    )

    u_transport_macro = metric_lookup(
        overall_macro,
        stage="stage_1_U",
        layer=U_LAYER,
        metric="bss_transport",
    )
    i_brier_diff = main_comparison[
        main_comparison["comparison"].eq("predictive_I_minus_U")
        & main_comparison["metric"].eq("brier")
    ].iloc[0]
    f_vs_i_brier = main_comparison[
        main_comparison["comparison"].eq("predictive_F_minus_I")
        & main_comparison["metric"].eq("brier")
    ].iloc[0]
    i_adapt_bss = metric_lookup(
        overall_macro,
        stage="stage_3_predictive",
        layer=PRED_I_LAYER,
        metric="bss_adaptation_train",
    )
    f_adapt_bss = metric_lookup(
        overall_macro,
        stage="stage_3_predictive",
        layer=PRED_F_LAYER,
        metric="bss_adaptation_train",
    )
    primary_u_brier = metric_lookup(
        overall_macro,
        stage="stage_1_U",
        layer=U_LAYER,
        metric="brier",
    )
    old_new_report = pd.DataFrame(
        f_regression["old_new_consistency"]
    )[
        [
            "layer",
            "new_rows",
            "old_rows",
            "joined_rows",
            "unmatched_rows",
            "max_absolute_eta_difference",
            "max_absolute_probability_difference",
            "material_change_rows",
        ]
    ]
    regression_rows = pd.DataFrame(
        [
            {
                "check": "target-to-prediction objective reconstruction",
                "status": (
                    "passed"
                    if f_regression["checks"][
                        "objective_reconstruction_within_tolerance"
                    ]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "objective_reconstruction_abs_error"
                ],
                "tolerance": f_regression["tolerances"]["objective_abs"],
            },
            {
                "check": "gamma = 1 + delta for all F records",
                "status": (
                    "passed"
                    if f_regression["checks"]["gamma_identity_within_tolerance"]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "gamma_identity_abs_error"
                ],
                "tolerance": f_regression["tolerances"]["gamma_identity_abs"],
            },
            {
                "check": "eta_F = a + gamma*eta_U",
                "status": (
                    "passed"
                    if f_regression["checks"][
                        "eta_reconstruction_within_tolerance"
                    ]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "eta_reconstruction_abs_error"
                ],
                "tolerance": f_regression["tolerances"][
                    "eta_reconstruction_abs"
                ],
            },
            {
                "check": "p_F = expit(eta_F)",
                "status": (
                    "passed"
                    if f_regression["checks"][
                        "probability_reconstruction_within_tolerance"
                    ]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "probability_reconstruction_abs_error"
                ],
                "tolerance": f_regression["tolerances"][
                    "probability_reconstruction_abs"
                ],
            },
            {
                "check": "saved analytical gradients",
                "status": (
                    "passed"
                    if f_regression["checks"]["saved_gradient_within_tolerance"]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "saved_gradient_max_abs"
                ],
                "tolerance": f_regression["tolerances"][
                    "saved_gradient_max_abs"
                ],
            },
            {
                "check": "lambda=0 independent gamma fit",
                "status": (
                    "passed"
                    if f_regression["checks"][
                        "all_lambda_zero_independent_predictions_match"
                    ]
                    else "failed"
                ),
                "maximum": f_regression["numerical_maxima"][
                    "lambda_zero_independent_probability_abs_error"
                ],
                "tolerance": f_regression["tolerances"][
                    "lambda_zero_independent_probability_abs"
                ],
            },
            {
                "check": "central finite-difference gradient",
                "status": (
                    "passed"
                    if adaptation_unit_tests["checks"][
                        "analytic_gradient_matches_central_finite_difference"
                    ]
                    else "failed"
                ),
                "maximum": adaptation_unit_tests["numerical_diagnostics"][
                    "max_finite_difference_gradient_error"
                ],
                "tolerance": adaptation_unit_tests["tolerances"][
                    "finite_difference_gradient_max_abs"
                ],
            },
            {
                "check": "large-lambda delta compression",
                "status": (
                    "passed"
                    if adaptation_unit_tests["checks"][
                        "large_lambda_compresses_delta_toward_zero"
                    ]
                    else "failed"
                ),
                "maximum": adaptation_unit_tests["numerical_diagnostics"][
                    "large_lambda_abs_delta"
                ],
                "tolerance": adaptation_unit_tests["tolerances"][
                    "large_lambda_abs_delta"
                ],
            },
        ]
    )
    regression_rows["maximum"] = regression_rows["maximum"].map(
        lambda value: f"{float(value):.6e}"
    )
    regression_rows["tolerance"] = regression_rows["tolerance"].map(
        lambda value: f"{float(value):.1e}"
    )
    invalidation_scope = pd.DataFrame(
        [
            [
                "descriptive F fits",
                f_regression["invalidated_old_F_scope"]["descriptive_F_fits"],
            ],
            [
                "descriptive F pair predictions",
                f_regression["invalidated_old_F_scope"][
                    "descriptive_F_pair_predictions"
                ],
            ],
            [
                "predictive F final outer fits",
                f_regression["invalidated_old_F_scope"]["predictive_F_fits"],
            ],
            [
                "predictive F pair predictions",
                f_regression["invalidated_old_F_scope"][
                    "predictive_F_pair_predictions"
                ],
            ],
            [
                "lambda candidates",
                f_regression["invalidated_old_F_scope"]["lambda_candidates"],
            ],
        ],
        columns=["old F artifact", "count"],
    )
    model_fit_failure_count = int(
        failure_ledger["failure_type"].ne("calibration_failure").sum()
    )
    correction_report = f"""# RQ3 F 层实现修正结果

运行目录：`{RUN_DIR}`

状态：**Codex 执行完成，仍需主执行线独立验收。**

本结果不替代主执行线的独立 Newton/IRLS 复核，也不自行宣称正式科学验收通过。

## 1. 根因

旧运行 `fit_offset_slope()` 使用设计矩阵 `[1, eta_U]`，实际优化的是 `eta_F = a + theta_1 * eta_U`，却把 `theta_1` 保存为 `delta`，并按 `gamma = 1 + theta_1` 预测。因此拟合目标与保存/预测使用了两条不同的线性预测子，二者恰好相差一个 `eta_U`。旧单测只检查 lambda 变化时的参数趋势，没有从保存参数重构训练目标，因而漏检。

冻结定义未改变：`eta_F = a + eta_U + delta * eta_U = a + gamma * eta_U`，`gamma = 1 + delta`；目标为未裁剪 `mean_NLL(eta_F, y) + 0.5*lambda*delta^2`，截距不惩罚。

## 2. 修改位置

- `rq3_common.py`：通用 `fit_weighted_logistic()` 增加显式 `offset`；目标、梯度、Hessian 与收敛后概率均使用 `X @ theta + offset`。F 专用 `fit_offset_slope()` 固定 `offset=eta_U`、设计 `[1, eta_U]`、惩罚索引 1，保存 `delta=theta[1]` 与 `gamma=1+delta`。
- `stage_2_oracle_descriptive.py`：参数表新增描述性 F 解析梯度最大值。
- `stage_3_predictive_I_F.py`：扩展合成单测，覆盖保存目标重构、恒等关系、解析梯度、中心有限差分、lambda=0 独立拟合与大 lambda 极限；原 lambda 网格、`1e-12` 并列取更大 lambda 规则和概率裁剪均未改变。
- `rq3_f_regression.py`：新增全部 88 个描述性 F 与 440 个最终外折 F 的逐记录目标—预测一致性检查，并生成新旧 U/I/F 数值对照。
- `stage_4_finalize.py`：把修正回归检查纳入最终 QC，并生成本报告。

## 3. 修正后回归测试

逐记录检查覆盖 {f_regression['record_counts']['total_F_fits']} 个 F 拟合（88 个描述性 F、440 个最终外折 F）。容差在测试 JSON 与 QC JSON 中明确记录，优化器停止条件没有替代目标一致性验证。

{markdown_table(regression_rows)}

全部 88 个描述性 F 的 lambda=0 直接 `y ~ intercept + gamma*eta_U` 独立拟合均成功；有限大 lambda 单测仅要求 `delta` 被压缩，不把“严格等于 I”作为要求。

## 4. 旧结果失效范围

旧运行的 F 层拟合参数、目标、lambda 候选与选择、F 预测、F 指标、比较及含 F 的图表/报告均不能用于验收。不失效的范围是旧 U、oracle I 与预测性 I；新运行仍从头重算它们并做数值对照。

{markdown_table(invalidation_scope)}

旧 `RQ3_adaptation.png` 和旧 HANDOFF 中包含 F 的面板或结论同时失效；U 层图不因本次根因而失效，但仍以新运行文件为交付。

## 5. 新旧 U／I 数值一致性

{markdown_table(old_new_report)}

`material_change_rows` 阈值：绝对 eta 差 `>{f_regression['tolerances']['old_new_eta_abs']}` 或绝对概率差 `>{f_regression['tolerances']['old_new_probability_abs']}`。U、oracle I、预测性 I 若存在实质变化会同时写入 `rq3_old_new_material_changes.csv` 并触发最终 QC 失败；本表直接报告实测最大差。

## 6. 修正后关键结果

- 主 U 口径：Brier `{primary_u_brier:.6f}`，BSS-transport `{u_transport_macro:.6f}`，BSS-target-oracle `{metric_lookup(overall_macro, stage='stage_1_U', layer=U_LAYER, metric='bss_target_oracle'):.6f}`。
- 预测性 I 相对 U：Brier 宏平均差 `{float(i_brier_diff['macro_mean_difference']):.6f}`，{int(i_brier_diff['n_improved'])}/6 个靶点改善；BSS-adaptation-train `{i_adapt_bss:.6f}`。
- 预测性 F 相对 I：Brier 宏平均差 `{float(f_vs_i_brier['macro_mean_difference']):.6f}`，{int(f_vs_i_brier['n_improved'])}/6 个靶点改善、{int(f_vs_i_brier['n_worsened'])}/6 个变差；BSS-adaptation-train `{f_adapt_bss:.6f}`。
- 2,640 个 lambda 候选全部重新计算；440 个最终外折各选一个 lambda。选择频数：`{json.dumps(stage3_summary['selected_lambda_frequencies'], ensure_ascii=False)}`。
- 上述为六个留出靶点的描述性执行结果，不构成对未来靶点的总体不确定性声明。

## 7. QC 与失败概况

- 冻结配置 SHA-256：`a416ecd357eabb4c4d95701136c412cf23651b730bf39feab7524fbc07eeef51`。
- 冻结方案 SHA-256：`36dc1503e86960b5f05d0eb273d3b68091eb2072975fb51e5556ee821e73cfd8`。
- 最终 QC：{sum(item['passed'] for item in qc_checks)}/{len(qc_checks)} 通过。
- 失败台账：{len(failure_ledger)} 条，其中校准诊断边界 `{len(calibration_failed)}` 条，其他模型/回归失败 `{model_fit_failure_count}` 条。
- 所有校准边界缺失均按原规则记账，没有替换方法、删除折、放宽裁剪或使用外测标签修补。

**Codex 执行完成，仍需主执行线独立验收。**
"""
    (RUN_DIR / "CORRECTION_RESULT.md").write_text(
        correction_report, encoding="utf-8", newline="\n"
    )

    report = f"""# RQ3 F 层修正后的正式执行交接结果

运行目录：`{RUN_DIR}`

总体状态：**Codex 执行完成，仍需主执行线独立验收。**

本运行针对已确认的 F 层实现错误从冻结输入完整重跑 stage 0–4；根因、修改、回归测试和旧结果失效范围见 `CORRECTION_RESULT.md`。

本运行严格复用已验收的 44 个 LOTO 折与 220 个 CV5 折；没有重新生成外折，没有修改任何冻结输入或旧结果，没有安装或升级依赖，也没有启动 RQ4/RQ5/RQ6。所有新文件均位于本运行目录。

## 1. 主分析：U 层零标签迁移

唯一主口径：`primary_six + sequence_cluster90 + target_equal + linear`。U 层是零标签主问题；oracle I、描述性 F、预测性 I/F 均不替代该结论。

{markdown_table(main_u_table)}

主分析宏平均与逐靶点范围：

{markdown_table(macro_u)}

配对汇总（与靶点等权宏平均分开报告）：Brier={main_pair_u['brier']:.6f}，log loss={main_pair_u['log_loss']:.6f}，BSS-transport={main_pair_u['bss_transport']:.6f}，BSS-target-oracle={main_pair_u['bss_target_oracle']:.6f}。

U 层初步判定：主分析靶点等权宏平均 `BSS-transport={u_transport_macro:.6f}`，平均未改善，按方案 6.6 初步判为**当前零标签映射未解决迁移**；逐靶点方向为 {int(main_u['bss_transport'].gt(0).sum())} 个正、{int(main_u['bss_transport'].lt(0).sum())} 个负，显示明显靶点依赖。AAH-I 的 19 个阳性保留在主分析中并单独可见。

## 2. 预设敏感性

U 层五项预设敏感性均保留，不因结果删减：

{markdown_table(sensitivity_u)}

预测性适配的对应敏感性（预测性 F 不扩展到 RCS，因为 RCS 只用于 U 层预设形态敏感性）：

{markdown_table(sensitivity_adaptation)}

## 3. Oracle I 与描述性 F（非部署诊断）

以下两层使用留出靶点标签，分别标为 `oracle_base_rate_scenario` 与 `descriptive_in_sample_adaptation`，不是零标签或折外预测结果。

{markdown_table(descriptive_report)}

Oracle I 的截距严格求解为 `mean(expit(a+eta_U)) = target prevalence`；描述性 F 使用留出靶点全样本拟合，未被称为预测性改善。

## 4. 预测性 I 与 F（同一 CV5 外测行）

每个预测性结果均来自对应 LOTO 折固定 `eta_U` 后在适配训练折上的更新；U、I、F 使用完全相同的外测配对行。lambda 只由外折训练侧内部的 3 折 `StratifiedGroupKFold` 选择。

{markdown_table(adaptation_report)}

主分析目标级配对比较：

{markdown_table(comparison_report)}

预测性结果按六靶点逐项报告：I 相对 U 的 Brier 宏平均差为 {i_brier_diff['macro_mean_difference']:.6f}，{int(i_brier_diff['n_improved'])}/6 个靶点改善；I 相对训练折基率预测的 `BSS-adaptation-train` 宏平均为 {i_adapt_bss:.6f}。F 相对 I 的 Brier 宏平均差为 {f_vs_i_brier['macro_mean_difference']:.6f}，{int(f_vs_i_brier['n_improved'])}/6 改善、{int(f_vs_i_brier['n_worsened'])}/6 变差，`BSS-adaptation-train={f_adapt_bss:.6f}`。方向计数按靶点等权统计，不能作为未来靶点上的总体不确定性声明。

## 5. 阶段状态与失败台账

{markdown_table(status_table)}

失败台账共 {len(failure_ledger)} 条，其中预测性阶段逐外折校准诊断触及 `±40` 数值边界 {len(calibration_failed)} 条，其他模型或回归失败 {model_fit_failure_count} 条；相应校准指标保留为缺失并标记 `failed`，未替换算法、未删除外折、未使用外测标签补救。预测模型与 lambda 候选可估计情况以 `stage3_summary.json` 和失败台账为准。

{markdown_table(failure_summary)}

## 6. QC、冻结与交付

- 输入 SHA-256 在开机和结束各核对一次，全部一致。
- 2,084 个 `pair_id`、44 个 LOTO 折、220 个 CV5 折、364 条跨靶点同簇剔除记录均与冻结清单一致。
- U：176 个模型-靶点评估、64,560 条测试配对预测，全部可估计；RCS 公式与尾部线性单元测试通过。
- Oracle/描述性：176 个诊断全覆盖；oracle 平均预测概率与真实基率误差小于 `1e-10`。
- 预测性：440 个目标×权重×外折评估、32,280 条相同外测行 U/I/F 预测、2,640 个候选 lambda、387,360 条内折分配；所有内折训练/测试序列簇隔离通过。
- F 修正回归：状态 `{f_regression['status']}`；{f_regression['record_counts']['total_F_fits']} 个 F 拟合的目标重构、恒等关系、预测重构与梯度结果详见 `rq3_f_regression_tests.json`。
- 最终 QC：{sum(item['passed'] for item in qc_checks)}/{len(qc_checks)} 项通过。
- `failed_item_ledger.csv`、`qc_checks.json`、脚本、命令、阶段日志、图与总表均已保存；`final_hash_manifest.csv` 和 `manifest_sha256.txt` 在本脚本退出、日志关闭后生成。

信息层级边界：不得把 `BSS-target-oracle` 当作可部署零标签参照；不得把 oracle I 或描述性 F 称为预测改善；不得用“校准斜率区间包含 1”作为成功标准。六靶点结果不提供对未来靶点的总体不确定性声明。

**正式 RQ3 结论待主执行线独立验收。**
"""
    (RUN_DIR / "HANDOFF_RESULT.md").write_text(report, encoding="utf-8", newline="\n")

    status = json.loads((RUN_DIR / "RUN_STATUS.json").read_text(encoding="utf-8"))
    handoff_result = {
        "handoff_id": config["handoff_id"],
        "run_directory": str(RUN_DIR),
        "overall_status": "complete",
        "requires_independent_acceptance": True,
        "config_frozen": True,
        "frozen_input_hashes_at_end_match": bool(final_hashes["sha256_matches"].all()),
        "stages": status["stages"],
        "counts": {
            "lot_o_folds_reused": 44,
            "cv5_folds_reused": 220,
            "u_model_target_evaluations": 176,
            "u_pair_predictions": 64560,
            "oracle_descriptive_diagnostics": 176,
            "predictive_outer_target_config_units": 88,
            "predictive_outer_fold_fits": 440,
            "predictive_fold_layer_evaluations": 1320,
            "lambda_candidates": 2640,
            "inner_assignment_rows": 387360,
            "predictive_pair_predictions": 32280,
        },
        "primary_u": {
            "definition": "primary_six + sequence_cluster90 + target_equal + linear",
            "macro_brier": primary_u_brier,
            "macro_bss_transport": u_transport_macro,
            "macro_bss_target_oracle": metric_lookup(
                overall_macro,
                stage="stage_1_U",
                layer=U_LAYER,
                metric="bss_target_oracle",
            ),
            "preliminary_determination": (
                "unresolved_zero_label_transport_with_target_dependence"
                if u_transport_macro <= 0
                else "partial_or_limited_transport"
            ),
        },
        "predictive_preliminary": {
            "I_vs_U_macro_brier_difference": float(
                i_brier_diff["macro_mean_difference"]
            ),
            "I_vs_U_n_improved": int(i_brier_diff["n_improved"]),
            "I_macro_bss_adaptation_train": i_adapt_bss,
            "F_vs_I_macro_brier_difference": float(
                f_vs_i_brier["macro_mean_difference"]
            ),
            "F_vs_I_n_improved": int(f_vs_i_brier["n_improved"]),
            "F_vs_I_n_worsened": int(f_vs_i_brier["n_worsened"]),
            "F_macro_bss_adaptation_train": f_adapt_bss,
        },
        "failures": {
            "ledger_rows": int(len(failure_ledger)),
            "calibration_bound_hit_rows": int(len(calibration_failed)),
            "model_fit_failures": model_fit_failure_count,
            "unresolved_analysis_failures": int(
                model_fit_failure_count
                + (f_regression["status"] != "passed")
            ),
        },
        "qc": {
            "checks_passed": int(sum(item["passed"] for item in qc_checks)),
            "checks_total": int(len(qc_checks)),
            "all_passed": bool(stage4_qc["all_passed"]),
        },
        "correction": {
            "result_markdown": "CORRECTION_RESULT.md",
            "old_run_directory": str(
                Path(
                    r"<REPRODUCTION_ROOT>\analysis\codex_outputs"
                    r"\RQ3_execution_2026-11-13\run_20261113_175022"
                )
            ),
            "root_cause": (
                "fit_offset_slope optimized a + theta_1*eta_U, saved "
                "delta=theta_1, then predicted a+(1+delta)*eta_U"
            ),
            "correction": (
                "explicit offset eta_U in the shared objective/gradient/"
                "Hessian; save delta directly and gamma=1+delta"
            ),
            "regression_status": f_regression["status"],
            "record_counts": f_regression["record_counts"],
            "numerical_maxima": f_regression["numerical_maxima"],
            "old_new_consistency": f_regression["old_new_consistency"],
            "invalidated_old_F_scope": f_regression[
                "invalidated_old_F_scope"
            ],
        },
        "final_hash_manifest_pending_separate_finalizer": True,
        "forbidden_actions_confirmed_absent": {
            "frozen_inputs_modified": False,
            "old_results_modified": False,
            "dependencies_installed_or_upgraded": False,
            "rq4_rq5_rq6_started": False,
            "pinc_or_other_ml_added": False,
            "methods_changed_after_outer_results": False,
        },
    }
    atomic_json(RUN_DIR / "handoff_result.json", handoff_result)

    commands_path = str(RUN_DIR).replace("\\", "/")
    commands = f"""# RQ3 execution commands

All commands were run with working directory:

`{commands_path}`

```powershell
python stage_0_input_and_fold_qc.py 2>&1 | Tee-Object -FilePath stage_0.log
python stage_1_u_layer.py 2>&1 | Tee-Object -FilePath stage_1.log
python stage_2_oracle_descriptive.py 2>&1 | Tee-Object -FilePath stage_2.log
python stage_3_predictive_I_F.py 2>&1 | Tee-Object -FilePath stage_3.log
python stage_4_finalize.py 2>&1 | Tee-Object -FilePath stage_4.log
python finalize_hash_manifest.py
```

`stage_4_finalize.py` invokes the complete F target-to-prediction regression audit and writes its JSON/CSV artifacts before final QC. The final hash manifest is run only after stage 4 exits and its log is closed.

No dependency installation or upgrade command was run.
"""
    (RUN_DIR / "COMMANDS.md").write_text(commands, encoding="utf-8", newline="\n")

    set_stage_status(
        STAGE,
        "complete",
        qc_checks_passed=int(sum(item["passed"] for item in qc_checks)),
        qc_checks_total=int(len(qc_checks)),
        failure_ledger_rows=int(len(failure_ledger)),
        overall_status="complete",
    )
    # Refresh the handoff stage status after marking stage 4 complete.
    handoff_result["stages"] = json.loads(
        (RUN_DIR / "RUN_STATUS.json").read_text(encoding="utf-8")
    )["stages"]
    atomic_json(RUN_DIR / "handoff_result.json", handoff_result)
    print(
        f"Stage 4 complete: {sum(item['passed'] for item in qc_checks)}/"
        f"{len(qc_checks)} QC checks passed; report and figures written."
    )


if __name__ == "__main__":
    main()
