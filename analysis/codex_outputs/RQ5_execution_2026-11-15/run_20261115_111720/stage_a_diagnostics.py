"""RQ5 Stage A: same-sample Pinc versus ipTM diagnostics."""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rq5_common import (
    ANALYSIS_SETS,
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    CLUSTER_RULE,
    RELIABILITY_EDGES,
    RUN_DIR,
    SCORES,
    analysis_frame,
    append_failure,
    atomic_csv,
    atomic_json,
    brier,
    load_config,
    load_core_inputs,
    log_loss,
    performance_point_metrics,
    reliability_frame,
    safe_ap,
    safe_auc,
    save_figure,
    score_point_metrics,
    set_stage_status,
    summarize_target_metrics,
)

STAGE = "stage_a_same_sample_diagnostics"


def bootstrap_target_pair(
    frame: pd.DataFrame,
    score_columns: dict[str, str],
    clip: float,
) -> tuple[list[dict], dict]:
    y = frame["outcome_value"].to_numpy(dtype=float)
    clusters = sorted(frame["sequence_cluster90"].astype(str).unique())
    indices_by_cluster = [
        np.flatnonzero(frame["sequence_cluster90"].astype(str).to_numpy() == cluster)
        for cluster in clusters
    ]
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = {
        metric: {"auc": [], "ap": [], "brier": []}
        for metric in ("ipTM", "Pinc")
    }
    paired = {"auc": [], "ap": [], "brier": []}
    for _ in range(BOOTSTRAP_REPLICATES):
        sampled = rng.choice(len(clusters), size=len(clusters), replace=True)
        chosen = np.concatenate([indices_by_cluster[index] for index in sampled])
        yy = y[chosen]
        values = {}
        for score_label, column in score_columns.items():
            ss = frame[column].to_numpy(dtype=float)[chosen]
            values[score_label] = {
                "auc": safe_auc(yy, ss),
                "ap": safe_ap(yy, ss),
                "brier": brier(yy, ss),
            }
        for metric in paired:
            left = values["ipTM"][metric]
            right = values["Pinc"][metric]
            if np.isfinite(left):
                draws["ipTM"][metric].append(float(left))
            if np.isfinite(right):
                draws["Pinc"][metric].append(float(right))
            if np.isfinite(left) and np.isfinite(right):
                paired[metric].append(float(right - left))

    observed = {}
    for score_label, column in score_columns.items():
        ss = frame[column].to_numpy(dtype=float)
        observed[score_label] = {
            "auc": safe_auc(y, ss),
            "ap": safe_ap(y, ss),
            "brier": brier(y, ss),
        }
    rows = []
    for metric in ("auc", "ap", "brier"):
        differences = np.asarray(paired[metric], dtype=float)
        low, high = (
            np.percentile(differences, [2.5, 97.5])
            if differences.size
            else (np.nan, np.nan)
        )
        rows.append(
            {
                "metric": metric,
                "iptm_value": observed["ipTM"][metric],
                "pinc_value": observed["Pinc"][metric],
                "paired_difference_pinc_minus_iptm": observed["Pinc"][metric] - observed["ipTM"][metric],
                "bootstrap_mean_difference": float(differences.mean()) if differences.size else np.nan,
                "bootstrap_95_low": float(low),
                "bootstrap_95_high": float(high),
                "n_valid_boot": int(differences.size),
                "n_requested_boot": BOOTSTRAP_REPLICATES,
                "n_clusters": int(len(clusters)),
                "seed": BOOTSTRAP_SEED,
                "same_resampled_clusters_for_both_scores": True,
                "interval_type": "percentile_95",
            }
        )
    audit = {
        "n_clusters": int(len(clusters)),
        "n_boot": BOOTSTRAP_REPLICATES,
        "seed": BOOTSTRAP_SEED,
        "valid": {metric: int(len(values)) for metric, values in paired.items()},
    }
    return rows, audit


def main() -> None:
    set_stage_status(STAGE, "running")
    config = load_config()
    _, _, _, joined = load_core_inputs(config)
    clip = float(config["probability_clip"])

    target_rows = []
    reliability_rows = []
    bootstrap_rows = []
    bootstrap_audits = []
    for analysis_set in ANALYSIS_SETS:
        outcome_name = config["analysis_sets"][analysis_set]["outcome"]
        frame = analysis_frame(config, joined, analysis_set)
        pooled_scores = {}
        for target, target_frame in frame.groupby("target", sort=True):
            target_frame = target_frame.sort_values("pair_id").reset_index(drop=True)
            y = target_frame["outcome_value"].to_numpy(dtype=float)
            target_metrics = {}
            for score_label, (raw_column, _) in SCORES.items():
                raw_score = target_frame[raw_column].to_numpy(dtype=float)
                point = score_point_metrics(y, raw_score, clip=clip)
                target_metrics[score_label] = point
                pooled_scores.setdefault(score_label, []).append(
                    pd.DataFrame(
                        {
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "score": score_label,
                            "target": target,
                            "outcome_value": y,
                            "raw_score": raw_score,
                        }
                    )
                )
                target_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": CLUSTER_RULE,
                        "target": target,
                        "score": score_label,
                        "raw_score_is_probability": "diagnostic_only",
                        **point,
                    }
                )
                for reliability in reliability_frame(y, raw_score, RELIABILITY_EDGES):
                    reliability_rows.append(
                        {
                            "analysis_set": analysis_set,
                            "outcome_name": outcome_name,
                            "cluster_rule": CLUSTER_RULE,
                            "target": target,
                            "score": score_label,
                            "raw_score_is_probability": "diagnostic_only",
                            **reliability,
                        }
                    )
            paired_rows, audit = bootstrap_target_pair(
                target_frame,
                {"ipTM": "boltz2_iptm", "Pinc": "pinc"},
                clip,
            )
            for row in paired_rows:
                row.update(
                    {
                        "analysis_set": analysis_set,
                        "outcome_name": outcome_name,
                        "cluster_rule": CLUSTER_RULE,
                        "target": target,
                        "comparison": "Pinc_minus_ipTM",
                        "brier_difference_negative_is_better": True,
                        "auc_difference_positive_is_better": True,
                        "ap_difference_positive_is_better": True,
                    }
                )
            bootstrap_rows.extend(paired_rows)
            bootstrap_audits.append(
                {
                    "analysis_set": analysis_set,
                    "target": target,
                    "outcome_name": outcome_name,
                    **audit,
                }
            )

        for score_label, frames in pooled_scores.items():
            pooled = pd.concat(frames, ignore_index=True)
            # Pooled raw-score diagnostics are shown separately from target-equal macro.
            if score_label not in {"ipTM", "Pinc"}:
                continue
        if len(pooled_scores.get("ipTM", [])) != len(pooled_scores.get("Pinc", [])):
            append_failure(STAGE, analysis_set, "pooled_score_row_mismatch", "ipTM and Pinc pooled row sets differ")

    target_metrics = pd.DataFrame(target_rows)
    atomic_csv(RUN_DIR / "rq5_target_score_metrics.csv", target_metrics)
    reliability = pd.DataFrame(reliability_rows)
    atomic_csv(RUN_DIR / "rq5_score_reliability.csv", reliability)
    bootstrap = pd.DataFrame(bootstrap_rows)
    atomic_csv(RUN_DIR / "rq5_score_bootstrap_comparisons.csv", bootstrap)
    atomic_json(RUN_DIR / "stage_a_bootstrap_audit.json", {"items": bootstrap_audits, "all_requested_draws": True})

    metric_columns = [
        "roc_auc",
        "ap",
        "brier",
        "log_loss",
        "bss_target_constant",
        "calibration_intercept",
        "calibration_slope",
    ]
    macro = summarize_target_metrics(
        target_metrics,
        metric_columns,
        ["analysis_set", "outcome_name", "cluster_rule", "score"],
    )
    macro["baseline_dependency_note"] = (
        "AP and Brier depend on target prevalence; compare targets as diagnostics only."
    )

    pooled_rows = []
    for analysis_set in ANALYSIS_SETS:
        frame = analysis_frame(config, joined, analysis_set)
        outcome_name = config["analysis_sets"][analysis_set]["outcome"]
        for score_label, (raw_column, _) in SCORES.items():
            y = frame["outcome_value"].to_numpy(dtype=float)
            score = frame["raw_score"] if "raw_score" in frame else frame[raw_column].to_numpy(dtype=float)
            result = score_point_metrics(y, np.asarray(score, dtype=float), clip=clip)
            pooled_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": CLUSTER_RULE,
                    "score": score_label,
                    "summary_type": "pair_pooled",
                    **result,
                    "baseline_dependency_note": "Pair-pooled AP/Brier reflect the pooled target mix, not target-equal performance.",
                }
            )
    pooled = pd.DataFrame(pooled_rows)
    # Pair-pooled rows are stored in a wide companion table; macro remains long.
    atomic_csv(RUN_DIR / "rq5_score_pair_summary.csv", pooled)
    atomic_csv(RUN_DIR / "rq5_score_macro_summary.csv", macro)

    # Figure 1: per-target AUC and AP for both raw scores.
    plot_sets = ["primary_six", "primary_five", "strict_six"]
    fig, axes = plt.subplots(len(plot_sets), 2, figsize=(13, 12), constrained_layout=True)
    for row_index, analysis_set in enumerate(plot_sets):
        subset = target_metrics[target_metrics["analysis_set"].eq(analysis_set)]
        targets = sorted(subset["target"].unique())
        positions = np.arange(len(targets))
        width = 0.36
        for metric, axis in zip(("roc_auc", "ap"), axes[row_index]):
            for offset, (score_label, color) in zip((-width / 2, width / 2), (("ipTM", "#2468A2"), ("Pinc", "#D97941"))):
                values = []
                for target in targets:
                    value = subset[
                        subset["target"].eq(target) & subset["score"].eq(score_label)
                    ][metric]
                    values.append(float(value.iloc[0]) if len(value) else np.nan)
                axis.barh(positions + offset, values, height=width, label=score_label, color=color, alpha=0.85)
            axis.set_yticks(positions)
            axis.set_yticklabels([str(target) for target in targets])
            axis.set_xlim(0, 1.02)
            axis.set_xlabel(metric.replace("_", " ").upper())
            axis.set_title(f"{analysis_set}: {metric.replace('_', ' ').upper()}")
            axis.grid(axis="x", alpha=0.2)
            axis.legend(frameon=False, loc="lower right")
    fig.suptitle("RQ5 exploratory same-sample diagnostics: raw-score discrimination", fontsize=14)
    save_figure(RUN_DIR / "RQ5_score_discrimination.png", fig)

    # Figure 2: observed rate against mean raw score, with target-equal macro panels.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    for axis, score_label in zip(axes, ("ipTM", "Pinc")):
        subset = reliability[
            reliability["analysis_set"].eq("primary_six")
            & reliability["score"].eq(score_label)
            & reliability["n"].gt(0)
        ].copy()
        for target, group in subset.groupby("target"):
            axis.plot(
                group["mean_prediction"],
                group["observed_rate"],
                marker="o",
                linewidth=1.2,
                label=str(target),
            )
        axis.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.9, label="identity")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel(f"Mean raw {score_label} in fixed bin")
        axis.set_ylabel("Observed positive rate")
        axis.set_title(f"{score_label}: fixed-bin raw-score reliability")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=7, frameon=False, loc="best")
    fig.suptitle("RQ5 exploratory raw-score reliability (primary_six)", fontsize=14)
    save_figure(RUN_DIR / "RQ5_raw_reliability.png", fig)

    audit = {
        "stage": STAGE,
        "target_metric_rows": int(len(target_metrics)),
        "bootstrap_comparison_rows": int(len(bootstrap)),
        "bootstrap_replicates_each_pair": BOOTSTRAP_REPLICATES,
        "reliability_rows_including_empty_bins": int(len(reliability)),
        "macro_rows": int(len(macro)),
        "pair_summary_rows": int(len(pooled)),
        "figures": ["RQ5_score_discrimination.png", "RQ5_raw_reliability.png"],
    }
    atomic_json(RUN_DIR / "stage_a_summary.json", audit)
    set_stage_status(STAGE, "complete", **{k: v for k, v in audit.items() if k != "stage"})
    print(f"Stage A complete: {len(target_metrics)} target-score metric rows; {len(bootstrap)} paired comparisons.")


if __name__ == "__main__":
    main()
