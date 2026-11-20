"""Run the frozen RQ1/RQ2 analysis from pair_master_v3.csv and audited clusters."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
ANALYSIS = ROOT / "analysis"
OUT = ANALYSIS / "rq12_outputs"
sys.path.insert(0, str(ANALYSIS))
from rq12_core import (
    bootstrap_summary,
    clipped_logit,
    raw_metrics,
    reliability_rows,
)
from rq12_adaptive import AdaptiveHierarchicalLogistic as HierarchicalLogistic

CONFIG = json.loads((ANALYSIS / "rq12_analysis_config.json").read_text(encoding="utf-8"))
MASTER = ANALYSIS / "v3_outputs" / "pair_master_v3.csv"
CLUSTERS = OUT / "sequence_clusters_rq12.csv"
B = int(CONFIG["bootstrap_replicates"])
SEED = int(CONFIG["random_seed"])
CLIP = float(CONFIG["probability_clip"])
EDGES = np.asarray(CONFIG["reliability_bin_edges"], dtype=float)
TARGETS = ["aah-i", "human-vegf110", "nipah-glycoprotein-g", "pd-l1", "sars-cov-2-rbd", "trem2"]
LABELS = {
    "aah-i": "AAH-I (SPR)*",
    "human-vegf110": "VEGF110",
    "nipah-glycoprotein-g": "Nipah G",
    "pd-l1": "PD-L1",
    "sars-cov-2-rbd": "SARS-CoV-2 RBD",
    "trem2": "TREM2",
}
COLORS = dict(zip(TARGETS, ["#c44e52", "#4c72b0", "#55a868", "#8172b2", "#ccb974", "#64b5cd"]))
METRICS = [
    "roc_auc", "average_precision", "positive_rate", "brier", "bss_target_prevalence",
    "log_loss", "ece", "calibration_intercept", "calibration_slope",
]
SUMMARY_METRICS = ["roc_auc", "average_precision", "positive_rate", "brier", "bss_target_prevalence", "log_loss", "ece"]
EXPECTED = {
    "primary_six": (2084, 748),
    "primary_five": (1951, 729),
    "strict_six": (2084, 678),
    "strict_five": (1951, 659),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def analysis_frame(df, name):
    spec = CONFIG["analysis_sets"][name]
    selected = df.loc[df.target != "aah-i"].copy() if spec["exclude_aah_i"] else df.copy()
    selected["outcome"] = selected[spec["outcome"]].astype(int)
    return selected


def targets_for(name):
    return [t for t in TARGETS if not (CONFIG["analysis_sets"][name]["exclude_aah_i"] and t == "aah-i")]


def make_target_row(name, target, point, boot):
    return {
        "analysis_set": name,
        "target": target,
        "assay": "spr" if target == "aah-i" else "bli",
        "metric": None,
        "estimate": None,
        **boot,
    }


def make_plots(target_metrics, reliability, hierarchical_targets):
    main = target_metrics[target_metrics.analysis_set == "primary_six"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    for ax, metric, title in zip(axes, ["roc_auc", "average_precision"], ["ROC-AUC", "Average precision"]):
        y_positions = np.arange(len(TARGETS))[::-1]
        for y_pos, target in zip(y_positions, TARGETS):
            row = main[(main.target == target) & (main.metric == metric)].iloc[0]
            estimate = row.estimate
            ax.errorbar(
                estimate, y_pos,
                xerr=[[estimate - row.ci_low], [row.ci_high - estimate]],
                fmt="o", color=COLORS[target], capsize=4, markersize=7,
            )
            if metric == "average_precision":
                baseline = main[(main.target == target) & (main.metric == "positive_rate")].iloc[0].estimate
                ax.scatter(baseline, y_pos, marker="d", color="gray", s=35, zorder=3)
        ax.set_yticks(y_positions, [LABELS[t] for t in TARGETS])
        ax.set_xlim(0, 1.02)
        ax.set_xlabel(title)
        ax.grid(axis="x", alpha=0.25)
        if metric == "roc_auc":
            ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
        else:
            ax.scatter([], [], marker="d", color="gray", label="Target prevalence")
            ax.legend(loc="lower left")
    fig.suptitle("RQ1 | Within-target discrimination (90% sequence-cluster bootstrap, 95% CI)")
    fig.text(0.01, 0.01, "* AAH-I uses SPR; its 19 positives include 18 selectively BLI-confirmed designs.", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    fig.savefig(OUT / "RQ1_target_discrimination.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, target in zip(axes.flat, TARGETS):
        subset = df_all[df_all.target == target]
        negative = subset.loc[subset.primary_outcome == 0, "boltz2_iptm"]
        positive = subset.loc[subset.primary_outcome == 1, "boltz2_iptm"]
        bins = np.linspace(0, 1, 26)
        ax.hist(negative, bins=bins, density=True, alpha=0.55, color="#4c72b0", label=f"Negative (n={len(negative)})")
        ax.hist(positive, bins=bins, density=True, alpha=0.55, color="#c44e52", label=f"Positive (n={len(positive)})")
        ax.set_title(LABELS[target])
        ax.set_xlim(0, 1)
        ax.set_xlabel("Raw Boltz-2 ipTM")
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)
    fig.suptitle("RQ1 | Raw ipTM distributions by observed outcome")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "RQ1_ipTM_distributions.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, target in zip(axes.flat, TARGETS):
        subset = reliability[(reliability.analysis_set == "primary_six") & (reliability.target == target) & (reliability.n > 0)]
        ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=1, label="Ideal")
        ax.errorbar(
            subset.mean_iptm, subset.observed_rate,
            yerr=[subset.observed_rate - subset.wilson_low, subset.wilson_high - subset.observed_rate],
            fmt="o-", color=COLORS[target], capsize=3, label="Fixed-bin observed rate",
        )
        for _, row in subset.iterrows():
            y_offset = [7, 17, -14][int(row.bin_index) % 3]
            ax.annotate(f"n={int(row.n)}", (row.mean_iptm, row.observed_rate), textcoords="offset points", xytext=(3, y_offset), fontsize=7)
        ax.set_title(LABELS[target])
        ax.set_xlim(0, 1.02)
        ax.set_ylim(-0.02, 1.05)
        ax.set_xlabel("Mean raw ipTM in bin")
        ax.set_ylabel("Observed binder fraction")
        ax.grid(alpha=0.2)
    fig.suptitle("RQ2 | Raw ipTM reliability; Wilson 95% intervals")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "RQ2_raw_reliability.png", dpi=180)
    plt.close(fig)

    main_h = hierarchical_targets[hierarchical_targets.analysis_set == "primary_six"].set_index("target")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    for ax, unpooled, pooled, title, ideal in [
        (axes[0], "alpha_unpooled", "alpha_partial", "Calibration intercept", 0.0),
        (axes[1], "beta_unpooled", "beta_partial", "Calibration slope", 1.0),
    ]:
        y_positions = np.arange(len(TARGETS))[::-1]
        for y_pos, target in zip(y_positions, TARGETS):
            row = main_h.loc[target]
            ax.annotate("", xy=(row[pooled], y_pos), xytext=(row[unpooled], y_pos), arrowprops={"arrowstyle": "->", "color": COLORS[target], "alpha": 0.65})
            ax.scatter(row[unpooled], y_pos, facecolor="white", edgecolor=COLORS[target], s=65, label="Unpooled" if target == TARGETS[0] else None)
            ax.scatter(row[pooled], y_pos, color=COLORS[target], s=65, label="Partial pooling" if target == TARGETS[0] else None)
        ax.set_yticks(y_positions, [LABELS[t] for t in TARGETS])
        ax.axvline(ideal, color="black", linestyle="--", linewidth=1, label="Ideal raw calibration")
        ax.set_xlabel(title)
        ax.grid(axis="x", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3)
    fig.suptitle("RQ2 | Target-specific estimates before and after partial pooling")
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    fig.savefig(OUT / "RQ2_partial_pooling.png", dpi=180)
    plt.close(fig)


started = time.time()
OUT.mkdir(parents=True, exist_ok=True)
df_all = pd.read_csv(MASTER)
clusters = pd.read_csv(CLUSTERS)
if len(df_all) != 2084 or df_all.pair_id.nunique() != 2084:
    raise RuntimeError("Frozen master quality gate failed")
if len(clusters) != 2084 or clusters.pair_id.nunique() != 2084:
    raise RuntimeError("Cluster mapping quality gate failed")
df_all = df_all.merge(clusters[["pair_id", "sequence_cluster90", "sequence_cluster80"]], on="pair_id", how="left", validate="one_to_one")
if df_all[["boltz2_iptm", "primary_outcome", "strict_outcome", "sequence_cluster90"]].isna().any().any():
    raise RuntimeError("Missing scores, outcomes, or cluster assignments")
for name, (expected_n, expected_positive) in EXPECTED.items():
    selected = analysis_frame(df_all, name)
    if len(selected) != expected_n or int(selected.outcome.sum()) != expected_positive:
        raise RuntimeError(f"Analysis-set gate failed for {name}: {len(selected)}/{selected.outcome.sum()}")

point_cache = {}
boot_cache = {}
bootstrap_audit = []
for target_index, target in enumerate(TARGETS):
    target_df = df_all[df_all.target == target].copy()
    scores = target_df.boltz2_iptm.to_numpy(dtype=float)
    outcomes = {
        "primary_outcome": target_df.primary_outcome.to_numpy(dtype=int),
        "strict_outcome": target_df.strict_outcome.to_numpy(dtype=int),
    }
    cluster_names = target_df.sequence_cluster90.drop_duplicates().tolist()
    cluster_rows = [np.flatnonzero(target_df.sequence_cluster90.to_numpy() == cluster) for cluster in cluster_names]
    rng = np.random.default_rng(SEED + target_index * 100003)
    for outcome_name, y in outcomes.items():
        point = raw_metrics(scores, y, EDGES, CLIP)
        if point is None:
            raise RuntimeError(f"Point fit failed: {target}/{outcome_name}")
        point_cache[(target, outcome_name)] = point
        boot_cache[(target, outcome_name)] = {metric: np.full(B, np.nan) for metric in point}
    failures = {"primary_outcome": 0, "strict_outcome": 0}
    valid_b = 0
    attempts = 0
    while valid_b < B:
        attempts += 1
        if attempts > B * 20:
            raise RuntimeError(f"Bootstrap attempt cap reached for {target}: {attempts}")
        drawn = rng.integers(0, len(cluster_names), size=len(cluster_names))
        sampled = np.concatenate([cluster_rows[c] for c in drawn])
        results = {}
        for outcome_name, y in outcomes.items():
            result = raw_metrics(scores[sampled], y[sampled], EDGES, CLIP)
            results[outcome_name] = result
            if result is None:
                failures[outcome_name] += 1
        if any(result is None for result in results.values()):
            continue
        for outcome_name, result in results.items():
            for metric, value in result.items():
                boot_cache[(target, outcome_name)][metric][valid_b] = value
        valid_b += 1
    audit = {
        "target": target,
        "clusters": len(cluster_names),
        "requested_bootstrap": B,
        "valid_bootstrap": valid_b,
        "attempts": attempts,
        **failures,
    }
    bootstrap_audit.append(audit)
    print(f"bootstrap_target={target} clusters={len(cluster_names)} valid={B} elapsed={time.time()-started:.1f}s", flush=True)

target_rows = []
reliability_output = []
for name in CONFIG["analysis_sets"]:
    outcome_name = CONFIG["analysis_sets"][name]["outcome"]
    selected = analysis_frame(df_all, name)
    for target in targets_for(name):
        target_df = selected[selected.target == target]
        point = point_cache[(target, outcome_name)]
        scores = target_df.boltz2_iptm.to_numpy(dtype=float)
        y = target_df.outcome.to_numpy(dtype=int)
        unique_scores, counts = np.unique(scores, return_counts=True)
        for metric in METRICS:
            boot = bootstrap_summary(boot_cache[(target, outcome_name)][metric])
            se = point.get(metric + "_se", np.nan)
            row = make_target_row(name, target, point, boot)
            row.update({
                "metric": metric,
                "estimate": point[metric],
                "model_standard_error": se,
                "model_wald_low": point[metric] - 1.95996398454 * se if np.isfinite(se) else np.nan,
                "model_wald_high": point[metric] + 1.95996398454 * se if np.isfinite(se) else np.nan,
                "n": len(target_df),
                "positives": int(y.sum()),
                "negatives": int(len(y) - y.sum()),
                "iptm_unique_values": len(unique_scores),
                "iptm_max_tie": int(counts.max()),
                "iptm_positive_median": float(np.median(scores[y == 1])),
                "iptm_negative_median": float(np.median(scores[y == 0])),
                "sequence_clusters": int(target_df.sequence_cluster90.nunique()),
            })
            target_rows.append(row)
        for bin_row in reliability_rows(scores, y, EDGES):
            reliability_output.append({"analysis_set": name, "target": target, **bin_row})

target_metrics = pd.DataFrame(target_rows)
target_metrics.to_csv(OUT / "rq12_target_metrics.csv", index=False)
target_metrics[target_metrics.metric.isin(["roc_auc", "average_precision", "positive_rate"])].to_csv(OUT / "rq1_target_metrics.csv", index=False)
target_metrics[target_metrics.metric.isin(["brier", "bss_target_prevalence", "log_loss", "ece", "calibration_intercept", "calibration_slope"])].to_csv(OUT / "rq2_target_metrics.csv", index=False)
reliability = pd.DataFrame(reliability_output)
reliability.to_csv(OUT / "rq2_reliability_bins.csv", index=False)

summary_rows = []
for name in CONFIG["analysis_sets"]:
    outcome_name = CONFIG["analysis_sets"][name]["outcome"]
    selected = analysis_frame(df_all, name)
    target_list = targets_for(name)
    for weighting in ["pair_weighted", "target_equal"]:
        weights = np.array([len(selected[selected.target == t]) if weighting == "pair_weighted" else 1.0 for t in target_list], dtype=float)
        for metric in SUMMARY_METRICS:
            if metric == "bss_target_prevalence":
                point_brier = np.average([point_cache[(t, outcome_name)]["brier"] for t in target_list], weights=weights)
                point_base = np.average([point_cache[(t, outcome_name)]["baseline_brier"] for t in target_list], weights=weights)
                estimate = 1.0 - point_brier / point_base
                boot_brier = np.average(np.array([boot_cache[(t, outcome_name)]["brier"] for t in target_list]), axis=0, weights=weights)
                boot_base = np.average(np.array([boot_cache[(t, outcome_name)]["baseline_brier"] for t in target_list]), axis=0, weights=weights)
                values = 1.0 - boot_brier / boot_base
            else:
                estimate = np.average([point_cache[(t, outcome_name)][metric] for t in target_list], weights=weights)
                values = np.average(np.array([boot_cache[(t, outcome_name)][metric] for t in target_list]), axis=0, weights=weights)
            summary_rows.append({
                "analysis_set": name,
                "weighting": weighting,
                "metric": metric,
                "estimate": float(estimate),
                **bootstrap_summary(values),
            })
summary_metrics = pd.DataFrame(summary_rows)
summary_metrics.to_csv(OUT / "rq12_summary_metrics.csv", index=False)

hierarchical_results = {}
hierarchical_target_rows = []
hierarchical_hyper_rows = []
stability_result = {"status": "not_run"}
for name in CONFIG["analysis_sets"]:
    outcome_name = CONFIG["analysis_sets"][name]["outcome"]
    selected = analysis_frame(df_all, name)
    target_list = targets_for(name)
    groups = []
    for target in target_list:
        target_df = selected[selected.target == target]
        groups.append((clipped_logit(target_df.boltz2_iptm, CLIP), target_df.outcome.to_numpy(dtype=int)))
    model = HierarchicalLogistic(groups, CONFIG["hierarchical_model"])
    result = model.fit()
    hierarchical_results[name] = result
    for target, posterior in zip(target_list, result["targets"]):
        unpooled = point_cache[(target, outcome_name)]
        hierarchical_target_rows.append({
            "analysis_set": name,
            "target": target,
            "model_type": result["model_type"],
            "alpha_unpooled": unpooled["calibration_intercept"],
            "beta_unpooled": unpooled["calibration_slope"],
            **posterior,
            "alpha_shift_from_unpooled": posterior["alpha_partial"] - unpooled["calibration_intercept"],
            "beta_shift_from_unpooled": posterior["beta_partial"] - unpooled["calibration_slope"],
        })
    hierarchical_hyper_rows.append({
        "analysis_set": name,
        "model_type": result["model_type"],
        "fallback_reason": result["fallback_reason"],
        "converged": result["converged"],
        "optimizer_message": result["optimizer_message"],
        "information_positive_definite": result["information_positive_definite"],
        "information_min_eigenvalue": result["information_min_eigenvalue"],
        "quadrature_nodes": result["quadrature_nodes"],
        **result["hyperparameters"],
    })
    print(f"hierarchical={name} type={result['model_type']} converged={result['converged']} pd={result['information_positive_definite']} elapsed={time.time()-started:.1f}s", flush=True)
    if name == "primary_six":
        stable_model = HierarchicalLogistic(groups, CONFIG["hierarchical_model"], nodes=CONFIG["hierarchical_model"]["quadrature_stability_nodes"])
        if result["model_type"] == "random_intercept_and_slope":
            stable = stable_model.fit(initial=result["parameters"])
        else:
            stable = stable_model._fit_random_intercept(result["parameters"][:2])
        parameter_keys = ["mu_intercept", "mu_slope", "tau_intercept", "tau_slope"]
        differences = {
            key: abs(result["hyperparameters"][key] - stable["hyperparameters"][key])
            for key in parameter_keys
            if result["hyperparameters"][key] is not None and stable["hyperparameters"][key] is not None
        }
        target_differences = []
        for a, b in zip(result["targets"], stable["targets"]):
            target_differences.extend([abs(a["alpha_partial"] - b["alpha_partial"]), abs(a["beta_partial"] - b["beta_partial"])])
        stability_result = {
            "status": "complete",
            "base_nodes": result["quadrature_nodes"],
            "check_nodes": stable["quadrature_nodes"],
            "same_model_type": result["model_type"] == stable["model_type"],
            "check_converged": stable["converged"],
            "check_information_positive_definite": stable["information_positive_definite"],
            "hyperparameter_absolute_differences": differences,
            "max_target_estimate_absolute_difference": max(target_differences),
            "tolerance": 0.02,
            "stable_result": stable,
        }

hierarchical_targets = pd.DataFrame(hierarchical_target_rows)
hierarchical_targets.to_csv(OUT / "rq2_hierarchical_targets.csv", index=False)
hierarchical_hyper = pd.DataFrame(hierarchical_hyper_rows)
hierarchical_hyper.to_csv(OUT / "rq2_hierarchical_hyperparameters.csv", index=False)
(OUT / "rq2_hierarchical_full_results.json").write_text(json.dumps({"models": hierarchical_results, "primary_six_stability": stability_result}, indent=2), encoding="utf-8")

make_plots(target_metrics, reliability, hierarchical_targets)

main_hierarchical = hierarchical_results["primary_six"]
stability_pass = (
    stability_result.get("status") == "complete"
    and stability_result["same_model_type"]
    and stability_result["check_converged"]
    and stability_result["check_information_positive_definite"]
    and stability_result["max_target_estimate_absolute_difference"] < stability_result["tolerance"]
    and max(stability_result["hyperparameter_absolute_differences"].values()) < stability_result["tolerance"]
)
quality_gates = {
    "frozen_master_2084_748_678": len(df_all) == 2084 and int(df_all.primary_outcome.sum()) == 748 and int(df_all.strict_outcome.sum()) == 678,
    "all_requested_bootstraps_valid": all(a["valid_bootstrap"] == B for a in bootstrap_audit),
    "summary_bootstraps_valid": bool((summary_metrics.valid_bootstrap == B).all()),
    "main_hierarchical_converged": bool(main_hierarchical["converged"]),
    "main_hierarchical_information_pd": bool(main_hierarchical["information_positive_definite"]),
    "all_hierarchical_models_converged_pd": bool(hierarchical_hyper.converged.all() and hierarchical_hyper.information_positive_definite.all()),
    "primary_six_quadrature_stability": bool(stability_pass),
    "figures_created": all((OUT / name).exists() for name in ["RQ1_target_discrimination.png", "RQ1_ipTM_distributions.png", "RQ2_raw_reliability.png", "RQ2_partial_pooling.png"]),
}
manifest = {
    "status": "complete" if all(quality_gates.values()) else "failed_quality_gate",
    "started_unix": started,
    "runtime_seconds": time.time() - started,
    "input_master": str(MASTER),
    "input_master_sha256": sha256(MASTER),
    "cluster_mapping": str(CLUSTERS),
    "cluster_mapping_sha256": sha256(CLUSTERS),
    "config": CONFIG,
    "quality_gates": quality_gates,
    "bootstrap_audit": bootstrap_audit,
    "hierarchical_stability": {k: v for k, v in stability_result.items() if k != "stable_result"},
    "software": {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
    },
    "interpretation_boundary": "Raw ipTM metrics are in-sample diagnostics of the frozen score, not out-of-fold prediction performance. Target posterior SDs are conditional on estimated hyperparameters; hyperparameter uncertainty is reported separately.",
}
(OUT / "rq12_run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps({"status": manifest["status"], "runtime_seconds": manifest["runtime_seconds"], "quality_gates": quality_gates}, indent=2), flush=True)
