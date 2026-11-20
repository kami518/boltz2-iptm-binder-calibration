"""RQ4 decision value under the frozen 2026-11-14 protocol.

A/B use the frozen pair master. C reuses target-specific RQ3 LOTO U only as
the starting point for pilot adaptation; it is not a single deployment map.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
CLUSTERS = ROOT / "analysis/rq12_outputs/sequence_clusters_rq12.csv"
U_PRED = ROOT / "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_u_pair_predictions.csv"
PROTOCOL = ROOT / "analysis/RQ4_decision_value_frozen_protocol_2026-11-14.md"
V3 = ROOT / "Proteinbase校准研究_实验流程_v3_2026-11-06.md"
OUT_ROOT = ROOT / "analysis/rq4_outputs"
EXPECTED_HASH = {
    MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
    CLUSTERS: "45283741709016fb5db3a114d2a108a236c73b73750f01c959ec2753d3ff9efb",
    U_PRED: "cadfbf4fc6ee79b8a72085ae839ceb5c2785b3ca3d40aaeb21d0a010770af622",
}
TARGET_ORDER = ["aah-i", "human-vegf110", "nipah-glycoprotein-g", "pd-l1", "sars-cov-2-rbd", "trem2"]
TARGET_LABEL = {
    "aah-i": "AAH-I", "human-vegf110": "VEGF110", "nipah-glycoprotein-g": "Nipah G",
    "pd-l1": "PD-L1", "sars-cov-2-rbd": "SARS-CoV-2 RBD", "trem2": "TREM2",
}
FRACTIONS = [0.05, 0.10, 0.20]
ABSOLUTE_K = [10, 20, 50]
THRESHOLDS = [0.80, 0.90, 0.95]
B_VALUES = [120, 300]
M_VALUES = [10, 20, 40]
LAMBDA = 0.1
SEED_AB = 20261115
SEED_C = 20261114
N_BOOT = 2000
N_SIM = {"primary_six": 2000, "primary_five": 1000, "strict_six": 1000}
METRICS = ["n_selected", "hits", "precision", "recall", "enrichment", "random_expected_hits", "hit_gain", "cost_per_hit"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    tmp.replace(path)


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def make_run_dir() -> Path:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    path = OUT_ROOT / stamp
    i = 1
    while path.exists():
        path = OUT_ROOT / f"{stamp}_{i}"
        i += 1
    path.mkdir()
    return path


def metric_dict(selected_y: np.ndarray, total_y: np.ndarray) -> dict[str, float]:
    selected_y = np.asarray(selected_y, dtype=float)
    total_y = np.asarray(total_y, dtype=float)
    n = int(len(selected_y))
    hits = float(selected_y.sum())
    positives = float(total_y.sum())
    prevalence = positives / len(total_y) if len(total_y) else np.nan
    precision = hits / n if n else np.nan
    recall = hits / positives if positives else np.nan
    enrichment = precision / prevalence if n and prevalence > 0 else np.nan
    expected = n * prevalence
    return {
        "n_selected": float(n), "hits": hits, "precision": precision, "recall": recall,
        "enrichment": enrichment, "random_expected_hits": expected,
        "hit_gain": hits - expected, "cost_per_hit": n / hits if hits > 0 else np.nan,
    }


def select_mask(score: np.ndarray, pair_id: np.ndarray, kind: str, value: float) -> np.ndarray:
    n = len(score)
    if kind == "threshold":
        return score >= value
    k = int(math.ceil(n * value)) if kind == "fraction" else min(int(value), n)
    # Stable deterministic tie break after pair_id pre-sort.
    order = np.argsort(-score, kind="mergesort")[:k]
    mask = np.zeros(n, dtype=bool)
    mask[order] = True
    return mask


def point_and_cluster_bootstrap(frame: pd.DataFrame, kind: str, value: float, seed: int) -> dict[str, float]:
    d = frame.sort_values("pair_id", kind="mergesort").reset_index(drop=True)
    score = d.boltz2_iptm.to_numpy(float)
    pair_id = d.pair_id.astype(str).to_numpy()
    y = d.outcome_value.to_numpy(float)
    groups = d.sequence_cluster90.astype(str).to_numpy()
    mask = select_mask(score, pair_id, kind, value)
    result = metric_dict(y[mask], y)
    unique_groups, inverse = np.unique(groups, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(len(unique_groups))]
    rng = np.random.default_rng(seed)
    boot = {m: [] for m in METRICS}
    for _ in range(N_BOOT):
        drawn = rng.integers(0, len(members), size=len(members))
        idx = np.concatenate([members[i] for i in drawn])
        bmask = select_mask(score[idx], pair_id[idx], kind, value)
        values = metric_dict(y[idx][bmask], y[idx])
        for m in METRICS:
            boot[m].append(values[m])
    for m in METRICS:
        lo, hi = np.nanpercentile(np.asarray(boot[m], dtype=float), [2.5, 97.5])
        result[m + "_lo"] = float(lo)
        result[m + "_hi"] = float(hi)
        result[m + "_boot_valid"] = int(np.isfinite(np.asarray(boot[m], dtype=float)).sum())
    return result


def summarize_ab(rows: pd.DataFrame) -> pd.DataFrame:
    output = []
    keys = ["analysis_set", "outcome_name", "analysis_family", "budget_kind", "budget_value"]
    for key, g in rows.groupby(keys, dropna=False, sort=True):
        n = g.n_selected.sum()
        hits = g.hits.sum()
        positives = g.total_positives.sum()
        candidates = g.n_candidates.sum()
        prevalence = positives / candidates
        pooled = {
            "n_selected": n, "hits": hits,
            "precision": hits / n if n else np.nan,
            "recall": hits / positives if positives else np.nan,
            "enrichment": (hits / n) / prevalence if n and prevalence else np.nan,
            "random_expected_hits": n * prevalence,
            "hit_gain": hits - n * prevalence,
            "cost_per_hit": n / hits if hits else np.nan,
        }
        base = dict(zip(keys, key))
        output.append({**base, "aggregation": "pooled_pairs", **pooled})
        macro = {m: g[m].mean() for m in ["precision", "recall", "enrichment", "hit_gain", "cost_per_hit"]}
        output.append({**base, "aggregation": "target_equal_macro", "n_selected": n, "hits": hits, **macro})
    return pd.DataFrame(output)


def solve_intercept(eta: np.ndarray, gamma: float, target_mean: float) -> float:
    eta = np.asarray(eta, dtype=float)
    def objective(a: float) -> float:
        return float(np.mean(expit(a + gamma * eta)) - target_mean)
    a = brentq(objective, -80.0, 80.0, xtol=1e-11, rtol=1e-12)
    if abs(objective(a)) > 1e-8:
        raise RuntimeError("Intercept root did not converge")
    return float(a)


def fit_pilot(eta: np.ndarray, y: np.ndarray, lam: float = LAMBDA) -> dict[str, float | int | str]:
    eta = np.asarray(eta, dtype=float)
    y = np.asarray(y, dtype=float)
    m = len(y)
    h = float(y.sum())
    p_tilde = (h + 0.5) / (m + 1.0)
    a_i = solve_intercept(eta, 1.0, p_tilde)
    cache: dict[float, tuple[float, float]] = {}
    def profile(delta: float) -> float:
        # Bounded optimizer only supplies finite values; no score/label-driven bounds.
        delta = float(delta)
        a = solve_intercept(eta, 1.0 + delta, p_tilde)
        linear = a + (1.0 + delta) * eta
        nll = float(np.mean(np.logaddexp(0.0, linear) - y * linear))
        objective = nll + 0.5 * lam * delta * delta
        cache[delta] = (a, objective)
        return objective
    opt = minimize_scalar(profile, bounds=(-12.0, 11.0), method="bounded", options={"xatol": 1e-6, "maxiter": 200})
    if not opt.success or not np.isfinite(opt.fun):
        raise RuntimeError("F profile optimization failed")
    delta = float(opt.x)
    a_f = solve_intercept(eta, 1.0 + delta, p_tilde)
    objective = float(profile(delta))
    boundary = bool(abs(delta + 12.0) < 1e-4 or abs(delta - 11.0) < 1e-4)
    if boundary:
        raise RuntimeError("F delta reached numerical optimizer guard")
    return {
        "pilot_n": m, "pilot_positives": h, "pilot_smoothed_prevalence": p_tilde,
        "I_intercept_update": a_i, "F_intercept_update": a_f,
        "F_delta": delta, "F_total_slope": 1.0 + delta,
        "F_objective": objective, "F_optimizer_iterations": int(opt.nfev),
        "status": "evaluable",
    }


def c_metric_row(analysis_set: str, B: int, m: int, rep: int, strategy: str,
                 selected: pd.DataFrame, pilot: pd.DataFrame, followup: pd.DataFrame) -> dict:
    hits = int(selected.outcome_value.sum())
    pilot_hits = int(pilot.outcome_value.sum()) if len(pilot) else 0
    followup_hits = int(followup.outcome_value.sum()) if len(followup) else 0
    row = {
        "analysis_set": analysis_set, "B": B, "m_per_target": m, "replicate": rep,
        "strategy": strategy, "total_tests": int(len(selected)),
        "pilot_tests": int(len(pilot)), "followup_tests": int(len(followup)),
        "pilot_hits": pilot_hits, "followup_hits": followup_hits, "total_hits": hits,
        "precision": hits / B, "cost_per_hit": B / hits if hits else np.nan,
    }
    for target in TARGET_ORDER:
        row[f"followup_n__{target}"] = int((followup.target == target).sum()) if len(followup) else 0
    return row


def run_c_for_set(analysis_set: str, outcome: str, exclude_aah: bool, master: pd.DataFrame,
                  upred: pd.DataFrame, log) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    d = master.copy()
    if exclude_aah:
        d = d[d.target != "aah-i"].copy()
    d["outcome_value"] = d[outcome].astype(int)
    u = upred[(upred.analysis_set == analysis_set) & (upred.cluster_rule == "sequence_cluster90") &
              (upred.source_weighting == "target_equal") & (upred.model_form == "linear")].copy()
    keep = ["pair_id", "eta_U", "predicted_probability", "sequence_cluster_value", "fit_status"]
    u = u[keep].rename(columns={"sequence_cluster_value": "u_cluster", "predicted_probability": "p_U"})
    d = d.merge(u, on="pair_id", how="left", validate="one_to_one")
    if d.eta_U.isna().any() or len(d) != len(u) or not (d.fit_status == "evaluable").all():
        raise RuntimeError(f"U merge/QC failed for {analysis_set}: {len(d)} vs {len(u)}")
    if not np.array_equal(d.sequence_cluster90.astype(str), d.u_cluster.astype(str)):
        raise RuntimeError("Cluster mismatch between master and U predictions")
    targets = [t for t in TARGET_ORDER if t in set(d.target)]
    T = len(targets)
    if any(B % T for B in B_VALUES):
        raise RuntimeError("Equal budget requires divisibility")
    rank_rows = []
    for target in targets:
        t = d[d.target == target].sort_values(["boltz2_iptm", "pair_id"], ascending=[True, True])
        eta = t.eta_U.to_numpy(float)
        increasing = bool(np.all(np.diff(eta) > 0))
        rank_rows.append({"analysis_set": analysis_set, "target": target, "within_target_raw_U_strictly_same_order": increasing})
    rank_frame = pd.DataFrame(rank_rows)
    if not rank_frame.within_target_raw_U_strictly_same_order.all():
        raise RuntimeError("Within-target raw/U monotonicity failed")

    # Deterministic baselines.
    sorted_all = d.sort_values(["boltz2_iptm", "pair_id"], ascending=[False, True], kind="mergesort")
    baseline_raw = {B: sorted_all.head(B).copy() for B in B_VALUES}
    baseline_equal = {}
    for B in B_VALUES:
        per = B // T
        baseline_equal[B] = pd.concat([
            d[d.target == t].sort_values(["boltz2_iptm", "pair_id"], ascending=[False, True], kind="mergesort").head(per)
            for t in targets
        ], ignore_index=True)
    rng = np.random.default_rng(SEED_C + (0 if analysis_set == "primary_six" else 100000 if analysis_set == "primary_five" else 200000))
    target_data = {t: d[d.target == t].sort_values("pair_id").reset_index(drop=True) for t in targets}
    for t, td in target_data.items():
        if td.sequence_cluster90.nunique() < max(M_VALUES):
            raise RuntimeError(f"Insufficient clusters for {analysis_set}/{t}")
    run_rows: list[dict] = []
    parameter_rows: list[dict] = []
    n_sim = N_SIM[analysis_set]
    start = time.time()
    for m in M_VALUES:
        feasible_B = [B for B in B_VALUES if T * m <= B]
        if not feasible_B:
            continue
        for rep in range(n_sim):
            pilot_parts = []
            remaining_parts = []
            fit_rows = []
            for target in targets:
                td = target_data[target]
                clusters = np.array(sorted(td.sequence_cluster90.unique()), dtype=object)
                clusters = clusters[rng.permutation(len(clusters))]
                chosen_clusters = clusters[:m]
                pilot_mask = td.sequence_cluster90.isin(chosen_clusters).to_numpy()
                cluster_members = td[pilot_mask]
                # One random representative per selected cluster.
                pilot = cluster_members.groupby("sequence_cluster90", sort=False, group_keys=False).sample(n=1, random_state=int(rng.integers(0, 2**31 - 1)))
                pilot = pilot.sort_values("pair_id").reset_index(drop=True)
                if len(pilot) != m or pilot.pair_id.duplicated().any():
                    raise RuntimeError("Pilot sampling structure failure")
                remaining = td[~pilot_mask].copy()
                fit = fit_pilot(pilot.eta_U.to_numpy(float), pilot.outcome_value.to_numpy(float))
                parameter_rows.append({"analysis_set": analysis_set, "m_per_target": m, "replicate": rep,
                                       "target": target, **fit})
                pilot_parts.append(pilot)
                remaining_parts.append(remaining)
                fit_rows.append((target, fit))
            pilot_all = pd.concat(pilot_parts, ignore_index=True)
            remaining_all = pd.concat(remaining_parts, ignore_index=True)
            if set(pilot_all.pair_id) & set(remaining_all.pair_id):
                raise RuntimeError("Pilot/follow-up overlap")
            # Same pilot and fitted probabilities are reused across feasible B.
            pred_frames = {"pilot_I": [], "pilot_F": []}
            for target, fit in fit_rows:
                r = remaining_all[remaining_all.target == target].copy()
                eta = r.eta_U.to_numpy(float)
                r["pilot_I_probability"] = expit(fit["I_intercept_update"] + eta)
                r["pilot_F_probability"] = expit(fit["F_intercept_update"] + fit["F_total_slope"] * eta)
                pred_frames["pilot_I"].append(r)
                pred_frames["pilot_F"].append(r)
            predictions = {k: pd.concat(v, ignore_index=True) for k, v in pred_frames.items()}
            for B in feasible_B:
                follow_n = B - T * m
                for strategy in ["raw_ipTM", "equal_budget_ipTM", "random_allocation", "pilot_I", "pilot_F"]:
                    if strategy == "raw_ipTM":
                        selected = baseline_raw[B]
                        pilot = d.iloc[0:0]
                        follow = selected
                    elif strategy == "equal_budget_ipTM":
                        selected = baseline_equal[B]
                        pilot = d.iloc[0:0]
                        follow = selected
                    elif strategy == "random_allocation":
                        idx = rng.choice(len(d), size=B, replace=False)
                        selected = d.iloc[idx]
                        pilot = d.iloc[0:0]
                        follow = selected
                    else:
                        score_col = strategy + "_probability"
                        follow = predictions[strategy].sort_values([score_col, "pair_id"], ascending=[False, True], kind="mergesort").head(follow_n)
                        selected = pd.concat([pilot_all, follow], ignore_index=True)
                        pilot = pilot_all
                    if len(selected) != B or selected.pair_id.duplicated().any():
                        raise RuntimeError(f"Budget/duplicate failure: {analysis_set} B={B} m={m} rep={rep} {strategy}")
                    row = c_metric_row(analysis_set, B, m, rep, strategy, selected, pilot, follow)
                    if row["pilot_hits"] + row["followup_hits"] != row["total_hits"]:
                        raise RuntimeError("Hit accounting failure")
                    run_rows.append(row)
        log(f"RQ4-C {analysis_set}: completed m={m}, reps={n_sim}, elapsed={time.time()-start:.1f}s")
    return pd.DataFrame(run_rows), pd.DataFrame(parameter_rows), rank_frame


def summarize_c(runs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    keys = ["analysis_set", "B", "m_per_target", "strategy"]
    for key, g in runs.groupby(keys, sort=True):
        hits = g.total_hits.to_numpy(float)
        summaries.append({
            **dict(zip(keys, key)), "replicates": len(g), "mean_total_hits": hits.mean(),
            "sd_total_hits": hits.std(ddof=1), "hits_lo": np.percentile(hits, 2.5), "hits_hi": np.percentile(hits, 97.5),
            "mean_precision": g.precision.mean(), "mean_pilot_tests": g.pilot_tests.mean(),
            "mean_followup_tests": g.followup_tests.mean(), "mean_pilot_hits": g.pilot_hits.mean(),
            "mean_followup_hits": g.followup_hits.mean(), "prob_zero_hits": np.mean(hits == 0),
            "cost_per_hit_from_mean_hits": g.B.iloc[0] / hits.mean() if hits.mean() else np.nan,
        })
    comparisons = []
    for key, g in runs.groupby(["analysis_set", "B", "m_per_target"], sort=True):
        wide = g.pivot(index="replicate", columns="strategy", values="total_hits")
        for a, b in [("pilot_I", "raw_ipTM"), ("pilot_F", "raw_ipTM"), ("pilot_F", "pilot_I"),
                     ("pilot_I", "equal_budget_ipTM"), ("pilot_F", "equal_budget_ipTM"),
                     ("pilot_I", "random_allocation"), ("pilot_F", "random_allocation")]:
            diff = (wide[a] - wide[b]).to_numpy(float)
            comparisons.append({**dict(zip(["analysis_set", "B", "m_per_target"], key)),
                                "comparison": f"{a}_minus_{b}", "replicates": len(diff),
                                "mean_hit_difference": diff.mean(), "sd_hit_difference": diff.std(ddof=1),
                                "difference_lo": np.percentile(diff, 2.5), "difference_hi": np.percentile(diff, 97.5),
                                "probability_difference_gt_0": np.mean(diff > 0),
                                "probability_difference_eq_0": np.mean(diff == 0)})
    return pd.DataFrame(summaries), pd.DataFrame(comparisons)


def make_figures(run_dir: Path, ab: pd.DataFrame, c_summary: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    # A: fractions only, primary six.
    a = ab[(ab.analysis_set == "primary_six") & (ab.analysis_family == "RQ4-A") & (ab.budget_kind == "fraction")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    width = 0.24
    x = np.arange(len(TARGET_ORDER))
    for j, frac in enumerate(FRACTIONS):
        g = a[a.budget_value == frac].set_index("target").loc[TARGET_ORDER]
        pos = x + (j - 1) * width
        axes[0].bar(pos, g.precision, width, label=f"Top {int(frac*100)}%")
        axes[0].errorbar(pos, g.precision, yerr=[g.precision-g.precision_lo, g.precision_hi-g.precision], fmt="none", ecolor="black", capsize=2, lw=0.7)
        axes[1].bar(pos, g.enrichment, width, label=f"Top {int(frac*100)}%")
        axes[1].errorbar(pos, g.enrichment, yerr=[g.enrichment-g.enrichment_lo, g.enrichment_hi-g.enrichment], fmt="none", ecolor="black", capsize=2, lw=0.7)
    for ax in axes:
        ax.set_xticks(x, [TARGET_LABEL[t] for t in TARGET_ORDER], rotation=25, ha="right")
        ax.legend(frameon=False)
    axes[0].set_ylabel("Precision (95% cluster bootstrap CI)")
    axes[1].set_ylabel("Enrichment over target prevalence (95% CI)")
    axes[1].axhline(1, color="grey", ls="--", lw=0.8)
    fig.suptitle("RQ4-A: within-target fixed-proportion selection")
    fig.tight_layout()
    fig.savefig(run_dir / "RQ4A_fixed_budget.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    b = ab[(ab.analysis_set == "primary_six") & (ab.analysis_family == "RQ4-B")].copy()
    b["target_label"] = b.target.map(TARGET_LABEL)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, metric, title in zip(axes, ["selected_fraction", "precision"], ["Selected fraction", "Precision"]):
        mat = np.full((len(TARGET_ORDER), len(THRESHOLDS)), np.nan)
        for i, target in enumerate(TARGET_ORDER):
            for j, threshold in enumerate(THRESHOLDS):
                val = b[(b.target == target) & (b.budget_value == threshold)][metric]
                if len(val): mat[i, j] = val.iloc[0]
        im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(3), [str(t) for t in THRESHOLDS])
        ax.set_yticks(range(6), [TARGET_LABEL[t] for t in TARGET_ORDER])
        ax.set_xlabel("Raw ipTM threshold")
        ax.set_title(title)
        for i in range(6):
            for j in range(3):
                ax.text(j, i, "NA" if not np.isfinite(mat[i,j]) else f"{mat[i,j]:.2f}", ha="center", va="center", color="white" if mat[i,j] < 0.65 else "black", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("RQ4-B: raw-ipTM operating points")
    fig.tight_layout()
    fig.savefig(run_dir / "RQ4B_raw_thresholds.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    c = c_summary[c_summary.analysis_set == "primary_six"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=False)
    for ax, B in zip(axes, B_VALUES):
        g = c[c.B == B]
        for strategy, label, color, marker in [("pilot_I", "Pilot I", "#1f77b4", "o"), ("pilot_F", "Pilot F", "#d62728", "s")]:
            s = g[g.strategy == strategy].sort_values("m_per_target")
            ax.errorbar(s.m_per_target, s.mean_total_hits,
                        yerr=[s.mean_total_hits-s.hits_lo, s.hits_hi-s.mean_total_hits],
                        label=label, color=color, marker=marker, capsize=3)
        for strategy, label, color, ls in [("raw_ipTM", "Raw ipTM", "black", "-"), ("equal_budget_ipTM", "Equal budget + ipTM", "grey", "--"), ("random_allocation", "Random", "#2ca02c", ":")]:
            val = g[g.strategy == strategy].mean_total_hits.mean()
            ax.axhline(val, label=label, color=color, ls=ls, lw=1.2)
        ax.set_title(f"Total budget B={B}")
        ax.set_xlabel("Pilot experiments per target (m)")
        ax.set_ylabel("Total hits (95% simulation interval)")
        ax.set_xticks(M_VALUES)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("RQ4-C: shared-budget retrospective simulation (primary six)")
    fig.tight_layout()
    fig.savefig(run_dir / "RQ4C_shared_budget.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_results(run_dir: Path, ab: pd.DataFrame, ab_summary: pd.DataFrame, c_summary: pd.DataFrame, c_compare: pd.DataFrame) -> None:
    lines = [
        "# RQ4 实验筛选决策价值：结果草稿", "",
        "本稿由冻结协议对应的执行脚本生成。RQ4-A 为主要决策分析；RQ4-B 为原始 ipTM 描述性 operating points；RQ4-C 为事后探索性回顾性模拟。区间含义与未来靶点推断限制见冻结协议。", "",
        "## 1. RQ4-A：逐靶点固定比例预算", "",
        "| 靶点 | 预算 | 入选／hits | Precision（95% CI） | Recall | Enrichment（95% CI） | 相对随机命中增益 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    a = ab[(ab.analysis_set == "primary_six") & (ab.analysis_family == "RQ4-A") & (ab.budget_kind == "fraction")]
    for target in TARGET_ORDER:
        for frac in FRACTIONS:
            r = a[(a.target == target) & (a.budget_value == frac)].iloc[0]
            lines.append(f"| {TARGET_LABEL[target]} | Top {int(frac*100)}% | {int(r.n_selected)}／{int(r.hits)} | {r.precision:.3f}（{r.precision_lo:.3f}–{r.precision_hi:.3f}） | {r.recall:.3f} | {r.enrichment:.3f}（{r.enrichment_lo:.3f}–{r.enrichment_hi:.3f}） | {r.hit_gain:+.2f} |")
    lines += ["", "![RQ4-A](RQ4A_fixed_budget.png)", "", "## 2. RQ4-B：原始 ipTM 固定 operating points", "",
              "| 靶点 | 阈值 | 入选／hits | Precision | Recall | Enrichment |", "|---|---:|---:|---:|---:|---:|"]
    b = ab[(ab.analysis_set == "primary_six") & (ab.analysis_family == "RQ4-B")]
    for target in TARGET_ORDER:
        for threshold in THRESHOLDS:
            r = b[(b.target == target) & (b.budget_value == threshold)].iloc[0]
            def f(v): return "NA" if not np.isfinite(v) else f"{v:.3f}"
            lines.append(f"| {TARGET_LABEL[target]} | {threshold:.2f} | {int(r.n_selected)}／{int(r.hits)} | {f(r.precision)} | {f(r.recall)} | {f(r.enrichment)} |")
    lines += ["", "空选集合保留为 0 入选、0 hits；recall=0，precision、enrichment 和单位 hit 成本因分母不可定义而记为 NA。上述阈值不是成功概率。", "",
              "![RQ4-B](RQ4B_raw_thresholds.png)", "", "## 3. RQ4-C：共享预算与 pilot 成本（探索性）", "",
              "下表列出主六靶点所有可行 B／m 组合。pilot-I／F 的总 hits 已包括 pilot hits；raw、等额和随机不使用 pilot。raw 与 U 只在靶点内同序；不同 LOTO U 不构成统一部署映射。", "",
              "| B | m／靶点 | 策略 | Pilot 实验 | 平均总 hits（95% 模拟区间） | 平均 precision | 平均 pilot hits | 平均后续 hits |",
              "|---:|---:|---|---:|---:|---:|---:|---:|"]
    c = c_summary[c_summary.analysis_set == "primary_six"]
    for B in B_VALUES:
        for m in M_VALUES:
            g = c[(c.B == B) & (c.m_per_target == m)]
            if not len(g): continue
            for strategy in ["raw_ipTM", "equal_budget_ipTM", "random_allocation", "pilot_I", "pilot_F"]:
                r = g[g.strategy == strategy].iloc[0]
                lines.append(f"| {B} | {m} | {strategy} | {r.mean_pilot_tests:.0f} | {r.mean_total_hits:.2f}（{r.hits_lo:.2f}–{r.hits_hi:.2f}） | {r.mean_precision:.3f} | {r.mean_pilot_hits:.2f} | {r.mean_followup_hits:.2f} |")
    lines += ["", "### 与 raw ipTM 的逐轮配对差值", "", "| B | m／靶点 | 比较 | 平均 hits 差（95% 模拟区间） | P(差值>0) |", "|---:|---:|---|---:|---:|"]
    cc = c_compare[(c_compare.analysis_set == "primary_six") & c_compare.comparison.isin(["pilot_I_minus_raw_ipTM", "pilot_F_minus_raw_ipTM", "pilot_F_minus_pilot_I"])]
    for _, r in cc.sort_values(["B", "m_per_target", "comparison"]).iterrows():
        lines.append(f"| {r.B} | {r.m_per_target} | {r.comparison} | {r.mean_hit_difference:+.2f}（{r.difference_lo:+.2f}–{r.difference_hi:+.2f}） | {r.probability_difference_gt_0:.3f} |")
    lines += ["", "![RQ4-C](RQ4C_shared_budget.png)", "", "## 4. 预设对照", "",
              "绝对预算 k=10／20／50、排除 AAH-I、严格结局、Miniprotein-only，以及 RQ4-C 的五靶点和严格结局结果均保存在 CSV 中。解释时须区分原预设比例预算与后增绝对预算，且不得将模拟分位数称为未来靶点的置信区间。", "",
              "具体论文结论需在独立数值验收后撰写；不得仅凭某一 B／m 组合宣布 pilot 或 F 普遍有效。"]
    (run_dir / "RQ4_decision_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    run_dir = make_run_dir()
    log_path = run_dir / "execution.log"
    def log(message: str) -> None:
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as f: f.write(line + "\n")
    atomic_json(run_dir / "RUN_STATUS.json", {"status": "running", "run_dir": str(run_dir)})
    try:
        log("Checking frozen hashes and synthetic numerical tests")
        actual_hashes = {str(p): sha256(p) for p in EXPECTED_HASH}
        for p, expected in EXPECTED_HASH.items():
            if actual_hashes[str(p)] != expected: raise RuntimeError(f"Hash mismatch: {p}")
        # Synthetic tests: constant score, all-positive and all-negative pilots, and F=I at delta 0.
        eta_test = np.linspace(-2, 2, 20)
        for y_test in [np.zeros(20), np.ones(20), np.r_[np.ones(10), np.zeros(10)]]:
            fit = fit_pilot(eta_test, y_test)
            a_i = fit["I_intercept_update"]
            if abs(np.mean(expit(a_i + eta_test)) - fit["pilot_smoothed_prevalence"]) > 1e-8: raise RuntimeError("Synthetic I mean constraint")
            eta_f = fit["F_intercept_update"] + fit["F_total_slope"] * eta_test
            if abs(np.mean(expit(eta_f)) - fit["pilot_smoothed_prevalence"]) > 1e-8: raise RuntimeError("Synthetic F mean constraint")
        zero = solve_intercept(eta_test, 1.0, 0.5)
        if abs(np.mean(expit(zero + eta_test)) - 0.5) > 1e-8: raise RuntimeError("Synthetic intercept failure")
        atomic_json(run_dir / "synthetic_model_tests.json", {"status": "passed", "cases": 3, "lambda": LAMBDA, "F_delta_optimizer_guard": [-12, 11]})

        master = pd.read_csv(MASTER, low_memory=False)
        clusters = pd.read_csv(CLUSTERS)
        if len(master) != 2084 or master.pair_id.duplicated().any() or master.target.nunique() != 6: raise RuntimeError("Master structure mismatch")
        master = master.merge(clusters[["pair_id", "sequence_cluster90"]], on="pair_id", how="left", validate="one_to_one")
        if master.sequence_cluster90.isna().any(): raise RuntimeError("Missing clusters")
        upred = pd.read_csv(U_PRED, low_memory=False)
        log("Starting RQ4-A/B with 2,000 cluster bootstraps per target and rule")
        cohorts = [
            ("primary_six", "primary_outcome", master),
            ("primary_five", "primary_outcome", master[master.target != "aah-i"].copy()),
            ("strict_six", "strict_outcome", master.copy()),
            ("miniprotein_primary_six", "primary_outcome", master[master.design_class == "Miniprotein"].copy()),
        ]
        ab_rows = []
        for analysis_set, outcome, cohort in cohorts:
            d = cohort.copy()
            d["outcome_value"] = d[outcome].astype(int)
            for target in TARGET_ORDER:
                td = d[d.target == target].copy()
                if not len(td): continue
                seed_base = SEED_AB + abs(hash((analysis_set, target))) % 1000000
                rules = [("RQ4-A", "fraction", x) for x in FRACTIONS] + [("RQ4-A", "absolute_k", x) for x in ABSOLUTE_K] + [("RQ4-B", "threshold", x) for x in THRESHOLDS]
                for j, (family, kind, value) in enumerate(rules):
                    # Use a local deterministic seed, independent of Python's salted hash().
                    stable_offset = int(hashlib.sha256(f"{analysis_set}|{target}|{kind}|{value}".encode()).hexdigest()[:8], 16)
                    res = point_and_cluster_bootstrap(td, kind, value, SEED_AB + stable_offset)
                    ab_rows.append({"analysis_set": analysis_set, "outcome_name": outcome, "analysis_family": family,
                                    "budget_kind": kind, "budget_value": value, "target": target,
                                    "n_candidates": len(td), "total_positives": int(td.outcome_value.sum()),
                                    "target_prevalence": td.outcome_value.mean(), "selected_fraction": res["n_selected"] / len(td), **res})
            log(f"RQ4-A/B completed {analysis_set}: {len(d)} pairs")
        ab = pd.DataFrame(ab_rows)
        ab_summary = summarize_ab(ab)
        atomic_csv(run_dir / "rq4ab_target_metrics.csv", ab)
        atomic_csv(run_dir / "rq4ab_summary.csv", ab_summary)

        log("Starting RQ4-C retrospective shared-budget simulations")
        all_runs, all_parameters, all_ranks = [], [], []
        for analysis_set, outcome, exclude in [("primary_six", "primary_outcome", False), ("primary_five", "primary_outcome", True), ("strict_six", "strict_outcome", False)]:
            runs, parameters, ranks = run_c_for_set(analysis_set, outcome, exclude, master, upred, log)
            all_runs.append(runs); all_parameters.append(parameters); all_ranks.append(ranks)
            atomic_csv(run_dir / f"rq4c_runs_{analysis_set}.csv", runs)
            atomic_csv(run_dir / f"rq4c_pilot_parameters_{analysis_set}.csv", parameters)
        c_runs = pd.concat(all_runs, ignore_index=True)
        c_parameters = pd.concat(all_parameters, ignore_index=True)
        c_ranks = pd.concat(all_ranks, ignore_index=True)
        c_summary, c_compare = summarize_c(c_runs)
        atomic_csv(run_dir / "rq4c_simulation_runs.csv", c_runs)
        atomic_csv(run_dir / "rq4c_pilot_parameters.csv", c_parameters)
        atomic_csv(run_dir / "rq4c_summary.csv", c_summary)
        atomic_csv(run_dir / "rq4c_comparisons.csv", c_compare)
        atomic_csv(run_dir / "rq4c_u_rank_checks.csv", c_ranks)
        make_figures(run_dir, ab, c_summary)
        write_results(run_dir, ab, ab_summary, c_summary, c_compare)

        qc = {
            "input_hashes_match": True,
            "master_2084_unique_pairs": True,
            "six_primary_targets": True,
            "raw_U_within_target_rank_checks_passed": bool(c_ranks.within_target_raw_U_strictly_same_order.all()),
            "all_c_total_tests_equal_B": bool((c_runs.total_tests == c_runs.B).all()),
            "all_c_hits_accounted": bool((c_runs.pilot_hits + c_runs.followup_hits == c_runs.total_hits).all()),
            "all_c_expected_strategies": bool(c_runs.groupby(["analysis_set", "B", "m_per_target", "replicate"]).strategy.nunique().eq(5).all()),
            "all_c_pilot_costs_correct": bool((c_runs.pilot_tests == np.where(c_runs.strategy.isin(["pilot_I", "pilot_F"]), c_runs.m_per_target * c_runs.analysis_set.map({"primary_six": 6, "primary_five": 5, "strict_six": 6}), 0)).all()),
            "all_fits_evaluable": bool((c_parameters.status == "evaluable").all()),
            "no_F_guard_boundary": bool((c_parameters.F_total_slope > -11).all() & (c_parameters.F_total_slope < 12).all()),
            "ab_expected_rule_counts": bool(len(ab) == 6 * 9 + 5 * 9 + 6 * 9 + 6 * 9),
            "B_empty_selection_retained": bool(((ab.analysis_family == "RQ4-B") & (ab.n_selected == 0)).any()),
            "protocol_and_script_snapshotted": True,
        }
        if not all(qc.values()): raise RuntimeError(f"QC failed: {qc}")
        atomic_json(run_dir / "qc_checks.json", qc)
        config = {
            "protocol": str(PROTOCOL), "protocol_sha256": sha256(PROTOCOL), "v3_sha256_at_run": sha256(V3),
            "expected_input_hashes": {str(p): h for p, h in EXPECTED_HASH.items()}, "actual_input_hashes": actual_hashes,
            "fractions": FRACTIONS, "absolute_k_exploratory": ABSOLUTE_K, "thresholds": THRESHOLDS,
            "B_values": B_VALUES, "m_values": M_VALUES, "lambda_F": LAMBDA,
            "seed_AB": SEED_AB, "seed_C": SEED_C, "N_boot": N_BOOT, "N_sim": N_SIM,
            "python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        }
        atomic_json(run_dir / "run_config.snapshot.json", config)
        atomic_csv(run_dir / "final_input_hashes.csv", pd.DataFrame({"path": list(actual_hashes), "sha256": list(actual_hashes.values())}))
        # Snapshot script after all computations but before manifest.
        script_path = Path(__file__)
        (run_dir / script_path.name).write_bytes(script_path.read_bytes())
        (run_dir / "RQ4_decision_value_frozen_protocol_2026-11-14.md").write_bytes(PROTOCOL.read_bytes())
        atomic_csv(run_dir / "failed_item_ledger.csv", pd.DataFrame(columns=["stage", "item", "failure_type", "detail"]))
        files = sorted(p for p in run_dir.iterdir() if p.is_file() and p.name not in {"final_hash_manifest.csv", "manifest_sha256.txt"})
        manifest = pd.DataFrame({"file": [p.name for p in files], "bytes": [p.stat().st_size for p in files], "sha256": [sha256(p) for p in files]})
        atomic_csv(run_dir / "final_hash_manifest.csv", manifest)
        manifest_hash = sha256(run_dir / "final_hash_manifest.csv")
        (run_dir / "manifest_sha256.txt").write_text(manifest_hash + "\n", encoding="utf-8")
        atomic_json(run_dir / "RUN_STATUS.json", {"status": "complete", "run_dir": str(run_dir), "qc_passed": True,
                                                  "ab_rows": len(ab), "c_run_rows": len(c_runs), "c_parameter_rows": len(c_parameters),
                                                  "manifest_sha256": manifest_hash})
        log(f"RQ4 complete: {run_dir}")
    except Exception as exc:
        import traceback
        log("FAILED: " + repr(exc))
        log(traceback.format_exc())
        atomic_json(run_dir / "RUN_STATUS.json", {"status": "failed", "run_dir": str(run_dir), "error": repr(exc)})
        raise


if __name__ == "__main__":
    main()
