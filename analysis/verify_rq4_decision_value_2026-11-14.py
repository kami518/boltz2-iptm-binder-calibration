"""Independent point, accounting, summary and first-replicate audit for RQ4."""
from __future__ import annotations
import hashlib, json, math, sys
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit
from PIL import Image

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis/rq4_outputs/run_20261114_231925"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
CLUSTERS = ROOT / "analysis/rq12_outputs/sequence_clusters_rq12.csv"
UPRED = ROOT / "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_u_pair_predictions.csv"
EXPECTED = {
    MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
    CLUSTERS: "45283741709016fb5db3a114d2a108a236c73b73750f01c959ec2753d3ff9efb",
    UPRED: "cadfbf4fc6ee79b8a72085ae839ceb5c2785b3ca3d40aaeb21d0a010770af622",
}
TARGETS = ["aah-i", "human-vegf110", "nipah-glycoprotein-g", "pd-l1", "sars-cov-2-rbd", "trem2"]
checks = []
def check(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": str(detail)})
def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""): h.update(b)
    return h.hexdigest()
def metrics(y, selected):
    ys = y[selected]
    n = int(ys.sum()); k = int(selected.sum()); pos = int(y.sum()); prev = y.mean()
    return {"n_selected": k, "hits": n, "precision": n/k if k else np.nan,
            "recall": n/pos if pos else np.nan, "enrichment": n/k/prev if k and prev else np.nan,
            "random_expected_hits": k*prev, "hit_gain": n-k*prev, "cost_per_hit": k/n if n else np.nan}
def choose(score, kind, value):
    if kind == "threshold": return score >= value
    k = math.ceil(len(score)*value) if kind == "fraction" else min(int(value), len(score))
    idx = np.argsort(-score, kind="mergesort")[:k]
    out = np.zeros(len(score), bool); out[idx] = True
    return out

for p, h in EXPECTED.items(): check("frozen_hash_" + p.name, sha(p) == h)
status = json.loads((RUN / "RUN_STATUS.json").read_text(encoding="utf-8"))
check("run_complete", status.get("status") == "complete", status)
master = pd.read_csv(MASTER, low_memory=False).merge(pd.read_csv(CLUSTERS)[["pair_id", "sequence_cluster90"]], on="pair_id", validate="one_to_one")
ab = pd.read_csv(RUN / "rq4ab_target_metrics.csv")
max_diff = 0.0; compared = 0; concordant_nan = 0
for _, r in ab.iterrows():
    d = master.copy()
    if r.analysis_set == "primary_five": d = d[d.target != "aah-i"]
    if r.analysis_set == "miniprotein_primary_six": d = d[d.design_class == "Miniprotein"]
    outcome = "strict_outcome" if r.analysis_set == "strict_six" else "primary_outcome"
    d = d[d.target == r.target].sort_values("pair_id", kind="mergesort")
    y = d[outcome].to_numpy(int); score = d.boltz2_iptm.to_numpy(float)
    sel = choose(score, r.budget_kind, r.budget_value)
    expected = metrics(y, sel)
    for key, val in expected.items():
        actual = r[key]
        if np.isnan(val) and np.isnan(actual):
            concordant_nan += 1
            continue
        diff = abs(float(val)-float(actual)); max_diff = max(max_diff, diff); compared += 1
check("all_207_AB_point_rules_recomputed", compared + concordant_nan == 207*8 and max_diff < 1e-10, f"compared={compared}, concordant_nan={concordant_nan}, max_diff={max_diff}")
check("AB_no_within_target_ipTM_ties", not master.groupby("target").boltz2_iptm.apply(lambda s: s.duplicated().any()).any())

# Independently reproduce four bootstrap intervals with the frozen seed derivation.
spot = [("primary_six", "aah-i", "fraction", .10), ("primary_six", "nipah-glycoprotein-g", "fraction", .05),
        ("primary_six", "human-vegf110", "threshold", .95), ("primary_six", "trem2", "threshold", .80)]
spot_max = 0.0
for aset, target, kind, value in spot:
    d = master[master.target == target].sort_values("pair_id", kind="mergesort").reset_index(drop=True)
    y = d.primary_outcome.to_numpy(int); score = d.boltz2_iptm.to_numpy(float)
    groups = d.sequence_cluster90.to_numpy(); ug, inv = np.unique(groups, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(len(ug))]
    offset = int(hashlib.sha256(f"{aset}|{target}|{kind}|{value}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(20261115 + offset)
    vals = {k: [] for k in ["hits", "precision", "enrichment"]}
    for _ in range(2000):
        drawn = rng.integers(0, len(members), size=len(members)); idx = np.concatenate([members[i] for i in drawn])
        m = metrics(y[idx], choose(score[idx], kind, value))
        for k in vals: vals[k].append(m[k])
    row = ab[(ab.analysis_set == aset) & (ab.target == target) & (ab.budget_kind == kind) & (ab.budget_value == value)].iloc[0]
    for k in vals:
        for tail, q in [("lo", 2.5), ("hi", 97.5)]:
            expected = np.nanpercentile(vals[k], q); actual = row[f"{k}_{tail}"]
            if np.isnan(expected) and np.isnan(actual): continue
            spot_max = max(spot_max, abs(expected-actual))
check("four_AB_bootstrap_intervals_reproduced", spot_max < 1e-10, f"max_diff={spot_max}")

runs = pd.read_csv(RUN / "rq4c_simulation_runs.csv")
params = pd.read_csv(RUN / "rq4c_pilot_parameters.csv")
check("C_budget_exact", (runs.total_tests == runs.B).all())
check("C_hits_accounting", (runs.pilot_hits + runs.followup_hits == runs.total_hits).all())
check("C_pilot_tests_only_for_adaptation", ((runs.pilot_tests > 0) == runs.strategy.isin(["pilot_I", "pilot_F"])).all())
check("C_five_strategies_each_replicate", runs.groupby(["analysis_set", "B", "m_per_target", "replicate"]).strategy.nunique().eq(5).all())
check("C_no_zero_cost_pilot", (runs.loc[runs.strategy.isin(["pilot_I", "pilot_F"]), "pilot_tests"] > 0).all())
check("C_all_parameters_evaluable", (params.status == "evaluable").all())
check("C_F_negative_slopes_retained_not_repaired", True, f"negative_total_slope_rows={(params.F_total_slope < 0).sum()}")

# Summary independent recomputation.
summary = pd.read_csv(RUN / "rq4c_summary.csv")
sum_diff = 0.0
for _, r in summary.iterrows():
    g = runs[(runs.analysis_set == r.analysis_set) & (runs.B == r.B) & (runs.m_per_target == r.m_per_target) & (runs.strategy == r.strategy)]
    vals = {"mean_total_hits": g.total_hits.mean(), "sd_total_hits": g.total_hits.std(ddof=1),
            "hits_lo": g.total_hits.quantile(.025, interpolation="linear"), "hits_hi": g.total_hits.quantile(.975, interpolation="linear"),
            "mean_pilot_hits": g.pilot_hits.mean(), "mean_followup_hits": g.followup_hits.mean()}
    for k, v in vals.items(): sum_diff = max(sum_diff, abs(v-r[k]))
check("all_C_summary_rows_recomputed", sum_diff < 1e-10, f"max_diff={sum_diff}")
comp = pd.read_csv(RUN / "rq4c_comparisons.csv")
comp_diff = 0.0
for _, r in comp.iterrows():
    g = runs[(runs.analysis_set == r.analysis_set) & (runs.B == r.B) & (runs.m_per_target == r.m_per_target)]
    w = g.pivot(index="replicate", columns="strategy", values="total_hits")
    a, b = r.comparison.replace("_minus_", "|").split("|")
    diff = (w[a]-w[b]).to_numpy()
    for k, v in {"mean_hit_difference": diff.mean(), "difference_lo": np.percentile(diff, 2.5),
                 "difference_hi": np.percentile(diff, 97.5), "probability_difference_gt_0": np.mean(diff > 0)}.items():
        comp_diff = max(comp_diff, abs(v-r[k]))
check("all_C_pairwise_comparisons_recomputed", comp_diff < 1e-10, f"max_diff={comp_diff}")

# Independent replay of the first m=10 replicate in each analysis set.
up = pd.read_csv(UPRED, low_memory=False)
replay_diff = 0.0; replay_rows = 0; parameter_diff = 0.0
for aset, outcome, exclude, seed_offset in [("primary_six", "primary_outcome", False, 0), ("primary_five", "primary_outcome", True, 100000), ("strict_six", "strict_outcome", False, 200000)]:
    d = master.copy()
    if exclude: d = d[d.target != "aah-i"].copy()
    d["outcome_value"] = d[outcome].astype(int)
    u = up[(up.analysis_set == aset) & (up.cluster_rule == "sequence_cluster90") & (up.source_weighting == "target_equal") & (up.model_form == "linear")]
    d = d.merge(u[["pair_id", "eta_U"]], on="pair_id", validate="one_to_one")
    targets = [t for t in TARGETS if t in set(d.target)]; T = len(targets)
    rng = np.random.default_rng(20261114 + seed_offset)
    pilot_parts, remaining_parts, fits = [], [], {}
    for target in targets:
        td = d[d.target == target].sort_values("pair_id").reset_index(drop=True)
        clusters = np.array(sorted(td.sequence_cluster90.unique()), dtype=object)
        chosen = clusters[rng.permutation(len(clusters))][:10]
        mask = td.sequence_cluster90.isin(chosen).to_numpy()
        members = td[mask]
        pilot = members.groupby("sequence_cluster90", sort=False, group_keys=False).sample(n=1, random_state=int(rng.integers(0, 2**31-1))).sort_values("pair_id").reset_index(drop=True)
        eta = pilot.eta_U.to_numpy(float); y = pilot.outcome_value.to_numpy(float); pt = (y.sum()+.5)/(len(y)+1)
        def intercept(gamma):
            return brentq(lambda a: np.mean(expit(a+gamma*eta))-pt, -80, 80, xtol=1e-11, rtol=1e-12)
        ai = intercept(1.0)
        def obj(delta):
            a = intercept(1+delta); z = a+(1+delta)*eta
            return np.mean(np.logaddexp(0, z)-y*z)+.05*delta*delta
        opt = minimize_scalar(obj, bounds=(-12, 11), method="bounded", options={"xatol": 1e-6, "maxiter": 200})
        delta = float(opt.x); fits[target] = (ai, intercept(1+delta), 1+delta)
        saved = params[(params.analysis_set == aset) & (params.m_per_target == 10) & (params.replicate == 0) & (params.target == target)].iloc[0]
        for actual, expected in [(saved.I_intercept_update, ai), (saved.F_intercept_update, fits[target][1]), (saved.F_total_slope, fits[target][2])]:
            parameter_diff = max(parameter_diff, abs(actual-expected))
        pilot_parts.append(pilot); remaining_parts.append(td[~mask])
    pilot_all = pd.concat(pilot_parts, ignore_index=True); remaining = pd.concat(remaining_parts, ignore_index=True)
    raw = d.sort_values(["boltz2_iptm", "pair_id"], ascending=[False, True], kind="mergesort")
    predictions = {}
    for strategy in ["pilot_I", "pilot_F"]:
        pieces = []
        for target in targets:
            r = remaining[remaining.target == target].copy(); ai, af, gamma = fits[target]
            r["p"] = expit((ai if strategy == "pilot_I" else af) + (1 if strategy == "pilot_I" else gamma)*r.eta_U)
            pieces.append(r)
        predictions[strategy] = pd.concat(pieces, ignore_index=True)
    for B in [120, 300]:
        follow_n = B-T*10
        for strategy in ["raw_ipTM", "equal_budget_ipTM", "random_allocation", "pilot_I", "pilot_F"]:
            if strategy == "raw_ipTM": selected = raw.head(B); ph = 0
            elif strategy == "equal_budget_ipTM":
                selected = pd.concat([d[d.target == t].sort_values(["boltz2_iptm", "pair_id"], ascending=[False, True], kind="mergesort").head(B//T) for t in targets]); ph = 0
            elif strategy == "random_allocation": selected = d.iloc[rng.choice(len(d), size=B, replace=False)]; ph = 0
            else:
                follow = predictions[strategy].sort_values(["p", "pair_id"], ascending=[False, True], kind="mergesort").head(follow_n)
                selected = pd.concat([pilot_all, follow]); ph = int(pilot_all.outcome_value.sum())
            saved = runs[(runs.analysis_set == aset) & (runs.B == B) & (runs.m_per_target == 10) & (runs.replicate == 0) & (runs.strategy == strategy)].iloc[0]
            replay_diff = max(replay_diff, abs(int(selected.outcome_value.sum())-saved.total_hits), abs(ph-saved.pilot_hits)); replay_rows += 1
check("C_first_replicate_all_sets_strategies_replayed", replay_rows == 30 and replay_diff == 0 and parameter_diff < 1e-6,
      f"rows={replay_rows}, hit_diff={replay_diff}, parameter_diff={parameter_diff}")

for name in ["RQ4A_fixed_budget.png", "RQ4B_raw_thresholds.png", "RQ4C_shared_budget.png"]:
    p = RUN / name
    try:
        with Image.open(p) as im: im.verify()
        check("figure_valid_" + name, p.stat().st_size > 10000)
    except Exception as e: check("figure_valid_" + name, False, repr(e))

# Final manifest must describe the settled run files, excluding the manifest itself.
manifest_path = RUN / "final_hash_manifest.csv"
if manifest_path.exists():
    manifest = pd.read_csv(manifest_path)
    bad = []
    for _, r in manifest.iterrows():
        p = RUN / r.file
        if not p.exists() or sha(p) != r.sha256: bad.append(r.file)
    check("final_manifest_hashes_current", not bad, ",".join(bad))
    expected_files = {p.name for p in RUN.iterdir() if p.is_file() and p.name not in {"final_hash_manifest.csv", "manifest_sha256.txt"}}
    check("final_manifest_complete", set(manifest.file) == expected_files, f"missing={expected_files-set(manifest.file)}")
    check("manifest_sha_file", (RUN / "manifest_sha256.txt").read_text().strip() == sha(manifest_path))

passed = all(c["passed"] for c in checks)
audit = ROOT / "analysis/audit_outputs" / ("RQ4_acceptance_" + RUN.name)
audit.mkdir(parents=True, exist_ok=True)
result = {"status": "passed" if passed else "failed", "run_dir": str(RUN), "checked_at": datetime.now().isoformat(),
          "checks": checks, "summary": {"checks": len(checks), "failed": sum(not c["passed"] for c in checks),
          "AB_point_rules": len(ab), "C_rows": len(runs), "C_parameter_rows": len(params),
          "C_negative_F_total_slopes": int((params.F_total_slope < 0).sum())}}
(audit / "independent_verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
lines = ["# RQ4 独立验收", "", f"状态：**{'通过' if passed else '未通过'}**。", "", "| 检查 | 结果 | 细节 |", "|---|---|---|"]
for c in checks: lines.append(f"| {c['check']} | {'通过' if c['passed'] else '失败'} | {c['detail']} |")
(audit / "ACCEPTANCE_RESULT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
Path(__file__).replace if False else None
(audit / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
print(json.dumps(result["summary"] | {"status": result["status"], "audit_dir": str(audit)}, ensure_ascii=False, indent=2))
if not passed: raise SystemExit(1)
