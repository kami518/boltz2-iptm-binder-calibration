from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
INPUT = ROOT / "data" / "proteinbase_latest.csv"
OUT = ROOT / "analysis" / "audit_outputs"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(value):
    if value is True or value == 1:
        return 1
    if value is False or value == 0:
        return 0
    return np.nan


def auc(y, score):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = rankdata(score, method="average")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def q(series, p):
    return float(series.quantile(p)) if len(series) else np.nan


print("Reading CSV...", flush=True)
df = pd.read_csv(INPUT)
print(f"Design rows: {len(df)}", flush=True)

required = ["id", "name", "sequence", "author", "designMethod", "evaluations"]
missing_columns = [c for c in required if c not in df.columns]
if missing_columns:
    raise ValueError(f"Missing columns: {missing_columns}")

records = []
parse_errors = []
for row_index, row in df.iterrows():
    try:
        evaluations = json.loads(row["evaluations"]) if pd.notna(row["evaluations"]) else []
        if not isinstance(evaluations, list):
            raise TypeError("evaluations is not a list")
    except Exception as exc:
        parse_errors.append({"row": row_index, "id": row["id"], "error": str(exc)})
        continue
    for eval_index, ev in enumerate(evaluations):
        value = ev.get("value")
        records.append({
            "row_index": row_index,
            "eval_index": eval_index,
            "design_id": row["id"],
            "design_name": row["name"],
            "sequence": row["sequence"],
            "author": row["author"],
            "designMethod": row["designMethod"],
            "type": ev.get("type"),
            "metric": ev.get("metric"),
            "method": ev.get("method"),
            "target": ev.get("target"),
            "unit": ev.get("unit"),
            "valueType": ev.get("valueType"),
            "measuredAt": ev.get("measuredAt"),
            "value": value if isinstance(value, (str, int, float, bool)) or value is None else json.dumps(value, sort_keys=True),
        })

print(f"Evaluation records: {len(records)}", flush=True)
ev = pd.DataFrame(records)
ev["target_norm"] = ev["target"].fillna("").astype(str).str.strip().str.lower()
ev.loc[ev["target_norm"] == "", "target_norm"] = np.nan

# One row per distinct metric value per design/target. This does not assert independence.
numeric = ev[(ev["type"] == "computational") & (ev["valueType"] == "numeric")].copy()
numeric["numeric_value"] = pd.to_numeric(numeric["value"], errors="coerce")

# Design class and other design-level annotations.
design_class = (
    numeric[numeric["metric"] == "design_class"]
    if False
    else ev[ev["metric"] == "design_class"].copy()
)
class_conflicts = design_class.groupby("design_id")["value"].nunique(dropna=True)
class_map = design_class.dropna(subset=["value"]).drop_duplicates(["design_id", "value"]).groupby("design_id")["value"].first()

binding = ev[(ev["type"] == "experimental") & (ev["metric"] == "binding")].copy()
binding["y_value"] = binding["value"].map(as_bool)
binding = binding.dropna(subset=["design_id", "target_norm", "y_value"])
binding["y_value"] = binding["y_value"].astype(int)
binding["method_norm"] = binding["method"].fillna("unknown").astype(str).str.lower()

# Exact duplicate rows: same design, target, assay, value and timestamp.
dedup_keys = ["design_id", "target_norm", "method_norm", "measuredAt", "y_value"]
binding_unique = binding.drop_duplicates(dedup_keys).copy()

pair_rows = []
for (design_id, target), raw_group in binding.groupby(["design_id", "target_norm"], dropna=False):
    unique_group = binding_unique[
        (binding_unique["design_id"] == design_id) & (binding_unique["target_norm"] == target)
    ]
    methods = sorted(raw_group["method_norm"].unique())
    timestamps = raw_group["measuredAt"].dropna().unique()
    raw_y = raw_group["y_value"].to_numpy()
    unique_y = unique_group["y_value"].to_numpy()
    pair_rows.append({
        "design_id": design_id,
        "target": target,
        "author": raw_group["author"].iloc[0],
        "designMethod": raw_group["designMethod"].iloc[0],
        "sequence": raw_group["sequence"].iloc[0],
        "design_class": class_map.get(design_id, np.nan),
        "methods": ";".join(methods),
        "n_methods": len(methods),
        "n_binding_records": len(raw_group),
        "n_binding_records_deduplicated": len(unique_group),
        "n_unique_timestamps": len(timestamps),
        "positive_records_raw": int(raw_y.sum()),
        "positive_fraction_raw": float(raw_y.mean()),
        "positive_fraction_deduplicated": float(unique_y.mean()) if len(unique_y) else np.nan,
        "y_raw_majority": int(raw_y.mean() >= 0.5),
        "y_deduplicated_majority": int(unique_y.mean() >= 0.5) if len(unique_y) else np.nan,
        "y_raw_all": int(raw_y.all()),
        "y_raw_any": int(raw_y.any()),
        "label_discordant_raw": int(raw_y.min() != raw_y.max()),
        "first_measuredAt": raw_group["measuredAt"].min(),
        "last_measuredAt": raw_group["measuredAt"].max(),
    })

pairs = pd.DataFrame(pair_rows)

# Attach one unique numeric value per design/target/metric; audit multiplicities separately.
score_metrics = [
    "boltz2_iptm", "boltz2_complex_plddt", "boltz2_ipsae", "boltz2_min_ipsae",
    "boltz2_lis", "boltz2_complex_pde", "boltz2_pdockq", "boltz2_pdockq2",
    "esmfold_plddt", "proteinmpnn_score", "shape_complimentarity_boltz2_binder_ss",
]
score_multiplicity = []
for metric in score_metrics:
    sm = numeric[(numeric["metric"] == metric) & numeric["numeric_value"].notna()].copy()
    key = ["design_id", "target_norm"]
    unique_values = sm.drop_duplicates(key + ["numeric_value"])
    counts = unique_values.groupby(key, dropna=False)["numeric_value"].nunique()
    score_multiplicity.append({
        "metric": metric,
        "numeric_records": len(sm),
        "unique_design_target_metric_values": len(unique_values),
        "pairs_with_multiple_distinct_values": int((counts > 1).sum()),
        "max_distinct_values_per_pair": int(counts.max()) if len(counts) else 0,
    })
    collapsed = unique_values.groupby(key, dropna=False)["numeric_value"].mean().rename(metric).reset_index()
    if metric in ["esmfold_plddt", "proteinmpnn_score"]:
        # These are design-level metrics: target is expected to be missing.
        design_values = sm.drop_duplicates(["design_id", "numeric_value"]).groupby("design_id")["numeric_value"].mean()
        pairs[metric] = pairs["design_id"].map(design_values)
    else:
        pairs = pairs.merge(collapsed.rename(columns={"target_norm": "target"}), on=["design_id", "target"], how="left")

pairs["iptm_covered"] = pairs["boltz2_iptm"].notna()

# Summaries by target, both for all binders and miniproteins only.
def target_summary(frame, suffix=""):
    rows = []
    for target, g in frame.groupby("target", dropna=False):
        scored = g[g["boltz2_iptm"].notna()].copy()
        y = scored["y_deduplicated_majority"].to_numpy()
        s = scored["boltz2_iptm"].to_numpy()
        pos = int(np.nansum(y)) if len(y) else 0
        row = {
            "target": target,
            "n_pairs": len(g),
            "n_designs": g["design_id"].nunique(),
            "n_positive_pairs_raw_majority": int(g["y_raw_majority"].sum()),
            "n_positive_pairs_deduplicated_majority": int(g["y_deduplicated_majority"].sum()),
            "hit_rate_deduplicated": float(g["y_deduplicated_majority"].mean()),
            "n_methods": g["methods"].nunique(),
            "methods": ";".join(sorted(g["methods"].unique())),
            "n_authors": g["author"].nunique(),
            "authors": ";".join(sorted(g["author"].dropna().astype(str).unique())),
            "n_design_classes": g["design_class"].nunique(),
            "design_classes": ";".join(sorted(g["design_class"].dropna().astype(str).unique())),
            "n_designMethods_nonmissing": int(g["designMethod"].notna().sum()),
            "n_with_iptm": len(scored),
            "iptm_coverage": float(len(scored) / len(g)) if len(g) else np.nan,
            "n_positive_with_iptm": pos,
            "hit_rate_with_iptm": float(np.nanmean(y)) if len(y) else np.nan,
            "iptm_min": float(scored["boltz2_iptm"].min()) if len(scored) else np.nan,
            "iptm_q05": q(scored["boltz2_iptm"], 0.05),
            "iptm_q25": q(scored["boltz2_iptm"], 0.25),
            "iptm_median": q(scored["boltz2_iptm"], 0.50),
            "iptm_q75": q(scored["boltz2_iptm"], 0.75),
            "iptm_q95": q(scored["boltz2_iptm"], 0.95),
            "iptm_max": float(scored["boltz2_iptm"].max()) if len(scored) else np.nan,
            "iptm_unique_values": int(scored["boltz2_iptm"].nunique()),
            "auc_iptm_descriptive": auc(y, s) if len(y) else np.nan,
            "raw_iptm_brier_descriptive": float(np.nanmean((s - y) ** 2)) if len(y) else np.nan,
        }
        for threshold in [0.8, 0.9, 0.95]:
            selected = scored[scored["boltz2_iptm"] >= threshold]
            key = str(threshold).replace(".", "p")
            row[f"n_ge_{key}"] = len(selected)
            row[f"precision_ge_{key}"] = float(selected["y_deduplicated_majority"].mean()) if len(selected) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_with_iptm", "n_pairs"], ascending=False)

targets_all = target_summary(pairs)
targets_mini = target_summary(pairs[pairs["design_class"] == "Miniprotein"])

# Metric coverage and value-range audit on the paired analysis table.
metric_audit = []
for metric in score_metrics:
    s = pd.to_numeric(pairs[metric], errors="coerce")
    present = s.dropna()
    metric_audit.append({
        "metric": metric,
        "n_pairs_present": len(present),
        "coverage_of_all_pairs": float(len(present) / len(pairs)),
        "min": float(present.min()) if len(present) else np.nan,
        "q01": q(present, 0.01),
        "q05": q(present, 0.05),
        "median": q(present, 0.50),
        "q95": q(present, 0.95),
        "max": float(present.max()) if len(present) else np.nan,
        "n_between_0_and_1": int(present.between(0, 1).sum()),
        "n_above_1": int((present > 1).sum()),
        "n_unique": int(present.nunique()),
    })
metric_audit = pd.DataFrame(metric_audit)

# Exact sequence redundancy and multi-target structure.
seq_counts = df.groupby("sequence")["id"].nunique().sort_values(ascending=False)
multi_target_counts = pairs.groupby("design_id")["target"].nunique()

# Target-method and target-author contingency tables.
target_method = (
    pairs.assign(n=1)
    .pivot_table(index="target", columns="methods", values="n", aggfunc="sum", fill_value=0)
)
target_author = (
    pairs.assign(n=1)
    .pivot_table(index="target", columns="author", values="n", aggfunc="sum", fill_value=0)
)

# Binding strength is audited separately because it may be derived rather than independent.
strength = ev[(ev["type"] == "experimental") & (ev["metric"] == "binding_strength")]

summary = {
    "input": {
        "path": str(INPUT),
        "bytes": INPUT.stat().st_size,
        "sha256": sha256(INPUT),
        "design_rows": len(df),
        "unique_id": int(df["id"].nunique()),
        "unique_name": int(df["name"].nunique()),
        "unique_sequence": int(df["sequence"].nunique()),
        "authors": df["author"].value_counts(dropna=False).to_dict(),
        "designMethod_missing": int(df["designMethod"].isna().sum()),
        "designMethod_values": df["designMethod"].value_counts(dropna=False).to_dict(),
        "evaluation_parse_errors": parse_errors,
        "evaluation_records": len(ev),
        "evaluation_records_by_type": ev["type"].value_counts(dropna=False).to_dict(),
        "experimental_records_by_method": ev[ev["type"] == "experimental"]["method"].value_counts(dropna=False).to_dict(),
    },
    "binding": {
        "raw_records": len(binding),
        "deduplicated_records": len(binding_unique),
        "exact_duplicate_records_removed": len(binding) - len(binding_unique),
        "pairs": len(pairs),
        "positive_pairs_raw_majority": int(pairs["y_raw_majority"].sum()),
        "positive_pairs_deduplicated_majority": int(pairs["y_deduplicated_majority"].sum()),
        "pairs_changed_by_exact_deduplication": int((pairs["y_raw_majority"] != pairs["y_deduplicated_majority"]).sum()),
        "pairs_with_discordant_raw_labels": int(pairs["label_discordant_raw"].sum()),
        "pairs_with_only_one_unique_timestamp": int((pairs["n_unique_timestamps"] == 1).sum()),
        "pairs_with_more_records_than_unique_timestamps": int((pairs["n_binding_records"] > pairs["n_unique_timestamps"]).sum()),
        "binding_strength_records": len(strength),
        "binding_strength_values": strength["value"].value_counts(dropna=False).to_dict(),
        "n_targets": int(pairs["target"].nunique()),
        "multi_target_designs": int((multi_target_counts > 1).sum()),
        "max_targets_per_design": int(multi_target_counts.max()),
        "sequences_shared_by_multiple_ids": int((seq_counts > 1).sum()),
        "ids_on_sequences_shared_by_multiple_ids": int(seq_counts[seq_counts > 1].sum()),
    },
    "design_class": {
        "pairs_by_class": pairs["design_class"].value_counts(dropna=False).to_dict(),
        "designs_by_class": pairs.drop_duplicates("design_id")["design_class"].value_counts(dropna=False).to_dict(),
        "designs_with_conflicting_classes": int((class_conflicts > 1).sum()),
        "miniprotein_pairs": int((pairs["design_class"] == "Miniprotein").sum()),
        "miniprotein_designs": int(pairs.loc[pairs["design_class"] == "Miniprotein", "design_id"].nunique()),
        "miniprotein_targets": int(pairs.loc[pairs["design_class"] == "Miniprotein", "target"].nunique()),
    },
    "confounding": {
        "targets_with_one_method_combination": int((targets_all["n_methods"] == 1).sum()),
        "targets_with_one_author": int((targets_all["n_authors"] == 1).sum()),
        "targets_with_one_design_class": int((targets_all["n_design_classes"] == 1).sum()),
        "targets_with_iptm_and_at_least_50_pairs_and_10_positives": int(((targets_all["n_with_iptm"] >= 50) & (targets_all["n_positive_with_iptm"] >= 10)).sum()),
        "miniprotein_targets_with_iptm_and_at_least_50_pairs_and_10_positives": int(((targets_mini["n_with_iptm"] >= 50) & (targets_mini["n_positive_with_iptm"] >= 10)).sum()),
    },
}

# JSON-safe conversion for numpy scalar dictionary values.
def json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if np.isnan(obj) else float(obj)
    if isinstance(obj, float) and np.isnan(obj):
        return None
    raise TypeError(f"Not JSON serializable: {type(obj)}")

pairs.to_csv(OUT / "pair_audit.csv", index=False)
targets_all.to_csv(OUT / "target_summary_all.csv", index=False)
targets_mini.to_csv(OUT / "target_summary_miniprotein.csv", index=False)
metric_audit.to_csv(OUT / "metric_audit.csv", index=False)
pd.DataFrame(score_multiplicity).to_csv(OUT / "score_multiplicity.csv", index=False)
target_method.to_csv(OUT / "target_method_crosstab.csv")
target_author.to_csv(OUT / "target_author_crosstab.csv")
seq_counts.rename("n_design_ids").to_csv(OUT / "sequence_redundancy.csv")
(OUT / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=json_default), encoding="utf-8")

print(json.dumps(summary, ensure_ascii=False, indent=2, default=json_default), flush=True)
print("\nTop targets (all classes):", flush=True)
print(targets_all.head(20).to_string(index=False), flush=True)
print("\nMetric audit:", flush=True)
print(metric_audit.to_string(index=False), flush=True)
print(f"\nOutputs: {OUT}", flush=True)
