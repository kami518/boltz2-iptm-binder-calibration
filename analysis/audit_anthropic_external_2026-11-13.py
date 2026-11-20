"""Audit external Anthropic design cohort without modifying the frozen RQ1–RQ3 tables."""
from pathlib import Path
import hashlib
import json
import duckdb
import numpy as np
import pandas as pd

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
DATA = ROOT / "data" / "external_target_audit_2026-11-13"
DESIGN = DATA / "anthropic_design_summary.parquet"
WETLAB = DATA / "anthropic_wetlab_summary.parquet"
COFOLD = DATA / "anthropic_cofold_predictions.parquet"
PB = ROOT / "data" / "proteinbase_latest.csv"
OVERATH = ROOT / "data" / "overath_hf.parquet"
con = duckdb.connect()

def readpq(path):
    return con.execute(f"select * from '{path.as_posix()}'").df()

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def normseq(s):
    return s.astype("string").str.strip().str.upper()

d = readpq(DESIGN)
w = readpq(WETLAB)
assert len(d) == len(w) == 1440
assert d.uuid.is_unique and w.uuid.is_unique
keep_w = ["uuid", "assessment_binder", "assessment_binder_confidence", "adaptyv_tested",
          "adaptyv_binder", "adaptyv_binding_strength", "twist_measured", "twist_hit_any_human",
          "twist_hit_any_human_qc_pass", "tested_at_both_vendors", "vendor_agreement_class"]
x = d.merge(w[keep_w], on="uuid", how="left", validate="one_to_one", suffixes=("", "_wet"))
x["sequence"] = normseq(x.sequence)
assert x.sequence.notna().all()
assert not x.sequence.str.contains("[^ACDEFGHIKLMNPQRSTVWY]", regex=True).any()
# Boolean comparisons deliberately retain missing values rather than treating them as negative.
x["integrated_valid"] = x.binder_final.notna()
x["integrated_positive"] = x.binder_final == True
x["integrated_high_medium_positive"] = x.integrated_positive & x.assessment_binder_confidence.isin(["high", "medium"])
x["adaptyv_valid"] = x.adaptyv_tested == True
x["adaptyv_main_positive"] = x.adaptyv_valid & x.adaptyv_binding_strength.str.lower().isin(["weak", "medium", "strong"])
x["adaptyv_strict_positive"] = x.adaptyv_valid & x.adaptyv_binding_strength.str.lower().isin(["medium", "strong"])
x["twist_valid"] = x.twist_measured == True
x["twist_permissive_positive"] = x.twist_valid & (x.twist_hit_any_human == True)
x["twist_qc_positive"] = x.twist_valid & (x.twist_hit_any_human_qc_pass == True)

c = readpq(COFOLD)
model_coverage = (c.groupby("cofolding_model").agg(score_rows=("uuid", "size"), designs=("uuid", "nunique"),
                  iptm_nonmissing=("iptm_pae", lambda s: int(s.notna().sum()))).reset_index())
model_coverage["iptm_design_coverage"] = model_coverage.cofolding_model.map(
    c[c.iptm_pae.notna()].groupby("cofolding_model").uuid.nunique()).fillna(0).astype(int)
model_coverage.to_csv(DATA / "anthropic_model_coverage.csv", index=False)

# Design-target, design-time stoichiometries only; no alternate-antigen or N:N counterfolds.
b = c[(c.cofolding_model == "boltz2") & (c.target_form.isna() | (c.target_form == "")) &
      c.stoichiometry.isin(["1to1", "1to2", "1to3", "rnp"])].copy()
assert not b.duplicated(["uuid", "seed"]).any()
# A single frozen predictor summary: ipTM at the seed selected by max ipSAE, ties by seed.
b["seed_text"] = b.seed.astype(str)
b = b.sort_values(["uuid", "ipsae_min", "seed_text"], ascending=[True, False, True])
best = b.drop_duplicates("uuid", keep="first")[["uuid", "iptm_pae", "ipsae_min", "seed", "stoichiometry"]].rename(
    columns={"iptm_pae":"boltz2_iptm_ipsae_best_seed", "ipsae_min":"boltz2_ipsae_selected_seed",
             "seed":"boltz2_selected_seed", "stoichiometry":"boltz2_stoichiometry"})
agg = b.groupby("uuid").iptm_pae.agg(boltz2_iptm_seed_mean="mean", boltz2_iptm_seed_median="median",
                                    boltz2_iptm_seed_max="max", boltz2_n_seeds="size").reset_index()
x = x.merge(best, on="uuid", how="left", validate="one_to_one").merge(agg, on="uuid", how="left", validate="one_to_one")
x["boltz2_iptm_available"] = x.boltz2_iptm_ipsae_best_seed.notna()

# Exact overlaps against both published reference sets.
pb = pd.read_csv(PB, usecols=["id", "name", "sequence"])
pb["sequence"] = normseq(pb.sequence)
overath = readpq(OVERATH)
# Verified against the full local CSV: HF SeqA is the unique binder; SeqB is target (19 unique).
overath["sequence"] = normseq(overath.SeqA)
assert overath.sequence.nunique() == len(overath) == 3669
pb_seqs = set(pb.sequence.dropna())
ov_seqs = set(overath.sequence.dropna())
x["exact_overlap_proteinbase"] = x.sequence.isin(pb_seqs)
x["exact_overlap_overath_hf"] = x.sequence.isin(ov_seqs)
pb_matches = x[x.exact_overlap_proteinbase].merge(pb, on="sequence", how="left", suffixes=("", "_proteinbase"))
pb_matches.to_csv(DATA / "anthropic_proteinbase_exact_overlap.csv", index=False)
ov_matches = x[x.exact_overlap_overath_hf][["uuid", "full_name", "target", "sequence", "binder_final"]]
ov_matches.to_csv(DATA / "anthropic_overath_exact_overlap.csv", index=False)

schemes = {
    "integrated": ("integrated_valid", "integrated_positive"),
    "integrated_high_medium": ("integrated_valid", "integrated_high_medium_positive"),
    "adaptyv_main_weak_plus": ("adaptyv_valid", "adaptyv_main_positive"),
    "adaptyv_strict_medium_plus": ("adaptyv_valid", "adaptyv_strict_positive"),
    "twist_permissive": ("twist_valid", "twist_permissive_positive"),
    "twist_qc": ("twist_valid", "twist_qc_positive"),
}
rows = []
for target, g in x.groupby("target", sort=True):
    row = {"target": target, "all_designs": len(g), "unique_sequences": g.sequence.nunique(),
           "boltz2_iptm_n": int(g.boltz2_iptm_available.sum()),
           "design_models": ";".join(sorted(g.design_model.unique())),
           "campaigns": ";".join(sorted(g.campaign.unique()))}
    for label, (valid, positive) in schemes.items():
        use = g[g[valid] & g.boltz2_iptm_available]
        row[f"{label}_n"] = len(use)
        row[f"{label}_positives"] = int(use[positive].sum())
        row[f"{label}_prevalence"] = float(use[positive].mean()) if len(use) else np.nan
    row["integrated_qualifies"] = bool(row["integrated_n"] >= 80 and row["integrated_positives"] >= 15)
    row["adaptyv_main_qualifies"] = bool(row["adaptyv_main_weak_plus_n"] >= 80 and row["adaptyv_main_weak_plus_positives"] >= 15)
    row["overlap_proteinbase_n"] = int(g.exact_overlap_proteinbase.sum())
    row["overlap_overath_n"] = int(g.exact_overlap_overath_hf.sum())
    rows.append(row)
target_counts = pd.DataFrame(rows)
current_six = {"Nipah-G", "PD-L1", "VEGF-A", "TREM2"}  # two current names are mapped below.
target_counts["in_current_six"] = target_counts.target.isin(current_six)
target_counts["new_target_vs_current_six"] = ~target_counts.in_current_six
target_counts = target_counts.sort_values(["integrated_qualifies", "integrated_n"], ascending=[False, False])
target_counts.to_csv(DATA / "anthropic_target_qualification.csv", index=False)

campaign_rows = []
for (target, design_model, campaign), g in x.groupby(["target", "design_model", "campaign"], sort=True):
    use = g[g.integrated_valid & g.boltz2_iptm_available]
    a = g[g.adaptyv_valid & g.boltz2_iptm_available]
    t = g[g.twist_valid & g.boltz2_iptm_available]
    campaign_rows.append({"target":target, "design_model":design_model, "campaign":campaign,
        "integrated_n":len(use), "integrated_positives":int(use.integrated_positive.sum()),
        "adaptyv_n":len(a), "adaptyv_positives":int(a.adaptyv_main_positive.sum()),
        "twist_n":len(t), "twist_qc_positives":int(t.twist_qc_positive.sum())})
campaign_counts = pd.DataFrame(campaign_rows)
campaign_counts.to_csv(DATA / "anthropic_campaign_counts.csv", index=False)
agreement = (x[x.tested_at_both_vendors == True].groupby(["target", "vendor_agreement_class"]).size()
             .unstack(fill_value=0).reset_index())
agreement.to_csv(DATA / "anthropic_vendor_agreement_by_target.csv", index=False)

x.to_csv(DATA / "anthropic_audited_pair_table.csv", index=False)
qual = target_counts[target_counts.integrated_qualifies]
summary = {
    "source_files": {p.name: {"bytes":p.stat().st_size, "sha256":sha(p)} for p in [DESIGN, WETLAB, COFOLD]},
    "n_designs":len(x), "n_targets":int(x.target.nunique()), "n_unique_sequences":int(x.sequence.nunique()),
    "integrated_valid_n":int(x.integrated_valid.sum()), "integrated_positive_n":int(x.integrated_positive.sum()),
    "boltz2_iptm_n":int(x.boltz2_iptm_available.sum()),
    "boltz2_design_target_rows":len(b),
    "primary_iptm_rule":"Design-target, design stoichiometry (1to1/1to2/1to3/rnp); ipTM from highest-ipSAE seed, ties by seed label.",
    "n_targets_integrated_80_15":len(qual),
    "qualified_targets_integrated":qual.target.tolist(),
    "qualified_new_targets_vs_current_six":qual[qual.new_target_vs_current_six].target.tolist(),
    "known_target_mapping_note":"VEGF-A maps to VEGF110 and Nipah-G to Nipah; PD-L1 and TREM2 are shared. Other Anthropic targets, including EGFR/IL-7Ra, are not among the frozen six-target RQ3 cohort.",
    "exact_overlap_proteinbase_designs":int(x.exact_overlap_proteinbase.sum()),
    "exact_overlap_proteinbase_sequences":int(x.loc[x.exact_overlap_proteinbase, "sequence"].nunique()),
    "exact_overlap_overath_designs":int(x.exact_overlap_overath_hf.sum()),
    "exact_overlap_overath_sequences":int(x.loc[x.exact_overlap_overath_hf, "sequence"].nunique()),
    "cautions":["Exact sequence identity is not a near-duplicate/cluster audit.",
                "Integrated labels combine vendors and an author/Claude rubric; do not treat as a third independent assay.",
                "Do not equate five-seed best-ipSAE ipTM with frozen Proteinbase ipTM before a configuration comparability check.",
                "Anthropic data have appeared in a 2026 calibration dataset; external to the frozen cohort does not mean unused by prior calibration literature.",
                "Campaign counts are descriptive; combining campaigns requires an a-priori target-level outcome rule."],
}
(DATA / "anthropic_audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
print("\nTARGET COUNTS\n", target_counts[["target","integrated_n","integrated_positives","adaptyv_main_weak_plus_n",
    "adaptyv_main_weak_plus_positives","twist_qc_n","twist_qc_positives","integrated_qualifies",
    "overlap_proteinbase_n","overlap_overath_n"]].to_string(index=False))
print("\nMODEL COVERAGE\n", model_coverage.to_string(index=False))
