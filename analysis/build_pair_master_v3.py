"""Build the locked v3 design-target master table from the raw Proteinbase CSV.

Primary assay:
  - AAH-I: SPR for every design (full-coverage assay; BLI is confirmation only).
  - Other eligible targets: the single observed assay; multi-assay pairs fail closed.

Primary outcome:
  - Majority binding after exact deduplication within design x target x primary assay.

Strict outcome:
  - A strict subset of primary positives.
  - Primary positive + strongest deduplicated primary-assay strength Medium/Strong = 1.
  - Primary positive + strongest strength None/Weak = 0.
  - Primary positive + no strength record = missing (never forced to 0).
  - Primary negative = 0; any Medium/Strong record is separately flagged as a conflict.
"""

from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
INPUT = ROOT / "data" / "proteinbase_latest.csv"
OUT = ROOT / "analysis" / "v3_outputs"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_SHA256 = "8eac08bfc931409e0703eda904cdd1066c9d4a2ceb46af43763f94c74dabb9cc"
MAIN_TARGETS = [
    "nipah-glycoprotein-g",
    "pd-l1",
    "human-vegf110",
    "sars-cov-2-rbd",
    "aah-i",
    "trem2",
]
PRIMARY_ASSAY_OVERRIDES = {"aah-i": "spr"}
STRENGTH_RANK = {"None": 0, "Weak": 1, "Medium": 2, "Strong": 3}
RANK_STRENGTH = {value: key for key, value in STRENGTH_RANK.items()}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_target(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    target = str(value).strip().lower()
    return target or None


def parse_binding(value):
    if isinstance(value, (bool, np.bool_)):
        return int(value)
    if value in (0, 1):
        return int(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "positive"}:
        return 1
    if text in {"false", "0", "no", "negative"}:
        return 0
    raise ValueError(f"Unrecognized binding value: {value!r}")


def parse_strength(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    canonical = {name.lower(): name for name in STRENGTH_RANK}
    if text.lower() not in canonical:
        raise ValueError(f"Unrecognized binding_strength value: {value!r}")
    return canonical[text.lower()]


def majority(series):
    if len(series) == 0:
        raise ValueError("Cannot take a majority of an empty group")
    return int(series.mean() >= 0.5)


def json_dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


input_hash = sha256(INPUT)
if input_hash != EXPECTED_SHA256:
    raise RuntimeError(f"Input hash changed: {input_hash}")

raw = pd.read_csv(
    INPUT,
    usecols=["id", "name", "sequence", "author", "designMethod", "evaluations"],
)
if raw["id"].nunique() != len(raw):
    raise RuntimeError("Design id is not unique in the raw input")

metadata = {}
binding_rows = []
strength_rows = []
score_rows = []
class_values = {}

for row_index, row in enumerate(raw.itertuples(index=False), start=2):
    metadata[row.id] = {
        "design_id": row.id,
        "design_name": row.name,
        "sequence": row.sequence,
        "author": row.author,
        "designMethod": row.designMethod,
        "source_row": row_index,
    }
    evaluations = json.loads(row.evaluations) if pd.notna(row.evaluations) else []
    if not isinstance(evaluations, list):
        raise TypeError(f"evaluations is not a list at source row {row_index}")
    for eval_index, evaluation in enumerate(evaluations):
        metric = evaluation.get("metric")
        if metric == "design_class":
            value = evaluation.get("value")
            if pd.notna(value):
                class_values.setdefault(row.id, set()).add(str(value))
        target = normalize_target(evaluation.get("target"))
        if target not in MAIN_TARGETS:
            continue
        base = {
            "design_id": row.id,
            "target": target,
            "assay": str(evaluation.get("method") or "unknown").strip().lower(),
            "measuredAt": evaluation.get("measuredAt"),
            "source_row": row_index,
            "eval_index": eval_index,
        }
        if evaluation.get("type") == "experimental" and metric == "binding":
            binding_rows.append({**base, "y": parse_binding(evaluation.get("value"))})
        elif evaluation.get("type") == "experimental" and metric == "binding_strength":
            strength = parse_strength(evaluation.get("value"))
            if strength is not None:
                strength_rows.append({**base, "strength": strength})
        elif evaluation.get("type") == "computational" and metric == "boltz2_iptm":
            score = float(evaluation.get("value"))
            if not np.isfinite(score) or not 0 <= score <= 1:
                raise ValueError(f"Invalid Boltz-2 ipTM at source row {row_index}: {score}")
            score_rows.append({"design_id": row.id, "target": target, "boltz2_iptm": score})

binding_raw = pd.DataFrame(binding_rows)
strength_raw = pd.DataFrame(strength_rows)
scores_raw = pd.DataFrame(score_rows)
if binding_raw.empty or scores_raw.empty:
    raise RuntimeError("No main-target binding or ipTM records were parsed")

binding_unique = binding_raw.drop_duplicates(
    ["design_id", "target", "assay", "measuredAt", "y"]
).copy()
strength_unique = strength_raw.drop_duplicates(
    ["design_id", "target", "assay", "measuredAt", "strength"]
).copy()

score_counts = scores_raw.groupby(["design_id", "target"])["boltz2_iptm"].nunique()
if (score_counts != 1).any():
    bad = score_counts.loc[score_counts != 1]
    raise RuntimeError(f"Expected one distinct ipTM per pair; found {len(bad)} exceptions")
scores = scores_raw.groupby(["design_id", "target"], as_index=False).agg(
    boltz2_iptm=("boltz2_iptm", "first"),
    iptm_n_records=("boltz2_iptm", "size"),
)

pair_keys = binding_raw[["design_id", "target"]].drop_duplicates()
pair_keys = pair_keys.merge(scores, on=["design_id", "target"], how="inner", validate="one_to_one")

master_rows = []
assay_audit_rows = []
for key in pair_keys.itertuples(index=False):
    design_id, target = key.design_id, key.target
    raw_pair = binding_raw.loc[(binding_raw["design_id"] == design_id) & (binding_raw["target"] == target)]
    unique_pair = binding_unique.loc[
        (binding_unique["design_id"] == design_id) & (binding_unique["target"] == target)
    ]
    assays = sorted(raw_pair["assay"].unique())
    if target in PRIMARY_ASSAY_OVERRIDES:
        primary_assay = PRIMARY_ASSAY_OVERRIDES[target]
        if primary_assay not in assays:
            raise RuntimeError(f"Missing required {primary_assay} for {design_id} × {target}")
        assay_rule = "full_coverage_spr_for_aah_i"
    else:
        if len(assays) != 1:
            raise RuntimeError(f"Multi-assay non-AAH-I pair requires a rule: {design_id} × {target}: {assays}")
        primary_assay = assays[0]
        assay_rule = "single_observed_assay"

    assay_labels = {}
    for assay in assays:
        raw_assay = raw_pair.loc[raw_pair["assay"] == assay]
        unique_assay = unique_pair.loc[unique_pair["assay"] == assay]
        raw_y = majority(raw_assay["y"])
        dedup_y = majority(unique_assay["y"])
        assay_labels[assay] = {
            "raw_n": int(len(raw_assay)),
            "deduplicated_n": int(len(unique_assay)),
            "raw_positive_n": int(raw_assay["y"].sum()),
            "deduplicated_positive_n": int(unique_assay["y"].sum()),
            "raw_majority": raw_y,
            "deduplicated_majority": dedup_y,
        }
        assay_audit_rows.append({
            "design_id": design_id,
            "target": target,
            "assay": assay,
            "is_primary_assay": int(assay == primary_assay),
            **assay_labels[assay],
        })

    primary_raw = raw_pair.loc[raw_pair["assay"] == primary_assay]
    primary_unique = unique_pair.loc[unique_pair["assay"] == primary_assay]
    y_raw = majority(primary_raw["y"])
    y_primary = majority(primary_unique["y"])
    cross_assay_y = majority(unique_pair["y"])

    strength_pair_raw = strength_raw.loc[
        (strength_raw["design_id"] == design_id)
        & (strength_raw["target"] == target)
        & (strength_raw["assay"] == primary_assay)
    ]
    strength_pair_unique = strength_unique.loc[
        (strength_unique["design_id"] == design_id)
        & (strength_unique["target"] == target)
        & (strength_unique["assay"] == primary_assay)
    ]
    observed_strengths = sorted(
        strength_pair_unique["strength"].unique(), key=lambda item: STRENGTH_RANK[item]
    )
    strongest_rank = max((STRENGTH_RANK[item] for item in observed_strengths), default=None)
    strongest = RANK_STRENGTH.get(strongest_rank) if strongest_rank is not None else pd.NA

    strict_missing_reason = ""
    if y_primary == 0:
        y_strict = 0
    elif strongest_rank is None:
        y_strict = pd.NA
        strict_missing_reason = "primary_positive_missing_primary_assay_strength"
    else:
        y_strict = int(strongest_rank >= STRENGTH_RANK["Medium"])

    negative_with_positive_strength = int(y_primary == 0 and strongest_rank is not None and strongest_rank >= 2)
    strict_downgrade = int(y_primary == 1 and pd.notna(y_strict) and y_strict == 0)
    confirmation = assay_labels.get("bli", {}) if target == "aah-i" and primary_assay == "spr" else {}

    classes = sorted(class_values.get(design_id, set()))
    if len(classes) > 1:
        raise RuntimeError(f"Conflicting design classes for {design_id}: {classes}")
    master_rows.append({
        "pair_id": f"{design_id}|{target}",
        **metadata[design_id],
        "design_class": classes[0] if classes else pd.NA,
        "target": target,
        "primary_assay": primary_assay,
        "primary_assay_rule": assay_rule,
        "assays_observed": ";".join(assays),
        "n_assays_observed": len(assays),
        "assay_label_audit_json": json.dumps(assay_labels, ensure_ascii=False, sort_keys=True),
        "boltz2_iptm": key.boltz2_iptm,
        "iptm_n_records": key.iptm_n_records,
        "primary_binding_raw_n": len(primary_raw),
        "primary_binding_deduplicated_n": len(primary_unique),
        "primary_binding_raw_positive_n": int(primary_raw["y"].sum()),
        "primary_binding_deduplicated_positive_n": int(primary_unique["y"].sum()),
        "primary_binding_raw_majority": y_raw,
        "primary_outcome": y_primary,
        "primary_outcome_all_deduplicated": int(primary_unique["y"].all()),
        "primary_outcome_any_deduplicated": int(primary_unique["y"].any()),
        "primary_outcome_if_concordant": float(y_primary) if primary_unique["y"].nunique() == 1 else np.nan,
        "primary_label_changed_by_exact_dedup": int(y_raw != y_primary),
        "primary_label_discordant_raw": int(primary_raw["y"].nunique() > 1),
        "primary_label_discordant_deduplicated": int(primary_unique["y"].nunique() > 1),
        "primary_outcome_equals_previous_cross_assay_label": int(y_primary == cross_assay_y),
        "primary_strength_raw_n": len(strength_pair_raw),
        "primary_strength_deduplicated_n": len(strength_pair_unique),
        "primary_strength_labels_observed": ";".join(observed_strengths),
        "primary_strength_strongest": strongest,
        "strict_outcome": y_strict,
        "strict_missing_reason": strict_missing_reason,
        "strict_downgrade_from_primary_positive": strict_downgrade,
        "primary_negative_but_medium_or_strong": negative_with_positive_strength,
        "aah_i_confirmation_bli_present": int(bool(confirmation)),
        "aah_i_confirmation_bli_outcome": confirmation.get("deduplicated_majority", pd.NA),
    })

master = pd.DataFrame(master_rows)
master["strict_outcome"] = master["strict_outcome"].astype("Int64")
master["aah_i_confirmation_bli_outcome"] = master["aah_i_confirmation_bli_outcome"].astype("Int64")
master = master.sort_values(["target", "design_id"]).reset_index(drop=True)
assay_audit = pd.DataFrame(assay_audit_rows).sort_values(["target", "design_id", "assay"])

# Reconcile against the earlier pair-level audit without using it to construct labels.
previous = pd.read_csv(
    ROOT / "analysis" / "audit_outputs" / "pair_audit.csv",
    usecols=["design_id", "target", "boltz2_iptm", "y_deduplicated_majority"],
)
previous = previous.loc[previous["target"].isin(MAIN_TARGETS) & previous["boltz2_iptm"].notna()]
reconciled = master.merge(
    previous,
    on=["design_id", "target"],
    how="outer",
    suffixes=("_v3", "_previous"),
    indicator=True,
    validate="one_to_one",
)
if not (reconciled["_merge"] == "both").all():
    raise RuntimeError("The v3 main-pair keys do not match the earlier covered-pair audit")
if not np.allclose(reconciled["boltz2_iptm_v3"], reconciled["boltz2_iptm_previous"], rtol=0, atol=1e-12):
    raise RuntimeError("ipTM values differ from the earlier audit")
previous_label_differences = int(
    (reconciled["primary_outcome"] != reconciled["y_deduplicated_majority"]).sum()
)

main = master
five = master.loc[master["target"] != "aah-i"]
aah_spr = assay_audit.loc[(assay_audit["target"] == "aah-i") & (assay_audit["assay"] == "spr")]
aah_bli = assay_audit.loc[(assay_audit["target"] == "aah-i") & (assay_audit["assay"] == "bli")]
aah_both = aah_bli.merge(
    aah_spr[["design_id", "deduplicated_majority"]],
    on="design_id",
    suffixes=("_bli", "_spr"),
    validate="one_to_one",
)

checks = {
    "main_pairs_equal_2084": len(main) == 2084,
    "unique_designs_equal_2084": main["design_id"].nunique() == 2084,
    "one_row_per_pair": not main.duplicated(["design_id", "target"]).any(),
    "no_design_on_multiple_main_targets": int((main.groupby("design_id")["target"].nunique() > 1).sum()) == 0,
    "primary_outcome_positives_equal_748": int(main["primary_outcome"].sum()) == 748,
    "five_target_pairs_equal_1951": len(five) == 1951,
    "five_target_positives_equal_729": int(five["primary_outcome"].sum()) == 729,
    "aah_i_spr_designs_equal_133": len(aah_spr) == 133,
    "aah_i_spr_positives_equal_19": int(aah_spr["deduplicated_majority"].sum()) == 19,
    "aah_i_bli_designs_equal_18": len(aah_bli) == 18,
    "aah_i_bli_all_positive": int(aah_bli["deduplicated_majority"].sum()) == 18,
    "aah_i_bli_spr_concordant": bool((aah_both["deduplicated_majority_bli"] == aah_both["deduplicated_majority_spr"]).all()),
    "iptm_complete": bool(main["boltz2_iptm"].notna().all()),
    "strict_positive_is_subset_of_primary_positive": bool(
        ((main["strict_outcome"] == 1) & (main["primary_outcome"] != 1)).sum() == 0
    ),
    "strict_missing_only_for_primary_positive": bool(
        (main["strict_outcome"].isna() & (main["primary_outcome"] != 1)).sum() == 0
    ),
}
failed = [name for name, passed in checks.items() if not passed]
if failed:
    raise RuntimeError(f"Critical v3 quality gates failed: {failed}")


def summarize_group(group):
    positive = group.loc[group["primary_outcome"] == 1]
    strength_counts = positive["primary_strength_strongest"].value_counts(dropna=False)
    return {
        "pairs": int(len(group)),
        "unique_designs": int(group["design_id"].nunique()),
        "authors": int(group["author"].nunique(dropna=True)),
        "assays": sorted(group["primary_assay"].unique().tolist()),
        "primary_positive": int(group["primary_outcome"].sum()),
        "primary_positive_rate": float(group["primary_outcome"].mean()),
        "strict_available": int(group["strict_outcome"].notna().sum()),
        "strict_positive": int(group["strict_outcome"].sum(skipna=True)),
        "strict_positive_rate_in_available": float(group["strict_outcome"].mean(skipna=True)),
        "strict_missing": int(group["strict_outcome"].isna().sum()),
        "strict_downgraded_from_primary_positive": int(group["strict_downgrade_from_primary_positive"].sum()),
        "primary_positive_strength": {
            str(key): int(value) for key, value in strength_counts.items()
        },
        "primary_negative_but_medium_or_strong": int(group["primary_negative_but_medium_or_strong"].sum()),
        "iptm_unique_values": int(group["boltz2_iptm"].nunique()),
        "iptm_max_tie_n": int(group["boltz2_iptm"].value_counts().max()),
        "primary_labels_changed_by_exact_dedup": int(group["primary_label_changed_by_exact_dedup"].sum()),
        "multi_assay_pairs": int((group["n_assays_observed"] > 1).sum()),
    }


target_summary = {target: summarize_group(group) for target, group in main.groupby("target", sort=True)}
overall = summarize_group(main)
five_target = summarize_group(five)
summary = {
    "protocol": "v3",
    "input": {
        "path": str(INPUT),
        "bytes": INPUT.stat().st_size,
        "sha256": input_hash,
        "raw_design_rows": int(len(raw)),
    },
    "definitions": {
        "analysis_unit": "design x target with one locked primary assay",
        "primary_assay": "AAH-I uses SPR; other targets use their single observed assay",
        "primary_outcome": "majority binding after exact deduplication within primary assay",
        "strength_aggregation": "strongest observed strength after exact deduplication within primary assay",
        "strict_outcome": "subset of primary positives; Medium/Strong=1, None/Weak=0, missing strength=missing for primary positives; primary negatives=0",
    },
    "overall": overall,
    "five_target_sensitivity": five_target,
    "targets": target_summary,
    "previous_cross_assay_label_differences": previous_label_differences,
    "checks": checks,
}

master_path = OUT / "pair_master_v3.csv"
assay_path = OUT / "assay_unit_audit_v3.csv"
summary_path = OUT / "reconciliation_v3.json"
master.to_csv(master_path, index=False)
assay_audit.to_csv(assay_path, index=False)
json_dump(summary_path, summary)

conflicts = master.loc[master["primary_negative_but_medium_or_strong"] == 1]
if len(conflicts):
    conflicts.to_csv(OUT / "strict_strength_conflicts_v3.csv", index=False)

lines = [
    "# v3 主表口径核对",
    "",
    f"- 输入 SHA-256：`{input_hash}`",
    f"- 主表：`{master_path}`",
    f"- 主表行数／唯一设计数：{len(main)}／{main['design_id'].nunique()}",
    f"- 主结局阳性：{int(main['primary_outcome'].sum())}",
    f"- 严格结局可用／阳性：{int(main['strict_outcome'].notna().sum())}／{int(main['strict_outcome'].sum(skipna=True))}",
    f"- 五靶点敏感性行数／主结局阳性：{len(five)}／{int(five['primary_outcome'].sum())}",
    f"- 与旧跨 assay 聚合标签不同的主标签数：{previous_label_differences}",
    "",
    "| 靶点 | 主 assay | 配对 | 主阳性 | 严格可用 | 严格阳性 | Weak | Medium | Strong | 严格缺失 | ipTM 唯一值／最大并列 |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for target, values in target_summary.items():
    strength = values["primary_positive_strength"]
    lines.append(
        f"| {target} | {';'.join(values['assays'])} | {values['pairs']} | {values['primary_positive']} | "
        f"{values['strict_available']} | {values['strict_positive']} | {strength.get('Weak', 0)} | "
        f"{strength.get('Medium', 0)} | {strength.get('Strong', 0)} | {values['strict_missing']} | "
        f"{values['iptm_unique_values']}／{values['iptm_max_tie_n']} |"
    )
lines.extend([
    "",
    "## 质量门",
    "",
])
for name, passed in checks.items():
    lines.append(f"- {'PASS' if passed else 'FAIL'}：`{name}`")
lines.extend([
    "",
    "严格结局始终是主结局阳性的子集；主结局阴性但出现 Medium／Strong 的记录单独标记，不反向改写成严格阳性。",
    "",
])
(OUT / "reconciliation_v3.md").write_text("\n".join(lines), encoding="utf-8")

print(pd.DataFrame(target_summary).T[["pairs", "primary_positive", "strict_available", "strict_positive"]].to_string())
print("\nOverall:", {key: overall[key] for key in ["pairs", "unique_designs", "primary_positive", "strict_available", "strict_positive", "strict_missing"]})
print("Five-target:", {key: five_target[key] for key in ["pairs", "primary_positive", "strict_available", "strict_positive"]})
print("Previous cross-assay label differences:", previous_label_differences)
print("All critical checks passed:", all(checks.values()))
print("Master table:", master_path)
