from pathlib import Path
import hashlib
import json
import duckdb
import pandas as pd

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
OUT = ROOT / "data" / "external_target_audit_2026-11-13"
csv_path = ROOT / "data" / "proteinbase_latest.csv"
pq_path = OUT / "proteinbase_current.parquet"

frozen = pd.read_csv(csv_path)
current = duckdb.sql(f"SELECT * FROM read_parquet('{pq_path.as_posix()}')").df()

result = {
    "frozen_csv": {
        "path": str(csv_path),
        "sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "shape": list(frozen.shape),
        "columns": list(frozen.columns),
    },
    "current_parquet": {
        "path": str(pq_path),
        "sha256": hashlib.sha256(pq_path.read_bytes()).hexdigest(),
        "shape": list(current.shape),
        "columns": list(current.columns),
    },
}

for key in ["id", "name", "sequence"]:
    if key in frozen.columns and key in current.columns:
        f = frozen[key].fillna("").astype(str)
        c = current[key].fillna("").astype(str)
        result[f"{key}_comparison"] = {
            "frozen_nunique": int(f.nunique()),
            "current_nunique": int(c.nunique()),
            "frozen_only": sorted(set(f) - set(c))[:20],
            "current_only": sorted(set(c) - set(f))[:20],
            "frozen_only_count": len(set(f) - set(c)),
            "current_only_count": len(set(c) - set(f)),
        }

if "id" in frozen.columns and "id" in current.columns:
    joined = frozen[["id", "sequence", "evaluations"]].merge(
        current[["id", "sequence", "evaluations"]], on="id", how="outer", suffixes=("_frozen", "_current"), indicator=True
    )
    result["join_status"] = joined["_merge"].value_counts().to_dict()
    both = joined[joined["_merge"] == "both"].copy()
    result["sequence_mismatch_same_id"] = int((both["sequence_frozen"].fillna("") != both["sequence_current"].fillna("")).sum())

    def normalize_value(value):
        if value is None:
            return "None"
        if isinstance(value, bool):
            return "true" if value else "false"
        return str(value)

    def canonical_evaluations(value, frozen_format):
        records = json.loads(value) if frozen_format else list(value)
        normalized = []
        fields = ["type", "metric", "valueType", "unit", "value", "target", "method", "measuredAt"]
        for record in records:
            record = dict(record)
            normalized_record = {key: normalize_value(record.get(key)) for key in fields}
            if record.get("valueType") == "json":
                try:
                    raw_json = record.get("value")
                    parsed_json = json.loads(raw_json) if isinstance(raw_json, str) else raw_json
                    normalized_record["value"] = json.dumps(
                        parsed_json, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                    )
                except (TypeError, json.JSONDecodeError):
                    pass
            normalized.append(normalized_record)
        return sorted(json.dumps(record, sort_keys=True, ensure_ascii=False) for record in normalized)

    evaluation_mismatches = []
    evaluation_count_mismatches = []
    first_evaluation_difference = None
    for row in both.itertuples(index=False):
        frozen_evals = canonical_evaluations(row.evaluations_frozen, True)
        current_evals = canonical_evaluations(row.evaluations_current, False)
        if len(frozen_evals) != len(current_evals):
            evaluation_count_mismatches.append(row.id)
        if frozen_evals != current_evals:
            evaluation_mismatches.append(row.id)
            if first_evaluation_difference is None:
                first_difference = next(
                    ({"index": i, "frozen": f, "current": c} for i, (f, c) in enumerate(zip(frozen_evals, current_evals)) if f != c),
                    {"note": "only list length differed"},
                )
                first_evaluation_difference = {"id": row.id, "difference": first_difference}
    result["evaluation_count_mismatches"] = evaluation_count_mismatches
    result["evaluation_content_mismatch_count"] = len(evaluation_mismatches)
    result["evaluation_content_mismatch_ids_first20"] = evaluation_mismatches[:20]
    result["first_evaluation_difference"] = first_evaluation_difference

(OUT / "proteinbase_snapshot_comparison.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
