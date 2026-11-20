from pathlib import Path
import hashlib
import json
import math
import re

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import range_boundaries

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
MAIN = ROOT / "manuscript/PEDS_Research_Article_English_draft_v2_2026-11-17.md"
GUIDE = ROOT / "manuscript/supplement/Supplementary_Materials_PEDS_2026-11-16.md"
WORKBOOK = ROOT / "manuscript/supplement/Supplementary_Tables_S1_S16.xlsx"
ORIGINAL = ROOT / "manuscript/supplement/Supplementary_Tables_S1_S15.xlsx"
BUILD = ROOT / "manuscript/supplement/supplement_build_manifest.json"
SOURCES = ROOT / "manuscript/supplement/supplement_source_manifest_2026-11-16.md"
RUN = ROOT / "analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855"

checks = []


def check(name, condition, detail=""):
    checks.append((name, bool(condition), str(detail)))
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(value):
    return f"{value:.4f}".replace("-", "−")


def signed_fmt(value):
    return f"{value:+.4f}".replace("-", "−")


def missing(value):
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


main = MAIN.read_text(encoding="utf-8")
guide = GUIDE.read_text(encoding="utf-8")
source_manifest = SOURCES.read_text(encoding="utf-8")
build = json.loads(BUILD.read_text(encoding="utf-8"))

check("main_has_rf_methods_once", main.count("### Random-forest model-complexity sensitivity") == 1)
check("main_has_rf_results_once", main.count("### Random forests alter probability error but do not consistently outperform") == 1)
check("original_pinc_results_retained", "0.2949 versus 0.3139" in main)
check("main_plain_language", "zero-label" not in main)
check("guide_plain_language", "zero-label" not in guide)
check("main_links_s16_workbook", "supplement/Supplementary_Tables_S1_S16.xlsx" in main)
check("guide_has_s16", guide.count("### Table S16.") == 1)
check("main_has_correct_table_count", "75 Excel tables; S16 contains six logical panels" in main)
check("no_prediction_panel_overclaim", "candidate-level tuning, predictions, and provenance" not in main)

# Validate that all relative Markdown links resolve.
for document in [MAIN, GUIDE]:
    text = document.read_text(encoding="utf-8")
    targets = re.findall(r"\]\(([^)]+)\)", text)
    unresolved = []
    for target in targets:
        if target.startswith(("http://", "https://", "#", "mailto:")):
            continue
        clean = target.split("#", 1)[0]
        if not (document.parent / clean).exists():
            unresolved.append(target)
    check(f"relative_links_resolve_{document.name}", not unresolved, unresolved)

six = pd.read_csv(RUN / "rq5rf_six_configuration_summary.csv")
macro = pd.read_csv(RUN / "rq5rf_macro_metrics.csv")
target = pd.read_csv(RUN / "rq5rf_target_metrics.csv")
main_six = six[(six.analysis_set == "primary_six") & (six.source_weighting == "target_equal")].iloc[0]
main_u = macro[
    (macro.analysis_set == "primary_six")
    & (macro.source_weighting == "target_equal")
    & (macro.model_layer_label == "U_transport_zero_label")
].set_index("model_exposure")

required_values = []
for model in ["pooled_logistic_ipTM", "rf_ipTM", "rf_combo"]:
    for metric in ["brier", "log_loss", "bss_transport", "roc_auc", "ap"]:
        if model == "rf_combo" and metric in ["roc_auc", "ap"]:
            continue
        required_values.append(fmt(main_u.loc[model, metric]))
required_values.append(fmt(main_u.loc["rf_ipTM", "bss_target_oracle"]))
check("main_rf_absolute_values_match_source", all(v in main for v in required_values), [v for v in required_values if v not in main])

five = six[(six.analysis_set == "primary_five") & (six.source_weighting == "target_equal")].iloc[0]
five_deltas = [
    five["rf_combo_minus_pooled_logistic_combo_brier_difference"],
    five["rf_combo_minus_pooled_logistic_combo_log_loss_difference"],
    five["rf_combo_minus_pooled_logistic_combo_bss_transport_difference"],
]
check("five_target_reversal_values_match", all(signed_fmt(v) in main for v in five_deltas), [signed_fmt(v) for v in five_deltas])
check("five_target_reversal_direction", five_deltas[0] > 0 and five_deltas[1] > 0 and five_deltas[2] < 0, five_deltas)

check(
    "pinc_within_rf_worse_all_six",
    (six.rf_combo_minus_rf_ipTM_brier_difference > 0).all()
    and (six.rf_combo_minus_rf_ipTM_log_loss_difference > 0).all()
    and (six.rf_combo_minus_rf_ipTM_bss_transport_difference < 0).all(),
)

u_target = target[
    (target.analysis_set == "primary_six")
    & (target.source_weighting == "target_equal")
    & (target.model_layer_label == "U_transport_zero_label")
]
pivot = u_target.pivot(index="target", columns="model_exposure", values="brier")
check("rf_iptm_brier_lower_in_five_of_six", int((pivot.rf_ipTM < pivot.pooled_logistic_ipTM).sum()) == 5)
check("rbd_difference_negligible", abs(pivot.loc["sars-cov-2-rbd", "rf_ipTM"] - pivot.loc["sars-cov-2-rbd", "pooled_logistic_ipTM"]) < 0.0001)

f_combo = macro[
    (macro.analysis_set == "primary_six")
    & (macro.source_weighting == "target_equal")
    & (macro.model_layer_label == "predictive_F")
    & (macro.model_exposure == "rf_combo")
].iloc[0]
check("local_adaptation_values_match", fmt(f_combo.brier) in main and fmt(f_combo.bss_target_oracle) in main)

# Workbook structure and exact transcription checks for native CSV panels.
wb = load_workbook(WORKBOOK, read_only=False, data_only=False)
check("workbook_sheet_set", wb.sheetnames == build["expected_sheets"])
check("workbook_has_s16", "S16" in wb.sheetnames)
check("build_table_count", len(build["tables"]) == 75, len(build["tables"]))
check("s16_table_count", len(wb["S16"].tables) == 11)
check("cover_s5_title_restored", wb["Cover"]["B25"].value == "Main-cohort LOTO U/I/F and sensitivity analyses")
check("cover_s10_title_restored", wb["Cover"]["B30"].value == "Post hoc exploratory external I/F adaptation")
legacy_cells = [
    (ws.title, cell.coordinate)
    for ws in wb.worksheets
    for row in ws.iter_rows()
    for cell in row
    if isinstance(cell.value, str) and "zero-label" in cell.value
]
check("workbook_plain_language", not legacy_cells, legacy_cells)
check("workbook_has_no_formulas", sum(1 for ws in wb.worksheets for row in ws.iter_rows() for cell in row if cell.data_type == "f") == 0)

panel_data = {
    "B.": macro[macro.model_layer_label == "U_transport_zero_label"].copy(),
    "C1.": target[target.model_layer_label == "U_transport_zero_label"].copy(),
    "C2.": pd.read_csv(RUN / "rq5rf_pairwise_target_differences.csv").query("model_layer_label == 'U_transport_zero_label'"),
    "D1.": macro[macro.model_layer_label.isin(["predictive_I", "predictive_F"])].copy(),
    "D2.": target[target.model_layer_label.isin(["predictive_I", "predictive_F"])].copy(),
    "E1.": pd.read_csv(RUN / "rq5rf_reliability.csv"),
    "E2.": pd.read_csv(RUN / "failed_item_ledger.csv"),
    "F1.": pd.read_csv(RUN / "rf_inner_candidate_selection.csv"),
    "F2.": pd.read_csv(RUN / "rf_final_candidate_fit_audit.csv"),
}
for entry in build["tables"]:
    if entry.get("sheet") != "S16":
        continue
    prefix = next((p for p in panel_data if entry["panel"].startswith(p)), None)
    if prefix is None:
        continue
    min_col, min_row, max_col, max_row = range_boundaries(entry["range"])
    values = list(wb["S16"].iter_rows(min_row=min_row, max_row=max_row, min_col=min_col, max_col=max_col, values_only=True))
    actual = pd.DataFrame(values[1:], columns=values[0])
    expected = panel_data[prefix].reset_index(drop=True)
    expected.columns = expected.columns.astype(str)
    check(f"panel_shape_{entry['panel']}", actual.shape == expected.shape, (actual.shape, expected.shape))
    check(f"panel_columns_{entry['panel']}", list(actual.columns) == list(expected.columns))
    mismatch = None
    for col in expected.columns:
        for idx, exp in enumerate(expected[col]):
            act = actual.iloc[idx][col]
            if missing(exp) and missing(act):
                continue
            if isinstance(exp, (int, float, np.integer, np.floating)) and not isinstance(exp, (bool, np.bool_)):
                try:
                    equal = math.isclose(float(act), float(exp), rel_tol=1e-12, abs_tol=1e-12)
                except (TypeError, ValueError):
                    equal = False
            else:
                equal = act == exp
            if not equal:
                mismatch = (col, idx, act, exp)
                break
        if mismatch:
            break
    check(f"panel_values_{entry['panel']}", mismatch is None, mismatch)

# Generated-file and source hashes must describe the actual delivered files.
check("workbook_hash_in_manifest", f"| `{WORKBOOK.name}` | {WORKBOOK.stat().st_size} | `{sha256(WORKBOOK)}` |" in source_manifest)
check("guide_hash_in_manifest", f"| `{GUIDE.name}` | {GUIDE.stat().st_size} | `{sha256(GUIDE)}` |" in source_manifest)
for path in RUN.iterdir():
    if path.name in {
        "analysis_set_audit.csv", "rq5rf_macro_metrics.csv", "rq5rf_target_metrics.csv",
        "rq5rf_pairwise_target_differences.csv", "rq5rf_reliability.csv", "failed_item_ledger.csv",
        "rf_inner_candidate_selection.csv", "rf_final_candidate_fit_audit.csv", "qc_checks.json",
        "RUN_STATUS.json", "independent_verification.json", "direction_consistency_10_2.json",
    }:
        relative = path.relative_to(ROOT).as_posix()
        check(f"source_hash_{path.name}", f"| `{relative}` | {path.stat().st_size} | `{sha256(path)}` |" in source_manifest)
check("original_workbook_preserved", sha256(ORIGINAL) == "3b2cc26f664e06abb1be51bbc8d22a78487d2dbef23cd42d0df6bfb339318bd4")

print(f"All {len(checks)} RQ5-RF integration checks passed.")
print(f"S16: 6 logical panels, 11 Excel tables, {wb['S16'].max_row} rows; total workbook tables: {len(build['tables'])}")
print(f"Workbook SHA-256: {sha256(WORKBOOK)}")
print("Native CSV panels were compared element by element against the frozen outputs.")
