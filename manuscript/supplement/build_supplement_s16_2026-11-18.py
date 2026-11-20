from pathlib import Path
from copy import copy
import hashlib
import json
import math

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
SUPP_DIR = ROOT / "manuscript/supplement"
SOURCE_WB = SUPP_DIR / "Supplementary_Tables_S1_S15.xlsx"
OUTPUT_WB = SUPP_DIR / "Supplementary_Tables_S1_S16.xlsx"
RUN = ROOT / "analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855"
PROTOCOL = ROOT / "analysis/RQ5RF_model_complexity_frozen_protocol_2026-11-17.md"
BUILD_MANIFEST = SUPP_DIR / "supplement_build_manifest.json"
SOURCE_MANIFEST = SUPP_DIR / "supplement_source_manifest_2026-11-16.md"
GUIDE = SUPP_DIR / "Supplementary_Materials_PEDS_2026-11-16.md"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def read_csv(name):
    return pd.read_csv(RUN / name)


def clean_value(value):
    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None if math.isnan(value) else str(value)
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            value = value.item()
        except Exception:
            pass
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, str) and value.startswith("="):
        return "'" + value
    return value


original_hash = sha256(SOURCE_WB)
wb = load_workbook(SOURCE_WB)
if "S16" in wb.sheetnames:
    raise RuntimeError("S16 already exists; refusing to overwrite an integrated workbook.")

# Update the cover while preserving its existing layout.
cover = wb["Cover"]
cover["B4"] = "Supplementary Tables S1–S16: context-dependent ranking, calibration, transport, and experimental-selection value of Boltz-2 ipTM"
cover["B5"] = "2026-11-18"
cover["B7"] = "Cover, Sources, and S1–S16. Multiple panels are stacked vertically and named as Excel tables."
cover["B10"] = "Main-cohort RQ1/RQ2/RQ3, external confirmatory U without destination-target labels, and prespecified sensitivity analyses."
cover["B11"] = "External post hoc I/F, retrospective shared-budget simulation, Pinc increments, random-forest model-complexity sensitivity, and paired predictor-pipeline descriptions."
cover["B25"] = "Main-cohort LOTO U/I/F and sensitivity analyses"
cover["C25"] = "Test prediction without destination-target labels separately from label-dependent local adaptation and identify what local measurements correct."
cover["D25"] = "U uses no destination-target labels. Predictive I/F use local labels and out-of-fold evaluation; descriptive/oracle fits are not deployable evidence."
cover["B30"] = "Post hoc exploratory external I/F adaptation"
cover["C30"] = "Assess what representative local labels can correct in the external context without changing the confirmatory result obtained without destination-target labels."
cover["D30"] = "Added after external U was examined. Five-fold outer predictions and nested tuning use local labels; these results are exploratory and do not validate prediction without destination-target labels."
wb["S5"]["A2"] = "Purpose: Test prediction without destination-target labels separately from label-dependent local adaptation and identify what local measurements correct."
wb["S10"]["A2"] = "Purpose: Assess what representative local labels can correct in the external context without changing the confirmatory result obtained without destination-target labels."
wb["S10"]["A3"] = "Evidence boundary: Added after external U was examined. Five-fold outer predictions and nested tuning use local labels; these results are exploratory and do not validate prediction without destination-target labels."
wb["S15"]["C441"] = "constant p = actual held-out-target prevalence; diagnostic only, never deployable without destination-target labels"
wb["S15"]["C460"] = "for each LOTO target, use all target labels only after U is fixed; solve a so mean(expit(a+eta_U)) equals full-target prevalence; label oracle_base_rate_scenario and never describe it as prediction without destination-target labels"
for col in range(1, 6):
    src = cover.cell(35, col)
    dst = cover.cell(36, col)
    if src.has_style:
        dst._style = copy(src._style)
    dst.alignment = copy(src.alignment)
cover.row_dimensions[36].height = cover.row_dimensions[35].height
cover.cell(36, 1).value = "S16"
cover.cell(36, 2).value = "Random-forest model-complexity sensitivity on common complete cases"
cover.cell(36, 3).value = "Test whether a more flexible functional form changes held-out-target probability error or the incremental-Pinc conclusion under identical rows, folds, weights, and adaptation rules."
cover.cell(36, 4).value = "Post hoc and exploratory; one random-forest grid only. Earlier target-intercept models are descriptive anchors, not pure model-family comparators. No independent external Pinc cohort or pure-interaction identification."
cover.cell(36, 5).value = 6

# Assemble the six logical panels. Subtables retain their native source schemas.
analysis = read_csv("analysis_set_audit.csv")
model_rows = []
for exposure, family, features in [
    ("pooled_logistic_ipTM", "pooled ridge logistic", "logit(ipTM)"),
    ("pooled_logistic_combo", "pooled ridge logistic", "logit(ipTM), logit(Pinc)"),
    ("rf_ipTM", "random forest", "logit(ipTM)"),
    ("rf_combo", "random forest", "logit(ipTM), logit(Pinc)"),
]:
    model_rows.append({
        "record_type": "model_definition",
        "model_exposure": exposure,
        "model_family": family,
        "features": features,
        "source_target_intercepts": "none; one global intercept for pooled logistic",
        "status": "frozen",
    })
analysis.insert(0, "record_type", "analysis_set_audit")
analysis["model_exposure"] = None
analysis["model_family"] = None
analysis["features"] = None
analysis["source_target_intercepts"] = None
panel_a = pd.concat([analysis, pd.DataFrame(model_rows)], ignore_index=True)

macro = read_csv("rq5rf_macro_metrics.csv")
target = read_csv("rq5rf_target_metrics.csv")
pairwise_target = read_csv("rq5rf_pairwise_target_differences.csv")
panel_b = macro[macro.model_layer_label == "U_transport_zero_label"].copy()
panel_c1 = target[target.model_layer_label == "U_transport_zero_label"].copy()
panel_c2 = pairwise_target[pairwise_target.model_layer_label == "U_transport_zero_label"].copy()
panel_d1 = macro[macro.model_layer_label.isin(["predictive_I", "predictive_F"])].copy()
panel_d2 = target[target.model_layer_label.isin(["predictive_I", "predictive_F"])].copy()
panel_e1 = read_csv("rq5rf_reliability.csv")
panel_e2 = read_csv("failed_item_ledger.csv")
panel_f1 = read_csv("rf_inner_candidate_selection.csv")
panel_f2 = read_csv("rf_final_candidate_fit_audit.csv")

source_names = [
    "analysis_set_audit.csv",
    "rq5rf_macro_metrics.csv",
    "rq5rf_target_metrics.csv",
    "rq5rf_pairwise_target_differences.csv",
    "rq5rf_reliability.csv",
    "failed_item_ledger.csv",
    "rf_inner_candidate_selection.csv",
    "rf_final_candidate_fit_audit.csv",
    "rq5rf_six_configuration_summary.csv",
    "rq5rf_pairwise_macro_comparisons.csv",
    "final_hash_manifest.csv",
    "run_config.snapshot.json",
    "environment.json",
    "qc_checks.json",
    "RUN_STATUS.json",
    "independent_verification.json",
    "direction_consistency_10_2.json",
    "rq5rf_unit_tests.json",
]
source_paths = [RUN / name for name in source_names] + [PROTOCOL]
for path in source_paths:
    if not path.is_file():
        raise FileNotFoundError(path)

qc = json.loads((RUN / "qc_checks.json").read_text(encoding="utf-8"))
run_status = json.loads((RUN / "RUN_STATUS.json").read_text(encoding="utf-8"))
direction = json.loads((RUN / "direction_consistency_10_2.json").read_text(encoding="utf-8"))
independent = json.loads((RUN / "independent_verification.json").read_text(encoding="utf-8"))
qc_rows = []
for item in qc["checks"]:
    qc_rows.append({
        "record_type": "qc_check",
        "check": item.get("check"),
        "passed": item.get("passed"),
        "details_json": json.dumps(item, ensure_ascii=False, sort_keys=True),
    })
for name, obj in [
    ("run_status", run_status),
    ("direction_consistency_10_2", direction),
    ("independent_verification", independent),
]:
    qc_rows.append({
        "record_type": "summary_record",
        "check": name,
        "passed": obj.get("all_passed", obj.get("all_checks_passed")) if isinstance(obj, dict) else None,
        "details_json": json.dumps(obj, ensure_ascii=False, sort_keys=True),
    })
for path in source_paths:
    qc_rows.append({
        "record_type": "source_provenance",
        "relative_path": rel(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    })
panel_f3 = pd.DataFrame(qc_rows)

# Add new sources to the existing Sources sheet.
sources_ws = wb["Sources"]
existing_sources = {sources_ws.cell(r, 1).value for r in range(7, sources_ws.max_row + 1)}
style_row = sources_ws.max_row
for path in source_paths:
    if rel(path) in existing_sources:
        continue
    r = sources_ws.max_row + 1
    for c in range(1, 4):
        if sources_ws.cell(style_row, c).has_style:
            sources_ws.cell(r, c)._style = copy(sources_ws.cell(style_row, c)._style)
    sources_ws.cell(r, 1).value = rel(path)
    sources_ws.cell(r, 2).value = path.stat().st_size
    sources_ws.cell(r, 3).value = sha256(path)

s16 = wb.create_sheet("S16")
s16.sheet_view.showGridLines = False
s16.freeze_panes = "A4"
for col in range(1, 4):
    s16.cell(1, col).value = None
s16["A1"] = "S16. Random-forest model-complexity sensitivity on common complete cases"
s16["A1"]._style = copy(wb["S13"]["A1"]._style)
s16["A2"] = "Purpose: Determine whether a more flexible functional form changes held-out-target probability error or the incremental-Pinc conclusion when samples, folds, feature sets, source weights, and local-adaptation rules are held fixed."
s16["A2"]._style = copy(wb["S13"]["A2"]._style)
s16["A3"] = "Evidence boundary: Post hoc exploratory analysis of the evaluated random-forest grid. Earlier target-intercept models are descriptive anchors, not pure model-family comparators. No independent external Pinc cohort; combined-feature contrasts do not identify pure interactions."
s16["A3"]._style = copy(wb["S13"]["A3"]._style)
for row in (1, 2, 3):
    s16.merge_cells(start_row=row, start_column=1, end_row=row, end_column=10)
    s16.row_dimensions[row].height = 30 if row == 1 else 42
    s16.cell(row, 1).alignment = Alignment(wrap_text=True, vertical="center")

panel_counter = 0
manifest_entries = []
current_row = 5


def add_panel(label, title, note, tables):
    global panel_counter, current_row
    s16.cell(current_row, 1).value = f"{label}. {title}"
    s16.cell(current_row, 1)._style = copy(wb["S13"]["A5"]._style)
    s16.merge_cells(start_row=current_row, start_column=1, end_row=current_row, end_column=10)
    s16.row_dimensions[current_row].height = 24
    current_row += 1
    s16.cell(current_row, 1).value = note
    s16.cell(current_row, 1)._style = copy(wb["S13"]["A6"]._style)
    s16.cell(current_row, 1).alignment = Alignment(wrap_text=True, vertical="center")
    s16.merge_cells(start_row=current_row, start_column=1, end_row=current_row, end_column=10)
    s16.row_dimensions[current_row].height = 34
    current_row += 2
    for sublabel, subtitle, df, source_list in tables:
        if df.empty:
            raise RuntimeError(f"Empty S16 table: {label}{sublabel}")
        s16.cell(current_row, 1).value = f"{label}{sublabel}. {subtitle}"
        s16.cell(current_row, 1)._style = copy(wb["S13"]["A6"]._style)
        s16.merge_cells(start_row=current_row, start_column=1, end_row=current_row, end_column=10)
        current_row += 1
        header_row = current_row
        headers = [str(c) for c in df.columns]
        if len(headers) != len(set(headers)) or any(not h for h in headers):
            raise RuntimeError(f"Invalid headers in {label}{sublabel}")
        for j, header in enumerate(headers, 1):
            cell = s16.cell(header_row, j, header)
            cell._style = copy(wb["S13"]["A8"]._style)
            cell.alignment = Alignment(wrap_text=True, vertical="center")
        s16.row_dimensions[header_row].height = 32
        for row_values in df.itertuples(index=False, name=None):
            current_row += 1
            for j, value in enumerate(row_values, 1):
                cell = s16.cell(current_row, j, clean_value(value))
                cell._style = copy(wb["S13"]["A9"]._style)
                if isinstance(cell.value, float):
                    cell.number_format = "0.###############"
        end_row = current_row
        end_col = len(headers)
        panel_counter += 1
        table = Table(displayName=f"S16_Panel_{panel_counter}", ref=f"A{header_row}:{get_column_letter(end_col)}{end_row}")
        table.tableStyleInfo = TableStyleInfo(name="TableStyleMedium2", showFirstColumn=False, showLastColumn=False, showRowStripes=True, showColumnStripes=False)
        s16.add_table(table)
        manifest_entries.append({
            "sheet": "S16",
            "panel": f"{label}{sublabel}. {subtitle}",
            "logical_panel": label,
            "rows": int(len(df)),
            "columns": int(end_col),
            "sources": source_list,
            "excel_table": table.name,
            "range": table.ref,
        })
        current_row += 2


add_panel("A", "Analysis-set, sample, and model-definition audit", "Three independently fitted analysis sets and the frozen two-by-two model definitions.", [
    ("", "Analysis-set and model-definition audit", panel_a, [rel(RUN / "analysis_set_audit.csv"), rel(PROTOCOL)]),
])
add_panel("B", "Complete six-configuration U macro-averaged metrics", "U applies each source-fitted mapping without labels from the held-out target; all expected targets must be evaluable for a macro summary.", [
    ("", "U macro-averaged metrics", panel_b, [rel(RUN / "rq5rf_macro_metrics.csv")]),
])
add_panel("C", "Target-level U metrics and prespecified pairwise differences", "Differences are left minus right. Negative Brier/log-loss differences and positive BSS/AUC/AP differences favour the left model.", [
    ("1", "Target-level U metrics", panel_c1, [rel(RUN / "rq5rf_target_metrics.csv")]),
    ("2", "Target-level prespecified pairwise differences", panel_c2, [rel(RUN / "rq5rf_pairwise_target_differences.csv")]),
])
add_panel("D", "Local-adaptation I/F results on identical outer-test rows", "I updates only the intercept; F updates the intercept and one calibration slope. These results use destination-target labels and remain separate from U.", [
    ("1", "I/F macro-averaged metrics", panel_d1, [rel(RUN / "rq5rf_macro_metrics.csv")]),
    ("2", "I/F target-level metrics", panel_d2, [rel(RUN / "rq5rf_target_metrics.csv")]),
])
add_panel("E", "Reliability summaries and non-evaluable diagnostics", "Blank numerical cells retain source-level non-evaluable quantities. The ledger records calibration diagnostics; unresolved failures remain zero.", [
    ("1", "Fixed-bin reliability summaries", panel_e1, [rel(RUN / "rq5rf_reliability.csv")]),
    ("2", "Non-evaluable and failure-item ledger", panel_e2, [rel(RUN / "failed_item_ledger.csv")]),
])
add_panel("F", "Candidate selection, implementation provenance, and QC", "All candidate selection and final-fit audits are retained. QC includes 3,492 random-forest fits, 16 unit tests, independent verification, and zero unresolved failures.", [
    ("1", "Inner-candidate metrics and selection", panel_f1, [rel(RUN / "rf_inner_candidate_selection.csv")]),
    ("2", "Final-candidate fit audit", panel_f2, [rel(RUN / "rf_final_candidate_fit_audit.csv")]),
    ("3", "QC, interpretation screen, and source hashes", panel_f3, [rel(RUN / "qc_checks.json"), rel(RUN / "independent_verification.json"), rel(RUN / "direction_consistency_10_2.json")]),
])

# Widths are capped for long provenance/JSON fields; full values remain in the cells.
for col in range(1, s16.max_column + 1):
    letter = get_column_letter(col)
    max_len = 0
    for row in range(1, s16.max_row + 1):
        value = s16.cell(row, col).value
        if value is not None:
            max_len = max(max_len, min(len(str(value)), 42))
    s16.column_dimensions[letter].width = max(11, min(42, max_len + 2))

wb.save(OUTPUT_WB)
if sha256(SOURCE_WB) != original_hash:
    raise RuntimeError("Original S1–S15 workbook changed unexpectedly.")

# Reopen to ensure the new workbook and named tables are structurally readable.
check = load_workbook(OUTPUT_WB, read_only=False, data_only=False)
assert check.sheetnames[-1] == "S16"
assert len(check["S16"].tables) == panel_counter
assert check["S16"].max_row > 5000
assert sum(1 for ws in check.worksheets for row in ws.iter_rows() for cell in row if cell.data_type == "f") == 0

# Keep the Markdown supplement guide synchronized with the integrated workbook.
guide = GUIDE.read_text(encoding="utf-8")
guide = guide.replace(
    "**Contents:** Supplementary Tables S1–S15 in one Excel workbook",
    "**Contents:** Supplementary Tables S1–S16 in one Excel workbook",
)
guide = guide.replace(
    "[Supplementary Tables S1–S15 (Excel workbook)](Supplementary_Tables_S1_S15.xlsx) contains Cover, Sources, and S1–S15 sheets.",
    "[Supplementary Tables S1–S16 (Excel workbook)](Supplementary_Tables_S1_S16.xlsx) contains Cover, Sources, and S1–S16 sheets.",
)
s16_section = """### Table S16. Random-forest model-complexity sensitivity on common complete cases

**Purpose.** Determine whether a more flexible functional form changes the cross-target probability-error and incremental-Pinc conclusions when samples, folds, feature sets, source weights, and local-adaptation rules are held fixed.

**Design.** A post hoc exploratory two-by-two comparison of pooled ridge logistic regression and random forests, each using ipTM alone or ipTM plus Pinc. Random-forest hyperparameters were selected within source training data using leave-one-source-target-out log loss. The six configurations comprised primary-six, primary-five, and strict-six analysis sets under target-equal and source-pair weighting. The primary-six, primary-five, and strict-six sets contained 2,074 pairs and 740 positives, 1,941 pairs and 721 positives, and 2,074 pairs and 674 positives, respectively. U used no destination-target labels; I and F used the corresponding frozen local-adaptation procedures.

**Evidence boundary.** This was not an algorithm search or a confirmatory model ranking. Only the prespecified random-forest family was added. The pooled logistic models omitted source-target intercepts; consequently, the earlier target-intercept RQ5 models are descriptive anchors rather than pure model-family comparators. Differences between ipTM-plus-Pinc and ipTM-only random forests combine possible nonlinear main effects and interactions and do not identify a pure interaction. No independent external Pinc cohort was available, and local-label I/F results cannot establish transport to held-out targets without their experimental labels.

**Workbook panels.** A. Analysis-set, sample, and model-definition audit; B. Complete six-configuration U macro-averaged metrics; C. Target-level U metrics and prespecified pairwise differences; D. Local-adaptation I/F results on identical outer-test rows; E. Reliability summaries and non-evaluable calibration diagnostics; F. Random-forest candidate selection, implementation provenance, and QC summary.

**Interpretation.** Negative Brier or log-loss differences and positive BSS, AUC, or AP differences favour the left-hand model in each prespecified contrast. The prespecified directional-consistency screen was not met. Results therefore support a statement of local or inconsistent changes under the evaluated models, not a claim that all nonlinear methods are ineffective or that random forests resolve calibration transport.

"""
if "### Table S16." not in guide:
    marker = "### Additional implementation notes for Table S13"
    if marker not in guide:
        raise RuntimeError("Could not locate the insertion point for Table S16 in the guide.")
    guide = guide.replace(marker, s16_section + marker, 1)
GUIDE.write_text(guide, encoding="utf-8")

build = json.loads(BUILD_MANIFEST.read_text(encoding="utf-8"))
build["generated"] = "2026-11-18"
build["workbook"] = OUTPUT_WB.name
build["expected_sheets"] = check.sheetnames
build["tables"] = [t for t in build["tables"] if t.get("sheet") != "S16"] + manifest_entries
build["s16_logical_panels"] = 6
build["s16_excel_tables"] = panel_counter
BUILD_MANIFEST.write_text(json.dumps(build, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

# Update source and generated-file hashes in the human-readable manifest.
manifest = SOURCE_MANIFEST.read_text(encoding="utf-8")
manifest = manifest.replace("Generated 16 November 2026.", "Updated 18 November 2026.")
new_source_rows = []
for path in source_paths:
    row = f"| `{rel(path)}` | {path.stat().st_size} | `{sha256(path)}` |"
    if f"| `{rel(path)}` |" not in manifest:
        new_source_rows.append(row)
manifest = manifest.replace("\n\n## Generated supplement files", "\n" + "\n".join(new_source_rows) + "\n\n## Generated supplement files")
lines = []
for line in manifest.splitlines():
    if line.startswith("| `Supplementary_Tables_S1_S15.xlsx` |") or line.startswith("| `Supplementary_Tables_S1_S16.xlsx` |"):
        line = f"| `{OUTPUT_WB.name}` | {OUTPUT_WB.stat().st_size} | `{sha256(OUTPUT_WB)}` |"
    elif line.startswith("| `Supplementary_Materials_PEDS_2026-11-16.md` |"):
        line = f"| `{GUIDE.name}` | {GUIDE.stat().st_size} | `{sha256(GUIDE)}` |"
    lines.append(line)
SOURCE_MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"Created {OUTPUT_WB}")
print(f"S16 rows={check['S16'].max_row}, columns={check['S16'].max_column}, logical panels=6, Excel tables={panel_counter}")
print(f"Workbook bytes={OUTPUT_WB.stat().st_size}, SHA-256={sha256(OUTPUT_WB)}")
print(f"Original workbook unchanged: {original_hash}")
print(f"Updated {BUILD_MANIFEST} and {SOURCE_MANIFEST}")
