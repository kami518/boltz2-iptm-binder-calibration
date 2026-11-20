from pathlib import Path

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
MAIN = ROOT / "manuscript/PEDS_Research_Article_English_draft_v2_2026-11-17.md"
REVIEW = ROOT / "manuscript/review/RQ5RF_results_and_supplement_note_for_review_2026-11-18.md"
SUPP = ROOT / "manuscript/supplement/Supplementary_Materials_PEDS_2026-11-16.md"


def replace_exact(text, old, new, expected=1):
    count = text.count(old)
    if count != expected:
        raise RuntimeError(f"Expected {expected} occurrence(s), found {count}: {old[:140]}")
    return text.replace(old, new)


main = MAIN.read_text(encoding="utf-8")
review = REVIEW.read_text(encoding="utf-8")
supp = SUPP.read_text(encoding="utf-8")

if "Random forests alter probability error but do not consistently outperform" in main:
    raise RuntimeError("RQ5-RF appears to be integrated already; stopping without changes.")

# Study-level provenance must acknowledge the explicitly post hoc model fits.
main = replace_exact(
    main,
    "No new wet-lab experiments were performed, and no model was refitted during manuscript preparation.",
    "No new wet-lab experiments were performed. Models for explicitly labelled post hoc extensions were fitted under separately frozen protocols; the original analysis outputs were not overwritten.",
)

methods = """### Random-forest model-complexity sensitivity

The random-forest analysis was added after the original Pinc results had been examined and was therefore post hoc and exploratory. A two-by-two comparison crossed pooled ridge logistic regression and random forests with ipTM alone or ipTM plus Pinc. All four models were refitted independently for the primary-six, primary-five, and strict-six analysis sets under target-equal and source-pair weighting, using identical complete-case rows and frozen outer folds. Source weights were recalculated within each current training row set after complete-case filtering and sequence-cluster exclusions. Unlike the original target-intercept RQ5 models, the pooled logistic comparators used one global intercept, allowing matched model-family comparisons without simultaneously changing the source-target intercept structure. The earlier target-intercept models were therefore retained as descriptive anchors rather than pure model-family comparators.

Random forests used scikit-learn 1.9.0 with 500 trees and a prespecified grid over maximum depth {2, 4, unrestricted}, minimum samples per leaf {10, 30}, and maximum features {1, 2}, with duplicate one-feature candidates removed. Hyperparameters were selected using leave-one-source-target-out logarithmic loss within each outer source-training set, with sequence-cluster exclusions and deterministic seeds. U used no destination-target labels. I and F applied the same frozen local intercept or intercept-plus-slope updates to each base model on identical outer-test rows. A prespecified directional-consistency screen required agreement across metrics, targets, analysis sets, and weighting schemes; it was a descriptive screen, not a significance test. No additional algorithm families were tried, and no independent external Pinc cohort was available. Complete parameters, candidate-level tuning, predictions, and provenance are assigned to Supplementary Table S16.

"""
main = replace_exact(main, "### Paired cross-predictor description", methods + "### Paired cross-predictor description")

result_start = review.index("### Random forests alter probability error")
result_end = review.index("## Proposed supplementary-table entry")
results = review[result_start:result_end].strip() + "\n\n"
main = replace_exact(main, "### Predictor pipelines differ in both score scale and design ordering", results + "### Predictor pipelines differ in both score scale and design ordering")

main = replace_exact(
    main,
    "Adding Pinc did not yield a stable increment across complete-case weighting schemes, despite the rationale for comparing interface-confidence measures (Badonyi and Toth-Petroczy, 2026).",
    "Adding Pinc did not yield a stable increment across complete-case weighting schemes, despite the rationale for comparing interface-confidence measures (Badonyi and Toth-Petroczy, 2026). The post hoc random-forest analysis produced metric- and target-dependent changes, but no consistent improvement over pooled logistic regression when predicting held-out targets without their experimental labels and no stable incremental contribution from Pinc.",
)
main = replace_exact(
    main,
    "Pinc analyses were complete-case, with missingness concentrated in PD-L1.",
    "Pinc analyses were complete-case, with missingness concentrated in PD-L1. The random-forest extension was post hoc, limited to the evaluated parameter grid, and lacked independent external Pinc validation; its combined-feature contrasts did not identify pure interaction effects.",
)

# Replace shorthand with language that states the information available at prediction time.
main_replacements = [
    ("a separately frozen external zero-label evaluation", "a separately frozen external evaluation without destination-target labels", 1),
    ("a descriptive zero-label stress test", "a descriptive stress test without destination-target labels", 1),
    ("the confirmatory zero-label evidence", "the confirmatory evidence obtained without destination-target labels", 1),
    ("the confirmatory zero-label result", "the confirmatory result obtained without destination-target labels", 1),
    ("are not zero-label evidence", "do not constitute evidence for prediction without destination-target labels", 1),
    ("external confirmatory zero-label U", "external confirmatory U without destination-target labels", 1),
    ("Zero-label transport and local adaptation", "Prediction without destination-target labels and local adaptation", 2),
    ("linear versus restricted-cubic-spline zero-label transport", "linear versus restricted-cubic-spline prediction without destination-target labels", 1),
]
for old, new, expected in main_replacements:
    main = replace_exact(main, old, new, expected)

main = replace_exact(main, "[Supplementary Tables S1–S15 workbook](supplement/Supplementary_Tables_S1_S15.xlsx)", "[Supplementary Tables S1–S16 workbook](supplement/Supplementary_Tables_S1_S16.xlsx)")
main = replace_exact(main, "15 supplementary-table sheets with 64 data panels", "16 supplementary-table sheets with 70 data panels")
main = replace_exact(
    main,
    "full predictor × target results; and reproducibility metadata with analysis-to-display mappings.",
    "full predictor × target results; reproducibility metadata with analysis-to-display mappings; and random-forest model-complexity sensitivity on common complete cases.",
)
main = replace_exact(main, "Supplementary Tables S1–S15 contain, respectively:", "Supplementary Tables S1–S16 contain, respectively:")
if "zero-label" in main:
    raise RuntimeError("Unreplaced zero-label wording remains in main manuscript.")

# Supplement guide: preserve existing numbering and add S16 after S15.
supp = replace_exact(supp, "**Scientific review draft — 16 November 2026**", "**Scientific review draft — updated 18 November 2026**")
supp = replace_exact(supp, "Supplementary Tables S1–S15 in one Excel workbook", "Supplementary Tables S1–S16 in one Excel workbook")
supp = replace_exact(supp, "[Supplementary Tables S1–S15 (Excel workbook)](Supplementary_Tables_S1_S15.xlsx) contains Cover, Sources, and S1–S15 sheets.", "[Supplementary Tables S1–S16 (Excel workbook)](Supplementary_Tables_S1_S16.xlsx) contains Cover, Sources, and S1–S16 sheets.")

s16_start = review.index("### Table S16.")
s16_end = review.index("## Numerical basis")
s16 = review[s16_start:s16_end].strip() + "\n\n"
supp = replace_exact(supp, "### Additional implementation notes for Table S13", s16 + "### Additional implementation notes for Table S13")

supp_replacements = [
    ("Ranking, absolute probability error, zero-label transport, label-dependent adaptation", "Ranking, absolute probability error, prediction without destination-target labels, label-dependent adaptation", 1),
    ("This is the zero-label transport comparator, not target-oracle calibration.", "This comparator applies when no destination-target labels are used; it is not target-oracle calibration.", 1),
    ("destination prevalence is unavailable in zero-label deployment", "destination prevalence is unavailable when deploying without destination-target labels", 1),
    ("without substituting a deployable zero-label reference", "without substituting a reference deployable without destination-target labels", 1),
    ("Test zero-label transport separately from label-dependent local adaptation", "Test prediction without destination-target labels separately from label-dependent local adaptation", 1),
    ("without changing the confirmatory zero-label result", "without changing the confirmatory result obtained without destination-target labels", 1),
    ("these results are exploratory, not zero-label validation", "these results are exploratory and do not validate prediction without destination-target labels", 1),
    ("Linear versus restricted-cubic-spline zero-label transport", "Linear versus restricted-cubic-spline prediction without destination-target labels", 2),
    ("the confirmatory zero-label U result", "the confirmatory external U result obtained without destination-target labels", 1),
    ("They do not establish zero-label success", "They do not establish successful prediction without destination-target labels", 1),
]
for old, new, expected in supp_replacements:
    supp = replace_exact(supp, old, new, expected)
if "zero-label" in supp:
    raise RuntimeError("Unreplaced zero-label wording remains in supplement guide.")

MAIN.write_text(main, encoding="utf-8")
SUPP.write_text(supp, encoding="utf-8")
print(f"Updated {MAIN}")
print(f"Updated {SUPP}")
print("Integrated Methods, Results, Discussion, limitations, Table S16 guide, and plain-language terminology.")
