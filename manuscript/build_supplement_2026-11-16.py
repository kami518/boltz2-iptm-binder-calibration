from pathlib import Path
import hashlib
import json
import math
import shutil
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo
from PIL import Image

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
OUT = ROOT / 'manuscript/supplement'
FIG = OUT / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
XLSX = OUT / 'Supplementary_Tables_S1_S15.xlsx'
DOC = OUT / 'Supplementary_Materials_PEDS_2026-11-16.md'
BUILD_MANIFEST = OUT / 'supplement_build_manifest.json'
SOURCE_MANIFEST = OUT / 'supplement_source_manifest_2026-11-16.md'

RQ3 = 'analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311'
RQ3E = 'analysis/rq3e_outputs/run_20261114_164950'
LOW = 'analysis/rq3e_outputs/run_20261114_162656'
RQ3EIF = 'analysis/rq3e_if_outputs/run_20261114_174403'
RQ4 = 'analysis/rq4_outputs/run_20261114_231925'
RQ5 = 'analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720'
RQ6 = 'analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147'
OVERATH = 'analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301'
PINC = 'analysis/pinc_batch/run_20261110_181500'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.titlesize': 10.5, 'axes.labelsize': 9.5,
    'figure.titlesize': 14, 'axes.spines.top': False,
    'axes.spines.right': False, 'savefig.facecolor': 'white',
})
NAVY = '#183B56'
BLUE = '#2878B5'
TEAL = '#1B9E91'
ORANGE = '#E07A3F'
RED = '#C44E52'
GREY = '#687684'
LIGHT = '#EAF0F4'
PURPLE = '#8064A2'
TARGETS = ['aah-i', 'human-vegf110', 'nipah-glycoprotein-g', 'pd-l1', 'sars-cov-2-rbd', 'trem2']
TARGET_LABELS = {
    'aah-i': 'AAH-I', 'human-vegf110': 'VEGF110',
    'nipah-glycoprotein-g': 'Nipah G', 'pd-l1': 'PD-L1',
    'sars-cov-2-rbd': 'SARS-CoV-2 RBD', 'trem2': 'TREM2'
}
SET_LABELS = {
    'primary_six': 'Primary, six targets',
    'primary_five': 'Primary, five targets\n(excluding AAH-I)',
    'strict_six': 'Strict, six targets',
    'strict_five': 'Strict, five targets\n(excluding AAH-I)',
}

used_sources = {}
table_records = []
figure_records = []


def abspath(rel):
    return ROOT / rel


def read_csv(rel, **kwargs):
    path = abspath(rel)
    if not path.is_file():
        raise FileNotFoundError(path)
    used_sources[str(path)] = rel
    return pd.read_csv(path, **kwargs)


def read_json(rel):
    path = abspath(rel)
    if not path.is_file():
        raise FileNotFoundError(path)
    used_sources[str(path)] = rel
    return json.loads(path.read_text(encoding='utf-8'))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def flatten_json(obj, prefix=''):
    rows = []
    if isinstance(obj, dict):
        if not obj:
            rows.append({'field': prefix, 'value': '{}'})
        for key, value in obj.items():
            path = f'{prefix}.{key}' if prefix else str(key)
            rows.extend(flatten_json(value, path))
    elif isinstance(obj, list):
        if not obj:
            rows.append({'field': prefix, 'value': '[]'})
        for idx, value in enumerate(obj):
            rows.extend(flatten_json(value, f'{prefix}[{idx}]'))
    else:
        if obj is None:
            value = ''
        elif isinstance(obj, bool):
            value = str(obj).lower()
        else:
            value = obj
        rows.append({'field': prefix, 'value': value})
    return rows


def json_frame(rel):
    return pd.DataFrame(flatten_json(read_json(rel)))


def rule_frame(rows):
    return pd.DataFrame(rows, columns=['Item', 'Operational definition or scope', 'Interpretation boundary'])


sheets = {}


def register_sheet(sheet, title, purpose, boundary):
    sheets[sheet] = {'title': title, 'purpose': purpose, 'boundary': boundary, 'panels': []}


def panel(sheet, name, df, note='', sources=None):
    if df is None or len(df) == 0:
        raise ValueError(f'Empty panel: {sheet} {name}')
    df = df.copy()
    sheets[sheet]['panels'].append({'name': name, 'df': df, 'note': note})
    table_records.append({
        'sheet': sheet, 'panel': name, 'rows': int(len(df)),
        'columns': int(len(df.columns)), 'sources': sources or []
    })


# S1: outcome definitions and full pair-level label audit.
pair = read_csv('analysis/v3_outputs/pair_master_v3.csv')
assert len(pair) == 2084
assert int(pair.primary_outcome.sum()) == 748
assert int(pair.strict_outcome.sum()) == 678
register_sheet(
    'S1', 'Outcome rules and alternative-label audit',
    'Make the binding endpoint and its denominator auditable before comparing score performance.',
    'Alternative labels expose definition sensitivity; they are not parallel models or independent experimental replicates.'
)
panel('S1', 'A. Frozen operational outcome rules', rule_frame([
    ('Analysis unit', 'One design–target pair with one prespecified primary assay.', 'Every design occurs once in the primary six-target table; no design-level random intercept is identifiable.'),
    ('Exact-record deduplication', 'Deduplicate exact records within design × target × primary assay.', 'Distinct timestamps are not assumed to represent independent experimental replicates.'),
    ('Primary outcome', 'Majority binding after exact-record deduplication, including Weak calls.', 'The estimand is an experimental binding call, not continuous affinity.'),
    ('Strict outcome', 'A primary positive whose highest observed deduplicated strength is Medium or Strong.', 'Primary negatives cannot be promoted to strict positives.'),
    ('All-deduplicated alternative', 'Positive only when all deduplicated primary-assay binding records are positive.', 'Definition-sensitivity audit only.'),
    ('Any-deduplicated alternative', 'Positive when any deduplicated primary-assay binding record is positive.', 'Definition-sensitivity audit only.'),
    ('If-concordant alternative', 'Retains positive calls only under the frozen concordance rule in the pair-level table.', 'Definition-sensitivity audit only; not an independent assay replication.'),
    ('Pooled total', '2,084 pairs, 748 primary positives, and 678 strict positives across six targets.', 'Counts describe this frozen retrospective cohort, not a future target population.'),
]), 'Definitions reproduce the frozen Methods; no endpoint was selected by prediction performance.', ['analysis/v3_outputs/pair_master_v3.csv'])
count_rows = []
for target, group in pair.groupby('target', sort=True):
    count_rows.append({
        'target': target, 'n_pairs': len(group),
        'primary_positives': int(group.primary_outcome.sum()),
        'strict_positives': int(group.strict_outcome.sum()),
        'all_deduplicated_positives': int(group.primary_outcome_all_deduplicated.sum()),
        'any_deduplicated_positives': int(group.primary_outcome_any_deduplicated.sum()),
        'if_concordant_positives': int(group.primary_outcome_if_concordant.sum()),
        'raw_discordant_pairs': int(group.primary_label_discordant_raw.sum()),
        'deduplicated_discordant_pairs': int(group.primary_label_discordant_deduplicated.sum()),
        'strict_downgrades_from_primary_positive': int(group.strict_downgrade_from_primary_positive.sum()),
    })
panel('S1', 'B. Target-level outcome audit', pd.DataFrame(count_rows),
      'Arithmetic tabulation of the frozen pair-level table; no model was refitted.', ['analysis/v3_outputs/pair_master_v3.csv'])
panel('S1', 'C. Complete pair-level outcome and discordance fields', pair,
      'All 39 frozen fields are retained, including primary, strict, alternative, strength, and discordance fields. Redistribution remains subject to source-data licensing.', ['analysis/v3_outputs/pair_master_v3.csv'])

# S2: assay composition and selective confirmation.
assay_audit = read_csv('analysis/v3_outputs/assay_unit_audit_v3.csv')
assay_counts = read_csv('analysis/audit_outputs/target_actual_assay_design_counts.csv')
register_sheet(
    'S2', 'Assay composition and AAH-I selective confirmation',
    'Show what was measured, on which denominator, and why AAH-I cannot support an assay-robustness claim.',
    'Primary assays were selected by measurement coverage, not by predictive strength. Target, assay, and confirmation structure remain associated.'
)
panel('S2', 'A. Assay-unit audit', assay_audit,
      'Complete frozen assay-unit audit.', ['analysis/v3_outputs/assay_unit_audit_v3.csv'])
panel('S2', 'B. Actual assay and design counts by target', assay_counts,
      'Source counts are retained without collapsing distinct assay contexts.', ['analysis/audit_outputs/target_actual_assay_design_counts.csv'])
panel('S2', 'C. AAH-I confirmation structure', rule_frame([
    ('Primary SPR denominator', 'All 133 AAH-I designs had SPR outcomes; SPR is the primary assay.', 'Supports a complete primary-outcome denominator.'),
    ('Selective BLI subset', '18 selectively acquired BLI confirmations.', 'Not a random or complete assay comparison subset.'),
    ('Joint calls', 'All 18 BLI-positive confirmations were also SPR-positive.', 'The 19 primary positives are not equivalent to a simple SPR-only discovery structure.'),
    ('Required sensitivity analysis', 'Report AAH-I separately and repeat target-level analyses after complete refitting without AAH-I.', 'Exclusion is a sensitivity analysis, not a new primary cohort.'),
]), 'The selective confirmation structure motivates separate reporting; it does not establish which instrument is more accurate.', ['analysis/v3_outputs/pair_master_v3.csv', 'analysis/v3_outputs/assay_unit_audit_v3.csv'])

# S3: complete main-cohort discrimination and raw calibration.
rq12_target = read_csv('analysis/rq12_outputs/rq12_target_metrics.csv')
rq12_summary = read_csv('analysis/rq12_outputs/rq12_summary_metrics.csv')
rq12_reliability = read_csv('analysis/rq12_outputs/rq2_reliability_bins.csv')
register_sheet(
    'S3', 'Main-cohort discrimination and raw probability error',
    'Separate within-target ordering from the stronger claim that raw ipTM is an experimental-success probability.',
    'AP and Brier-based measures depend on prevalence. Raw-ipTM probability interpretation is a stress test, not a claim about the model developers’ intended meaning.'
)
panel('S3', 'A. Complete target-level metrics and bootstrap intervals', rq12_target,
      'Primary/strict and six-/five-target configurations are retained at source precision.', ['analysis/rq12_outputs/rq12_target_metrics.csv'])
panel('S3', 'B. Target-equal and pair-weighted summaries', rq12_summary,
      'Intervals quantify conditional precision in observed targets, not uncertainty over future targets.', ['analysis/rq12_outputs/rq12_summary_metrics.csv'])
panel('S3', 'C. Frozen common-bin reliability results', rq12_reliability,
      'Common bins support visual comparison; binned curves are diagnostics rather than a fitted probability mapping.', ['analysis/rq12_outputs/rq2_reliability_bins.csv'])

# S4: hierarchical heterogeneity and numerical diagnostics.
register_sheet(
    'S4', 'Hierarchical calibration heterogeneity and numerical diagnostics',
    'Quantify target-dependent intercept/slope heterogeneity while separating description from transport validation.',
    'The hierarchy is descriptive and is not a validated mapping for an unseen target. Six contexts cannot identify independent target-protein covariate effects.'
)
panel('S4', 'A. Target-specific unpooled and partially pooled estimates', read_csv('analysis/rq12_outputs/rq2_hierarchical_targets.csv'),
      'All primary/strict and six-/five-target models are retained.', ['analysis/rq12_outputs/rq2_hierarchical_targets.csv'])
panel('S4', 'B. Hierarchical hyperparameters', read_csv('analysis/rq12_outputs/rq2_hierarchical_hyperparameters.csv'),
      'Between-target estimates describe heterogeneity in the fitted models.', ['analysis/rq12_outputs/rq2_hierarchical_hyperparameters.csv'])
panel('S4', 'C. Full hierarchical model output', json_frame('analysis/rq12_outputs/rq2_hierarchical_full_results.json'),
      'Nested JSON is flattened to field/value rows without changing values.', ['analysis/rq12_outputs/rq2_hierarchical_full_results.json'])
panel('S4', 'D. Adaptive quadrature convergence diagnostic', json_frame('analysis/rq12_outputs/rq12_quadrature_convergence_diagnostic.json'),
      'Adaptive Gauss–Hermite 9/11-node check; maximum target-estimate difference 2.71×10⁻⁶.', ['analysis/rq12_outputs/rq12_quadrature_convergence_diagnostic.json'])
panel('S4', 'E. Frozen RQ1/RQ2 analysis configuration', json_frame('analysis/rq12_outputs/rq12_analysis_config.json'),
      'Includes clipping, reliability bins, seeds, bootstrap rules, priors, and fallback criteria.', ['analysis/rq12_outputs/rq12_analysis_config.json'])

# S5: main-cohort target-blocked transport and local adaptation.
register_sheet(
    'S5', 'Main-cohort LOTO U/I/F and sensitivity analyses',
    'Test zero-label transport separately from label-dependent local adaptation and identify what local measurements correct.',
    'U uses no destination-target labels. Predictive I/F use local labels and out-of-fold evaluation; descriptive/oracle fits are not deployable evidence.'
)
for name, filename, note in [
    ('A. Complete LOTO target-level results', 'rq3_loto_metrics.csv', 'Includes linear and restricted-cubic-spline models and all frozen sensitivity configurations.'),
    ('B. LOTO macro summaries', 'rq3_loto_macro_summary.csv', 'BSS references are named separately; target-oracle is diagnostic, not deployable.'),
    ('C. Target-level local-adaptation results', 'rq3_target_adaptation_metrics.csv', 'Predictive I/F use identical outer rows within each target.'),
    ('D. Paired adaptation comparisons', 'rq3_adaptation_comparisons.csv', 'Negative Brier/log-loss differences favor the named left layer.'),
    ('E. Adaptation macro summaries', 'rq3_adaptation_macro_summary.csv', 'Counts of improving targets are descriptive, not significance tests.'),
    ('F. Complete sensitivity summary', 'rq3_sensitivity_summary.csv', 'Includes strict outcome, five-target refitting, source weighting, nonlinear form, and cluster sensitivity where prespecified.'),
]:
    panel('S5', name, read_csv(f'{RQ3}/{filename}'), note, [f'{RQ3}/{filename}'])

# S6: exact and near-sequence audits.
register_sheet(
    'S6', 'Overath exact-sequence and near-sequence audits',
    'Distinguish sequence non-overlap from temporal, training-data, structural, or epitope independence.',
    'The accessible mirror has 3,669 records versus 3,766 in the article. Short local alignments can have high identity; identity must be interpreted with bidirectional coverage.'
)
panel('S6', 'A. Exact-sequence overlap audit', read_csv(f'{OVERATH}/overath_exact_sequence_overlap.csv'),
      'Full-source audit: 345/4,280 designs and 342/4,943 pairs; frozen Boltz-2 main cohort: zero exact matches.', [f'{OVERATH}/overath_exact_sequence_overlap.csv'])
panel('S6', 'B. Maximum similarity per unique sequence', read_csv(f'{OVERATH}/overath_sequence_similarity_unique.csv'),
      'Retains maximum-identity distributions and coverage information instead of only a binary threshold indicator.', [f'{OVERATH}/overath_sequence_similarity_unique.csv'])
panel('S6', 'C. Sequence-similarity run summary', json_frame(f'{OVERATH}/overath_sequence_similarity_summary.json'),
      '7,646,196 comparisons; no 80% or 90% identity links with query and reference coverage both ≥80%.', [f'{OVERATH}/overath_sequence_similarity_summary.json'])
panel('S6', 'D. Mirror/article version reconciliation', json_frame(f'{OVERATH}/overath_version_reconciliation.json'),
      'The 97-record discrepancy remains explicit; observed percentages are not called lower bounds for the complete original dataset.', [f'{OVERATH}/overath_version_reconciliation.json'])

# S7: source/campaign/generator descriptions.
register_sheet(
    'S7', 'Campaign, source, and generator composition',
    'Make changed deployment contexts visible rather than attributing all heterogeneity to target biology.',
    'Generator, target, campaign, design model, assay, and score protocol are associated; independent causal effects are not identifiable.'
)
panel('S7', 'A. Main-cohort pair-level source fields', pair[['pair_id', 'design_id', 'author', 'designMethod', 'design_class', 'target', 'source_row']],
      'Original source fields are supplied without inventing a new causal grouping.', ['analysis/v3_outputs/pair_master_v3.csv'])
panel('S7', 'B. External design-model × generator composition', read_csv(f'{RQ3E}/rq3e_generator_composition.csv'),
      'Mythos Preview: 480 designs; Opus 4.8: 240. Nine generators are represented in the 720-design external cohort.', [f'{RQ3E}/rq3e_generator_composition.csv'])
panel('S7', 'C. External campaign/source descriptions', read_csv(f'{RQ3E}/rq3e_campaign_metrics.csv'),
      'Campaign/source descriptions are not independent replicate cohorts or causal effect estimates.', [f'{RQ3E}/rq3e_campaign_metrics.csv'])

# S8: external U and low-event stress test.
rq3e_target = read_csv(f'{RQ3E}/rq3e_target_metrics.csv')
rq3e_group = read_csv(f'{RQ3E}/rq3e_group_metrics.csv')
low_metrics = read_csv(f'{LOW}/rq3e_low_event_stress_test_metrics.csv')
assert len(low_metrics) == 6 and int(low_metrics.n.sum()) == 521
register_sheet(
    'S8', 'Frozen external U and nonconfirmatory low-event stress test',
    'Evaluate a frozen released-score-plus-mapping pipeline without external-label refitting and expose the prevalence range actually tested.',
    'External target/campaign and score-generation conditions change together. Low-event results are descriptive and do not extend confirmatory coverage.'
)
panel('S8', 'A. Primary external U target-level results', rq3e_target[rq3e_target.analysis_role.eq('primary')],
      'Primary Adaptyv Weak+, best-ipSAE-seed score rule, target-equal source weighting.', [f'{RQ3E}/rq3e_target_metrics.csv'])
panel('S8', 'B. Primary external U grouped results', rq3e_group[rq3e_group.analysis_role.eq('primary')],
      'Includes all eight targets, four new targets, four same-target new campaigns, and pair-pooled diagnostics.', [f'{RQ3E}/rq3e_group_metrics.csv'])
panel('S8', 'C. Low-event stress-test target metrics', low_metrics,
      'Six targets, 521 pairs. Zero/one-event targets do not support substantive slope inference; non-evaluable cells remain blank with source status fields.', [f'{LOW}/rq3e_low_event_stress_test_metrics.csv'])
panel('S8', 'D. Low-event stress-test reliability bins', read_csv(f'{LOW}/rq3e_low_event_stress_test_reliability.csv'),
      'Final low-event run 20261114_162656; earlier-run images are not used.', [f'{LOW}/rq3e_low_event_stress_test_reliability.csv'])

# S9: external sensitivity configurations.
register_sheet(
    'S9', 'External outcome, score-rule, and source-weight sensitivity',
    'Test whether plausible outcome, score-summary, seed, or source-weight choices remove context dependence.',
    'Sensitivity configurations are not competing results from which to select a preferred calibration claim. Each outcome/provider view has its own experimental meaning.'
)
panel('S9', 'A. All nonprimary group-level configurations', rq3e_group[~rq3e_group.analysis_role.eq('primary')],
      'analysis_role, outcome_name, score_rule, source_weighting, aggregation, and stratum identify every configuration.', [f'{RQ3E}/rq3e_group_metrics.csv'])
panel('S9', 'B. All nonprimary target-level configurations', rq3e_target[~rq3e_target.analysis_role.eq('primary')],
      'Target-level values preserve where sensitivity analyses help or fail; non-evaluable fields retain their status and message.', [f'{RQ3E}/rq3e_target_metrics.csv'])
panel('S9', 'C. Interpretation guide', rule_frame([
    ('Outcome sensitivity', 'Adaptyv Medium+, integrated, Twist QC, and Twist permissive outcomes are separate outcome views.', 'Alternative providers are not independent experimental replicates.'),
    ('Score sensitivity', 'Five-seed mean/median and fixed seed 0 are compared with the best-ipSAE-seed rule.', 'The primary rule selects maximum ipSAE_min and uses that seed’s iptm_pae; it is not maximum ipTM.'),
    ('Source weighting', 'Target-equal and source-pair weighting examine sensitivity to the large Nipah G contribution.', 'A favorable BSS-transport under one weighting does not establish calibration or positive BSS-target-oracle.'),
]), 'This panel records the interpretation rules; it does not add a new analysis.', [f'{RQ3E}/run_config.snapshot.json'])

# S10: post hoc external I/F.
register_sheet(
    'S10', 'Post hoc exploratory external I/F adaptation',
    'Assess what representative local labels can correct in the external context without changing the confirmatory zero-label result.',
    'Added after external U was examined. Five-fold outer predictions and nested tuning use local labels; these results are exploratory, not zero-label validation.'
)
for name, filename, note in [
    ('A. Grouped U/I/F metrics', 'group_metrics.csv', 'All eight targets and transport strata; U/I/F are evaluated on identical outer rows.'),
    ('B. Target-level paired comparisons', 'target_comparisons.csv', 'Negative Brier/log-loss differences favor the named left layer.'),
    ('C. Target-level U/I/F metrics', 'target_metrics.csv', 'Small target samples limit slope precision; direction counts are descriptive.'),
    ('D. Outer-fold metrics', 'fold_metrics.csv', 'Outer-fold values expose heterogeneity hidden by target-level averages.'),
    ('E. Nested regularization selection', 'lambda_selection.csv', 'Lambda selection is confined to adaptation-training folds.'),
]:
    panel('S10', name, read_csv(f'{RQ3EIF}/{filename}'), note, [f'{RQ3EIF}/{filename}'])

# S11: fixed budgets and raw thresholds.
register_sheet(
    'S11', 'Fixed-budget selection and raw-ipTM operating points',
    'Translate ranking into experimental selections and show why a common numerical score is not a common operating point.',
    'The 0.80/0.90/0.95 thresholds are descriptive operating points, not calibrated probabilities or recommended universal cutoffs.'
)
panel('S11', 'A. Complete target-level fixed-budget and threshold metrics', read_csv(f'{RQ4}/rq4ab_target_metrics.csv'),
      'Includes proportional/absolute budgets, empty selections, strict/five-target and miniprotein-scope configurations.', [f'{RQ4}/rq4ab_target_metrics.csv'])
panel('S11', 'B. Pooled and target-equal summaries', read_csv(f'{RQ4}/rq4ab_summary.csv'),
      'Precision, recall, enrichment, hit gain, and cost per hit are distinct outcomes. Non-evaluable precision/enrichment remain blank.', [f'{RQ4}/rq4ab_summary.csv'])

# S12: shared-budget retrospective simulations.
register_sheet(
    'S12', 'Post hoc exploratory retrospective shared-budget simulations',
    'Charge pilot assays to the same experimental budget and test when local information changes total confirmed hits.',
    'Observed-pool retrospective simulation only. Intervals are simulation ranges, not future-target confidence intervals; no universal pilot size or prospective saving is established.'
)
panel('S12', 'A. Complete strategy-level B/m results', read_csv(f'{RQ4}/rq4c_summary.csv'),
      'Covers primary_six, primary_five, and strict_six. Primary six-target analyses use 2,000 replicates; sensitivity analyses use 1,000.', [f'{RQ4}/rq4c_summary.csv'])
panel('S12', 'B. Complete paired strategy comparisons', read_csv(f'{RQ4}/rq4c_comparisons.csv'),
      'Pilot assays and pilot-identified hits count toward budget and hits; already measured pairs are not selected again.', [f'{RQ4}/rq4c_comparisons.csv'])
panel('S12', 'C. Pilot parameter audit', read_csv(f'{RQ4}/rq4c_pilot_parameters.csv'),
      'Small-pilot models use prespecified smoothing/regularization and are distinct from larger-fold transport I/F models.', [f'{RQ4}/rq4c_pilot_parameters.csv'])

# S13: Pinc complete-case and implementation provenance.
register_sheet(
    'S13', 'Pinc complete-case analyses and implementation provenance',
    'Test incremental value on identical rows and folds without converting a context-specific comparison into a universal metric ranking.',
    '2,074 complete-case pairs and 740 primary positives; ten missing Pinc values are all in PD-L1 and are not imputed. Chain-specific Pinc1/Pinc2 are not interpreted.'
)
for name, filename, note in [
    ('A. Same-row target score metrics', 'rq5_target_score_metrics.csv', 'Raw ipTM and Pinc are compared on common complete cases.'),
    ('B. Same-row score macro summaries', 'rq5_score_macro_summary.csv', 'Lower raw Brier does not establish calibration.'),
    ('C. LOTO U target metrics', 'rq5_u_target_metrics.csv', 'ipTM-only, Pinc-only, and ipTM+Pinc are refitted on identical rows and folds.'),
    ('D. LOTO U macro summaries', 'rq5_u_macro_summary.csv', 'Source weighting is retained as a separate dimension.'),
    ('E. Complete pairwise comparisons', 'rq5_pairwise_comparisons.csv', 'Comparisons follow the stored direction convention; no universal superiority claim is made.'),
    ('F. Local-adaptation pair summaries', 'rq5_adaptation_pair_summary.csv', 'Only local calibration parameters are updated; relative feature weights are not silently re-estimated.'),
    ('G. Complete Pinc sensitivity summary', 'rq5_sensitivity_summary.csv', 'Includes analysis set, outcome, cluster rule, source weighting, summary type, and paired contrast fields.'),
]:
    panel('S13', name, read_csv(f'{RQ5}/{filename}'), note, [f'{RQ5}/{filename}'])
panel('S13', 'H. Pinc pair-level extension', read_csv(f'{PINC}/pinc_pair_extension.csv'),
      'Frozen main table was not modified; Pinc is supplied as a separate extension.', [f'{PINC}/pinc_pair_extension.csv'])
panel('S13', 'I. Pair-processing manifest', read_csv(f'{PINC}/pair_processing_manifest.csv'),
      'Retains success/missing-input status and implementation provenance for every processed pair.', [f'{PINC}/pair_processing_manifest.csv'])
panel('S13', 'J. Batch summary', json_frame(f'{PINC}/pinc_batch_summary.json'),
      'Official legacy R was accepted for batch throughput after small-sample agreement with compiled official C; source/script hashes are preserved.', [f'{PINC}/pinc_batch_summary.json'])
CLOSURE = 'analysis/codex_outputs/C03_official_pinc_handoff_2026-11-09/run_20261109_220449/closure/closure_result.json'
panel('S13', 'K. Official C/R closure and PAE headline status', json_frame(CLOSURE),
      'Small-sample official C/R agreement supports the accepted R batch toolchain. Only symmetric Pinc is used; chain-specific interpretation remains open.', [CLOSURE])

# S14: paired predictor-pipeline description.
rq6_summary = read_csv(f'{RQ6}/rq6_target_equal_summary.csv')
assert len(rq6_summary) == 180
register_sheet(
    'S14', 'Paired description across ten predictor pipelines',
    'Determine whether alternative predictor outputs merely shift the numerical scale or also change within-target design ordering.',
    'Same 720 external designs and outcomes; no main-cohort cross-model validation and no U/I/F transport across predictors. Complete pipelines differ in more than weights.'
)
for name, filename, note in [
    ('A. Target-equal predictor × outcome summaries', 'rq6_target_equal_summary.csv', 'Means, medians, and ranges are descriptive; ranks are point-estimate descriptions, not superiority tests.'),
    ('B. Complete predictor × target metrics', 'rq6_metrics_by_target.csv', 'Both Adaptyv Weak+ and Medium+ outcomes are retained.'),
    ('C. Score shifts and within-target agreement', 'rq6_score_shift_by_target.csv', 'Spearman agreement tests whether candidate ordering changes as well as score scale.'),
    ('D. Paired bootstrap intervals', 'rq6_bootstrap_intervals.csv', 'Design-level paired bootstrap within target; intervals are conditional on observed targets.'),
    ('E. Selected run/model audit', 'rq6_selected_run_audit.csv', 'Ten complete selected pipelines per design. af3of3 is AF3 architecture code with OpenFold3 weights.'),
]:
    panel('S14', name, read_csv(f'{RQ6}/{filename}'), note, [f'{RQ6}/{filename}'])

# S15: reproducibility metadata and analysis-to-display mapping.
register_sheet(
    'S15', 'Reproducibility metadata and analysis-to-display mapping',
    'Connect every scientific display to its frozen numerical source, run configuration, and file hashes.',
    'This is a local reproducibility draft. Source-release versions, citations, and redistribution terms are verified. The public repository/DOI and author metadata remain unresolved.'
)
manifest_frames = []
manifest_paths = [
    f'{RQ3}/final_hash_manifest.csv', f'{RQ3E}/final_hash_manifest.csv',
    f'{RQ3EIF}/final_hash_manifest.csv', f'{RQ4}/final_hash_manifest.csv',
    f'{RQ5}/final_hash_manifest.csv', f'{RQ6}/final_hash_manifest.csv',
]
for rel in manifest_paths:
    frame = read_csv(rel)
    frame.insert(0, 'run_manifest_source', rel)
    manifest_frames.append(frame)
manifest_all = pd.concat(manifest_frames, ignore_index=True, sort=False)
panel('S15', 'A. Final hash manifests from accepted runs', manifest_all,
      'Run-specific columns are preserved in a combined long table.', manifest_paths)
config_frames = []
config_paths = [
    'analysis/rq12_outputs/rq12_run_manifest.json',
    f'{RQ3}/run_config.snapshot.json', f'{RQ3E}/run_config.snapshot.json',
    f'{RQ3EIF}/run_config.snapshot.json', f'{RQ4}/run_config.snapshot.json',
    f'{RQ5}/run_config.snapshot.json', f'{RQ6}/run_manifest.json',
]
for rel in config_paths:
    frame = json_frame(rel)
    frame.insert(0, 'configuration_source', rel)
    config_frames.append(frame)
panel('S15', 'B. Run manifests and configuration snapshots', pd.concat(config_frames, ignore_index=True),
      'Nested configurations are flattened to source/field/value rows; no missing public link is invented.', config_paths)
display_rows = [
    ('Main Table 1', 'Main/external pair and composition tables', 'S1, S2, S7', 'Cohort denominators and assay/outcome contexts.'),
    ('Main Table 2', 'analysis/rq12_outputs/rq12_target_metrics.csv; rq12_summary_metrics.csv', 'S3; Supplementary Figure S1', 'Discrimination and raw probability error.'),
    ('Main Table 3A', f'{RQ3}/rq3_sensitivity_summary.csv and adaptation tables', 'S5; Supplementary Figure S2', 'Main LOTO U and predictive local I/F.'),
    ('Main Table 3B', f'{RQ3E}/rq3e_group_metrics.csv; {RQ3EIF}/group_metrics.csv', 'S8–S10; Supplementary Figures S3–S6', 'External confirmatory U and post hoc I/F remain separate.'),
    ('Main Table 4', f'{RQ4}/rq4c_summary.csv; rq4c_comparisons.csv', 'S12; Supplementary Figure S7', 'Post hoc shared-budget simulation.'),
    ('Main Figure 1', 'manuscript/figures/figure_source_manifest_2026-11-15.md', 'S1, S7–S8, S15', 'Study design and cohort flow.'),
    ('Main Figure 2', 'analysis/rq12_outputs/RQ1_ipTM_distributions.png; RQ1_target_discrimination.png', 'S3', 'ipTM distributions and discrimination.'),
    ('Main Figure 3', 'analysis/rq12_outputs/RQ2_raw_reliability.png; RQ2_partial_pooling.png', 'S3–S4', 'Raw calibration and descriptive hierarchy.'),
    ('Main Figure 4', f'{RQ3}; {RQ3E}', 'S5, S8–S10', 'Main/external U and main predictive I/F; post hoc external I/F is supplementary.'),
    ('Main Figure 5', f'{RQ4}', 'S11–S12', 'Experimental selection and shared-budget consequences.'),
    ('Main Figure 6', f'{RQ5}', 'S13; Supplementary Figures S8–S9', 'Pinc complete-case analyses.'),
    ('Main Figure 7', f'{RQ6}; analysis/audit_outputs/RQ6_acceptance_run_20261115_151147', 'S14; Supplementary Figure S10', 'Paired complete-pipeline description.'),
    ('Supplementary Figure S1', 'analysis/rq12_outputs/rq12_target_metrics.csv; rq12_summary_metrics.csv', 'S3', 'Strict/five-target sensitivity.'),
    ('Supplementary Figure S2', f'{RQ3}/rq3_loto_metrics.csv; rq3_loto_macro_summary.csv', 'S5', 'Linear versus restricted-cubic-spline U.'),
    ('Supplementary Figures S3–S4', f'{RQ3EIF}/RQ3E_IF_OOF_metrics.png; RQ3E_IF_OOF_reliability.png', 'S10', 'Post hoc external I/F.'),
    ('Supplementary Figure S5', f'{LOW}/RQ3E_U_low_event_stress_test.png', 'S8', 'Final low-event stress-test run only.'),
    ('Supplementary Figure S6', f'{RQ3E}/RQ3E_U_external_sensitivity.png', 'S9', 'External outcome/score/source-weight sensitivity.'),
    ('Supplementary Figure S7', f'{RQ4}/rq4c_comparisons.csv', 'S12', 'Shared-budget sensitivity.'),
    ('Supplementary Figure S8', f'{RQ5}/RQ5_raw_reliability.png', 'S13', 'Raw Pinc and ipTM reliability.'),
    ('Supplementary Figure S9', f'{RQ5}/rq5_sensitivity_summary.csv', 'S13', 'Pinc increment by analysis set and source weighting.'),
    ('Supplementary Figure S10', f'{RQ6}/rq6_target_equal_summary.csv', 'S14', 'Ten-pipeline target-equal descriptive summary.'),
]
panel('S15', 'C. Analysis-to-display mapping', pd.DataFrame(display_rows, columns=['Display', 'Frozen source or run', 'Related supplement', 'Evidence purpose']),
      'This mapping is a navigation aid; the numerical sources and accepted run manifests remain authoritative.', [])
panel('S15', 'D. Source-release and submission metadata status', rule_frame([
    ('Proteinbase version', 'v2026-11-01; Hugging Face revision 2d2825070da5ba77e0c7e1a0942146e8307e9282.', 'Local CSV saved 5 November 2026; analysis access recorded 7 November 2026.'),
    ('External release citation', 'Claude protein binder design data release v1.0; revision 9e1b81696da46835e9e9cde9a3da976e0abc92ab.', 'Official citation author: The protein binder design data release authors.'),
    ('Public code/data repository', 'Public repository and DOI remain to be supplied.', 'Local manifests are not a public reproducibility deposit.'),
    ('Redistribution', 'Proteinbase: ODC-BY with required attribution. External release data/documentation: CC BY 4.0; included scripts: MIT.', 'Third-party structure outputs, sequences, vendor images, and other third-party material retain their own terms.'),
]), 'Verified source metadata is separated from the remaining public-repository and authorship items.', [])

assert set(sheets) == {f'S{i}' for i in range(1, 16)}

# Build styled workbook with multiple named Excel tables per sheet.
if XLSX.exists():
    XLSX.unlink()
with pd.ExcelWriter(XLSX, engine='openpyxl') as writer:
    cover = pd.DataFrame([
        ('Title', 'Supplementary Tables S1–S15: context-dependent ranking, calibration, transport, and experimental-selection value of Boltz-2 ipTM'),
        ('Draft date', '2026-11-16'),
        ('Status', 'Complete scientific supplement draft for review; not yet a submission-formatted publisher file.'),
        ('Workbook structure', 'Cover, Sources, and S1–S15. Multiple panels are stacked vertically and named as Excel tables.'),
        ('Numerical policy', 'Frozen accepted outputs are transcribed without scientific reanalysis. Arithmetic count tabulations are explicitly labelled.'),
        ('Blank cells', 'A blank numeric cell can denote a source-level non-evaluable quantity. Adjacent status/failure/message fields and panel notes explain why.'),
        ('Primary evidence', 'Main-cohort RQ1/RQ2/RQ3, external confirmatory zero-label U, and prespecified sensitivity analyses.'),
        ('Exploratory evidence', 'External post hoc I/F, retrospective shared-budget simulation, Pinc increments, and paired predictor-pipeline descriptions.'),
        ('Main cohort', '2,084 design–target pairs; 748 primary positives; 678 strict positives; six targets.'),
        ('External confirmatory cohort', '720 pairs; 287 Adaptyv Weak+ positives; 254 Medium+ positives; eight targets.'),
        ('Low-event stress test', '521 pairs across six targets; descriptive and nonconfirmatory.'),
        ('Pinc complete cases', '2,074 pairs and 740 primary positives; ten missing Pinc values in PD-L1.'),
        ('Unresolved before submission', 'Authors/declarations, repository/DOI, journal sizing.'),
    ], columns=['Field', 'Value'])
    cover.to_excel(writer, sheet_name='Cover', index=False, startrow=2)
    toc_rows = []
    for sheet, info in sheets.items():
        toc_rows.append({'Sheet': sheet, 'Title': info['title'], 'Purpose': info['purpose'], 'Evidence boundary': info['boundary'], 'Panels': len(info['panels'])})
    pd.DataFrame(toc_rows).to_excel(writer, sheet_name='Cover', index=False, startrow=19)

    # Register accepted figure panels before writing Sources; copy_figure later reuses these records.
    for rel in [
        f'{RQ3EIF}/RQ3E_IF_OOF_metrics.png', f'{RQ3EIF}/RQ3E_IF_OOF_reliability.png',
        f'{LOW}/RQ3E_U_low_event_stress_test.png', f'{RQ3E}/RQ3E_U_external_sensitivity.png',
        f'{RQ5}/RQ5_raw_reliability.png',
    ]:
        path = abspath(rel)
        if not path.is_file():
            raise FileNotFoundError(path)
        used_sources[str(path)] = rel
    source_rows = []
    for abs_path, rel in sorted(used_sources.items(), key=lambda item: item[1]):
        path = Path(abs_path)
        source_rows.append({'Relative path': rel, 'Bytes': path.stat().st_size, 'SHA-256': sha256(path)})
    pd.DataFrame([{'Field': 'Source policy', 'Value': 'Local frozen sources used to construct this supplement. Proteinbase is ODC-BY with required attribution; external-release data/documentation are CC BY 4.0 and included scripts are MIT. Third-party material retains its own terms.'}]).to_excel(writer, sheet_name='Sources', index=False, startrow=1)
    pd.DataFrame(source_rows).to_excel(writer, sheet_name='Sources', index=False, startrow=5)

    for sheet, info in sheets.items():
        pd.DataFrame().to_excel(writer, sheet_name=sheet, index=False)
        ws = writer.book[sheet]
        ws['A1'] = f'{sheet}. {info["title"]}'
        ws['A2'] = f'Purpose: {info["purpose"]}'
        ws['A3'] = f'Evidence boundary: {info["boundary"]}'
        row = 5
        for idx, item in enumerate(info['panels'], start=1):
            df = item['df']
            ws.cell(row=row, column=1, value=item['name'])
            ws.cell(row=row + 1, column=1, value=item['note'])
            header_row = row + 3
            df.to_excel(writer, sheet_name=sheet, index=False, startrow=header_row - 1)
            end_row = header_row + len(df)
            end_col = max(1, len(df.columns))
            ref = f'A{header_row}:{get_column_letter(end_col)}{end_row}'
            table_name = f'{sheet}_Panel_{idx}'
            tab = Table(displayName=table_name, ref=ref)
            tab.tableStyleInfo = TableStyleInfo(name='TableStyleMedium2', showFirstColumn=False, showLastColumn=False, showRowStripes=True, showColumnStripes=False)
            ws.add_table(tab)
            item['table_name'] = table_name
            item['header_row'] = header_row
            item['end_row'] = end_row
            item['end_col'] = end_col
            row = end_row + 3

# Workbook styling and navigation.
wb = load_workbook(XLSX)
navy_fill = PatternFill('solid', fgColor='183B56')
light_fill = PatternFill('solid', fgColor='EAF0F4')
subtle_fill = PatternFill('solid', fgColor='F4F7F9')
white_font = Font(color='FFFFFF', bold=True, size=11)
title_font = Font(color='183B56', bold=True, size=15)
section_font = Font(color='183B56', bold=True, size=11)
small_grey = Font(color='536471', italic=True, size=9)
thin_grey = Side(style='thin', color='D8E1E7')

for ws in wb.worksheets:
    ws.sheet_view.showGridLines = False
    ws.freeze_panes = 'A4'
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0
    ws.sheet_properties.outlinePr.summaryBelow = True
    max_col = min(max(ws.max_column, 2), 10)
    if ws.title not in {'Cover', 'Sources'}:
        for row in (1, 2, 3):
            ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=max_col)
        ws['A1'].font = title_font
        ws['A1'].fill = light_fill
        ws['A1'].alignment = Alignment(wrap_text=True, vertical='center')
        ws.row_dimensions[1].height = 30
        for row in (2, 3):
            ws.cell(row=row, column=1).font = small_grey
            ws.cell(row=row, column=1).alignment = Alignment(wrap_text=True, vertical='top')
            ws.row_dimensions[row].height = 30
        for item in sheets[ws.title]['panels']:
            heading_row = item['header_row'] - 3
            note_row = item['header_row'] - 2
            ws.merge_cells(start_row=heading_row, start_column=1, end_row=heading_row, end_column=max_col)
            ws.merge_cells(start_row=note_row, start_column=1, end_row=note_row, end_column=max_col)
            ws.cell(row=heading_row, column=1).font = section_font
            ws.cell(row=heading_row, column=1).fill = subtle_fill
            ws.cell(row=note_row, column=1).font = small_grey
            ws.cell(row=note_row, column=1).alignment = Alignment(wrap_text=True, vertical='top')
            ws.row_dimensions[note_row].height = 28
            ws.cell(row=heading_row, column=1).alignment = Alignment(wrap_text=True, vertical='center')
            ws.row_dimensions[heading_row].height = 24
            for cell in ws[item['header_row']]:
                if cell.column <= item['end_col']:
                    cell.font = white_font
                    cell.fill = navy_fill
                    cell.alignment = Alignment(wrap_text=True, vertical='center')
            ws.row_dimensions[item['header_row']].height = 32
    for col_idx in range(1, ws.max_column + 1):
        values = []
        for row_idx in range(1, min(ws.max_row, 350) + 1):
            value = ws.cell(row=row_idx, column=col_idx).value
            if value is not None:
                values.append(str(value))
        max_len = max([len(v) for v in values] or [8])
        header_values = [str(ws.cell(row=r, column=col_idx).value or '') for r in range(1, min(ws.max_row, 50) + 1)]
        is_long = any(any(k in v.lower() for k in ['definition', 'boundary', 'message', 'path', 'source', 'note', 'purpose', 'value', 'json', 'sequence']) for v in header_values)
        width = min(55 if is_long else 32, max(11, max_len * 0.9 + 2))
        ws.column_dimensions[get_column_letter(col_idx)].width = width
    for row in ws.iter_rows():
        for cell in row:
            if cell.row > 3 and cell.value is not None and cell.alignment == Alignment():
                cell.alignment = Alignment(vertical='top')
            if isinstance(cell.value, str) and len(cell.value) > 100:
                cell.alignment = Alignment(wrap_text=True, vertical='top')

# Cover styling.
cover_ws = wb['Cover']
cover_ws['A1'] = 'Supplementary materials | Scientific review draft'
cover_ws['A1'].font = Font(size=17, bold=True, color='FFFFFF')
cover_ws['A1'].fill = navy_fill
cover_ws.merge_cells('A1:D1')
cover_ws.row_dimensions[1].height = 32
for row in (3, 20):
    for cell in cover_ws[row]:
        if cell.value is not None:
            cell.font = white_font
            cell.fill = navy_fill
            cell.alignment = Alignment(wrap_text=True)
cover_ws.column_dimensions['A'].width = 30
cover_ws.column_dimensions['B'].width = 85
cover_ws.column_dimensions['C'].width = 65
cover_ws.column_dimensions['D'].width = 70
cover_ws.column_dimensions['E'].width = 12
for row in range(4, 36):
    for cell in cover_ws[row]:
        cell.alignment = Alignment(wrap_text=True, vertical='top')
    cover_ws.row_dimensions[row].height = 34
cover_ws.freeze_panes = 'A20'

source_ws = wb['Sources']
source_ws['A1'] = 'Frozen input sources and SHA-256 hashes'
source_ws['A1'].font = title_font
source_ws.merge_cells('A1:C1')
for cell in source_ws[6]:
    if cell.value is not None:
        cell.font = white_font
        cell.fill = navy_fill
source_ws.freeze_panes = 'A7'
source_ws.auto_filter.ref = f'A6:C{source_ws.max_row}'
source_ws.column_dimensions['A'].width = 120
source_ws.column_dimensions['B'].width = 18
source_ws.column_dimensions['C'].width = 70
wb.save(XLSX)

# Reopen to verify workbook integrity and table geometry.
check_wb = load_workbook(XLSX, data_only=False, read_only=False)
assert check_wb.sheetnames == ['Cover', 'Sources'] + [f'S{i}' for i in range(1, 16)]
for rec in table_records:
    ws = check_wb[rec['sheet']]
    assert ws.max_row >= rec['rows']
    assert len(ws.tables) >= 1
assert all(check_wb[f'S{i}'].max_row > 10 for i in range(1, 16))
check_wb.close()


def finish_figure(fig, number, filename, title, caption, sources, footnote=None):
    if footnote:
        fig.text(0.01, 0.008, footnote, fontsize=7.4, color=GREY, va='bottom')
    path = FIG / filename
    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    with Image.open(path) as image:
        size = image.size
    figure_records.append({'number': number, 'filename': filename, 'title': title, 'caption': caption, 'sources': sources, 'pixels': list(size), 'bytes': path.stat().st_size})


def copy_figure(number, source_rel, filename, title, caption):
    source = abspath(source_rel)
    if not source.is_file():
        raise FileNotFoundError(source)
    used_sources[str(source)] = source_rel
    target = FIG / filename
    shutil.copy2(source, target)
    with Image.open(target) as image:
        size = image.size
    figure_records.append({'number': number, 'filename': filename, 'title': title, 'caption': caption, 'sources': [source_rel], 'pixels': list(size), 'bytes': target.stat().st_size})


# S1: outcome and cohort sensitivity.
sets = ['primary_six', 'primary_five', 'strict_six', 'strict_five']
short_set_labels = ['Primary\n6 targets', 'Primary\n5 targets\n(no AAH-I)', 'Strict\n6 targets', 'Strict\n5 targets\n(no AAH-I)']
fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4))
fig.suptitle('Supplementary Figure S1 | Outcome and cohort sensitivity', x=0.02, ha='left', fontweight='bold', color=NAVY)
summary_auc = rq12_summary[(rq12_summary.weighting == 'target_equal') & (rq12_summary.metric == 'roc_auc')].set_index('analysis_set')
x = np.arange(4)
ax = axes[0, 0]
vals = [summary_auc.loc[s, 'estimate'] for s in sets]
ax.bar(x, vals, color=[BLUE, TEAL, ORANGE, PURPLE], width=0.68)
ax.set_xticks(x, short_set_labels, fontsize=7.6)
ax.set_ylim(0, 0.82)
ax.set_ylabel('Target-equal mean AUC')
ax.set_title('A | Discrimination is sensitive to AAH-I, not to Weak inclusion', loc='left', fontweight='bold')
for xi, value in zip(x, vals):
    ax.text(xi, value + 0.015, f'{value:.3f}', ha='center', fontsize=8)
ax.grid(axis='y', alpha=0.2)

ax = axes[0, 1]
primary_target = rq12_target[(rq12_target.analysis_set == 'primary_six') & (rq12_target.metric == 'roc_auc')].set_index('target')
strict_target = rq12_target[(rq12_target.analysis_set == 'strict_six') & (rq12_target.metric == 'roc_auc')].set_index('target')
xx = np.arange(len(TARGETS))
for offset, frame, label, color in [(-0.17, primary_target, 'Primary (Weak+)', BLUE), (0.17, strict_target, 'Strict (Medium/Strong)', ORANGE)]:
    values = frame.loc[TARGETS, 'estimate'].to_numpy(float)
    low = frame.loc[TARGETS, 'ci_low'].to_numpy(float)
    high = frame.loc[TARGETS, 'ci_high'].to_numpy(float)
    ax.errorbar(xx + offset, values, yerr=[values-low, high-values], fmt='o', capsize=3, markersize=4.5, linewidth=1, label=label, color=color)
ax.axhline(0.5, color=GREY, linestyle='--', linewidth=0.8)
ax.set_xticks(xx, [TARGET_LABELS[t] for t in TARGETS], rotation=32, ha='right', fontsize=7.5)
ax.set_ylim(0.35, 1.03)
ax.set_ylabel('Within-target AUC (95% interval)')
ax.set_title('B | Target-level primary and strict outcomes', loc='left', fontweight='bold')
ax.legend(frameon=False, fontsize=7.5, loc='lower left')
ax.grid(axis='y', alpha=0.2)

for ax, metric, title, ylabel in [
    (axes[1, 0], 'brier', 'C | Raw probability error', 'Target-equal mean Brier score'),
    (axes[1, 1], 'bss_target_prevalence', 'D | Raw score versus target-constant prevalence', 'Target-equal mean BSS'),
]:
    frame = rq12_summary[(rq12_summary.weighting == 'target_equal') & (rq12_summary.metric == metric)].set_index('analysis_set')
    vals = [frame.loc[s, 'estimate'] for s in sets]
    colors = [BLUE, TEAL, ORANGE, PURPLE]
    ax.bar(x, vals, color=colors, width=0.68)
    if metric.startswith('bss'):
        ax.axhline(0, color=NAVY, linewidth=0.9)
    ax.set_xticks(x, short_set_labels, fontsize=7.6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left', fontweight='bold')
    for xi, value in zip(x, vals):
        ax.text(xi, value + (0.025 if metric == 'brier' else 0.07), f'{value:.3f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=8)
    ax.grid(axis='y', alpha=0.2)
fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.16, hspace=0.55, wspace=0.28)
finish_figure(fig, 1, 'Supplementary_Figure_S1.png', 'Outcome and cohort sensitivity of discrimination and raw probability error',
              'Target-equal AUC, Brier score, and Brier skill score under primary and strict outcomes, with and without AAH-I. Target-level intervals are 95% within-target sequence-cluster bootstrap intervals. Lower Brier is better; BSS=0 equals the target-constant prevalence reference. These sensitivity analyses do not create a new primary cohort.',
              ['analysis/rq12_outputs/rq12_target_metrics.csv', 'analysis/rq12_outputs/rq12_summary_metrics.csv'],
              'Primary includes Weak calls; strict requires Medium/Strong. Five-target analyses exclude AAH-I. Intervals describe observed-target precision, not future-target uncertainty.')

# S2: nonlinear U sensitivity.
loto_macro = read_csv(f'{RQ3}/rq3_loto_macro_summary.csv')
loto_metrics = read_csv(f'{RQ3}/rq3_loto_metrics.csv')
base = loto_macro[(loto_macro.cluster_rule == 'sequence_cluster90') & (loto_macro.source_weighting == 'target_equal') & (loto_macro.summary_type == 'target_equal_macro') & (loto_macro.model_layer_label == 'U_transport_zero_label') & (loto_macro.metric == 'bss_transport')]
assert len(base) >= 6
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9), gridspec_kw={'width_ratios': [1, 1.4]})
fig.suptitle('Supplementary Figure S2 | Nonlinear score relationships do not resolve zero-label transport', x=0.02, ha='left', fontweight='bold', color=NAVY)
ax = axes[0]
sets3 = ['primary_six', 'primary_five', 'strict_six']
xx = np.arange(3)
for offset, form, label, color in [(-0.18, 'linear', 'Linear', BLUE), (0.18, 'restricted_cubic_spline', 'Restricted cubic spline', ORANGE)]:
    vals = [base[(base.analysis_set == s) & (base.model_form == form)]['mean'].iloc[0] for s in sets3]
    ax.bar(xx + offset, vals, width=0.34, label=label, color=color)
    for xi, val in zip(xx + offset, vals):
        ax.text(xi, val + 0.004 if val < 0 else val - 0.004, f'{val:.3f}', ha='center', va='bottom' if val < 0 else 'top', fontsize=7, color='white', fontweight='bold')
ax.axhline(0, color=NAVY, linewidth=0.9)
ax.set_xticks(xx, ['Primary\n6 targets', 'Primary\n5 targets\n(no AAH-I)', 'Strict\n6 targets'])
ax.set_ylim(-0.105, 0.012)
ax.set_ylabel('BSS-transport (target-equal mean)')
ax.set_title('A | Cohort-level U sensitivity', loc='left', fontweight='bold')
ax.legend(frameon=False, fontsize=8, loc='upper right')
ax.grid(axis='y', alpha=0.2)

ax = axes[1]
required = {'analysis_set', 'cluster_rule', 'source_weighting', 'model_form', 'model_layer_label', 'held_out_target', 'bss_transport'}
missing = required - set(loto_metrics.columns)
assert not missing, missing
per = loto_metrics[(loto_metrics.analysis_set == 'primary_six') & (loto_metrics.cluster_rule == 'sequence_cluster90') & (loto_metrics.source_weighting == 'target_equal') & (loto_metrics.model_layer_label == 'U_transport_zero_label')]
for form, label, color, marker in [('linear', 'Linear', BLUE, 'o'), ('restricted_cubic_spline', 'Restricted cubic spline', ORANGE, 's')]:
    frame = per[per.model_form == form].set_index('held_out_target').loc[TARGETS]
    ax.plot(np.arange(6), frame.bss_transport.to_numpy(float), marker=marker, linewidth=1.2, markersize=5, color=color, label=label)
ax.axhline(0, color=NAVY, linewidth=0.9)
ax.set_xticks(np.arange(6), [TARGET_LABELS[t] for t in TARGETS], rotation=30, ha='right', fontsize=8)
ax.set_ylabel('BSS-transport')
ax.set_title('B | Held-out targets in the primary six-target analysis', loc='left', fontweight='bold')
ax.legend(frameon=False, fontsize=8)
ax.grid(axis='y', alpha=0.2)
fig.subplots_adjust(left=0.08, right=0.98, top=0.83, bottom=0.25, wspace=0.30)
finish_figure(fig, 2, 'Supplementary_Figure_S2.png', 'Linear versus restricted-cubic-spline zero-label transport',
              'The sole prespecified nonlinear sensitivity analysis uses a three-knot restricted cubic spline on clipped logit-ipTM, with knots learned only from source training data. Primary six-target BSS-transport improves from −0.0809 to −0.0240 but remains negative. No further model search was triggered.',
              [f'{RQ3}/rq3_loto_macro_summary.csv', f'{RQ3}/rq3_loto_metrics.csv'],
              'BSS-transport uses source-prevalence references. Positive means improvement over that reference; it does not establish target-oracle calibration.')

copy_figure(3, f'{RQ3EIF}/RQ3E_IF_OOF_metrics.png', 'Supplementary_Figure_S3.png', 'Post hoc external U/I/F out-of-fold metrics',
            'External local intercept and intercept-plus-slope adaptation, evaluated out of fold on identical rows as U. Intercept adaptation reduces Brier in seven of eight targets; F improves on I only in PD-L1 and TREM2. This analysis was added after external U was examined and is exploratory.')
copy_figure(4, f'{RQ3EIF}/RQ3E_IF_OOF_reliability.png', 'Supplementary_Figure_S4.png', 'Post hoc external U/I/F reliability',
            'Out-of-fold reliability for external U, I, and F. Local labels are used during adaptation; the display does not alter or rescue the confirmatory zero-label U result.')
copy_figure(5, f'{LOW}/RQ3E_U_low_event_stress_test.png', 'Supplementary_Figure_S5.png', 'Nonconfirmatory low-event stress test',
            'Six low-event targets (521 pairs) evaluated with the unchanged frozen U mapping. Mean predictions exceed observed prevalence in all six. MBP has zero positives; BBF-14 and 15-PGDH each have one, limiting or precluding slope and target-oracle interpretation.')
copy_figure(6, f'{RQ3E}/RQ3E_U_external_sensitivity.png', 'Supplementary_Figure_S6.png', 'External outcome, score-rule, and source-weight sensitivity',
            'Prespecified external sensitivity analyses across outcome definitions, score summaries/seeds, and source weighting. Changes in Brier or BSS-transport do not establish a stable probability mapping; BSS-target-oracle remains a diagnostic reference.')

# S7: complete shared-budget sensitivity grid.
rq4c = read_csv(f'{RQ4}/rq4c_comparisons.csv')
pilot = rq4c[rq4c.comparison == 'pilot_I_minus_raw_ipTM'].copy()
fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.5), sharey=True)
fig.suptitle('Supplementary Figure S7 | Pilot value depends on budget, pilot cost, and analysis set', x=0.02, ha='left', fontweight='bold', color=NAVY)
for ax, analysis_set in zip(axes, ['primary_six', 'primary_five', 'strict_six']):
    frame = pilot[pilot.analysis_set == analysis_set]
    for budget, color in [(120, BLUE), (300, ORANGE)]:
        sub = frame[frame.B == budget].sort_values('m_per_target')
        xx = sub.m_per_target.to_numpy(float)
        yy = sub.mean_hit_difference.to_numpy(float)
        lo = sub.difference_lo.to_numpy(float)
        hi = sub.difference_hi.to_numpy(float)
        ax.errorbar(xx, yy, yerr=[yy-lo, hi-yy], fmt='o-', capsize=3, markersize=4.5, linewidth=1.2, color=color, label=f'B={budget}')
        for xval, yval in zip(xx, yy):
            ax.annotate(f'{yval:+.2f}', (xval, yval), xytext=(0, 6), textcoords='offset points', ha='center', fontsize=7, color=color)
    ax.axhline(0, color=NAVY, linestyle='--', linewidth=0.8)
    ax.set_xticks([10, 20, 40])
    ax.set_xlabel('Pilot size m per target')
    ax.set_title(SET_LABELS[analysis_set].replace('\n', ' '), fontsize=9, fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
axes[0].set_ylabel('Mean hit difference: pilot I − raw ipTM')
axes[-1].legend(frameon=False, fontsize=8, loc='lower left')
fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.22, wspace=0.18)
finish_figure(fig, 7, 'Supplementary_Figure_S7.png', 'Shared-budget pilot sensitivity across analysis sets',
              'Paired retrospective differences between pilot-I and raw-ipTM allocation. Error bars are central 95% simulation ranges. Pilot assays count against B and pilot-identified hits count toward total hits. Primary six-target scenarios use 2,000 replicates; five-target and strict scenarios use 1,000.',
              [f'{RQ4}/rq4c_comparisons.csv'],
              'Post hoc exploratory observed-pool simulation. Ranges are not confidence intervals for future targets; no universal pilot size is implied.')

copy_figure(8, f'{RQ5}/RQ5_raw_reliability.png', 'Supplementary_Figure_S8.png', 'Raw Pinc and ipTM reliability on common complete cases',
            'Raw-score reliability for Pinc and ipTM on the same 2,074 complete-case pairs and 740 primary positives. Both scores have negative target-constant BSS in every target; lower raw error is not achieved calibration.')

# S9: Pinc increment by source weighting and analysis set.
rq5_sens = read_csv(f'{RQ5}/rq5_sensitivity_summary.csv')
inc = rq5_sens[(rq5_sens.analysis_set.isin(['primary_six', 'primary_five', 'strict_six'])) &
               (rq5_sens.outcome_name.isin(['primary_outcome', 'strict_outcome'])) &
               (rq5_sens.cluster_rule == 'sequence_cluster90') &
               (rq5_sens.comparison_family == 'preset_U_exposure') &
               (rq5_sens.comparison == 'ipTM_plus_Pinc_minus_ipTM_only') &
               (rq5_sens.layer == 'U_transport_zero_label') &
               (rq5_sens.metric == 'brier') &
               (rq5_sens.summary_type == 'target_equal_macro') &
               (rq5_sens.source_weighting.isin(['target_equal', 'pair']))].copy()
assert len(inc) == 6, inc[['analysis_set', 'source_weighting', 'mean_difference']]
heat = np.empty((3, 2))
for i, analysis_set in enumerate(['primary_six', 'primary_five', 'strict_six']):
    for j, weighting in enumerate(['target_equal', 'pair']):
        rows = inc[(inc.analysis_set == analysis_set) & (inc.source_weighting == weighting)]
        assert len(rows) == 1, (analysis_set, weighting, len(rows))
        heat[i, j] = rows.iloc[0].mean_difference
fig, ax = plt.subplots(figsize=(8.8, 5.2))
limit = max(abs(heat.min()), abs(heat.max())) * 1.1
im = ax.imshow(heat, cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect='auto')
ax.set_xticks([0, 1], ['Target-equal source weighting', 'Source-pair weighting'])
ax.set_yticks([0, 1, 2], ['Primary, six targets', 'Primary, five targets\n(excluding AAH-I)', 'Strict, six targets'])
for i in range(3):
    for j in range(2):
        ax.text(j, i, f'{heat[i, j]:+.4f}', ha='center', va='center', fontsize=12, fontweight='bold', color='white' if abs(heat[i, j]) > limit * 0.55 else NAVY)
ax.set_title('Supplementary Figure S9 | Pinc increment changes sign with source weighting', loc='left', fontweight='bold', color=NAVY, pad=16)
ax.text(0, 1.02, 'ΔBrier = ipTM + Pinc − ipTM-only under LOTO U; negative favors the combined model', transform=ax.transAxes, fontsize=9, color=GREY)
cbar = fig.colorbar(im, ax=ax, shrink=0.78, pad=0.03)
cbar.set_label('Target-equal macro ΔBrier')
for spine in ax.spines.values():
    spine.set_visible(False)
fig.subplots_adjust(left=0.26, right=0.90, top=0.79, bottom=0.16)
finish_figure(fig, 9, 'Supplementary_Figure_S9.png', 'Pinc increment across analysis sets and source weighting',
              'Complete-case LOTO U comparison of ipTM+Pinc versus ipTM-only on identical rows and folds. The summary_type filter is target_equal_macro; source_weighting varies independently between target_equal and pair. Negative ΔBrier favors the combined model. The sign reversal under source-pair weighting precludes a stable general increment claim.',
              [f'{RQ5}/rq5_sensitivity_summary.csv'],
              'Exploratory complete-case analysis; 2,074 pairs and 740 primary positives. This is not a universal comparison of Pinc and ipTM.')

# S10: full target-equal predictor summaries with observed target ranges.
metric_specs = [('auc', 'AUC ↑'), ('ap', 'AP ↑'), ('brier', 'Raw Brier ↓'), ('log_loss', 'Raw logarithmic loss ↓')]
predictors = ['boltz2', 'afm3', 'chai1', 'ef2fast', 'ef2full', 'ptxv2', 'of3', 'odde', 'rf3', 'af3of3']
outcomes = [('adaptyv_weak_plus', 'Adaptyv Weak+', BLUE), ('adaptyv_medium_plus', 'Adaptyv Medium+', ORANGE)]
fig, axes = plt.subplots(2, 2, figsize=(12, 8.3))
fig.suptitle('Supplementary Figure S10 | Ten predictor pipelines on the same 720 external designs', x=0.02, ha='left', fontweight='bold', color=NAVY)
for ax, (metric, ylabel) in zip(axes.flat, metric_specs):
    xx = np.arange(len(predictors))
    for offset, (outcome, label, color) in zip([-0.16, 0.16], outcomes):
        frame = rq6_summary[(rq6_summary.metric == metric) & (rq6_summary.outcome == outcome)].set_index('predictor').loc[predictors]
        mean = frame.target_equal_mean.to_numpy(float)
        low = frame.target_equal_min.to_numpy(float)
        high = frame.target_equal_max.to_numpy(float)
        ax.errorbar(xx + offset, mean, yerr=[mean-low, high-mean], fmt='o', color=color, markersize=4.2, linewidth=0.9, capsize=2, label=label)
    ax.set_xticks(xx, predictors, rotation=35, ha='right', fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(f'{ylabel}: target-equal mean and observed target range', loc='left', fontsize=9.5, fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    if metric == 'auc':
        ax.axhline(0.5, color=GREY, linestyle='--', linewidth=0.8)
axes[0, 0].legend(frameon=False, fontsize=8, loc='lower left')
fig.subplots_adjust(left=0.07, right=0.98, top=0.91, bottom=0.17, hspace=0.55, wspace=0.24)
finish_figure(fig, 10, 'Supplementary_Figure_S10.png', 'Target-equal description across ten predictor pipelines',
              'Points are target-equal means across eight targets; lines extend from the minimum to maximum target-specific value and are not confidence intervals. Ten complete pipelines are compared on the same 720 external designs under two outcomes. No U/I/F mapping is fitted across predictors and no superiority test is performed.',
              [f'{RQ6}/rq6_target_equal_summary.csv'],
              'Complete-pipeline description only. af3of3 uses AF3 architecture code with OpenFold3 weights. Lower raw Brier/log loss does not establish calibration.')

assert [r['number'] for r in figure_records] == list(range(1, 11))

# Build concise but complete English supplementary guide.
def md_link(text, target):
    return f'[{text}]({target})'

lines = []
lines += [
    '# Supplementary Materials',
    '',
    '## Context-dependent ranking, calibration, transport, and experimental-selection value of Boltz-2 ipTM in computational binder design',
    '',
    '**Scientific review draft — updated 18 November 2026**  ',
    '**Contents:** Supplementary Tables S1–S15 in one Excel workbook and Supplementary Figures S1–S10.  ',
    '**Status:** Complete scientific supplement draft. Journal-specific file conversion, public repository/DOI, and author metadata remain submission-stage tasks. Source versions, formal citations, and redistribution terms for the two study datasets have been verified.',
    '',
    'The supplement is organized around the evidentiary question rather than around software runs. Ranking, absolute probability error, zero-label transport, label-dependent adaptation, and experimental decision value are different claims. Full configurations are retained so that a favorable average cannot conceal unfavorable targets, alternative outcomes, source weighting, or non-evaluable settings.',
    '',
    '## Files and reading guide',
    '',
    f'- {md_link("Supplementary Tables S1–S15 (Excel workbook)", "Supplementary_Tables_S1_S15.xlsx")} contains Cover, Sources, and S1–S15 sheets. Panels are stacked vertically and named as Excel tables.',
    '- The workbook preserves source precision and configuration columns. Blank numerical cells may represent source-level non-evaluable quantities; consult the adjacent status/failure/message fields and the panel note.',
    '- Figure files are in `figures/`. Captions below distinguish confirmatory, sensitivity, descriptive, and post hoc exploratory evidence.',
    '- `supplement_source_manifest_2026-11-16.md` records source and generated-file hashes. `supplement_build_manifest.json` records panel dimensions and figure sources for independent quality control.',
    '',
    '### Core notation and interpretation',
    '',
    '| Term | Meaning | Why it matters |',
    '|---|---|---|',
    '| AUC / AP | Within-target ranking measures. AP must be read against target prevalence. | Ranking information does not establish a probability scale. |',
    '| Brier / logarithmic loss | Absolute probability-error measures; lower is better. | Values should not be compared naively across cohorts with different outcomes and prevalences. |',
    '| BSS-transport | Skill relative to the source-prevalence transport reference. | This is the zero-label transport comparator, not target-oracle calibration. |',
    '| BSS-target-oracle | Skill relative to held-out-target observed prevalence. | Diagnostic only: destination prevalence is unavailable in zero-label deployment. |',
    '| BSS-adaptation-train | Skill relative to prevalence in the corresponding adaptation-training folds. | Evaluates label-dependent adaptation without substituting a deployable zero-label reference. |',
    '| U | Frozen source mapping with no destination-target labels. | Tests transport of the released-score-plus-mapping pipeline. |',
    '| I | Local intercept update using U as an offset. | Isolates correction of baseline-level mismatch. |',
    '| F | Local intercept plus one calibration-slope update. | Tests whether changing the score–outcome slope adds stable value beyond I. |',
    '',
    '## Supplementary tables',
    '',
]
for i in range(1, 16):
    key = f'S{i}'
    info = sheets[key]
    lines += [
        f'### Table {key}. {info["title"]}',
        '',
        f'**Purpose.** {info["purpose"]}',
        '',
        f'**Evidence boundary.** {info["boundary"]}',
        '',
        '**Workbook panels.** ' + '; '.join(item['name'] for item in info['panels']) + '.',
        '',
    ]
lines += [
    '### Additional implementation notes for Table S13',
    '',
    'Pinc analyses use common complete cases (2,074 pairs and 740 primary positives), with ten PD-L1 pairs lacking Pinc left un-imputed. The frozen main table was not overwritten; Pinc was added through an extension table. Batch processing used official legacy R for throughput after small-sample agreement with the compiled official C implementation. This accepted toolchain deviation is recorded with source and script hashes. The analysis uses the symmetric Pinc score; chain-specific Pinc1/Pinc2 interpretation remains open because released PAE conversion metadata do not close chain orientation.',
    '',
    '### Additional scope note for Table S14',
    '',
    'The predictor analysis is a paired description on the same 720 external designs, not main-cohort cross-model validation. Implementations, MSA construction, recycles, and sampling can differ between pipelines. `af3of3` denotes AlphaFold 3 architecture code with OpenFold3 weights, not AlphaFold 3 weights. No cross-predictor U/I/F mapping or superiority test was performed.',
    '',
    '## Supplementary figures',
    '',
]
for rec in figure_records:
    lines += [
        f'### Supplementary Figure S{rec["number"]}. {rec["title"]}',
        '',
        f'![Supplementary Figure S{rec["number"]}. {rec["title"]}](figures/{rec["filename"]})',
        '',
        rec['caption'],
        '',
        '*Numerical or accepted panel source:* ' + '; '.join(f'`{source}`' for source in rec['sources']) + '.',
        '',
    ]
lines += [
    '## Reproducibility and release boundaries',
    '',
    '1. The numerical sources are accepted frozen outputs. Workbook construction did not rerun the scientific analyses, refit mappings, or choose configurations after inspecting performance.',
    '2. Main-cohort and external transport conditions differ in score-generation protocol as well as target/campaign. External U therefore evaluates the released score plus frozen mapping, not an isolated target effect under an identical protocol.',
    '3. The external confirmatory cohort covers primary prevalences of 21.1%–72.2%. The 521-pair low-event stress test is nonconfirmatory; the main cohort’s 84.3% high-prevalence extreme remains uncovered externally.',
    '4. The accessible Overath mirror has 3,669 records versus 3,766 reported in the article. Exact and near-sequence audits do not establish temporal, training-data, structural, or epitope independence.',
    '5. External I/F and shared-budget analyses are post hoc exploratory. They do not establish zero-label success, prospective assay savings, or a universal pilot size.',
    '6. Local source paths and hashes make this review package auditable, but they do not replace a public reproducibility deposit. Proteinbase is pinned to version `v2026-11-01` and revision `2d2825070da5ba77e0c7e1a0942146e8307e9282` under ODC-BY. The Claude protein binder design data release is pinned to v1.0 and revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`; its data and documentation use CC BY 4.0, included scripts use MIT, and third-party material retains its own terms. The public repository/DOI and author metadata must still be resolved before submission.',
    '',
]
DOC.write_text('\n'.join(lines), encoding='utf-8')

# Hash manifest is deliberately external to the workbook so its own hash is stable.
source_rows = []
for abs_path, rel in sorted(used_sources.items(), key=lambda item: item[1]):
    path = Path(abs_path)
    source_rows.append((rel, path.stat().st_size, sha256(path)))
manifest_lines = [
    '# Supplement source and output manifest', '',
    'Generated 16 November 2026. SHA-256 hashes support local traceability; they do not imply permission to redistribute source data.', '',
    '## Frozen input and accepted panel sources', '',
    '| Relative path | Bytes | SHA-256 |', '|---|---:|---|'
]
for rel, size, digest in source_rows:
    manifest_lines.append(f'| `{rel}` | {size} | `{digest}` |')
manifest_lines += ['', '## Generated supplement files', '', '| File | Bytes | SHA-256 |', '|---|---:|---|']
for path in [XLSX, DOC] + [FIG / rec['filename'] for rec in figure_records]:
    manifest_lines.append(f'| `{path.relative_to(OUT).as_posix()}` | {path.stat().st_size} | `{sha256(path)}` |')
manifest_lines += ['', '## Figure-source mapping', '', '| Figure | Source(s) |', '|---|---|']
for rec in figure_records:
    manifest_lines.append(f'| S{rec["number"]} | ' + '<br>'.join(f'`{source}`' for source in rec['sources']) + ' |')
SOURCE_MANIFEST.write_text('\n'.join(manifest_lines) + '\n', encoding='utf-8')

build = {
    'generated': '2026-11-16',
    'workbook': XLSX.name,
    'document': DOC.name,
    'expected_sheets': ['Cover', 'Sources'] + [f'S{i}' for i in range(1, 16)],
    'tables': table_records,
    'figures': figure_records,
    'source_hashes': {rel: sha256(abspath(rel)) for rel in sorted(set(used_sources.values()))},
    'output_hashes': {path.name: sha256(path) for path in [XLSX, DOC, SOURCE_MANIFEST]},
    'anchors': {
        'main_pairs': 2084, 'main_primary_positives': 748, 'main_strict_positives': 678,
        'external_pairs': 720, 'external_weak_positives': 287, 'external_medium_positives': 254,
        'low_event_pairs': 521, 'low_event_targets': 6,
        'pinc_complete_pairs': 2074, 'pinc_complete_positives': 740,
        'rq6_summary_rows': 180,
    },
}
BUILD_MANIFEST.write_text(json.dumps(build, ensure_ascii=False, indent=2), encoding='utf-8')

print(f'Workbook: {XLSX}')
print(f'Document: {DOC}')
print(f'Figures: {len(figure_records)}')
print(f'Sheets: {len(sheets)}; data panels: {len(table_records)}')
print(f'Workbook bytes: {XLSX.stat().st_size}')
for rec in figure_records:
    print(f'S{rec["number"]}: {rec["pixels"][0]}x{rec["pixels"][1]}, {rec["bytes"]} bytes, {rec["filename"]}')
print('BUILD COMPLETE')
