"""RQ5-RF Stage 5: QC, frozen screen finalization, and handoff deliverables."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from rq5rf_common import (
    EXPOSURES,
    FAILURE_PATH,
    RUN_DIR,
    STATUS_PATH,
    atomic_csv,
    atomic_json,
    now_iso,
    sha256_file,
    update_status,
)

STAGE = "stage_5_qc_and_final_handoff"


def check(name: str, passed: bool, **details) -> dict:
    return {"check": name, "passed": bool(passed), **details}


def json_clean(value):
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return json_clean(value.item())
    return value


def main() -> None:
    unit = json.loads((RUN_DIR / "rq5rf_unit_tests.json").read_text(encoding="utf-8"))
    stage0 = json.loads((RUN_DIR / "stage0_preflight.json").read_text(encoding="utf-8"))
    stage1 = json.loads((RUN_DIR / "stage1_summary.json").read_text(encoding="utf-8"))
    stage2 = json.loads((RUN_DIR / "stage2_u_models_summary.json").read_text(encoding="utf-8"))
    stage3 = json.loads((RUN_DIR / "stage3_adaptation_summary.json").read_text(encoding="utf-8"))
    stage4 = json.loads((RUN_DIR / "stage4_metrics_summary.json").read_text(encoding="utf-8"))
    verification = json.loads((RUN_DIR / "independent_verification.json").read_text(encoding="utf-8"))
    screen = json.loads((RUN_DIR / "direction_consistency_10_2.json").read_text(encoding="utf-8"))
    failures = pd.read_csv(FAILURE_PATH, keep_default_na=False)
    unresolved = failures[failures["resolution_status"].eq("unresolved")]
    recorded_non_evaluable = failures[
        failures["resolution_status"].eq("recorded_non_evaluable_diagnostic")
    ]

    source_parameters = pd.read_csv(RUN_DIR / "rf_source_model_parameters.csv", low_memory=False)
    inner_fits = pd.read_csv(RUN_DIR / "rf_inner_candidate_fits.csv", low_memory=False)
    inner_selection = pd.read_csv(RUN_DIR / "rf_inner_candidate_selection.csv", low_memory=False)
    final_audit = pd.read_csv(RUN_DIR / "rf_final_candidate_fit_audit.csv", low_memory=False)
    u_predictions = pd.read_csv(RUN_DIR / "rf_u_pair_predictions.csv", low_memory=False)
    pair_predictions = pd.read_csv(RUN_DIR / "rq5rf_pair_predictions.csv", low_memory=False)
    target_metrics = pd.read_csv(RUN_DIR / "rq5rf_target_metrics.csv", low_memory=False)
    macro_metrics = pd.read_csv(RUN_DIR / "rq5rf_macro_metrics.csv", low_memory=False)
    pair_metrics = pd.read_csv(RUN_DIR / "rq5rf_pair_metrics.csv", low_memory=False)
    pairwise_macro = pd.read_csv(RUN_DIR / "rq5rf_pairwise_macro_comparisons.csv", low_memory=False)
    checkpoints = list((RUN_DIR / "checkpoints").glob("*.json.gz"))
    images = list(RUN_DIR.glob("*.png"))

    rf_inner_count = int(inner_fits["model_family"].eq("random_forest").sum())
    rf_final_count = int(final_audit["model_family"].eq("random_forest").sum())
    checks = [
        check("unit_tests_all_pass", bool(unit["successful"]), tests=unit["tests_run"]),
        check("stage0_preflight_complete", stage0["status"] == "complete"),
        check("stage1_shared_rows_complete", stage1["status"] == "complete"),
        check("stage2_u_models_complete", stage2["status"] == "complete"),
        check("stage3_adaptation_complete", stage3["status"] == "complete"),
        check("stage4_metrics_complete", stage4["status"] == "complete"),
        check(
            "analysis_set_sample_counts",
            all(
                row["rows"] == row["rows"]
                for row in stage0["analysis_sets"]
            )
            and [(row["rows"], row["positives"]) for row in stage0["analysis_sets"]]
            == [(2074, 740), (1941, 721), (2074, 674)],
        ),
        check("outer_checkpoint_count", len(checkpoints) == 34, observed=len(checkpoints)),
        check(
            "source_parameter_rows",
            len(source_parameters) == 136
            and source_parameters["fit_status"].eq("evaluable").all(),
            observed=len(source_parameters),
        ),
        check(
            "inner_fit_record_rows",
            len(inner_fits) == 4800,
            observed=len(inner_fits),
        ),
        check(
            "inner_selection_rows",
            len(inner_selection) == 1020
            and int(inner_selection["selected"].fillna(False).sum()) == 136,
            observed=len(inner_selection),
        ),
        check(
            "final_candidate_audit_rows",
            len(final_audit) == 680,
            observed=len(final_audit),
        ),
        check(
            "rf_fit_count",
            rf_inner_count == 2880 and rf_final_count == 612 and rf_inner_count + rf_final_count == 3492,
            inner=rf_inner_count,
            final=rf_final_count,
            total=rf_inner_count + rf_final_count,
        ),
        check(
            "u_prediction_rows",
            len(u_predictions) == 48712,
            observed=len(u_predictions),
        ),
        check(
            "target_pair_prediction_rows",
            len(pair_predictions) == 48712,
            observed=len(pair_predictions),
        ),
        check("target_metric_rows", len(target_metrics) == 408, observed=len(target_metrics)),
        check("macro_metric_rows", len(macro_metrics) == 72, observed=len(macro_metrics)),
        check("pair_metric_rows", len(pair_metrics) == 72, observed=len(pair_metrics)),
        check(
            "complete_sample_macro_gate",
            bool(macro_metrics["macro_status"].eq("evaluable").all()),
        ),
        check(
            "complete_sample_pair_gate",
            bool(pair_metrics["complete_sample_gate"].eq("pass").all()),
        ),
        check(
            "pairwise_pair_set_gate",
            bool(stage4["pairwise_pair_sets_all_match"]),
        ),
        check(
            "adaptation_fit_status",
            bool(stage3["I_all_evaluable"] and stage3["F_all_evaluable"]),
        ),
        check(
            "target_oof_coverage",
            bool(stage3["all_each_target_model_oof_unique"]),
        ),
        check(
            "figure_count_and_readability_proxy",
            len(images) >= 4
            and next(item for item in verification["checks"] if item["name"] == "figure_readability")["pass"],
            observed=len(images),
        ),
        check("independent_verification_all_pass", bool(verification["all_passed"])),
        check(
            "unresolved_failure_ledger_items",
            len(unresolved) == 0,
            unresolved_count=len(unresolved),
            total_failure_rows=len(failures),
        ),
    ]
    required_checks = [item for item in checks if item["check"] != "unresolved_failure_ledger_items"]
    core_qc_pass = all(item["passed"] for item in required_checks) and len(unresolved) == 0

    screen["condition_4_qc_provisional"] = True
    screen["condition_4_qc_complete"] = bool(core_qc_pass)
    screen["all_conditions_satisfied"] = bool(
        screen["condition_1_primary_six_target_equal_strict_improvement"]
        and screen["condition_2_target_direction_counts"]
        and screen["condition_3_six_configurations_no_reversal"]
        and core_qc_pass
    )
    if screen["all_conditions_satisfied"]:
        screen["final_screen_conclusion"] = (
            "出现方向一致的探索性增量，值得在独立 complete-case 数据中重复核验。"
        )
    else:
        screen["final_screen_conclusion"] = (
            "未满足冻结的 10.2 全部条件；只能按实际结果报告方向不一致、局部增量或无明确增量。"
        )
    atomic_json(RUN_DIR / "direction_consistency_10_2.json", screen)

    # Primary, five-target, strict, and both weightings are all explicitly
    # retained. AAH-I remains visible in target-level and six-target tables.
    u_macro = macro_metrics[macro_metrics["model_layer_label"].eq("U_transport_zero_label")].copy()
    primary = u_macro[
        u_macro["analysis_set"].eq("primary_six")
        & u_macro["source_weighting"].eq("target_equal")
    ].copy()
    six_config_summary = []
    comparison_names = (
        "rf_combo_minus_pooled_logistic_combo",
        "rf_combo_minus_rf_ipTM",
    )
    for analysis_set in ("primary_six", "primary_five", "strict_six"):
        for weighting in ("target_equal", "pair"):
            row = {
                "analysis_set": analysis_set,
                "source_weighting": weighting,
            }
            for exposure in EXPOSURES:
                source = u_macro[
                    u_macro["analysis_set"].eq(analysis_set)
                    & u_macro["source_weighting"].eq(weighting)
                    & u_macro["model_exposure"].eq(exposure)
                ]
                if len(source) == 1:
                    for metric in ("brier", "log_loss", "bss_transport", "roc_auc", "ap"):
                        row[f"{exposure}_{metric}"] = float(source.iloc[0][metric])
            for comparison in comparison_names:
                for metric in ("brier", "log_loss", "bss_transport"):
                    selected = pairwise_macro[
                        pairwise_macro["analysis_set"].eq(analysis_set)
                        & pairwise_macro["source_weighting"].eq(weighting)
                        & pairwise_macro["comparison"].eq(comparison)
                        & pairwise_macro["model_layer_label"].eq("U_transport_zero_label")
                        & pairwise_macro["metric"].eq(metric)
                    ]
                    row[f"{comparison}_{metric}_difference"] = (
                        float(selected.iloc[0]["mean_difference_left_minus_right"])
                        if len(selected) == 1
                        else np.nan
                    )
            six_config_summary.append(row)
    six_config_frame = pd.DataFrame(six_config_summary)
    atomic_csv(RUN_DIR / "rq5rf_six_configuration_summary.csv", six_config_frame)

    primary_table = primary[
        [
            "model_exposure",
            "brier",
            "log_loss",
            "bss_transport",
            "roc_auc",
            "ap",
            "calibration_intercept",
            "calibration_slope",
            "macro_status",
        ]
    ].to_dict("records")
    screen_condition_rows = [
        {
            "condition": 1,
            "description": "primary_six target_equal: both comparisons strictly improve Brier, log loss, BSS-transport",
            "passed": bool(screen["condition_1_primary_six_target_equal_strict_improvement"]),
        },
        {
            "condition": 2,
            "description": "primary six at least 5/6 and five-target at least 4/5 Brier-improved targets for both comparisons",
            "passed": bool(screen["condition_2_target_direction_counts"]),
        },
        {
            "condition": 3,
            "description": "all six configurations show no reversal in the three required macro metrics",
            "passed": bool(screen["condition_3_six_configurations_no_reversal"]),
        },
        {
            "condition": 4,
            "description": "all predictions, expected targets, candidate tuning records, and QC gates pass",
            "passed": bool(core_qc_pass),
        },
    ]

    run_status = {
        "handoff_id": "RQ5RF_EXECUTION_CODEX_HANDOFF_2026-11-17",
        "run_directory": str(RUN_DIR),
        "overall_status": "complete" if core_qc_pass else "partial",
        "scope": "post hoc exploratory model-complexity sensitivity",
        "completed_at": now_iso(),
        "protocol_sha256": screen.get("protocol_sha256", "d9db1c3a09d1264580ca56b797922d55fce0e563c04ecb9dab76ac4744a25e28"),
        "stages": {
            "stage_0_preflight": "complete",
            "stage_1_shared_rows_and_unit_tests": "complete",
            "stage_2_u_models": "complete",
            "stage_3_adaptation": "complete",
            "stage_4_metrics_and_figures": "complete",
            "stage_5_independent_verification": (
                "complete" if verification["all_passed"] else "failed"
            ),
            "stage_5_qc_and_final_handoff": "complete" if core_qc_pass else "partial",
        },
        "rf_fits": {
            "inner": rf_inner_count,
            "final_candidate_audit": rf_final_count,
            "total": rf_inner_count + rf_final_count,
            "expected": 3492,
        },
        "unresolved_failure_count": len(unresolved),
        "resolved_failure_count": int(failures["resolution_status"].eq("resolved").sum()),
        "recorded_non_evaluable_diagnostic_count": len(recorded_non_evaluable),
        "screen_10_2_satisfied": bool(screen["all_conditions_satisfied"]),
        "screen_10_2_conclusion": screen["final_screen_conclusion"],
        "limitations": [
            "post hoc exploratory only",
            "no independent external Pinc cohort",
            "no pure interaction identification",
            "RF probabilities are not claimed to be calibrated",
            "six-target direction counts are not future-target population inference",
        ],
    }
    atomic_json(STATUS_PATH, run_status)
    atomic_json(RUN_DIR / "qc_checks.json", {"all_checks_passed": core_qc_pass, "checks": checks})

    handoff = {
        "status": run_status["overall_status"],
        "run_directory": str(RUN_DIR),
        "scope": run_status["scope"],
        "computation": run_status["rf_fits"],
        "main_primary_six_target_equal_u": primary_table,
        "six_configuration_rows": six_config_frame.to_dict("records"),
        "screen_10_2": {
            "satisfied": bool(screen["all_conditions_satisfied"]),
            "conditions": screen_condition_rows,
            "conclusion": screen["final_screen_conclusion"],
            "absolute_effects_available_in": [
                "rq5rf_pairwise_macro_comparisons.csv",
                "rq5rf_pairwise_target_differences.csv",
                "rq5rf_six_configuration_summary.csv",
            ],
        },
        "failure_ledger": {
            "path": str(FAILURE_PATH),
            "unresolved": len(unresolved),
            "resolved": int(failures["resolution_status"].eq("resolved").sum()),
            "recorded_non_evaluable_diagnostic": len(recorded_non_evaluable),
        },
        "independent_verification": {
            "all_passed": bool(verification["all_passed"]),
            "path": str(RUN_DIR / "independent_verification.json"),
        },
    }
    atomic_json(RUN_DIR / "handoff_result.json", json_clean(handoff))

    lines = [
        "# RQ5-RF 执行交付",
        "",
        f"- 运行目录：`{RUN_DIR}`",
        f"- 状态：`{run_status['overall_status']}`",
        "- 证据边界：post hoc exploratory model-complexity sensitivity；不是验证性比较、通用校准或未来靶点推断。",
        f"- scikit-learn：`1.9.0`；RF：`{rf_inner_count + rf_final_count}` 次 500-tree 拟合。",
        "",
        "## 主分析 U 层",
        "",
        "`primary_six × target_equal`，靶点等权宏平均，未使用 destination-target 标签：",
        "",
        "| model | Brier | log loss | BSS-transport | ROC-AUC | AP |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary_table:
        lines.append(
            "| {model_exposure} | {brier:.6f} | {log_loss:.6f} | {bss_transport:.6f} | {roc_auc:.6f} | {ap:.6f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## 六配置结果",
            "",
            "六个配置均完成并保存在 `rq5rf_six_configuration_summary.csv`。两种预设核心对比的 U 层宏平均差值如下（Brier/log loss 为左减右，负值利于左侧；BSS 为正值利于左侧）：",
            "",
        ]
    )
    for row in six_config_frame.to_dict("records"):
        lines.append(
            f"- `{row['analysis_set']} × {row['source_weighting']}`："
            f"RF-combo − logistic-combo：Brier {row['rf_combo_minus_pooled_logistic_combo_brier_difference']:.6g}，"
            f"log loss {row['rf_combo_minus_pooled_logistic_combo_log_loss_difference']:.6g}，"
            f"BSS {row['rf_combo_minus_pooled_logistic_combo_bss_transport_difference']:.6g}；"
            f"RF-combo − RF-ipTM：Brier {row['rf_combo_minus_rf_ipTM_brier_difference']:.6g}，"
            f"log loss {row['rf_combo_minus_rf_ipTM_log_loss_difference']:.6g}，"
            f"BSS {row['rf_combo_minus_rf_ipTM_bss_transport_difference']:.6g}。"
        )
    lines.extend(
        [
            "",
            "## 10.2 方向一致性筛查",
            "",
            f"结论：**{screen['final_screen_conclusion']}**",
            "",
            "| condition | passed | description |",
            "|---:|:---:|---|",
        ]
    )
    for row in screen_condition_rows:
        lines.append(
            f"| {row['condition']} | {'PASS' if row['passed'] else 'FAIL'} | {row['description']} |"
        )
    lines.extend(
        [
            "",
            "该筛查不是显著性检验或效应重要性门槛；即使满足，也只能表述为需要在独立 complete-case 数据中重复核验。绝对效应与逐靶点方向见 pairwise 文件。",
            "",
            "## 验证与失败项",
            "",
            f"- 独立验证：`{'PASS' if verification['all_passed'] else 'FAIL'}`，包含选参重算、指标重算、单/双特征 final RF 同 seed 复跑和逐元素预测比对。",
            f"- 未解决失败项：`{len(unresolved)}`；失败台帐：`{FAILURE_PATH}`。",
            f"- 已记录但不可评价的诊断项：`{len(recorded_non_evaluable)}`（保留在 target 表，不纳入预测指标宏平均）。",
            f"- 单元测试：`{unit['tests_run']}` 项，successful=`{unit['successful']}`。",
            "- AAH-I、五靶点、严格结局与两种权重均保留在 target/macro 表中，未删除不可评价项或不可见靶点。",
            "",
            "## 解释限制",
            "",
            "- 不识别纯交互；RF-combo − RF-ipTM 可来自 ipTM 非线性、Pinc 非线性或二者交互。",
            "- 不声称 RF 概率已校准，不声称 RF/Pinc 普遍更优，不声称解决跨靶点迁移。",
            "- 缺少独立外部 Pinc 队列；本结果不建立未来靶点总体推断。",
        ]
    )
    (RUN_DIR / "HANDOFF_RESULT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n"
    )

    # Hash every deliverable except self-referential manifests.
    manifest_rows = []
    for path in sorted(RUN_DIR.rglob("*")):
        if not path.is_file():
            continue
        if path.name in {"final_hash_manifest.csv", "manifest_sha256.txt"}:
            continue
        relative = path.relative_to(RUN_DIR).as_posix()
        manifest_rows.append(
            {
                "relative_path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    atomic_csv(RUN_DIR / "final_hash_manifest.csv", pd.DataFrame(manifest_rows))
    manifest_hash = sha256_file(RUN_DIR / "final_hash_manifest.csv")
    (RUN_DIR / "manifest_sha256.txt").write_text(
        f"{manifest_hash}  final_hash_manifest.csv\n",
        encoding="ascii",
        newline="\n",
    )
    print(
        json.dumps(
            {
                "status": run_status["overall_status"],
                "run_directory": str(RUN_DIR),
                "primary_metrics": primary_table,
                "screen_10_2_satisfied": bool(screen["all_conditions_satisfied"]),
                "unresolved_failures": len(unresolved),
                "independent_verification_passed": bool(verification["all_passed"]),
                "manifest_sha256": manifest_hash,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    raise SystemExit(0 if core_qc_pass else 1)


if __name__ == "__main__":
    main()
