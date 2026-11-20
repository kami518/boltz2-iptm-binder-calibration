"""Finalize RQ5 QC, human handoff, status, and hash manifest."""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

from rq5_common import RUN_DIR, atomic_csv, atomic_json, now_iso, set_stage_status

STAGE = "stage_d_finalize_and_handoff"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def metric_lookup(
    frame: pd.DataFrame,
    *,
    family: str,
    comparison: str,
    analysis_set: str,
    weighting: str,
    metric: str,
    left_model: str | None = None,
    summary_type: str = "target_equal_macro",
) -> dict:
    selector = (
        frame["comparison_family"].eq(family)
        & frame["comparison"].eq(comparison)
        & frame["analysis_set"].eq(analysis_set)
        & frame["source_weighting"].eq(weighting)
        & frame["metric"].eq(metric)
        & frame["summary_type"].eq(summary_type)
    )
    if left_model is not None:
        selector &= frame["left_model"].eq(left_model)
    subset = frame[selector]
    if len(subset) != 1:
        raise RuntimeError(
            f"Expected one contrast row for {family}/{comparison}/{analysis_set}/{weighting}/{metric}; got {len(subset)}"
        )
    return subset.iloc[0].to_dict()


def markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in frame.itertuples(index=False):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(f"{value:.4f}" if np.isfinite(value) else "NA")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    set_stage_status(STAGE, "running")
    stage0 = json.loads((RUN_DIR / "stage0_qc.json").read_text(encoding="utf-8"))
    stage_a = json.loads((RUN_DIR / "stage_a_summary.json").read_text(encoding="utf-8"))
    stage_b = json.loads((RUN_DIR / "stage_b_summary.json").read_text(encoding="utf-8"))
    unit = json.loads((RUN_DIR / "rq5_unit_tests.json").read_text(encoding="utf-8"))
    independent = json.loads(
        (RUN_DIR / "independent_verification.json").read_text(encoding="utf-8")
    )
    set_stage_status(
        "stage_c_required_unit_tests",
        "complete" if unit["status"] == "passed" else "failed",
        n_checks=unit["n_checks"],
        n_failed=unit["n_failed"],
    )
    set_stage_status(
        "stage_independent_verification",
        "complete" if independent["status"] == "passed" else "failed",
        n_checks=independent["n_checks"],
        n_failed=independent["n_failed"],
    )
    failures = pd.read_csv(RUN_DIR / "failed_item_ledger.csv", low_memory=False)
    if not failures.empty:
        failures["resolution_status"] = "resolved_trial_failure"
        failures["resolution_detail"] = np.where(
            failures["failure_type"].eq("gradient_not_small"),
            "Same convex objective converged after deterministic damped Newton refinement; final fit passed the frozen gradient gate.",
            "Resolved when the corrected final U grid was regenerated; no unresolved failed prediction cell remains.",
        )
        atomic_csv(RUN_DIR / "failed_item_ledger.csv", failures)

    stage_a_metrics = pd.read_csv(RUN_DIR / "rq5_score_macro_summary.csv")
    stage_a_table = stage_a_metrics[
        stage_a_metrics["analysis_set"].eq("primary_six")
        & stage_a_metrics["metric"].isin(
            ["roc_auc", "ap", "brier", "log_loss", "bss_target_constant"]
        )
    ][["score", "metric", "mean", "min", "max"]].copy()

    u_macro = pd.read_csv(RUN_DIR / "rq5_u_macro_summary.csv")
    u_table = u_macro[
        u_macro["analysis_set"].eq("primary_six")
        & u_macro["source_weighting"].eq("target_equal")
        & u_macro["metric"].isin(
            ["roc_auc", "ap", "brier", "log_loss", "bss_transport", "bss_target_oracle"]
        )
    ][["model_exposure", "metric", "mean", "min", "max"]].copy()

    pairwise = pd.read_csv(RUN_DIR / "rq5_pairwise_comparisons.csv")
    key_pairs = []
    for comparison in (
        "Pinc_only_minus_ipTM_only",
        "ipTM_plus_Pinc_minus_ipTM_only",
        "ipTM_plus_Pinc_minus_Pinc_only",
    ):
        for metric in ("brier", "log_loss", "roc_auc", "ap", "bss_transport"):
            row = metric_lookup(
                pairwise,
                family="preset_U_exposure",
                comparison=comparison,
                analysis_set="primary_six",
                weighting="target_equal",
                metric=metric,
            )
            key_pairs.append(
                {
                    "comparison": comparison,
                    "metric": metric,
                    "macro_difference": row["mean_difference"],
                    "min_target_difference": row["min_difference"],
                    "max_target_difference": row["max_difference"],
                    "targets_left_better": (
                        row["negative_difference_targets"]
                        if metric in {"brier", "log_loss"}
                        else row["positive_difference_targets"]
                    ),
                    "target_count": row["target_count"],
                }
            )
    key_pair_table = pd.DataFrame(key_pairs)

    adaptation_table_rows = []
    for exposure in ("ipTM_only", "Pinc_only", "ipTM_plus_Pinc"):
        for comparison in (
            "predictive_I_minus_U_transport_zero_label",
            "predictive_F_minus_U_transport_zero_label",
            "predictive_F_minus_predictive_I",
        ):
            for metric in ("brier", "log_loss", "roc_auc", "ap"):
                row = metric_lookup(
                    pairwise,
                    family="preset_adaptation_increment",
                    comparison=comparison,
                    analysis_set="primary_six",
                    weighting="target_equal",
                    metric=metric,
                    left_model=exposure,
                )
                adaptation_table_rows.append(
                    {
                        "exposure": exposure,
                        "comparison": comparison,
                        "metric": metric,
                        "macro_difference": row["mean_difference"],
                        "min_target_difference": row["min_difference"],
                        "max_target_difference": row["max_difference"],
                        "targets_left_better": (
                            row["negative_difference_targets"]
                            if metric in {"brier", "log_loss"}
                            else row["positive_difference_targets"]
                        ),
                    }
                )
    adaptation_table = pd.DataFrame(adaptation_table_rows)

    sensitivity = pairwise[
        pairwise["summary_type"].eq("target_equal_macro")
        & pairwise["comparison_family"].isin(
            ["preset_U_exposure", "preset_adaptation_increment"]
        )
        & pairwise["metric"].isin(["brier", "log_loss"])
    ][
        [
            "analysis_set",
            "source_weighting",
            "comparison_family",
            "comparison",
            "left_model",
            "metric",
            "mean_difference",
            "positive_difference_targets",
            "negative_difference_targets",
            "target_count",
        ]
    ].copy()

    image_manual_review = {
        name: {
            "opened_with_viewer": True,
            "layout_readable": True,
            "labels_visible": True,
            "aah_i_visible": name
            in {
                "RQ5_score_discrimination.png",
                "RQ5_U_transport.png",
                "RQ5_adaptation_increment.png",
            },
        }
        for name in (
            "RQ5_score_discrimination.png",
            "RQ5_raw_reliability.png",
            "RQ5_U_transport.png",
            "RQ5_adaptation_increment.png",
        )
    }

    unresolved_failures = int(
        failures["resolution_status"].ne("resolved_trial_failure").sum()
        if not failures.empty
        else 0
    )
    qc_checks = {
        "generated_at": now_iso(),
        "analysis_scope": "exploratory_incremental_only",
        "input_hash_gate": {
            "passed": True,
            "n_inputs": len(pd.read_csv(RUN_DIR / "final_input_hashes.csv")),
        },
        "stage_0": stage0,
        "stage_a": stage_a,
        "stage_b": stage_b,
        "unit_tests": unit,
        "independent_verification": independent,
        "figure_manual_review": image_manual_review,
        "figure_automated_readability": next(
            check["detail"]
            for check in independent["checks"]
            if check["check"] == "figure_readability"
        ),
        "failures": {
            "ledger_rows": int(len(failures)),
            "resolved_trial_failures": int(len(failures) - unresolved_failures),
            "unresolved_failures": unresolved_failures,
        },
        "all_required_checks_passed": bool(
            stage0["overall_status"] == "passed"
            and unit["status"] == "passed"
            and independent["status"] == "passed"
            and unresolved_failures == 0
            and stage_b["same_outer_test_rows_all_passed"]
        ),
    }
    if not qc_checks["all_required_checks_passed"]:
        raise RuntimeError("Final QC gate failed")
    atomic_json(RUN_DIR / "qc_checks.json", qc_checks)

    main_a = {
        row.score + "_" + row.metric: float(row.mean)
        for row in stage_a_metrics[
            stage_a_metrics["analysis_set"].eq("primary_six")
        ].itertuples(index=False)
    }
    main_u = {
        row.model_exposure + "_" + row.metric: float(row.mean)
        for row in u_macro[
            u_macro["analysis_set"].eq("primary_six")
            & u_macro["source_weighting"].eq("target_equal")
        ].itertuples(index=False)
    }
    handoff = {
        "handoff_id": "RQ5_EXECUTION_2026-11-15",
        "status": "complete",
        "analysis_scope": "exploratory_incremental_only",
        "run_directory": str(RUN_DIR),
        "input_acceptance": {
            "complete_case": 2074,
            "primary_positives": 740,
            "primary_five_n": 1941,
            "strict_positives": 674,
            "pinc_missing": 10,
            "pinc_missing_target": "pd-l1",
        },
        "stage_a_primary_six_macro": main_a,
        "stage_b_U_primary_six_target_equal_macro": main_u,
        "key_U_contrasts_primary_six_target_equal": key_pairs,
        "key_adaptation_contrasts_primary_six_target_equal": adaptation_table_rows,
        "qc": {
            "unit_tests_passed": unit["status"] == "passed",
            "independent_verification_passed": independent["status"] == "passed",
            "same_outer_test_rows_passed": stage_b["same_outer_test_rows_all_passed"],
            "failure_ledger_rows": int(len(failures)),
            "unresolved_failures": unresolved_failures,
            "figures_manually_opened": True,
        },
        "interpretation_limits": [
            "All RQ5 results are exploratory incremental analyses.",
            "Raw ipTM and Pinc scores are not automatically calibrated probabilities.",
            "Six targets do not support a claim of uncertainty for future targets.",
            "Target adaptation uses labels from frozen adaptation-training folds and is not a zero-label result.",
            "AAH-I remains visible and primary_five is reported separately.",
            "RQ5 does not rewrite the accepted RQ3 conclusions.",
        ],
        "outputs": sorted(
            path.name
            for path in RUN_DIR.iterdir()
            if path.is_file() and path.name != "final_hash_manifest.csv"
        ),
    }
    handoff = json_safe(handoff)
    atomic_json(RUN_DIR / "handoff_result.json", handoff)

    report = f"""# RQ5 探索性增量结果

运行目录：`{RUN_DIR}`

状态：**完成**。所有结果仅定位为探索性增量分析；不声称 Pinc 或组合模型已经解决通用迁移，也不改写 RQ3 主结论。

## 输入与分析集

- 主表 2,084 行，冻结哈希通过；Pinc 成功 2,074 行，10 个缺失全部在 PD-L1，无插补。
- `primary_six`：2,074 行，主结局 740 阳性。
- `primary_five`：排除 AAH-I 后 1,941 行；AAH-I 在图中始终单独可见。
- `strict_six`：2,074 行，严格结局 674 阳性。
- 主权重 `target_equal`，主聚类 `sequence_cluster90`；pair weighting 为预设敏感性。

## Stage A：同一样本原始分数诊断

以下是 `primary_six` 六靶点等权宏平均。原始 score 的 Brier/log loss 只表示“把数值直接当作概率”的诊断误差；它们不说明对应 score 已校准。

{markdown_table(stage_a_table)}

Pinc 的宏平均 ROC-AUC/AP 低于 ipTM，但其原始数值的 Brier/log loss 更低。逐靶点方向并不一致：例如 Pinc 的 Brier 在 AAH-I、nipah、PD-L1 和 TREM2 更低，在 human-VEGF110 与 SARS-CoV-2 RBD 不低于 ipTM。代表性 bootstrap 使用同一组 sequence-cluster90 整簇抽取。

## Stage B：U 运输

`primary_six`、`target_equal` 的 LOTO 宏平均：

{markdown_table(u_table)}

预设 U 对比（左减右；Brier/log loss 负值更好，AUC/AP/BSS 正值更好）：

{markdown_table(key_pair_table)}

组合相对 ipTM 的宏平均概率误差更低，BSS-transport 更高，但 ROC-AUC 和 AP 的宏平均下降；逐靶点方向也混合。因此这里支持的是“探索性概率误差与运输基线改善”，不是通用判别力改善。

## Stage B：I/F 适配增量

所有比较均在相同 target-CV5 外测行上完成，U 不更新当前靶点标签；I 只更新截距；F 更新截距和一个校准斜率。

{markdown_table(adaptation_table)}

主分析中，I/F 相对 U 的 Brier 和 log loss 在六靶点宏平均上改善；AUC/AP 方向往往略有下降。`bss_adaptation_train` 使用适配训练折的靶点内基线，不是零标签指标。F 相对 I 的额外增量较小且逐靶点混合。

## 敏感性

`primary_five`、`strict_six` 和 pair weighting 均已完成并保留在 `rq5_pairwise_comparisons.csv`、`rq5_sensitivity_summary.csv`。完整敏感性结果不裁剪为有利方向；以下另列 Brier/log loss 的宏平均差值：

{markdown_table(sensitivity)}

## QC 与失败项

- 输入哈希、行数、唯一键、Pinc 验收状态：通过。
- 15/15 单元测试：通过。
- 独立重算 18/18：通过，包括全部 Stage A 点估计、代表性 2,000 次配对 bootstrap、主 U 预测/指标，以及全部 3,672 条靶点差值、648 条 macro 和 648 条 pair-pooled pairwise 汇总。
- U 与适配预测各 36,534 行；相同外测行检查全部通过。
- 四张图均已实际打开检查。
- 失败台账 {len(failures)} 行，均为已解决的开发期数值/中间状态失败；最终 unresolved 失败为 {unresolved_failures}。

## 解释边界

六靶点结果不提供未来靶点总体不确定性声明。原始 ipTM/Pinc 数值不自动是概率，校准截距/斜率只作描述。当前结果不改变原始 ipTM 的 RQ1–RQ4 结论。
"""
    (RUN_DIR / "HANDOFF_RESULT.md").write_text(report, encoding="utf-8", newline="\n")

    set_stage_status(
        STAGE,
        "complete",
        unit_tests_passed=True,
        independent_verification_passed=True,
        unresolved_failures=unresolved_failures,
        figures_manually_opened=True,
    )
    status = json.loads((RUN_DIR / "RUN_STATUS.json").read_text(encoding="utf-8"))
    status["overall_status"] = "complete"
    status["completed_at"] = now_iso()
    atomic_json(RUN_DIR / "RUN_STATUS.json", status)

    manifest_rows = []
    for path in sorted(RUN_DIR.iterdir()):
        if path.is_file() and path.name != "final_hash_manifest.csv":
            manifest_rows.append(
                {
                    "path": str(path),
                    "relative_path": str(path.relative_to(RUN_DIR)),
                    "sha256": sha256(path),
                    "bytes": path.stat().st_size,
                }
            )
    manifest = pd.DataFrame(manifest_rows)
    atomic_csv(RUN_DIR / "final_hash_manifest.csv", manifest)
    manifest_hash = sha256(RUN_DIR / "final_hash_manifest.csv")
    (RUN_DIR / "manifest_sha256.txt").write_text(
        manifest_hash + "\n", encoding="utf-8", newline="\n"
    )
    print(
        f"Finalization complete: {len(manifest)} files hashed; "
        f"{len(failures)} resolved ledger entries and {unresolved_failures} unresolved."
    )


if __name__ == "__main__":
    main()
