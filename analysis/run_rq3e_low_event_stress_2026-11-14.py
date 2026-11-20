"""RQ3-E low-event descriptive stress test only.

Fits deployment U on the frozen six-target Proteinbase primary cohort, then
applies it unchanged to the six prespecified low-event Anthropic targets.
Does not evaluate the 720-pair confirmatory cohort and does not fit I/F.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import shutil
import sys
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.special import expit

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
PROTOCOL = ROOT / "analysis/RQ3E_external_transport_frozen_protocol_2026-11-14.md"
COMMON = ROOT / "analysis/codex_outputs/RQ3_execution_2026-11-13/run_20261113_175022/rq3_common.py"
PREP = ROOT / "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/phase_ab.py"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
EXTERNAL = ROOT / "data/external_target_audit_2026-11-13/anthropic_audited_pair_table.csv"
OUTPUT_ROOT = ROOT / "analysis/rq3e_outputs"
CLIP = 1e-6
EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]
STRESS_ORDER = ["TNFa", "Cas9", "EGFR", "BBF-14", "15-PGDH", "MBP"]
DISPLAY = {"TNFa": "TNFα", **{name: name for name in STRESS_ORDER if name != "TNFa"}}
EXPECTED_STRESS = {
    "TNFa": (150, 12),
    "Cas9": (82, 10),
    "EGFR": (90, 10),
    "BBF-14": (79, 1),
    "15-PGDH": (30, 1),
    "MBP": (90, 0),
}
EXPECTED_HASHES = {
    MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
    ROOT / "analysis/codex_tasks/RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13/rq3_execution_config.json": "a416ecd357eabb4c4d95701136c412cf23651b730bf39feab7524fbc07eeef51",
    COMMON: "276f46ad6d068ecead50a9236d1dcf6953428e415eef79da2424ad850c471f3e",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_design_summary.parquet": "3e796869b8dc0c90d7ab0e12daf25308c0cf272d04c6a5e27281f903ecd18744",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_wetlab_summary.parquet": "f80c3cfce67cb7da53576a020d476c5b74efed44e61c6fceaec9ee7de5722f88",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_cofold_predictions.parquet": "ba48010610daa9b9d184aa1ddb419beebdf728ccf745f5a8f3b48e4041e1a87d",
    EXTERNAL: "1f2a4ceab9de06cb75f4aae9a889fc5164e144232101417653944d03ee33d0c6",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_target_qualification.csv": "bb3ac9e2a0f3bc253d7b500dcc4909f1aa2291026e56d2960f028f8911e50579",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_audit_summary.json": "df377fc5b755e41a5c2585b14d87f5bf73325b0c8b2d8f94c17381b59c066b81",
    ROOT / "data/external_target_audit_2026-11-13/anthropic_sequence_independence_summary.json": "8549f59b535357dbe0003a9aa1fa7ffcdadee423bed6de3d9dcf18d0beaf1316",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def finite_or_nan(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def fmt(value: object, digits: int = 4) -> str:
    number = finite_or_nan(value)
    return "NA" if not math.isfinite(number) else f"{number:.{digits}f}"


def make_run_dir() -> Path:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    candidate = OUTPUT_ROOT / stamp
    suffix = 1
    while candidate.exists():
        candidate = OUTPUT_ROOT / f"{stamp}_{suffix}"
        suffix += 1
    candidate.mkdir()
    return candidate


def main() -> None:
    run_dir = make_run_dir()
    log_path = run_dir / "execution.log"

    def log(message: str) -> None:
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    atomic_json(
        run_dir / "RUN_STATUS.json",
        {
            "status": "running",
            "scope": "low_event_stress_test_only",
            "confirmatory_720_evaluated": False,
            "run_dir": str(run_dir),
        },
    )

    try:
        log("Checking frozen input hashes.")
        actual_hashes = {str(path): sha256(path) for path in EXPECTED_HASHES}
        mismatches = {
            str(path): {"expected": expected, "actual": actual_hashes[str(path)]}
            for path, expected in EXPECTED_HASHES.items()
            if actual_hashes[str(path)] != expected
        }
        if mismatches:
            raise RuntimeError(f"Frozen input hash mismatch: {mismatches}")

        spec = importlib.util.spec_from_file_location("rq3_common_frozen", COMMON)
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot import frozen rq3_common.py")
        common = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(common)

        log("Loading and checking the six-target source cohort.")
        source = pd.read_csv(MASTER, low_memory=False)
        required = {"pair_id", "target", "primary_outcome", "boltz2_iptm"}
        if not required.issubset(source.columns):
            raise RuntimeError(f"Missing source columns: {required - set(source.columns)}")
        if len(source) != 2084 or source.pair_id.duplicated().any():
            raise RuntimeError("Source cohort row count or pair_id uniqueness failed")
        if source.target.nunique() != 6:
            raise RuntimeError("Source cohort does not contain exactly six targets")
        y_source = pd.to_numeric(source.primary_outcome, errors="raise").to_numpy(dtype=float)
        s_source = pd.to_numeric(source.boltz2_iptm, errors="raise").to_numpy(dtype=float)
        if not np.isin(y_source, [0.0, 1.0]).all():
            raise RuntimeError("Source outcome is not binary")
        if not np.isfinite(s_source).all():
            raise RuntimeError("Source ipTM is not finite")
        x_source = common.logit_clip(s_source, CLIP)
        source_targets = source.target.astype(str).to_numpy()
        counts = source.groupby("target", sort=True).size().astype(float)
        raw_map = 1.0 / counts
        raw_weights = source.target.astype(str).map(raw_map).to_numpy(dtype=float)
        weights = raw_weights / raw_weights.mean()
        raw_sums = pd.Series(raw_weights).groupby(source.target.astype(str)).sum()
        if not np.allclose(raw_sums.to_numpy(), 1.0, atol=1e-12, rtol=0.0):
            raise RuntimeError("Target-equal raw weights do not sum to one per target")
        if not math.isclose(float(weights.mean()), 1.0, abs_tol=1e-12):
            raise RuntimeError("Normalized target-equal weights do not have mean one")

        log("Fitting deployment U using source rows only.")
        fit = common.fit_source_linear(
            x_source,
            y_source,
            source_targets,
            weights,
            "RQ3E_deployment_U_all_six_target_equal",
        )
        if not fit.get("estimable", False):
            raise RuntimeError(f"Deployment U failed QC: {fit.get('failure_type')} / {fit.get('message')}")
        levels = list(fit["levels"])
        if len(levels) != 6:
            raise RuntimeError("Deployment U does not have six source intercepts")
        theta = np.asarray(fit["theta"], dtype=float)
        covariance = np.asarray(fit["covariance"], dtype=float)
        alphas = theta[:6]
        beta = float(theta[6])
        transport_intercept = float(alphas.mean())
        contrast = np.r_[np.full(6, 1.0 / 6.0), 0.0]
        transport_se = math.sqrt(max(float(contrast @ covariance @ contrast), 0.0))
        source_prevalence = source.groupby("target", sort=True).primary_outcome.mean()
        transport_reference = float(source_prevalence.loc[levels].mean())

        parameter_rows = []
        for index, target in enumerate(levels):
            parameter_rows.append(
                {
                    "parameter_name": f"alpha[{target}]",
                    "source_target": target,
                    "estimate": float(alphas[index]),
                    "standard_error": math.sqrt(max(float(covariance[index, index]), 0.0)),
                    "n_source": int(counts.loc[target]),
                    "positive_source": int(source.loc[source.target.eq(target), "primary_outcome"].sum()),
                    "source_prevalence": float(source_prevalence.loc[target]),
                    "fit_status": fit["status"],
                }
            )
        parameter_rows.extend(
            [
                {
                    "parameter_name": "common_logit_iptm_slope",
                    "source_target": "",
                    "estimate": beta,
                    "standard_error": math.sqrt(max(float(covariance[6, 6]), 0.0)),
                    "n_source": len(source),
                    "positive_source": int(y_source.sum()),
                    "source_prevalence": float(y_source.mean()),
                    "fit_status": fit["status"],
                },
                {
                    "parameter_name": "transport_intercept_equal_mean",
                    "source_target": "",
                    "estimate": transport_intercept,
                    "standard_error": transport_se,
                    "n_source": len(source),
                    "positive_source": int(y_source.sum()),
                    "source_prevalence": transport_reference,
                    "fit_status": fit["status"],
                },
            ]
        )
        atomic_csv(run_dir / "rq3e_source_model_parameters.csv", pd.DataFrame(parameter_rows))
        atomic_json(
            run_dir / "rq3e_source_model.json",
            {
                "model": "unpenalized weighted logistic; six target intercepts and one common logit-ipTM slope",
                "levels": levels,
                "alpha": [float(value) for value in alphas],
                "beta_common": beta,
                "transport_intercept": transport_intercept,
                "transport_intercept_se": transport_se,
                "transport_reference": transport_reference,
                "probability_clip": CLIP,
                "fit_status": fit["status"],
                "gradient_max_abs": float(fit["gradient_max_abs"]),
                "objective": float(fit["objective"]),
                "iterations": int(fit["iterations"]),
                "information_qc_message": fit["information_qc_message"],
                "external_labels_used": False,
                "source_rows": len(source),
            },
        )
        log(f"Deployment U frozen: intercept={transport_intercept:.8f}, beta={beta:.8f}, reference={transport_reference:.8f}.")

        log("Selecting only the six prespecified low-event external targets.")
        external = pd.read_csv(EXTERNAL, low_memory=False)
        valid = bool_series(external.adaptyv_valid)
        score = pd.to_numeric(external.boltz2_iptm_ipsae_best_seed, errors="coerce")
        stress = external.loc[external.target.isin(STRESS_ORDER) & valid & score.notna()].copy()
        stress["score"] = score.loc[stress.index].astype(float)
        stress["outcome_value"] = bool_series(stress.adaptyv_main_positive).astype(int)
        if len(stress) != 521 or stress.uuid.duplicated().any():
            raise RuntimeError("Stress cohort row count or uuid uniqueness failed")
        if set(stress.target.astype(str)) != set(STRESS_ORDER):
            raise RuntimeError("Stress cohort target set mismatch")
        for target, (expected_n, expected_positive) in EXPECTED_STRESS.items():
            group = stress.loc[stress.target.eq(target)]
            if len(group) != expected_n or int(group.outcome_value.sum()) != expected_positive:
                raise RuntimeError(f"Stress target count mismatch: {target}")
        if int(stress.outcome_value.sum()) != 34:
            raise RuntimeError("Stress cohort total positives mismatch")
        if not np.isfinite(stress.score.to_numpy(dtype=float)).all():
            raise RuntimeError("Nonfinite external stress scores")

        x_external = common.logit_clip(stress.score.to_numpy(dtype=float), CLIP)
        eta = transport_intercept + beta * x_external
        probability = expit(eta)
        if not np.isfinite(probability).all() or (probability <= 0).any() or (probability >= 1).any():
            raise RuntimeError("External U predictions are nonfinite or outside (0,1)")
        stress["x"] = x_external
        stress["eta_U"] = eta
        stress["predicted_probability"] = probability
        stress["bss_transport_reference"] = transport_reference
        stress["bss_target_oracle_reference"] = stress.groupby("target").outcome_value.transform("mean")

        prediction_columns = [
            "uuid", "full_name", "target", "design_model", "campaign", "generator",
            "score", "x", "eta_U", "predicted_probability", "bss_transport_reference",
            "bss_target_oracle_reference",
        ]
        frozen_predictions = stress[prediction_columns].sort_values(["target", "uuid"])
        atomic_csv(run_dir / "rq3e_low_event_stress_test_predictions_frozen.csv", frozen_predictions)
        log("External predictions saved and frozen before metric calculation.")

        pair_columns = prediction_columns + ["outcome_value", "adaptyv_binding_strength", "adaptyv_expression"]
        atomic_csv(
            run_dir / "rq3e_low_event_stress_test_pairs.csv",
            stress[pair_columns].sort_values(["target", "uuid"]),
        )

        metric_rows = []
        reliability_frames = []
        for target in STRESS_ORDER:
            group = stress.loc[stress.target.eq(target)].copy()
            y = group.outcome_value.to_numpy(dtype=float)
            p = group.predicted_probability.to_numpy(dtype=float)
            transport_ref = np.full(len(group), transport_reference, dtype=float)
            oracle_ref = np.full(len(group), float(y.mean()), dtype=float)
            metrics = common.compute_performance_metrics(
                y,
                p,
                clip=CLIP,
                transport_reference=transport_ref,
                target_oracle_reference=oracle_ref,
            )
            prevalence_low, prevalence_high = common.wilson_interval(int(y.sum()), len(y))
            metrics.update(
                {
                    "target": target,
                    "target_display": DISPLAY[target],
                    "analysis_layer": "nonconfirmatory_low_event_stress_test",
                    "mean_prediction": float(p.mean()),
                    "median_prediction": float(np.median(p)),
                    "prevalence_wilson_95_low": prevalence_low,
                    "prevalence_wilson_95_high": prevalence_high,
                    "transport_reference_value": transport_reference,
                    "target_oracle_reference_value": float(y.mean()),
                    "brier_transport_reference": common.brier_score(y, transport_ref),
                    "brier_target_oracle_reference": common.brier_score(y, oracle_ref),
                }
            )
            if metrics["calibration_status"] == "evaluable":
                metrics["calibration_intercept_95_low"] = metrics["calibration_intercept"] - 1.959963984540054 * metrics["calibration_intercept_se"]
                metrics["calibration_intercept_95_high"] = metrics["calibration_intercept"] + 1.959963984540054 * metrics["calibration_intercept_se"]
                metrics["calibration_slope_95_low"] = metrics["calibration_slope"] - 1.959963984540054 * metrics["calibration_slope_se"]
                metrics["calibration_slope_95_high"] = metrics["calibration_slope"] + 1.959963984540054 * metrics["calibration_slope_se"]
            else:
                for name in ["calibration_intercept_95_low", "calibration_intercept_95_high", "calibration_slope_95_low", "calibration_slope_95_high"]:
                    metrics[name] = float("nan")
            if target == "MBP":
                metrics["interpretation_warning"] = "zero_positive_no_calibration_slope_oracle_BSS"
            elif target in {"BBF-14", "15-PGDH"}:
                metrics["interpretation_warning"] = "single_positive_descriptive_only_no_slope_inference"
            elif target == "Cas9":
                metrics["interpretation_warning"] = "low_event_RNP_apo_assay_complexity"
            else:
                metrics["interpretation_warning"] = "low_event_descriptive_only"
            metric_rows.append(metrics)
            reliability = pd.DataFrame(common.reliability_rows(group, "predicted_probability", EDGES))
            reliability.insert(0, "target", target)
            reliability.insert(1, "analysis_layer", "nonconfirmatory_low_event_stress_test")
            reliability_frames.append(reliability)

        metrics_frame = pd.DataFrame(metric_rows)
        metric_columns = [
            "target", "target_display", "analysis_layer", "n", "positives", "prevalence",
            "prevalence_wilson_95_low", "prevalence_wilson_95_high", "mean_prediction",
            "median_prediction", "brier", "log_loss", "bss_transport", "bss_target_oracle",
            "transport_reference_value", "target_oracle_reference_value",
            "brier_transport_reference", "brier_target_oracle_reference",
            "calibration_intercept", "calibration_intercept_se", "calibration_intercept_95_low",
            "calibration_intercept_95_high", "calibration_slope", "calibration_slope_se",
            "calibration_slope_95_low", "calibration_slope_95_high", "calibration_status",
            "calibration_failure_type", "calibration_message", "interpretation_warning",
        ]
        metrics_frame = metrics_frame.loc[:, metric_columns]
        reliability_frame = pd.concat(reliability_frames, ignore_index=True)
        atomic_csv(run_dir / "rq3e_low_event_stress_test_metrics.csv", metrics_frame)
        atomic_csv(run_dir / "rq3e_low_event_stress_test_reliability.csv", reliability_frame)
        if len(metrics_frame) != 6 or len(reliability_frame) != 42:
            raise RuntimeError("Metric or reliability output row count mismatch")

        log("Drawing the frozen-bin reliability figure.")
        fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True, constrained_layout=True)
        for axis, target in zip(axes.ravel(), STRESS_ORDER):
            rel = reliability_frame.loc[reliability_frame.target.eq(target) & reliability_frame.n.gt(0)]
            met = metrics_frame.loc[metrics_frame.target.eq(target)].iloc[0]
            axis.plot([0, 1], [0, 1], linestyle="--", color="0.65", linewidth=1)
            lower_error = np.maximum((rel.observed_rate - rel.wilson_95_low).to_numpy(dtype=float), 0.0)
            upper_error = np.maximum((rel.wilson_95_high - rel.observed_rate).to_numpy(dtype=float), 0.0)
            axis.errorbar(
                rel.mean_prediction,
                rel.observed_rate,
                yerr=[lower_error, upper_error],
                fmt="o-",
                markersize=4,
                linewidth=1.2,
                capsize=2,
                color="#2166ac",
            )
            axis.set_xlim(0, 1)
            axis.set_ylim(0, 1)
            axis.set_title(
                f"{DISPLAY[target]}: n={int(met.n)}, pos={int(met.positives)}\n"
                f"Brier={met.brier:.3f}; BSS-transport={met.bss_transport:.3f}",
                fontsize=10,
            )
            axis.grid(alpha=0.15)
        for axis in axes[-1, :]:
            axis.set_xlabel("Mean predicted probability")
        for axis in axes[:, 0]:
            axis.set_ylabel("Observed positive rate")
        fig.suptitle("RQ3-E low-event stress test: fixed U, frozen reliability bins", fontsize=14)
        fig.savefig(run_dir / "RQ3E_U_low_event_stress_test.png", dpi=220)
        plt.close(fig)

        master_hash_after = sha256(MASTER)
        if master_hash_after != EXPECTED_HASHES[MASTER]:
            raise RuntimeError("Frozen source master changed during execution")
        qc = {
            "status": "passed",
            "scope": "low_event_stress_test_only",
            "confirmatory_720_evaluated": False,
            "input_hash_gate": "passed",
            "source_rows": len(source),
            "source_targets": len(levels),
            "source_weight_raw_sums": {str(key): float(value) for key, value in raw_sums.items()},
            "source_weight_normalized_mean": float(weights.mean()),
            "source_model_status": fit["status"],
            "source_gradient_max_abs": float(fit["gradient_max_abs"]),
            "stress_rows": len(stress),
            "stress_unique_uuids": int(stress.uuid.nunique()),
            "stress_positives": int(stress.outcome_value.sum()),
            "stress_target_counts": {
                target: {
                    "n": int((stress.target == target).sum()),
                    "positives": int(stress.loc[stress.target == target, "outcome_value"].sum()),
                }
                for target in STRESS_ORDER
            },
            "all_predictions_finite_and_inside_open_unit_interval": True,
            "metric_rows": len(metrics_frame),
            "reliability_rows": len(reliability_frame),
            "master_hash_before": EXPECTED_HASHES[MASTER],
            "master_hash_after": master_hash_after,
            "I_or_F_fitted": False,
            "result_driven_changes": False,
        }
        atomic_json(run_dir / "qc_checks.json", qc)
        atomic_csv(run_dir / "failed_items.csv", pd.DataFrame(columns=["stage", "item", "failure_type", "message"]))

        config_snapshot = {
            "protocol_path": str(PROTOCOL),
            "protocol_sha256": sha256(PROTOCOL),
            "execution_script_sha256": sha256(Path(__file__)),
            "frozen_common_sha256": sha256(COMMON),
            "original_weight_rule_reference": str(PREP),
            "original_weight_rule_sha256": sha256(PREP),
            "scope": "six-target deployment U plus six-target / 521-pair low-event stress test only",
            "confirmatory_720_evaluated": False,
            "source": {
                "path": str(MASTER),
                "rows": 2084,
                "outcome": "primary_outcome",
                "score": "boltz2_iptm",
                "weight_rule": "raw=1/n_source_target; normalized=raw/mean(raw over 2084 source rows)",
            },
            "external": {
                "path": str(EXTERNAL),
                "targets": STRESS_ORDER,
                "rows": 521,
                "outcome": "Adaptyv Weak+ / adaptyv_main_positive",
                "score": "boltz2_iptm_ipsae_best_seed",
                "labels_used_for_training_or_tuning": False,
            },
            "formula": "eta_U = mean(six source alpha_j) + beta_common * logit(clip(ipTM,1e-6,1-1e-6))",
            "references": {
                "transport": "unweighted mean of six source-target prevalences",
                "target_oracle": "external target prevalence; diagnostic only",
            },
            "reliability_edges": EDGES,
            "probability_clip": CLIP,
            "input_hashes": actual_hashes,
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__,
            },
        }
        atomic_json(run_dir / "run_config.snapshot.json", config_snapshot)

        result_lines = [
            "# RQ3-E 低事件描述性压力测试结果",
            "",
            "**范围：仅原六靶点部署 U＋六靶点／521 行低事件压力测试；未评价 720 行确认性主队列，未拟合 I／F。**",
            "",
            f"部署 U 的运输截距为 **{transport_intercept:.6f}**，共享 logit-ipTM 斜率为 **{beta:.6f}**，源靶点等权基率参照为 **{transport_reference:.6f}**。源模型数值 QC 通过，外部标签未参与训练或调参。",
            "",
            "| 靶点 | n／阳性 | 基率（95% Wilson） | 平均预测 | Brier | Log loss | BSS-transport | BSS-target-oracle | 校准截距 | 校准斜率 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in metric_rows:
            result_lines.append(
                f"| {DISPLAY[row['target']]} | {int(row['n'])}／{int(row['positives'])} | "
                f"{fmt(row['prevalence'])}（{fmt(row['prevalence_wilson_95_low'])}–{fmt(row['prevalence_wilson_95_high'])}） | "
                f"{fmt(row['mean_prediction'])} | {fmt(row['brier'])} | {fmt(row['log_loss'])} | "
                f"{fmt(row['bss_transport'])} | {fmt(row['bss_target_oracle'])} | "
                f"{fmt(row['calibration_intercept'])} | {fmt(row['calibration_slope'])} |"
            )
        result_lines.extend(
            [
                "",
                "## 强制解释边界",
                "",
                "- 本表是非确认性低事件压力测试，不与 720 行主队列合并，不设置通过门槛。",
                "- MBP 零阳性：校准截距／斜率与 BSS-target-oracle 不可评价。",
                "- BBF-14、15-PGDH 各仅一个阳性：即使校准模型数值收敛，也不作斜率推断。",
                "- 压力测试不能证明极低基率校准可靠，也没有覆盖原 RQ3 的高基率端。",
                "- 新旧队列 score-generation protocol 不完全相同，结论是已发布 score＋冻结 U 的运输表现。",
                "",
                "逐行预测、固定分箱可靠性、QC、运行配置与哈希清单见同目录文件。",
            ]
        )
        (run_dir / "RQ3E_low_event_stress_test_results.md").write_text("\n".join(result_lines) + "\n", encoding="utf-8")

        shutil.copy2(__file__, run_dir / "rq3e_low_event_stress_test.py")
        shutil.copy2(COMMON, run_dir / "rq3_common_snapshot.py")
        shutil.copy2(PREP, run_dir / "phase_ab_weight_rule_snapshot.py")
        atomic_json(
            run_dir / "RUN_STATUS.json",
            {
                "status": "complete",
                "scope": "low_event_stress_test_only",
                "confirmatory_720_evaluated": False,
                "qc": "passed",
                "run_dir": str(run_dir),
                "source_rows": 2084,
                "stress_rows": 521,
                "stress_positives": 34,
            },
        )
        log(f"Complete with QC passed: {run_dir}")
        manifest_rows = []
        for path in sorted(run_dir.iterdir(), key=lambda item: item.name):
            if path.is_file() and path.name != "final_hash_manifest.csv":
                manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
        atomic_csv(run_dir / "final_hash_manifest.csv", pd.DataFrame(manifest_rows))
    except Exception as error:
        log(f"FAILED: {error}")
        (run_dir / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        atomic_json(
            run_dir / "RUN_STATUS.json",
            {
                "status": "failed",
                "scope": "low_event_stress_test_only",
                "confirmatory_720_evaluated": False,
                "error": str(error),
                "run_dir": str(run_dir),
            },
        )
        raise


if __name__ == "__main__":
    main()
