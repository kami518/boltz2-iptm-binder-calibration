"""RQ3-E confirmatory 720-pair external transport evaluation.

Reuses the previously frozen target-equal deployment U unchanged. Also runs
only the prespecified outcome, seed, and source-weight sensitivity analyses.
No external labels are used for model fitting or tuning; no I/F is fitted.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import shutil
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
DATA = ROOT / "data/external_target_audit_2026-11-13"
EXTERNAL = DATA / "anthropic_audited_pair_table.csv"
COFOLD = DATA / "anthropic_cofold_predictions.parquet"
UPSTREAM_RUN = ROOT / "analysis/rq3e_outputs/run_20261114_162656"
UPSTREAM_MODEL = UPSTREAM_RUN / "rq3e_source_model.json"
OUTPUT_ROOT = ROOT / "analysis/rq3e_outputs"
UPSTREAM_MODEL_SHA = "c3a22f05199729828d9b96051ea8759a64d9f4c632a6558ed6ae84f0f902aabd"
CLIP = 1e-6
EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]
TARGET_ORDER = ["BHRF1", "IL-7Ra", "RBX1", "TrkA", "Nipah-G", "PD-L1", "TREM2", "VEGF-A"]
NEW_TARGETS = {"BHRF1", "IL-7Ra", "RBX1", "TrkA"}
SAME_TARGETS = {"Nipah-G", "PD-L1", "TREM2", "VEGF-A"}
EXPECTED_MAIN = {
    "BHRF1": (90, 20), "IL-7Ra": (90, 47), "RBX1": (90, 23), "TrkA": (90, 19),
    "Nipah-G": (90, 19), "PD-L1": (90, 39), "TREM2": (90, 65), "VEGF-A": (90, 55),
}
OUTCOMES = {
    "adaptyv_weak_plus": ("adaptyv_valid", "adaptyv_main_positive", 720, 287),
    "adaptyv_medium_plus": ("adaptyv_valid", "adaptyv_strict_positive", 720, 254),
    "integrated": ("integrated_valid", "integrated_positive", 720, 304),
    "twist_qc": ("twist_valid", "twist_qc_positive", 718, 247),
    "twist_permissive": ("twist_valid", "twist_permissive_positive", 718, 289),
}
SCORE_COLUMNS = {
    "best_ipsae_seed": "boltz2_iptm_ipsae_best_seed",
    "seed_mean": "boltz2_iptm_seed_mean",
    "seed_median": "boltz2_iptm_seed_median",
    "seed0": "boltz2_iptm_seed0",
}
EXPECTED_HASHES = {
    MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
    ROOT / "analysis/codex_tasks/RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13/rq3_execution_config.json": "a416ecd357eabb4c4d95701136c412cf23651b730bf39feab7524fbc07eeef51",
    COMMON: "276f46ad6d068ecead50a9236d1dcf6953428e415eef79da2424ad850c471f3e",
    DATA / "anthropic_design_summary.parquet": "3e796869b8dc0c90d7ab0e12daf25308c0cf272d04c6a5e27281f903ecd18744",
    DATA / "anthropic_wetlab_summary.parquet": "f80c3cfce67cb7da53576a020d476c5b74efed44e61c6fceaec9ee7de5722f88",
    COFOLD: "ba48010610daa9b9d184aa1ddb419beebdf728ccf745f5a8f3b48e4041e1a87d",
    EXTERNAL: "1f2a4ceab9de06cb75f4aae9a889fc5164e144232101417653944d03ee33d0c6",
    DATA / "anthropic_target_qualification.csv": "bb3ac9e2a0f3bc253d7b500dcc4909f1aa2291026e56d2960f028f8911e50579",
    DATA / "anthropic_audit_summary.json": "df377fc5b755e41a5c2585b14d87f5bf73325b0c8b2d8f94c17381b59c066b81",
    DATA / "anthropic_sequence_independence_summary.json": "8549f59b535357dbe0003a9aa1fa7ffcdadee423bed6de3d9dcf18d0beaf1316",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def fmt(value: object, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
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

    atomic_json(run_dir / "RUN_STATUS.json", {"status": "running", "scope": "confirmatory_720_and_prespecified_sensitivities", "run_dir": str(run_dir)})

    try:
        log("Checking frozen input and upstream deployment-U hashes.")
        actual_hashes = {str(path): sha256(path) for path in EXPECTED_HASHES}
        mismatches = {
            str(path): {"expected": expected, "actual": actual_hashes[str(path)]}
            for path, expected in EXPECTED_HASHES.items() if actual_hashes[str(path)] != expected
        }
        if sha256(UPSTREAM_MODEL) != UPSTREAM_MODEL_SHA:
            mismatches[str(UPSTREAM_MODEL)] = {"expected": UPSTREAM_MODEL_SHA, "actual": sha256(UPSTREAM_MODEL)}
        if mismatches:
            raise RuntimeError(f"Frozen hash mismatch: {mismatches}")

        spec = importlib.util.spec_from_file_location("rq3_common_frozen", COMMON)
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot import frozen rq3_common.py")
        common = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(common)

        primary_model = json.loads(UPSTREAM_MODEL.read_text(encoding="utf-8"))
        if primary_model["fit_status"] != "evaluable" or primary_model["external_labels_used"] is not False:
            raise RuntimeError("Upstream deployment U failed provenance checks")
        primary_intercept = float(primary_model["transport_intercept"])
        primary_beta = float(primary_model["beta_common"])
        primary_reference = float(primary_model["transport_reference"])
        levels = list(primary_model["levels"])
        if len(levels) != 6:
            raise RuntimeError("Upstream U does not have six source targets")

        log("Fitting only the prespecified source-pair-weight sensitivity model.")
        source = pd.read_csv(MASTER, low_memory=False)
        if len(source) != 2084 or source.pair_id.duplicated().any() or source.target.nunique() != 6:
            raise RuntimeError("Frozen source cohort structure mismatch")
        source_targets = source.target.astype(str).to_numpy()
        y_source = source.primary_outcome.to_numpy(dtype=float)
        x_source = common.logit_clip(source.boltz2_iptm.to_numpy(dtype=float), CLIP)
        pair_fit = common.fit_source_linear(
            x_source, y_source, source_targets, np.ones(len(source)), "RQ3E_pair_weight_sensitivity"
        )
        if not pair_fit.get("estimable", False):
            raise RuntimeError(f"Source pair-weight model failed QC: {pair_fit.get('message')}")
        pair_levels = list(pair_fit["levels"])
        pair_theta = np.asarray(pair_fit["theta"], dtype=float)
        source_counts = source.groupby("target", sort=True).size().astype(float)
        ordered_counts = source_counts.loc[pair_levels].to_numpy(dtype=float)
        pair_intercept = float(np.dot(ordered_counts, pair_theta[:6]) / ordered_counts.sum())
        pair_beta = float(pair_theta[6])
        pair_reference = float(y_source.sum() / len(y_source))

        models = {
            "target_equal": {
                "intercept": primary_intercept, "beta": primary_beta, "reference": primary_reference,
                "fit_status": primary_model["fit_status"], "external_labels_used": False,
            },
            "pair": {
                "intercept": pair_intercept, "beta": pair_beta, "reference": pair_reference,
                "fit_status": pair_fit["status"], "external_labels_used": False,
            },
        }
        shutil.copy2(UPSTREAM_MODEL, run_dir / "rq3e_source_model_target_equal_reused.json")
        atomic_json(
            run_dir / "rq3e_source_model_pair_sensitivity.json",
            {
                "model": "unpenalized weighted logistic; six target intercepts and one common logit-ipTM slope",
                "levels": pair_levels,
                "alpha": [float(value) for value in pair_theta[:6]],
                "beta_common": pair_beta,
                "transport_intercept": pair_intercept,
                "transport_reference": pair_reference,
                "fit_status": pair_fit["status"],
                "gradient_max_abs": float(pair_fit["gradient_max_abs"]),
                "source_rows": len(source),
                "external_labels_used": False,
            },
        )
        primary_parameters = pd.read_csv(UPSTREAM_RUN / "rq3e_source_model_parameters.csv")
        primary_parameters.insert(0, "source_weighting", "target_equal")
        primary_parameters.insert(1, "model_role", "primary_reused_unchanged")
        # The frozen common fitter stores inv(mean Hessian), so its SE needs 1/sqrt(N).
        # Correct only the tabulated SE here; frozen coefficients and upstream JSON stay unchanged.
        primary_parameters["standard_error"] = pd.to_numeric(primary_parameters["standard_error"], errors="coerce") / math.sqrt(len(source))
        pair_parameter_rows = []
        pair_covariance = np.asarray(pair_fit["covariance"], dtype=float) / len(source)
        source_prevalence = source.groupby("target", sort=True).primary_outcome.mean()
        for index, target in enumerate(pair_levels):
            pair_parameter_rows.append({
                "source_weighting": "pair", "model_role": "source_only_sensitivity",
                "parameter_name": f"alpha[{target}]", "source_target": target,
                "estimate": float(pair_theta[index]), "standard_error": math.sqrt(max(float(pair_covariance[index, index]), 0.0)),
                "n_source": int(source_counts.loc[target]), "positive_source": int(source.loc[source.target.eq(target), "primary_outcome"].sum()),
                "source_prevalence": float(source_prevalence.loc[target]), "fit_status": pair_fit["status"],
            })
        pair_parameter_rows.extend([
            {"source_weighting": "pair", "model_role": "source_only_sensitivity", "parameter_name": "common_logit_iptm_slope", "source_target": "", "estimate": pair_beta, "standard_error": math.sqrt(max(float(pair_covariance[6, 6]), 0.0)), "n_source": len(source), "positive_source": int(y_source.sum()), "source_prevalence": pair_reference, "fit_status": pair_fit["status"]},
            {"source_weighting": "pair", "model_role": "source_only_sensitivity", "parameter_name": "transport_intercept_pair_count_weighted", "source_target": "", "estimate": pair_intercept, "standard_error": float("nan"), "n_source": len(source), "positive_source": int(y_source.sum()), "source_prevalence": pair_reference, "fit_status": pair_fit["status"]},
        ])
        atomic_csv(run_dir / "rq3e_source_model_parameters.csv", pd.concat([primary_parameters, pd.DataFrame(pair_parameter_rows)], ignore_index=True))

        log("Building the frozen 720-pair external cohort and seed-0 score.")
        external = pd.read_csv(EXTERNAL, low_memory=False)
        integrated_valid = bool_series(external.integrated_valid)
        score_available = pd.to_numeric(external.boltz2_iptm_ipsae_best_seed, errors="coerce").notna()
        cohort = external.loc[external.target.isin(TARGET_ORDER) & integrated_valid & score_available].copy()
        if len(cohort) != 720 or cohort.uuid.duplicated().any() or set(cohort.target) != set(TARGET_ORDER):
            raise RuntimeError("External confirmatory cohort structure mismatch")
        if not bool_series(cohort.adaptyv_valid).all():
            raise RuntimeError("Primary Adaptyv outcome is missing within the 720-pair cohort")

        import duckdb
        cofold = duckdb.connect().execute("SELECT * FROM read_parquet(?)", [str(COFOLD)]).df()
        boltz = cofold.loc[
            cofold.cofolding_model.eq("boltz2")
            & (cofold.target_form.isna() | cofold.target_form.eq(""))
            & cofold.stoichiometry.isin(["1to1", "1to2", "1to3", "rnp"])
        ].copy()
        if boltz.duplicated(["uuid", "seed"]).any():
            raise RuntimeError("Boltz design-target seed rows are not unique")
        seed0 = boltz.loc[pd.to_numeric(boltz.seed, errors="raise").eq(0), ["uuid", "iptm_pae"]].rename(columns={"iptm_pae": "boltz2_iptm_seed0"})
        if seed0.uuid.duplicated().any():
            raise RuntimeError("Seed-0 values are not unique per uuid")
        cohort = cohort.merge(seed0, on="uuid", how="left", validate="one_to_one")
        for score_column in SCORE_COLUMNS.values():
            cohort[score_column] = pd.to_numeric(cohort[score_column], errors="raise").astype(float)
            if not np.isfinite(cohort[score_column]).all():
                raise RuntimeError(f"Missing or nonfinite score: {score_column}")
        for outcome_name, (valid_column, positive_column, _, _) in OUTCOMES.items():
            cohort[f"valid__{outcome_name}"] = bool_series(cohort[valid_column])
            cohort[f"outcome__{outcome_name}"] = bool_series(cohort[positive_column]).astype(int)
        cohort["stratum"] = np.where(cohort.target.isin(NEW_TARGETS), "new_target", "same_target_new_campaign")
        main_n = int(cohort.valid__adaptyv_weak_plus.sum())
        main_pos = int(cohort.loc[cohort.valid__adaptyv_weak_plus, "outcome__adaptyv_weak_plus"].sum())
        if main_n != 720 or main_pos != 287:
            raise RuntimeError("Primary external outcome totals mismatch")
        for target, (expected_n, expected_pos) in EXPECTED_MAIN.items():
            group = cohort.loc[cohort.target.eq(target) & cohort.valid__adaptyv_weak_plus]
            if len(group) != expected_n or int(group.outcome__adaptyv_weak_plus.sum()) != expected_pos:
                raise RuntimeError(f"Primary target count mismatch: {target}")
        for outcome_name, (_, _, expected_n, expected_pos) in OUTCOMES.items():
            valid = cohort[f"valid__{outcome_name}"]
            if int(valid.sum()) != expected_n or int(cohort.loc[valid, f"outcome__{outcome_name}"].sum()) != expected_pos:
                raise RuntimeError(f"Outcome coverage mismatch: {outcome_name}")
        if int(cohort.target.isin(NEW_TARGETS).sum()) != 360 or int(cohort.target.isin(SAME_TARGETS).sum()) != 360:
            raise RuntimeError("External stratum sizes mismatch")

        pair_columns = [
            "uuid", "full_name", "target", "stratum", "design_model", "campaign", "generator", "sequence_design_method",
            *SCORE_COLUMNS.values(),
            *[f"valid__{name}" for name in OUTCOMES], *[f"outcome__{name}" for name in OUTCOMES],
        ]
        atomic_csv(run_dir / "rq3e_external_pairs.csv", cohort[pair_columns].sort_values(["target", "uuid"]))

        configurations = [
            {"config_id": "primary_target_equal_weak_plus_best_ipsae", "analysis_role": "primary", "outcome": "adaptyv_weak_plus", "score": "best_ipsae_seed", "source_weighting": "target_equal"},
            {"config_id": "outcome_medium_plus", "analysis_role": "outcome_sensitivity", "outcome": "adaptyv_medium_plus", "score": "best_ipsae_seed", "source_weighting": "target_equal"},
            {"config_id": "outcome_integrated", "analysis_role": "outcome_sensitivity", "outcome": "integrated", "score": "best_ipsae_seed", "source_weighting": "target_equal"},
            {"config_id": "outcome_twist_qc", "analysis_role": "outcome_sensitivity", "outcome": "twist_qc", "score": "best_ipsae_seed", "source_weighting": "target_equal"},
            {"config_id": "outcome_twist_permissive", "analysis_role": "outcome_sensitivity", "outcome": "twist_permissive", "score": "best_ipsae_seed", "source_weighting": "target_equal"},
            {"config_id": "score_seed_mean", "analysis_role": "score_sensitivity", "outcome": "adaptyv_weak_plus", "score": "seed_mean", "source_weighting": "target_equal"},
            {"config_id": "score_seed_median", "analysis_role": "score_sensitivity", "outcome": "adaptyv_weak_plus", "score": "seed_median", "source_weighting": "target_equal"},
            {"config_id": "score_seed0", "analysis_role": "score_sensitivity", "outcome": "adaptyv_weak_plus", "score": "seed0", "source_weighting": "target_equal"},
            {"config_id": "source_pair_weight", "analysis_role": "source_weight_sensitivity", "outcome": "adaptyv_weak_plus", "score": "best_ipsae_seed", "source_weighting": "pair"},
        ]

        log("Generating and freezing all nine prespecified external prediction sets.")
        prediction_frames = []
        for config in configurations:
            valid = cohort[f"valid__{config['outcome']}"]
            frame = cohort.loc[valid].copy()
            model = models[config["source_weighting"]]
            score = frame[SCORE_COLUMNS[config["score"]]].to_numpy(dtype=float)
            x = common.logit_clip(score, CLIP)
            eta = model["intercept"] + model["beta"] * x
            probability = expit(eta)
            if not np.isfinite(probability).all() or (probability <= 0).any() or (probability >= 1).any():
                raise RuntimeError(f"Invalid probabilities: {config['config_id']}")
            oracle = frame.groupby("target")[f"outcome__{config['outcome']}"].transform("mean").to_numpy(dtype=float)
            predictions = pd.DataFrame({
                "uuid": frame.uuid.astype(str).to_numpy(),
                "target": frame.target.astype(str).to_numpy(),
                "stratum": frame.stratum.astype(str).to_numpy(),
                "config_id": config["config_id"],
                "analysis_role": config["analysis_role"],
                "outcome_name": config["outcome"],
                "score_rule": config["score"],
                "source_weighting": config["source_weighting"],
                "score_value": score,
                "x": x,
                "eta_U": eta,
                "predicted_probability": probability,
                "bss_transport_reference": model["reference"],
                "bss_target_oracle_reference": oracle,
            })
            prediction_frames.append(predictions)
        all_predictions = pd.concat(prediction_frames, ignore_index=True)
        atomic_csv(run_dir / "rq3e_u_predictions_frozen.csv", all_predictions.sort_values(["config_id", "target", "uuid"]))
        log("Predictions frozen; now joining outcomes only for evaluation.")

        outcome_lookup = {}
        for outcome_name in OUTCOMES:
            valid = cohort[f"valid__{outcome_name}"]
            outcome_lookup[outcome_name] = cohort.loc[valid, ["uuid", f"outcome__{outcome_name}"]].rename(columns={f"outcome__{outcome_name}": "outcome_value"})
        evaluated_frames = []
        for config in configurations:
            pred = all_predictions.loc[all_predictions.config_id.eq(config["config_id"])].copy()
            evaluated = pred.merge(outcome_lookup[config["outcome"]], on="uuid", how="left", validate="many_to_one")
            if len(evaluated) != len(pred) or evaluated.outcome_value.isna().any():
                raise RuntimeError(f"Outcome join failed: {config['config_id']}")
            evaluated_frames.append(evaluated)
        evaluated_all = pd.concat(evaluated_frames, ignore_index=True)
        atomic_csv(run_dir / "rq3e_u_pair_predictions.csv", evaluated_all.sort_values(["config_id", "target", "uuid"]))

        metric_rows = []
        reliability_frames = []
        for config in configurations:
            evaluated = evaluated_all.loc[evaluated_all.config_id.eq(config["config_id"])]
            for target in TARGET_ORDER:
                group = evaluated.loc[evaluated.target.eq(target)]
                y = group.outcome_value.to_numpy(dtype=float)
                p = group.predicted_probability.to_numpy(dtype=float)
                transport = group.bss_transport_reference.to_numpy(dtype=float)
                oracle = group.bss_target_oracle_reference.to_numpy(dtype=float)
                metrics = common.compute_performance_metrics(y, p, clip=CLIP, transport_reference=transport, target_oracle_reference=oracle)
                # Unit weights: convert inv(mean Hessian) SE to inv(sum Hessian) SE.
                for se_name in ["calibration_intercept_se", "calibration_slope_se"]:
                    metrics[se_name] = float(metrics[se_name] / math.sqrt(len(y))) if np.isfinite(metrics[se_name]) else float("nan")
                low, high = common.wilson_interval(int(y.sum()), len(y))
                metrics.update({
                    "config_id": config["config_id"], "analysis_role": config["analysis_role"],
                    "outcome_name": config["outcome"], "score_rule": config["score"], "source_weighting": config["source_weighting"],
                    "target": target, "stratum": "new_target" if target in NEW_TARGETS else "same_target_new_campaign",
                    "prevalence_wilson_95_low": low, "prevalence_wilson_95_high": high,
                    "mean_prediction": float(p.mean()), "median_prediction": float(np.median(p)),
                    "transport_reference_value": float(transport[0]), "target_oracle_reference_value": float(oracle[0]),
                    "brier_transport_reference": common.brier_score(y, transport),
                    "brier_target_oracle_reference": common.brier_score(y, oracle),
                })
                metric_rows.append(metrics)
                reliability = pd.DataFrame(common.reliability_rows(group, "predicted_probability", EDGES))
                reliability.insert(0, "target", target)
                reliability.insert(0, "config_id", config["config_id"])
                reliability_frames.append(reliability)
        target_metrics = pd.DataFrame(metric_rows)
        target_metric_columns = [
            "config_id", "analysis_role", "outcome_name", "score_rule", "source_weighting", "target", "stratum",
            "n", "positives", "prevalence", "prevalence_wilson_95_low", "prevalence_wilson_95_high",
            "mean_prediction", "median_prediction", "brier", "log_loss", "bss_transport", "bss_target_oracle",
            "transport_reference_value", "target_oracle_reference_value", "brier_transport_reference", "brier_target_oracle_reference",
            "calibration_intercept", "calibration_intercept_se", "calibration_slope", "calibration_slope_se",
            "calibration_status", "calibration_failure_type", "calibration_message",
        ]
        target_metrics = target_metrics[target_metric_columns]
        reliability_all = pd.concat(reliability_frames, ignore_index=True)
        atomic_csv(run_dir / "rq3e_target_metrics.csv", target_metrics)
        atomic_csv(run_dir / "rq3e_reliability.csv", reliability_all)

        log("Computing target-equal, stratum, and pair-pooled summaries.")
        group_rows = []
        metric_names = ["brier", "log_loss", "bss_transport", "bss_target_oracle", "calibration_intercept", "calibration_slope"]
        for config in configurations:
            config_metrics = target_metrics.loc[target_metrics.config_id.eq(config["config_id"])]
            config_pairs = evaluated_all.loc[evaluated_all.config_id.eq(config["config_id"])]
            for aggregation, stratum, targets in [
                ("target_equal_macro", "all_8", set(TARGET_ORDER)),
                ("target_equal_macro", "new_target_4", NEW_TARGETS),
                ("target_equal_macro", "same_target_campaign_4", SAME_TARGETS),
            ]:
                subset = config_metrics.loc[config_metrics.target.isin(targets)]
                row = {
                    "config_id": config["config_id"], "analysis_role": config["analysis_role"],
                    "outcome_name": config["outcome"], "score_rule": config["score"], "source_weighting": config["source_weighting"],
                    "aggregation": aggregation, "stratum": stratum, "n_targets": len(subset),
                    "n_pairs": int(subset.n.sum()), "positives": int(subset.positives.sum()),
                    "prevalence": float(subset.prevalence.mean()),
                    "calibration_evaluable_targets": int(subset.calibration_status.eq("evaluable").sum()),
                    "bss_transport_positive_targets": int(subset.bss_transport.gt(0).sum()),
                    "bss_transport_negative_targets": int(subset.bss_transport.lt(0).sum()),
                    "bss_target_oracle_positive_targets": int(subset.bss_target_oracle.gt(0).sum()),
                    "bss_target_oracle_negative_targets": int(subset.bss_target_oracle.lt(0).sum()),
                }
                for name in metric_names:
                    values = pd.to_numeric(subset[name], errors="coerce").dropna().to_numpy(dtype=float)
                    row[name] = float(values.mean()) if len(values) else float("nan")
                    row[f"{name}_min"] = float(values.min()) if len(values) else float("nan")
                    row[f"{name}_max"] = float(values.max()) if len(values) else float("nan")
                group_rows.append(row)
            y = config_pairs.outcome_value.to_numpy(dtype=float)
            p = config_pairs.predicted_probability.to_numpy(dtype=float)
            transport = config_pairs.bss_transport_reference.to_numpy(dtype=float)
            oracle = config_pairs.bss_target_oracle_reference.to_numpy(dtype=float)
            pooled = common.compute_performance_metrics(y, p, clip=CLIP, transport_reference=transport, target_oracle_reference=oracle)
            for se_name in ["calibration_intercept_se", "calibration_slope_se"]:
                pooled[se_name] = float(pooled[se_name] / math.sqrt(len(y))) if np.isfinite(pooled[se_name]) else float("nan")
            pooled.update({
                "config_id": config["config_id"], "analysis_role": config["analysis_role"],
                "outcome_name": config["outcome"], "score_rule": config["score"], "source_weighting": config["source_weighting"],
                "aggregation": "pair_pooled", "stratum": "all_available_pairs", "n_targets": int(config_pairs.target.nunique()),
                "n_pairs": len(config_pairs), "positives": int(y.sum()),
                "calibration_evaluable_targets": int(pooled["calibration_status"] == "evaluable"),
                "bss_transport_positive_targets": int(pooled["bss_transport"] > 0),
                "bss_transport_negative_targets": int(pooled["bss_transport"] < 0),
                "bss_target_oracle_positive_targets": int(pooled["bss_target_oracle"] > 0),
                "bss_target_oracle_negative_targets": int(pooled["bss_target_oracle"] < 0),
            })
            for name in metric_names:
                pooled[f"{name}_min"] = float("nan")
                pooled[f"{name}_max"] = float("nan")
            group_rows.append(pooled)
        group_metrics = pd.DataFrame(group_rows)
        atomic_csv(run_dir / "rq3e_group_metrics.csv", group_metrics)

        log("Computing campaign and generator descriptive summaries for the primary analysis.")
        primary_pairs = evaluated_all.loc[evaluated_all.config_id.eq("primary_target_equal_weak_plus_best_ipsae")].copy()
        primary_with_meta = primary_pairs.merge(
            cohort[["uuid", "design_model", "campaign", "generator", "sequence_design_method"]],
            on="uuid", how="left", validate="many_to_one",
        )
        campaign_rows = []
        for (target, design_model, campaign), group in primary_with_meta.groupby(["target", "design_model", "campaign"], sort=True):
            y = group.outcome_value.to_numpy(dtype=float)
            p = group.predicted_probability.to_numpy(dtype=float)
            campaign_rows.append({
                "target": target, "design_model": design_model, "campaign": campaign,
                "n": len(group), "positives": int(y.sum()), "observed_prevalence": float(y.mean()),
                "mean_prediction": float(p.mean()), "brier": common.brier_score(y, p),
                "log_loss": common.log_loss_score(y, p, CLIP),
                "generators_observed": ";".join(sorted(group.generator.astype(str).unique())),
                "sequence_design_methods": ";".join(sorted(group.sequence_design_method.astype(str).unique())),
            })
        atomic_csv(run_dir / "rq3e_campaign_metrics.csv", pd.DataFrame(campaign_rows))
        generator_rows = []
        for (design_model, generator), group in primary_with_meta.groupby(["design_model", "generator"], sort=True):
            generator_rows.append({
                "design_model": design_model, "generator": generator, "n": len(group),
                "positives": int(group.outcome_value.sum()), "observed_prevalence": float(group.outcome_value.mean()),
                "mean_prediction": float(group.predicted_probability.mean()),
            })
        atomic_csv(run_dir / "rq3e_generator_composition.csv", pd.DataFrame(generator_rows))

        log("Drawing primary target, reliability, and sensitivity figures.")
        primary_metrics = target_metrics.loc[target_metrics.config_id.eq("primary_target_equal_weak_plus_best_ipsae")].set_index("target").loc[TARGET_ORDER]
        fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True, sharey=True, constrained_layout=True)
        for axis, target in zip(axes.ravel(), TARGET_ORDER):
            rel = reliability_all.loc[reliability_all.config_id.eq("primary_target_equal_weak_plus_best_ipsae") & reliability_all.target.eq(target) & reliability_all.n.gt(0)]
            met = primary_metrics.loc[target]
            lower = np.maximum((rel.observed_rate - rel.wilson_95_low).to_numpy(dtype=float), 0.0)
            upper = np.maximum((rel.wilson_95_high - rel.observed_rate).to_numpy(dtype=float), 0.0)
            axis.plot([0, 1], [0, 1], "--", color="0.65", linewidth=1)
            axis.errorbar(rel.mean_prediction, rel.observed_rate, yerr=[lower, upper], fmt="o-", markersize=4, linewidth=1.1, capsize=2, color="#2166ac")
            axis.set_xlim(0, 1); axis.set_ylim(0, 1); axis.grid(alpha=0.15)
            axis.set_title(f"{target}: n={int(met.n)}, pos={int(met.positives)}\nBrier={met.brier:.3f}; BSS-T={met.bss_transport:.3f}", fontsize=10)
        fig.suptitle("RQ3-E confirmatory cohort: frozen U reliability", fontsize=14)
        fig.savefig(run_dir / "RQ3E_U_external_reliability.png", dpi=220)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
        x = np.arange(len(TARGET_ORDER))
        colors = ["#1b9e77" if target in NEW_TARGETS else "#7570b3" for target in TARGET_ORDER]
        axes[0].bar(x, primary_metrics.bss_transport, color=colors)
        axes[0].axhline(0, color="black", linewidth=0.8)
        axes[0].set_xticks(x, TARGET_ORDER, rotation=35, ha="right")
        axes[0].set_ylabel("BSS-transport")
        axes[0].set_title("Improvement over source transport reference")
        axes[1].errorbar(primary_metrics.mean_prediction, primary_metrics.prevalence, fmt="o", color="#2166ac")
        for _, row in primary_metrics.iterrows():
            axes[1].annotate(row.name, (row.mean_prediction, row.prevalence), xytext=(4, 4), textcoords="offset points", fontsize=8)
        axes[1].plot([0, 1], [0, 1], "--", color="0.65")
        axes[1].set_xlim(0, 1); axes[1].set_ylim(0, 1)
        axes[1].set_xlabel("Mean predicted probability"); axes[1].set_ylabel("Observed prevalence")
        axes[1].set_title("Calibration in the large")
        fig.savefig(run_dir / "RQ3E_U_external_target_metrics.png", dpi=220)
        plt.close(fig)

        sensitivity = group_metrics.loc[group_metrics.aggregation.eq("target_equal_macro") & group_metrics.stratum.eq("all_8")].copy()
        sensitivity["label"] = sensitivity.config_id.str.replace("primary_target_equal_weak_plus_best_ipsae", "PRIMARY", regex=False)
        fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
        ypos = np.arange(len(sensitivity))
        axes[0].barh(ypos, sensitivity.brier, color="#4c78a8")
        axes[0].set_yticks(ypos, sensitivity.label); axes[0].invert_yaxis(); axes[0].set_xlabel("Target-equal macro Brier")
        axes[1].barh(ypos, sensitivity.bss_transport, color=["#1b9e77" if role == "primary" else "#9e9e9e" for role in sensitivity.analysis_role])
        axes[1].axvline(0, color="black", linewidth=0.8); axes[1].set_yticks(ypos, sensitivity.label); axes[1].invert_yaxis(); axes[1].set_xlabel("Target-equal macro BSS-transport")
        fig.savefig(run_dir / "RQ3E_U_external_sensitivity.png", dpi=220)
        plt.close(fig)

        log("Writing QC, frozen configuration, and result tables.")
        if sha256(MASTER) != EXPECTED_HASHES[MASTER] or sha256(UPSTREAM_MODEL) != UPSTREAM_MODEL_SHA:
            raise RuntimeError("Frozen source model or master changed during execution")
        qc = {
            "status": "passed",
            "scope": "confirmatory_720_and_prespecified_sensitivities",
            "input_hash_gate": "passed",
            "upstream_target_equal_model_sha256": sha256(UPSTREAM_MODEL),
            "target_equal_model_reused_without_refit": True,
            "pair_sensitivity_fit_used_source_rows_only": True,
            "external_rows": len(cohort),
            "external_unique_uuids": int(cohort.uuid.nunique()),
            "external_targets": int(cohort.target.nunique()),
            "external_main_positives": main_pos,
            "new_target_pairs": int(cohort.target.isin(NEW_TARGETS).sum()),
            "same_target_campaign_pairs": int(cohort.target.isin(SAME_TARGETS).sum()),
            "outcome_coverage": {
                name: {"n": int(cohort[f"valid__{name}"].sum()), "positives": int(cohort.loc[cohort[f"valid__{name}"], f"outcome__{name}"].sum())}
                for name in OUTCOMES
            },
            "score_rules": list(SCORE_COLUMNS),
            "seed0_from_raw_design_target_stoichiometry_rows": True,
            "configuration_count": len(configurations),
            "prediction_rows": len(all_predictions),
            "target_metric_rows": len(target_metrics),
            "group_metric_rows": len(group_metrics),
            "reliability_rows": len(reliability_all),
            "all_predictions_finite_and_inside_open_unit_interval": True,
            "external_design_models": int(cohort.design_model.nunique()),
            "external_generators": int(cohort.generator.nunique()),
            "I_or_F_fitted": False,
            "Pinc_or_cross_predictor_or_ML_run": False,
            "master_hash_after": sha256(MASTER),
            "result_driven_changes": False,
            "standard_error_correction": "inv(mean Hessian) SE divided by sqrt(sum of unit weights); point estimates unchanged",
        }
        if len(configurations) != 9 or len(target_metrics) != 72 or len(reliability_all) != 504:
            raise RuntimeError("Prespecified analysis output dimensions mismatch")
        atomic_json(run_dir / "qc_checks.json", qc)
        atomic_csv(run_dir / "failed_items.csv", pd.DataFrame(columns=["stage", "item", "failure_type", "message"]))

        config_snapshot = {
            "protocol_path": str(PROTOCOL), "protocol_sha256": sha256(PROTOCOL),
            "execution_script_sha256": sha256(Path(__file__)), "frozen_common_sha256": sha256(COMMON),
            "original_weight_rule_reference": str(PREP), "original_weight_rule_sha256": sha256(PREP),
            "upstream_target_equal_model": str(UPSTREAM_MODEL), "upstream_target_equal_model_sha256": sha256(UPSTREAM_MODEL),
            "target_equal_model_reused_without_refit": True,
            "scope": "8-target / 720-pair external primary analysis plus prespecified sensitivities",
            "primary_config": configurations[0], "configurations": configurations,
            "models": models, "targets": TARGET_ORDER, "new_targets": sorted(NEW_TARGETS), "same_target_campaigns": sorted(SAME_TARGETS),
            "outcomes": {name: {"valid_column": value[0], "positive_column": value[1], "expected_n": value[2], "expected_positives": value[3]} for name, value in OUTCOMES.items()},
            "score_columns": SCORE_COLUMNS, "reliability_edges": EDGES, "probability_clip": CLIP,
            "input_hashes": actual_hashes,
            "external_labels_used_for_fitting_tuning_or_selection": False,
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        }
        atomic_json(run_dir / "run_config.snapshot.json", config_snapshot)

        primary_group = group_metrics.loc[group_metrics.config_id.eq("primary_target_equal_weak_plus_best_ipsae")]
        result_lines = [
            "# RQ3-E 720 行确认性外部迁移结果",
            "",
            "**范围：8 靶点／720 行主外部队列及预定义敏感性分析；target-equal U 从上一运行原样复用，未使用外部标签训练或调参，未拟合 I／F。**",
            "",
            f"主 U：运输截距 **{primary_intercept:.6f}**，共享 logit-ipTM 斜率 **{primary_beta:.6f}**，源靶点等权基率参照 **{primary_reference:.6f}**。",
            "",
            "## 主分析：逐靶点",
            "",
            "| 靶点 | 分层 | n／阳性 | 基率 | 平均预测 | Brier | Log loss | BSS-transport | BSS-target-oracle | 校准截距 | 校准斜率 |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for target, row in primary_metrics.iterrows():
            result_lines.append(
                f"| {target} | {row.stratum} | {int(row.n)}／{int(row.positives)} | {fmt(row.prevalence)} | {fmt(row.mean_prediction)} | "
                f"{fmt(row.brier)} | {fmt(row.log_loss)} | {fmt(row.bss_transport)} | {fmt(row.bss_target_oracle)} | "
                f"{fmt(row.calibration_intercept)} | {fmt(row.calibration_slope)} |"
            )
        result_lines.extend([
            "", "## 主分析：预定义汇总", "",
            "| 汇总 | n 靶点／配对 | 阳性 | Brier | Log loss | BSS-transport | BSS-target-oracle | 校准截距 | 校准斜率 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for _, row in primary_group.iterrows():
            result_lines.append(
                f"| {row.aggregation}／{row.stratum} | {int(row.n_targets)}／{int(row.n_pairs)} | {int(row.positives)} | "
                f"{fmt(row.brier)} | {fmt(row.log_loss)} | {fmt(row.bss_transport)} | {fmt(row.bss_target_oracle)} | {fmt(row.calibration_intercept)} | {fmt(row.calibration_slope)} |"
            )
        result_lines.extend([
            "", "## 敏感性分析：8 靶点等权宏平均", "",
            "| 配置 | 类型 | n 配对 | Brier | Log loss | BSS-transport | BSS-target-oracle | BSS-T 为正靶点数 |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ])
        for _, row in sensitivity.iterrows():
            result_lines.append(
                f"| {row.config_id} | {row.analysis_role} | {int(row.n_pairs)} | {fmt(row.brier)} | {fmt(row.log_loss)} | "
                f"{fmt(row.bss_transport)} | {fmt(row.bss_target_oracle)} | {int(row.bss_transport_positive_targets)}／8 |"
            )
        result_lines.extend([
            "", "## 强制解释边界", "",
            "- 确认性结论仅适用于 Adaptyv Weak+ 基率约 21.1%–72.2% 的 8 个达标靶点；低事件压力测试另行报告，高基率端仍未覆盖。",
            "- 本研究评价已发布 score＋冻结 U 的运输，不是完全同一 score-generation protocol 下的纯靶点验证。",
            "- 外部队列有两种高层 design model、九种底层 generator；不识别独立 generator 因果效应。",
            "- 靶点等权宏平均、4 新增／4 同靶点分层与配对加权汇总不得混用；不同结局的绝对 Brier 不作同基率直接比较。",
            "- BSS>0 只表示优于对应参照；校准斜率置信区间包含 1 不能单独作为校准成功证据。",
            "- 未运行 Pinc、跨预测器、其他 ML 或任何外部标签适配。",
        ])
        (run_dir / "RQ3E_confirmatory_720_results.md").write_text("\n".join(result_lines) + "\n", encoding="utf-8")

        shutil.copy2(__file__, run_dir / "rq3e_confirmatory_720.py")
        shutil.copy2(COMMON, run_dir / "rq3_common_snapshot.py")
        shutil.copy2(PREP, run_dir / "phase_ab_weight_rule_snapshot.py")
        atomic_json(run_dir / "RUN_STATUS.json", {
            "status": "complete", "scope": "confirmatory_720_and_prespecified_sensitivities", "qc": "passed",
            "run_dir": str(run_dir), "external_rows": 720, "external_main_positives": 287,
            "configuration_count": 9, "target_equal_model_reused_without_refit": True,
        })
        log(f"Complete with QC passed: {run_dir}")
        manifest_rows = []
        for path in sorted(run_dir.iterdir(), key=lambda item: item.name):
            if path.is_file() and path.name != "final_hash_manifest.csv":
                manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
        atomic_csv(run_dir / "final_hash_manifest.csv", pd.DataFrame(manifest_rows))
    except Exception as error:
        log(f"FAILED: {error}")
        (run_dir / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        atomic_json(run_dir / "RUN_STATUS.json", {"status": "failed", "error": str(error), "run_dir": str(run_dir)})
        raise


if __name__ == "__main__":
    main()
