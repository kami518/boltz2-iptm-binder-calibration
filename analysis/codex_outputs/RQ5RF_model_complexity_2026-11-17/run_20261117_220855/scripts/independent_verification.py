"""Independent RQ5-RF verification; intentionally avoids production summaries."""
from __future__ import annotations

import hashlib
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = ROOT / "analysis" / "codex_outputs" / "RQ5RF_model_complexity_2026-11-17" / "run_20261117_220855"
CONFIG = json.loads((RUN / "run_config.snapshot.json").read_text(encoding="utf-8"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def own_logit(probability: np.ndarray, clip: float = 1e-6) -> np.ndarray:
    p = np.asarray(probability, dtype=float)
    q = np.clip(p, clip, 1.0 - clip)
    return np.log(q) - np.log1p(-q)


def own_brier(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean((np.asarray(p, dtype=float) - np.asarray(y, dtype=float)) ** 2))


def own_log_loss(y: np.ndarray, p: np.ndarray, clip: float) -> float:
    yy = np.asarray(y, dtype=float)
    pp = np.clip(np.asarray(p, dtype=float), clip, 1.0 - clip)
    return float(-np.mean(yy * np.log(pp) + (1.0 - yy) * np.log1p(-pp)))


def own_bss(y: np.ndarray, p: np.ndarray, reference: np.ndarray) -> float:
    yy = np.asarray(y, dtype=float)
    denominator = float(np.sum((np.asarray(reference, dtype=float) - yy) ** 2))
    if denominator <= 0.0 or not np.isfinite(denominator):
        return float("nan")
    return float(1.0 - np.sum((np.asarray(p, dtype=float) - yy) ** 2) / denominator)


def own_calibration(y: np.ndarray, p: np.ndarray, clip: float) -> tuple[float, float]:
    yy = np.asarray(y, dtype=float)
    q = own_logit(p, clip)
    design = np.column_stack([np.ones(len(q)), q])

    def objective(theta: np.ndarray) -> float:
        eta = design @ theta
        return float(np.mean(np.logaddexp(0.0, eta) - yy * eta))

    def gradient(theta: np.ndarray) -> np.ndarray:
        return design.T @ (expit(design @ theta) - yy) / len(yy)

    result = minimize(
        objective,
        np.array([0.0, 1.0]),
        method="BFGS",
        jac=gradient,
        options={"gtol": 1e-10, "maxiter": 5000},
    )
    if not result.success:
        return float("nan"), float("nan")
    return float(result.x[0]), float(result.x[1])


def own_metrics(
    y: np.ndarray,
    p: np.ndarray,
    reference_transport: np.ndarray,
    reference_oracle: np.ndarray,
    reference_adaptation: np.ndarray | None,
) -> dict:
    yy = np.asarray(y, dtype=float)
    pp = np.asarray(p, dtype=float)
    calibration_intercept, calibration_slope = own_calibration(yy, pp, 1e-6)
    return {
        "brier": own_brier(yy, pp),
        "log_loss": own_log_loss(yy, pp, 1e-6),
        "bss_transport": own_bss(yy, pp, reference_transport),
        "bss_target_oracle": own_bss(yy, pp, reference_oracle),
        "bss_adaptation_train": (
            own_bss(yy, pp, reference_adaptation)
            if reference_adaptation is not None
            else float("nan")
        ),
        "roc_auc": (
            float(roc_auc_score(yy, pp)) if np.unique(yy).size == 2 else float("nan")
        ),
        "ap": float(average_precision_score(yy, pp)) if np.sum(yy) > 0 else float("nan"),
        "calibration_intercept": calibration_intercept,
        "calibration_slope": calibration_slope,
    }


def expected_targets(analysis_set: str) -> set[str]:
    values = {
        "aah-i",
        "human-vegf110",
        "nipah-glycoprotein-g",
        "pd-l1",
        "sars-cov-2-rbd",
        "trem2",
    }
    if analysis_set == "primary_five":
        values.remove("aah-i")
    return values


def analysis_frames() -> dict[str, pd.DataFrame]:
    master = pd.read_csv(CONFIG["inputs"]["pair_master"]["path"], low_memory=False)
    pinc = pd.read_csv(CONFIG["inputs"]["pinc_extension"]["path"], low_memory=False)
    clusters = pd.read_csv(CONFIG["inputs"]["sequence_clusters"]["path"], low_memory=False)
    joined = master.merge(
        pinc[["pair_id", "pinc", "pinc_status"]],
        on="pair_id",
        how="left",
        validate="one_to_one",
    ).merge(
        clusters[["pair_id", "sequence_cluster90"]],
        on="pair_id",
        how="left",
        validate="one_to_one",
    )
    joined["x_ipTM"] = np.nan
    joined["x_Pinc"] = np.nan
    success = joined["pinc_status"].eq("success")
    joined.loc[success, "x_ipTM"] = own_logit(
        joined.loc[success, "boltz2_iptm"].to_numpy(dtype=float)
    )
    joined.loc[success, "x_Pinc"] = own_logit(
        joined.loc[success, "pinc"].to_numpy(dtype=float)
    )
    joined["sequence_cluster_value"] = joined["sequence_cluster90"].astype(str)
    frames = {}
    for analysis_set in CONFIG["analysis_sets"]:
        spec = CONFIG["analysis_sets"][analysis_set]
        frame = joined[success].copy()
        if spec["exclude_target"] is not None:
            frame = frame[frame["target"].astype(str).ne(spec["exclude_target"])].copy()
        frame = frame[frame[spec["outcome"]].notna()].copy()
        frame["outcome_value"] = frame[spec["outcome"]].astype(int)
        frames[analysis_set] = frame.sort_values("pair_id", kind="mergesort").reset_index(drop=True)
    return frames


def selection_checks() -> dict:
    fits = pd.read_csv(RUN / "rf_inner_candidate_fits.csv", low_memory=False)
    selection = pd.read_csv(RUN / "rf_inner_candidate_selection.csv", low_memory=False)
    manifest = pd.read_csv(RUN / "shared_inner_task_manifest.csv", low_memory=False)
    expected_counts = manifest.groupby(
        ["analysis_set", "source_weighting", "held_out_target"]
    ).size()
    checks = []
    selected_ids = []
    for keys, group in fits.groupby(
        [
            "analysis_set",
            "source_weighting",
            "held_out_target",
            "model_exposure",
            "candidate_id",
        ],
        sort=True,
    ):
        analysis_set, weighting, held_target, exposure, candidate_id = keys
        expected = int(expected_counts.loc[(analysis_set, weighting, held_target)])
        evaluable = group["fit_status"].eq("evaluable").all()
        loss = (
            float(pd.to_numeric(group["validation_log_loss"], errors="coerce").mean())
            if evaluable and len(group) == expected
            else float("nan")
        )
        recorded = selection[
            selection["analysis_set"].eq(analysis_set)
            & selection["source_weighting"].eq(weighting)
            & selection["held_out_target"].eq(held_target)
            & selection["model_exposure"].eq(exposure)
            & selection["candidate_id"].astype(str).eq(str(candidate_id))
        ]
        if len(recorded) != 1:
            checks.append({"pass": False, "reason": "recorded candidate count", "keys": keys})
            continue
        recorded_value = float(recorded.iloc[0]["macro_mean_inner_log_loss"])
        checks.append(
            {
                "pass": bool(
                    len(group) == expected
                    and np.isfinite(loss)
                    and np.isfinite(recorded_value)
                    and abs(loss - recorded_value) <= 1e-12
                ),
                "family": str(group["model_family"].iloc[0]),
                "keys": keys,
                "recomputed": loss,
                "recorded": recorded_value,
            }
        )
    for (analysis_set, weighting, held_target, exposure), group in selection.groupby(
        ["analysis_set", "source_weighting", "held_out_target", "model_exposure"],
        sort=True,
    ):
        evaluable = group[group["candidate_status"].eq("evaluable")]
        best = float(evaluable["macro_mean_inner_log_loss"].min())
        ties = evaluable[
            evaluable["macro_mean_inner_log_loss"].astype(float) - best <= 1e-12
        ]
        family = str(group["model_family"].iloc[0])
        if family == "pooled_logistic":
            expected_id = str(ties.sort_values("lambda").iloc[-1]["candidate_id"])
        else:
            parameters = ties["parameters_json"].map(json.loads)
            tie_frame = ties.copy()
            tie_frame["depth_rank"] = parameters.map(
                lambda value: {2: 0, 4: 1, None: 2}[value["max_depth"]]
            )
            tie_frame["leaf"] = parameters.map(lambda value: value["min_samples_leaf"])
            tie_frame["features"] = parameters.map(lambda value: value["max_features"])
            tie_frame = tie_frame.sort_values(
                ["depth_rank", "leaf", "features"],
                ascending=[True, False, True],
                kind="mergesort",
            )
            expected_id = str(tie_frame.iloc[0]["candidate_id"])
        recorded_id = group[group["selected"].astype(bool)]["candidate_id"].astype(str).tolist()
        selected_ids.append((analysis_set, weighting, held_target, exposure, expected_id))
        checks.append(
            {
                "pass": bool(recorded_id == [expected_id]),
                "reason": "selection tie rule",
                "keys": (analysis_set, weighting, held_target, exposure),
                "expected_id": expected_id,
                "recorded_id": recorded_id,
            }
        )
    all_pass = bool(
        len(selected_ids) == 136
        and all(item["pass"] for item in checks)
        and len(fits) == 4800
    )
    return {
        "pass": all_pass,
        "candidate_rows": len(fits),
        "selected_exposure_tasks": len(selected_ids),
        "checks": len(checks),
        "failed_checks": [item for item in checks if not item["pass"]][:20],
    }


def rerun_final_predictions(frames: dict[str, pd.DataFrame]) -> dict:
    task_id = "primary_six|target_equal|aah-i"
    analysis_set, weighting, held_target = task_id.split("|")
    fold_id = f"loto|{analysis_set}|{CONFIG['cluster_rule']}|{held_target}"
    loto = pd.read_csv(CONFIG["inputs"]["loto_assignments"]["path"], low_memory=False)
    exclusions = pd.read_csv(CONFIG["inputs"]["cross_target_exclusions"]["path"], low_memory=False)
    valid = frames[analysis_set]
    fold = loto[loto["fold_id"].eq(fold_id)].copy()
    fold = fold[fold["pair_id"].astype(str).isin(set(valid["pair_id"].astype(str)))]
    train_ids = set(fold.loc[fold["split_status"].eq("train"), "pair_id"].astype(str))
    train_ids -= set(exclusions.loc[exclusions["fold_id"].eq(fold_id), "pair_id"].astype(str))
    test_ids = set(fold.loc[fold["split_status"].eq("test"), "pair_id"].astype(str))
    features = ["pair_id", "target", "outcome_value", "x_ipTM", "x_Pinc", "sequence_cluster_value"]
    train = valid[valid["pair_id"].astype(str).isin(train_ids)][features].sort_values(
        "pair_id", kind="mergesort"
    ).reset_index(drop=True)
    test = valid[valid["pair_id"].astype(str).isin(test_ids)][features].sort_values(
        "pair_id", kind="mergesort"
    ).reset_index(drop=True)
    counts = train["target"].astype(str).value_counts().to_dict()
    target_count = len(counts)
    weights = np.asarray(
        [1.0 / (target_count * counts[target]) for target in train["target"].astype(str)],
        dtype=float,
    )
    canonical = (
        f"RQ5RF|v1|{analysis_set}|{weighting}|{fold_id}|none|final"
    )
    seed = int.from_bytes(hashlib.sha256(canonical.encode("utf-8")).digest()[:4], "big")
    source = pd.read_csv(RUN / "rf_source_model_parameters.csv", low_memory=False)
    checkpoint_path = RUN / "checkpoints" / "primary_six__target_equal__aah-i.json.gz"
    with gzip.open(checkpoint_path, "rt", encoding="utf-8") as handle:
        checkpoint = json.load(handle)
    checkpoint_predictions = pd.DataFrame(checkpoint["u_predictions"])
    checks = []
    for exposure in ("rf_ipTM", "rf_combo"):
        source_row = source[
            source["analysis_set"].eq(analysis_set)
            & source["source_weighting"].eq(weighting)
            & source["held_out_target"].eq(held_target)
            & source["model_exposure"].eq(exposure)
        ]
        if len(source_row) != 1:
            checks.append({"exposure": exposure, "pass": False, "reason": "source row"})
            continue
        parameters = json.loads(source_row.iloc[0]["parameters_json"])
        columns = ["x_ipTM"] if exposure == "rf_ipTM" else ["x_ipTM", "x_Pinc"]
        estimator = RandomForestClassifier(
            n_estimators=500,
            criterion="gini",
            max_depth=parameters["max_depth"],
            min_samples_split=2,
            min_samples_leaf=parameters["min_samples_leaf"],
            min_weight_fraction_leaf=0.0,
            max_features=parameters["max_features"],
            max_leaf_nodes=None,
            min_impurity_decrease=0.0,
            bootstrap=True,
            oob_score=False,
            n_jobs=1,
            random_state=seed,
            verbose=0,
            warm_start=False,
            class_weight=None,
            max_samples=None,
            ccp_alpha=0.0,
        )
        estimator.fit(
            train[columns].to_numpy(dtype=float),
            train["outcome_value"].to_numpy(dtype=int),
            sample_weight=weights,
        )
        if not np.array_equal(estimator.classes_, np.array([0, 1])):
            checks.append({"exposure": exposure, "pass": False, "reason": "classes"})
            continue
        positive_index = int(np.flatnonzero(estimator.classes_ == 1)[0])
        recomputed = estimator.predict_proba(test[columns].to_numpy(dtype=float))[:, positive_index]
        stored = checkpoint_predictions[
            checkpoint_predictions["analysis_set"].eq(analysis_set)
            & checkpoint_predictions["source_weighting"].eq(weighting)
            & checkpoint_predictions["held_out_target"].eq(held_target)
            & checkpoint_predictions["model_exposure"].eq(exposure)
        ].sort_values("pair_id", kind="mergesort")
        stored = stored[stored["pair_id"].astype(str).isin(set(test["pair_id"].astype(str)))]
        checks.append(
            {
                "exposure": exposure,
                "features": columns,
                "seed": seed,
                "n_train": len(train),
                "n_test": len(test),
                "pair_ids_identical": stored["pair_id"].astype(str).tolist()
                == test["pair_id"].astype(str).tolist(),
                "max_abs_probability_difference": float(
                    np.max(np.abs(recomputed - stored["p_U"].to_numpy(dtype=float)))
                )
                if len(stored) == len(test)
                else float("nan"),
                "elementwise_identical": bool(
                    len(stored) == len(test)
                    and stored["pair_id"].astype(str).tolist()
                    == test["pair_id"].astype(str).tolist()
                    and np.array_equal(recomputed, stored["p_U"].to_numpy(dtype=float))
                ),
            }
        )
        checks[-1]["pass"] = bool(checks[-1]["elementwise_identical"])
    return {
        "pass": bool(len(checks) == 2 and all(item["pass"] for item in checks)),
        "task_id": task_id,
        "checks": checks,
    }


def recompute_target_and_macro(frames: dict[str, pd.DataFrame]) -> dict:
    predictions = pd.read_csv(RUN / "rq5rf_pair_predictions.csv", low_memory=False)
    recorded_target = pd.read_csv(RUN / "rq5rf_target_metrics.csv", low_memory=False)
    recorded_macro = pd.read_csv(RUN / "rq5rf_macro_metrics.csv", low_memory=False)
    recorded_pair = pd.read_csv(RUN / "rq5rf_pair_metrics.csv", low_memory=False)
    probability_columns = {
        "U_transport_zero_label": "p_U",
        "predictive_I": "p_I",
        "predictive_F": "p_F",
    }
    predictions = pd.concat(
        [
            predictions.assign(
                model_layer_label=layer,
                selected_probability=predictions[probability_column],
            )
            for layer, probability_column in probability_columns.items()
        ],
        ignore_index=True,
    )
    rows = []
    for keys, group in predictions.groupby(
        [
            "analysis_set",
            "outcome_name",
            "cluster_rule",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
            "target",
        ],
        sort=True,
    ):
        analysis_set, outcome_name, cluster_rule, weighting, exposure, layer, target_name = keys
        group = group.sort_values("pair_id", kind="mergesort")
        adaptation = (
            group["bss_adaptation_train_reference"].to_numpy(dtype=float)
            if layer != "U_transport_zero_label"
            else None
        )
        metrics = own_metrics(
            group["outcome_value"].to_numpy(dtype=float),
            group["selected_probability"].to_numpy(dtype=float),
            group["bss_transport_reference"].to_numpy(dtype=float),
            group["bss_target_oracle_reference"].to_numpy(dtype=float),
            adaptation,
        )
        rows.append(
            {
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "cluster_rule": cluster_rule,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "model_layer_label": layer,
                "target": target_name,
                **metrics,
            }
        )
    independent = pd.DataFrame(rows)
    merge_keys = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "model_exposure",
        "model_layer_label",
        "target",
    ]
    aligned_target = independent.merge(
        recorded_target, on=merge_keys, how="outer", validate="one_to_one", suffixes=("__independent", "__recorded"), indicator=True
    )
    target_differences = {}
    target_pass = True
    for metric in (
        "brier",
        "log_loss",
        "bss_transport",
        "bss_target_oracle",
        "roc_auc",
        "ap",
        "calibration_intercept",
        "calibration_slope",
    ):
        left = pd.to_numeric(aligned_target[f"{metric}__independent"], errors="coerce")
        right = pd.to_numeric(aligned_target[f"{metric}__recorded"], errors="coerce")
        mask = left.notna() & right.notna()
        difference = float(np.max(np.abs(left[mask] - right[mask]))) if mask.any() else float("nan")
        target_differences[metric] = difference
        if not metric.startswith("calibration_"):
            target_pass &= bool(
                aligned_target["_merge"].eq("both").all()
                and mask.all()
                and difference <= 1e-9
            )

    macro_rows = []
    for keys, group in independent.groupby(
        [
            "analysis_set",
            "outcome_name",
            "cluster_rule",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
        ],
        sort=True,
    ):
        analysis_set = keys[0]
        complete = set(group["target"].astype(str)) == expected_targets(analysis_set)
        row = dict(zip(
            [
                "analysis_set",
                "outcome_name",
                "cluster_rule",
                "source_weighting",
                "model_exposure",
                "model_layer_label",
            ],
            keys,
        ))
        for metric in (
            "brier",
            "log_loss",
            "bss_transport",
            "bss_target_oracle",
            "roc_auc",
            "ap",
            "calibration_intercept",
            "calibration_slope",
        ):
            values = pd.to_numeric(group[metric], errors="coerce")
            row[metric] = float(values.mean()) if complete and values.notna().all() else np.nan
        if keys[5] != "U_transport_zero_label":
            values = pd.to_numeric(group["bss_adaptation_train"], errors="coerce")
            row["bss_adaptation_train"] = float(values.mean()) if complete and values.notna().all() else np.nan
        else:
            row["bss_adaptation_train"] = np.nan
        macro_rows.append(row)
    macro_independent = pd.DataFrame(macro_rows)
    macro_aligned = macro_independent.merge(
        recorded_macro,
        on=[
            "analysis_set",
            "outcome_name",
            "cluster_rule",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
        ],
        how="outer",
        validate="one_to_one",
        suffixes=("__independent", "__recorded"),
        indicator=True,
    )
    macro_differences = {}
    macro_pass = True
    for metric in (
        "brier",
        "log_loss",
        "bss_transport",
        "bss_target_oracle",
        "roc_auc",
        "ap",
        "calibration_intercept",
        "calibration_slope",
    ):
        left = pd.to_numeric(macro_aligned[f"{metric}__independent"], errors="coerce")
        right = pd.to_numeric(macro_aligned[f"{metric}__recorded"], errors="coerce")
        mask = left.notna() & right.notna()
        difference = float(np.max(np.abs(left[mask] - right[mask]))) if mask.any() else float("nan")
        macro_differences[metric] = difference
        if not metric.startswith("calibration_"):
            macro_pass &= bool(
                macro_aligned["_merge"].eq("both").all()
                and mask.all()
                and difference <= 1e-9
            )
    pair_rows = []
    for keys, group in predictions.groupby(
        [
            "analysis_set",
            "outcome_name",
            "cluster_rule",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
        ],
        sort=True,
    ):
        analysis_set, outcome_name, cluster_rule, weighting, exposure, layer = keys
        expected_ids = set(frames[analysis_set]["pair_id"].astype(str))
        complete = (
            set(group["pair_id"].astype(str)) == expected_ids
            and not group["pair_id"].astype(str).duplicated().any()
            and len(group) == len(expected_ids)
        )
        if not complete:
            pair_rows.append(
                {
                    "analysis_set": analysis_set,
                    "outcome_name": outcome_name,
                    "cluster_rule": cluster_rule,
                    "source_weighting": weighting,
                    "model_exposure": exposure,
                    "model_layer_label": layer,
                    "complete": False,
                }
            )
            continue
        adaptation = (
            group["bss_adaptation_train_reference"].to_numpy(dtype=float)
            if layer != "U_transport_zero_label"
            else None
        )
        pair_rows.append(
            {
                "analysis_set": analysis_set,
                "outcome_name": outcome_name,
                "cluster_rule": cluster_rule,
                "source_weighting": weighting,
                "model_exposure": exposure,
                "model_layer_label": layer,
                "complete": True,
                **own_metrics(
                    group["outcome_value"].to_numpy(dtype=float),
                    group["selected_probability"].to_numpy(dtype=float),
                    group["bss_transport_reference"].to_numpy(dtype=float),
                    group["bss_target_oracle_reference"].to_numpy(dtype=float),
                    adaptation,
                ),
            }
        )
    independent_pair = pd.DataFrame(pair_rows)
    pair_aligned = independent_pair.merge(
        recorded_pair,
        on=[
            "analysis_set",
            "outcome_name",
            "cluster_rule",
            "source_weighting",
            "model_exposure",
            "model_layer_label",
        ],
        how="outer",
        validate="one_to_one",
        suffixes=("__independent", "__recorded"),
        indicator=True,
    )
    pair_differences = {}
    pair_pass = bool(
        pair_aligned["_merge"].eq("both").all()
        and pair_aligned["complete"].all()
        and pair_aligned["complete_sample_gate"].eq("pass").all()
    )
    for metric in (
        "brier",
        "log_loss",
        "bss_transport",
        "bss_target_oracle",
        "roc_auc",
        "ap",
        "calibration_intercept",
        "calibration_slope",
    ):
        left = pd.to_numeric(pair_aligned[f"{metric}__independent"], errors="coerce")
        right = pd.to_numeric(pair_aligned[f"{metric}__recorded"], errors="coerce")
        mask = left.notna() & right.notna()
        difference = float(np.max(np.abs(left[mask] - right[mask]))) if mask.any() else float("nan")
        pair_differences[metric] = difference
        if not metric.startswith("calibration_"):
            pair_pass &= bool(mask.all() and difference <= 1e-9)
    return {
        "pass": bool(target_pass and macro_pass and pair_pass),
        "target_rows_independent": len(independent),
        "target_rows_recorded": len(recorded_target),
        "target_max_absolute_differences": target_differences,
        "macro_rows_independent": len(macro_independent),
        "macro_rows_recorded": len(recorded_macro),
        "macro_max_absolute_differences": macro_differences,
        "target_recompute_pass": bool(target_pass),
        "macro_recompute_pass": bool(macro_pass),
        "pair_rows_independent": len(independent_pair),
        "pair_rows_recorded": len(recorded_pair),
        "pair_max_absolute_differences": pair_differences,
        "pair_recompute_pass": bool(pair_pass),
        "independent_target_metrics": independent,
    }


def pairwise_and_screen(independent_target: pd.DataFrame) -> dict:
    pairs = [
        ("pooled_logistic_combo", "pooled_logistic_ipTM"),
        ("rf_ipTM", "pooled_logistic_ipTM"),
        ("rf_combo", "pooled_logistic_combo"),
        ("rf_combo", "rf_ipTM"),
        ("rf_combo", "pooled_logistic_ipTM"),
    ]
    rows = []
    for left_model, right_model in pairs:
        comparison = f"{left_model}_minus_{right_model}"
        for layer in ("U_transport_zero_label", "predictive_I", "predictive_F"):
            for keys, group in independent_target[
                independent_target["model_layer_label"].eq(layer)
            ].groupby(
                ["analysis_set", "source_weighting", "target"], sort=True
            ):
                analysis_set, weighting, target_name = keys
                left = group[group["model_exposure"].eq(left_model)]
                right = group[group["model_exposure"].eq(right_model)]
                if len(left) != 1 or len(right) != 1:
                    continue
                for metric in ("brier", "log_loss", "bss_transport"):
                    rows.append(
                        {
                            "analysis_set": analysis_set,
                            "source_weighting": weighting,
                            "target": target_name,
                            "comparison": comparison,
                            "layer": layer,
                            "metric": metric,
                            "difference": float(left.iloc[0][metric]) - float(right.iloc[0][metric]),
                        }
                    )
    differences = pd.DataFrame(rows)
    macro = differences.groupby(
        ["analysis_set", "source_weighting", "comparison", "layer", "metric"],
        sort=True,
    )["difference"].mean().reset_index()

    def value(analysis_set: str, weighting: str, comparison: str, layer: str, metric: str) -> float:
        subset = macro[
            macro["analysis_set"].eq(analysis_set)
            & macro["source_weighting"].eq(weighting)
            & macro["comparison"].eq(comparison)
            & macro["layer"].eq(layer)
            & macro["metric"].eq(metric)
        ]
        return float(subset.iloc[0]["difference"]) if len(subset) == 1 else np.nan

    def improve(metric: str, difference: float) -> bool:
        if not np.isfinite(difference):
            return False
        return difference < -1e-12 if metric in {"brier", "log_loss"} else difference > 1e-12

    def reversal(metric: str, difference: float) -> bool:
        if not np.isfinite(difference):
            return True
        return difference > 1e-12 if metric in {"brier", "log_loss"} else difference < -1e-12

    comparisons = ("rf_combo_minus_pooled_logistic_combo", "rf_combo_minus_rf_ipTM")
    condition_1 = all(
        improve(
            metric,
            value("primary_six", "target_equal", comparison, "U_transport_zero_label", metric),
        )
        for comparison in comparisons
        for metric in ("brier", "log_loss", "bss_transport")
    )
    condition_2 = True
    for analysis_set, required in (("primary_six", 5), ("primary_five", 4)):
        for weighting in ("target_equal", "pair"):
            for comparison in comparisons:
                subset = differences[
                    differences["analysis_set"].eq(analysis_set)
                    & differences["source_weighting"].eq(weighting)
                    & differences["comparison"].eq(comparison)
                    & differences["layer"].eq("U_transport_zero_label")
                    & differences["metric"].eq("brier")
                ]
                condition_2 &= len(subset) == len(expected_targets(analysis_set)) and int((subset["difference"] < -1e-12).sum()) >= required
    condition_3 = True
    for analysis_set in CONFIG["analysis_sets"]:
        for weighting in ("target_equal", "pair"):
            for comparison in comparisons:
                for metric in ("brier", "log_loss", "bss_transport"):
                    condition_3 &= not reversal(
                        metric,
                        value(analysis_set, weighting, comparison, "U_transport_zero_label", metric),
                    )

    recorded = json.loads((RUN / "direction_consistency_10_2.json").read_text(encoding="utf-8"))
    recorded_condition_1 = bool(recorded["condition_1_primary_six_target_equal_strict_improvement"])
    recorded_condition_2 = bool(recorded["condition_2_target_direction_counts"])
    recorded_condition_3 = bool(recorded["condition_3_six_configurations_no_reversal"])
    return {
        "pass": bool(
            condition_1 == recorded_condition_1
            and condition_2 == recorded_condition_2
            and condition_3 == recorded_condition_3
        ),
        "independent_condition_1": bool(condition_1),
        "independent_condition_2": bool(condition_2),
        "independent_condition_3": bool(condition_3),
        "recorded_conditions": [
            recorded_condition_1,
            recorded_condition_2,
            recorded_condition_3,
        ],
        "independent_pairwise_macro": macro,
    }


def main() -> None:
    checks = []
    config = CONFIG
    input_manifest = pd.read_csv(
        RUN / "final_input_hashes.csv", keep_default_na=False
    )
    for _, row in input_manifest.iterrows():
        path = Path(str(row["path"]))
        actual = sha(path) if path.exists() else ""
        expected_move = str(row["expected_sha256"]).strip().lower()
        checks.append(
            {
                "name": f"input_hash:{row['role']}",
                "pass": bool(
                    path.exists()
                    and (
                        actual == expected_move
                        if expected_move
                        else actual == str(row["actual_sha256"]).lower()
                    )
                ),
                "actual": actual,
            }
        )
    frames = analysis_frames()
    for analysis_set, frame in frames.items():
        expected = config["analysis_sets"][analysis_set]
        checks.append(
            {
                "name": f"sample:{analysis_set}",
                "pass": bool(
                    len(frame) == expected["expected_rows"]
                    and int(frame["outcome_value"].sum()) == expected["expected_positives"]
                    and frame["pair_id"].is_unique
                ),
                "rows": len(frame),
                "positives": int(frame["outcome_value"].sum()),
            }
        )
    selection = selection_checks()
    checks.append({"name": "selection_recomputation", **selection})
    rerun = rerun_final_predictions(frames)
    checks.append({"name": "single_and_dual_final_rf_rerun", **rerun})
    metrics = recompute_target_and_macro(frames)
    independent_target = metrics.pop("independent_target_metrics")
    checks.append({"name": "target_and_macro_recomputation", **metrics})
    pairwise = pairwise_and_screen(independent_target)
    independent_pairwise = pairwise.pop("independent_pairwise_macro")
    checks.append({"name": "pairwise_and_10_2_recomputation", **pairwise})
    independent_target.to_csv(RUN / "independent_target_metrics.csv", index=False)
    independent_pairwise.to_csv(RUN / "independent_pairwise_macro.csv", index=False)

    images = sorted(RUN.glob("*.png"))
    image_checks = []
    for path in images:
        raw = path.read_bytes()
        image_checks.append(
            {
                "file": path.name,
                "bytes": len(raw),
                "png_signature": raw[:8] == b"\x89PNG\r\n\x1a\n",
                "pass": len(raw) > 10000 and raw[:8] == b"\x89PNG\r\n\x1a\n",
            }
        )
    checks.append(
        {
            "name": "figure_readability",
            "pass": bool(len(image_checks) >= 4 and all(item["pass"] for item in image_checks)),
            "images": image_checks,
        }
    )
    checkpoint_files = sorted((RUN / "checkpoints").glob("*.json.gz"))
    checks.append(
        {
            "name": "outer_checkpoint_count",
            "pass": len(checkpoint_files) == 34,
            "observed": len(checkpoint_files),
        }
    )
    result = {
        "verification_scope": "independent of production fitting, selection, aggregation, and screening functions",
        "all_passed": all(item["pass"] for item in checks),
        "checks": checks,
        "independent_target_rows": len(independent_target),
        "independent_pairwise_rows": len(independent_pairwise),
    }
    (RUN / "independent_verification.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({key: value for key, value in result.items() if key != "checks"}, indent=2))
    for check in checks:
        print(f"{check['name']}: {'PASS' if check['pass'] else 'FAIL'}")
    raise SystemExit(0 if result["all_passed"] else 1)


if __name__ == "__main__":
    main()
