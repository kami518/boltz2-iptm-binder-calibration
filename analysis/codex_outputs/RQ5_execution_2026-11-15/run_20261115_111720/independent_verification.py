"""Independent RQ5 verification from frozen inputs and saved predictions.

This script deliberately does not import rq5_common or call the production
aggregation functions.  It rebuilds the analysis frame, Stage A estimates,
one complete paired cluster bootstrap, the main ipTM-only U model, and every
pairwise macro/pair-pooled summary with separate metric code.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score

RUN_DIR = Path(__file__).resolve().parent
ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
CLIP = 1e-6
EDGES = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0)
ANALYSIS_SETS = {
    "primary_six": {"outcome": "primary_outcome", "exclude": None},
    "primary_five": {"outcome": "primary_outcome", "exclude": "aah-i"},
    "strict_six": {"outcome": "strict_outcome", "exclude": None},
}
COMPARISON_METRICS = (
    "roc_auc",
    "ap",
    "brier",
    "log_loss",
    "bss_transport",
    "bss_target_oracle",
    "bss_adaptation_train",
    "calibration_intercept",
    "calibration_slope",
)
CHECKS: list[dict] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_nll(eta: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.logaddexp(0.0, eta) - y * eta


def logit(probability: np.ndarray) -> np.ndarray:
    values = np.asarray(probability, dtype=float)
    values = np.clip(values, CLIP, 1.0 - CLIP)
    return np.log(values) - np.log1p(-values)


def auc(y: np.ndarray, score: np.ndarray) -> float:
    if np.unique(y).size < 2 or not np.isfinite(score).all():
        return np.nan
    return float(roc_auc_score(y, score))


def ap(y: np.ndarray, score: np.ndarray) -> float:
    if np.sum(y) <= 0 or not np.isfinite(score).all():
        return np.nan
    return float(average_precision_score(y, score))


def brier(y: np.ndarray, probability: np.ndarray) -> float:
    return float(np.mean((probability - y) ** 2))


def log_loss(y: np.ndarray, probability: np.ndarray) -> float:
    probability = np.clip(probability, CLIP, 1.0 - CLIP)
    return float(-np.mean(y * np.log(probability) + (1.0 - y) * np.log1p(-probability)))


def bss(y: np.ndarray, probability: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.sum((reference - y) ** 2))
    if denominator <= 0 or not math.isfinite(denominator):
        return np.nan
    return float(1.0 - np.sum((probability - y) ** 2) / denominator)


def newton_logistic(
    design: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    x = np.asarray(design, dtype=float)
    yy = np.asarray(y, dtype=float)
    ww = np.ones(len(yy)) if weights is None else np.asarray(weights, dtype=float)
    beta = np.zeros(x.shape[1], dtype=float)
    for column in range(min(x.shape[1], 8)):
        if np.allclose(x[:, column], 1.0):
            prevalence = np.clip(np.average(yy, weights=ww), 1e-6, 1.0 - 1e-6)
            beta[column] = math.log(prevalence / (1.0 - prevalence))
    weight_sum = float(ww.sum())
    for _ in range(200):
        eta = x @ beta
        p = expit(eta)
        gradient = x.T @ (ww * (p - yy)) / weight_sum
        if float(np.max(np.abs(gradient))) < 1e-11:
            return beta
        hessian = (x.T * (ww * p * (1.0 - p))) @ x / weight_sum
        step = np.linalg.solve(hessian, -gradient)
        current = float(np.sum(ww * stable_nll(eta, yy)) / weight_sum)
        scale = 1.0
        accepted = False
        while scale > 1e-12:
            candidate = beta + scale * step
            candidate_objective = float(
                np.sum(ww * stable_nll(x @ candidate, yy)) / weight_sum
            )
            if candidate_objective <= current + 1e-15:
                beta = candidate
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            break
    return beta


def calibration(y: np.ndarray, probability: np.ndarray) -> tuple[float, float]:
    q = logit(probability)
    design = np.column_stack([np.ones(len(q)), q])
    beta = newton_logistic(design, y)
    return float(beta[0]), float(beta[1])


def point_metrics(
    y: np.ndarray,
    probability: np.ndarray,
    *,
    transport_reference: np.ndarray,
    oracle_reference: np.ndarray,
    adaptation_reference: np.ndarray | None = None,
) -> dict:
    result = {
        "n": int(len(y)),
        "positives": int(np.sum(y)),
        "prevalence": float(np.mean(y)),
        "roc_auc": auc(y, probability),
        "ap": ap(y, probability),
        "brier": brier(y, probability),
        "log_loss": log_loss(y, probability),
        "bss_transport": bss(y, probability, transport_reference),
        "bss_target_oracle": bss(y, probability, oracle_reference),
    }
    if np.unique(y).size == 2:
        result["calibration_intercept"], result["calibration_slope"] = calibration(
            y, probability
        )
    else:
        result["calibration_intercept"] = np.nan
        result["calibration_slope"] = np.nan
    if adaptation_reference is not None:
        result["bss_adaptation_train"] = bss(
            y, probability, adaptation_reference
        )
    return result


def record_check(name: str, passed: bool, detail) -> None:
    CHECKS.append({"check": name, "passed": bool(passed), "detail": detail})


def main() -> None:
    global CHECKS
    CHECKS = []
    with (RUN_DIR / "run_config.snapshot.json").open("r", encoding="utf-8") as handle:
        config = json.load(handle)

    hash_table = pd.read_csv(RUN_DIR / "final_input_hashes.csv")
    hash_mismatches = []
    for row in hash_table.itertuples(index=False):
        actual = sha256(Path(row.path))
        if actual.lower() != str(row.expected_sha256).lower() or actual.lower() != str(row.actual_sha256).lower():
            hash_mismatches.append(row.input_key)
    record_check(
        "frozen_input_hashes",
        not hash_mismatches,
        {"mismatches": hash_mismatches, "n_checked": len(hash_table)},
    )

    master = pd.read_csv(config["inputs"]["pair_master"]["path"], low_memory=False)
    pinc = pd.read_csv(config["inputs"]["pinc_extension"]["path"], low_memory=False)
    clusters = pd.read_csv(
        config["inputs"]["sequence_clusters"]["path"], low_memory=False
    )
    record_check(
        "core_input_rows",
        len(master) == 2084 and len(pinc) == 2084 and len(clusters) == 2084,
        {"master": len(master), "pinc": len(pinc), "clusters": len(clusters)},
    )
    record_check(
        "core_unique_pair_ids",
        all(frame["pair_id"].is_unique for frame in (master, pinc, clusters)),
        {
            "master": int(master["pair_id"].nunique()),
            "pinc": int(pinc["pair_id"].nunique()),
            "clusters": int(clusters["pair_id"].nunique()),
        },
    )
    joined = master.merge(
        pinc[["pair_id", "target", "pinc", "pinc_status"]],
        on="pair_id",
        validate="one_to_one",
        suffixes=("", "_pinc"),
    )
    if not joined["target"].astype(str).eq(joined["target_pinc"].astype(str)).all():
        raise RuntimeError("Independent target join disagreement")
    joined = joined.drop(columns=["target_pinc"]).merge(
        clusters[["pair_id", "sequence_cluster90"]],
        on="pair_id",
        validate="one_to_one",
    )
    joined["x_ipTM"] = logit(joined["boltz2_iptm"].to_numpy(dtype=float))
    missing = joined[joined["pinc_status"].ne("success")]
    record_check(
        "pinc_acceptance",
        len(missing) == 10
        and set(missing["target"].astype(str)) == {"pd-l1"}
        and int(missing["primary_outcome"].sum()) == 8
        and int(missing["strict_outcome"].sum()) == 4,
        {
            "success": int(joined["pinc_status"].eq("success").sum()),
            "missing": int(len(missing)),
            "missing_primary_positive": int(missing["primary_outcome"].sum()),
            "missing_strict_positive": int(missing["strict_outcome"].sum()),
        },
    )

    frames = {}
    for set_name, spec in ANALYSIS_SETS.items():
        frame = joined[joined["pinc_status"].eq("success")].copy()
        if spec["exclude"] is not None:
            frame = frame[frame["target"].astype(str).ne(spec["exclude"])]
        frame = frame[frame[spec["outcome"]].notna()].copy()
        frame["outcome_value"] = frame[spec["outcome"]].astype(int)
        frames[set_name] = frame
    record_check(
        "analysis_set_counts",
        len(frames["primary_six"]) == 2074
        and int(frames["primary_six"]["outcome_value"].sum()) == 740
        and len(frames["primary_five"]) == 1941
        and int(frames["strict_six"]["outcome_value"].sum()) == 674,
        {
            name: {
                "n": int(len(frame)),
                "positives": int(frame["outcome_value"].sum()),
            }
            for name, frame in frames.items()
        },
    )

    stage_a = pd.read_csv(RUN_DIR / "rq5_target_score_metrics.csv")
    stage_a_errors = []
    for row in stage_a.itertuples(index=False):
        frame = frames[row.analysis_set]
        target = frame[frame["target"].astype(str).eq(str(row.target))].sort_values(
            "pair_id"
        )
        y = target["outcome_value"].to_numpy(dtype=float)
        score = target["boltz2_iptm" if row.score == "ipTM" else "pinc"].to_numpy(
            dtype=float
        )
        recalculated = {
            "n": len(y),
            "positives": int(np.sum(y)),
            "roc_auc": auc(y, score),
            "ap": ap(y, score),
            "brier": brier(y, score),
            "log_loss": log_loss(y, score),
            "bss_target_constant": bss(
                y, score, np.full(len(y), float(np.mean(y)))
            ),
        }
        cal_intercept, cal_slope = calibration(y, score)
        recalculated["calibration_intercept"] = cal_intercept
        recalculated["calibration_slope"] = cal_slope
        for metric, value in recalculated.items():
            observed = getattr(row, metric)
            if not (
                (pd.isna(value) and pd.isna(observed))
                or np.isclose(value, observed, atol=2e-7, rtol=1e-8)
            ):
                stage_a_errors.append(
                    {
                        "analysis_set": row.analysis_set,
                        "target": row.target,
                        "score": row.score,
                        "metric": metric,
                        "expected_saved": observed,
                        "independently_recomputed": value,
                    }
                )
    record_check(
        "stage_a_all_point_estimates",
        not stage_a_errors,
        {"n_rows": len(stage_a), "errors": stage_a_errors[:10]},
    )

    # Reproduce a complete paired cluster bootstrap for one representative item.
    representative = stage_a[
        stage_a["analysis_set"].eq("primary_six")
        & stage_a["target"].eq("aah-i")
    ]["score"].tolist()
    bootstrap_saved = pd.read_csv(RUN_DIR / "rq5_score_bootstrap_comparisons.csv")
    bootstrap_saved = bootstrap_saved[
        bootstrap_saved["analysis_set"].eq("primary_six")
        & bootstrap_saved["target"].eq("aah-i")
    ]
    boot_frame = frames["primary_six"]
    boot_frame = boot_frame[boot_frame["target"].astype(str).eq("aah-i")].sort_values(
        "pair_id"
    )
    cluster_values = sorted(boot_frame["sequence_cluster90"].astype(str).unique())
    cluster_indices = [
        np.flatnonzero(
            boot_frame["sequence_cluster90"].astype(str).to_numpy() == cluster
        )
        for cluster in cluster_values
    ]
    rng = np.random.default_rng(20261115)
    differences = {"auc": [], "ap": [], "brier": []}
    y_boot = boot_frame["outcome_value"].to_numpy(dtype=float)
    score_boot = {
        "ipTM": boot_frame["boltz2_iptm"].to_numpy(dtype=float),
        "Pinc": boot_frame["pinc"].to_numpy(dtype=float),
    }
    for _ in range(2000):
        sampled = rng.choice(len(cluster_values), size=len(cluster_values), replace=True)
        selected = np.concatenate([cluster_indices[index] for index in sampled])
        yy = y_boot[selected]
        values = {}
        for score_label in ("ipTM", "Pinc"):
            ss = score_boot[score_label][selected]
            values[score_label] = {
                "auc": auc(yy, ss),
                "ap": ap(yy, ss),
                "brier": brier(yy, ss),
            }
        for metric in differences:
            left = values["ipTM"][metric]
            right = values["Pinc"][metric]
            if np.isfinite(left) and np.isfinite(right):
                differences[metric].append(float(right - left))
    bootstrap_errors = []
    for metric, values in differences.items():
        saved = bootstrap_saved[bootstrap_saved["metric"].eq(metric)].iloc[0]
        low, high = np.percentile(values, [2.5, 97.5])
        recalculated = {
            "bootstrap_mean_difference": float(np.mean(values)),
            "bootstrap_95_low": float(low),
            "bootstrap_95_high": float(high),
            "n_valid_boot": len(values),
        }
        for key, value in recalculated.items():
            if not np.isclose(value, float(saved[key]), atol=1e-12, rtol=1e-12):
                bootstrap_errors.append(
                    {
                        "metric": metric,
                        "field": key,
                        "saved": float(saved[key]),
                        "recomputed": value,
                    }
                )
    record_check(
        "representative_paired_cluster_bootstrap",
        not bootstrap_errors,
        {
            "analysis_set": "primary_six",
            "target": "aah-i",
            "n_boot": 2000,
            "seed": 20261115,
            "errors": bootstrap_errors,
        },
    )

    loto = pd.read_csv(config["inputs"]["loto_assignments"]["path"], low_memory=False)
    loto = loto[
        loto["analysis_set"].isin(ANALYSIS_SETS)
        & loto["cluster_rule"].eq("sequence_cluster90")
    ]
    fold_errors = []
    for fold_id, fold in loto.groupby("fold_id"):
        train = fold[fold["split_status"].eq("train")]
        test = fold[fold["split_status"].eq("test")]
        held_target = str(test["held_out_target"].iloc[0])
        if train["target"].astype(str).eq(held_target).any():
            fold_errors.append((fold_id, "held_target_in_train"))
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            fold_errors.append((fold_id, "cluster_overlap"))
    record_check("frozen_loto_isolation", not fold_errors, fold_errors[:10])

    # Independently reconstruct primary_six/target_equal/ipTM_only U.
    u_saved = pd.read_csv(RUN_DIR / "rq5_u_pair_predictions.csv", low_memory=False)
    main_u = u_saved[
        u_saved["analysis_set"].eq("primary_six")
        & u_saved["source_weighting"].eq("target_equal")
        & u_saved["model_exposure"].eq("ipTM_only")
    ].copy()
    joined_index = joined.set_index("pair_id")
    prediction_errors = []
    for held_target in sorted(main_u["held_out_target"].unique()):
        fold_id = f"loto|primary_six|sequence_cluster90|{held_target}"
        fold = loto[loto["fold_id"].eq(fold_id)]
        train = fold[fold["split_status"].eq("train")].merge(
            joined_index[["boltz2_iptm", "x_ipTM", "pinc_status"]],
            left_on="pair_id",
            right_index=True,
            how="left",
        )
        test = fold[fold["split_status"].eq("test")].merge(
            joined_index[["boltz2_iptm", "x_ipTM", "pinc_status"]],
            left_on="pair_id",
            right_index=True,
            how="left",
        )
        train = train[train["pinc_status"].eq("success")]
        test = test[test["pinc_status"].eq("success")].sort_values("pair_id")
        targets = sorted(train["target"].astype(str).unique())
        design = np.column_stack(
            [
                (train["target"].astype(str).to_numpy()[:, None] == np.array(targets)[None, :]).astype(float),
                train["x_ipTM"].to_numpy(dtype=float),
            ]
        )
        beta = newton_logistic(
            design,
            train["outcome_value"].to_numpy(dtype=float),
            train["weight_target_equal"].to_numpy(dtype=float),
        )
        transport = float(np.mean(beta[: len(targets)]))
        eta = transport + beta[-1] * test["x_ipTM"].to_numpy(dtype=float)
        saved = main_u[main_u["held_out_target"].eq(held_target)].sort_values("pair_id")
        if not np.array_equal(saved["pair_id"].astype(str).to_numpy(), test["pair_id"].astype(str).to_numpy()):
            prediction_errors.append((held_target, "row_order"))
        max_eta_error = float(np.max(np.abs(saved["eta_U"].to_numpy(dtype=float) - eta)))
        if max_eta_error > 5e-6:
            prediction_errors.append((held_target, "eta", max_eta_error))
    record_check(
        "main_U_prediction_reconstruction",
        not prediction_errors,
        {"errors": prediction_errors[:10], "analysis": "primary_six/target_equal/ipTM_only"},
    )

    # Recompute target-level U metrics from reconstructed saved predictions.
    u_target = pd.read_csv(RUN_DIR / "rq5_u_target_metrics.csv")
    main_u_target = u_target[
        u_target["analysis_set"].eq("primary_six")
        & u_target["source_weighting"].eq("target_equal")
        & u_target["model_exposure"].eq("ipTM_only")
    ]
    metric_errors = []
    for held_target in sorted(main_u["held_out_target"].unique()):
        group = main_u[main_u["held_out_target"].eq(held_target)].sort_values("pair_id")
        recalculated = point_metrics(
            group["outcome_value"].to_numpy(dtype=float),
            group["p_U"].to_numpy(dtype=float),
            transport_reference=group["bss_transport_reference"].to_numpy(dtype=float),
            oracle_reference=group["bss_target_oracle_reference"].to_numpy(dtype=float),
        )
        saved = main_u_target[main_u_target["held_out_target"].eq(held_target)].iloc[0]
        for metric, value in recalculated.items():
            if not np.isclose(value, float(saved[metric]), atol=1e-5, rtol=1e-8):
                metric_errors.append(
                    {
                        "target": held_target,
                        "metric": metric,
                        "saved": float(saved[metric]),
                        "recomputed": value,
                    }
                )
    record_check(
        "main_U_target_metrics_reconstruction",
        not metric_errors,
        {"errors": metric_errors[:10]},
    )

    # All target-level pairwise differences and macro summaries.
    adaptation = pd.read_csv(
        RUN_DIR / "rq5_adaptation_pair_predictions.csv", low_memory=False
    )
    adaptation_target = pd.read_csv(RUN_DIR / "rq5_adaptation_target_metrics.csv")
    saved_target_differences = pd.read_csv(
        RUN_DIR / "rq5_target_pairwise_differences.csv"
    )
    saved_pairwise = pd.read_csv(RUN_DIR / "rq5_pairwise_comparisons.csv")
    grouping = ["analysis_set", "outcome_name", "cluster_rule", "source_weighting"]
    independent_target_differences = []
    for family in ("preset_U_exposure", "preset_adaptation_increment"):
        if family == "preset_U_exposure":
            comparisons = [
                ("Pinc_only", "ipTM_only"),
                ("ipTM_plus_Pinc", "ipTM_only"),
                ("ipTM_plus_Pinc", "Pinc_only"),
            ]
            layer = "U_transport_zero_label"
            for left_model, right_model in comparisons:
                for keys, left_group in u_target[
                    u_target["model_exposure"].eq(left_model)
                ].groupby(grouping + ["held_out_target"], sort=True):
                    *base_keys, target = keys
                    base = dict(zip(grouping + ["held_out_target"], keys))
                    right_group = u_target[
                        u_target["analysis_set"].eq(base["analysis_set"])
                        & u_target["outcome_name"].eq(base["outcome_name"])
                        & u_target["cluster_rule"].eq(base["cluster_rule"])
                        & u_target["source_weighting"].eq(base["source_weighting"])
                        & u_target["model_exposure"].eq(right_model)
                        & u_target["held_out_target"].eq(target)
                    ]
                    if len(left_group) != 1 or len(right_group) != 1:
                        continue
                    left_row = left_group.iloc[0]
                    right_row = right_group.iloc[0]
                    for metric in COMPARISON_METRICS:
                        independent_target_differences.append(
                            {
                                **base,
                                "comparison_family": family,
                                "comparison": f"{left_model}_minus_{right_model}",
                                "layer": layer,
                                "left_model": left_model,
                                "right_model": right_model,
                                "metric": metric,
                                "difference_left_minus_right": (
                                    float(left_row[metric]) - float(right_row[metric])
                                    if np.isfinite(left_row[metric]) and np.isfinite(right_row[metric])
                                    else np.nan
                                ),
                            }
                        )
        else:
            for exposure in ("ipTM_only", "Pinc_only", "ipTM_plus_Pinc"):
                subset = adaptation_target[
                    adaptation_target["model_exposure"].eq(exposure)
                ]
                for left_layer, right_layer in (
                    ("predictive_I", "U_transport_zero_label"),
                    ("predictive_F", "U_transport_zero_label"),
                    ("predictive_F", "predictive_I"),
                ):
                    for keys, left_group in subset[
                        subset["model_layer_label"].eq(left_layer)
                    ].groupby(grouping + ["held_out_target"], sort=True):
                        base = dict(zip(grouping + ["held_out_target"], keys))
                        right_group = subset[
                            subset["analysis_set"].eq(base["analysis_set"])
                            & subset["outcome_name"].eq(base["outcome_name"])
                            & subset["cluster_rule"].eq(base["cluster_rule"])
                            & subset["source_weighting"].eq(base["source_weighting"])
                            & subset["model_layer_label"].eq(right_layer)
                            & subset["held_out_target"].eq(base["held_out_target"])
                        ]
                        if len(left_group) != 1 or len(right_group) != 1:
                            continue
                        left_row = left_group.iloc[0]
                        right_row = right_group.iloc[0]
                        for metric in COMPARISON_METRICS:
                            independent_target_differences.append(
                                {
                                    **base,
                                    "comparison_family": family,
                                    "comparison": f"{left_layer}_minus_{right_layer}",
                                    "layer": f"{left_layer}|{right_layer}",
                                    "left_model": exposure,
                                    "right_model": exposure,
                                    "metric": metric,
                                    "difference_left_minus_right": (
                                        float(left_row[metric]) - float(right_row[metric])
                                        if np.isfinite(left_row[metric]) and np.isfinite(right_row[metric])
                                        else np.nan
                                    ),
                                }
                            )
    independent_target_frame = pd.DataFrame(independent_target_differences)
    saved_keys = [
        "analysis_set",
        "outcome_name",
        "cluster_rule",
        "source_weighting",
        "comparison_family",
        "comparison",
        "layer",
        "left_model",
        "right_model",
        "held_out_target",
        "metric",
    ]
    compare_target = saved_target_differences[saved_keys + ["difference_left_minus_right"]].merge(
        independent_target_frame[saved_keys + ["difference_left_minus_right"]],
        on=saved_keys,
        how="outer",
        suffixes=("__saved", "__ind"),
        indicator=True,
    )
    target_diff_errors = compare_target[
        (compare_target["_merge"].ne("both"))
        | (
            compare_target["difference_left_minus_right__saved"].notna()
            ^ compare_target["difference_left_minus_right__ind"].notna()
        )
        | (
            np.abs(
                compare_target["difference_left_minus_right__saved"]
                - compare_target["difference_left_minus_right__ind"]
            )
            > 2e-8
        )
    ]
    record_check(
        "all_target_pairwise_differences",
        len(target_diff_errors) == 0,
        {
            "saved_rows": len(saved_target_differences),
            "independent_rows": len(independent_target_frame),
            "mismatch_rows": len(target_diff_errors),
        },
    )

    macro_keys = grouping + [
        "comparison_family",
        "comparison",
        "layer",
        "left_model",
        "right_model",
        "metric",
    ]
    independent_macro = (
        independent_target_frame.groupby(macro_keys, dropna=False)[
            "difference_left_minus_right"
        ]
        .mean()
        .reset_index()
    )
    saved_macro = saved_pairwise[
        saved_pairwise["summary_type"].eq("target_equal_macro")
    ][macro_keys + ["mean_difference"]]
    compare_macro = saved_macro.merge(
        independent_macro,
        on=macro_keys,
        how="outer",
        suffixes=("__saved", "__ind"),
        indicator=True,
    )
    macro_errors = compare_macro[
        (compare_macro["_merge"].ne("both"))
        | (
            np.abs(
                (
                    compare_macro["mean_difference__saved"]
                    if "mean_difference__saved" in compare_macro.columns
                    else compare_macro["mean_difference"]
                )
                - compare_macro["difference_left_minus_right"]
            )
            > 2e-8
        )
    ]
    record_check(
        "all_pairwise_macro_summaries",
        len(macro_errors) == 0,
        {
            "saved_rows": len(saved_macro),
            "independent_rows": len(independent_macro),
            "mismatch_rows": len(macro_errors),
        },
    )

    # Independently recompute every pair-pooled pairwise difference.
    pair_rows = []
    pair_saved = saved_pairwise[saved_pairwise["summary_type"].eq("pair_pooled")]
    for saved in pair_saved.itertuples(index=False):
        base_mask = (
            adaptation["analysis_set"].eq(saved.analysis_set)
            & adaptation["source_weighting"].eq(saved.source_weighting)
            & adaptation["model_exposure"].eq(saved.left_model)
        )
        if saved.comparison_family == "preset_U_exposure":
            left_frame = u_saved[
                u_saved["analysis_set"].eq(saved.analysis_set)
                & u_saved["source_weighting"].eq(saved.source_weighting)
                & u_saved["model_exposure"].eq(saved.left_model)
            ]
            right_frame = u_saved[
                u_saved["analysis_set"].eq(saved.analysis_set)
                & u_saved["source_weighting"].eq(saved.source_weighting)
                & u_saved["model_exposure"].eq(saved.right_model)
            ]
            left_probability = "p_U"
            right_probability = "p_U"
        else:
            left_frame = adaptation[base_mask]
            right_frame = left_frame
            left_probability = {
                "predictive_I": "p_I",
                "predictive_F": "p_F",
            }[saved.comparison.split("_minus_")[0]]
            right_probability = {
                "U_transport_zero_label": "p_U",
                "predictive_I": "p_I",
            }[saved.comparison.split("_minus_")[1]]
        aligned = left_frame.merge(
            right_frame,
            on="pair_id",
            how="inner",
            suffixes=("__left", "__right"),
            validate="one_to_one",
        )
        if saved.comparison_family == "preset_U_exposure":
            left_prob = aligned["p_U__left"].to_numpy(dtype=float)
            right_prob = aligned["p_U__right"].to_numpy(dtype=float)
            transport = aligned["bss_transport_reference__left"].to_numpy(dtype=float)
            oracle = aligned["bss_target_oracle_reference__left"].to_numpy(dtype=float)
            adaptation_reference = None
        else:
            left_prob = aligned[f"{left_probability}__left"].to_numpy(dtype=float)
            right_prob = aligned[f"{right_probability}__right"].to_numpy(dtype=float)
            transport = aligned["bss_transport_reference__left"].to_numpy(dtype=float)
            oracle = aligned["bss_target_oracle_reference__left"].to_numpy(dtype=float)
            adaptation_reference = aligned[
                "bss_adaptation_train_reference__left"
            ].to_numpy(dtype=float)
        y = aligned["outcome_value__left"].to_numpy(dtype=float)
        left_metrics = point_metrics(
            y,
            left_prob,
            transport_reference=transport,
            oracle_reference=oracle,
            adaptation_reference=adaptation_reference,
        )
        right_metrics = point_metrics(
            y,
            right_prob,
            transport_reference=transport,
            oracle_reference=oracle,
            adaptation_reference=adaptation_reference,
        )
        for metric in COMPARISON_METRICS:
            pair_rows.append(
                {
                    **{
                        key: getattr(saved, key)
                        for key in macro_keys
                    },
                    "summary_type": "pair_pooled",
                    "metric": metric,
                    "independent_difference": (
                        left_metrics.get(metric, np.nan)
                        - right_metrics.get(metric, np.nan)
                    ),
                }
            )
    pair_independent = pd.DataFrame(pair_rows)
    saved_pair_metric = pair_saved.melt(
        id_vars=macro_keys + ["summary_type"],
        value_vars=[f"difference_{metric}" for metric in COMPARISON_METRICS],
        var_name="metric_column",
        value_name="saved_difference",
    )
    saved_pair_metric["metric"] = saved_pair_metric["metric_column"].str.replace(
        "difference_", "", regex=False
    )
    pair_base_keys = [key for key in macro_keys if key != "metric"]
    compare_pair = saved_pair_metric[
        pair_base_keys + ["summary_type", "metric", "saved_difference"]
    ].merge(
        pair_independent[
            pair_base_keys + ["summary_type", "metric", "independent_difference"]
        ],
        on=pair_base_keys + ["summary_type", "metric"],
        how="outer",
        indicator=True,
    )
    pair_errors = compare_pair[
        (compare_pair["_merge"].ne("both"))
        | (
            compare_pair["saved_difference"].notna()
            ^ compare_pair["independent_difference"].notna()
        )
        | (
            np.abs(
                compare_pair["saved_difference"]
                - compare_pair["independent_difference"]
            )
            > 3e-7
        )
    ]
    record_check(
        "all_pair_pooled_pairwise_summaries",
        len(pair_errors) == 0,
        {
            "saved_rows": len(saved_pair_metric),
            "independent_rows": len(pair_independent),
            "mismatch_rows": len(pair_errors),
        },
    )

    # Row counts and equal-row checks.
    expected_rows = sum(len(frame) for frame in frames.values()) * 2 * 3
    record_check(
        "prediction_row_counts",
        len(u_saved) == expected_rows and len(adaptation) == expected_rows,
        {
            "u": len(u_saved),
            "adaptation": len(adaptation),
            "expected": expected_rows,
        },
    )
    u_row_failures = []
    for keys, group in u_saved.groupby(
        ["analysis_set", "source_weighting", "outer_fold_id"]
    ):
        sets = group.groupby("model_exposure")["pair_id"].apply(
            lambda values: frozenset(values.astype(str))
        )
        if sets.nunique() != 1:
            u_row_failures.append(keys)
    record_check("same_U_outer_rows", not u_row_failures, u_row_failures[:10])
    adaptation_row_failures = []
    for keys, group in adaptation.groupby(
        ["analysis_set", "source_weighting", "held_out_target", "outer_fold_id"]
    ):
        sets = group.groupby("model_exposure")["pair_id"].apply(
            lambda values: frozenset(values.astype(str))
        )
        if sets.nunique() != 1:
            adaptation_row_failures.append(keys)
    record_check(
        "same_adaptation_outer_rows",
        not adaptation_row_failures,
        adaptation_row_failures[:10],
    )
    inner = pd.read_csv(
        RUN_DIR / "rq5_adaptation_inner_assignments.csv", low_memory=False
    )
    inner_leaks = []
    for keys, group in inner.groupby(
        [
            "outer_fold_id",
            "source_weighting",
            "model_exposure",
            "inner_fold_number",
        ]
    ):
        train = group[group["split"].eq("train")]
        test = group[group["split"].eq("test")]
        if set(train["sequence_cluster_value"].astype(str)) & set(
            test["sequence_cluster_value"].astype(str)
        ):
            inner_leaks.append(keys)
    record_check("adaptation_inner_isolation", not inner_leaks, inner_leaks[:10])

    image_results = {}
    image_errors = []
    for name in (
        "RQ5_score_discrimination.png",
        "RQ5_raw_reliability.png",
        "RQ5_U_transport.png",
        "RQ5_adaptation_increment.png",
    ):
        path = RUN_DIR / name
        image = mpimg.imread(path)
        finite = bool(np.isfinite(image).all())
        readable = (
            path.stat().st_size > 10000
            and image.ndim >= 2
            and min(image.shape[:2]) >= 300
            and finite
        )
        image_results[name] = {
            "shape": list(image.shape),
            "bytes": path.stat().st_size,
            "finite_pixels": finite,
            "readable": readable,
        }
        if not readable:
            image_errors.append(name)
    record_check("figure_readability", not image_errors, image_results)

    result = {
        "status": "passed" if all(check["passed"] for check in CHECKS) else "failed",
        "independent_of_production_summary_functions": True,
        "n_checks": len(CHECKS),
        "n_failed": int(sum(not check["passed"] for check in CHECKS)),
        "checks": CHECKS,
    }
    def json_safe(value):
        if isinstance(value, (np.generic,)):
            value = value.item()
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {str(key): json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(item) for item in value]
        return value

    result = json_safe(result)
    with (RUN_DIR / "independent_verification.json").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps(result, ensure_ascii=False, indent=2)[:12000])
    if result["status"] != "passed":
        raise RuntimeError("Independent verification failed")


if __name__ == "__main__":
    main()
