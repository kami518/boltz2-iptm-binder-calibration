"""C01B: minimally specified ipTM + ipSAE logistic complementarity test."""
from __future__ import annotations

import importlib.util
import json
import os
import platform
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis" / "postreview_revision_2026-11-21"
COMBO_OUT = OUT / "ipSAE_combo"
BASE_SCRIPT = OUT / "scripts" / "run_c01_c02_postreview.py"
BASE_PAIR_EXTENSION = OUT / "ipSAE" / "ipSAE_pair_extension.csv"
BASELINE_TARGET_METRICS = OUT / "ipSAE" / "ipSAE_target_metrics.csv"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = load_module("postreview_base", BASE_SCRIPT)


def combo_design(frame: pd.DataFrame, levels: list[str]) -> np.ndarray:
    return np.column_stack(
        [
            base.rq3.target_one_hot(frame["target"].astype(str).to_numpy(), levels),
            frame["x_ipTM"].to_numpy(dtype=float),
            frame["x_ipSAE"].to_numpy(dtype=float),
        ]
    )


def fit_combo_source(train: pd.DataFrame, weighting: str, label: str) -> dict:
    levels = sorted(train["target"].astype(str).unique())
    design = combo_design(train, levels)
    outcome = train["outcome_value"].to_numpy(dtype=float)
    weights = base.source_weights(train, weighting)
    prevalence_start = np.zeros(design.shape[1])
    for index, target in enumerate(levels):
        mask = train["target"].astype(str).eq(target).to_numpy()
        prevalence = np.clip(
            float(np.average(outcome[mask], weights=weights[mask])),
            1e-6,
            1.0 - 1e-6,
        )
        prevalence_start[index] = np.log(prevalence / (1.0 - prevalence))
    starts = []
    for slope_ipTM, slope_ipSAE in (
        (0.0, 0.0),
        (0.25, 0.0),
        (-0.25, 0.0),
        (0.0, 0.25),
        (0.0, -0.25),
        (0.25, 0.25),
        (0.25, -0.25),
        (-0.25, 0.25),
        (-0.25, -0.25),
    ):
        start = prevalence_start.copy()
        start[-2] = slope_ipTM
        start[-1] = slope_ipSAE
        starts.append(start)
    fit = base.rq3.fit_weighted_logistic(
        design,
        outcome,
        weights,
        label=label,
        starts=starts,
    )
    fit["levels"] = levels
    if fit.get("estimable", False):
        intercepts = fit["theta"][: len(levels)]
        fit["transport_intercept"] = base.transport_intercept(
            levels,
            intercepts,
            train["target"].astype(str).to_numpy(),
            weighting,
        )
        fit["beta_ipTM"] = float(fit["theta"][-2])
        fit["beta_ipSAE"] = float(fit["theta"][-1])
    else:
        fit["transport_intercept"] = np.nan
        fit["beta_ipTM"] = np.nan
        fit["beta_ipSAE"] = np.nan
    return fit


def build_u_predictions(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    loto = pd.read_csv(base.LOTO, low_memory=False)
    complete_ids = set(data["pair_id"].astype(str))
    fold_ids = (
        loto[
            loto["analysis_set"].eq("primary_six")
            & loto["cluster_rule"].eq("sequence_cluster90")
        ]["fold_id"]
        .drop_duplicates()
        .tolist()
    )
    if len(fold_ids) != 6:
        raise RuntimeError("C01B requires exactly six frozen LOTO folds")
    prediction_frames = []
    parameter_rows = []
    columns = [
        "pair_id",
        "boltz2_iptm",
        "official_ipsae_min",
        "x_ipTM",
        "x_ipSAE",
    ]
    for fold_id in fold_ids:
        fold = loto[loto["fold_id"].eq(fold_id)].copy()
        fold = fold[fold["pair_id"].astype(str).isin(complete_ids)]
        train = fold[fold["split_status"].eq("train")].merge(
            data[columns],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        test = fold[fold["split_status"].eq("test")].merge(
            data[columns + ["sequence_cluster90"]],
            on="pair_id",
            how="left",
            validate="one_to_one",
        )
        if train[["x_ipTM", "x_ipSAE"]].isna().any().any():
            raise RuntimeError("Missing combo source feature")
        if test[["x_ipTM", "x_ipSAE"]].isna().any().any():
            raise RuntimeError("Missing combo held-out feature")
        held_target = str(test["target"].iloc[0])
        ref, _ = base.source_reference(train, "target_equal")
        ref_pair, _ = base.source_reference(train, "pair")
        for weighting in base.WEIGHTINGS:
            fit = fit_combo_source(
                train,
                weighting,
                f"{fold_id}|{weighting}|ipTM_plus_ipSAE",
            )
            if not fit.get("estimable", False):
                raise RuntimeError(
                    f"Combination U fit failed for {fold_id} {weighting}: "
                    f"{fit.get('failure_type')} {fit.get('message')}"
                )
            x_test = test[["x_ipTM", "x_ipSAE"]].to_numpy(dtype=float)
            eta = (
                fit["transport_intercept"]
                + x_test @ np.array([fit["beta_ipTM"], fit["beta_ipSAE"]])
            )
            prediction_frames.append(
                pd.DataFrame(
                    {
                        "pair_id": test["pair_id"].astype(str),
                        "target": test["target"].astype(str),
                        "sequence_cluster_value": test[
                            "sequence_cluster90"
                        ].astype(str),
                        "outer_fold_id": fold_id,
                        "model": "ipTM_plus_ipSAE",
                        "source_weighting": weighting,
                        "outcome_value": test["outcome_value"].astype(int),
                        "raw_score": np.nan,
                        "eta_U": eta,
                        "p_U": expit(eta),
                        "bss_transport_reference": (
                            ref if weighting == "target_equal" else ref_pair
                        ),
                        "bss_target_oracle_reference": float(
                            test["outcome_value"].mean()
                        ),
                    }
                )
            )
            levels = list(fit["levels"])
            intercepts = fit["theta"][: len(levels)]
            parameter_rows.append(
                {
                    "fold_id": fold_id,
                    "held_out_target": held_target,
                    "source_weighting": weighting,
                    "fit_status": fit["status"],
                    "failure_type": fit.get("failure_type", ""),
                    "transport_intercept": fit["transport_intercept"],
                    "beta_ipTM": fit["beta_ipTM"],
                    "beta_ipSAE": fit["beta_ipSAE"],
                    "source_intercepts_json": json.dumps(
                        {
                            level: float(value)
                            for level, value in zip(levels, intercepts)
                        },
                        sort_keys=True,
                    ),
                    "gradient_max_abs": fit.get("gradient_max_abs"),
                    "source_n": len(train),
                    "source_targets": len(levels),
                }
            )
    return pd.concat(prediction_frames, ignore_index=True), pd.DataFrame(
        parameter_rows
    )


def build_combo_metrics(u_predictions: pd.DataFrame):
    cv5 = pd.read_csv(base.CV5, low_memory=False)
    adapted_frames = []
    adaptation_parameters = []
    lambda_frames = []
    for (weighting, target), target_u in u_predictions.groupby(
        ["source_weighting", "target"], sort=True
    ):
        cv_target = cv5[
            cv5["analysis_set"].eq("primary_six")
            & cv5["outcome_name"].eq("primary_outcome")
            & cv5["cluster_rule"].eq("sequence_cluster90")
            & cv5["target"].astype(str).eq(str(target))
        ].copy()
        adapted, parameters, lambdas = base.fit_adaptation(
            target_u,
            cv_target,
            score_model="ipTM_plus_ipSAE",
            weighting=weighting,
        )
        adapted_frames.append(adapted)
        adaptation_parameters.append(parameters)
        lambda_frames.append(lambdas)
    adapted = pd.concat(adapted_frames, ignore_index=True)
    if len(adapted) != 2074 * 2:
        raise RuntimeError("Combination adaptation row count mismatch")
    target_rows = []
    for (weighting, target), group in adapted.groupby(
        ["source_weighting", "target"], sort=True
    ):
        y = group["outcome_value"].to_numpy(dtype=float)
        for layer in base.LAYERS:
            target_rows.append(
                {
                    "model": "ipTM_plus_ipSAE",
                    "source_weighting": weighting,
                    "target": target,
                    "layer": layer,
                    **base.metric_set(
                        y,
                        group[f"p_{layer}"].to_numpy(dtype=float),
                        transport_reference=group[
                            "bss_transport_reference"
                        ].to_numpy(dtype=float),
                        target_reference=group[
                            "bss_target_oracle_reference"
                        ].to_numpy(dtype=float),
                    ),
                    "bss_adaptation_train": (
                        base.bss(
                            y,
                            group[f"p_{layer}"].to_numpy(dtype=float),
                            group["adaptation_train_reference"].to_numpy(dtype=float),
                        )
                        if layer in {"I", "F"}
                        else np.nan
                    ),
                }
            )
    return (
        adapted,
        pd.concat(adaptation_parameters, ignore_index=True),
        pd.concat(lambda_frames, ignore_index=True),
        pd.DataFrame(target_rows),
    )


def build_comparisons(target_metrics: pd.DataFrame) -> pd.DataFrame:
    baseline = pd.read_csv(BASELINE_TARGET_METRICS, low_memory=False)
    baseline = baseline[
        baseline["model"].isin(["ipTM_only", "ipSAE_only"])
        & baseline["layer"].isin(base.LAYERS)
    ][
        [
            "model",
            "source_weighting",
            "target",
            "layer",
            "roc_auc",
            "ap",
            "brier",
            "log_loss",
            "bss_transport",
        ]
    ].copy()
    combo = target_metrics[
        [
            "model",
            "source_weighting",
            "target",
            "layer",
            "roc_auc",
            "ap",
            "brier",
            "log_loss",
            "bss_transport",
        ]
    ].copy()
    combined = pd.concat([baseline, combo], ignore_index=True)
    if combined.groupby(["source_weighting", "target", "layer"]).size().ne(3).any():
        raise RuntimeError("C01B comparison table is incomplete")
    rows = []
    for (weighting, target, layer), group in combined.groupby(
        ["source_weighting", "target", "layer"], sort=True
    ):
        by_model = group.set_index("model")
        for comparator in ("ipTM_only", "ipSAE_only"):
            rows.append(
                {
                    "source_weighting": weighting,
                    "target": target,
                    "layer": layer,
                    "comparison": f"ipTM_plus_ipSAE_vs_{comparator}",
                    "auc_combination_minus_comparator": by_model.loc[
                        "ipTM_plus_ipSAE", "roc_auc"
                    ]
                    - by_model.loc[comparator, "roc_auc"],
                    "ap_combination_minus_comparator": by_model.loc[
                        "ipTM_plus_ipSAE", "ap"
                    ]
                    - by_model.loc[comparator, "ap"],
                    "brier_combination_minus_comparator": by_model.loc[
                        "ipTM_plus_ipSAE", "brier"
                    ]
                    - by_model.loc[comparator, "brier"],
                    "log_loss_combination_minus_comparator": by_model.loc[
                        "ipTM_plus_ipSAE", "log_loss"
                    ]
                    - by_model.loc[comparator, "log_loss"],
                    "bss_transport_combination_minus_comparator": by_model.loc[
                        "ipTM_plus_ipSAE", "bss_transport"
                    ]
                    - by_model.loc[comparator, "bss_transport"],
                }
            )
    per_target = pd.DataFrame(rows)
    macro_rows = []
    for (weighting, layer, comparison), group in per_target.groupby(
        ["source_weighting", "layer", "comparison"], sort=True
    ):
        row = {
            "source_weighting": weighting,
            "target": "__macro_target_equal__",
            "layer": layer,
            "comparison": comparison,
            "target_count": len(group),
        }
        for metric in (
            "auc_combination_minus_comparator",
            "ap_combination_minus_comparator",
            "brier_combination_minus_comparator",
            "log_loss_combination_minus_comparator",
            "bss_transport_combination_minus_comparator",
        ):
            row[metric] = float(group[metric].mean())
            row[f"{metric}_improved_targets"] = int(
                (group[metric] > 0).sum()
                if metric.startswith(("auc", "ap", "bss"))
                else (group[metric] < 0).sum()
            )
        macro_rows.append(row)
    return pd.concat(
        [
            per_target.assign(summary_type="per_target"),
            pd.DataFrame(macro_rows).assign(summary_type="target_equal_macro"),
        ],
        ignore_index=True,
        sort=False,
    )


def main() -> None:
    COMBO_OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(BASE_PAIR_EXTENSION, low_memory=False)
    if len(data) != 2074 or int(data["primary_outcome"].sum()) != 740:
        raise RuntimeError("C01B complete-case source mismatch")
    data["outcome_value"] = data["primary_outcome"].astype(int)
    data["x_ipTM"] = base.logit_clip(data["boltz2_iptm"])
    data["x_ipSAE"] = base.logit_clip(data["official_ipsae_min"])
    u_predictions, source_parameters = build_u_predictions(data)
    (
        adapted,
        adaptation_parameters,
        lambda_selection,
        target_metrics,
    ) = build_combo_metrics(u_predictions)
    macro_rows = []
    for (weighting, layer), group in target_metrics.groupby(
        ["source_weighting", "layer"], sort=True
    ):
        row = {
            "model": "ipTM_plus_ipSAE",
            "source_weighting": weighting,
            "layer": layer,
            "summary_type": "target_equal_macro",
            "target_count": len(group),
            "n_pairs": int(group["n"].sum()),
            "positives": int(group["positives"].sum()),
        }
        for metric in (
            "prevalence",
            "roc_auc",
            "ap",
            "brier",
            "log_loss",
            "bss_transport",
            "bss_target_oracle",
            "bss_adaptation_train",
        ):
            row[metric] = float(group[metric].mean())
        macro_rows.append(row)
    macro_metrics = pd.DataFrame(macro_rows)
    comparisons = build_comparisons(target_metrics)

    base.write_csv(COMBO_OUT / "combo_UIF_pair_predictions.csv", adapted)
    base.write_csv(COMBO_OUT / "combo_source_parameters.csv", source_parameters)
    base.write_csv(
        COMBO_OUT / "combo_adaptation_parameters.csv",
        adaptation_parameters,
    )
    base.write_csv(COMBO_OUT / "combo_lambda_selection.csv", lambda_selection)
    base.write_csv(COMBO_OUT / "combo_target_metrics.csv", target_metrics)
    base.write_csv(COMBO_OUT / "combo_macro_metrics.csv", macro_metrics)
    base.write_csv(COMBO_OUT / "combo_paired_comparisons.csv", comparisons)

    macro = comparisons[
        comparisons["summary_type"].eq("target_equal_macro")
        & comparisons["source_weighting"].eq("target_equal")
    ]
    target_equal_metrics = macro_metrics[
        macro_metrics["source_weighting"].eq("target_equal")
    ].set_index("layer")
    report = [
        "# C01B: ipTM + ipSAE complementarity",
        "",
        "## Design",
        "",
        "- Cohort: the same 2,074 complete-case pairs, including 740 primary positives and six targets.",
        "- Models compared: `ipTM_only`, `ipSAE_only`, and `ipTM_plus_ipSAE`.",
        "- Evidence level: post hoc sensitivity extension added after the single-score C01 results had been examined. The combination execution rules were frozen before the combination result was inspected.",
        "- U used one intercept per source target and two jointly estimated common slopes on `logit(ipTM)` and `logit(ipSAE_min)`.",
        "- I updated only an overall intercept; F updated an overall intercept plus one scalar calibration slope.",
        "- The relative ipTM/ipSAE weights were not re-estimated locally.",
        "- Frozen primary-six LOTO and target-CV5 assignments, 90% sequence clusters, and both existing source-weighting schemes were reused.",
        "- No feature selection, penalty tuning, Pinc combination, or random forest was run.",
        "",
        "## Target-equal combination metrics",
        "",
        "| Layer | AUC | AP | Brier | Log loss | BSS-transport |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for layer in base.LAYERS:
        row = target_equal_metrics.loc[layer]
        report.append(
            f"| {layer} | {row.roc_auc:.4f} | {row.ap:.4f} | "
            f"{row.brier:.4f} | {row.log_loss:.4f} | {row.bss_transport:.4f} |"
        )
    report += [
        "",
        "## Target-equal paired comparisons",
        "",
        "| Comparison/layer | AUC difference | AP difference | Brier difference | Log-loss difference | BSS difference |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in macro.sort_values(["comparison", "layer"]).itertuples(
        index=False
    ):
        report.append(
            f"| {row.comparison}/{row.layer} | "
            f"{row.auc_combination_minus_comparator:.4f} | "
            f"{row.ap_combination_minus_comparator:.4f} | "
            f"{row.brier_combination_minus_comparator:.4f} | "
            f"{row.log_loss_combination_minus_comparator:.4f} | "
            f"{row.bss_transport_combination_minus_comparator:.4f} |"
        )
    target_vs_iptm = comparisons[
        comparisons["summary_type"].eq("per_target")
        & comparisons["comparison"].eq("ipTM_plus_ipSAE_vs_ipTM_only")
    ]
    consistency_lines = []
    for weighting in base.WEIGHTINGS:
        for layer in base.LAYERS:
            group = target_vs_iptm[
                target_vs_iptm["source_weighting"].eq(weighting)
                & target_vs_iptm["layer"].eq(layer)
            ]
            consistency_lines.append(
                f"| {weighting} | {layer} | "
                f"{int((group['auc_combination_minus_comparator'] > 0).sum())}/6 | "
                f"{int((group['ap_combination_minus_comparator'] > 0).sum())}/6 | "
                f"{int((group['brier_combination_minus_comparator'] < 0).sum())}/6 | "
                f"{int((group['log_loss_combination_minus_comparator'] < 0).sum())}/6 | "
                f"{int((group['bss_transport_combination_minus_comparator'] > 0).sum())}/6 |"
            )
    slope_summary = (
        source_parameters.groupby("source_weighting")
        .agg(
            beta_ipTM_min=("beta_ipTM", "min"),
            beta_ipTM_max=("beta_ipTM", "max"),
            beta_ipSAE_min=("beta_ipSAE", "min"),
            beta_ipSAE_max=("beta_ipSAE", "max"),
        )
        .reset_index()
    )
    report += [
        "",
        "Negative Brier/log-loss differences favor the combination; positive AUC/AP/BSS differences favor it.",
        "These are point estimates from the existing folds, not a new model-selection exercise.",
        "",
        "## Direction consistency versus ipTM-only",
        "",
        "| Weighting | Layer | AUC improved | AP improved | Brier improved | Log loss improved | BSS improved |",
        "|---|---|---:|---:|---:|---:|---:|",
        *consistency_lines,
        "",
        "## Conditional source slopes",
        "",
        "| Weighting | ipTM slope range | ipSAE slope range |",
        "|---|---:|---:|",
    ]
    for row in slope_summary.itertuples(index=False):
        report.append(
            f"| {row.source_weighting} | {row.beta_ipTM_min:.4f} to "
            f"{row.beta_ipTM_max:.4f} | {row.beta_ipSAE_min:.4f} to "
            f"{row.beta_ipSAE_max:.4f} |"
        )
    report += [
        "",
        "The conditional ipSAE slope was negative in every source fold and under both weighting schemes. "
        "Because ipTM and ipSAE are strongly correlated, this sign should not be interpreted as a causal "
        "or mechanistic effect. It does show that adding ipSAE did not provide a stable positive marginal "
        "signal after ipTM.",
        "",
        "## Interpretation",
        "",
        "The combination improved at least one probability-error metric in some configurations, but ranking "
        "was not consistently improved and the gains were small and target-dependent. The appropriate "
        "manuscript statement is that ipSAE, either alone or jointly with ipTM, did not show a stable "
        "incremental improvement in this evaluated cohort; this does not establish that ipSAE is generally "
        "uninformative or that no other calibration map could work.",
        "",
    ]
    (COMBO_OUT / "C01B_RESULTS.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )
    manifest = {
        "status": "complete",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "analysis": "C01B",
        "analysis_role": "post_hoc_sensitivity_extension",
        "complete_cases": len(data),
        "positives": int(data["outcome_value"].sum()),
        "models": ["ipTM_only", "ipSAE_only", "ipTM_plus_ipSAE"],
        "combination_model": {
            "U": "target-specific intercept + unpenalized ipTM slope + unpenalized ipSAE slope",
            "I": "overall intercept update on the fixed combination linear predictor",
            "F": "overall intercept plus one ridge-shrunk scalar slope update; relative feature weights fixed",
            "lambda_grid": list(base.LAMBDAS),
        },
        "random_forest_run": False,
        "other_features_added": False,
        "outputs": {
            str(path.relative_to(ROOT)): base.sha256(path)
            for path in (
                COMBO_OUT / "combo_UIF_pair_predictions.csv",
                COMBO_OUT / "combo_source_parameters.csv",
                COMBO_OUT / "combo_adaptation_parameters.csv",
                COMBO_OUT / "combo_lambda_selection.csv",
                COMBO_OUT / "combo_target_metrics.csv",
                COMBO_OUT / "combo_macro_metrics.csv",
                COMBO_OUT / "combo_paired_comparisons.csv",
                COMBO_OUT / "C01B_RESULTS.md",
            )
        },
        "inputs": {
            str(path.relative_to(ROOT)): base.sha256(path)
            for path in (BASE_PAIR_EXTENSION, BASELINE_TARGET_METRICS, base.LOTO, base.CV5)
        },
        "software": {
            "python": platform.python_version(),
            "numpy": base.np.__version__,
            "pandas": base.pd.__version__,
            "scipy": base.scipy.__version__,
            "scikit_learn": base.sklearn.__version__,
        },
        "script_sha256": base.sha256(Path(__file__)),
    }
    base.write_json(COMBO_OUT / "run_manifest.json", manifest)


if __name__ == "__main__":
    main()
