"""Frozen minimal RQ6 execution.

This script rebuilds the 720 x 10 paired predictor score table directly from
the released parquet, computes the frozen descriptive metrics for both
prespecified Adaptyv outcomes, and runs the target-stratified paired bootstrap.
It does not fit a model or change any frozen input.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import traceback
from datetime import datetime
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score


ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN_DIR = Path(__file__).resolve().parent

PROTOCOL = ROOT / "analysis/RQ6_minimal_descriptive_protocol_2026-11-15.md"
START = ROOT / "analysis/codex_tasks/RQ6_EXECUTION_CODEX_2026-11-15/START_HERE_CODEX.md"
AUDIT = ROOT / "analysis/audit_outputs/RQ6_availability_schema_2026-11-15/AVAILABILITY_AUDIT.md"
ACCEPTANCE = ROOT / "analysis/audit_outputs/RQ6_availability_acceptance_2026-11-15/ACCEPTANCE_RESULT.md"

DATA_DIR = ROOT / "data/external_target_audit_2026-11-13"
COFOLD = DATA_DIR / "anthropic_cofold_predictions.parquet"
PAIR_TABLE = DATA_DIR / "anthropic_audited_pair_table.csv"
LOOKUP = DATA_DIR / "anthropic_LOOKUP_TABLES.md"
INSILICO = DATA_DIR / "anthropic_INSILICO.md"
WETLAB = DATA_DIR / "anthropic_WETLAB.md"
AVAILABILITY = (
    ROOT
    / "analysis/audit_outputs/RQ6_availability_schema_2026-11-15"
    / "anthropic_selected_pair_score_availability.csv"
)
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"

TARGETS = [
    "BHRF1",
    "IL-7Ra",
    "RBX1",
    "TrkA",
    "Nipah-G",
    "PD-L1",
    "TREM2",
    "VEGF-A",
]
PREDICTORS = [
    "boltz2",
    "afm3",
    "chai1",
    "ef2fast",
    "ef2full",
    "ptxv2",
    "of3",
    "odde",
    "rf3",
    "af3of3",
]
OUTCOMES = {
    "adaptyv_weak_plus": "adaptyv_main_positive",
    "adaptyv_medium_plus": "adaptyv_strict_positive",
}
EXPECTED_POSITIVES = {
    "BHRF1": (20, 18),
    "IL-7Ra": (47, 36),
    "RBX1": (23, 20),
    "TrkA": (19, 19),
    "Nipah-G": (19, 17),
    "PD-L1": (39, 35),
    "TREM2": (65, 62),
    "VEGF-A": (55, 47),
}

EXPECTED_HASHES = {
    COFOLD: "ba48010610daa9b9d184aa1ddb419beebdf728ccf745f5a8f3b48e4041e1a87d",
    PAIR_TABLE: "1f2a4ceab9de06cb75f4aae9a889fc5164e144232101417653944d03ee33d0c6",
    MASTER: "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91",
}
BOOTSTRAP_SEED = 20261115
BOOTSTRAP_REPLICATES = 2000
LOG_LOSS_CLIP = 1e-6
PAIR_TOLERANCE = 1e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def write_text(path: Path, text: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def as_bool(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    return series.astype(str).str.strip().str.lower().eq("true").to_numpy(dtype=bool)


def auc_mann_whitney(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    positives = int(y.sum())
    negatives = int(len(y) - positives)
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = rankdata(score, method="average")
    rank_sum = float(ranks[y == 1].sum())
    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def brier(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    return float(np.mean((score - y) ** 2))


def log_loss(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    p = np.clip(np.asarray(score, dtype=float), LOG_LOSS_CLIP, 1.0 - LOG_LOSS_CLIP)
    return float(-np.mean(y * np.log(p) + (1.0 - y) * np.log(1.0 - p)))


def ap_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    if int(y.sum()) == 0 or int(y.sum()) == len(y):
        return float("nan")
    return float(average_precision_score(y, score))


def correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    if method == "spearman":
        return float(spearmanr(x, y).statistic)
    return float(np.corrcoef(x, y)[0, 1])


def fmt(value: object, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(number):
        return "NA"
    return f"{number:.{digits}f}"


def metric_row(
    target: str,
    predictor: str,
    outcome_name: str,
    y: np.ndarray,
    score: np.ndarray,
    boltz_score: np.ndarray,
    boltz_metrics: dict[str, float],
) -> dict[str, object]:
    metrics = {
        "auc": auc_mann_whitney(y, score),
        "ap": ap_score(y, score),
        "brier": brier(y, score),
        "log_loss": log_loss(y, score),
    }
    return {
        "target": target,
        "predictor": predictor,
        "outcome": outcome_name,
        "n": len(y),
        "positives": int(np.sum(y)),
        "observed_prevalence": float(np.mean(y)),
        "mean_iptm": float(np.mean(score)),
        "auc": metrics["auc"],
        "ap": metrics["ap"],
        "brier": metrics["brier"],
        "log_loss": metrics["log_loss"],
        "calibration_in_large": float(np.mean(score) - np.mean(y)),
        "delta_auc_vs_boltz2": metrics["auc"] - boltz_metrics["auc"],
        "delta_ap_vs_boltz2": metrics["ap"] - boltz_metrics["ap"],
        "delta_brier_vs_boltz2": metrics["brier"] - boltz_metrics["brier"],
        "delta_log_loss_vs_boltz2": metrics["log_loss"] - boltz_metrics["log_loss"],
    }


def bootstrap_interval(values: list[float], estimate: float) -> dict[str, object]:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if len(finite):
        low, high = np.percentile(finite, [2.5, 97.5])
        median = float(np.median(finite))
    else:
        low = high = median = float("nan")
    return {
        "estimate": float(estimate),
        "bootstrap_median": median,
        "ci_2_5_percentile": float(low),
        "ci_97_5_percentile": float(high),
        "valid_iterations": int(len(finite)),
        "requested_iterations": BOOTSTRAP_REPLICATES,
    }


def rank_winners(
    metrics: pd.DataFrame, outcome: str, metric: str, higher_is_better: bool
) -> pd.DataFrame:
    subset = metrics.loc[metrics.outcome.eq(outcome)]
    if higher_is_better:
        winners = subset.loc[subset.groupby("target")[metric].transform("max").eq(subset[metric])]
    else:
        winners = subset.loc[subset.groupby("target")[metric].transform("min").eq(subset[metric])]
    return winners[["target", "predictor", metric]].sort_values(["target", "predictor"])


def make_heatmap(
    axis: plt.Axes,
    matrix: pd.DataFrame,
    title: str,
    cmap: str,
    center: float | None = None,
    value_format: str = ".3f",
) -> object:
    values = matrix.to_numpy(dtype=float)
    if center is None:
        image = axis.imshow(values, aspect="auto", cmap=cmap)
    else:
        limit = max(float(np.nanmax(np.abs(values - center))), 1e-6)
        image = axis.imshow(values, aspect="auto", cmap=cmap, vmin=center - limit, vmax=center + limit)
    axis.set_xticks(np.arange(matrix.shape[1]), matrix.columns, rotation=40, ha="right")
    axis.set_yticks(np.arange(matrix.shape[0]), matrix.index)
    axis.set_title(title, fontsize=11)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            axis.text(
                column,
                row,
                format(values[row, column], value_format),
                ha="center",
                va="center",
                fontsize=7,
                color="black",
            )
    return image


def main() -> None:
    if not RUN_DIR.name.startswith("run_"):
        raise RuntimeError(f"Formal output must be inside a run_* directory: {RUN_DIR}")
    log_path = RUN_DIR / "execution.log"

    def log(message: str) -> None:
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    write_json(
        RUN_DIR / "RUN_STATUS.json",
        {"status": "running", "scope": "frozen_minimal_RQ6", "run_dir": str(RUN_DIR)},
    )
    log("Starting frozen minimal RQ6 execution.")
    log("Checking required frozen inputs and SHA-256 gates.")

    input_paths = [
        START,
        PROTOCOL,
        AUDIT,
        ACCEPTANCE,
        COFOLD,
        PAIR_TABLE,
        AVAILABILITY,
        LOOKUP,
        INSILICO,
        WETLAB,
        MASTER,
    ]
    missing = [str(path) for path in input_paths if not path.exists()]
    if missing:
        raise RuntimeError(f"Required input missing: {missing}")
    input_hashes_before = {str(path): sha256(path) for path in input_paths}
    mismatches = {
        str(path): {"expected": expected, "actual": input_hashes_before[str(path)]}
        for path, expected in EXPECTED_HASHES.items()
        if input_hashes_before[str(path)] != expected
    }
    if mismatches:
        raise RuntimeError(f"Frozen SHA-256 mismatch: {mismatches}")
    log("Frozen core input hashes match the execution contract.")

    design = pd.read_csv(PAIR_TABLE, low_memory=False)
    availability = pd.read_csv(AVAILABILITY, low_memory=False)
    if len(design) != 1440 or design.uuid.duplicated().any():
        raise RuntimeError("Released audited pair table structure mismatch.")
    if len(availability) != 1440 or availability.uuid.duplicated().any():
        raise RuntimeError("Availability connection table structure mismatch.")

    cohort = design.loc[design.target.isin(TARGETS)].copy()
    if len(cohort) != 720 or cohort.uuid.nunique() != 720:
        raise RuntimeError("Frozen RQ6 cohort is not exactly 720 unique designs.")
    if set(cohort.target) != set(TARGETS):
        raise RuntimeError("Frozen 8-target set mismatch.")
    if not cohort.groupby("target").size().eq(90).all():
        raise RuntimeError("Frozen cohort does not have 90 designs in every target.")
    if not as_bool(cohort.adaptyv_valid).all():
        raise RuntimeError("The primary outcome is not valid for all 720 designs.")

    for index, outcome_name in enumerate(OUTCOMES):
        if outcome_name == "adaptyv_weak_plus":
            positive_column = "adaptyv_main_positive"
            expected_total = 287
        else:
            positive_column = "adaptyv_strict_positive"
            expected_total = 254
        positives = as_bool(cohort[positive_column]).astype(int)
        if int(positives.sum()) != expected_total:
            raise RuntimeError(f"{positive_column} total mismatch.")
        for target, counts in EXPECTED_POSITIVES.items():
            actual = int(positives[cohort.target.eq(target).to_numpy()].sum())
            if actual != counts[index]:
                raise RuntimeError(f"{positive_column} positive count mismatch for {target}.")

    log("Cohort/outcome QC passed: 720 designs; Weak+ 287; Medium+ 254.")
    log("Rebuilding predictor scores directly from the parquet.")

    connection = duckdb.connect()
    raw_candidates = connection.execute(
        """
        SELECT uuid, target, cofolding_model, stoichiometry, target_form, seed,
               ipsae_min, iptm_pae
        FROM read_parquet(?)
        WHERE target_form IS NULL OR target_form = ''
        """,
        [str(COFOLD)],
    ).df()
    connection.close()
    raw_candidates = raw_candidates.loc[raw_candidates.cofolding_model.isin(PREDICTORS)].copy()
    raw_candidates = raw_candidates.merge(
        cohort[["uuid", "target", "boltz2_stoichiometry"]],
        on="uuid",
        how="inner",
        validate="many_to_one",
    )
    if not raw_candidates.target_x.eq(raw_candidates.target_y).all():
        raise RuntimeError("Target disagreement after joining parquet to the design table.")
    raw_candidates = raw_candidates.loc[
        raw_candidates.stoichiometry.eq(raw_candidates.boltz2_stoichiometry)
    ].copy()
    raw_candidates["seed_numeric"] = pd.to_numeric(raw_candidates.seed, errors="raise").astype(int)
    if raw_candidates.duplicated(["uuid", "cofolding_model", "seed_numeric"]).any():
        raise RuntimeError("Duplicated uuid x predictor x numeric run/model label.")

    run_counts = raw_candidates.groupby(["uuid", "cofolding_model"]).size()
    if len(run_counts) != 7200 or not run_counts.eq(5).all():
        raise RuntimeError("Every uuid x predictor must have exactly five released runs/models.")

    selected = raw_candidates.sort_values(
        ["uuid", "cofolding_model", "ipsae_min", "seed_numeric"],
        ascending=[True, True, False, True],
        kind="mergesort",
    ).drop_duplicates(["uuid", "cofolding_model"], keep="first")
    if len(selected) != 7200:
        raise RuntimeError("Selected predictor score table is not 7,200 rows.")
    if not np.isfinite(pd.to_numeric(selected.iptm_pae, errors="raise").to_numpy(dtype=float)).all():
        raise RuntimeError("A selected ipTM is missing or non-finite.")

    score_wide = selected.pivot(
        index="uuid", columns="cofolding_model", values="iptm_pae"
    ).reindex(columns=PREDICTORS)
    score_wide.columns = [f"iptm__{predictor}" for predictor in PREDICTORS]
    score_wide = score_wide.reindex(cohort.uuid)
    if score_wide.shape != (720, 10) or score_wide.isna().any().any():
        raise RuntimeError("The rebuilt 720 x 10 score table is incomplete.")

    availability_indexed = availability.set_index("uuid").reindex(score_wide.index)
    availability_columns = [f"iptm__{predictor}" for predictor in PREDICTORS]
    availability_diff = np.abs(
        score_wide[availability_columns].to_numpy(dtype=float)
        - availability_indexed[availability_columns].to_numpy(dtype=float)
    )
    availability_max_diff = float(np.nanmax(availability_diff))
    if not np.isfinite(availability_max_diff) or availability_max_diff > PAIR_TOLERANCE:
        raise RuntimeError("Rebuilt scores disagree with the availability connection table.")

    design_indexed = cohort.set_index("uuid").reindex(score_wide.index)
    boltz_rebuild_diff = np.abs(
        score_wide["iptm__boltz2"].to_numpy(dtype=float)
        - pd.to_numeric(
            design_indexed["boltz2_iptm_ipsae_best_seed"], errors="raise"
        ).to_numpy(dtype=float)
    )
    boltz_max_diff = float(np.max(boltz_rebuild_diff))
    if boltz_max_diff > PAIR_TOLERANCE:
        raise RuntimeError("Rebuilt Boltz-2 scores disagree with RQ3-E published scores.")
    log(
        f"Score reconstruction QC passed; max availability diff={availability_max_diff:.3e}, "
        f"max Boltz-2 audit diff={boltz_max_diff:.3e}."
    )

    wide = pd.DataFrame(
        {
            "uuid": score_wide.index.astype(str),
            "target": design_indexed.target.to_numpy(dtype=str),
            "adaptyv_valid": as_bool(design_indexed.adaptyv_valid),
            "adaptyv_main_positive": as_bool(design_indexed.adaptyv_main_positive).astype(int),
            "adaptyv_strict_positive": as_bool(design_indexed.adaptyv_strict_positive).astype(int),
            "n_predictors": 10,
        }
    )
    for column in availability_columns:
        wide[column] = score_wide[column].to_numpy(dtype=float)
    wide = wide.sort_values(["target", "uuid"], kind="mergesort").reset_index(drop=True)
    if (
        len(wide) != 720
        or wide.uuid.nunique() != 720
        or wide[[f"iptm__{predictor}" for predictor in PREDICTORS]].isna().any().any()
    ):
        raise RuntimeError("Wide deliverable failed final structural QC.")

    long = wide.melt(
        id_vars=[
            "uuid",
            "target",
            "adaptyv_valid",
            "adaptyv_main_positive",
            "adaptyv_strict_positive",
        ],
        value_vars=availability_columns,
        var_name="score_column",
        value_name="iptm_pae",
    )
    long["predictor"] = long.score_column.str.removeprefix("iptm__")
    long = long.drop(columns="score_column")[
        [
            "uuid",
            "target",
            "predictor",
            "iptm_pae",
            "adaptyv_valid",
            "adaptyv_main_positive",
            "adaptyv_strict_positive",
        ]
    ].sort_values(["target", "uuid", "predictor"], kind="mergesort")
    if len(long) != 7200 or long.duplicated(["uuid", "predictor"]).any():
        raise RuntimeError("Long deliverable is not 7,200 unique uuid x predictor rows.")

    selected_audit = selected[
        ["uuid", "cofolding_model", "stoichiometry", "seed", "seed_numeric", "ipsae_min", "iptm_pae"]
    ].rename(columns={"cofolding_model": "predictor", "seed": "selected_run_or_model_label"})
    selected_audit = selected_audit.sort_values(["uuid", "predictor"], kind="mergesort")

    write_csv(RUN_DIR / "rq6_pair_scores_wide.csv", wide)
    write_csv(RUN_DIR / "rq6_pair_scores_long.csv", long)
    write_csv(RUN_DIR / "rq6_selected_run_audit.csv", selected_audit)

    target_score = {
        target: wide.loc[wide.target.eq(target), availability_columns].to_numpy(dtype=float)
        for target in TARGETS
    }
    target_y = {
        outcome: {
            target: wide.loc[
                wide.target.eq(target), positive_column
            ].to_numpy(dtype=int)
            for target in TARGETS
        }
        for outcome, positive_column in OUTCOMES.items()
    }

    log("Computing shared score shifts and predictor correlations.")
    score_rows = []
    for target in TARGETS:
        scores = target_score[target]
        boltz = scores[:, 0]
        for predictor_index, predictor in enumerate(PREDICTORS):
            values = scores[:, predictor_index]
            difference = values - boltz
            score_rows.append(
                {
                    "summary_level": "target",
                    "target": target,
                    "predictor": predictor,
                    "outcome_scope": "shared_main_and_strict",
                    "n": len(values),
                    "score_mean": float(np.mean(values)),
                    "score_sd": float(np.std(values, ddof=1)),
                    "score_median": float(np.median(values)),
                    "score_q1": float(np.percentile(values, 25)),
                    "score_q3": float(np.percentile(values, 75)),
                    "score_min": float(np.min(values)),
                    "score_max": float(np.max(values)),
                    "difference_mean": float(np.mean(difference)),
                    "difference_median": float(np.median(difference)),
                    "difference_q1": float(np.percentile(difference, 25)),
                    "difference_q3": float(np.percentile(difference, 75)),
                    "mean_absolute_difference": float(np.mean(np.abs(difference))),
                    "spearman_vs_boltz2": correlation(values, boltz, "spearman"),
                    "pearson_vs_boltz2": correlation(values, boltz, "pearson"),
                    "target_equal_aggregation": "",
                }
            )
    score_by_target = pd.DataFrame(score_rows)

    target_equal_score_rows = []
    for predictor in PREDICTORS:
        subset = score_by_target.loc[score_by_target.predictor.eq(predictor)]
        row = {
            "summary_level": "target_equal",
            "target": "ALL_8_EQUAL_WEIGHT",
            "predictor": predictor,
            "outcome_scope": "shared_main_and_strict",
            "n": 720,
            "score_mean": float(subset.score_mean.mean()),
            "score_sd": float("nan"),
            "score_median": float(subset.score_median.median()),
            "score_q1": float(subset.score_q1.median()),
            "score_q3": float(subset.score_q3.median()),
            "score_min": float(subset.score_min.min()),
            "score_max": float(subset.score_max.max()),
            "difference_mean": float(subset.difference_mean.mean()),
            "difference_median": float(subset.difference_median.median()),
            "difference_q1": float(subset.difference_q1.median()),
            "difference_q3": float(subset.difference_q3.median()),
            "mean_absolute_difference": float(subset.mean_absolute_difference.mean()),
            "spearman_vs_boltz2": float(subset.spearman_vs_boltz2.mean()),
            "pearson_vs_boltz2": float(subset.pearson_vs_boltz2.mean()),
            "target_equal_aggregation": (
                "mean_of_target_means; median_of_target_medians; "
                "median_of_target_Q1s; median_of_target_Q3s; "
                "min_of_target_mins; max_of_target_maxs; "
                "mean_target_mean_absolute_difference; mean_target_correlation"
            ),
        }
        target_equal_score_rows.append(row)
    score_shift = pd.concat(
        [score_by_target, pd.DataFrame(target_equal_score_rows)], ignore_index=True
    )
    target_range_rows = []
    for predictor in PREDICTORS:
        subset = score_by_target.loc[score_by_target.predictor.eq(predictor)]
        target_range_rows.append(
            {
                "predictor": predictor,
                "median_score_range_across_targets": float(
                    subset.score_median.max() - subset.score_median.min()
                ),
                "median_difference_range_across_targets": float(
                    subset.difference_median.max() - subset.difference_median.min()
                ),
                "spearman_min_across_targets": float(subset.spearman_vs_boltz2.min()),
                "spearman_max_across_targets": float(subset.spearman_vs_boltz2.max()),
            }
        )
    score_ranges = pd.DataFrame(target_range_rows)
    write_csv(RUN_DIR / "rq6_score_shift_by_target.csv", score_shift)
    write_csv(RUN_DIR / "rq6_score_shift_target_ranges.csv", score_ranges)

    log("Computing target x predictor x outcome AUC, AP, Brier, and log loss.")
    metric_rows = []
    point_metrics: dict[tuple[str, str], dict[str, float]] = {}
    for outcome, _ in OUTCOMES.items():
        for target in TARGETS:
            y = target_y[outcome][target]
            scores = target_score[target]
            boltz = scores[:, 0]
            boltz_metrics = {
                "auc": auc_mann_whitney(y, boltz),
                "ap": ap_score(y, boltz),
                "brier": brier(y, boltz),
                "log_loss": log_loss(y, boltz),
            }
            for predictor_index, predictor in enumerate(PREDICTORS):
                row = metric_row(
                    target,
                    predictor,
                    outcome,
                    y,
                    scores[:, predictor_index],
                    boltz,
                    boltz_metrics,
                )
                metric_rows.append(row)
                point_metrics[(outcome, target, predictor)] = row
    metrics = pd.DataFrame(metric_rows)
    if len(metrics) != 160:
        raise RuntimeError("Metrics table is not 160 rows (8 targets x 10 predictors x 2 outcomes).")
    write_csv(RUN_DIR / "rq6_metrics_by_target.csv", metrics)

    log("Building target-equal metric summaries and predictor ranks.")
    metric_names = [
        "auc",
        "ap",
        "brier",
        "log_loss",
        "calibration_in_large",
        "delta_auc_vs_boltz2",
        "delta_ap_vs_boltz2",
        "delta_brier_vs_boltz2",
        "delta_log_loss_vs_boltz2",
    ]
    target_equal_rows = []
    for outcome in OUTCOMES:
        for metric in metric_names:
            outcome_predictor_means = []
            for predictor in PREDICTORS:
                values = metrics.loc[
                    metrics.outcome.eq(outcome) & metrics.predictor.eq(predictor), metric
                ].to_numpy(dtype=float)
                row = {
                    "outcome": outcome,
                    "predictor": predictor,
                    "metric": metric,
                    "n_targets": int(np.isfinite(values).sum()),
                    "target_equal_mean": float(np.nanmean(values)),
                    "target_equal_median": float(np.nanmedian(values)),
                    "target_equal_min": float(np.nanmin(values)),
                    "target_equal_max": float(np.nanmax(values)),
                    "target_equal_range": float(np.nanmax(values) - np.nanmin(values)),
                    "rank_by_target_equal_mean": float("nan"),
                }
                target_equal_rows.append(row)
                outcome_predictor_means.append((predictor, row["target_equal_mean"]))
            if metric in {"auc", "ap"}:
                sorted_predictors = sorted(
                    outcome_predictor_means, key=lambda item: (-item[1], PREDICTORS.index(item[0]))
                )
                direction = "higher_is_better"
            elif metric in {"brier", "log_loss"}:
                sorted_predictors = sorted(
                    outcome_predictor_means, key=lambda item: (item[1], PREDICTORS.index(item[0]))
                )
                direction = "lower_is_better"
            else:
                sorted_predictors = []
                direction = "not_ranked"
            rank_map = {predictor: rank for rank, (predictor, _) in enumerate(sorted_predictors, 1)}
            for row in target_equal_rows:
                if row["outcome"] == outcome and row["metric"] == metric:
                    row["rank_by_target_equal_mean"] = rank_map.get(row["predictor"], float("nan"))
                    row["rank_direction"] = direction
    target_equal = pd.DataFrame(target_equal_rows)
    write_csv(RUN_DIR / "rq6_target_equal_summary.csv", target_equal)

    log("Running 2,000 target-stratified paired bootstrap replicates.")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_values: dict[tuple[str, str, str, str], list[float]] = {}
    score_bootstrap: dict[tuple[str, str, str], list[float]] = {}

    def add_boot(key: tuple[str, ...], value: float) -> None:
        bootstrap_values.setdefault(key, []).append(float(value))

    for target in TARGETS:
        scores = target_score[target]
        boltz = scores[:, 0]
        n = len(scores)
        for _ in range(BOOTSTRAP_REPLICATES):
            draw = rng.integers(0, n, size=n)
            sampled_scores = scores[draw, :]
            sampled_boltz = boltz[draw]
            for predictor_index, predictor in enumerate(PREDICTORS):
                sampled_predictor = sampled_scores[:, predictor_index]
                difference = sampled_predictor - sampled_boltz
                score_values = {
                    "score_difference_mean": float(np.mean(difference)),
                    "score_difference_median": float(np.median(difference)),
                    "score_difference_q1": float(np.percentile(difference, 25)),
                    "score_difference_q3": float(np.percentile(difference, 75)),
                    "score_mean_absolute_difference": float(np.mean(np.abs(difference))),
                    "spearman_vs_boltz2": correlation(
                        sampled_predictor, sampled_boltz, "spearman"
                    ),
                    "pearson_vs_boltz2": correlation(
                        sampled_predictor, sampled_boltz, "pearson"
                    ),
                }
                for metric_name, value in score_values.items():
                    score_bootstrap.setdefault(
                        (target, predictor, metric_name), []
                    ).append(value)

                outcome_metrics: dict[str, dict[str, float]] = {}
                for outcome in OUTCOMES:
                    sampled_y = target_y[outcome][target][draw]
                    outcome_metrics[outcome] = {
                        "auc": auc_mann_whitney(sampled_y, sampled_predictor),
                        "ap": ap_score(sampled_y, sampled_predictor),
                        "brier": brier(sampled_y, sampled_predictor),
                        "log_loss": log_loss(sampled_y, sampled_predictor),
                    }
                boltz_outcome_metrics: dict[str, dict[str, float]] = {}
                for outcome in OUTCOMES:
                    sampled_y = target_y[outcome][target][draw]
                    boltz_outcome_metrics[outcome] = {
                        "auc": auc_mann_whitney(sampled_y, sampled_boltz),
                        "ap": ap_score(sampled_y, sampled_boltz),
                        "brier": brier(sampled_y, sampled_boltz),
                        "log_loss": log_loss(sampled_y, sampled_boltz),
                    }
                for outcome in OUTCOMES:
                    for metric_name in ["auc", "ap", "brier", "log_loss"]:
                        value = (
                            outcome_metrics[outcome][metric_name]
                            - boltz_outcome_metrics[outcome][metric_name]
                        )
                        add_boot(
                            (
                                "outcome",
                                outcome,
                                target,
                                predictor,
                                f"delta_{metric_name}_vs_boltz2",
                            ),
                            value,
                        )

    bootstrap_rows = []
    for target in TARGETS:
        for predictor in PREDICTORS:
            score_specs = {
                "score_difference_mean": float(
                    np.mean(target_score[target][:, PREDICTORS.index(predictor)]
                           - target_score[target][:, 0])
                ),
                "score_difference_median": float(
                    np.median(target_score[target][:, PREDICTORS.index(predictor)]
                             - target_score[target][:, 0])
                ),
                "score_difference_q1": float(
                    np.percentile(
                        target_score[target][:, PREDICTORS.index(predictor)]
                        - target_score[target][:, 0],
                        25,
                    )
                ),
                "score_difference_q3": float(
                    np.percentile(
                        target_score[target][:, PREDICTORS.index(predictor)]
                        - target_score[target][:, 0],
                        75,
                    )
                ),
                "score_mean_absolute_difference": float(
                    np.mean(
                        np.abs(
                            target_score[target][:, PREDICTORS.index(predictor)]
                            - target_score[target][:, 0]
                        )
                    )
                ),
                "spearman_vs_boltz2": correlation(
                    target_score[target][:, PREDICTORS.index(predictor)],
                    target_score[target][:, 0],
                    "spearman",
                ),
                "pearson_vs_boltz2": correlation(
                    target_score[target][:, PREDICTORS.index(predictor)],
                    target_score[target][:, 0],
                    "pearson",
                ),
            }
            for metric_name, estimate in score_specs.items():
                interval = bootstrap_interval(
                    score_bootstrap[(target, predictor, metric_name)], estimate
                )
                bootstrap_rows.append(
                    {
                        "outcome": "shared_main_and_strict",
                        "target": target,
                        "predictor": predictor,
                        "metric": metric_name,
                        **interval,
                        "validity_note": "correlation undefined when a resample is constant",
                    }
                )
            for outcome in OUTCOMES:
                for metric_name in ["auc", "ap", "brier", "log_loss"]:
                    delta_name = f"delta_{metric_name}_vs_boltz2"
                    estimate = float(
                        point_metrics[(outcome, target, predictor)][delta_name]
                    )
                    interval = bootstrap_interval(
                        bootstrap_values[(("outcome", outcome, target, predictor, delta_name))],
                        estimate,
                    )
                    bootstrap_rows.append(
                        {
                            "outcome": outcome,
                            "target": target,
                            "predictor": predictor,
                            "metric": delta_name,
                            **interval,
                            "validity_note": (
                                "AUC/AP undefined when the resample has only one outcome class; "
                                "those iterations are excluded, not set to zero"
                            )
                            if metric_name in {"auc", "ap"}
                            else "all paired resamples are defined",
                        }
                    )
    bootstrap = pd.DataFrame(bootstrap_rows)
    write_csv(RUN_DIR / "rq6_bootstrap_intervals.csv", bootstrap)

    log("Drawing score-shift and discrimination/error figures.")
    shift_target = score_by_target.pivot(
        index="target", columns="predictor", values="difference_median"
    ).reindex(index=TARGETS, columns=PREDICTORS)
    spearman_target = score_by_target.pivot(
        index="target", columns="predictor", values="spearman_vs_boltz2"
    ).reindex(index=TARGETS, columns=PREDICTORS)
    figure, axes = plt.subplots(1, 2, figsize=(18, 7.5), constrained_layout=True)
    image_shift = make_heatmap(
        axes[0],
        shift_target,
        "Median paired ipTM change vs Boltz-2",
        "RdBu_r",
        center=0.0,
        value_format=".3f",
    )
    figure.colorbar(image_shift, ax=axes[0], label="Predictor ipTM minus Boltz-2 ipTM")
    image_corr = make_heatmap(
        axes[1],
        spearman_target,
        "Within-target Spearman correlation with Boltz-2",
        "viridis",
        value_format=".3f",
    )
    figure.colorbar(image_corr, ax=axes[1], label="Spearman rho")
    figure.suptitle(
        "RQ6 score shifts on 720 common Anthropic designs (scores shared by both outcomes)",
        fontsize=14,
    )
    figure.savefig(RUN_DIR / "RQ6_figure1_score_shift.png", dpi=240)
    figure.savefig(RUN_DIR / "RQ6_figure1_score_shift.svg")
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(17, 11), constrained_layout=True)
    for row_index, outcome in enumerate(OUTCOMES):
        for column_index, metric_name in enumerate(["auc", "brier"]):
            matrix = (
                metrics.loc[metrics.outcome.eq(outcome)]
                .pivot(index="target", columns="predictor", values=metric_name)
                .reindex(index=TARGETS, columns=PREDICTORS)
            )
            axis = axes[row_index, column_index]
            image = make_heatmap(
                axis,
                matrix,
                f"{outcome}: {metric_name.upper() if metric_name == 'auc' else 'Raw Brier'}",
                "viridis_r" if metric_name == "brier" else "viridis",
                value_format=".3f",
            )
            figure.colorbar(image, ax=axis, label=metric_name)
    figure.suptitle(
        "RQ6 discrimination and raw probability error by target (same 90 designs per cell)",
        fontsize=14,
    )
    figure.savefig(RUN_DIR / "RQ6_figure2_discrimination_error.png", dpi=240)
    figure.savefig(RUN_DIR / "RQ6_figure2_discrimination_error.svg")
    plt.close(figure)

    log("Computing report findings and writing self-check documents.")
    score_equal = score_shift.loc[score_shift.summary_level.eq("target_equal")].copy()
    score_equal["absolute_target_equal_median_shift"] = score_equal.difference_median.abs()
    top_shift = score_equal.sort_values(
        ["absolute_target_equal_median_shift", "predictor"], ascending=[False, True]
    ).head(5)

    best_auc_main = rank_winners(metrics, "adaptyv_weak_plus", "auc", True)
    best_auc_strict = rank_winners(metrics, "adaptyv_medium_plus", "auc", True)
    best_brier_main = rank_winners(metrics, "adaptyv_weak_plus", "brier", False)
    best_brier_strict = rank_winners(metrics, "adaptyv_medium_plus", "brier", False)

    auc_compare = best_auc_main.merge(
        best_auc_strict,
        on="target",
        suffixes=("_weak_plus", "_medium_plus"),
        how="outer",
        validate="one_to_one",
    )
    auc_same = int(
        (auc_compare.predictor_weak_plus == auc_compare.predictor_medium_plus).sum()
    )
    brier_compare = best_brier_main.merge(
        best_brier_strict,
        on="target",
        suffixes=("_weak_plus", "_medium_plus"),
        how="outer",
        validate="one_to_one",
    )
    brier_same = int(
        (brier_compare.predictor_weak_plus == brier_compare.predictor_medium_plus).sum()
    )

    auc_unique_main = int(best_auc_main.predictor.nunique())
    auc_unique_strict = int(best_auc_strict.predictor.nunique())
    brier_unique_main = int(best_brier_main.predictor.nunique())
    brier_unique_strict = int(best_brier_strict.predictor.nunique())
    brier_equal = target_equal.loc[target_equal.metric.eq("brier")]
    brier_ranges = {
        outcome: (
            float(brier_equal.loc[brier_equal.outcome.eq(outcome), "target_equal_mean"].min()),
            float(brier_equal.loc[brier_equal.outcome.eq(outcome), "target_equal_mean"].max()),
        )
        for outcome in OUTCOMES
    }
    auc_equal = target_equal.loc[target_equal.metric.eq("auc")]
    auc_ranges = {
        outcome: (
            float(auc_equal.loc[auc_equal.outcome.eq(outcome), "target_equal_mean"].min()),
            float(auc_equal.loc[auc_equal.outcome.eq(outcome), "target_equal_mean"].max()),
        )
        for outcome in OUTCOMES
    }

    report_lines = [
        "# RQ6 最小描述性跨预测器结果",
        "",
        "**定位：Anthropic 外部队列内的同配对、探索性描述性补充。**本结果不是 Proteinbase 主队列内的跨预测器验证，也不检验固定映射迁移。",
        "",
        "## 方法与执行",
        "",
        "- 队列：BHRF1、IL-7Ra、RBX1、TrkA、Nipah-G、PD-L1、TREM2、VEGF-A；每靶点 90 个设计，共 720 个唯一 `uuid`。",
        "- 主结局 Adaptyv Weak+：720 个有效设计，287 个阳性；严格结局 Adaptyv Medium+：720 个有效设计，254 个阳性。",
        "- 10 个预测器均从原始 parquet 独立重建。仅取 `target_form` 为空的 design-target 记录，并要求化学计量等于该设计的 `boltz2_stoichiometry`。",
        "- 每个 `uuid × predictor` 从 5 个发布 runs/models 中选择 `ipsae_min` 最大者；并列取数值最小标签；读取该行的 `iptm_pae`。这些 runs/models 不视为独立预测器。",
        f"- 重建分数与可用性连接表最大绝对差为 `{availability_max_diff:.3e}`，与已审计 Boltz-2 分数最大绝对差为 `{boltz_max_diff:.3e}`，均不超过 `1e-6`。",
        f"- 不确定性：固定种子 `{BOOTSTRAP_SEED}`，每靶点 `{BOOTSTRAP_REPLICATES}` 次设计层配对 bootstrap；同一抽样索引同时用于 10 个分数与两个结局。区间为 percentile 95%；单一结局类别的 AUC/AP 迭代无效而非填零。",
        "- 未拟合 U／I／F，未做校准模型搜索、显著性筛选、联网、下载、安装依赖或新结构预测。",
        "",
        "## 主要发现",
        "",
        "### 1. 分数水平与设计排序",
        "",
        "下表为 8 靶点等权的预测器汇总；分数统计在主、严格结局间共用，没有重复计算为两套分数。",
        "",
        "| 预测器 | 等权平均分数 | 中位分数偏移 | 平均绝对差 | 平均 Spearman | 平均 Pearson | 跨靶点 Spearman 范围 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in score_equal.iterrows():
        predictor_range = score_ranges.loc[score_ranges.predictor.eq(row.predictor)].iloc[0]
        report_lines.append(
            f"| {row.predictor} | {fmt(row.score_mean)} | {fmt(row.difference_median)} | "
            f"{fmt(row.mean_absolute_difference)} | {fmt(row.spearman_vs_boltz2)} | "
            f"{fmt(row.pearson_vs_boltz2)} | "
            f"{fmt(predictor_range.spearman_min_across_targets)}–"
            f"{fmt(predictor_range.spearman_max_across_targets)} |"
        )
    report_lines.extend(
        [
            "",
            "按靶点等权中位偏移的绝对值，偏离 Boltz-2 最大的五个预测器为："
            + "、".join(
                f"`{row.predictor}`（{fmt(row.difference_median)}）"
                for _, row in top_shift.iterrows()
            )
            + "。正负号表示预测器 `iptm_pae` 高于或低于同一设计上的 Boltz-2，不代表真值方向。",
            "",
            "### 2. 靶点内区分能力",
            "",
            f"Weak+ 的 8 个靶点最佳 AUC 共涉及 {auc_unique_main} 个不同预测器；Medium+ 共涉及 {auc_unique_strict} 个。"
            f"等权平均 AUC 的预测器范围分别为 Weak+ `{fmt(auc_ranges['adaptyv_weak_plus'][0])}`–"
            f"`{fmt(auc_ranges['adaptyv_weak_plus'][1])}`、Medium+ `{fmt(auc_ranges['adaptyv_medium_plus'][0])}`–"
            f"`{fmt(auc_ranges['adaptyv_medium_plus'][1])}`。",
            "",
            "| 靶点 | Weak+ 最佳 AUC 预测器 | AUC | Medium+ 最佳 AUC 预测器 | AUC |",
            "|---|---|---:|---|---:|",
        ]
    )
    for _, row in auc_compare.iterrows():
        report_lines.append(
            f"| {row.target} | {row.predictor_weak_plus} | {fmt(row.auc_weak_plus)} | "
            f"{row.predictor_medium_plus} | {fmt(row.auc_medium_plus)} |"
        )
    report_lines.extend(
        [
            "",
            "最佳预测器身份随靶点改变；另外，AP 的逐靶点结果与数值区间见 `rq6_metrics_by_target.csv` 和 `rq6_bootstrap_intervals.csv`。",
            "",
            "### 3. 原始概率误差与结局一致性",
            "",
            f"Weak+ 的等权平均 Brier 在预测器间为 `{fmt(brier_ranges['adaptyv_weak_plus'][0])}`–"
            f"`{fmt(brier_ranges['adaptyv_weak_plus'][1])}`；Medium+ 为 "
            f"`{fmt(brier_ranges['adaptyv_medium_plus'][0])}`–"
            f"`{fmt(brier_ranges['adaptyv_medium_plus'][1])}`。"
            "这些是未经重新校准的 ipTM 原始误差，不应解释为校准能力。",
            "",
            f"Weak+ 的 8 个靶点最低 Brier 共涉及 {brier_unique_main} 个不同预测器；Medium+ 共涉及 {brier_unique_strict} 个。"
            f"两种结局的最佳 AUC 预测器在 {auc_same}/8 个靶点完全相同，最佳 Brier 预测器在 {brier_same}/8 个靶点完全相同。",
            "",
            "| 靶点 | Weak+ 最低 Brier 预测器 | Brier | Medium+ 最低 Brier 预测器 | Brier |",
            "|---|---|---:|---|---:|",
        ]
    )
    for _, row in brier_compare.iterrows():
        report_lines.append(
            f"| {row.target} | {row.predictor_weak_plus} | {fmt(row.brier_weak_plus)} | "
            f"{row.predictor_medium_plus} | {fmt(row.brier_medium_plus)} |"
        )
    report_lines.extend(
        [
            "",
            "### 4. 不确定性",
            "",
            "`rq6_bootstrap_intervals.csv` 对每个靶点×预测器报告 ΔAUC、ΔAP、ΔBrier、Δlog loss 和与 Boltz-2 的 Spearman／Pearson 相关，以及分数差值的均值、中位数、四分位与平均绝对差 95% 区间。"
            "AUC/AP 的有效次数逐单元格记录；本次完整 720 行点估计均有双结局类别，但单个 bootstrap 重采样可能出现单一类别。",
            "",
            "## 图表",
            "",
            "- **图 1**：`RQ6_figure1_score_shift.png` / `.svg`。逐靶点中位分数偏移及 Spearman 相关热图。",
            "- **图 2**：`RQ6_figure2_discrimination_error.png` / `.svg`。两种结局的逐靶点 AUC 与原始 Brier 摘要热图。",
            "",
            "## 解释边界",
            "",
            "- 结论仅适用于 Anthropic 8 个靶点、经上游计算筛选的同一批 720 个设计及其 Adaptyv 结局。",
            "- 不同代码代表完整共折叠预测器管线；差异不能仅归因于模型权重，MSA、recycles、采样与实现也共同变化。",
            "- 本分析没有证明 Proteinbase 主队列中的跨模型迁移，没有解决 RQ3 的跨靶点失败，也不识别“普遍最佳”预测器。",
            "- `af3of3` 是 AlphaFold 3 架构代码配合 OpenFold3 权重，不是 AlphaFold 3 权重结果。",
            "- 主结局与严格结局的绝对概率误差对应不同基率，跨结局绝对数值不作直接优劣比较。",
            "",
            "正式运行目录：`" + str(RUN_DIR) + "`。",
        ]
    )
    write_text(RUN_DIR / "RQ6_results.md", "\n".join(report_lines) + "\n")

    self_check_lines = [
        "# RQ6 正式运行自检报告",
        "",
        "| 核查项 | 结果 | 证据 |",
        "|---|---|---|",
        "| 协议与执行入口 | 通过 | 三份冻结规范已完整读取；脚本按冻结范围执行 |",
        f"| 核心输入 SHA-256 | 通过 | parquet `{EXPECTED_HASHES[COFOLD]}`；pair table `{EXPECTED_HASHES[PAIR_TABLE]}`；主表只核验 `{EXPECTED_HASHES[MASTER]}` |",
        "| 样本 | 通过 | 720 个唯一 `uuid`；8 靶点各 90 |",
        "| 主结局 | 通过 | Adaptyv Weak+：720 有效，287 阳性 |",
        "| 严格结局 | 通过 | Adaptyv Medium+：720 有效，254 阳性 |",
        "| 预测器 | 通过 | 固定 10 个；无缺失、无额外预测器 |",
        "| design-target 记录 | 通过 | 只取 `target_form` 为空的记录 |",
        "| 化学计量 | 通过 | 每条记录匹配置计表中该 `uuid` 的 `boltz2_stoichiometry` |",
        "| 运行/模型选择 | 通过 | 每 `uuid × predictor` 恰有 5 行；`ipsae_min` 最大、并列取最小数值标签 |",
        "| 分数取值 | 通过 | 读取所选行 `iptm_pae`；未使用 max-ipTM、seed 均值或中位数 |",
        "| 长/宽表完整性 | 通过 | 720×10 wide；7,200 行 long；无缺、无重复 |",
        f"| 可用性逐值 QC | 通过 | 最大绝对差 `{availability_max_diff:.12e}` ≤ `1e-6` |",
        f"| Boltz-2 审计 QC | 通过 | 最大绝对差 `{boltz_max_diff:.12e}` ≤ `1e-6` |",
        "| 指标定义 | 通过 | AUC 为 Mann–Whitney/平均秩等价；AP 用 sklearn；Brier、裁剪 `[1e-6,1-1e-6]` log loss |",
        "| 配对差值 | 通过 | 所有 Δ 都为 predictor 减同一 `uuid` 的 Boltz-2 |",
        "| Bootstrap | 通过 | 每靶点 2,000 次；固定种子 20261115；每个 draw 同步复用于 10 分数和 2 结局 |",
        "| 单类重采样 | 通过 | AUC/AP 单一类别迭代不填 0，记录 `valid_iterations` |",
        "| 主/严格结局 | 通过 | 两结局全部预定义指标均计算；分数统计共享且未伪造成两套 |",
        "| 靶点等权摘要 | 通过 | 每个结局×预测器×指标给出 8 靶点 mean、median、min、max |",
        "| 图表 | 通过 | 两张图均交付 PNG 和 SVG |",
        "| 解释边界 | 通过 | 报告明确外部队列、上游筛选、管线差异与不可外推性 |",
        "| 禁止操作 | 通过 | 未联网、下载、安装、新预测、模型拟合、映射迁移、显著性筛选或修改既有结果 |",
        "| 冻结主表稳定 | 通过 | 执行结束前再次核验 SHA-256 |",
        "",
        "最终文件哈希见 `final_hash_manifest.csv`；该清单按惯例不包含自身哈希，`run_manifest.json` 记录输入与运行元数据。",
    ]
    write_text(RUN_DIR / "RQ6_self_check.md", "\n".join(self_check_lines) + "\n")

    environment_lines = [
        "# RQ6 执行环境",
        "",
        f"- Python: `{platform.python_version()}`",
        f"- NumPy: `{np.__version__}`",
        f"- pandas: `{pd.__version__}`",
        f"- SciPy: `{scipy.__version__}`",
        f"- scikit-learn: `{__import__('sklearn').__version__}`",
        f"- DuckDB: `{duckdb.__version__}`",
        f"- Matplotlib: `{matplotlib.__version__}`",
        "- Parquet 读取：现有 DuckDB；未安装 `pyarrow` 或其他依赖。",
        "- 网络：未使用。",
        "",
        "所有依赖均为运行前环境中已有版本。",
    ]
    write_text(RUN_DIR / "requirements_or_environment.md", "\n".join(environment_lines) + "\n")

    command_lines = [
        "# RQ6 正式执行命令",
        "",
        "工作目录：`E:\\diyipian`",
        "",
        "```powershell",
        "python \"analysis/codex_outputs/RQ6_execution_2026-11-15/"
        + RUN_DIR.name
        + "/rq6_analysis.py\"",
        "```",
        "",
        "脚本使用冻结种子与完整样本；不接受结果驱动参数，不联网，不安装依赖。",
    ]
    write_text(RUN_DIR / "COMMANDS.md", "\n".join(command_lines) + "\n")

    input_hashes_after = {str(path): sha256(path) for path in input_paths}
    changed_inputs = {
        path: {"before": input_hashes_before[path], "after": input_hashes_after[path]}
        for path in input_hashes_before
        if input_hashes_before[path] != input_hashes_after[path]
    }
    if changed_inputs:
        raise RuntimeError(f"A frozen input changed during execution: {changed_inputs}")

    qc = {
        "status": "passed",
        "pair_rows": len(wide),
        "unique_uuids": int(wide.uuid.nunique()),
        "targets": TARGETS,
        "target_sizes": {target: int((wide.target == target).sum()) for target in TARGETS},
        "predictors": PREDICTORS,
        "long_rows": len(long),
        "selected_run_rows": len(selected_audit),
        "five_runs_or_models_per_uuid_predictor": True,
        "outcome_totals": {
            "adaptyv_weak_plus": {"n": 720, "positives": 287},
            "adaptyv_medium_plus": {"n": 720, "positives": 254},
        },
        "outcome_by_target_positives": {
            target: {
                "adaptyv_weak_plus": EXPECTED_POSITIVES[target][0],
                "adaptyv_medium_plus": EXPECTED_POSITIVES[target][1],
            }
            for target in TARGETS
        },
        "availability_max_abs_diff": availability_max_diff,
        "boltz2_published_max_abs_diff": boltz_max_diff,
        "tolerance": PAIR_TOLERANCE,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates_per_target": BOOTSTRAP_REPLICATES,
        "bootstrap_rng": "numpy.random.default_rng(20261115); TARGETS order; 2,000 draws each",
        "metrics_rows": len(metrics),
        "bootstrap_interval_rows": len(bootstrap),
        "single_class_auc_ap_filled_with_zero": False,
        "score_statistics_duplicated_by_outcome": False,
        "fitted_models": [],
        "network_used": False,
        "downloads_or_installs": False,
        "new_structure_predictions": False,
        "protected_inputs_unchanged_during_run": True,
        "input_hashes_before": input_hashes_before,
        "input_hashes_after": input_hashes_after,
        "software": platform.platform(),
    }
    write_json(RUN_DIR / "rq6_qc.json", qc)

    manifest = {
        "status": "complete",
        "scope": "frozen_minimal_RQ6_descriptive_cross_predictor_analysis",
        "run_directory": str(RUN_DIR),
        "protocol": str(PROTOCOL),
        "protocol_sha256": input_hashes_before[str(PROTOCOL)],
        "execution_script": str(Path(__file__).resolve()),
        "execution_script_sha256": sha256(Path(__file__).resolve()),
        "inputs": input_hashes_before,
        "targets": TARGETS,
        "predictors": PREDICTORS,
        "outcomes": OUTCOMES,
        "bootstrap": {
            "seed": BOOTSTRAP_SEED,
            "replicates_per_target": BOOTSTRAP_REPLICATES,
            "design_layer_paired": True,
            "shared_draw_across_predictors_and_outcomes": True,
        },
        "qc": "passed",
        "interpretation_boundary": (
            "Anthropic 8-target, upstream-filtered 720-design paired descriptive "
            "comparison only; not Proteinbase transport or general model validation"
        ),
    }
    write_json(RUN_DIR / "run_manifest.json", manifest)

    handoff = {
        "status": "complete",
        "run_directory": str(RUN_DIR),
        "formal_result": "RQ6_results.md",
        "pair_wide": "rq6_pair_scores_wide.csv",
        "pair_long": "rq6_pair_scores_long.csv",
        "metrics": "rq6_metrics_by_target.csv",
        "bootstrap": "rq6_bootstrap_intervals.csv",
        "target_equal_summary": "rq6_target_equal_summary.csv",
        "figures": [
            "RQ6_figure1_score_shift.png",
            "RQ6_figure1_score_shift.svg",
            "RQ6_figure2_discrimination_error.png",
            "RQ6_figure2_discrimination_error.svg",
        ],
        "key_counts": {
            "designs": 720,
            "targets": 8,
            "predictors": 10,
            "weak_plus_positives": 287,
            "medium_plus_positives": 254,
            "bootstrap_replicates_per_target": 2000,
        },
        "protected_inputs_unchanged": True,
    }
    write_json(RUN_DIR / "handoff_result.json", handoff)

    log("All calculations and documents completed with QC passed.")
    write_json(
        RUN_DIR / "RUN_STATUS.json",
        {
            "status": "complete",
            "qc": "passed",
            "run_dir": str(RUN_DIR),
            "scope": "frozen_minimal_RQ6",
        },
    )
    log("Writing final SHA-256 manifest.")
    manifest_rows = []
    for path in sorted(RUN_DIR.iterdir(), key=lambda item: item.name):
        if path.is_file() and path.name != "final_hash_manifest.csv":
            manifest_rows.append(
                {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            )
    write_csv(RUN_DIR / "final_hash_manifest.csv", pd.DataFrame(manifest_rows))
    print(f"Complete with QC passed: {RUN_DIR}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        traceback_text = traceback.format_exc()
        (RUN_DIR / "failure_traceback.txt").write_text(traceback_text, encoding="utf-8")
        try:
            write_json(
                RUN_DIR / "RUN_STATUS.json",
                {"status": "failed", "error": str(error), "run_dir": str(RUN_DIR)},
            )
        except Exception:
            pass
        raise
