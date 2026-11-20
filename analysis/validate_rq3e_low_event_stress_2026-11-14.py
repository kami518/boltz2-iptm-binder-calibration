"""Independent numerical validation of the completed RQ3-E stress-test run."""

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = ROOT / "analysis/rq3e_outputs/run_20261114_162656"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
EXPECTED_MASTER = "581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91"
CLIP = 1e-6


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), CLIP, 1.0 - CLIP)
    return np.log(p) - np.log1p(-p)


def bss(y: np.ndarray, p: np.ndarray, r: np.ndarray) -> float:
    denominator = float(np.sum((r - y) ** 2))
    return float("nan") if denominator == 0 else float(1.0 - np.sum((p - y) ** 2) / denominator)


def main() -> None:
    source = pd.read_csv(MASTER, low_memory=False)
    model = json.loads((RUN / "rq3e_source_model.json").read_text(encoding="utf-8"))
    levels = model["levels"]
    targets = source.target.astype(str)
    one_hot = np.column_stack([(targets == level).to_numpy(dtype=float) for level in levels])
    x = logit(source.boltz2_iptm.to_numpy(dtype=float))
    design = np.column_stack([one_hot, x])
    counts = targets.value_counts()
    raw_w = targets.map(1.0 / counts).to_numpy(dtype=float)
    w = raw_w / raw_w.mean()
    independent_source = LogisticRegression(
        penalty=None, fit_intercept=False, solver="lbfgs", max_iter=10000, tol=1e-12
    )
    independent_source.fit(design, source.primary_outcome.to_numpy(dtype=int), sample_weight=w)
    expected_theta = np.asarray(model["alpha"] + [model["beta_common"]], dtype=float)
    source_difference = np.abs(independent_source.coef_[0] - expected_theta)
    assert float(source_difference.max()) < 1e-5

    pairs = pd.read_csv(RUN / "rq3e_low_event_stress_test_pairs.csv", low_memory=False)
    assert len(pairs) == 521 and pairs.uuid.nunique() == 521
    assert int(pairs.outcome_value.sum()) == 34
    p_recomputed = expit(model["transport_intercept"] + model["beta_common"] * logit(pairs.score.to_numpy(dtype=float)))
    prediction_difference = float(np.max(np.abs(p_recomputed - pairs.predicted_probability.to_numpy(dtype=float))))
    assert prediction_difference < 1e-12

    metrics = pd.read_csv(RUN / "rq3e_low_event_stress_test_metrics.csv").set_index("target")
    max_metric_difference = 0.0
    max_calibration_difference = 0.0
    calibration_differences = {}
    for target, group in pairs.groupby("target"):
        y = group.outcome_value.to_numpy(dtype=float)
        p = group.predicted_probability.to_numpy(dtype=float)
        transport = np.full(len(group), model["transport_reference"], dtype=float)
        oracle = np.full(len(group), float(y.mean()), dtype=float)
        recomputed = {
            "brier": float(np.mean((p - y) ** 2)),
            "log_loss": float(-np.mean(y * np.log(np.clip(p, CLIP, 1 - CLIP)) + (1 - y) * np.log1p(-np.clip(p, CLIP, 1 - CLIP)))),
            "bss_transport": bss(y, p, transport),
            "bss_target_oracle": bss(y, p, oracle),
        }
        for name, value in recomputed.items():
            stored = float(metrics.loc[target, name])
            if np.isnan(value):
                assert np.isnan(stored)
            else:
                difference = abs(value - stored)
                max_metric_difference = max(max_metric_difference, difference)
                assert difference < 1e-12
        if len(np.unique(y)) == 2:
            calibration = LogisticRegression(
                penalty=None, fit_intercept=True, solver="lbfgs", max_iter=10000, tol=1e-12
            )
            calibration.fit(logit(p).reshape(-1, 1), y.astype(int))
            intercept = float(calibration.intercept_[0])
            slope = float(calibration.coef_[0, 0])
            intercept_difference = abs(intercept - float(metrics.loc[target, "calibration_intercept"]))
            slope_difference = abs(slope - float(metrics.loc[target, "calibration_slope"]))
            calibration_differences[target] = {"intercept": intercept_difference, "slope": slope_difference}
            max_calibration_difference = max(max_calibration_difference, intercept_difference, slope_difference)
            assert intercept_difference < 1e-4 and slope_difference < 1e-4
        else:
            assert np.isnan(float(metrics.loc[target, "calibration_intercept"]))
            assert np.isnan(float(metrics.loc[target, "calibration_slope"]))

    reliability = pd.read_csv(RUN / "rq3e_low_event_stress_test_reliability.csv")
    assert len(reliability) == 42
    for target, group in pairs.groupby("target"):
        assert int(reliability.loc[reliability.target == target, "n"].sum()) == len(group)
        assert int(reliability.loc[reliability.target == target, "positives"].sum()) == int(group.outcome_value.sum())

    manifest = pd.read_csv(RUN / "final_hash_manifest.csv")
    manifest_mismatches = []
    for row in manifest.itertuples(index=False):
        path = RUN / row.file
        if not path.exists() or sha(path) != row.sha256:
            manifest_mismatches.append(row.file)
    assert not manifest_mismatches
    assert sha(MASTER) == EXPECTED_MASTER
    status = json.loads((RUN / "RUN_STATUS.json").read_text(encoding="utf-8"))
    assert status["status"] == "complete" and status["confirmatory_720_evaluated"] is False

    report = {
        "status": "passed",
        "independent_source_implementation": "sklearn unpenalized weighted logistic regression",
        "source_max_absolute_coefficient_difference": float(source_difference.max()),
        "external_max_absolute_prediction_difference": prediction_difference,
        "max_absolute_brier_logloss_bss_difference": max_metric_difference,
        "max_absolute_calibration_intercept_or_slope_difference": max_calibration_difference,
        "calibration_differences_by_target": calibration_differences,
        "pairs": len(pairs),
        "positives": int(pairs.outcome_value.sum()),
        "reliability_rows": len(reliability),
        "final_manifest_files_verified": len(manifest),
        "master_hash_unchanged": True,
        "confirmatory_720_evaluated": False,
    }
    (RUN / "independent_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    shutil.copy2(__file__, RUN / "independent_validation.py")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
