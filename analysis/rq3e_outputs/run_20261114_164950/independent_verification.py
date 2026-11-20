"""Independent numerical verification of the frozen RQ3-E 720-pair evaluation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = ROOT / "analysis/rq3e_outputs/run_20261114_164950"
RAW = ROOT / "data/external_target_audit_2026-11-13/anthropic_audited_pair_table.csv"
COFOLD = ROOT / "data/external_target_audit_2026-11-13/anthropic_cofold_predictions.parquet"
MASTER = ROOT / "analysis/v3_outputs/pair_master_v3.csv"
UPSTREAM = ROOT / "analysis/rq3e_outputs/run_20261114_162656/rq3e_source_model.json"
CLIP = 1e-6
EDGES = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
NEW_TARGETS = {"BHRF1", "IL-7Ra", "RBX1", "TrkA"}
DIFFS: dict[str, float] = {}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def close(name: str, actual, expected, rtol=2e-7, atol=2e-9) -> None:
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    valid = np.isfinite(a) | np.isfinite(e)
    if valid.any():
        diff = float(np.nanmax(np.abs(a[valid] - e[valid])))
        DIFFS[name] = max(DIFFS.get(name, 0.0), diff)
    np.testing.assert_allclose(a, e, rtol=rtol, atol=atol, equal_nan=True, err_msg=name)


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1"])


def logit(p):
    p = np.asarray(p, dtype=float)
    return np.log(p) - np.log1p(-p)


def fit_logistic(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Independently implemented unpenalized logistic Newton/IRLS fit."""
    theta = np.zeros(X.shape[1], dtype=float)
    theta[0] = float(logit(np.clip(np.mean(y), 1e-8, 1 - 1e-8)))

    def objective(t):
        eta = X @ t
        return float(np.sum(y * eta - np.logaddexp(0.0, eta)))

    current = objective(theta)
    converged = False
    for _ in range(200):
        p = expit(X @ theta)
        gradient = X.T @ (y - p)
        information = X.T @ ((p * (1 - p))[:, None] * X)
        step = np.linalg.solve(information, gradient)
        scale = 1.0
        accepted = False
        for _ in range(60):
            candidate = theta + scale * step
            candidate_objective = objective(candidate)
            if np.isfinite(candidate_objective) and candidate_objective >= current - 1e-12:
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            raise RuntimeError("Logistic line search failed")
        change = float(np.max(np.abs(candidate - theta)))
        theta = candidate
        current = candidate_objective
        if change < 2e-11:
            converged = True
            break
    p = expit(X @ theta)
    information = X.T @ ((p * (1 - p))[:, None] * X)
    covariance = np.linalg.inv(information)
    gradient_max = float(np.max(np.abs(X.T @ (y - p))))
    if not converged or gradient_max > 2e-5:
        raise RuntimeError(f"Logistic fit did not meet QC: converged={converged}, gradient={gradient_max}")
    return theta, covariance, gradient_max


def calculate(g: pd.DataFrame) -> dict[str, float]:
    y = g.outcome_value.to_numpy(dtype=float)
    p = g.predicted_probability.to_numpy(dtype=float)
    x = logit(p)
    n = len(g)
    positives = int(y.sum())
    prevalence = float(y.mean())
    ref_t = float(g.bss_transport_reference.iloc[0])
    ref_o_values = g.bss_target_oracle_reference.to_numpy(dtype=float)
    ref_o = float(ref_o_values[0])
    brier = float(np.mean((p - y) ** 2))
    log_loss = float(-np.mean(y * np.log(p) + (1 - y) * np.log1p(-p)))
    brier_t = float(np.mean((ref_t - y) ** 2))
    # Pooled summaries retain each target's own oracle prevalence, not one pooled constant.
    brier_o = float(np.mean((ref_o_values - y) ** 2))
    X = np.column_stack([np.ones(n), x])
    theta, covariance, gradient = fit_logistic(X, y)
    return {
        "n": n,
        "positives": positives,
        "prevalence": prevalence,
        "mean_prediction": float(np.mean(p)),
        "median_prediction": float(np.median(p)),
        "brier": brier,
        "log_loss": log_loss,
        "bss_transport": 1 - brier / brier_t,
        "bss_target_oracle": 1 - brier / brier_o,
        "transport_reference_value": ref_t,
        "target_oracle_reference_value": ref_o,
        "brier_transport_reference": brier_t,
        "brier_target_oracle_reference": brier_o,
        "calibration_intercept": float(theta[0]),
        "calibration_intercept_se": float(np.sqrt(covariance[0, 0])),
        "calibration_slope": float(theta[1]),
        "calibration_slope_se": float(np.sqrt(covariance[1, 1])),
        "calibration_gradient_max": gradient,
    }


def main() -> None:
    config = json.loads((RUN / "run_config.snapshot.json").read_text(encoding="utf-8"))
    target_model = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    pair_model = json.loads((RUN / "rq3e_source_model_pair_sensitivity.json").read_text(encoding="utf-8"))
    assert sha256(UPSTREAM) == config["upstream_target_equal_model_sha256"]
    assert sha256(MASTER) == config["input_hashes"][str(MASTER)]
    assert sha256(ROOT / "analysis/RQ3E_external_transport_frozen_protocol_2026-11-14.md") == config["protocol_sha256"]
    assert sha256(RUN / "rq3e_source_model_target_equal_reused.json") == config["upstream_target_equal_model_sha256"]

    pairs = pd.read_csv(RUN / "rq3e_external_pairs.csv")
    assert len(pairs) == 720 and pairs.uuid.nunique() == 720
    raw = pd.read_csv(RAW, low_memory=False)
    selected = raw.loc[
        raw.target.isin(config["targets"])
        & as_bool(raw.integrated_valid)
        & pd.to_numeric(raw.boltz2_iptm_ipsae_best_seed, errors="coerce").notna()
    ].copy()
    assert len(selected) == 720 and set(selected.uuid) == set(pairs.uuid)
    raw_by_id = selected.set_index("uuid").loc[pairs.uuid]
    for score_rule, column in config["score_columns"].items():
        if score_rule != "seed0":
            close(f"raw_score_{score_rule}", pairs[column], raw_by_id[column].to_numpy(), rtol=0, atol=0)
    for outcome, spec in config["outcomes"].items():
        assert np.array_equal(as_bool(pairs[f"valid__{outcome}"]), as_bool(raw_by_id[spec["valid_column"]]).to_numpy())
        assert np.array_equal(as_bool(pairs[f"outcome__{outcome}"]), as_bool(raw_by_id[spec["positive_column"]]).to_numpy())

    seed0 = duckdb.connect().execute(
        """SELECT uuid, iptm_pae FROM read_parquet(?)
           WHERE cofolding_model = 'boltz2'
             AND (target_form IS NULL OR target_form = '')
             AND stoichiometry IN ('1to1', '1to2', '1to3', 'rnp')
             AND CAST(seed AS INTEGER) = 0""",
        [str(COFOLD)],
    ).df()
    assert not seed0.uuid.duplicated().any()
    expected_seed0 = pairs[["uuid"]].merge(seed0, on="uuid", how="left", validate="one_to_one")
    close("raw_seed0", pairs.boltz2_iptm_seed0, expected_seed0.iptm_pae, rtol=0, atol=2e-15)

    pred = pd.read_csv(RUN / "rq3e_u_predictions_frozen.csv")
    scored = pd.read_csv(RUN / "rq3e_u_pair_predictions.csv")
    assert len(pred) == len(scored) == 6476
    common = [c for c in pred.columns if c in scored.columns]
    pd.testing.assert_frame_equal(pred[common], scored[common], check_exact=True)
    pairs_by_id = pairs.set_index("uuid")
    for cfg in config["configurations"]:
        cid = cfg["config_id"]
        idx = pred.config_id.eq(cid)
        block = pred.loc[idx]
        source = pairs_by_id.loc[block.uuid]
        assert np.array_equal(block.target.to_numpy(), source.target.to_numpy())
        assert as_bool(source[f"valid__{cfg['outcome']}"]).all()
        score = source[config["score_columns"][cfg["score"]]].to_numpy(dtype=float)
        model = target_model if cfg["source_weighting"] == "target_equal" else pair_model
        expected_x = logit(score)
        expected_eta = model["transport_intercept"] + model["beta_common"] * expected_x
        expected_p = np.clip(expit(expected_eta), CLIP, 1 - CLIP)
        close(f"score_{cid}", block.score_value, score, rtol=0, atol=0)
        close(f"x_{cid}", block.x, expected_x)
        close(f"eta_{cid}", block.eta_U, expected_eta)
        close(f"probability_{cid}", block.predicted_probability, expected_p)
        close(f"transport_reference_{cid}", block.bss_transport_reference, np.full(len(block), model["transport_reference"]), rtol=0, atol=2e-15)
        expected_y = as_bool(source[f"outcome__{cfg['outcome']}"]).to_numpy(dtype=int)
        assert np.array_equal(scored.loc[idx, "outcome_value"].to_numpy(dtype=int), expected_y)
        oracle = pd.Series(expected_y).groupby(block.target.to_numpy()).mean()
        close(f"oracle_{cid}", block.bss_target_oracle_reference, block.target.map(oracle), rtol=0, atol=2e-15)

    recalculated = []
    for (cid, target), g in scored.groupby(["config_id", "target"], sort=False):
        row = {"config_id": cid, "target": target, "stratum": g.stratum.iloc[0], **calculate(g)}
        recalculated.append(row)
    recalc = pd.DataFrame(recalculated)
    assert len(recalc) == 72
    observed = pd.read_csv(RUN / "rq3e_target_metrics.csv")
    merged = observed.merge(recalc, on=["config_id", "target"], suffixes=("_observed", "_expected"), validate="one_to_one")
    assert len(merged) == 72
    metric_columns = [
        "n", "positives", "prevalence", "mean_prediction", "median_prediction", "brier", "log_loss",
        "bss_transport", "bss_target_oracle", "transport_reference_value", "target_oracle_reference_value",
        "brier_transport_reference", "brier_target_oracle_reference", "calibration_intercept",
        "calibration_intercept_se", "calibration_slope", "calibration_slope_se",
    ]
    for column in metric_columns:
        close(f"target_{column}", merged[f"{column}_observed"], merged[f"{column}_expected"], rtol=5e-6, atol=2e-8)

    recalc_groups = []
    for cfg in config["configurations"]:
        cid = cfg["config_id"]
        for stratum in ["all_8", "new_target_4", "same_target_campaign_4"]:
            if stratum == "all_8":
                subset = recalc.loc[recalc.config_id.eq(cid)]
            elif stratum == "new_target_4":
                subset = recalc.loc[recalc.config_id.eq(cid) & recalc.target.isin(NEW_TARGETS)]
            else:
                subset = recalc.loc[recalc.config_id.eq(cid) & ~recalc.target.isin(NEW_TARGETS)]
            row = {
                "config_id": cid, "aggregation": "target_equal_macro", "stratum": stratum,
                "n_targets": len(subset), "n_pairs": int(subset.n.sum()), "positives": int(subset.positives.sum()),
            }
            for column in ["prevalence", "brier", "log_loss", "bss_transport", "bss_target_oracle", "calibration_intercept", "calibration_slope"]:
                row[column] = float(subset[column].mean())
                row[f"{column}_min"] = float(subset[column].min())
                row[f"{column}_max"] = float(subset[column].max())
            recalc_groups.append(row)
        pooled = calculate(scored.loc[scored.config_id.eq(cid)])
        recalc_groups.append({
            "config_id": cid, "aggregation": "pair_pooled", "stratum": "all_available_pairs",
            "n_targets": int(scored.loc[scored.config_id.eq(cid), "target"].nunique()),
            "n_pairs": pooled["n"], "positives": pooled["positives"], **pooled,
        })
    recalc_group = pd.DataFrame(recalc_groups)
    observed_group = pd.read_csv(RUN / "rq3e_group_metrics.csv")
    gm = observed_group.merge(
        recalc_group, on=["config_id", "aggregation", "stratum"], suffixes=("_observed", "_expected"), validate="one_to_one"
    )
    assert len(gm) == 36
    group_columns = [
        "n_targets", "n_pairs", "positives", "prevalence", "brier", "brier_min", "brier_max", "log_loss",
        "log_loss_min", "log_loss_max", "bss_transport", "bss_transport_min", "bss_transport_max",
        "bss_target_oracle", "bss_target_oracle_min", "bss_target_oracle_max", "calibration_intercept",
        "calibration_intercept_min", "calibration_intercept_max", "calibration_slope", "calibration_slope_min",
        "calibration_slope_max", "calibration_intercept_se", "calibration_slope_se",
    ]
    for column in group_columns:
        close(f"group_{column}", gm[f"{column}_observed"], gm[f"{column}_expected"], rtol=5e-6, atol=2e-8)

    reliability = pd.read_csv(RUN / "rq3e_reliability.csv")
    expected_reliability = []
    z = float(norm.ppf(0.975))
    for (cid, target), g in scored.groupby(["config_id", "target"], sort=False):
        p = g.predicted_probability.to_numpy(dtype=float)
        y = g.outcome_value.to_numpy(dtype=int)
        bins = np.minimum(np.searchsorted(EDGES, p, side="right") - 1, len(EDGES) - 2)
        for b in range(len(EDGES) - 1):
            take = bins == b
            n = int(take.sum())
            pos = int(y[take].sum())
            if n:
                phat = pos / n
                denominator = 1 + z * z / n
                center = (phat + z * z / (2 * n)) / denominator
                half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denominator
                mean_p = float(p[take].mean())
                low, high = float(center - half), float(center + half)
                observed_rate = float(phat)
            else:
                mean_p = observed_rate = low = high = np.nan
            expected_reliability.append({
                "config_id": cid, "target": target, "bin_index": b, "bin_lower": EDGES[b], "bin_upper": EDGES[b + 1],
                "n": n, "positives": pos, "mean_prediction": mean_p, "observed_rate": observed_rate,
                "wilson_95_low": low, "wilson_95_high": high,
            })
    er = pd.DataFrame(expected_reliability)
    rm = reliability.merge(er, on=["config_id", "target", "bin_index"], suffixes=("_observed", "_expected"), validate="one_to_one")
    assert len(rm) == 504
    for column in ["bin_lower", "bin_upper", "n", "positives", "mean_prediction", "observed_rate", "wilson_95_low", "wilson_95_high"]:
        close(f"reliability_{column}", rm[f"{column}_observed"], rm[f"{column}_expected"], rtol=2e-7, atol=2e-9)

    source = pd.read_csv(MASTER)
    levels = pair_model["levels"]
    X = np.column_stack([
        *[source.target.eq(level).to_numpy(dtype=float) for level in levels],
        logit(source.boltz2_iptm.to_numpy(dtype=float)),
    ])
    y = source.primary_outcome.to_numpy(dtype=float)
    theta, covariance, gradient = fit_logistic(X, y)
    counts = source.target.value_counts()
    transport = float(np.dot(theta[:6], counts.loc[levels].to_numpy(dtype=float)) / len(source))
    close("source_pair_alpha", theta[:6], pair_model["alpha"], rtol=2e-7, atol=1e-6)
    close("source_pair_beta", theta[6], pair_model["beta_common"], rtol=2e-7, atol=2e-9)
    close("source_pair_transport", transport, pair_model["transport_intercept"], rtol=2e-7, atol=2e-9)
    close("source_pair_reference", y.mean(), pair_model["transport_reference"], rtol=0, atol=0)
    close("source_pair_se", np.sqrt(np.diag(covariance)), pair_model["alpha_se"] + [pair_model["beta_common_se"]], rtol=5e-6, atol=2e-8) if "alpha_se" in pair_model else None

    manifest = pd.read_csv(RUN / "final_hash_manifest.csv")
    for row in manifest.itertuples(index=False):
        path = RUN / row.file
        assert path.stat().st_size == row.bytes
        assert sha256(path) == row.sha256

    report = {
        "status": "passed",
        "scope": "independent_recalculation_without_using_original_metric_functions",
        "external_pairs_verified": len(pairs),
        "prediction_rows_verified": len(pred),
        "target_metric_rows_verified": len(merged),
        "group_metric_rows_verified": len(gm),
        "reliability_rows_verified": len(rm),
        "manifest_files_verified": len(manifest),
        "source_pair_rows": len(source),
        "source_pair_gradient_max_abs": gradient,
        "maximum_absolute_differences": DIFFS,
        "target_equal_model_reused_unchanged": True,
        "external_labels_used_for_fitting": False,
        "master_sha256": sha256(MASTER),
        "protocol_sha256": sha256(ROOT / "analysis/RQ3E_external_transport_frozen_protocol_2026-11-14.md"),
    }
    (RUN / "independent_verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
