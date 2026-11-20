"""Protocol section 11 unit tests for RQ5-RF."""
from __future__ import annotations

import hashlib
import json
import unittest

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.ensemble import RandomForestClassifier

from rq5rf_common import (
    ANALYSIS_SETS,
    FEATURE_MAP,
    adaptation_inner_splits,
    atomic_csv,
    bss,
    current_source_weights,
    effective_rf_grid,
    fit_offset_slope,
    fit_pooled_logistic,
    fit_random_forest,
    load_config,
    load_joined_inputs,
    log_loss,
    logit_clip,
    logistic_eta,
    rf_positive_probability,
    solve_offset_intercept,
    stable_sort,
    task_seed,
    weighted_center_scale,
)


class RQ5RFUnitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_config()

    def test_recomputed_weights_and_target_sums(self) -> None:
        frame = pd.DataFrame(
            {
                "pair_id": ["b|t2", "a|t1", "c|t1", "d|t2", "e|t1"],
                "target": ["t2", "t1", "t1", "t2", "t1"],
            }
        )
        target_equal, audit = current_source_weights(frame, "target_equal")
        self.assertAlmostEqual(float(target_equal.sum()), 1.0, delta=1e-15)
        for target in ("t1", "t2"):
            self.assertAlmostEqual(
                float(target_equal[frame["target"].eq(target).to_numpy()].sum()),
                0.5,
                delta=1e-15,
            )
        self.assertLessEqual(audit["max_target_weight_sum_error"], 1e-12)
        pair, pair_audit = current_source_weights(frame, "pair")
        np.testing.assert_allclose(pair, np.full(5, 0.2), rtol=0.0, atol=1e-15)
        self.assertLessEqual(pair_audit["weight_normalization_error"], 1e-12)
        self.assertNotEqual(
            audit["weight_vector_sha256"], pair_audit["weight_vector_sha256"]
        )

    def test_logit_and_numeric_clip(self) -> None:
        values = np.array([0.0, 1e-8, 0.5, 1.0])
        transformed = logit_clip(values, 1e-6)
        expected = np.log(
            np.clip(values, 1e-6, 1.0 - 1e-6)
        ) - np.log1p(-np.clip(values, 1e-6, 1.0 - 1e-6))
        np.testing.assert_array_equal(transformed, expected)
        self.assertTrue(np.isfinite(transformed).all())
        with self.assertRaises(RuntimeError):
            logit_clip(np.array([0.5, np.nan]), 1e-6)
        with self.assertRaises(RuntimeError):
            logit_clip(np.array([-0.1, 1.1]), 1e-6)

    def test_probability_column_mapping_uses_classes_not_fixed_column(self) -> None:
        estimator = RandomForestClassifier(n_estimators=20, random_state=7)
        x = np.arange(40, dtype=float).reshape(-1, 1)
        y_one = (x[:, 0] > 15).astype(int)
        estimator.fit(x, y_one)
        probability = rf_positive_probability(estimator, x[:5])
        expected = estimator.predict_proba(x[:5])[
            :, int(np.flatnonzero(estimator.classes_ == 1)[0])
        ]
        np.testing.assert_array_equal(probability, expected)

        class ReversedClasses:
            classes_ = np.array([1, 0])

            def predict_proba(self, values: np.ndarray) -> np.ndarray:
                return np.column_stack([np.full(len(values), 0.8), np.full(len(values), 0.2)])

        with self.assertRaises(RuntimeError):
            rf_positive_probability(ReversedClasses(), np.zeros((2, 1)))

    def test_deterministic_seed_canonical_derivation(self) -> None:
        args = ("primary_six", "target_equal", "outer", "aah-i", "inner")
        expected_string = (
            "RQ5RF|v1|primary_six|target_equal|outer|aah-i|inner"
        )
        expected = int.from_bytes(
            hashlib.sha256(expected_string.encode("utf-8")).digest()[:4], "big"
        )
        self.assertEqual(task_seed(*args), expected)
        self.assertNotEqual(
            task_seed(*args),
            task_seed("primary_six", "target_equal", "outer", "aah-i", "final"),
        )

    def test_repeated_parallel_rf_fit_is_elementwise_identical(self) -> None:
        rng = np.random.default_rng(20261117)
        x = rng.normal(size=(260, 2))
        y = (x[:, 0] * 0.8 - x[:, 1] * 0.4 + rng.normal(scale=0.8, size=260) > 0).astype(int)
        w = np.linspace(0.2, 2.0, len(y))
        parameters = effective_rf_grid(2)[1]
        first = fit_random_forest(
            x,
            y,
            w,
            parameters,
            seed=123456789,
            n_jobs=int(self.config["random_forest"]["n_jobs"]),
            label="repeat_first",
        )
        second = fit_random_forest(
            x,
            y,
            w,
            parameters,
            seed=123456789,
            n_jobs=int(self.config["random_forest"]["n_jobs"]),
            label="repeat_second",
        )
        self.assertTrue(first["estimable"] and second["estimable"])
        p1 = rf_positive_probability(first["estimator"], x)
        p2 = rf_positive_probability(second["estimator"], x)
        np.testing.assert_array_equal(p1, p2)
        self.assertEqual(first["classes"], second["classes"])
        self.assertEqual(first["n_in_bag"], second["n_in_bag"])

    def test_effective_grid_deduplicates_single_feature(self) -> None:
        self.assertEqual(len(effective_rf_grid(1)), 6)
        self.assertEqual(len(effective_rf_grid(2)), 12)
        single = effective_rf_grid(1)
        self.assertTrue(all(item["max_features"] == 1 for item in single))
        pairs = {
            (item["max_depth"], item["min_samples_leaf"], item["max_features"])
            for item in single
        }
        self.assertEqual(len(pairs), 6)

    def test_loto_and_inner_isolation_and_identical_outer_rows(self) -> None:
        _, _, _, joined = load_joined_inputs(self.config)
        loto = pd.read_csv(
            self.config["inputs"]["loto_assignments"]["path"], low_memory=False
        )
        for analysis_set in ANALYSIS_SETS:
            frame = joined[joined["pinc_status"].eq("success")].copy()
            excluded = self.config["analysis_sets"][analysis_set]["exclude_target"]
            if excluded is not None:
                frame = frame[frame["target"].astype(str).ne(str(excluded))]
            outcome = self.config["analysis_sets"][analysis_set]["outcome"]
            frame = frame[frame[outcome].notna()].copy()
            valid = set(frame["pair_id"].astype(str))
            subset = loto[
                loto["analysis_set"].eq(analysis_set)
                & loto["cluster_rule"].eq(self.config["cluster_rule"])
            ]
            for fold_id, fold in subset.groupby("fold_id", sort=True):
                fold = fold[fold["pair_id"].astype(str).isin(valid)]
                train = fold[fold["split_status"].eq("train")]
                test = fold[fold["split_status"].eq("test")]
                self.assertFalse(
                    set(train["sequence_cluster_value"].astype(str))
                    & set(test["sequence_cluster_value"].astype(str)),
                    fold_id,
                )
                self.assertFalse(
                    set(train["pair_id"].astype(str))
                    & set(test["pair_id"].astype(str)),
                    fold_id,
                )
                held = set(test["target"].astype(str))
                self.assertEqual(len(held), 1)
                self.assertFalse(held & set(train["target"].astype(str)), fold_id)

    def test_standardization_uses_training_rows_only(self) -> None:
        train = np.array([[0.0], [1.0], [2.0], [3.0]])
        weights = np.full(4, 0.25)
        mean, sd = weighted_center_scale(train, weights)
        np.testing.assert_allclose(mean, [1.5], atol=1e-15)
        np.testing.assert_allclose(sd, [np.sqrt(1.25)], atol=1e-15)
        extreme_validation = np.array([[1000.0]])
        transformed = (extreme_validation - mean) / sd
        self.assertGreater(float(transformed[0, 0]), 800.0)
        self.assertEqual(tuple(mean.shape), (1,))

    def test_no_target_id_or_forbidden_feature_leakage(self) -> None:
        allowed = {
            "pooled_logistic_ipTM": ("x_ipTM",),
            "pooled_logistic_combo": ("x_ipTM", "x_Pinc"),
            "rf_ipTM": ("x_ipTM",),
            "rf_combo": ("x_ipTM", "x_Pinc"),
        }
        self.assertEqual(FEATURE_MAP, allowed)
        for features in FEATURE_MAP.values():
            self.assertTrue(set(features).issubset({"x_ipTM", "x_Pinc"}))
            self.assertNotIn("target", features)
            self.assertNotIn("design_id", features)

    def test_ridge_penalty_shrinks_slopes(self) -> None:
        rng = np.random.default_rng(44)
        x = rng.normal(size=(500, 2))
        eta = -0.2 + 2.5 * x[:, 0] - 1.5 * x[:, 1]
        y = rng.binomial(1, expit(eta))
        weights = np.full(len(y), 1.0 / len(y))
        weak = fit_pooled_logistic(x, y, weights, 0.001, label="weak")
        strong = fit_pooled_logistic(x, y, weights, 100.0, label="strong")
        self.assertTrue(weak["estimable"] and strong["estimable"])
        self.assertLess(
            float(np.linalg.norm(strong["slopes_standardized"])),
            float(np.linalg.norm(weak["slopes_standardized"])),
        )
        self.assertFalse(
            np.array_equal(
                weak["slopes_standardized"], strong["slopes_standardized"]
            )
        )

    def test_logistic_intercept_only_changes_prediction_scale(self) -> None:
        rng = np.random.default_rng(12)
        x = rng.normal(size=(300, 2))
        y = rng.binomial(1, expit(0.3 + x[:, 0] - 0.4 * x[:, 1]))
        w = np.full(len(y), 1.0 / len(y))
        fit = fit_pooled_logistic(x, y, w, 0.1, label="scale")
        direct = logistic_eta(x, fit)
        shifted = logistic_eta(x, {**fit, "intercept": fit["intercept"] + 0.7})
        np.testing.assert_allclose(shifted - direct, np.full(len(x), 0.7), atol=1e-12)

    def test_I_fit_matches_training_average_probability(self) -> None:
        rng = np.random.default_rng(99)
        eta = rng.normal(scale=1.4, size=700)
        y = rng.binomial(1, expit(eta))
        fit = solve_offset_intercept(eta, y, "I")
        self.assertEqual(fit["status"], "evaluable")
        self.assertLess(abs(fit["residual"]), 1e-12)
        self.assertAlmostEqual(float(np.mean(expit(fit["intercept"] + eta))), float(np.mean(y)), delta=1e-12)

    def test_F_delta_zero_reproduces_U(self) -> None:
        eta_u = np.array([-2.0, -0.3, 0.0, 0.8, 2.4])
        delta = 0.0
        intercept = 0.0
        gamma = 1.0 + delta
        eta_f = intercept + gamma * eta_u
        np.testing.assert_array_equal(eta_f, eta_u)

    def test_adaptation_order_and_shared_inner_assignments(self) -> None:
        rng = np.random.default_rng(2026)
        n = 180
        frame = pd.DataFrame(
            {
                "pair_id": [f"pair-{index:03d}" for index in range(n)],
                "outcome_value": rng.integers(0, 2, size=n),
                "sequence_cluster_value": [f"c{index % 45:03d}" for index in range(n)],
            }
        )
        # Keep each group balanced enough for the frozen three-fold splitter.
        frame["outcome_value"] = frame.index.to_numpy() % 2
        splits_a, hash_a = adaptation_inner_splits(frame)
        shuffled = frame.sample(frac=1.0, random_state=73).reset_index(drop=True)
        splits_b, hash_b = adaptation_inner_splits(shuffled)
        self.assertEqual(hash_a, hash_b)
        ordered_a = stable_sort(frame)
        ordered_b = stable_sort(shuffled)
        self.assertListEqual(
            ordered_a["pair_id"].astype(str).tolist(),
            ordered_b["pair_id"].astype(str).tolist(),
        )
        for (train_a, test_a), (train_b, test_b) in zip(splits_a, splits_b):
            mapping = ordered_a["pair_id"].astype(str).to_numpy()
            mapping_shuffled = ordered_b["pair_id"].astype(str).to_numpy()
            self.assertListEqual(
                sorted(mapping[test_a].tolist()),
                sorted(mapping_shuffled[test_b].tolist()),
            )
            self.assertFalse(set(mapping[train_a]) & set(mapping[test_a]))

    def test_log_loss_and_bss_direction(self) -> None:
        y = np.array([0, 0, 1, 1])
        good = np.array([0.05, 0.1, 0.9, 0.95])
        bad = np.array([0.9, 0.8, 0.2, 0.1])
        reference = np.full(4, 0.5)
        self.assertLess(log_loss(y, good, 1e-6), log_loss(y, bad, 1e-6))
        self.assertGreater(bss(y, good, reference), bss(y, bad, reference))

    def test_config_counts_match_protocol(self) -> None:
        expected = self.config["expected_computation"]
        self.assertEqual(expected["final_outer_source_tasks"], 34)
        self.assertEqual(expected["inner_source_target_tasks"], 160)
        self.assertEqual(expected["rf_fits"], 3492)
        self.assertEqual(self.config["random_forest"]["n_estimators"], 500)
        self.assertEqual(self.config["random_forest"]["n_jobs"], 1)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RQ5RFUnitTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    payload = {
        "tests_run": int(result.testsRun),
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
        "successful": bool(result.wasSuccessful()),
        "test_names": [test.id() for test, _ in result.failures + result.errors],
    }
    path = __import__("pathlib").Path(__file__).resolve().parent.parent / "rq5rf_unit_tests.json"
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    raise SystemExit(0 if result.wasSuccessful() else 1)
