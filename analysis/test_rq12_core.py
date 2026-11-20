"""Small deterministic checks for the RQ1/RQ2 numerical core."""
from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
sys.path.insert(0, str(ROOT / "analysis"))
sys.path.insert(0, str(ROOT / "analysis" / "rq12_lib"))

import parasail
from rq12_core import HierarchicalLogistic, fit_logistic, raw_metrics, reliability_rows
from rq12_adaptive import AdaptiveHierarchicalLogistic

spec = importlib.util.spec_from_file_location(
    "rq12_clustering",
    ROOT / "analysis" / "build_sequence_clusters_rq12_2026-11-09.py",
)
clustering = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clustering)
CONFIG = json.loads((ROOT / "analysis" / "rq12_analysis_config.json").read_text(encoding="utf-8"))["hierarchical_model"]


class RQ12CoreTests(unittest.TestCase):
    def test_logistic_recovers_signal(self):
        rng = np.random.default_rng(17)
        scores = rng.uniform(0.05, 0.95, 1000)
        probability = 1 / (1 + np.exp(-(-2 + 4 * scores)))
        outcomes = rng.binomial(1, probability)
        fit = fit_logistic(scores, outcomes)
        self.assertTrue(fit["converged"])
        self.assertTrue(np.isfinite(fit["intercept_se"]))
        self.assertGreater(fit["slope"], 0)

    def test_raw_brier_and_reliability_counts(self):
        scores = np.array([0.1, 0.2, 0.85, 0.95])
        outcomes = np.array([0, 1, 0, 1])
        edges = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
        result = raw_metrics(scores, outcomes, edges)
        self.assertAlmostEqual(result["brier"], np.mean((scores - outcomes) ** 2))
        bins = reliability_rows(scores, outcomes, edges)
        self.assertEqual(sum(row["n"] for row in bins), len(scores))
        self.assertEqual(sum(row["positives"] for row in bins), int(outcomes.sum()))

    def test_alignment_cigar_identity_and_coverage(self):
        result = parasail.sw_trace_striped_16("ACDEFG", "ACDXXG", 11, 1, parasail.blosum62)
        identity, query_coverage, target_coverage = clustering.alignment_values(result, 6, 6)
        self.assertAlmostEqual(identity, 4 / 6)
        self.assertAlmostEqual(query_coverage, 1.0)
        self.assertAlmostEqual(target_coverage, 1.0)

    def test_quadrature_weights_and_degenerate_marginal(self):
        groups = [(np.array([-1.0, 0.0, 1.0]), np.array([0.0, 0.0, 1.0]))]
        model = AdaptiveHierarchicalLogistic(groups, CONFIG, nodes=11)
        self.assertAlmostEqual(float(np.exp(model.log_w2).sum()), 1.0)
        parameters = [0.2, 0.8, np.log(1e-6), np.log(1e-6)]
        terms, _, _ = model._group_log_terms(*groups[0], *parameters)
        eta = parameters[0] + parameters[1] * groups[0][0]
        expected = float((eta * groups[0][1]).sum() - np.logaddexp(0, eta).sum())
        self.assertAlmostEqual(float(logsumexp(terms)), expected, places=6)

    def test_adaptive_quadrature_against_independent_dense_grid(self):
        rng = np.random.default_rng(43)
        x = rng.normal(0.0, 0.6, 36)
        y = rng.binomial(1, 1 / (1 + np.exp(-(-0.4 + 1.1 * x))))
        parameters = [0.15, -0.2, np.log(0.7), np.log(0.9)]
        model = AdaptiveHierarchicalLogistic([(x, y)], CONFIG, nodes=21)
        terms, _, _ = model._group_log_terms(x, y, *parameters)
        adaptive = float(logsumexp(terms))

        # Independent brute-force integral of likelihood times both normal densities.
        alpha = np.linspace(-5.0, 5.0, 401)
        beta = np.linspace(-5.0, 5.0, 401)
        eta = alpha[:, None, None] + beta[None, :, None] * x[None, None, :]
        log_terms = (eta * y).sum(axis=2) - np.logaddexp(0.0, eta).sum(axis=2)
        mu_a, mu_b, log_tau_a, log_tau_b = parameters
        log_terms += -0.5 * ((alpha[:, None] - mu_a) / np.exp(log_tau_a)) ** 2 - log_tau_a
        log_terms += -0.5 * ((beta[None, :] - mu_b) / np.exp(log_tau_b)) ** 2 - log_tau_b
        log_terms -= np.log(2.0 * np.pi)
        peak = float(log_terms.max())
        weights = np.exp(log_terms - peak)
        integral = np.trapezoid(np.trapezoid(weights, beta, axis=1), alpha)
        direct = peak + float(np.log(integral))
        self.assertAlmostEqual(adaptive, direct, delta=2e-5)

    def test_hierarchical_fit_on_simulated_groups(self):
        rng = np.random.default_rng(29)
        groups = []
        for alpha, beta in [(-1.0, 0.8), (0.0, 1.2), (0.8, 0.6)]:
            x = rng.normal(size=160)
            outcomes = rng.binomial(1, 1 / (1 + np.exp(-(alpha + beta * x))))
            groups.append((x, outcomes))
        model = AdaptiveHierarchicalLogistic(groups, CONFIG, nodes=11)
        result = model.fit()
        self.assertTrue(result["converged"])
        self.assertEqual(result["model_type"], "random_intercept_and_slope", result["fallback_reason"])
        self.assertTrue(result["information_positive_definite"])
        self.assertEqual(len(result["targets"]), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
