"""Numerical core for the frozen Proteinbase RQ1/RQ2 analysis."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit, logsumexp
from scipy.stats import gaussian_kde
from sklearn.metrics import average_precision_score, roc_auc_score


def clipped_logit(scores, clip):
    return logit(np.clip(np.asarray(scores, dtype=float), clip, 1.0 - clip))


def fit_logistic(scores, outcomes, clip=1e-6):
    x = clipped_logit(scores, clip)
    y = np.asarray(outcomes, dtype=float)
    X = np.column_stack([np.ones(len(x)), x])

    def objective(beta):
        eta = X @ beta
        nll = np.logaddexp(0.0, eta).sum() - y @ eta
        gradient = X.T @ (expit(eta) - y)
        return nll, gradient

    start = np.array([logit(np.clip(y.mean(), 1e-4, 1 - 1e-4)), 0.0])
    result = minimize(
        objective,
        start,
        jac=True,
        method="L-BFGS-B",
        bounds=[(-40.0, 40.0), (-40.0, 40.0)],
        options={"maxiter": 300, "ftol": 1e-12, "gtol": 1e-7, "maxls": 50},
    )
    eta = X @ result.x
    weights = expit(eta) * (1.0 - expit(eta))
    information = X.T @ (weights[:, None] * X)
    eigenvalues = np.linalg.eigvalsh(information)
    covariance = np.linalg.pinv(information)
    return {
        "intercept": float(result.x[0]),
        "slope": float(result.x[1]),
        "intercept_se": float(np.sqrt(max(covariance[0, 0], 0.0))),
        "slope_se": float(np.sqrt(max(covariance[1, 1], 0.0))),
        "converged": bool(result.success),
        "message": str(result.message),
        "information_min_eigenvalue": float(eigenvalues.min()),
        "x": x,
    }


def raw_metrics(scores, outcomes, edges, clip=1e-6):
    s = np.asarray(scores, dtype=float)
    y = np.asarray(outcomes, dtype=int)
    if len(np.unique(y)) < 2:
        return None
    fit = fit_logistic(s, y, clip)
    if not fit["converged"]:
        return None
    probability = np.clip(s, clip, 1.0 - clip)
    brier = float(np.mean((s - y) ** 2))
    prevalence = float(y.mean())
    baseline = float(np.mean((prevalence - y) ** 2))
    bin_index = np.digitize(s, edges[1:-1], right=False)
    ece = 0.0
    for b in range(len(edges) - 1):
        mask = bin_index == b
        if mask.any():
            ece += mask.mean() * abs(y[mask].mean() - s[mask].mean())
    return {
        "roc_auc": float(roc_auc_score(y, s)),
        "average_precision": float(average_precision_score(y, s)),
        "positive_rate": prevalence,
        "brier": brier,
        "baseline_brier": baseline,
        "bss_target_prevalence": float(1.0 - brier / baseline) if baseline > 0 else np.nan,
        "log_loss": float(-np.mean(y * np.log(probability) + (1 - y) * np.log(1 - probability))),
        "ece": float(ece),
        "calibration_intercept": fit["intercept"],
        "calibration_slope": fit["slope"],
        "calibration_intercept_se": fit["intercept_se"],
        "calibration_slope_se": fit["slope_se"],
    }


def reliability_rows(scores, outcomes, edges):
    s = np.asarray(scores, dtype=float)
    y = np.asarray(outcomes, dtype=int)
    indices = np.digitize(s, edges[1:-1], right=False)
    rows = []
    z = 1.959963984540054
    for b in range(len(edges) - 1):
        mask = indices == b
        n = int(mask.sum())
        positives = int(y[mask].sum())
        observed = positives / n if n else np.nan
        if n:
            denominator = 1.0 + z * z / n
            center = (observed + z * z / (2 * n)) / denominator
            half = z * np.sqrt(observed * (1 - observed) / n + z * z / (4 * n * n)) / denominator
            low, high = max(0.0, center - half), min(1.0, center + half)
        else:
            low = high = np.nan
        rows.append({
            "bin_index": b,
            "bin_left": edges[b],
            "bin_right": edges[b + 1],
            "n": n,
            "positives": positives,
            "mean_iptm": float(s[mask].mean()) if n else np.nan,
            "observed_rate": float(observed) if n else np.nan,
            "wilson_low": float(low) if n else np.nan,
            "wilson_high": float(high) if n else np.nan,
        })
    return rows


def quantile_mcse(values, quantile, estimate):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 20 or np.unique(values).size < 3:
        return np.nan
    try:
        density = float(gaussian_kde(values)([estimate])[0])
        if not np.isfinite(density) or density <= 1e-12:
            return np.nan
        return float(np.sqrt(quantile * (1 - quantile) / len(values)) / density)
    except Exception:
        return np.nan


def bootstrap_summary(values):
    values = np.asarray(values, dtype=float)
    valid = values[np.isfinite(values)]
    if not len(valid):
        return {"valid_bootstrap": 0}
    low, high = np.quantile(valid, [0.025, 0.975])
    mean = float(valid.mean())
    sd = float(valid.std(ddof=1)) if len(valid) > 1 else 0.0
    return {
        "bootstrap_mean": mean,
        "bootstrap_sd": sd,
        "ci_low": float(low),
        "ci_high": float(high),
        "bootstrap_mean_mcse": float(sd / np.sqrt(len(valid))),
        "ci_low_mcse": quantile_mcse(valid, 0.025, low),
        "ci_high_mcse": quantile_mcse(valid, 0.975, high),
        "valid_bootstrap": int(len(valid)),
    }


class HierarchicalLogistic:
    """Independent random target intercepts/slopes integrated by Gauss-Hermite quadrature."""

    def __init__(self, groups, config, nodes=None):
        self.groups = [(np.asarray(x, dtype=float), np.asarray(y, dtype=float)) for x, y in groups]
        self.config = config
        self.nodes = int(nodes or config["quadrature_nodes"])
        self.tau_low, self.tau_high = config["tau_bounds"]
        self.fit_result = None
        self.model_type = "random_intercept_and_slope"
        self.fallback_reason = None
        self._set_nodes(self.nodes)

    def _set_nodes(self, nodes):
        self.nodes = int(nodes)
        self.gh_x, self.gh_w = np.polynomial.hermite.hermgauss(self.nodes)
        self.log_w2 = np.log(self.gh_w[:, None] * self.gh_w[None, :] / np.pi)

    def _hyperprior(self, mu_a, mu_b, log_taus):
        prior = self.config
        value = -0.5 * (mu_a / prior["mu_intercept_prior"]["sd"]) ** 2
        value += -0.5 * (mu_b / prior["mu_slope_prior"]["sd"]) ** 2
        value += -0.5 * np.sum((np.asarray(log_taus) - prior["log_tau_prior"]["mean_log_tau"]) / prior["log_tau_prior"]["sd"]) ** 2
        return float(value)

    def _group_log_terms(self, x, y, mu_a, mu_b, log_tau_a, log_tau_b):
        tau_a, tau_b = np.exp(log_tau_a), np.exp(log_tau_b)
        alpha = mu_a + np.sqrt(2.0) * tau_a * self.gh_x[:, None]
        beta = mu_b + np.sqrt(2.0) * tau_b * self.gh_x[None, :]
        eta = alpha[:, :, None] + beta[:, :, None] * x[None, None, :]
        log_likelihood = (eta * y).sum(axis=2) - np.logaddexp(0.0, eta).sum(axis=2)
        return self.log_w2 + log_likelihood, alpha, beta

    def objective(self, parameters):
        mu_a, mu_b, log_tau_a, log_tau_b = parameters
        value = self._hyperprior(mu_a, mu_b, [log_tau_a, log_tau_b])
        for x, y in self.groups:
            terms, _, _ = self._group_log_terms(x, y, mu_a, mu_b, log_tau_a, log_tau_b)
            value += logsumexp(terms)
        return float(-value)

    def _initial_parameters(self):
        x = np.concatenate([g[0] for g in self.groups])
        y = np.concatenate([g[1] for g in self.groups])
        X = np.column_stack([np.ones(len(x)), x])

        def objective(beta):
            eta = X @ beta
            return np.logaddexp(0, eta).sum() - y @ eta, X.T @ (expit(eta) - y)

        result = minimize(objective, [0.0, 0.0], jac=True, method="BFGS")
        return np.array([result.x[0], result.x[1], 0.0, 0.0])

    @staticmethod
    def observed_information(objective, parameters):
        parameters = np.asarray(parameters, dtype=float)
        n = len(parameters)
        hessian = np.zeros((n, n))
        steps = np.array([1e-4 * max(1.0, abs(v)) for v in parameters])
        for i in range(n):
            ei = np.zeros(n)
            ei[i] = steps[i]
            hessian[i, i] = (objective(parameters + ei) - 2 * objective(parameters) + objective(parameters - ei)) / steps[i] ** 2
            for j in range(i + 1, n):
                ej = np.zeros(n)
                ej[j] = steps[j]
                value = (
                    objective(parameters + ei + ej)
                    - objective(parameters + ei - ej)
                    - objective(parameters - ei + ej)
                    + objective(parameters - ei - ej)
                ) / (4 * steps[i] * steps[j])
                hessian[i, j] = hessian[j, i] = value
        # objective is negative log posterior, so its Hessian is observed information.
        information = hessian
        eigenvalues = np.linalg.eigvalsh(information)
        positive_definite = bool(eigenvalues.min() > 1e-7)
        covariance = np.linalg.inv(information) if positive_definite else np.linalg.pinv(information)
        return covariance, positive_definite, eigenvalues

    def fit(self, initial=None):
        initial = self._initial_parameters() if initial is None else np.asarray(initial, dtype=float)
        bounds = [(-30.0, 30.0), (-30.0, 30.0), (np.log(self.tau_low), np.log(self.tau_high)), (np.log(self.tau_low), np.log(self.tau_high))]
        starts = [initial]
        if initial[2] != -0.5:
            starts.append(np.array([initial[0], initial[1], -0.5, -0.5]))
        results = [minimize(self.objective, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6, "maxls": 50}) for start in starts]
        result = min(results, key=lambda r: r.fun if np.isfinite(r.fun) else np.inf)
        self.fit_result = result
        tau_b = float(np.exp(result.x[3]))
        covariance, positive_definite, eigenvalues = self.observed_information(self.objective, result.x)
        identifiable = (
            result.success
            and positive_definite
            and self.tau_low * 1.001 < tau_b < self.tau_high * 0.999
            and tau_b >= self.config["random_slope_collapse_threshold"]
        )
        if not identifiable:
            self.fallback_reason = (
                "random-slope convergence/identifiability gate: "
                f"success={result.success}, information_pd={positive_definite}, "
                f"tau_beta={tau_b:.6g}, message={result.message}"
            )
            return self._fit_random_intercept(initial[:2])
        targets = self._posterior_targets(result.x)
        return self._format_result(result, covariance, positive_definite, eigenvalues, targets)

    def _posterior_targets(self, parameters):
        mu_a, mu_b, log_tau_a, log_tau_b = parameters
        targets = []
        for x, y in self.groups:
            terms, alpha, beta = self._group_log_terms(x, y, mu_a, mu_b, log_tau_a, log_tau_b)
            posterior = np.exp(terms - logsumexp(terms))
            mean_a = float((posterior * alpha).sum())
            mean_b = float((posterior * beta).sum())
            var_a = float((posterior * (alpha - mean_a) ** 2).sum())
            var_b = float((posterior * (beta - mean_b) ** 2).sum())
            covariance = float((posterior * (alpha - mean_a) * (beta - mean_b)).sum())
            targets.append({
                "alpha_partial": mean_a,
                "beta_partial": mean_b,
                "alpha_posterior_sd": np.sqrt(max(var_a, 0.0)),
                "beta_posterior_sd": np.sqrt(max(var_b, 0.0)),
                "alpha_beta_posterior_correlation": covariance / np.sqrt(max(var_a * var_b, 1e-30)),
            })
        return targets

    def _format_result(self, result, covariance, positive_definite, eigenvalues, targets):
        se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        return {
            "model_type": self.model_type,
            "fallback_reason": self.fallback_reason,
            "converged": bool(result.success),
            "optimizer_message": str(result.message),
            "objective": float(result.fun),
            "quadrature_nodes": self.nodes,
            "information_positive_definite": positive_definite,
            "information_min_eigenvalue": float(eigenvalues.min()),
            "hyperparameters": {
                "mu_intercept": float(result.x[0]),
                "mu_slope": float(result.x[1]),
                "tau_intercept": float(np.exp(result.x[2])),
                "tau_slope": float(np.exp(result.x[3])) if self.model_type == "random_intercept_and_slope" else None,
                "mu_intercept_se": float(se[0]),
                "mu_slope_se": float(se[1]),
                "log_tau_intercept_se": float(se[2]),
                "log_tau_slope_se": float(se[-1]) if self.model_type == "random_intercept_and_slope" else None,
            },
            "targets": targets,
            "parameters": [float(v) for v in result.x],
        }

    def _fit_random_intercept(self, initial_fixed):
        self.model_type = "random_intercept_common_slope"

        def objective(parameters):
            mu_a, beta_common, log_tau_a = parameters
            value = self._hyperprior(mu_a, beta_common, [log_tau_a])
            alpha = mu_a + np.sqrt(2.0) * np.exp(log_tau_a) * self.gh_x
            log_weights = np.log(self.gh_w / np.sqrt(np.pi))
            for x, y in self.groups:
                eta = alpha[:, None] + beta_common * x[None, :]
                log_likelihood = (eta * y).sum(axis=1) - np.logaddexp(0.0, eta).sum(axis=1)
                value += logsumexp(log_weights + log_likelihood)
            return float(-value)

        start = np.array([initial_fixed[0], initial_fixed[1], 0.0])
        bounds = [(-30.0, 30.0), (-30.0, 30.0), (np.log(self.tau_low), np.log(self.tau_high))]
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6, "maxls": 50})
        self.fit_result = result
        covariance, positive_definite, eigenvalues = self.observed_information(objective, result.x)
        targets = []
        mu_a, beta_common, log_tau_a = result.x
        alpha = mu_a + np.sqrt(2.0) * np.exp(log_tau_a) * self.gh_x
        log_weights = np.log(self.gh_w / np.sqrt(np.pi))
        for x, y in self.groups:
            eta = alpha[:, None] + beta_common * x[None, :]
            terms = log_weights + (eta * y).sum(axis=1) - np.logaddexp(0.0, eta).sum(axis=1)
            posterior = np.exp(terms - logsumexp(terms))
            mean_a = float((posterior * alpha).sum())
            var_a = float((posterior * (alpha - mean_a) ** 2).sum())
            targets.append({
                "alpha_partial": mean_a,
                "beta_partial": float(beta_common),
                "alpha_posterior_sd": np.sqrt(max(var_a, 0.0)),
                "beta_posterior_sd": np.nan,
                "alpha_beta_posterior_correlation": np.nan,
            })
        return self._format_result(result, covariance, positive_definite, eigenvalues, targets)


def finite_hessian_check(result):
    values = []
    for key, value in result["hyperparameters"].items():
        if value is not None:
            values.append(np.isfinite(value))
    return bool(all(values) and result["converged"] and result["information_positive_definite"])
