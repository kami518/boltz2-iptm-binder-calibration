"""Adaptive Gauss-Hermite integration for the RQ2 hierarchical logistic model."""
from __future__ import annotations

import numpy as np
from scipy.special import expit, logsumexp
from scipy.optimize import minimize

from rq12_core import HierarchicalLogistic


class AdaptiveHierarchicalLogistic(HierarchicalLogistic):
    """Use mode- and curvature-centred quadrature separately within each target."""

    @staticmethod
    def _mode_2d(x, y, mu_a, mu_b, log_tau_a, log_tau_b):
        mu = np.array([mu_a, mu_b], dtype=float)
        theta = mu.copy()
        precision = np.exp(-2.0 * np.array([log_tau_a, log_tau_b]))

        def log_posterior(value):
            eta = value[0] + value[1] * x
            log_likelihood = float((eta * y).sum() - np.logaddexp(0.0, eta).sum())
            return log_likelihood - 0.5 * float(np.sum(precision * (value - mu) ** 2))

        converged = False
        information = None
        for _ in range(100):
            eta = theta[0] + theta[1] * x
            probability = expit(eta)
            residual = y - probability
            gradient = np.array([residual.sum(), (residual * x).sum()]) - precision * (theta - mu)
            weights = probability * (1.0 - probability)
            information = np.array([
                [weights.sum() + precision[0], (weights * x).sum()],
                [(weights * x).sum(), (weights * x * x).sum() + precision[1]],
            ])
            step = np.linalg.solve(information, gradient)
            # Judge convergence by parameter displacement, not a fixed residual
            # sum: at n~800 the latter can plateau near 1e-6 in floating point.
            if np.linalg.norm(step) < 1e-7 * (1.0 + np.linalg.norm(theta)):
                converged = True
                break
            current = log_posterior(theta)
            scale = 1.0
            improved = False
            for _ in range(60):
                candidate = theta + scale * step
                if log_posterior(candidate) >= current:
                    theta = candidate
                    improved = True
                    break
                scale *= 0.5
            if not improved:
                if np.linalg.norm(step) < 1e-7 * (1.0 + np.linalg.norm(theta)):
                    converged = True
                    break
                raise RuntimeError("Adaptive-quadrature mode line search failed")
        if not converged:
            gradient_norm = np.linalg.norm(gradient)
            raise RuntimeError(
                f"Adaptive-quadrature mode did not converge; |gradient|={gradient_norm}, "
                f"|step|={np.linalg.norm(step)}, theta={theta.tolist()}"
            )
        eta = theta[0] + theta[1] * x
        probability = expit(eta)
        weights = probability * (1.0 - probability)
        information = np.array([
            [weights.sum() + precision[0], (weights * x).sum()],
            [(weights * x).sum(), (weights * x * x).sum() + precision[1]],
        ])
        return theta, information

    def _group_log_terms(self, x, y, mu_a, mu_b, log_tau_a, log_tau_b):
        mode, information = self._mode_2d(x, y, mu_a, mu_b, log_tau_a, log_tau_b)
        sign, logdet_information = np.linalg.slogdet(information)
        if sign <= 0:
            raise RuntimeError("Adaptive-quadrature curvature is not positive definite")
        chol = np.linalg.cholesky(np.linalg.inv(information))
        z1 = self.gh_x[:, None]
        z2 = self.gh_x[None, :]
        alpha = mode[0] + np.sqrt(2.0) * chol[0, 0] * z1
        beta = mode[1] + np.sqrt(2.0) * (chol[1, 0] * z1 + chol[1, 1] * z2)
        eta = alpha[:, :, None] + beta[:, :, None] * x[None, None, :]
        log_likelihood = (eta * y).sum(axis=2) - np.logaddexp(0.0, eta).sum(axis=2)
        precision_a, precision_b = np.exp(-2.0 * log_tau_a), np.exp(-2.0 * log_tau_b)
        log_prior = -0.5 * (precision_a * (alpha - mu_a) ** 2 + precision_b * (beta - mu_b) ** 2)
        # Two-dimensional Jacobian 2*det(chol) times the joint normal density,
        # relative to normalized Gauss-Hermite weights w1*w2/pi.
        adjustment = -0.5 * logdet_information - log_tau_a - log_tau_b
        terms = self.log_w2 + log_likelihood + log_prior + z1**2 + z2**2 + adjustment
        return terms, alpha, beta

    def _fit_random_intercept(self, initial_fixed):
        self.model_type = "random_intercept_common_slope"

        def mode_1d(x, y, mu_a, beta_common, log_tau_a):
            precision = np.exp(-2.0 * log_tau_a)
            theta = float(mu_a)
            converged = False
            for _ in range(100):
                eta = theta + beta_common * x
                probability = expit(eta)
                gradient = float((y - probability).sum() - precision * (theta - mu_a))
                weights = probability * (1.0 - probability)
                information = float(weights.sum() + precision)
                step = gradient / information
                if abs(step) < 1e-10 * (1.0 + abs(theta)):
                    converged = True
                    break
                theta += step
            if not converged:
                raise RuntimeError("Adaptive random-intercept mode did not converge")
            probability = expit(theta + beta_common * x)
            weights = probability * (1.0 - probability)
            information = float(weights.sum() + precision)
            return theta, information

        def group_terms(x, y, mu_a, beta_common, log_tau_a):
            mode, information = mode_1d(x, y, mu_a, beta_common, log_tau_a)
            alpha = mode + np.sqrt(2.0 / information) * self.gh_x
            eta = alpha[:, None] + beta_common * x[None, :]
            log_likelihood = (eta * y).sum(axis=1) - np.logaddexp(0.0, eta).sum(axis=1)
            log_prior = -0.5 * np.exp(-2.0 * log_tau_a) * (alpha - mu_a) ** 2
            # One-dimensional Jacobian and normal density, relative to w/sqrt(pi).
            adjustment = -0.5 * np.log(information) - log_tau_a
            terms = np.log(self.gh_w / np.sqrt(np.pi)) + log_likelihood + log_prior + self.gh_x**2 + adjustment
            return terms, alpha

        def objective(parameters):
            mu_a, beta_common, log_tau_a = parameters
            value = self._hyperprior(mu_a, beta_common, [log_tau_a])
            for x, y in self.groups:
                terms, _ = group_terms(x, y, mu_a, beta_common, log_tau_a)
                value += logsumexp(terms)
            return float(-value)

        start = np.array([initial_fixed[0], initial_fixed[1], 0.0])
        bounds = [(-30.0, 30.0), (-30.0, 30.0), (np.log(self.tau_low), np.log(self.tau_high))]
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6, "maxls": 50})
        self.fit_result = result
        covariance, positive_definite, eigenvalues = self.observed_information(objective, result.x)
        targets = []
        mu_a, beta_common, log_tau_a = result.x
        for x, y in self.groups:
            terms, alpha = group_terms(x, y, mu_a, beta_common, log_tau_a)
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
