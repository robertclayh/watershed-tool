"""Numerical helpers shared by the estimators."""

from __future__ import annotations

import numpy as np
from scipy.optimize import nnls


def project_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection onto the probability simplex (Duchi et al., 2008)."""
    n = v.size
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    ind = np.arange(1, n + 1)
    cond = u - css / ind > 0
    rho = ind[cond][-1]
    theta = css[cond][-1] / rho
    return np.maximum(v - theta, 0.0)


def _simplex_fista(A, b, ridge, max_iter, tol, w0=None):
    AtA = A.T @ A
    Atb = A.T @ b
    J = A.shape[1]
    L = 2.0 * (np.linalg.eigvalsh(AtA)[-1] + ridge)
    if L <= 0:
        return np.full(J, 1.0 / J)
    w = np.full(J, 1.0 / J) if w0 is None else project_simplex(np.asarray(w0, float))
    z = w.copy()
    t = 1.0
    for _ in range(max_iter):
        grad = 2.0 * (AtA @ z - Atb) + 2.0 * ridge * z
        w_new = project_simplex(z - grad / L)
        t_new = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        z = w_new + ((t - 1.0) / t_new) * (w_new - w)
        if np.max(np.abs(w_new - w)) < tol:
            return w_new
        w, t = w_new, t_new
    return w


def simplex_lstsq(
    A: np.ndarray,
    b: np.ndarray,
    ridge: float = 0.0,
    max_iter: int = 20000,
    tol: float = 1e-10,
) -> np.ndarray:
    """Solve ``min ||A w - b||^2 + ridge ||w||^2`` subject to ``w >= 0, sum(w) = 1``.

    The sum-to-one constraint is enforced through a heavily weighted extra row,
    and the result is solved exactly by Lawson-Hanson NNLS. Synthetic-control
    problems are typically rank-deficient (more donors than pre-periods), where
    first-order methods converge slowly. FISTA (accelerated projected gradient)
    is the fallback if NNLS fails to converge.
    """
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    J = A.shape[1]
    if J == 1:
        return np.ones(1)
    rows, rhs = [A], [b]
    if ridge > 0:
        rows.append(np.sqrt(ridge) * np.eye(J))
        rhs.append(np.zeros(J))
    M = 1e3 * max(np.abs(A).max(), np.abs(b).max(), 1.0)
    rows.append(np.full((1, J), M))
    rhs.append(np.array([M]))
    try:
        w, _ = nnls(np.vstack(rows), np.concatenate(rhs), maxiter=50 * J + 100)
        if w.sum() > 0:
            return w / w.sum()
    except RuntimeError:
        pass
    return _simplex_fista(A, b, ridge, max_iter, tol)


def twoway_demean(M: np.ndarray) -> np.ndarray:
    """Two-way (unit and time) within transformation for a balanced panel.

    ``M`` has shape ``(..., N, T)``; leading axes are treated as a batch, so this
    also demeans every bootstrap draw or every regressor at once.
    """
    return (
        M
        - M.mean(axis=-1, keepdims=True)
        - M.mean(axis=-2, keepdims=True)
        + M.mean(axis=(-2, -1), keepdims=True)
    )


def ridge_loo_lambda(X: np.ndarray, Y: np.ndarray, grid: np.ndarray | None = None) -> float:
    """Choose a ridge penalty by closed-form leave-one-out CV.

    ``X`` is (n, p) centred predictors, ``Y`` is (n, q) centred responses.
    """
    n = X.shape[0]
    if grid is None:
        s = np.linalg.svd(X, compute_uv=False)
        top = max(s[0] ** 2, 1e-12)
        grid = top * np.logspace(-4, 1, 30)
    best, best_err = grid[0], np.inf
    for lam in grid:
        H = X @ np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T)
        resid = Y - H @ Y
        denom = 1.0 - np.clip(np.diag(H), None, 1 - 1e-9)
        err = np.mean((resid / denom[:, None]) ** 2) if n > 1 else np.inf
        if err < best_err:
            best, best_err = lam, err
    return float(best)
