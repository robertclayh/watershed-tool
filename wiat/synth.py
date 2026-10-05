"""Synthetic control for a single treated basin.

The counterfactual for the treated basin is a convex combination of donor basins
chosen to reproduce its pre-intervention trajectory (Abadie, Diamond & Hainmueller,
2010). Three variants:

* ``"standard"``: weights on the simplex, fitted on pre-period outcome levels
  (plus optional standardised static covariates).
* ``"demeaned"``: fits deviations from pre-period means and adds back the level
  difference (Ferman & Pinto, 2021). Use it when the treated basin is wetter or drier
  than every donor, which is common for headwater basins.
* ``"augmented"``: ridge-augmented SCM (Ben-Michael, Feller & Rothstein, 2021),
  built on the de-meaned fit. It corrects the bias that remains when the treated
  basin's pre-period path lies outside the donors' convex hull, at the cost of
  extrapolation.

Inference:

* **Placebo-in-space** (Abadie et al., 2010): every donor is treated as if it were
  the treated basin. The p-value is the treated basin's rank among the placebos on
  ``|ATT| / pre-RMSPE``; Abadie's post/pre RMSPE ratio is also reported. Inverting
  the test gives a confidence interval (Firpo & Possebom, 2018). The smallest
  attainable p-value is ``1 / (J + 1)``. With 8 donors, no effect can ever be
  significant at 10%, and the interval is the whole real line.
* **Conformal** (Chernozhukov, Wüthrich & Zhu, 2021): a constant-effect confidence
  interval from time-permutation of residuals. Its resolution is set by the number of
  periods, not donors.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ._linalg import ridge_loo_lambda, simplex_lstsq
from .panel import TIME, UNIT


@dataclass
class SCMFit:
    weights: np.ndarray
    counterfactual: np.ndarray  # (T,)
    pre_rmspe: float
    post_rmspe: float


def _fit(
    y1: np.ndarray,
    Y0: np.ndarray,
    pre: np.ndarray,
    method: str = "standard",
    z1: np.ndarray | None = None,
    Z0: np.ndarray | None = None,
    covariate_weight: float = 1.0,
    ridge: float = 0.0,
) -> SCMFit:
    """Fit SCM weights. ``y1`` (T,), ``Y0`` (T, J), ``pre`` boolean (T,)."""
    if method not in {"standard", "demeaned", "augmented"}:
        raise ValueError(f"unknown method {method!r}")
    yp, Yp = y1[pre], Y0[pre]
    centered = method in {"demeaned", "augmented"}
    if centered:
        A, b = Yp - Yp.mean(axis=0), yp - yp.mean()
    else:
        A, b = Yp, yp
    if z1 is not None and Z0 is not None and covariate_weight > 0:
        # Rescale covariate rows so they carry `covariate_weight` times the
        # influence of an average pre-period outcome row.
        scale = np.sqrt(covariate_weight * np.mean(A**2) / max(np.mean(Z0**2), 1e-12))
        A = np.vstack([A, scale * Z0])
        b = np.concatenate([b, scale * z1])
    w = simplex_lstsq(A, b, ridge=ridge)
    cf = Y0 @ w
    if centered:
        cf = cf + (yp.mean() - Yp.mean(axis=0) @ w)
    if method == "augmented" and Y0.shape[1] > 2 and (~pre).any():
        # Ridge-augmented SCM on de-meaned series: regress each donor's post-period
        # deviations (from its own pre-mean) on its de-meaned pre-period path. Then
        # correct the counterfactual by the predicted effect of the treated basin's
        # remaining pre-period imbalance.
        mu_d = Yp.mean(axis=0)
        Xd = (Yp - mu_d).T  # (J, T0)
        Yd = (Y0[~pre] - mu_d).T  # (J, T1)
        Xc = Xd - Xd.mean(axis=0)
        Yc = Yd - Yd.mean(axis=0)
        lam = ridge_loo_lambda(Xc, Yc)
        eta = np.linalg.solve(Xc.T @ Xc + lam * np.eye(Xc.shape[1]), Xc.T @ Yc)
        imbalance = (yp - yp.mean()) - Xd.T @ w
        cf = cf.copy()
        # Pre-period counterfactual is left as the SCM fit, so pre-RMSPE stays
        # comparable across methods and reflects true donor-pool adequacy.
        cf[~pre] = cf[~pre] + imbalance @ eta
    gap = y1 - cf
    pre_rmspe = float(np.sqrt(np.mean(gap[pre] ** 2)))
    post_rmspe = float(np.sqrt(np.mean(gap[~pre] ** 2))) if (~pre).any() else np.nan
    return SCMFit(w, cf, pre_rmspe, post_rmspe)


@dataclass
class SCMResult:
    treated: object
    times: list
    donors: list
    weights: pd.Series
    y_treated: pd.Series
    y_synthetic: pd.Series
    pre_mask: np.ndarray
    att: float
    pre_rmspe: float
    post_rmspe: float
    method: str
    alpha: float
    placebo: pd.DataFrame | None = None
    placebo_gaps: pd.DataFrame | None = None
    p_value: float = np.nan  # placebo test on |ATT| / pre-RMSPE
    p_value_rmspe: float = np.nan  # Abadie et al. post/pre RMSPE-ratio test
    min_attainable_p: float = np.nan
    ci: tuple[float, float] = (-np.inf, np.inf)
    conformal: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)

    @property
    def gap(self) -> pd.Series:
        return self.y_treated - self.y_synthetic

    @property
    def significant(self) -> bool:
        return bool(self.p_value <= self.alpha)

    @property
    def effective_donors(self) -> float:
        w = self.weights.values
        return float(1.0 / np.sum(w**2)) if w.sum() > 0 else 0.0


def synthetic_control(
    df: pd.DataFrame,
    outcome: str,
    treated,
    treat_start,
    donors: Sequence | None = None,
    unit: str = UNIT,
    time: str = TIME,
    method: str = "demeaned",
    static: pd.DataFrame | None = None,
    covariates: Sequence[str] = (),
    covariate_weight: float = 0.5,
    alpha: float = 0.1,
    placebo: bool = True,
    placebo_fit_cutoff: float | None = 5.0,
    conformal: bool = False,
    conformal_grid: int = 81,
    conformal_permutation: str = "moving_block",
) -> SCMResult:
    """Estimate the effect of an intervention in one basin by synthetic control.

    Parameters
    ----------
    donors
        Candidate control basins. Defaults to every other basin in ``df``. Pass the
        output of :func:`wiat.matching.select_donors` to restrict to suitable ones.
    static, covariates
        Optional static basin attributes (indexed by ``unit``) to balance alongside
        the pre-period outcome path. Each is standardised across basins first.
    placebo_fit_cutoff
        Placebo runs whose pre-period RMSPE exceeds this multiple of the treated
        basin's are discarded, as in Abadie et al. (2010). ``None`` keeps all.
    """
    wide = df.pivot(index=time, columns=unit, values=outcome).sort_index()
    times = wide.index.tolist()
    if donors is None:
        donors = [u for u in wide.columns if u != treated]
    donors = [d for d in donors if d != treated]
    if len(donors) < 1:
        raise ValueError("synthetic control needs at least one donor basin")
    sub = wide[[treated, *donors]]
    if sub.isna().any().any():
        bad = sub.columns[sub.isna().any()].tolist()
        raise ValueError(f"missing outcome values for basins {bad[:10]}; balance the panel first")
    start = pd.Timestamp(treat_start) if isinstance(times[0], pd.Timestamp) else treat_start
    pre = np.array([t < start for t in times])
    if pre.sum() < 2 or (~pre).sum() < 1:
        raise ValueError("need at least two pre-treatment and one post-treatment period")
    y1 = sub[treated].to_numpy(float)
    Y0 = sub[donors].to_numpy(float)

    Zall = None
    if static is not None and covariates:
        Zs = static.loc[[treated, *donors], list(covariates)].astype(float)
        Zs = (Zs - Zs.mean()) / Zs.std(ddof=0).replace(0, 1)
        Zall = Zs.to_numpy()

    def fit_one(y, Y, z1=None, Z0=None):
        return _fit(y, Y, pre, method, z1, Z0, covariate_weight)

    z1 = Zall[0] if Zall is not None else None
    Z0 = Zall[1:].T if Zall is not None else None
    main = fit_one(y1, Y0, z1, Z0)
    gap = y1 - main.counterfactual
    att = float(gap[~pre].mean())

    res = SCMResult(
        treated=treated,
        times=times,
        donors=list(donors),
        weights=pd.Series(main.weights, index=donors, name="weight"),
        y_treated=pd.Series(y1, index=times, name=treated),
        y_synthetic=pd.Series(main.counterfactual, index=times, name="synthetic"),
        pre_mask=pre,
        att=att,
        pre_rmspe=main.pre_rmspe,
        post_rmspe=main.post_rmspe,
        method=method,
        alpha=alpha,
    )

    # Fit diagnostics: how well can the donor pool reproduce the treated basin at all?
    sd_pre = float(np.std(y1[pre], ddof=1)) if pre.sum() > 1 else np.nan
    # Share of pre-periods where the treated basin lies outside the donors' range.
    # No convex combination can reproduce such periods. For de-meaned methods only
    # the de-meaned range matters, but level extremity is still worth reporting.
    lo, hi = Y0[pre].min(axis=1), Y0[pre].max(axis=1)
    outside_level = float(np.mean((y1[pre] < lo) | (y1[pre] > hi)))
    yc = y1[pre] - y1[pre].mean()
    Yc = Y0[pre] - Y0[pre].mean(axis=0)
    outside_dm = float(np.mean((yc < Yc.min(axis=1)) | (yc > Yc.max(axis=1))))
    outside = outside_dm if method in {"demeaned", "augmented"} else outside_level
    res.diagnostics = {
        "n_donors": len(donors),
        "effective_donors": res.effective_donors,
        "pre_rmspe_normalized": main.pre_rmspe / sd_pre if sd_pre and sd_pre > 0 else np.nan,
        "pre_r2": 1 - main.pre_rmspe**2 / np.var(y1[pre]) if np.var(y1[pre]) > 0 else np.nan,
        "frac_pre_outside_donor_range": outside,
        "frac_pre_outside_donor_levels": outside_level,
        "max_weight": float(main.weights.max()),
    }

    if placebo and len(donors) >= 2:
        rows, gaps = [], {}
        for j, d in enumerate(donors):
            others = [k for k in range(len(donors)) if k != j]
            zj = Zall[1 + j] if Zall is not None else None
            Zo = Zall[1:][others].T if Zall is not None else None
            f = fit_one(Y0[:, j], Y0[:, others], zj, Zo)
            g = Y0[:, j] - f.counterfactual
            gaps[d] = g
            rows.append(
                {
                    unit: d,
                    "pre_rmspe": f.pre_rmspe,
                    "post_rmspe": f.post_rmspe,
                    "ratio": f.post_rmspe / max(f.pre_rmspe, 1e-12),
                    "att": float(g[~pre].mean()),
                }
            )
        pl = pd.DataFrame(rows).set_index(unit)
        pl["kept"] = True
        if placebo_fit_cutoff is not None:
            pl["kept"] = pl.pre_rmspe <= placebo_fit_cutoff * max(main.pre_rmspe, 1e-12)
        kept = pl[pl.kept]
        ratio = main.post_rmspe / max(main.pre_rmspe, 1e-12)
        pl["att_stat"] = pl.att.abs() / pl.pre_rmspe.clip(lower=1e-12)
        kept = pl[pl.kept]
        stat = abs(att) / max(main.pre_rmspe, 1e-12)
        res.p_value = float((1 + np.sum(kept.att_stat >= stat)) / (1 + len(kept)))
        res.p_value_rmspe = float((1 + np.sum(kept.ratio >= ratio)) / (1 + len(kept)))
        res.min_attainable_p = 1.0 / (1 + len(kept))
        res.diagnostics["placebos_discarded"] = int((~pl.kept).sum())
        res.ci = placebo_inverted_interval(att, main.pre_rmspe, kept.att_stat.to_numpy(), alpha)
        res.placebo = pl
        res.placebo_gaps = pd.DataFrame(gaps, index=times)
    else:
        res.min_attainable_p = 1.0

    if conformal:
        res.conformal = conformal_interval(
            y1, Y0, pre, method=method if method != "augmented" else "standard",
            alpha=alpha, n_grid=conformal_grid, center=att, permutation=conformal_permutation,
            scale=max(main.pre_rmspe, np.std(gap[~pre]) if (~pre).sum() > 1 else 0, 1e-9),
        )
    return res


def placebo_inverted_interval(
    att: float, pre_rmspe: float, placebo_stats: np.ndarray, alpha: float
) -> tuple[float, float]:
    """Confidence interval for the ATT by inverting the placebo-in-space test.

    Follows Firpo & Possebom (2018) with the standardised statistic
    ``|mean post gap| / pre-RMSPE``. The null ``ATT = theta`` is accepted when the
    treated basin's statistic, after subtracting ``theta`` from its post-period gap,
    is not among the ``alpha`` share of most extreme placebo statistics. Placebo fits
    do not depend on ``theta``, so the interval is ``att ± c * pre_rmspe``, where ``c``
    is the corresponding order statistic of the placebo distribution.

    The interval is ``(-inf, inf)`` when ``alpha * (J + 1) <= 1``: with fewer than
    ``1/alpha - 1`` usable donors, the data cannot rule out *any* effect size.
    """
    J = len(placebo_stats)
    m = int(np.floor(alpha * (J + 1) - 1 + 1e-9)) + 1
    if m <= 0:
        return (-np.inf, np.inf)
    c = float(np.sort(np.asarray(placebo_stats))[::-1][m - 1])
    return (att - c * pre_rmspe, att + c * pre_rmspe)


def _conformal_p(y1, Y0, pre, theta, method, permutation="moving_block", n_perm=999, seed=0):
    y = y1.copy()
    y[~pre] = y[~pre] - theta
    all_t = np.ones_like(pre, dtype=bool)
    f = _fit(y, Y0, all_t, method)
    u = y - f.counterfactual
    T = u.size
    post_idx = np.where(~pre)[0]
    T1 = post_idx.size
    stat = np.sum(np.abs(u[post_idx])) / np.sqrt(T1)
    if permutation == "moving_block":
        # Cyclic shifts preserve serial dependence; resolution is limited to 1/T.
        shifted = np.array([np.roll(u, -s)[post_idx] for s in range(T)])
    else:
        # i.i.d. permutations: finer resolution, but assume exchangeable residuals.
        rng = np.random.default_rng(seed)
        idx = np.argsort(rng.random((n_perm, T)), axis=1)[:, :T1]
        shifted = np.vstack([u[post_idx][None], u[idx]])
    perm = np.sum(np.abs(shifted), axis=1) / np.sqrt(T1)
    return float(np.mean(perm >= stat - 1e-12))


def conformal_interval(
    y1: np.ndarray,
    Y0: np.ndarray,
    pre: np.ndarray,
    method: str = "demeaned",
    alpha: float = 0.1,
    n_grid: int = 81,
    center: float = 0.0,
    scale: float = 1.0,
    width: float = 6.0,
    permutation: str = "moving_block",
) -> dict:
    """Conformal confidence set for a constant post-period effect (CWZ, 2021).

    Post-period outcomes require the treated basin to be untreated after the null
    effect is removed. All periods are then refit, and the null is rejected when the
    post-period residuals are atypically large relative to their permutations.
    ``permutation="moving_block"`` (cyclic shifts) is robust to serial correlation, but
    its smallest attainable p-value is ``1/T``, and it loses power when the post period
    is a large share of the record. ``"iid"`` permutes freely and assumes
    exchangeable residuals.
    """
    grid = center + scale * np.linspace(-width, width, n_grid)
    grid = np.unique(np.append(grid, 0.0))
    ps = np.array([_conformal_p(y1, Y0, pre, th, method, permutation) for th in grid])
    acc = grid[ps > alpha]
    p0 = float(ps[grid == 0.0][0])
    if acc.size == 0:
        lo = hi = np.nan
    else:
        lo, hi = float(acc.min()), float(acc.max())
        if lo <= grid.min():
            lo = -np.inf
        if hi >= grid.max():
            hi = np.inf
    min_p = 1.0 / len(pre) if permutation == "moving_block" else 1.0 / 1000
    return {"ci_low": lo, "ci_high": hi, "p_value": p0, "min_attainable_p": min_p, "permutation": permutation}
