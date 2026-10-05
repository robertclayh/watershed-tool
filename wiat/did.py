"""Difference-in-differences estimators for basin panels.

Three things matter for watershed attribution and drive the design here:

1. **Few treated basins.** Interventions are usually evaluated on one basin or a
   handful. Cluster-robust standard errors are unreliable with few treated clusters,
   and the wild cluster bootstrap is unreliable when only one cluster is treated
   (MacKinnon & Webb, 2017). The default inference is therefore the placebo/permutation
   approach of Conley & Taber (2011). Its resolution is set *mechanically* by the number
   of control basins, which makes the cost of a thin donor pool explicit.
2. **Climate confounding.** Precipitation is time-varying and basin-specific, so it
   enters as a covariate rather than being absorbed by fixed effects.
3. **Staggered adoption.** With multiple start dates, two-way fixed effects can be
   badly biased. :func:`staggered_att` provides the Callaway & Sant'Anna (2021)
   group-time estimator with never-treated comparisons.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from ._linalg import twoway_demean
from .panel import TIME, UNIT


# --------------------------------------------------------------------------- helpers
@dataclass
class _PanelArrays:
    units: list
    times: list
    Y: np.ndarray  # (N, T)
    X: np.ndarray  # (k, N, T)
    start_idx: dict  # unit -> index of first treated period (treated units only)

    @property
    def treated(self) -> list:
        return [u for u in self.units if u in self.start_idx]

    @property
    def controls(self) -> list:
        return [u for u in self.units if u not in self.start_idx]


def _start_index(times: list, start) -> int:
    if isinstance(times[0], pd.Timestamp):
        start = pd.Timestamp(start)
    for k, t in enumerate(times):
        if t >= start:
            return k
    return len(times)


def _arrays(
    df: pd.DataFrame,
    outcome: str,
    treat_start: Mapping,
    covariates: Sequence[str] = (),
    unit: str = UNIT,
    time: str = TIME,
) -> _PanelArrays:
    cols = [outcome, *covariates]
    sub = df[[unit, time, *cols]]
    wide = sub.pivot(index=unit, columns=time, values=cols)
    if wide.isna().any().any():
        bad = wide.index[wide.isna().any(axis=1)].tolist()
        raise ValueError(
            f"panel is unbalanced or has missing values for units {bad[:10]}; "
            "use wiat.panel.balance() or impute first"
        )
    units = wide.index.tolist()
    times = sorted(sub[time].unique().tolist())
    if isinstance(times[0], np.datetime64):
        times = [pd.Timestamp(t) for t in times]
    Y = wide[outcome].reindex(columns=times).to_numpy(float)
    X = (
        np.stack([wide[c].reindex(columns=times).to_numpy(float) for c in covariates])
        if covariates
        else np.zeros((0, *Y.shape))
    )
    start_idx = {}
    for u, s in treat_start.items():
        if u not in units:
            raise KeyError(f"treated unit {u!r} not found in panel")
        k = _start_index(times, s)
        if k == 0 or k >= len(times):
            raise ValueError(f"treated unit {u!r} needs at least one pre and one post period")
        start_idx[u] = k
    return _PanelArrays(units, times, Y, X, start_idx)


def _partial_out(y: np.ndarray, R: np.ndarray, X: np.ndarray):
    """Two-way demean then partial controls ``X`` out of ``y`` (NT,) and ``R`` (NT, p)."""
    if X.shape[0] == 0:
        return y, R, None
    x = twoway_demean(X).reshape(X.shape[0], -1).T
    Q = np.linalg.pinv(x)
    return y - x @ (Q @ y), R - x @ (Q @ R), (x, Q)


def _cr1(r: np.ndarray, e: np.ndarray, N: int, T: int, k_extra: int):
    """Cluster-robust (by unit) CR1 covariance for OLS of e-residualised data on r."""
    n, p = r.shape
    bread = np.linalg.pinv(r.T @ r)
    scores = (r * e[:, None]).reshape(N, T, p).sum(axis=1)
    meat = scores.T @ scores
    c = (N / (N - 1)) * ((n - 1) / max(n - p - k_extra, 1))
    return c * bread @ meat @ bread


def _webb_weights(rng, B: int, G: int) -> np.ndarray:
    vals = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5), math.sqrt(0.5), 1.0, math.sqrt(1.5)])
    return vals[rng.integers(0, 6, size=(B, G))]


def _bootstrap_weights(rng, B: int, G: int) -> tuple[np.ndarray, str]:
    if 2**G <= B:
        W = np.array(list(itertools.product([-1.0, 1.0], repeat=G)))
        return W, f"rademacher-enumerated ({2**G} draws)"
    return _webb_weights(rng, B, G), f"webb ({B} draws)"


def placebo_interval(estimate: float, placebo: np.ndarray, alpha: float):
    """Permutation p-value and the confidence interval obtained by inverting it.

    Placebo estimates come from assigning the treatment to untreated basins
    (Conley & Taber, 2011). The test rejects ``effect = theta`` when
    ``|estimate - theta|`` is among the ``alpha`` share of most extreme absolute
    placebo effects. Inverting it gives ``estimate ± c``, where ``c`` (the
    *detection threshold*) is the smallest absolute effect the design could
    declare significant. By construction, the interval excludes zero exactly when
    ``p <= alpha``.

    With fewer than ``1/alpha - 1`` placebos (9 at ``alpha = 0.10``, 19 at 0.05),
    nothing can be significant: ``c = inf`` and the interval is the whole real line.
    Just above that limit, ``c`` is the single most extreme placebo, so a lone
    unusual control basin sets the whole study's resolution.

    Returns ``(ci_low, ci_high, p_value, detection_threshold)``.
    """
    e = np.asarray(placebo, float)
    M = e.size
    p = (1 + np.sum(np.abs(e) >= abs(estimate))) / (M + 1)
    m = int(math.floor(alpha * (M + 1) + 1e-9)) - 1
    if m < 0:
        return -np.inf, np.inf, float(p), np.inf
    thresh = float(np.sort(np.abs(e))[::-1][m]) if m < M else 0.0
    return estimate - thresh, estimate + thresh, float(p), thresh


# --------------------------------------------------------------------------- 2x2 / TWFE
@dataclass
class DiDResult:
    """Result of :func:`did`. ``ci``/``p_value`` follow the primary ``inference`` method."""

    estimate: float
    se: float
    ci: tuple[float, float]
    p_value: float
    inference: str
    alpha: float
    n_treated: int
    n_control: int
    n_periods: int
    covariates: tuple = ()
    methods: dict = field(default_factory=dict)
    placebo_estimates: np.ndarray | None = None
    detection_threshold: float = np.nan

    @property
    def significant(self) -> bool:
        return bool(self.p_value <= self.alpha)

    def summary(self) -> pd.DataFrame:
        rows = []
        for name, m in self.methods.items():
            rows.append({"method": name, **{k: v for k, v in m.items() if k != "note"}})
        return pd.DataFrame(rows).set_index("method")


def did(
    df: pd.DataFrame,
    outcome: str,
    treat_start: Mapping,
    covariates: Sequence[str] = (),
    unit: str = UNIT,
    time: str = TIME,
    controls: Sequence | None = None,
    inference: str = "auto",
    alpha: float = 0.1,
    n_boot: int = 999,
    n_placebo: int = 2000,
    seed: int | None = 0,
    compute: Sequence[str] = ("cluster", "wild", "placebo"),
) -> DiDResult:
    """Two-way fixed-effects difference-in-differences.

    ``y_it = a_i + l_t + tau * D_it + X_it' b + e_it``.

    Parameters
    ----------
    treat_start
        Mapping ``{basin_id: first treated period}``. Units not in the mapping are
        controls (optionally restricted to ``controls``).
    inference
        ``"placebo"`` (Conley-Taber permutation), ``"wild"`` (restricted wild cluster
        bootstrap, inverted for a CI), ``"cluster"`` (CR1, t with G-1 df) or ``"auto"``,
        which uses placebo inference when fewer than 5 basins are treated and the wild
        bootstrap otherwise. All applicable methods are computed and kept in
        ``result.methods``.
    compute
        Which inference methods to compute. The wild bootstrap's test inversion is
        the expensive one; drop it for Monte Carlo work.
    """
    if controls is not None:
        keep = set(controls) | set(treat_start)
        df = df[df[unit].isin(keep)]
    P = _arrays(df, outcome, treat_start, covariates, unit, time)
    N, T = P.Y.shape
    D = np.zeros((N, T))
    for i, u in enumerate(P.units):
        if u in P.start_idx:
            D[i, P.start_idx[u]:] = 1.0
    if len(P.controls) < 1:
        raise ValueError("DiD needs at least one control basin")

    y, r, _ = _partial_out(twoway_demean(P.Y).ravel(), twoway_demean(D).reshape(-1, 1), P.X)
    d = r[:, 0]
    dd = d @ d
    tau = float(d @ y / dd)
    e = y - tau * d
    k = P.X.shape[0]
    V = _cr1(r, e, N, T, k)
    se = float(np.sqrt(V[0, 0]))
    methods: dict = {}

    tcrit = stats.t.ppf(1 - alpha / 2, N - 1)
    methods["cluster"] = {
        "estimate": tau,
        "se": se,
        "ci_low": tau - tcrit * se,
        "ci_high": tau + tcrit * se,
        "p_value": float(2 * stats.t.sf(abs(tau / se), N - 1)) if se > 0 else np.nan,
    }

    rng = np.random.default_rng(seed)
    if "wild" in compute:
        # Restricted wild cluster bootstrap (WCR), inverted for a confidence set.
        W, wdesc = _bootstrap_weights(rng, n_boot, N)
        xs = None
        if k:
            xs = twoway_demean(P.X).reshape(k, -1).T
            Qx = np.linalg.pinv(xs)
        c = (N / (N - 1)) * ((N * T - 1) / max(N * T - 1 - k, 1))

        def wild_p(tau0: float) -> float:
            u = (y - tau0 * d).reshape(N, T)
            Uw = twoway_demean(u[None] * W[:, :, None]).reshape(W.shape[0], -1)
            if xs is not None:
                Uw = Uw - (Uw @ Qx.T) @ xs.T
            tau_b = tau0 + Uw @ d / dd
            e_b = Uw - (tau_b - tau0)[:, None] * d[None]
            s = (e_b * d[None]).reshape(W.shape[0], N, T).sum(axis=2)
            se_b = np.sqrt(c * (s**2).sum(axis=1)) / dd
            t_obs = (tau - tau0) / se if se > 0 else 0.0
            with np.errstate(divide="ignore", invalid="ignore"):
                t_b = (tau_b - tau0) / se_b
            return float(np.mean(np.abs(t_b) >= abs(t_obs) - 1e-12))

        def invert(direction: int) -> float:
            step = max(se, 1e-8 * max(abs(tau), 1.0))
            inner, outer = tau, tau + direction * step
            n = 0
            while wild_p(outer) > alpha:
                inner, outer = outer, outer + direction * step * (2**n)
                n += 1
                if n > 30:
                    return direction * np.inf
            for _ in range(40):
                mid = 0.5 * (inner + outer)
                if wild_p(mid) > alpha:
                    inner = mid
                else:
                    outer = mid
            return 0.5 * (inner + outer)

        methods["wild"] = {
            "estimate": tau,
            "se": se,
            "ci_low": invert(-1),
            "ci_high": invert(+1),
            "p_value": wild_p(0.0),
            "note": wdesc,
        }

    # Placebo / permutation inference over control basins (Conley & Taber, 2011).
    ctrl_rows = [i for i, u in enumerate(P.units) if u not in P.start_idx]
    starts = [P.start_idx[u] for u in P.treated]
    n_tr = len(starts)
    placebo = None
    if "placebo" in compute and len(ctrl_rows) > n_tr:
        Yc = P.Y[ctrl_rows]
        Xc = P.X[:, ctrl_rows]
        J = len(ctrl_rows)
        if math.comb(J, n_tr) <= n_placebo:
            combos = [list(c) for c in itertools.combinations(range(J), n_tr)]
        else:
            combos = [list(rng.choice(J, n_tr, replace=False)) for _ in range(n_placebo)]
        Dp = np.zeros((len(combos), J, T))
        for b, combo in enumerate(combos):
            for j, s in zip(combo, starts):
                Dp[b, j, s:] = 1.0
        yc, _, xq = _partial_out(twoway_demean(Yc).ravel(), np.zeros((J * T, 1)), Xc)
        dpt = twoway_demean(Dp).reshape(len(combos), -1)
        if xq is not None:
            xx, Q = xq
            dpt = dpt - (dpt @ Q.T) @ xx.T
        placebo = (dpt @ yc) / np.einsum("ij,ij->i", dpt, dpt)
        lo, hi, p, thresh = placebo_interval(tau, placebo, alpha)
        methods["placebo"] = {
            "estimate": tau,
            "se": float(np.std(placebo, ddof=1)) if placebo.size > 1 else np.nan,
            "ci_low": lo,
            "ci_high": hi,
            "p_value": p,
            "detection_threshold": thresh,
            "n_placebo": int(placebo.size),
        }

    if inference == "auto":
        inference = "placebo" if (n_tr < 5 and "placebo" in methods) else "wild"
    if inference not in methods:
        raise ValueError(f"inference method {inference!r} unavailable (have {list(methods)})")
    m = methods[inference]
    return DiDResult(
        estimate=tau,
        se=m["se"],
        ci=(m["ci_low"], m["ci_high"]),
        p_value=m["p_value"],
        inference=inference,
        alpha=alpha,
        n_treated=n_tr,
        n_control=len(P.controls),
        n_periods=T,
        covariates=tuple(covariates),
        methods=methods,
        placebo_estimates=placebo,
        detection_threshold=methods.get("placebo", {}).get("detection_threshold", np.nan),
    )


# --------------------------------------------------------------------------- event study
@dataclass
class EventStudyResult:
    coefs: pd.DataFrame  # index rel_time; columns coef, se, ci_low, ci_high, [pl_low, pl_high]
    pretrend_f: float
    pretrend_p: float
    pretrend_placebo_p: float
    alpha: float

    @property
    def pretrend_ok(self) -> bool:
        p = self.pretrend_placebo_p if np.isfinite(self.pretrend_placebo_p) else self.pretrend_p
        return bool(np.isnan(p) or p > self.alpha)


def _event_regressors(N, T, rows_starts, bins):
    lo, hi = bins[0], bins[-1]
    R = np.zeros((len(bins), N, T))
    tt = np.arange(T)
    for i, s in rows_starts:
        rel = np.clip(tt - s, lo, hi)
        for b, e in enumerate(bins):
            R[b, i] = rel == e
    return R


def event_study(
    df: pd.DataFrame,
    outcome: str,
    treat_start: Mapping,
    window: tuple[int, int] = (-5, 5),
    covariates: Sequence[str] = (),
    unit: str = UNIT,
    time: str = TIME,
    controls: Sequence | None = None,
    alpha: float = 0.1,
    n_placebo: int = 2000,
    seed: int | None = 0,
) -> EventStudyResult:
    """Dynamic DiD with binned endpoints, reference period -1.

    Reports cluster-robust coefficients and a joint pre-trend Wald test. When fewer than
    5 basins are treated it also reports placebo bands and a placebo-based pre-trend
    test, built by assigning the treated start dates to control basins.
    """
    if controls is not None:
        keep = set(controls) | set(treat_start)
        df = df[df[unit].isin(keep)]
    P = _arrays(df, outcome, treat_start, covariates, unit, time)
    N, T = P.Y.shape
    lo, hi = window
    bins = [e for e in range(lo, hi + 1) if e != -1]
    rows_starts = [(i, P.start_idx[u]) for i, u in enumerate(P.units) if u in P.start_idx]
    R = _event_regressors(N, T, rows_starts, bins)
    # Drop bins with no support (e.g., window longer than the panel).
    support = R.reshape(len(bins), -1).sum(axis=1) > 0
    bins = [b for b, s in zip(bins, support) if s]
    R = R[support]
    y, r, _ = _partial_out(twoway_demean(P.Y).ravel(), twoway_demean(R).reshape(len(bins), -1).T, P.X)
    beta, *_ = np.linalg.lstsq(r, y, rcond=None)
    e = y - r @ beta
    V = _cr1(r, e, N, T, P.X.shape[0])
    se = np.sqrt(np.clip(np.diag(V), 0, None))
    tcrit = stats.t.ppf(1 - alpha / 2, N - 1)
    coefs = pd.DataFrame(
        {"coef": beta, "se": se, "ci_low": beta - tcrit * se, "ci_high": beta + tcrit * se},
        index=pd.Index(bins, name="rel_time"),
    )
    pre = [k for k, b in enumerate(bins) if b < -1]
    f_stat = p_f = np.nan
    if pre:
        bp = beta[pre]
        Vp = V[np.ix_(pre, pre)]
        f_stat = float(bp @ np.linalg.pinv(Vp) @ bp / len(pre))
        p_f = float(stats.f.sf(f_stat, len(pre), max(N - 1, 1)))

    p_pl = np.nan
    n_tr = len(rows_starts)
    ctrl_rows = [i for i, u in enumerate(P.units) if u not in P.start_idx]
    if n_tr < 5 and len(ctrl_rows) > n_tr:
        rng = np.random.default_rng(seed)
        J = len(ctrl_rows)
        starts = [s for _, s in rows_starts]
        if math.comb(J, n_tr) <= n_placebo:
            combos = [list(c) for c in itertools.combinations(range(J), n_tr)]
        else:
            combos = [list(rng.choice(J, n_tr, replace=False)) for _ in range(n_placebo)]
        Yc = P.Y[ctrl_rows]
        Xc = P.X[:, ctrl_rows]
        yc, _, xq = _partial_out(twoway_demean(Yc).ravel(), np.zeros((J * T, 1)), Xc)
        pl = []
        for combo in combos:
            Rp = _event_regressors(J, T, list(zip(combo, starts)), bins)
            rp = twoway_demean(Rp).reshape(len(bins), -1).T
            if xq is not None:
                xx, Q = xq
                rp = rp - xx @ (Q @ rp)
            bp, *_ = np.linalg.lstsq(rp, yc, rcond=None)
            pl.append(bp)
        pl = np.array(pl)
        # Placebo bands are centred on zero: the range of "effects" seen in untreated basins.
        coefs["pl_low"] = np.quantile(pl, alpha / 2, axis=0)
        coefs["pl_high"] = np.quantile(pl, 1 - alpha / 2, axis=0)
        if pre:
            # Normalise each lead by its placebo spread so no single lead dominates.
            scale = pl[:, pre].std(axis=0, ddof=1)
            scale[scale == 0] = 1.0
            s_obs = np.sum((beta[pre] / scale) ** 2)
            s_pl = np.sum((pl[:, pre] / scale) ** 2, axis=1)
            p_pl = float((1 + np.sum(s_pl >= s_obs)) / (len(s_pl) + 1))
    return EventStudyResult(coefs, f_stat, p_f, p_pl, alpha)


# --------------------------------------------------------------------------- staggered
@dataclass
class StaggeredResult:
    att: float
    se: float
    ci: tuple[float, float]
    group_time: pd.DataFrame
    dynamic: pd.DataFrame
    n_boot: int


def staggered_att(
    df: pd.DataFrame,
    outcome: str,
    treat_start: Mapping,
    unit: str = UNIT,
    time: str = TIME,
    alpha: float = 0.1,
    n_boot: int = 999,
    seed: int | None = 0,
) -> StaggeredResult:
    """Callaway & Sant'Anna (2021) ATT(g, t) with never-treated comparison basins.

    Uses a universal base period (g - 1) for every cohort. The overall ATT averages
    post-treatment cells weighted by cohort size. Inference is a unit bootstrap,
    stratified by cohort.
    """
    P = _arrays(df, outcome, treat_start, (), unit, time)
    N, T = P.Y.shape
    cohort = np.full(N, -1)
    for i, u in enumerate(P.units):
        if u in P.start_idx:
            cohort[i] = P.start_idx[u]
    never = np.where(cohort < 0)[0]
    if never.size < 2:
        raise ValueError("staggered_att needs at least two never-treated basins")
    groups = sorted(set(cohort[cohort >= 0]))

    def estimate(idx: np.ndarray):
        Y = P.Y[idx]
        coh = cohort[idx]
        nev = coh < 0
        cells = []
        for g in groups:
            m = coh == g
            if m.sum() == 0 or nev.sum() == 0:
                continue
            base = g - 1
            dy_t = Y[m] - Y[m][:, [base]]
            dy_c = Y[nev] - Y[nev][:, [base]]
            att_gt = dy_t.mean(axis=0) - dy_c.mean(axis=0)
            for t in range(T):
                if t == base:
                    continue
                cells.append((g, t, t - g, att_gt[t], m.sum()))
        cdf = pd.DataFrame(cells, columns=["group", "t", "rel_time", "att", "n"])
        post = cdf[cdf.rel_time >= 0]
        overall = np.average(post.att, weights=post.n) if len(post) else np.nan
        dyn = cdf.groupby("rel_time").apply(lambda g: np.average(g.att, weights=g.n), include_groups=False)
        return overall, cdf, dyn

    att, cdf, dyn = estimate(np.arange(N))
    rng = np.random.default_rng(seed)
    strata = [np.where(cohort == g)[0] for g in groups] + [never]
    boots, dyn_boots = [], []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(s, s.size, replace=True) for s in strata])
        a, _, d = estimate(idx)
        boots.append(a)
        dyn_boots.append(d.reindex(dyn.index))
    boots = np.array(boots)
    dyn_b = np.array(dyn_boots)
    se = float(np.nanstd(boots, ddof=1))
    ci = (float(np.nanquantile(boots, alpha / 2)), float(np.nanquantile(boots, 1 - alpha / 2)))
    dyn_df = pd.DataFrame(
        {
            "att": dyn.values,
            "se": np.nanstd(dyn_b, axis=0, ddof=1),
            "ci_low": np.nanquantile(dyn_b, alpha / 2, axis=0),
            "ci_high": np.nanquantile(dyn_b, 1 - alpha / 2, axis=0),
        },
        index=pd.Index(dyn.index, name="rel_time"),
    )
    cdf["time"] = [P.times[t] for t in cdf.t]
    cdf["group"] = [P.times[g] for g in cdf.group]
    return StaggeredResult(float(att), se, ci, cdf.drop(columns="t"), dyn_df, n_boot)
