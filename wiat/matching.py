"""Donor-basin matching: who is a credible counterfactual for the treated basin?

Matching runs on *basin signatures*: pre-intervention climate and vegetation
summaries, plus static physiography such as elevation, slope and forest cover.
A donor is "suitable" when two conditions hold: its root-mean-square standardised
difference from the treated basin is within ``caliper``, and no single covariate
differs by more than ``hard_caliper``. The number of suitable donors is the key input to the scarcity analysis.
A basin that is wetter, steeper or more degraded than anything in the region has
few suitable donors, and attribution there is correspondingly weak.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .panel import TIME, UNIT


def basin_signatures(
    df: pd.DataFrame,
    baseline_end,
    variables: Sequence[str] = ("precip_mm", "et_mm", "ndvi"),
    trend_variables: Sequence[str] = ("ndvi",),
    seasonality_variables: Sequence[str] = ("precip_mm",),
    unit: str = UNIT,
    time: str = TIME,
) -> pd.DataFrame:
    """Summarise each basin's pre-intervention behaviour.

    * ``<var>_mean``: mean over the baseline
    * ``<var>_trend``: OLS slope per period over the baseline (pre-trends matter for DiD)
    * ``<var>_seasonality``: coefficient of variation of the month-of-year climatology
      (monthly panels only)
    * ``aridity``: mean ET / mean P (evaporative index), when both are present
    """
    d = df.copy()
    is_dt = not np.issubdtype(d[time].dtype, np.number)
    tt = pd.to_datetime(d[time]) if is_dt else d[time]
    cutoff = pd.Timestamp(baseline_end) if is_dt else baseline_end
    base = d[tt < cutoff].copy()
    base["_t"] = (pd.to_datetime(base[time]).dt.year + pd.to_datetime(base[time]).dt.month / 12) if is_dt else base[time]
    out = {}
    g = base.groupby(unit)
    for v in variables:
        if v in base:
            out[f"{v}_mean"] = g[v].mean()
    for v in trend_variables:
        if v in base:
            out[f"{v}_trend"] = g.apply(
                lambda x: np.polyfit(x["_t"], x[v], 1)[0] if x[v].notna().sum() > 2 else np.nan,
                include_groups=False,
            )
    if is_dt:
        base["_m"] = pd.to_datetime(base[time]).dt.month
        for v in seasonality_variables:
            if v in base:
                clim = base.groupby([unit, "_m"])[v].mean().unstack()
                out[f"{v}_seasonality"] = clim.std(axis=1) / clim.mean(axis=1).replace(0, np.nan)
    if {"precip_mm", "et_mm"} <= set(base.columns):
        out["aridity"] = g["et_mm"].mean() / g["precip_mm"].mean()
    return pd.DataFrame(out)


@dataclass
class DonorSelection:
    treated: list
    candidates: list
    suitable: list
    distances: pd.DataFrame  # Mahalanobis distance: rows candidates, cols treated
    std_diffs: dict  # treated -> DataFrame (candidates x covariates)
    caliper: float
    hard_caliper: float
    covariates: list
    excluded: dict

    def _ok(self, t) -> pd.Series:
        sd = self.std_diffs[t]
        rms = np.sqrt((sd**2).mean(axis=1))
        return (rms <= self.caliper) & (sd.abs() <= self.hard_caliper).all(axis=1)

    @property
    def n_suitable(self) -> int:
        return len(self.suitable)

    def nearest(self, k: int, treated=None) -> list:
        """The ``k`` closest candidates (by Mahalanobis distance) to ``treated``."""
        col = self.distances[treated if treated is not None else self.treated[0]]
        return col.sort_values().index[:k].tolist()

    def summary(self) -> pd.DataFrame:
        rows = []
        for t in self.treated:
            sd = self.std_diffs[t]
            ok = self._ok(t)
            rows.append(
                {
                    "treated": t,
                    "candidates": len(self.candidates),
                    "suitable": int(ok.sum()),
                    "median_rms_diff": float(np.sqrt((sd**2).mean(axis=1)).median()),
                    # Covariate on which the treated basin is most atypical.
                    "most_atypical_covariate": sd.abs().median().idxmax(),
                }
            )
        return pd.DataFrame(rows).set_index("treated")


def select_donors(
    features: pd.DataFrame,
    treated: Iterable,
    covariates: Sequence[str] | None = None,
    candidates: Iterable | None = None,
    caliper: float = 1.0,
    hard_caliper: float = 2.5,
    exclude: Mapping | None = None,
    max_donors: int | None = None,
) -> DonorSelection:
    """Select suitable donor basins by calipered covariate matching.

    Parameters
    ----------
    features
        One row per basin (index = basin id), e.g. :func:`basin_signatures` joined
        with static attributes.
    caliper
        Maximum root-mean-square standardised difference (in pooled SDs across all
        basins) between donor and treated basin. Scale-free in the number of covariates.
    hard_caliper
        Maximum absolute standardised difference allowed on *any single* covariate.
    exclude
        ``{basin_id: reason}`` for basins that must not serve as donors, e.g. basins
        hydrologically connected to the treated basin (spillovers), or basins with
        their own interventions (contamination).
    max_donors
        Keep at most this many suitable donors, nearest by Mahalanobis distance.
    """
    treated = list(treated)
    covariates = list(covariates or features.columns)
    F = features[covariates].astype(float)
    exclude = dict(exclude or {})
    pool = list(candidates) if candidates is not None else [u for u in F.index if u not in treated]
    excluded = {u: exclude[u] for u in pool if u in exclude}
    missing = F.loc[pool].isna().any(axis=1)
    for u in missing[missing].index:
        excluded.setdefault(u, "missing covariates")
    pool = [u for u in pool if u not in excluded and u not in treated]

    sd = F.std(ddof=1).replace(0, 1.0)
    Z = (F - F.mean()) / sd
    cov = np.atleast_2d(np.cov(Z.dropna().to_numpy().T))
    inv = np.linalg.pinv(cov)
    dist = {}
    std_diffs = {}
    for t in treated:
        diff = Z.loc[pool] - Z.loc[t]
        std_diffs[t] = diff
        dist[t] = pd.Series(np.sqrt(np.einsum("ij,jk,ik->i", diff.values, inv, diff.values)), index=pool)
    distances = pd.DataFrame(dist)
    sel = DonorSelection(treated, pool, [], distances, std_diffs, caliper, hard_caliper, covariates, excluded)
    suitable_set: list = []
    for t in treated:
        ok = sel._ok(t)
        chosen = distances.loc[ok[ok].index, t].sort_values()
        if max_donors is not None:
            chosen = chosen.iloc[:max_donors]
        suitable_set.extend(u for u in chosen.index if u not in suitable_set)
    sel.suitable = suitable_set
    return sel


def balance_table(
    features: pd.DataFrame,
    treated,
    donors: Sequence,
    weights: pd.Series | None = None,
    covariates: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Standardised mean differences between the treated basin and (weighted) donors.

    |SMD| < 0.25 is a common rule of thumb for adequate balance.
    """
    covariates = list(covariates or features.columns)
    F = features[covariates].astype(float)
    sd = F.std(ddof=1).replace(0, 1.0)
    t_vals = F.loc[treated]
    pool = F.loc[list(donors)]
    unweighted = pool.mean()
    out = pd.DataFrame(
        {
            "treated": t_vals,
            "donor_mean": unweighted,
            "smd_unweighted": (t_vals - unweighted) / sd,
            "pool_min": pool.min(),
            "pool_max": pool.max(),
        }
    )
    if weights is not None:
        w = weights.reindex(pool.index).fillna(0.0)
        wm = (pool.mul(w, axis=0)).sum() / max(w.sum(), 1e-12)
        out["synthetic"] = wm
        out["smd_weighted"] = (t_vals - wm) / sd
    out["outside_pool_range"] = (t_vals < out.pool_min) | (t_vals > out.pool_max)
    return out


def hydrologically_connected(
    basins: pd.DataFrame, target, id_col: str = "HYBAS_ID", next_col: str = "NEXT_DOWN"
) -> set:
    """Basins upstream or downstream of ``target`` in a HydroBASINS-style table.

    Such basins share water, sediment, and often the intervention itself, so they
    violate SUTVA as controls and should be excluded from the donor pool.
    """
    nxt = dict(zip(basins[id_col], basins[next_col]))
    down = set()
    cur = nxt.get(target)
    while cur and cur in nxt and cur not in down:
        down.add(cur)
        cur = nxt.get(cur)
    up_map: dict = {}
    for b, n in nxt.items():
        up_map.setdefault(n, []).append(b)
    up, stack = set(), [target]
    while stack:
        for b in up_map.get(stack.pop(), []):
            if b not in up:
                up.add(b)
                stack.append(b)
    return (up | down) - {target}
