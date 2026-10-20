"""Validation on real panels: inject a known effect, see whether the design recovers it.

The Monte Carlo in :mod:`wiat.scarcity` uses simulated basins, so its conclusions are
only as good as the simulator. This module repeats the experiment on an observed
panel. Each untreated basin in turn plays the treated basin. Its donors are its
``k`` best-matched real neighbours, ranked exactly as :class:`AttributionStudy`
ranks them, and a known step effect is added to its outcome after the real
programme's start year. Everything else (climate shocks, sensor noise, cloud gaps,
ET model artefacts, the real spread of basin types) is the region's own.

Two worlds are estimated on the same data, as in the simulation:

* ``no_effect``: the panel as observed. Untreated basins have no intervention, so
  any rejection is a false positive (``size``).
* ``effect``: the outcome raised by ``effect`` times the basin's own pre-period
  mean from the start year on (``power``, ``coverage``, ``confident_correct``).

Rows match :func:`wiat.scarcity.run_scarcity_experiment`, so
:func:`wiat.scarcity.summarize_scarcity` and :func:`wiat.report.plot_scarcity`
work unchanged. ``extremity`` holds the pseudo-treated basin's atypicality group
(0 = most typical third of the region, 2 = least typical third), measured as the
RMS standardised distance of its matching covariates from the regional mean.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from .attribution import AttributionStudy, _to_period
from .matching import select_donors
from .panel import aggregate, balance
from .scarcity import _estimate_all

ATYPICALITY_LABELS = {0.0: "typical third", 1.0: "middle third", 2.0: "atypical third"}


def prepare_region(
    panel: pd.DataFrame,
    outcome: str,
    treat_start: Mapping,
    static: pd.DataFrame | None = None,
    exclude: Mapping | None = None,
    covariates: Sequence[str] = ("precip_mm",),
) -> dict:
    """Annual balanced panel, matching features and eligible basins for one region.

    Eligible basins are those with complete records that are neither treated nor in
    ``exclude`` (hydrologically connected to, or contaminated by, the real programme).
    """
    study = AttributionStudy(panel, outcome, treat_start, static=static)
    feats = study._features()
    cov = study._match_covariates(feats)
    annual = aggregate(panel, "year")
    annual, _ = balance(annual, [outcome, *covariates])
    blocked = set(treat_start) | set(exclude or {})
    eligible = [u for u in annual.basin_id.unique()
                if u not in blocked and u in feats.index and feats.loc[u, cov].notna().all()]
    F = feats.loc[eligible, cov].astype(float)
    Z = (F - F.mean()) / F.std(ddof=1).replace(0, 1.0)
    atypicality = np.sqrt((Z**2).mean(axis=1))
    start = _to_period(min(treat_start.values()), "year", 1)
    return {"annual": annual[annual.basin_id.isin(eligible)].reset_index(drop=True), "features": feats,
            "covariates": cov, "eligible": eligible, "atypicality": atypicality, "start": start,
            "outcome": outcome, "did_covariates": tuple(covariates)}


def _one_basin(args) -> list[dict]:
    region, basin, group, donor_counts, effects, alpha, estimators = args
    annual, outcome, start = region["annual"], region["outcome"], region["start"]
    pool = [u for u in region["eligible"] if u != basin]
    sel = select_donors(region["features"], [basin], covariates=region["covariates"], candidates=pool,
                        caliper=np.inf, hard_caliper=np.inf)
    ranked = sel.distances[basin].sort_values().index.tolist()
    own = annual[annual.basin_id == basin]
    base = float(own.loc[own.date < start, outcome].mean())
    post = (annual.basin_id == basin) & (annual.date >= start)
    rows = []
    for k in donor_counts:
        if k > len(ranked):
            continue
        keep = annual.basin_id.isin([basin, *ranked[:k]])
        worlds = [("no_effect", annual, 0.0, 0.0)]
        for e in effects:
            shifted = annual.copy()
            shifted.loc[post, outcome] += e * base
            worlds.append(("effect", shifted, e * base, e))
        for world, data, truth, rel in worlds:
            res = _estimate_all(data[keep], outcome, basin, start, alpha, region["did_covariates"], estimators)
            for est, vals in res.items():
                e_, lo, hi, p = vals[:4]
                row = {"n_donors": k, "extremity": group, "rep": basin, "world": world, "effect": rel,
                       "estimator": est, "truth": truth, "estimate": e_, "ci_low": lo, "ci_high": hi,
                       "p_value": p, "atypicality": float(region["atypicality"][basin])}
                if len(vals) > 4:
                    row["pre_rmspe_norm"], row["outside_range"] = vals[4], vals[5]
                rows.append(row)
    return rows


def injection_experiment(
    region: dict,
    donor_counts: Sequence[int] = (4, 6, 9, 14, 19, 29, 49),
    effects: Sequence[float] = (0.05, 0.10),
    n_basins: int | None = None,
    alpha: float = 0.1,
    estimators: Sequence[str] = ("did_placebo", "did_cluster", "scm"),
    seed: int = 0,
    n_jobs: int | None = None,
) -> pd.DataFrame:
    """Run the injected-effect experiment over (a sample of) a region's eligible basins.

    ``effects`` are fractions of each pseudo-treated basin's own pre-period mean
    outcome. Returns one row per (basin, donor count, world, effect, estimator).
    Summarise one effect size at a time with :func:`summarize_injection`.
    """
    basins = list(region["eligible"])
    if n_basins is not None and n_basins < len(basins):
        basins = list(np.random.default_rng(seed).choice(basins, n_basins, replace=False))
    groups = pd.qcut(region["atypicality"].rank(method="first"), 3, labels=False).astype(float)
    jobs = [(region, b, float(groups[b]), tuple(donor_counts), tuple(effects), alpha, tuple(estimators))
            for b in basins]
    n_jobs = n_jobs or max(1, (os.cpu_count() or 1))
    if n_jobs == 1:
        rows = [r for j in jobs for r in _one_basin(j)]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            rows = [r for res in ex.map(_one_basin, jobs) for r in res]
    df = pd.DataFrame(rows)
    df.attrs.update({"outcome": region["outcome"], "alpha": alpha})
    return df


def summarize_injection(results: pd.DataFrame, effect: float, alpha: float = 0.1) -> pd.DataFrame:
    """:func:`summarize_scarcity` for one injected effect size (the null world is shared)."""
    from .scarcity import summarize_scarcity

    sub = results[(results.world == "no_effect") | np.isclose(results.effect, effect)]
    out = summarize_scarcity(sub, alpha=alpha)
    out.insert(1, "effect", effect)
    return out
