"""Quantifying where attribution fails: confidence versus donor-basin scarcity.

Comparative designs need comparison basins. In practice the basins where
nature-based investment is most needed are the ones with the fewest credible
controls:

* Small, steep headwater regions have few basins of the right size.
* Neighbouring basins share the intervention through hydrologic connectivity, so they
  are excluded to avoid spillovers.
* Other basins carry their own programmes (contamination).
* Persistent cloud cover leaves basins with gaps in the optical record.
* Degraded, high-rainfall basins are atypical: they sit at the edge of the regional
  distribution, outside the donors' convex hull.

Two tools are provided:

* :func:`run_scarcity_experiment` is a Monte Carlo over the *size* of the donor pool
  and the *atypicality* (``extremity``) of the treated basin, with known truth. It
  measures bias, RMSE, interval coverage, interval finiteness, power, size, and the
  probability of a "confident, correct attribution".
* :func:`donor_thinning` is an empirical check on a real panel. It re-estimates the
  effect while restricting the donor pool to the k best-matched (or k random)
  basins, showing how much a study's conclusion depends on donor availability.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from .did import did
from .panel import aggregate
from .simulate import simulate_watersheds
from .synth import synthetic_control

ESTIMATORS = ("did_placebo", "did_cluster", "scm", "ascm")


def _estimate_all(annual: pd.DataFrame, outcome: str, treated: str, start: int, alpha: float,
                  covariates: Sequence[str], estimators: Sequence[str]) -> dict:
    out = {}
    if {"did_placebo", "did_cluster"} & set(estimators):
        r = did(annual, outcome, {treated: start}, covariates=covariates, alpha=alpha,
                inference="cluster", compute=("cluster", "placebo"))
        if "did_cluster" in estimators:
            m = r.methods["cluster"]
            out["did_cluster"] = (r.estimate, m["ci_low"], m["ci_high"], m["p_value"])
        if "did_placebo" in estimators:
            m = r.methods.get("placebo")
            if m is None:
                out["did_placebo"] = (r.estimate, -np.inf, np.inf, 1.0)
            else:
                out["did_placebo"] = (r.estimate, m["ci_low"], m["ci_high"], m["p_value"])
    for name, method in (("scm", "demeaned"), ("ascm", "augmented")):
        if name in estimators:
            s = synthetic_control(annual, outcome, treated, start, method=method, alpha=alpha,
                                  placebo_fit_cutoff=None)
            ci = s.ci if s.placebo is not None else (-np.inf, np.inf)
            p = s.p_value if s.placebo is not None else 1.0
            out[name] = (s.att, ci[0], ci[1], p, s.diagnostics["pre_rmspe_normalized"],
                         s.diagnostics["frac_pre_outside_donor_range"])
    return out


def _one_rep(args) -> list[dict]:
    (n_donors, extremity, rep, seed, outcome, alpha, n_years, start_offset,
     ndvi_effect, covariates, estimators) = args
    rs = np.random.SeedSequence([seed, n_donors, int(extremity * 100), rep]).generate_state(1)[0]
    start = 2001 + start_offset
    sim = simulate_watersheds(
        n_basins=n_donors + 1, n_years=n_years, treat_start={"B000": start},
        ndvi_effect=ndvi_effect, treated_extremity=extremity, seed=int(rs),
    )
    truth = sim.true_att(outcome)
    rows = []
    for world, panel, true in (("effect", sim.panel, truth), ("no_effect", sim.counterfactual, 0.0)):
        annual = aggregate(panel, "year")
        res = _estimate_all(annual, outcome, "B000", start, alpha, covariates, estimators)
        for est, vals in res.items():
            e, lo, hi, p = vals[:4]
            row = {
                "n_donors": n_donors, "extremity": extremity, "rep": rep, "world": world,
                "estimator": est, "truth": true, "estimate": e, "ci_low": lo, "ci_high": hi,
                "p_value": p,
            }
            if len(vals) > 4:
                row["pre_rmspe_norm"], row["outside_range"] = vals[4], vals[5]
            rows.append(row)
    return rows


def run_scarcity_experiment(
    donor_counts: Sequence[int] = (2, 4, 6, 9, 14, 19, 29, 49),
    extremities: Sequence[float] = (0.0, 1.0, 2.0),
    n_reps: int = 200,
    outcome: str = "et_mm",
    alpha: float = 0.1,
    n_years: int = 24,
    start_offset: int = 14,
    ndvi_effect: float = 0.08,
    covariates: Sequence[str] = ("precip_mm",),
    estimators: Sequence[str] = ESTIMATORS,
    seed: int = 20260101,
    n_jobs: int | None = None,
) -> pd.DataFrame:
    """Monte Carlo of attribution quality versus donor-pool size and treated atypicality.

    Each replication simulates a region of ``n_donors + 1`` basins (see
    :func:`wiat.simulate.simulate_watersheds`). The region is observed twice with
    identical draws: once with the intervention (``world="effect"``) and once without
    (``world="no_effect"``). Running every estimator in both worlds measures power and size
    on the same data. Returns one row per (replication, world, estimator). Pass the
    result to :func:`summarize_scarcity`.
    """
    jobs = [
        (J, x, r, seed, outcome, alpha, n_years, start_offset, ndvi_effect, tuple(covariates), tuple(estimators))
        for J in donor_counts for x in extremities for r in range(n_reps)
    ]
    n_jobs = n_jobs or max(1, (os.cpu_count() or 1))
    if n_jobs == 1:
        rows = [row for j in jobs for row in _one_rep(j)]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            rows = [row for res in ex.map(_one_rep, jobs, chunksize=8) for row in res]
    df = pd.DataFrame(rows)
    df.attrs.update({"outcome": outcome, "alpha": alpha})
    return df


def summarize_scarcity(results: pd.DataFrame, alpha: float | None = None) -> pd.DataFrame:
    """Aggregate Monte Carlo rows into degradation curves.

    Columns per (estimator, extremity, n_donors):

    * ``bias``, ``rmse``, ``rel_rmse``: point-estimate error (relative to mean |truth|)
    * ``coverage``: share of intervals containing the truth (target ``1 - alpha``)
    * ``ci_bounded``: share of intervals that are finite at all
    * ``median_rel_width``: median interval width / mean |truth| (bounded intervals only)
    * ``power``: rejection rate of "no effect" when the effect is real
    * ``size``: rejection rate when there is no effect (target ``alpha``)
    * ``confident_correct``: share of replications reaching a bounded interval, a
      significant result and the correct sign. This is the practical probability
      that a study can make the attribution claim at all.
    * ``min_attainable_p``: ``1 / (J + 1)`` for placebo-based inference
    """
    alpha = alpha if alpha is not None else results.attrs.get("alpha", 0.1)
    keys = ["estimator", "extremity", "n_donors"]
    eff = results[results.world == "effect"].copy()
    nul = results[results.world == "no_effect"].copy()
    eff["err"] = eff.estimate - eff.truth
    eff["covered"] = (eff.ci_low <= eff.truth) & (eff.truth <= eff.ci_high)
    eff["bounded"] = np.isfinite(eff.ci_low) & np.isfinite(eff.ci_high)
    eff["width"] = np.where(eff.bounded, eff.ci_high - eff.ci_low, np.nan)
    eff["reject"] = eff.p_value <= alpha
    eff["correct"] = eff.reject & eff.bounded & (np.sign(eff.estimate) == np.sign(eff.truth))
    nul["reject"] = nul.p_value <= alpha
    g = eff.groupby(keys)
    scale = g.truth.apply(lambda x: np.mean(np.abs(x)))
    out = pd.DataFrame(
        {
            "truth": g.truth.mean(),
            "bias": g.err.mean(),
            "rmse": g.err.apply(lambda x: float(np.sqrt(np.mean(x**2)))),
            "coverage": g.covered.mean(),
            "ci_bounded": g.bounded.mean(),
            "median_rel_width": g.width.median() / scale,
            "power": g.reject.mean(),
            "confident_correct": g.correct.mean(),
            "n_reps": g.size(),
        }
    )
    out["rel_rmse"] = out.rmse / scale
    out["size"] = nul.groupby(keys).reject.mean()
    if "pre_rmspe_norm" in eff:
        out["pre_fit_rmspe_norm"] = g.pre_rmspe_norm.median()
        out["outside_donor_range"] = g.outside_range.mean()
    out["min_attainable_p"] = 1.0 / (out.index.get_level_values("n_donors") + 1)
    return out.reset_index()


def donors_needed(alpha: float) -> dict:
    """Donor counts below which placebo-based inference cannot succeed, by construction."""
    return {
        # 1/(J+1) <= alpha: below this, placebo p-values cannot reach alpha and
        # placebo-inverted intervals are unbounded.
        "significance_possible": int(np.ceil(1 / alpha - 1 - 1e-9)),
        # Below this, the detection threshold is set by one or two extreme placebos.
        "stable_threshold": int(np.ceil(2 / alpha - 1 - 1e-9)),
    }


def donor_thinning(
    df: pd.DataFrame,
    outcome: str,
    treated,
    treat_start,
    ranked_donors: Sequence,
    sizes: Sequence[int] = (3, 5, 8, 12, 19, 29),
    mode: str = "nearest",
    n_draws: int = 50,
    alpha: float = 0.1,
    covariates: Sequence[str] = (),
    unit: str = "basin_id",
    time: str = "date",
    seed: int = 0,
) -> pd.DataFrame:
    """Re-estimate an observed study with progressively smaller donor pools.

    Parameters
    ----------
    ranked_donors
        Donor basins ordered from best to worst match, e.g.
        ``selection.nearest(len(selection.candidates))``.
    mode
        ``"nearest"`` keeps the ``k`` best-matched donors (one estimate per ``k``).
        ``"random"`` draws ``n_draws`` random ``k``-subsets, showing how much the
        answer depends on *which* basins happen to be available.
    """
    rng = np.random.default_rng(seed)
    ranked = list(ranked_donors)
    rows = []
    for k in sizes:
        if k > len(ranked):
            continue
        subsets = [ranked[:k]] if mode == "nearest" else [
            list(rng.choice(ranked, k, replace=False)) for _ in range(n_draws)
        ]
        for draw, pool in enumerate(subsets):
            sub = df[df[unit].isin([treated, *pool])]
            r = did(sub, outcome, {treated: treat_start}, covariates=covariates, unit=unit,
                    time=time, alpha=alpha, inference="placebo" if k > 1 else "cluster",
                    compute=("cluster", "placebo"))
            s = synthetic_control(sub, outcome, treated, treat_start, donors=pool, unit=unit,
                                  time=time, alpha=alpha)
            rows.append({
                "k": k, "draw": draw,
                "did_estimate": r.estimate, "did_ci_low": r.ci[0], "did_ci_high": r.ci[1],
                "did_p": r.p_value, "did_detection_threshold": r.detection_threshold,
                "scm_estimate": s.att, "scm_ci_low": s.ci[0], "scm_ci_high": s.ci[1],
                "scm_p": s.p_value, "scm_pre_rmspe_norm": s.diagnostics["pre_rmspe_normalized"],
            })
    return pd.DataFrame(rows)
