"""Socio Bosque (Ecuador): matching, staggered ATT, event study, dose response.

Runs on the outputs of examples/socio_bosque.py (prepare, then extract/assemble):

    python examples/socio_bosque_analysis.py --data data/socio_bosque --out docs/results --fig-dir docs/figures

Design (docs/socio_bosque_plan.md): HydroBASINS level-12 basins; treated when at least
20% of the basin is under agreement (first full year after crossing); never-enrolled
controls (<= 1%) matched within region on 2001-07 baseline characteristics;
Callaway-Sant'Anna ATT with a basin bootstrap stratified by cohort.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from wiat.did import staggered_att
from wiat.matching import select_donors
from wiat.report import BLUE, GRID, MUTED, ORANGE, PLACEBO, _plt

BASELINE = (2001, 2007)
MATCH_ON = ["forest_2001", "loss_rate_base", "loss_trend_base", "ndvi_base", "et_base", "precip_base",
            "elevation_m", "slope_deg"]
OUTCOMES = {  # column: (label, unit)
    "forest_loss_pct": ("Forest loss", "% of 2000 forest / yr"),
    "forest_cover_pct": ("Forest cover", "% of basin"),
    "ndvi": ("NDVI", "annual median"),
    "et_mm": ("ET", "mm / yr"),
    "evaporative_index": ("ET / P", "ratio"),
    "water_yield_mm": ("P − ET", "mm / yr"),
    "precip_mm": ("Precipitation (placebo)", "mm / yr"),
}


def load(d: Path):
    panel = pd.read_csv(d / "panel.csv", parse_dates=["date"])
    forest = pd.read_csv(d / "forest.csv", parse_dates=["date"])
    p = panel.merge(forest, on=["basin_id", "date"], how="left")
    p["year"] = p.date.dt.year
    p["forest_loss_pct"] = 100 * p.forest_loss_frac
    p["forest_cover_pct"] = 100 * p.forest_cover_frac
    static = pd.read_csv(d / "static.csv", index_col=0)
    regions = pd.read_csv(d / "regions.csv", index_col=0)
    treat = pd.read_csv(d / "treatments.csv")
    ctrl = pd.read_csv(d / "controls.csv").basin_id.tolist()
    expo = pd.read_csv(d / "exposure.csv", index_col=0)
    expo.columns = expo.columns.astype(int)
    return p, static, regions, dict(zip(treat.basin_id, treat.start)), ctrl, expo


def baseline_features(p: pd.DataFrame, static: pd.DataFrame) -> pd.DataFrame:
    b = p[p.year.between(*BASELINE)]
    g = b.groupby("basin_id")

    def slope(x):
        x = x.dropna(subset=["forest_loss_pct"])
        return np.polyfit(x.year, x.forest_loss_pct, 1)[0] if len(x) > 2 else np.nan

    f = pd.DataFrame({
        "forest_2001": p[p.year == 2001].set_index("basin_id").forest_cover_pct,
        # Loss is heavily right-skewed (most treated basins lose ~0); asinh keeps the
        # Mahalanobis distance from being driven by a few high-loss basins.
        "loss_rate_base": np.arcsinh(g.forest_loss_pct.mean() * 10),
        "loss_trend_base": np.arcsinh(g.apply(slope, include_groups=False) * 100),
        "ndvi_base": g.ndvi.mean(),
        "et_base": g.et_mm.mean(),
        "precip_base": g.precip_mm.mean(),
    })
    return f.join(static[["elevation_m", "slope_deg"]])


def match(feats: pd.DataFrame, regions: pd.Series, treated: list, controls: list, k: int = 5) -> tuple[list, pd.DataFrame]:
    """k nearest never-enrolled basins per treated basin, exact on region (with replacement)."""
    chosen, rows = [], []
    for reg in sorted(regions.loc[treated].unique()):
        tr = [t for t in treated if regions[t] == reg and feats.loc[t, MATCH_ON].notna().all()]
        cand = [c for c in controls if regions.get(c) == reg and feats.loc[c, MATCH_ON].notna().all()]
        if not tr or len(cand) < k:
            continue
        sel = select_donors(feats.loc[tr + cand], tr, covariates=MATCH_ON, candidates=cand,
                            caliper=np.inf, hard_caliper=np.inf)
        for t in tr:
            nn = sel.nearest(k, t)
            chosen += nn
            rows.append({"treated": t, "region": reg, "nearest": nn})
    return sorted(set(chosen)), pd.DataFrame(rows)


def balance(feats: pd.DataFrame, treated: list, controls: list) -> pd.DataFrame:
    t, c = feats.loc[treated, MATCH_ON], feats.loc[controls, MATCH_ON]
    sd = np.sqrt((t.var() + c.var()) / 2).replace(0, np.nan)
    return pd.DataFrame({"treated_mean": t.mean(), "control_mean": c.mean(), "smd": (t.mean() - c.mean()) / sd})


def estimate(p: pd.DataFrame, outcome: str, treat: dict, controls: list, alpha: float, n_boot: int, seed: int = 0):
    units = list(treat) + list(controls)
    sub = p[p.basin_id.isin(units)][["basin_id", "year", outcome]]
    complete = sub.groupby("basin_id")[outcome].apply(lambda x: x.notna().all() and len(x) == p.year.nunique())
    keep = complete[complete].index
    sub = sub[sub.basin_id.isin(keep)]
    tr = {u: s for u, s in treat.items() if u in keep}
    r = staggered_att(sub, outcome, tr, unit="basin_id", time="year", alpha=alpha, n_boot=n_boot, seed=seed)
    pre = r.dynamic[(r.dynamic.index < -1) & (r.dynamic.index >= -6)]
    base_mean = sub[(sub.basin_id.isin(tr)) & (sub.year < min(tr.values()))][outcome].mean()
    return r, {"outcome": outcome, "att": r.att, "ci_low": r.ci[0], "ci_high": r.ci[1], "se": r.se,
               "n_treated": len(tr), "n_controls": int(sum(u not in tr for u in keep)),
               "pre_mean_att": pre.att.mean(), "pre_max_abs_t": float((pre.att / pre.se).abs().max()),
               "treated_baseline_mean": base_mean}


def placebo_validation(p: pd.DataFrame, cohorts: list, pool: list, outcomes, reps: int, alpha: float,
                       n_boot: int, seed: int = 0) -> pd.DataFrame:
    """Random enrollment among never-enrolled basins, with the real cohort years.

    No basin in the pool is treated, so every rejection is a false positive and the
    interval half-width is the design's resolution for each outcome.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(reps):
        fake = rng.choice(pool, len(cohorts), replace=False)
        tr = dict(zip(fake, rng.permutation(cohorts)))
        rest = [c for c in pool if c not in tr]
        for o in outcomes:
            _, row = estimate(p, o, tr, rest, alpha, n_boot, seed=rep)
            rows.append({"rep": rep, **row})
    r = pd.DataFrame(rows)
    r["reject"] = (r.ci_low > 0) | (r.ci_high < 0)
    r["half_width"] = (r.ci_high - r.ci_low) / 2
    return r.groupby("outcome").agg(false_positive_rate=("reject", "mean"), median_half_width=("half_width", "median"),
                                    sd_estimate=("att", "std"), reps=("rep", "nunique")).reset_index()


def plot_event(results: dict, path: Path, window=(-6, 10), title: str = ""):
    plt = _plt()
    keys = list(results)
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 3.2 * nrow), squeeze=False)
    for ax, key in zip(axes.ravel(), keys):
        dyn = results[key].dynamic
        d = dyn[(dyn.index >= window[0]) & (dyn.index <= window[1])]
        pre = d.index < 0
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.axvline(-0.5, color=GRID, lw=6)
        ax.fill_between(d.index, d.ci_low, d.ci_high, color=PLACEBO, alpha=0.5, lw=0)
        ax.plot(d.index[pre], d.att[pre], "o-", color=MUTED, ms=3.5)
        ax.plot(d.index[~pre], d.att[~pre], "o-", color=ORANGE if key == "precip_mm" else BLUE, ms=3.5)
        label, unit = OUTCOMES[key]
        ax.set_title(label)
        ax.set_ylabel(unit)
        ax.set_xlabel("years since enrollment")
    for ax in axes.ravel()[len(keys):]:
        ax.set_visible(False)
    fig.suptitle(title or "Socio Bosque: event-study ATT, matched never-enrolled controls (90% bootstrap band)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_summary(summary: pd.DataFrame, path: Path):
    plt = _plt()
    specs = ["unmatched", "matched", "matched, Amazon", "dose 20-50%", "dose >=50%"]
    keys = [k for k in OUTCOMES if k in set(summary.outcome)]
    fig, axes = plt.subplots(1, len(keys), figsize=(2.3 * len(keys), 3.6), sharey=True)
    for ax, key in zip(np.atleast_1d(axes), keys):
        s = summary[summary.outcome == key].set_index("spec").reindex(specs)
        y = np.arange(len(specs))[::-1]
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.hlines(y, s.ci_low, s.ci_high, color=BLUE, lw=2)
        ax.plot(s.att, y, "o", color=BLUE, ms=5)
        ax.set_title(OUTCOMES[key][0], fontsize=9)
        ax.set_xlabel(OUTCOMES[key][1], fontsize=8)
    np.atleast_1d(axes)[0].set_yticks(np.arange(len(specs))[::-1], specs)
    fig.suptitle("Socio Bosque: overall ATT by specification (90% bootstrap CI)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/socio_bosque")
    ap.add_argument("--out", default="docs/results")
    ap.add_argument("--fig-dir", default="docs/figures")
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--n-boot", type=int, default=499)
    ap.add_argument("--k", type=int, default=1, help="matched controls per treated basin (main spec)")
    ap.add_argument("--placebo-reps", type=int, default=200)
    a = ap.parse_args(argv)
    d, out, figdir = Path(a.data), Path(a.out), Path(a.fig_dir)
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    p, static, regions, treat, controls, expo = load(d)
    region = regions.region
    feats = baseline_features(p, static)
    treated = [t for t in treat if t in feats.index]
    matched, pairs = match(feats, region, treated, controls, k=a.k)
    treated = sorted(set(pairs.treated))
    treat = {t: treat[t] for t in treated}
    all_ctrl = [c for c in controls if c in feats.index]

    bal = pd.concat({"unmatched": balance(feats, treated, all_ctrl), "matched": balance(feats, treated, matched)},
                    axis=1)
    bal.to_csv(out / "socio_bosque_balance.csv")
    print(f"treated {len(treated)} | controls: all {len(all_ctrl)}, matched {len(matched)}")
    print(bal.round(3).to_string())

    peak = expo.reindex(treated).max(axis=1)
    amazon = [t for t in treated if region[t] == "Amazon"]
    amazon_ctrl = [c for c in matched if region[c] == "Amazon"]
    specs = {
        "unmatched": (treat, all_ctrl),
        "matched": (treat, matched),
        "matched, Amazon": ({t: treat[t] for t in amazon}, amazon_ctrl),
        "dose 20-50%": ({t: treat[t] for t in treated if peak[t] < 0.5}, matched),
        "dose >=50%": ({t: treat[t] for t in treated if peak[t] >= 0.5}, matched),
    }
    rows, dyn_main, dyn_amazon, dyn_rows = [], {}, {}, []
    for spec, (tr, ct) in specs.items():
        for outcome in OUTCOMES:
            r, row = estimate(p, outcome, tr, ct, a.alpha, a.n_boot)
            rows.append({"spec": spec, **row})
            dyn_rows.append(r.dynamic.assign(spec=spec, outcome=outcome).reset_index())
            if spec == "matched":
                dyn_main[outcome] = r
            if spec == "matched, Amazon":
                dyn_amazon[outcome] = r
        print(f"done: {spec}")
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "socio_bosque_att.csv", index=False)
    pd.concat(dyn_rows).to_csv(out / "socio_bosque_event_study.csv", index=False)
    plot_event(dyn_amazon, figdir / "socio_bosque_event_study.png",
               title=f"Socio Bosque, Amazon: event-study ATT vs matched never-enrolled basins (k={a.k}, 90% band)")
    plot_event(dyn_main, figdir / "socio_bosque_event_study_all_regions.png",
               title=f"Socio Bosque, all regions: event-study ATT vs matched never-enrolled basins (k={a.k}, 90% band)")

    # Sensitivity to the number of matched controls per treated basin.
    sens = []
    for k in (1, 2, 3, 5):
        m_k, pairs_k = match(feats, region, treated, controls, k=k)
        b_k = balance(feats, treated, m_k)
        for spec, (tr, ct) in {"matched": (treat, m_k),
                               "matched, Amazon": ({t: treat[t] for t in amazon}, [c for c in m_k if region[c] == "Amazon"])}.items():
            for outcome in OUTCOMES:
                _, row = estimate(p, outcome, tr, ct, a.alpha, a.n_boot)
                sens.append({"k": k, "max_abs_smd": b_k.smd.abs().max(), "spec": spec, **row})
    pd.DataFrame(sens).to_csv(out / "socio_bosque_sensitivity_k.csv", index=False)

    # Placebo validation among never-enrolled Amazon basins.
    pool = [c for c in all_ctrl if region[c] == "Amazon"]
    val = placebo_validation(p, [treat[t] for t in amazon], pool, OUTCOMES, a.placebo_reps, a.alpha, 99)
    val.to_csv(out / "socio_bosque_placebo_validation.csv", index=False)
    print(val.round(4).to_string(index=False))
    plot_summary(summary, figdir / "socio_bosque_att.png")
    cols = ["spec", "outcome", "att", "ci_low", "ci_high", "n_treated", "n_controls", "pre_mean_att", "pre_max_abs_t",
            "treated_baseline_mean"]
    print(summary[cols].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
