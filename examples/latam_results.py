"""Results on the real Latin American panels (docs/latam_results.md).

Needs the panels from ``wiat extract --config examples/latam_water_funds.toml``:

    python examples/latam_results.py --data data --out docs/results --fig-dir docs/figures

Three analyses per study:

1. the attribution study itself, for ET, the evaporative index, NDVI and EVI;
2. donor thinning: the same study re-estimated with the k best-matched donors and
   with random k-subsets;
3. injected-effect validation: every eligible untreated basin plays the treated
   basin, with a known effect added after the real start year
   (:mod:`wiat.validation`).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from wiat import AttributionStudy, aggregate, donor_thinning
from wiat.report import _plt, plot_donor_thinning, plot_estimator_comparison, plot_scarcity
from wiat.validation import ATYPICALITY_LABELS, injection_experiment, prepare_region, summarize_injection

STUDIES = {"extrema": "Extrema (MG, Brazil)", "camboriu": "Camboriú (SC, Brazil)",
           "quito_paramo": "Quito páramo (Ecuador)"}
OUTCOMES = ("et_mm", "evaporative_index", "ndvi", "evi")
EFFECTS = (0.025, 0.05, 0.10)
DONOR_COUNTS = (4, 6, 9, 14, 19, 29, 49)


def load(root: Path, name: str):
    d = root / name
    p = pd.read_csv(d / "panel.csv", parse_dates=["date"])
    p["basin_id"] = p.basin_id.astype(str)
    s = pd.read_csv(d / "static.csv", index_col=0)
    s.index = s.index.astype(str)
    tr = pd.read_csv(d / "treatments.csv")
    ex = pd.read_csv(d / "exclusions.csv")
    return p, s, dict(zip(tr.basin_id.astype(str), tr.start)), dict(zip(ex.basin_id.astype(str), ex.reason))


def observed(panel, static, treat, excl, alpha):
    rows = []
    for outcome in OUTCOMES:
        for caliper in (1.0, 1.5):
            try:
                rep = AttributionStudy(panel, outcome, treat, static=static, exclude=excl, caliper=caliper,
                                       did_covariates=["precip_mm"], alpha=alpha).run()
            except ValueError as e:  # treated basin with gaps in this outcome
                rows.append({"outcome": outcome, "caliper": caliper, "verdict": f"not run: {e}"})
                continue
            base = rep.did.estimate if rep.did is not None else np.nan
            row = {"outcome": outcome, "caliper": caliper, "verdict": rep.verdict, "n_donors": len(rep.donors),
                   "flags": "; ".join(f.code for f in rep.flags if f.severity in ("fatal", "major"))}
            if rep.did is not None:
                row.update(did_estimate=base, did_ci_low=rep.did.ci[0], did_ci_high=rep.did.ci[1],
                           did_p=rep.did.p_value, did_detection_threshold=rep.did.detection_threshold)
            for s in rep.scm.values():
                row.update(scm_estimate=s.att, scm_ci_low=s.ci[0], scm_ci_high=s.ci[1], scm_p=s.p_value,
                           scm_pre_r2=s.diagnostics["pre_r2"])
            rows.append(row)
    return pd.DataFrame(rows)


def thinning(panel, static, treat, excl, alpha, outcome="et_mm"):
    study = AttributionStudy(panel, outcome, treat, static=static, exclude=excl)
    feats = study._features()
    cov = study._match_covariates(feats)
    from wiat.matching import select_donors
    from wiat.panel import balance

    annual, _ = balance(aggregate(panel, "year"), [outcome, "precip_mm"])
    t = next(iter(treat))
    pool = [u for u in annual.basin_id.unique() if u != t and u not in excl and u in feats.index]
    sel = select_donors(feats, [t], covariates=cov, candidates=pool, caliper=np.inf, hard_caliper=np.inf)
    ranked = sel.nearest(len(sel.candidates))
    start = study.treat_start_raw[t]
    from wiat.attribution import _to_period

    start = _to_period(start, "year", 1)
    sizes = [k for k in DONOR_COUNTS if k <= len(ranked)]
    near = donor_thinning(annual, outcome, t, start, ranked, sizes=sizes, mode="nearest",
                          covariates=["precip_mm"], alpha=alpha)
    rand = donor_thinning(annual, outcome, t, start, ranked, sizes=sizes, mode="random", n_draws=40,
                          covariates=["precip_mm"], alpha=alpha)
    return near.assign(mode="nearest"), rand.assign(mode="random")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="docs/results")
    ap.add_argument("--fig-dir", default="docs/figures")
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--jobs", type=int, default=None)
    a = ap.parse_args(argv)
    root, out, figdir = Path(a.data), Path(a.out), Path(a.fig_dir)
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    plt = _plt()

    obs, inj, thin = [], [], []
    for name, label in STUDIES.items():
        if not (root / name / "panel.csv").exists():
            print(f"{name}: no panel, skipped")
            continue
        panel, static, treat, excl = load(root, name)
        print(f"== {label}: {panel.basin_id.nunique()} basins")
        obs.append(observed(panel, static, treat, excl, a.alpha).assign(study=name))

        near, rand = thinning(panel, static, treat, excl, a.alpha)
        thin += [near.assign(study=name), rand.assign(study=name)]
        fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
        plot_donor_thinning(near, ax=axes[0])
        axes[0].set_title(f"{label}: k best-matched donors")
        plot_donor_thinning(rand, ax=axes[1])
        axes[1].set_title(f"{label}: random k donors (median of 40 draws)")
        for ax in axes:
            ax.set_ylabel("ET effect (mm/yr)")
        fig.tight_layout()
        fig.savefig(figdir / f"latam_thinning_{name}.png", bbox_inches="tight")
        plt.close(fig)

        for outcome in ("et_mm", "ndvi"):
            reg = prepare_region(panel, outcome, treat, static, excl)
            r = injection_experiment(reg, donor_counts=DONOR_COUNTS, effects=EFFECTS, alpha=a.alpha,
                                     n_jobs=a.jobs)
            inj.append(r.assign(study=name, outcome=outcome, n_pre=reg["start"] - int(reg["annual"].date.min())))
            print(f"   injection {outcome}: {len(reg['eligible'])} pseudo-treated basins, {len(r)} rows")

    obs = pd.concat(obs, ignore_index=True)
    obs.to_csv(out / "latam_observed.csv", index=False)
    thin = pd.concat(thin, ignore_index=True)
    thin.to_csv(out / "latam_thinning.csv", index=False)
    inj = pd.concat(inj, ignore_index=True)
    inj.to_csv(out / "latam_injection_raw.csv.gz", index=False)

    summ = []
    for (study, outcome), g in inj.groupby(["study", "outcome"]):
        for e in EFFECTS:
            summ.append(summarize_injection(g, e, a.alpha).assign(study=study, outcome=outcome))
        pooled = g.assign(extremity=-1.0)
        for e in EFFECTS:
            summ.append(summarize_injection(pooled, e, a.alpha).assign(study=study, outcome=outcome))
    summ = pd.concat(summ, ignore_index=True)
    summ.to_csv(out / "latam_injection_summary.csv", index=False)

    # Figures: degradation curves (5% effect on ET) by atypicality, one per study and estimator.
    for study in summ.study.unique():
        s = summ[(summ.study == study) & (summ.outcome == "et_mm") & np.isclose(summ.effect, 0.05)
                 & (summ.extremity >= 0)]
        for est in ("scm", "did_placebo"):
            fig = plot_scarcity(s, est, alpha=a.alpha, labels=ATYPICALITY_LABELS,
                                title=f"{STUDIES[study]}: real basins, injected 5% ET effect — {est}")
            fig.savefig(figdir / f"latam_injection_{study}_{est}.png", bbox_inches="tight")
            plt.close(fig)
    # Estimator comparison, all basins pooled.
    studies = list(summ.study.unique())
    fig, axes = plt.subplots(len(studies), 3, figsize=(13, 3.4 * len(studies)), squeeze=False)
    for row, study in zip(axes, studies):
        s = summ[(summ.study == study) & (summ.outcome == "et_mm") & np.isclose(summ.effect, 0.05)]
        plot_estimator_comparison(s, "size", -1.0, estimators=("did_cluster", "did_placebo", "scm"),
                                  ylabel="false-positive rate", ref=a.alpha, ax=row[0])
        plot_estimator_comparison(s, "coverage", -1.0, estimators=("did_cluster", "did_placebo", "scm"),
                                  ylabel="interval coverage", ref=1 - a.alpha, ax=row[1])
        plot_estimator_comparison(s, "confident_correct", -1.0, estimators=("did_cluster", "did_placebo", "scm"),
                                  ylabel="P(confident, correct)", ax=row[2])
        row[0].set_title(f"{STUDIES[study]}: no effect")
        row[1].set_title("injected 5% ET effect: coverage")
        row[2].set_title("injected 5% ET effect: usable result")
    fig.tight_layout()
    fig.savefig(figdir / "latam_estimators_et_mm.png", bbox_inches="tight")
    plt.close(fig)

    cols = ["study", "outcome", "caliper", "verdict", "n_donors", "did_estimate", "did_ci_low", "did_ci_high",
            "scm_estimate", "scm_ci_low", "scm_ci_high"]
    print(obs[[c for c in cols if c in obs]].round(4).to_string(index=False))
    show = ["study", "outcome", "estimator", "extremity", "n_donors", "effect", "size", "coverage", "power",
            "confident_correct", "median_rel_width"]
    print(summ[show].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
