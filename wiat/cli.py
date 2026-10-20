"""Command-line interface.

    wiat simulate   --out data/sim            # synthetic study with known truth
    wiat extract    --config study.toml       # Earth Engine panel for real basins
    wiat analyze    --panel ... --treatments ... --outcome et_mm --out reports/x
    wiat scarcity   --reps 200 --out docs     # Monte Carlo: confidence vs donor scarcity
"""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from pathlib import Path

import pandas as pd


def _cmd_simulate(a) -> int:
    from .simulate import simulate_watersheds

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    s = simulate_watersheds(n_basins=a.basins, n_years=a.years, treated_extremity=a.extremity,
                            ndvi_effect=a.effect, seed=a.seed)
    s.panel.to_csv(out / "panel.csv", index=False)
    s.static.to_csv(out / "static.csv")
    pd.DataFrame({"basin_id": list(s.treat_start), "start": [f"{v}-01-01" for v in s.treat_start.values()]}).to_csv(
        out / "treatments.csv", index=False)
    truth = {o: s.true_att(o) for o in ("et_mm", "water_yield_mm", "evaporative_index", "ndvi")}
    (out / "truth.json").write_text(json.dumps(truth, indent=2))
    print(f"wrote {out}/panel.csv, static.csv, treatments.csv, truth.json")
    print("true annual ATT:", {k: round(v, 4) for k, v in truth.items()})
    return 0


def _cmd_analyze(a) -> int:
    from .attribution import AttributionStudy
    from .report import write_report

    panel = pd.read_csv(a.panel)
    if "date" in panel and not pd.api.types.is_numeric_dtype(panel["date"]):
        panel["date"] = pd.to_datetime(panel["date"])
    static = pd.read_csv(a.static, index_col=0) if a.static else None
    tr = pd.read_csv(a.treatments)
    panel["basin_id"] = panel["basin_id"].astype(str)
    tr["basin_id"] = tr["basin_id"].astype(str)
    if static is not None:
        static.index = static.index.astype(str)
    excl = {}
    if a.exclusions:
        ex = pd.read_csv(a.exclusions)
        excl = dict(zip(ex.basin_id.astype(str), ex.reason))
    study = AttributionStudy(
        panel, a.outcome, dict(zip(tr.basin_id, tr.start)), static=static,
        freq=None if a.freq == "none" else a.freq, water_year_start_month=a.water_year_start,
        did_covariates=a.covariates, caliper=a.caliper, exclude=excl, alpha=a.alpha,
        scm_method=a.scm_method, min_donors=a.min_donors,
    )
    rep = study.run()
    path = write_report(rep, a.out)
    print(f"verdict: {rep.verdict}")
    for f in rep.flags:
        print(f"  {f}")
    print(f"report: {path}")
    return 0


def _cmd_scarcity(a) -> int:
    from .report import _plt, plot_estimator_comparison, plot_scarcity
    from .scarcity import run_scarcity_experiment, summarize_scarcity

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.from_csv:
        raw = pd.read_csv(a.from_csv)
    else:
        raw = run_scarcity_experiment(
            donor_counts=a.donors, extremities=a.extremities, n_reps=a.reps, outcome=a.outcome,
            alpha=a.alpha, n_jobs=a.jobs, seed=a.seed,
        )
        raw.to_csv(out / f"scarcity_raw_{a.outcome}.csv.gz", index=False)
    summ = summarize_scarcity(raw, alpha=a.alpha)
    summ.to_csv(out / f"scarcity_summary_{a.outcome}.csv", index=False)
    plt = _plt()
    figdir = Path(a.fig_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    for est in summ.estimator.unique():
        fig = plot_scarcity(summ, est, alpha=a.alpha)
        fig.savefig(figdir / f"scarcity_{a.outcome}_{est}.png", bbox_inches="tight")
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    plot_estimator_comparison(summ, "size", 0.0, ylabel="false-positive rate (no effect)", ref=a.alpha, ax=axes[0])
    plot_estimator_comparison(summ, "coverage", 1.0, ylabel="interval coverage", ref=1 - a.alpha, ax=axes[1])
    plot_estimator_comparison(summ, "confident_correct", 1.0, ylabel="P(confident, correct)", ax=axes[2])
    axes[0].set_title("Typical basin: size")
    axes[1].set_title("Atypical basin (1 SD): coverage")
    axes[2].set_title("Atypical basin (1 SD): usable result")
    fig.tight_layout()
    fig.savefig(figdir / f"scarcity_{a.outcome}_estimators.png", bbox_inches="tight")
    plt.close(fig)
    cols = ["estimator", "extremity", "n_donors", "bias", "rel_rmse", "coverage", "ci_bounded",
            "power", "size", "confident_correct"]
    print(summ[cols].round(3).to_string(index=False))
    return 0


def _cmd_extract(a) -> int:
    from .data import gee
    from .matching import hydrologically_connected

    cfg = tomllib.loads(Path(a.config).read_text())
    gee.initialize(project=cfg.get("project") or a.project)
    out_root = Path(a.out)
    for st in cfg["studies"]:
        name = st["name"]
        out = out_root / name
        out.mkdir(parents=True, exist_ok=True)
        level = st.get("hydrobasins_level", 8)
        fc = gee.basins(level, region=st.get("region"), country=st.get("country"))
        table = gee.basin_table(fc)
        treat_rows, excl = [], {}
        for iv in st["interventions"]:
            b = gee.basin_at(iv["lon"], iv["lat"], level)
            bid = b.get("HYBAS_ID").getInfo()
            treat_rows.append({"basin_id": bid, "start": iv["start"], "name": iv.get("name", "")})
            for c in hydrologically_connected(table, bid):
                excl[c] = f"hydrologically connected to {bid}"
        for bid in st.get("exclude", []):
            excl[bid] = "listed in config (other intervention / known disturbance)"
        spec = gee.PanelSpec(start=st.get("start", "2001-01-01"), end=st.get("end", "2025-01-01"),
                             freq=st.get("freq", "month"), scale=st.get("scale", 500),
                             sensors=tuple(st.get("sensors", ("L5", "L7", "L8", "L9", "S2"))))
        pd.DataFrame(treat_rows).to_csv(out / "treatments.csv", index=False)
        pd.DataFrame({"basin_id": list(excl), "reason": list(excl.values())}).to_csv(out / "exclusions.csv", index=False)
        gee.static_covariates(fc).to_csv(out / "static.csv")
        if a.export:
            task = gee.export_panel(fc, spec, description=f"wiat_{name}")
            print(f"{name}: started export task {task.id}. When it finishes, download the CSV and run\n"
                  f"  python -c \"from wiat.data.gee import load_export; load_export('<csv>').to_csv('{out}/panel.csv', index=False)\"")
        else:
            panel = gee.quality_mask(gee.fetch_panel(fc, spec))
            from .panel import add_hydrologic_outcomes

            add_hydrologic_outcomes(panel).to_csv(out / "panel.csv", index=False)
            print(f"{name}: wrote {out}/panel.csv ({panel.basin_id.nunique()} basins)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="wiat", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="write a synthetic study with known effect")
    s.add_argument("--out", default="data/sim")
    s.add_argument("--basins", type=int, default=40)
    s.add_argument("--years", type=int, default=24)
    s.add_argument("--extremity", type=float, default=0.0)
    s.add_argument("--effect", type=float, default=0.08, help="full-ramp NDVI increase")
    s.add_argument("--seed", type=int, default=4)
    s.set_defaults(func=_cmd_simulate)

    an = sub.add_parser("analyze", help="run an attribution study on a panel")
    an.add_argument("--panel", required=True)
    an.add_argument("--treatments", required=True, help="CSV with basin_id,start")
    an.add_argument("--static")
    an.add_argument("--exclusions", help="CSV with basin_id,reason")
    an.add_argument("--outcome", default="et_mm")
    an.add_argument("--freq", default="year", choices=["year", "water_year", "none"])
    an.add_argument("--water-year-start", type=int, default=1)
    an.add_argument("--covariates", nargs="*", default=["precip_mm"])
    an.add_argument("--caliper", type=float, default=1.0)
    an.add_argument("--min-donors", type=int, default=5)
    an.add_argument("--alpha", type=float, default=0.1)
    an.add_argument("--scm-method", default="demeaned", choices=["standard", "demeaned", "augmented"])
    an.add_argument("--out", default="reports/study")
    an.set_defaults(func=_cmd_analyze)

    sc = sub.add_parser("scarcity", help="Monte Carlo of attribution quality vs donor scarcity")
    sc.add_argument("--donors", type=int, nargs="+", default=[2, 4, 6, 9, 14, 19, 29, 49])
    sc.add_argument("--extremities", type=float, nargs="+", default=[0.0, 1.0, 2.0])
    sc.add_argument("--reps", type=int, default=200)
    sc.add_argument("--outcome", default="et_mm")
    sc.add_argument("--alpha", type=float, default=0.1)
    sc.add_argument("--jobs", type=int, default=None)
    sc.add_argument("--seed", type=int, default=20260101)
    sc.add_argument("--from-csv", help="re-summarise an existing raw results CSV")
    sc.add_argument("--out", default="docs/results")
    sc.add_argument("--fig-dir", default="docs/figures")
    sc.set_defaults(func=_cmd_scarcity)

    ex = sub.add_parser("extract", help="extract a basin panel from Google Earth Engine")
    ex.add_argument("--config", required=True)
    ex.add_argument("--project", help="Earth Engine cloud project")
    ex.add_argument("--export", action="store_true", help="batch export to Drive instead of getInfo")
    ex.add_argument("--out", default="data")
    ex.set_defaults(func=_cmd_extract)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
