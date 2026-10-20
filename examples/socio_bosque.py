"""Socio Bosque (Ecuador) staggered study: data preparation and extraction.

See docs/socio_bosque_plan.md. Needs the ``gee`` and ``geo`` extras.

    python examples/socio_bosque.py prepare --project <ee-project>
    python examples/socio_bosque.py extract --project <ee-project> --pilot 20    # interactive, small
    python examples/socio_bosque.py extract --project <ee-project> --export      # full, to Google Drive
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from wiat.data import gee
from wiat.data.interventions import basin_exposure, control_pool, socio_bosque, summarize_exposure, treatment_years

LEVEL = 12
COUNTRY = "Ecuador"


def _basins_gdf(out: Path):
    import geopandas as gpd

    path = out / f"hybas{LEVEL}_ecuador.geojson"
    if not path.exists():
        fc = gee.basins(LEVEL, country=COUNTRY).select(["HYBAS_ID", "NEXT_DOWN", "SUB_AREA"])
        urllib.request.urlretrieve(fc.getDownloadURL(filetype="geojson"), path)
    return gpd.read_file(path)


def prepare(a) -> int:
    import geopandas as gpd

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    gee.initialize(a.project)
    agr_path = out / "agreements.gpkg"
    if agr_path.exists():
        agr = gpd.read_file(agr_path)
    else:
        agr = socio_bosque()
        agr.to_file(agr_path)
    print(f"agreements: {len(agr)}, {agr.to_crs(32717).area.sum() / 1e6:,.0f} km2, signed {agr.year.min()}-{agr.year.max()}")

    basins = _basins_gdf(out)
    basins["HYBAS_ID"] = basins.HYBAS_ID.astype(np.int64)
    basins[["HYBAS_ID", "NEXT_DOWN", "SUB_AREA"]].to_csv(out / "basins.csv", index=False)
    share = basin_exposure(agr, basins)
    share.to_csv(out / "exposure.csv")
    print(f"basins: {len(basins)}; touched by an agreement: {len(share)}")
    print(summarize_exposure(share).to_string())

    # Galápagos (inside Ecuador's boundary; no Socio Bosque; CHIRPS unreliable over small
    # islands) and tiny coastal fragments are outside the study population.
    lon = basins.to_crs(32717).centroid.to_crs(4326).x.to_numpy()
    keep = (lon > -84) & (basins.SUB_AREA.to_numpy() >= a.min_area)
    print(f"study basins: {int(keep.sum())} (dropped {int((lon <= -84).sum())} Galapagos, "
          f"{int(((lon > -84) & (basins.SUB_AREA < a.min_area)).sum())} under {a.min_area:g} km2)")
    basins = basins[keep]
    basins[["HYBAS_ID"]].rename(columns={"HYBAS_ID": "basin_id"}).to_csv(out / "study_basins.csv", index=False)
    share = share[share.index.isin(basins.HYBAS_ID)]
    treat = treatment_years(share, threshold=a.threshold, max_year=a.max_cohort)
    pd.DataFrame({"basin_id": list(treat), "start": list(treat.values())}).to_csv(out / "treatments.csv", index=False)
    ctrl = control_pool(share, basins.HYBAS_ID.tolist(), max_share=a.control_max_share)
    pd.DataFrame({"basin_id": ctrl}).to_csv(out / "controls.csv", index=False)
    print(f"treated (share >= {a.threshold:.0%}, first treated year <= {a.max_cohort}): {len(treat)}; "
          f"never-enrolled controls (share <= {a.control_max_share:.0%}): {len(ctrl)}")

    if not (out / "static.csv").exists():
        gee.static_covariates(gee.basins(LEVEL, country=COUNTRY)).to_csv(out / "static.csv")
    if not (out / "forest.csv").exists():
        gee.forest_loss_panel(gee.basins(LEVEL, country=COUNTRY), 2001, a.end_year - 1).to_csv(
            out / "forest.csv", index=False)
    print(f"wrote {out}/: exposure, treatments, controls, static, forest")
    return 0


def extract(a) -> int:
    out = Path(a.out)
    gee.initialize(a.project)
    treat = pd.read_csv(out / "treatments.csv")
    ctrl = pd.read_csv(out / "controls.csv")
    if a.pilot:
        rng = np.random.default_rng(0)
        ids = list(rng.choice(treat.basin_id, min(a.pilot, len(treat)), replace=False)) + \
            list(rng.choice(ctrl.basin_id, min(a.pilot, len(ctrl)), replace=False))
        fc = gee.basins_by_id(ids, LEVEL)
        spec = gee.PanelSpec(start=a.pilot_start, end=f"{a.end_year}-01-01", freq="year", scale=a.scale)
        panel = gee.quality_mask(gee.fetch_panel(fc, spec, verbose=True), min_clear=a.min_clear)
        from wiat.panel import add_hydrologic_outcomes

        add_hydrologic_outcomes(panel).to_csv(out / "panel_pilot.csv", index=False)
        print(f"pilot: {panel.basin_id.nunique()} basins x {panel.date.nunique()} years -> {out}/panel_pilot.csv")
        return 0
    if not a.export:
        raise SystemExit("pass --pilot N for a small interactive run, or --export for the full panel")
    fc = gee.basins_by_id(pd.read_csv(out / "study_basins.csv").basin_id, LEVEL)
    for y0 in range(a.start_year, a.end_year, a.block):
        y1 = min(y0 + a.block, a.end_year)
        spec = gee.PanelSpec(start=f"{y0}-01-01", end=f"{y1}-01-01", freq="year", scale=a.scale)
        task = gee.export_panel(fc, spec, description=f"wiat_socio_bosque_{y0}_{y1 - 1}", folder=a.folder)
        print(f"started {task.id}: {y0}-{y1 - 1}")
    print(f"When the tasks finish, download the CSVs from Drive folder '{a.folder}' into {out}/export/ and run\n"
          f"  python examples/socio_bosque.py assemble")
    return 0


def assemble(a) -> int:
    out = Path(a.out)
    parts = sorted((out / "export").glob("*.csv"))
    if not parts:
        raise SystemExit(f"no CSVs in {out}/export/")
    panel = pd.concat([gee.load_export(p, min_clear=a.min_clear) for p in parts], ignore_index=True)
    panel = panel.drop_duplicates(["basin_id", "date"]).sort_values(["basin_id", "date"])
    panel.to_csv(out / "panel.csv", index=False)
    print(f"panel: {panel.basin_id.nunique()} basins x {panel.date.nunique()} years -> {out}/panel.csv")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["prepare", "extract", "assemble"])
    ap.add_argument("--project", help="Earth Engine cloud project")
    ap.add_argument("--out", default="data/socio_bosque")
    ap.add_argument("--threshold", type=float, default=0.2, help="basin share enrolled that counts as treated")
    ap.add_argument("--control-max-share", type=float, default=0.01)
    ap.add_argument("--min-area", type=float, default=10, help="km2; smaller basins are dropped")
    ap.add_argument("--max-cohort", type=int, default=2016, help="last first-treated year kept")
    ap.add_argument("--start-year", type=int, default=2001)
    ap.add_argument("--end-year", type=int, default=2025, help="exclusive")
    ap.add_argument("--block", type=int, default=4, help="years per export task")
    ap.add_argument("--scale", type=float, default=250)
    ap.add_argument("--min-clear", type=float, default=0.3)
    ap.add_argument("--pilot", type=int, default=0, help="treated and control basins each, interactive")
    ap.add_argument("--pilot-start", default="2019-01-01")
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--folder", default="wiat_socio_bosque")
    a = ap.parse_args(argv)
    return {"prepare": prepare, "extract": extract, "assemble": assemble}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
