"""Socio Bosque parcel-level first stage: enrolled forest pixels vs matched non-enrolled pixels.

The basin-level study (docs/socio_bosque_results.md) could not separate avoided
deforestation from selection, and basin averaging dilutes any effect inside
agreements. This script tests the first stage where the programme acts: on forest
pixels inside agreements, against matched forest pixels elsewhere.

    python examples/socio_bosque_parcels.py sample
    python examples/socio_bosque_parcels.py extract --project <ee-project>
    python examples/socio_bosque_parcels.py analyze

Design:
* Treated: random points inside agreements (60 m inside the boundary), sampled in
  proportion to area. Controls: random points in mainland Ecuador at least 1 km from
  any agreement. Protected areas (WDPA) and non-forest (tree cover 2000 < 30%) are
  dropped from both.
* Matching covariates: travel time to cities, distance to forest edge (2000), slope,
  elevation, tree cover 2000, and forest loss in 2001-04 within 1 km and 5 km (local
  deforestation pressure). Matched 1:1 with replacement, exact on region and
  enrollment cohort; a control must still be forest at the treated pixel's reference
  year.
* Outcome: annual forest-loss hazard (events per pixel-year at risk) after enrollment.
  Falsification: the same comparison in the 4 years *before* enrollment, among
  pixels forested 4 years earlier, matched only on 2001-04 information.
* Inference: bootstrap over agreements (treated clusters), carrying their matched pairs.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

COVARIATES = ["log_access", "log_edge", "slope", "elevation", "treecover", "loss_1km", "loss_5km"]
END_YEAR = 2024
PRE_YEARS = 4
COHORTS = range(2009, 2016)  # agreement signing years with a clean 2001-04 matching window


def sample(a) -> int:
    import geopandas as gpd

    d = Path(a.data)
    agr = gpd.read_file(d / "agreements.gpkg").to_crs(32717)
    agr = agr[agr.tpcnv != "MANGLAR"]  # mangrove concessions: different ecosystem and threat
    inner = agr.assign(geometry=agr.buffer(-60)).loc[lambda x: ~x.is_empty]
    rng = np.random.default_rng(a.seed)
    w = inner.area / inner.area.sum()
    counts = rng.multinomial(a.n_treated, w.to_numpy())
    pts = inner.geometry.sample_points(size=counts, rng=a.seed)
    tr = gpd.GeoDataFrame({"codcnv": inner.codcnv.to_numpy(), "year": inner.year.to_numpy(),
                           "tpcnv": inner.tpcnv.to_numpy(), "geometry": pts}, crs=32717).explode(index_parts=False)
    tr = tr[~tr.is_empty].assign(sample="area")
    # Area weighting leaves most small individual agreements without a point; add a fixed
    # number per individual agreement so they can be analysed on their own.
    ind = inner[inner.tpcnv == "INDIVIDUAL"]
    extra = ind.geometry.sample_points(size=a.per_individual, rng=a.seed + 2)
    extra = gpd.GeoDataFrame({"codcnv": ind.codcnv.to_numpy(), "year": ind.year.to_numpy(), "tpcnv": "INDIVIDUAL",
                              "geometry": extra}, crs=32717).explode(index_parts=False)
    tr = pd.concat([tr, extra[~extra.is_empty].assign(sample="individual")], ignore_index=True)

    basins = gpd.read_file(d / f"hybas12_ecuador.geojson").to_crs(32717)
    lon = basins.centroid.to_crs(4326).x
    mainland = basins[lon > -84].union_all()
    blocked = agr.buffer(1000).union_all()
    area = mainland.difference(blocked)
    ctl = gpd.GeoSeries([area], crs=32717).sample_points(size=a.n_control, rng=a.seed + 1).explode(index_parts=False)
    ctl = gpd.GeoDataFrame({"codcnv": None, "year": np.nan, "tpcnv": None, "sample": "control",
                            "geometry": ctl.to_numpy()}, crs=32717)

    pts = pd.concat([tr.assign(treated=1), ctl.assign(treated=0)], ignore_index=True)
    pts = gpd.GeoDataFrame(pts, geometry="geometry", crs=32717)
    regions = pd.read_csv(d / "regions.csv", index_col=0).region
    j = gpd.sjoin(pts, basins[["HYBAS_ID", "geometry"]], predicate="within", how="left")
    pts["basin_id"] = j.HYBAS_ID.to_numpy()
    pts["region"] = regions.reindex(pts.basin_id).to_numpy()
    pts = pts.dropna(subset=["region"]).to_crs(4326)
    pts["pid"] = np.arange(len(pts))
    pts["lon"], pts["lat"] = pts.geometry.x, pts.geometry.y
    pts.drop(columns="geometry").to_csv(d / "parcel_points.csv", index=False)
    print(f"points: {int(pts.treated.sum())} treated in {pts.codcnv.nunique()} agreements, "
          f"{int((pts.treated == 0).sum())} control candidates")
    print(pts.groupby(["region", "treated"]).size().unstack())
    return 0


def _image():
    import ee

    gfc = ee.Image("UMD/hansen/global_forest_change_2025_v1_13")
    cover = gfc.select("treecover2000")
    forest = cover.gte(30)
    lossyear = gfc.select("lossyear").unmask(0)  # masked where no loss
    early = lossyear.gte(1).And(lossyear.lte(4)).And(forest)
    loss_1km = early.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(1000, "meters"))
    coarse = early.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject("EPSG:4326", None, 300)
    loss_5km = coarse.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(5000, "meters"))
    edge = forest.Not().fastDistanceTransform(256).sqrt().multiply(30)
    dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation")
    access = ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0").select("accessibility")
    pa = ee.Image(0).paint(ee.FeatureCollection("WCMC/WDPA/current/polygons").filter(ee.Filter.eq("ISO3", "ECU")), 1)
    return (cover.rename("treecover").addBands(lossyear.rename("lossyear"))
            .addBands(loss_1km.rename("loss_1km")).addBands(loss_5km.rename("loss_5km"))
            .addBands(edge.rename("edge_m")).addBands(dem.rename("elevation"))
            .addBands(ee.Terrain.slope(dem).rename("slope")).addBands(access.rename("access_min"))
            .addBands(pa.rename("protected")))


def extract(a) -> int:
    """Run the extraction as batch tasks writing Earth Engine table assets, then download them.

    Interactive requests time out on the neighbourhood and distance bands, so each chunk
    of points is a batch task (no time limit). Re-running downloads finished tables and
    starts only missing tasks.
    """
    import time

    import ee

    from wiat.data import gee

    gee.initialize(a.project)
    d = Path(a.data)
    pts = pd.read_csv(d / "parcel_points.csv", low_memory=False)
    out = d / "parcel_extract"
    out.mkdir(exist_ok=True)
    root = f"projects/{a.project}/assets/wiat_parcels"
    try:
        ee.data.getAsset(root)
    except ee.EEException:
        ee.data.createFolder(root)
    img = _image()
    bands = ["treecover", "lossyear", "loss_1km", "loss_5km", "edge_m", "elevation", "slope", "access_min", "protected"]
    chunks = list(range(0, len(pts), a.chunk))
    active = {t["description"] for t in ee.data.getTaskList() if t["state"] in ("READY", "RUNNING")}
    for i in chunks:
        name = f"chunk_{i // a.chunk:03d}"
        if (out / f"{name}.csv").exists():
            continue
        try:
            ee.data.getAsset(f"{root}/{name}")
            continue
        except ee.EEException:
            pass
        if f"wiat_parcels_{name}" in active:
            continue
        c = pts.iloc[i:i + a.chunk]
        fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([x, y]), {"pid": int(p)})
                                   for x, y, p in zip(c.lon, c.lat, c.pid)])
        # Table assets need geometries, so the points are kept (and dropped again on download).
        s = img.sampleRegions(collection=fc, properties=["pid"], scale=30, tileScale=8, geometries=True)
        task = ee.batch.Export.table.toAsset(collection=s.select(["pid", *bands]), description=f"wiat_parcels_{name}",
                                             assetId=f"{root}/{name}")
        task.start()
        print(f"started {name}: {len(c)} points")
    while True:
        pending = []
        for i in chunks:
            name = f"chunk_{i // a.chunk:03d}"
            path = out / f"{name}.csv"
            if path.exists():
                continue
            try:
                ee.data.getAsset(f"{root}/{name}")
            except ee.EEException:
                pending.append(name)
                continue
            url = ee.FeatureCollection(f"{root}/{name}").getDownloadURL(filetype="csv", selectors=["pid", *bands])
            urllib.request.urlretrieve(url, path)
            print(f"downloaded {name}")
        if not pending or not a.wait:
            break
        print(f"waiting for {len(pending)} task(s)...")
        time.sleep(60)
    print(f"{len(list(out.glob('chunk_*.csv')))} of {len(chunks)} chunks downloaded")
    return 0


def _prepare(d: Path) -> pd.DataFrame:
    pts = pd.read_csv(d / "parcel_points.csv", low_memory=False)
    x = pd.concat([pd.read_csv(p) for p in sorted((d / "parcel_extract").glob("chunk_*.csv"))])
    df = pts.merge(x, on="pid")
    df = df[(df.treecover >= 30) & (df.protected == 0)].copy()
    df["loss_year"] = np.where(df.lossyear > 0, df.lossyear + 2000, np.inf)
    df["log_access"] = np.log1p(df.access_min)
    # fastDistanceTransform searches 256 pixels (~10.9 km diagonal); beyond that it returns a
    # fixed sentinel (~1,390 km). Cap so "deep interior" is one value, not a fake distance.
    df["log_edge"] = np.log1p(df.edge_m.clip(upper=10_000))
    return df.dropna(subset=COVARIATES)


def _match(tr: pd.DataFrame, ctl: pd.DataFrame) -> np.ndarray:
    """Nearest control (Mahalanobis on COVARIATES) for each treated row; returns positions in ctl."""
    from scipy.spatial import cKDTree

    pool = pd.concat([tr[COVARIATES], ctl[COVARIATES]])
    mu, cov = pool.mean().to_numpy(), np.cov(pool.to_numpy().T)
    w = np.linalg.cholesky(np.linalg.pinv(cov))
    _, idx = cKDTree((ctl[COVARIATES].to_numpy() - mu) @ w).query((tr[COVARIATES].to_numpy() - mu) @ w, k=1)
    return idx


def matched_pairs(df: pd.DataFrame, window: str) -> pd.DataFrame:
    """Treated/control pairs and their outcome windows.

    ``post``: at risk from the signing year + 1 to END_YEAR, if forested at the end of the
    signing year. ``pre``: at risk in the PRE_YEARS before signing, if forested
    PRE_YEARS + 1 years before it.
    """
    rows = []
    tr_all, ctl_all = df[df.treated == 1], df[df.treated == 0]
    for reg in sorted(tr_all.region.unique()):
        for e in COHORTS:
            t0, t1 = (e + 1, END_YEAR) if window == "post" else (e - PRE_YEARS + 1, e)
            tr = tr_all[(tr_all.region == reg) & (tr_all.year == e) & (tr_all.loss_year >= t0)]
            ctl = ctl_all[(ctl_all.region == reg) & (ctl_all.loss_year >= t0)]
            if len(tr) == 0 or len(ctl) < 50:
                continue
            idx = _match(tr, ctl)
            c = ctl.iloc[idx]
            for side, g in (("treated", tr), ("control", c)):
                lost = (g.loss_year >= t0) & (g.loss_year <= t1)
                years = np.minimum(g.loss_year, t1) - t0 + 1  # pixel-years at risk
                rows.append(pd.DataFrame({"side": side, "codcnv": tr.codcnv.to_numpy(), "region": reg, "cohort": e,
                                          "tpcnv": tr.tpcnv.to_numpy(), "sample": tr["sample"].to_numpy(),
                                          "lost": lost.to_numpy().astype(int), "years": years.to_numpy(),
                                          **{c_: g[c_].to_numpy() for c_ in COVARIATES}}))
    return pd.concat(rows, ignore_index=True)


def hazard_diff(pairs: pd.DataFrame, n_boot: int, alpha: float, seed: int = 0) -> dict:
    """Treated minus control loss hazard (%/yr), bootstrapping agreements."""
    agg = pairs.groupby(["codcnv", "side"])[["lost", "years"]].sum().unstack("side").fillna(0)

    def est(a):
        t = a[("lost", "treated")].sum() / a[("years", "treated")].sum()
        c = a[("lost", "control")].sum() / a[("years", "control")].sum()
        return 100 * t, 100 * c

    t, c = est(agg)
    rng = np.random.default_rng(seed)
    boots = np.array([np.subtract(*est(agg.iloc[rng.integers(0, len(agg), len(agg))])) for _ in range(n_boot)])
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return {"treated_hazard": t, "control_hazard": c, "diff": t - c, "ci_low": lo, "ci_high": hi,
            "ratio": t / c if c > 0 else np.nan, "n_pairs": int((pairs.side == "treated").sum()),
            "n_agreements": len(agg)}


def change_in_gap(post: pd.DataFrame, pre: pd.DataFrame, n_boot: int, alpha: float, seed: int = 0) -> dict:
    """(post treated-control hazard gap) - (pre gap), %/yr, bootstrapping agreements jointly.

    The pre-enrollment gap measures how differently matched treated and control pixels
    were already behaving; subtracting it removes time-invariant differences in
    threat that the covariates miss (a difference-in-differences on hazards).
    """
    def table(pairs):
        return pairs.groupby(["codcnv", "side"])[["lost", "years"]].sum().unstack("side").fillna(0)

    tp, tq = table(post), table(pre)
    ids = np.array(sorted(set(tp.index) | set(tq.index)))
    tp, tq = tp.reindex(ids, fill_value=0), tq.reindex(ids, fill_value=0)

    def gap(t):
        return 100 * (t[("lost", "treated")].sum() / t[("years", "treated")].sum()
                      - t[("lost", "control")].sum() / t[("years", "control")].sum())

    def ratio(t):
        return (t[("lost", "treated")].sum() / t[("years", "treated")].sum()) / max(
            t[("lost", "control")].sum() / t[("years", "control")].sum(), 1e-12)

    est = gap(tp) - gap(tq)
    # Same comparison on the ratio scale: proportional changes in clearing cancel out.
    rr = ratio(tp) / ratio(tq) if ratio(tq) > 0 else np.nan
    rng = np.random.default_rng(seed)
    boots, rboots = [], []
    for _ in range(n_boot):
        i = rng.integers(0, len(ids), len(ids))
        boots.append(gap(tp.iloc[i]) - gap(tq.iloc[i]))
        q = ratio(tq.iloc[i])
        rboots.append(ratio(tp.iloc[i]) / q if q > 0 else np.nan)
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    rlo, rhi = np.nanquantile(rboots, [alpha / 2, 1 - alpha / 2])
    return {"post_gap": gap(tp), "pre_gap": gap(tq), "change_in_gap": est, "ci_low": lo, "ci_high": hi,
            "post_ratio": ratio(tp), "pre_ratio": ratio(tq), "ratio_of_ratios": rr, "rr_ci_low": rlo,
            "rr_ci_high": rhi, "rr_undefined_share": float(np.mean(np.isnan(rboots)))}


def balance(pairs: pd.DataFrame) -> pd.DataFrame:
    t, c = pairs[pairs.side == "treated"][COVARIATES], pairs[pairs.side == "control"][COVARIATES]
    sd = np.sqrt((t.var() + c.var()) / 2)
    return pd.DataFrame({"treated": t.mean(), "control": c.mean(), "smd": (t.mean() - c.mean()) / sd})


def analyze(a) -> int:
    d, out = Path(a.data), Path(a.out)
    df = _prepare(d)
    print(f"forest, unprotected: {int(df.treated.sum())} treated, {int((df.treated == 0).sum())} controls")
    unm = df[df.treated == 0]
    raw = pd.DataFrame({"treated": df[df.treated == 1][COVARIATES].mean(), "control": unm[COVARIATES].mean()})
    rows, bals, pairs_by = [], {}, {}

    def scopes(pairs):
        area = pairs[pairs["sample"] == "area"]
        return {"all regions": area, **{r: area[area.region == r] for r in ("Amazon", "Andes", "Coast")},
                "collective": area[area.tpcnv == "COLECTIVO"],
                "individual": pairs[pairs.tpcnv == "INDIVIDUAL"],
                "individual, Amazon": pairs[(pairs.tpcnv == "INDIVIDUAL") & (pairs.region == "Amazon")],
                "individual, Andes": pairs[(pairs.tpcnv == "INDIVIDUAL") & (pairs.region == "Andes")]}

    for window in ("post", "pre"):
        pairs = matched_pairs(df, window)
        pairs_by[window] = scopes(pairs)
        bals[window] = balance(pairs[pairs["sample"] == "area"])
        bals[f"{window}, individual"] = balance(pairs[pairs.tpcnv == "INDIVIDUAL"])
        for scope, sub in pairs_by[window].items():
            if (sub.side == "treated").sum() < 200:
                continue
            rows.append({"window": window, "scope": scope, **hazard_diff(sub, a.n_boot, a.alpha)})
        if window == "post":
            area = pairs[pairs["sample"] == "area"]
            by_cohort = [{"cohort": e, **hazard_diff(area[area.cohort == e], a.n_boot, a.alpha)}
                         for e in COHORTS if (area.cohort == e).sum() > 400]
            pd.DataFrame(by_cohort).to_csv(out / "socio_bosque_parcel_by_cohort.csv", index=False)
    res = pd.DataFrame(rows)
    res.to_csv(out / "socio_bosque_parcel.csv", index=False)
    did = pd.DataFrame([{"scope": sc, "n_agreements": pairs_by["post"][sc].codcnv.nunique(),
                         **change_in_gap(pairs_by["post"][sc], pairs_by["pre"][sc], a.n_boot, a.alpha)}
                        for sc in pairs_by["post"] if (pairs_by["post"][sc].side == "treated").sum() >= 200])
    did.to_csv(out / "socio_bosque_parcel_did.csv", index=False)
    bal = pd.concat({"unmatched": raw.assign(smd=(raw.treated - raw.control) / np.sqrt(
        (df[df.treated == 1][COVARIATES].var() + unm[COVARIATES].var()) / 2)), **bals}, axis=1)
    bal.to_csv(out / "socio_bosque_parcel_balance.csv")
    print(bal.round(3).to_string())
    print(res.round(4).to_string(index=False))
    print(did.round(4).to_string(index=False))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["sample", "extract", "analyze"])
    ap.add_argument("--data", default="data/socio_bosque")
    ap.add_argument("--out", default="docs/results")
    ap.add_argument("--project")
    ap.add_argument("--n-treated", type=int, default=20000)
    ap.add_argument("--n-control", type=int, default=60000)
    ap.add_argument("--per-individual", type=int, default=3, help="extra points per individual agreement")
    ap.add_argument("--chunk", type=int, default=10000)
    ap.add_argument("--wait", action="store_true", help="poll until all extraction tasks finish")
    ap.add_argument("--n-boot", type=int, default=999)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    return {"sample": sample, "extract": extract, "analyze": analyze}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
