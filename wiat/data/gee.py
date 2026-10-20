"""Google Earth Engine extraction of basin panels.

Produces the long panel the rest of the toolkit consumes: one row per
(HydroBASINS unit, period) with

* ``ndvi``, ``evi``, ``ndmi``: median composites of cloud-masked Landsat 5/7/8/9
  Collection 2 Level 2 and Sentinel-2 L2A surface reflectance, harmonised to OLI.
  TM/ETM+ use the Roy et al. (2016) OLS transformation; MSI uses the HLS v1.4
  bandpass adjustment. ``veg_clear_frac`` is the share of the basin with at least
  one clear observation in the period.
* ``precip_mm``: CHIRPS v2 daily, summed over the period.
* ``et_mm``, ``pet_mm``: MODIS MOD16A2GF v061. Each 8-day sum is converted to a
  daily rate (the year's last composite covers only 5 or 6 days) and multiplied by
  the days it overlaps the period, so months receive exactly their share.

and a static covariate table (NASADEM terrain, Hansen forest cover, basin area and
latitude) for donor matching.

Nothing here touches Earth Engine at import time. The pure helpers (period bounds,
composite lengths, harmonisation coefficients, table conversion, quality masking)
work without ``earthengine-api`` installed.
"""

from __future__ import annotations

import calendar
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

try:  # optional dependency: pip install "wiat[gee]"
    import ee
except ImportError:  # pragma: no cover - exercised only without the extra
    ee = None

# --- Assets -------------------------------------------------------------------

HYDROBASINS = "WWF/HydroSHEDS/v1/Basins/hybas_{level}"
COUNTRIES = "USDOS/LSIB_SIMPLE/2017"
CHIRPS = "UCSB-CHG/CHIRPS/DAILY"
MOD16 = "MODIS/061/MOD16A2GF"
NASADEM = "NASA/NASADEM_HGT/001"
HANSEN = "UMD/hansen/global_forest_change_2025_v1_13"
S2 = "COPERNICUS/S2_SR_HARMONIZED"
LANDSAT = {
    "L5": "LANDSAT/LT05/C02/T1_L2",
    "L7": "LANDSAT/LE07/C02/T1_L2",
    "L8": "LANDSAT/LC08/C02/T1_L2",
    "L9": "LANDSAT/LC09/C02/T1_L2",
}

# --- Harmonisation --------------------------------------------------------------

COMMON_BANDS = ("blue", "green", "red", "nir", "swir1", "swir2")
LANDSAT_BANDS = {
    "L5": ("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"),
    "L7": ("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"),
    "L8": ("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"),
    "L9": ("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"),
}
S2_BANDS = ("B2", "B3", "B4", "B8A", "B11", "B12")  # B8A matches OLI's narrow NIR

# Roy et al. (2016), Remote Sens. Environ. 185, Table 2: OLS ETM+ -> OLI surface reflectance.
ETM_TO_OLI_SLOPE = (0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071)
ETM_TO_OLI_INTERCEPT = (0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172)
# HLS v1.4 user guide: MSI -> OLI bandpass adjustment (Sentinel-2A coefficients).
MSI_TO_OLI_SLOPE = (0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.0030)
MSI_TO_OLI_INTERCEPT = (-0.0040, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012)

LANDSAT_SCALE, LANDSAT_OFFSET = 2.75e-05, -0.2
S2_SCALE = 1e-4
S2_CLEAR_SCL = (4, 5, 6)  # vegetation, bare soil, water; drops cloud, shadow, cirrus, snow
MOD16_SCALE = 0.1  # kg m-2 per composite == mm
MOD16_FILL_MIN = 32700
MODIS_COMPOSITE = 8

VEG_COLUMNS = ("ndvi", "evi", "ndmi")
PANEL_COLUMNS = (*VEG_COLUMNS, "veg_clear_frac", "precip_mm", "et_mm", "pet_mm")


@dataclass(frozen=True)
class PanelSpec:
    """What to extract: period range (end exclusive), period length and reduction scale."""

    start: str = "2001-01-01"
    end: str = "2025-01-01"
    freq: str = "month"
    scale: float = 500
    sensors: tuple = ("L5", "L7", "L8", "L9", "S2")
    max_cloud: float = 80.0  # scene-level pre-filter; pixel masks do the real work
    tile_scale: int = 4


# --- Pure helpers -----------------------------------------------------------------


def period_bounds(start, end, freq: str = "month") -> list[tuple[str, str]]:
    """Half-open ``[start, end)`` periods as ISO date strings.

    ``start`` is floored to the beginning of its month (or year); periods that would
    extend past ``end`` are dropped.
    """
    if freq not in {"month", "year"}:
        raise ValueError("freq must be 'month' or 'year'")
    s, e = pd.Timestamp(start), pd.Timestamp(end)
    s = s.replace(day=1) if freq == "month" else s.replace(month=1, day=1)
    step = pd.DateOffset(months=1) if freq == "month" else pd.DateOffset(years=1)
    out = []
    while s + step <= e:
        out.append((s.strftime("%Y-%m-%d"), (s + step).strftime("%Y-%m-%d")))
        s = s + step
    return out


def modis_composite_days(doy0: int, year: int) -> int:
    """Days covered by the MODIS 8-day composite starting on 0-based day ``doy0``.

    Composites restart on 1 January, so the last one of the year is truncated to
    5 days (6 in leap years).
    """
    days_in_year = 366 if calendar.isleap(year) else 365
    return int(min(MODIS_COMPOSITE, days_in_year - doy0))


def composite_overlap_days(composite_start, period_start, period_end) -> int:
    """Days a MODIS composite beginning on ``composite_start`` shares with a period."""
    c0 = pd.Timestamp(composite_start)
    c1 = c0 + pd.Timedelta(days=modis_composite_days(c0.dayofyear - 1, c0.year))
    lo, hi = max(c0, pd.Timestamp(period_start)), min(c1, pd.Timestamp(period_end))
    return max(0, (hi - lo).days)


def features_to_frame(info: dict, id_col: str = "HYBAS_ID") -> pd.DataFrame:
    """Turn ``FeatureCollection.getInfo()`` output into a sorted long panel."""
    df = pd.DataFrame([f["properties"] for f in info.get("features", [])])
    if df.empty:
        return df
    df = df.rename(columns={id_col: "basin_id"})
    if "date" in df:
        df["date"] = pd.to_datetime(df["date"])
    keys = [c for c in ("basin_id", "date") if c in df]
    return df.sort_values(keys).reset_index(drop=True)


def quality_mask(df: pd.DataFrame, min_clear: float = 0.3) -> pd.DataFrame:
    """Blank vegetation indices in periods where too little of the basin was cloud-free.

    A median over a few clear pixels in one corner of a basin is not a basin NDVI.
    Blanked periods become gaps that ``wiat.panel.balance`` later drops, which thins
    the donor pool; that cost is reported, not hidden.
    """
    out = df.copy()
    if "veg_clear_frac" not in out:
        return out
    bad = out["veg_clear_frac"].fillna(0) < min_clear
    for c in VEG_COLUMNS:
        if c in out:
            out.loc[bad, c] = np.nan
    return out


# --- Earth Engine -------------------------------------------------------------------


def _ee():
    if ee is None:
        raise ImportError('Earth Engine extraction needs the optional extra: pip install "wiat[gee]"')
    return ee


def initialize(project: str | None = None) -> None:
    """Initialise Earth Engine (run ``earthengine authenticate`` once beforehand)."""
    _ee().Initialize(project=project)


def basins(level: int = 8, region: Sequence[float] | None = None, country: str | None = None):
    """HydroBASINS units at ``level`` (1-12) intersecting a bbox and/or a country."""
    e = _ee()
    if not 1 <= int(level) <= 12:
        raise ValueError("HydroBASINS level must be 1-12")
    fc = e.FeatureCollection(HYDROBASINS.format(level=int(level)))
    if region is not None:
        fc = fc.filterBounds(e.Geometry.Rectangle(list(region)))
    if country is not None:
        geom = e.FeatureCollection(COUNTRIES).filter(e.Filter.eq("country_na", country)).geometry()
        fc = fc.filterBounds(geom)
    return fc


def basin_at(lon: float, lat: float, level: int = 8):
    """The HydroBASINS unit containing a point."""
    e = _ee()
    fc = e.FeatureCollection(HYDROBASINS.format(level=int(level)))
    return fc.filterBounds(e.Geometry.Point([lon, lat])).first()


def basin_table(fc) -> pd.DataFrame:
    """Attribute table (id, downstream id, area) used for connectivity exclusion."""
    cols = ["HYBAS_ID", "NEXT_DOWN", "SUB_AREA"]
    info = fc.select(cols, None, False).getInfo()
    return pd.DataFrame([f["properties"] for f in info["features"]], columns=cols)


def _harmonise(img, slope, intercept):
    e = _ee()
    return img.multiply(e.Image.constant(list(slope))).add(e.Image.constant(list(intercept)))


def _indices(sr):
    """NDVI, EVI, NDMI from a 6-band OLI-equivalent reflectance image."""
    ndvi = sr.normalizedDifference(["nir", "red"]).rename("ndvi")
    evi = sr.expression(
        "2.5 * (n - r) / (n + 6 * r - 7.5 * b + 1)",
        {"n": sr.select("nir"), "r": sr.select("red"), "b": sr.select("blue")},
    ).rename("evi")
    ndmi = sr.normalizedDifference(["nir", "swir1"]).rename("ndmi")
    return _ee().Image(ndvi.addBands([evi, ndmi]).toFloat().copyProperties(sr, ["system:time_start"]))


def _landsat(sensor: str, geom, start, end, max_cloud: float):
    e = _ee()
    bands, oli = list(LANDSAT_BANDS[sensor]), sensor in ("L8", "L9")

    def prep(img):
        qa = img.select("QA_PIXEL")
        # bits: 1 dilated cloud, 2 cirrus, 3 cloud, 4 shadow, 5 snow
        clear = qa.bitwiseAnd(0b111110).eq(0).And(img.select("QA_RADSAT").eq(0))
        sr = img.select(bands, list(COMMON_BANDS)).multiply(LANDSAT_SCALE).add(LANDSAT_OFFSET)
        if not oli:
            sr = _harmonise(sr, ETM_TO_OLI_SLOPE, ETM_TO_OLI_INTERCEPT)
        sr = sr.updateMask(clear).updateMask(sr.reduce(e.Reducer.min()).gt(0))
        return _indices(e.Image(sr.copyProperties(img, ["system:time_start"])))

    return (e.ImageCollection(LANDSAT[sensor]).filterBounds(geom).filterDate(start, end)
            .filter(e.Filter.lt("CLOUD_COVER", max_cloud)).map(prep))


def _sentinel2(geom, start, end, max_cloud: float):
    e = _ee()

    def prep(img):
        scl = img.select("SCL")
        clear = scl.eq(S2_CLEAR_SCL[0])
        for c in S2_CLEAR_SCL[1:]:
            clear = clear.Or(scl.eq(c))
        sr = img.select(list(S2_BANDS), list(COMMON_BANDS)).multiply(S2_SCALE)
        sr = _harmonise(sr, MSI_TO_OLI_SLOPE, MSI_TO_OLI_INTERCEPT).updateMask(clear)
        return _indices(e.Image(sr.copyProperties(img, ["system:time_start"])))

    return (e.ImageCollection(S2).filterBounds(geom).filterDate(start, end)
            .filter(e.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", max_cloud)).map(prep))


def vegetation_composite(geom, start, end, sensors=PanelSpec.sensors, max_cloud: float = 80.0):
    """Median NDVI/EVI/NDMI over ``[start, end)`` plus a clear-observation fraction band."""
    e = _ee()
    col = None
    for s in sensors:
        c = _sentinel2(geom, start, end, max_cloud) if s == "S2" else _landsat(s, geom, start, end, max_cloud)
        col = c if col is None else col.merge(c)
    empty = e.Image.constant([0, 0, 0]).toFloat().rename(list(VEG_COLUMNS)).updateMask(0)
    comp = e.Image(e.Algorithms.If(col.size().gt(0), col.median(), empty)).select(list(VEG_COLUMNS))
    clear = comp.select("ndvi").mask().gt(0).unmask(0).rename("veg_clear_frac")
    return comp.addBands(clear)


def precipitation(start, end):
    """CHIRPS precipitation summed over ``[start, end)`` (mm)."""
    e = _ee()
    return e.ImageCollection(CHIRPS).filterDate(start, end).select("precipitation").sum().rename("precip_mm")


def evapotranspiration(start, end):
    """MOD16 ET and PET apportioned to ``[start, end)`` by day overlap (mm).

    Server-side counterpart of :func:`composite_overlap_days`.
    """
    e = _ee()
    s, t = e.Date(start), e.Date(end)

    def weight(img):
        c0 = e.Date(img.get("system:time_start"))
        year_end = e.Date.fromYMD(c0.get("year").add(1), 1, 1)
        length = year_end.difference(c0, "day").min(MODIS_COMPOSITE)
        c1 = c0.advance(length, "day")
        lo = e.Number(c0.millis()).max(s.millis())
        hi = e.Number(c1.millis()).min(t.millis())
        overlap = hi.subtract(lo).divide(86_400_000).max(0)
        vals = img.select(["ET", "PET"])
        vals = vals.updateMask(vals.lt(MOD16_FILL_MIN))
        # cast: zero-overlap composites otherwise get a degenerate Float<0, 0> type
        return vals.multiply(MOD16_SCALE).divide(length).multiply(overlap).toFloat()

    col = e.ImageCollection(MOD16).filterDate(s.advance(-MODIS_COMPOSITE, "day"), t).map(weight)
    return col.sum().rename(["et_mm", "pet_mm"])


def period_image(geom, start, end, spec: PanelSpec):
    """All panel variables for one period as a single multi-band image."""
    return (vegetation_composite(geom, start, end, spec.sensors, spec.max_cloud)
            .addBands(precipitation(start, end))
            .addBands(evapotranspiration(start, end)))


def _period_table(fc, start, end, spec: PanelSpec):
    e = _ee()
    img = period_image(fc.geometry().bounds(1000), start, end, spec)
    red = img.reduceRegions(collection=fc.select(["HYBAS_ID"]), reducer=e.Reducer.mean(),
                            scale=spec.scale, tileScale=spec.tile_scale)
    return red.map(lambda f: e.Feature(None, f.toDictionary()).set("date", start))


def fetch_panel(fc, spec: PanelSpec, verbose: bool = True) -> pd.DataFrame:
    """Interactive extraction, one ``getInfo`` per period.

    Fine for a few hundred basins over a study region; use :func:`export_panel` for
    anything larger, where interactive requests hit Earth Engine's memory and time limits.
    """
    frames = []
    periods = period_bounds(spec.start, spec.end, spec.freq)
    for i, (s, t) in enumerate(periods):
        frames.append(features_to_frame(_period_table(fc, s, t, spec).getInfo()))
        if verbose and (i + 1) % 12 == 0:
            print(f"  {i + 1}/{len(periods)} periods")
    df = pd.concat(frames, ignore_index=True)
    for c in PANEL_COLUMNS:
        if c not in df:
            df[c] = np.nan
    return df[["basin_id", "date", *PANEL_COLUMNS]].sort_values(["basin_id", "date"]).reset_index(drop=True)


def export_panel(fc, spec: PanelSpec, description: str = "wiat_panel", folder: str | None = None):
    """Batch-export the panel to Google Drive as CSV; returns the started task.

    Turn the downloaded CSV into an analysis-ready panel with :func:`load_export`.
    """
    e = _ee()
    tables = [_period_table(fc, s, t, spec) for s, t in period_bounds(spec.start, spec.end, spec.freq)]
    merged = e.FeatureCollection(tables).flatten()
    task = e.batch.Export.table.toDrive(
        collection=merged, description=description, folder=folder, fileFormat="CSV",
        selectors=["HYBAS_ID", "date", *PANEL_COLUMNS],
    )
    task.start()
    return task


def load_export(path, min_clear: float = 0.3) -> pd.DataFrame:
    """Read a CSV written by :func:`export_panel`, quality-mask it and add outcomes."""
    from ..panel import add_hydrologic_outcomes

    df = pd.read_csv(path).rename(columns={"HYBAS_ID": "basin_id"})
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values(["basin_id", "date"]).reset_index(drop=True)
    return add_hydrologic_outcomes(quality_mask(df, min_clear))


def static_covariates(fc, scale: float = 250, forest_threshold: float = 30) -> pd.DataFrame:
    """Basin physiography for matching, indexed by ``basin_id``.

    ``forest_frac`` is the share of the basin with >= ``forest_threshold``% tree cover
    in 2000; ``forest_loss_frac`` is the share of that forest lost since.
    """
    e = _ee()
    dem = e.Image(NASADEM).select("elevation")
    gfc = e.Image(HANSEN)
    forest = gfc.select("treecover2000").gte(forest_threshold)
    img = (dem.rename("elevation_m")
           .addBands(e.Terrain.slope(dem).rename("slope_deg"))
           .addBands(forest.rename("forest_frac"))
           .addBands(gfc.select("loss").And(forest).rename("forest_loss_px")))
    red = img.reduceRegions(collection=fc.select(["HYBAS_ID", "SUB_AREA"]), reducer=e.Reducer.mean(),
                            scale=scale, tileScale=4)
    red = red.map(lambda f: e.Feature(None, f.toDictionary()).set(
        "latitude", f.geometry().centroid(100).coordinates().get(1)))
    df = features_to_frame(red.getInfo()).rename(columns={"SUB_AREA": "area_km2"})
    df["forest_loss_frac"] = df["forest_loss_px"] / df["forest_frac"].replace(0, np.nan)
    cols = ["latitude", "elevation_m", "slope_deg", "area_km2", "forest_frac", "forest_loss_frac"]
    return df.set_index("basin_id")[cols]
