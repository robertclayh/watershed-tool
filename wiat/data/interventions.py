"""Programme footprints: parcel polygons -> basin-level treatment exposure.

Single-basin studies use one reference point per programme. Programmes that publish
parcel polygons with enrollment dates (see ``docs/intervention_data_sources.md``)
support a staggered design instead: each basin's *share* of area under the
programme, cumulative by year, defines whether and when it is treated.

Needs the optional ``geo`` extra (geopandas).
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from collections.abc import Mapping

import pandas as pd

SOCIO_BOSQUE_WFS = "http://ide.ambiente.gob.ec:8080/geoserver/ows"
SOCIO_BOSQUE_LAYER = "mae_ide:v_hc005_area_bajo_conservacion_a_anonimo"
SOCIO_BOSQUE_CRS = 32717  # UTM 17S, as served


def _gpd():
    try:
        import geopandas
    except ImportError as e:  # pragma: no cover - exercised only without the extra
        raise ImportError('programme footprints need the optional extra: pip install "wiat[geo]"') from e
    return geopandas


def fetch_wfs(url: str, layer: str, page_size: int = 1000, crs: int | None = None, sort_by: str | None = None,
              timeout: int = 600):
    """Download every feature of a WFS layer as a GeoDataFrame, ``page_size`` at a time.

    Paging needs a stable order, so pass ``sort_by`` (any unique attribute) for
    servers that do not order features by default.
    """
    gpd = _gpd()
    frames, start, total, crs_name = [], 0, None, None
    while total is None or start < total:
        q = {"service": "WFS", "version": "1.0.0", "request": "GetFeature", "typeName": layer,
             "outputFormat": "application/json", "maxFeatures": page_size, "startIndex": start}
        if sort_by:
            q["sortBy"] = sort_by
        with urllib.request.urlopen(f"{url}?{urllib.parse.urlencode(q)}", timeout=timeout) as r:
            page = json.load(r)
        feats = page.get("features", [])
        total = page.get("totalFeatures", len(feats)) if total is None else total
        crs_name = crs_name or (page.get("crs") or {}).get("properties", {}).get("name")
        if not feats:
            break
        frames.append(gpd.GeoDataFrame.from_features(feats))
        start += len(feats)
    out = pd.concat(frames, ignore_index=True)
    epsg = crs or (int(crs_name.rsplit(":", 1)[-1]) if crs_name else 4326)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=f"EPSG:{epsg}")


def socio_bosque(timeout: int = 600):
    """Ecuador's Socio Bosque conservation agreements, with a ``year`` of signing."""
    gdf = fetch_wfs(SOCIO_BOSQUE_WFS, SOCIO_BOSQUE_LAYER, crs=SOCIO_BOSQUE_CRS, sort_by="codcnv", timeout=timeout)
    gdf["year"] = pd.to_datetime(gdf["fec"].astype(str).str[:10]).dt.year
    return gdf


def basin_exposure(parcels, basins, year_col: str = "year", id_col: str = "HYBAS_ID",
                   equal_area_crs: int | None = None) -> pd.DataFrame:
    """Share of each basin's area under the programme, cumulative by year.

    Returns a (basin x year) DataFrame for every basin that intersects a parcel,
    with years as columns. Overlapping parcels are dissolved within each year so
    their area is not double counted. Areas are computed in ``equal_area_crs``
    (default: the parcels' projected CRS, or World Cylindrical Equal Area).
    """
    gpd = _gpd()
    crs = equal_area_crs or (parcels.crs if parcels.crs and parcels.crs.is_projected else 6933)
    p = parcels[[year_col, "geometry"]].to_crs(crs)
    p = p.assign(geometry=p.geometry.buffer(0))
    b = basins[[id_col, "geometry"]].to_crs(crs)
    area = pd.Series(b.area.to_numpy(), index=b[id_col])
    years = sorted(p[year_col].dropna().unique())
    out = {}
    for y in years:
        upto = gpd.GeoDataFrame(geometry=[p[p[year_col] <= y].union_all()], crs=crs)
        ov = gpd.overlay(b, upto, how="intersection", keep_geom_type=True)
        out[int(y)] = pd.Series(ov.area.to_numpy(), index=ov[id_col]).groupby(level=0).sum()
    share = pd.DataFrame(out).fillna(0.0)
    share = share.div(area.reindex(share.index), axis=0).clip(upper=1.0)
    share.index.name = id_col
    return share.cummax(axis=1)


def treatment_years(share: pd.DataFrame, threshold: float = 0.2, lag: int = 1,
                    max_year: int | None = None) -> dict:
    """``{basin: first treated year}`` for basins whose share reaches ``threshold``.

    The first treated year is the crossing year plus ``lag`` (default 1: agreements
    are signed through the year, so the first *fully* enrolled year is the next).
    Basins crossing after ``max_year`` are left out.
    """
    hit = share.ge(threshold)
    first = hit.idxmax(axis=1)[hit.any(axis=1)].astype(int) + lag
    if max_year is not None:
        first = first[first <= max_year]
    return first.to_dict()


def control_pool(share: pd.DataFrame, all_basins, max_share: float = 0.01) -> list:
    """Basins never above ``max_share`` enrolled: candidate never-treated controls."""
    peak = share.max(axis=1) if len(share.columns) else pd.Series(dtype=float)
    exposed = set(peak[peak > max_share].index)
    return [b for b in all_basins if b not in exposed]


def summarize_exposure(share: pd.DataFrame, thresholds: Mapping | list = (0.1, 0.2, 0.5)) -> pd.DataFrame:
    """Number of basins first reaching each threshold, by year."""
    rows = {}
    for th in thresholds:
        hit = share.ge(th)
        rows[f">= {th:.0%}"] = hit.idxmax(axis=1)[hit.any(axis=1)].value_counts().sort_index()
    return pd.DataFrame(rows).fillna(0).astype(int).rename_axis("year")


__all__ = ["basin_exposure", "control_pool", "fetch_wfs", "socio_bosque", "summarize_exposure",
           "treatment_years"]
