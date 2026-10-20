"""Panel utilities: validation, hydrologic outcome derivation, temporal aggregation.

All estimators in the toolkit consume a *long* panel: one row per (basin, period),
with a unit column (default ``basin_id``), a time column (default ``date``) and one
column per variable. Fluxes are in millimetres over the period (``*_mm``); indices
(NDVI, EVI, ratios) are period means.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

UNIT = "basin_id"
TIME = "date"


def validate_panel(df: pd.DataFrame, unit: str = UNIT, time: str = TIME) -> pd.DataFrame:
    """Return a sorted copy of ``df`` after checking structural assumptions."""
    missing = {unit, time} - set(df.columns)
    if missing:
        raise KeyError(f"panel is missing required column(s): {sorted(missing)}")
    if df.duplicated([unit, time]).any():
        dups = df.loc[df.duplicated([unit, time], keep=False), [unit, time]].head()
        raise ValueError(f"panel has duplicate (unit, time) rows, e.g.:\n{dups}")
    out = df.copy()
    if not np.issubdtype(out[time].dtype, np.number):
        out[time] = pd.to_datetime(out[time])
    return out.sort_values([unit, time]).reset_index(drop=True)


def add_hydrologic_outcomes(
    df: pd.DataFrame, precip: str = "precip_mm", et: str = "et_mm"
) -> pd.DataFrame:
    """Add water-balance outcomes derived from precipitation and evapotranspiration.

    * ``water_yield_mm`` = P - ET: the long-run water-balance proxy for runoff plus
      storage change. Over annual or longer periods storage change is small, so this
      approximates the water a basin releases downstream.
    * ``evaporative_index`` = ET / P: the Budyko evaporative index. Interventions that
      change vegetation water use shift this ratio. Normalising by P removes much of the
      interannual climate signal that otherwise dominates basin-level fluxes.
    """
    out = df.copy()
    if precip in out and et in out:
        out["water_yield_mm"] = out[precip] - out[et]
        with np.errstate(divide="ignore", invalid="ignore"):
            out["evaporative_index"] = np.where(out[precip] > 0, out[et] / out[precip], np.nan)
    return out


def _default_how(columns: Iterable[str]) -> dict[str, str]:
    return {c: ("sum" if c.endswith("_mm") else "mean") for c in columns}


def aggregate(
    df: pd.DataFrame,
    freq: str = "year",
    how: Mapping[str, str] | None = None,
    unit: str = UNIT,
    time: str = TIME,
    water_year_start_month: int = 1,
    min_months: int = 12,
    min_months_mean: int | None = 4,
) -> pd.DataFrame:
    """Aggregate a monthly panel to annual (or water-year) periods.

    Fluxes (columns ending in ``_mm``) are summed; everything else is averaged, unless
    ``how`` overrides it. A period with fewer than ``min_months`` valid months gets NaN
    rather than a biased partial sum. Averaged columns (vegetation indices) need only
    ``min_months_mean`` valid months, because optical records have cloud gaps in most
    wet seasons and a mean over the clear months is still an unbiased-in-expectation
    level, whereas a sum over them is not. Pass ``None`` to apply ``min_months`` to
    every column. Derived ratios such as ``evaporative_index`` are recomputed from the
    aggregated fluxes so they stay internally consistent.

    ``water_year_start_month`` lets hydrologic years follow regional seasonality
    (e.g. 10 for an October-September water year). The result's integer ``time``
    column is the calendar year in which the water year *ends*.
    """
    if freq not in {"year", "water_year"}:
        raise ValueError("freq must be 'year' or 'water_year'")
    d = validate_panel(df, unit, time)
    dates = pd.to_datetime(d[time])
    shift = 0 if freq == "year" else (12 - water_year_start_month + 1) % 12
    period = (dates + pd.DateOffset(months=shift)).dt.year
    value_cols = [c for c in d.columns if c not in (unit, time) and np.issubdtype(d[c].dtype, np.number)]
    rules = _default_how(value_cols)
    if how:
        rules.update(how)
    rules = {c: r for c, r in rules.items() if c in value_cols}
    d = d.assign(_period=period.values)
    grouped = d.groupby([unit, "_period"])
    out = grouped[list(rules)].agg(rules)
    counts = grouped[list(rules)].count()
    need = pd.Series({c: min_months if r == "sum" or min_months_mean is None else min(min_months_mean, min_months)
                      for c, r in rules.items()})
    out = out.where(counts.ge(need, axis=1))
    out = out.reset_index().rename(columns={"_period": time})
    if "evaporative_index" in out and {"precip_mm", "et_mm"} <= set(out.columns):
        out["evaporative_index"] = out["et_mm"] / out["precip_mm"]
    return out


def deseasonalize(
    df: pd.DataFrame,
    columns: Iterable[str],
    baseline_end,
    unit: str = UNIT,
    time: str = TIME,
) -> pd.DataFrame:
    """Subtract each basin's month-of-year climatology computed over the baseline.

    The climatology is estimated only from periods before ``baseline_end`` so that the
    intervention does not leak into the reference.
    """
    d = validate_panel(df, unit, time)
    dates = pd.to_datetime(d[time])
    month = dates.dt.month
    base = dates < pd.Timestamp(baseline_end)
    for c in columns:
        clim = d.loc[base].groupby([d.loc[base, unit], month[base]])[c].mean()
        key = pd.MultiIndex.from_arrays([d[unit], month])
        d[f"{c}_anom"] = d[c].values - clim.reindex(key).values
    return d


def to_wide(df: pd.DataFrame, value: str, unit: str = UNIT, time: str = TIME) -> pd.DataFrame:
    """Pivot a long panel to a (time x unit) matrix."""
    return df.pivot(index=time, columns=unit, values=value).sort_index()


def balance(
    df: pd.DataFrame, columns: Iterable[str], unit: str = UNIT, time: str = TIME
) -> tuple[pd.DataFrame, list]:
    """Keep only units observed (non-NaN in ``columns``) in every period.

    Returns the balanced panel and the list of dropped units.
    """
    d = validate_panel(df, unit, time)
    cols = list(columns)
    n_periods = d[time].nunique()
    ok = d.dropna(subset=cols).groupby(unit)[time].nunique() == n_periods
    keep = ok[ok].index
    dropped = sorted(set(d[unit]) - set(keep))
    return d[d[unit].isin(keep)].reset_index(drop=True), dropped


def period_index(times: Iterable) -> dict:
    """Map sorted unique period labels to consecutive integers."""
    uniq = sorted(pd.unique(pd.Series(list(times))))
    return {t: i for i, t in enumerate(uniq)}


def treatment_matrix(
    units: list, times: list, treat_start: Mapping, time_is_ordered: bool = True
) -> np.ndarray:
    """Binary (unit x time) matrix that is 1 from each treated unit's start period on."""
    D = np.zeros((len(units), len(times)))
    tarr = np.asarray(times)
    for i, u in enumerate(units):
        if u in treat_start and treat_start[u] is not None:
            start = treat_start[u]
            if isinstance(tarr[0], (pd.Timestamp, np.datetime64)):
                start = pd.Timestamp(start)
            D[i] = (tarr >= start).astype(float)
    return D
