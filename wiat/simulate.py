"""Synthetic watershed panels with a known intervention effect.

Real attribution never observes the counterfactual, so the toolkit is validated on
simulated basins. The generator mimics the remote-sensing panel that
:mod:`wiat.data.gee` extracts: monthly CHIRPS-like precipitation, MODIS-like ET and
Landsat/Sentinel-like NDVI/EVI, plus static physiography.

Precipitation follows an interactive fixed-effects (factor) model. Regional modes
such as ENSO load differently on each basin. Basins with different loadings
therefore diverge over time even without an intervention, which is exactly the
violation of parallel trends that synthetic control is meant to handle. Vegetation
responds to lagged water availability, and ET follows vegetation, energy and water
supply.

The intervention (e.g. riparian restoration, páramo protection, reforestation) raises
NDVI along a linear ramp. ET rises through the vegetation coupling, and water yield
(P - ET) falls as a consequence. The true effect is computed by simulating the same
basins with identical random draws but without the intervention.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .panel import aggregate, add_hydrologic_outcomes


@dataclass
class SimulatedStudy:
    panel: pd.DataFrame  # monthly, observed (with intervention)
    counterfactual: pd.DataFrame  # monthly, same draws, no intervention
    static: pd.DataFrame  # one row per basin
    treat_start: dict  # basin_id -> first treated year
    loadings: pd.DataFrame  # latent climate loadings (unobservable in practice)

    def annual(self) -> pd.DataFrame:
        return aggregate(self.panel, "year")

    def true_effect(self, outcome: str, freq: str = "year") -> pd.DataFrame:
        """True per-period effect (observed - counterfactual) for each treated basin."""
        obs = aggregate(self.panel, freq) if freq else self.panel
        cf = aggregate(self.counterfactual, freq) if freq else self.counterfactual
        key = ["basin_id", "date"]
        m = obs[key + [outcome]].merge(cf[key + [outcome]], on=key, suffixes=("", "_cf"))
        m = m[m.basin_id.isin(self.treat_start)]
        m["effect"] = m[outcome] - m[f"{outcome}_cf"]
        return m.pivot(index="date", columns="basin_id", values="effect")

    def true_att(self, outcome: str) -> float:
        """Average annual effect over treated basin-years after each start."""
        eff = self.true_effect(outcome)
        vals = [eff.loc[eff.index >= s, b].mean() for b, s in self.treat_start.items()]
        return float(np.mean(vals))


def simulate_watersheds(
    n_basins: int = 40,
    start_year: int = 2001,
    n_years: int = 24,
    treat_start: Mapping | None = None,
    ndvi_effect: float = 0.08,
    ramp_years: float = 4.0,
    treated_extremity: float = 0.0,
    loading_sd: float = 0.35,
    n_factors: int = 2,
    seed: int | None = 0,
) -> SimulatedStudy:
    """Simulate a regional set of basins with one or more intervention basins.

    Parameters
    ----------
    treat_start
        ``{basin_id: start_year}``. Defaults to ``{"B000": start_year + n_years // 2}``.
        Basin ids are ``"B000"``, ``"B001"``, ...
    ndvi_effect
        Full-ramp NDVI increase caused by the intervention.
    treated_extremity
        Shifts the treated basins' latent loadings, elevation, precipitation and
        forest cover by this many SDs toward the edge of the regional distribution.
        Larger values mean fewer natural matches. Basins that most need investment
        (degraded, steep, high-rainfall headwaters) usually sit at that edge.
    """
    rng = np.random.default_rng(seed)
    ids = [f"B{i:03d}" for i in range(n_basins)]
    if treat_start is None:
        treat_start = {ids[0]: start_year + n_years // 2}
    treated_mask = np.array([b in treat_start for b in ids])
    T = n_years * 12
    dates = pd.date_range(f"{start_year}-01-01", periods=T, freq="MS")
    month = dates.month.to_numpy()
    years = dates.year.to_numpy()

    # Static physiography (correlated, Andean-foothill-like ranges).
    shift = treated_extremity * treated_mask
    lat = rng.uniform(-18, 8, n_basins)
    elev = np.clip(np.exp(rng.normal(np.log(1500), 0.6, n_basins) + 0.6 * shift), 100, 4800)
    slope = np.clip(4 + 0.006 * elev + rng.normal(0, 3, n_basins) + 3 * shift, 0.5, 45)
    area = np.exp(rng.normal(np.log(400), 0.8, n_basins))
    forest = np.clip(rng.beta(4, 3, n_basins) - 0.15 * shift, 0.02, 0.98)
    map_mm = np.exp(rng.normal(np.log(1400), 0.35, n_basins) + 0.35 * shift) * (1 - 0.00005 * elev)

    # Latent climate loadings: partly predictable from location (so matching helps),
    # partly idiosyncratic (so it never fully succeeds).
    L = np.column_stack(
        [0.4 * np.sin(np.deg2rad(lat) * 4) + rng.normal(0, loading_sd, n_basins) for _ in range(n_factors)]
    )
    L = L + 1.5 * loading_sd * shift[:, None]

    # Regional factors: AR(1) anomalies (e.g. ENSO, Atlantic modes).
    F = np.zeros((T, n_factors))
    for k in range(n_factors):
        rho = 0.92 - 0.1 * k
        e = rng.normal(0, np.sqrt(1 - rho**2) * 0.45, T)
        for t in range(1, T):
            F[t, k] = rho * F[t - 1, k] + e[t]

    amp = np.clip(0.3 + 0.5 * np.abs(lat) / 18 + rng.normal(0, 0.1, n_basins), 0.1, 0.95)
    phase = np.where(lat < 0, 1.0, 7.0)
    seas = 1 + amp[:, None] * np.cos(2 * np.pi * (month[None] - phase[:, None]) / 12)
    anomaly = (L @ F.T) + 0.25 * rng.normal(size=(n_basins, T))
    precip = (map_mm / 12)[:, None] * np.clip(seas, 0.03, None) * np.exp(anomaly - 0.5 * 0.25**2)

    pet_annual = np.clip(1500 - 0.22 * elev + rng.normal(0, 60, n_basins), 500, None)
    pet = (pet_annual / 12)[:, None] * (1 + 0.18 * np.cos(2 * np.pi * (month[None] - phase[:, None] - 1) / 12))

    # Water availability: exponential moving average of precipitation.
    S = np.zeros_like(precip)
    S[:, 0] = precip[:, 0]
    for t in range(1, T):
        S[:, t] = 0.6 * S[:, t - 1] + 0.4 * precip[:, t]

    ndvi0 = 0.30 + 0.45 * forest
    trend = rng.normal(0.0015, 0.002, n_basins)
    ndvi_noise = 0.02 * rng.normal(size=(n_basins, T))
    et_noise = rng.normal(0, 0.05, size=(n_basins, T))
    evi_noise = 0.015 * rng.normal(size=(n_basins, T))
    tyears = (years - start_year)[None].astype(float)

    ramp = np.zeros((n_basins, T))
    for i, b in enumerate(ids):
        if b in treat_start:
            elapsed = (dates - pd.Timestamp(f"{int(treat_start[b])}-01-01")).days.to_numpy() / 365.25
            ramp[i] = np.clip(elapsed / ramp_years, 0, 1) * (elapsed >= 0)

    def realise(intervention: np.ndarray) -> pd.DataFrame:
        wet = np.tanh((S - (map_mm / 12)[:, None]) / (map_mm / 12)[:, None])
        ndvi = np.clip(ndvi0[:, None] + 0.12 * wet + trend[:, None] * tyears + ndvi_effect * intervention + ndvi_noise, 0.05, 0.95)
        demand = pet * (0.15 + 0.95 * ndvi) * (1 - np.exp(-S / (0.6 * pet)))
        et = np.clip(np.minimum(demand, precip + 0.5 * S) * np.exp(et_noise), 0, None)
        evi = np.clip(0.62 * ndvi - 0.03 + evi_noise, 0, 1)
        return pd.DataFrame(
            {
                "basin_id": np.repeat(ids, T),
                "date": np.tile(dates, n_basins),
                "precip_mm": precip.ravel(),
                "et_mm": et.ravel(),
                "pet_mm": pet.ravel(),
                "ndvi": ndvi.ravel(),
                "evi": evi.ravel(),
            }
        )

    obs = add_hydrologic_outcomes(realise(ramp))
    cf = add_hydrologic_outcomes(realise(np.zeros_like(ramp)))
    static = pd.DataFrame(
        {
            "latitude": lat,
            "elevation_m": elev,
            "slope_deg": slope,
            "area_km2": area,
            "forest_frac": forest,
        },
        index=pd.Index(ids, name="basin_id"),
    )
    loadings = pd.DataFrame(L, index=static.index, columns=[f"lambda_{k}" for k in range(n_factors)])
    return SimulatedStudy(obs, cf, static, dict(treat_start), loadings)
