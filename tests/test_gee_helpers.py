"""Earth Engine helpers that do not need an authenticated session."""

import pandas as pd

from wiat.data import gee


def test_period_bounds_month_and_year():
    m = gee.period_bounds("2020-01-15", "2020-04-01", "month")
    assert m == [("2020-01-01", "2020-02-01"), ("2020-02-01", "2020-03-01"), ("2020-03-01", "2020-04-01")]
    y = gee.period_bounds("2001-01-01", "2003-01-01", "year")
    assert y == [("2001-01-01", "2002-01-01"), ("2002-01-01", "2003-01-01")]


def test_modis_last_composite_is_short():
    assert gee.modis_composite_days(0, 2021) == 8
    assert gee.modis_composite_days(360, 2021) == 5  # DOY 361 (1-based), non-leap
    assert gee.modis_composite_days(360, 2020) == 6  # leap year


def test_harmonisation_coefficients_shape():
    for arr in (gee.ETM_TO_OLI_SLOPE, gee.ETM_TO_OLI_INTERCEPT, gee.MSI_TO_OLI_SLOPE, gee.MSI_TO_OLI_INTERCEPT):
        assert len(arr) == len(gee.COMMON_BANDS)
    for bands in gee.LANDSAT_BANDS.values():
        assert len(bands) == len(gee.COMMON_BANDS)


def test_features_to_frame():
    info = {"features": [
        {"properties": {"HYBAS_ID": 2, "date": "2001-02-01", "ndvi": 0.5}},
        {"properties": {"HYBAS_ID": 1, "date": "2001-01-01", "ndvi": 0.4}},
    ]}
    df = gee.features_to_frame(info)
    assert df.basin_id.tolist() == [1, 2]
    assert pd.api.types.is_datetime64_any_dtype(df.date)


def test_quality_mask():
    df = pd.DataFrame({"ndvi": [0.5, 0.6], "evi": [0.3, 0.3], "veg_clear_frac": [0.1, 0.9]})
    out = gee.quality_mask(df, 0.3)
    assert out.ndvi.isna().tolist() == [True, False]


def test_composite_overlap_apportions_whole_year():
    # Every day of a year is counted exactly once across its composites and months.
    for year in (2020, 2021):
        starts = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="8D")
        months = gee.period_bounds(f"{year}-01-01", f"{year + 1}-01-01", "month")
        total = sum(gee.composite_overlap_days(c, s, e) for c in starts for s, e in months)
        assert total == (366 if year == 2020 else 365)
    assert gee.composite_overlap_days("2021-01-25", "2021-02-01", "2021-03-01") == 1
    assert gee.composite_overlap_days("2021-12-27", "2021-12-01", "2022-01-01") == 5


def test_load_export(tmp_path):
    path = tmp_path / "export.csv"
    pd.DataFrame({
        "HYBAS_ID": [7, 7], "date": ["2001-02-01", "2001-01-01"], "ndvi": [0.5, 0.6], "evi": [0.3, 0.3],
        "ndmi": [0.1, 0.1], "veg_clear_frac": [0.9, 0.1], "precip_mm": [100.0, 80.0], "et_mm": [60.0, 40.0],
        "pet_mm": [90.0, 90.0],
    }).to_csv(path, index=False)
    df = gee.load_export(path)
    assert df.date.is_monotonic_increasing
    assert df.ndvi.isna().tolist() == [True, False]
    assert df.water_yield_mm.tolist() == [40.0, 40.0]
