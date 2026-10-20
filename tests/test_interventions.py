"""Programme footprints -> basin exposure (needs the geo extra)."""

import numpy as np
import pandas as pd
import pytest

gpd = pytest.importorskip("geopandas")
from shapely.geometry import box  # noqa: E402

from wiat.data.gee import forest_frame  # noqa: E402
from wiat.data.interventions import basin_exposure, control_pool, treatment_years  # noqa: E402


def _basins():
    # Two 10 km x 10 km basins side by side, in a projected CRS (metres).
    return gpd.GeoDataFrame({"HYBAS_ID": [1, 2]}, geometry=[box(0, 0, 10_000, 10_000), box(10_000, 0, 20_000, 10_000)],
                            crs="EPSG:32717")


def test_basin_exposure_cumulative_and_dissolved():
    parcels = gpd.GeoDataFrame(
        {"year": [2009, 2011, 2011]},
        geometry=[box(0, 0, 5_000, 10_000),          # half of basin 1 in 2009
                  box(2_500, 0, 7_500, 10_000),      # overlaps the first parcel: only 2,500 m new
                  box(10_000, 0, 11_000, 10_000)],   # 10% of basin 2 in 2011
        crs="EPSG:32717")
    share = basin_exposure(parcels, _basins())
    assert share.loc[1, 2009] == pytest.approx(0.5)
    assert share.loc[1, 2011] == pytest.approx(0.75)  # overlap not double counted
    assert share.loc[2, 2009] == 0 and share.loc[2, 2011] == pytest.approx(0.1)
    assert treatment_years(share, 0.2) == {1: 2010}
    assert treatment_years(share, 0.05, lag=0) == {1: 2009, 2: 2011}
    assert treatment_years(share, 0.2, max_year=2009) == {}
    assert control_pool(share, [1, 2, 3], max_share=0.2) == [2, 3]


def test_forest_frame_loss_and_cover():
    wide = pd.DataFrame({"basin_id": [7], "forest2000": [0.8], "loss_2001": [0.08], "loss_2002": [0.0]})
    f = forest_frame(wide, [2001, 2002]).set_index("date")
    assert f.forest_loss_frac.tolist() == pytest.approx([0.1, 0.0])
    assert f.forest_cover_frac.tolist() == pytest.approx([0.72, 0.72])
    assert np.isnan(forest_frame(wide.assign(forest2000=0.0), [2001]).forest_loss_frac.iloc[0])
