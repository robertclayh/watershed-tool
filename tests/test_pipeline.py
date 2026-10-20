import numpy as np
import pandas as pd
import pytest

from wiat import AttributionStudy, aggregate, select_donors, simulate_watersheds
from wiat.cli import main
from wiat.matching import basin_signatures, hydrologically_connected
from wiat.panel import add_hydrologic_outcomes, balance, deseasonalize
from wiat.scarcity import _one_rep, donor_thinning, donors_needed, summarize_scarcity


@pytest.fixture(scope="module")
def sim():
    return simulate_watersheds(n_basins=40, seed=4)


def test_simulator_truth_is_counterfactual_difference(sim):
    eff = sim.true_effect("ndvi")
    assert (eff.loc[eff.index < 2013].abs() < 1e-12).all().all()
    assert eff.loc[2020, "B000"] == pytest.approx(0.08, abs=0.02)
    assert sim.true_att("water_yield_mm") == pytest.approx(-sim.true_att("et_mm"))


def test_aggregate_sums_fluxes_and_means_indices(sim):
    m = sim.panel[sim.panel.basin_id == "B001"]
    a = aggregate(m, "year")
    y = a.set_index("date").loc[2005]
    mm = m[pd.to_datetime(m.date).dt.year == 2005]
    assert y.precip_mm == pytest.approx(mm.precip_mm.sum())
    assert y.ndvi == pytest.approx(mm.ndvi.mean())
    assert y.evaporative_index == pytest.approx(mm.et_mm.sum() / mm.precip_mm.sum())


def test_aggregate_tolerates_index_gaps_but_not_flux_gaps():
    dates = pd.date_range("2001-01-01", periods=24, freq="MS")
    ndvi = np.r_[[np.nan] * 8, np.full(4, 0.5), [np.nan] * 9, np.full(3, 0.7)]
    precip = np.r_[np.full(12, 10.0), [np.nan], np.full(11, 10.0)]
    df = pd.DataFrame({"basin_id": "a", "date": dates, "ndvi": ndvi, "precip_mm": precip})
    a = aggregate(df, "year").set_index("date")
    assert a.loc[2001, "ndvi"] == pytest.approx(0.5) and np.isnan(a.loc[2002, "ndvi"])  # 4 vs 3 clear months
    assert a.loc[2001, "precip_mm"] == 120.0 and np.isnan(a.loc[2002, "precip_mm"])
    strict = aggregate(df, "year", min_months_mean=None).set_index("date")
    assert np.isnan(strict.loc[2001, "ndvi"])


def test_water_year_aggregation_labels_by_end_year():
    dates = pd.date_range("2001-01-01", periods=36, freq="MS")
    df = pd.DataFrame({"basin_id": "a", "date": dates, "precip_mm": 1.0})
    a = aggregate(df, "water_year", water_year_start_month=10)
    # Oct 2001 - Sep 2002 is water year 2002 and complete.
    assert a.set_index("date").loc[2002, "precip_mm"] == 12
    assert np.isnan(a.set_index("date").loc[2001, "precip_mm"])  # Jan-Sep 2001 only


def test_deseasonalize_and_balance(sim):
    d = deseasonalize(sim.panel, ["ndvi"], "2013-01-01")
    base = d[pd.to_datetime(d.date) < "2013-01-01"]
    assert abs(base.groupby("basin_id").ndvi_anom.mean()).max() < 1e-10
    p = sim.panel.copy()
    p.loc[p.index[5], "et_mm"] = np.nan
    bal, dropped = balance(p, ["et_mm"])
    assert dropped == ["B000"] and "B000" not in set(bal.basin_id)


def test_hydrologic_outcomes():
    df = pd.DataFrame({"precip_mm": [100.0, 0.0], "et_mm": [60.0, 5.0]})
    out = add_hydrologic_outcomes(df)
    assert out.water_yield_mm.tolist() == [40.0, -5.0]
    assert out.evaporative_index.iloc[0] == 0.6 and np.isnan(out.evaporative_index.iloc[1])


def test_matching_prefers_similar_basins(sim):
    feats = basin_signatures(sim.panel, "2013-01-01").join(sim.static)
    sel = select_donors(feats, ["B000"], caliper=1.2,
                        covariates=["precip_mm_mean", "aridity", "elevation_m", "forest_frac"],
                        exclude={"B001": "connected"})
    assert "B001" not in sel.candidates and sel.excluded["B001"] == "connected"
    assert 0 < sel.n_suitable < len(sel.candidates)
    near = sel.nearest(3)
    assert sel.distances.loc[near, "B000"].is_monotonic_increasing


def test_hydrologically_connected():
    t = pd.DataFrame({"HYBAS_ID": [1, 2, 3, 4, 5], "NEXT_DOWN": [3, 3, 4, 0, 0]})
    assert hydrologically_connected(t, 3) == {1, 2, 4}
    assert hydrologically_connected(t, 5) == set()


def test_attribution_study_typical_basin_attributes(sim):
    rep = AttributionStudy(sim.panel, "et_mm", {"B000": "2013-01-01"}, static=sim.static).run()
    assert rep.verdict == "ATTRIBUTED"
    truth = sim.true_att("et_mm")
    assert rep.did.ci[0] <= truth <= rep.did.ci[1]
    assert rep.scm["B000"].ci[0] <= truth <= rep.scm["B000"].ci[1]
    assert "Verdict: ATTRIBUTED" in rep.to_markdown()


def test_attribution_study_flags_scarce_and_extreme():
    s = simulate_watersheds(n_basins=40, seed=3, treated_extremity=2.0)
    rep = AttributionStudy(s.panel, "et_mm", {"B000": "2013-01-01"}, static=s.static).run()
    codes = {f.code for f in rep.flags}
    assert rep.verdict == "NOT ATTRIBUTABLE (DESIGN)"
    assert {"DONOR_SCARCITY", "OUTSIDE_CONVEX_HULL", "POOR_PRE_FIT"} <= codes


def test_midyear_start_maps_to_next_full_year(sim):
    rep = AttributionStudy(sim.panel, "ndvi", {"B000": "2012-06-15"}, static=sim.static).run()
    assert rep.treat_start == {"B000": 2013}


def test_donors_needed():
    assert donors_needed(0.1) == {"significance_possible": 9, "stable_threshold": 19}
    assert donors_needed(0.05)["significance_possible"] == 19


def test_scarcity_rep_and_summary():
    rows = []
    for rep in range(2):
        rows += _one_rep((4, 0.0, rep, 1, "et_mm", 0.1, 16, 9, 0.08, ("precip_mm",), ("did_placebo", "scm")))
    df = pd.DataFrame(rows)
    assert set(df.world) == {"effect", "no_effect"}
    s = summarize_scarcity(df, alpha=0.1)
    # 4 donors at alpha=0.1: placebo inference can neither bound nor reject.
    assert (s.ci_bounded == 0).all() and (s.power == 0).all() and (s.confident_correct == 0).all()


def test_donor_thinning(sim):
    annual = aggregate(sim.panel, "year")
    ranked = [b for b in annual.basin_id.unique() if b != "B000"]
    out = donor_thinning(annual, "et_mm", "B000", 2013, ranked, sizes=(4, 9, 19))
    assert list(out.k) == [4, 9, 19]
    assert np.isinf(out.loc[out.k == 4, ["did_ci_low", "scm_ci_low"]]).all().all()
    # DiD placebo intervals become bounded at 9 donors. SCM may need more, because
    # placebos with poor pre-fit are discarded (Abadie et al.).
    assert np.isfinite(out.loc[out.k == 9, "did_ci_low"]).all()
    assert np.isfinite(out.loc[out.k == 19, ["did_ci_low", "scm_ci_low"]]).all().all()


def test_cli_simulate_and_analyze(tmp_path):
    assert main(["simulate", "--out", str(tmp_path / "sim"), "--basins", "25"]) == 0
    assert main([
        "analyze", "--panel", str(tmp_path / "sim/panel.csv"),
        "--treatments", str(tmp_path / "sim/treatments.csv"),
        "--static", str(tmp_path / "sim/static.csv"), "--out", str(tmp_path / "rep"),
    ]) == 0
    assert (tmp_path / "rep/report.md").exists()
    assert (tmp_path / "rep/scm_B000.png").exists()


def test_injection_experiment_recovers_effect_and_holds_size(sim):
    from wiat.validation import injection_experiment, prepare_region, summarize_injection

    region = prepare_region(sim.counterfactual, "et_mm", {"B000": "2013-01-01"}, static=sim.static)
    assert "B000" not in region["eligible"]
    r = injection_experiment(region, donor_counts=(4, 12), effects=(0.10,), n_basins=12, n_jobs=1)
    assert set(r.world) == {"effect", "no_effect"}
    eff = r[(r.world == "effect") & (r.estimator == "did_placebo")]
    assert (eff.truth > 0).all()
    s = summarize_injection(r.assign(extremity=-1.0), 0.10)
    did4 = s[(s.estimator == "did_placebo") & (s.n_donors == 4)].iloc[0]
    assert did4.ci_bounded == 0 and did4.power == 0  # below 9 donors nothing is attributable
    did12 = s[(s.estimator == "did_placebo") & (s.n_donors == 12)].iloc[0]
    assert did12.ci_bounded == 1 and did12.coverage >= 0.75


def test_scm_overfit_flag_with_short_baseline():
    s = simulate_watersheds(n_basins=40, n_years=12, treat_start={"B000": 2004}, seed=4)
    rep = AttributionStudy(s.panel, "et_mm", {"B000": "2004-01-01"}, static=s.static, caliper=np.inf,
                           hard_caliper=np.inf).run()
    codes = [f.code for f in rep.flags]
    assert len(rep.donors) > 9 and "SCM_OVERFIT" in codes
    assert not np.isfinite(rep.scm["B000"].ci).all()
