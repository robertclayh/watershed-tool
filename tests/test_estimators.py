import numpy as np
import pandas as pd
import pytest

from wiat._linalg import project_simplex, simplex_lstsq
from wiat.did import did, event_study, placebo_interval, staggered_att
from wiat.synth import placebo_inverted_interval, synthetic_control


def factor_panel(n_units=25, n_periods=20, effect=2.0, start=12, seed=0, treated_loading=None, loading_sd=0.3):
    """y_it = a_i + l_t + lambda_i f_t + effect * D_it + e."""
    rng = np.random.default_rng(seed)
    a = rng.normal(0, 3, n_units)
    lt = np.cumsum(rng.normal(0, 1, n_periods))
    f = np.cumsum(rng.normal(0, 1, n_periods))
    lam = rng.normal(1, loading_sd, n_units)
    if treated_loading is not None:
        lam[0] = treated_loading
    y = a[:, None] + lt[None] + lam[:, None] * f[None] + rng.normal(0, 0.3, (n_units, n_periods))
    y[0, start:] += effect
    rows = [
        {"basin_id": f"U{i:02d}", "date": 2000 + t, "y": y[i, t], "x": rng.normal()}
        for i in range(n_units) for t in range(n_periods)
    ]
    return pd.DataFrame(rows), {"U00": 2000 + start}


def test_project_simplex():
    w = project_simplex(np.array([0.5, 2.0, -1.0]))
    assert np.isclose(w.sum(), 1) and (w >= 0).all()
    assert np.allclose(w, [0, 1, 0])


def test_simplex_lstsq_recovers_convex_combination():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(30, 8))
    w_true = np.array([0.5, 0.3, 0.2, 0, 0, 0, 0, 0])
    w = simplex_lstsq(A, A @ w_true)
    assert np.allclose(w, w_true, atol=1e-6)


def test_simplex_lstsq_rank_deficient():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(5, 40))
    w = simplex_lstsq(A, A[:, 3])
    assert np.isclose(w.sum(), 1) and (w >= -1e-12).all()
    assert np.linalg.norm(A @ w - A[:, 3]) < 1e-6


def test_did_recovers_effect_parallel_trends():
    df, ts = factor_panel(effect=2.0, seed=3, loading_sd=0.0)  # exact parallel trends
    r = did(df, "y", ts, alpha=0.1, compute=("cluster", "placebo"))
    assert abs(r.estimate - 2.0) < 0.8
    assert r.inference == "placebo"
    assert r.n_control == 24 and r.methods["placebo"]["n_placebo"] == 24


def test_did_matches_manual_2x2():
    # Balanced 2-period, no covariates: TWFE equals the difference of mean changes.
    df = pd.DataFrame({
        "basin_id": ["a", "a", "b", "b", "c", "c"],
        "date": [0, 1, 0, 1, 0, 1],
        "y": [1.0, 4.0, 2.0, 3.0, 0.0, 2.0],
    })
    r = did(df, "y", {"a": 1}, compute=("cluster",), inference="cluster")
    manual = (4 - 1) - np.mean([3 - 2, 2 - 0])
    assert np.isclose(r.estimate, manual)


def test_did_with_covariate_runs_and_wild_ci_contains_estimate():
    df, _ = factor_panel(n_units=12, seed=4)
    ts = {"U00": 2012, "U01": 2012, "U02": 2012}
    r = did(df, "y", ts, covariates=["x"], n_boot=199, inference="wild")
    lo, hi = r.methods["wild"]["ci_low"], r.methods["wild"]["ci_high"]
    assert lo <= r.estimate <= hi


def test_placebo_interval_resolution():
    # 9 placebos at alpha=0.1: min attainable p = 0.1; threshold is the largest |placebo|.
    lo, hi, p, thresh = placebo_interval(10.0, np.linspace(-1, 1, 9), alpha=0.1)
    assert np.isclose(p, 0.1) and thresh == 1.0 and (lo, hi) == (9.0, 11.0)
    # Interval excludes zero exactly when significant.
    for est in (0.5, 0.95, 1.05, 3.0):
        lo, hi, p, thresh = placebo_interval(est, np.linspace(-1, 1, 39), alpha=0.1)
        assert (p <= 0.1) == (not lo <= 0 <= hi)
    # Fewer than 1/alpha - 1 placebos: nothing is ever significant.
    lo, hi, p, thresh = placebo_interval(10.0, np.linspace(-1, 1, 5), alpha=0.1)
    assert np.isinf(thresh) and np.isinf(lo) and p > 0.1


def test_event_study_flat_pretrends_and_effect():
    df, ts = factor_panel(effect=3.0, seed=5, n_units=30, loading_sd=0.0)
    es = event_study(df, "y", ts, window=(-4, 4))
    post = es.coefs.loc[es.coefs.index >= 0, "coef"].mean()
    assert abs(post - 3.0) < 1.5
    assert -1 not in es.coefs.index
    assert "pl_low" in es.coefs


def test_scm_exact_fit_recovers_effect():
    rng = np.random.default_rng(6)
    T, J = 20, 10
    Y0 = np.cumsum(rng.normal(size=(T, J)), axis=0) + rng.normal(0, 5, J)
    w = np.zeros(J)
    w[[1, 4]] = [0.6, 0.4]
    y1 = Y0 @ w
    y1[14:] += 5.0
    rows = [{"basin_id": "T", "date": t, "y": y1[t]} for t in range(T)]
    rows += [{"basin_id": f"D{j}", "date": t, "y": Y0[t, j]} for j in range(J) for t in range(T)]
    res = synthetic_control(pd.DataFrame(rows), "y", "T", 14, method="standard", placebo_fit_cutoff=None)
    assert abs(res.att - 5.0) < 1e-3
    assert res.weights.loc["D1"] == pytest.approx(0.6, abs=1e-3)
    assert res.min_attainable_p == pytest.approx(1 / 11)
    assert res.p_value == pytest.approx(1 / 11)


@pytest.mark.parametrize("method", ["standard", "demeaned", "augmented"])
def test_scm_methods_run(method):
    df, ts = factor_panel(effect=2.0, seed=7)
    res = synthetic_control(df, "y", "U00", ts["U00"], method=method)
    assert np.isfinite(res.att)
    assert res.placebo is not None and len(res.placebo) == 24


def test_demeaned_scm_handles_level_shift():
    df, ts = factor_panel(effect=2.0, seed=8)
    df.loc[df.basin_id == "U00", "y"] += 100.0  # far outside donor levels
    std = synthetic_control(df, "y", "U00", ts["U00"], method="standard")
    dm = synthetic_control(df, "y", "U00", ts["U00"], method="demeaned")
    assert abs(dm.att - 2.0) < abs(std.att - 2.0)
    assert std.diagnostics["frac_pre_outside_donor_levels"] == 1.0


def test_placebo_inverted_interval():
    assert placebo_inverted_interval(1.0, 1.0, np.ones(5), 0.1) == (-np.inf, np.inf)
    lo, hi = placebo_inverted_interval(1.0, 0.5, np.linspace(0, 2, 19), 0.1)
    assert lo < 1.0 < hi and np.isclose(1.0 - lo, hi - 1.0)


def test_staggered_att():
    df, _ = factor_panel(effect=2.0, seed=9, n_units=30)
    ts = {"U00": 2012}
    # add a second cohort with the same effect
    m = (df.basin_id.isin(["U01", "U02"])) & (df.date >= 2008)
    df.loc[m, "y"] += 2.0
    ts.update({"U01": 2008, "U02": 2008})
    r = staggered_att(df, "y", ts, n_boot=99)
    assert abs(r.att - 2.0) < 1.0
    assert set(r.dynamic.index) >= {0, 1}
