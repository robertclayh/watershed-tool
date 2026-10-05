"""End-to-end walk-through on simulated basins with a known effect.

Two treated basins in the same simulated region:

1. a *typical* basin, with plenty of hydro-climatically similar donors;
2. an *atypical* basin (2 SD wetter, steeper, more degraded), the kind where
   nature-based investment is usually targeted.

The same intervention (NDVI +0.08 over four years, raising ET) is applied in both.
The toolkit attributes the first and refuses to attribute the second. It also shows
how the first study's conclusion would change if fewer donors were available.

    python examples/quickstart_synthetic.py            # writes to reports/quickstart
"""

from pathlib import Path

from wiat import AttributionStudy, aggregate, donor_thinning, simulate_watersheds
from wiat.report import _plt, plot_donor_thinning, write_report

OUT = Path("reports/quickstart")

for label, extremity, seed in (("typical", 0.0, 4), ("atypical", 2.0, 3)):
    sim = simulate_watersheds(n_basins=40, seed=seed, treated_extremity=extremity)
    truth = sim.true_att("et_mm")
    study = AttributionStudy(sim.panel, "et_mm", {"B000": "2013-01-01"}, static=sim.static)
    rep = study.run()
    write_report(rep, OUT / label)
    print(f"\n=== {label} basin: true ATT on annual ET = {truth:.1f} mm ===")
    print(f"verdict: {rep.verdict}  ({len(rep.donors)} donors)")
    print(rep.estimates()[["estimator", "estimate", "ci_low", "ci_high", "p_value"]].round(3).to_string(index=False))
    for f in rep.flags:
        print("  ", f)

    if label == "typical":
        # How would this study have fared with fewer donors available?
        annual = aggregate(sim.panel, "year")
        ranked = rep.selection.nearest(len(rep.selection.candidates))
        thin = donor_thinning(annual, "et_mm", "B000", 2013, ranked, sizes=(3, 5, 8, 12, 19, 29, 39),
                              covariates=["precip_mm"])
        plt = _plt()
        fig = plot_donor_thinning(thin, truth=truth)
        fig.savefig(OUT / "donor_thinning.png", bbox_inches="tight")
        plt.close(fig)
        print(thin[["k", "did_estimate", "did_ci_low", "did_ci_high", "scm_estimate", "scm_ci_low", "scm_ci_high"]]
              .round(1).to_string(index=False))

print(f"\nreports written to {OUT}/")
