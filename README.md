# Watershed Intervention Attribution Toolkit (`wiat`)

An open-source toolkit for **attributing hydrologic change to nature-based
interventions** (riparian restoration, páramo protection, reforestation, payments
for ecosystem services) using remote sensing and causal inference. It also
**documents where that attribution fails**.

- **Data:** Google Earth Engine extraction of basin panels. Landsat 5/7/8/9 and
  Sentinel-2 vegetation indices (harmonised), CHIRPS precipitation, MODIS MOD16
  evapotranspiration, NASADEM terrain, Hansen forest cover, and HydroBASINS units and
  topology.
- **Designs:** difference-in-differences (TWFE, event study, Callaway–Sant'Anna for
  staggered adoption) and synthetic control (standard, de-meaned, ridge-augmented),
  with donor basins chosen by calipered matching on hydro-climatic signatures.
- **Inference for the realistic case of one or a few treated basins:** placebo
  permutation (Conley–Taber; Abadie et al.), test-inverted confidence intervals, a
  detection threshold, plus wild-bootstrap and cluster-robust results for comparison.
- **Failure diagnostics:** every study ends in an explicit verdict: `ATTRIBUTED`,
  `NO DETECTABLE EFFECT`, `INCONCLUSIVE` or `NOT ATTRIBUTABLE (DESIGN)`, with the
  flags that produced it.
- **Scarcity analysis:** a Monte Carlo and an empirical donor-thinning tool that
  quantify how confidence degrades as suitable control basins become scarce.

## Key result

Attribution confidence does not degrade gracefully. It falls off a cliff, and
conventional methods hide the cliff.

![Attribution quality vs donor scarcity](docs/figures/scarcity_et_mm_scm.png)

- **Fewer than 9 suitable donor basins (α = 0.10): attribution is impossible by
  construction.** No effect of any size can reach significance, and intervals are
  unbounded.
- **Conventional cluster-robust DiD looks fine but is wrong.** With one treated basin,
  its false-positive rate is 44–78% where 10% is promised, and worse with *more*
  donors.
- **Typical basins with ≥ 9 donors:** honest methods keep their nominal error rate,
  but a real 5–9% ET change is confidently attributed only 38–58% of the time.
- **Atypical basins, the ones where investment is usually targeted, are the danger
  zone.** Adding poorly matched donors narrows intervals around a biased estimate.
  At 2 SD atypicality and 49 donors, SCM's false-positive rate reaches 40%, even
  though its pre-period fit is the best in the experiment.

Full write-up: [docs/where_attribution_fails.md](docs/where_attribution_fails.md).

### On real basins (Extrema and Camboriú, Brazil)

The same experiment on 293 real HydroBASINS units confirms the simulation:

- **Camboriú: a bounded null.** The programme did not change basin ET by more than
  about 10 mm/yr (≈ 0.9%) or NDVI by more than about 0.03.
- **Extrema: not attributable.** The programme basin is more atypical than 96% of
  its region (elevation +2.7 SD, forest +2.3 SD), so it has 6 suitable donors. A
  typical basin in the same region has 92.
- **Injecting known effects into every untreated basin:** cluster-robust DiD gives
  38–83% false positives, and placebo methods hold 10%. Nothing is detectable below
  9 donors. With 9 or more, a 5% ET change is found 87–100% of the time.
- **Synthetic control with a 6-year baseline loses all power as donors are added**,
  because placebos fit the short pre-period exactly.

Write-up: [docs/latam_results.md](docs/latam_results.md).
Methods: [docs/methodology.md](docs/methodology.md).

## Install

```bash
pip install -e .            # core (numpy, pandas, scipy, matplotlib)
pip install -e ".[gee]"     # + Google Earth Engine extraction
pip install -e ".[dev]"     # + pytest
```

Python ≥ 3.11.

## Quick start (no Earth Engine account needed)

```bash
python examples/quickstart_synthetic.py
```

The script simulates a region of 40 basins with a known intervention effect and runs
two studies:

```
=== typical basin: true ATT on annual ET = 55.6 mm ===
verdict: ATTRIBUTED  (21 donors)
     estimator  estimate  ci_low  ci_high  p_value
 DiD (placebo)    50.168   0.508   99.828    0.091
SCM (demeaned)    56.850  25.538   88.163    0.050

=== atypical basin: true ATT on annual ET = 61.6 mm ===
verdict: NOT ATTRIBUTABLE (DESIGN)  (5 donors)
   [FATAL] DONOR_SCARCITY: 5 usable donor basin(s). Placebo inference can never reach p <= 0.1 ...
   [MAJOR] POOR_PRE_FIT: B000: synthetic control misses the pre-period by 1.13 SD ...
   [MAJOR] OUTSIDE_CONVEX_HULL: B000: 75% of pre-periods lie outside the donor range ...
```

### Python API

```python
from wiat import AttributionStudy, simulate_watersheds

sim = simulate_watersheds(n_basins=40, seed=4)          # or a panel from Earth Engine
study = AttributionStudy(
    sim.panel,                       # long monthly panel: basin_id, date, precip_mm, et_mm, ndvi, ...
    outcome="et_mm",                 # or water_yield_mm, evaporative_index, ndvi, ...
    treat_start={"B000": "2013-01-01"},
    static=sim.static,               # elevation, slope, forest cover, ... for matching
    alpha=0.1,
)
report = study.run()
print(report.verdict)
for flag in report.flags:
    print(flag)
print(report.estimates())

from wiat.report import write_report
write_report(report, "reports/my_study")   # markdown + synthetic-control and event-study figures
```

The building blocks are usable on their own: `did`, `event_study`, `staggered_att`,
`synthetic_control`, `select_donors`, `donor_thinning`, `run_scarcity_experiment`.

## Real basins with Google Earth Engine

```bash
earthengine authenticate
wiat extract --config examples/latam_water_funds.toml --project <ee-project> --out data
wiat analyze --panel data/extrema/panel.csv --treatments data/extrema/treatments.csv \
             --static data/extrema/static.csv --exclusions data/extrema/exclusions.csv \
             --outcome et_mm --out reports/extrema
```

`examples/latam_water_funds.toml` defines illustrative studies for Latin American
watershed-investment programmes: Extrema (Brazil), Camboriú (Brazil) and Quito's
FONAG páramo catchments (Ecuador). For each, the extractor:

1. finds the HydroBASINS unit containing each intervention point;
2. excludes basins upstream or downstream of it (spillovers) and any listed as
   contaminated;
3. builds a monthly panel of NDVI/EVI/NDMI (cloud-masked, Landsat–Sentinel
   harmonised), CHIRPS precipitation, MODIS ET/PET and a cloud-free coverage fraction;
4. extracts static covariates (elevation, slope, tree cover, forest loss).

Use `--export` for regional extractions, which runs as a batch export to Drive
instead of interactive `getInfo`. The coordinates and start years in the example
config are approximate reference points. Replace them with actual intervention
footprints and the year land management changed before drawing conclusions.

## Reproducing the scarcity analysis

```bash
wiat scarcity --reps 200 --out docs/results --fig-dir docs/figures   # ~7 min on 4 cores
```

## Repository layout

```
wiat/
  data/gee.py      Earth Engine extraction (Landsat/S2, CHIRPS, MODIS ET, HydroBASINS, static)
  panel.py         validation, aggregation (calendar/water year), hydrologic outcomes
  matching.py      basin signatures, calipered donor selection, balance, connectivity
  did.py           TWFE DiD (cluster / wild bootstrap / placebo), event study, Callaway-Sant'Anna
  synth.py         synthetic control (standard / de-meaned / augmented), placebo & conformal inference
  attribution.py   AttributionStudy: pipeline, failure flags, verdict
  scarcity.py      Monte Carlo over donor count x atypicality; empirical donor thinning
  validation.py    injected-effect validation on real panels
  simulate.py      synthetic watershed panels with known counterfactuals
  report.py        figures and markdown reports
  cli.py           `wiat simulate | extract | analyze | scarcity`
docs/              methodology, failure analysis, figures, raw Monte Carlo results
examples/          synthetic quick start, Latin American Earth Engine config, latam_results.py
tests/             pytest suite
```

## Testing

```bash
pytest -q
```

The tests check estimator correctness on data with known answers: exact
synthetic-control recovery, the manual 2×2 DiD, Callaway–Sant'Anna on staggered
cohorts, and placebo-interval resolution limits. They also cover the matching and
aggregation logic, end-to-end verdicts on simulated typical and atypical basins, and
the CLI. The Earth Engine module is tested for its pure helpers. Live extraction
needs an authenticated account and is not run in CI.

## License

MIT
