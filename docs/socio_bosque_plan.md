# Study plan: Socio Bosque and catchment hydrology (Ecuador)

Status: **done**; results in [socio_bosque_results.md](socio_bosque_results.md),
including the spillover and threshold checks.

## Question

Did enrolling forest in Socio Bosque change catchment-scale vegetation, forest loss,
evapotranspiration and the water balance, relative to comparable never-enrolled
catchments? And how does the answer depend on the share of each catchment enrolled?

This differs from the three single-basin studies in
[latam_results.md](latam_results.md). There, one treated basin per programme meant
placebo inference, at most a few dozen donors, and the "below 9 donors" cliff. Here,
hundreds of catchments enter treatment in different years. That allows a staggered
design (Callaway–Sant'Anna) with basin-bootstrap inference, which is valid with many
treated units, and a dose-response test.

## Data

| | Source | Notes |
|---|---|---|
| Treatment | MAATE Socio Bosque agreements (WFS, `docs/intervention_data_sources.md`) | 2,569 active agreements, 18,300 km², signing date per agreement |
| Units | HydroBASINS level 12 within Ecuador | 2,032 basins, median 138 km². Levels 10–12 are almost identical in size here, so level 12 is the finest available. |
| Panel | Earth Engine, **annual** periods 2001–2024 | Annual composites (12× fewer requests than monthly, and an annual median NDVI is far more robust to Andean and Amazonian cloud than monthly values) |

## Unit of analysis and treatment

**Unit:** HydroBASINS level-12 basin. **Exposure:** the share of basin area under
Socio Bosque agreement, cumulative by signing year (`basin_share_by_year`).

Measured on the real overlay:

| Basin share enrolled reaches | Basins | First-crossing years |
|---|---|---|
| any | 985 | |
| ≥ 10% | 280 | 2008: 13, 2009: 40, 2010: 32, 2011: 60, 2012: 27, 2013: 11, 2014: 38, later: 56 |
| ≥ 20% | 201 | 2009: 28, 2010: 25, 2011: 45, 2012: 21, 2014: 27, … |
| ≥ 50% | 97 | |
| 0 (never) | 1,047 | |

Collective agreements (communities and Indigenous territories) make up 89% of the
enrolled area within basins (16,300 km²), individual agreements 8%, and mangrove
agreements 3%.

- **Treated:** share ≥ 20% (main specification). Treatment year = the first full year
  after the share crosses the threshold (signing year + 1). Sensitivity: 10% and 50%.
- **Cohorts:** 2009–2015 for the main analysis. Cohorts after 2015 have too few
  basins and too short a post-period, so they are dropped from treatment *and*
  from controls.
- **Dose:** share bins (20–50%, ≥ 50%) estimated separately, plus a continuous
  share interaction as a check.
- **Excluded from controls:** any basin with share > 1% at any time (partially
  treated) and mangrove agreements' coastal basins, unless they are matched as their
  own stratum.
- **Outside the study population:** 37 Galápagos basins (inside Ecuador's boundary,
  but with no Socio Bosque, and CHIRPS reports ~100 mm/yr there against 1,660 mm of
  MODIS ET) and 68 coastal fragments under 10 km².

**Study population as built:** 1,927 basins; **162 treated** (share ≥ 20%, first
treated 2008–2016: 2009: 7, 2010: 28, 2011: 25, 2012: 45, 2013: 21, 2015: 26,
2016: 7, others: 3); **1,278 never-enrolled controls** (share ≤ 1%).

## Outcomes

In order of the causal chain (each is a check on the next):

1. **Forest loss** (Hansen GFC v1.13): annual loss as a share of the basin's 2000
   forest, plus cumulative forest cover. This is the first stage: Socio Bosque pays
   to keep forest standing, so if loss doesn't fall, hydrology has no channel to change.
2. **NDVI / EVI**: annual median composites (Landsat 5/7/8/9 + Sentinel-2).
3. **ET** (MOD16) and **ET/P**.
4. **P − ET** as the annual water-yield proxy.

Precipitation (CHIRPS) is a covariate, not an outcome.

## Control selection

The forest panel already shows why matching is essential. Before enrollment,
treated basins lost forest **about 6½× more slowly** than never-enrolled ones:

| Pooled annual forest loss (% of 2000 forest) | 2001–07 | 2008–15 | 2016–24 |
|---|---|---|---|
| Never-enrolled (1,278) | 0.26 | 0.33 | 0.33 |
| Treated (162) | 0.04 | 0.06 | 0.12 |

Treated basins are flatter (median slope 1.5° vs 5.4°) and almost fully forested
(median 100% vs 93%): remote lowland forest under collective title. An unmatched
comparison would credit Socio Bosque with differences that existed before it. The
treated basins' loss also *rises* over time, which a naive DiD would read as the
programme increasing deforestation. Matching on the pre-period level and trend of
loss, and the event study, have to sort this out.

Never-enrolled basins (share ≤ 1% throughout) are matched to treated basins on
pre-2008 characteristics:

- **Region stratum:** Coast / Andes / Amazon. Exact match, because the regions differ
  in climate, deforestation drivers and enrollment (collective agreements are
  concentrated in Amazonian Indigenous territories).
- **Baseline forest cover** (2000) and **pre-period forest-loss rate and trend**
  (2001–07). Deforestation pressure is the main selection channel: Socio Bosque
  enrolls remote, low-threat forest, which would lose little forest anyway.
- Elevation, slope, mean precipitation, baseline ET and NDVI.

Matching uses the existing `select_donors` signatures (Mahalanobis distance with
calipers) at the pooled level, reporting balance tables before and after. The pool is
large (about 1,000 basins), so the "donor scarcity" failure should not apply. That is
itself worth stating against the single-basin results.

## Estimation and inference

- `staggered_att` (Callaway–Sant'Anna, never-treated matched controls, universal base
  period, basin bootstrap stratified by cohort) for each outcome.
- Event-study ATT(e) by years since enrollment, from −6 to +10, as the pre-trend check.
- TWFE DiD only as a comparison, flagged as biased under staggered adoption.
- **Validation on the same panel:** the injected-effect experiment from
  `wiat.validation`, adapted to many pseudo-treated basins drawn at random from the
  never-enrolled pool with the real cohort-size distribution. It reports size,
  coverage and the minimum detectable effect for each outcome.

## Extraction size

| Item | Size |
|---|---|
| Basins | 1,927 study basins (162 treated, 1,278 controls, the rest partially enrolled and kept for dose checks) |
| Periods | 24 annual (2001–2024) |
| Panel rows | ~49,000 |
| Earth Engine route | `--export` to Google Drive: six batch tasks (4-year blocks), run in parallel on Google's servers. Interactive `getInfo` is used only for the pilot (40 basins × 6 years in 3.5 minutes; no NDVI gaps after 2019 with annual medians). |
| Forest panel | Done: Hansen loss by year per basin for all 2,032 basins (one reduction, ~1 minute) |

## Steps

1. ✅ Download agreements (`wiat.data.interventions.fetch_wfs`).
2. ✅ Overlay with level-12 basins: exposure by year (`basin_exposure`).
3. ✅ Annual extraction and Hansen forest-loss panel in `wiat.data.gee`.
4. ✅ Pilot on a small subset of basins; Galápagos and fragments excluded as a result.
5. ✅ Full extraction: Drive export tasks completed. The CSVs were fetched with Earth
   Engine's direct table download (`getDownloadURL`, one year per request, ~25 s
   each), which runs the same computation without a Drive round trip. Then
   `examples/socio_bosque.py assemble`.
6. ✅ Region strata, matching and balance tables (k = 1 main; k = 2, 3, 5 sensitivity).
7. ✅ Staggered ATT, event study and dose response for each outcome.
8. ✅ Placebo validation: random enrollment among never-enrolled Amazon basins.
9. ✅ Write-up: [socio_bosque_results.md](socio_bosque_results.md).
10. ✅ Spillover checks (buffer, strict buffer, leakage) and 10% / 50% thresholds
    (`socio_bosque_analysis.py --checks`). They overturned the first-stage finding: see
    results §2 and §6.
11. ⬜ Parcel-level first stage: agreement polygons against matched non-enrolled
    forest pixels.

## Known threats

- **Survivorship:** only active agreements are published. Basins whose agreements
  ended (for example after deforestation) are missing from the treated group, which
  biases toward finding avoided loss.
- **Selection on low threat:** handled by matching on pre-period loss, but the
  residual risk remains. The event study's pre-period is the check.
- **Spillovers / leakage:** deforestation displaced to neighbouring non-enrolled
  basins would *inflate* the estimate. A check uses controls at least one basin away
  from any treated basin.
- **MODIS ET at 500 m** is coarse for 138 km² basins with partial enrollment, but
  adequate at basin scale. It is still model output that depends on vegetation
  indices.
- **Cloud:** annual medians reduce gaps but do not remove them in the wettest Andean
  slopes. `veg_clear_frac` is carried through, and basins with no clear observation
  in a year get NaN.
