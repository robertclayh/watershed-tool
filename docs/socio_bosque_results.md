# Socio Bosque and catchment hydrology: results

The design is in [socio_bosque_plan.md](socio_bosque_plan.md). In short: 1,927 HydroBASINS
level-12 basins in mainland Ecuador. A basin is treated from the first full year after
≥ 20% of its area is under a Socio Bosque agreement (162 basins, cohorts 2008–2016).
Each treated basin is matched to its nearest never-enrolled basin (≤ 1% enrolled), within
its region (Amazon / Andes / Coast), on 2001–07 forest cover, forest-loss level and
trend, NDVI, ET, precipitation, elevation and slope. Effects are Callaway–Sant'Anna ATTs
with a basin bootstrap stratified by cohort (499 draws), reported as 90% intervals.
Annual panels run 2001–2024.

Reproduce with `examples/socio_bosque.py` (prepare / extract / assemble) and
`examples/socio_bosque_analysis.py` (add `--checks` for §6), and
`examples/socio_bosque_parcels.py` for §7. Results are in
`results/socio_bosque_*`.

## Summary

1. **Matching is essential, and even matching is not enough.** Treated basins are
   remote, flat, almost fully forested Amazonian land under collective title. Before
   enrollment they lost forest about 6½× more slowly than never-enrolled basins.
   Unmatched comparisons say Socio Bosque *increased* deforestation and cut rainfall by
   163 mm/yr, both artefacts of selection. After matching, Amazon controls still lose
   about twice as much forest at baseline as the treated basins.
2. **No robust evidence of avoided deforestation.** The main matched estimate for the
   Amazon is −0.038 percentage points of forest per year (90% CI −0.086 to +0.002). It
   is significant only in specifications with looser balance, and it shrinks toward
   zero as controls are matched more closely on baseline deforestation:
   - controls at least one basin away from any enrolled land: −0.018
   - 10% threshold: −0.015
   - 50% threshold: −0.001, a tight null
   - log-type scale, best-balanced match: −0.037 [−0.105, +0.027]

   That pattern is what residual selection on deforestation pressure produces.
3. **No detectable hydrological effect.** NDVI, ET, ET/P and P − ET show no change in
   the Amazon in any design. ET is bounded within a few mm/yr (about ±0.4%). Even the
   largest forest-loss estimate would save too little forest to move basin ET beyond
   what the design can resolve.
4. **Pooled across regions, a rainfall placebo fails in every design**, so basin
   hydrology estimates outside the Amazon (27 treated basins) are confounded by climate.
5. **No detectable leakage** of deforestation into neighbouring basins, but the test
   is imprecise.
6. **At the pixel level, inside agreements (§7), the same story holds with sharper
   evidence.** Matching on travel time, forest-edge distance, terrain and local
   deforestation pressure balances every covariate to within 0.02 SD. Yet treated
   pixels were already losing forest more slowly than their matches *before*
   enrollment: 86% more slowly for individual agreements. Net of that pre-existing
   gap, avoided loss is −0.017 pp/yr [−0.043, +0.004] on enrolled forest, not
   significant on either the difference or the ratio scale. The binding constraint is
   that forest loss is rare: 22 losses in treated pixels against 65 in controls.

This is the toolkit's central finding in a design that should favour attribution:
135 treated Amazonian basins, hundreds of candidate controls, and a 7–15-year baseline,
with hundreds of treated units making the inference close to nominal. Precision is no
longer the problem. Selection is: the programme enrolled forest that was under little
threat, and no observable comparison basin fully reproduces that. Where the comparison
is cleanest (nearly fully enrolled basins against controls with the same near-zero
baseline loss), there is no effect to attribute, on forest or on water.

## 1. Selection and balance

![Overall ATT by specification](figures/socio_bosque_att.png)

| Standardised mean difference, treated vs controls | Unmatched (1,278) | Matched (92 controls, k = 1) |
|---|---|---|
| Forest cover 2001 | 1.05 | 0.12 |
| Baseline forest-loss rate (asinh) | −1.31 | −0.23 |
| Baseline loss trend | −0.25 | −0.02 |
| Baseline NDVI | 1.20 | 0.13 |
| Baseline ET | 1.11 | 0.19 |
| Baseline precipitation | 1.07 | 0.12 |
| Elevation | −0.57 | −0.13 |
| Slope | −0.68 | −0.15 |

With one control per treated basin, every covariate is within the usual 0.25 threshold.
Using more controls per basin (k = 2, 3, 5) narrows intervals but loosens balance
(maximum |SMD| 0.30, 0.38, 0.43).

The unmatched comparison shows what selection does:

| Unmatched ATT | Estimate [90% CI] |
|---|---|
| Forest loss (pp of 2000 forest / yr) | **+0.080** [+0.054, +0.109] |
| Precipitation, placebo (mm/yr) | **−163** [−209, −123] |
| ET (mm/yr) | +22.8 [+19.5, +26.1] |

A forest-protection programme cannot cut rainfall by 163 mm/yr. The unmatched controls
are drier, more deforested basins on different trajectories.

An SMD of −0.23 hides a large *proportional* gap. In the Amazon, the treated basins'
2001–07 loss averaged 0.031 %/yr and their matched controls' 0.065 %/yr, about twice as
much. §2 shows why that matters.

## 2. Forest loss: an apparent effect that does not survive tighter matching

![Event study, Amazon](figures/socio_bosque_event_study.png)

The Amazon event study (k = 1) looks like avoided deforestation. Loss rates are
indistinguishable before enrollment (leads within ±0.02 pp, one of five excluding
zero). They stay near zero for three years, then run 0.05–0.10 pp/yr lower in most
years from year 4 on. Across matching specifications:

| Amazon forest loss (pp of 2000 forest / yr) | Estimate [90% CI] | Max \|SMD\| | Baseline loss, treated vs controls (%/yr) |
|---|---|---|---|
| k = 5 | −0.061 [−0.093, −0.030] | 0.43 | 0.031 vs 0.082 |
| k = 3 | −0.056 [−0.096, −0.020] | 0.38 | 0.031 vs 0.082 |
| k = 2 | −0.053 [−0.095, −0.011] | 0.30 | 0.031 vs 0.072 |
| **k = 1 (main)** | **−0.038 [−0.086, +0.002]** | 0.23 | 0.031 vs 0.065 |
| Strict buffer (§6) | −0.018 [−0.065, +0.024] | 0.32 | 0.031 vs 0.047 |
| 10% threshold (§6) | −0.015 [−0.055, +0.024] | 0.16 | 0.034 vs 0.063 |
| 50% threshold (§6) | −0.001 [−0.010, +0.009] | 0.15 | 0.005 vs 0.006 |

**The estimate tracks the baseline gap.** The more the controls out-deforest the treated
basins before enrollment, the larger the apparent effect. Where baseline loss is equal
(the 50% design), the effect is zero with a tight interval.

**The likely mechanism is proportional growth.** After 2010, Amazonian deforestation
roughly tripled in both groups: treated basins went from 0.031 to 0.094 %/yr, and
matched controls from 0.065 to 0.182 %/yr. If clearing grows in proportion to its
existing level, a difference-in-differences in absolute percentage points reads that
common growth as a programme effect. On a log-type scale (asinh of loss), where
proportional growth counts as no effect, the estimate is not significant in any
well-balanced design:

| Amazon, asinh(10 × loss %) | Estimate [90% CI] |
|---|---|
| k = 1 | −0.037 [−0.105, +0.027] |
| k = 3 | −0.057 [−0.113, −0.003] |
| k = 5 | −0.073 [−0.119, −0.021] |
| Strict buffer | +0.014 [−0.044, +0.077] |
| 10% threshold | −0.006 [−0.068, +0.054] |

**The earlier dose result was an artefact.** In the main run, basins ≥ 50% enrolled
showed −0.047 [−0.086, −0.011] against the controls matched to *all* treated basins.
Those nearly fully enrolled basins lose almost no forest (0.006 %/yr at baseline),
far less than that shared control set. Matched to their own controls (the 50% threshold
design), the effect disappears.

Cumulative forest cover adds nothing independent: it integrates the same rate
difference, and its leads drift by about 0.06 pp/yr before enrollment.

**Conclusion on the first stage:** the data cannot distinguish a small avoided-loss
effect from residual selection on deforestation pressure. Survivorship (only active
agreements are published) would bias *toward* finding an effect, and the best-matched
comparisons still find none.

## 3. Hydrology: bounded nulls in the Amazon

| Amazon, matched (k = 1) | ATT [90% CI] | Treated baseline mean | Placebo-validation half-width |
|---|---|---|---|
| NDVI | −0.0015 [−0.0047, +0.0019] | 0.81 | ±0.0032 |
| ET (mm/yr) | +1.8 [−2.1, +5.7] | 1,480 | ±3.5 |
| ET / P | +0.0025 [−0.0021, +0.0070] | 0.46 | ±0.0041 |
| P − ET (mm/yr) | −20 [−82, +39] | 1,827 | ±45 |
| Precipitation, placebo (mm/yr) | −18 [−80, +42] ✓ | 3,306 | ±45 |

Nothing moves, in this or any robustness design (§6: Amazon ET between +0.7 and +2.8
mm/yr, every interval including zero). The precipitation placebo passes in the Amazon,
borderline under the strict buffer (−59 [−115, +1]). ET is tightly bounded: the main
interval excludes changes larger than about 6 mm/yr (0.4%). The ET event study has
pre-period leads of −2 to −7 mm, so the comparison is not clean at the ±5 mm level, but
no specification shows a hydrologically meaningful change.

**Even the largest forest effect would be invisible in hydrology.** Take the largest
estimate (k = 5): its cumulative cover gap reaches about 1% of a 138 km² basin after a
decade. If tropical forest converted to pasture loses 20–30% of its ET (an illustrative
figure, not one measured here), that moves basin ET by 0.2–0.3%, or about 3–4 mm/yr.
That is at the design's ±3.5 mm resolution. P − ET is noisier still: its ±45 mm
resolution is set by rainfall variability.

## 4. Validation: what this design can see

The validation runs 200 times. Each time, 135 never-enrolled Amazonian basins are
assigned the real cohort years at random, and the rest serve as controls. No basin is
treated, so every rejection is a false positive.

| Outcome | False-positive rate (target 10%) | Median 90% half-width | SD of placebo estimates |
|---|---|---|---|
| Forest loss (pp/yr) | 13% | 0.069 | 0.046 |
| Forest cover (pp) | 14% | 0.77 | 0.51 |
| NDVI | 14% | 0.0032 | 0.0019 |
| ET (mm/yr) | 17% | 3.5 | 2.5 |
| ET / P | 12% | 0.0041 | 0.0026 |
| P − ET (mm/yr) | 13% | 45 | 28 |
| Precipitation (mm/yr) | 13% | 45 | 28 |

Random assignment among never-enrolled basins removes selection, so this measures noise
only. The intervals run slightly narrow (12–17% false positives at a nominal 10%), and
the resolution is fine: ±3.5 mm/yr for ET, against ±12–30 mm/yr for one basin with a
full donor pool ([latam_results.md](latam_results.md)). What validation by random
assignment cannot measure is the selection bias in §2. Real enrollment was not random,
and that, not noise, limits this study.

## 5. Outside the Amazon

Pooled across regions (162 treated), the precipitation placebo gives −46 mm/yr
[−99, +5] at k = 1 and excludes zero at k = 2, 3 and 5 and in every §6 design except the
50% threshold. Treated Andean and Coastal basins (27) follow different rainfall
trajectories from their matches, so pooled ET (+8.5 [+3.6, +13.5] mm/yr) and P − ET
estimates reflect climate, not the programme. With 14 treated Andean and 13 Coastal
basins, separate regional estimates would be too noisy to be useful.

## 6. Robustness: spillovers and enrollment thresholds

Each design is re-matched at k = 1 (`socio_bosque_analysis.py --checks`;
`results/socio_bosque_robustness.csv`). HydroBASINS units are adjacent if their
boundaries are within 100 m (median 6 neighbours).

| Design | Treated (Amazon) | Matched controls (Amazon) | Max \|SMD\| |
|---|---|---|---|
| **Buffer:** controls not adjacent to any treated basin | 162 (135) | 78 (56) | 0.34 |
| **Strict buffer:** controls not adjacent to any basin > 1% enrolled | 162 (135) | 71 (52) | 0.32 |
| **Leakage:** never-enrolled neighbours of treated basins as "treated" | 130 (98) | 75 (55) | 0.12 |
| **10% threshold** | 231 (164) | 137 (89) | 0.16 |
| **50% threshold** | 78 (74) | 45 (41) | 0.15 |

| Amazon ATT [90% CI] | Forest loss (pp/yr) | NDVI | ET (mm/yr) | Precipitation, placebo (mm/yr) |
|---|---|---|---|---|
| Main (k = 1) | −0.038 [−0.086, +0.002] | −0.0015 [−0.0047, +0.0019] | +1.8 [−2.1, +5.7] | −18 [−80, +42] |
| Buffer | −0.036 [−0.080, +0.003] | −0.0025 [−0.0053, +0.0008] | +2.7 [−0.9, +7.0] | −36 [−98, +19] |
| Strict buffer | −0.018 [−0.065, +0.024] | −0.0035 [−0.0064, −0.0001] | +1.3 [−2.6, +5.6] | −59 [−115, +1] |
| Leakage | +0.025 [−0.064, +0.110] | −0.0054 [−0.0099, −0.0007] | +0.8 [−3.6, +4.9] | −15 [−93, +50] |
| 10% threshold | −0.015 [−0.055, +0.024] | −0.0019 [−0.0047, +0.0008] | +2.8 [−0.4, +6.2] | −34 [−91, +23] |
| 50% threshold | −0.001 [−0.010, +0.009] | −0.0027 [−0.0064, +0.0005] | +0.7 [−3.5, +5.5] | −32 [−110, +46] |

- **Spillovers do not explain the result.** Excluding controls next to treated basins
  barely changes the estimate (−0.036). Excluding controls next to *any* enrolled
  land halves it (−0.018), mostly because those distant controls have lower
  baseline loss (0.047 vs 0.065 %/yr), not because of displaced clearing.
- **No detectable leakage.** Never-enrolled neighbours of treated basins show
  +0.025 pp/yr more forest loss than distant controls, the sign leakage would produce,
  but the interval (−0.064 to +0.110) is too wide to conclude anything.
- **Thresholds:** at 10% the effect is −0.015, and at 50% it is zero (§2).
- **NDVI** is slightly negative in the buffered and leakage designs (−0.003 to −0.005,
  about 0.5% of baseline). That is small, inconsistent in sign with avoided
  deforestation, and one of many outcomes tested, so it is not a finding.

## 7. Parcel level: enrolled forest pixels against matched pixels

Basin averaging dilutes any effect confined to agreement land, so the first stage was
also tested where the programme acts (`examples/socio_bosque_parcels.py`;
`results/socio_bosque_parcel*.csv`).

**Design.**

- **Treated:** 19,988 random points inside agreements (60 m inside the boundary),
  area-weighted, plus 3 points in each of 2,208 individual agreements, which area
  weighting otherwise misses. Mangrove agreements are excluded.
- **Controls:** 59,916 random points in mainland Ecuador at least 1 km from any
  agreement. Protected areas and non-forest (tree cover in 2000 < 30%) are dropped
  from both groups.
- **Matching:** 1:1 nearest neighbour (Mahalanobis), exact on region and signing
  cohort (2009–2015), on travel time to cities (2015), distance to forest edge (2000),
  slope, elevation, tree cover, and 2001–04 forest loss within 1 km and 5 km.
- **Outcome:** the annual Hansen loss hazard after signing, among pixels still
  forested at signing.
- **Falsification:** the same comparison in the 4 years *before* signing, re-matched
  among pixels forested 4 years earlier. The matching uses only 2001–04 information,
  so it cannot have absorbed the pre-window.
- **Inference:** bootstrap over agreements (999 draws).

**Balance is near perfect.** Before matching, treated pixels were far more remote,
deeper inside forest and under less local clearing pressure (SMD 0.5–1.2). After
matching, every covariate is within 0.01 SD for the area sample and 0.02 SD for
individual agreements.

**Results** (loss hazard in % of pixels per year; the gap is treated minus control):

| Scope | Agreements | Post gap [90% CI] | Pre gap [90% CI] | Change in gap [90% CI] | Ratio of ratios [90% CI] |
|---|---|---|---|---|---|
| **All regions** (area sample) | 565 | −0.027 [−0.048, −0.012] | −0.010 [−0.032, +0.006] | **−0.017 [−0.043, +0.004]** | 0.53 [0.16, 1.43] |
| Amazon | 142 | −0.021 [−0.044, −0.007] | −0.009 [−0.030, +0.005] | −0.012 [−0.037, +0.009] | 0.79 [0.17, 3.39] |
| Collective agreements | 111 | −0.018 [−0.038, −0.006] | −0.004 [−0.022, +0.010] | −0.014 [−0.038, +0.006] | 0.55 [0.13, 1.69] |
| **Individual agreements** | 1,896 | −0.253 [−0.294, −0.213] | −0.216 [−0.283, −0.159] | **−0.037 [−0.109, +0.037]** | 1.02 [0.52, 2.53] |
| Individual, Andes | 986 | −0.089 [−0.129, −0.048] | −0.038 [−0.099, +0.023] | −0.051 [−0.123, +0.024] | 0.59 [0.17, 1.98] |

Treated hazards in the area sample are tiny: 0.014%/yr after signing against 0.041%/yr
for matched controls. Coast, Andes (area sample) and Amazonian individual agreements
have too few loss events for usable estimates.

**What this shows.**

- **A naive post-enrollment comparison would claim large avoided deforestation:**
  treated pixels lose forest at a third of the matched-control rate overall, and at
  15% of it in individual agreements. Almost all of that gap was already there before
  enrollment. For individual agreements, matched controls on the same terrain, at the
  same distance from cities and forest edge, and under the same local clearing
  pressure lost forest about 7× faster *before* signing. Landholders enroll forest
  they were not going to clear, for reasons no satellite covariate captures.
- **Net of the pre-existing gap, the point estimates are negative but not
  significant** on the difference or the ratio scale, in every scope. The bound is
  informative: on enrolled forest, avoided loss is at most about 0.04 pp/yr (lower
  limit of the main interval). Over 15 years that is at most about 0.6% of enrolled
  forest.
- **Precision is limited by rare events, not sample size.** The main comparison rests
  on 22 losses in treated pixels against 65 in controls after signing, and 9 against 14
  in the 4-year pre-window. Ten times more points would narrow the intervals about
  threefold. That is feasible with the same batch extraction, but the pre-window
  would still contain few events.

**Hydrological implication.** Even at the bound, saving 0.6% of enrolled forest over 15
years, in a basin that is half enrolled, changes basin forest cover by about 0.3%. That
is far below anything the basin panel can detect in ET (±3.5 mm/yr resolution) or water
yield. The parcel analysis closes the loop: the programme's land-cover effect, if any,
is too small to have a measurable hydrological consequence at catchment scale.

## Caveats

- **Survivorship:** only currently active agreements are published. Basins whose
  agreements ended, possibly after clearing, are missing from the treated group. That
  biases toward finding avoided loss, which the best-matched designs still do not find.
- **Treatment timing** uses signing dates. Payments and on-the-ground enforcement may lag.
- **Basin scale:** effects inside agreement boundaries could be diluted at 138 km² basin
  scale. The parcel analysis (§7) addresses this and reaches the same conclusion.
- **Parcel covariates:** travel time is a 2015 map, so it partly reflects roads built
  after many agreements were signed. Edge distances beyond ~10 km are capped (the
  distance transform's search limit).
- **Hansen loss** counts only loss of 2000 forest; regrowth is not tracked.
- **MODIS ET** is model output partly driven by vegetation indices, and coarse (500 m).
- **EVI** was dropped: annual composites give unstable values where its denominator
  approaches zero.
