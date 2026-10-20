# Methodology

This note describes what `wiat` computes and why. For results on when the methods
fail, see [where_attribution_fails.md](where_attribution_fails.md).

## 1. The attribution problem

A nature-based intervention (riparian restoration, páramo protection, reforestation,
payments for ecosystem services) is implemented in basin *i* from period *g*. We want
the average effect on a hydrologic outcome *Y* after *g*:

$$\tau = \frac{1}{|\text{post}|}\sum_{t \ge g} \big(Y_{it}(1) - Y_{it}(0)\big)$$

$Y_{it}(0)$, the basin's trajectory had nothing been done, is never observed. Every
estimator here builds it from **other basins** (donors) that were not treated. The
quality of the attribution is bounded by the quality and number of those basins.

## 2. Data layer (Google Earth Engine)

`wiat.data.gee` extracts a monthly basin panel:

| Variable | Source | Processing |
|---|---|---|
| `ndvi`, `evi`, `ndmi` | Landsat 5/7/8/9 C2 L2, Sentinel-2 L2A | QA_PIXEL / SCL cloud, shadow and snow masks. ETM+/TM→OLI (Roy et al. 2016) and MSI→OLI (HLS v1.4) harmonisation. Per-period median composite. `veg_clear_frac` records cloud-free basin share. |
| `precip_mm` | CHIRPS v2 daily | Sum over the period. |
| `et_mm`, `pet_mm` | MODIS MOD16A2GF v061 | 8-day sums → daily rates (short year-end composite handled) × days of overlap with the period. Fill values masked. |
| `elevation_m`, `slope_deg` | NASADEM | Basin mean. |
| `forest_frac`, `forest_loss_frac` | Hansen GFC v1.13 | Tree cover 2000, loss 2001–2025. |
| basins, topology | HydroSHEDS HydroBASINS | Level 1–12 polygons; `NEXT_DOWN` used to exclude connected basins. |

Derived outcomes (`wiat.panel.add_hydrologic_outcomes`):

- `water_yield_mm = P − ET`: over annual or longer periods, the water a basin
  releases downstream (runoff + recharge).
- `evaporative_index = ET / P`: the Budyko evaporative index. Vegetation
  interventions shift it, and normalising by P removes much interannual climate
  noise.

Monthly panels are aggregated to calendar or water years (`aggregate`). Fluxes are
summed, indices averaged, and incomplete years set to missing rather than biased.
Intervention dates are mapped to the first *fully* treated year.

### Caveats specific to remote sensing

- MODIS ET is itself a model driven partly by vegetation indices. An NDVI increase
  mechanically raises MOD16 ET, so ET effects are not independent of NDVI effects.
  Treat them as one line of evidence, not two.
- CHIRPS is coarse (0.05°) and blends station data that are sparse in the Andes, so
  precipitation errors in small headwater basins are large. Using precipitation as
  a DiD covariate attenuates, but does not remove, climate confounding.
- Persistent cloud (Andes, Amazon, Atlantic Forest) leaves gaps in the optical
  record. `quality_mask` blanks low-coverage months, and `balance` then drops
  basins with gaps. Both silently thin the donor pool.

## 3. Donor selection (`wiat.matching`)

Each basin is summarised by pre-intervention **signatures**: mean P, P seasonality
(CV of the monthly climatology), aridity (ET/P), baseline outcome mean and trend.
These are joined with static physiography. Covariates are standardised across the
region. A candidate is *suitable* if:

- its root-mean-square standardised difference from the treated basin is ≤ `caliper`
  (default 1.0 SD), and
- no single covariate differs by more than `hard_caliper` (default 2.5 SD).

Candidates are excluded if they are hydrologically connected to the treated basin
(SUTVA: shared water and often a shared programme), carry their own interventions,
or have missing data. The number of suitable donors, `J`, is the central quantity
in the rest of the toolkit.

## 4. Difference-in-differences (`wiat.did`)

Two-way fixed effects with time-varying covariates (default: precipitation):

$$Y_{it} = \alpha_i + \lambda_t + \tau D_{it} + X_{it}'\beta + \varepsilon_{it}$$

estimated by the within transformation on a balanced panel. Inference options, all
reported side by side:

| Method | Valid when | Behaviour with one treated basin |
|---|---|---|
| Cluster-robust (CR1, *t*<sub>G−1</sub>) | many treated *and* control clusters | Severely over-rejects (see results) |
| Wild cluster restricted bootstrap (Rademacher enumeration or Webb weights), inverted for a CI | moderate number of treated clusters | Degenerate. Intervals are extremely wide or uninformative (MacKinnon & Webb 2017) |
| **Placebo / permutation** (Conley & Taber 2011) | few treated, exchangeable controls | Valid. Resolution is 1/(J+1) |

The default (`inference="auto"`) is placebo inference when fewer than five basins are
treated. Treatment is assigned in turn to every control (or to random sets of
controls with the same start dates), giving a distribution of placebo effects. The
test rejects when |τ̂ − θ| is in the top α share of |placebo| values. Inverting it
gives `τ̂ ± c`, where `c` is the **detection threshold**: the smallest effect the design
can distinguish from untreated-basin variability. The interval excludes zero exactly
when *p* ≤ α.

**Event study.** Leads and lags with binned endpoints and reference period −1. It
reports cluster-robust coefficients, a joint pre-trend F-test, and, with few treated
basins, placebo bands plus a placebo-based pre-trend test. That test standardises
each lead by its placebo spread and ranks the sum of squares.

**Staggered adoption.** TWFE is biased when start dates vary and effects evolve.
`staggered_att` implements Callaway & Sant'Anna (2021) group-time ATTs with
never-treated comparisons, a universal base period, cohort-size weighting and a
cohort-stratified unit bootstrap.

## 5. Synthetic control (`wiat.synth`)

Weights *w* on the simplex minimise the pre-period discrepancy between the treated
basin and the weighted donors (optionally plus standardised static covariates). The
problem is solved exactly as non-negative least squares with a heavily weighted
sum-to-one row (Lawson–Hanson). FISTA is the fallback, because SCM problems are
usually rank-deficient (J > T₀) and first-order methods converge slowly there.

| `method` | Counterfactual | Use when |
|---|---|---|
| `standard` | Σ w<sub>j</sub> Y<sub>jt</sub> | the treated basin sits inside the donors' range in levels |
| `demeaned` (default) | Σ w<sub>j</sub> Y<sub>jt</sub> + (Ȳ<sub>1,pre</sub> − Σ w<sub>j</sub> Ȳ<sub>j,pre</sub>) | the treated basin is wetter or drier than all donors (Ferman & Pinto 2021) |
| `augmented` | de-meaned SCM + ridge correction for remaining pre-period imbalance | pre-fit is imperfect. It extrapolates, so read it with the hull diagnostics (Ben-Michael et al. 2021) |

**Inference.** Placebo-in-space: each donor is treated as if it were the treated
basin, with the other donors as its pool. The test statistic is |ATT| / pre-RMSPE.
Abadie's post/pre RMSPE ratio is also reported. Placebos whose pre-fit is more than
5× worse than the treated basin's are discarded, following Abadie et al. (2010).
Inverting the test (Firpo & Possebom 2018) gives `ATT ± c · pre-RMSPE`, unbounded
when α(J+1) < 1. A conformal interval (Chernozhukov, Wüthrich & Zhu 2021) is
available (`conformal=True`). It is only informative when the pre-period is much
longer than the post-period, which is rarely true for annual hydrologic records.

**Fit diagnostics.**

- `pre_rmspe_normalized`: pre-RMSPE / SD of the treated basin's pre-series. Above
  ~0.5, the synthetic basin does not track the real one.
- `frac_pre_outside_donor_range`: share of pre-periods where the treated basin lies
  outside every donor, so no convex combination can match it.
- `effective_donors = 1/Σw²`.

## 6. Failure diagnostics and verdicts (`wiat.attribution`)

`AttributionStudy.run()` evaluates each condition below and reports all of them. It
never silently drops a study.

| Code | Severity | Trigger |
|---|---|---|
| `DONOR_SCARCITY` | fatal | J < 1/α − 1 (9 at α = 0.10): significance is impossible by construction |
| `UNMATCHED_DONORS_USED` | major | fewer than `min_donors` pass the caliper; nearest unmatched basins were added |
| `PRE_TREND` | major | event-study pre-trend test rejects |
| `POOR_PRE_FIT` | major | SCM normalised pre-RMSPE > 0.5 |
| `OUTSIDE_CONVEX_HULL` | major | > 25% of pre-periods outside the donor range |
| `PLACEBOS_DISCARDED` | major | after discarding poor-fit placebos, SCM can no longer reach α |
| `ESTIMATOR_DISAGREEMENT` | major | DiD and SCM differ in sign, or SCM lies outside the DiD interval |
| `THIN_PLACEBO_DISTRIBUTION` | minor | J < 2/α − 1: the threshold rests on one or two placebos |
| `TREATED_OUTSIDE_DONOR_SUPPORT` | minor | treated basin outside the donor range on a matching covariate |
| `CONCENTRATED_WEIGHTS` | minor | SCM effective donors < 2 |
| `STAGGERED_ADOPTION` | minor | multiple start dates (CS estimator also reported) |
| `SHORT_BASELINE` | minor | fewer than 5 pre-periods |
| `SINGLE_TREATED_BASIN`, `BELOW_DETECTION_THRESHOLD` | info | context for reading the estimates |

Verdict: any fatal → **NOT ATTRIBUTABLE (DESIGN)**. Otherwise any major →
**INCONCLUSIVE**. Otherwise DiD and SCM significant with the same sign →
**ATTRIBUTED**. Otherwise neither significant → **NO DETECTABLE EFFECT**, read
together with the detection threshold.

## References

- Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies. *JASA*, 105(490).
- Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The augmented synthetic control method. *JASA*, 116(536).
- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *J. Econometrics*, 225(2).
- Chernozhukov, V., Wüthrich, K., & Zhu, Y. (2021). An exact and robust conformal inference method for counterfactual and synthetic controls. *JASA*, 116(536).
- Claverie, M., et al. (2018). The Harmonized Landsat and Sentinel-2 surface reflectance data set. *Remote Sens. Environ.*, 219.
- Conley, T. G., & Taber, C. R. (2011). Inference with "difference in differences" with a small number of policy changes. *REStat*, 93(1).
- Ferman, B., & Pinto, C. (2021). Synthetic controls with imperfect pretreatment fit. *Quantitative Economics*, 12(4).
- Firpo, S., & Possebom, V. (2018). Synthetic control method: Inference, sensitivity analysis and confidence sets. *J. Causal Inference*, 6(2).
- Funk, C., et al. (2015). The climate hazards infrared precipitation with stations (CHIRPS). *Scientific Data*, 2.
- Lehner, B., & Grill, G. (2013). Global river hydrography and network routing (HydroSHEDS/HydroBASINS). *Hydrol. Process.*, 27(15).
- MacKinnon, J. G., & Webb, M. D. (2017). Wild bootstrap inference for wildly different cluster sizes. *J. Applied Econometrics*, 32(2).
- Roy, D. P., et al. (2016). Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity. *Remote Sens. Environ.*, 185.
- Running, S. W., Mu, Q., Zhao, M., & Moreno, A. (2019). MODIS Global Terrestrial Evapotranspiration (MOD16A2GF) User's Guide.
