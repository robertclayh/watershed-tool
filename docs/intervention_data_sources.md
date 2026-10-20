# Intervention data sources (checked 2026-10-20)

Parcel-level records of payment-for-ecosystem-services (PES) and watershed programmes,
checked for public access, coverage and fit with `wiat`. "Verified" means the
polygons were actually downloaded, not just listed in a catalogue.

## Summary

| Source | Public polygons? | Treated units | Start dates | Years | Fit |
|---|---|---|---|---|---|
| **Ecuador Socio Bosque** (MAATE) | **Yes, verified** | 2,569 agreements, 1.83 Mha | Signing date per agreement | 2008–2025, mostly 2009–2014 | **Best** |
| **Mexico PSA / PSAH** (CONAFOR) | **Yes, verified, 2018+ only** | ~200–490 projects/yr, 0.17–0.42 Mha/yr | Approval year | 2018–2025 nationally; 2003–2020 Chihuahua only | Good, short post-period |
| Ecuador restoration actions (MAATE) | Yes, verified | 5,831 polygons, ~25 kha | Intervention year | 2019–2025 | Too recent |
| Costa Rica PSA (FONAFIFO) | Not found | – | – | 1997– | Needs a data request |
| Brazil Produtor de Água (ANA) | Not found | 76 projects, 1,283 producers | – | 2006– | Needs a data request |
| Extrema, Conservador das Águas | Not found | Sub-basins in a fixed order | Order only, not years | 2007– | Needs a data request |
| Peru MERESE (SUNASS/MINAM) | Not found | – | – | 2014– | Registry only |

## Ecuador: Socio Bosque conservation agreements

- WFS: `http://ide.ambiente.gob.ec:8080/geoserver/ows`, layer
  `mae_ide:v_hc005_area_bajo_conservacion_a_anonimo` (the `https://ide.` and
  `geoportal.` hosts list the layer but return empty responses).
- Fields: `codcnv` (agreement code, e.g. `MAE-PSB-II-2011-I-204`), `fec` (signing
  date), `tpcnv` (INDIVIDUAL 2,317 / COLECTIVO 206 / MANGLAR 46), `abc` (hectares),
  `estejc` (status; all ACTIVO), geometry (MultiPolygon, projected metres, likely
  UTM 17S).
- Enrollment by signing year (agreements, kha): 2008: 49, 72 · 2009: 309, 293 ·
  2010: 493, 238 · 2011: 499, 433 · 2012: 406, 196 · 2013: 360, 47 · 2014: 389, 189 ·
  2015–2025: 62 agreements, ~360 kha.
- Parcel size: median 40 ha, 90th percentile 312 ha. Collective agreements reach
  10,000–20,000 ha.
- Caveats: only *active* agreements are published, so terminated agreements are
  missing (survivorship). Socio Bosque pays to keep standing forest, so any
  hydrological effect is avoided loss rather than restoration. Enrollment is
  voluntary and targeted, which calls for matching on deforestation pressure.

## Mexico: CONAFOR PSA

- WFS: `https://idefor.cnf.gob.mx/geoserver/ows`. Layers `geonode:psa_2018_geo_1`,
  `psa_2019_geo_1`, `psa_2020_geo`, `psa_2021_geo`, `psa_2022_geo`, `psa_2023_geo`,
  `psa_2024_geo`, `psa_2025_geo`. Large layers download slowly (minutes); the 2021
  file arrived truncated once and 2025 timed out, so retry or page with
  `startIndex`.
- Fields: `FOL_PROG`, `ESTADO`, `MUNICIPIO`, `SOLICITANT`, `REG_PROP`, `MODALIDAD`,
  `SUPERFICIE` (ha), `EJERCICIO` (year), `ESQUEMA`.
- Verified counts: 2018: 420 projects, 316 kha (173 hydrological) · 2019: 489,
  417 kha (252 hydrological) · 2020: 242, 187 kha · 2022: 206, 174 kha · 2023: 217,
  237 kha · 2024: 256, 265 kha. From 2020 the hydrological and biodiversity
  modalities are merged into one PSA modality.
- Also public: `mlpsa_fc_2018_geo` … `mlpsa_fc_2025_geo` (local matched-fund PES
  mechanisms, the closest Mexican analogue to water funds),
  `comp_iv_nacional_rop2024_4326` (national PSA *eligible areas*, 2024: a natural
  pool of eligible-but-not-enrolled controls), and `psa_tle_2003_2020_geo`
  (Chihuahua only, 432 projects, 2003–2020).
- Not public: national polygons for 2003–2017, the period used in the main
  evaluations of PSAH (e.g. Alix-Garcia, Sims and Yañez-Pagans 2015). These would
  need a transparency request to CONAFOR.
