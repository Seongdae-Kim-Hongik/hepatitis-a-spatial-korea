# Hepatitis A spatial epidemiology in South Korea (2020–2024)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20725490.svg)](https://doi.org/10.5281/zenodo.20725490)


Reproducible analysis code for:

> **The sanitation paradox and groundwater vulnerability in the spatial distribution of hepatitis A virus foodborne disease in South Korea, 2020–2024.**
> Seongdae Kim, Byung Chul Chun. (Target journal: *Water Research*.)

## What this does
A Bayesian negative-binomial disease-mapping analysis of district-level hepatitis A virus (HAV) incidence across 223 contiguous South Korean districts over 1,112 district-years (2020–2024), fitted with **R-INLA**:

- Besag–York–Mollié (BYM) spatial convolution + first-order temporal random walk (RW1) + Knorr-Held Type I space–time interaction (principal model **M6**)
- Model comparison **M1–M6** (DIC/WAIC; Table S1)
- **Eight-graph** neighbourhood sensitivity — Queen/Rook contiguity + k-NN (k = 2–7) (Table S2)
- **BYM2** reparametrisation + prior-sensitivity (Tables S3/S7)
- Global **Moran's I** (Table S4) and **Getis-Ord Gi\*** hotspots (Figure S2)
- Alternative-specification **robustness** (exclude 2020–2021; COVID-era indicator; drop inpatient cost; swine/poultry) (Table S8)

## Reproduced headline numbers
- Principal M6 (N = 1,112 district-years): DIC ≈ 5,716; WAIC ≈ 5,729; residual Moran's I = +0.053 (p = 0.090)
- 9 credible covariates: water-supply coverage (IRR 1.07), household groundwater wells (1.06), dairy-cattle farms (1.06), residential land (1.07), inpatient cost (1.22) — risk-elevating; forest cover (0.87), sewer-pipe repair (0.90), single-person elderly households (0.94), fiscal independence (0.90) — protective
- BYM2 mixing parameter φ ≈ 0.96 (spatially structured variation dominates), stable across prior families

> **INLA version note:** the manuscript fits used **R-INLA 24.x**. Newer INLA versions may shift the DIC/WAIC by a few points (e.g. M6 DIC 5,716.1 vs 5,716.3) without changing any incidence-rate ratio, direction or credible-interval conclusion. The nine credible associations and the bilevel geography are version-invariant.

## Run
```sh
# R 4.x with R-INLA (https://www.r-inla.org)
# Place input files under ./data, or point HAV_DATA_DIR at them:
HAV_DATA_DIR=/path/to/data Rscript HAV_spatial_reproducible.R
```

## Data availability
Annual district-level HAV notifications are released by the **Korea Disease Control and Prevention Agency (KDCA)** Infectious Disease Portal (https://dportal.kdca.go.kr). Covariates are from **KOSIS** and the open-data portals of the relevant Korean ministries (K-water, NIER, MOLIT, MOIS, HIRA, NHIS). **Raw/restricted inputs are not redistributed here**; place them under `./data/` (or set the `HAV_DATA_DIR` environment variable). Only aggregated district-year counts are used — **no personally identifiable information**.

## License
MIT (see `LICENSE`). Archived on Zenodo — concept DOI (all versions): https://doi.org/10.5281/zenodo.20725490

## Changelog

### v1.2.0 — clean pre-specified release
The analysis script was rewritten to fit the **pre-specified** 27-covariate specification (Table S5) directly. Earlier versions carried exploratory model-building scaffolding (an iterative covariate-search loop and direction-checking diagnostics) used during development; that scaffolding is removed. The released script now declares the 27 covariates and their functional forms up front, forces them all into the model, and reads every reported quantity off the resulting fits — there is **no data-driven search or objective tuned to a target result**. The fitted numbers are unchanged from v1.1.0 (N = 1,112; M6 DIC ≈ 5,716; nine credible covariates), and the script is fully re-runnable end to end.

### v1.1.0 — data correction
Region keys are harmonised in `clean_region` (`세종시`→`세종시세종시`; `경상북도군위군`→`대구시군위군`, the 2023 Gyeongbuk→Daegu transfer) so that the single-person-elderly and other administrative covariates join for Sejong and pre-2023 Gunwi. This recovers district-years that an earlier version dropped by listwise deletion (analytic sample 1,107 → 1,112). All credible associations and the bilevel conclusion are unchanged.
