# Hepatitis A spatial epidemiology in South Korea (2020–2024)

Reproducible analysis code for:

> **The sanitation paradox and groundwater vulnerability in the spatial distribution of hepatitis A virus foodborne disease in South Korea, 2020–2024.**
> Seongdae Kim, Byung Chul Chun. (Target journal: *Water Research*.)

## What this does
A Bayesian negative-binomial disease-mapping analysis of district-level hepatitis A virus (HAV) incidence across 223 contiguous South Korean districts over 1,107 district-years (2020–2024), fitted with **R-INLA**:

- Besag–York–Mollié (BYM) spatial convolution + first-order temporal random walk (RW1) + Knorr-Held Type I space–time interaction (principal model **M6**)
- Model comparison **M1–M6** (DIC/WAIC; Table S1)
- **Eight-graph** neighbourhood sensitivity — Queen/Rook contiguity + k-NN (k = 2–7) (Table S2)
- **BYM2** reparametrisation + prior-sensitivity (Tables S3/S7)
- Global **Moran's I** (Table S4) and **Getis-Ord Gi\*** hotspots (Figure S2)
- Alternative-specification **robustness** (exclude 2020–2021; COVID-era indicator; drop inpatient cost; swine/poultry) (Table S8)

## Reproduced headline numbers
- Principal M6: DIC ≈ 5,699; residual Moran's I = +0.05 (p ≈ 0.09)
- 9 credible covariates incl. water-supply coverage (IRR 1.07), household groundwater wells (1.06), forest (0.87), inpatient cost (1.21)

> **INLA version note:** the manuscript used **R-INLA 23.12.16**. Pin this version to match the DIC exactly; newer INLA may shift DIC by ~3 points (e.g. 5,696 vs 5,699) without changing any incidence rate ratio, direction or credibility.

## Run
```r
# R 4.x with R-INLA (https://www.r-inla.org)
Rscript HAV_spatial_reproducible.R
```

## Data availability
Annual district-level HAV notifications are released by the **Korea Disease Control and Prevention Agency (KDCA)** Infectious Disease Portal (https://dportal.kdca.go.kr). Covariates are from **KOSIS** and the open-data portals of the relevant Korean ministries (K-water, NIER, MOLIT, MOIS, HIRA, NHIS). **Raw/restricted inputs are not redistributed here**; place them under `./data/` and set `BASE_IV`. Only aggregated district-year counts are used — **no personally identifiable information**.

## License
MIT (see `LICENSE`). Archived on Zenodo (concept DOI to be added on release).
