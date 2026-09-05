# Hepatitis A spatial epidemiology in South Korea (2020–2024)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20725490.svg)](https://doi.org/10.5281/zenodo.20725490)

Reproducibility package for:

> **Piped water, private wells, and forest cover shape hepatitis A geography in South Korea.**
> Seongdae Kim, Byung Chul Chun. (Submitted to *Frontiers in Ecology and the Environment* as a Research Communication.)

Earlier releases of this repository accompanied the same analysis under previous working titles; the concept DOI above always resolves to the latest version.

## What this does
A Bayesian negative-binomial disease-mapping analysis of district-level hepatitis A virus (HAV) incidence across 223 contiguous South Korean districts over 1,112 district-years (2020–2024), fitted with **R-INLA**:

- Besag–York–Mollié (BYM) spatial convolution + first-order temporal random walk (RW1) + Knorr-Held Type I space–time interaction (principal model **M6**), log(population + 1) offset, PC priors (PC.prec(0.5, 0.01))
- 27 district-level covariates (water infrastructure, groundwater use and testing, land cover, livestock, shellfish, demographic, fiscal and health-care measures); biennially available administrative measures are carried to adjacent years
- Model comparison **M1–M6** (DIC/WAIC), crude and residual global **Moran's I**, CPO/PIT diagnostics, VIF audit
- Optional (`RUN_EXTENDED=true`): eight-graph neighbourhood sensitivity (Queen/Rook + k-NN, k = 2–7), Getis-Ord Gi\* local clustering, alternative-specification robustness checks

**Interpretation.** The final 27-covariate specification was developed through exploratory model building. This package reproduces the selected final model; it does not convert the analysis into a prospectively specified confirmatory study. The reported credible associations are hypothesis-generating.

## Reproduced headline numbers (principal model M6, N = 1,112)
- DIC ≈ 5,716; WAIC ≈ 5,729; residual Moran's I = +0.053 (p = 0.090)
- 9 credible covariates: piped-water coverage, household groundwater wells, dairy-cattle farms, residential land area, inpatient cost — risk-elevating; forest cover, sewer-pipe repair, single-person elderly households, municipal fiscal independence — protective
- Seven of the nine associations persist under all eight neighbourhood definitions; residential land area is the only purely local (adjacent-district) signal

## Contents
| Path | What it is |
|---|---|
| `HAV_spatial_reproducible.R` | Corrected, directly runnable implementation of the final model (v2.0.0) |
| `results/analysis_dataset_compiled.csv` | The exact 1,112-row district-year analytic table used by the verified run |
| `results/table2_principal_IRR.csv` | 27 incidence-rate ratios with 95% credible intervals (Table 1 / Appendix S1: Table S2) |
| `results/tableS1_model_comparison.csv` | M1–M6 DIC/WAIC (Appendix S1: Table S3a) |
| `results/core_diagnostics.csv` | Fit criteria and crude/residual Moran's I |
| `results/cpo_pit_diagnostics.csv` | CPO and PIT values, all district-years |
| `results/vif_audit.csv` | Variance-inflation audit of the 27 covariates |
| `results/fast_principal_run.log`, `results/core_run.log` | Logs of the verified runs |
| `results/sessionInfo.txt` | R 4.6.0 / INLA 25.10.19 session information |
| `DATA_PROVENANCE.md` | Source domains and redistribution note for the compiled file |
| `CITATION.cff`, `.zenodo.json`, `LICENSE` | Citation and archive metadata; MIT license |

## Run
```sh
# R 4.x with R-INLA (https://www.r-inla.org)
# Raw inputs (KDCA surveillance extract, KOSIS and ministry files, district shapefile)
# go under ./data, or point HAV_DATA_DIR at them. Outputs go to ./results (HAV_OUTPUT_DIR).
HAV_DATA_DIR=/path/to/data FAST_PRINCIPAL=true Rscript HAV_spatial_reproducible.R
# Full core run (principal model + diagnostics):
HAV_DATA_DIR=/path/to/data Rscript HAV_spatial_reproducible.R
# Add neighbourhood-graph, Gi* and alternative-specification sections:
HAV_DATA_DIR=/path/to/data RUN_EXTENDED=true Rscript HAV_spatial_reproducible.R
```
`FAST_PRINCIPAL=true` is the verified path and exits after the principal model and core diagnostics. Minor DIC/WAIC differences across INLA versions are expected; compare effect directions and credible-interval conclusions as well as fit criteria. The groundwater-quality variable is a testing count, not a pass rate.

## Data availability
Annual district-level HAV notifications are released by the **Korea Disease Control and Prevention Agency (KDCA)** Infectious Disease Portal (https://dportal.kdca.go.kr). Covariates are from the **Korean Statistical Information Service (KOSIS)** and the open-data portals of the relevant Korean ministries and agencies. Raw source extracts are not redistributed here; the compiled district-year analytic table is provided in `results/` (see `DATA_PROVENANCE.md`). Only aggregated district-year counts are used — **no personally identifiable information**.

## License
MIT (see `LICENSE`). Archived on Zenodo — concept DOI (all versions): https://doi.org/10.5281/zenodo.20725490

## Changelog

### v2.0.0 — corrected reproducibility release (2026-09-05)
- The analysis script is replaced by a corrected implementation. Fixes relative to v1.x: the `nb2INLA` argument order in the graph export; the precision prior is now PC.prec(0.5, 0.01) as reported; covariates enter as district-year values with biennial carry-forward rather than five-year means; the groundwater-quality variable is documented as a testing count; quantile-coded covariates are documented as ordinal scores (1-SD effects).
- The description of the covariate set is corrected: it was developed exploratorily, not pre-specified. The associated manuscript is framed as hypothesis-generating accordingly.
- BYM2 prior-sensitivity results were withdrawn after audit (the archived fits used 1,107 rather than 1,112 district-years, and refits were numerically unstable); CPO-based outlier counts are no longer reported.
- The compiled 1,112-row analytic dataset, principal-model outputs, diagnostics, logs and session information are now committed under `results/`.
- The reproduced point estimates are unchanged: nine credible associations, M6 DIC ≈ 5,716, residual Moran's I = +0.053.

### v1.2.0 — clean release
Analysis script rewritten without development scaffolding (an iterative covariate-search loop and direction-checking diagnostics). Superseded by v2.0.0, which also corrects the pre-specification claim made in this release.

### v1.1.0 — data correction
Region keys harmonised (`세종시`, `경상북도군위군`) so administrative covariates join for Sejong and pre-2023 Gunwi; analytic sample 1,107 → 1,112.
