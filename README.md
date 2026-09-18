# Hepatitis A spatial epidemiology in South Korea (2020–2024)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20725490.svg)](https://doi.org/10.5281/zenodo.20725490)

Reproducibility package for:

> **Spatial Clustering of Hepatitis A in South Korea, 2020-2024: Nationwide Bayesian Analysis of Groundwater, Land Cover, and Socioeconomic Gradients**
> Seongdae Kim, Byung Chul Chun. (Prepared for submission to *JMIR Public Health and Surveillance* as an Original Paper.)

The concept DOI above always resolves to the latest version. Releases up to v2.0.6 accompanied an earlier version of this analysis under the title *Sanitation Infrastructure and Environmental Vulnerability of Hepatitis A Transmission in South Korea, 2020-2024*.

## Important: v2.1 corrects a data-quality defect, and v2.2.0 adopts a further data-quality rule found afterwards
An audit of the compiled dataset shipped with v2.0.x found that several administrative source files store **unavailable values as 0 instead of missing** and that some series were **redefined in particular years**. Under binary and quantile coding, these values made covariates encode *which years or regions had data* rather than their true level:

- piped water-supply coverage was 0 in every study year for the 76 autonomous districts of the metropolitan cities and Jeju (77 from 2023), because it has been published only for each metropolitan city as a whole since 2018;
- whole survey years were 0 for the food-security rate (2020–2023), basic livelihood security recipients and wastewater-discharging facilities (2021–2023; 2024 absent), older adults living alone (2021) and groundwater-quality tests (every year before 2023; 2024 absent);
- inpatient medical cost was 0 in 85 districts in 2023 (not released; 4 in each of the 3 preceding years) and absent for 2024 in all but 4 districts, where it was stored as 0 (`results/inpatient_cost_zero_audit.csv`), and its 2019 level was about three times that of any other year;
- the series for older adults living alone was replaced by a different indicator in 2024 (national median 23.6 → 59.5);
- the urban-area population proportion exceeded 100% in 2020 in 10 cities with non-autonomous wards (ward values summed).

v2.1 recodes these values as missing, fills them with the most recent earlier available year of the same district (the nearest later year when no earlier value exists), and logs every change in `results/data_repair_log.csv` (rules R1–R6 in the script, section `[2b]`, and in `DATA_PROVENANCE.md`). **The v2.0.x headline associations for piped water-supply coverage and inpatient medical cost, and the v2.0.x claim of nine credible covariates, were artefacts of this defect and are withdrawn.** A second artefact (an apparent protective association of older adults living alone, created by the 2024 redefinition) appeared after the first repair and was removed by the series-break rule; it is documented in the manuscript. Two simple checks are needed to expose such problems: the share of zeros (region-wide and partial zero coding) and the share of each covariate's variance explained by calendar year (`qa/year_r2_screen.py`; whole-year zero coding and redefined series; `results/year_r2_before_after.csv`). Household groundwater wells and small-scale water-supply facilities are available only through 2022 and dairy-cattle farms through 2023; later years carry the last available value. Because only one year of data exists, three covariates (basic livelihood security recipients, wastewater-discharging facilities, groundwater-quality tests) are constant over the study years.

v2.1 also fixes a descriptive-statistics bug: the Getis-Ord Gi\* input treated the three districts without covariate data as having zero incidence, which created spurious cold spots. Gi\* now uses the observed incidence of all 223 graph districts.

**v2.2.0.** Inspecting the v2.1 repaired analysis with a post hoc extreme-value screen (values of rate, percentage and cost covariates lying more than 5 interquartile ranges beyond the quartiles) found 17 further apparent aggregation errors in 4 covariates, mostly in cities with non-autonomous wards (eg, a sex ratio near 50 males per 100 females; basic livelihood security recipients above 30%; `results/extreme_value_flags.csv`). Two of these, Suwon-si and Seongnam-si, had every study-year value of basic livelihood security recipients flagged and could not be filled, so they leave the analysis under the complete-case rule. v2.2.0 adds this screen as rule R7 of the principal analysis (`DATA_PROVENANCE.md`), reducing the analytic frame to **218 districts, 1,090 district-years**, and adds sewer-pipe repair sites as a fifth credible covariate. The v2.1.x claim of four credible covariates and the v2.1.x analytic frame (220 districts, 1,100 district-years) are superseded; the R7-free analysis is retained as a sensitivity check (`results/sens_extreme_values.csv`, specification `flagged_values_retained`).

## What this does
A Bayesian negative-binomial disease-mapping analysis of district-level hepatitis A virus (HAV) incidence, fitted with **R-INLA**:

- contiguity graph of 223 districts (229 national districts minus 6 islands without a land neighbour); **218 districts × 5 years = 1,090 district-years** have complete, non-flagged covariates and enter the likelihood (Jeju-si and Seogwipo-si have no fiscal statistics of their own; Jeju-si and Busanjin-gu have no non-zero sewerage-coverage value; Suwon-si and Seongnam-si have every study-year value of basic livelihood security recipients flagged by rule R7)
- Besag–York–Mollié (BYM) spatial convolution + first-order temporal random walk (RW1) + Knorr-Held Type I space–time interaction (principal model **M6**), log(population + 1) offset, PC priors (PC.prec(0.5, 0.01))
- 27 district-level covariates (water infrastructure, groundwater use and testing, land cover, livestock, shellfish, demographic, fiscal and health-care measures); applied transformations, cut points and descriptive statistics are written to `results/covariate_dictionary.csv`. Net migration enters as log(1 + max(x, 0)), ie, net out-migration is set to 0, as in the specification developed during model building
- model comparison **M1–M6** (DIC, WAIC, effective parameters), crude and residual global **Moran's I** (count and Pearson residuals), CPO/PIT diagnostics, VIF audit
- `RUN_EXTENDED=true`: eight-graph neighbourhood sensitivity (Queen/Rook + symmetrised k-NN, k = 2–7), Getis-Ord Gi\*, alternative specifications (excluding 2020–2021, COVID-19 indicator, without inpatient cost, with swine and poultry farms, two treatments of oyster production, three sensitivity analyses of the handling of extreme values — R7-flagged values retained, sewer-pipe repair sites log transformed, or both), covariate estimates under nested random-effect structures, and M6 without covariates
- `RUN_STABILITY=true`: eight INLA approximation/prior/parameterisation settings (all 27 covariates monitored), a glmmTMB cross-fit without spatial structure, and repeated fits of every nested model

**Interpretation.** The 27-covariate specification was developed through exploratory model building and is not prospectively specified. The reported credible associations are hypothesis-generating.

## Headline numbers (principal model M6, N = 1,090; `results/`)
- M6: DIC 5,618.36; WAIC 5,629.59; effective parameters 148.2. INLA fits of this model are not bit-reproducible: over the archived fit and 20 repeated fits with identical settings the DIC ranged 5,618.02–5,619.97 and the WAIC 5,629.43–5,631.48, the incidence-rate ratios differed in the fourth decimal place, and the credible set (the same 5 covariates) was identical in every fit (`results/repeated_fits.csv`)
- **5 credible covariates** (IRR = exponentiated posterior mean of the coefficient, per 1 SD of the transformed covariate; 95% CrI): household groundwater wells 1.077 (1.001–1.158); dairy-cattle farms 1.042 (1.003–1.083); sewer-pipe repair sites 0.907 (0.826–0.997); forest area 0.836 (0.757–0.925); fiscal independence 0.861 (0.789–0.940). Forest area and fiscal independence are credible under all 8 neighbourhood graphs and all 8 INLA settings; dairy-cattle farms under 7 of 8 graphs and all 8 INLA settings; sewer-pipe repair sites under all 8 graphs but only 7 of 8 INLA settings; household groundwater wells under 5 of 8 graphs and 6 of 8 INLA settings — the last 2 are therefore borderline
- not credible: piped water-supply coverage 1.018 (0.937–1.106); inpatient medical cost 1.062 (0.974–1.158); older adults living alone 0.962 (0.908–1.018); small-scale water-supply facilities 1.045 (0.983–1.112); food-security rate 1.031 (0.990–1.074); basic livelihood security recipients 1.005 (0.950–1.063); sex ratio 1.061 (0.991–1.136); no covariate outside the 5 was credible under any neighbourhood graph or INLA setting
- excluding 2020–2021, only household groundwater wells and forest area remain credible; dairy-cattle farms, sewer-pipe repair sites and fiscal independence keep their direction with intervals including 1, and the food-security rate and the urban-area population proportion, neither credible in the principal model, have intervals excluding 1 (`results/alternative_specifications.csv`)
- oyster production: the 5 associations are maintained when the annual score is replaced by a time-invariant producer indicator (under which the indicator itself, 0.948 (0.901–0.997), is also credible) or removed (`results/sens_oyster_coverage.csv`)
- sensitivity to the handling of extreme values (`results/sens_extreme_values.csv`): dairy-cattle farms, forest area and fiscal independence remain credible in all 3 analyses. Without rule R7 (the 17 flagged values retained, N = 1,100), sewer-pipe repair sites is not credible (0.911 (0.829–1.000)) but household groundwater wells is (1.079 (1.006–1.157)). With sewer-pipe repair sites log transformed (N = 1,090 with R7, N = 1,100 without), it is not credible either way (1.014 (0.977–1.051); 1.019 (0.983–1.056)), and household groundwater wells is not credible either way (1.067 (0.991–1.148); 1.067 (0.995–1.144))
- crude Moran's I = +0.735 (223 districts); residual Moran's I = +0.212 (count residuals, p < 0.001) and +0.101 (Pearson residuals, p = 0.008); fitted values include the random effects
- the covariates take up only part of the spatial structure: the posterior median SD of the structured spatial component is 0.29 without and 0.25 with the 27 covariates, the covariate-free model has a DIC (5,620.61) and WAIC (5,632.91) within 3 points of M6's, which is of the order of the variation between repeated fits, and the number of districts with credibly elevated/reduced effects is 78/58 without and 74/56 with covariates (218 analysed districts; `results/spatial_structure_with_without_covariates.csv`)
- Gi\* (223 districts, 95% level): 51 hot spots, 56 cold spots
- M4 (interaction without a temporal main effect) is numerically unstable — effective parameters 650–718 of 1,090, DIC 5,645.20–5,859.04 across fits — and is not used for model selection

## Contents
| Path | What it is |
|---|---|
| `HAV_spatial_reproducible.R` | Single script: data repair, model fits, all sensitivity and stability analyses |
| `make_figures.R` | Figure 2 and the supplementary model figures from the script outputs |
| `qa/year_r2_screen.py`, `qa/year_r2_before_after.py` | Year-R² screen for covariates that encode data availability instead of a district characteristic |
| `qa/inpatient_cost_zero_audit.py`, `results/inpatient_cost_zero_audit.csv` | Year-by-year counts of 0 and missing values of inpatient medical cost in the source file; reconciles the 96 recoded values and the 321 missing values of `results/data_repair_log.csv` |
| `results/analysis_dataset_compiled.csv` | The exact 1,090-row district-year analytic table (repaired, rule R7 applied) |
| `results/data_repair_log.csv`, `results/year_r2_before_after.csv` | Per-variable counts of values recoded as missing and of missing values before/after filling; share of covariate variance explained by calendar year before (v2.0.6 dataset), after recovery of the zero-coded values alone, and after the full repair |
| `results/covariate_dictionary.csv` | Applied transformation, cut points, descriptive statistics, number of districts constant over years |
| `results/table2_principal_IRR.csv` | 27 incidence-rate ratios with 95% credible intervals (manuscript Table 2; unrounded) |
| `results/model_comparison.csv`, `results/repeated_fits.csv` | M1–M6 DIC/WAIC/effective parameters; run-to-run ranges |
| `results/core_diagnostics.csv`, `results/descriptive_223.csv`, `results/m6_hyperparameter_sd.csv`, `results/spatial_structure_with_without_covariates.csv` | Fit criteria, Moran's I, risk classification, hyperparameters on the SD scale, M6 with and without covariates |
| `results/m6_outputs.rds` | Fitted values, district effects and CPO/PIT of the archived M6 (input of `make_figures.R`) |
| `results/district_effects.csv` | Posterior combined district effects (the 5 districts not analysed — Jeju-si, Seogwipo-si, Busanjin-gu, Suwon-si, Seongnam-si — carry interpolated effects and are not classified) |
| `results/graph_sensitivity.csv`, `results/graph_counts_analysed220.csv` | Eight-graph sensitivity |
| `results/alternative_specifications.csv`, `results/sens_oyster_coverage.csv`, `results/sens_extreme_values.csv`, `results/extreme_value_flags.csv`, `results/nested_model_covariates.csv` | Alternative specifications, sensitivity to the handling of extreme values (rule R7) and random-effect structures (all 27 covariates) |
| `results/stability_inla_variants.csv`, `results/stability_hyperparameters.csv`, `results/stability_glmmTMB.csv` | Numerical stability (all 27 covariates) |
| `results/getis_ord_gi.csv`, `results/cpo_pit_diagnostics.csv`, `results/vif_audit.csv` | Gi\* z scores, CPO/PIT values, VIF audit |
| `results/full_run.log`, `results/sessionInfo.txt` | Log of the archived run; R 4.6.0 / INLA 25.10.19 session information |
| `DATA_PROVENANCE.md` | Source domains, repair rules and redistribution note |
| `CITATION.cff`, `.zenodo.json`, `LICENSE` | Citation and archive metadata; MIT license |

Output file names carry no supplementary-table numbers, because the manuscript numbers its supplementary tables in order of first citation.

## Run
```sh
# R 4.x with R-INLA (https://www.r-inla.org); glmmTMB for the cross-fit
# Raw inputs (KDCA surveillance extract, KOSIS and ministry files, district shapefile)
# go under ./data, or point HAV_DATA_DIR at them. Outputs go to ./results (HAV_OUTPUT_DIR).
HAV_DATA_DIR=/path/to/data FAST_PRINCIPAL=true Rscript HAV_spatial_reproducible.R      # principal model only
HAV_DATA_DIR=/path/to/data Rscript HAV_spatial_reproducible.R                          # + model comparison and Moran's I
HAV_DATA_DIR=/path/to/data RUN_EXTENDED=true RUN_STABILITY=true Rscript HAV_spatial_reproducible.R   # everything (archived run)
HAV_DATA_DIR=/path/to/data Rscript make_figures.R
python3 qa/year_r2_screen.py results/analysis_dataset_compiled.csv
```
Compare effect directions and credible-interval conclusions rather than the second decimal of DIC/WAIC (see the run-to-run ranges above). The groundwater-quality variable is a testing count, not a pass rate.

## Data availability
Annual district-level HAV notifications are released by the **Korea Disease Control and Prevention Agency (KDCA)** Infectious Disease Portal (https://dportal.kdca.go.kr). Covariates are from the **Korean Statistical Information Service (KOSIS)** and the open-data portals of the relevant Korean ministries and agencies; part of the livestock series was obtained through an information-disclosure request. Raw source extracts are not redistributed here; the compiled district-year analytic table is provided in `results/` (see `DATA_PROVENANCE.md`). Only aggregated district-year counts are used — **no personally identifiable information**.

## License
MIT (see `LICENSE`). Archived on Zenodo — concept DOI (all versions): https://doi.org/10.5281/zenodo.20725490

## Changelog

### v2.2.0 — extreme-value rule adopted in the principal analysis (2026-09-19)
A post hoc extreme-value screen of the v2.1.x repaired analysis (`qa`-style check of rate, percentage and cost covariates lying more than 5 interquartile ranges beyond the quartiles of the study years) flagged 17 further district-years in 4 covariates as likely aggregation errors, mostly in cities with non-autonomous wards. This is now rule R7 of the principal analysis (`DATA_PROVENANCE.md`), applied in section `[2b]` before the design matrix is built. Two districts, Suwon-si and Seongnam-si, had every study-year value of basic livelihood security recipients flagged and could not be filled, so they leave the analysis under the complete-case rule.
- **Analytic frame changed.** 220 districts / 1,100 district-years → **218 districts / 1,090 district-years** (`results/analysis_dataset_compiled.csv`).
- **Results changed.** Credible covariates 4 → **5**: household groundwater wells, dairy-cattle farms, forest area and fiscal independence remain, and sewer-pipe repair sites becomes credible (0.907, 0.826–0.997). Only forest area is credible in every sensitivity analysis; household groundwater wells and sewer-pipe repair sites are borderline (lost under some neighbourhood graphs, some INLA settings, or when the extreme-value handling changes); dairy-cattle farms, sewer-pipe repair sites and fiscal independence are lost when 2020-2021 is excluded. Residual Moran's I is +0.212 (previously +0.171).
- **The v2.1.x analysis without rule R7 is retained as a sensitivity check**, not withdrawn: it is specification `flagged_values_retained` of `results/sens_extreme_values.csv` (N = 1,100), in which sewer-pipe repair sites is not credible (0.911, 0.829–1.000) and the other 4 associations match the v2.1.x principal model.
- **New file.** `results/extreme_value_flags.csv` (17 rows) now also drives the principal repair log; the R7 counts are not yet added to `results/data_repair_log.csv`, which continues to report R1–R4 only (R7 is applied to the already-repaired frame in section `[2b]`, after the R1–R6 loop).
- Figure 2 and Figures S2–S4 are regenerated from the new principal fit; Figure 1 and Figure S1 (which describe the full 223-district incidence, not the modelled frame) are unchanged.
- Script, `make_figures.R`, `README.md`, `DATA_PROVENANCE.md`, `CITATION.cff` and `.zenodo.json` updated accordingly. No further code-logic change beyond adding rule R7 and the `EXTREME_RULE` environment-variable switch (`EXTREME_RULE=false` reproduces the v2.1.x frame without it).

### v2.1.2 — documentation correction (2026-09-18)
Documentation-only. The description of inpatient medical cost is corrected: the source file stores 0 in 85 districts in 2023 (not 81) and in 4 districts in each of 2020–2022, and in 2024 it stores 0 in 4 districts and is empty in the other 225. One district (Yangyang-gun) meets the structural-zero rule R1b, so 84 + 3 × 4 = 96 values were recoded as missing and 96 + 225 = 321 were missing before filling, as already logged in `results/data_repair_log.csv`. `qa/inpatient_cost_zero_audit.py` and `results/inpatient_cost_zero_audit.csv` are added to document this. `DATA_PROVENANCE.md` now names every rule counted in the repair log. No change to code logic, data or results.

### v2.1.1 — traceability (2026-09-18)
Documentation-only. `results/year_r2_before_after.csv` gains the column `year_r2_after_zero_recovery_only` (the year-R² of each covariate after the zero-coded values had been recovered but before the series-break rule R4 existed; 0.79 for older adults living alone, the value quoted in the manuscript), produced by an optional fourth argument of `qa/year_r2_before_after.py`. Author affiliations completed in the archive metadata. No change to code logic, data or results.

### v2.1.0 — data-repair release (2026-09-18)
- **Data repair.** Values stored as 0 when unavailable, percentages above 100 and redefined series-years are recoded as missing and filled with the most recent earlier available year of the same district (section `[2b]`; `results/data_repair_log.csv`). Functional-form coding no longer turns missing values into 0. Analytic sample 1,112 → 1,100 district-years (220 districts).
- **Results changed.** Credible covariates 9 → 4. The associations of piped water-supply coverage, inpatient medical cost, residential land area, single-person households aged 80–84 years and sewer-pipe repair are no longer credible; household groundwater wells, dairy-cattle farms, forest area and fiscal independence remain. Residual Moran's I is +0.171 (previously +0.053). The earlier results are withdrawn.
- **Gi\* fix.** The three districts without covariate data were entered with zero incidence; Gi\* now uses the observed incidence of all 223 districts (51 hot spots, 56 cold spots).
- **Risk classification** counts only the 220 analysed districts (80 high, 59 low).
- **New analyses** in the single script: post hoc extreme-value screen and sensitivity, M6 without covariates, Pearson-residual Moran's I, oyster-coverage sensitivity, nested random-effect structures, eight INLA settings, glmmTMB cross-fit, repeated fits; hyperparameters reported on the SD scale; covariate dictionary with cut points. The swine/poultry sensitivity now stops instead of mean-imputing missing values silently. `table2_principal_IRR.csv` is written unrounded so that the two-decimal values of the paper are rounded only once.
- Title, target manuscript and archive metadata updated; result files renamed without supplementary-table numbers. `results/core_run.log` and `results/fast_principal_run.log` (v2.0.x runs) are removed from the tree; they remain available in the v2.0.6 archive.

### v2.0.6 — data-provenance alignment (2026-09-18)
Documentation-only. `DATA_PROVENANCE.md` now states that some district-level livestock series were obtained through information-disclosure requests rather than downloaded from a public portal, matching the Data Availability statement and Table 1 of the submitted manuscript. No change to code logic, data or results.

### v2.0.5 — JMIR copyediting alignment (2026-09-18)
Documentation-only release marking the state of the archive when the JMIR package was prepared (that package was not submitted; see v2.1.0). The two archived principal-model fits and the values reported in the manuscript are now listed side by side in this README, matching the reproducibility disclosure added to the manuscript Methods. A script comment referred to "Table S4 **of** Multimedia Appendix 1"; the submitted manuscript uses "**in** Multimedia Appendix N" throughout, and the comment is aligned. No change to code logic, data or results.

### v2.0.4 — JMIR Public Health and Surveillance submission (2026-09-17)
Manuscript title (JMIR title case, year moved before the colon) and target journal updated to the version prepared for *JMIR Public Health and Surveillance*; supplementary material is now referenced as Multimedia Appendix 1 (tables) and Multimedia Appendix 2 (figures). No change to code logic, data or results.

### v2.0.3 — variable-name alignment (2026-09-12)
Reproduced-numbers list in this README now uses the manuscript Table 2 covariate names (water-supply coverage, forest area, fiscal independence, inpatient medical cost) instead of looser paraphrases, and states the seven-of-eight result for fiscal independence explicitly. A script comment described the swine and poultry specificity check as farm *density*; the variables are farm **counts** (source columns `농가수(호)`), and the comment is corrected. No change to code logic, data or results.

### v2.0.2 — One Health title alignment (2026-09-12)
Manuscript title, target journal and archive metadata updated to the version submitted to *One Health*. Supplementary table references in the script header and section comments renumbered to the One Health supplement (model comparison S1, neighbourhood sensitivity S2, Moran's I S3, variable dictionary S4, alternative specifications S6). No change to code logic, data or results.

### v2.0.1 — title alignment (2026-09-06)
Manuscript title and archive metadata updated to the final submitted wording; no change to code, data, or results.

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
