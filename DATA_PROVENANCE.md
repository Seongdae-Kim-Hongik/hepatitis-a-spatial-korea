# Data provenance, repair rules and redistribution note

The analysis uses aggregated district-year counts and administrative indicators; it contains no individual-level records or direct personal identifiers.

## Input domains

- Hepatitis A notifications: Korea Disease Control and Prevention Agency infectious-disease surveillance extracts; district population denominators: resident-registration population from the Korean Statistical Information Service, merged into the surveillance file
- Community-health, demographic, fiscal and healthcare indicators: Korean Statistical Information Service and associated public administrative statistics
- Water supply, sewerage, groundwater use and groundwater-quality testing: Korean environmental and water-related public-data sources
- Land use and district boundaries: Korean land and administrative-boundary public-data sources
- Livestock and shellfish production: Korean agricultural and fisheries statistics (Statistics Korea fishery production, livestock-rearing farm and census-of-agriculture series); some district-level livestock series were obtained through information-disclosure requests rather than downloaded from a public portal

The loader documents the Korean source filenames expected from the archived working directory.

## Known defects of the source files and the repair rules (R0-R7)

Several source files store unavailable values as 0 rather than as missing, and some series change definition, unit or aggregation in particular years (listed in `README.md`). The repair is applied to the full 2008-2024 history of each series before the study years are selected:

| Rule | Action |
|---|---|
| R0 | Busanjin-gu sewerage values, lost when the source extracts were compiled (the district name "부산진구" had its leading "부산" stripped as a province abbreviation), are restored from an earlier compilation of the same Ministry of Environment statistics (`data_patches/busanjin_sewerage_restore.csv`; checked against the 2023 sewerage yearbook and KOSIS DT_1YL20751). Applied before R1-R7 and only where the source holds no non-zero value |
| R1 | A value of 0 in a rate, percentage or cost covariate is recoded as missing |
| R1b | Inpatient medical cost is kept as a structural 0 in districts with at most one positive year in the whole history (no inpatient facility; 1 district, Yangyang-gun); the single stray positive value is set to 0 |
| R2 | In count covariates, a year is recoded as missing when at least 90% of the districts with data are 0 in that year and the other years are mostly non-zero; genuinely sparse covariates (eg, oyster production) are left unchanged |
| R3 | Percentages above 100 are recoded as missing |
| R4 | For rate, percentage and cost covariates, a year whose national median departs by a factor of more than 1.6 from the median of the (up to four) nearest years is treated as non-comparable and recoded as missing (older adults living alone 2024; inpatient medical cost 2019) |
| R5 | Missing values are filled with the value of the most recent earlier available year of the same district, or of the nearest later year when no earlier value exists |
| R6 | Binary, tertile and quartile coding preserves missing values (no silent conversion to 0 or to the lowest class) |
| R7 | Among the rate, percentage and cost covariates of the study years (2020-2024), a value lying more than 5 interquartile ranges below the first or above the third quartile of the district-years of the contiguity graph is treated as an aggregation error, recoded as missing and filled by R5; count covariates, which are genuinely skewed, are not screened. Added after the first corrected (R1-R6) analysis was inspected, so the analysis without R7 (R0-R6 applied) is retained as a sensitivity check (`results/sens_extreme_values.csv`, specification `flagged_values_retained`) |

`results/data_repair_log.csv` gives, per source variable, the number of study-year values recoded as missing under R1 or R2 (one column; R1 applies to rate, percentage and cost covariates and R2 to count covariates), R3 and R4 and the number of missing values before and after filling; `results/extreme_value_flags.csv` gives the 17 district-year values recoded under R7. Districts that still have missing values after filling (Jeju-si and Seogwipo-si: no fiscal statistics of their own; Jeju-si also has no sewerage values in the compiled files; Suwon-si and Seongnam-si: every study-year value of basic livelihood security recipients was flagged as an aggregation error by R7 and could not be filled) are excluded under a complete-case rule, leaving 219 of the 223 graph districts (1,095 district-years) in the principal analysis.

Because only one year of data exists, basic livelihood security recipients and wastewater-discharging facilities (2020) and groundwater-quality tests (2023) are constant over the study years after filling. A further limitation is handled by sensitivity analysis: oyster production is available only for 2020, 2022 and 2023, and the 2020 file omits the major producing districts Tongyeong, Goseong and Yeosu (`results/sens_oyster_coverage.csv`).

## Input corrections made at revision R1 (v2.4.0)

Three further problems were found while the covariates were being rescaled for the revision; they are corrected in `revision_R1/data/analysis_dataset_compiled_v2.4.0.csv`, which the revision script reads, and they do not alter `results/` (v2.3.3), which is kept as the reviewed analysis.

| Problem | Correction |
|---|---|
| The district-name truncation that lost Busanjin-gu's sewerage values (rule R0) had also affected the land-use (forest, paddy-field, residential area) and household groundwater-well files, so Busanjin-gu carried 0 in those series; the land-use values of Cheonan-si were 0 in every year | Newer compilations of the same source files supply the values (35 district-year cells in 2 districts); the sewerage values of Busanjin-gu are still restored by R0 |
| The area field carried in the compiled frame (`총면적(㎢)`, from the sewerage statistics) is about twice the cadastral area in most districts (eg, Hongcheon-gun 3,641 vs 1,815 km²); it was never a covariate | District land area for every share and density is computed from the official district boundary file (Korea 2000 Unified Coordinate System, EPSG:5179); `revision_R1/results/area_check_sewer_vs_shapefile.csv` lists both |
| In the reviewed frame the notifications of Michuhol-gu (Incheon Nam-gu) and Namdong-gu were attached to each other's polygon through a label inconsistency in an earlier copy of the boundary file | The district index is derived from the boundary file by district code and geometry; the current boundary file is consistent (code 23030, 24.6 km², centroid 126.665E 37.453N = Michuhol-gu; code 23050, 57.8 km² = Namdong-gu) |

Refitting the reviewed specification on the corrected inputs (arm A0 of `revision_R1/HAV_revision_R1.R`) reproduces its 5 credible associations to the second decimal.

## Age-sex-standardized expected counts (v2.4.0, sex added in v2.4.1)

`revision_R1/data/ST_HepatitisA_agesexcase.csv` holds the KDCA Infectious Disease Portal counts of notified hepatitis A cases by district, year (2020-2024), sex and 10-year age band (0-9 to 70+; the age-only file `ST_HepatitisA_agecase.csv` is its sum over sex); `revision_R1/data/ST_agesexpop_2020_2024.csv` the resident-registration population by district, year, sex and the same bands (December of each year, built from the Ministry of the Interior and Safety monthly registration file by `build_agesexpop.py`; its sum over sex equals `ST_agepop.csv` cell by cell). The portal's own population field was not used because it repeats the 2022 population of Changwon-si for 2023. The national rate of each sex-age stratum and year is the sum of district counts divided by the sum of district populations over all 229 districts; the expected count of a district-year is the sum over the 16 strata of its population multiplied by the national rate (`revision_R1/results/national_agesex_rates_2020_2024.csv`; `national_age_rates_2020_2024.csv` gives the both-sexes rates used for the age-only comparison). Cases with unknown age are not in the age-band cross-tabulation, so the age-band totals are slightly below the district totals. `revision_R1/data/physicians_count.csv` (KOSIS health care resources) supplies the physicians-per-1000 ascertainment proxy.

## Compiled analytic file

`results/analysis_dataset_compiled.csv` is the exact 1,095-row district-year table (219 districts, 2020-2024) used by the archived run, after repair (rules R0 and R7 applied). It is included to support computational reproduction of the reported model. The compiled file shipped with v2.0.x (1,112 rows) contained the unrepaired values and must not be used.

The compiled file contains only district-year aggregates derived from publicly released statistics. Users who redistribute it should check the terms of the source portals.

## Interpretation

The final 27-covariate specification was developed through exploratory model building. The compiled dataset and code reproduce that selected model; they do not convert the analysis into a prospectively specified confirmatory study.
