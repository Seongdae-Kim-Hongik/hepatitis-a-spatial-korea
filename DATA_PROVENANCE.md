# Data provenance, repair rules and redistribution note

The analysis uses aggregated district-year counts and administrative indicators; it contains no individual-level records or direct personal identifiers.

## Input domains

- Hepatitis A notifications: Korea Disease Control and Prevention Agency infectious-disease surveillance extracts; district population denominators: resident-registration population from the Korean Statistical Information Service, merged into the surveillance file
- Community-health, demographic, fiscal and healthcare indicators: Korean Statistical Information Service and associated public administrative statistics
- Water supply, sewerage, groundwater use and groundwater-quality testing: Korean environmental and water-related public-data sources
- Land use and district boundaries: Korean land and administrative-boundary public-data sources
- Livestock and shellfish production: Korean agricultural and fisheries statistics (Statistics Korea fishery production, livestock-rearing farm and census-of-agriculture series); some district-level livestock series were obtained through information-disclosure requests rather than downloaded from a public portal

The loader documents the Korean source filenames expected from the archived working directory.

## Known defects of the source files and the v2.1 repair rules

Several source files store unavailable values as 0 rather than as missing, and some series change definition, unit or aggregation in particular years (listed in `README.md`). The repair is applied to the full 2008-2024 history of each series before the study years are selected:

| Rule | Action |
|---|---|
| R1 | A value of 0 in a rate, percentage or cost covariate is recoded as missing |
| R1b | Inpatient medical cost is kept as a structural 0 in districts with at most one positive year in the whole history (no inpatient facility; 1 district, Yangyang-gun); the single stray positive value is set to 0 |
| R2 | In count covariates, a year is recoded as missing when at least 90% of the districts with data are 0 in that year and the other years are mostly non-zero; genuinely sparse covariates (eg, oyster production) are left unchanged |
| R3 | Percentages above 100 are recoded as missing |
| R4 | For rate, percentage and cost covariates, a year whose national median departs by a factor of more than 1.6 from the median of the (up to four) nearest years is treated as non-comparable and recoded as missing (older adults living alone 2024; inpatient medical cost 2019) |
| R5 | Missing values are filled with the value of the most recent earlier available year of the same district, or of the nearest later year when no earlier value exists |
| R6 | Binary, tertile and quartile coding preserves missing values (no silent conversion to 0 or to the lowest class) |
| R7 | Among the rate, percentage and cost covariates of the study years (2020-2024), a value lying more than 5 interquartile ranges below the first or above the third quartile of the district-years of the contiguity graph is treated as an aggregation error, recoded as missing and filled by R5; count covariates, which are genuinely skewed, are not screened. Added after the first corrected (R1-R6) analysis was inspected, so that analysis is retained as a sensitivity check (`results/sens_extreme_values.csv`, specification `flagged_values_retained`) |

`results/data_repair_log.csv` gives, per source variable, the number of study-year values recoded as missing under R1 or R2 (one column; R1 applies to rate, percentage and cost covariates and R2 to count covariates), R3 and R4 and the number of missing values before and after filling; `results/extreme_value_flags.csv` gives the 17 district-year values recoded under R7. Districts that still have missing values after filling (Jeju-si and Seogwipo-si: no fiscal statistics of their own; Jeju-si and Busanjin-gu: their records in the compiled sewerage files contain no non-zero value in any field or year, which indicates a linkage failure during compilation of the source rather than true absence, so the original values could not be recovered; Suwon-si and Seongnam-si: every study-year value of basic livelihood security recipients was flagged as an aggregation error by R7 and could not be filled) are excluded under a complete-case rule, leaving 218 of the 223 graph districts (1,090 district-years) in the principal analysis.

Because only one year of data exists, basic livelihood security recipients and wastewater-discharging facilities (2020) and groundwater-quality tests (2023) are constant over the study years after filling. A further limitation is handled by sensitivity analysis: oyster production is available only for 2020, 2022 and 2023, and the 2020 file omits the major producing districts Tongyeong, Goseong and Yeosu (`results/sens_oyster_coverage.csv`).

## Compiled analytic file

`results/analysis_dataset_compiled.csv` is the exact 1,090-row district-year table (218 districts, 2020-2024) used by the archived run, after repair (rule R7 applied). It is included to support computational reproduction of the reported model. The compiled file shipped with v2.0.x (1,112 rows) contained the unrepaired values and must not be used.

The compiled file contains only district-year aggregates derived from publicly released statistics. Users who redistribute it should check the terms of the source portals.

## Interpretation

The final 27-covariate specification was developed through exploratory model building. The compiled dataset and code reproduce that selected model; they do not convert the analysis into a prospectively specified confirmatory study.
