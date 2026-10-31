# revision_R1 (v2.4.1)

Additional analyses for revision R1 of the manuscript (JMIR Public Health and Surveillance #112448). See the top section of the repository README for the summary and the run commands.

| Path | Content |
|---|---|
| `HAV_revision_R1.R` | All arms A0-A16 (A0 reviewed specification on corrected inputs; A1 offset without +1; A2 age-standardized offset; A3 proportion aged 20-49; A4 rescaled + density; **A5 revised principal model**; A6/A7 uniform forms; A8 parsimonious P1/P2; A9 restricted spatial regression; A10 zero-inflated NB; A11 nested M1-M6 with residual Moran's I; A12 eight graphs; A13 spatially varying coefficients; A14 physicians per 1000; A15 2022-2024; A16 repeated fits). `HAV_FRAME`, `HAV_ARMS`, `HAV_OUT` select an alternative input frame, a subset of arms, and the output folder |
| `figures_R1.R` | Figure 2 (zero-boundary classes), Figures S3 and S4, and the district-effect table from the A5 fit |
| `data/` | `analysis_dataset_compiled_v2.4.0.csv` (1,095-row frame of the corrected pipeline, 219 districts), `covariate_dictionary_v2.4.0.csv`, `ST_HepatitisA_agesexcase.csv` and `ST_agesexpop_2020_2024.csv` (KDCA district x sex x age-band notifications and resident-registration populations, 2020-2024; `build_agesexpop.py` rebuilds the population file), `ST_HepatitisA_agecase.csv` and `ST_agepop.csv` (age-only, used when `HAV_STD=age`), `physicians_count.csv` |
| `results/` | `coefficients_all_arms.csv` (IRR, 95% CrI, credible flag for every term of every arm), `fit_all_arms.csv` (DIC, WAIC, pD, residual Moran's I on count and Poisson-scaled residuals, two-sided P, fixed-Queen reference, high/low-risk counts), `credible_set_persistence.csv`, `A5_district_effects.csv`, `A5_hyperparameters.csv`, `A5_cpo_pit.csv`, `A13_SVC_*_local_slopes.csv` (total local slopes as INLA linear combinations), `repeated_fits_A0_A5.csv`, `national_agesex_rates_2020_2024.csv`, `national_age_rates_2020_2024.csv`, `rescaled_covariate_dictionary_full.csv` (units, forms, cut points), `analysis_frame_R1.csv` (frame with the rescaled covariates and expected counts), `area_check_sewer_vs_shapefile.csv`, `run.log` |
| `thresholds/` | `HAV_spatial_sens_thresholds.R` = the main pipeline with the rule-2 share (`ZERO_SHARE`), rule-4 factor (`MEDIAN_FACTOR`) and fill method (`FILL_METHOD=nearest|interp`) as environment variables; `run_thresholds.sh` |
| `results_thresholds/` | `submitted_table2_*.csv` and `submitted_table_S11_thresholds_interp.csv` (reviewed specification under each threshold), `pipeline_<run>_data_repair_log.csv` / `_extreme_value_flags.csv` / `_table2_principal_IRR.csv` (pipeline outputs per run), `A5_<run>_coefficients_all_arms.csv` / `_fit_all_arms.csv` (revised principal model refitted on each run's frame) |
| `figures/` | `Figure2.png`, `FigureS3.png`, `FigureS4.png` |

The A5 fit object (`A5_fit.rds`, 7 MB) is not committed; `HAV_revision_R1.R` writes it to `results/` when run.
