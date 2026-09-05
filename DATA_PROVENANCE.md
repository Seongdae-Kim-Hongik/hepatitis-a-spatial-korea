# Data provenance and redistribution note

The analysis uses aggregated district-year counts and administrative indicators; it contains no individual-level records or direct personal identifiers.

## Input domains

- Hepatitis A notifications and district population denominators: Korea Disease Control and Prevention Agency infectious-disease surveillance extracts
- Community-health, demographic, fiscal and healthcare indicators: Korean Statistical Information Service and associated public administrative statistics
- Water supply, sewerage, groundwater use and groundwater-quality testing: Korean environmental and water-related public-data sources
- Land use and district boundaries: Korean land and administrative-boundary public-data sources
- Livestock and shellfish production: Korean agricultural and fisheries public statistics

The corrected loader documents the Korean source filenames expected from the archived working directory. Biennially available administrative measures are carried to adjacent years as described in the manuscript.

## Compiled analytic file

`results/analysis_dataset_compiled.csv` is the exact 1,112-row district-year table used by the verified principal-model run. It is included to support computational reproduction of the reported model.

Before public release, the depositing authors must verify that redistribution of each derived field complies with the terms of its source portal. If any source terms prohibit redistribution, remove the compiled file from the public release and retain the acquisition/assembly code, field dictionary, source URLs and access dates so authorised users can reconstruct it.

## Interpretation

The final 27-covariate specification was developed through exploratory model building. The compiled dataset and code reproduce that selected model; they do not convert the analysis into a prospectively specified confirmatory study.
