#!/bin/zsh
# Reruns the complete data pipeline (HAV_spatial_sens_thresholds.R = HAV_spatial_reproducible.R with the rule-2 share,
# rule-4 factor and fill method exposed as environment variables) under alternative screening thresholds.
# Raw inputs (KDCA extract, KOSIS and ministry files, boundary file) must be under $HAV_DATA_DIR (not redistributed).
cd "$(dirname "$0")"
run() { name=$1; shift; echo "=== $name $(date +%H:%M:%S)"; env "$@" HAV_OUTPUT_DIR="$PWD/out_$name" FAST_PRINCIPAL=true Rscript HAV_spatial_sens_thresholds.R > "log_$name.txt" 2>&1; echo "exit $?"; }
run base   ZERO_SHARE=0.9  MEDIAN_FACTOR=1.6
run zs080  ZERO_SHARE=0.8  MEDIAN_FACTOR=1.6
run zs095  ZERO_SHARE=0.95 MEDIAN_FACTOR=1.6
run mf140  ZERO_SHARE=0.9  MEDIAN_FACTOR=1.4
run mf200  ZERO_SHARE=0.9  MEDIAN_FACTOR=2.0
run interp ZERO_SHARE=0.9  MEDIAN_FACTOR=1.6 FILL_METHOD=interp
