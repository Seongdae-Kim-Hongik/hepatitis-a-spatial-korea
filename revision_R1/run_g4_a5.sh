#!/bin/zsh
# Refits the revised principal model (A5) on the compiled frames produced by thresholds/run_thresholds.sh.
cd "$(dirname "$0")/.."
for r in zs080 zs095 mf140 mf200 interp; do
  echo "=== $r"; env HAV_FRAME="$PWD/revision_R1/thresholds/out_$r/analysis_dataset_compiled.csv" HAV_OUT="$PWD/revision_R1/results_thresholds/A5_$r" HAV_ARMS=A5 Rscript revision_R1/HAV_revision_R1.R > "revision_R1/results_thresholds/A5_$r.log" 2>&1; echo "exit $?"
done
