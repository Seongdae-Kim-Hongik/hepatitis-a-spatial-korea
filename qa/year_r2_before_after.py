#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Year-R2 of every model covariate before and after the v2.1 data repair.

before = compiled dataset archived with v2.0.6 (unrepaired; https://doi.org/10.5281/zenodo.22823390)
after  = results/analysis_dataset_compiled.csv of this release
intermediate (optional 4th argument) = dataset after the zero-coded values had been recovered but before the
         series-break rule R4 was added (not archived); it shows that the 2024 redefinition of the series for
         older adults living alone was still visible at that stage.
The statistic is the share of the variance of the transformed, standardized covariate (the *_z column that
enters the model) explained by calendar year. Writes results/year_r2_before_after.csv.

Usage: python3 qa/year_r2_before_after.py <v2.0.6 analysis_dataset_compiled.csv> results/analysis_dataset_compiled.csv results/year_r2_before_after.csv [<intermediate dataset>]
"""
import collections, csv, sys

def year_r2(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig"))); out = {}
    for col in rows[0]:
        if not col.endswith("_z"): continue
        by = collections.defaultdict(list)
        for r in rows:
            try: by[r["year"]].append(float(r[col]))
            except ValueError: pass
        vals = [v for vs in by.values() for v in vs]; m = sum(vals) / len(vals)
        sst = sum((v - m) ** 2 for v in vals)
        ssb = sum(len(vs) * (sum(vs) / len(vs) - m) ** 2 for vs in by.values())
        out[col[:-2]] = ssb / sst if sst else 0.0
    return out

before, after = year_r2(sys.argv[1]), year_r2(sys.argv[2])
inter = year_r2(sys.argv[4]) if len(sys.argv) > 4 else {}
with open(sys.argv[3], "w", newline="") as f:
    w = csv.writer(f); w.writerow(["covariate", "year_r2_before_repair", "year_r2_after_zero_recovery_only", "year_r2_after_repair"])
    for c in after: w.writerow([c, f"{before.get(c, float('nan')):.4f}", f"{inter[c]:.4f}" if c in inter else "", f"{after[c]:.4f}"])
for c in sorted(after, key=lambda k: -before.get(k, 0))[:10]: print(f"{c:18s} before {before.get(c, float('nan')):.3f}  after {after[c]:.3f}")
