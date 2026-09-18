#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Year-by-year audit of the 0 and missing values of in-district inpatient medical cost in the source file.

Reconciles results/data_repair_log.csv for this covariate: the number of 0 values recoded as missing (rule R1)
equals the 0 values in the study years minus those of the structural-zero districts (rule R1b), and the number
missing before filling equals that count plus the district-years that are absent from the source file.
Writes results/inpatient_cost_zero_audit.csv.

Usage: python3 qa/inpatient_cost_zero_audit.py <health_indicators.parquet> results/inpatient_cost_zero_audit.csv
"""
import sys
import pandas as pd

VAR = "관내진료비_입원"
src, out = sys.argv[1], sys.argv[2]
d = pd.read_parquet(src)
d["region"] = d["region"].astype(str).str.replace(r"\s+", "", regex=True)
d[VAR] = pd.to_numeric(d[VAR], errors="coerce")
g = d.groupby(["region", "year"])[VAR].mean().reset_index()
npos = g.assign(p=g[VAR] > 0).groupby("region")["p"].sum()
struct = set(npos[npos <= 1].index)                       # rule R1b: at most 1 positive year in the whole history
rows = []
for y, s in g.groupby("year"):
    z = s[s[VAR] == 0]
    rows.append(dict(year=int(y), districts=len(s), missing_in_source=int(s[VAR].isna().sum()), zero_in_source=len(z),
                     zero_structural_kept=int(z["region"].isin(struct).sum()),
                     zero_recoded_missing=int((~z["region"].isin(struct)).sum()),
                     median_positive=float(s.loc[s[VAR] > 0, VAR].median()) if (s[VAR] > 0).any() else float("nan")))
r = pd.DataFrame(rows); r["structural_zero_districts"] = len(struct)
r.to_csv(out, index=False)
st = r[(r.year >= 2020) & (r.year <= 2024)]
print(r.to_string(index=False))
print("study years: zeros recoded =", st.zero_recoded_missing.sum(), "| missing before filling =",
      st.zero_recoded_missing.sum() + st.missing_in_source.sum())
