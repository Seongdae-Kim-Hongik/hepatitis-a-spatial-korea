#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Year-R2 screen for district-year covariate panels.

For every numeric column, reports the share of variance explained by calendar year
(between-year sum of squares / total sum of squares). A covariate that mainly encodes
WHICH YEARS HAVE DATA (missing values stored as 0, a survey run only in some years, or a
series redefined in one year) has a high value; a genuine district characteristic has a
low value. In this project the screen exposed both the zero-coded missing years and the
2024 redefinition of the "older adults living alone" series, which a check of the share
of zeros alone had missed.

Usage:  python3 qa/year_r2_screen.py results/analysis_dataset_compiled.csv [--year year] [--top 15] [--flag 0.30]
No third-party packages are required.
"""
import argparse, collections, csv

ap = argparse.ArgumentParser()
ap.add_argument("panel"); ap.add_argument("--year", default="year"); ap.add_argument("--top", type=int, default=15)
ap.add_argument("--flag", type=float, default=0.30, help="flag covariates whose year R2 exceeds this value")
ap.add_argument("--skip", default="cases,population,idarea,idtime,idarea_time,rate_100k", help="comma-separated columns to skip")
a = ap.parse_args()
rows = list(csv.DictReader(open(a.panel, encoding="utf-8-sig")))
skip = set(a.skip.split(",")) | {a.year}

def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None

out = []
for col in rows[0]:
    if col in skip: continue
    by = collections.defaultdict(list)
    for r in rows:
        v = num(r[col])
        if v is not None: by[r[a.year]].append(v)
    vals = [v for vs in by.values() for v in vs]
    if len(vals) < 20 or len(by) < 2: continue
    m = sum(vals) / len(vals); sst = sum((v - m) ** 2 for v in vals)
    if sst == 0: continue
    ssb = sum(len(vs) * (sum(vs) / len(vs) - m) ** 2 for vs in by.values())
    out.append((ssb / sst, col, {y: sum(vs) / len(vs) for y, vs in sorted(by.items())}))

out.sort(reverse=True)
print(f"{len(rows):,} rows | {len(out)} numeric columns screened | flag threshold {a.flag}")
for r2, col, means in out[:a.top]:
    flag = "FLAG" if r2 > a.flag else "    "
    print(f"{flag} {r2:6.3f}  {col:32s} " + "  ".join(f"{y}:{v:.3g}" for y, v in means.items()))
n = sum(1 for r2, _, _ in out if r2 > a.flag)
print(f"{n} column(s) above the threshold" + ("" if n else " - no covariate is dominated by calendar year"))
