#!/usr/bin/env python3
"""Readings of the corrected-variance test made after its results. Post hoc: they enter no rule. Plain Python, from saved results.
Usage: python3 -I tools/variance-mi-posthoc.py [OUTPUT_FILE]"""
import csv
import gzip
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
OUT = []


def say(text=""):
    OUT.append(text)


groups = {}
with gzip.open(ROOT/"validation/variance-mi-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            groups.setdefault((r["analysis"], float(r["fraction"]), r["imputer"]), []).append(r)


def mean(key, field):
    return statistics.mean(float(r[field]) for r in groups[key])


say("Post hoc: made after the results, enters no rule")
say()
say("IMPUTED-COVARIATE CELLS: relative bias of the slope, percent; RMSE; RMSE ratio against linear MI; coverage under cv3")
for errors in ("covariate_shift", "covariate_under", "covariate_row_only"):
    for f in (.25, .5):
        rows = {}
        for imp in ("oracle", "linear_mi", "random_intercept"):
            e = [float(r["estimate"])-.5 for r in groups[(errors, f, imp)]]
            rows[imp] = (100*statistics.mean(e)/.5, (sum(x*x for x in e)/len(e))**.5, 100*mean((errors, f, imp), "cover_cv3"))
        say(f"  {errors:20s} {int(f*100)}%: " + "; ".join(f"{imp} bias {rows[imp][0]:+.1f}, RMSE {rows[imp][1]:.4f}, coverage {rows[imp][2]:.1f}" for imp in rows) + f"; RMSE ratio {rows['random_intercept'][1]/rows['linear_mi'][1]:.3f}")
say()
say("REAL ITEMS AT 25%, 2,000 COHORTS: items at 92% or more, of 11; lowest coverage; mean width ratio against cr1s")
for imp in ("oracle", "linear_mi", "random_intercept"):
    for v in ("cr1s", "cr2", "cv3"):
        values = {a: mean((a, .25, imp), f"cover_{v}") for a in ITEMS}
        ratios = [mean((a, .25, imp), f"width_{v}")/mean((a, .25, imp), "width_cr1s") for a in ITEMS]
        say(f"  {imp:18s} {v}: {sum(x >= .92 for x in values.values())} items; lowest {100*min(values.values()):.1f}% on {min(values, key=values.get)}; width ratio {100*(min(ratios)-1):+.0f}% to {100*(max(ratios)-1):+.0f}%")
say()
say("MONTE CARLO ERROR OF COVERAGE, points: at 92% with 2,000 cohorts " + f"{100*(.92*.08/2000)**.5:.2f}" + ", at 92% with 500 cohorts " + f"{100*(.92*.08/500)**.5:.2f}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
