#!/usr/bin/env python3
"""Description only for validation/WITHIN_BETWEEN_PLAN.md: the error's intraclass correlation in SNAP, the share of reviews with an error, and cohort-bootstrap intervals for the differences behind the rule.
Usage: python3 -I tools/within-between-posthoc.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import random
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/within-between-study"
DRAWS = 500
OUT = []


def say(text=""):
    OUT.append(text)


with (ROOT/"validation/within-between/snap-pool.csv").open(newline="", encoding="utf-8") as stream:
    pool = list(csv.DictReader(stream))
errors = {}
for r in pool:
    errors.setdefault(r["state"], []).append(float(r["recorded"])-float(r["reference"]))
every = [e for v in errors.values() for e in v]
n, g = len(every), len(errors)
grand = statistics.mean(every)
msb = sum(len(v)*(statistics.mean(v)-grand)**2 for v in errors.values())/(g-1)
msw = sum((e-statistics.mean(v))**2 for v in errors.values() for e in v)/(n-g)
k = (n-sum(len(v)**2 for v in errors.values())/n)/(g-1)
say("SNAP POOL")
say(f"  reviews {n}; states {g}; share with an error {100*sum(e != 0 for e in every)/n:.1f}%; mean error {100*grand:.2f} dollars; standard deviation of error {100*statistics.pstdev(every):.2f} dollars")
say(f"  intraclass correlation of the error by state, by moments: {(msb-msw)/(msb+(k-1)*msw):.3f}")
say()

cells = {}
with gzip.open(STUDY/"deciding.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["method"] in ("random_intercept", "within_between") and r["status"] == "ok":
            cells.setdefault((r["source"], r["role"], r["analysis"], r["fraction"]), {}).setdefault(int(r["replicate"]), {})[r["method"]] = (float(r["estimate"])-float(r["target"]), float(r["cover"]))


def measure(picked):
    """The difference in distance, and the log RMSE ratio, for one cell's cohorts."""
    a2 = sum(p["within_between"][0]**2 for p in picked)
    b2 = sum(p["random_intercept"][0]**2 for p in picked)
    ca = 100*sum(p["within_between"][1] for p in picked)/len(picked)
    cb = 100*sum(p["random_intercept"][1] for p in picked)/len(picked)
    return abs(ca-95)-abs(cb-95), .5*math.log(a2/b2)


random.seed(20261102)
say(f"COHORT BOOTSTRAP, {DRAWS} draws: difference in distance, candidate minus option, then RMSE ratio, each with a 95% interval")
for source in ("snap", "fresh"):
    for role in ("outcome", "predictor"):
        chosen = [[p for p in v.values() if len(p) == 2] for k, v in cells.items() if k[0] == source and k[1] == role]
        if not chosen:
            continue
        point = [measure(c) for c in chosen]
        distance, ratio = [], []
        for _ in range(DRAWS):
            drawn = [measure(random.choices(c, k=len(c))) for c in chosen]
            distance.append(statistics.mean(d[0] for d in drawn))
            ratio.append(math.exp(statistics.mean(d[1] for d in drawn)))
        distance.sort(); ratio.sort()
        low, high = int(.025*DRAWS), int(.975*DRAWS)-1
        say(f"  {source} {role}: distance {statistics.mean(p[0] for p in point):+.2f} ({distance[low]:+.2f} to {distance[high]:+.2f}); ratio {math.exp(statistics.mean(p[1] for p in point)):.3f} ({ratio[low]:.3f} to {ratio[high]:.3f})")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
