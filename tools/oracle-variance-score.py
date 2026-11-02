#!/usr/bin/env python3
"""Score the reading and predictions in validation/VARIANCE_CHECK_PLAN.md from saved results, in plain Python, and check the control against the study.
Usage: python3 -I tools/oracle-variance-score.py [OUTPUT_FILE] [DIRECTORY]"""
import csv
import gzip
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
FOLDER = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/"validation/variance-check"
OUT = []


def say(text=""):
    OUT.append(text)


def truthy(v):
    return v in ("TRUE", "True", "1")


with gzip.open(FOLDER/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    records = list(csv.DictReader(stream))
groups = {}
for r in records:
    groups.setdefault((r["level"], r["analysis"], float(r["fraction"])), []).append(r)
with gzip.open(ROOT/"validation/clustered-fix-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    saved = {(r["analysis"], float(r["fraction"]), r["replicate"]): r for r in csv.DictReader(stream) if r["method"] == "oracle"}
checked = 0
largest = 0.0
for r in records:
    s = saved[(r["analysis"], float(r["fraction"]), r["replicate"])]
    largest = max(largest, abs(float(r["estimate"])-float(s["estimate"])), abs(float(r["se_cr1s"])-float(s["se"])))
    checked += 1
assert largest < 1e-8, largest
METHODS = ["cr1s", "cr2", "cv3", "wcr"]


def cover(g, m):
    return sum(float(r["wcr_p"]) > .05 for r in g)/len(g) if m == "wcr" else sum(truthy(r[f"cover_{m}"]) for r in g)/len(g)


def width(g, m):
    return sum(float(r[f"width_{m}"]) for r in g)/len(g)


def se_ratio(g, m):
    return (sum(float(r[f"se_{m}"]) for r in g)/len(g))/statistics.stdev(float(r["estimate"]) for r in g)


say("Scored from validation/variance-check by tools/oracle-variance-score.py")
say(f"Requests: {len(records)}. Control: {checked} estimates and cr1s standard errors match the study's oracle rows, largest difference {largest:.2e}.")
say()
real = {k: g for k, g in groups.items() if k[0] != "simulation"}
sim = {k: g for k, g in groups.items() if k[0] == "simulation"}
assert len(real) == 11 and len(sim) == 38
say("REAL ITEMS: coverage, percent; mean width; mean standard error over empirical SD; Satterthwaite df; largest partial leverage")
for k, g in real.items():
    say(f"  {k[1]:12s} coverage " + ", ".join(f"{m} {100*cover(g, m):.1f}" for m in METHODS) + "; width " + ", ".join(f"{m} {width(g, m):.3f}" for m in METHODS[:3]) + "; SE ratio " + ", ".join(f"{m} {se_ratio(g, m):.2f}" for m in METHODS[:3]) + f"; df {statistics.mean(float(r['cr2_df']) for r in g):.1f} of {g[0]['clusters']} clusters; leverage {statistics.mean(float(r['leverage']) for r in g):.2f}")
say()
say("SIMULATED CELLS: coverage, percent; SE ratio; Satterthwaite df")
for k, g in sim.items():
    say(f"  {k[1]:22s} {int(k[2]*100)}% coverage " + ", ".join(f"{m} {100*cover(g, m):.1f}" for m in METHODS) + "; SE ratio " + ", ".join(f"{m} {se_ratio(g, m):.2f}" for m in METHODS[:3]) + f"; df {statistics.mean(float(r['cr2_df']) for r in g):.1f}")
say()
say("MEANS BY SETTING: coverage, percent")
SETTINGS = [("clustered", lambda a: a.endswith(("shift_low", "shift", "under", "outlier"))), ("rare counts", lambda a: a.endswith("rare")), ("row-only errors", lambda a: a.endswith("row_only") or a == "control_one_row")]
for name, test in SETTINGS:
    keys = [k for k in sim if test(k[1])]
    say(f"  {name:16s} {len(keys):2d} cells: " + ", ".join(f"{m} {100*statistics.mean(cover(sim[k], m) for k in keys):.1f}" for m in METHODS))
say("  Real items, mean distance from 95%, points: " + ", ".join(f"{m} {100*statistics.mean(abs(cover(g, m)-.95) for g in real.values()):.2f}" for m in METHODS))
say("  Real items, lowest coverage: " + ", ".join(f"{m} {100*min(cover(g, m) for g in real.values()):.1f}" for m in METHODS))
say()
say("READING FIXED IN ADVANCE")
qualifying = []
for m in METHODS:
    low = sum(cover(g, m) >= .92 for g in real.values())
    high = all(cover(g, m) <= .98 for g in real.values())
    inside = sum(.92 <= cover(g, m) <= .98 for g in sim.values())
    ok = low >= 9 and high and inside >= .8*38
    if ok:
        qualifying.append(m)
    say(f"  {m}: real items at least 92% on {low} of 11, at most 98% on all {high}; simulated cells within limits {inside} of 38; candidate {ok}")
say("  CANDIDATES: " + (", ".join(qualifying) if qualifying else "none"))
say()
say("PREDICTIONS")
lbw = next(g for k, g in real.items() if k[1] == "lbw")
say(f"  1. cr2 and cv3 raise lbw to at least 92%: {cover(lbw, 'cr2') >= .92 and cover(lbw, 'cv3') >= .92} (cr2 {100*cover(lbw, 'cr2'):.1f}, cv3 {100*cover(lbw, 'cv3'):.1f})")
rare = [g for k, g in sim.items() if k[1].endswith("rare")]
means = {m: statistics.mean(cover(g, m) for g in rare) for m in METHODS}
say(f"  2. No variance reaches a mean of 92% on rare cells: {all(v < .92 for v in means.values())} (" + ", ".join(f"{m} {100*v:.1f}" for m, v in means.items()) + ")")
order = sum(width(g, "cv3") > width(g, "cr2") > width(g, "cr1s") for g in real.values())
say(f"  3. cv3 wider than cr2 wider than cr1s on all 11: {order == 11} ({order} of 11)")
near = sum(abs(cover(g, "wcr")-cover(g, "cv3")) <= .02 for g in real.values())
say(f"  4. wcr within 2 points of cv3 on at least nine: {near >= 9} ({near} of 11)")
under = sum(se_ratio(g, "cr1s") < 1 for g in real.values())
say(f"  5. cr1s SE ratio below 1 on at least nine: {under >= 9} ({under} of 11)")
say(f"  6. At least one candidate: {bool(qualifying)} ({len(qualifying)})")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
