#!/usr/bin/env python3
"""Readings of the test of independently written setups and of the rare-count test, made after their results. Post hoc: they enter no rule. Plain Python, from saved results.
Usage: python3 -I tools/independent-setups-posthoc.py [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = []


def say(text=""):
    OUT.append(text)


ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
setups = {}
with gzip.open(ROOT/"validation/independent-setups-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            setups.setdefault((r["analysis"], float(r["fraction"]), r["method"]), {})[int(r["replicate"])] = (float(r["estimate"])-float(r["target"]), float(r["cover"]), float(r["width"]))
earlier = {}
with gzip.open(ROOT/"validation/variance-mi-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            earlier.setdefault((r["analysis"], float(r["fraction"]), r["imputer"]), {})[int(r["replicate"])] = (float(r["estimate"])-float(r["target"]), float(r["cover_cv3"]), float(r["width_cv3"]))
get = lambda a, f, m: setups[(a, f, m)] if m.startswith("ind_") else earlier[(a, f, m)]
cover = lambda a, f, m: 100*statistics.mean(v[1] for v in get(a, f, m).values())
width = lambda a, f, m: statistics.mean(v[2] for v in get(a, f, m).values())
say("Post hoc: made after the results, enters no rule")
say()
say("ind_mice ABOVE 98 ON REAL ITEMS, with the recount's coverage; this alone failed N3 and N4")
for f in (.25, .5):
    above = [(a, cover(a, f, "ind_mice"), cover(a, f, "oracle")) for a in ITEMS if cover(a, f, "ind_mice") > 98]
    say(f"  {int(f*100)}%: "+"; ".join(f"{a} {x:.2f}, recount {y:.2f}" for a, x, y in above))
say()
say("MEAN WIDTH OVER THE OPTION'S, simulated cells: rare cells, other audited-outcome cells")
rare = [(c, f) for c in ("rwanda_rare", "pairs_rare", "none_rare") for f in (.25, .5)]
other = [(f"{d}_{e}", f) for d in ("rwanda", "pairs", "none") for e in ("shift_low", "shift", "under", "outlier", "row_only") for f in (.25, .5)]
for m in ("ind_mice", "ind_survey", "ind_design"):
    say(f"  {m}: {statistics.mean(width(a, f, m)/width(a, f, 'random_intercept') for a, f in rare):.3f}, {statistics.mean(width(a, f, m)/width(a, f, 'random_intercept') for a, f in other):.3f}")
say()
say("MEAN COVERAGE ON THE SIX RARE CELLS")
for m in ("oracle", "linear_mi", "random_intercept", "ind_mice", "ind_survey", "ind_design"):
    say(f"  {m}: {statistics.mean(cover(a, f, m) for a, f in rare):.2f}")
say()
say("nb_mi IN THE RARE-COUNT TEST, median absolute error over linear MI's, and the largest absolute error, to separate typical behavior from blow-ups")
nb = {}
with gzip.open(ROOT/"validation/rare-counts-study/part2.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            nb.setdefault((r["analysis"], float(r["fraction"]), r["method"]), {})[int(r["replicate"])] = float(r["estimate"])-float(r["target"])
for a, f, m in sorted(nb):
    if m != "nb_mi":
        continue
    base = earlier[(a, f, "linear_mi")]
    both = [r for r in nb[(a, f, m)] if r in base]
    say(f"  {a} {int(f*100)}%: {statistics.median(abs(nb[(a, f, m)][r]) for r in both)/statistics.median(abs(base[r][0]) for r in both):.3f}; {max(abs(nb[(a, f, m)][r]) for r in both):.3g}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
