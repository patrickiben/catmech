#!/usr/bin/env python3
"""Post hoc, written after the results of validation/WITHIN_BETWEEN_PLAN.md were read, and entering no rule: prediction 4 on the cells where ind_mice gave intervals, and the deciding classes split by audit size.
Usage: python3 -I tools/within-between-shared-cells.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/within-between-study"
OPTION, CAND = "random_intercept", "within_between"
OUT = []


def say(text=""):
    OUT.append(text)


store = {}
with gzip.open(STUDY/"deciding.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["method"] in (OPTION, CAND, "ind_mice") and r["status"] == "ok":
            store.setdefault((r["source"], r["role"], r["analysis"], r["fraction"]), {}).setdefault(r["method"], {})[int(r["replicate"])] = (float(r["estimate"])-float(r["target"]), float(r["target"]))


def bias(cell, method):
    values = list(cell[method].values())
    return 100*abs(statistics.mean(v[0] for v in values)/values[0][1])


def ratio(cell):
    shared = cell[CAND].keys() & cell[OPTION].keys()
    return math.sqrt(sum(cell[CAND][i][0]**2 for i in shared)/sum(cell[OPTION][i][0]**2 for i in shared))


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


fresh_predictor = [v for k, v in store.items() if k[0] == "fresh" and k[1] == "predictor"]
shared = [v for v in fresh_predictor if "ind_mice" in v]
say("POST HOC. Written after the results were read. It enters no rule.")
say()
say("PREDICTION 4: new simulation, predictor class, mean absolute relative bias %")
say(f"  all {len(fresh_predictor)} cells: candidate {statistics.mean(bias(v, CAND) for v in fresh_predictor):.2f}, option {statistics.mean(bias(v, OPTION) for v in fresh_predictor):.2f}, ind_mice refused on {len(fresh_predictor)-len(shared)}")
say(f"  the {len(shared)} cells with ind_mice: candidate {statistics.mean(bias(v, CAND) for v in shared):.2f}, ind_mice {statistics.mean(bias(v, 'ind_mice') for v in shared):.2f}, option {statistics.mean(bias(v, OPTION) for v in shared):.2f}")
say()
say("RMSE RATIO, candidate over option, geometric mean across cells, by audit size")
for source in ("snap", "fresh"):
    for role in ("outcome", "predictor"):
        chosen = {k: v for k, v in store.items() if k[0] == source and k[1] == role}
        parts = "; ".join(f"{100*float(f):.0f}% {geometric(ratio(v) for k, v in chosen.items() if k[3] == f):.3f}" for f in sorted({k[3] for k in chosen}))
        say(f"  {source} {role}: {parts}")

guard = {}
with gzip.open(ROOT/"validation/variance-mi-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["imputer"] == OPTION and r["status"] != "ok":
            guard[(r["analysis"], r["fraction"])] = guard.get((r["analysis"], r["fraction"]), 0)+1
with gzip.open(STUDY/"guard.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    missing = {}
    for r in csv.DictReader(stream):
        if r["status"] != "ok":
            missing[(r["analysis"], r["fraction"])] = missing.get((r["analysis"], r["fraction"]), 0)+1
say()
say("GUARD REQUESTS WITH NO INTERVAL, candidate then the option's saved count")
for k in sorted(set(missing) | set(guard)):
    say(f"  {k[0]} {100*float(k[1]):.0f}%: {missing.get(k, 0)}, {guard.get(k, 0)}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
