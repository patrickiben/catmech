#!/usr/bin/env python3
"""Scores the comparison with outside methods against the rule in validation/OUTSIDE_METHODS_PLAN.md. Plain Python, from saved results.
Usage: python3 -I tools/outside-methods-score.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/outside-methods-study"
EARLIER = ROOT/"validation/variance-mi-study/replicates.csv.gz"
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
JUDGED = ["audited_only", "difference", "raking", "mice_2l_lmer", "mice_2l_pan"]
OPTION = "random_intercept"
ALL_AUDIT = JUDGED+[OPTION, "linear_mi"]
OUT = []


def say(text=""):
    OUT.append(text)


# method -> (analysis, fraction) -> replicate -> (estimate-target or None, cover, width)
store = {}
cell_level = {}
seconds = {}
logged = {}
fitted = {}
counts = {}
messages = {}


def add(method, level, analysis, fraction, replicate, status, estimate, target, cover, width):
    key = (analysis, fraction)
    cell_level[key] = level
    good = status == "ok"
    store.setdefault(method, {}).setdefault(key, {})[replicate] = (float(estimate)-float(target), float(cover), float(width)) if good else None


with gzip.open(STUDY/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        fraction = float(r["fraction"])
        replicate = int(r["replicate"])
        add(r["method"], r["level"], r["analysis"], fraction, replicate, r["status"], r["estimate"], r["target"], r["cover"], r["width"])
        counts[(r["method"], r["status"])] = counts.get((r["method"], r["status"]), 0)+1
        if r["status"] != "ok":
            messages[(r["method"], r["message"])] = messages.get((r["method"], r["message"]), 0)+1
        seconds.setdefault(r["method"], []).append(float(r["seconds"]))
        logged[r["method"]] = logged.get(r["method"], 0)+(int(r["logged"]) > 0)
        if r["method"] == "difference" and r["fitted_estimate"] not in ("", "NA"):
            fitted.setdefault((r["analysis"], fraction), {})[replicate] = float(r["fitted_estimate"])-float(r["target"])
with gzip.open(EARLIER, "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        add(r["imputer"], r["level"], r["analysis"], float(r["fraction"]), int(r["replicate"]), r["status"], r["estimate"], r["target"], r["cover_cv3"], r["width_cv3"])


def rows(method, key):
    return [v for v in store[method][key].values() if v is not None]


def rmse(method, key):
    values = rows(method, key)
    return math.sqrt(sum(v[0]**2 for v in values)/len(values)) if values else float("nan")


def coverage(method, key):
    values = rows(method, key)
    return 100*statistics.mean(v[1] for v in values) if values else float("nan")


def mean_width(method, key):
    values = rows(method, key)
    return statistics.mean(v[2] for v in values) if values else float("nan")


def bias(method, key):
    values = rows(method, key)
    return statistics.mean(v[0] for v in values) if values else float("nan")


def paired_ratio(method, key):
    a, b = store[method][key], store[OPTION][key]
    both = [r for r in a if a[r] is not None and b.get(r) is not None]
    if not both:
        return float("nan")
    return math.sqrt(sum(a[r][0]**2 for r in both)/sum(b[r][0]**2 for r in both))


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def available(method, keys):
    requested = sum(len(store[method][k]) for k in keys)
    return 100*sum(1 for k in keys for v in store[method][k].values() if v is not None)/requested


real25 = [(a, .25) for a in ITEMS]
real50 = [(a, .5) for a in ITEMS]
outcome_cells = sorted(k for k, level in cell_level.items() if level == "simulation" and k[0] != "control_one_row")
control_cells = sorted(k for k, level in cell_level.items() if level == "simulation" and k[0] == "control_one_row")
predictor_cells = sorted(k for k, level in cell_level.items() if level == "imputed_covariate")


def distance(method, keys):
    return statistics.mean(abs(coverage(method, k)-95) for k in keys)


def at_least(method, keys, limit=92):
    return sum(coverage(method, k) >= limit for k in keys)


def at_most(method, keys, limit=98):
    return sum(coverage(method, k) <= limit for k in keys)


def in_limits(method, keys):
    return sum(coverage(method, k) <= 98 and (coverage(method, k) >= 92 or coverage(method, k) >= coverage("oracle", k)) for k in keys)


say("Scored against validation/OUTSIDE_METHODS_PLAN.md")
say()
say("REQUESTS AND STATUSES")
for method in ["recorded"]+JUDGED:
    parts = ", ".join(f"{s} {counts.get((method, s), 0)}" for s in ("ok", "refused", "failed") if counts.get((method, s), 0))
    say(f"  {method}: {parts}; requests with a logged event from mice {logged.get(method, 0)}")
for (method, message), n in sorted(messages.items()):
    say(f"    {method}: {n} x {message}")
say()
say(f"OPTION, the comparison point: real 25% at least 92 on {at_least(OPTION, real25)} of 11, distance {distance(OPTION, real25):.2f}; real 50% at least 92 on {at_least(OPTION, real50)} of 11, distance {distance(OPTION, real50):.2f}; simulated cells in limits {in_limits(OPTION, outcome_cells)} of {len(outcome_cells)}")
say()
say("DECISION RULE")
passing = []
results = {}
for c in JUDGED:
    n1 = geometric([paired_ratio(c, k) for k in real25])
    n2 = geometric([paired_ratio(c, k) for k in real50])
    n3 = (at_least(c, real25), at_most(c, real25), distance(c, real25))
    n4 = (at_least(c, real50), at_most(c, real50), distance(c, real50))
    n5 = (geometric([paired_ratio(c, k) for k in outcome_cells]), in_limits(c, outcome_cells))
    n6 = (available(c, real25), available(c, real50), available(c, outcome_cells))
    passes = {
        "N1": n1 <= 1.05,
        "N2": n2 <= 1.05,
        "N3": n3[0] >= at_least(OPTION, real25) and n3[1] >= 10 and n3[2] <= distance(OPTION, real25)+.25,
        "N4": n4[0] >= at_least(OPTION, real50)-1 and n4[1] >= 10 and n4[2] <= distance(OPTION, real50)+.5,
        "N5": n5[0] <= 1.05 and n5[1] >= in_limits(OPTION, outcome_cells)-2,
        "N6": min(n6) >= 99,
    }
    results[c] = (n1, n3[2])
    say(f"  {c}:")
    say(f"    N1 RMSE ratio, real 25%: {n1:.3f} ({'pass' if passes['N1'] else 'FAIL'})")
    say(f"    N2 RMSE ratio, real 50%: {n2:.3f} ({'pass' if passes['N2'] else 'FAIL'})")
    say(f"    N3 real 25%: at least 92 on {n3[0]}, at most 98 on {n3[1]}, distance {n3[2]:.2f} ({'pass' if passes['N3'] else 'FAIL'})")
    say(f"    N4 real 50%: at least 92 on {n4[0]}, at most 98 on {n4[1]}, distance {n4[2]:.2f} ({'pass' if passes['N4'] else 'FAIL'})")
    say(f"    N5 simulations: RMSE ratio {n5[0]:.3f}, cells in limits {n5[1]} of {len(outcome_cells)} ({'pass' if passes['N5'] else 'FAIL'})")
    say(f"    N6 available: {n6[0]:.2f}%, {n6[1]:.2f}%, {n6[2]:.2f}% ({'pass' if passes['N6'] else 'FAIL'})")
    say(f"    No worse: {'YES' if all(passes.values()) else 'no'}")
    if all(passes.values()):
        passing.append(c)
say()
if passing:
    named = min(passing, key=lambda c: (results[c][0], results[c][1]))
    say(f"OUTCOME: RETIRE. No worse: {', '.join(passing)}. Named in the deprecation warning: {named}.")
else:
    say("OUTCOME: KEEP. No comparator was no worse than the option.")
say()
say("REAL ITEMS AT 25%: coverage, then RMSE ratio against the option")
methods = ["oracle", "linear_mi", OPTION, "recorded"]+JUDGED
for k in real25:
    say(f"  {k[0]:12s} "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in methods))
say()
say("REAL ITEMS AT 50%: coverage, then RMSE ratio against the option")
for k in real50:
    say(f"  {k[0]:12s} "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in methods))
say()
say("MEAN WIDTH RATIO AGAINST THE OPTION, real items, range over items")
for f, keys in ((25, real25), (50, real50)):
    for m in ["recorded"]+JUDGED:
        ratios = [mean_width(m, k)/mean_width(OPTION, k) for k in keys]
        say(f"  {f}% {m}: {100*(min(ratios)-1):+.0f}% to {100*(max(ratios)-1):+.0f}%")
say()
say("DESCRIPTION: real items where recorded has lower RMSE than every audit-using method")
for f, keys in ((25, real25), (50, real50)):
    wins = [k[0] for k in keys if all(rmse("recorded", k) < rmse(m, k) for m in ALL_AUDIT)]
    say(f"  {f}%: {len(wins)} of 11: {', '.join(wins)}")
say()
say("SIMULATED CELLS, AUDITED VALUE AS OUTCOME: coverage, then RMSE ratio against the option")
for k in outcome_cells+control_cells:
    say(f"  {k[0]:24s} {int(k[1]*100)}% "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in methods))
say()
say("AUDITED VALUE AS PREDICTOR, no rule: bias %, coverage, RMSE, RMSE ratio against linear MI")
for k in predictor_cells:
    a = rmse("linear_mi", k)
    say(f"  {k[0]:20s} {int(k[1]*100)}% "+"; ".join(f"{m} {100*bias(m, k)/.5:+.1f}, {coverage(m, k):.1f}, {rmse(m, k):.4f}, {rmse(m, k)/a:.3f}" for m in methods))
say()
say("FITTED_VALUE, exploratory: RMSE ratio against linear MI and the option, real items")
for keys, label in ((real25, "25%"), (real50, "50%")):
    a, b = [], []
    for k in keys:
        e = fitted[k]
        for name, target in (("linear_mi", a), (OPTION, b)):
            other = store[name][k]
            both = [r for r in e if other.get(r) is not None]
            target.append(math.sqrt(sum(e[r]**2 for r in both)/sum(other[r][0]**2 for r in both)))
    say(f"  {label}: against linear MI, geometric mean {geometric(a):.3f}; against the option {geometric(b):.3f}")
say()
say("SECONDS PER REQUEST, on one core of a busy machine, mean")
for m in ["recorded"]+JUDGED:
    say(f"  {m}: {statistics.mean(seconds[m]):.3f}")
control = STUDY/"control.csv"
if control.exists():
    with control.open(newline="", encoding="utf-8") as stream:
        values = [float(r["seconds"]) for r in csv.DictReader(stream)]
    say(f"  {OPTION}, from the control runs: {statistics.mean(values):.3f}")
say()
say("PREDICTIONS")
ratio = {c: geometric([paired_ratio(c, k) for k in real25]) for c in JUDGED}
wins = sum(all(rmse("recorded", k) < rmse(m, k) for m in ALL_AUDIT) for k in real25)
predictions = [
    ("The outcome is RETIRE, through mice_2l_lmer", "mice_2l_lmer" in passing),
    ("raking fails N1", ratio["raking"] > 1.05),
    ("difference reaches 92% on at least the option's count at 25%, and fails N1", at_least("difference", real25) >= at_least(OPTION, real25) and ratio["difference"] > 1.05),
    ("audited_only has an N1 ratio above 1.20", ratio["audited_only"] > 1.20),
    ("mice_2l_pan's distance at 25% is within 0.50 points of mice_2l_lmer's", abs(distance("mice_2l_pan", real25)-distance("mice_2l_lmer", real25)) <= .5),
    ("recorded has lower RMSE than every audit-using method on at least three items at 25%", wins >= 3),
]
for i, (text, held) in enumerate(predictions, 1):
    say(f"  {i}. {text}: {'HELD' if held else 'FAILED'}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
