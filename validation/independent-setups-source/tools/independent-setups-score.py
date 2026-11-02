#!/usr/bin/env python3
"""Scores comparators set up by fresh agents with the rule of the comparison with outside methods, as validation/INDEPENDENT_SETUPS_PLAN.md requires.
The conditions are the same as tools/outside-methods-score.py, which is frozen. As a control this file also reproduces that study's saved values.
Usage: python3 -I tools/independent-setups-score.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/independent-setups-study"
EARLIER = ROOT/"validation/variance-mi-study/replicates.csv.gz"
FIRST_RUN = ROOT/"validation/outside-methods-study/replicates.csv.gz"
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
JUDGED = ["ind_mice", "ind_survey", "ind_design"]
OPTION = "random_intercept"
OUT = []


def say(text=""):
    OUT.append(text)


store = {}
cell_level = {}
counts = {}
messages = {}
seconds = {}


def add(method, level, analysis, fraction, replicate, ok, estimate, target, cover, width):
    key = (analysis, fraction)
    cell_level[key] = level
    store.setdefault(method, {}).setdefault(key, {})[replicate] = (float(estimate)-float(target), float(cover), float(width)) if ok else None


with gzip.open(STUDY/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        ok = r["status"] == "ok"
        add(r["method"], r["level"], r["analysis"], float(r["fraction"]), int(r["replicate"]), ok, r["estimate"], r["target"], r["cover"], r["width"])
        counts[(r["method"], r["status"])] = counts.get((r["method"], r["status"]), 0)+1
        if not ok:
            messages[(r["method"], r["message"])] = messages.get((r["method"], r["message"]), 0)+1
        seconds.setdefault(r["method"], []).append(float(r["seconds"]))
with gzip.open(EARLIER, "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        add(r["imputer"], r["level"], r["analysis"], float(r["fraction"]), int(r["replicate"]), r["status"] == "ok", r["estimate"], r["target"], r["cover_cv3"], r["width_cv3"])
with gzip.open(FIRST_RUN, "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["method"] in ("mice_2l_lmer", "audited_only", "recorded"):
            add(r["method"], r["level"], r["analysis"], float(r["fraction"]), int(r["replicate"]), r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])


def rows(method, key):
    return [v for v in store[method][key].values() if v is not None]


def coverage(method, key):
    values = rows(method, key)
    return 100*statistics.mean(v[1] for v in values) if values else float("nan")


def rmse(method, key):
    values = rows(method, key)
    return math.sqrt(sum(v[0]**2 for v in values)/len(values)) if values else float("nan")


def bias(method, key):
    values = rows(method, key)
    return statistics.mean(v[0] for v in values) if values else float("nan")


def mean_width(method, key):
    values = rows(method, key)
    return statistics.mean(v[2] for v in values) if values else float("nan")


def paired_ratio(method, key):
    a, b = store[method][key], store[OPTION][key]
    both = [r for r in a if a[r] is not None and b.get(r) is not None]
    return math.sqrt(sum(a[r][0]**2 for r in both)/sum(b[r][0]**2 for r in both)) if both else float("nan")


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def available(method, keys):
    return 100*sum(1 for k in keys for v in store[method][k].values() if v is not None)/sum(len(store[method][k]) for k in keys)


real25 = [(a, .25) for a in ITEMS]
real50 = [(a, .5) for a in ITEMS]
outcome_cells = sorted(k for k, level in cell_level.items() if level == "simulation" and k[0] != "control_one_row")
predictor_cells = sorted(k for k, level in cell_level.items() if level == "imputed_covariate")


def distance(method, keys):
    return statistics.mean(abs(coverage(method, k)-95) for k in keys)


def at_least(method, keys, limit=92):
    return sum(coverage(method, k) >= limit for k in keys)


def at_most(method, keys, limit=98):
    return sum(coverage(method, k) <= limit for k in keys)


def in_limits(method, keys):
    return sum(coverage(method, k) <= 98 and (coverage(method, k) >= 92 or coverage(method, k) >= coverage("oracle", k)) for k in keys)


def judge(c):
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
    return n1, n2, n3, n4, n5, n6, passes


saved = (ROOT/"validation/outside-methods-decisions.txt").read_text(encoding="utf-8")
control = judge("mice_2l_lmer")
wanted = [f"N1 RMSE ratio, real 25%: {control[0]:.3f}", f"N3 real 25%: at least 92 on {control[2][0]}, at most 98 on {control[2][1]}, distance {control[2][2]:.2f}", f"N5 simulations: RMSE ratio {control[4][0]:.3f}, cells in limits {control[4][1]} of {len(outcome_cells)}"]
block = saved[saved.index("  mice_2l_lmer:"):saved.index("  mice_2l_pan:")]
if len(rows("ind_mice", real25[0])) >= 2000:
    assert all(w in block for w in wanted), "This scorer does not reproduce the earlier study's values"
    control_text = "Control: this scorer reproduces the earlier study's N1, N3 and N5 for mice_2l_lmer."
else:
    control_text = "Control skipped on a pilot."
say("Scored against validation/INDEPENDENT_SETUPS_PLAN.md, with the rule of validation/OUTSIDE_METHODS_PLAN.md")
say(control_text)
say()
say("REQUESTS AND STATUSES")
for m in JUDGED:
    say(f"  {m}: "+", ".join(f"{s} {counts.get((m, s), 0)}" for s in ("ok", "failed") if counts.get((m, s), 0)))
for (m, text), n in sorted(messages.items()):
    say(f"    {m}: {n} x {text}")
say()
say(f"OPTION: real 25% at least 92 on {at_least(OPTION, real25)} of 11, distance {distance(OPTION, real25):.2f}; real 50% at least 92 on {at_least(OPTION, real50)} of 11, distance {distance(OPTION, real50):.2f}; simulated cells in limits {in_limits(OPTION, outcome_cells)} of {len(outcome_cells)}")
say()
say("DECISION RULE")
passing = []
results = {}
for c in JUDGED:
    n1, n2, n3, n4, n5, n6, passes = judge(c)
    results[c] = (n1, n3[2])
    mark = lambda name: "pass" if passes[name] else "FAIL"
    say(f"  {c}:")
    say(f"    N1 RMSE ratio, real 25%: {n1:.3f} ({mark('N1')})")
    say(f"    N2 RMSE ratio, real 50%: {n2:.3f} ({mark('N2')})")
    say(f"    N3 real 25%: at least 92 on {n3[0]}, at most 98 on {n3[1]}, distance {n3[2]:.2f} ({mark('N3')})")
    say(f"    N4 real 50%: at least 92 on {n4[0]}, at most 98 on {n4[1]}, distance {n4[2]:.2f} ({mark('N4')})")
    say(f"    N5 simulations: RMSE ratio {n5[0]:.3f}, cells in limits {n5[1]} of {len(outcome_cells)} ({mark('N5')})")
    say(f"    N6 available: {n6[0]:.2f}%, {n6[1]:.2f}%, {n6[2]:.2f}% ({mark('N6')})")
    say(f"    No worse: {'YES' if all(passes.values()) else 'no'}")
    if all(passes.values()):
        passing.append(c)
say()
if passing:
    named = min(passing, key=lambda c: results[c])
    say(f"OUTCOME: RETIRE. No worse: {', '.join(passing)}. Named in the deprecation warning: {named}.")
else:
    say("OUTCOME: KEEP. No setup was no worse than the option.")
say()
methods = ["oracle", "linear_mi", OPTION]+JUDGED
for keys, title in ((real25, "REAL ITEMS AT 25%"), (real50, "REAL ITEMS AT 50%")):
    say(f"{title}: coverage, then RMSE ratio against the option")
    for k in keys:
        say(f"  {k[0]:12s} "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in methods))
    say("  width over the option, range: "+"; ".join(f"{m} {100*(min(mean_width(m, k)/mean_width(OPTION, k) for k in keys)-1):+.0f}% to {100*(max(mean_width(m, k)/mean_width(OPTION, k) for k in keys)-1):+.0f}%" for m in JUDGED))
    say()
say("SIMULATED CELLS, AUDITED VALUE AS OUTCOME: coverage, then RMSE ratio against the option")
for k in outcome_cells:
    say(f"  {k[0]:24s} {int(k[1]*100)}% "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in methods))
say()
say("AUDITED VALUE AS PREDICTOR, no rule: bias %, coverage, RMSE, RMSE ratio against linear MI")
for k in predictor_cells:
    a = rmse("linear_mi", k)
    say(f"  {k[0]:20s} {int(k[1]*100)}% "+"; ".join(f"{m} {100*bias(m, k)/.5:+.1f}, {coverage(m, k):.1f}, {rmse(m, k):.4f}, {rmse(m, k)/a:.3f}" for m in methods))
say()
say("SECOND HEALTH SYSTEM, whole-cluster audits, description only: coverage and RMSE")
zambia = {}
with gzip.open(STUDY/"zambia.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        zambia.setdefault((r["audit"], r["item"], r["method"]), []).append(r)
for audit in ("district_33", "district_50"):
    say(f"  {audit}:")
    for m in JUDGED:
        parts = []
        dists = []
        ok_total = n_total = 0
        for item in sorted({k[1] for k in zambia}):
            got = zambia[(audit, item, m)]
            good = [g for g in got if g["status"] == "ok"]
            ok_total += len(good); n_total += len(got)
            if good:
                cov = 100*statistics.mean(float(g["cover"]) for g in good)
                dists.append(abs(cov-95))
                parts.append(f"{item} {cov:.1f}")
        say(f"    {m}: available {100*ok_total/n_total:.1f}%; mean distance {statistics.mean(dists) if dists else float('nan'):.2f}; "+", ".join(parts))
say()
say("SECONDS PER REQUEST, mean")
for m in JUDGED:
    say(f"  {m}: {statistics.mean(seconds[m]):.3f}")
say()
say("PREDICTIONS")
predictions = [
    ("The outcome is KEEP", not passing),
    ("ind_mice has an N1 ratio at most 1.05, and fails N3", judge("ind_mice")[0] <= 1.05 and not judge("ind_mice")[6]["N3"]),
    ("ind_survey and ind_design fail N1", not judge("ind_survey")[6]["N1"] and not judge("ind_design")[6]["N1"]),
    ("ind_survey reaches 92% on more real items at 25% than the earlier raking recipe under approx, which reached none", at_least("ind_survey", real25) > 0),
]
for i, (text, held) in enumerate(predictions, 1):
    say(f"  {i}. {text}: {'HELD' if held else 'FAILED'}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
