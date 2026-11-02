#!/usr/bin/env python3
"""Scores the second-system test against the rule in validation/SECOND_SYSTEM_PLAN.md. Plain Python, from saved results.
Usage: python3 -I tools/second-system-score.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/second-system-study"
ITEMS = ["opd", "deliveries", "anc", "postnatal", "immunization", "ipt", "fpmethods", "counselling", "exposedbabies"]
AUDITS = ["clinic_25", "clinic_50", "district_33", "district_50"]
METHODS = ["oracle", "recorded", "audited_only", "linear_mi", "random_intercept", "mice_2l_lmer"]
OUT = []


def say(text=""):
    OUT.append(text)


store = {}
counts = {}
messages = {}
seconds = {}
with gzip.open(STUDY/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        key = (r["run"], r["item"], r["audit"], r["method"])
        good = r["status"] == "ok"
        store.setdefault(key, {})[int(r["replicate"])] = dict(error=float(r["estimate"])-float(r["target"]), cover=float(r["cover"]), width=float(r["width"]), kept=float(r["kept"]), clusters=int(r["audited_clusters"])) if good else None
        counts[(r["run"], r["method"], r["status"])] = counts.get((r["run"], r["method"], r["status"]), 0)+1
        if not good:
            messages[(r["run"], r["method"], r["message"])] = messages.get((r["run"], r["method"], r["message"]), 0)+1
        seconds.setdefault(r["method"], []).append(float(r["seconds"]))


def good(run, item, audit, method):
    return [v for v in store[(run, item, audit, method)].values() if v is not None]


def coverage(run, item, audit, method):
    values = good(run, item, audit, method)
    return 100*statistics.mean(v["cover"] for v in values) if values else float("nan")


def distance(run, audit, method):
    return statistics.mean(abs(coverage(run, i, audit, method)-95) for i in ITEMS)


def rmse(values):
    return math.sqrt(sum(v*v for v in values)/len(values))


def ratio(run, item, audit, top, bottom):
    a, b = store[(run, item, audit, top)], store[(run, item, audit, bottom)]
    both = [r for r in a if a[r] is not None and b.get(r) is not None]
    return rmse([a[r]["error"] for r in both])/rmse([b[r]["error"] for r in both]) if both else float("nan")


def geometric(run, audit, top, bottom):
    return math.exp(statistics.mean(math.log(ratio(run, i, audit, top, bottom)) for i in ITEMS))


def available(run, audit, method):
    total = sum(len(store[(run, i, audit, method)]) for i in ITEMS)
    return 100*sum(len(good(run, i, audit, method)) for i in ITEMS)/total


def at_least(run, audit, method, limit=92):
    return sum(coverage(run, i, audit, method) >= limit for i in ITEMS)


say("Scored against validation/SECOND_SYSTEM_PLAN.md")
say()
say("REQUESTS AND STATUSES")
for run in ("primary", "sensitivity"):
    for m in METHODS:
        parts = ", ".join(f"{s} {counts[(run, m, s)]}" for s in ("ok", "failed") if (run, m, s) in counts)
        if parts:
            say(f"  {run} {m}: {parts}")
for (run, m, text), n in sorted(messages.items()):
    say(f"    {run} {m}: {n} x {text}")
say()
say("DECISION RULE: D is mean distance from 95 over nine items; R is the geometric mean RMSE of the option over linear MI")
verdicts = {}
for audit in ("clinic_25", "district_50"):
    d_option, d_linear = distance("primary", audit, "random_intercept"), distance("primary", audit, "linear_mi")
    r = geometric("primary", audit, "random_intercept", "linear_mi")
    av = available("primary", audit, "random_intercept")
    supports = d_option <= d_linear-.5 and r <= 1.02 and av >= 95
    against = not d_option < d_linear or r > 1.05 or av < 95
    verdicts[audit] = "supports" if supports else "goes against" if against else "neither"
    say(f"  {audit}: D option {d_option:.2f}, D linear MI {d_linear:.2f}, R {r:.3f}, available {av:.2f}%: {verdicts[audit]}")
values = list(verdicts.values())
outcome = "REPLICATES" if values == ["supports", "supports"] else "CONTRADICTS" if values == ["goes against", "goes against"] else "MIXED"
say()
say(f"OUTCOME: {outcome}")
for audit in AUDITS:
    say()
    say(f"AUDIT {audit}: coverage, then RMSE ratio against linear MI; last row is mean distance from 95, items at 92 or more, geometric mean ratio")
    for item in ITEMS:
        say(f"  {item:14s} "+"; ".join(f"{m} {coverage('primary', item, audit, m):.1f}/{ratio('primary', item, audit, m, 'linear_mi'):.3f}" for m in METHODS))
    say("  summary        "+"; ".join(f"{m} {distance('primary', audit, m):.2f}, {at_least('primary', audit, m)}, {geometric('primary', audit, m, 'linear_mi'):.3f}" for m in METHODS))
    say("  available      "+"; ".join(f"{m} {available('primary', audit, m):.2f}%" for m in METHODS))
    kept = {i: statistics.mean(v["kept"] for v in good("primary", i, audit, "linear_mi")) for i in ITEMS}
    say("  recorded kept, linear MI, share of requests: "+", ".join(f"{i} {100*kept[i]:.1f}%" for i in ITEMS if kept[i] > 0))
    say("  mean width against the option: "+"; ".join(f"{m} {statistics.mean(statistics.mean(v['width'] for v in good('primary', i, audit, m))/statistics.mean(v['width'] for v in good('primary', i, audit, 'random_intercept')) for i in ITEMS):.3f}" for m in METHODS))
say()
say("DESCRIPTION: items where recorded has lower RMSE than linear MI, and than every audit-using method")
for audit in AUDITS:
    against_linear = [i for i in ITEMS if ratio("primary", i, audit, "recorded", "linear_mi") < 1]
    against_all = [i for i in ITEMS if all(ratio("primary", i, audit, "recorded", m) < 1 for m in ("audited_only", "linear_mi", "random_intercept", "mice_2l_lmer"))]
    say(f"  {audit}: {len(against_linear)} and {len(against_all)} of 9")
say()
say("SENSITIVITY, both-zero rows kept, 500 cohorts: D for oracle, linear MI, option; R option over linear MI")
for audit in AUDITS:
    say(f"  {audit}: {distance('sensitivity', audit, 'oracle'):.2f}, {distance('sensitivity', audit, 'linear_mi'):.2f}, {distance('sensitivity', audit, 'random_intercept'):.2f}; {geometric('sensitivity', audit, 'random_intercept', 'linear_mi'):.3f}")
say()
say("ERROR INTRACLASS CORRELATION BY DISTRICT, description only")
with (STUDY/"error-icc.csv").open(newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        say(f"  {r['item']:14s} clinics {r['clinics']}, disagree {100*float(r['share_disagree']):.1f}%, mean error {float(r['mean_error']):.2f}, error SD {float(r['sd_error']):.2f}, ICC {float(r['error_icc']):.3f}, mean verified {float(r['outcome_mean']):.1f}")
say()
say("SECONDS PER REQUEST, mean")
for m in METHODS:
    say(f"  {m}: {statistics.mean(seconds[m]):.3f}")
say()
say("PREDICTIONS")
predictions = [
    ("The outcome is MIXED", outcome == "MIXED"),
    ("The oracle reaches 92% on at most six of nine items at district_50", at_least("primary", "district_50", "oracle") <= 6),
    ("The option gives an interval on at least 95% of requests in every audit", all(available("primary", a, "random_intercept") >= 95 for a in AUDITS)),
    ("recorded has a lower RMSE than linear_mi on at least five of nine items at clinic_25", sum(ratio("primary", i, "clinic_25", "recorded", "linear_mi") < 1 for i in ITEMS) >= 5),
    ("Linear MI reaches 92% on at least seven of nine items at clinic_50", at_least("primary", "clinic_50", "linear_mi") >= 7),
]
for i, (text, held) in enumerate(predictions, 1):
    say(f"  {i}. {text}: {'HELD' if held else 'FAILED'}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
