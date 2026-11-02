#!/usr/bin/env python3
"""Scores the test of within-cluster and between-cluster slopes, as validation/WITHIN_BETWEEN_PLAN.md requires. As a control it reproduces the option's saved values.
Usage: python3 -I tools/within-between-score.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/within-between-study"
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
CAND = "within_between"
OPTION = "random_intercept"
OUT = []


def say(text=""):
    OUT.append(text)


def reader(path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        yield from csv.DictReader(stream)


def entry(ok, estimate, target, cover, width):
    return (float(estimate)-float(target), float(cover), float(width), float(target)) if ok else None


store = {}
roles = {}
fallback = {}
seconds = {}
statuses = {}
messages = {}
for r in reader(STUDY/"deciding.csv.gz"):
    key = (r["source"], r["analysis"], float(r["fraction"]))
    roles[key] = r["role"]
    store.setdefault(r["method"], {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])
    statuses[(r["source"], r["method"], r["status"])] = statuses.get((r["source"], r["method"], r["status"]), 0)+1
    seconds.setdefault((r["source"], r["method"]), []).append(float(r["seconds"]))
    if r["status"] != "ok":
        messages[(r["source"], r["method"], r["message"])] = messages.get((r["source"], r["method"], r["message"]), 0)+1
    if r["method"] == CAND and r["status"] == "ok":
        fb = fallback.setdefault(r["source"], [0, 0])
        fb[0] += float(r["fallback"]) == 1
        fb[1] += 1
levels = {}
for r in reader(STUDY/"guard.csv.gz"):
    key = ("guard", r["analysis"], float(r["fraction"]))
    levels[key] = r["level"]
    store.setdefault(CAND, {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])
    statuses[("guard", CAND, r["status"])] = statuses.get(("guard", CAND, r["status"]), 0)+1
    if r["status"] != "ok":
        messages[("guard", CAND, r["message"])] = messages.get(("guard", CAND, r["message"]), 0)+1
    if r["status"] == "ok":
        fb = fallback.setdefault("guard", [0, 0])
        fb[0] += float(r["fallback"]) == 1
        fb[1] += 1
for r in reader(ROOT/"validation/variance-mi-study/replicates.csv.gz"):
    key = ("guard", r["analysis"], float(r["fraction"]))
    store.setdefault(r["imputer"], {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover_cv3"], r["width_cv3"])
for r in reader(ROOT/"validation/independent-setups-study/replicates.csv.gz"):
    if r["method"] == "ind_mice":
        key = ("guard", r["analysis"], float(r["fraction"]))
        store.setdefault("ind_mice", {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])
for r in reader(STUDY/"zambia.csv.gz"):
    key = ("zambia", r["item"], r["audit"])
    store.setdefault(CAND, {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])
    statuses[("zambia", CAND, r["status"])] = statuses.get(("zambia", CAND, r["status"]), 0)+1
    if r["status"] != "ok":
        messages[("zambia", CAND, r["message"])] = messages.get(("zambia", CAND, r["message"]), 0)+1
    if r["status"] == "ok":
        fb = fallback.setdefault("zambia", [0, 0])
        fb[0] += float(r["fallback"]) == 1
        fb[1] += 1
for r in reader(ROOT/"validation/second-system-study/replicates.csv.gz"):
    if r["run"] == "primary" and r["method"] == OPTION:
        key = ("zambia", r["item"], r["audit"])
        store.setdefault(OPTION, {}).setdefault(key, {})[int(r["replicate"])] = entry(r["status"] == "ok", r["estimate"], r["target"], r["cover"], r["width"])


def rows(method, key):
    return [v for v in store[method][key].values() if v is not None]


def coverage(method, key):
    values = rows(method, key)
    return 100*statistics.mean(v[1] for v in values) if values else float("nan")


def rmse(method, key):
    values = rows(method, key)
    return math.sqrt(sum(v[0]**2 for v in values)/len(values)) if values else float("nan")


def relative_bias(method, key):
    values = rows(method, key)
    return abs(statistics.mean(v[0] for v in values)/values[0][3]) if values else float("nan")


def mean_width(method, key):
    values = rows(method, key)
    return statistics.mean(v[2] for v in values) if values else float("nan")


def paired_ratio(method, key, other=OPTION):
    a, b = store[method][key], store[other][key]
    both = [r for r in a if a[r] is not None and b.get(r) is not None]
    return math.sqrt(sum(a[r][0]**2 for r in both)/sum(b[r][0]**2 for r in both)) if both else float("nan")


def geometric(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def available(method, keys):
    return 100*sum(1 for k in keys for v in store[method][k].values() if v is not None)/sum(len(store[method][k]) for k in keys)


def ratio(method, keys):
    return geometric([paired_ratio(method, k) for k in keys])


def distance(method, keys):
    return statistics.mean(abs(coverage(method, k)-95) for k in keys)


def width_ratio(method, keys):
    return geometric([mean_width(method, k)/mean_width(OPTION, k) for k in keys])


def bias(method, keys):
    return statistics.mean(relative_bias(method, k) for k in keys)


def at_least(method, keys, limit=92):
    return sum(coverage(method, k) >= limit for k in keys)


def at_most(method, keys, limit=98):
    return sum(coverage(method, k) <= limit for k in keys)


def in_limits(method, keys):
    return sum(coverage(method, k) <= 98 and (coverage(method, k) >= 92 or coverage(method, k) >= coverage("oracle", k)) for k in keys)


def judge_class(keys, predictor):
    measures = {"ratio": ratio(CAND, keys), "distance": distance(CAND, keys), "option distance": distance(OPTION, keys), "width": width_ratio(CAND, keys), "bias": bias(CAND, keys), "option bias": bias(OPTION, keys), "available": available(CAND, keys)}
    no_worse = measures["ratio"] <= 1.02 and measures["distance"] <= measures["option distance"]+.25 and measures["available"] >= 99
    better = no_worse and measures["width"] <= 1.25 and (measures["distance"] <= measures["option distance"]-.25 or measures["ratio"] <= .98)
    if predictor:
        better = better and measures["bias"] <= measures["option bias"]
    return measures, no_worse, better


def show_class(title, result):
    m, no_worse, better = result
    say(f"    {title}: ratio {m['ratio']:.3f}; distance {m['distance']:.2f} against {m['option distance']:.2f}; width {m['width']:.3f}; bias {100*m['bias']:.2f}% against {100*m['option bias']:.2f}%; available {m['available']:.2f}%; no worse {'YES' if no_worse else 'no'}; better {'YES' if better else 'no'}")


real25 = [("guard", a, .25) for a in ITEMS]
real50 = [("guard", a, .5) for a in ITEMS]
old_outcome = sorted(k for k, level in levels.items() if level == "simulation" and k[1] != "control_one_row")
old_predictor = sorted(k for k, level in levels.items() if level == "imputed_covariate")
AUDITS = ["clinic_25", "clinic_50", "district_33", "district_50"]
zambia_items = sorted({k[1] for k in store[CAND] if k[0] == "zambia"})

full = len(rows(OPTION, real25[0])) >= 2000 and len(rows(CAND, real25[0])) >= 2000
if full:
    saved = (ROOT/"validation/independent-setups-decisions.txt").read_text(encoding="utf-8")
    wanted = f"OPTION: real 25% at least 92 on {at_least(OPTION, real25)} of 11, distance {distance(OPTION, real25):.2f}; real 50% at least 92 on {at_least(OPTION, real50)} of 11, distance {distance(OPTION, real50):.2f}; simulated cells in limits {in_limits(OPTION, old_outcome)} of {len(old_outcome)}"
    assert wanted in saved, "This scorer does not reproduce the saved values for the option"
    control = "Control: this scorer reproduces the saved real-item and simulation values for the option."
else:
    control = "Control skipped on a pilot."
say("Scored against validation/WITHIN_BETWEEN_PLAN.md")
say(control)
say()
say("REQUESTS AND STATUSES")
for (source, method, status), n in sorted(statuses.items()):
    say(f"  {source} {method} {status}: {n}")
for (source, method, text), n in sorted(messages.items()):
    say(f"    {source} {method}: {n} x {text}")
for source, (fell, total) in sorted(fallback.items()):
    say(f"  {source}: requests that fell back to the option's model: {fell} of {total}")
say()
say("DECISION RULE")
passes = {}
for source in ("snap", "fresh"):
    outcome = sorted(k for k, role in roles.items() if k[0] == source and role == "outcome")
    predictor = sorted(k for k, role in roles.items() if k[0] == source and role == "predictor")
    a, b = judge_class(outcome, False), judge_class(predictor, True)
    say(f"  {source}: {len(outcome)} outcome cells, {len(predictor)} predictor cells")
    show_class("outcome", a)
    show_class("predictor", b)
    passes[source] = a[1] and b[1] and (a[2] or b[2])
    say(f"    Source passes: {'YES' if passes[source] else 'no'}")
guard = {}
for name, keys in (("G1 real 25%", real25), ("G2 real 50%", real50)):
    r, d, up, cap = ratio(CAND, keys), distance(CAND, keys), at_least(CAND, keys), at_most(CAND, keys)
    guard[name] = r <= 1.02 and up >= at_least(OPTION, keys) and cap >= 10 and d <= distance(OPTION, keys)+.25
    say(f"  {name}: ratio {r:.3f}; at least 92 on {up} (option {at_least(OPTION, keys)}); at most 98 on {cap}; distance {d:.2f} (option {distance(OPTION, keys):.2f}); {'pass' if guard[name] else 'FAIL'}")
r, lim = ratio(CAND, old_outcome), in_limits(CAND, old_outcome)
guard["G3"] = r <= 1.02 and lim >= in_limits(OPTION, old_outcome)
say(f"  G3 existing outcome cells: ratio {r:.3f}; in limits {lim} of {len(old_outcome)} (option {in_limits(OPTION, old_outcome)}); {'pass' if guard['G3'] else 'FAIL'}")
r, b = ratio(CAND, old_predictor), bias(CAND, old_predictor)
guard["G4"] = r <= 1.02 and b <= bias(OPTION, old_predictor)
say(f"  G4 existing predictor cells: ratio {r:.3f}; bias {100*b:.2f}% (option {100*bias(OPTION, old_predictor):.2f}%); {'pass' if guard['G4'] else 'FAIL'}")
for audit in AUDITS:
    keys = [("zambia", item, audit) for item in zambia_items]
    r, d, av = ratio(CAND, keys), distance(CAND, keys), available(CAND, keys)
    guard[f"G5 {audit}"] = r <= 1.02 and d <= distance(OPTION, keys)+.5 and av >= available(OPTION, keys)-1
    say(f"  G5 Zambia {audit}: ratio {r:.3f}; distance {d:.2f} (option {distance(OPTION, keys):.2f}); available {av:.2f}% (option {available(OPTION, keys):.2f}%); {'pass' if guard[f'G5 {audit}'] else 'FAIL'}")
holds = all(guard.values())
say(f"  Guard holds: {'YES' if holds else 'no'}")
say()
adopt = passes["snap"] and passes["fresh"] and holds
say(f"OUTCOME: {'ADOPT' if adopt else 'NOT ADOPTED'}. SNAP passes: {'yes' if passes['snap'] else 'no'}. New simulation passes: {'yes' if passes['fresh'] else 'no'}. Guard holds: {'yes' if holds else 'no'}.")
say()
methods = ["oracle", "recorded", "linear_mi", OPTION, CAND, "ind_mice"]
for source, title in (("snap", "SNAP"), ("fresh", "NEW SIMULATION")):
    say(f"{title}, each cell: coverage, RMSE ratio against the option, bias %")
    for k in sorted(k for k in roles if k[0] == source):
        say(f"  {k[1]:32s} {roles[k][:4]} {int(k[2]*100)}% "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}/{100*relative_bias(m, k):.1f}" for m in methods))
    say()
say("DESCRIPTION: each class against ind_mice, no rule")
for source in ("snap", "fresh"):
    for role in ("outcome", "predictor"):
        keys = sorted(k for k, value in roles.items() if k[0] == source and value == role)
        say(f"  {source} {role}: "+"; ".join(f"{m} distance {distance(m, keys):.2f}, RMSE ratio {ratio(m, keys):.3f}, width {width_ratio(m, keys):.3f}, bias {100*bias(m, keys):.2f}%, at 92 or more {at_least(m, keys)} of {len(keys)}" for m in ("linear_mi", OPTION, CAND, "ind_mice")))
say()
say("GUARD DETAIL: coverage, then RMSE ratio against the option")
for keys, title in ((real25, "real 25%"), (real50, "real 50%")):
    say(f"  {title}")
    for k in keys:
        say(f"    {k[1]:12s} "+"; ".join(f"{m} {coverage(m, k):.1f}/{paired_ratio(m, k):.3f}" for m in ("oracle", OPTION, CAND, "ind_mice")))
say("  existing predictor cells: bias %, coverage")
for k in old_predictor:
    say(f"    {k[1]:20s} {int(k[2]*100)}% "+"; ".join(f"{m} {100*relative_bias(m, k):.1f}, {coverage(m, k):.1f}" for m in ("linear_mi", OPTION, CAND, "ind_mice")))
say()
say("SECONDS PER REQUEST, mean")
for (source, method), values in sorted(seconds.items()):
    say(f"  {source} {method}: {statistics.mean(values):.3f}")
say()
say("PREDICTIONS")
snap_outcome = sorted(k for k, role in roles.items() if k[0] == "snap" and role == "outcome")
fresh_predictor = sorted(k for k, role in roles.items() if k[0] == "fresh" and role == "predictor")
predictions = [
    ("The outcome is NOT ADOPTED", not adopt),
    ("On SNAP, the outcome class has a ratio between 0.98 and 1.02", .98 <= ratio(CAND, snap_outcome) <= 1.02),
    ("On the new simulation, the candidate's predictor class has bias below the option's", bias(CAND, fresh_predictor) < bias(OPTION, fresh_predictor)),
    ("On the new simulation, ind_mice has a smaller predictor-class bias than the candidate", bias("ind_mice", fresh_predictor) < bias(CAND, fresh_predictor)),
]
for i, (text, held) in enumerate(predictions, 1):
    say(f"  {i}. {text}: {'HELD' if held else 'FAILED'}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
