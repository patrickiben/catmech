#!/usr/bin/env python3
"""Scores the rare-counts test against the rules in validation/RARE_COUNTS_PLAN.md. Plain Python, from saved results.
Usage: python3 -I tools/rare-counts-score.py [STUDY_DIRECTORY] [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"validation/rare-counts-study"
EARLIER = ROOT/"validation/variance-mi-study/replicates.csv.gz"
OUT = []


def say(text=""):
    OUT.append(text)


def rows(name):
    with gzip.open(STUDY/name, "rt", newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


num = lambda v: float(v) if v not in ("", "NA") else None
part1 = {}
for r in rows("part1.csv.gz"):
    part1.setdefault(r["label"], []).append(r)
DESIGNS, FRAILTY, LEVELS = ("rwanda", "pairs", "none"), ("none", "shape_2", "shape_05"), ("rare", "common")
labels = [f"{d}_{f}_{l}" for l in LEVELS for f in FRAILTY for d in DESIGNS]
assert set(labels) == set(part1)


def cover(label, name):
    values = [num(r[f"cover_{name}"]) for r in part1[label]]
    values = [v for v in values if v is not None]
    return 100*statistics.mean(values) if values else float("nan")


def width(label, name):
    values = [num(r[f"width_{name}"]) for r in part1[label]]
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else float("nan")


say("Scored against validation/RARE_COUNTS_PLAN.md")
say()
say("PART 1, THE RECOUNT: coverage by reference distribution")
NAMES = ("cv3", "cr2", "wcr", "boot", "poisson")
for label in labels:
    values = [statistics.mean(float(r["D"]) for r in part1[label] if r["D"] not in ("", "NA"))]
    say(f"  {label:28s} "+"; ".join(f"{n} {cover(label, n):.1f}" for n in NAMES)+f"; mean D {values[0]:.3f}; flagged {100*statistics.mean(float(r['D']) > .5 for r in part1[label] if r['D'] not in ('', 'NA')):.1f}%")
say()
rare_cells = lambda f: [f"{d}_{f}_rare" for d in DESIGNS]
m0, m2, m5 = (statistics.mean(cover(c, "cv3") for c in rare_cells(f)) for f in FRAILTY)
outcome_p1 = "HEAVY TAILS" if m0 >= 93 and m5 <= 90 else "RARITY" if m0 < 92 else "MIXED"
common = lambda f: [f"{d}_{f}_common" for d in DESIGNS]
say(f"P1: M0 {m0:.2f}, shape 2 {m2:.2f}, M5 {m5:.2f}. Common level under cv3: no frailty {statistics.mean(cover(c, 'cv3') for c in common('none')):.2f}, shape 2 {statistics.mean(cover(c, 'cv3') for c in common('shape_2')):.2f}, shape 0.5 {statistics.mean(cover(c, 'cv3') for c in common('shape_05')):.2f}")
say(f"OUTCOME P1: {outcome_p1}")
say()
say("P2: mean over rare shape 0.5 cells; at most 98 in all 18; mean over the other 15, and cv3's")
heavy = rare_cells("shape_05")
others = [c for c in labels if c not in heavy]
base = statistics.mean(cover(c, "cv3") for c in others)
eligible = []
for n in ("cr2", "wcr", "boot", "poisson"):
    a = statistics.mean(cover(c, n) for c in heavy)
    top = sum(cover(c, n) <= 98 for c in labels)
    b = statistics.mean(cover(c, n) for c in others)
    ok = a >= 92 and top == 18 and b >= base-.5
    if ok:
        eligible.append(n)
    extra = "" if n == "wcr" else f", width over cv3 in those cells {statistics.mean(width(c, n)/width(c, 'cv3') for c in heavy) if n != 'poisson' else float('nan'):.2f}"
    say(f"  {n}: {a:.2f}; {top} of 18; {b:.2f} against {base:.2f}{extra}: {'ELIGIBLE' if ok else 'not eligible'}")
say(f"OUTCOME P2: {', '.join(eligible) if eligible else 'none eligible'}")
say()

saved = {}
saved_target = {}
with gzip.open(EARLIER, "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            saved[(r["analysis"], float(r["fraction"]), int(r["replicate"]), r["imputer"])] = (float(r["estimate"])-float(r["target"]), float(r["cover_cv3"]), float(r["width_cv3"]))
            saved_target[(r["analysis"], float(r["fraction"]), int(r["replicate"]))] = float(r["target"])
new = {}
states = {}
messages = {}
fallback = {}
gap = 0.0
for r in rows("part2.csv.gz"):
    key = (r["analysis"], float(r["fraction"]), r["method"])
    gap = max(gap, abs(float(r["target"])-saved_target.get((r["analysis"], float(r["fraction"]), int(r["replicate"])), float(r["target"]))))
    states.setdefault(key, []).append(r["status"] == "ok")
    if r["status"] == "ok":
        new.setdefault(key, {})[int(r["replicate"])] = (float(r["estimate"])-float(r["target"]), float(r["cover"]), float(r["width"]))
        fallback.setdefault(r["method"], []).append(r["fallback"] == "TRUE")
    else:
        messages[(r["method"], r["message"])] = messages.get((r["method"], r["message"]), 0)+1
cells = sorted({k[:2] for k in states})
assert len(cells) == 12, cells


def get(method, cell):
    if method in ("nb_mi", "pmm_recorded"):
        return new.get(cell+(method,), {})
    out = {}
    for (a, f, r, m), v in saved.items():
        if (a, f) == cell and m == method and r in reps[cell]:
            out[r] = v
    return out


reps = {c: {r for m in ("nb_mi", "pmm_recorded") for r in new.get(c+(m,), {})} for c in cells}
table = {m: {c: get(m, c) for c in cells} for m in ("oracle", "linear_mi", "random_intercept", "nb_mi", "pmm_recorded")}
cov = lambda m, c: 100*statistics.mean(v[1] for v in table[m][c].values())
dist = lambda m: statistics.mean(abs(cov(m, c)-95) for c in cells)


def rmse_ratio(m, c, base_method="linear_mi"):
    a, b = table[m][c], table[base_method][c]
    both = [r for r in a if r in b]
    return math.sqrt(sum(a[r][0]**2 for r in both)/sum(b[r][0]**2 for r in both))


say(f"PART 2, COUNT-AWARE IMPUTATION: coverage, then RMSE ratio against linear MI. Largest gap between these targets and the earlier study's: {gap:.1e}")
for c in cells:
    say(f"  {c[0]:16s} {int(c[1]*100)}% "+"; ".join(f"{m} {cov(m, c):.1f}/{rmse_ratio(m, c):.3f}" for m in table))
say()
say("P3")
eligible3 = []
for m in ("nb_mi", "pmm_recorded"):
    g = math.exp(statistics.mean(math.log(rmse_ratio(m, c)) for c in cells))
    top = sum(cov(m, c) <= 98 for c in cells)
    total = sum(len(states[c+(m,)]) for c in cells)
    av = 100*sum(sum(states[c+(m,)]) for c in cells)/total
    ok = dist(m) <= dist("linear_mi")-.5 and g <= 1 and top >= 10 and av >= 99
    if ok:
        eligible3.append(m)
    width_ratio = statistics.mean(statistics.mean(v[2] for v in table[m][c].values())/statistics.mean(v[2] for v in table["linear_mi"][c].values()) for c in cells)
    say(f"  {m}: distance {dist(m):.2f} against linear MI {dist('linear_mi'):.2f}, option {dist('random_intercept'):.2f}, oracle {dist('oracle'):.2f}; RMSE ratio {g:.3f}; at most 98 on {top} of 12; available {av:.2f}%; width over linear MI {width_ratio:.3f}: {'ELIGIBLE' if ok else 'not eligible'}")
for (m, text), n in sorted(messages.items()):
    say(f"    {m}: {n} x {text}")
say(f"  nb_mi fell back to Poisson on {100*statistics.mean(fallback.get('nb_mi', [False])):.1f}% of requests")
say(f"OUTCOME P3: {', '.join(eligible3) if eligible3 else 'none eligible'}")
say()
pooled = [r for label in labels for r in part1[label] if r["D"] not in ("", "NA") and r["cover_cv3"] not in ("", "NA")]
flag = [float(r["cover_cv3"]) for r in pooled if float(r["D"]) > .5]
rest = [float(r["cover_cv3"]) for r in pooled if float(r["D"]) <= .5]
heavy_rows = [r for label in heavy for r in part1[label] if r["D"] not in ("", "NA")]
share_heavy = 100*statistics.mean(float(r["D"]) > .5 for r in heavy_rows)
part3 = rows("part3.csv.gz")
false_alarm = 100*statistics.mean(float(r["D"]) > .5 for r in part3)
real3 = [r for r in part3 if r["kind"] == "real"]
sim3 = [r for r in part3 if r["kind"] == "simulation"]
cover_flag = 100*statistics.mean(flag) if flag else float("nan")
cover_rest = 100*statistics.mean(rest)
adopted = bool(flag) and cover_flag <= cover_rest-5 and share_heavy >= 25 and false_alarm <= 5
say(f"P4: flagged requests in Part 1 {len(flag)} of {len(pooled)}; coverage under cv3 flagged {cover_flag:.1f}, unflagged {cover_rest:.1f}; flagged in rare shape 0.5 cells {share_heavy:.1f}%; flagged in Part 3 {false_alarm:.2f}% (real items {100*statistics.mean(float(r['D']) > .5 for r in real3):.2f}%, simulated {100*statistics.mean(float(r['D']) > .5 for r in sim3):.2f}%)")
say(f"OUTCOME P4: {'ADOPTED' if adopted else 'NOT ADOPTED'}")
say()
say("PREDICTIONS")
predictions = [
    ("P1 is HEAVY TAILS", outcome_p1 == "HEAVY TAILS"),
    ("No method in P2 is eligible, with boot_t the closest", not eligible and statistics.mean(cover(c, "boot") for c in heavy) >= max(statistics.mean(cover(c, n) for c in heavy) for n in ("cr2", "wcr", "poisson"))),
    ("Neither candidate in P3 is eligible", not eligible3),
    ("The warning in P4 is adopted", adopted),
]
for i, (text, held) in enumerate(predictions, 1):
    say(f"  {i}. {text}: {'HELD' if held else 'FAILED'}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
print(text)
