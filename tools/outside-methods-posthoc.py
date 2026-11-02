#!/usr/bin/env python3
"""Readings of the comparison with outside methods made after its results. Post hoc: they enter no rule. Plain Python, from saved results.
Usage: python3 -I tools/outside-methods-posthoc.py [OUTPUT_FILE]"""
import csv
import gzip
import math
from pathlib import Path
import random
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
SEED = 20261102
RESAMPLES = 2000
OUT = []


def say(text=""):
    OUT.append(text)


def load(path, method_field, cover, width, se):
    out = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        for r in csv.DictReader(stream):
            if r["status"] != "ok":
                continue
            out[(r[method_field], r["analysis"], float(r["fraction"]), int(r["replicate"]))] = {
                "error": float(r["estimate"])-float(r["target"]), "estimate": float(r["estimate"]), "cover": float(r[cover]), "width": float(r[width]), "se": float(r[se]), "logged": int(r.get("logged") or 0), "between": float(r.get("between") or 0)}
    return out


def failed_rows(path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        return {(r["method"], r["analysis"], float(r["fraction"]), int(r["replicate"])) for r in csv.DictReader(stream) if r["status"] != "ok"}


study = load(ROOT/"validation/outside-methods-study/replicates.csv.gz", "method", "cover", "width", "se")
rerun = load(ROOT/"validation/outside-methods-sensitivity/replicates.csv.gz", "method", "cover", "width", "se")
earlier = load(ROOT/"validation/variance-mi-study/replicates.csv.gz", "imputer", "cover_cv3", "width_cv3", "se_cv3")
rerun_failed = failed_rows(ROOT/"validation/outside-methods-sensitivity/replicates.csv.gz")


def with_rerun(method, analysis, fraction, replicate):
    key = (method, analysis, fraction, replicate)
    if key in rerun_failed:
        return None
    return rerun.get(key, study.get(key))


def rmse(values):
    return math.sqrt(sum(v*v for v in values)/len(values))


say("Post hoc: made after the results, enters no rule")
say()
say("THE MICE DEFAULT, mice_2l_lmer at 25%: requests where mice dropped a predictor; RMSE and coverage on those cohorts under the default, with eps = 0, and for the option")
for item in ITEMS:
    dropped = [r for r in range(1, 2001) if study.get(("mice_2l_lmer", item, .25, r), {}).get("logged", 0) > 0]
    if len(dropped) < 50:
        say(f"  {item}: {len(dropped)} of 2,000")
        continue
    usable = [r for r in dropped if ("mice_2l_lmer", item, .25, r) in rerun and ("random_intercept", item, .25, r) in earlier]
    default = [study[("mice_2l_lmer", item, .25, r)] for r in usable]
    fixed = [rerun[("mice_2l_lmer", item, .25, r)] for r in usable]
    option = [earlier[("random_intercept", item, .25, r)] for r in usable]
    line = "; ".join(f"{name} RMSE {rmse([v['error'] for v in group]):.4f}, coverage {100*statistics.mean(v['cover'] for v in group):.1f}" for name, group in (("default", default), ("eps = 0", fixed), ("option", option)))
    say(f"  {item}: {len(dropped)} of 2,000, compared on {len(usable)}: {line}")
say()
say("OPTION AGAINST MICE_2L_LMER WITH EPS = 0, real items at 25%: coverage difference in points with Monte Carlo error, mean width ratio, RMSE ratio")
pairs = {}
for item in ITEMS:
    both = [r for r in range(1, 2001) if ("random_intercept", item, .25, r) in earlier and with_rerun("mice_2l_lmer", item, .25, r) is not None]
    a = [earlier[("random_intercept", item, .25, r)] for r in both]
    b = [with_rerun("mice_2l_lmer", item, .25, r) for r in both]
    pairs[item] = (a, b)
    gaps = [x["cover"]-y["cover"] for x, y in zip(a, b)]
    say(f"  {item}: {100*statistics.mean(gaps):+.2f} (MC error {100*statistics.stdev(gaps)/math.sqrt(len(gaps)):.2f}); option width over mice {statistics.mean(x['width'] for x in a)/statistics.mean(y['width'] for y in b):.3f}; mice RMSE over option {rmse([y['error'] for y in b])/rmse([x['error'] for x in a]):.3f}")


def distance(group):
    return abs(100*statistics.mean(v["cover"] for v in group)-95)


observed = statistics.mean(distance(pairs[i][1]) for i in ITEMS)-statistics.mean(distance(pairs[i][0]) for i in ITEMS)
rng = random.Random(SEED)
draws = []
for _ in range(RESAMPLES):
    gap = 0
    for item in ITEMS:
        a, b = pairs[item]
        pick = [rng.randrange(len(a)) for _ in a]
        gap += distance([b[j] for j in pick])-distance([a[j] for j in pick])
    draws.append(gap/len(ITEMS))
draws.sort()
say(f"  Mean distance from 95, mice minus option: {observed:.2f} points; cohort resampling, {RESAMPLES} draws, 2.5% to 97.5%: {draws[int(.025*RESAMPLES)]:.2f} to {draws[int(.975*RESAMPLES)]:.2f}. This treats the 11 items as fixed.")
say()
with gzip.open(ROOT/"validation/outside-methods-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    lmer = [int(r["singular"]) > 0 for r in csv.DictReader(stream) if r["method"] == "mice_2l_lmer"]
say(f"MICE_2L_LMER, first run: requests with a singular fit in at least one imputation, {100*sum(lmer)/len(lmer):.1f}%")
say()
lost = {}
for key in rerun_failed:
    if key[0] == "mice_2l_lmer":
        lost[(key[1], key[2])] = lost.get((key[1], key[2]), 0)+1
say("MICE_2L_LMER WITH EPS = 0, requests with no interval: "+"; ".join(f"{a} {int(f*100)}% {n}" for (a, f), n in sorted(lost.items())))
say()
say("OPTION AGAINST MICE_2L_LMER, AUDITED VALUE AS PREDICTOR, first run, paired: bias difference in percent of the target, with Monte Carlo error, and RMSE ratio option over mice")
for cell in ("covariate_shift", "covariate_under", "covariate_row_only"):
    for fraction in (.25, .5):
        both = [r for r in range(1, 501) if ("random_intercept", cell, fraction, r) in earlier and ("mice_2l_lmer", cell, fraction, r) in study]
        a = [earlier[("random_intercept", cell, fraction, r)]["error"] for r in both]
        b = [study[("mice_2l_lmer", cell, fraction, r)]["error"] for r in both]
        gaps = [100*(x-y)/.5 for x, y in zip(a, b)]
        say(f"  {cell} {int(fraction*100)}%: {statistics.mean(gaps):+.1f} (MC error {statistics.stdev(gaps)/math.sqrt(len(gaps)):.1f}); {rmse(a)/rmse(b):.3f}")
say()
say("EVERY AUDITED ROW AGREED, real items at 25%, so recorded values were kept: share of requests, and coverage on them against the rest. Linear MI and the option under cv3")
for item in ITEMS:
    parts = []
    for name in ("linear_mi", "random_intercept"):
        rows_all = [earlier[(name, item, .25, r)] for r in range(1, 2001) if (name, item, .25, r) in earlier]
        kept = [v for v in rows_all if v["between"] == 0]
        rest = [v for v in rows_all if v["between"] != 0]
        parts.append(f"{name} {100*len(kept)/len(rows_all):.1f}%"+(f", coverage {100*statistics.mean(v['cover'] for v in kept):.1f} against {100*statistics.mean(v['cover'] for v in rest):.1f}" if kept else ""))
    say(f"  {item}: "+"; ".join(parts))
say()
say("RECORDED OVER THE BEST AUDIT-USING METHOD, real items at 25%, RMSE ratio with a cohort-resampling range, 2.5% to 97.5%")
AUDIT = [("audited_only", study), ("difference", study), ("raking", study), ("mice_2l_lmer", study), ("mice_2l_pan", study), ("random_intercept", earlier), ("linear_mi", earlier)]
for item in ITEMS:
    reps = list(range(1, 2001))
    errs = {name: {r: source[(name, item, .25, r)]["error"] for r in reps if (name, item, .25, r) in source} for name, source in AUDIT}
    errs["recorded"] = {r: study[("recorded", item, .25, r)]["error"] for r in reps}
    def ratio(sample):
        best = min(rmse([errs[name][r] for r in sample if r in errs[name]]) for name, _ in AUDIT)
        return rmse([errs["recorded"][r] for r in sample])/best
    resampled = sorted(ratio([reps[rng.randrange(len(reps))] for _ in reps]) for _ in range(400))
    say(f"  {item}: {ratio(reps):.3f} ({resampled[10]:.3f} to {resampled[389]:.3f})")
say()
say("RECORDED AGAINST THE OPTION, real items: RMSE ratio and recorded coverage")
for fraction in (.25, .5):
    say(f"  {int(fraction*100)}%:")
    for item in ITEMS:
        reps = [r for r in range(1, 2001) if ("recorded", item, fraction, r) in study and ("random_intercept", item, fraction, r) in earlier]
        ratio = rmse([study[("recorded", item, fraction, r)]["error"] for r in reps])/rmse([earlier[("random_intercept", item, fraction, r)]["error"] for r in reps])
        say(f"    {item}: {ratio:.3f}, coverage {100*statistics.mean(study[('recorded', item, fraction, r)]['cover'] for r in reps):.1f}")
say()
say("STANDARD ERROR OVER THE SPREAD OF THE ESTIMATE, simulated shift cells")
for cell in ("none_shift", "rwanda_shift", "pairs_shift"):
    for fraction in (.25, .5):
        parts = []
        for name, source in (("difference", study), ("raking, full", study), ("raking, approx", rerun), ("audited_only", study)):
            method = name.split(",")[0]
            values = [source[(method, cell, fraction, r)] for r in range(1, 501) if (method, cell, fraction, r) in source]
            parts.append(f"{name} {statistics.mean(v['se'] for v in values)/statistics.stdev(v['estimate'] for v in values):.2f}")
        say(f"  {cell} {int(fraction*100)}%: "+"; ".join(parts))
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
