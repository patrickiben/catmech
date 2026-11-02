#!/usr/bin/env python3
"""Score the decision rule and predictions in validation/VARIANCE_MI_PLAN.md from saved results, in plain Python, and check the controls.
Usage: python3 -I tools/variance-mi-score.py [OUTPUT_FILE] [DIRECTORY]"""
import csv
import gzip
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
FOLDER = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/"validation/variance-mi-study"
STUDY = FOLDER.name == "variance-mi-study"
ITEMS = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early", "pnc", "syphilis", "iron", "anc4", "ancnewreg", "livebirths"]
VARIANCES = ["cr1s", "cr2", "cv3"]
IMPUTERS = ["oracle", "linear_mi", "random_intercept"]
OUT = []


def say(text=""):
    OUT.append(text)


with gzip.open(FOLDER/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    records = list(csv.DictReader(stream))
groups = {}
for r in records:
    groups.setdefault((r["level"], r["analysis"], float(r["fraction"]), r["imputer"]), []).append(r)
good = {k: [r for r in g if r["status"] == "ok"] for k, g in groups.items()}


def real_key(a, f, imp):
    level = "held_out" if a in ITEMS[:5] else "development" if a in ITEMS[5:9] else "unused"
    return (level, a, f, imp)


def cover(key, v):
    g = good[key]
    if v == "wcr":
        return sum(float(r["wcr_p"]) > .05 for r in g)/len(g)
    return sum(float(r[f"cover_{v}"]) for r in g)/len(g)


def width(key, v):
    g = good[key]
    return sum(float(r[f"width_{v}"]) for r in g)/len(g)


def mean_field(key, field):
    return statistics.mean(float(r[field]) for r in good[key])


say("Scored from validation/variance-mi-study by tools/variance-mi-score.py")
statuses = {}
for r in records:
    statuses[r["status"]] = statuses.get(r["status"], 0)+1
say(f"Requests: {len(records)}; " + ", ".join(f"{k} {v}" for k, v in sorted(statuses.items())))
if STUDY:
    with gzip.open(ROOT/"validation/clustered-fix-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        saved = {(r["analysis"], float(r["fraction"]), r["method"], r["replicate"]): r for r in csv.DictReader(stream) if r["method"] in IMPUTERS}
    largest, count = 0.0, 0
    for r in records:
        s = saved.get((r["analysis"], float(r["fraction"]), r["imputer"], r["replicate"]))
        if s is None or r["level"] == "imputed_covariate":
            continue
        assert s["status"] == r["status"], (r["analysis"], r["replicate"], r["imputer"])
        count += 1
        if r["status"] == "ok":
            largest = max(largest, abs(float(r["estimate"])-float(s["estimate"])), abs(float(r["se_cr1s"])-float(s["se"])), abs(float(r["df_cr1s"])-float(s["df"])))
    assert largest < 1e-8, largest
    say(f"Control: {count} requests reproduce the replacement study under the current variance, largest difference {largest:.2e}.")
    with gzip.open(ROOT/"validation/variance-check/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        earlier = {(r["analysis"], float(r["fraction"]), r["replicate"]): r for r in csv.DictReader(stream)}
    largest, count = 0.0, 0
    for r in records:
        s = earlier.get((r["analysis"], float(r["fraction"]), r["replicate"]))
        if s is None or r["imputer"] != "oracle" or r["level"] == "imputed_covariate":
            continue
        count += 1
        largest = max(largest, abs(float(r["se_cr2"])-float(s["se_cr2"])), abs(float(r["se_cv3"])-float(s["se_cv3"])), abs(float(r["df_cr2"])-float(s["cr2_df"])))
    assert largest < 1e-8, largest
    say(f"Control: {count} corrected recount requests reproduce the earlier scoring, largest difference {largest:.2e}.")
say()
say("REAL ITEMS AT 25%, RANDOM INTERCEPT AND LINEAR MI: coverage, percent; mean width; pooled df under cr2")
for a in ITEMS:
    ri, li, orc = (real_key(a, .25, m) for m in ["random_intercept", "linear_mi", "oracle"])
    say(f"  {a:12s} n {len(good[ri])}: random intercept " + ", ".join(f"{v} {100*cover(ri, v):.1f}" for v in VARIANCES+["wcr"]) + "; linear MI " + ", ".join(f"{v} {100*cover(li, v):.1f}" for v in VARIANCES+["wcr"]) + "; oracle " + ", ".join(f"{v} {100*cover(orc, v):.1f}" for v in VARIANCES+["wcr"]))
    say(f"  {'':12s} width, random intercept " + ", ".join(f"{v} {width(ri, v):.3f}" for v in VARIANCES) + "; linear MI " + ", ".join(f"{v} {width(li, v):.3f}" for v in VARIANCES) + f"; df cr2 {mean_field(ri, 'df_cr2'):.1f}, cr1s {mean_field(ri, 'df_cr1s'):.1f}")
say()
say("REAL ITEMS AT 50%: coverage, percent")
for a in ITEMS:
    ri, li = (real_key(a, .5, m) for m in ["random_intercept", "linear_mi"])
    say(f"  {a:12s} random intercept " + ", ".join(f"{v} {100*cover(ri, v):.1f}" for v in VARIANCES+["wcr"]) + "; linear MI " + ", ".join(f"{v} {100*cover(li, v):.1f}" for v in VARIANCES+["wcr"]))
say()
SIM = sorted({(k[0], k[1], k[2]) for k in good if k[0] in ("simulation", "imputed_covariate")})
assert len(SIM) == 44
say("SIMULATED CELLS: coverage, percent, under cr1s, cr2, cv3 and wcr")
for level, a, f in SIM:
    say(f"  {a:22s} {int(f*100)}% " + "; ".join(f"{imp} " + ", ".join(f"{100*cover((level, a, f, imp), v):.1f}" for v in VARIANCES+(["wcr"] if level == "simulation" else [])) for imp in IMPUTERS))
say()
say("MEANS BY SETTING: coverage, percent, cr1s, cr2, cv3")
SETTINGS = [("clustered", lambda k: k[0] == "simulation" and k[1].endswith(("shift_low", "shift", "under", "outlier"))), ("rare counts", lambda k: k[0] == "simulation" and k[1].endswith("rare")),
            ("row-only", lambda k: k[0] == "simulation" and (k[1].endswith("row_only") or k[1] == "control_one_row")), ("imputed covariate", lambda k: k[0] == "imputed_covariate")]
for name, test in SETTINGS:
    keys = [k for k in SIM if test(k)]
    say(f"  {name:18s} {len(keys):2d} cells: " + "; ".join(f"{imp} " + ", ".join(f"{100*statistics.mean(cover((k[0], k[1], k[2], imp), v) for k in keys):.1f}" for v in VARIANCES) for imp in IMPUTERS))
say()


def distance(v, f=.25, imp="random_intercept"):
    return statistics.mean(abs(cover(real_key(a, f, imp), v)-.95) for a in ITEMS)


say("MEAN DISTANCE FROM 95% ACROSS THE 11 ITEMS, points")
for imp in IMPUTERS:
    say(f"  {imp:18s} 25%: " + ", ".join(f"{v} {100*distance(v, .25, imp):.2f}" for v in VARIANCES+["wcr"]) + "; 50%: " + ", ".join(f"{v} {100*distance(v, .5, imp):.2f}" for v in VARIANCES+["wcr"]))
say("  Mean recovery is unchanged by the variance, so accuracy results from the replacement study stand.")
say()
say("DECISION RULE, using random_intercept unless stated")
eligible = []
for v in ["cr2", "cv3"]:
    c25 = {a: cover(real_key(a, .25, "random_intercept"), v) for a in ITEMS}
    c50 = {a: cover(real_key(a, .5, "random_intercept"), v) for a in ITEMS}
    r1 = sum(x >= .92 for x in c25.values()) >= 9 and sum(x <= .98 for x in c25.values()) >= 10
    r2 = sum(x >= .92 for x in c50.values()) >= 8 and sum(x <= .98 for x in c50.values()) >= 10
    r3 = distance("cr1s")-distance(v) >= .005
    within = sum((cover((k[0], k[1], k[2], "random_intercept"), v) >= .92 or cover((k[0], k[1], k[2], "random_intercept"), v) >= cover((k[0], k[1], k[2], "oracle"), v)) and cover((k[0], k[1], k[2], "random_intercept"), v) <= .98 for k in SIM)
    s1 = within >= .8*44
    linear_real = sum(cover(real_key(a, .25, "linear_mi"), v) <= .98 for a in ITEMS)
    linear_sim = sum(cover((k[0], k[1], k[2], "linear_mi"), v) <= .98 for k in SIM)
    l1 = linear_real >= 10 and linear_sim >= .9*44
    wider = {a: width(real_key(a, .25, "random_intercept"), v)/width(real_key(a, .25, "random_intercept"), "cr1s") for a in ITEMS}
    w1 = all(x <= 1.30 for x in wider.values())
    ok = r1 and r2 and r3 and s1 and l1 and w1
    if ok:
        eligible.append(v)
    say(f"  {v}:")
    say(f"    R1 at 25%: at least 92% on {sum(x >= .92 for x in c25.values())} of 11, at most 98% on {sum(x <= .98 for x in c25.values())}; holds {r1}; lowest {100*min(c25.values()):.1f}% on {min(c25, key=c25.get)}, highest {100*max(c25.values()):.1f}% on {max(c25, key=c25.get)}")
    say(f"    R2 at 50%: at least 92% on {sum(x >= .92 for x in c50.values())} of 11, at most 98% on {sum(x <= .98 for x in c50.values())}; holds {r2}")
    say(f"    R3: mean distance {100*distance(v):.2f} against {100*distance('cr1s'):.2f}; holds {r3}")
    say(f"    S1: within limits in {within} of 44; holds {s1}")
    say(f"    L1: linear MI at most 98% on {linear_real} of 11 and in {linear_sim} of 44; holds {l1}")
    say(f"    W1: width ratio {100*(min(wider.values())-1):+.0f}% to {100*(max(wider.values())-1):+.0f}%, widest on {max(wider, key=wider.get)}; holds {w1}")
    say(f"    Eligible: {ok}")
if eligible == ["cr2", "cv3"]:
    adopted = "cv3" if distance("cr2")-distance("cv3") > .005 else "cr2"
elif eligible:
    adopted = eligible[0]
else:
    adopted = None
say(f"OUTCOME: {adopted if adopted else 'none'}")
current = [k for k in SIM if (cover((k[0], k[1], k[2], "random_intercept"), "cr1s") >= .92 or cover((k[0], k[1], k[2], "random_intercept"), "cr1s") >= cover((k[0], k[1], k[2], "oracle"), "cr1s")) and cover((k[0], k[1], k[2], "random_intercept"), "cr1s") <= .98]
c25 = {a: cover(real_key(a, .25, "random_intercept"), "cr1s") for a in ITEMS}
say(f"  For comparison, the current variance with random_intercept: at least 92% on {sum(x >= .92 for x in c25.values())} of 11 at 25%, simulation limits in {len(current)} of 44")
say()
say("IMPUTED-COVARIATE CELLS, random_intercept: coverage under cr2 and cv3; under 92% in how many of 6")
cells = [k for k in SIM if k[0] == "imputed_covariate"]
for v in ["cr2", "cv3"]:
    say(f"  {v}: under 92% in {sum(cover((k[0], k[1], k[2], 'random_intercept'), v) < .92 for k in cells)} of 6")
say()
say("PREDICTIONS")
say(f"  1. cr2 is eligible: {'cr2' in eligible}")
wider_cv3 = {a: width(real_key(a, .25, "random_intercept"), "cv3")/width(real_key(a, .25, "random_intercept"), "cr1s") for a in ITEMS}
say(f"  2. cv3 fails W1 on at least one item: {any(x > 1.30 for x in wider_cv3.values())} (largest {100*(max(wider_cv3.values())-1):+.0f}%)")
sy = cover(real_key("syphilis", .25, "random_intercept"), "cr2")
say(f"  3. syphilis below 92% under cr2 at 25%: {sy < .92} ({100*sy:.1f}%)")
re_cv3 = cover(real_key("resus", .25, "random_intercept"), "cv3")
say(f"  4. resus above 98% under cv3 at 25%: {re_cv3 > .98} ({100*re_cv3:.1f}%)")
gain = {imp: statistics.mean(cover(real_key(a, .25, imp), "cr2")-cover(real_key(a, .25, imp), "cr1s") for a in ITEMS) for imp in ["random_intercept", "oracle"]}
say(f"  5. cr2 gains less through imputation than on the recount: {gain['random_intercept'] < gain['oracle']} (random intercept {100*gain['random_intercept']:+.2f}, oracle {100*gain['oracle']:+.2f})")
dfs = {a: mean_field(real_key(a, .25, "random_intercept"), "df_cr2") for a in ITEMS}
say(f"  6. Pooled df under cr2 below 15 on all 11: {all(x < 15 for x in dfs.values())} (range {min(dfs.values()):.1f} to {max(dfs.values()):.1f})")
say()
say("VERY LOW DEGREES OF FREEDOM, random_intercept across all requests: share under 5")
everything = [r for k, g in good.items() if k[3] == "random_intercept" for r in g]
for v in VARIANCES:
    say(f"  {v}: {100*sum(float(r[f'df_{v}']) < 5 for r in everything)/len(everything):.2f}% of {len(everything)}; smallest {min(float(r[f'df_{v}']) for r in everything):.2f}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
