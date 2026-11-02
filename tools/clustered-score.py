#!/usr/bin/env python3
"""Score the decision rule and predictions in validation/CLUSTERED_PLAN.md from saved results, in plain Python.
Usage: python3 -I tools/clustered-score.py [OUTPUT_FILE] [STUDY_DIRECTORY]"""
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/"validation/clustered-study"
HELD = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early"]
DEV = ["pnc", "syphilis", "iron", "anc4"]
METHODS = ["oracle", "recorded", "linear_mi", "linear_mi_clustered"]
OUT = []


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def say(text=""):
    OUT.append(text)


summary = {(r["level"], r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"summary.csv")}
ratios = {(r["level"], r["analysis"], float(r["fraction"])): r for r in rows(STUDY/"ratios.csv")}
pools = {(r["level"], r["analysis"]): r for r in rows(STUDY/"pool-targets.csv")}
grid = rows(STUDY/"simulation-grid.csv")


def s(level, analysis, fraction, method, field):
    return float(summary[(level, analysis, fraction, method)][field])


say("Scored from validation/clustered-study by tools/clustered-score.py")
say(f"Requests: {sum(int(r['requested']) for r in summary.values())}; returned: {sum(int(r['returned_intervals']) for r in summary.values())}; refused: {sum(int(r['refused']) for r in summary.values())}; other errors: {sum(int(r['other_errors']) for r in summary.values())}")
say()
say("HELD-OUT RWANDA ITEMS AT A 25% FACILITY AUDIT")
closer, within, accurate = 0, 0, 0
for a in HELD:
    oracle, clustered, plain = (s("held_out", a, .25, m, "coverage") for m in ["oracle", "linear_mi_clustered", "linear_mi"])
    r = ratios[("held_out", a, .25)]
    is_closer = abs(clustered-oracle) < abs(plain-oracle)
    is_within = abs(clustered-oracle) <= .05
    is_accurate = float(r["recovery_upper"]) <= 1.05
    closer += is_closer
    within += is_within
    accurate += is_accurate
    say(f"  {a:12s} coverage: clustered {100*clustered:.1f}%, package {100*plain:.1f}%, oracle {100*oracle:.1f}%, closer {is_closer}, within 5 points {is_within}; paired difference {100*float(r['coverage_difference']):+.1f} (MCSE {100*float(r['coverage_difference_mcse']):.1f}); recovery ratio {float(r['recovery_ratio']):.3f} [{float(r['recovery_lower']):.3f}, {float(r['recovery_upper']):.3f}], at most 1.05 {is_accurate}; pool ratio {float(r['pool_ratio']):.3f} [{float(r['pool_lower']):.3f}, {float(r['pool_upper']):.3f}]")
coverage_holds = closer >= 4 and within >= 4
accuracy_holds = accurate >= 4
say(f"  Coverage holds: {coverage_holds} (closer {closer} of 5, within 5 points {within} of 5, need 4 of 5 each)")
say(f"  Accuracy holds: {accuracy_holds} ({accurate} of 5, need 4 of 5)")
say()
say("SIMULATIONS WITH 70 OR MORE CLUSTERS, 25% AND 50%: coverage at least 92% and within 3 points of the oracle")
eligible, passing = [], 0
for g in grid:
    if g["rows"] != "3" or int(g["clusters"]) < 70:
        continue
    label = f"g{g['clusters']}_icc{round(100*float(g['icc'])):02d}_{g['errors']}"
    fraction = float(g["fraction"])
    clustered, oracle = s("simulation", label, fraction, "linear_mi_clustered", "coverage"), s("simulation", label, fraction, "oracle", "coverage")
    ok = clustered >= .92 and abs(clustered-oracle) <= .03
    passing += ok
    eligible.append((label, fraction, clustered, oracle, s("simulation", label, fraction, "linear_mi", "coverage"), ok))
simulation_holds = passing >= .8*len(eligible)
say(f"  Passing {passing} of {len(eligible)} cells ({100*passing/len(eligible):.0f}%, need 80%). Simulation holds: {simulation_holds}")
verdict = "GO" if coverage_holds and accuracy_holds and simulation_holds else "PARTIAL" if coverage_holds else "NO-GO"
say(f"VERDICT: {verdict}")
say()
say("PREDICTIONS")
say(f"  1. Closer to the oracle in at least 4 of 5: {closer >= 4} ({closer} of 5)")
widths = {a: s("held_out", a, .25, "linear_mi_clustered", "mean_width")/s("held_out", a, .25, "linear_mi", "mean_width") for a in HELD}
say("  2. Width 15% to 60% wider in all five: " + str(all(1.15 <= w <= 1.60 for w in widths.values())) + "; " + ", ".join(f"{a} {100*(w-1):+.0f}%" for a, w in widths.items()))
say(f"  3. Recovery RMSE no worse than 5% in at least 4 of 5: {accuracy_holds} ({accurate} of 5)")
few = [g for g in grid if g["rows"] == "3" and g["clusters"] == "30" and g["fraction"] == "0.25"]
below = []
for g in few:
    label = f"g30_icc{round(100*float(g['icc'])):02d}_{g['errors']}"
    below.append((label, s("simulation", label, .25, "linear_mi_clustered", "coverage")))
say(f"  4. With 30 clusters at 25%, coverage below 92% in most cells: {sum(v < .92 for _, v in below) > len(below)/2} ({sum(v < .92 for _, v in below)} of {len(below)}); range {100*min(v for _, v in below):.1f}% to {100*max(v for _, v in below):.1f}%")
controls = {f: s("simulation", "control_one_row", f, "linear_mi_clustered", "mean_width")/s("simulation", "control_one_row", f, "linear_mi", "mean_width") for f in (.25, .5)}
say("  5. Control width within 5%: " + str(all(abs(v-1) <= .05 for v in controls.values())) + "; " + ", ".join(f"{int(f*100)}% {100*(v-1):+.1f}%" for f, v in controls.items()))
gaps = {a: (abs(s("held_out", a, .25, "linear_mi_clustered", "coverage")-s("held_out", a, .25, "linear_mi", "coverage")), abs(s("held_out", a, .5, "linear_mi_clustered", "coverage")-s("held_out", a, .5, "linear_mi", "coverage"))) for a in HELD}
say("  6. Coverage gap between methods shrinks from 25% to 50%: " + f"{sum(v[1] < v[0]-1e-9 for v in gaps.values())} of 5, tied on {sum(abs(v[1]-v[0]) <= 1e-9 for v in gaps.values())}; " + ", ".join(f"{a} {100*v[0]:.1f} to {100*v[1]:.1f}" for a, v in gaps.items()))
say()
say("TABLES: coverage, mean width, recovery RMSE, bias, for oracle, recorded, package linear MI, clustered linear MI")
for level, names in [("held_out", HELD), ("development", DEV)]:
    say(f"  {level}")
    for a in names:
        p = pools[(level, a)]
        say(f"    {a}: rows {p['rows']}, facilities {p['units']}, wrong {p['wrong_rows']}, target {float(p['target_coefficient']):.3f}, reported {float(p['recorded_coefficient']):.3f}")
        for f in (.25, .5):
            cells = "; ".join(f"{m} {100*s(level, a, f, m, 'coverage'):.1f}%, {s(level, a, f, m, 'mean_width'):.3f}, {s(level, a, f, m, 'recovery_rmse'):.3f}, {s(level, a, f, m, 'bias'):+.3f}" for m in METHODS)
            say(f"      {int(f*100)}%: {cells}; returned {summary[(level, a, f, 'linear_mi_clustered')]['returned_intervals']}")
say("  simulation: label, fraction, audited clusters, coverage oracle, recorded, package, clustered; width package, clustered; recovery RMSE package, clustered; recorded bias; returned")
sim_keys = sorted({(k[1], k[2]) for k in summary if k[0] == "simulation"})
for label, f in sim_keys:
    v = lambda m, field: s("simulation", label, f, m, field)
    say(f"    {label:24s} {int(f*100)}% {v('linear_mi_clustered', 'mean_audited_clusters'):5.1f}: {100*v('oracle', 'coverage'):.1f}, {100*v('recorded', 'coverage'):.1f}, {100*v('linear_mi', 'coverage'):.1f}, {100*v('linear_mi_clustered', 'coverage'):.1f}; {v('linear_mi', 'mean_width'):.3f}, {v('linear_mi_clustered', 'mean_width'):.3f}; {v('linear_mi', 'recovery_rmse'):.4f}, {v('linear_mi_clustered', 'recovery_rmse'):.4f}; {v('recorded', 'bias'):+.3f}; {summary[('simulation', label, f, 'linear_mi_clustered')]['returned_intervals']}")
say()
say("SIMULATION SUMMARY BY FACTOR: mean coverage of clustered against package")
for name, test in [("30 clusters", lambda k: k[0].startswith("g30_")), ("70 clusters", lambda k: k[0].startswith("g70_")), ("150 clusters", lambda k: k[0].startswith("g150_")),
                   ("iid errors", lambda k: k[0].endswith("iid") and k[0].startswith("g")), ("persistent errors", lambda k: k[0].endswith("persistent")), ("tracking errors", lambda k: k[0].endswith("tracking"))]:
    keys = [k for k in sim_keys if test(k)]
    mean = lambda m: sum(s("simulation", k[0], k[1], m, "coverage") for k in keys)/len(keys)
    say(f"  {name:18s} cells {len(keys):2d}: oracle {100*mean('oracle'):.1f}%, package {100*mean('linear_mi'):.1f}%, clustered {100*mean('linear_mi_clustered'):.1f}%")
say()
say("DEVELOPMENT RECORD AT 25% AND 50%, enters no rule: recovery ratio, clustered against package; paired coverage difference")
for a in DEV:
    for f in (.25, .5):
        r = ratios[("development", a, f)]
        say(f"  {a:9s} {int(f*100)}%: {float(r['recovery_ratio']):.3f} [{float(r['recovery_lower']):.3f}, {float(r['recovery_upper']):.3f}]; {100*float(r['coverage_difference']):+.1f} (MCSE {100*float(r['coverage_difference_mcse']):.1f})")
say()
say("HELD-OUT AT 50%: recovery ratio, clustered against package; paired coverage difference")
for a in HELD:
    r = ratios[("held_out", a, .5)]
    say(f"  {a:12s} {float(r['recovery_ratio']):.3f} [{float(r['recovery_lower']):.3f}, {float(r['recovery_upper']):.3f}]; {100*float(r['coverage_difference']):+.1f} (MCSE {100*float(r['coverage_difference_mcse']):.1f})")
say()
say("EXTREME WIDTHS: where mean width is far from the median, by method and level. Width over 5 times the cell median; smallest degrees of freedom; largest width")
detail = {}
for r in rows(STUDY/"replicates.csv"):
    if r["status"] == "ok" and r["method"] in ("linear_mi", "linear_mi_clustered"):
        detail.setdefault((r["level"], r["analysis"], float(r["fraction"]), r["method"]), []).append((float(r["width"]), float(r["df"])))
totals = {}
for key, values in sorted(detail.items()):
    widths = sorted(v[0] for v in values)
    median = widths[len(widths)//2]
    extreme = sum(v > 5*median for v in widths)
    t = totals.setdefault((key[0], key[3]), [0, 0, 0.0, 1e9])
    t[0] += extreme
    t[1] += len(widths)
    t[2] = max(t[2], widths[-1]/median)
    t[3] = min(t[3], min(v[1] for v in values))
    if extreme:
        say(f"  {key[1]:24s} {int(key[2]*100)}% {key[3]:20s} median {median:.3f}, over 5x median {extreme} of {len(widths)}, largest {widths[-1]:.1f}, smallest df {min(v[1] for v in values):.2f}")
for key, t in sorted(totals.items()):
    say(f"  {key[0]:12s} {key[1]:20s} over 5x median {t[0]} of {t[1]}; largest over median {t[2]:.1f}; smallest df {t[3]:.2f}")
wider = {}
for a in HELD:
    wider[a] = sorted(v[0] for v in detail[("held_out", a, .25, "linear_mi_clustered")])[250]/sorted(v[0] for v in detail[("held_out", a, .25, "linear_mi")])[250]
say("  Median width, clustered over package, held-out at 25%: " + ", ".join(f"{a} {v:.3f}" for a, v in wider.items()))
sim30 = [k for k in sorted({k[1] for k in detail if k[0] == "simulation" and k[1].startswith("g30_")})]
medians = {m: [sorted(v[0] for v in detail[("simulation", k, .25, m)])[250] for k in sim30] for m in ("linear_mi", "linear_mi_clustered")}
say("  Median width over the 30-cluster cells at 25%, clustered over package: " + ", ".join(f"{b/a:.3f}" for a, b in zip(medians["linear_mi"], medians["linear_mi_clustered"])))
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
