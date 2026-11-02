#!/usr/bin/env python3
"""Select the candidate and score the decision rule and predictions in validation/CLUSTERED_FIX_PLAN.md from saved results, in plain Python.
Usage: python3 -I tools/clustered-fix-score.py [OUTPUT_FILE] [STUDY_DIRECTORY]"""
import csv
import gzip
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/"validation/clustered-fix-study"
FIVE = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early"]
DEV = ["pnc", "syphilis", "iron", "anc4"]
UNUSED = ["ancnewreg", "livebirths"]
NINE = FIVE+DEV
RARE = ["lbw", "resus", "stillbirths"]
COMMON = ["deliveries", "anc1_early"]+DEV
CANDIDATES = ["jackknife", "random_intercept"]
METHODS = ["oracle", "linear_mi", "linear_mi_clustered", "no_strata", "rao_wu", "jackknife", "random_intercept"]
SHORT = {"oracle": "oracle", "linear_mi": "existing", "linear_mi_clustered": "first", "no_strata": "no strata", "rao_wu": "rao-wu", "jackknife": "jackknife", "random_intercept": "random int"}
OUT = []


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def say(text=""):
    OUT.append(text)


summary = {(r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"summary.csv")}
ratios = {(r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"ratios.csv")}
pools = {r["analysis"]: r for r in rows(STUDY/"pool-targets.csv")}
grid = rows(STUDY/"simulation-grid.csv")
LEVEL = {a: "held_out" for a in FIVE} | {a: "development" for a in DEV} | {a: "unused" for a in UNUSED}


def s(analysis, fraction, method, field="coverage"):
    return float(summary[(analysis, fraction, method)][field])


def ratio(analysis, fraction, method, field="recovery_ratio"):
    return float(ratios[(analysis, fraction, method)][field])


CLUSTERED = [(g["label"], float(g["fraction"])) for g in grid if g["errors"] in ("shift_low", "shift", "under", "outlier")]
ALL_CELLS = [(g["label"], float(g["fraction"])) for g in grid]
assert len(CLUSTERED) == 24 and len(ALL_CELLS) == 38


def distance(method):
    return sum(abs(s(a, f, method)-.95) for a, f in CLUSTERED)/len(CLUSTERED)


def mean_ratio(method):
    return sum(ratio(a, f, method) for a, f in ALL_CELLS)/len(ALL_CELLS)


def refused_share(method):
    requested = sum(int(summary[(a, f, method)]["requested"]) for a, f in ALL_CELLS)
    returned = sum(int(summary[(a, f, method)]["returned_intervals"]) for a, f in ALL_CELLS)
    return (requested-returned)/requested


say("Scored from validation/clustered-fix-study by tools/clustered-fix-score.py")
say(f"Requests: {sum(int(r['requested']) for r in summary.values())}; returned: {sum(int(r['returned_intervals']) for r in summary.values())}; refused: {sum(int(r['refused']) for r in summary.values())}; other errors: {sum(int(r['other_errors']) for r in summary.values())}")
for m in METHODS[1:]:
    say(f"  {m}: returned {sum(int(r['returned_intervals']) for k, r in summary.items() if k[2] == m)} of {sum(int(r['requested']) for k, r in summary.items() if k[2] == m)}")
say()
say("SELECTION, SIMULATIONS ONLY")
eligible = [m for m in CANDIDATES if refused_share(m) <= .01]
for m in CANDIDATES:
    say(f"  {m}: not returned {100*refused_share(m):.2f}% of simulated requests; mean distance from 95% over 24 clustered cells {100*distance(m):.2f} points; mean recovery ratio over 38 cells {mean_ratio(m):.3f}")
for m in ["linear_mi", "linear_mi_clustered", "no_strata", "rao_wu"]:
    say(f"  {m}: mean distance from 95% over 24 clustered cells {100*distance(m):.2f} points" + ("" if m == "linear_mi" else f"; mean recovery ratio over 38 cells {mean_ratio(m):.3f}"))
if len(eligible) == 2:
    a, b = eligible
    if abs(distance(a)-distance(b)) <= .005 + 1e-12:
        selected = min(eligible, key=mean_ratio)
        how = "within 0.5 points, so the lower mean recovery ratio decides"
    else:
        selected = min(eligible, key=distance)
        how = "smaller mean distance from 95%"
elif eligible:
    selected, how = eligible[0], "the other refused more than 1%"
else:
    selected, how = None, "both refused more than 1%"
say(f"  SELECTED: {selected} ({how})")
say()
if selected is None:
    say("VERDICT: NO-GO")
    print("\n".join(OUT))
    sys.exit(0)

say(f"CONDITIONS FOR {selected} AT A 25% AUDIT ON REAL ITEMS")
passes = {}
for a in NINE:
    c, o = s(a, .25, selected), s(a, .25, "oracle")
    passes[a] = (c >= .90 or c >= o) and c <= .98
    say(f"  {a:12s} coverage {100*c:.1f}%, oracle {100*o:.1f}%, existing {100*s(a, .25, 'linear_mi'):.1f}%, first option {100*s(a, .25, 'linear_mi_clustered'):.1f}%; within limits {passes[a]}; recovery ratio {ratio(a, .25, selected):.3f} [{ratio(a, .25, selected, 'recovery_lower'):.3f}, {ratio(a, .25, selected, 'recovery_upper'):.3f}]; first option {ratio(a, .25, 'linear_mi_clustered'):.3f}")
held = sum(passes.values())
coverage_holds = held >= 8 and passes["anc1_early"]
say(f"  Coverage holds: {coverage_holds} ({held} of 9 within limits, need 8 with anc1_early among them)")
accurate = sum(ratio(a, .25, selected, "recovery_upper") <= 1.05 for a in FIVE)
accuracy_holds = accurate >= 4
say(f"  Accuracy holds: {accuracy_holds} ({accurate} of 5 with upper limit at most 1.05, need 4)")
unused_each = {}
for a in UNUSED:
    near = abs(s(a, .25, selected)-s(a, .25, "oracle")) <= .05
    light = ratio(a, .25, selected, "recovery_upper") <= 1.05
    unused_each[a] = near and light
    say(f"  {a:12s} unused item: rows {pools[a]['rows']}, wrong {pools[a]['wrong_rows']}; coverage {100*s(a, .25, selected):.1f}%, oracle {100*s(a, .25, 'oracle'):.1f}%, existing {100*s(a, .25, 'linear_mi'):.1f}%; within 5 points {near}; recovery ratio {ratio(a, .25, selected):.3f} [{ratio(a, .25, selected, 'recovery_lower'):.3f}, {ratio(a, .25, selected, 'recovery_upper'):.3f}], at most 1.05 {light}")
unused_holds = all(unused_each.values())
say(f"  Unused items hold: {unused_holds}")
say()
say("CONDITIONS IN THE SIMULATIONS, 25% AND 50%")
ok_cells = 0
for a, f in CLUSTERED:
    c, o = s(a, f, selected), s(a, f, "oracle")
    good = (c >= .92 or c >= o) and c <= .98
    ok_cells += good
    if not good:
        say(f"  Outside limits: {a} {int(f*100)}%: {100*c:.1f}%, oracle {100*o:.1f}%")
simulation_coverage = ok_cells >= .8*len(CLUSTERED)
say(f"  Simulation coverage holds: {simulation_coverage} ({ok_cells} of {len(CLUSTERED)} clustered cells, need 80%)")
ok_ratio = 0
for a, f in ALL_CELLS:
    upper = ratio(a, f, selected, "recovery_upper")
    ok_ratio += upper <= 1.05
    if upper > 1.05:
        say(f"  Upper limit above 1.05: {a} {int(f*100)}%: {ratio(a, f, selected):.3f} [{ratio(a, f, selected, 'recovery_lower'):.3f}, {upper:.3f}]")
simulation_accuracy = ok_ratio >= .8*len(ALL_CELLS)
say(f"  Simulation accuracy holds: {simulation_accuracy} ({ok_ratio} of {len(ALL_CELLS)} cells, need 80%)")
verdict = "GO" if coverage_holds and accuracy_holds and unused_holds and simulation_coverage and simulation_accuracy else "PARTIAL" if coverage_holds else "NO-GO"
say(f"VERDICT: {verdict}")
say()
say("THE FIRST TEST'S COVERAGE CONDITION, FOR COMPARISON ONLY, ENTERS NO VERDICT")
for m in ["linear_mi_clustered", "no_strata", "rao_wu", "jackknife", "random_intercept"]:
    closer = sum(abs(s(a, .25, m)-s(a, .25, "oracle")) < abs(s(a, .25, "linear_mi")-s(a, .25, "oracle")) for a in FIVE)
    within = sum(abs(s(a, .25, m)-s(a, .25, "oracle")) <= .05 for a in FIVE)
    limits = sum((s(a, .25, m) >= .90 or s(a, .25, m) >= s(a, .25, "oracle")) and s(a, .25, m) <= .98 for a in NINE)
    upper = sum(ratio(a, .25, m, "recovery_upper") <= 1.05 for a in FIVE)
    say(f"  {m:20s} closer to oracle {closer} of 5; within 5 points of oracle {within} of 5; this plan's limits {limits} of 9; accuracy upper limit {upper} of 5")
say()
say("PREDICTIONS")
say(f"  1. Simulations select random_intercept: {selected == 'random_intercept'} ({selected})")
stratified = [(a, f) for a, f in CLUSTERED if f == .25 and not a.startswith("none")]
over = sum(s(a, f, "rao_wu") > .96 for a, f in stratified)
say(f"  2. rao_wu above 96% in most clustered cells of the stratified designs at 25%: {over > len(stratified)/2} ({over} of {len(stratified)}); range {100*min(s(a, f, 'rao_wu') for a, f in stratified):.1f}% to {100*max(s(a, f, 'rao_wu') for a, f in stratified):.1f}%")
say(f"  3. Both candidates closer to 95% than the first option on average: {all(distance(m) < distance('linear_mi_clustered') for m in CANDIDATES)} (" + ", ".join(f"{m} {100*distance(m):.2f}" for m in CANDIDATES+["linear_mi_clustered"]) + ")")
say(f"  4. anc1_early at 25% at least 90%: {s('anc1_early', .25, selected) >= .90} ({100*s('anc1_early', .25, selected):.1f}%)")
say(f"  5. Rare items, recovery ratio below the first option's on all three: {all(ratio(a, .25, selected) < ratio(a, .25, 'linear_mi_clustered') for a in RARE)} (" + ", ".join(f"{a} {ratio(a, .25, selected):.3f} against {ratio(a, .25, 'linear_mi_clustered'):.3f}" for a in RARE) + ")")
say(f"  6. Common items, recovery ratio at most 1.02 on all six: {all(ratio(a, .25, selected) <= 1.02 for a in COMMON)} (" + ", ".join(f"{a} {ratio(a, .25, selected):.3f}" for a in COMMON) + ")")
controls = {f: s("control_one_row", f, selected, "mean_width")/s("control_one_row", f, "linear_mi", "mean_width") for f in (.25, .5)}
say(f"  7. One-row controls, mean width within 5%: {all(abs(v-1) <= .05 for v in controls.values())} (" + ", ".join(f"{int(f*100)}% {100*(v-1):+.1f}%" for f, v in controls.items()) + ")")
say(f"  8. Unused items at 25%, coverage within 3 points of linear_mi: {all(abs(s(a, .25, selected)-s(a, .25, 'linear_mi')) <= .03 for a in UNUSED)} (" + ", ".join(f"{a} {100*(s(a, .25, selected)-s(a, .25, 'linear_mi')):+.1f}" for a in UNUSED) + ")")
wider = {a: s(a, .25, selected, "mean_width")/s(a, .25, "linear_mi", "mean_width") for a in FIVE}
say(f"  9. Five items at 25%, mean width at most 40% wider: {all(v <= 1.40 for v in wider.values())} (" + ", ".join(f"{a} {100*(v-1):+.0f}%" for a, v in wider.items()) + ")")
wide = {m: sum(int(float(r["wide_requests"])) for k, r in summary.items() if k[2] == m and r["wide_requests"] not in ("", "NA")) for m in METHODS}
say(f"  10. Fewer widths over five times the cell median than the first option: {wide[selected] < wide['linear_mi_clustered']} ({selected} {wide[selected]}, first option {wide['linear_mi_clustered']})")
say()
say("REAL ITEMS: coverage, percent; mean width; recovery RMSE ratio against existing with 95% limits")
for level, names in [("five", FIVE), ("development", DEV), ("unused", UNUSED)]:
    say(f"  {level}")
    for a in names:
        p = pools[a]
        say(f"    {a}: rows {p['rows']}, facilities {p['units']}, wrong {p['wrong_rows']}, target {float(p['target_coefficient']):.3f}, reported {float(p['recorded_coefficient']):.3f}")
        for f in (.25, .5):
            say(f"      {int(f*100)}% coverage: " + "; ".join(f"{SHORT[m]} {100*s(a, f, m):.1f}" for m in METHODS))
            say(f"      {int(f*100)}% width: " + "; ".join(f"{SHORT[m]} {s(a, f, m, 'mean_width'):.3f}" for m in METHODS))
            say(f"      {int(f*100)}% ratio: " + "; ".join(f"{SHORT[m]} {ratio(a, f, m):.3f} [{ratio(a, f, m, 'recovery_lower'):.3f}, {ratio(a, f, m, 'recovery_upper'):.3f}]" for m in METHODS[2:]))
            say(f"      {int(f*100)}% returned: " + "; ".join(f"{SHORT[m]} {summary[(a, f, m)]['returned_intervals']}" for m in METHODS) + f"; audited clusters {s(a, f, 'linear_mi', 'mean_audited_clusters'):.1f}")
say()
say("SIMULATIONS: coverage, percent, then recovery ratio against existing, then mean width")
for a, f in ALL_CELLS:
    say(f"  {a:22s} {int(f*100)}% clusters {s(a, f, 'linear_mi', 'mean_audited_clusters'):4.1f}: " + "; ".join(f"{SHORT[m]} {100*s(a, f, m):.1f}" for m in METHODS))
    say(f"  {'':22s}     ratio: " + "; ".join(f"{SHORT[m]} {ratio(a, f, m):.3f}" for m in METHODS[2:]))
    say(f"  {'':22s}     width: " + "; ".join(f"{SHORT[m]} {s(a, f, m, 'mean_width'):.3f}" for m in METHODS))
say()
say("MEANS BY GROUP: coverage over cells, percent")
GROUPS = [("clustered, rwanda design", lambda a: a.startswith("rwanda_") and not a.endswith(("rare", "row_only"))), ("clustered, pairs design", lambda a: a.startswith("pairs_") and not a.endswith(("rare", "row_only"))),
          ("clustered, no strata", lambda a: a.startswith("none_") and not a.endswith(("rare", "row_only"))), ("rare counts", lambda a: a.endswith("rare")), ("row-only errors", lambda a: a.endswith("row_only")), ("one-row controls", lambda a: a == "control_one_row")]
for name, test in GROUPS:
    for f in (.25, .5):
        keys = [(a, g) for a, g in ALL_CELLS if test(a) and g == f]
        say(f"  {name:26s} {int(f*100)}%: " + "; ".join(f"{SHORT[m]} {100*sum(s(a, g, m) for a, g in keys)/len(keys):.1f}" for m in METHODS))
say()
say("BETWEEN-IMPUTATION VARIANCE: median over the cell, divided by the squared recovery RMSE, on clustered cells at 25% and 50%")
for a, f in CLUSTERED:
    say(f"  {a:22s} {int(f*100)}%: " + "; ".join(f"{SHORT[m]} {s(a, f, m, 'median_between')/s(a, f, m, 'recovery_rmse')**2:.2f}" for m in METHODS[1:]))
say()
say("EXTREME WIDTHS: requests over five times their cell median, and smallest degrees of freedom")
detail = {}
with gzip.open(STUDY/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
    for r in csv.DictReader(stream):
        if r["status"] == "ok":
            d = detail.setdefault((r["level"] == "simulation", r["method"]), [0, 0, 1e18])
            d[1] += 1
            d[2] = min(d[2], float(r["df"]))
for (simulated, m), d in sorted(detail.items()):
    d[0] = sum(int(float(r["wide_requests"])) for k, r in summary.items() if k[2] == m and (r["level"] == "simulation") == simulated)
    say(f"  {'simulation' if simulated else 'real':10s} {m:20s} over 5x median {d[0]} of {d[1]}; smallest df {d[2]:.2f}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
