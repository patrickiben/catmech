#!/usr/bin/env python3
"""Readings of the clustered fix study made after its results. Post hoc: they enter no rule. Plain Python, from saved results.
Usage: python3 -I tools/clustered-fix-posthoc.py [OUTPUT_FILE]"""
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT/"validation/clustered-fix-study"
FIVE = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early"]
DEV = ["pnc", "syphilis", "iron", "anc4"]
NINE = FIVE+DEV
METHODS = ["linear_mi", "linear_mi_clustered", "no_strata", "rao_wu", "jackknife", "random_intercept"]
OUT = []


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def say(text=""):
    OUT.append(text)


summary = {(r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"summary.csv")}
ratios = {(r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"ratios.csv")}
grid = rows(STUDY/"simulation-grid.csv")
CELLS = [(g["label"], float(g["fraction"])) for g in grid]
CLUSTERED = [(g["label"], float(g["fraction"])) for g in grid if g["errors"] in ("shift_low", "shift", "under", "outlier")]


def cover(a, f, m):
    return float(summary[(a, f, m)]["coverage"])


def ratio(a, f, m, field="recovery_ratio"):
    return float(ratios[(a, f, m)][field])


say("Post hoc: made after the results, enters no rule")
say()
say("NINE ITEMS AT 25%: mean distance of coverage from 95%, points; lowest coverage; widest mean width against the existing option")
say(f"  {'oracle':20s} {100*sum(abs(cover(a, .25, 'oracle')-.95) for a in NINE)/9:.2f}; lowest {100*min(cover(a, .25, 'oracle') for a in NINE):.1f}%")
for m in METHODS:
    distance = sum(abs(cover(a, .25, m)-.95) for a in NINE)/9
    lowest = min(NINE, key=lambda a: cover(a, .25, m))
    widths = {a: float(summary[(a, .25, m)]["mean_width"])/float(summary[(a, .25, "linear_mi")]["mean_width"]) for a in NINE}
    widest = max(widths, key=widths.get)
    say(f"  {m:20s} {100*distance:.2f}; lowest {100*cover(lowest, .25, m):.1f}% on {lowest}; widest {100*(widths[widest]-1):+.0f}% on {widest}, narrowest {100*(min(widths.values())-1):+.0f}%")
say()
say("NINE ITEMS AT 25%: the replacement against the first option")
faster = [a for a in NINE if ratio(a, .25, "random_intercept") < ratio(a, .25, "linear_mi_clustered")-.0015]
slower = [a for a in NINE if ratio(a, .25, "random_intercept") > ratio(a, .25, "linear_mi_clustered")+.0015]
say(f"  Lower recovery ratio on {len(faster)}: {', '.join(faster)}; higher on {len(slower)}: {', '.join(slower)}; the rest within 0.0015")
closer = [a for a in NINE if abs(cover(a, .25, "random_intercept")-.95) < abs(cover(a, .25, "linear_mi_clustered")-.95)-1e-9]
farther = [a for a in NINE if abs(cover(a, .25, "random_intercept")-.95) > abs(cover(a, .25, "linear_mi_clustered")-.95)+1e-9]
say(f"  Closer to 95% on {len(closer)}: {', '.join(closer)}; farther on {len(farther)}: {', '.join(farther)}")
narrower = [a for a in NINE if float(summary[(a, .25, "random_intercept")]["mean_width"]) < float(summary[(a, .25, "linear_mi_clustered")]["mean_width"])]
say(f"  Narrower mean width on {len(narrower)}: {', '.join(narrower)}")
for a in NINE:
    say(f"  {a:12s} width: existing {float(summary[(a, .25, 'linear_mi')]['mean_width']):.3f}, first {float(summary[(a, .25, 'linear_mi_clustered')]['mean_width']):.3f}, replacement {float(summary[(a, .25, 'random_intercept')]['mean_width']):.3f}, jackknife {float(summary[(a, .25, 'jackknife')]['mean_width']):.3f}")
say()
say("NINE ITEMS AT 50%: replacement and first option, mean distance from 95%, points")
for m in ["linear_mi", "linear_mi_clustered", "random_intercept"]:
    say(f"  {m:20s} {100*sum(abs(cover(a, .5, m)-.95) for a in NINE)/9:.2f}; lowest {100*min(cover(a, .5, m) for a in NINE):.1f}%")
say()
say("SIMULATIONS, 38 CELLS: recovery ratio against the existing option")
for m in METHODS[1:]:
    points = {c: ratio(c[0], c[1], m) for c in CELLS}
    uppers = sum(ratio(c[0], c[1], m, "recovery_upper") <= 1.05 for c in CELLS)
    worst = max(points, key=points.get)
    best = min(points, key=points.get)
    say(f"  {m:20s} point ratio at most 1.05 in {sum(v <= 1.05 for v in points.values())}; upper limit at most 1.05 in {uppers}; highest {points[worst]:.3f} ({worst[0]} {int(worst[1]*100)}%); lowest {points[best]:.3f} ({best[0]} {int(best[1]*100)}%)")
say()
say("SIMULATIONS, 24 CLUSTERED CELLS: lowest and highest coverage; cells under 92%; cells over 98%")
for m in ["oracle"]+METHODS:
    values = [cover(a, f, m) for a, f in CLUSTERED]
    say(f"  {m:20s} {100*min(values):.1f}% to {100*max(values):.1f}%; under 92% in {sum(v < .92 for v in values)}; over 98% in {sum(v > .98 for v in values)}")
say()
say("ALL REQUESTS: widths over five times their cell median; smallest degrees of freedom")
for m in METHODS:
    group = [r for k, r in summary.items() if k[2] == m]
    say(f"  {m:20s} {sum(int(float(r['wide_requests'])) for r in group)} of {sum(int(r['returned_intervals']) for r in group)}; smallest df {min(float(r['smallest_df']) for r in group):.2f}")
say()
say(f"MONTE CARLO ERROR: standard error of coverage near 90% with 500 cohorts, {100*(.9*.1/500)**.5:.2f} points")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
