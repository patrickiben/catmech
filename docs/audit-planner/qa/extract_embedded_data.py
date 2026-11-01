#!/usr/bin/env python3
"""Derive the page's embedded data from the study files (read-only) and, with --check, compare it
to the JSON embedded in ../index.html. Run from anywhere: python3 extract_embedded_data.py [--check]"""
import csv, json, math, pathlib, re, sys
HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STUDY = ROOT / "validation/audit-plan-study"
csv.field_size_limit(10**8)

def rows(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))

def q7(values, p):
    v = sorted(values)
    h = (len(v) - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])

def sig(x, digits=8):
    return float(f"{x:.{digits}g}")

planner = json.load(open(STUDY / "planner-data.json"))
summary = rows(STUDY / "summary.csv")
targets = rows(STUDY / "targets.csv")
plans = rows(STUDY / "plans.csv")
frozen = rows(ROOT / "validation/real-data-study/summary.csv")
# plans.csv repeats each plan once per target fraction, so keep one row per plan.
seen, unique_plans = set(), []
for r in plans:
    key = (r["cell"], r["replicate"], r["pilot_records"])
    if key not in seen:
        seen.add(key)
        unique_plans.append(r)
SHORT = {"ctgov": "ClinicalTrials.gov", "nhanes": "NHANES"}

data = {"verdict": planner["verdict"], "package_version": planner["package_version"], "datasets": {}, "accuracy": [], "methods": {}, "boundary": {}}
for name, block in planner["datasets"].items():
    single = []
    for f in (0.05, 0.1, 0.25, 0.5):
        ratios = [float(r["ratio"]) for r in targets if r["cell"] == name and float(r["target_fraction"]) == f and r["status"] == "ok"]
        single.append({"f": f, "records": round(f * 600), "ok": len(ratios), "q05": sig(q7(ratios, .05), 6), "q25": sig(q7(ratios, .25), 6), "q50": sig(q7(ratios, .5), 6),
                       "q75": sig(q7(ratios, .75), 6), "q95": sig(q7(ratios, .95), 6), "mean": sig(sum(ratios) / len(ratios), 6)})
    answers = {}
    for pilot in ("60", "150"):
        vals = [float(r["answer_records"]) for r in unique_plans if r["cell"] == name and r["pilot_records"] == pilot and r["status"] == "ok"]
        answers[pilot] = {"med": q7(vals, .5), "q05": q7(vals, .05), "q95": q7(vals, .95), "plans": len(vals)}
    data["datasets"][name] = {
        "short": SHORT[name], "label": block["label"], "pool_rows": block["pool_rows"], "pool_k": block["pool_k"], "pool_discard_k": block["pool_discard_k"],
        "pool_classical": block["pool_classical"], "pool_discard_classical": block["pool_discard_classical"],
        "frozen": [{"f": x["fraction"], "mi": sig(x["linear_mi_mean_width"]), "ref": sig(x["reference_normal_mean_width"]), "oracle": sig(x["oracle_mean_width"]), "cohorts": x["cohorts"]} for x in block["frozen_study"]["by_fraction"]],
        "pilots": {p: {k: (sig(v) if isinstance(v, float) else v) for k, v in block["pilot_estimates"][p].items()} for p in ("60", "150")},
        "answers": answers, "single": single}
for r in summary:
    if r["cell"] in SHORT:
        data["accuracy"].append({"d": r["cell"], "pilot": int(r["pilot_records"]), "f": float(r["target_fraction"]), "compared": int(r["compared"]),
            "med": sig(float(r["median_error"]), 6), "lo": sig(float(r["median_error_lower"]), 6), "hi": sig(float(r["median_error_upper"]), 6), "p90": sig(float(r["p90_error"]), 6),
            "signed": sig(float(r["median_signed_error"]), 6), "floor": sig(float(r["floor_median_error"]), 6), "retest": sig(float(r["retest_median_error"]), 6),
            "kcov": sig(float(r["interval_coverage"]), 6), "pcov": sig(float(r["predictive_coverage"]), 6)})
METHODS = ["cole_mi", "reference_normal", "contamination", "raking", "naive_recorded", "oracle"]
for name in SHORT:
    data["methods"][name] = {}
    for r in frozen:
        if r["family"] == name and r["method"] in METHODS and (r["scenario"].startswith("audit_") or r["scenario"] == "targeted_audit"):
            data["methods"][name].setdefault(r["scenario"], {})[r["method"]] = {"bias": sig(float(r["bias"]), 5), "cov": sig(float(r["coverage_conditional"]), 5), "width": sig(float(r["mean_width"]), 5), "rmse": sig(float(r["rmse"]), 5), "ret": int(r["returned_intervals"]), "req": int(r["requested"])}
gross = {}
for r in unique_plans:
    if r["status"] == "ok":
        g = gross.setdefault(r["cell"], [0, 0])
        g[1] += 1
        g[0] += "gross_errors" in r["flags"]
coverage = {}
for r in targets:
    if r["cell"] == "error_gross_and_small" and r["status"] == "ok":
        c = coverage.setdefault(r["target_fraction"], [0, 0])
        c[1] += 1
        c[0] += float(r["lower"]) <= 0.8 <= float(r["upper"])
data["boundary"] = {"gross_warning": gross, "gross_cell_slope_coverage": {f: sig(c[0] / c[1], 5) for f, c in sorted(coverage.items(), key=lambda kv: float(kv[0]))},
    "error": [{"cell": r["cell"], "pilot": int(r["pilot_records"]), "f": float(r["target_fraction"]), "med": sig(float(r["median_error"]), 6), "p90": sig(float(r["p90_error"]), 6)} for r in summary if r["cell"] not in SHORT]}
data["vectors"] = planner["js_test_vectors"]

if "--check" in sys.argv:
    html = (HERE.parent / "index.html").read_text()
    embedded = json.loads(re.search(r'<script id="study-data" type="application/json">(.*?)</script>', html, re.S).group(1))
    same = embedded == json.loads(json.dumps(data))
    print("embedded data equals a fresh derivation from the study files:", same)
    sys.exit(0 if same else 1)
print(json.dumps(data, separators=(",", ":")))
