#!/usr/bin/env python3
"""Post hoc look at where the clustered option helped, written after the results were seen. It enters no rule.

Usage: python3 -I tools/clustered-posthoc.py
Reads the frozen pool and the saved study, and prints the table behind the post-hoc paragraphs of
validation/CLUSTERED_SYNTHESIS.md. The saved copy is validation/clustered-posthoc.txt.
"""
import csv
import pathlib
import statistics as st

ROOT = pathlib.Path(__file__).resolve().parents[1]
ITEMS = {
    "lbw": ("held_out", ("lbw_notpremature_h", "lbw_premature_h"), ("lbw_r",), "hf_deliveries_h"),
    "resus": ("held_out", ("newborns_nocryres_succ_h",), ("nocry_res_success_r",), "hf_deliveries_h"),
    "stillbirths": ("held_out", ("sbmacerated_h", "sbfresh_h"), ("sbmacerated_r", "sbfresh_r"), "hf_deliveries_h"),
    "deliveries": ("held_out", ("hf_deliveries_h",), ("hf_deliveries_r",), "ancnewreg_h"),
    "anc1_early": ("held_out", ("anc1vis1trim_h",), ("anc1vis_16weeks_r",), "ancnewreg_h"),
    "pnc": ("development", ("pnc_newborn_24h_h",), ("pnc_newborn_24h_r",), "hf_deliveries_h"),
    "syphilis": ("development", ("anc_syphilis_test_h",), ("anc_syphilis_test_r",), "ancnewreg_h"),
    "iron": ("development", ("womenreciron_h",), ("womenreciron_r",), "ancnewreg_h"),
    "anc4": ("development", ("anc_4visits_h",), ("anc_4visits_r",), "ancnewreg_h"),
}


def rows(path):
    return list(csv.DictReader((ROOT/path).open(newline="")))


def residuals(y, x):
    """Least squares with an intercept, by normal equations."""
    k = len(x[0])+1
    a, b = [[0.0]*k for _ in range(k)], [0.0]*k
    for yi, xi in zip(y, x):
        v = [1.0, *xi]
        for i in range(k):
            b[i] += v[i]*yi
            for j in range(k):
                a[i][j] += v[i]*v[j]
    for i in range(k):
        p = max(range(i, k), key=lambda q: abs(a[q][i]))
        a[i], a[p], b[i], b[p] = a[p], a[i], b[p], b[i]
        for q in range(i+1, k):
            f = a[q][i]/a[i][i]
            a[q] = [u-f*w for u, w in zip(a[q], a[i])]
            b[q] -= f*b[i]
    beta = [0.0]*k
    for i in reversed(range(k)):
        beta[i] = (b[i]-sum(a[i][j]*beta[j] for j in range(i+1, k)))/a[i][i]
    return [yi-sum(c*v for c, v in zip(beta, [1.0, *xi])) for yi, xi in zip(y, x)]


def main():
    pool = [r for r in rows("validation/real-audit/rwanda-pool.csv") if r["facility_type"] == "health_center"]
    saved = {r["analysis"]: r for r in rows("validation/clustered-study/pool-targets.csv")}
    ratios = {(r["analysis"], float(r["fraction"])): r for r in rows("validation/clustered-study/ratios.csv") if r["analysis"] in ITEMS}
    summary = {(r["analysis"], float(r["fraction"]), r["method"]): r for r in rows("validation/clustered-study/summary.csv") if r["analysis"] in ITEMS}
    number = lambda r, names: None if any(r[n] == "" for n in names) else sum(float(r[n]) for n in names)
    table = []
    for name, (role, w, reference, covariate) in ITEMS.items():
        data = [(number(r, w), number(r, reference), float(r[covariate])) for r in pool]
        data = [d for d in data if d[0] is not None and d[1] is not None]
        recorded, truth, x = zip(*data)
        wrong = sum(abs(a-b) > 1e-8 for a, b in zip(recorded, truth))
        assert int(saved[name]["rows"]) == len(data) and int(saved[name]["wrong_rows"]) == wrong, name
        spread = st.correlation([abs(e) for e in residuals(truth, list(zip(recorded, x)))], x)
        cell = lambda fraction, method, col="coverage": float(summary[(name, fraction, method)][col])
        table.append(dict(name=name, role=role, n=len(data), recorded=st.mean(recorded), truth=st.mean(truth), zeros=sum(v == 0 for v in truth)/len(data),
                          wrong=wrong/len(data), spread=spread, oracle=cell(.25, "oracle"), existing=cell(.25, "linear_mi"), clustered=cell(.25, "linear_mi_clustered"),
                          existing_50=cell(.5, "linear_mi"), rmse_25=cell(.25, "linear_mi", "recovery_rmse"), rmse_50=cell(.5, "linear_mi", "recovery_rmse"),
                          ratio_25=float(ratios[(name, .25)]["recovery_ratio"]), ratio_50=float(ratios[(name, .5)]["recovery_ratio"]),
                          diff_25=float(ratios[(name, .25)]["coverage_difference"]), diff_50=float(ratios[(name, .5)]["coverage_difference"])))
    table.sort(key=lambda t: t["recorded"])
    print("Post hoc, written after the results were seen. It enters no rule.")
    print("Nine Rwanda items, sorted by mean recorded count per facility-month.")
    print()
    print("item, role, rows: mean recorded count, mean recount, zeros in the recount, rows that disagree, correlation of absolute residual with the covariate;")
    print("  coverage at 25% for oracle, existing, clustered; recovery RMSE ratio at 25%, 50%; coverage difference at 25%, 50%")
    for t in table:
        print(f"{t['name']}, {t['role']}, {t['n']}: {t['recorded']:.1f}, {t['truth']:.1f}, {100*t['zeros']:.0f}%, {100*t['wrong']:.0f}%, {t['spread']:.2f};")
        print(f"  {100*t['oracle']:.1f}, {100*t['existing']:.1f}, {100*t['clustered']:.1f}; {t['ratio_25']:.3f}, {t['ratio_50']:.3f}; {100*t['diff_25']:+.1f}, {100*t['diff_50']:+.1f}")
    cost = [t for t in table if t["ratio_25"] > 1]
    helped = [t for t in table if t["ratio_25"] < 1]
    print()
    print("WHICH ITEMS")
    print(f"Recovery RMSE ratio above 1 at 25%: {', '.join(t['name'] for t in cost)}. Largest mean recorded count among them: {max(t['recorded'] for t in cost):.1f}. Zeros in the recount: {100*min(t['zeros'] for t in cost):.0f}% to {100*max(t['zeros'] for t in cost):.0f}%.")
    print(f"Recovery RMSE ratio below 1 at 25%: {', '.join(t['name'] for t in helped)}. Mean recorded count: {min(t['recorded'] for t in helped):.1f} to {max(t['recorded'] for t in helped):.1f}. Zeros in the recount: {100*min(t['zeros'] for t in helped):.0f}% to {100*max(t['zeros'] for t in helped):.0f}%.")
    print(f"At 25%, coverage difference: {100*min(t['diff_25'] for t in cost):+.1f} to {100*max(t['diff_25'] for t in cost):+.1f} on the first group, {100*min(t['diff_25'] for t in helped):+.1f} to {100*max(t['diff_25'] for t in helped):+.1f} on the second.")
    print(f"At 50%, recovery RMSE ratio: {min(t['ratio_50'] for t in cost):.3f} to {max(t['ratio_50'] for t in cost):.3f} on the first group, {min(t['ratio_50'] for t in helped):.3f} to {max(t['ratio_50'] for t in helped):.3f} on the second.")
    print(f"At 50%, coverage difference: {100*min(t['diff_50'] for t in cost):+.1f} to {100*max(t['diff_50'] for t in cost):+.1f} on the first group, {100*min(t['diff_50'] for t in helped):+.1f} to {100*max(t['diff_50'] for t in helped):+.1f} on the second.")
    print(f"Correlation of absolute residual with the covariate: {min(t['spread'] for t in cost):.2f} to {max(t['spread'] for t in cost):.2f} on the first group, {min(t['spread'] for t in helped):.2f} to {max(t['spread'] for t in helped):.2f} on the second.")
    print()
    print("THE OPTION AGAINST DOUBLING THE AUDIT")
    doubled = [t["rmse_50"]/t["rmse_25"]-1 for t in table]
    option = [t["ratio_25"]-1 for t in table]
    print(f"Doubling the audit from 25% to 50%, existing option, change in recovery RMSE: {100*min(doubled):+.0f}% to {100*max(doubled):+.0f}% across all nine.")
    print(f"The clustered option at 25%, change in recovery RMSE: {100*min(option):+.0f}% to {100*max(option):+.0f}% across all nine.")
    print("Coverage gain in points on the second group, the clustered option at 25% against doubling the audit with the existing option:")
    for t in helped:
        print(f"  {t['name']}: {100*t['diff_25']:+.1f} against {100*(t['existing_50']-t['existing']):+.1f}")
    print(f"The clustered option gained more on {sum(t['diff_25'] > t['existing_50']-t['existing'] for t in helped)} of {len(helped)}.")

    study = {(r["analysis"], float(r["fraction"]), r["method"], r["replicate"]): r for r in rows("validation/clustered-study/replicates.csv")}
    good = lambda r: r is not None and r["status"] == "ok"
    print("Paired coverage difference, the clustered option at 25% minus the existing option at 50%, with MCSE; change in mean width from the option, and from doubling the audit:")
    for t in helped:
        pairs = [(study.get((t["name"], .25, "linear_mi_clustered", str(i))), study.get((t["name"], .5, "linear_mi", str(i)))) for i in range(1, 501)]
        d = [float(a["coverage"])-float(b["coverage"]) for a, b in pairs if good(a) and good(b)]
        width = lambda fraction, method: float(summary[(t["name"], fraction, method)]["mean_width"])
        print(f"  {t['name']}: {100*st.mean(d):+.1f} ({100*st.stdev(d)/len(d)**.5:.1f}); {100*(width(.25, 'linear_mi_clustered')/width(.25, 'linear_mi')-1):+.0f}%, {100*(width(.5, 'linear_mi')/width(.25, 'linear_mi')-1):+.0f}%")

    print()
    print("THE EXISTING OPTION AGAINST THE SAME CONDITIONS")
    full = {(r["analysis"], float(r["fraction"]), r["method"]): float(r["coverage"]) for r in rows("validation/clustered-study/summary.csv")}
    held = [t for t in table if t["role"] == "held_out"]
    print("Existing minus oracle coverage on the held-out items at 25%, points: " + ", ".join(f"{t['name']} {100*(t['existing']-t['oracle']):+.1f}" for t in held) + f". Within 5 points: {sum(abs(t['existing']-t['oracle']) <= .05 for t in held)} of {len(held)}.")
    cells = sorted({(a, f) for a, f, _ in full if a.startswith(("g70", "g150"))})
    meets = lambda method: sum(full[(a, f, method)] >= .92 and abs(full[(a, f, method)]-full[(a, f, "oracle")]) <= .03 for a, f in cells)
    print(f"Simulation condition, cells with 70 or more clusters: existing {meets('linear_mi')} of {len(cells)}, clustered {meets('linear_mi_clustered')} of {len(cells)}.")
    oracle_30 = [v for (a, _, m), v in full.items() if a.startswith("g30") and m == "oracle"]
    print(f"Oracle coverage, the cluster-robust variance on complete data, with 30 clusters: mean {100*st.mean(oracle_30):.1f}%, lowest {100*min(oracle_30):.1f}%, {sum(v < .92 for v in oracle_30)} of {len(oracle_30)} cells under 92%.")

    print()
    print("SIMULATIONS")
    simulated = [r for r in rows("validation/clustered-study/ratios.csv") if r["level"] == "simulation" and not r["analysis"].startswith("control")]
    above = [r for r in simulated if float(r["recovery_lower"]) > 1]
    below = [r for r in simulated if float(r["recovery_upper"]) < 1]
    print(f"Of {len(simulated)} cells, recovery RMSE ratio clearly above 1 in {len(above)}: " + ", ".join(f"{r['analysis']} at {100*float(r['fraction']):.0f}% {float(r['recovery_ratio']):.3f}" for r in above) + ".")
    print(f"Clearly below 1 in {len(below)}, of which {sum('tracking' in r['analysis'] for r in below)} have tracking errors.")
    sim_icc = rows("validation/clustered-review/simulation-icc.csv")
    for errors in ("iid", "persistent", "tracking"):
        v = [float(r["residual_icc"]) for r in sim_icc if r["errors"] == errors]
        print(f"Within-cluster correlation of the imputation model's residuals, {errors} errors: {min(v):.2f} to {max(v):.2f}.")

    print()
    print("VERY WIDE INTERVALS, 30 clusters at 25%")
    for method in ("linear_mi", "linear_mi_clustered"):
        wide, usual = [], []
        for a in {a for a, _, _ in full if a.startswith("g30")}:
            cell = [r for i in range(1, 501) if good(r := study.get((a, .25, method, str(i))))]
            median = st.median(float(r["width"]) for r in cell)
            for r in cell:
                (wide if float(r["width"]) > 5*median else usual).append(float(r["recovery_error"]))
        print(f"{method}: {len(wide)} requests over 5 times the cell median width, recovery RMSE {(st.mean(e*e for e in wide))**.5:.3f}; others {(st.mean(e*e for e in usual))**.5:.3f}.")
    low_df = {m: [float(r["df"]) < 5 for i in range(1, 501) if good(r := study.get(("iron", .25, m, str(i))))] for m in ("linear_mi", "linear_mi_clustered")}
    print(f"Iron at 25%, requests under 5 degrees of freedom: existing {100*st.mean(low_df['linear_mi']):.1f}%, clustered {100*st.mean(low_df['linear_mi_clustered']):.1f}%.")

    print()
    print("WIDTH RATIO WITHIN AN ITEM, 25%: coverage in the third of audits where the clustered interval is widest relative to the existing one, against the other two thirds")
    worse = 0
    for t in table:
        pairs = [(study.get((t["name"], .25, "linear_mi", str(i))), study.get((t["name"], .25, "linear_mi_clustered", str(i)))) for i in range(1, 501)]
        pairs = sorted(((float(b["width"])/float(a["width"]), float(a["coverage"]), float(b["coverage"])) for a, b in pairs if good(a) and good(b)))
        k = len(pairs)//3
        top, rest = pairs[-k:], pairs[:-k]
        mean = lambda part, j: 100*st.mean(p[j] for p in part)
        worse += mean(top, 1) < mean(rest, 1)
        print(f"  {t['name']}: existing {mean(top, 1):.1f} against {mean(rest, 1):.1f}; clustered {mean(top, 2):.1f} against {mean(rest, 2):.1f}")
    print(f"The existing option covered less in that third on {worse} of {len(table)} items.")

    print()
    print("WITHIN-FACILITY CORRELATION OF THE IMPUTATION MODEL'S RESIDUALS (recount on recorded value, covariate and catchment)")
    pool_icc = {r["analysis"]: r for r in rows("validation/clustered-review/pool-icc.csv")}
    print("Whole pool, with the same quantity leaving out the three facilities with the largest residuals:")
    for t in table:
        print(f"  {t['name']}: {float(pool_icc[t['name']]['icc']):+.2f}, {float(pool_icc[t['name']]['icc_without_three_facilities']):+.2f}")
    print(f"First group: {min(float(pool_icc[t['name']]['icc']) for t in cost):+.2f} to {max(float(pool_icc[t['name']]['icc']) for t in cost):+.2f}. Second group: {min(float(pool_icc[t['name']]['icc']) for t in helped):+.2f} to {max(float(pool_icc[t['name']]['icc']) for t in helped):+.2f}.")
    variants = rows("validation/clustered-review/variants.csv")
    audited = {(r["analysis"], r["replicate"]): r["audited_icc"] for r in variants if r["variant"] == "control" and float(r["fraction"]) == .25}
    print("Same quantity on each audit's own rows at 25%. Coverage gain in points, and recovery RMSE ratio, where it is at least .15 against below .15:")
    pooled = []
    for t in table:
        hi, lo = [], []
        for i in range(1, 501):
            a, b, icc = study.get((t["name"], .25, "linear_mi", str(i))), study.get((t["name"], .25, "linear_mi_clustered", str(i))), audited.get((t["name"], str(i)), "NA")
            if good(a) and good(b) and icc not in ("NA", ""):
                row = (float(b["coverage"])-float(a["coverage"]), float(a["recovery_error"]), float(b["recovery_error"]), float(icc))
                (hi if row[3] >= .15 else lo).append(row)
                if t in helped:
                    pooled.append(row)
        ratio = lambda part: (sum(r[2]**2 for r in part)/sum(r[1]**2 for r in part))**.5
        show = lambda part: f"{len(part)} audits {100*st.mean(r[0] for r in part):+.1f}, {ratio(part):.2f}" if part else "no audits"
        print(f"  {t['name']}: {show(hi)}; {show(lo)}")
    print("The six items of the second group pooled, by cutoff: " + "; ".join(f"{cut:.2f}: {100*st.mean(r[0] for r in pooled if r[3] >= cut):+.1f} against {100*st.mean(r[0] for r in pooled if r[3] < cut):+.1f}" for cut in (0, .1, .15, .2, .3)) + ".")

    print()
    print("RESAMPLING VARIANTS, unplanned, on the study's cohorts, audits and seeds")
    chosen = {(r["analysis"], float(r["fraction"]), r["variant"], r["replicate"]): r for r in variants}
    gaps = [abs(float(r["estimate"])-float(study[(r["analysis"], float(r["fraction"]), "linear_mi_clustered", r["replicate"])]["estimate"])) for r in variants if r["variant"] == "control" and r["status"] == "ok"]
    mismatch = sum((r["status"] == "ok") != good(study[(r["analysis"], float(r["fraction"]), "linear_mi_clustered", r["replicate"])]) for r in variants if r["variant"] == "control")
    print(f"Control: the package's stratified cluster resampling, rerun, against the study's saved records: {len(gaps)} requests, largest difference in estimate {max(gaps):.1e}, {mismatch} status differences.")
    stratum = {}
    for r in pool:
        stratum.setdefault(r["hosp_id"], set()).add(r["facility_id"])
    sizes = sorted(len(v) for v in stratum.values())
    print(f"Facilities per catchment: {sizes}. Audited at 25%: {[max(1, round(.25*n)) for n in sizes]}, {sum(max(1, round(.25*n)) for n in sizes)} in all.")
    print("At each audit size: coverage for existing, stratified clusters (the study), rows within strata, clusters without strata, rescaled (not a package feature), oracle; recovery RMSE ratio against existing for the last four methods; requests refused")
    for fraction in (.25, .5):
        print(f"  {100*fraction:.0f}%")
        for t in table:
            coverage, ratio, refused = [], [], []
            for variant in ("control", "rows", "no_strata", "rescaled"):
                pairs = [(study.get((t["name"], fraction, "linear_mi", str(i))), chosen.get((t["name"], fraction, variant, str(i))), study.get((t["name"], fraction, "oracle", str(i)))) for i in range(1, 501)]
                returned = [(a, b, o) for a, b, o in pairs if b["status"] == "ok"]
                target = float(pairs[0][2]["truth"])
                coverage.append(100*st.mean(float(b["lower"]) <= target <= float(b["upper"]) for _, b, _ in returned))
                both = [(a, b, o) for a, b, o in returned if good(a) and good(o)]
                ratio.append((sum((float(b["estimate"])-float(o["estimate"]))**2 for _, b, o in both)/sum((float(a["estimate"])-float(o["estimate"]))**2 for a, _, o in both))**.5)
                refused.append(500-len(returned))
            print(f"    {t['name']}: {100*float(summary[(t['name'], fraction, 'linear_mi')]['coverage']):.1f}, " + ", ".join(f"{c:.1f}" for c in coverage) + f", {100*float(summary[(t['name'], fraction, 'oracle')]['coverage']):.1f}; " + ", ".join(f"{v:.3f}" for v in ratio) + f"; {refused}")
    row_gain = []
    for t_ in helped:
        pairs = [(study.get((t_["name"], .25, "linear_mi", str(i))), chosen.get((t_["name"], .25, "rows", str(i)))) for i in range(1, 501)]
        target = float(pairs[0][0]["truth"])
        row_gain.append(100*(st.mean(float(b["lower"]) <= target <= float(b["upper"]) for _, b in pairs if b["status"] == "ok")-float(summary[(t_["name"], .25, "linear_mi")]["coverage"])))
    print(f"Rows within strata minus existing, coverage at 25% on the second group: {min(row_gain):+.1f} to {max(row_gain):+.1f}.")
    early = next(t_ for t_ in table if t_["name"] == "anc1_early")
    print(f"anc1_early at 25%, oracle minus the study's stratified clusters: {100*(early['oracle']-early['clustered']):.1f} points.")
    print("Median between-imputation variance, rescaled over the study's stratified clusters, at 25%:")
    for t in table:
        v = [float(chosen[(t["name"], .25, "rescaled", str(i))]["between_variance"])/float(chosen[(t["name"], .25, "control", str(i))]["between_variance"]) for i in range(1, 501)
             if chosen[(t["name"], .25, "rescaled", str(i))]["status"] == "ok" and chosen[(t["name"], .25, "control", str(i))]["status"] == "ok" and float(chosen[(t["name"], .25, "control", str(i))]["between_variance"]) > 0]
        print(f"  {t['name']}: {st.median(v):.2f}")


if __name__ == "__main__":
    main()
