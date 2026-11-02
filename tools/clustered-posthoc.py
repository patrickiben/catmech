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


if __name__ == "__main__":
    main()
