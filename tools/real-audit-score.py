#!/usr/bin/env python3
"""Score the decision rules and predictions in validation/REAL_AUDIT_PLAN.md from saved results, in plain Python.
Usage: python3 -I tools/real-audit-score.py [OUTPUT_FILE] [STUDY_DIRECTORY]
Without a study directory it scores the frozen validation/real-audit-study. The post-hoc sections always read the frozen post-hoc tables."""
import csv
import importlib.util
import math
from pathlib import Path
import random
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT/"validation/real-audit-study"
spec = importlib.util.spec_from_file_location("verify_artifacts", ROOT/"tools/verify-artifacts.py")
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)
OUT = []


def say(text=""):
    OUT.append(text)


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


summary = {(r["level"], r["analysis"], float(r["fraction"]), r["method"]): r for r in rows(STUDY/"summary.csv")}
ratios = {(r["level"], r["analysis"], float(r["fraction"]), r["method"], r["comparator"]): r for r in rows(STUDY/"recovery-ratios.csv")}
pools = {(r["level"], r["analysis"]): r for r in rows(STUDY/"pool-targets.csv")}
scores = {}
for r in rows(STUDY/"corruption-scores.csv"):
    scores.setdefault((r["level"], r["analysis"], .25 if r["scenario"].endswith("25pct") else .5), []).append(r)
posthoc_path = ROOT/"validation/real-audit-posthoc/swap-in-summary.csv"
posthoc = {(r["level"], r["analysis"], float(r["fraction"])): r for r in rows(posthoc_path)} if posthoc_path.exists() else {}
MONTH = ["pnc", "syphilis", "iron", "anc4"]
EXACT = ["pnc", "syphilis", "iron"]
METHODS = ["recorded", "audited_only", "linear_mi", "raking", "contamination"]


def s(level, analysis, fraction, method, field):
    return float(summary[(level, analysis, fraction, method)][field])


def ratio(level, analysis, fraction, method, comparator):
    r = ratios[(level, analysis, fraction, method, comparator)]
    return float(r["recovery_rmse_ratio"]), float(r["lower"]), float(r["upper"])


say("Scored from validation/real-audit-study by tools/real-audit-score.py")
say(f"Requests: {sum(int(r['requested']) for r in summary.values())}; returned: {sum(int(r['returned_intervals']) for r in summary.values())}; "
    f"convergence refusals: {sum(int(r['convergence_failures']) for r in summary.values())}; refused: {sum(int(r['refused']) for r in summary.values())}; other errors: {sum(int(r['other_errors']) for r in summary.values())}")
say()
say("RULE 1. Does a partial audit pay? linear_mi against recorded, facility-months, recovery RMSE ratio")
rule1 = {}
for analysis in MONTH:
    for fraction in (.25, .5):
        value, low, high = ratio("rwanda_month", analysis, fraction, "linear_mi", "recorded")
        rule1[(analysis, fraction)] = high < 1
        say(f"  {analysis:9s} {int(fraction*100)}%: {value:.3f} [{low:.3f}, {high:.3f}] -> {'pays' if high < 1 else ('costs' if low > 1 else 'tie')}")
say()
say("RULE 2. Does the contamination model add anything? contamination against linear_mi, facility-months, 25%")
passes = []
contradicted = 0
for analysis in EXACT:
    value, low, high = ratio("rwanda_month", analysis, .25, "contamination", "linear_mi")
    returned = s("rwanda_month", analysis, .25, "contamination", "interval_return_rate")
    coverage = s("rwanda_month", analysis, .25, "contamination", "coverage_cluster")
    oracle = s("rwanda_month", analysis, .25, "oracle", "coverage_cluster")
    quarter, qlow, qhigh = ratio("rwanda_quarter", analysis, .25, "contamination", "linear_mi")
    lower_error = high < 1
    conditions = returned >= .90 and coverage >= oracle-.03 and quarter < 1
    passes.append(lower_error and conditions)
    contradicted += low > 1
    say(f"  {analysis:9s} ratio {value:.3f} [{low:.3f}, {high:.3f}]; returned {returned:.3f} (need 0.90); clustered coverage {coverage:.3f} vs oracle {oracle:.3f} (need within 0.03); quarter ratio {quarter:.3f} -> lower error {lower_error}, conditions met {conditions}")
outcome = "SUPPORTED" if sum(passes) >= 2 else "CONTRADICTED" if contradicted >= 2 else "INCONCLUSIVE"
say(f"  Outcome: {outcome}")
value, low, high = ratio("rwanda_month", "anc4", .25, "contamination", "linear_mi")
say(f"  anc4, stress case outside the rule: {value:.3f} [{low:.3f}, {high:.3f}]")
value, low, high = ratio("rwanda_month", "anc4", .5, "contamination", "linear_mi")
say(f"  anc4 at 50%: {value:.3f} [{low:.3f}, {high:.3f}]")
say()
say("RULE 3. Control. Recovery RMSE ratio against recorded, effect sizes")
for analysis in ["effect_smd", "effect_r"]:
    for fraction in (.25, .5):
        parts = []
        for method in ["linear_mi", "raking", "audited_only", "contamination"]:
            v, l, h = ratio("effect_size", analysis, fraction, method, "recorded")
            parts.append(f"{method} {v:.3f} [{l:.3f}, {h:.3f}]")
        say(f"  {analysis:10s} {int(fraction*100)}%: "+"; ".join(parts))
say()
say("PREDICTIONS")
anc4 = ratio("rwanda_month", "anc4", .25, "linear_mi", "recorded")
pnc = ratio("rwanda_month", "pnc", .25, "linear_mi", "recorded")
say(f"  1. anc4 linear_mi beats recorded: {anc4[2] < 1}; pnc does not: {not pnc[2] < 1}")
parts = []
for analysis in EXACT:
    v, _, _ = ratio("rwanda_month", analysis, .25, "contamination", "linear_mi")
    parts.append(f"{analysis} lower: {v < 1}")
v, _, _ = ratio("rwanda_month", "anc4", .25, "contamination", "linear_mi")
say(f"  2. {'; '.join(parts)}; anc4 not lower: {not v < 1}")
undercover = all(s("rwanda_month", a, f, m, "coverage") < .95 for a in MONTH for f in (.25, .5) for m in ["oracle"]+METHODS)
closer = all(abs(s("rwanda_month", a, f, m, "coverage_cluster")-.95) < abs(s("rwanda_month", a, f, m, "coverage")-.95) for a in MONTH for f in (.25, .5) for m in ["oracle"]+METHODS)
within = {a: s("rwanda_month", a, .25, "oracle", "coverage_cluster")-s("rwanda_month", a, .25, "contamination", "coverage_cluster") for a in MONTH}
say(f"  3. default undercovers for every method and fraction: {undercover}; clustered closer to 95% for every method and fraction: {closer}; contamination gap to oracle in points: " + ", ".join(f"{a} {100*g:.1f}" for a, g in within.items()))
for analysis in ["effect_smd", "effect_r"]:
    values = {m: s("effect_size", analysis, .25, m, "recovery_rmse") for m in METHODS}
    say(f"  4. {analysis}: lowest is {min(values, key=values.get)}: " + ", ".join(f"{m} {v:.4f}" for m, v in values.items()))
for level, names in [("rwanda_month", MONTH), ("rwanda_quarter", MONTH), ("effect_size", ["effect_smd", "effect_r"])]:
    default = {a: s(level, a, .25, "raking", "mean_width") > s(level, a, .25, "linear_mi", "mean_width") for a in names}
    clustered = {a: s(level, a, .25, "raking", "mean_width_cluster") > s(level, a, .25, "linear_mi", "mean_width_cluster") for a in names}
    say(f"  5. {level}: raking wider, default {sum(default.values())}/{len(names)}, clustered {sum(clustered.values())}/{len(names)}")
auc = {}
for (level, analysis, fraction), group in scores.items():
    values = [float(r["auc"]) for r in group if r["auc"] not in ("", "NA")]
    auc[(level, analysis, fraction)] = sum(values)/len(values)
say("  6. mean AUC at facility-months: " + ", ".join(f"{a} {int(f*100)}% {auc[('rwanda_month', a, f)]:.3f}" for a in MONTH for f in (.25, .5)))
failures = []
for (level, analysis, fraction, method), r in summary.items():
    if fraction == .25 and method in ["audited_only", "linear_mi", "raking", "contamination"]:
        if not s(level, analysis, .5, method, "recovery_rmse") < float(r["recovery_rmse"]):
            failures.append((level, analysis, method))
say(f"  7. audit-using methods better at 50% than 25% in every cell: {not failures} {failures}")
say()
say("TABLES: recovery RMSE, bias, clustered coverage")
for level, names in [("rwanda_month", MONTH), ("rwanda_quarter", MONTH), ("effect_size", ["effect_smd", "effect_r"])]:
    for analysis in names:
        for fraction in (.25, .5):
            cells = []
            for method in METHODS:
                cells.append(f"{s(level, analysis, fraction, method, 'recovery_rmse'):.3f}, {s(level, analysis, fraction, method, 'bias'):+.3f}, {100*s(level, analysis, fraction, method, 'coverage_cluster'):.1f}%")
            r = summary[(level, analysis, fraction, "contamination")]
            say(f"| {level} | {analysis} | {int(fraction*100)}% | {r['returned_intervals']}/{r['requested']} | " + " | ".join(cells) + f" | {s(level, analysis, fraction, 'oracle', 'bias'):+.3f}, {100*s(level, analysis, fraction, 'oracle', 'coverage_cluster'):.1f}% |")
say()
say("DEFAULT COVERAGE AND WIDTHS, facility-months")
for analysis in MONTH:
    for fraction in (.25, .5):
        say(f"  {analysis:9s} {int(fraction*100)}%: " + "; ".join(f"{m} default {100*s('rwanda_month', analysis, fraction, m, 'coverage'):.1f}% width {s('rwanda_month', analysis, fraction, m, 'mean_width'):.3f} clustered width {s('rwanda_month', analysis, fraction, m, 'mean_width_cluster'):.3f}" for m in ["oracle", "linear_mi", "contamination"]))
say()
say("COVERING INTERVALS PER REQUEST, clustered, facility-months")
for analysis in MONTH:
    for fraction in (.25, .5):
        say(f"  {analysis:9s} {int(fraction*100)}%: " + "; ".join(f"{m} {100*s('rwanda_month', analysis, fraction, m, 'coverage_cluster')*s('rwanda_month', analysis, fraction, m, 'interval_return_rate'):.1f}%" for m in ["linear_mi", "raking", "contamination"]))
say()
say("PER-RECORD SCORES, mean across cohorts: AUC; Brier vs prevalence; actual vs expected errors; posterior vs recorded RMSE")
for (level, analysis, fraction), group in scores.items():
    mean = lambda field: sum(float(r[field]) for r in group)/len(group)
    say(f"  {level:15s}{analysis:11s}{int(fraction*100)}%: AUC {auc[(level, analysis, fraction)]:.3f}; Brier {mean('brier'):.3f} vs {mean('brier_prevalence'):.3f}; errors {mean('actual_errors'):.1f} vs {mean('expected_errors'):.1f}; RMSE {mean('completion_rmse'):.3f} vs {mean('recorded_rmse'):.3f}")
say()
say("RMSE AGAINST THE POOL SLOPE, as the plan listed: method values, then paired ratios with 95% limits from 2,000 resamples")
records = {}
for r in rows(STUDY/"replicates.csv"):
    if r["status"] == "ok":
        records.setdefault((r["level"], r["analysis"], float(r["fraction"]), r["method"]), {})[r["replicate"]] = float(r["error"])


def pool_ratio(level, analysis, fraction, first, second, seed):
    a, b = records[(level, analysis, fraction, first)], records[(level, analysis, fraction, second)]
    shared = sorted(set(a) & set(b))
    x, y = [a[k]**2 for k in shared], [b[k]**2 for k in shared]
    generator = random.Random(seed)
    draws = []
    for _ in range(2000):
        picks = generator.choices(range(len(x)), k=len(x))
        draws.append(math.sqrt(sum(x[i] for i in picks)/sum(y[i] for i in picks)))
    draws.sort()
    return math.sqrt(sum(x)/sum(y)), draws[49], draws[1949]


seed = 20261102
for level, names in [("rwanda_month", MONTH), ("effect_size", ["effect_smd", "effect_r"])]:
    for analysis in names:
        for fraction in (.25, .5):
            values = ", ".join(f"{m} {s(level, analysis, fraction, m, 'rmse'):.3f}" for m in ["oracle"]+METHODS)
            seed += 1
            first = pool_ratio(level, analysis, fraction, "linear_mi", "recorded", seed)
            second = pool_ratio(level, analysis, fraction, "contamination", "linear_mi", seed+1000)
            say(f"  {analysis:10s} {int(fraction*100)}%: {values}; linear_mi/recorded {first[0]:.3f} [{first[1]:.3f}, {first[2]:.3f}]; contamination/linear_mi {second[0]:.3f} [{second[1]:.3f}, {second[2]:.3f}]")
say()
say("COUNTS ACROSS THE 20 CELLS (non-oracle planned methods)")
cells = [(level, analysis, fraction) for level, names in [("rwanda_month", MONTH), ("rwanda_quarter", MONTH), ("effect_size", ["effect_smd", "effect_r"])] for analysis in names for fraction in (.25, .5)]
for field, label in [("recovery_rmse", "lowest recovery RMSE"), ("rmse", "lowest RMSE against the pool")]:
    tally = {m: 0 for m in METHODS}
    for cell in cells:
        tally[min(METHODS, key=lambda m: s(*cell, m, field))] += 1
    say(f"  {label}: " + ", ".join(f"{m} {n}" for m, n in tally.items()))
tally = {m: 0 for m in METHODS}
for cell in cells:
    tally[min(METHODS, key=lambda m: abs(s(*cell, m, "coverage_cluster")-.95))] += 1
say("  clustered coverage closest to 95%: " + ", ".join(f"{m} {n}" for m, n in tally.items()))
say("  worst clustered coverage, four facility-month analyses at 25%: " + ", ".join(f"{m} {100*min(s('rwanda_month', a, .25, m, 'coverage_cluster') for a in MONTH):.1f}%" for m in METHODS))
say()
say("REFUSALS, facility-months at 25%: other methods on refused against returned cohorts, and falling back to linear MI when the model refuses")
everything = {}
for r in rows(STUDY/"replicates.csv"):
    if r["level"] == "rwanda_month" and r["fraction"] == "0.25":
        everything.setdefault((r["analysis"], r["method"]), {})[r["replicate"]] = r
for analysis in MONTH:
    model = everything[(analysis, "contamination")]
    returned = {k for k, v in model.items() if v["status"] == "ok"}
    refused = set(model)-returned
    parts = []
    for m in ["recorded", "audited_only", "linear_mi"]:
        group = everything[(analysis, m)]
        rmse_of = lambda keys: math.sqrt(sum(float(group[k]["recovery_error"])**2 for k in keys if group[k]["status"] == "ok")/sum(group[k]["status"] == "ok" for k in keys))
        parts.append(f"{m} {rmse_of(refused)/rmse_of(returned):.3f}")
    fallback = {k: (model[k] if k in returned else everything[(analysis, "linear_mi")][k]) for k in model if k in returned or everything[(analysis, "linear_mi")][k]["status"] == "ok"}
    recovery = math.sqrt(sum(float(v["recovery_error"])**2 for v in fallback.values())/len(fallback))
    coverage = sum(float(v["coverage_cluster"]) for v in fallback.values())/len(fallback)
    say(f"  {analysis:9s} refused {len(refused)}; refused/returned recovery RMSE: " + ", ".join(parts) + f"; with fallback recovery RMSE {recovery:.3f}, clustered coverage {100*coverage:.1f}%")
say()
say("SHARE OF REPORTED BIAS KEPT, facility-months: by pool bias, then by recovery bias")
for analysis in ["pnc", "iron", "anc4"]:
    for fraction in (.25, .5):
        by_pool = s("rwanda_month", analysis, fraction, "contamination", "bias")/s("rwanda_month", analysis, fraction, "recorded", "bias")
        by_recovery = s("rwanda_month", analysis, fraction, "contamination", "recovery_bias")/s("rwanda_month", analysis, fraction, "recorded", "recovery_bias")
        lin_pool = s("rwanda_month", analysis, fraction, "linear_mi", "bias")/s("rwanda_month", analysis, fraction, "recorded", "bias")
        say(f"  {analysis:9s} {int(fraction*100)}%: contamination {100*by_pool:.0f}% and {100*by_recovery:.0f}%; linear_mi {100*lin_pool:.0f}%")
say()
say("POOL DESCRIPTION")
frames = verify.real_audit_pools()
for key, frame in frames.items():
    r = pools[key]
    errors = [row[1]-row[2] for row in frame]
    xs = [row[3] for row in frame]
    correlation = statistics.correlation(errors, xs)
    say(f"  {key[0]:15s}{key[1]:11s} rows {r['rows']:>4s} units {r['units']:>4s} target {float(r['target_coefficient']):+.4f} reported {float(r['recorded_coefficient']):+.4f} gap {float(r['recorded_minus_target']):+.4f} wrong {r['wrong_rows']:>4s} above {r['recorded_above_reference']:>4s} below {r['recorded_below_reference']:>4s} error/resid {float(r['error_to_residual']):.2f} corr(error, covariate) {correlation:+.3f}")
say()
say("INFLUENCE: slope using reference values when each facility is left out, facility-months")
centers = [r for r in rows(ROOT/"validation/real-audit/rwanda-pool.csv") if r["facility_type"] == "health_center"]
for name, (w, reference, x) in verify.RWANDA_COLUMNS.items():
    kept = [r for r in centers if r[w] != "" and r[reference] != "" and r[x] != ""]
    facilities = sorted({r["facility_id"] for r in kept})
    full = verify.real_audit_slope([(r["hosp_id"], float(r[w]), float(r[reference]), float(r[x])) for r in kept], 2)
    values = []
    for facility in facilities:
        frame = [(r["hosp_id"], float(r[w]), float(r[reference]), float(r[x])) for r in kept if r["facility_id"] != facility]
        values.append(verify.real_audit_slope(frame, 2))
    say(f"  {name:9s} full {full:.3f}; leave-one-out range {min(values):.3f} to {max(values):.3f}")
if posthoc:
    say()
    say("POST HOC, not in the plan: swap-in uses reference values on audited rows and recorded values elsewhere")
    for level, names in [("rwanda_month", MONTH), ("rwanda_quarter", MONTH), ("effect_size", ["effect_smd", "effect_r"])]:
        for analysis in names:
            for fraction in (.25, .5):
                r = posthoc[(level, analysis, fraction)]
                recorded_bias, cont, swap = float(r["recorded_bias"]), float(r["contamination_bias"]), float(r["swap_in_bias_same_cohorts"])
                share = f"contamination keeps {100*cont/recorded_bias:.0f}% of reported bias, swap-in {100*swap/recorded_bias:.0f}%" if abs(recorded_bias) > .02 else "reported bias near zero"
                say(f"  {level:15s}{analysis:11s}{int(fraction*100)}%: recovery RMSE recorded {float(r['recorded_recovery_rmse']):.4f}, swap-in {float(r['swap_in_recovery_rmse']):.4f}, linear_mi {float(r['linear_mi_recovery_rmse']):.4f}; same cohorts contamination {float(r['contamination_recovery_rmse']):.4f} vs swap-in {float(r['swap_in_recovery_rmse_same_cohorts']):.4f}; cor {float(r['correlation_contamination_swap_in']):.3f}, gap {float(r['rms_gap_contamination_swap_in']):.4f}; {share}")
review = ROOT/"validation/real-audit-review"
if review.exists():
    say()
    say("POST HOC, AFTER INDEPENDENT REVIEW, not in the plan")
    say("  Slopes: reported; one-item target; both sides recounted; reported minus target; reported minus both-sides recount")
    for r in rows(review/"estimand-slopes.csv"):
        say(f"  {r['analysis']:9s} {float(r['reported']):.3f}; {float(r['one_item_target']):.3f}; {float(r['both_sides_recount']):.3f}; {float(r['gap_one_item']):+.3f} [{float(r['gap_one_item_lower']):.3f}, {float(r['gap_one_item_upper']):.3f}] ({float(r['gap_one_item_percent']):+.0f}%); {float(r['gap_both_sides']):+.3f} [{float(r['gap_both_sides_lower']):.3f}, {float(r['gap_both_sides_upper']):.3f}] ({float(r['gap_both_sides_percent']):+.0f}%); covariate wrong on {r['covariate_rows_wrong']} of {r['rows']}")
    say("  Gap structure: reported equals covariate; wrong rows copying it; gap and gap without three facilities; within-catchment correlation; recounts of zero and slopes without them")
    for r in rows(review/"gap-structure.csv"):
        say(f"  {r['analysis']:9s} {r['reported_equals_covariate']} of {r['rows']}; {r['wrong_rows_copying_covariate']} of {r['wrong_rows']}; {float(r['gap']):.3f} to {float(r['gap_without_top3']):.3f} ({r['top3']}); {float(r['within_catchment_correlation']):.3f}; {r['zero_recount_rows']} rows in {r['zero_recount_facilities']} facilities, target {float(r['one_item_without_zero']):.3f}, reported {float(r['reported_without_zero']):.3f}, gap {float(r['gap_without_zero']):.3f}")
    say("  Mean error when no key facility is audited, then when any is")
    for r in rows(review/"conditional-on-key-facilities.csv"):
        say(f"  {r['analysis']:9s} {r['method']:14s} {float(r['mean_error_none_audited']):+.3f}, {float(r['mean_error_any_audited']):+.3f}")
    say("  Rule 1 on fresh cohorts: recovery ratio; ratio against the pool slope")
    for r in rows(review/"rule1-extension.csv"):
        say(f"  {r['analysis']:30s} {float(r['fraction']):.2f} cohorts {r['cohorts']}: {float(r['recovery_ratio']):.3f} [{float(r['recovery_lower']):.3f}, {float(r['recovery_upper']):.3f}]; {float(r['pool_ratio']):.3f} [{float(r['pool_lower']):.3f}, {float(r['pool_upper']):.3f}]")
    say("  Checking errors against the covariate on audited rows, 25%: share significant; share with the wrong sign")
    for r in rows(review/"audited-rows-detection.csv"):
        say(f"  {r['analysis']:9s} {100*float(r['share_significant']):.1f}%; {100*float(r['share_wrong_sign']):.1f}%")
    say("  Oracle coverage: default; allowing unequal variance; clustered")
    for r in rows(review/"oracle-interval-types.csv"):
        say(f"  {r['analysis']:9s} {100*float(r['default']):.1f}%; {100*float(r['hc1_unequal_variance']):.1f}%; {100*float(r['cr1_clustered']):.1f}%")
    say("  Errors shuffled within catchment, 25%: oracle coverage; model returned, bias, coverage, SE over SD; linear MI coverage")
    for r in rows(review/"mechanism-shuffled-errors.csv"):
        say(f"  {r['analysis']:9s} {100*float(r['oracle_coverage']):.1f}%; {r['model_returned']} of {r['cohorts']}, {float(r['model_bias']):+.3f}, {100*float(r['model_coverage']):.1f}%, {float(r['model_se_over_sd']):.3f}; {100*float(r['linear_coverage']):.1f}%")
    say("  Error means by facility-size group, paired cohorts: strata then default then linear MI, for recovery RMSE, bias and coverage")
    for r in rows(review/"error-strata-by-size.csv"):
        say(f"  {r['analysis']:9s} paired {r['paired']}: {float(r['strata_recovery_rmse']):.3f}, {float(r['default_recovery_rmse']):.3f}, {float(r['linear_recovery_rmse']):.3f}; {float(r['strata_bias']):+.3f}, {float(r['default_bias']):+.3f}, {float(r['linear_bias']):+.3f}; {100*float(r['strata_coverage']):.1f}%, {100*float(r['default_coverage']):.1f}%, {100*float(r['linear_coverage']):.1f}%")
    say("  Imputation resampling audited facilities, first 200 cohorts at 25%: coverage scratch, package, oracle; recovery RMSE scratch, package")
    for r in rows(review/"linear-mi-facility-bootstrap.csv"):
        say(f"  {r['analysis']:9s} {100*float(r['facility_bootstrap_coverage']):.1f}%, {100*float(r['package_coverage']):.1f}%, {100*float(r['oracle_coverage']):.1f}%; {float(r['facility_bootstrap_recovery_rmse']):.3f}, {float(r['package_recovery_rmse']):.3f}")
    say("  Raking with the small-sample factor: clustered coverage saved, then with the factor")
    for r in rows(review/"raking-like-with-like.csv"):
        say(f"  {r['scenario']:40s} {100*float(r['coverage_saved']):.1f}%, {100*float(r['coverage_with_factor']):.1f}%")
    difference = rows(review/"difference-estimator.csv")
    say(f"  Difference estimator: below linear MI in {sum(float(r['difference_recovery_rmse']) < float(r['linear_mi_recovery_rmse']) for r in difference)} of {len(difference)} cells; largest ratio to linear MI {max(float(r['difference_recovery_rmse'])/float(r['linear_mi_recovery_rmse']) for r in difference):.3f}; raking above it in {sum(float(r['raking_recovery_rmse']) > float(r['difference_recovery_rmse']) for r in difference)}")
    fixed = rows(review/"fixed-pool-swap-in.csv")
    say(f"  Fixed pool: swap-in beats recorded in {sum(float(r['swap_over_recorded']) < 1 for r in fixed)} of {len(fixed)} cells")
    for r in rows(review/"fixed-pool-methods.csv"):
        say(f"  {r['analysis']:9s} fixed pool, 25%: linear MI over recorded {float(r['linear_over_recorded']):.3f}; contamination over linear MI {float(r['model_over_linear']):.3f}")
text = "\n".join(OUT)+"\n"
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
