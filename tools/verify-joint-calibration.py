#!/usr/bin/env python3
"""Verify the staged joint-calibration study independently of its R summaries.

Scientific calculations use Python's standard library. Only bootstrap index
generation delegates to base R, to replay the protocol's specified R RNG stream;
bootstrap coverage distances and percentiles are recomputed here. No sampler or
study-runner code is loaded. This checks reporting, not inferential validity.
"""
import argparse
import array
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tempfile


CELLS = ["gaussian_balanced", "gaussian_unbalanced", "near_zero_icc",
         "observed_mar", "null_slope", "student_t5", "contextual_x"]
PRIMARY = ["gaussian_balanced", "gaussian_unbalanced", "observed_mar", "null_slope"]
PRIOR_CELLS = ["gaussian_balanced", "near_zero_icc", "observed_mar"]
MAIN = ["old_cv3", "new_cv3", "direct", "oracle_cv3"]
PRIOR_METHODS = ["new_cv3_prior_quarter", "direct_prior_quarter",
                 "new_cv3_prior_four", "direct_prior_four"]
OLD_CONTROL = dict(chains=4, burnin=1000, iterations=2000, max_iterations=8000, prior_scale=1)
NEW_CONTROL = dict(chains=4, burnin=1000, iterations=4000, max_iterations=32000, prior_scale=1)
MISSING = {"", "NA", "NaN", "null"}
Z975 = statistics.NormalDist().inv_cdf(.975)


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row, field):
    value = row[field]
    return math.nan if value in MISSING else float(value)


def close(actual, expected, label):
    if math.isnan(expected):
        require(math.isnan(actual), f"{label}: expected absent, got {actual}")
    else:
        require(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=2e-10, abs_tol=2e-12),
                f"{label}: {actual} != {expected}")


def mean(values):
    return statistics.mean(values) if values else math.nan


def sd(values):
    return statistics.stdev(values) if len(values) > 1 else math.nan


def mcse(values):
    return sd(values) / math.sqrt(len(values)) if values else math.nan


def quantile(values, probability):
    values = sorted(values)
    if not values:
        return math.nan
    index = (len(values) - 1) * probability
    low = math.floor(index)
    high = math.ceil(index)
    return values[low] + (index - low) * (values[high] - values[low])


def wilson(successes, n):
    if not n:
        return math.nan, math.nan
    p = successes / n
    den = 1 + Z975**2 / n
    center = (p + Z975**2 / (2*n)) / den
    half = Z975 * math.sqrt(p*(1-p)/n + Z975**2/(4*n*n)) / den
    return max(0, center-half), min(1, center+half)


def beta_cf(a, b, x):
    """Continued fraction for the regularized incomplete beta function."""
    tiny = 1e-300
    c = 1.0
    d = 1 - (a+b)*x/(a+1)
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1001):
        aa = m*(b-m)*x/((a+2*m-1)*(a+2*m))
        d = 1+aa*d
        c = 1+aa/c
        d = 1/(d if abs(d) > tiny else tiny)
        c = c if abs(c) > tiny else tiny
        h *= d*c
        aa = -(a+m)*(a+b+m)*x/((a+2*m)*(a+2*m+1))
        d = 1+aa*d
        c = 1+aa/c
        d = 1/(d if abs(d) > tiny else tiny)
        c = c if abs(c) > tiny else tiny
        delta = d*c
        h *= delta
        if abs(delta-1) < 3e-14:
            return h
    raise AssertionError("Incomplete-beta calculation did not converge")


def regularized_beta(x, a, b):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    leading = math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)
                       + a*math.log(x) + b*math.log1p(-x))
    if x < (a+1)/(a+b+2):
        return leading * beta_cf(a, b, x) / a
    return 1 - leading * beta_cf(b, a, 1-x) / b


def t975(df):
    require(math.isfinite(df) and df > 0, "Invalid t degrees of freedom")
    def upper_tail(x):
        return .5 * regularized_beta(df/(df+x*x), df/2, .5)
    low, high = 0.0, 2.0
    while upper_tail(high) > .025:
        high *= 2
        require(math.isfinite(high), "Unbounded t quantile")
    for _ in range(80):
        middle = (low+high)/2
        if upper_tail(middle) > .025:
            low = middle
        else:
            high = middle
    return (low+high)/2


def index_rows(rows, fields, label):
    out = {tuple(r[f] for f in fields): r for r in rows}
    require(len(out) == len(rows), f"Duplicate {label} keys")
    return out


def key(row, kind="method"):
    return row["cell"], int(row["replicate"]), row[kind]


def stage_settings(mode):
    require(mode in ("pilot", "development", "confirmation"), "Unknown stage")
    return {
        "pilot": (1, 1, 1, 1590000000, 1600000000, None),
        "development": (30, 30, 30, 1610000000, 1620000000, 1650000000),
        "confirmation": (400, 200, 200, 1630000000, 1640000000, 1660000000),
    }[mode]


def expected_grid(mode):
    ordinary, contextual, prior, *_ = stage_settings(mode)
    analyses, fits = set(), set()
    for cell in CELLS:
        for rep in range(1, (contextual if cell == "contextual_x" else ordinary)+1):
            methods = MAIN + (PRIOR_METHODS if cell in PRIOR_CELLS and rep <= prior else [])
            analyses.update((cell, rep, method) for method in methods)
            ids = ["old", "new", "oracle"] + (["new_prior_quarter", "new_prior_four"]
                   if cell in PRIOR_CELLS and rep <= prior else [])
            fits.update((cell, rep, name) for name in ids)
    return analyses, fits


def verify_sources(out, root, settings, completed):
    before, after = read_csv(out/"source-hashes-before.csv"), read_csv(out/"source-hashes-after.csv")
    require(before == after, "Sources changed during the stage")
    expected = {"DESCRIPTION", "NAMESPACE", "validation/JOINT_CALIBRATION_PLAN.md",
                "validation/run-joint-calibration.R", "validation/v0.8.0-source/DESCRIPTION",
                "validation/v0.8.0-source/NAMESPACE"}
    expected.update(str(p.relative_to(root)) for p in (root/"R").glob("*.R"))
    expected.update(str(p.relative_to(root)) for p in (root/"validation/v0.8.0-source/R").glob("*.R"))
    require(len(before) == len(expected) and {r["file"] for r in before} == expected,
            "Hashed source set differs from candidate/baseline protocol set")
    for row in before:
        path = root/row["file"]
        require(path.is_file(), f"Missing source: {path}")
        content = path.read_bytes()
        require(len(content) == int(row["bytes"]), f"Source length changed: {row['file']}")
        require(hashlib.sha256(content).hexdigest() == row["sha256"], f"Source hash changed: {row['file']}")
    require(settings["source_fingerprint"] == completed["source_fingerprint"], "Completion source fingerprint differs")


def verify_records(out, root):
    settings = json.loads((out/"settings.json").read_text())
    completed = json.loads((out/"COMPLETE.json").read_text())
    mode = settings["mode"]
    ordinary, contextual, prior, data_base, imp_base, boot_seed = stage_settings(mode)
    require(settings["cells"] == CELLS and settings["primary_cells"] == PRIMARY
            and settings["prior_cells"] == PRIOR_CELLS, "Cell order/roles changed")
    for field, expected in dict(primary_reps=ordinary, contextual_reps=contextual, prior_reps=prior,
        data_seed_base=data_base, imputation_seed_base=imp_base, bootstrap_seed=boot_seed,
        bootstrap_reps=2000, m=40, old_control=OLD_CONTROL, new_control=NEW_CONTROL,
        baseline_source_commit="b236293", pan_version="2.0", accuracy_exported=mode != "pilot").items():
        require(settings[field] == expected, f"Protocol setting changed: {field}")
    require(1 <= settings["workers"] <= 16, "Unsupported worker count")
    require(completed["status"] == "complete" and completed["mode"] == mode, "Stage is not completed")
    verify_sources(out, root, settings, completed)
    records, fits = read_csv(out/"records.csv"), read_csv(out/"fit-records.csv")
    record_map = {key(r): r for r in records}
    fit_map = {key(r, "fit_id"): r for r in fits}
    expected_records, expected_fits = expected_grid(mode)
    require(len(records) == len(record_map) and set(record_map) == expected_records, "Analysis request grid mismatch")
    require(len(fits) == len(fit_map) and set(fit_map) == expected_fits, "MCMC/oracle fit grid mismatch")
    cohorts = {}
    metadata = ("cell_index", "data_seed", "imputation_seed", "n", "sites", "gold_n", "paired_n",
                "missing_recorded_n", "wholly_unaudited_sites", "target", "structural_slope",
                "missing_logit_intercept", "audit_logit_intercept", "input_sha256")
    for row in records+fits:
        cell, rep = row["cell"], int(row["replicate"])
        cell_index = CELLS.index(cell)+1
        require(int(row["cell_index"]) == cell_index, "Cell index mismatch")
        require(int(row["data_seed"]) == data_base+cell_index*100000+rep, "Data seed mismatch")
        require(int(row["imputation_seed"]) == imp_base+cell_index*100000+rep, "Imputation seed mismatch")
        require(int(row["n"]) == 240 and int(row["sites"]) == 24, "Cohort size changed")
        require(6 <= int(row["wholly_unaudited_sites"]) <= 24, "Unaudited-site count invalid")
        require(0 <= int(row["paired_n"]) <= int(row["gold_n"]) <= 240, "Reference counts invalid")
        if cell != "observed_mar":
            require(int(row["gold_n"]) == 60, "Fixed audit size changed")
            require(int(row["missing_recorded_n"]) == (84 if cell == "gaussian_unbalanced" else 72),
                    "Recorded-value masking changed")
        close(number(row, "structural_slope"), 0 if cell == "null_slope" else .7, "Structural slope")
        close(number(row, "target"), .9 if cell == "contextual_x" else 0 if cell == "null_slope" else .7, "Estimand")
        previous = cohorts.setdefault((cell, rep), tuple(row[f] for f in metadata))
        require(previous == tuple(row[f] for f in metadata), f"Unpaired cohort information: {cell}/{rep}")
        require(row["source_fingerprint"] == settings["source_fingerprint"], "Worker source fingerprint mismatch")
    statistical = ("estimate", "spread", "lower", "upper", "width", "df", "within_variance",
                   "between_variance", "error", "interval_score", "posterior_rhat", "posterior_ess_bulk",
                   "posterior_ess_tail", "mcse_mean", "mcse_lower", "mcse_upper", "posterior_draws")
    for k, row in record_map.items():
        method = k[2]
        direct = method.startswith("direct")
        fit_id = "new"+method[6:] if direct else method.replace("_cv3", "", 1)
        parent = fit_map[k[:2]+(fit_id,)]
        require(row["fit_id"] == fit_id and row["analysis"] == ("direct" if direct else "cv3"), f"Wrong shared fit: {k}")
        require(row["spread_type"] == ("posterior_sd" if direct else "pooled_se"), f"Mislabeled uncertainty: {k}")
        scale = .25 if method.endswith("quarter") else 4 if method.endswith("four") else 1
        close(number(row, "prior_scale"), scale, f"Prior scale: {k}")
        require(row["status"] in ("ok", "error"), f"Unknown analysis status: {k}")
        if parent["fit_status"] == "error":
            require(row["status"] == "error" and row["failure_stage"] == "fit", f"Failed fit reused: {k}")
            require(row["error_class"] == parent["error_class"] and row["message"] == parent["message"], f"Fit error lost: {k}")
        elif row["status"] == "error":
            require(row["failure_stage"] == "analysis", f"Analysis failure mislabeled: {k}")
        else:
            require(not row["failure_stage"] and not row["error_class"] and not row["message"], f"Successful error metadata: {k}")
        if row["status"] == "error":
            require(bool(row["error_class"]) and bool(row["message"]), f"Missing refusal reason: {k}")
        if mode == "pilot" or row["status"] == "error":
            for field in statistical:
                require(math.isnan(number(row, field)), f"Pilot/refusal exports accuracy: {k}/{field}")
            close(number(row, "covered"), math.nan if mode == "pilot" else 0, f"Failure coverage: {k}")
            continue
        estimate, spread, lower, upper, target = (number(row, f) for f in ("estimate", "spread", "lower", "upper", "target"))
        require(all(math.isfinite(v) for v in (estimate, spread, lower, upper)) and spread > 0 and upper > lower,
                f"Invalid interval: {k}")
        close(number(row, "error"), estimate-target, f"Raw error: {k}")
        close(number(row, "width"), upper-lower, f"Raw width: {k}")
        close(number(row, "covered"), float(lower <= target <= upper), f"Raw coverage: {k}")
        close(number(row, "interval_score"), upper-lower+40*max(0, lower-target)+40*max(0, target-upper), f"Interval score: {k}")
        if direct:
            require(row["interpretation"] == "conditional_joint_model_credible_interval", f"Direct interpretation: {k}")
            for field in ("df", "within_variance", "between_variance"):
                require(math.isnan(number(row, field)), f"Direct interval has pooling metadata: {k}")
            require(number(row, "posterior_draws") == 4*number(parent, "iterations_used"), f"Not all retained draws used: {k}")
            require(math.isfinite(number(row, "posterior_rhat")) and number(row, "posterior_rhat") <= 1.05
                    and number(row, "posterior_ess_bulk") >= 400 and number(row, "posterior_ess_tail") >= 400,
                    f"Direct interval failed diagnostics: {k}")
            for field in ("mcse_mean", "mcse_lower", "mcse_upper"):
                require(math.isfinite(number(row, field)) and number(row, field) >= 0, f"Invalid MCSE: {k}/{field}")
        else:
            within, between, df = (number(row, f) for f in ("within_variance", "between_variance", "df"))
            require(within >= 0 and between >= 0 and df > 0, f"Invalid pooling quantities: {k}")
            total = within+1.025*between
            close(spread*spread, total, f"Rubin variance: {k}")
            fraction = 1.025*between/total
            tmp = (1-fraction)*24*23
            expected_df = 23 if fit_id == "oracle" else 39*tmp/(26*39+fraction*fraction*tmp)
            close(df, expected_df, f"Pooled df: {k}")
            half = t975(df)*spread
            close(lower, estimate-half, f"t lower endpoint: {k}")
            close(upper, estimate+half, f"t upper endpoint: {k}")
            for field in ("posterior_rhat", "posterior_ess_bulk", "posterior_ess_tail", "mcse_mean", "mcse_lower", "mcse_upper", "posterior_draws"):
                require(math.isnan(number(row, field)), f"CV3 has posterior metadata: {k}/{field}")
    return settings, completed, records, fits, record_map, fit_map


def expected_summary(rows):
    ok = [r for r in rows if r["status"] == "ok"]
    n, nk = len(rows), len(ok)
    covers = sum(number(r, "covered") for r in rows)
    pc, cc = covers/n, covers/nk if nk else math.nan
    failure = (n-nk)/n
    pw, cw, fw = wilson(covers, n), wilson(covers, nk), wilson(n-nk, n)
    errors = [number(r, "error") for r in ok]
    squares = [e*e for e in errors]
    rmse = math.sqrt(mean(squares)) if squares else math.nan
    ans = dict(requests=n, intervals=nk, failures=n-nk,
        coverage_per_request=pc, coverage_per_request_mcse=math.sqrt(pc*(1-pc)/n),
        coverage_per_request_lower=pw[0], coverage_per_request_upper=pw[1],
        coverage_conditional=cc, coverage_conditional_mcse=math.sqrt(cc*(1-cc)/nk) if nk else math.nan,
        coverage_conditional_lower=cw[0], coverage_conditional_upper=cw[1], failure_rate=failure,
        failure_rate_mcse=math.sqrt(failure*(1-failure)/n), failure_rate_lower=fw[0], failure_rate_upper=fw[1],
        bias=mean(errors), bias_mcse=mcse(errors), rmse=rmse,
        rmse_mcse=mcse(squares)/(2*rmse) if rmse > 0 else math.nan,
        empirical_sd=sd([number(r, "estimate") for r in ok]), mean_spread=mean([number(r, "spread") for r in ok]),
        rms_spread=math.sqrt(mean([number(r, "spread")**2 for r in ok])) if nk else math.nan)
    for field in ("width", "interval_score", "within_variance", "between_variance", "df", "mcse_mean", "mcse_lower", "mcse_upper"):
        values = [number(r, field) for r in ok if math.isfinite(number(r, field))]
        ans["mean_"+field] = mean(values)
        ans["mean_"+field+"_mcse"] = mcse(values)
    return ans


def verify_summaries(out, records):
    groups = {}
    for row in records:
        groups.setdefault((row["cell"], row["method"]), []).append(row)
    actual = index_rows(read_csv(out/"summary.csv"), ("cell", "method"), "summary")
    require(set(actual) == set(groups), "Summary grid mismatch")
    checked = 0
    expected = {}
    for k, rows in groups.items():
        expected[k] = expected_summary(rows)
        for field, value in expected[k].items():
            close(number(actual[k], field), value, f"Summary {k}/{field}")
            checked += 1
    return groups, expected, checked


def boolean(value, label):
    require(value in ("TRUE", "FALSE"), f"Invalid Boolean: {label}")
    return value == "TRUE"


def verify_diagnostics(out, settings, fits, fit_map):
    rows = read_csv(out/"diagnostics.csv")
    indexed = index_rows(rows, ("cell", "replicate", "fit_id", "parameter"), "diagnostic")
    groups = {}
    for row in indexed.values():
        k = key(row, "fit_id")
        require(k in fit_map and k[2] != "oracle", f"Unexpected diagnostic fit: {k}")
        groups.setdefault(k, []).append(row)
    base_parameters = {"beta_T:(Intercept)", "beta_T:x", "beta_W:(Intercept)", "beta_W:x",
        "sigma_TT", "sigma_TW", "sigma_WW", "psi_TT", "psi_TW", "psi_WW",
        "conditional_target_mean_recorded", "conditional_target_mean_absent"}
    for k, fit in fit_map.items():
        require(fit["fit_status"] in ("ok", "error"), f"Unknown fit status: {k}")
        scale = .25 if k[2].endswith("quarter") else 4 if k[2].endswith("four") else 1
        close(number(fit, "prior_scale"), scale, f"Fit prior: {k}")
        require(math.isfinite(number(fit, "elapsed_seconds")) and number(fit, "elapsed_seconds") >= 0,
                f"Invalid fit timing: {k}")
        if fit["fit_status"] == "error":
            require(fit["error_class"] and fit["message"], f"Fit refusal lost: {k}")
        else:
            require(not fit["error_class"] and not fit["message"], f"Successful fit has error: {k}")
        if k[2] == "oracle":
            require(k not in groups and not fit["attempts"], f"Oracle has MCMC diagnostics: {k}")
            for field in ("max_rhat", "min_ess_bulk", "min_ess_tail", "iterations_used"):
                require(math.isnan(number(fit, field)), f"Oracle has MCMC metadata: {k}/{field}")
            continue
        if k not in groups:
            require(fit["fit_status"] == "error" and fit["error_class"] != "catmech_convergence_error",
                    f"Fit diagnostics missing: {k}")
            require(not fit["attempts"] and math.isnan(number(fit, "iterations_used")), f"Partial diagnostics missing: {k}")
            continue
        diag = groups[k]
        expected = base_parameters | ({"beta_T_original:(Intercept)"} if k[2] != "old" else set())
        require({r["parameter"] for r in diag} == expected, f"Diagnostic parameter set changed: {k}")
        for field, source, op in (("max_rhat", "rhat", max), ("min_ess_bulk", "ess_bulk", min),
                                  ("min_ess_tail", "ess_tail", min)):
            values = [number(r, source) for r in diag if math.isfinite(number(r, source))]
            close(number(fit, field), op(values) if values else math.nan, f"Diagnostic extremum: {k}/{field}")
        passes = all(math.isfinite(number(r, "rhat")) and number(r, "rhat") <= 1.05
            and math.isfinite(number(r, "ess_bulk")) and number(r, "ess_bulk") >= 400
            and math.isfinite(number(r, "ess_tail")) and number(r, "ess_tail") >= 400 for r in diag)
        attempts = json.loads(fit["attempts"])
        require(isinstance(attempts, list) and attempts, f"Attempt history missing: {k}")
        initial, maximum = (2000, 8000) if k[2] == "old" else (4000, 32000)
        schedule = [initial*2**i for i in range(len(attempts))]
        require(schedule[-1] <= maximum and [a["iterations"] for a in attempts] == schedule,
                f"Unplanned retry schedule: {k}")
        close(number(fit, "iterations_used"), schedule[-1], f"Final iteration count: {k}")
        require(all(not a["convergence_pass"] for a in attempts[:-1]), f"Fit rerun after convergence: {k}")
        require(attempts[-1]["convergence_pass"] == passes, f"Final diagnostic decision changed: {k}")
        if fit["fit_status"] == "ok":
            require(passes, f"Returned fit failed convergence: {k}")
        elif fit["error_class"] == "catmech_convergence_error":
            require(not passes and schedule[-1] == maximum, f"Premature convergence refusal: {k}")
    return len(rows)


def verify_completion(out, settings, completed, records, fits):
    counts = dict(cohorts=len({(r["cell"], r["replicate"]) for r in records}),
        analysis_requests=len(records), fitted_workflow_requests=len(fits),
        mcmc_fit_requests=sum(f["fit_id"] != "oracle" for f in fits),
        returned_intervals=sum(r["status"] == "ok" for r in records),
        analysis_failures=sum(r["status"] != "ok" for r in records),
        fit_failures=sum(f["fit_status"] != "ok" for f in fits))
    for field, value in counts.items():
        require(completed[field] == value, f"Completion count changed: {field}")
    require(completed["elapsed_seconds"] >= 0 and completed["started_utc"] and completed["finished_utc"],
            "Completion timing missing")
    fields = ("cell", "replicate", "method", "fit_id", "failure_stage", "error_class", "message", "warnings")
    errors = index_rows(read_csv(out/"errors.csv"), ("cell", "replicate", "method"), "error")
    expected = {tuple(r[f] for f in ("cell", "replicate", "method")): {f: r[f] for f in fields}
                for r in records if r["status"] != "ok"}
    require(errors == expected, "Error export omits or changes refusals")
    require((out/"session-info.txt").is_file() and (out/"session-info.txt").stat().st_size > 0,
            "Session/dependency record missing")
    if settings["mode"] == "pilot":
        for name in ("summary.csv", "paired.csv", "prior-triples.csv", "prior-summary.csv",
                     "calibration-bootstrap.csv", "availability.csv", "adoption.json"):
            require(not (out/name).exists(), f"Pilot exports accuracy: {name}")
    return counts


def verify_paired(out, records, record_map):
    actual = index_rows(read_csv(out/"paired.csv"), ("cell", "method", "baseline", "metric"), "paired")
    expected = {}
    for cell in CELLS:
        comparisons = [("new_cv3", "old_cv3"), ("direct", "new_cv3"), ("direct", "old_cv3")]
        if cell in PRIOR_CELLS:
            comparisons += [(a+s, a) for a in ("new_cv3", "direct") for s in ("_prior_quarter", "_prior_four")]
        for method, baseline in comparisons:
            reps = sorted(k[1] for k in record_map if k[0] == cell and k[2] == method)
            pairs = [(record_map[cell, r, method], record_map[cell, r, baseline]) for r in reps]
            complete = [(a, b) for a, b in pairs if a["status"] == b["status"] == "ok"]
            for metric in ("coverage_per_request", "coverage_complete_pairs", "error", "squared_error", "width", "interval_score"):
                chosen = pairs if metric == "coverage_per_request" else complete
                field = "covered" if metric.startswith("coverage") else "error" if metric == "squared_error" else metric
                power = 2 if metric == "squared_error" else 1
                delta = [number(a, field)**power-number(b, field)**power for a, b in chosen]
                expected[cell, method, baseline, metric] = dict(pairs=len(chosen), difference=mean(delta), mcse=mcse(delta))
    require(set(actual) == set(expected), "Paired comparison grid mismatch")
    for k, values in expected.items():
        for field, value in values.items():
            close(number(actual[k], field), value, f"Paired contrast {k}/{field}")
    return len(expected)


def verify_prior(out, settings, record_map):
    triples = index_rows(read_csv(out/"prior-triples.csv"), ("cell", "replicate"), "prior triple")
    expected_keys = {(c, str(r)) for c in PRIOR_CELLS for r in range(1, settings["prior_reps"]+1)}
    require(set(triples) == expected_keys, "Prior triple grid mismatch")
    computed = {}
    for c, r in expected_keys:
        row = triples[c, r]
        answer = {}
        for prefix, method in (("direct", "direct"), ("cv3", "new_cv3")):
            selected = [record_map[c, int(r), method+s] for s in ("", "_prior_quarter", "_prior_four")]
            complete = all(x["status"] == "ok" for x in selected)
            answer[prefix+"_complete"] = complete
            require(boolean(row[prefix+"_complete"], f"{c}/{r}/{prefix}") == complete, "Prior completeness changed")
            sensitivity = mean([abs(math.log(number(x, "width")/number(selected[0], "width"))) for x in selected[1:]]) if complete else math.nan
            shift = mean([abs(number(x, "estimate")-number(selected[0], "estimate")) for x in selected[1:]]) if complete else math.nan
            answer[prefix+"_sensitivity"], answer[prefix+"_estimate_shift"] = sensitivity, shift
            close(number(row, prefix+"_sensitivity"), sensitivity, f"Prior sensitivity {c}/{r}/{prefix}")
            close(number(row, prefix+"_estimate_shift"), shift, f"Prior estimate shift {c}/{r}/{prefix}")
        answer["difference"] = answer["direct_sensitivity"]-answer["cv3_sensitivity"]
        close(number(row, "difference"), answer["difference"], f"Prior paired difference {c}/{r}")
        computed[c, r] = answer
    summaries = index_rows(read_csv(out/"prior-summary.csv"), ("cell",), "prior summary")
    require(set(summaries) == {(c,) for c in PRIOR_CELLS}, "Prior summary grid mismatch")
    means, ses = [], []
    for c in PRIOR_CELLS:
        rows = [v for k, v in computed.items() if k[0] == c]
        good = [v for v in rows if v["direct_complete"] and v["cv3_complete"]]
        diff = [v["difference"] for v in good]
        values = dict(requested_triples=len(rows), direct_triples=sum(v["direct_complete"] for v in rows),
            cv3_triples=sum(v["cv3_complete"] for v in rows), complete_pairs=len(good),
            mean_direct_sensitivity=mean([v["direct_sensitivity"] for v in good]),
            mean_cv3_sensitivity=mean([v["cv3_sensitivity"] for v in good]), difference=mean(diff), mcse=mcse(diff))
        for field, value in values.items():
            close(number(summaries[c,], field), value, f"Prior summary {c}/{field}")
        means.append(values["difference"]); ses.append(values["mcse"])
    delta, se = mean(means), math.sqrt(sum(s*s for s in ses))/3
    no_positive = all(math.isfinite(d) and d <= 0 for d in means)
    return dict(difference=delta, mcse=se, lower=delta-Z975*se, upper=delta+Z975*se,
                no_positive_cell=no_positive, pass_=math.isfinite(delta) and math.isfinite(se) and delta+Z975*se < 0 and no_positive)


def bootstrap_indices(seed, sizes, replicates):
    """Replay only R's specified sample.int stream; never load study code."""
    code = '''args <- commandArgs(trailingOnly = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(as.integer(args[1]))
sizes <- as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]])
con <- file(args[4], "wb")
for (b in seq_len(as.integer(args[3]))) for (n in sizes)
  writeBin(as.integer(sample.int(n, n, replace = TRUE) - 1L), con, size = 4L, endian = "little")
close(con)
'''
    with tempfile.TemporaryDirectory(prefix="catmech-bootstrap-rng-") as temp:
        script, binary = Path(temp)/"indices.R", Path(temp)/"indices.bin"
        script.write_text(code)
        subprocess.run(["Rscript", "--vanilla", str(script), str(seed), ",".join(map(str, sizes)), str(replicates), str(binary)],
                       check=True, capture_output=True, text=True)
        indices = array.array("i", binary.read_bytes())
    require(indices.itemsize == 4, "Platform integer width cannot replay bootstrap")
    if sys.byteorder != "little":
        indices.byteswap()
    require(len(indices) == replicates*sum(sizes), "Bootstrap index stream truncated")
    return indices


def verify_calibration(out, settings, record_map, summaries):
    methods = ["old_cv3", "new_cv3", "direct"]
    cond = {m: {c: summaries[c, m]["coverage_conditional"] for c in PRIMARY} for m in methods}
    distances = {m: mean([abs(v-.95) for v in cond[m].values()]) for m in methods}
    biases = {c: summaries[c, "direct"]["bias"] for c in PRIMARY}
    matrices = []
    for c in PRIMARY:
        reps = sorted(k[1] for k in record_map if k[0] == c and k[2] == "direct")
        matrices.append([[number(record_map[c, r, m], "covered") if record_map[c, r, m]["status"] == "ok" else math.nan
                          for m in methods] for r in reps])
    indices = bootstrap_indices(settings["bootstrap_seed"], [len(z) for z in matrices], settings["bootstrap_reps"])
    actual = index_rows(read_csv(out/"calibration-bootstrap.csv"), ("bootstrap",), "bootstrap")
    require(set(actual) == {(str(b),) for b in range(1, settings["bootstrap_reps"]+1)}, "Bootstrap grid mismatch")
    values = [[], []]; offset = 0
    for b in range(1, settings["bootstrap_reps"]+1):
        by_cell = []
        for matrix in matrices:
            selected = indices[offset:offset+len(matrix)]; offset += len(matrix)
            require(all(0 <= j < len(matrix) for j in selected), "Bootstrap index outside cohort")
            by_cell.append([abs(mean([matrix[j][m] for j in selected if math.isfinite(matrix[j][m])])-.95) for m in range(3)])
        scores = [mean([c[m] for c in by_cell]) for m in range(3)]
        for j, field in enumerate(("direct_minus_old", "direct_minus_new_cv3")):
            delta = scores[2]-scores[j]
            values[j].append(delta)
            close(number(actual[str(b),], field), delta, f"Paired bootstrap {b}/{field}")
    comparisons = {}
    for j, field in enumerate(("direct_minus_old", "direct_minus_new_cv3")):
        finite = all(math.isfinite(v) for v in values[j])
        comparisons[field] = dict(difference=distances["direct"]-distances[methods[j]],
            lower=quantile(values[j], .025) if finite else math.nan,
            upper=quantile(values[j], .975) if finite else math.nan)
    coverage_ok = all(math.isfinite(v) and .925 <= v <= .975 for v in cond["direct"].values())
    bias_ok = all(math.isfinite(v) and abs(v) <= .02 for v in biases.values())
    finite = all(math.isfinite(v) for z in cond.values() for v in z.values()) and all(math.isfinite(v) for v in biases.values())
    passed = finite and coverage_ok and bias_ok and all(distances["direct"] < distances[m] for m in methods[:2]) and all(
        math.isfinite(v["upper"]) and v["upper"] < 0 for v in comparisons.values())
    return dict(conditional_coverage=cond, methods=methods, mean_absolute_distance=distances, direct_bias=biases,
                coverage_band_pass=coverage_ok, bias_bound_pass=bias_ok, **comparisons, pass_=passed)


def verify_availability(out, settings, fit_map):
    actual = index_rows(read_csv(out/"availability.csv"), ("cell",), "availability")
    require(set(actual) == {(c,) for c in CELLS[:-1]}, "Availability grid mismatch")
    deltas, ses, rate_ok, no_increase = [], [], [], []
    for c in CELLS[:-1]:
        reps = range(1, settings["primary_reps"]+1)
        old = [fit_map[c, r, "old"]["fit_status"] != "ok" for r in reps]
        new = [fit_map[c, r, "new"]["fit_status"] != "ok" for r in reps]
        differences = [int(a)-int(b) for a, b in zip(new, old)]
        expected = dict(requests=len(old), old_failures=sum(old), new_failures=sum(new), old_failure_rate=mean(old),
                        new_failure_rate=mean(new), difference=mean(differences), mcse=mcse(differences))
        for field, value in expected.items():
            close(number(actual[c,], field), value, f"Fit availability {c}/{field}")
        deltas.append(expected["difference"]); ses.append(expected["mcse"])
        rate_ok.append(expected["new_failure_rate"] <= .01)
        no_increase.append(expected["new_failure_rate"] <= expected["old_failure_rate"])
    delta, se = mean(deltas), math.sqrt(sum(s*s for s in ses))/6
    return dict(difference=delta, mcse=se, lower=delta-Z975*se, upper=delta+Z975*se,
        every_return_rate_at_least_99=all(rate_ok), no_observed_failure_increase=all(no_increase),
        pass_=all(rate_ok) and all(no_increase) and math.isfinite(se) and delta+Z975*se < 0)


def check_json(actual, expected, path="adoption"):
    if isinstance(expected, dict):
        expected = {"pass" if k == "pass_" else k: v for k, v in expected.items()}
        require(isinstance(actual, dict) and set(actual) == set(expected), f"JSON fields changed: {path}")
        for name, value in expected.items():
            check_json(actual[name], value, path+"/"+name)
    elif isinstance(expected, (str, bool, list)):
        require(actual == expected, f"JSON value changed: {path}")
    elif math.isnan(expected):
        require(actual is None, f"Expected absent JSON value: {path}")
    else:
        close(float(actual) if actual is not None else math.nan, expected, path)


def verify(out, root):
    settings, completed, records, fits, record_map, fit_map = verify_records(out, root)
    counts = verify_completion(out, settings, completed, records, fits)
    diagnostics = verify_diagnostics(out, settings, fits, fit_map)
    checks = dict(stage=settings["mode"], **counts, diagnostic_rows=diagnostics,
                  source_files=len(read_csv(out/"source-hashes-before.csv")))
    if settings["mode"] != "pilot":
        _, summaries, checked = verify_summaries(out, records)
        checks["summary_values_checked"] = checked
        checks["paired_contrasts_checked"] = verify_paired(out, records, record_map)
        prior = verify_prior(out, settings, record_map)
        calibration = verify_calibration(out, settings, record_map, summaries)
        availability = verify_availability(out, settings, fit_map)
        adoption = dict(stage=settings["mode"], calibration=calibration, prior_sensitivity=prior, availability=availability,
            activate_prior_development_branch=settings["mode"] == "development" and not prior["pass_"],
            interpretation="development_only_no_confirmation_claim" if settings["mode"] == "development" else
                "separate_prespecified_screens_no_universal_validity")
        check_json(json.loads((out/"adoption.json").read_text()), adoption)
        checks["bootstrap_replicates_recomputed"] = settings["bootstrap_reps"]
        checks["adoption_screens_recomputed"] = 3
    return checks


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def tamper_tests(out, root):
    """Alter temporary copies only. Each corruption must cause a hard refusal."""
    tests = []
    mode = json.loads((out/"settings.json").read_text())["mode"]
    names = ["missing_entire_cell", "missing_fit", "changed_seed", "changed_source_hash"]
    if mode != "pilot":
        names += ["missing_summary", "changed_summary", "flipped_raw_covered", "changed_prior", "changed_adoption"]
    for name in names:
        with tempfile.TemporaryDirectory(prefix="catmech-verifier-tamper-") as temp:
            dest = Path(temp)/"copy"; shutil.copytree(out, dest)
            filename = "records.csv"
            if name == "missing_fit": filename = "fit-records.csv"
            if name in ("missing_summary", "changed_summary"): filename = "summary.csv"
            if name == "changed_prior": filename = "prior-triples.csv"
            if name == "changed_source_hash": filename = "source-hashes-before.csv"
            if name == "changed_adoption":
                p = dest/"adoption.json"; obj = json.loads(p.read_text())
                obj["calibration"]["pass"] = not obj["calibration"]["pass"]
                p.write_text(json.dumps(obj))
            else:
                rows = read_csv(dest/filename); fields = list(rows[0])
                if name == "missing_entire_cell": rows = [r for r in rows if r["cell"] != CELLS[0]]
                elif name in ("missing_fit", "missing_summary"): rows.pop(0)
                elif name == "changed_seed": rows[0]["data_seed"] = str(int(rows[0]["data_seed"])+1)
                elif name == "changed_source_hash": rows[0]["sha256"] = "0"*64
                elif name == "changed_summary": rows[0]["coverage_conditional"] = "0.1234567"
                elif name == "flipped_raw_covered":
                    row = next(r for r in rows if r["status"] == "ok")
                    row["covered"] = str(1-int(row["covered"]))
                elif name == "changed_prior":
                    row = next(r for r in rows if r["direct_complete"] == "TRUE")
                    row["direct_sensitivity"] = str(number(row, "direct_sensitivity")+.1234567)
                write_csv(dest/filename, rows, fields)
            try:
                verify(dest, root)
            except (AssertionError, KeyError, ValueError, FileNotFoundError) as error:
                tests.append(dict(test=name, rejected=True, reason=str(error)))
            else:
                raise AssertionError(f"Tampering escaped detection: {name}")
    return tests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Completed stage output directory")
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1]/"validation/v0.9.0-source",
                        help="Frozen v0.9.0 source tree matching the stage's before/after hashes")
    parser.add_argument("--tamper-tests", action="store_true", help="Confirm altered temporary copies are rejected")
    parser.add_argument("--report", type=Path, help="Write JSON verification report (does not overwrite study evidence)")
    args = parser.parse_args()
    try:
        # Numerical controls independent of the study outcomes.
        close(t975(1), 12.706204736432095, "Cauchy t quantile control")
        close(t975(10), 2.2281388519649385, "t10 quantile control")
        close(quantile([0, 1, 4, 9], .25), .75, "Type-7 quantile control")
        result = dict(status="PASS", output=str(args.output.resolve()), source_root=str(args.source_root.resolve()),
                      checks=verify(args.output.resolve(), args.source_root.resolve()),
                      scope="Independent reporting/grid/arithmetic audit; not proof of model validity or sampler correctness.")
        if args.tamper_tests:
            result["tamper_tests"] = tamper_tests(args.output.resolve(), args.source_root.resolve())
        text = json.dumps(result, indent=2, allow_nan=False)+"\n"
        if args.report:
            require(not args.report.exists(), "Refusing to overwrite a verifier report")
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(text)
        print(text, end="")
    except (AssertionError, KeyError, ValueError, FileNotFoundError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
