#!/usr/bin/env python3
"""Independently verify the fixed joint-cluster study grid and scalar summaries.

Only Python's standard library is used. This does not validate the imputer's
mathematics; it checks completeness, provenance and reported Monte Carlo math.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

CELLS = ["gaussian_balanced", "gaussian_unbalanced", "near_zero_icc",
         "observed_mar", "null_slope", "student_t5"]
SENSITIVITY = {"gaussian_balanced", "near_zero_icc", "observed_mar"}
MAIN = ["joint", "plain_recorded", "oracle"]
PRIOR = ["joint_prior_quarter", "joint_prior_four"]


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row, field):
    value = row[field]
    return float(value) if value not in ("", "NA", "NaN") else math.nan


def close(actual, expected, label):
    if math.isnan(expected):
        require(math.isnan(actual), f"{label}: expected absent, got {actual}")
    else:
        require(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=2e-11, abs_tol=2e-12),
                f"{label}: {actual} != {expected}")


def mean(values):
    return statistics.mean(values) if values else math.nan


def sd(values):
    return statistics.stdev(values) if len(values) > 1 else math.nan


def mcse(values):
    return sd(values) / math.sqrt(len(values)) if values else math.nan


def wilson(successes, n):
    if not n:
        return math.nan, math.nan
    p = successes / n
    z = statistics.NormalDist().inv_cdf(.975)
    denominator = 1 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return max(0, center-half), min(1, center+half)


def verify(out, source_root):
    config = json.loads((out / "settings.json").read_text())
    complete = json.loads((out / "COMPLETE.json").read_text())
    mode = config["mode"]
    require(mode in ("pilot", "heldout"), "Unexpected mode")
    pilot = mode == "pilot"
    reps = 3 if pilot else 200
    prior_reps = 3 if pilot else 100
    require(config["reps"] == reps and config["m"] == 40, "Protocol m/repetitions changed")
    require(config["cells"] == CELLS, "Cell list/order changed")
    require(set(config["sensitivity_cells"]) == SENSITIVITY, "Prior sensitivity cells changed")
    require(config["sensitivity_reps"] == prior_reps, "Prior sensitivity repetitions changed")
    require(config["control"] == dict(chains=4, burnin=1000, iterations=2000,
                                      max_iterations=8000, prior_scale=1), "Chain controls changed")
    require(config["pan_version"] == "2.0", "Unexpected pan version; inspect provenance")
    data_base, imp_base = (1510000000, 1520000000) if pilot else (1530000000, 1540000000)
    require(config["data_seed_base"] == data_base and config["imputation_seed_base"] == imp_base,
            "Seed bases changed")
    expected = set()
    for cell in CELLS:
        for rep in range(1, reps+1):
            methods = MAIN + (PRIOR if cell in SENSITIVITY and rep <= prior_reps else [])
            expected.update((cell, rep, method) for method in methods)
    rows = read_csv(out / "records.csv")
    keyed = {(r["cell"], int(r["replicate"]), r["method"]): r for r in rows}
    require(len(keyed) == len(rows), "Duplicate request keys")
    require(set(keyed) == expected,
            f"Request grid mismatch; missing {len(expected-set(keyed))}, extra {len(set(keyed)-expected)}")
    require(complete["status"] == "complete" and complete["requests"] == len(rows), "Completion count mismatch")
    require(complete["mode"] == mode, "Completion mode mismatch")
    fingerprints = {r["source_fingerprint"] for r in rows}
    require(fingerprints == {config["source_fingerprint"]} == {complete["source_fingerprint"]},
            "Worker/source fingerprints differ")
    before = read_csv(out / "source-hashes-before.csv")
    after = read_csv(out / "source-hashes-after.csv")
    require(before == after, "Sources changed during run")
    expected_sources = {"DESCRIPTION", "NAMESPACE", "validation/JOINT_CLUSTER_PLAN.md",
                        "validation/run-joint-cluster.R"}
    expected_sources.update(str(f.relative_to(source_root)) for f in (source_root / "R").glob("*.R"))
    require({r["file"] for r in before} == expected_sources, "Hashed source set differs from expected production set")
    for r in before:
        path = source_root / r["file"]
        require(path.is_file(), f"Missing source: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == r["sha256"], f"Source changed: {r['file']}")
    input_hashes = {}
    nerrors = 0
    for key, r in keyed.items():
        cell, rep, method = key
        index = CELLS.index(cell)+1
        require(int(r["cell_index"]) == index, f"Cell index: {key}")
        require(int(r["data_seed"]) == data_base + index*100000 + rep, f"Data seed: {key}")
        require(int(r["imputation_seed"]) == imp_base + index*100000 + rep, f"Imputation seed: {key}")
        require(int(r["n"]) == 240 and int(r["sites"]) == 24, f"Cohort size: {key}")
        require(int(r["wholly_unaudited_sites"]) >= 6, f"Missing unaudited sites: {key}")
        if cell != "observed_mar":
            require(int(r["gold_n"]) == 60, f"Audit size: {key}")
            require(int(r["mixed_stochastic_sites"]) >= 6, f"Mixed site support: {key}")
        truth = 0 if cell == "null_slope" else .7
        close(number(r, "truth"), truth, f"Truth: {key}")
        cohort_key = cell, rep
        old = input_hashes.setdefault(cohort_key, r["input_sha256"])
        require(old == r["input_sha256"], f"Different method inputs: {key}")
        require(r["status"] in ("ok", "error"), f"Unknown status: {key}")
        failed = r["status"] == "error"
        nerrors += failed
        if failed:
            require(bool(r["error_class"]) and bool(r["message"]), f"Failure details missing: {key}")
        if pilot:
            for field in ("estimate", "se", "df", "lower", "upper", "width", "error", "covered"):
                require(math.isnan(number(r, field)), f"Pilot accuracy exported: {key}/{field}")
        elif failed:
            close(number(r, "covered"), 0, f"Failure denominator: {key}")
            for field in ("estimate", "se", "df", "lower", "upper", "width", "error"):
                require(math.isnan(number(r, field)), f"Failure has result: {key}/{field}")
        else:
            estimate, lower, upper = (number(r, f) for f in ("estimate", "lower", "upper"))
            require(all(math.isfinite(x) for x in (estimate, lower, upper)), f"Nonfinite interval: {key}")
            require(upper > lower and number(r, "se") > 0 and number(r, "df") > 0, f"Degenerate interval: {key}")
            close(number(r, "error"), estimate-truth, f"Error: {key}")
            close(number(r, "width"), upper-lower, f"Width: {key}")
            close(number(r, "covered"), float(lower <= truth <= upper), f"Coverage: {key}")
    require(complete["errors"] == nerrors and complete["intervals"] == len(rows)-nerrors,
            "Completion status counts differ")
    errors = read_csv(out / "errors.csv")
    require({(r["cell"], int(r["replicate"]), r["method"]) for r in errors} ==
            {k for k, r in keyed.items() if r["status"] == "error"}, "Error table incomplete")
    if pilot:
        require(not (out / "summary.csv").exists() and not (out / "paired.csv").exists(), "Pilot accuracy summary exists")
        return len(rows), nerrors, 0

    groups = {}
    for r in rows:
        groups.setdefault((r["cell"], r["method"]), []).append(r)
    summaries = read_csv(out / "summary.csv")
    by_group = {(r["cell"], r["method"]): r for r in summaries}
    require(len(by_group) == len(summaries) and set(by_group) == set(groups), "Summary grid mismatch")
    checked = 0
    for key, z in groups.items():
        s = by_group[key]
        ok = [r for r in z if r["status"] == "ok"]
        n, nk = len(z), len(ok)
        covered = sum(number(r, "covered") for r in z)
        pc = covered/n
        cc = covered/nk if nk else math.nan
        failure = (n-nk)/n
        pw, cw, fw = wilson(covered, n), wilson(covered, nk), wilson(n-nk, n)
        errors = [number(r, "error") for r in ok]
        sq = [e*e for e in errors]
        rmse = math.sqrt(mean(sq)) if sq else math.nan
        widths = [number(r, "width") for r in ok]
        expected_summary = dict(requests=n, intervals=nk, failures=n-nk,
            coverage_per_request=pc, coverage_per_request_mcse=math.sqrt(pc*(1-pc)/n),
            coverage_per_request_lower=pw[0], coverage_per_request_upper=pw[1],
            coverage_conditional=cc, coverage_conditional_mcse=math.sqrt(cc*(1-cc)/nk) if nk else math.nan,
            coverage_conditional_lower=cw[0], coverage_conditional_upper=cw[1],
            failure_rate=failure, failure_rate_mcse=math.sqrt(failure*(1-failure)/n),
            failure_rate_lower=fw[0], failure_rate_upper=fw[1],
            bias=mean(errors), bias_mcse=mcse(errors), rmse=rmse,
            rmse_mcse=mcse(sq)/(2*rmse) if rmse > 0 else math.nan,
            empirical_sd=sd([number(r, "estimate") for r in ok]),
            mean_se=mean([number(r, "se") for r in ok]), mean_width=mean(widths), mean_width_mcse=mcse(widths),
            total_seconds=sum(number(r, "elapsed_seconds") for r in z),
            median_seconds=statistics.median(number(r, "elapsed_seconds") for r in z))
        for field, value in expected_summary.items():
            close(number(s, field), value, f"Summary {key}/{field}")
            checked += 1

    pairs = read_csv(out / "paired.csv")
    expected_pairs = {}
    for cell in CELLS:
        methods = {method for c, method in groups if c == cell}
        for method in methods - {"plain_recorded"}:
            baseline = "joint" if method in PRIOR else "plain_recorded"
            a = {int(r["replicate"]): r for r in groups[cell, method]}
            b = {int(r["replicate"]): r for r in groups[cell, baseline]}
            shared = sorted(set(a) & set(b))
            complete_pairs = [r for r in shared if a[r]["status"] == b[r]["status"] == "ok"]
            for metric in ("coverage_per_request", "coverage_complete_pairs", "error", "squared_error", "width"):
                selected = shared if metric == "coverage_per_request" else complete_pairs
                field = "covered" if metric.startswith("coverage_") else "error" if metric == "squared_error" else metric
                values = []
                for rep in selected:
                    av, bv = number(a[rep], field), number(b[rep], field)
                    values.append(av*av-bv*bv if metric == "squared_error" else av-bv)
                expected_pairs[cell, method, baseline, metric] = len(selected), mean(values), mcse(values)
    actual_pairs = {(r["cell"], r["method"], r["baseline"], r["metric"]): r for r in pairs}
    require(len(actual_pairs) == len(pairs) and set(actual_pairs) == set(expected_pairs), "Paired table grid mismatch")
    for key, expected_pair in expected_pairs.items():
        for field, value in zip(("pairs", "difference", "mcse"), expected_pair):
            close(number(actual_pairs[key], field), value, f"Paired {key}/{field}")
            checked += 1
    return len(rows), nerrors, checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    requests, errors, summaries = verify(args.output.resolve(), args.source_root.resolve())
    print(f"PASS: {requests} planned requests; {errors} recorded failures; {summaries} independently recomputed scalar summaries; source hashes match.")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, KeyError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
