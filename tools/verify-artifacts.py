#!/usr/bin/env python3
"""Check saved evidence against code and replicate-level data; no validity claim."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def verify_summary(summary, records):
    for row in summary:
        group = [r for r in records if r["scenario"] == row["scenario"]
                 and r["method"] == row["method"]]
        good = [r for r in group if r["status"] == "ok"]
        assert len(group) == int(row["requested"])
        assert len(good) == int(row["successful"])
        assert len(group)-len(good) == int(row["failures"])
        error = [float(r["estimate"])-float(row["estimand_truth"]) for r in good]
        cover = [float(r["lower"]) <= float(row["estimand_truth"]) <= float(r["upper"])
                 for r in good]
        calculated = {
            "bias": statistics.mean(error),
            "bias_mcse": statistics.stdev(error)/math.sqrt(len(error)),
            "rmse": math.sqrt(statistics.mean(e*e for e in error)),
            "coverage": statistics.mean(cover),
            "mean_se": statistics.mean(float(r["se"]) for r in good),
            "mean_interval_width": statistics.mean(float(r["upper"])-float(r["lower"]) for r in good),
        }
        calculated["coverage_mcse"] = math.sqrt(calculated["coverage"]*(1-calculated["coverage"])/len(good))
        for name, value in calculated.items():
            assert math.isclose(float(row[name]), value, rel_tol=1e-11, abs_tol=1e-12), (row["scenario"], row["method"], name)


def check_hashes(manifest, base):
    for row in rows(manifest):
        filename = row.get("file", row.get("path"))
        actual = hashlib.sha256((base/filename).read_bytes()).hexdigest()
        assert actual == row["sha256"], filename


def close(saved, calculated, context):
    assert math.isclose(float(saved), calculated, rel_tol=1e-10, abs_tol=1e-11), context


def mean_se(values):
    return statistics.stdev(values) / math.sqrt(len(values))


def wilson(k, n):
    z = statistics.NormalDist().inv_cdf(.975)
    p = k/n
    denominator = 1 + z*z/n
    centre = (p + z*z/(2*n))/denominator
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/denominator
    return max(0, centre-half), min(1, centre+half)


def verify_coverage(folder, source_root=None):
    check_hashes(folder/"source-hashes.csv", source_root or ROOT)
    records = rows(folder/"replicates.csv")
    if folder.name.startswith("extension-"):
        verify_extension_plan(folder, records)
    keys = [(r["family"], r["scenario"], r["replicate"], r["method"]) for r in records]
    assert len(keys) == len(set(keys)), "Duplicate coverage request"
    groups = {}
    for r in records:
        groups.setdefault((r["family"], r["scenario"], r["method"]), []).append(r)
        if r["status"] == "ok":
            close(r["error"], float(r["estimate"])-float(r["truth"]), "record error")
            close(r["width"], float(r["upper"])-float(r["lower"]), "record width")
            assert int(r["coverage"]) == (float(r["lower"]) <= float(r["truth"]) <= float(r["upper"]))
    summary = rows(folder/"summary.csv")
    assert len(summary) == len(groups)
    for r in summary:
        key = r["family"], r["scenario"], r["method"]
        group = groups[key]
        good = [g for g in group if g["status"] == "ok"]
        n = len(good)
        assert int(r["requested"]) == len(group)
        assert int(r["returned_intervals"]) == n
        for field, status in [("refused", "refused"), ("convergence_failures", "convergence_failure"), ("other_errors", "error")]:
            assert int(r.get(field, 0)) == sum(g["status"] == status for g in group)
        close(r["interval_return_rate"], n/len(group), key)
        if not n:
            assert r["coverage_conditional"] == "NA"
            assert r["bias"] == "NA"
            continue
        error = [float(g["error"]) for g in good]
        cover = [int(g["coverage"]) for g in good]
        se = [float(g["se"]) for g in good]
        width = [float(g["width"]) for g in good]
        bias = statistics.mean(error)
        rmse = math.sqrt(statistics.mean(e*e for e in error))
        coverage = statistics.mean(cover)
        ci = wilson(sum(cover), n)
        diagnostic = statistics.mean(float(g["lower"])-bias <= float(g["truth"]) <= float(g["upper"])-bias for g in good)
        values = dict(bias=bias, bias_mcse=mean_se(error), rmse=rmse,
            rmse_mcse=mean_se([e*e for e in error])/(2*rmse),
            empirical_sd=statistics.stdev(float(g["estimate"]) for g in good),
            mean_se=statistics.mean(se), mean_se_mcse=mean_se(se),
            mean_width=statistics.mean(width), width_mcse=mean_se(width),
            coverage_conditional=coverage, coverage_mcse=math.sqrt(coverage*(1-coverage)/n),
            coverage_wilson_lower=ci[0], coverage_wilson_upper=ci[1],
            covering_interval_per_request=sum(cover)/len(group),
            bias_centered_coverage_diagnostic=diagnostic)
        for field, value in values.items():
            close(r[field], value, (key, field))
        if "primary_criterion_applicable" in r:
            primary = (r["family"] == "distribution" and r["scenario"] in
                ["missing_only", "perfect_flags", "false_positive_flags", "category_error", "null_effect", "strong_mar"]
                and r["method"] in ["normal_no_cat", "mice_normal_no_cat"]) or (
                r["family"] == "reference" and r["method"] in ["reference_normal", "mice_reference_normal"])
            assert (r["primary_criterion_applicable"] == "TRUE") == (primary and len(group) >= 1000)
            if primary and n >= 1000:
                assert (r["approximate_nominal_criterion"] == "TRUE") == (ci[0] >= .925 and ci[1] <= .975)
            else:
                assert r["approximate_nominal_criterion"] == "NA"
        else:
            assert (r["approximate_nominal_criterion"] == "TRUE") == (ci[0] >= .925 and ci[1] <= .975)
    for r in rows(folder/"request-status.csv"):
        group = groups[r["family"], r["scenario"], r["method"]]
        k = sum(g["status"] == r["status"] for g in group)
        assert int(r["count"]) == k and int(r["requested"]) == len(group)
        ci = wilson(k, len(group))
        close(r["rate"], k/len(group), "status rate")
        close(r["lower"], ci[0], "status Wilson lower")
        close(r["upper"], ci[1], "status Wilson upper")
    for r in rows(folder/"paired-comparisons.csv"):
        a = {g["replicate"]: g for g in groups[r["family"], r["scenario"], r["method"]] if g["status"] == "ok"}
        b = {g["replicate"]: g for g in groups[r["family"], r["scenario"], r["comparator"]] if g["status"] == "ok"}
        pairs = sorted(a.keys() & b.keys())
        assert int(r["complete_pairs"]) == len(pairs)
        for saved, source, squared, mcse_field in [
            ("error_difference", "error", False, "error_mcse"),
            ("mse_difference", "error", True, "mse_mcse"),
            ("coverage_difference", "coverage", False, "coverage_difference_mcse"),
            ("width_difference", "width", False, "width_difference_mcse")]:
            differences = [float(a[k][source])**(2 if squared else 1)-float(b[k][source])**(2 if squared else 1) for k in pairs]
            close(r[saved], statistics.mean(differences), saved)
            close(r[mcse_field], mean_se(differences), mcse_field)
    decomposition = rows(folder/"bias-decomposition.csv")
    for r in decomposition:
        close(r["selection_difference"], float(r["oracle_retained"])-float(r["oracle_full"]), "selection identity")
        close(r["residual_error_difference"], float(r["recorded_retained"])-float(r["oracle_retained"]), "residual-error identity")
        close(float(r["selection_difference"])+float(r["residual_error_difference"]),
              float(r["recorded_retained"])-float(r["oracle_full"]), "sum identity")
    for r in rows(folder/"bias-decomposition-summary.csv"):
        group = [g for g in decomposition if g["family"] == r["family"] and g["scenario"] == r["scenario"]]
        assert int(r["replicates"]) == len(group)
        for field, mcse_field in [("selection_difference", "selection_mcse"), ("residual_error_difference", "residual_error_mcse")]:
            values = [float(g[field]) for g in group]
            close(r[field], statistics.mean(values), field)
            close(r[mcse_field], mean_se(values), mcse_field)
    print(f"{folder.name}: {len(records)} unique requests, {len(summary)} recomputed groups, statuses and paired MCSEs: PASS")


def verify_extension_plan(folder, records):
    """Independently require every protocol cell, not merely every saved group."""
    settings = json.loads((folder/"settings.json").read_text())
    mode = settings["mode"]
    assert settings["m"] == 20 and settings["donors"] == 5
    assert settings["synthetic_n"] == 300 and settings["reference_n"] == 600
    assert settings["fixture_master_seed"] == 390173
    distribution = ["missing_only", "perfect_flags", "false_positive_flags", "category_error",
                    "null_effect", "heavy_tails", "heavy_category_error", "strong_mar"]
    ref = ["clean", "flagged_errors", "missed_flags", "site_shift", "rounded_records"]
    if mode == "distribution":
        cells = {("distribution", s): 1000 for s in distribution}
    elif mode == "reference":
        cells = {("reference", s): 1000 for s in ref}
    elif mode == "stress":
        cells = {("support", f"{s}_gold{fraction}"): 250
                 for s in ["missed_flags", "site_shift"] for fraction in [10, 50]}
        cells.update({("stress", "nonignorable_selection"): 1000,
                      ("stress", "noisy_reference"): 500})
    elif mode == "controls":
        cells = {("control", s): 20 for s in ["all_reference", "no_reference"]}
    else:
        raise AssertionError("Only final held-out extension modes are release evidence")
    methods = {
        "distribution": ["oracle", "trusted_observed", "normal_no_cat", "normal_cat",
                         "pmm_no_cat", "pmm_cat", "mice_normal_no_cat", "mice_normal_cat",
                         "mice_pmm_no_cat", "mice_pmm_cat"],
        "control": ["reference_normal", "reference_pmm"],
        "reference": ["oracle", "v02_mask_normal", "gold_only", "reference_normal",
                      "reference_pmm", "mice_reference_normal", "mice_reference_pmm"]}
    expected = {(family, scenario, method, str(rep))
                for (family, scenario), reps in cells.items()
                for method in methods.get(family, methods["reference"])
                for rep in range(1, reps+1)}
    actual = {(r["family"], r["scenario"], r["method"], r["replicate"]) for r in records}
    assert actual == expected and len(actual) == len(records), "Missing, duplicate or unplanned protocol requests"
    for r in records:
        rep = int(r["replicate"])
        if mode == "distribution":
            seed = 60000000 + 100000*(distribution.index(r["scenario"])+1) + rep
            assert int(r["data_seed"]) == seed and int(r["imputation_seed"]) == seed+1000000
            assert int(r["n"]) == 300
        else:
            assert int(r["data_seed"]) == 80100000+rep and int(r["imputation_seed"]) == 90000000+rep
            assert r["data_stream"] == f"master390173:heldout:{rep}"
            assert int(r["n"]) == 600
        if r["scenario"] == "no_reference":
            assert r["status"] == "refused" and int(r["available_n"]) == 0
        elif r["scenario"] == "all_reference":
            assert r["status"] == "ok" and int(r["available_n"]) == 600
            assert float(r["control_max_abs_difference"]) <= 1e-10
    if mode == "reference":
        baseline = {(r["method"], r["replicate"]): r for r in records if r["scenario"] == "clean"}
        for r in records:
            if r["method"] == "v02_mask_normal":
                continue
            first = baseline[r["method"], r["replicate"]]
            for field in ["available_n", "status", "estimate", "se", "df", "lower", "upper"]:
                assert r[field] == first[field], ("Same-reference invariance", field, r["method"])
    print(f"{folder.name}: all {len(expected)} prespecified requests, seeds and endpoint controls accounted for: PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path)
    args = parser.parse_args()
    benchmark = ROOT/"validation/benchmark"
    check_hashes(benchmark/"source-hashes.csv", ROOT/"validation/v0.1.0-source")
    records = rows(benchmark/"replicates.csv")
    summary = rows(benchmark/"summary.csv")
    assert len({(r["scenario"],r["replicate"],r["method"]) for r in records}) == len(records)
    verify_summary(summary,records)
    wrong = [dict(r) for r in summary]
    wrong[0]["bias"] = str(float(wrong[0]["bias"])+0.2)
    try:
        verify_summary(wrong,records)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: corrupted summary went undetected")
    if args.original:
        check_hashes(ROOT/"validation/original-sha256.csv",args.original)
        print("Original draft file hashes unchanged: PASS")
    print(f"Frozen v0.1 source hashes and independently recomputed {len(summary)} summary rows: PASS")
    print("Known-wrong bias mutation detected: PASS")
    for name in ["coverage-confirmatory", "coverage-clinical"]:
        verify_coverage(ROOT/"validation"/name, ROOT/"validation/v0.2.0-source")
    check_hashes(ROOT/"validation/model-oracle-final/source-hashes.csv", ROOT/"validation/v0.2.0-source")
    for r in rows(ROOT/"validation/model-oracle-final/posterior-moment-checks.csv"):
        assert r["pass"] == "TRUE", r
    print("Independent posterior oracle source hashes and moment checks: PASS")
    for name in ["extension-distribution", "extension-reference", "extension-stress", "extension-controls"]:
        verify_coverage(ROOT/"validation"/name, ROOT/"validation/v0.3.0-source")
    check_hashes(ROOT/"validation/reference-oracle/source-hashes.csv", ROOT)
    for filename in ["posterior-predictive-moments.csv", "data-contract-controls.csv", "negative-controls.csv"]:
        for r in rows(ROOT/"validation/reference-oracle"/filename):
            assert r.get("pass", r.get("detected")) == "TRUE", (filename, r)
    print("Independent reference predictive oracle, contracts and negative controls: PASS")
    check_hashes(ROOT/"validation/pooling-guards/source-hashes.csv", ROOT)
    for filename in ["scalar-numerical-equivalence.csv", "guard-controls.csv"]:
        checks = rows(ROOT/"validation/pooling-guards"/filename)
        assert len(checks) == (54 if filename.startswith("scalar") else 15)
        for r in checks:
            assert r.get("pass", r.get("bit_identical")) == "TRUE", (filename, r)
            if "max_absolute_difference" in r:
                assert float(r["max_absolute_difference"]) == 0
    print("Current pooling guards and frozen-v0.3 scalar numerical equivalence: PASS")
    extension = ROOT/"validation/extension-controls"
    incomplete = [r for r in rows(extension/"replicates.csv") if r["scenario"] != "no_reference"]
    try:
        verify_extension_plan(extension, incomplete)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: omitted planned scenario went undetected")
    print("Known-wrong deletion of an entire planned scenario detected: PASS")
    print("This verifies correspondence and repeatability, not inferential validity.")


if __name__ == "__main__":
    main()
