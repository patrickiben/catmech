#!/usr/bin/env python3
"""Check saved evidence against code and replicate-level data; no validity claim."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
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


def verify_study_bundle(folder, records):
    """Recompute every summary, status rate and paired contrast from replicate records."""
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
        assert int(r["requested"]) == len(group) and int(r["returned_intervals"]) == n
        for field, status in [("refused", "refused"), ("convergence_failures", "convergence_failure"), ("other_errors", "error")]:
            assert int(r[field]) == sum(g["status"] == status for g in group)
        close(r["interval_return_rate"], n/len(group), key)
        if not n:
            assert r["coverage_conditional"] == "NA" and r["bias"] == "NA"
            continue
        error = [float(g["error"]) for g in good]
        cover = [int(g["coverage"]) for g in good]
        se = [float(g["se"]) for g in good]
        width = [float(g["width"]) for g in good]
        bias = statistics.mean(error)
        rmse = math.sqrt(statistics.mean(e*e for e in error))
        coverage = statistics.mean(cover)
        ci = wilson(sum(cover), n)
        values = dict(bias=bias, bias_mcse=mean_se(error), rmse=rmse, empirical_sd=statistics.stdev(float(g["estimate"]) for g in good),
            mean_se=statistics.mean(se), mean_se_mcse=mean_se(se), mean_width=statistics.mean(width), width_mcse=mean_se(width),
            coverage_conditional=coverage, coverage_mcse=math.sqrt(coverage*(1-coverage)/n),
            coverage_wilson_lower=ci[0], coverage_wilson_upper=ci[1], covering_interval_per_request=sum(cover)/len(group))
        for field, value in values.items():
            close(r[field], value, (key, field))
        if "primary_criterion_applicable" in r:
            applies = r["family"] == "primary" and r["method"] == "contamination"
            assert (r["primary_criterion_applicable"] == "TRUE") == (applies and len(group) >= 1000)
            if applies and n >= 1000:
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
    return summary


def contamination_cells():
    """The frozen protocol, written out independently of the R scenario list."""
    primary = ["clean", "flagged_errors", "missed_flags", "site_shift_declared"]
    stress = ["reference_05pct", "reference_10pct", "reference_50pct", "rounded_records", "site_shift_undeclared",
              "targeted_audit", "nonignorable_clean", "nonignorable_missed", "noisy_declared", "noisy_wrong_sd",
              "noisy_undeclared", "heavy_tail_target", "nonlinear_target", "heteroskedastic_target", "error_bimodal",
              "error_gross_and_small", "error_multiplicative", "differential_outcome", "value_dependent_error",
              "flag_from_recorded", "assay_noise_tolerance", "assay_noise_no_tolerance", "recorded_missing_mnar"]
    alone = {"noisy_wrong_sd", "noisy_undeclared", "assay_noise_no_tolerance"}
    matching = {"missed_flags", "reference_05pct", "reference_10pct", "reference_50pct", "heavy_tail_target", "nonlinear_target"}
    every = ["oracle", "gold_only", "v02_mask_normal", "reference_normal", "cole_mi", "raking", "contamination"]
    expected = set()
    for scenario, reps in [(s, 1000) for s in primary] + [(s, 150) for s in stress]:
        methods = ["contamination"] if scenario in alone else list(every)
        if scenario in matching:
            methods += ["reference_pmm", "reference_pmm_abb"]
        expected |= {(scenario, method, str(rep)) for method in methods for rep in range(1, reps+1)}
    return expected


def verify_contamination_records(records):
    actual = {(r["scenario"], r["method"], r["replicate"]) for r in records}
    assert actual == contamination_cells() and len(actual) == len(records), "Missing, duplicate or unplanned requests"
    for r in records:
        rep = int(r["replicate"])
        assert int(r["data_seed"]) == 180100000+rep and int(r["imputation_seed"]) == 190000000+rep
        assert r["data_stream"] == f"master410273:heldout:{rep}"
        assert int(r["n"]) == 600 and float(r["truth"]) == 0.8


def verify_contamination_study(folder, source_root):
    check_hashes(folder/"source-hashes.csv", source_root)
    settings = json.loads((folder/"settings.json").read_text())
    assert settings["mode"] == "study" and settings["m"] == 20 and settings["fixture_master_seed"] == 410273
    assert settings["phase"] == "heldout" and not (folder/"checkpoint.rds").exists()
    records = rows(folder/"replicates.csv")
    verify_contamination_records(records)
    summary = verify_study_bundle(folder, records)
    scores = rows(folder/"corruption-scores.csv")
    for r in rows(folder/"corruption-scores-summary.csv"):
        group = [g for g in scores if g["scenario"] == r["scenario"]]
        assert int(r["replicates"]) == len(group)
        close(r["mean_actual_errors"], statistics.mean(float(g["actual_errors"]) for g in group), "actual errors")
        close(r["mean_expected_errors"], statistics.mean(float(g["expected_errors"]) for g in group), "expected errors")
        close(r["brier"], statistics.mean(float(g["brier"]) for g in group), "brier")
    holdout = rows(folder/"holdout-checks.csv")
    for r in rows(folder/"holdout-summary.csv"):
        group = [g for g in holdout if g["scenario"] == r["scenario"]]
        assert int(r["replicates"]) == len(group)
        assert int(r["detected"]) == sum(g["status"] == "miscalibration_detected" for g in group)
    wrong = [dict(r) for r in records if r["scenario"] != "missed_flags"]
    try:
        verify_contamination_records(wrong)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: omitted planned cell went undetected")
    print(f"{folder.name}: all {len(records)} planned requests, {len(summary)} recomputed groups, statuses, paired MCSEs and deletion control: PASS")


def least_squares_last(design, y):
    """Last coefficient of an ordinary regression, by Gaussian elimination in plain Python."""
    p = len(design[0])
    a = [[0.0]*(p+1) for _ in range(p)]
    for row, value in zip(design, y):
        for i in range(p):
            a[i][p] += row[i]*value
            for j in range(i, p):
                a[i][j] += row[i]*row[j]
    for i in range(p):
        for j in range(i):
            a[i][j] = a[j][i]
    for i in range(p):
        pivot = max(range(i, p), key=lambda k: abs(a[k][i]))
        a[i], a[pivot] = a[pivot], a[i]
        for k in range(i+1, p):
            factor = a[k][i]/a[i][i]
            for j in range(i, p+1):
                a[k][j] -= factor*a[i][j]
    beta = [0.0]*p
    for i in reversed(range(p)):
        beta[i] = (a[i][p]-sum(a[i][j]*beta[j] for j in range(i+1, p)))/a[i][i]
    return beta[-1]


def real_data_targets():
    """Recompute both target coefficients from the frozen pools, independently of R."""
    pool = [r for r in rows(ROOT/"validation/real-data/ctgov-pool.csv") if r["split"] == "evaluation"]
    design = [[1.0, r["phase"] == "early", r["phase"] == "phase2", r["phase"] == "phase3_4", r["sponsor"] == "industry",
               math.log(float(r["sites"])), float(r["start_year"])-2010, math.log(float(r["reference_n"]))] for r in pool]
    design = [[float(v) for v in row] for row in design]
    ctgov = least_squares_last(design, [math.log(float(r["duration_months"])) for r in pool]), pool
    pool = [r for r in rows(ROOT/"validation/real-data/nhanes-pool.csv") if r["split"] == "evaluation"]
    levels = ["black", "mexican_american", "other_hispanic", "other"]
    design = [[1.0, float(r["age"]), float(r["sex"] == "male")] + [float(r["race"] == v) for v in levels] +
              [float(r["height"]), float(r["reference_kg"])] for r in pool]
    nhanes = least_squares_last(design, [float(r["sbp"]) for r in pool]), pool
    return {"ctgov": ctgov, "nhanes": nhanes}


def real_data_cells():
    every = ["oracle", "naive_recorded", "gold_only", "reference_normal", "cole_mi", "raking", "contamination", "contamination_strata"]
    expected = set()
    for family, cells in [("ctgov", ["audit_25pct", "audit_10pct", "audit_50pct", "audit_05pct", "targeted_audit"]),
                          ("nhanes", ["audit_25pct", "audit_10pct", "audit_50pct", "targeted_audit"])]:
        for cell in cells:
            methods = list(every)
            if cell in ["audit_25pct", "audit_10pct", "audit_50pct"]:
                methods += ["reference_pmm", "reference_pmm_abb"]
            if family == "nhanes" and cell == "audit_25pct":
                methods.append("contamination_no_tol")
            reps = 1000 if cell == "audit_25pct" else 200
            expected |= {(family, cell, method, str(rep)) for method in methods for rep in range(1, reps+1)}
    return expected


def verify_real_data_records(records):
    actual = {(r["family"], r["scenario"], r["method"], r["replicate"]) for r in records}
    assert actual == real_data_cells() and len(actual) == len(records), "Missing, duplicate or unplanned requests"
    for r in records:
        rep = int(r["replicate"])
        offset = 5000000 if r["family"] == "nhanes" else 0
        assert int(r["data_seed"]) == 210100000+offset+rep and int(r["imputation_seed"]) == 220000000+offset+rep
        assert r["data_stream"] == f'{r["family"]}:evaluation:{rep}' and int(r["n"]) == 600


def verify_real_data_study(folder, source_root):
    for row in rows(folder/"source-hashes.csv"):
        candidate = source_root/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    manifest = json.loads((ROOT/"validation/real-data/fetch-manifest.json").read_text())
    for name, entry in [("ctgov-pool.csv", manifest["ctgov"]), ("nhanes-pool.csv", manifest["nhanes_pool"])]:
        assert hashlib.sha256((ROOT/"validation/real-data"/name).read_bytes()).hexdigest() == entry["pool_sha256"], name
    settings = json.loads((folder/"settings.json").read_text())
    assert settings["mode"] == "study" and settings["split"] == "evaluation" and settings["m"] == 20
    assert not (folder/"checkpoint.rds").exists()
    records = rows(folder/"replicates.csv")
    verify_real_data_records(records)
    summary = verify_study_bundle(folder, records)
    targets = real_data_targets()
    saved = {r["dataset"]: r for r in rows(folder/"pool-targets.csv")}
    for name, (value, pool) in targets.items():
        assert int(saved[name]["pool_rows"]) == len(pool)
        assert math.isclose(float(saved[name]["target_coefficient"]), value, rel_tol=1e-7), name
        assert all(math.isclose(float(r["truth"]), value, rel_tol=1e-7) for r in records if r["family"] == name)
    ctgov = targets["ctgov"][1]
    disagree = sum(r["recorded_n"] != r["reference_n"] for r in ctgov)
    overall = {r["dataset"]: r for r in rows(folder/"pool-audit-overall.csv")}
    assert int(overall["ctgov"]["disagreements"]) == disagree and int(overall["ctgov"]["calibration_rows"]) == len(ctgov)
    nhanes = [r for r in targets["nhanes"][1] if r["recorded_kg"] != "NA"]
    assert int(overall["nhanes"]["calibration_rows"]) == len(nhanes)
    assert int(overall["nhanes"]["disagreements"]) == sum(abs(float(r["recorded_kg"])-float(r["reference_kg"])) > 1+1e-8*max(1, float(r["reference_kg"])) for r in nhanes)
    wrong = [r for r in records if not (r["family"] == "nhanes" and r["scenario"] == "targeted_audit")]
    try:
        verify_real_data_records(wrong)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: omitted planned cell went undetected")
    print(f"{folder.name}: all {len(records)} planned requests, {len(summary)} recomputed groups, independent pool targets and audit counts, deletion control: PASS")


def verify_contamination_supplement(folder):
    scores = rows(folder/"corruption-scores.csv")
    main = rows(ROOT/"validation/contamination-study/replicates.csv")
    returned = {(r["scenario"], r["replicate"]) for r in main if r["method"] == "contamination" and r["status"] == "ok"}
    recovered = {(r["scenario"], r["replicate"]) for r in scores}
    missing = {"missed_flags", "reference_05pct", "reference_10pct", "reference_50pct", "heavy_tail_target", "nonlinear_target"}
    assert recovered == {key for key in returned if key[0] in missing}, "Recovered records do not match returned main-study fits"
    print(f"{folder.name}: {len(scores)} recovered records match every returned main-study fit in the six affected cells: PASS")


def verify_merged_audit(folder, source_root):
    """Recompute the go/no-go contrasts and verdict from replicate records."""
    for row in rows(folder/"source-hashes.csv"):
        candidate = source_root/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    records = rows(folder/"replicates.csv")
    assert len(records) == 4*400*9*2, "Missing or unplanned go/no-go requests"
    effect = math.log(.5)
    table = {}
    for r in records:
        if r["analysis"] == "pooled" and r["status"] == "ok":
            table.setdefault((r["regime"], float(r["direction"]), r["policy"]), {})[r["replicate"]] = (float(r["estimate"])-effect)**2

    def contrast(regime, direction, new, old):
        a, b = table[regime, direction, new], table[regime, direction, old]
        keys = a.keys() & b.keys()
        d = [a[k]-b[k] for k in keys]
        return statistics.mean(a[k] for k in keys)/statistics.mean(b[k] for k in keys), statistics.mean(d), statistics.mean(d)/mean_se(d)
    ratio, _, z = contrast("moderate", 1.0, "full_audit", "recorded_only")
    matter = ratio <= .9 and z <= -2.5
    r1, _, z1 = contrast("moderate", 1.0, "transmitter", "reactive")
    r2, _, z2 = contrast("moderate", 1.0, "transmitter", "random")
    network = r1 <= .9 and z1 <= -2.5 and r2 < 1 and z2 <= -2.5
    _, full_gap, z3 = contrast("moderate", 1.0, "transmitter", "undirected")
    _, half_gap, _ = contrast("moderate", .5, "transmitter", "undirected")
    direction = z3 <= -2 and full_gap < half_gap
    verdict = "NO-GO" if not (matter and network) else ("GO" if direction else "PARTIAL")
    saved = rows(folder/"verdict.csv")[0]
    assert saved["verdict"] == verdict and verdict == "NO-GO", (saved, verdict)
    zero = [table["moderate", 0.0, "transmitter"][k] == table["moderate", 0.0, "undirected"][k] for k in table["moderate", 0.0, "transmitter"]]
    assert all(zero), "Direction 0 must make the two leverage policies identical"
    print(f"{folder.name}: {len(records)} requests, independently recomputed {verdict}, direction-0 identity control: PASS")


def solve_normal(matrix, targets):
    """Inverse-free solve of a symmetric system by Gauss-Jordan; returns solution rows for each right side."""
    p = len(matrix)
    a = [list(matrix[i])+[t[i] for t in targets] for i in range(p)]
    for i in range(p):
        pivot = max(range(i, p), key=lambda k: abs(a[k][i]))
        a[i], a[pivot] = a[pivot], a[i]
        lead = a[i][i]
        a[i] = [v/lead for v in a[i]]
        for k in range(p):
            if k != i and a[k][i] != 0:
                factor = a[k][i]
                a[k] = [x-factor*y for x, y in zip(a[k], a[i])]
    return [[a[i][p+t] for i in range(p)] for t in range(len(targets))]


def cross(design, weights=None):
    p = len(design[0])
    out = [[0.0]*p for _ in range(p)]
    for index, row in enumerate(design):
        w = 1.0 if weights is None else weights[index]
        for i in range(p):
            value = row[i]*w
            if value != 0:
                target = out[i]
                for j in range(i, p):
                    target[j] += value*row[j]
    for i in range(p):
        for j in range(i):
            out[i][j] = out[j][i]
    return out


def ols(design, y):
    """Coefficients, residuals and residual variance of an ordinary regression, in plain Python."""
    p = len(design[0])
    right = [sum(row[i]*v for row, v in zip(design, y)) for i in range(p)]
    beta = solve_normal(cross(design), [right])[0]
    residuals = [v-sum(b*x for b, x in zip(beta, row)) for row, v in zip(design, y)]
    return beta, residuals, sum(e*e for e in residuals)/(len(y)-p)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def audit_imputation(rows, plain, full, gold, observed):
    """Imputation mean and variance for every record, from reference rows. full is None when discarding."""
    n = len(plain)
    mu, tau2 = [0.0]*n, [0.0]*n
    if full is None or any(o is None for o in observed):
        beta, _, variance = ols([plain[i] for i in rows], [gold[i] for i in rows])
        for i in range(n):
            mu[i], tau2[i] = dot(beta, plain[i]), variance
    if full is not None:
        use = [i for i in rows if observed[i] is not None]
        if all(abs(observed[i]-gold[i]) <= 1e-8*max(1.0, abs(gold[i])) for i in use):
            for i in range(n):
                if observed[i] is not None:
                    mu[i], tau2[i] = observed[i], 0.0
        else:
            beta, _, variance = ols([full[i] for i in use], [gold[i] for i in use])
            for i in range(n):
                if observed[i] is not None:
                    mu[i], tau2[i] = dot(beta, full[i]), variance
    return mu, tau2


def audit_k(covariates, outcome, known, gold, mu, tau2):
    """k under the working models for outcome ~ value + covariates, with the term being the value."""
    n = len(outcome)
    y = [gold[i] if known[i] else mu[i] for i in range(n)]
    design = [[1.0, y[i]]+covariates[i] for i in range(n)]
    p = len(design[0])
    matrix = cross(design)
    matrix[1][1] += sum(tau2[i] for i in range(n) if not known[i])
    right = [sum(row[i]*o for row, o in zip(design, outcome)) for i in range(p)]
    unit = [1.0 if i == 1 else 0.0 for i in range(p)]
    beta, row = solve_normal(matrix, [right, unit])
    sse = sum((o-dot(beta, d))**2 for d, o in zip(design, outcome))+sum(tau2[i]*beta[1]**2 for i in range(n) if not known[i])
    s2 = sse/(n-p)
    weight = [n*v for v in row]
    total = 0.0
    a1, e1 = weight[1], -beta[1]
    for i in range(n):
        base = [1.0, mu[i]]+covariates[i]
        u, v = dot(base, weight), outcome[i]-dot(base, beta)
        total += tau2[i]*(u*e1+v*a1)**2+2*(a1*e1)**2*tau2[i]**2
    return total/n/(s2*weight[1])


def audit_classical(covariates, outcome, gold, mu, tau2, plain, full, observed):
    """The empirical ratio on data where every row has a reference, with leave-one-out imputation means."""
    n = len(outcome)
    design = [[1.0, gold[i]]+covariates[i] for i in range(n)]
    p = len(design[0])
    matrix = cross(design)
    right = [sum(row[i]*o for row, o in zip(design, outcome)) for i in range(p)]
    unit = [1.0 if i == 1 else 0.0 for i in range(p)]
    beta, row = solve_normal(matrix, [right, unit])
    weight = [n*v for v in row]
    a1, e1 = weight[1], -beta[1]
    loo = [None]*n
    for source, rows in ([(full, [i for i in range(n) if observed[i] is not None]),
                          (plain, [i for i in range(n) if observed[i] is None])] if full is not None else [(plain, list(range(n)))]):
        if not rows:
            continue
        fit_rows = rows if source is full or full is None else list(range(n))
        x = [source[i] for i in fit_rows]
        q = len(x[0])
        fitted = solve_normal(cross(x), [[sum(r[j]*gold[i] for r, i in zip(x, fit_rows)) for j in range(q)]]+[[float(j == k) for j in range(q)] for k in range(q)])
        coefficients, inverse = fitted[0], fitted[1:]
        for i in rows:
            r = source[i]
            h = sum(r[a]*dot(inverse[a], r) for a in range(q))
            e = gold[i]-dot(coefficients, r)
            loo[i] = gold[i]-e/max(1-h, 1e-8)
    top, bottom = 0.0, 0.0
    for i in range(n):
        psi = dot(design[i], weight)*(outcome[i]-dot(design[i], beta))
        base = [1.0, loo[i]]+covariates[i]
        expected = dot(base, weight)*(outcome[i]-dot(base, beta))+a1*e1*tau2[i]
        top += (psi-expected)**2
        bottom += psi*psi
    return top/bottom


def audit_pool_inputs(name, split="evaluation"):
    pool = [r for r in rows(ROOT/f"validation/real-data/{name}-pool.csv") if r["split"] == split]
    if name == "ctgov":
        covariates = [[float(r["phase"] == "early"), float(r["phase"] == "phase2"), float(r["phase"] == "phase3_4"), float(r["sponsor"] == "industry"),
                       math.log(float(r["sites"])), float(r["start_year"])-2010] for r in pool]
        outcome = [math.log(float(r["duration_months"])) for r in pool]
        gold = [math.log(float(r["reference_n"])) for r in pool]
        recorded = [math.log(float(r["recorded_n"])) for r in pool]
    else:
        levels = ["black", "mexican_american", "other_hispanic", "other"]
        covariates = [[float(r["age"]), float(r["sex"] == "male")]+[float(r["race"] == v) for v in levels]+[float(r["height"])] for r in pool]
        outcome = [float(r["sbp"]) for r in pool]
        gold = [float(r["reference_kg"]) for r in pool]
        recorded = [None if r["recorded_kg"] == "NA" else float(r["recorded_kg"]) for r in pool]
    return covariates, outcome, gold, recorded


def audit_values(covariates, outcome, gold, recorded, known):
    n = len(outcome)
    plain = [[1.0]+covariates[i]+[outcome[i]] for i in range(n)]
    full = [plain[i]+[recorded[i] if recorded[i] is not None else 0.0] for i in range(n)]
    rows_known = [i for i in range(n) if known[i]]
    mu, tau2 = audit_imputation(rows_known, plain, full, gold, recorded)
    mu0, tau02 = audit_imputation(rows_known, plain, None, gold, recorded)
    filled = [g if g is not None else 0.0 for g in gold]
    return {"k": audit_k(covariates, outcome, known, filled, mu, tau2),
            "discard_k": audit_k(covariates, outcome, known, filled, mu0, tau02),
            "mu": mu, "tau2": tau2, "mu0": mu0, "tau02": tau02, "plain": plain, "full": full}


def order_interval(values):
    x = sorted(values)
    n = len(x)
    low = max(1, math.floor(n/2-1.96*math.sqrt(n)/2))
    high = min(n, math.ceil(n/2+1+1.96*math.sqrt(n)/2))
    return x[low-1], x[high-1]


def quantile7(values, p):
    x = sorted(values)
    h = (len(x)-1)*p
    low = math.floor(h)
    return x[low]+(h-low)*(x[min(low+1, len(x)-1)]-x[low])


def audit_verdict(summary):
    def holds(size):
        for name in ["ctgov", "nhanes"]:
            group = [r for r in summary if r["cell"] == name and int(r["pilot_records"]) == size]
            if len(group) != 3 or any(float(r["median_error"]) > .10 for r in group):
                return False
        return True
    return "GO" if holds(60) else "PARTIAL" if holds(150) else "NO-GO"


def audit_expected_keys(settings):
    fractions = [size/settings["n"] for size in settings["target_records"]]
    reps = {"ctgov": settings["real_reps"], "nhanes": settings["real_reps"],
            "missed_flags": settings["synthetic_reps"], "error_gross_and_small": settings["synthetic_reps"]}
    return fractions, reps


def verify_audit_records(targets, plans, settings):
    fractions, reps = audit_expected_keys(settings)
    expected_targets = {(cell, rep, round(f, 6)) for cell, count in reps.items() for rep in range(1, count+1) for f in fractions}
    actual = {(r["cell"], int(r["replicate"]), round(float(r["target_fraction"]), 6)) for r in targets}
    assert actual == expected_targets and len(targets) == len(expected_targets), "Missing, duplicate or unplanned target fits"
    expected_plans = {(key[0], key[1], size, key[2]) for key in expected_targets for size in settings["pilot_records"]}
    actual = {(r["cell"], int(r["replicate"]), int(r["pilot_records"]), round(float(r["target_fraction"]), 6)) for r in plans}
    assert actual == expected_plans and len(plans) == len(expected_plans), "Missing, duplicate or unplanned plans"
    for r in plans:
        base = settings["plan_10_base"] if int(r["pilot_records"]) == 60 else settings["plan_25_base"]
        assert int(r["plan_seed"]) == base+settings["offsets"][r["cell"]]+4*int(r["replicate"])


def recompute_audit_summary(targets, plans):
    """Recomputes every summary field from replicate rows, independently of R."""
    actual = {(r["cell"], r["replicate"], round(float(r["target_fraction"]), 6)): r for r in targets}
    groups = {}
    for r in plans:
        fraction = round(float(r["target_fraction"]), 6)
        if abs(fraction-float(r["pilot_fraction"])) > 1e-9:
            groups.setdefault((r["cell"], int(r["pilot_records"]), fraction), []).append((r, actual[(r["cell"], r["replicate"], fraction)]))
    out = {}
    for key, pairs in groups.items():
        usable = [(p, t) for p, t in pairs if t["status"] == "ok"]
        errors, signed, widths, ratios = [], [], [], []
        interval, predictive, predicted_ratio = [], [], []
        for p, t in usable:
            width = float(t["width"])
            widths.append(width)
            ratio = width/float(t["full_width"])
            ratios.append(ratio)
            if p["status"] == "ok":
                gap = (float(p["predicted_width"])-width)/width
                errors.append(abs(gap))
                signed.append(gap)
                interval.append(float(p["relative_lower"]) <= ratio <= float(p["relative_upper"]))
                predictive.append(float(p["predictive_lower"]) <= ratio <= float(p["predictive_upper"]))
                predicted_ratio.append(float(p["relative_width"]))
            else:
                errors.append(math.inf)
        center = statistics.median(widths)
        retest = [abs(float(t["retest_width"])-float(t["width"]))/float(t["width"]) for _, t in usable if t["retest_width"] != "NA"]
        low, high = order_interval(errors)
        out[key] = {"requested": len(pairs), "target_failures": len(pairs)-len(usable),
                    "plan_failures": sum(p["status"] != "ok" for p, _ in usable), "compared": len(usable),
                    "median_error": statistics.median(errors), "median_error_lower": low, "median_error_upper": high,
                    "p90_error": quantile7(errors, .9), "median_signed_error": statistics.median(signed),
                    "interval_coverage": statistics.mean(interval), "predictive_coverage": statistics.mean(predictive),
                    "floor_median_error": statistics.median(abs(center-w)/w for w in widths),
                    "retest_median_error": statistics.median(retest) if retest else None,
                    "mean_actual_ratio": statistics.mean(ratios), "median_actual_ratio": statistics.median(ratios),
                    "mean_predicted_ratio": statistics.mean(predicted_ratio)}
    return out


def check_audit_summary(summary, calculated):
    assert len(summary) == len(calculated), "Summary has missing or unplanned cells"
    for row in summary:
        key = (row["cell"], int(row["pilot_records"]), round(float(row["target_fraction"]), 6))
        for name, value in calculated[key].items():
            if value is None:
                assert row[name] == "NA", (key, name)
            elif math.isinf(value):
                assert row[name] == "Inf", (key, name)
            else:
                assert math.isclose(float(row[name]), value, rel_tol=1e-9, abs_tol=1e-11), (key, name, row[name], value)


def verify_audit_picks(targets, plans, saved, ratio):
    truth = {}
    for r in targets:
        key = (r["cell"], r["replicate"])
        truth.setdefault(key, [])
        if r["status"] == "ok" and float(r["width"])/float(r["full_width"]) <= ratio:
            truth[key].append(float(r["target_fraction"]))
    seen, groups = set(), {}
    for r in plans:
        key = (r["cell"], r["replicate"], r["pilot_records"])
        if key in seen:
            continue
        seen.add(key)
        true_pick = min(truth[(r["cell"], r["replicate"])], default=None)
        pick = None if r["smallest_grid_fraction"] == "NA" else float(r["smallest_grid_fraction"])
        same = r["status"] == "ok" and ((true_pick is None and pick is None) or
                                        (true_pick is not None and pick is not None and abs(true_pick-pick) < 1e-9))
        groups.setdefault((r["cell"], int(r["pilot_records"])), []).append(same)
    assert len(saved) == len(groups)
    for row in saved:
        same = groups[(row["cell"], int(row["pilot_records"]))]
        assert int(row["cohorts"]) == len(same)
        assert math.isclose(float(row["same_pick"]), statistics.mean(same), rel_tol=1e-11), row["cell"]


def verify_audit_k_summary(plans, saved):
    seen, groups = set(), {}
    for r in plans:
        key = (r["cell"], r["replicate"], r["pilot_records"])
        if key not in seen:
            seen.add(key)
            if r["status"] == "ok":
                groups.setdefault((r["cell"], int(r["pilot_records"])), []).append(r)
    for row in saved:
        group = groups[(row["cell"], int(row["pilot_records"]))]
        k = [float(r["k"]) for r in group]
        assert int(row["plans_returned"]) == len(group)
        for name, value in [("k_mean", statistics.mean(k)), ("k_median", statistics.median(k)), ("k_q05", quantile7(k, .05)), ("k_q95", quantile7(k, .95))]:
            assert math.isclose(float(row[name]), value, rel_tol=1e-9), (row["cell"], name)
        if row["pool_k"] != "NA":
            pool = float(row["pool_k"])
            cover = statistics.mean(float(r["k_lower"]) <= pool <= float(r["k_upper"]) for r in group)
            assert math.isclose(float(row["interval_covers_pool"]), cover, rel_tol=1e-9), row["cell"]
        assert math.isclose(float(row["flag_gross_errors"]), statistics.mean("gross_errors" in r["flags"] for r in group), abs_tol=1e-12)


def verify_audit_vectors(folder):
    """Recompute k from the browser test vectors in plain Python, which also tests the written description."""
    data = json.loads((folder/"planner-data.json").read_text())
    vectors = data["js_test_vectors"]
    assert len(vectors) >= 3 and data["k_calculation"]
    for vector in vectors:
        recorded = vector["recorded"]
        gold = vector["reference"]
        known = [g is not None for g in gold]
        values = audit_values(vector["covariates"], vector["outcome"], gold, recorded, known)
        assert math.isclose(values["k"], vector["expected"]["k"], rel_tol=1e-8, abs_tol=1e-12), vector["name"]
        assert math.isclose(values["discard_k"], vector["expected"]["discard_k"], rel_tol=1e-8), vector["name"]
        for row in vector["expected"]["table"]:
            f = row["fraction"]
            assert math.isclose(row["relative_width"], math.sqrt(1+values["k"]*(1-f)/f), rel_tol=1e-8), vector["name"]
            assert math.isclose(row["discard_relative_width"], math.sqrt(1+values["discard_k"]*(1-f)/f), rel_tol=1e-8), vector["name"]
            assert row["records"] == min(len(gold), math.ceil(f*len(gold)-1e-8)), vector["name"]
    return data, len(vectors)


def verify_audit_plan_study(folder, source_root):
    for row in rows(folder/"source-hashes.csv"):
        candidate = source_root/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    manifest = json.loads((ROOT/"validation/real-data/fetch-manifest.json").read_text())
    for name, entry in [("ctgov-pool.csv", manifest["ctgov"]), ("nhanes-pool.csv", manifest["nhanes_pool"])]:
        assert hashlib.sha256((ROOT/"validation/real-data"/name).read_bytes()).hexdigest() == entry["pool_sha256"], name
    settings = json.loads((folder/"settings.json").read_text())
    assert settings["mode"] == "study" and settings["split"] == "evaluation" and settings["m"] == 20 and settings["n"] == 600
    assert settings["target_records"] == [30, 60, 150, 300] and settings["pilot_records"] == [60, 150]
    assert settings["cohort_base"] > 505100000, "Seed streams must be disjoint from earlier studies"
    assert not (folder/"checkpoint.rds").exists()
    targets = rows(folder/"targets.csv")
    plans = rows(folder/"plans.csv")
    verify_audit_records(targets, plans, settings)
    summary = rows(folder/"summary.csv")
    calculated = recompute_audit_summary(targets, plans)
    check_audit_summary(summary, calculated)
    verdict = rows(folder/"verdict.csv")[0]["verdict"]
    assert audit_verdict(summary) == verdict
    verify_audit_picks(targets, plans, rows(folder/"fraction-picks.csv"), settings["target_ratio"])
    verify_audit_k_summary(plans, rows(folder/"k-summary.csv"))
    pool_saved = {r["dataset"]: r for r in rows(folder/"pool-k.csv")}
    for name in ["ctgov", "nhanes"]:
        covariates, outcome, gold, recorded = audit_pool_inputs(name)
        known = [True]*len(gold)
        values = audit_values(covariates, outcome, gold, recorded, known)
        saved = pool_saved[name]
        assert int(saved["pool_rows"]) == len(gold)
        assert math.isclose(float(saved["k"]), values["k"], rel_tol=1e-7), name
        assert math.isclose(float(saved["discard_k"]), values["discard_k"], rel_tol=1e-7), name
        classical = audit_classical(covariates, outcome, gold, values["mu"], values["tau2"], values["plain"], values["full"], recorded)
        discard = audit_classical(covariates, outcome, gold, values["mu0"], values["tau02"], values["plain"], None, recorded)
        assert math.isclose(float(saved["classical"]), classical, rel_tol=1e-6), name
        assert math.isclose(float(saved["discard_classical"]), discard, rel_tol=1e-6), name
    wrong = [dict(r) for r in summary]
    try:
        check_audit_summary([dict(wrong[0], median_error=str(float(wrong[0]["median_error"])+.01))]+wrong[1:], calculated)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: altered median error went undetected")
    try:
        verify_audit_records(targets[:-1], plans, settings)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: omitted target fit went undetected")
    flipped = [dict(r, median_error="0.05") for r in summary]
    assert audit_verdict(flipped) == "GO"
    flipped = [dict(r, median_error="0.05" if int(r["pilot_records"]) == 150 else "0.5") for r in summary]
    assert audit_verdict(flipped) == "PARTIAL"
    assert audit_verdict([dict(r, median_error="0.5") for r in summary]) == "NO-GO"
    data, count = verify_audit_vectors(folder)
    assert data["verdict"] == verdict
    for name in ["ctgov", "nhanes"]:
        block = data["datasets"][name]
        assert math.isclose(block["pool_k"], float(pool_saved[name]["k"]), rel_tol=1e-12), name
    frozen = {(r["family"], r["scenario"], r["method"]): float(r["mean_width"]) for r in rows(ROOT/"validation/real-data-study/summary.csv")}
    for name in ["ctgov", "nhanes"]:
        block = data["datasets"][name]["frozen_study"]
        for fraction, scenario in [(.05, "audit_05pct"), (.10, "audit_10pct"), (.25, "audit_25pct"), (.50, "audit_50pct")]:
            entry = [e for e in block["by_fraction"] if abs(e["fraction"]-fraction) < 1e-9]
            if (name, scenario, "cole_mi") in frozen:
                assert math.isclose(entry[0]["linear_mi_mean_width"], frozen[(name, scenario, "cole_mi")], rel_tol=1e-12)
                assert math.isclose(entry[0]["reference_normal_mean_width"], frozen[(name, scenario, "reference_normal")], rel_tol=1e-12)
                assert math.isclose(entry[0]["oracle_mean_width"], frozen[(name, scenario, "oracle")], rel_tol=1e-12)
            else:
                assert not entry
    frozen_rows = {(r["family"], r["scenario"], r["method"]): r for r in rows(ROOT/"validation/real-data-study/summary.csv")}
    scenarios = {.05: "audit_05pct", .10: "audit_10pct", .25: "audit_25pct", .50: "audit_50pct"}
    comparison = rows(folder/"frozen-comparison.csv")
    assert len(comparison) == 7
    for row in comparison:
        fraction = float(row["fraction"])
        scenario = scenarios[fraction]
        odds = (1-fraction)/fraction
        get = lambda method, field: float(frozen_rows[(row["dataset"], scenario, method)][field])
        k_width = ((get("cole_mi", "mean_width")/get("oracle", "mean_width"))**2-1)/odds
        k_rmse = ((get("cole_mi", "rmse")/get("oracle", "rmse"))**2-1)/odds
        ratio = get("reference_normal", "mean_width")/get("oracle", "mean_width")
        planned = math.sqrt(1+float(pool_saved[row["dataset"]]["discard_k"])*odds)
        for name, value in [("k_implied_by_width", k_width), ("k_implied_by_rmse", k_rmse), ("discard_understatement", 1-planned/ratio)]:
            assert math.isclose(float(row[name]), value, rel_tol=1e-9), (row["dataset"], fraction, name)
    print(f"{folder.name}: all {len(targets)} target fits and {len(plans)} plan rows, {len(summary)} recomputed cells, verdict {verdict}, "
          f"independent pool k and {count} browser test vectors recomputed: PASS")


def sponsor_trials(split):
    mapping = {r["nct_id"]: r for r in rows(ROOT/"validation/sponsor-clustering/sponsor-map.csv")}
    pool = rows(ROOT/"validation/real-data/ctgov-pool.csv")
    assert set(mapping) == {r["nct_id"] for r in pool}
    trials = []
    for r in pool:
        named = mapping[r["nct_id"]]
        assert named["sponsor_key"] == named["lead_sponsor_name"].strip().casefold()
        if r["split"] == split:
            year, recorded, reference = int(r["start_year"]), int(r["recorded_n"]), int(r["reference_n"])
            trials.append(dict(id=r["nct_id"], key=named["sponsor_key"], phase=r["phase"], cls=r["sponsor"],
                               band=0 if year <= 2013 else 1 if year <= 2015 else 2 if year <= 2017 else 3, n=recorded,
                               D=int(recorded != reference), D2=int(abs(math.log(recorded/reference)) > .25)))
    return trials


def sponsor_counts(trials):
    counts = {}
    for t in trials:
        counts[t["key"]] = counts.get(t["key"], 0)+1
    return counts


def sponsor_strata(trials, outcome):
    counts = sponsor_counts(trials)
    chosen = [t for t in trials if counts[t["key"]] >= 3]
    cuts = [quantile7([t["n"] for t in chosen], p) for p in (.25, .5, .75)]
    strata = {}
    for t in chosen:
        strata.setdefault((t["phase"], t["cls"], t["band"], sum(t["n"] > c for c in cuts)), []).append(t)
    expected, observed = {}, {}
    for members in strata.values():
        rate = sum(t[outcome] for t in members)/len(members)
        for t in members:
            expected[t["key"]] = expected.get(t["key"], 0)+rate
            observed[t["key"]] = observed.get(t["key"], 0)+t[outcome]
    return chosen, strata, expected, observed, cuts


def sponsor_terms(expected, observed):
    return {k: (observed[k]-e)**2/e if e > 0 else 0.0 for k, e in expected.items()}


def sponsor_shuffled(strata, expected, outcome, count, seed):
    stream = random.Random(seed)
    results = []
    for _ in range(count):
        observed = dict.fromkeys(expected, 0)
        for members in strata.values():
            values = [t[outcome] for t in members]
            stream.shuffle(values)
            for t, value in zip(members, values):
                observed[t["key"]] += value
        results.append(sum(sponsor_terms(expected, observed).values()))
    return results


def sponsor_features(t, history=None):
    x = [1.0, float(t["phase"] == "phase2"), float(t["phase"] == "phase3_4"), float(t["phase"] == "not_applicable"),
         float(t["cls"] == "industry"), float(t["band"] == 1), float(t["band"] == 2), float(t["band"] == 3), math.log(t["n"])]
    return x if history is None else x+[history]


def sponsor_fit(xs, ys):
    p = len(xs[0])
    beta = [0.0]*p
    for _ in range(100):
        matrix = [[0.0]*p for _ in range(p)]
        gradient = [0.0]*p
        for x, y in zip(xs, ys):
            mu = 1/(1+math.exp(-dot(x, beta)))
            weight = mu*(1-mu)
            for i in range(p):
                gradient[i] += x[i]*(y-mu)
                row = matrix[i]
                for j in range(i, p):
                    row[j] += weight*x[i]*x[j]
        for i in range(p):
            for j in range(i):
                matrix[i][j] = matrix[j][i]
        delta = solve_normal(matrix, [gradient])[0]
        beta = [b+d for b, d in zip(beta, delta)]
        if max(abs(d) for d in delta) < 1e-10:
            return beta
    raise AssertionError("Plain-Python logistic fit did not converge")


def sponsor_score(x, beta):
    return 1/(1+math.exp(-dot(x, beta)))


def sponsor_auc(scores, ys):
    order = sorted(range(len(scores)), key=scores.__getitem__)
    ranks = [0.0]*len(scores)
    i = 0
    while i < len(order):
        j = i
        while j+1 < len(order) and scores[order[j+1]] == scores[order[i]]:
            j += 1
        for k in range(i, j+1):
            ranks[order[k]] = (i+j)/2+1
        i = j+1
    positives = sum(ys)
    negatives = len(ys)-positives
    return (sum(r for r, y in zip(ranks, ys) if y)-positives*(positives+1)/2)/(positives*negatives)


def sponsor_brier(scores, ys):
    return sum((s-y)**2 for s, y in zip(scores, ys))/len(ys)


def sponsor_audit(population, outcome):
    """population holds (trial, baseline score) in nct_id order; returns per-seed audit counts."""
    budget = (len(population)+5)//10
    first_size = budget//2
    second_size = budget-first_size
    ranked = sorted(range(len(population)), key=lambda i: (-population[i][1], hashlib.sha256(population[i][0]["id"].encode()).hexdigest()))
    ys = [p[0][outcome] for p in population]
    keys = [p[0]["key"] for p in population]
    risk = sum(ys[i] for i in ranked[:budget])
    seeds = []
    for seed in range(200):
        stream = random.Random(20261104+seed)
        drawn = sorted((stream.random(), i) for i in range(len(population)))[:budget]
        picked = [i for _, i in drawn]
        first = picked[:first_size]
        flagged = {keys[i] for i in first if ys[i]}
        seen = set(first)
        second = [i for i in ranked if i not in seen and keys[i] in flagged][:second_size]
        shortfall = second_size-len(second)
        second_found = sum(ys[i] for i in second)
        if shortfall:
            taken = seen | set(second)
            second_found += sum(ys[i] for i in [i for i in ranked if i not in taken][:shortfall])
        seeds.append(dict(found_random=sum(ys[i] for i in picked), found_adaptive=sum(ys[i] for i in first)+second_found,
                          found_second_stage=second_found, flagged_sponsors=len(flagged), second_stage_shortfall=shortfall))
    return dict(budget=budget, second_size=second_size, risk=risk, seeds=seeds, ranked=ranked, ys=ys, keys=keys)


def sponsor_audit_metrics(audit, sponsors, total):
    budget, seeds = audit["budget"], audit["seeds"]
    random_yield = 100*statistics.mean(s["found_random"] for s in seeds)/budget
    risk_yield = 100*audit["risk"]/budget
    adaptive_yield = 100*statistics.mean(s["found_adaptive"] for s in seeds)/budget
    return dict(population_trials=total, sponsors=sponsors, budget=budget, discrepancies=sum(audit["ys"]),
                yield_random=random_yield, yield_risk=risk_yield, yield_adaptive=adaptive_yield,
                yield_second_stage=100*statistics.mean(s["found_second_stage"] for s in seeds)/audit["second_size"],
                mean_flagged_sponsors=statistics.mean(s["flagged_sponsors"] for s in seeds),
                mean_second_stage_shortfall=statistics.mean(s["second_stage_shortfall"] for s in seeds),
                yield_adaptive_seed_se=100*statistics.stdev(s["found_adaptive"] for s in seeds)/budget/math.sqrt(200),
                adaptive_to_risk=adaptive_yield/risk_yield, adaptive_to_random=adaptive_yield/random_yield,
                risk_to_random=risk_yield/random_yield)


def sponsor_verdict(p, auc_difference, auc_lower, ratio, ratio_lower):
    if p >= .01:
        return "NO-CLUSTERING"
    if auc_difference >= .02 and auc_lower > 0 and ratio >= 1.2 and ratio_lower > 1:
        return "CLUSTERING-SUPPORTED"
    return "WEAK"


def sponsor_recompute(folder, trials):
    """Rebuild every summary number from the frozen pool, the sponsor map and the saved permutation and bootstrap draws."""
    calc = {}
    limits = {}
    counts = sponsor_counts(trials)
    for metric, value in [("trials_in_split", len(trials)), ("sponsors_in_split", len(counts)),
                          ("sponsors_at_least_3", sum(v >= 3 for v in counts.values())),
                          ("sponsors_at_least_4", sum(v >= 4 for v in counts.values()))]:
        calc["counts", "both", "all", metric] = value
    saved_sponsors = {r["sponsor_key"]: r for r in rows(folder/"sponsors.csv")}
    assert set(saved_sponsors) == {k for k, v in counts.items() if v >= 3}
    null = rows(folder/"test1-null.csv")
    assert len(null) == 2000
    eligible = sorted((k for k, v in counts.items() if v >= 4), key=lambda k: hashlib.sha256(k.encode()).hexdigest())
    first = (len(eligible)+1)//2
    train, test = eligible[:first], eligible[first:]
    role = {**{k: "train" for k in train}, **{k: "test" for k in test}}
    stream = random.Random(20261102)
    draws = {}
    for t in trials:
        if t["key"] in role:
            draws.setdefault(t["key"], []).append((stream.random(), t["id"]))
    revealed = {i for values in draws.values() for _, i in sorted(values)[:len(values)//2]}
    for k, r in saved_sponsors.items():
        assert r["role"] == role.get(k, "not_eligible") and int(r["trials"]) == counts[k]
        assert int(r["revealed"]) == sum(i in revealed for _, i in draws.get(k, [])), k
    stream = random.Random(20261103)
    boot_draws = [[int(stream.random()*len(test)) for _ in range(len(test))] for _ in range(25)]
    bootstrap = rows(folder/"bootstrap.csv")
    assert len(bootstrap) == 4000
    coefficients = rows(folder/"coefficients.csv")
    predictions = {r["nct_id"]: r for r in rows(folder/"predictions.csv")}
    seed_rows = rows(folder/"test3-seeds.csv")
    tested = set(test)
    trained = set(train)
    verdict_inputs = {}
    for outcome in ["D", "D2"]:
        calc["counts", outcome, "all", "discrepancies_in_split"] = sum(t[outcome] for t in trials)
        chosen, strata, expected, observed, cuts = sponsor_strata(trials, outcome)
        terms = sponsor_terms(expected, observed)
        kinds = {t["key"]: t["cls"] for t in chosen}
        for k, r in saved_sponsors.items():
            assert int(r[f"observed_{outcome}"]) == observed[k]
            close(r[f"expected_{outcome}"], expected[k], (k, "expected"))
            close(r[f"term_{outcome}"], terms[k], (k, "term"))
        calc["test1", "both", "all", "nonempty_strata"] = len(strata)
        for label, value in zip(["q1", "q2", "q3"], cuts):
            calc["test1", "both", "all", f"size_cut_{label}"] = value
        for subgroup, column in [("all", "statistic"), ("industry", "industry"), ("other", "other")]:
            keys = [k for k in terms if subgroup == "all" or kinds[k] == subgroup]
            part = [t for t in chosen if subgroup == "all" or t["cls"] == subgroup]
            statistic = sum(terms[k] for k in keys)
            values = [float(r[f"{column}_{outcome}"]) for r in null]
            where = ("test1", outcome, subgroup)
            calc[where+("sponsors",)] = len(keys)
            calc[where+("trials",)] = len(part)
            calc[where+("discrepancies",)] = sum(t[outcome] for t in part)
            calc[where+("statistic",)] = statistic
            calc[where+("null_mean",)] = statistics.fmean(values)
            limits[where+("null_mean",)] = (quantile7(values, .025), quantile7(values, .975))
            calc[where+("p_value",)] = (1+sum(v >= statistic for v in values))/2001
            calc[where+("observed_to_null_ratio",)] = statistic/statistics.fmean(values)
        hits, seen = {}, {}
        for t in trials:
            if t["id"] in revealed:
                hits[t["key"]] = hits.get(t["key"], 0)+t[outcome]
                seen[t["key"]] = seen.get(t["key"], 0)+1
        history = {k: math.log((hits[k]+.5)/(seen[k]-hits[k]+.5)) for k in seen}
        fit_rows = [t for t in trials if t["key"] in trained and t["id"] not in revealed]
        ys = [t[outcome] for t in fit_rows]
        beta_base = sponsor_fit([sponsor_features(t) for t in fit_rows], ys)
        beta_hist = sponsor_fit([sponsor_features(t, history[t["key"]]) for t in fit_rows], ys)
        for model, beta in [("baseline", beta_base), ("history", beta_hist)]:
            saved = [r for r in coefficients if r["outcome"] == outcome and r["model"] == model]
            assert len(saved) == len(beta)
            for r, value in zip(saved, beta):
                assert math.isclose(float(r["estimate"]), value, rel_tol=1e-7, abs_tol=1e-8), (outcome, model, r["term"])
        scored = [t for t in trials if t["key"] in tested and t["id"] not in revealed]
        assert len(scored) == len(predictions)
        base = [sponsor_score(sponsor_features(t), beta_base) for t in scored]
        hist = [sponsor_score(sponsor_features(t, history[t["key"]]), beta_hist) for t in scored]
        for t, b, h in zip(scored, base, hist):
            r = predictions[t["id"]]
            assert int(r[outcome]) == t[outcome]
            assert math.isclose(float(r[f"baseline_{outcome}"]), b, rel_tol=1e-8) and math.isclose(float(r[f"history_{outcome}"]), h, rel_tol=1e-8)
        calc["test2", outcome, "all", "train_sponsors"] = len(train)
        calc["test2", outcome, "all", "test_sponsors"] = len(test)
        calc["test2", outcome, "all", "train_trials_fit"] = len(fit_rows)
        calc["test2", outcome, "all", "train_discrepancies_fit"] = sum(ys)
        calc["test2", outcome, "all", "revealed_trials"] = len(revealed)
        boot = [r for r in bootstrap if r["outcome"] == outcome]
        assert [int(r["replicate"]) for r in boot] == list(range(1, 2001))
        for subgroup in ["all", "industry", "other"]:
            take = [i for i, t in enumerate(scored) if subgroup == "all" or t["cls"] == subgroup]
            y_part = [scored[i][outcome] for i in take]
            base_part, hist_part = [base[i] for i in take], [hist[i] for i in take]
            where = ("test2", outcome, subgroup)
            calc[where+("test_trials_scored",)] = len(take)
            calc[where+("test_discrepancies",)] = sum(y_part)
            auc_b, auc_h = sponsor_auc(base_part, y_part), sponsor_auc(hist_part, y_part)
            brier_b, brier_h = sponsor_brier(base_part, y_part), sponsor_brier(hist_part, y_part)
            for metric, value in [("auc_baseline", auc_b), ("auc_history", auc_h), ("auc_difference", auc_h-auc_b),
                                  ("brier_baseline", brier_b), ("brier_history", brier_h), ("brier_difference", brier_h-brier_b)]:
                calc[where+(metric,)] = value
                column = metric if subgroup == "all" else f"{metric}_{subgroup}"
                if column in boot[0]:
                    values = [float(r[column]) for r in boot if math.isfinite(float(r[column]))]
                    limits[where+(metric,)] = (quantile7(values, .025), quantile7(values, .975))
        score = {t["id"]: sponsor_score(sponsor_features(t), beta_base) for t in trials if t["key"] in tested}
        population = [(t, score[t["id"]]) for t in trials if t["key"] in tested]
        audit = sponsor_audit(population, outcome)
        metrics = sponsor_audit_metrics(audit, len(test), len(population))
        for metric, value in metrics.items():
            calc["test3", outcome, "all", metric] = value
        for metric in ["yield_random", "yield_risk", "yield_adaptive", "adaptive_to_risk", "adaptive_to_random", "risk_to_random"]:
            values = [float(r[metric]) for r in boot if math.isfinite(float(r[metric]))]
            limits["test3", outcome, "all", metric] = (quantile7(values, .025), quantile7(values, .975))
            assert abs(statistics.fmean(values)-metrics[metric]) < statistics.stdev(values), (outcome, metric)
        for subgroup in ["all", "industry", "other"]:
            if subgroup == "all":
                part_audit = audit
            else:
                part = [p for p in population if p[0]["cls"] == subgroup]
                part_audit = sponsor_audit(part, outcome)
                for metric, value in sponsor_audit_metrics(part_audit, len({p[0]["key"] for p in part}), len(part)).items():
                    calc["test3", outcome, subgroup, metric] = value
            saved = [r for r in seed_rows if r["outcome"] == outcome and r["subgroup"] == subgroup]
            assert len(saved) == 200
            for r, s in zip(saved, part_audit["seeds"]):
                assert int(r["budget"]) == part_audit["budget"] and int(r["found_risk"]) == part_audit["risk"]
                for name, value in s.items():
                    assert int(r[name]) == value, (outcome, subgroup, r["seed"], name)
        for b, draw in enumerate(boot_draws):
            multiplicity = {}
            for j in draw:
                multiplicity[test[j]] = multiplicity.get(test[j], 0)+1
            repeat = [i for i, t in enumerate(scored) for _ in range(multiplicity.get(t["key"], 0))]
            y_boot = [scored[i][outcome] for i in repeat]
            base_boot, hist_boot = [base[i] for i in repeat], [hist[i] for i in repeat]
            r = boot[b]
            close(r["auc_difference"], sponsor_auc(hist_boot, y_boot)-sponsor_auc(base_boot, y_boot), (outcome, b, "auc"))
            close(r["brier_difference"], sponsor_brier(hist_boot, y_boot)-sponsor_brier(base_boot, y_boot), (outcome, b, "brier"))
            copies = [i for i in audit["ranked"] for _ in range(multiplicity.get(audit["keys"][i], 0))]
            budget = (len(copies)+5)//10
            assert int(r["budget"]) == budget
            close(r["yield_risk"], 100*sum(audit["ys"][i] for i in copies[:budget])/budget, (outcome, b, "risk"))
        verdict_inputs[outcome] = (calc["test1", outcome, "all", "p_value"], calc["test2", outcome, "all", "auc_difference"],
                                   limits["test2", outcome, "all", "auc_difference"][0], metrics["adaptive_to_risk"],
                                   limits["test3", outcome, "all", "adaptive_to_risk"][0])
    p, auc_difference, auc_lower, ratio, ratio_lower = verdict_inputs["D"]
    for metric, value in [("test1_p_below_0.01", p < .01), ("auc_difference_at_least_0.02", auc_difference >= .02),
                          ("auc_lower_above_0", auc_lower > 0), ("ratio_at_least_1.2", ratio >= 1.2),
                          ("ratio_lower_above_1", ratio_lower > 1), ("verdict", sponsor_verdict(*verdict_inputs["D"]))]:
        calc["verdict", "D", "all", metric] = str(value)
    return calc, limits


def check_sponsor_summary(summary, calc, limits):
    saved = {(r["section"], r["outcome"], r["subgroup"], r["metric"]): r for r in summary}
    assert set(saved) == set(calc), sorted(set(saved) ^ set(calc))[:5]
    for key, value in calc.items():
        if isinstance(value, str):
            assert saved[key]["value"] == value, key
        else:
            close(saved[key]["value"], value, key)
        if key in limits:
            close(saved[key]["lower"], limits[key][0], key)
            close(saved[key]["upper"], limits[key][1], key)
        else:
            assert saved[key]["lower"] == "" and saved[key]["upper"] == "", key


def verify_sponsor_clustering(folder, split):
    base = ROOT/"validation/sponsor-clustering"
    manifest = json.loads((base/"fetch-manifest.json").read_text())
    hashes = {r["input"]: r for r in rows(folder/"input-hashes.csv")}
    assert set(hashes) == {"plan", "pool", "sponsor_map", "script"}
    for r in hashes.values():
        actual = hashlib.sha256((ROOT/r["path"]).read_bytes()).hexdigest()
        assert r["sha256_before"] == r["sha256_after"] == actual and r["unchanged"] == "True", r["input"]
    assert hashes["pool"]["sha256_after"] == manifest["pool_sha256"]
    assert hashes["sponsor_map"]["sha256_after"] == manifest["sponsor_map_sha256"]
    trials = sponsor_trials(split)
    assert manifest["match_by_split"][split] == dict(pool_rows=len(trials), matched=len(trials), rate=1.0)
    info = json.loads((folder/"run-info.json").read_text())
    assert info["split"] == split and info["permutations"] == 2000 and info["bootstraps"] == 2000 and info["audit_seeds"] == 200
    summary = rows(folder/"summary.csv")
    calc, limits = sponsor_recompute(folder, trials)
    check_sponsor_summary(summary, calc, limits)
    assert info["verdict"] == calc["verdict", "D", "all", "verdict"]
    for arguments, expected in [((.0005, .03, .01, 1.3, 1.1), "CLUSTERING-SUPPORTED"), ((.0099, .02, .001, 1.2, 1.001), "CLUSTERING-SUPPORTED"),
                                ((.0005, .03, .01, 1.1, .9), "WEAK"), ((.0005, .03, .01, 1.3, .9), "WEAK"), ((.0005, .019, .01, 1.3, 1.1), "WEAK"),
                                ((.0005, .03, -.001, 1.3, 1.1), "WEAK"), ((.01, .5, .4, 2, 1.5), "NO-CLUSTERING"), ((.02, .5, .4, 2, 1.5), "NO-CLUSTERING")]:
        assert sponsor_verdict(*arguments) == expected, arguments
    print(f"Sponsor clustering {split}: {len(calc)} summary values, model fits, audit seeds, intervals and verdict "
          f"({calc['verdict', 'D', 'all', 'verdict']}) recomputed in plain Python from raw files: PASS")
    wrong = [dict(r) for r in summary]
    target = next(r for r in wrong if r["section"] == "test1" and r["outcome"] == "D" and r["subgroup"] == "all" and r["metric"] == "statistic")
    target["value"] = str(float(target["value"])+1)
    try:
        check_sponsor_summary(wrong, calc, limits)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: wrong Test 1 statistic went undetected")
    _, strata, expected, _, _ = sponsor_strata(trials, "D")
    shuffled = sponsor_shuffled(strata, expected, "D", 200, 20261199)
    null = sorted(float(r["statistic_D"]) for r in rows(folder/"test1-null.csv"))
    null_mean = statistics.fmean(null)
    standard_error = statistics.stdev(null)/math.sqrt(200)
    assert abs(statistics.fmean(shuffled)-null_mean) < 4*standard_error, (statistics.fmean(shuffled), null_mean)
    extreme = sum((1+sum(v >= s for v in null))/2001 < .01 for s in shuffled)
    assert extreme <= 10, extreme
    print(f"Sponsor clustering {split} negative controls: wrong statistic detected; 200 within-strata shuffles of D give "
          f"mean ratio {statistics.fmean(shuffled)/null_mean:.3f} and {extreme} of 200 below p = 0.01: PASS")


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
    check_hashes(ROOT/"validation/reference-oracle/source-hashes.csv", ROOT/"validation/v0.3.1-source")
    for filename in ["posterior-predictive-moments.csv", "data-contract-controls.csv", "negative-controls.csv"]:
        for r in rows(ROOT/"validation/reference-oracle"/filename):
            assert r.get("pass", r.get("detected")) == "TRUE", (filename, r)
    print("Independent reference predictive oracle, contracts and negative controls: PASS")
    check_hashes(ROOT/"validation/pooling-guards/source-hashes.csv", ROOT/"validation/v0.3.1-source")
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
    verify_contamination_study(ROOT/"validation/contamination-study", ROOT/"validation/v0.4.0-source")
    check_hashes(ROOT/"validation/contamination-oracle/source-hashes.csv", ROOT/"validation/v0.4.0-source")
    for filename, field in [("posterior-moments.csv", "pass"), ("row-posteriors.csv", "pass"),
                            ("negative-controls.csv", "detected"), ("data-contract-controls.csv", "pass")]:
        for r in rows(ROOT/"validation/contamination-oracle"/filename):
            assert r[field] == "TRUE", (filename, r)
    print("Independent contamination oracle, contracts and negative controls: PASS")
    verify_contamination_supplement(ROOT/"validation/contamination-supplement")
    verify_real_data_study(ROOT/"validation/real-data-study", ROOT/"validation/v0.4.0-source")
    verify_merged_audit(ROOT/"validation/merged-audit-pilot/study", ROOT/"validation/v0.4.0-source")
    verify_audit_plan_study(ROOT/"validation/audit-plan-study", ROOT/"validation/v0.4.1-source")
    verify_sponsor_clustering(ROOT/"validation/sponsor-clustering/development", "development")
    verify_sponsor_clustering(ROOT/"validation/sponsor-clustering/evaluation", "evaluation")
    print("This verifies correspondence and repeatability, not inferential validity.")


if __name__ == "__main__":
    main()
