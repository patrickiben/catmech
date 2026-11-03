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


REAL_AUDIT_ANALYSES = [
    ("rwanda_month", "pnc"), ("rwanda_month", "syphilis"), ("rwanda_month", "iron"), ("rwanda_month", "anc4"),
    ("rwanda_quarter", "pnc"), ("rwanda_quarter", "syphilis"), ("rwanda_quarter", "iron"), ("rwanda_quarter", "anc4"),
    ("effect_size", "effect_smd"), ("effect_size", "effect_r")]
REAL_AUDIT_METHODS = ["oracle", "recorded", "audited_only", "linear_mi", "raking", "contamination"]
RWANDA_COLUMNS = {"pnc": ("pnc_newborn_24h_h", "pnc_newborn_24h_r", "hf_deliveries_h"),
                  "syphilis": ("anc_syphilis_test_h", "anc_syphilis_test_r", "ancnewreg_h"),
                  "iron": ("womenreciron_h", "womenreciron_r", "ancnewreg_h"),
                  "anc4": ("anc_4visits_h", "anc_4visits_r", "ancnewreg_h")}


def real_audit_pools():
    """Rebuild each analysis frame from the pool files, independently of R: (stratum, recorded, reference, covariate)."""
    centers = [r for r in rows(ROOT/"validation/real-audit/rwanda-pool.csv") if r["facility_type"] == "health_center"]
    frames = {}
    for name, (w, reference, x) in RWANDA_COLUMNS.items():
        kept = [r for r in centers if r[w] != "" and r[reference] != "" and r[x] != ""]
        frames[("rwanda_month", name)] = [(r["hosp_id"], float(r[w]), float(r[reference]), float(r[x])) for r in kept]
        months = {}
        for r in centers:
            months.setdefault(r["facility_id"], []).append(r)
        frames[("rwanda_quarter", name)] = [
            (group[0]["hosp_id"], sum(float(r[w]) for r in group), sum(float(r[reference]) for r in group), sum(float(r[x]) for r in group))
            for group in months.values() if all(r[w] != "" and r[reference] != "" and r[x] != "" for r in group)]
    effect = [r for r in rows(ROOT/"validation/real-audit/effect-size-pool.csv") if r["info"] in ("0", "1")]
    for name, family in [("effect_smd", "smd"), ("effect_r", "r")]:
        kept = [r for r in effect if r["family"] == family]
        counts = {}
        for r in kept:
            counts[r["meta"]] = counts.get(r["meta"], 0)+1
        frames[("effect_size", name)] = [(r["meta"], float(r["recorded"]), float(r["reference"]), math.log(float(r["n"])))
                                         for r in kept if counts[r["meta"]] >= 2]
    return frames


def real_audit_slope(frame, column):
    strata = sorted({row[0] for row in frame})
    design = [[1.0]+[float(row[0] == level) for level in strata[1:]]+[row[3]] for row in frame]
    return least_squares_last(design, [row[column] for row in frame])


def verify_rwanda_published_table():
    """The pool must reproduce the published facility table: facilities, exact matches and within-range counts."""
    pool = rows(ROOT/"validation/real-audit/rwanda-pool.csv")
    published = [("ancnewreg_r", "ancnewreg_h", 69, 49, 66), ("anc1vis_16weeks_r", "anc1vis1trim_h", 68, 5, 17),
                 ("anc_syphilis_test_r", "anc_syphilis_test_h", 68, 25, 46), ("womenreciron_r", "womenreciron_h", 60, 29, 47),
                 ("anc_4visits_r", "anc_4visits_h", 69, 2, 12), ("hf_deliveries_r", "hf_deliveries_h", 76, 42, 75),
                 ("livebirths_r", "livebirths_h", 76, 39, 72), ("pnc_newborn_24h_r", "pnc_newborn_24h_h", 65, 10, 53)]
    within_matches = 0
    for reference, recorded, n, exact, within in published:
        months = {}
        for r in pool:
            months.setdefault(r["facility_id"], []).append(r)
        ratios = []
        for group in months.values():
            if all(r[reference] != "" and r[recorded] != "" for r in group):
                a, b = sum(float(r[reference]) for r in group), sum(float(r[recorded]) for r in group)
                if a == 0 and b == 0:
                    ratios.append(1.0)
                elif b > 0:
                    ratios.append(a/b)
                else:
                    ratios.append(math.inf)
        assert len(ratios) == n, (reference, len(ratios))
        assert sum(value == 1 for value in ratios) == exact, (reference, "exact")
        count = sum(0.9 <= value <= 1.1 for value in ratios)
        assert abs(count-within) <= 1, (reference, count, within)
        within_matches += count == within
    assert within_matches == 6, within_matches
    return within_matches


def real_audit_planned(settings):
    reps = {"rwanda_month": settings["month_reps"], "rwanda_quarter": settings["quarter_reps"], "effect_size": settings["effect_reps"]}
    expected = set()
    for index, (level, analysis) in enumerate(REAL_AUDIT_ANALYSES, 1):
        for fraction_index, fraction in enumerate(settings["fractions"], 1):
            for replicate in range(1, reps[level]+1):
                for method in REAL_AUDIT_METHODS:
                    expected.add((level, analysis, float(fraction), method, replicate))
    return expected


def verify_real_audit_records(records, settings):
    actual = {(r["level"], r["analysis"], float(r["fraction"]), r["method"], int(r["replicate"])) for r in records}
    assert actual == real_audit_planned(settings) and len(actual) == len(records), "Missing, duplicate or unplanned requests"
    position = {pair: index for index, pair in enumerate(REAL_AUDIT_ANALYSES, 1)}
    for r in records:
        index = position[(r["level"], r["analysis"])]
        fraction_index = settings["fractions"].index(float(r["fraction"]))+1
        replicate = int(r["replicate"])
        assert int(r["selection_seed"]) == settings["selection_base"]+index*1000000+fraction_index*100000+replicate
        assert int(r["imputation_seed"]) == settings["imputation_base"]+index*1000000+fraction_index*100000+replicate


def real_audit_groups(records):
    groups = {}
    for r in records:
        groups.setdefault((r["scenario"], r["method"]), []).append(r)
    return groups


def verify_real_audit_summary(summary, groups):
    close = lambda a, b: math.isclose(float(a), b, rel_tol=1e-7, abs_tol=1e-9)
    for row in summary:
        group = groups[(row["scenario"], row["method"])]
        good = [r for r in group if r["status"] == "ok"]
        assert int(row["requested"]) == len(group) and int(row["returned_intervals"]) == len(good), row["scenario"]
        recovery = [float(r["recovery_error"]) for r in good]
        error = [float(r["error"]) for r in good]
        assert close(row["recovery_rmse"], math.sqrt(sum(v*v for v in recovery)/len(recovery))), ("recovery", row["scenario"], row["method"])
        assert close(row["rmse"], math.sqrt(sum(v*v for v in error)/len(error))), ("rmse", row["scenario"], row["method"])
        assert close(row["bias"], sum(error)/len(error)), ("bias", row["scenario"], row["method"])
        assert close(row["coverage"], sum(float(r["coverage"]) for r in good)/len(good)), ("coverage", row["scenario"], row["method"])
        assert close(row["coverage_cluster"], sum(float(r["coverage_cluster"]) for r in good)/len(good)), ("cluster", row["scenario"], row["method"])
        assert close(row["mean_width"], sum(float(r["width"]) for r in good)/len(good)), ("width", row["scenario"], row["method"])
        assert close(row["mean_width_cluster"], sum(float(r["width_cluster"]) for r in good)/len(good)), ("cluster width", row["scenario"], row["method"])
        assert close(row["interval_return_rate"], len(good)/len(group)), ("return rate", row["scenario"], row["method"])
        assert int(row["convergence_failures"]) == sum(r["status"] == "convergence_failure" for r in group), ("convergence", row["scenario"], row["method"])
        assert int(row["refused"]) == sum(r["status"] == "refused" for r in group) and int(row["other_errors"]) == sum(r["status"] == "error" for r in group), row["scenario"]


def verify_real_audit_writeup(summary):
    """Regenerate the scored output, tie the post-hoc table to the study, and check the write-up's numbers against both."""
    import re
    import subprocess
    import sys
    scored = subprocess.run([sys.executable, "-I", str(ROOT/"tools/real-audit-score.py")], capture_output=True, text=True, check=True).stdout
    saved = (ROOT/"validation/real-audit-decisions.txt").read_text(encoding="utf-8")
    assert scored.strip() == saved.strip(), "Scored output differs from the saved decisions"
    lookup = {(r["level"], r["analysis"], r["fraction"], r["method"]): float(r["recovery_rmse"]) for r in summary}
    posthoc = rows(ROOT/"validation/real-audit-posthoc/swap-in-summary.csv")
    assert len(posthoc) == 20
    for r in posthoc:
        key = (r["level"], r["analysis"], r["fraction"])
        assert math.isclose(float(r["recorded_recovery_rmse"]), lookup[key+("recorded",)], rel_tol=1e-6), key
        assert math.isclose(float(r["linear_mi_recovery_rmse"]), lookup[key+("linear_mi",)], rel_tol=1e-6), key
        assert math.isclose(float(r["contamination_recovery_rmse"]), lookup[key+("contamination",)], rel_tol=1e-6), key
        assert float(r["swap_in_recovery_rmse"]) < float(r["recorded_recovery_rmse"]), key
    text = (ROOT/"validation/REAL_AUDIT_SYNTHESIS.md").read_text(encoding="utf-8")
    triples = verify_real_audit_tables(text, summary)
    assert triples == 100, triples
    limits = re.findall(r"\d\.\d{3} \[\d\.\d{3}, \d\.\d{3}\]", text)
    assert len(limits) >= 20 and all(t in saved for t in limits), [t for t in limits if t not in saved][:3]
    fourths = re.findall(r"\b0\.\d{4}\b", text)
    table = {f"{float(r[f]):.4f}" for r in posthoc for f in ["recorded_recovery_rmse", "swap_in_recovery_rmse", "linear_mi_recovery_rmse", "contamination_recovery_rmse", "swap_in_recovery_rmse_same_cohorts"]}
    summary_values = {f"{value:.4f}" for value in lookup.values()}
    missing = [v for v in fourths if v not in table and v not in summary_values]
    assert not missing, missing
    outcome = real_audit_rule2(summary, rows(ROOT/"validation/real-audit-study/recovery-ratios.csv"))
    assert outcome == "INCONCLUSIVE" and f"Outcome: {outcome}" in saved and f"the result is **{outcome.lower()}**" in text, outcome
    return triples, len(limits)


def verify_real_audit_tables(text, summary):
    """Check each method-table cell by position: section, row label, then method column."""
    sections = {"**Rwanda facility-months**": "rwanda_month", "**Rwanda quarter totals**": "rwanda_quarter", "**Effect sizes, the control**": "effect_size"}
    labels = {"PNC": "pnc", "Syphilis": "syphilis", "Iron": "iron", "ANC 4th visit": "anc4", "SMD family": "effect_smd", "r family": "effect_r"}
    lookup = {(r["level"], r["analysis"], float(r["fraction"]), r["method"]): r for r in summary}
    level, checked = None, 0
    for line in text.splitlines():
        for marker, name in sections.items():
            if line.startswith(marker):
                level = name
        if line.startswith("**RMSE against the pool slope**"):
            level = None
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if level is None or len(cells) != 8 or ", " not in cells[0] or cells[0].split(", ")[0] not in labels:
            continue
        name, audit = cells[0].rsplit(", ", 1)
        fraction = float(audit.rstrip("%"))/100
        for method, cell in zip(["recorded", "audited_only", "linear_mi", "raking", "contamination"], cells[2:7]):
            row = lookup[(level, labels[name], fraction, method)]
            expected = f'{float(row["recovery_rmse"]):.3f}, {float(row["bias"]):+.3f}, {100*float(row["coverage_cluster"]):.1f}%'
            assert cell == expected, (level, name, audit, method, cell, expected)
            checked += 1
        model = lookup[(level, labels[name], fraction, "contamination")]
        assert cells[1] == f'{model["returned_intervals"]}/{model["requested"]}', (level, name, audit)
    return checked


def real_audit_rule2(summary, ratio_rows):
    """Re-derive Rule 2 of the plan from saved tables, separately from the scoring script."""
    lookup = {(r["level"], r["analysis"], float(r["fraction"]), r["method"]): r for r in summary}
    ratios = {(r["level"], r["analysis"], float(r["fraction"]), r["method"], r["comparator"]): r for r in ratio_rows}
    supported, contradicted = 0, 0
    for analysis in ["pnc", "syphilis", "iron"]:
        pair = ratios[("rwanda_month", analysis, .25, "contamination", "linear_mi")]
        model = lookup[("rwanda_month", analysis, .25, "contamination")]
        oracle = lookup[("rwanda_month", analysis, .25, "oracle")]
        quarter = ratios[("rwanda_quarter", analysis, .25, "contamination", "linear_mi")]
        if (float(pair["upper"]) < 1 and float(model["interval_return_rate"]) >= .90
                and float(model["coverage_cluster"]) >= float(oracle["coverage_cluster"])-.03 and float(quarter["recovery_rmse_ratio"]) < 1):
            supported += 1
        contradicted += float(pair["lower"]) > 1
    return "SUPPORTED" if supported >= 2 else "CONTRADICTED" if contradicted >= 2 else "INCONCLUSIVE"


def verify_real_audit_replay(records):
    """Saved results must match a replay from seeds, which catches relabeled or rescaled results."""
    lookup = {(r["scenario"], r["replicate"], r["method"]): r for r in records}
    replay = rows(ROOT/"validation/real-audit-replay/replay.csv")
    assert len(replay) == 600 and {r["method"] for r in replay} == set(REAL_AUDIT_METHODS)
    for row in replay:
        saved = lookup[(row["scenario"], row["replicate"], row["method"])]
        assert row["status"] == saved["status"], row
        if row["status"] == "ok":
            for field in ["estimate", "se", "se_cluster"]:
                assert abs(float(row[field])-float(saved[field])) < 1e-7, (row["scenario"], row["replicate"], row["method"], field)


def verify_real_audit_study(folder, require_study=True):
    if require_study:
        for row in rows(folder/"source-hashes.csv"):
            candidate = ROOT/"validation/v0.4.1-source"/row["file"]
            path = candidate if candidate.exists() else ROOT/row["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    manifest = json.loads((ROOT/"validation/real-audit/fetch-manifest.json").read_text())
    for name, entry in [("rwanda-pool.csv", manifest["rwanda"]), ("effect-size-pool.csv", manifest["effect_sizes"])]:
        assert hashlib.sha256((ROOT/"validation/real-audit"/name).read_bytes()).hexdigest() == entry["pool_sha256"], name
    for entry in [manifest["rwanda"], manifest["effect_sizes"]]:
        assert hashlib.sha256((ROOT/"validation/real-audit/raw"/entry["file"]).read_bytes()).hexdigest() == entry["sha256"], entry["file"]
    within_matches = verify_rwanda_published_table()
    effect = rows(ROOT/"validation/real-audit/effect-size-pool.csv")
    assert [sum(r["info"] == str(code) for r in effect) for code in range(4)] == [276, 74, 54, 96]
    settings = json.loads((folder/"settings.json").read_text())
    if require_study:
        assert settings["mode"] == "study" and settings["m"] == 20 and settings["fractions"] == [0.25, 0.5]
        assert settings["month_reps"] == 500 and settings["quarter_reps"] == 200 and settings["effect_reps"] == 500
        assert (settings["cohort_base"], settings["selection_base"], settings["imputation_base"]) == (400100000, 410100000, 420000000)
    assert not (folder/"checkpoint.rds").exists()
    records = rows(folder/"replicates.csv")
    verify_real_audit_records(records, settings)
    frames = real_audit_pools()
    saved = {(r["level"], r["analysis"]): r for r in rows(folder/"pool-targets.csv")}
    for key, frame in frames.items():
        target, recorded = real_audit_slope(frame, 2), real_audit_slope(frame, 1)
        tol = {"effect_smd": .049, "effect_r": .025}.get(key[1], 0)
        wrong = sum(abs(row[1]-row[2]) > tol+1e-8 for row in frame)
        assert int(saved[key]["rows"]) == len(frame) and int(saved[key]["wrong_rows"]) == wrong, key
        assert math.isclose(float(saved[key]["target_coefficient"]), target, rel_tol=1e-7), key
        assert math.isclose(float(saved[key]["recorded_coefficient"]), recorded, rel_tol=1e-7), key
        assert all(math.isclose(float(r["truth"]), target, rel_tol=1e-7) for r in records if (r["level"], r["analysis"]) == key), key
    groups = real_audit_groups(records)
    oracle = {(r["scenario"], r["replicate"]): float(r["estimate"]) for r in records if r["method"] == "oracle"}
    for r in records:
        if r["status"] != "ok":
            continue
        estimate, truth = float(r["estimate"]), float(r["truth"])
        assert math.isclose(float(r["error"]), estimate-truth, rel_tol=1e-6, abs_tol=1e-9)
        assert math.isclose(float(r["recovery_error"]), estimate-oracle[(r["scenario"], r["replicate"])], rel_tol=1e-6, abs_tol=1e-9)
        assert (float(r["lower"]) <= truth <= float(r["upper"])) == (float(r["coverage"]) == 1)
        assert (float(r["lower_cluster"]) <= truth <= float(r["upper_cluster"])) == (float(r["coverage_cluster"]) == 1)
    summary = rows(folder/"summary.csv")
    assert len(summary) == len(groups)
    verify_real_audit_summary(summary, groups)
    generator = random.Random(20261101)
    worst = 0.0
    for row in rows(folder/"recovery-ratios.csv"):
        first = {r["replicate"]: float(r["recovery_error"])**2 for r in groups[(row["scenario"], row["method"])] if r["status"] == "ok"}
        second = {r["replicate"]: float(r["recovery_error"])**2 for r in groups[(row["scenario"], row["comparator"])] if r["status"] == "ok"}
        shared = sorted(set(first) & set(second))
        a, b = [first[k] for k in shared], [second[k] for k in shared]
        assert int(row["complete_pairs"]) == len(shared)
        assert math.isclose(float(row["recovery_rmse_ratio"]), math.sqrt(sum(a)/sum(b)), rel_tol=1e-7), (row["scenario"], row["method"])
        draws = []
        for _ in range(400):
            picks = generator.choices(range(len(a)), k=len(a))
            denominator = sum(b[i] for i in picks)
            if denominator > 0:
                draws.append(math.sqrt(sum(a[i] for i in picks)/denominator))
        draws.sort()
        worst = max(worst, abs(draws[int(.975*len(draws))]-float(row["upper"]))/max(float(row["recovery_rmse_ratio"]), 1e-9),
                    abs(draws[int(.025*len(draws))]-float(row["lower"]))/max(float(row["recovery_rmse_ratio"]), 1e-9))
    if require_study:
        assert worst < .15, worst
        verify_real_audit_replay(records)
        swapped = [dict(r, method={"linear_mi": "contamination", "contamination": "linear_mi"}.get(r["method"], r["method"])) for r in records]
        rescaled = [dict(r, se_cluster=str(float(r["se_cluster"])*3)) if r["method"] == "contamination" and r["status"] == "ok" else r for r in records]
        for wrong, label in [(swapped, "relabeled methods"), (rescaled, "rescaled standard errors")]:
            try:
                verify_real_audit_replay(wrong)
            except AssertionError:
                pass
            else:
                raise AssertionError(f"Negative control: {label} went undetected")
    mutated = [dict(r) for r in summary]
    mutated[0]["recovery_rmse"] = str(float(mutated[0]["recovery_rmse"])+.2)
    try:
        verify_real_audit_summary(mutated, groups)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: corrupted real-audit summary went undetected")
    try:
        verify_real_audit_records([r for r in records if r["analysis"] != "anc4"], settings)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: omitted planned analysis went undetected")
    print(f"{folder.name}: all {len(records)} planned requests, {len(summary)} recomputed groups, independent pool targets, published Rwanda table ({within_matches} of 8 within-range counts exact, the other 2 off by one), published effect-size counts, ratio limits recomputed, two negative controls: PASS")
    if require_study:
        triples, limits = verify_real_audit_writeup(summary)
        print(f"Real-audit write-up: scored output regenerated identically, post-hoc table tied to the study, Rule 2 re-derived, {triples} table entries checked by position, {limits} ratio intervals found in the scored output, replay from seeds matches, relabeled and rescaled results caught: PASS")


def verify_clustered_summary(summary, groups):
    close = lambda a, b: math.isclose(float(a), b, rel_tol=1e-7, abs_tol=1e-9)
    for row in summary:
        group = groups[(row["scenario"], row["method"])]
        good = [r for r in group if r["status"] == "ok"]
        assert int(row["requested"]) == len(group) and int(row["returned_intervals"]) == len(good), row["scenario"]
        assert int(row["refused"]) == sum(r["status"] == "refused" for r in group) and int(row["other_errors"]) == sum(r["status"] == "error" for r in group), row["scenario"]
        assert close(row["coverage"], sum(float(r["coverage"]) for r in good)/len(good)), ("coverage", row["scenario"], row["method"])
        assert close(row["mean_width"], sum(float(r["width"]) for r in good)/len(good)), ("width", row["scenario"], row["method"])
        assert close(row["recovery_rmse"], math.sqrt(sum(float(r["recovery_error"])**2 for r in good)/len(good))), ("recovery", row["scenario"], row["method"])
        assert close(row["rmse"], math.sqrt(sum(float(r["error"])**2 for r in good)/len(good))), ("rmse", row["scenario"], row["method"])


def verify_clustered_study(folder):
    """Clustered-audit study: planned requests, seeds, targets, summaries, the earlier study's estimates, and the verdict."""
    import re
    import subprocess
    import sys
    for row in rows(folder/"source-hashes.csv"):
        candidate = ROOT/"validation/v0.5.0-source"/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    settings = json.loads((folder/"settings.json").read_text())
    assert settings["mode"] == "study" and settings["reps"] == 500 and settings["m"] == 20
    assert not (folder/"checkpoint.rds").exists()
    records = rows(folder/"replicates.csv")
    held = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early"]
    development = ["pnc", "syphilis", "iron", "anc4"]
    grid = rows(folder/"simulation-grid.csv")
    assert len(grid) == 56
    labels = [(("control_one_row" if g["rows"] == "1" else f"g{g['clusters']}_icc{round(100*float(g['icc'])):02d}_{g['errors']}"), float(g["fraction"])) for g in grid]
    expected = set()
    for level, names in [("held_out", held), ("development", development)]:
        for name in names:
            for fraction in [.25, .5]:
                expected |= {(level, name, fraction, m, str(r)) for m in ["oracle", "recorded", "linear_mi", "linear_mi_clustered"] for r in range(1, 501)}
    for label, fraction in labels:
        expected |= {("simulation", label, fraction, m, str(r)) for m in ["oracle", "recorded", "linear_mi", "linear_mi_clustered"] for r in range(1, 501)}
    actual = {(r["level"], r["analysis"], float(r["fraction"]), r["method"], r["replicate"]) for r in records}
    assert actual == expected and len(actual) == len(records), "Missing, duplicate or unplanned requests"
    seeds = {name: i for i, name in enumerate(held, 1)}
    for r in records:
        if r["level"] == "held_out":
            fraction_index = 1 if float(r["fraction"]) == .25 else 2
            assert int(r["imputation_seed"]) == 880000000+seeds[r["analysis"]]*1000000+fraction_index*100000+int(r["replicate"])
    centers = [r for r in rows(ROOT/"validation/real-audit/rwanda-pool.csv") if r["facility_type"] == "health_center"]
    number = lambda r, *names: None if any(r[n] == "" for n in names) else sum(float(r[n]) for n in names)
    definitions = {"lbw": (("lbw_notpremature_h", "lbw_premature_h"), ("lbw_r",), "hf_deliveries_h"), "resus": (("newborns_nocryres_succ_h",), ("nocry_res_success_r",), "hf_deliveries_h"),
                   "stillbirths": (("sbmacerated_h", "sbfresh_h"), ("sbmacerated_r", "sbfresh_r"), "hf_deliveries_h"), "deliveries": (("hf_deliveries_h",), ("hf_deliveries_r",), "ancnewreg_h"),
                   "anc1_early": (("anc1vis1trim_h",), ("anc1vis_16weeks_r",), "ancnewreg_h")}
    saved = {(r["level"], r["analysis"]): r for r in rows(folder/"pool-targets.csv")}
    for name, (w, reference, x) in definitions.items():
        frame = [(r["hosp_id"], number(r, *w), number(r, *reference), float(r[x])) for r in centers if number(r, *w) is not None and number(r, *reference) is not None and r[x] != ""]
        assert int(saved[("held_out", name)]["rows"]) == len(frame), name
        assert int(saved[("held_out", name)]["wrong_rows"]) == sum(abs(a[1]-a[2]) > 1e-8 for a in frame), name
        assert math.isclose(float(saved[("held_out", name)]["target_coefficient"]), real_audit_slope(frame, 2), rel_tol=1e-7), name
        assert math.isclose(float(saved[("held_out", name)]["recorded_coefficient"]), real_audit_slope(frame, 1), rel_tol=1e-7), name
    groups = {}
    for r in records:
        groups.setdefault((r["scenario"], r["method"]), []).append(r)
    summary = rows(folder/"summary.csv")
    assert len(summary) == len(groups)
    verify_clustered_summary(summary, groups)
    wrong = [dict(r) for r in summary]
    wrong[0]["coverage"] = str(float(wrong[0]["coverage"])+.1)
    try:
        verify_clustered_summary(wrong, groups)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: corrupted clustered summary went undetected")
    earlier = {}
    for r in rows(ROOT/"validation/real-audit-study/replicates.csv"):
        if r["level"] == "rwanda_month" and r["method"] in ("oracle", "recorded", "linear_mi"):
            earlier[(r["analysis"], float(r["fraction"]), r["replicate"], r["method"])] = r
    matched = 0
    for r in records:
        if r["level"] == "development" and r["method"] != "linear_mi_clustered":
            old = earlier[(r["analysis"], float(r["fraction"]), r["replicate"], r["method"])]
            assert r["status"] == old["status"], r["scenario"]
            if r["status"] == "ok":
                assert abs(float(r["estimate"])-float(old["estimate"])) < 1e-8 and abs(float(r["se"])-float(old["se_cluster"])) < 1e-8, (r["scenario"], r["replicate"], r["method"])
                matched += 1
    scored = subprocess.run([sys.executable, "-I", str(ROOT/"tools/clustered-score.py")], capture_output=True, text=True, check=True).stdout
    assert scored.strip() == (ROOT/"validation/clustered-decisions.txt").read_text(encoding="utf-8").strip(), "Scored output differs from the saved decisions"
    posthoc = subprocess.run([sys.executable, "-I", str(ROOT/"tools/clustered-posthoc.py")], capture_output=True, text=True, check=True).stdout
    assert posthoc.strip() == (ROOT/"validation/clustered-posthoc.txt").read_text(encoding="utf-8").strip(), "Post-hoc output differs from the saved table"
    verdict = re.search(r"VERDICT: (\S+)", scored).group(1)
    lookup = {(r["level"], r["analysis"], float(r["fraction"]), r["method"]): float(r["coverage"]) for r in summary}
    ratio = {r["analysis"]: float(r["recovery_upper"]) for r in rows(folder/"ratios.csv") if r["level"] == "held_out" and float(r["fraction"]) == .25}
    closer = sum(abs(lookup[("held_out", a, .25, "linear_mi_clustered")]-lookup[("held_out", a, .25, "oracle")]) < abs(lookup[("held_out", a, .25, "linear_mi")]-lookup[("held_out", a, .25, "oracle")]) for a in held)
    near = sum(abs(lookup[("held_out", a, .25, "linear_mi_clustered")]-lookup[("held_out", a, .25, "oracle")]) <= .05 for a in held)
    accurate = sum(ratio[a] <= 1.05 for a in held)
    eligible = [(g, float(g["fraction"])) for g in grid if g["rows"] == "3" and int(g["clusters"]) >= 70]
    passing = 0
    for g, fraction in eligible:
        label = f"g{g['clusters']}_icc{round(100*float(g['icc'])):02d}_{g['errors']}"
        c, o = lookup[("simulation", label, fraction, "linear_mi_clustered")], lookup[("simulation", label, fraction, "oracle")]
        passing += c >= .92 and abs(c-o) <= .03
    coverage_holds = closer >= 4 and near >= 4
    derived = "GO" if coverage_holds and accurate >= 4 and passing >= .8*len(eligible) else "PARTIAL" if coverage_holds else "NO-GO"
    assert derived == verdict, (derived, verdict)
    synthesis = (ROOT/"validation/CLUSTERED_SYNTHESIS.md").read_text(encoding="utf-8")
    assert f"**{verdict}**" in synthesis, "Write-up does not state the verdict"
    sources = scored+posthoc+(ROOT/"validation/CLUSTERED_PLAN.md").read_text(encoding="utf-8")
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    unsupported = sorted(decimals(re.sub(r"\[doi:[^\]]*\]\([^)]*\)", "", synthesis))-decimals(sources))
    assert not unsupported, ("Numbers in the write-up with no source", unsupported)
    assert decimals("A made-up 12.345 here.")-decimals(sources), "Negative control: a made-up number went undetected"
    print(f"{folder.name}: all {len(records)} planned requests, {len(summary)} recomputed groups, held-out targets recomputed, {matched} development requests match the earlier study, verdict re-derived as {verdict}, corrupted-summary control: PASS")


def verify_clustered_fix_summary(summary, groups):
    close = lambda a, b: math.isclose(float(a), b, rel_tol=1e-7, abs_tol=1e-9)
    for row in summary:
        group = groups[(row["scenario"], row["method"])]
        good = [r for r in group if r["status"] == "ok"]
        assert int(row["requested"]) == len(group) and int(row["returned_intervals"]) == len(good), row["scenario"]
        assert int(row["refused"]) == sum(r["status"] == "refused" for r in group) and int(row["other_errors"]) == sum(r["status"] == "error" for r in group), row["scenario"]
        assert close(row["coverage"], sum(float(r["coverage"]) for r in good)/len(good)), ("coverage", row["scenario"], row["method"])
        assert close(row["mean_width"], sum(float(r["width"]) for r in good)/len(good)), ("width", row["scenario"], row["method"])
        assert close(row["median_width"], statistics.median(float(r["width"]) for r in good)), ("median width", row["scenario"], row["method"])
        assert close(row["recovery_rmse"], math.sqrt(sum(float(r["recovery_error"])**2 for r in good)/len(good))), ("recovery", row["scenario"], row["method"])
        assert close(row["rmse"], math.sqrt(sum(float(r["error"])**2 for r in good)/len(good))), ("rmse", row["scenario"], row["method"])
        assert close(row["median_between"], statistics.median(float(r["between_variance"]) for r in good)), ("between", row["scenario"], row["method"])


def verify_clustered_fix_study(folder):
    """Test of replacements for the clustered resampling: planned requests, seeds, targets, summaries, controls, selection, verdict and replay."""
    import gzip
    import re
    import subprocess
    import sys
    for row in rows(folder/"source-hashes.csv"):
        candidate = ROOT/"validation/clustered-fix-source"/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    settings = json.loads((folder/"settings.json").read_text())
    assert settings["mode"] == "study" and settings["reps"] == 500 and settings["m"] == 20
    assert not (folder/"checkpoint.rds").exists()
    with gzip.open(folder/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    five = ["lbw", "resus", "stillbirths", "deliveries", "anc1_early"]
    development = ["pnc", "syphilis", "iron", "anc4"]
    unused = ["ancnewreg", "livebirths"]
    methods = ["oracle", "linear_mi", "linear_mi_clustered", "no_strata", "rao_wu", "jackknife", "random_intercept"]
    grid = rows(folder/"simulation-grid.csv")
    assert len(grid) == 38
    expected = {(level, name, fraction, m, str(r)) for level, names in [("held_out", five), ("development", development), ("unused", unused)] for name in names for fraction in [.25, .5] for m in methods for r in range(1, 501)}
    expected |= {("simulation", g["label"], float(g["fraction"]), m, str(r)) for g in grid for m in methods for r in range(1, 501)}
    actual = {(r["level"], r["analysis"], float(r["fraction"]), r["method"], r["replicate"]) for r in records}
    assert actual == expected and len(actual) == len(records) == 210000, "Missing, duplicate or unplanned requests"
    index = {name: i for i, name in enumerate(five+unused, 1)} | {name: i for i, name in enumerate(development, 1)}
    cells = {(g["label"], float(g["fraction"])): i for i, g in enumerate(grid, 1)}
    for r in records:
        fraction_index = 1 if float(r["fraction"]) == .25 else 2
        if r["level"] == "simulation":
            seed = 1310000000+cells[(r["analysis"], float(r["fraction"]))]*100000+int(r["replicate"])
        else:
            base = 420000000 if r["level"] == "development" else 880000000
            seed = base+index[r["analysis"]]*1000000+fraction_index*100000+int(r["replicate"])
        assert int(r["imputation_seed"]) == seed, (r["scenario"], r["replicate"])
    centers = [r for r in rows(ROOT/"validation/real-audit/rwanda-pool.csv") if r["facility_type"] == "health_center"]
    number = lambda r, *names: None if any(r[n] == "" for n in names) else sum(float(r[n]) for n in names)
    definitions = {"lbw": (("lbw_notpremature_h", "lbw_premature_h"), ("lbw_r",), "hf_deliveries_h"), "resus": (("newborns_nocryres_succ_h",), ("nocry_res_success_r",), "hf_deliveries_h"),
                   "stillbirths": (("sbmacerated_h", "sbfresh_h"), ("sbmacerated_r", "sbfresh_r"), "hf_deliveries_h"), "deliveries": (("hf_deliveries_h",), ("hf_deliveries_r",), "ancnewreg_h"),
                   "anc1_early": (("anc1vis1trim_h",), ("anc1vis_16weeks_r",), "ancnewreg_h"), "ancnewreg": (("ancnewreg_h",), ("ancnewreg_r",), "hf_deliveries_h"), "livebirths": (("livebirths_h",), ("livebirths_r",), "ancnewreg_h"),
                   "pnc": (("pnc_newborn_24h_h",), ("pnc_newborn_24h_r",), "hf_deliveries_h"), "syphilis": (("anc_syphilis_test_h",), ("anc_syphilis_test_r",), "ancnewreg_h"), "iron": (("womenreciron_h",), ("womenreciron_r",), "ancnewreg_h"), "anc4": (("anc_4visits_h",), ("anc_4visits_r",), "ancnewreg_h")}
    saved = {r["analysis"]: r for r in rows(folder/"pool-targets.csv")}
    assert set(saved) == set(definitions)
    for name, (w, reference, x) in definitions.items():
        frame = [(r["hosp_id"], number(r, *w), number(r, *reference), float(r[x])) for r in centers if number(r, *w) is not None and number(r, *reference) is not None and r[x] != ""]
        assert int(saved[name]["rows"]) == len(frame), name
        assert int(saved[name]["wrong_rows"]) == sum(abs(a[1]-a[2]) > 1e-8 for a in frame), name
        assert math.isclose(float(saved[name]["target_coefficient"]), real_audit_slope(frame, 2), rel_tol=1e-7), name
        assert math.isclose(float(saved[name]["recorded_coefficient"]), real_audit_slope(frame, 1), rel_tol=1e-7), name
    groups = {}
    for r in records:
        groups.setdefault((r["scenario"], r["method"]), []).append(r)
    summary = rows(folder/"summary.csv")
    assert len(summary) == len(groups) == 420
    verify_clustered_fix_summary(summary, groups)
    wrong = [dict(r) for r in summary]
    wrong[0]["coverage"] = str(float(wrong[0]["coverage"])+.1)
    try:
        verify_clustered_fix_summary(wrong, groups)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: corrupted summary went undetected")
    lookup_records = {(r["analysis"], float(r["fraction"]), r["method"], r["replicate"]): r for r in records}
    same = lambda a, b, fields: all((a[f] == "" and b[f] == "") or (a[f] not in ("", "NA") and b[f] not in ("", "NA") and abs(float(a[f])-float(b[f])) < 1e-10) or (a[f] in ("", "NA") and b[f] in ("", "NA")) for f in fields)
    matched_first = 0
    for r in rows(ROOT/"validation/clustered-study/replicates.csv"):
        if r["level"] in ("held_out", "development") and r["method"] in ("oracle", "linear_mi", "linear_mi_clustered"):
            now = lookup_records[(r["analysis"], float(r["fraction"]), r["method"], r["replicate"])]
            assert now["status"] == r["status"] and same(now, r, ["estimate", "se", "df", "lower", "upper"]), ("first study", r["scenario"], r["replicate"], r["method"])
            matched_first += 1
    assert matched_first == 27000
    matched_review = 0
    for r in rows(ROOT/"validation/clustered-review/variants.csv"):
        if r["variant"] in ("no_strata", "rescaled"):
            now = lookup_records[(r["analysis"], float(r["fraction"]), "no_strata" if r["variant"] == "no_strata" else "rao_wu", r["replicate"])]
            assert (now["status"] == "ok") == (r["status"] == "ok") and same(now, r, ["estimate", "se", "df", "lower", "upper", "between_variance"]), ("review", r["analysis"], r["replicate"], r["variant"])
            matched_review += 1
    assert matched_review == 18000
    scored = subprocess.run([sys.executable, "-I", str(ROOT/"tools/clustered-fix-score.py")], capture_output=True, text=True, check=True).stdout
    assert scored.strip() == (ROOT/"validation/clustered-fix-decisions.txt").read_text(encoding="utf-8").strip(), "Scored output differs from the saved decisions"
    posthoc = subprocess.run([sys.executable, "-I", str(ROOT/"tools/clustered-fix-posthoc.py")], capture_output=True, text=True, check=True).stdout
    assert posthoc.strip() == (ROOT/"validation/clustered-fix-posthoc.txt").read_text(encoding="utf-8").strip(), "Post-hoc output differs from the saved table"
    selected = re.search(r"SELECTED: (\S+)", scored).group(1)
    verdict = re.search(r"VERDICT: (\S+)", scored).group(1)
    cover = {(r["analysis"], float(r["fraction"]), r["method"]): float(r["coverage"]) for r in summary}
    upper = {(r["analysis"], float(r["fraction"]), r["method"]): float(r["recovery_upper"]) for r in rows(folder/"ratios.csv")}
    point = {(r["analysis"], float(r["fraction"]), r["method"]): float(r["recovery_ratio"]) for r in rows(folder/"ratios.csv")}
    clustered = [(g["label"], float(g["fraction"])) for g in grid if g["errors"] in ("shift_low", "shift", "under", "outlier")]
    everything = [(g["label"], float(g["fraction"])) for g in grid]
    distance = {m: sum(abs(cover[(a, f, m)]-.95) for a, f in clustered)/24 for m in ("jackknife", "random_intercept")}
    ratio_mean = {m: sum(point[(a, f, m)] for a, f in everything)/38 for m in distance}
    chosen = min(distance, key=ratio_mean.get) if abs(distance["jackknife"]-distance["random_intercept"]) <= .005+1e-12 else min(distance, key=distance.get)
    assert chosen == selected, (chosen, selected)
    limits = {a: (cover[(a, .25, selected)] >= .90 or cover[(a, .25, selected)] >= cover[(a, .25, "oracle")]) and cover[(a, .25, selected)] <= .98 for a in five+development}
    coverage_holds = sum(limits.values()) >= 8 and limits["anc1_early"]
    accuracy_holds = sum(upper[(a, .25, selected)] <= 1.05 for a in five) >= 4
    unused_hold = all(abs(cover[(a, .25, selected)]-cover[(a, .25, "oracle")]) <= .05 and upper[(a, .25, selected)] <= 1.05 for a in unused)
    simulated_cover = sum((cover[(a, f, selected)] >= .92 or cover[(a, f, selected)] >= cover[(a, f, "oracle")]) and cover[(a, f, selected)] <= .98 for a, f in clustered) >= .8*24
    simulated_error = sum(upper[(a, f, selected)] <= 1.05 for a, f in everything) >= .8*38
    derived = "GO" if coverage_holds and accuracy_holds and unused_hold and simulated_cover and simulated_error else "PARTIAL" if coverage_holds else "NO-GO"
    assert derived == verdict, (derived, verdict)
    replayed = 0
    if verdict != "NO-GO":
        with gzip.open(ROOT/"validation/clustered-fix-replay/replay.csv.gz", "rt", newline="", encoding="utf-8") as stream:
            for r in csv.DictReader(stream):
                now = lookup_records[(r["analysis"], float(r["fraction"]), selected, r["replicate"])]
                assert now["status"] == r["status"] and same(now, r, ["estimate", "se", "df", "lower", "upper"]), ("replay", r["analysis"], r["fraction"], r["replicate"])
                replayed += 1
        assert replayed == 30000, replayed
    synthesis = (ROOT/"validation/CLUSTERED_FIX_SYNTHESIS.md").read_text(encoding="utf-8")
    assert f"**{verdict}**" in synthesis, "Write-up does not state the verdict"
    sources = scored+posthoc+(ROOT/"validation/CLUSTERED_FIX_PLAN.md").read_text(encoding="utf-8")+(ROOT/"validation/clustered-fix-replay/README.md").read_text(encoding="utf-8")
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    unsupported = sorted(decimals(re.sub(r"\[doi:[^\]]*\]\([^)]*\)", "", synthesis))-decimals(sources))
    assert not unsupported, ("Numbers in the write-up with no source", unsupported)
    assert decimals("A made-up 12.345 here.")-decimals(sources), "Negative control: a made-up number went undetected"
    print(f"{folder.name}: all {len(records)} planned requests, {len(summary)} recomputed groups, seeds and pool targets recomputed, {matched_first} requests match the first study, {matched_review} match the review, selection and verdict re-derived as {selected} and {verdict}, {replayed} replayed through the package, corrupted-summary control: PASS")


def verify_variance_check(folder):
    """Recount under corrected variances: planned requests, the study's estimates, scored output, reading, write-up numbers and the algebra record."""
    import gzip
    import re
    import subprocess
    import sys
    with gzip.open(folder/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == 24500 and len({(r["analysis"], r["fraction"], r["replicate"]) for r in records}) == 24500, "Missing or duplicate requests"
    for row in open(folder/"plan-hashes.txt", encoding="utf-8").read().splitlines()[:3]:
        digest, name = row.split(None, 1)
        assert hashlib.sha256((ROOT/name.strip()).read_bytes()).hexdigest() == digest, name
    scored = subprocess.run([sys.executable, "-I", str(ROOT/"tools/oracle-variance-score.py")], capture_output=True, text=True, check=True).stdout
    assert scored.strip() == (ROOT/"validation/variance-check-decisions.txt").read_text(encoding="utf-8").strip(), "Scored output differs from the saved decisions"
    assert "CANDIDATES: cr2, cv3, wcr" in scored
    groups = {}
    for r in records:
        groups.setdefault((r["level"], r["analysis"], r["fraction"]), []).append(r)
    real = [g for k, g in groups.items() if k[0] != "simulation"]
    rate = lambda g, field: sum(r[field] == "TRUE" for r in g)/len(g)
    assert sum(rate(g, "cover_cr1s") >= .92 for g in real) == 7, "Reading of the current variance changed"
    results = (ROOT/"validation/VARIANCE_CHECK_RESULTS.md").read_text(encoding="utf-8")
    algebra = (ROOT/"validation/cluster-algebra-check.txt").read_text(encoding="utf-8")
    assert "Failures: 0" in algebra and algebra.count("\nPASS") == 50 and "All 50 scored checks passed" in results
    sources = scored+algebra+(ROOT/"validation/VARIANCE_CHECK_PLAN.md").read_text(encoding="utf-8")+(ROOT/"validation/CLUSTERED_FIX_SYNTHESIS.md").read_text(encoding="utf-8")
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    known = {"5.3"}  # Mean of one Zambia item, whose data are not in this repository.
    unsupported = sorted(decimals(re.sub(r"\]\([^)]*\)", "]", results))-decimals(sources)-known)
    assert not unsupported, ("Numbers in the results with no source", unsupported)
    print(f"{folder.name}: all {len(records)} planned requests, scored output regenerated identically, reading re-derived, numbers in the results traced, algebra record has 50 passes: PASS")


def verify_variance_mi_study(folder):
    """Corrected variances through imputation: planned requests, sources, scored output, the rule's outcome, the replay and the write-up's numbers."""
    import gzip
    import re
    import subprocess
    import sys
    for row in rows(folder/"source-hashes.csv"):
        candidate = ROOT/"validation/variance-mi-source"/row["file"]
        path = candidate if candidate.exists() else ROOT/row["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    with gzip.open(folder/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == 148500 and len({(r["level"], r["analysis"], r["fraction"], r["replicate"], r["imputer"]) for r in records}) == 148500, "Missing or duplicate requests"
    scored = subprocess.run([sys.executable, "-I", str(ROOT/"tools/variance-mi-score.py")], capture_output=True, text=True, check=True).stdout
    assert scored.strip() == (ROOT/"validation/variance-mi-decisions.txt").read_text(encoding="utf-8").strip(), "Scored output differs from the saved decisions"
    posthoc = subprocess.run([sys.executable, "-I", str(ROOT/"tools/variance-mi-posthoc.py")], capture_output=True, text=True, check=True).stdout
    assert posthoc.strip() == (ROOT/"validation/variance-mi-posthoc.txt").read_text(encoding="utf-8").strip(), "Post-hoc output differs from the saved table"
    groups = {}
    for r in records:
        if r["status"] == "ok":
            groups.setdefault((r["level"], r["analysis"], float(r["fraction"]), r["imputer"]), []).append(r)
    items = sorted({k[1] for k in groups if k[0] not in ("simulation", "imputed_covariate")})
    assert len(items) == 11
    cells = sorted({k[:3] for k in groups if k[0] in ("simulation", "imputed_covariate")})
    assert len(cells) == 44
    real = {(k[1], k[2], k[3]): g for k, g in groups.items() if k[0] not in ("simulation", "imputed_covariate")}
    cover = lambda g, v: sum(float(r[f"cover_{v}"]) for r in g)/len(g)
    width = lambda g, v: sum(float(r[f"width_{v}"]) for r in g)/len(g)
    distance = lambda v: sum(abs(cover(real[(a, .25, "random_intercept")], v)-.95) for a in items)/11
    eligible = {}
    for v in ("cr2", "cv3"):
        c25 = [cover(real[(a, .25, "random_intercept")], v) for a in items]
        c50 = [cover(real[(a, .5, "random_intercept")], v) for a in items]
        inside = sum((cover(groups[k+("random_intercept",)], v) >= .92 or cover(groups[k+("random_intercept",)], v) >= cover(groups[k+("oracle",)], v)) and cover(groups[k+("random_intercept",)], v) <= .98 for k in cells)
        linear = sum(cover(real[(a, .25, "linear_mi")], v) <= .98 for a in items) >= 10 and sum(cover(groups[k+("linear_mi",)], v) <= .98 for k in cells) >= .9*44
        wide = all(width(real[(a, .25, "random_intercept")], v)/width(real[(a, .25, "random_intercept")], "cr1s") <= 1.30 for a in items)
        eligible[v] = sum(x >= .92 for x in c25) >= 9 and sum(x <= .98 for x in c25) >= 10 and sum(x >= .92 for x in c50) >= 8 and sum(x <= .98 for x in c50) >= 10 and distance("cr1s")-distance(v) >= .005 and inside >= .8*44 and linear and wide
    assert eligible == {"cr2": False, "cv3": True}, eligible
    assert "OUTCOME: cv3" in scored
    saved = {(r["analysis"], float(r["fraction"]), r["replicate"], r["imputer"]): r for r in records}
    replayed = 0
    with gzip.open(ROOT/"validation/variance-mi-replay/replay.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        for r in csv.DictReader(stream):
            now = saved[(r["analysis"], float(r["fraction"]), r["replicate"], r["imputer"])]
            assert (r["status"] == "ok") == (now["status"] == "ok"), r["analysis"]
            if r["status"] == "ok":
                assert abs(float(r["estimate"])-float(now["estimate"])) < 1e-8 and abs(float(r["se"])-float(now["se_cv3"])) < 1e-8 and abs(float(r["df"])-float(now["df_cv3"])) < .01, (r["analysis"], r["replicate"])
            replayed += 1
    assert replayed == 148500, replayed
    synthesis = (ROOT/"validation/VARIANCE_MI_SYNTHESIS.md").read_text(encoding="utf-8")
    assert "the outcome is **CV3**" in synthesis, "Write-up does not state the outcome"
    sources = scored+posthoc+(ROOT/"validation/VARIANCE_MI_PLAN.md").read_text(encoding="utf-8")+(ROOT/"validation/variance-mi-replay/README.md").read_text(encoding="utf-8")
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    unsupported = sorted(decimals(re.sub(r"\]\([^)]*\)", "]", synthesis))-decimals(sources))
    assert not unsupported, ("Numbers in the write-up with no source", unsupported)
    assert decimals("A made-up 12.345 here.")-decimals(sources), "Negative control: a made-up number went undetected"
    print(f"{folder.name}: all {len(records)} planned requests, sources hashed, scored output and post-hoc table regenerated identically, the outcome re-derived as cv3, {replayed} replayed through the package: PASS")


def verify_outside_methods(folder):
    """Comparison with outside methods: sources, planned requests, scored output, the accuracy failures behind the outcome, the rerun and the write-up's numbers."""
    import gzip
    import math
    import re
    import subprocess
    import sys
    rerun = ROOT/"validation/outside-methods-sensitivity"
    for where in (folder, rerun):
        for row in rows(where/"source-hashes.csv"):
            candidate = ROOT/"validation/outside-methods-source"/row["file"]
            path = candidate if candidate.exists() else ROOT/row["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    with gzip.open(folder/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == 297000 and len({(r["analysis"], r["fraction"], r["replicate"], r["method"]) for r in records}) == 297000, "Missing or duplicate requests"
    with gzip.open(rerun/"replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        replaced = sum(1 for _ in csv.DictReader(stream))
    assert replaced == 58356, replaced
    run = lambda *args: subprocess.run([sys.executable, "-I", *map(str, args)], capture_output=True, text=True, check=True).stdout
    scored = run(ROOT/"tools/outside-methods-score.py", folder)
    assert scored.strip() == (ROOT/"validation/outside-methods-decisions.txt").read_text(encoding="utf-8").strip(), "Scored output differs from the saved decisions"
    assert "OUTCOME: KEEP" in scored
    second = run(ROOT/"tools/outside-methods-sensitivity-score.py")
    assert second.strip() == (ROOT/"validation/outside-methods-sensitivity-decisions.txt").read_text(encoding="utf-8").strip(), "The rerun's score differs from the saved file"
    assert second.startswith("POST HOC") and "OUTCOME: KEEP" in second
    posthoc = run(ROOT/"tools/outside-methods-posthoc.py")
    assert posthoc.strip() == (ROOT/"validation/outside-methods-posthoc.txt").read_text(encoding="utf-8").strip(), "Post-hoc output differs from the saved table"
    # Separately from the scorer: every judged method's accuracy ratio at 25% on real items exceeds the limit in the first run.
    with gzip.open(ROOT/"validation/variance-mi-study/replicates.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        option = {(r["analysis"], r["replicate"]): float(r["estimate"])-float(r["target"]) for r in csv.DictReader(stream) if r["imputer"] == "random_intercept" and r["status"] == "ok" and float(r["fraction"]) == .25 and r["level"] not in ("simulation", "imputed_covariate")}
    items = sorted({k[0] for k in option})
    assert len(items) == 11
    for method in ("audited_only", "difference", "raking", "mice_2l_lmer", "mice_2l_pan"):
        errors = {(r["analysis"], r["replicate"]): float(r["estimate"])-float(r["target"]) for r in records if r["method"] == method and r["status"] == "ok" and float(r["fraction"]) == .25 and r["analysis"] in items}
        logs = []
        for item in items:
            both = [k for k in errors if k[0] == item and k in option]
            logs.append(.5*math.log(sum(errors[k]**2 for k in both)/sum(option[k]**2 for k in both)))
        assert math.exp(sum(logs)/11) > 1.05, method
    control = (folder/"control.txt").read_text(encoding="utf-8")
    assert "660 requests" in control and "status mismatches 0" in control, "The control did not pass"
    assert "390 requests" in (rerun/"control.txt").read_text(encoding="utf-8"), "The rerun's control did not pass"
    synthesis = (ROOT/"validation/OUTSIDE_METHODS_SYNTHESIS.md").read_text(encoding="utf-8")
    assert "the outcome is **KEEP**" in synthesis, "Write-up does not state the outcome"
    sources = scored+second+posthoc+"".join((ROOT/"validation"/name).read_text(encoding="utf-8") for name in ("OUTSIDE_METHODS_PLAN.md", "OUTSIDE_METHODS_SENSITIVITY_PLAN.md", "outside-methods-check.txt", "outside-methods-twophase-check.txt"))
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    unsupported = sorted(decimals(re.sub(r"\]\([^)]*\)", "]", synthesis))-decimals(sources))
    assert not unsupported, ("Numbers in the write-up with no source", unsupported)
    assert decimals("A made-up 12.345 here.")-decimals(sources), "Negative control: a made-up number went undetected"
    print(f"{folder.name}: all {len(records)} planned requests and {replaced} rerun, sources hashed, three scored files regenerated identically, accuracy failures re-derived, controls passed: PASS")


def verify_later_studies():
    """Second health system, rare counts and independently written setups: sources, planned requests, scored output, posthoc numbers and the write-ups' numbers."""
    import gzip
    import re
    import subprocess
    import sys
    def hashes(folder, frozen=None):
        for row in rows(folder/"source-hashes.csv"):
            candidate = (ROOT/"validation"/frozen/row["file"]) if frozen else None
            path = candidate if candidate and candidate.exists() else ROOT/row["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["file"]
    def count(path):
        with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
            return sum(1 for _ in csv.DictReader(stream))
    run = lambda *args: subprocess.run([sys.executable, "-I", *map(str, args)], capture_output=True, text=True, check=True).stdout
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    read = lambda name: (ROOT/"validation"/name).read_text(encoding="utf-8")
    def trace(name, *sources):
        unsupported = sorted(decimals(re.sub(r"\]\([^)]*\)", "]", read(name)))-decimals("".join(read(s) for s in sources)))
        assert not unsupported, (name, "Numbers with no source", unsupported)
    second = ROOT/"validation/second-system-study"
    hashes(second)
    assert count(second/"replicates.csv.gz") == 270000
    scored = run(ROOT/"tools/second-system-score.py", second)
    assert scored.strip() == read("second-system-decisions.txt").strip() and "OUTCOME: MIXED" in scored
    assert "the outcome is **MIXED**" in read("SECOND_SYSTEM_SYNTHESIS.md")
    trace("SECOND_SYSTEM_SYNTHESIS.md", "second-system-decisions.txt", "second-system-posthoc.txt", "SECOND_SYSTEM_PLAN.md")
    rare = ROOT/"validation/rare-counts-study"
    hashes(rare)
    assert (count(rare/"part1.csv.gz"), count(rare/"part2.csv.gz"), count(rare/"part3.csv.gz")) == (18000, 21000, 20500)
    scored = run(ROOT/"tools/rare-counts-score.py", rare)
    assert scored.strip() == read("rare-counts-decisions.txt").strip()
    assert all(k in scored for k in ("OUTCOME P1: HEAVY TAILS", "OUTCOME P2: none eligible", "OUTCOME P3: none eligible", "OUTCOME P4: NOT ADOPTED"))
    trace("RARE_COUNTS_SYNTHESIS.md", "rare-counts-decisions.txt", "independent-setups-posthoc.txt", "RARE_COUNTS_PLAN.md", "second-system-posthoc.txt")
    setups = ROOT/"validation/independent-setups-study"
    hashes(setups, "independent-setups-source")
    assert count(setups/"replicates.csv.gz") == 148500 and count(setups/"zambia.csv.gz") == 54000
    scored = run(ROOT/"tools/independent-setups-score.py", setups)
    assert scored.strip() == read("independent-setups-decisions.txt").strip() and "OUTCOME: KEEP" in scored and "Control: this scorer reproduces" in scored
    assert run(ROOT/"tools/independent-setups-posthoc.py").strip() == read("independent-setups-posthoc.txt").strip()
    trace("INDEPENDENT_SETUPS_SYNTHESIS.md", "independent-setups-decisions.txt", "independent-setups-posthoc.txt", "INDEPENDENT_SETUPS_PLAN.md", "second-system-decisions.txt", "outside-methods-decisions.txt", "independent-setups/agent-reports.json")
    assert decimals("A made-up 12.345 here.")-decimals(read("independent-setups-decisions.txt")), "Negative control: a made-up number went undetected"
    slopes = ROOT/"validation/within-between-study"
    hashes(slopes)
    assert (count(slopes/"guard.csv.gz"), count(slopes/"zambia.csv.gz"), count(slopes/"deciding.csv.gz")) == (49500, 36000, 444000)
    scored = run(ROOT/"tools/within-between-score.py", slopes)
    assert scored.strip() == read("within-between-decisions.txt").strip() and "OUTCOME: NOT ADOPTED" in scored and "Control: this scorer reproduces" in scored
    assert "the outcome is **NOT ADOPTED**" in read("WITHIN_BETWEEN_SYNTHESIS.md")
    assert run(ROOT/"tools/within-between-posthoc.py", slopes).strip() == read("within-between-posthoc.txt").strip()
    assert run(ROOT/"tools/within-between-shared-cells.py", slopes).strip() == read("within-between-shared-cells.txt").strip()
    trace("WITHIN_BETWEEN_SYNTHESIS.md", "within-between-decisions.txt", "within-between-posthoc.txt", "within-between-shared-cells.txt", "WITHIN_BETWEEN_PLAN.md")
    print("Later studies: second health system, rare counts, independently written setups and separate slopes, all planned requests, sources hashed, scored files regenerated identically, numbers traced: PASS")


def verify_learn_pages():
    """The guide and simulator in docs/learn: every decimal in their visible text has a source, the guide's
    diagram matches docs/figures/how-it-works.svg, and, where node is installed, the simulator's engine still
    matches the R package's saved results."""
    import html
    import re
    import shutil
    import subprocess
    learn = ROOT/"docs/learn"
    def visible(path):
        text = path.read_text(encoding="utf-8")
        text = re.sub(r"<(script|style)\b.*?</\1>", " ", text, flags=re.S)
        text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
        return html.unescape(re.sub(r"<[^>]+>", " ", text))
    decimals = lambda text: set(re.findall(r"(?<![\w.])\d+\.\d+(?![\d])", text))
    read = lambda *paths: "".join((ROOT/p).read_text(encoding="utf-8") for p in paths)
    # The guide claims its decimals come from the README or the example's output, so only those count for it.
    for page, sources in (("index.html", read("README.md", "docs/learn/qa/guide-example.txt")),
                          ("simulator.html", read("README.md", "docs/learn/qa/cross-check.txt", "docs/learn/audit-sim.js"))):
        unsupported = sorted(decimals(visible(learn/page)) - decimals(sources))
        assert not unsupported, (page, "Numbers with no source", unsupported)
    assert decimals("a made-up 12.345")-decimals(read("README.md", "docs/learn/qa/guide-example.txt")), "Negative control: a made-up number went undetected"
    guide = (learn/"index.html").read_text(encoding="utf-8")
    inline = re.search(r"<!--HOW-FIGURE start[^>]*-->(.*?)<!--HOW-FIGURE end-->", guide, re.S).group(1)
    figure = re.sub(r"\s*<style>.*?</style>", "", (ROOT/"docs/figures/how-it-works.svg").read_text(encoding="utf-8"), flags=re.S)
    shapes = lambda svg: re.findall(r"<(?:rect|line|circle|path|text)\b[^>]*>(?:[^<]*)", re.sub(r'\s(?:id|aria-labelledby)="[^"]*"', "", svg))
    assert shapes(inline) == shapes(figure) and shapes(figure), "The guide's diagram differs from docs/figures/how-it-works.svg"
    example = (learn/"qa/guide-example.txt").read_text(encoding="utf-8")
    assert "0.126, 95% interval 0.097 to 0.155" in example, "The guide's example no longer matches the README"
    node = shutil.which("node")
    if node:
        result = subprocess.run([node, str(learn/"qa/compare.js")], capture_output=True, text=True)
        assert result.returncode == 0 and result.stdout == (learn/"qa/cross-check.txt").read_text(encoding="utf-8"), "Simulator engine check changed or failed"
        engine = "engine rechecked against the R package's saved results"
    else:
        engine = "engine check NOT RUN, node is not installed"
    print(f"Learning pages: every decimal traced, diagram matches the figure, {engine}: PASS")


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
    verify_real_audit_study(ROOT/"validation/real-audit-study")
    verify_clustered_study(ROOT/"validation/clustered-study")
    verify_clustered_fix_study(ROOT/"validation/clustered-fix-study")
    verify_variance_check(ROOT/"validation/variance-check")
    verify_variance_mi_study(ROOT/"validation/variance-mi-study")
    verify_outside_methods(ROOT/"validation/outside-methods-study")
    verify_later_studies()
    verify_learn_pages()
    print("This verifies correspondence and repeatability, not inferential validity.")


if __name__ == "__main__":
    main()
