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
    print("This verifies correspondence and repeatability, not inferential validity.")


if __name__ == "__main__":
    main()
