#!/usr/bin/env python3
"""Check every number in the README front page, case study, PDF and blurbs against its source, then write FACT_CHECK.md.

Usage: python3 docs/case-study/fact_check.py
Each claim recomputes its value from a source file, usually a frozen results table, which is only read. Any digit or number
word in a document that no claim covers fails the check, and so does a claim that no longer appears in its document.
Also checks relative links, the README example (by running it), the PDF page count and wording rules. Exit 0 means all passed.
"""
import csv, datetime, json, pathlib, re, statistics, subprocess, sys, tempfile

sys.dont_write_bytecode = True

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import build_case_study as build  # noqa: E402

README, CASE, BLURBS = ROOT / "README.md", HERE / "CASE_STUDY.md", HERE / "BLURBS.md"
RD, RP, RC, RS = [f"validation/real-data-study/{x}.csv" for x in ("summary", "paired-comparisons", "corruption-scores", "request-status")]
SY = "validation/contamination-study/summary.csv"
AP, AK, AV, AR = [f"validation/audit-plan-study/{x}.csv" for x in ("summary", "pool-k", "verdict", "request-status")]
MV = "validation/merged-audit-pilot/study/verdict.csv"
PT, RSET, ASET = "validation/real-data-study/pool-targets.csv", "validation/real-data-study/settings.json", "validation/audit-plan-study/settings.json"
SYN_R, SYN_A, VAL, ASSESS = "validation/REAL_DATA_SYNTHESIS.md", "validation/AUDIT_PLAN_SYNTHESIS.md", "validation/VALIDATION.md", "docs/ASSESSMENT.md"
RA, RR, RE, DEC, SYN_T = "validation/real-audit-study/summary.csv", "validation/real-audit-study/recovery-ratios.csv", "validation/real-audit-review/rule1-extension.csv", "validation/real-audit-decisions.txt", "validation/REAL_AUDIT_SYNTHESIS.md"
BIASED = ("pnc", "iron", "anc4")

_cache = {}


def rows(path):
    return _cache.setdefault(path, list(csv.DictReader((ROOT / path).open())))


def one(path, **kw):
    hits = [r for r in rows(path) if all(r[k] == str(v) for k, v in kw.items())]
    assert len(hits) == 1, (path, kw, len(hits))
    return hits[0]


def rd(family, scenario, method, col):
    return float(one(RD, family=family, scenario=scenario, method=method)[col])


def sy(family, scenario, method, col="coverage_conditional"):
    return float(one(SY, family=family, scenario=scenario, method=method)[col])


def has(path, regex, shown):
    return shown if re.search(regex, (ROOT / path).read_text()) else f"NOT FOUND: {regex}"


f3, f4, f2 = (lambda v: f"{v:.3f}"), (lambda v: f"{v:.4f}"), (lambda v: f"{v:.2f}")
s4 = lambda v: f"{v:+.4f}"
pct0, pct1 = (lambda v: f"{100 * v:.0f}%"), (lambda v: f"{100 * v:.1f}%")
GOOD = lambda family, name: ("cole_mi", family, name)
PRIMARY = ("clean", "flagged_errors", "missed_flags", "site_shift_declared")
REAL_CELLS = ("ctgov", "nhanes")


def narrower():
    return [1 - sy("primary", c, "contamination", "mean_width") / sy("primary", c, "reference_normal", "mean_width") for c in PRIMARY]


def plan_errors():
    return {(r["cell"], r["pilot_records"], float(r["target_fraction"])): float(r["median_error"]) for r in rows(AP) if r["cell"] in REAL_CELLS}


def auc(family):
    return statistics.mean(float(r["auc"]) for r in rows(RC) if r["family"] == family and r["scenario"] == "audit_25pct" and r["method"] == "contamination" and r["auc"] != "NA")


def screen_passes():
    n = 0
    for c in PRIMARY:
        r = one(SY, family="primary", scenario=c, method="contamination")
        n += float(r["coverage_wilson_lower"]) >= .925 and float(r["coverage_wilson_upper"]) <= .975
    return n


def model_vs_linear():
    """Paired MSE z of the contamination model against linear MI in every real-data cell. Negative favors the model."""
    return {(r["family"], r["scenario"]): float(r["mse_difference"]) / float(r["mse_mcse"]) for r in rows(RP)
            if r["method"] == "contamination" and r["comparator"] == "cole_mi"}


def worded(condition, shown, detail):
    return shown if condition else f"FALSE: {detail}"


def width_ratio(family, method):
    return rd(family, "audit_25pct", "reference_normal", "mean_width") / rd(family, "audit_25pct", method, "mean_width")


def audit_pct(family, scenario):
    return rd(family, scenario, "cole_mi", "mean_available_n") / json.loads((ROOT / RSET).read_text())["n"]


def ra(analysis, method, col="coverage_cluster"):
    """Rwanda facility-months at a 25% audit, from the real-audit study."""
    return float(one(RA, level="rwanda_month", analysis=analysis, fraction="0.25", method=method)[col])


def ratio_25(analysis, method, comparator):
    return one(RR, level="rwanda_month", analysis=analysis, fraction="0.25", method=method, comparator=comparator)


def mi_gaps():
    return [ra(a, "oracle") - ra(a, "linear_mi") for a in ("pnc", "syphilis", "iron", "anc4")]


EXAMPLE = {}


def run_example():
    """Run the README's second R block exactly, loading the package from this tree."""
    front = README.read_text().split("<!-- front-page-end")[0]
    code = re.findall(r"```r\n(.*?)```", front, re.S)[1].replace("library(catmech)", f'suppressMessages(pkgload::load_all("{ROOT}", quiet = TRUE))')
    with tempfile.NamedTemporaryFile("w", suffix=".R") as script:
        script.write(code)
        script.flush()
        out = subprocess.run(["Rscript", script.name], cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr[-800:]
    m = re.search(r"recorded_kg\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", out.stdout)
    EXAMPLE.update(estimate=float(m.group(1)), lower=float(m.group(4)), upper=float(m.group(5)))
    EXAMPLE["warnings"] = "Warning" in out.stderr
    return out.stdout


# Each claim: text as shown, documents using it (R README front page, C case study and PDF, B blurbs, L README lower sections
# changed in this rewrite), what it says, source file, and a function returning the text the source supports.
claims = [
    ("35%", "RC", "Real audit, 25%: lowest clustered coverage of the contamination model on slopes with a material gap", RA, lambda: pct0(min(ra(a, "contamination") for a in BIASED))),
    ("72%", "RC", "Real audit, 25%: highest clustered coverage of the contamination model on slopes with a material gap", RA, lambda: pct0(max(ra(a, "contamination") for a in BIASED))),
    ("90%", "RC", "Real audit: lowest clustered coverage of a one-item recount on those slopes", RA, lambda: pct0(min(ra(a, "oracle") for a in BIASED))),
    ("97%", "RC", "Real audit: highest clustered coverage of a one-item recount on those slopes", RA, lambda: pct0(max(ra(a, "oracle") for a in BIASED))),
    ("79.5%", "R", "Real audit, 25%: worst clustered coverage of linear MI across the four facility-month slopes", RA, lambda: pct1(min(ra(a, "linear_mi") for a in ("pnc", "syphilis", "iron", "anc4")))),
    ("1.4", "R", "Real audit, 25%: RMSE against the pool, linear MI over trusting reported values, syphilis", RA, lambda: f'{ra("syphilis", "linear_mi", "rmse") / ra("syphilis", "recorded", "rmse"):.1f}'),
    ("2.6", "R", "Real audit, 25%: recovery RMSE on 2,000 fresh cohorts, linear MI over trusting reported values, syphilis (2.4 in the study's 500)", RE, lambda: f'{float(one(RE, analysis="syphilis", fraction="0.25")["recovery_ratio"]):.1f}'),
    ("9", "R", "Real audit, 25%: smallest gap between linear MI and one-item recount clustered coverage, points", RA, lambda: f"{100 * min(mi_gaps()):.0f}"),
    ("11", "R", "Real audit, 25%: largest gap between linear MI and one-item recount clustered coverage, points", RA, lambda: f"{100 * max(mi_gaps()):.0f}"),
    ("two of three", "RC", "Real audit, 25%: exact-or-wrong slopes where the model's error was clearly below linear MI's", RR, lambda: worded(sum(float(ratio_25(a, "contamination", "linear_mi")["upper"]) < 1 for a in ("pnc", "syphilis", "iron")) == 2, "two of three", "not two slopes")),
    ("0.4.1", "RC", "Package version", "DESCRIPTION", lambda: re.search(r"^Version: (.+)$", (ROOT / "DESCRIPTION").read_text(), re.M).group(1)),
    ("November 1, 2026", "CB", "Status date", ASSESS, lambda: has(ASSESS, r"Assessed 2026-11-01", "November 1, 2026")),
    ("Version 0.4", "C", "Release that added the recording-error model", "NEWS.md", lambda: has("NEWS.md", r"# catmech 0\.4\.0\n\n- Add `catmech_contamination\(\)`", "Version 0.4")),
    ("25%", "RCB", "Audit fraction: mean referenced rows over cohort size, ClinicalTrials.gov cell", RD, lambda: pct0(audit_pct("ctgov", "audit_25pct"))),
    ("25%", "RC", "Audit fraction: mean referenced rows over cohort size, NHANES cell", RD, lambda: pct0(audit_pct("nhanes", "audit_25pct"))),
    ("0.0286", "RC", "RMSE, linear MI, ClinicalTrials.gov, 25%", RD, lambda: f4(rd("ctgov", "audit_25pct", "cole_mi", "rmse"))),
    ("0.0308", "RC", "RMSE, contamination model, ClinicalTrials.gov, 25%", RD, lambda: f4(rd("ctgov", "audit_25pct", "contamination", "rmse"))),
    ("0.0539", "C", "RMSE, discarding recorded values, ClinicalTrials.gov, 25%", RD, lambda: f4(rd("ctgov", "audit_25pct", "reference_normal", "rmse"))),
    ("0.0275", "C", "RMSE, all references, ClinicalTrials.gov, 25%", RD, lambda: f4(rd("ctgov", "audit_25pct", "oracle", "rmse"))),
    ("+0.0001", "R", "Bias, linear MI, ClinicalTrials.gov, 25%", RD, lambda: s4(rd("ctgov", "audit_25pct", "cole_mi", "bias"))),
    ("+0.0068", "R", "Bias, contamination model, ClinicalTrials.gov, 25%", RD, lambda: s4(rd("ctgov", "audit_25pct", "contamination", "bias"))),
    ("0.0409", "RC", "RMSE, linear MI, NHANES, 25%", RD, lambda: f4(rd("nhanes", "audit_25pct", "cole_mi", "rmse"))),
    ("0.0408", "RC", "RMSE, contamination model, NHANES, 25%", RD, lambda: f4(rd("nhanes", "audit_25pct", "contamination", "rmse"))),
    ("0.0785", "C", "RMSE, discarding recorded values, NHANES, 25%", RD, lambda: f4(rd("nhanes", "audit_25pct", "reference_normal", "rmse"))),
    ("0.0357", "C", "RMSE, all references, NHANES, 25%", RD, lambda: f4(rd("nhanes", "audit_25pct", "oracle", "rmse"))),
    ("10%", "R", "Audit fraction: mean referenced rows over cohort size, NHANES 10% cell", RD, lambda: pct0(audit_pct("nhanes", "audit_10pct"))),
    ("0.0435", "R", "RMSE, contamination model, NHANES, 10%", RD, lambda: f4(rd("nhanes", "audit_10pct", "contamination", "rmse"))),
    ("0.0490", "R", "RMSE, linear MI, NHANES, 10%", RD, lambda: f4(rd("nhanes", "audit_10pct", "cole_mi", "rmse"))),
    ("14 of 200", "R", "Model requests with no interval, NHANES, 10%: requested minus returned", RD,
     lambda: f'{int(rd("nhanes", "audit_10pct", "contamination", "requested") - rd("nhanes", "audit_10pct", "contamination", "returned_intervals"))} of {int(rd("nhanes", "audit_10pct", "contamination", "requested"))}'),
    ("0.54", "RC", "Mean AUC for finding wrong records, ClinicalTrials.gov, 25%", RC, lambda: f2(auc("ctgov"))),
    ("0.50", "RC", "Mean AUC for finding wrong records, NHANES, 25%", RC, lambda: f2(auc("nhanes"))),
    ("about half", "R", "Linear MI width over discard width, 25%: must be 0.45 to 0.55 in both datasets", RD,
     lambda: worded(all(.45 <= 1 / width_ratio(f, "cole_mi") <= .55 for f in REAL_CELLS), "about half", "width ratio")),
    ("about twice", "R", "Discard width over linear MI width, 25%: must be 1.8 to 2.2 in both datasets", RD,
     lambda: worded(all(1.8 <= width_ratio(f, "cole_mi") <= 2.2 for f in REAL_CELLS), "about twice", "width ratio")),
    ("0.102", "R", "Mean width, linear MI, ClinicalTrials.gov, 25%", RD, lambda: f3(rd("ctgov", "audit_25pct", "cole_mi", "mean_width"))),
    ("0.205", "R", "Mean width, discarding recorded values, ClinicalTrials.gov, 25%", RD, lambda: f3(rd("ctgov", "audit_25pct", "reference_normal", "mean_width"))),
    ("0.154", "R", "Mean width, linear MI, NHANES, 25%", RD, lambda: f3(rd("nhanes", "audit_25pct", "cole_mi", "mean_width"))),
    ("0.297", "R", "Mean width, discarding recorded values, NHANES, 25%", RD, lambda: f3(rd("nhanes", "audit_25pct", "reference_normal", "mean_width"))),
    ("91.8%", "RC", "Coverage, linear MI, ClinicalTrials.gov, 25%", RD, lambda: pct1(rd("ctgov", "audit_25pct", "cole_mi", "coverage_conditional"))),
    ("91.2%", "RC", "Coverage, all references, ClinicalTrials.gov, 25%", RD, lambda: pct1(rd("ctgov", "audit_25pct", "oracle", "coverage_conditional"))),
    ("91.2%", "C", "Coverage, contamination model, ClinicalTrials.gov, 25%", RD, lambda: pct1(rd("ctgov", "audit_25pct", "contamination", "coverage_conditional"))),
    ("93.9%", "RC", "Coverage, linear MI, NHANES, 25%", RD, lambda: pct1(rd("nhanes", "audit_25pct", "cole_mi", "coverage_conditional"))),
    ("94.5%", "C", "Coverage, contamination model, NHANES, 25%", RD, lambda: pct1(rd("nhanes", "audit_25pct", "contamination", "coverage_conditional"))),
    ("94.9%", "RC", "Coverage, all references, NHANES, 25%", RD, lambda: pct1(rd("nhanes", "audit_25pct", "oracle", "coverage_conditional"))),
    ("64.7%", "RC", "Coverage, linear MI, synthetic gross-and-small errors", SY, lambda: pct1(sy("wrong_model", "error_gross_and_small", "cole_mi"))),
    ("94.0%", "RC", "Coverage, contamination model, synthetic gross-and-small errors", SY, lambda: pct1(sy("wrong_model", "error_gross_and_small", "contamination"))),
    ("91.8%", "R", "Coverage, contamination model, synthetic noisy reference with declared SD", SY, lambda: pct1(sy("credibility", "noisy_declared", "contamination"))),
    ("87.3%", "R", "Coverage, contamination model, synthetic heavy-tailed target", SY, lambda: pct1(sy("wrong_model", "heavy_tail_target", "contamination"))),
    ("79.3%", "R", "Coverage, contamination model, synthetic selection on true value with errors", SY, lambda: pct1(sy("credibility", "nonignorable_missed", "contamination"))),
    ("7.3%", "R", "Coverage, reference normal, synthetic selective references", SY, lambda: pct1(sy("credibility", "nonignorable_missed", "reference_normal"))),
    ("52.7%", "R", "Coverage, reference normal, synthetic noisy references", SY, lambda: pct1(sy("credibility", "noisy_declared", "reference_normal"))),
    ("10%", "RC", "Planner's limit on median width error", SYN_A, lambda: has(SYN_A, r"The limit is 10%", "10%")),
    ("5%", "R", "Smallest target audit in the planner test", AP, lambda: pct0(min(k[2] for k in plan_errors()))),
    ("11.3%", "RC", "Lowest planner median error at the 5% target, of four real-data cells", AP, lambda: pct1(min(v for k, v in plan_errors().items() if k[2] == .05))),
    ("16.8%", "RC", "Highest planner median error at the 5% target, of four real-data cells", AP, lambda: pct1(max(v for k, v in plan_errors().items() if k[2] == .05))),
    ("7.8%", "RC", "Largest planner median error at any larger target, real-data cells", AP, lambda: pct1(max(v for k, v in plan_errors().items() if k[2] > .05))),
    ("30 of 600 records", "C", "Records at the 5% target, cohort size", ASET, lambda: (lambda c: f'{c["target_records"][0]} of {c["n"]} records')(json.loads((ROOT / ASET).read_text()))),
    ("NO-GO", "RC", "Planner verdict", AV, lambda: rows(AV)[0]["verdict"]),
    ("NO-GO", "RCB", "Merged-idea verdict", MV, lambda: one(MV, regime="moderate")["verdict"]),
    ("set.seed(1)", "R", "Example parameter, not a finding. Run by this check.", "README.md", lambda: "set.seed(1)"),
    ("nrow(d) / 4", "R", "Example parameter: a 25% audit. Run by this check.", "README.md", lambda: "nrow(d) / 4"),
    ("agreement_tol = 1", "R", "Example parameter: the 1 kg tolerance the frozen NHANES audit used", SYN_R, lambda: has(SYN_R, r"67\.3% beyond 1 kg", "agreement_tol = 1")),
    ("m = 20", "RC", "Imputations per fit in the real-data study", RSET, lambda: f'm = {json.loads((ROOT / RSET).read_text())["m"]}'),
    ("seed = 123", "R", "Example parameter, not a finding. Run by this check.", "README.md", lambda: "seed = 123"),
    ("0.126", "R", "Example output: recorded_kg estimate, one random audit", "README.md, by running it", lambda: f3(EXAMPLE["estimate"])),
    ("95% interval", "R", "Default level of catmech_pool_lm", "R/analysis.R", lambda: has("R/analysis.R", r"catmech_pool_lm <- function\(result, formula, conf_level = 0\.95", "95% interval")),
    ("0.097", "R", "Example output: lower limit", "README.md, by running it", lambda: f3(EXAMPLE["lower"])),
    ("0.155", "R", "Example output: upper limit", "README.md, by running it", lambda: f3(EXAMPLE["upper"])),
    ("0.127", "R", "NHANES evaluation pool's coefficient using every measured weight", PT, lambda: f3(float(one(PT, dataset="nhanes")["target_coefficient"]))),
    ("1,000", "RC", "Cohorts at a 25% audit", RD, lambda: f'{int(rd("ctgov", "audit_25pct", "cole_mi", "requested")):,}'),
    ("600", "RC", "Records per cohort", RSET, lambda: str(json.loads((ROOT / RSET).read_text())["n"])),
    ("v0.3", "R", "Figure label for the old default, which discards recorded values", ASSESS, lambda: has(ASSESS, r"Version 0\.3 treated recorded values as unavailable", "v0.3")),
    ("two public datasets", "RCB", "Real-data pools in the study", AK, lambda: worded(len(rows(AK)) == 2, "two public datasets", "pool count")),
    ("30%", "C", "Smallest interval reduction, model against old default, four primary synthetic cells", SY, lambda: pct0(min(narrower()))),
    ("38%", "C", "Largest interval reduction, model against old default, four primary synthetic cells", SY, lambda: pct0(max(narrower()))),
    ("three of four", "C", "Primary synthetic cells whose Wilson interval is inside the screen 0.925 to 0.975", SY, lambda: worded(screen_passes() == 3, "three of four", f"{screen_passes()} passes")),
    ("in the fourth", "C", "Site-shift cell is the one that failed the screen", SY, lambda: worded(float(one(SY, family="primary", scenario="site_shift_declared", method="contamination")["coverage_wilson_lower"]) < .925, "in the fourth", "screen")),
    ("nine real-data scenarios", "C", "Paired contrasts of the model against linear MI", RP, lambda: worded(len(model_vs_linear()) == 9, "nine real-data scenarios", "cell count")),
    ("never clearly better", "C", "Most favorable paired MSE z for the model across the nine cells must be above -2", RP, lambda: worded(min(model_vs_linear().values()) > -2, "never clearly better", "paired z")),
    ("within Monte Carlo noise", "R", "Paired MSE z for the model, NHANES 10%: must be between -2 and 0, and the most favorable cell", RP,
     lambda: worded(-2 < model_vs_linear()[("nhanes", "audit_10pct")] < 0 and min(model_vs_linear(), key=model_vs_linear().get) == ("nhanes", "audit_10pct"), "within Monte Carlo noise", "paired z")),
    ("nearly halved", "C", "RMSE reduction, linear MI against discarding, 25%: must be 40% to 50% in both datasets", RD,
     lambda: worded(all(.4 <= 1 - rd(f, "audit_25pct", "cole_mi", "rmse") / rd(f, "audit_25pct", "reference_normal", "rmse") < .5 for f in REAL_CELLS), "nearly halved", "RMSE ratio")),
    ("95% ranges", "C", "Nominal level of the intervals scored for coverage", "validation/CONTAMINATION_PLAN.md", lambda: has("validation/CONTAMINATION_PLAN.md", r"95%|0\.95", "95% ranges")),
    ("two of five", "C", "Advance real-data predictions that held", VAL, lambda: has(VAL, r"Two of five advance predictions held", "two of five")),
    ("Two more ideas", "C", "Merged idea and planner, both NO-GO", MV, lambda: worded(rows(AV)[0]["verdict"] == "NO-GO" and one(MV, regime="moderate")["verdict"] == "NO-GO", "Two more ideas", "verdicts")),
    ("Two other ideas", "B", "Merged idea and planner, both NO-GO", MV, lambda: worded(rows(AV)[0]["verdict"] == "NO-GO" and one(MV, regime="moderate")["verdict"] == "NO-GO", "Two other ideas", "verdicts")),
    ("27,208", "C", "ClinicalTrials.gov evaluation pool rows", AK, lambda: f'{int(one(AK, dataset="ctgov")["pool_rows"]):,}'),
    ("3,908", "C", "NHANES evaluation pool rows", AK, lambda: f'{int(one(AK, dataset="nhanes")["pool_rows"]):,}'),
    ("1.31e-4", "C", "Paired MSE difference, model minus linear MI, ClinicalTrials.gov, 25%", RP, lambda: f'{float(one(RP, family="ctgov", scenario="audit_25pct", method="contamination", comparator="cole_mi")["mse_difference"]):.2e}'.replace("e-04", "e-4")),
    ("2.2e-5", "C", "MCSE of that difference", RP, lambda: f'{float(one(RP, family="ctgov", scenario="audit_25pct", method="contamination", comparator="cole_mi")["mse_mcse"]):.1e}'.replace("e-05", "e-5")),
    ("86.7%", "L", "Moved text: synthetic standard reference PMM coverage, 25% references", SY, lambda: pct1(sy("primary", "missed_flags", "reference_pmm"))),
    ("95.9%", "L", "Moved text: synthetic bootstrap-donor PMM coverage, 25% references", SY, lambda: pct1(sy("primary", "missed_flags", "reference_pmm_abb"))),
    ("51%", "L", "Moved text: bootstrap-donor widths against reference normal on real data, largest increase. REAL_DATA_SYNTHESIS.md says 52%, from rounded table values", RD,
     lambda: pct0(max(rd(f, s, "reference_pmm_abb", "mean_width") / rd(f, s, "reference_normal", "mean_width") - 1 for f, s in (("ctgov", "audit_10pct"), ("ctgov", "audit_25pct"), ("ctgov", "audit_50pct"), ("nhanes", "audit_10pct"), ("nhanes", "audit_25pct"), ("nhanes", "audit_50pct"))))),
    ("about 380 of 600", "L", "New text: mean records by formula when discarding, both pools (383 and 378)", AK,
     lambda: worded(all(375 <= 600 * float(r["discard_k"]) / (1.25 ** 2 - 1 + float(r["discard_k"])) <= 385 for r in rows(AK)), "about 380 of 600", "discard records")),
    ("about 40 to 70", "L", "New text: records by formula when keeping recorded values (41 and 66)", AK,
     lambda: worded(all(40 <= 600 * float(r["k"]) / (1.25 ** 2 - 1 + float(r["k"])) <= 70 for r in rows(AK)), "about 40 to 70", "records")),
]

# Number-like text that states no finding. Listed so the coverage check can pass, and so a reader can see what was excluded.
not_findings = [
    ("R", "This one audit", "Describes the single illustrative run in the example"),
    ("R", "# 1.", "Step number in a code comment"), ("R", "# 2.", "Step number in a code comment"),
    ("C", "two poor options", "Framing of the problem, not a measurement"),
    ("C", "until one passed", "Ordinary phrase, not a count"),
    ("B", "a checked one", "Ordinary phrase, not a count"),
    ("B", "One-line", "Heading from the brief"), ("B", "Three-line", "Heading from the brief"), ("B", "Two-sentence", "Heading from the brief"),
]


BOUND_L, BOUND_R = r"(?<![\d.,])", r"(?!\d|\.\d)"


# Statements without a number, each checked by grep against the record that supports it.
statements = [
    ("The real-audit rule on the contamination model came out inconclusive", DEC, r"Outcome: INCONCLUSIVE"),
    ("The README says so", "README.md", r"came out inconclusive at a 25% audit"),
    ("The README says the model had lower error on two of three slopes, not the reverse", "README.md", r"model had lower error than linear MI on two of three slopes"),
    ("The README says the model's coverage was the lower one", "README.md", r"its intervals covered 35% to 72%, against 90% to 97% for the recount"),
    ("No dataset is a patient-level clinical source-data audit", "README.md", r"No real dataset here is a patient-level clinical source-data audit"),
    ("The write-up says what the target is", SYN_T, r"Every \"full audit\" here recounts one item and keeps the covariate as reported"),
    ("On ClinicalTrials.gov both values are sponsor-entered, so the reference is not independent truth", SYN_R, r"both are entered by sponsors[\s\S]*reference is not independent truth"),
    ("Simulation data are synthetic, and the primary cells favor the model", "validation/CONTAMINATION_SYNTHESIS.md", r"favor it by construction[\s\S]*All data are synthetic"),
    ("No new statistical method is claimed", ASSESS, r"does not support a novel-method"),
    ("Not on CRAN, with no outside peer review or clinical validation", "docs/GITHUB.md", r"A CRAN submission, DOI, independent scientific peer review and\s+clinical validation are separate future steps"),
    ("Linear MI is the recommended start and a supported option", ASSESS, r"`catmech_reference\(use_recorded = TRUE\)` and is the recommended starting point"),
    ("The model could not tell which records were wrong", SYN_R, r"cannot tell \*which\* record is wrong"),
    ("Every study had a written plan first, with earlier looks disclosed", VAL, r"Plans preceded their runs[\s\S]*Development runs and the changes they caused are disclosed[\s\S]*\[plan\]\(AUDIT_PLAN_PLAN\.md\) preceded its evaluation run"),
    ("Results were not retuned after viewing", SYN_A, r"nothing was retried or retuned after viewing"),
    ("Coverage was judged against the oracle", SYN_R, r"coverage is judged against the oracle"),
    ("Do not write the merged paper on that thesis", ASSESS, r"A merged paper on that thesis is not recommended"),
    ("Merged idea: audits helped, targeting did not beat the most suspicious records, direction added nothing", ASSESS, r"audits helped, but network targeting did not beat auditing the most suspicious records, and direction added nothing"),
    ("The planner is not validated", SYN_A, r"should not be described as validated"),
    ("The planner test was not fully blind", SYN_A, r"prospective for the planner's code, and not fully blind for the pools"),
    ("A bookkeeping bug and missed run-time estimates are on the record", VAL, r"bookkeeping bug[\s\S]*Run times exceeded my estimates"),
    ("Plain Python recomputes the reported summaries", VAL, r"Plain Python recomputes all summaries"),
    ("AI use and its limits", "CONTRIBUTING.md", r"OpenAI Codex[\s\S]*Anthropic Claude Code[\s\S]*No AI system is listed as an author[\s\S]*not external peer review"),
]


def self_test():
    """Negative controls: each deliberate fault must be caught, or the check itself has failed."""
    misses = []
    phrases = [c[0] for c in claims if "R" in c[1]]
    if not any("0.0296" in hit for hit in uncovered("RMSE 0.0296 here", phrases)):
        misses.append("a changed number")
    if not any("12.3%" in hit for hit in uncovered("a new 12.3% appears", phrases)):
        misses.append("an added number")
    if not uncovered("about seven records", phrases):
        misses.append("an added number word")
    if f4(rd("ctgov", "audit_25pct", "contamination", "rmse")) == "0.0296":
        misses.append("a wrong claim value")
    if not wording_rules({"test": "The package is published and has many users."}):
        misses.append("an unnegated publication claim")
    if not wording_rules({"test": "A third party used this."}):
        misses.append("a third-party reference")
    probe = HERE / "_link_probe.md"
    try:
        probe.write_text("[broken](no-such-file.md) and [anchor](../../README.md#no-such-heading)\n")
        if len(check_links([probe])[1]) != 2:
            misses.append("broken links")
    finally:
        probe.unlink(missing_ok=True)
    return misses


def mask(text, phrases):
    """Remove claimed phrases, longest first, respecting number boundaries, and return what is left."""
    for phrase in sorted(phrases, key=len, reverse=True):
        text = re.sub(BOUND_L + re.escape(phrase) + BOUND_R, " ", text, flags=re.I)
    return text


NUMBER_WORDS = r"\b(?:one|two|three|four|five|six|seven|eight|nine|ten|fourth|half|halved|twice)\b"


def uncovered(text, phrases):
    rest = mask(text, phrases)
    hits = re.findall(r"\S*\d\S*", rest) + [rest[max(0, m.start() - 20):m.end() + 20].replace("\n", " ") for m in re.finditer(NUMBER_WORDS, rest, re.I)]
    return hits


def clean(md):
    """Drop link targets and HTML comments. Link text, image alt text and code stay in."""
    return re.sub(r"\]\([^)]*\)", "]", re.sub(r"<!--.*?-->", " ", md, flags=re.S))


def slug(heading):
    return re.sub(r"[^\w\- ]", "", heading.lower().replace("`", "")).replace(" ", "-")


def check_links(paths):
    problems, count = [], 0
    for path in paths:
        text = re.sub(r"```.*?```", "", path.read_text(), flags=re.S)
        for target in re.findall(r"\]\(([^)\s]+)\)", text):
            if re.match(r"https?:|mailto:", target):
                continue
            count += 1
            file_part, _, anchor = target.partition("#")
            dest = (path.parent / file_part).resolve() if file_part else path
            if not dest.exists() and dest != HERE / "FACT_CHECK.md":  # This script writes that file.
                problems.append(f"{path.relative_to(ROOT)}: missing {target}")
            elif anchor:
                slugs = {slug(h) for h in re.findall(r"^#+ (.+)$", dest.read_text(), re.M)}
                if anchor not in slugs:
                    problems.append(f"{path.relative_to(ROOT)}: no heading for #{anchor}")
    reproduce = re.search(r"## Reproduce\n\n```sh\n(.*?)```", README.read_text(), re.S).group(1)
    for tool in re.findall(r"\btools/[\w.\-]+", reproduce):
        count += 1
        if not (ROOT / tool).exists():
            problems.append(f"README.md Reproduce: missing {tool}")
    return count, problems


def wording_rules(docs):
    """Negated forms only for publication, users, review and validation claims; no third party; American spelling."""
    problems = []
    # Specific names live in an untracked local file, so this public script never contains them.
    blocked = HERE / ".local-terms"
    extra = [re.escape(line.strip()) for line in blocked.read_text().splitlines() if line.strip()] if blocked.exists() else []
    pattern = "|".join(["third party", r"third.party"] + extra)
    for name, text in docs.items():
        if re.search(pattern, text, re.I):
            problems.append(f"{name}: names a third party")
        if re.search(r"\b(?:favour|colour|modelling|analyse[ds]?|behaviour)\b", text, re.I):
            problems.append(f"{name}: non-American spelling")
        if "—" in text:
            problems.append(f"{name}: em dash")
        for sentence in re.split(r"(?<=[.!?])\s+|\n", text):
            if re.search(r"\b(?:published|users|peer[- ]review|peer[- ]reviewed|clinical validation|clinically validated|validated)\b", sentence, re.I) and \
               not re.search(r"\b(?:not|no|neither|never|without)\b", sentence, re.I):
                problems.append(f"{name}: unnegated claim: {sentence.strip()[:110]}")
        if re.search(r"\b(?:novel|groundbreaking|breakthrough|state-of-the-art|cutting[- ]edge|revolutionary)\b", text, re.I):
            problems.append(f"{name}: hype or novelty word")
    return problems


def main():
    run_example()
    full = README.read_text()
    assert "<!-- front-page-end" in full
    front, lower = full.split("<!-- front-page-end", 1)
    pdf = HERE / "catmech-case-study.pdf"
    pdf_text = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True, check=True).stdout
    pdf_text = re.sub(r"\s+", " ", pdf_text)
    docs = {"R": clean(front), "C": clean(CASE.read_text()), "B": clean(BLURBS.read_text()), "L": clean(lower)}
    names = {"R": "README.md front page", "C": "docs/case-study/CASE_STUDY.md and PDF", "B": "docs/case-study/BLURBS.md", "L": "README.md lower sections, numbers added or moved in this rewrite"}

    results, failures = [], []
    for shown, where, meaning, source, fn in claims:
        try:
            got = fn()
        except Exception as error:  # A broken lookup is a failed claim, never a skipped one.
            got = f"ERROR: {error!r}"
        ok = got == shown
        results.append((shown, where, meaning, source, got, ok))
        if not ok:
            failures.append(f"{shown!r} ({meaning}): source supports {got!r}")
        for d in where:
            if not re.search(BOUND_L + re.escape(shown) + BOUND_R, docs[d], re.I):
                failures.append(f"{shown!r} is claimed for {names[d]} but does not appear there")

    exempt = {d: [p for w, p, _ in not_findings if w == d] for d in "RCB"}
    for d in "RCB":
        phrases = [c[0] for c in claims if d in c[1]] + exempt[d]
        for hit in uncovered(docs[d], phrases):
            failures.append(f"Uncovered in {names[d]}: {hit}")
    for hit in uncovered(pdf_text, [c[0] for c in claims if "C" in c[1]] + exempt["C"]):
        failures.append(f"Uncovered in PDF: {hit}")

    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
    pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    if pages != 1:
        failures.append(f"PDF has {pages} pages")
    expected = build.build_html(CASE.read_text())
    if (HERE / "catmech-case-study.html").read_text() != expected:
        failures.append("HTML is stale against CASE_STUDY.md: rerun build_case_study.py")
    if pdf.stat().st_mtime < (HERE / "catmech-case-study.html").stat().st_mtime:
        failures.append("PDF is older than the HTML: rerun build_case_study.py")
    if EXAMPLE.get("warnings"):
        failures.append("The README example raised a warning")

    link_count, link_problems = check_links([README, CASE, BLURBS] + ([HERE / "FACT_CHECK.md"] if (HERE / "FACT_CHECK.md").exists() else []))
    failures += link_problems
    failures += wording_rules({"README front page": front, "CASE_STUDY.md": CASE.read_text(), "BLURBS.md": BLURBS.read_text(), "PDF": pdf_text})
    z = model_vs_linear()
    statement_results = [(text, source, bool(re.search(rx, (ROOT / source).read_text()))) for text, source, rx in statements]
    failures += [f"Statement not supported: {text} ({source})" for text, source, ok in statement_results if not ok]
    remotes = subprocess.run(["git", "remote"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    if remotes:
        failures.append(f"A git remote now exists ({remotes}), so recheck every 'not published' statement")
    misses = self_test()
    failures += [f"Negative control missed: {m}" for m in misses]

    lines = ["# Fact check", "",
             f"Written by [`fact_check.py`](fact_check.py). Last run: {datetime.date.today().isoformat()}. Each row recomputes the number from its source, and most sources are frozen results tables that the script only reads. "
             "A document fails if it contains a digit or number word that no row covers, or if a row's text no longer appears in it. Rerun with `python3 docs/case-study/fact_check.py`.", "",
             f"**Result: {'all ' + str(len(results)) + ' claims match their sources' if not failures else str(len(failures)) + ' problems'}.** "
             f"PDF: {pages} page, generated from the current `CASE_STUDY.md`. Relative links and README tool paths checked: {link_count}, {len(link_problems)} broken. External DOI links were not checked, because this check stays offline. "
             f"Negative controls: a changed number, an added number, an added number word, a wrong claim value, an unnegated publication claim, a third-party reference and broken links {'were all caught' if not misses else 'WERE NOT ALL CAUGHT'}.", ""]
    for d in "RCBL":
        lines += [f"## {names[d]}", "", "| As shown | What it says | Source file | Result |", "|---|---|---|---|"]
        for shown, where, meaning, source, got, ok in results:
            if d in where:
                lines.append(f"| `{shown}` | {meaning} | `{source}` | {'Match' if ok else 'MISMATCH: ' + got} |")
        lines.append("")
    lines += ["## Statements checked by grep", "", "| Statement | Source file | Result |", "|---|---|---|"]
    for text, source, ok in statement_results:
        lines.append(f"| {text} | `{source}` | {'Supported' if ok else 'NOT SUPPORTED'} |")
    lines += ["| Not published: this repository has no git remote | `git remote` | " + ("Supported" if not remotes else "NOT SUPPORTED") + " |", ""]
    lines += ["## Excluded as not findings", "", "| Document | Text | Why |", "|---|---|---|"]
    for w, p, why in not_findings:
        lines.append(f"| {names[w]} | {p} | {why} |")
    lines += ["", "## The README example", "",
              f"The second R block on the front page ran as written, loading the package from this tree, with no warnings. `recorded_kg` estimate {EXAMPLE['estimate']:.4f}, interval {EXAMPLE['lower']:.4f} to {EXAMPLE['upper']:.4f}. The seeds were chosen before the first run and were not changed. "
              "The numbers are one illustrative audit, so another seed gives other values, and exact cross-platform equality is not promised.", "",
              "## Paired check behind \"never clearly better\"", "",
              "Paired MSE z, contamination model minus linear MI, from `" + RP + "`. Negative favors the model. Conditional on returned intervals. This reading of the frozen table is descriptive and was not a planned test.", "",
              "| Dataset | Audit | z |", "|---|---|---:|"]
    for (family, scenario), value in sorted(z.items()):
        lines.append(f"| {family} | {scenario} | {value:+.2f} |")
    lines += ["", "## Corrections from the fact check", ""] + [f"{i}. {c}" for i, c in enumerate(CORRECTIONS, 1)] + [""]
    if not failures:
        (HERE / "FACT_CHECK.md").write_text("\n".join(lines))
    print(f"{len(results)} claims, {sum(r[5] for r in results)} match; links {link_count}, broken {len(link_problems)}; PDF pages {pages}; problems {len(failures)}")
    for failure in failures:
        print("  -", failure)
    return 1 if failures else 0


CORRECTIONS = [
    "Round 1. README and case study said the model's \"one clear win\" was the 10% NHANES audit. The paired table shows that edge is within Monte Carlo noise (z = -1.13), so both now say it was not clear. The marginal RMSE values were right, and the reading of them was too strong.",
    "Round 1. The case study said the planner \"predicted precision within 7.8%\". That number is a median error, and the 90th percentile reached 20.6%. It now says median error.",
    "Round 1. The first draft of the README example used `simulate_catmech_lab()`, whose default errors are gross, which is the documented failure case for linear MI. It now uses the public NHANES pool, with the study's 1 kg tolerance in the audit summary. Its comment on `fit$audit` was also wrong and was fixed: that table has one row per record and imputation.",
    "Round 1. README text moved from the old version said bootstrap-donor intervals were up to 52% wider. Raw widths give 51.3% (0.7557 over 0.4994, NHANES 10%). The README now says 51% and names the comparator. `validation/REAL_DATA_SYNTHESIS.md` still says 52%, from rounded table values, because this task does not own it.",
    "Round 2. Direct grep of every distinct numeric claim in the written records (the syntheses, validation record, assessment and planning guide), a path independent of the table recomputation. Of 65, 55 are printed there. The other 10 are the date format, example parameters and output, two phrases containing 95%, the 30-of-600 phrase and the corrected 51%. Nothing further needed fixing.",
]

if __name__ == "__main__":
    sys.exit(main())
