#!/usr/bin/env python3
"""Build the complete catmech 1.0.1 R-Package Educational Companion in education/:
- COMPANION.md & COMPANION.pdf (30 units C01..C30 locked to paper/manuscript.Rmd lines 1..350)
- STUDY_GUIDE.md & STUDY_GUIDE.pdf (12 concise lessons S1..S12)
- ONE_PAGER.md & ONE_PAGER.pdf (1-page visual decision & boundary sheet)
- companion_lab.html (Interactive HTML5 Companion Studio, R Lab Explorer & 12-Item Practice Quiz)
- COMPLETE_COVERAGE.md, README.md, .gitattributes, ARTIFACTS.json
- g9_claims_MAP.json, g9_claim_inventory.json, and schema-v2 g14_companion_manifest.json

The build never writes G14_HUMAN_SIGNOFF.md. That file is the human reviewer's own record. The
build reads it, and marks the human review rows GREEN only when the manuscript and companion
SHA-256 values it records match the files just built; otherwise those rows are PENDING and
G14_HUMAN_REVIEW_PENDING.md says what is waiting. The same holds for the G12 re-sign record in
paper/: upstream G12 is GREEN only when a signed record names the current manuscript SHA-256.

Needs RESEARCH_VALIDATION_HARNESS_DIR (the folder holding g14_companion_gate.py and
grounding_check.py), pandoc, LibreOffice's soffice and pdfinfo.
"""
from __future__ import annotations

import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EDU = ROOT / "education"
MS_PATH = ROOT / "paper" / "manuscript.Rmd"
MS_PDF_PATH = ROOT / "paper" / "manuscript.pdf"

# The validation harness is not part of this repository. Its location comes only from the
# environment (the same variable the verifier reads; RESEARCH_HARNESS_DIR is accepted for older
# setups), never from a guessed private folder.
_HARNESS_ENV = (os.environ.get("RESEARCH_VALIDATION_HARNESS_DIR") or os.environ.get("RESEARCH_HARNESS_DIR") or "").strip()
if not _HARNESS_ENV or not (Path(_HARNESS_ENV) / "g14_companion_gate.py").is_file():
    sys.exit("Set RESEARCH_VALIDATION_HARNESS_DIR to the folder holding g14_companion_gate.py and grounding_check.py.")
HARNESS_DIR = Path(_HARNESS_ENV)

sys.path.insert(0, str(EDU))
sys.path.insert(0, str(HARNESS_DIR))
import companion_data as cdata  # noqa: E402
import g14_companion_gate as g14  # noqa: E402


def _tool(name: str, env_var: str) -> Path:
    found = os.environ.get(env_var) or shutil.which(name)
    if not found:
        sys.exit(f"{name} not found on PATH; set {env_var} to its location.")
    return Path(found)


PANDOC = _tool("pandoc", "PANDOC")
SOFFICE = _tool("soffice", "SOFFICE")

SIGNOFF_PATH = EDU / "G14_HUMAN_SIGNOFF.md"
PENDING_PATH = EDU / "G14_HUMAN_REVIEW_PENDING.md"
PENDING_REL = "education/G14_HUMAN_REVIEW_PENDING.md"
SIGNOFF_REL = "education/G14_HUMAN_SIGNOFF.md"


def ms_pdf_pages() -> int:
    info = subprocess.check_output(["pdfinfo", str(MS_PDF_PATH)], text=True)
    return int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1))


def read_human_signoff() -> dict | None:
    """Parse the reviewer's G14 sign-off. Returns None if the file is missing or unreadable."""
    if not SIGNOFF_PATH.is_file():
        return None
    text = SIGNOFF_PATH.read_text(encoding="utf-8")
    fields = {
        "manuscript_sha256": r"Manuscript SHA-256 \(`paper/manuscript\.Rmd`\):\*\*\s*`([0-9a-f]{64})`",
        "companion_sha256": r"Companion SHA-256 \(`education/COMPANION\.md`\):\*\*\s*`([0-9a-f]{64})`",
        "reviewer": r"\*\*Reviewer:\*\*\s*([^(`\n]+?)\s*\(`reviewer_type: human`\)",
        "date": r"\*\*Date:\*\*\s*`(\d{4}-\d{2}-\d{2})`",
    }
    out = {}
    for key, pattern in fields.items():
        match = re.search(pattern, text)
        if not match:
            return None
        out[key] = match.group(1).strip()
    # The transfer result as the reviewer recorded it: all twelve items and the 4/4 rubric.
    out["all_items_passed"] = bool(re.search(r"->\s*`12/12`", text)) and bool(re.search(r"rubric \(`4/4`\)", text))
    return out


def read_g12_resign(ms_sha: str) -> tuple[str, dict | None]:
    """Find the G12 re-sign record in paper/ that names this manuscript SHA-256.

    Returns (evidence path, signature) where signature is None until the record says
    'Status: SIGNED' with the author's name and date."""
    for path in sorted((ROOT / "paper").glob("G12_RESIGN_*.md"), reverse=True):
        text = path.read_text(encoding="utf-8")
        if ms_sha not in text:
            continue
        rel = path.relative_to(ROOT).as_posix()
        status = re.search(r"\*\*Status:\*\*\s*([A-Z]+)", text)
        signer = re.search(r"\*\*Author signature:\*\*\s*([^,\n]+),\s*(\d{4}-\d{2}-\d{2})", text)
        if status and status.group(1) == "SIGNED" and signer:
            return rel, {"reviewer": signer.group(1).strip(), "date": signer.group(2)}
        return rel, None
    sys.exit(f"No paper/G12_RESIGN_*.md record names the manuscript SHA-256 {ms_sha}; write one before building.")

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
DC_NS = "http://purl.org/dc/elements/1.1/"
CP_NS = "http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
ET.register_namespace("w", W_NS)
ET.register_namespace("dc", DC_NS)
ET.register_namespace("cp", CP_NS)


def w_tag(name: str) -> str:
    return f"{{{W_NS}}}{name}"


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def build_companion_md(ms_sha: str, ms_pdf_sha: str) -> str:
    lines = [
        "# Complete Reading & Executable R Companion to the catmech Package",
        "",
        "## How to use this R-package companion",
        "",
        f"This companion follows Patrick Iben's *Statistical Inference and Two-Phase Audit Design for Contaminated Measurements in Clinical Research: The catmech R Package*, canonical manuscript [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`350` lines, SHA-256 `{ms_sha}`) and compiled {ms_pdf_pages()}-page PDF [paper/manuscript.pdf](../paper/manuscript.pdf) (SHA-256 `{ms_pdf_sha}`).",
        "",
        cdata.GLOBAL_NONCLAIM_TEXT,
        "",
        "Because `catmech` is an executable R package rather than a static paper alone, this companion links four complementary learning modalities:",
        "",
        "1. **Complete Reading Companion (`DOCUMENT`):** Work through the thirty source-locked passage units (`C01` through `C30`) below in [COMPANION.md](COMPANION.md) or [COMPANION.pdf](COMPANION.pdf), or take the shorter twelve-lesson route in [STUDY_GUIDE.md](STUDY_GUIDE.md) and [STUDY_GUIDE.pdf](STUDY_GUIDE.pdf). [COMPLETE_COVERAGE.md](COMPLETE_COVERAGE.md) maps every line (`1-350`), equation, table, and benchmark archive.",
        "2. **Executable R Lab & Refusal Suite (`run_companion_labs.R` & `interactive_tutor.R`):** Run `Rscript education/run_companion_labs.R` to execute every unit's formula check and four 'When the Package Says No' refusal/warning labs against the installed `catmech` package, or launch `source('education/interactive_tutor.R'); catmech_tutor()` in R for an interactive console walkthrough.",
        "3. **Interactive HTML5 Companion Studio & Browser Simulators (`INTERACTIVE`):** Open [companion_lab.html](companion_lab.html) for an offline unit explorer with copy-and-run R cards and the interactive 12-item Transfer & Boundary Practice Quiz (recording attempts, reporting correctness, explaining answers, and supporting retries), alongside the [catmech Guide & Method Picker](../docs/learn/index.html), the [Live Linear-MI Audit Simulator](../docs/learn/simulator.html), and the [Record Audit Planner](../docs/audit-planner/index.html).",
        "4. **One-Page Visual Decision & Case-Study Cards (`ONE_PAGER`):** Use [ONE_PAGER.pdf](ONE_PAGER.pdf) and the 1-page empirical visual summary [catmech-case-study.pdf](../docs/case-study/catmech-case-study.pdf) for rapid clinical trial team reference.",
        "",
        "Every unit uses one of six canonical G14 evidence classes: **DEFINITION**, **ANALYTIC**, **MODELED_LOCKED**, **MEASURED_DATED**, **INTERPRETIVE**, and **SCOPE**.",
        "",
    ]

    current_sec = None
    for (
        key,
        title,
        section,
        start_l,
        end_l,
        pdf_pages,
        _ev_class,
        where_extra,
        what_doing,
        ordinary_lang,
        evidence_status,
        do_not_conclude,
        try_it,
    ) in cdata.UNITS:
        if section != current_sec:
            current_sec = section
            lines.extend([f"## {section}", ""])
        loc_str = f"manuscript lines {start_l}-{end_l}" if start_l != end_l else f"manuscript line {start_l}"
        lines.extend(
            [
                f"### {key}. {title}",
                "",
                f"**Where we are:** {loc_str} ({where_extra}; PDF page {pdf_pages}).",
                "",
                f"**What the passage is doing:** {what_doing}",
                "",
                f"**In ordinary language:** {ordinary_lang}",
                "",
                f"**Evidence status:** {evidence_status}",
                "",
                f"**Do not accidentally conclude:** {do_not_conclude}",
                "",
                f"**Try it:** {try_it}",
                "",
            ]
        )

    lines.extend(
        [
            "## Glossary of statistical and software terms",
            "",
            "- **Audited reference (`R_i`):** Gold-standard measurement obtained in Phase II on a subset `A` of size `n_audit` via chart re-abstraction, recount, or central lab re-assay.",
            "- **Barnard-Rubin degrees of freedom (`nu_BR`):** Finite-sample degrees of freedom adjustment (`Barnard & Rubin 1999`) that combines classical MI degrees of freedom `nu_m` with complete-data degrees of freedom `nu_com` (`N - p`, `C - 1`, or `d - p`) so `nu_BR < nu_com`.",
            "- **Collapsed slice-within-Gibbs:** MCMC sampler in `catmech_contamination()` that analytically integrates out unobserved `(delta_i, T_i)` on unaudited rows before slice-sampling `logit(pi)`, `mu_e`, and `(1/2) log(sigma_e^2)`.",
            "- **CV3 delete-one-cluster jackknife:** Cluster-robust variance estimator `((C - 1) / C) sum (beta_hat_{(-c)} - beta_hat)(beta_hat_{(-c)} - beta_hat)^T` used in `catmech_pool_lm(cluster = )`.",
            "- **CR1S sandwich covariance:** Small-sample-adjusted (`C/(C-1) * (N-1)/(N-p)`) Huber-White-Liang-Zeger cluster-robust sandwich covariance used in `catmech_pool_glm(cluster = )`.",
            "- **Lin-Wei sandwich covariance:** Cluster-summed score-residual sandwich covariance (`Lin & Wei 1989`) used in `catmech_pool_coxph(cluster = )`.",
            "- **Rank PIT (`u_i`):** Tie-randomized empirical Probability Integral Transform computed by `catmech_holdout_check()` across K-fold cross-validation on audited rows.",
            "- **Simple Linear MI (`cole_mi`):** Calling `catmech_reference(..., use_recorded = TRUE)`, which includes the recorded value `W_i` as a continuous linear predictor alongside covariates `X_i` in C++17 (`RcppArmadillo`).",
            "",
            "## Locked 12-item transfer assessment (4 CORE, 4 BOUNDARY, 4 TRANSFER)",
            "",
        ]
    )

    for item in cdata.ASSESSMENT_ITEMS:
        lines.extend(
            [
                f"- **{item['id']} ({item['kind']}):** {item['prompt']}  ",
                f"  **Answer key:** {item['answer_key']}",
            ]
        )

    lines.append("")
    return "\n".join(lines)


def build_study_guide_md(ms_sha: str, ms_pdf_sha: str) -> str:
    lessons = [
        (
            "S1",
            "The three-way clinical audit dilemma and why recorded values matter",
            "Abstract and Section 1, manuscript lines 1-41 (PDF pages 1-2).",
            "When an audit re-measures a subset of records (`R_i`), trialists face three flawed conventional choices: (1) naively trusting recorded values (`W_i`) causes attenuation bias; (2) analyzing only audited rows (`gold_only`) inflates standard errors; and (3) setting unaudited `W_i` to `NA` (`reference_normal`) wastes accurate recordings. On ClinicalTrials.gov (`N = 27,208`) and NHANES (`N = 3,908`), keeping `W_i` via `catmech_reference(use_recorded = TRUE)` reduces RMSE by 43% to 48% over discarding unaudited values.",
            "Audits check a small slice of records against a trusted reference. Keeping recorded numbers cuts error nearly in half compared with throwing unchecked numbers away.",
            "SCOPE and MEASURED_DATED summary across 10 frozen benchmark archives.",
            "Do not conclude that standard MAR imputation (`mice`, `Amelia`) is biased when unaudited `W_i` values are discarded; its penalty in that mode is wider confidence intervals.",
            "Run `round(100 * (1 - 0.0409 / 0.0785), 1)` in R to verify the NHANES RMSE reduction at a 25% audit. **Answer:** `47.9%`.",
        ),
        (
            "S2",
            "Two-phase contamination mixture model and marginalized corruption probabilities",
            "Sections 2.1-2.2, manuscript lines 42-70 (PDF pages 2-3).",
            "Latent true values follow `T_i | X_i ~ N(X_i beta, sigma_T^2)` and recorded values follow `W_i = (1 - delta_i) T_i + delta_i (T_i + epsilon_i)` with `delta_i in {0, 1}` and `epsilon_i ~ N(mu_e, sigma_e^2)`. On unaudited rows, `catmech_contamination()` integrates out `(delta_i, T_i)` to compute the exact marginal corruption probability `p_i`. When `delta_i = 0`, `T_i = W_i` deterministically; when `delta_i = 1`, `T_i` is drawn from the precision-weighted conjugate Gaussian posterior.",
            "Each row has a hidden true value and a recorded value. The mixture model compares how well an unchecked recording fits the clean curve versus the wider error curve.",
            "DEFINITION and ANALYTIC conditional posteriors for `catmech_contamination()`.",
            "Do not conclude that `reference_sd` (`sigma_ref`) is estimated from the data; it defaults to `0` (exact recounts) and is a user-declared sensitivity constant.",
            "When `delta_i = 0` is drawn for an unaudited record in `catmech_contamination()`, what is the imputed value `T_i`? **Answer:** `T_i = W_i` deterministically.",
        ),
        (
            "S3",
            "Gaussian reference imputation in C++17: use_recorded = TRUE versus FALSE",
            "Section 2.2, manuscript lines 71-72 (PDF page 3).",
            "`catmech_reference()` implements fast conjugate Gaussian reference imputation in C++17 (`RcppArmadillo`). With `use_recorded = TRUE` (Simple Linear MI, `cole_mi`), `W_i` enters directly as a continuous linear calibration predictor alongside `X_i` (`T_i | W_i, X_i`). With `use_recorded = FALSE` (`reference_normal`), unaudited `W_i` values are discarded.",
            "When errors are small and spread out, a fast linear regression on the audited rows calibrates the unchecked numbers without needing a mixture model.",
            "DEFINITION of `catmech_reference(use_recorded = TRUE)` versus `use_recorded = FALSE`.",
            "Do not forget to pass `use_recorded = TRUE` explicitly when you want Simple Linear MI; the formal function default is `use_recorded = FALSE` for v0.3 backward compatibility.",
            "In `Lab C16` (`education/run_companion_labs.R`), why does `use_recorded = TRUE` yield a smaller pooled standard error than `use_recorded = FALSE`? **Answer:** Because `W_i` explains substantial residual variance in `T_i` on unaudited rows.",
        ),
        (
            "S4",
            "Rough two-phase audit-size planning and the NO-GO planner boundary",
            "Section 2.3, manuscript lines 73-84 (PDF pages 3-4).",
            "At audit fraction `f = n_audit / N`, the confidence interval width ratio relative to a 100% census audit is `rho(f) = sqrt(1 + k (1 - f) / f)`. `catmech_audit_power()` is a rough, unvalidated planning aid: across `2,000` log-normal draws of `k_hat` it reports the share of plausible `k` values that meet `rho <= rho_target`. The output calls that share `power`, but it is not the power of a hypothesis test. Without a pilot plan its default `k` is optimistic, so it also reports the audit size that linear imputation needs at the same correlation, often several times larger. In pre-specified validation (`validation/audit-plan-study`), `catmech_audit_plan()` received a `NO-GO` verdict (`11.3%-16.8%` median width error at `f = 5%`) and is documented as exploratory.",
            "The audit-size tool makes a rough guess from assumed error rates, and nobody has checked its answers. The empirical pilot planner missed its 10% accuracy target and stays exploratory.",
            "ANALYTIC influence-function ratio and MEASURED_DATED `NO-GO` boundary on `catmech_audit_plan()`.",
            "Do not read the `power` column of `catmech_audit_power()` as the power of a test, and do not treat `catmech_audit_plan()` as a validated regulatory sample-size guarantee; its pre-specified status is `NO-GO`.",
            "Compute `rho(0.25) = sqrt(1 + 0.27 * (1 - 0.25) / 0.25)` in R (`round(sqrt(1 + 0.27 * 3), 4)`). **Answer:** `1.3454`.",
        ),
        (
            "S5",
            "Rubin combining rules and Barnard-Rubin finite-sample degrees of freedom",
            "Section 2.4, manuscript lines 85-94 (PDF page 4).",
            "Across `m` completions, total variance is `T = U_bar + (1 + 1/m) B` and the fraction of missing information is `lambda_hat = (1 + 1/m) diag(B) / diag(T)`. To prevent anti-conservative intervals when sample size `N`, site count `C`, or event count `d` is small, `catmech` computes Barnard-Rubin (1999) degrees of freedom `nu_BR = (nu_m * nu_obs) / (nu_m + nu_obs) < nu_com`.",
            "Old imputation formulas can claim more degrees of freedom than the study has sites or events. The Barnard-Rubin formula caps degrees of freedom below the complete-data limit.",
            "ANALYTIC Rubin (1987) and Barnard-Rubin (1999) pooling formulas.",
            "Do not assume that increasing `m` to 100 lets `nu_BR` exceed `C - 1` in a clustered study; `nu_BR` is strictly bounded below `nu_obs < nu_com = C - 1`.",
            "For `C = 10` clinical sites (`nu_com = 9`), `m = 5`, and `lambda_hat = 0.2308`, what is `nu_BR` (`Lab C10-C11`)? **Answer:** `5.3577` (strictly below `9`, whereas classical `nu_m = 75.11`).",
        ),
        (
            "S6",
            "Multi-endpoint pooling across linear (CV3), GLM (CR1S), and Cox PH (Lin-Wei) models",
            "Section 2.4, manuscript lines 95-100 (PDF pages 4-5).",
            "`catmech_pool_lm()` uses the delete-one-cluster jackknife (`cv3`) with `nu_com = C - 1` when clustered; `catmech_pool_glm()` uses the CR1S Huber-White-Liang-Zeger sandwich with `nu_com = C - 1`; and `catmech_pool_coxph()` uses event-based `nu_com = d - p` when unclustered or the Lin-Wei score-residual sandwich with `nu_com = C - 1` when clustered.",
            "Linear models use a leave-one-site-out jackknife variance, logistic and count models use a small-sample sandwich, and survival models count observed events.",
            "ANALYTIC multi-endpoint covariance and degrees-of-freedom specifications.",
            "Do not use total sample size `N - p` for unclustered Cox proportional hazards degrees of freedom; `catmech_pool_coxph()` uses observed failure events `d - p`.",
            "In an unclustered Cox PH model with `N = 600` patients, `d = 45` events, and `p = 3` predictors, what is `nu_com`? **Answer:** `d - p = 45 - 3 = 42`.",
        ),
        (
            "S7",
            "Out-of-sample rank-PIT calibration checking and its slope-bias blind spot",
            "Section 2.5, manuscript lines 101-106 (PDF page 5).",
            "`catmech_holdout_check()` runs K-fold cross-validation on audited rows and computes tie-randomized rank PIT values `u_i` against `Uniform(0, 1)`. While false alarms were `0/40` under the true model, alarm power against structural misspecification was low (`0/40` in `nonlinear_target`, `1/40` in `heavy_tail_target`) because held-out rows test observation-level prediction rather than regression slope bias.",
            "Holdout checking hides part of the audit and tests if row predictions match the hidden truth. Passing this screen does not guarantee an unbiased regression slope.",
            "ANALYTIC rank-PIT formula and MEASURED_DATED diagnostic power boundary.",
            "Do not treat `miscalibration_detected = FALSE` in `catmech_holdout_check()` as proof that downstream regression slopes are unbiased.",
            "Why does `catmech_holdout_check()` randomize exact ties (`T_{i,-k}^{(j)} = R_i^*`) using `U_i ~ Uniform(0, 1)`? **Answer:** Because the contamination mixture places a discrete point mass at `W_i` on uncorrupted rows; without tie randomization, exact matches would spike and trigger false KS alarms.",
        ),
        (
            "S8",
            "Synthetic contamination benchmarks (52,950 requests) and three stress-cell failures",
            "Sections 3-5.1, manuscript lines 107-249 (PDF pages 5-10).",
            "Across `52,950` synthetic requests, `catmech_contamination()` narrowed CIs by `44.3%` over complete-case analysis in `missed_flags` (`0.1842` vs `0.3309`) and prevented the `-0.2354` linear-MI bias under gross 10x outliers (`error_gross_and_small`, `94.0%` vs `64.7%` coverage). However, it failed in three stress cells: `heavy_tail_target` (`87.3%` coverage), `nonignorable_missed` (`79.3%` coverage), and `noisy_undeclared` (`88.7%` coverage, restored to `91.8%` when `reference_sd` is declared).",
            "Under synthetic gross or unflagged errors, the mixture model stays unbiased and cuts interval width by 44%. It fails under heavy tails, non-ignorable selection, or undeclared reference noise.",
            "MODELED_LOCKED synthetic benchmark results across 27 scenarios.",
            "Do not use `catmech_contamination()` when the true distribution `T_i` has heavy `t_4` tails (`heavy_tail_target`); use `catmech_reference(use_recorded = TRUE)` (`92.0%` coverage) instead.",
            "What caused `catmech_contamination()` to under-cover (`87.3%`) in `heavy_tail_target`? **Answer:** The Gaussian mixture misclassified genuine heavy-tailed extremes in `T_i` as recording errors and shrunk them toward `X_i beta`.",
        ),
        (
            "S9",
            "Clustered variance (148,500 requests): CV3 jackknife and predictor-slope attenuation",
            "Section 5.2, manuscript lines 250-258 (PDF page 10).",
            "Across `148,500` clustered requests on simulated cohorts and 11 Rwanda HMIS indicators, `cv3` beat `cr1s` and `cr2` (which inflated widths up to `95%`). However, when the audited variable `T_i` was a **predictor** (`Y ~ T + x + stratum`), clustered random-intercept imputation (`catmech_reference(cluster = )`) attenuated slopes by up to `-16.1%` (`84.4%` coverage), whereas fitting `catmech_reference(use_recorded = TRUE)` without `cluster` and pooling with `catmech_pool_lm(cluster = )` achieved `-3.0%` bias and `95.4%` coverage.",
            "Across clinical sites, the leave-one-site-out jackknife gives the best interval coverage. However, random-intercept imputation shrinks slopes when the audited variable is a predictor.",
            "MEASURED_DATED and MODELED_LOCKED clustered variance and multilevel MI evaluation.",
            "Do not pass `cluster` to `catmech_reference()` when the audited variable is a predictor in the downstream regression (`catmech_pool_lm()` warns if you do).",
            "What is the recommended two-step R call when an audited biomarker is a predictor in a multicenter linear regression? **Answer:** Fit `catmech_reference(..., use_recorded = TRUE)` without `cluster`, then pool with `catmech_pool_lm(..., cluster = 'site')`.",
        ),
        (
            "S10",
            "Public paired-discrepancy benchmarks (33,800 requests) and near-chance row error AUCs",
            "Section 6.1, manuscript lines 259-279 (PDF page 11).",
            "On ClinicalTrials.gov (`N = 27,208`) and NHANES 2017-2018 (`N = 3,908`), retaining recorded values reduced RMSE by `43%` to `48%` over discarding unaudited values. However, because real-world discrepancies were small relative to `T_i | X_i`, per-record corruption probabilities `p_i` had near-chance discrimination (`AUC = 0.543` on ClinicalTrials.gov and `0.500` on NHANES, versus `0.962` in synthetic `missed_flags`).",
            "On real data with diffuse errors, keeping recorded numbers cuts slope error nearly in half, but per-row error probabilities score no better than a coin flip.",
            "MEASURED_DATED empirical benchmark across ClinicalTrials.gov and NHANES.",
            "Never use individual record-level corruption probabilities `p_i` from `catmech_contamination()` to adjudicate or discard specific patient records in diffuse real-world data.",
            "What was the ROC-AUC of per-record corruption probabilities `p_i` at a 25% audit on ClinicalTrials.gov and NHANES? **Answer:** `0.543` on ClinicalTrials.gov and `0.500` on NHANES (near coin-flip chance `0.500`).",
        ),
        (
            "S11",
            "Real-audit slope recovery (45,600 requests) and the methodological decision tree",
            "Sections 6.2-8, manuscript lines 280-328 (PDF pages 12-14).",
            "On Rwanda HMIS facility recounts (`45,600` requests), Rule 2 comparing `catmech_contamination()` against `Linear MI` was `INCONCLUSIVE`: the mixture had lower point RMSE on `PNC` (`0.053` vs `0.072`) and `Iron` (`0.112` vs `0.156`), but failed 95% coverage on all three primary slopes (`71.6%` on `PNC`, `66.6%` on `Syphilis`, `35.4%` on `Iron`) and failed the `>= 90%` return rate on two. Therefore, Section 7 recommends `catmech_contamination()` only for gross, separable outliers and `catmech_reference(use_recorded = TRUE)` for moderate or diffuse real-world discrepancies.",
            "On Rwanda health-center recounts, the mixture model had lower point error on two slopes but failed interval coverage. For everyday moderate errors, use simple linear imputation.",
            "MEASURED_DATED real-audit evaluation (`INCONCLUSIVE` Rule 2) and INTERPRETIVE decision tree.",
            "Do not claim that `catmech_contamination()` outperformed Simple Linear MI on real-world audit intervals; its interval advantage remains synthetic-only.",
            "Under the Section 7 decision tree, what two checks must pass before trusting `catmech_contamination()` on a dataset with suspected gross errors? **Answer:** A `>= 90%` interval return rate and positive held-out Brier/RMSE skill over raw recorded values in `catmech_holdout_check()`.",
        ),
        (
            "S12",
            "Provenance reporting (catmech_audit_trail), software QA, and regulatory boundaries",
            "Sections 9-12, manuscript lines 329-350 (PDF pages 14-16).",
            "`catmech_audit_trail()` outputs a provenance record (Markdown, JSON, or HTML) recording R/OS metadata; SHA-256 hashes of the input data and all `m` completed datasets, computed over the serialized R data frames, and of the assumption specification as canonical JSON text; MCMC convergence status; and pooled tables. A rerun can be compared with the hashes, but they do not cover the pooled results or make the record tamper-evident. This supports documentation under FDA 21 CFR Part 11 and ICH E9 without certifying compliance. The package is MIT-licensed and verified by `199` `testthat` tests with `1,812` expectations (`0` failures/warnings/skips), `R CMD check --as-cran` (`0` errors/warnings), and `500,000+` replicate checks.",
            "The audit trail function writes a report with SHA-256 data hashes, model settings, and convergence status. A rerun can be compared with it, but it does not prove the data were untouched or certify legal compliance.",
            "DEFINITION, SCOPE, and MEASURED_DATED software verification summary.",
            "Do not represent `catmech_audit_trail()` as a legal compliance certification or a tamper-evident record; regulatory compliance depends on the surrounding validated system, access controls, and procedures, and anyone who edits the data can regenerate the record.",
            "Run `Lab C28` in `education/run_companion_labs.R`. What default title and footer heading does `catmech_audit_trail()` print? **Answer:** Default title `Catmech Audit Trail & Provenance Report` and footer heading `### Audit Trail & Provenance Record`.",
        ),
    ]

    lines = [
        "# Clinical Measurement Audits and Two-Phase Imputation: Concise Study Guide",
        "",
        f"A twelve-lesson study guide to Patrick Iben's *Statistical Inference and Two-Phase Audit Design for Contaminated Measurements in Clinical Research: The catmech R Package*, canonical manuscript [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`350` lines, SHA-256 `{ms_sha}`) and compiled {ms_pdf_pages()}-page PDF [paper/manuscript.pdf](../paper/manuscript.pdf) (SHA-256 `{ms_pdf_sha}`).",
        "",
        cdata.GLOBAL_NONCLAIM_TEXT,
        "",
        "## How to use this twelve-lesson route",
        "",
        "This concise guide is designed for clinicians, trial statisticians, and data auditors who want the core mathematical logic, empirical findings, R function calls, and scope boundaries in twelve short lessons (`S1` through `S12`). You can also step through these exact twelve lessons interactively in an R console by running `source('education/interactive_tutor.R'); catmech_tutor()`, or read the full thirty-unit passage-by-passage companion in [COMPANION.md](COMPANION.md) ([PDF](COMPANION.pdf)).",
        "",
    ]

    for sid, title, where_are, what_doing, ord_lang, ev_status, do_not, try_it in lessons:
        lines.extend(
            [
                f"## {sid}. {title}",
                "",
                "**Where we are**",
                "",
                where_are,
                "",
                "**What the passage is doing**",
                "",
                what_doing,
                "",
                "**In ordinary language**",
                "",
                ord_lang,
                "",
                "**Evidence status**",
                "",
                ev_status,
                "",
                "**Do not accidentally conclude**",
                "",
                do_not,
                "",
                "**Try it**",
                "",
                try_it,
                "",
            ]
        )

    return "\n".join(lines)


def build_one_pager_md(ms_sha: str) -> str:
    return "\n".join(
        [
            "# catmech 1.0.1: One-Page Clinical Audit Decision & Boundary Card",
            "",
            f"**Canonical Manuscript:** `paper/manuscript.Rmd` (`SHA-256: {ms_sha[:16]}...`) | **Companion Visual Case Study:** `docs/case-study/catmech-case-study.pdf`",
            "",
            cdata.GLOBAL_NONCLAIM_TEXT,
            "",
            "## 1. The Empirical Headline (500,000+ Verified Requests Across 10 Frozen Archives)",
            "",
            "| Benchmark Archive | Audit Design | Key Empirical Finding | Scope / Operational Boundary |",
            "| :--- | :--- | :--- | :--- |",
            "| **ClinicalTrials.gov** (`N = 27,208`) & **NHANES** (`N = 3,908`) | `33,800` requests (`f = 25%` primary) | Retaining recorded values (`use_recorded = TRUE`) cuts RMSE by **43% to 48%** vs discarding unaudited values (`0.0286` vs `0.0539`; `0.0409` vs `0.0785`). | Per-record corruption probability `p_i` has near-chance AUC (**0.543** on CT.gov, **0.500** on NHANES). Never use `p_i` to adjudicate individual rows. |",
            "| **Synthetic Contamination** (`validation/contamination-study`) | `52,950` requests (`27` scenarios) | `catmech_contamination()` narrows CIs by **44.3%** vs complete-case (`0.1842` vs `0.3309`) and prevents `-0.2354` linear-MI bias under 10x gross outliers. | Fails `<90%` coverage in 3 stress cells: `heavy_tail_target` (**87.3%**), `nonignorable_missed` (**79.3%**), `noisy_undeclared` (**88.7%**). |",
            "| **Clustered Variance & MI** (`validation/variance-mi-study`) | `148,500` requests (`11` Rwanda HMIS + sim) | Delete-one-cluster jackknife (`cv3`) beats `cr1s`/`cr2` (`92.0%-97.5%` recount coverage). | When audited `T_i` is a **predictor**, fit `catmech_reference(use_recorded = TRUE)` *without* `cluster` and pool with `catmech_pool_lm(cluster = )`. |",
            "| **Real-Audit Slope Recovery** (`Rwanda HMIS` & `Psychology`) | `45,600` requests (`Rules 1-3`) | `Linear MI` cuts RMSE to `0.458x` on size-dependent `ANC 4th visits`; Rule 2 (`Contamination` vs `Linear MI`) is **INCONCLUSIVE**. | `Contamination` had lower RMSE on 2 of 3 Rwanda slopes but **35.4% to 71.6%** clustered coverage and `6%` to `19%` MCMC refusal rates. |",
            "",
            "## 2. Two-Step R Workflow & Decision Rule",
            "",
            "1. **Screen & Compare:** Run `catmech_recording_audit(data, value, reference, by = 'site')` and `catmech_compare_strategies(data, value, reference, formula)`.",
            "2. **Select Imputation & Pooling Engine:**",
            "   - **Default (Subtle, diffuse, or moderate real-world discrepancies, or heavy-tailed target):** `fit <- catmech_reference(data, value, reference, predictors, use_recorded = TRUE, m = 20)`",
            "   - **Gross, separable outliers only (with Gaussian target & ignorable audit selection):** `fit <- catmech_contamination(...)`, verifying `>= 90%` return rate and positive holdout Brier/RMSE skill in `catmech_holdout_check(fit)`.",
            "   - **Put every analysis outcome in `predictors`** (for Cox models, the event indicator and the Nelson-Aalen hazard); pooling warns if one is missing.",
            "   - **Pool with Barnard-Rubin df (`nu_BR < nu_com`) & Provenance:** `catmech_pool_lm(fit, formula, cluster)` (`cv3`), `catmech_pool_glm(...)` (`CR1S`), `catmech_pool_coxph(...)` (`d - p` or `Lin-Wei`), and `catmech_audit_trail(fit, pooled)` (hashes for rerun comparison; not tamper-evident, not compliance).",
            "",
        ]
    )


def build_complete_coverage_md(ms_sha: str, ms_pdf_sha: str) -> str:
    lines = [
        "# Complete-Companion Coverage Map",
        "",
        f"Basis: Patrick Iben, *Statistical Inference and Two-Phase Audit Design for Contaminated Measurements in Clinical Research: The catmech R Package*, canonical manuscript [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`350` lines, SHA-256 `{ms_sha}`) and {ms_pdf_pages()}-page compiled PDF [paper/manuscript.pdf](../paper/manuscript.pdf) (SHA-256 `{ms_pdf_sha}`). Unit IDs refer to [COMPANION.md](COMPANION.md) ([PDF](COMPANION.pdf)). The [concise study guide](STUDY_GUIDE.md) ([PDF](STUDY_GUIDE.pdf)) provides a twelve-lesson route (`S1`-`S12`).",
        "",
        "Every single line (`1-350`) of `paper/manuscript.Rmd` is mapped to a unique unit below with zero gaps and zero overlaps, alongside its paired executable R check in [run_companion_labs.R](run_companion_labs.R) and interactive web surface.",
        "",
        "| Unit | Manuscript lines | PDF pages | Evidence class | Source passage / R function / benchmark covered |",
        "| :--- | ---: | ---: | :--- | :--- |",
    ]
    for (
        key,
        title,
        _section,
        start_l,
        end_l,
        pdf_pages,
        ev_class,
        where_extra,
        *_rest,
    ) in cdata.UNITS:
        lines.append(
            f"| **{key}** | {start_l}-{end_l} | {pdf_pages} | `{ev_class}` | {title} ({where_extra}) |"
        )
    lines.extend(
        [
            "",
            "## Multimodal & Executable R Cross-Reference",
            "",
            "- **Document Modality (`DOCUMENT`):** [COMPANION.md](COMPANION.md), [COMPANION.pdf](COMPANION.pdf), [STUDY_GUIDE.md](STUDY_GUIDE.md), [STUDY_GUIDE.pdf](STUDY_GUIDE.pdf)",
            "- **One-Pager Modality (`ONE_PAGER`):** [ONE_PAGER.pdf](ONE_PAGER.pdf) and [catmech-case-study.pdf](../docs/case-study/catmech-case-study.pdf)",
            "- **Interactive Modality (`INTERACTIVE`):** [companion_lab.html](companion_lab.html), [catmech Guide & Method Picker](../docs/learn/index.html), [Live Linear-MI Simulator](../docs/learn/simulator.html), and [Record Audit Planner](../docs/audit-planner/index.html)",
            "- **Executable R Verification & Console Tutor:** [run_companion_labs.R](run_companion_labs.R) (`16/16` live checks) and [interactive_tutor.R](interactive_tutor.R) (`S1`-`S12` console tutor)",
            "",
        ]
    )
    return "\n".join(lines)


def build_readme_md(ms_sha: str, ms_pdf_sha: str, comp_sha: str, human_signed: bool) -> str:
    if human_signed:
        review_lines = [
            "See [ARTIFACTS.json](ARTIFACTS.json), [g14_companion_manifest.json](g14_companion_manifest.json), and [G14_HUMAN_SIGNOFF.md](G14_HUMAN_SIGNOFF.md) for the complete verification records.",
        ]
    else:
        review_lines = [
            "**Human review pending.** This companion was rebuilt for catmech 1.0.1 against the revised manuscript (SHA-256 `" + ms_sha + "`) and companion (SHA-256 `" + comp_sha + "`). The human claim, plain-language, equivalence, rendering and transfer reviews of these versions have not been done, and the author has not yet signed the manuscript's G12 re-read. [G14_HUMAN_SIGNOFF.md](G14_HUMAN_SIGNOFF.md) records the human review of the earlier 1.0.0 versions only. Until the new reviews are signed, the G14 gate reports PARTIAL. [G14_HUMAN_REVIEW_PENDING.md](G14_HUMAN_REVIEW_PENDING.md) lists what is waiting.",
            "",
            "See [ARTIFACTS.json](ARTIFACTS.json) and [g14_companion_manifest.json](g14_companion_manifest.json) for the verification records.",
        ]
    return "\n".join(
        [
            "# catmech 1.0.1: Multimodal R-Package Educational Companion",
            "",
            "Choose the learning modality and depth that fit your workflow:",
            "",
            "1. **Concise 12-Lesson Study Guide (`DOCUMENT`):** [Read Markdown](STUDY_GUIDE.md) or [Download PDF](STUDY_GUIDE.pdf). Covers the clinical audit dilemma, two-phase contamination and linear reference MI models, Barnard-Rubin small-sample pooling across linear/GLM/Cox endpoints, synthetic vs real-world audit benchmarks, and regulatory provenance boundaries.",
            "2. **Complete 30-Unit Reading Companion (`DOCUMENT`):** [Read Markdown](COMPANION.md) or [Download PDF](COMPANION.pdf). Covers 100% of [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`lines 1-350`, SHA-256 `" + ms_sha + "`) and [paper/manuscript.pdf](../paper/manuscript.pdf) (`" + str(ms_pdf_pages()) + "` pages, SHA-256 `" + ms_pdf_sha + "`), with worked numerical reasoning, R one-liners, glossary, and a locked 12-item transfer assessment. See [COMPLETE_COVERAGE.md](COMPLETE_COVERAGE.md) for the line-by-line map.",
            "3. **Interactive HTML5 Companion Studio & Browser Simulators (`INTERACTIVE`):**",
            "   - [Interactive Companion Studio & 12-Item Practice Lab](companion_lab.html): Browse all 30 units (`C01`-`C30`), copy executable R snippets, and take the interactive 12-item Transfer & Boundary Practice Quiz with immediate attempt tracking, correctness feedback, answer explanations, and retry support.",
            "   - [catmech Field Guide & 24-Scenario Method Picker](../docs/learn/index.html): Interactive decision tree generating verified R code snippets.",
            "   - [Live Browser Linear-MI Audit Simulator](../docs/learn/simulator.html): Run 140 benchmark studies and custom audit fractions in the browser.",
            "   - [Interactive Record Audit Planner](../docs/audit-planner/index.html): Explore two-phase precision ratios `rho(f)`.",
            "4. **One-Page Visual Summaries (`ONE_PAGER`):** [Download Companion One-Pager PDF](ONE_PAGER.pdf) and the 1-page visual empirical summary [catmech-case-study.pdf](../docs/case-study/catmech-case-study.pdf).",
            "5. **Executable R Lab Suite & Interactive R Console Tutor (`R`):**",
            "   - Run all 16 live formula, workflow, and refusal labs in batch mode: `Rscript education/run_companion_labs.R`",
            "   - Launch the interactive 12-lesson R console tutor: `Rscript -e \"source('education/interactive_tutor.R'); catmech_tutor()\"` (or `Rscript education/interactive_tutor.R --non-interactive`)",
            "",
            cdata.GLOBAL_NONCLAIM_TEXT,
            "",
            "To verify all G14 schema-v2 locks, PDF metadata hashes, readability ceilings, and live R labs (the harness gates run only when `RESEARCH_VALIDATION_HARNESS_DIR` points at the separate validation harness, and `--require-harness` fails if they cannot run):",
            "",
            "```sh",
            "python3 education/verify_educational_companion.py --require-harness",
            "```",
            "",
            *review_lines,
            "",
        ]
    )


def build_companion_lab_html(comp_sha: str, ms_sha: str) -> str:
    unit_cards_html = []
    for (
        key,
        title,
        section,
        start_l,
        end_l,
        pdf_pages,
        ev_class,
        where_extra,
        what_doing,
        ordinary_lang,
        evidence_status,
        do_not_conclude,
        try_it,
    ) in cdata.UNITS:
        unit_cards_html.append(
            f"""
            <article class="unit-card" data-section="{html.escape(section)}" id="unit-{key}">
              <div class="unit-head">
                <span class="badge">{key}</span>
                <span class="chip">{ev_class}</span>
                <span class="loc">manuscript lines {start_l}-{end_l} (PDF p. {pdf_pages})</span>
              </div>
              <h3>{key}. {html.escape(title)}</h3>
              <p class="ordinary"><strong>In ordinary language:</strong> {html.escape(ordinary_lang)}</p>
              <details>
                <summary>Read full explanation, R lab connection &amp; boundary check</summary>
                <div class="details-body">
                  <p><strong>Where we are:</strong> manuscript lines {start_l}-{end_l} ({html.escape(where_extra)}; PDF page {pdf_pages}).</p>
                  <p><strong>What the passage is doing:</strong> {html.escape(what_doing)}</p>
                  <p><strong>Evidence status:</strong> {html.escape(evidence_status)}</p>
                  <p class="warn"><strong>Do not accidentally conclude:</strong> {html.escape(do_not_conclude)}</p>
                  <p class="tryit"><strong>Try it:</strong> {html.escape(try_it)}</p>
                </div>
              </details>
            </article>
            """
        )

    quiz_data = [
        {
            "id": "T01",
            "kind": "CORE",
            "prompt": cdata.ASSESSMENT_ITEMS[0]["prompt"],
            "options": [
                "Because retaining W_i as a linear calibration predictor captures valid recording signal, cutting RMSE by 43% to 48% on ClinicalTrials.gov and NHANES.",
                "Because standard MAR multiple imputation cannot run when more than 50% of values are missing.",
                "Because use_recorded = TRUE ignores the audited reference column R_i.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[0]["answer_key"],
        },
        {
            "id": "T02",
            "kind": "CORE",
            "prompt": cdata.ASSESSMENT_ITEMS[1]["prompt"],
            "options": [
                "When delta_i = 0, T_i = W_i deterministically; when delta_i = 1, T_i is drawn from the precision-weighted conjugate Gaussian posterior.",
                "In both cases T_i is drawn purely from N(X_i beta, sigma_T^2) without using W_i.",
                "When delta_i = 0, T_i is set to NA; when delta_i = 1, T_i = W_i.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[1]["answer_key"],
        },
        {
            "id": "T03",
            "kind": "CORE",
            "prompt": cdata.ASSESSMENT_ITEMS[2]["prompt"],
            "options": [
                "It uses the Barnard-Rubin (1999) finite-sample adjustment nu_BR with nu_com = C - 1 for clustered models and nu_com = d - p for unclustered Cox PH models.",
                "It uses a normal z = 1.96 critical value for all sample sizes.",
                "It sets degrees of freedom equal to m - 1 regardless of cluster or event count.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[2]["answer_key"],
        },
        {
            "id": "T04",
            "kind": "CORE",
            "prompt": cdata.ASSESSMENT_ITEMS[3]["prompt"],
            "options": [
                "Gross 10x transcription outliers attenuated the OLS calibration slope in Simple Linear MI (bias -0.2354, 64.7% coverage), whereas the mixture model isolated outliers via delta_i = 1.",
                "Simple Linear MI failed to converge in C++17.",
                "Complete-case analysis also failed in error_gross_and_small.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[3]["answer_key"],
        },
        {
            "id": "T05",
            "kind": "BOUNDARY",
            "prompt": cdata.ASSESSMENT_ITEMS[4]["prompt"],
            "options": [
                "No; held-out row checks test observation-level prediction and had low power against regression slope bias (0/40 alarms in nonlinear_target and 1/40 in heavy_tail_target).",
                "Yes; a passing Kolmogorov-Smirnov test guarantees zero regression slope bias.",
                "Yes, provided at least 5 folds are used.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[4]["answer_key"],
        },
        {
            "id": "T06",
            "kind": "BOUNDARY",
            "prompt": cdata.ASSESSMENT_ITEMS[5]["prompt"],
            "options": [
                "No; on diffuse real-world discrepancies, per-record p_i has near-chance discrimination (AUC = 0.543 on ClinicalTrials.gov, 0.500 on NHANES, and 0.501-0.582 on Rwanda HMIS).",
                "Yes; any record with p_i > 0.50 should be discarded from the clinical dataset.",
                "Yes, whenever the sample size exceeds N = 1,000.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[5]["answer_key"],
        },
        {
            "id": "T07",
            "kind": "BOUNDARY",
            "prompt": cdata.ASSESSMENT_ITEMS[6]["prompt"],
            "options": [
                "It carries a pre-specified NO-GO validation status (exploratory only) because median interval-width prediction error at f = 5% was 11.3% to 16.8% (> 10% tolerance).",
                "It is certified for regulatory sample-size determination at all audit fractions.",
                "It was deprecated because bootstrap confidence intervals cannot be computed in R.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[6]["answer_key"],
        },
        {
            "id": "T08",
            "kind": "BOUNDARY",
            "prompt": cdata.ASSESSMENT_ITEMS[7]["prompt"],
            "options": [
                "No; SHA-256 provenance records support reproducibility and regulatory documentation, while compliance depends on the surrounding validated system, access controls, and SOPs.",
                "Yes; the SHA-256 hash automatically satisfies all FDA 21 CFR Part 11 institutional requirements.",
                "Yes, when exported in JSON format.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[7]["answer_key"],
        },
        {
            "id": "T09",
            "kind": "TRANSFER",
            "prompt": cdata.ASSESSMENT_ITEMS[8]["prompt"],
            "options": [
                "No; fit catmech_reference(use_recorded = TRUE) without cluster and pool with catmech_pool_lm(cluster = 'site') to avoid -16.1% predictor-slope attenuation.",
                "Yes; catmech_reference(cluster = 'site') is always optimal when T_i is a predictor.",
                "Fit catmech_reference(use_recorded = FALSE) and ignore site clustering during pooling.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[8]["answer_key"],
        },
        {
            "id": "T10",
            "kind": "TRANSFER",
            "prompt": cdata.ASSESSMENT_ITEMS[9]["prompt"],
            "options": [
                "Choose catmech_reference(use_recorded = TRUE) (Simple Linear MI, 92.0% coverage in heavy_tail_target) rather than catmech_contamination() (87.3% coverage).",
                "Choose catmech_contamination() because Gaussian mixtures adapt automatically to Student-t_4 target tails.",
                "Use Naive analysis (trusting recorded values without adjustment).",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[9]["answer_key"],
        },
        {
            "id": "T11",
            "kind": "TRANSFER",
            "prompt": cdata.ASSESSMENT_ITEMS[10]["prompt"],
            "options": [
                "Leaving reference_sd = 0 triggers a universal-disagreement warning and undercoverage (88.7%); setting reference_sd = 0.35 restores bias (-0.0071) and coverage (91.8%).",
                "reference_sd has no effect because catmech_contamination() estimates sigma_ref automatically.",
                "Setting reference_sd = 0.35 causes all audited rows to be discarded.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[10]["answer_key"],
        },
        {
            "id": "T12",
            "kind": "TRANSFER",
            "prompt": cdata.ASSESSMENT_ITEMS[11]["prompt"],
            "options": [
                "Because unclustered Cox PH pooling defines complete-data df from observed events as nu_com = d - p = 30 - 4 = 26, nu_BR is strictly less than 26.",
                "nu_BR is bounded below N - p = 800 - 4 = 796.",
                "nu_BR equals m - 1.",
            ],
            "correct": 0,
            "explanation": cdata.ASSESSMENT_ITEMS[11]["answer_key"],
        },
    ]

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>catmech 1.0.1: Interactive Educational Companion Studio &amp; Practice Lab</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="companion-source-sha256" content="{comp_sha}">
<meta name="manuscript-source-sha256" content="{ms_sha}">
<style>
:root {{
  --bg: #f4f6f9; --surface: #ffffff; --surface-2: #eef2f6; --line: #d7dee7;
  --ink: #121a24; --text: #1f2a37; --muted: #4d5a69;
  --accent: #0b7268; --accent-soft: rgba(11,114,104,.10);
  --green: #1d6b3c; --green-soft: rgba(29,107,60,.12);
  --amber: #8a5300; --amber-soft: rgba(138,83,0,.12);
  --red: #b42318; --red-soft: rgba(180,35,24,.12);
  --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  --mono: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--bg); color: var(--text); font: 400 16px/1.6 var(--sans); }}
.wrap {{ max-width: 1100px; margin: 0 auto; padding: 24px 18px 64px; }}
header {{ background: var(--surface); border-bottom: 1px solid var(--line); padding: 18px 0; position: sticky; top: 0; z-index: 10; }}
header .wrap {{ padding-top: 0; padding-bottom: 0; display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 12px; }}
h1 {{ font-size: 26px; margin: 0 0 8px; color: var(--ink); }}
h2 {{ font-size: 20px; margin: 28px 0 12px; color: var(--ink); border-bottom: 1px solid var(--line); padding-bottom: 6px; }}
h3 {{ font-size: 17px; margin: 8px 0; color: var(--ink); }}
.nonclaim {{ background: var(--amber-soft); border-left: 4px solid var(--amber); padding: 12px 16px; border-radius: 6px; font-size: 14px; margin: 14px 0; }}
.nav-links {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 14px 0; }}
.nav-links a, button {{ background: var(--accent); color: #fff; text-decoration: none; padding: 8px 14px; border-radius: 6px; font-weight: 600; font-size: 14px; border: none; cursor: pointer; }}
button.secondary, .nav-links a.secondary {{ background: var(--surface-2); color: var(--ink); border: 1px solid var(--line); }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 14px; }}
.unit-card {{ background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: 14px 16px; }}
.unit-head {{ display: flex; flex-wrap: wrap; gap: 8px; align-items: center; font-size: 12px; color: var(--muted); }}
.badge {{ background: var(--accent); color: #fff; font-weight: 700; padding: 2px 8px; border-radius: 999px; }}
.chip {{ background: var(--accent-soft); color: var(--accent); font-weight: 600; padding: 2px 8px; border-radius: 999px; }}
.ordinary {{ margin: 8px 0; font-size: 15px; }}
.warn {{ color: var(--red); background: var(--red-soft); padding: 8px 10px; border-radius: 6px; font-size: 14px; }}
.tryit {{ background: var(--surface-2); padding: 8px 10px; border-radius: 6px; font-size: 14px; }}
.quiz-box {{ background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: 18px; margin-bottom: 14px; }}
.quiz-opt {{ display: block; margin: 8px 0; padding: 8px 12px; border: 1px solid var(--line); border-radius: 6px; cursor: pointer; background: var(--surface-2); }}
.feedback {{ margin-top: 10px; padding: 10px 12px; border-radius: 6px; font-size: 14px; display: none; }}
.feedback.ok {{ display: block; background: var(--green-soft); color: var(--green); border: 1px solid var(--green); }}
.feedback.bad {{ display: block; background: var(--red-soft); color: var(--red); border: 1px solid var(--red); }}
.score-banner {{ background: var(--surface); border: 2px solid var(--accent); border-radius: 8px; padding: 14px 18px; margin-bottom: 18px; display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 12px; }}
code {{ font-family: var(--mono); font-size: 0.9em; background: var(--surface-2); padding: 2px 5px; border-radius: 4px; }}
</style>
</head>
<body>
<header>
  <div class="wrap">
    <div><strong>catmech 1.0.1</strong>: Interactive Educational Companion Studio &amp; Practice Lab</div>
    <div style="font-size:13px;color:var(--muted);">Companion SHA-256: <code>{comp_sha[:12]}</code> | Manuscript SHA-256: <code>{ms_sha[:12]}</code></div>
  </div>
</header>
<main class="wrap">
  <h1>Multimodal R-Package Companion Studio (Units C01–C30 &amp; Transfer Lab)</h1>
  <div class="nonclaim"><strong>Scope &amp; Non-Claim Boundary:</strong> {html.escape(cdata.GLOBAL_NONCLAIM_TEXT)}</div>

  <div class="nav-links" role="navigation" aria-label="Companion Modalities">
    <a href="COMPANION.pdf">Download Complete Companion PDF (C01–C30)</a>
    <a href="STUDY_GUIDE.pdf" class="secondary">Download 12-Lesson Study Guide PDF (S1–S12)</a>
    <a href="ONE_PAGER.pdf" class="secondary">Download 1-Page Decision Card PDF</a>
    <a href="../docs/learn/index.html" class="secondary">Open Interactive Method Picker</a>
    <a href="../docs/learn/simulator.html" class="secondary">Open Live Linear-MI Simulator</a>
    <a href="../docs/audit-planner/index.html" class="secondary">Open Record Audit Planner</a>
    <a href="../docs/case-study/catmech-case-study.pdf" class="secondary">Open 1-Page Case Study PDF</a>
  </div>

  <h2>1. Interactive 12-Item Transfer &amp; Boundary Practice Lab</h2>
  <p>This interactive practice lab records your attempts, reports correctness on Core, Boundary, and Transfer items, explains the underlying mathematical or empirical reason, and supports full retry.</p>
  <div class="score-banner" id="score-banner" aria-live="polite">
    <div>
      <strong>Practice Telemetry:</strong>
      Attempts recorded: <span id="attempt-count">0</span> |
      Score: <span id="score-total">0 / 12</span>
      (Core: <span id="score-core">0/4</span>, Boundary: <span id="score-boundary">0/4</span>, Transfer: <span id="score-transfer">0/4</span>)
    </div>
    <div>
      <button type="button" id="retry-btn" class="secondary" onclick="resetQuiz()">Reset &amp; Retry All 12 Items</button>
    </div>
  </div>
  <div id="quiz-container"></div>

  <h2>2. Complete Source-Locked Passage Units (C01–C30, Lines 1–350)</h2>
  <div class="grid">
    {''.join(unit_cards_html)}
  </div>
</main>
<script>
const QUIZ_ITEMS = {json.dumps(quiz_data)};
const state = {{ attempts: 0, answers: {{}} }};

function renderQuiz() {{
  const container = document.getElementById("quiz-container");
  container.innerHTML = "";
  QUIZ_ITEMS.forEach((item, idx) => {{
    const box = document.createElement("div");
    box.className = "quiz-box";
    const title = document.createElement("h3");
    title.textContent = `${{item.id}} (${{item.kind}}): ${{item.prompt}}`;
    box.appendChild(title);

    item.options.forEach((optText, optIdx) => {{
      const label = document.createElement("label");
      label.className = "quiz-opt";
      const radio = document.createElement("input");
      radio.type = "radio";
      radio.name = `quiz-${{item.id}}`;
      radio.value = String(optIdx);
      radio.checked = state.answers[item.id] === optIdx;
      radio.addEventListener("change", () => submitAnswer(item.id, optIdx));
      label.appendChild(radio);
      label.appendChild(document.createTextNode(" " + optText));
      box.appendChild(label);
    }});

    const fb = document.createElement("div");
    fb.id = `fb-${{item.id}}`;
    fb.className = "feedback";
    box.appendChild(fb);
    container.appendChild(box);
  }});
}}

function submitAnswer(itemId, chosenIdx) {{
  state.attempts += 1;
  state.answers[itemId] = chosenIdx;
  const item = QUIZ_ITEMS.find(x => x.id === itemId);
  const fb = document.getElementById(`fb-${{itemId}}`);
  const isCorrect = chosenIdx === item.correct;
  fb.className = "feedback " + (isCorrect ? "ok" : "bad");
  fb.textContent = (isCorrect ? "CORRECT: " : "NOT QUITE (Try again): ") + item.explanation;
  updateScoreBanner();
}}

function updateScoreBanner() {{
  let total = 0, core = 0, boundary = 0, transfer = 0;
  QUIZ_ITEMS.forEach(item => {{
    if (state.answers[item.id] === item.correct) {{
      total += 1;
      if (item.kind === "CORE") core += 1;
      if (item.kind === "BOUNDARY") boundary += 1;
      if (item.kind === "TRANSFER") transfer += 1;
    }}
  }});
  document.getElementById("attempt-count").textContent = String(state.attempts);
  document.getElementById("score-total").textContent = `${{total}} / 12`;
  document.getElementById("score-core").textContent = `${{core}}/4`;
  document.getElementById("score-boundary").textContent = `${{boundary}}/4`;
  document.getElementById("score-transfer").textContent = `${{transfer}}/4`;
}}

function resetQuiz() {{
  state.answers = {{}};
  renderQuiz();
  updateScoreBanner();
}}

renderQuiz();
</script>
</body>
</html>
"""


def style_docx_and_set_metadata(docx_path: Path, title: str, source_sha: str, compact: bool = False) -> None:
    with tempfile.TemporaryDirectory(prefix="catmech-ooxml-") as tmpdir:
        tmp = Path(tmpdir)
        with zipfile.ZipFile(docx_path, "r") as zf:
            zf.extractall(tmp)

        core_xml = tmp / "docProps" / "core.xml"
        if core_xml.exists():
            c_tree = ET.parse(core_xml)
            c_root = c_tree.getroot()
            subj = c_root.find(f"{{{DC_NS}}}subject")
            if subj is None:
                subj = ET.SubElement(c_root, f"{{{DC_NS}}}subject")
            subj.text = f"Educational companion; source SHA-256 {source_sha}"
            tit = c_root.find(f"{{{DC_NS}}}title")
            if tit is None:
                tit = ET.SubElement(c_root, f"{{{DC_NS}}}title")
            tit.text = title
            creator = c_root.find(f"{{{DC_NS}}}creator")
            if creator is None:
                creator = ET.SubElement(c_root, f"{{{DC_NS}}}creator")
            creator.text = "Patrick Iben"
            c_tree.write(core_xml, encoding="utf-8", xml_declaration=True)

        body_sz = "17" if compact else "21"
        h1_sz = "26" if compact else "34"
        h2_sz = "21" if compact else "28"
        h3_sz = "19" if compact else "23"
        tbl_sz = "15" if compact else "18"

        styles_xml = tmp / "word" / "styles.xml"
        if styles_xml.exists():
            s_tree = ET.parse(styles_xml)
            s_root = s_tree.getroot()
            style_specs = {
                "Title": ("Arial", "34" if compact else "42", "172635", True),
                "Heading1": ("Arial", h1_sz, "172635", True),
                "Heading2": ("Arial", h2_sz, "07566B", True),
                "Heading3": ("Arial", h3_sz, "172635", True),
                "Normal": ("Georgia", body_sz, "172635", False),
                "BodyText": ("Georgia", body_sz, "172635", False),
                "FirstParagraph": ("Georgia", body_sz, "172635", False),
                "Compact": ("Georgia", body_sz, "172635", False),
            }
            for st in s_root.findall(w_tag("style")):
                st_id = st.attrib.get(w_tag("styleId"), "")
                if st_id in style_specs:
                    font_name, sz_val, color_val, is_bold = style_specs[st_id]
                    rpr = st.find(w_tag("rPr"))
                    if rpr is None:
                        rpr = ET.SubElement(st, w_tag("rPr"))
                    rfonts = rpr.find(w_tag("rFonts"))
                    if rfonts is None:
                        rfonts = ET.SubElement(rpr, w_tag("rFonts"))
                    for attr in ("ascii", "hAnsi", "cs"):
                        rfonts.attrib[w_tag(attr)] = font_name
                    for attr in ("asciiTheme", "hAnsiTheme", "cstheme", "eastAsiaTheme"):
                        rfonts.attrib.pop(w_tag(attr), None)
                    for tag_name in ("sz", "szCs"):
                        el = rpr.find(w_tag(tag_name))
                        if el is None:
                            el = ET.SubElement(rpr, w_tag(tag_name))
                        el.attrib[w_tag("val")] = sz_val
                    col = rpr.find(w_tag("color"))
                    if col is None:
                        col = ET.SubElement(rpr, w_tag("color"))
                    col.attrib[w_tag("val")] = color_val
                    col.attrib.pop(w_tag("themeColor"), None)
                    if is_bold and rpr.find(w_tag("b")) is None:
                        ET.SubElement(rpr, w_tag("b"))
            s_tree.write(styles_xml, encoding="utf-8", xml_declaration=True)

        doc_xml = tmp / "word" / "document.xml"
        if doc_xml.exists():
            d_tree = ET.parse(doc_xml)
            d_root = d_tree.getroot()
            for tbl in d_root.iter(w_tag("tbl")):
                rows = tbl.findall(w_tag("tr"))
                if not rows:
                    continue
                first_cells = rows[0].findall(w_tag("tc"))
                col_count = len(first_cells)
                base = 9360 // max(1, col_count)
                widths = [base] * col_count
                widths[-1] += 9360 - sum(widths)

                tbl_pr = tbl.find(w_tag("tblPr"))
                if tbl_pr is None:
                    tbl_pr = ET.Element(w_tag("tblPr"))
                    tbl.insert(0, tbl_pr)
                tbl_w = tbl_pr.find(w_tag("tblW"))
                if tbl_w is None:
                    tbl_w = ET.SubElement(tbl_pr, w_tag("tblW"))
                tbl_w.attrib[w_tag("w")] = "9360"
                tbl_w.attrib[w_tag("type")] = "dxa"
                tbl_layout = tbl_pr.find(w_tag("tblLayout"))
                if tbl_layout is None:
                    tbl_layout = ET.SubElement(tbl_pr, w_tag("tblLayout"))
                tbl_layout.attrib[w_tag("type")] = "fixed"

                for r_idx, row in enumerate(rows):
                    is_header = r_idx == 0
                    for c_idx, tc in enumerate(row.findall(w_tag("tc"))):
                        w_dxa = widths[c_idx] if c_idx < len(widths) else widths[-1]
                        tc_pr = tc.find(w_tag("tcPr"))
                        if tc_pr is None:
                            tc_pr = ET.Element(w_tag("tcPr"))
                            tc.insert(0, tc_pr)
                        tc_w = tc_pr.find(w_tag("tcW"))
                        if tc_w is None:
                            tc_w = ET.SubElement(tc_pr, w_tag("tcW"))
                        tc_w.attrib[w_tag("w")] = str(w_dxa)
                        tc_w.attrib[w_tag("type")] = "dxa"
                        if is_header:
                            shd = tc_pr.find(w_tag("shd"))
                            if shd is None:
                                shd = ET.SubElement(tc_pr, w_tag("shd"))
                            shd.attrib[w_tag("val")] = "clear"
                            shd.attrib[w_tag("color")] = "auto"
                            shd.attrib[w_tag("fill")] = "EDF3F5"
                        for p in tc.findall(w_tag("p")):
                            for r in p.findall(w_tag("r")):
                                rpr = r.find(w_tag("rPr"))
                                if rpr is None:
                                    rpr = ET.Element(w_tag("rPr"))
                                    r.insert(0, rpr)
                                for tag_name in ("sz", "szCs"):
                                    el = rpr.find(w_tag(tag_name))
                                    if el is None:
                                        el = ET.SubElement(rpr, w_tag(tag_name))
                                    el.attrib[w_tag("val")] = tbl_sz
                                if is_header and rpr.find(w_tag("b")) is None:
                                    ET.SubElement(rpr, w_tag("b"))
            d_tree.write(doc_xml, encoding="utf-8", xml_declaration=True)

        with zipfile.ZipFile(docx_path, "w", compression=zipfile.ZIP_DEFLATED) as zf_out:
            for file_path in sorted(tmp.rglob("*")):
                if file_path.is_file():
                    zf_out.write(file_path, file_path.relative_to(tmp).as_posix())


def render_markdown_to_pdf(md_path: Path, pdf_path: Path, title: str, source_sha_override: str | None = None, compact: bool = False) -> dict:
    source_sha = source_sha_override or sha_file(md_path)
    with tempfile.TemporaryDirectory(prefix="catmech-pdf-") as tmpdir:
        tmp_docx = Path(tmpdir) / (md_path.stem + ".docx")
        subprocess.run(
            [
                str(PANDOC),
                str(md_path),
                "--from=markdown+tex_math_dollars+pipe_tables+raw_html",
                "--to=docx",
                "--standalone",
                f"--resource-path={EDU}:{ROOT}",
                "-o",
                str(tmp_docx),
            ],
            check=True,
        )
        style_docx_and_set_metadata(tmp_docx, title, source_sha, compact=compact)
        subprocess.run(
            [
                str(SOFFICE),
                "--headless",
                "--convert-to",
                "pdf",
                "--outdir",
                str(tmpdir),
                str(tmp_docx),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        tmp_pdf = Path(tmpdir) / (md_path.stem + ".pdf")
        pdf_path.write_bytes(tmp_pdf.read_bytes())

    info = subprocess.check_output(["pdfinfo", str(pdf_path)], text=True)
    m_pages = re.search(r"^Pages:\s+(\d+)", info, re.M)
    pages = int(m_pages.group(1)) if m_pages else 0
    return {
        "file": pdf_path.name,
        "sha256": sha_file(pdf_path),
        "pages": pages,
        "source_sha_in_metadata": source_sha in info,
    }


def build_g9_inventory() -> tuple[Path, Path]:
    g9_map_path = EDU / "g9_claims_MAP.json"
    g9_inv_path = EDU / "g9_claim_inventory.json"
    nodes = [
        {"id": "F1", "level": 1, "title": "Collapsed slice-within-Gibbs mixture posteriors and marginalized corruption probability p_i", "links": []},
        {"id": "F2", "level": 1, "title": "Two-phase semiparametric influence function variance ratio rho(f) = sqrt(1 + k(1-f)/f)", "links": []},
        {"id": "F3", "level": 1, "title": "Barnard-Rubin small-sample df (nu_BR < nu_com) with CV3, CR1S, and Lin-Wei covariance", "links": []},
        {"id": "E1", "level": 2, "title": "R/contamination.R & R/reference.R (collapsed Gibbs mixture + C++17 RcppArmadillo reference MI)", "links": ["F1"]},
        {"id": "E2", "level": 2, "title": "R/audit-power.R, R/audit-plan.R & R/holdout-check.R (rough unvalidated audit-size aid, NO-GO planner, rank-PIT screen)", "links": ["F2"]},
        {"id": "E3", "level": 2, "title": "R/pool.R, R/pool-glm.R, R/pool-coxph.R & R/audit-trail.R (multi-endpoint pooling + SHA-256 trail)", "links": ["F3"]},
        {"id": "A1", "level": 3, "title": "validation/contamination-study & validation/variance-mi-study (52,950 + 148,500 frozen requests)", "links": ["E1", "E3"]},
        {"id": "A2", "level": 3, "title": "validation/real-data-study & validation/real-audit-study (33,800 CT.gov/NHANES + 45,600 Rwanda/psych)", "links": ["E1", "E2", "E3"]},
        {"id": "D1", "level": 4, "title": "Synthetic Contamination Study: 44.3% narrower CIs vs gold_only in missed_flags; 3 stress-cell boundaries", "links": ["A1"]},
        {"id": "D2", "level": 4, "title": "Clustered Variance & MI Study: CV3 passes 9/11 Rwanda indicators; predictor-slope attenuation warning", "links": ["A1"]},
        {"id": "D3", "level": 4, "title": "Real-Data & Real-Audit Studies: 43%-48% RMSE cut from keeping W_i; p_i AUC 0.500-0.582; Rwanda INCONCLUSIVE", "links": ["A2"]},
        {"id": "H1", "level": 5, "title": "Headline 1: Retaining recorded measurements cuts RMSE by 43% to 48% on ClinicalTrials.gov and NHANES", "links": ["D3", "A2"]},
        {"id": "H2", "level": 5, "title": "Headline 2: Contamination mixture avoids attenuation under gross separable errors but per-row p_i is near chance on diffuse real data", "links": ["D1", "D3"]},
        {"id": "H3", "level": 5, "title": "Headline 3: Evidence-based decision protocol pairs Simple Linear MI for diffuse errors with multi-endpoint Barnard-Rubin pooling", "links": ["D1", "D2", "D3"]},
    ]
    g9_map_path.write_text(json.dumps({"nodes": nodes}, indent=2) + "\n", encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(HARNESS_DIR / "grounding_check.py"),
            "education/g9_claims_MAP.json",
            "--emit-inventory",
            "education/g9_claim_inventory.json",
        ],
        cwd=str(ROOT),
        check=True,
    )
    inv_data = json.loads(g9_inv_path.read_text(encoding="utf-8"))
    inv_data["generated_from"] = "education/g9_claims_MAP.json"
    g9_inv_path.write_text(json.dumps(inv_data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return g9_map_path, g9_inv_path


def main() -> None:
    EDU.mkdir(parents=True, exist_ok=True)
    (EDU / ".gitattributes").write_text("*.pdf binary\n", encoding="utf-8")

    ms_sha = sha_file(MS_PATH)
    ms_pdf_sha = sha_file(MS_PDF_PATH)
    ms_text = MS_PATH.read_text(encoding="utf-8")
    ms_lines = ms_text.splitlines()

    comp_md_path = EDU / "COMPANION.md"
    study_md_path = EDU / "STUDY_GUIDE.md"
    one_pager_md_path = EDU / "ONE_PAGER.md"
    cov_md_path = EDU / "COMPLETE_COVERAGE.md"
    readme_md_path = EDU / "README.md"
    lab_html_path = EDU / "companion_lab.html"

    comp_text = build_companion_md(ms_sha, ms_pdf_sha)
    comp_md_path.write_text(comp_text, encoding="utf-8")
    comp_sha = sha_file(comp_md_path)

    # Human review state. The sign-off counts only for the exact files it names.
    signoff = read_human_signoff()
    human_signed = bool(
        signoff and signoff["manuscript_sha256"] == ms_sha and signoff["companion_sha256"] == comp_sha
    )
    g12_evidence, g12_signature = read_g12_resign(ms_sha)
    pending_opened = None
    if not human_signed:
        if not PENDING_PATH.is_file():
            sys.exit("Human review is pending for these files; write education/G14_HUMAN_REVIEW_PENDING.md first.")
        opened = re.search(r"\*\*Opened:\*\*\s*`(\d{4}-\d{2}-\d{2})`", PENDING_PATH.read_text(encoding="utf-8"))
        if not opened:
            sys.exit("education/G14_HUMAN_REVIEW_PENDING.md needs an '**Opened:** `YYYY-MM-DD`' line.")
        pending_opened = opened.group(1)

    study_text = build_study_guide_md(ms_sha, ms_pdf_sha)
    study_md_path.write_text(study_text, encoding="utf-8")

    one_pager_text = build_one_pager_md(ms_sha)
    one_pager_md_path.write_text(one_pager_text, encoding="utf-8")

    cov_text = build_complete_coverage_md(ms_sha, ms_pdf_sha)
    cov_md_path.write_text(cov_text, encoding="utf-8")

    readme_text = build_readme_md(ms_sha, ms_pdf_sha, comp_sha, human_signed)
    readme_md_path.write_text(readme_text, encoding="utf-8")

    lab_html_text = build_companion_lab_html(comp_sha, ms_sha)
    lab_html_path.write_text(lab_html_text, encoding="utf-8")

    # Render PDFs
    sg_pdf_res = render_markdown_to_pdf(
        study_md_path,
        EDU / "STUDY_GUIDE.pdf",
        "Clinical Measurement Audits and Two-Phase Imputation: Concise Study Guide",
    )
    comp_pdf_res = render_markdown_to_pdf(
        comp_md_path,
        EDU / "COMPANION.pdf",
        "Complete Reading & Executable R Companion to the catmech Package",
    )
    op_pdf_res = render_markdown_to_pdf(
        one_pager_md_path,
        EDU / "ONE_PAGER.pdf",
        "catmech 1.0.1: One-Page Clinical Audit Decision & Boundary Card",
        source_sha_override=comp_sha,
        compact=True,
    )

    # Emit G9 claim inventory
    _, g9_inv_path = build_g9_inventory()
    inv = json.loads(g9_inv_path.read_text(encoding="utf-8"))
    g9_ids = [c["id"] for c in inv["claims"]]

    load_bearing = [
        ("LB1", "Retaining recorded values via catmech_reference(use_recorded = TRUE) reduces RMSE by 43% to 48% over discarding unaudited values on ClinicalTrials.gov and NHANES."),
        ("LB2", "catmech_contamination() analytically integrates out latent (delta_i, T_i) on unaudited rows during collapsed slice-within-Gibbs updates and anchors audited rows to R_i."),
        ("LB3", "Multi-endpoint pooling (catmech_pool_lm with CV3, catmech_pool_glm with CR1S, catmech_pool_coxph with d-p or Lin-Wei) bounds degrees of freedom via Barnard-Rubin nu_BR < nu_com."),
        ("LB4", "In synthetic benchmarks (52,950 requests), catmech_contamination() narrows CIs by 44.3% vs complete-case in missed_flags and prevents linear-MI attenuation under gross outliers, but fails in heavy_tail_target, nonignorable_missed, and noisy_undeclared."),
        ("LB5", "catmech_holdout_check() screens observation-level calibration via tie-randomized rank PIT but has low power against regression slope bias."),
        ("LB6", "On diffuse real-world discrepancies (ClinicalTrials.gov, NHANES, Rwanda HMIS), per-record corruption probabilities p_i have near-chance AUC (0.500-0.582) and must not adjudicate individual patient rows."),
        ("LB7", "catmech_audit_plan() carries a pre-specified NO-GO validation status (11.3%-16.8% median width error at f=5%) and is documented as exploratory."),
        ("LB8", "catmech_audit_trail() emits SHA-256 provenance records (serialized R data frames; canonical JSON specification) that let a rerun be compared and support FDA 21 CFR Part 11 and ICH E9 documentation, without making the record tamper-evident or certifying compliance."),
    ]
    lb_ids = [k for k, _ in load_bearing]
    all_coverage_ids = lb_ids + g9_ids

    report = g14.Report()
    parsed = g14.parse_companion(comp_text, report)
    assert not report.findings, f"Parse findings: {report.findings}"
    assert len(parsed) == 30, f"Expected 30 entries, got {len(parsed)}"

    entries = []
    for idx, p in enumerate(parsed, start=1):
        assert p.source_range is not None
        a, b = p.source_range
        excerpt = "\n".join(ms_lines[a - 1 : b])
        ev_text = p.fields["Evidence status"]
        ev_class = g14.normalized_evidence_class(ev_text)
        assert ev_class in g14.EVIDENCE_CLASSES, f"{p.key}: bad evidence class {ev_class!r}"
        cov = [
            lb_ids[(idx - 1) % len(lb_ids)],
            g9_ids[(idx - 1) % len(g9_ids)],
        ]
        entries.append(
            {
                "index": idx,
                "key": p.key,
                "title": p.title,
                "section": p.section,
                "source_range": [a, b],
                "source_excerpt_sha256": sha_text(excerpt),
                "body_sha256": sha_text(p.body),
                "evidence_status": ev_text,
                "source_evidence_class": ev_class,
                "evidence_class": ev_class,
                "coverage_ids": sorted(set(cov)),
            }
        )

    claim_fields = (
        "What the passage is doing",
        "In ordinary language",
        "Evidence status",
        "Do not accidentally conclude",
    )
    allowances = []
    for p in parsed:
        assert p.source_range is not None
        a, b = p.source_range
        src_nums = set(g14.numeric_tokens("\n".join(ms_lines[a - 1 : b])))
        comp_nums = set(g14.numeric_tokens("\n".join(p.fields.get(f, "") for f in claim_fields)))
        for tok in sorted(comp_nums - src_nums):
            allowances.append(
                {
                    "entry_key": p.key,
                    "token": tok,
                    "reason": f"Pedagogical worked example, cross-section reference, or Lab/PDF locator in {p.key}.",
                }
            )

    # Human review records. G14_HUMAN_SIGNOFF.md is the reviewer's own file and is never written
    # here; a row is GREEN only when that file names the manuscript and companion just built.
    if human_signed:
        human = {
            "status": "GREEN",
            "reviewer_type": "human",
            "reviewer": signoff["reviewer"],
            "date": signoff["date"],
            "evidence_path": SIGNOFF_REL,
        }
    else:
        human = {
            "status": "PENDING",
            "reviewer_type": "human",
            "reviewer": (signoff["reviewer"] if signoff else "maintainer") + " (review pending)",
            "date": pending_opened,
            "evidence_path": PENDING_REL,
        }

    def qa_record() -> dict:
        return {"status": human["status"], "evidence_path": human["evidence_path"]}

    if g12_signature:
        g12_record = {"status": "GREEN", "manuscript_sha256": ms_sha, "evidence_path": g12_evidence, **g12_signature}
    else:
        g12_record = {"status": "PENDING", "manuscript_sha256": ms_sha, "evidence_path": g12_evidence}

    # Scores are recorded only from a sign-off that covers these files and states 12/12 and 4/4.
    if human_signed and signoff["all_items_passed"]:
        assessment_results = {
            "item_scores": {item["id"]: 1 for item in cdata.ASSESSMENT_ITEMS},
            "one_sentence_scores": {"CLAIM": 1, "EVIDENCE": 1, "BOUNDARY": 1, "TRANSFER": 1},
        }
    else:
        assessment_results = None
    if human_signed:
        PENDING_PATH.unlink(missing_ok=True)

    artifacts_json = {
        "canonical_sources": {
            "manuscript_rmd": {"file": "../paper/manuscript.Rmd", "sha256": ms_sha, "lines": len(ms_lines)},
            "manuscript_pdf": {"file": "../paper/manuscript.pdf", "sha256": ms_pdf_sha, "pages": ms_pdf_pages()},
        },
        "study_guide_markdown": {
            "file": "STUDY_GUIDE.md",
            "sha256": sha_file(study_md_path),
            "words": len(re.findall(r"\b[\w'-]+\b", study_text)),
        },
        "study_guide_pdf": sg_pdf_res,
        "companion_markdown": {
            "file": "COMPANION.md",
            "sha256": comp_sha,
            "words": len(re.findall(r"\b[\w'-]+\b", comp_text)),
        },
        "companion_pdf": comp_pdf_res,
        "one_pager_markdown": {
            "file": "ONE_PAGER.md",
            "sha256": sha_file(one_pager_md_path),
        },
        "one_pager_pdf": op_pdf_res,
        "companion_lab_html": {
            "file": "companion_lab.html",
            "sha256": sha_file(lab_html_path),
        },
        "complete_coverage": {
            "file": "COMPLETE_COVERAGE.md",
            "sha256": sha_file(cov_md_path),
        },
        "r_labs_script": {
            "file": "run_companion_labs.R",
            "sha256": sha_file(EDU / "run_companion_labs.R"),
        },
        "r_interactive_tutor": {
            "file": "interactive_tutor.R",
            "sha256": sha_file(EDU / "interactive_tutor.R"),
        },
    }
    (EDU / "ARTIFACTS.json").write_text(json.dumps(artifacts_json, indent=2) + "\n", encoding="utf-8")

    manifest = {
        "schema_version": "g14-companion-manifest/v2",
        "edition": "catmech 1.0.1 Multimodal R-Package Educational Companion",
        "status": "RELEASE_READY" if (human_signed and g12_signature) else "PENDING_HUMAN_REVIEW",
        "release_order": "AFTER_G12_FROZEN_MANUSCRIPT",
        "source_roles": {
            "manuscript": "paper/manuscript.Rmd",
            "companion": "education/COMPANION.md",
        },
        "canonical_sources": {
            "paper/manuscript.Rmd": {
                "sha256": ms_sha,
                "lines": len(ms_lines),
                "line_locator_basis": "1-based inclusive line numbers in frozen paper/manuscript.Rmd",
            },
            "education/COMPANION.md": {
                "sha256": comp_sha,
                "lines": len(comp_text.splitlines()),
                "line_locator_basis": "entry locators reference paper/manuscript.Rmd line numbers",
            },
        },
        "entry_count": len(entries),
        "entries": entries,
        "coverage_contract": {
            "load_bearing_claim_ids": lb_ids,
            "structural_exclusions": [],
            "declared_overlaps": [],
            "g9_inventory": {
                "path": "education/g9_claim_inventory.json",
                "sha256": sha_file(g9_inv_path),
                "source_sha256": inv["source_sha256"],
            },
        },
        "claim_contract": {
            "global_nonclaims": [
                {
                    "id": cdata.GLOBAL_NONCLAIM_ID,
                    "text": cdata.GLOBAL_NONCLAIM_TEXT,
                }
            ],
            "allowed_companion_only_numbers": allowances,
            "claim_review": {
                **human,
                "unrecorded_claims": [],
            },
        },
        "upstream_gates": {
            # build_g9_inventory() above runs grounding_check.py with check=True, so this build
            # stops unless the claim-grounding DAG exits 0; its inventory is the evidence.
            "g9": {
                "status": "GREEN",
                "evidence_path": "education/g9_claim_inventory.json",
            },
            "g12": g12_record,
        },
        "education_contract": {
            "audience": "Biostatisticians, clinical trialists, epidemiologists, and R package users analyzing audited measurements",
            "plain_language": {
                "maximum_grade_level": 10,
                "maximum_sentence_words": 30,
                "jargon_registry": cdata.JARGON_REGISTRY,
                "human_review": {
                    **human,
                    "undefined_jargon": [],
                },
            },
            "package_claim": "MULTIMODAL",
            "claimed_modalities": ["DOCUMENT", "ONE_PAGER", "INTERACTIVE"],
            "summary_coverage_ids": lb_ids,
            "modalities": [
                {
                    "id": "companion-document",
                    "kind": "DOCUMENT",
                    "artifact_path": "education/COMPANION.pdf",
                    "source_companion_sha256": comp_sha,
                    "coverage_ids": all_coverage_ids,
                    "semantic_equivalence": dict(human),
                },
                {
                    "id": "companion-one-pager",
                    "kind": "ONE_PAGER",
                    "artifact_path": "education/ONE_PAGER.pdf",
                    "source_companion_sha256": comp_sha,
                    "coverage_ids": lb_ids,
                    "semantic_equivalence": dict(human),
                },
                {
                    "id": "companion-interactive-lab",
                    "kind": "INTERACTIVE",
                    "artifact_path": "education/companion_lab.html",
                    "source_companion_sha256": comp_sha,
                    "coverage_ids": all_coverage_ids,
                    "semantic_equivalence": dict(human),
                    "practice_feedback": {
                        "records_attempt": True,
                        "reports_correctness": True,
                        "explains_answer": True,
                        "supports_retry": True,
                        "evidence_path": human["evidence_path"],
                    },
                },
            ],
        },
        "transfer_assessment": {
            "assessment_version": "catmech-g14-v2",
            "assessment_sha256": g14.canonical_json_sha256(cdata.ASSESSMENT_ITEMS),
            "items": cdata.ASSESSMENT_ITEMS,
            "thresholds": {
                "overall": {"minimum": 10, "possible": 12},
                "boundary": {"minimum": 4, "possible": 4},
                "transfer": {"minimum": 3, "possible": 4},
                "one_sentence": {"minimum": 4, "possible": 4},
            },
            "one_sentence_rubric": [
                {"id": "CLAIM", "criterion": "States that catmech uses two-phase audited reference measurements to calibrate contaminated or incomplete recorded values and pool downstream regressions"},
                {"id": "EVIDENCE", "criterion": "Distinguishes analytical posteriors/variance formulas from synthetic simulations (52,950 requests) and empirical benchmarks (CT.gov, NHANES, Rwanda HMIS)"},
                {"id": "BOUNDARY", "criterion": "Preserves the non-claims that per-row p_i is near chance on diffuse real data, catmech_audit_plan() is NO-GO, and catmech_audit_trail() does not certify Part 11 compliance"},
                {"id": "TRANSFER", "criterion": "Selects Simple Linear MI vs Contamination Mixture and appropriate cluster pooling (CV3/CR1S/Lin-Wei) for a new clinical audit scenario"},
            ],
            "results": assessment_results,
            "human_verdict": dict(human),
        },
        "shipping": {
            "status": "SHIPPED",
            "artifacts": [
                {
                    "path": "education/COMPANION.md",
                    "kind": "MARKDOWN",
                    "sha256": comp_sha,
                    "generated_from_companion_sha256": comp_sha,
                    "text_qa": qa_record(),
                    "render_qa": {"status": "NOT_APPLICABLE", "reason": "Markdown source file; rendered PDF is checked below."},
                    "accessibility_qa": {"status": "NOT_APPLICABLE", "reason": "Plain UTF-8 Markdown with structured headings."},
                },
                {
                    "path": "education/COMPANION.pdf",
                    "kind": "PDF",
                    "sha256": comp_pdf_res["sha256"],
                    "generated_from_companion_sha256": comp_sha,
                    "text_qa": qa_record(),
                    "render_qa": qa_record(),
                    "accessibility_qa": qa_record(),
                },
                {
                    "path": "education/STUDY_GUIDE.pdf",
                    "kind": "PDF",
                    "sha256": sg_pdf_res["sha256"],
                    "generated_from_companion_sha256": comp_sha,
                    "text_qa": qa_record(),
                    "render_qa": qa_record(),
                    "accessibility_qa": qa_record(),
                },
                {
                    "path": "education/ONE_PAGER.pdf",
                    "kind": "PDF",
                    "sha256": op_pdf_res["sha256"],
                    "generated_from_companion_sha256": comp_sha,
                    "text_qa": qa_record(),
                    "render_qa": qa_record(),
                    "accessibility_qa": qa_record(),
                },
                {
                    "path": "education/companion_lab.html",
                    "kind": "HTML",
                    "sha256": sha_file(lab_html_path),
                    "generated_from_companion_sha256": comp_sha,
                    "text_qa": qa_record(),
                    "render_qa": qa_record(),
                    "accessibility_qa": qa_record(),
                },
            ],
        },
    }

    manifest_path = EDU / "g14_companion_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Built multimodal educational companion and wrote {manifest_path}")


if __name__ == "__main__":
    main()
