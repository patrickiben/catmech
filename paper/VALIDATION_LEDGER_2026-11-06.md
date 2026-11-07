<!-- gate-basis: paper/manuscript.Rmd@c56d286f2042992f47d072c95c68c7b343ded835e610c5a72b753e59549bb8b2 -->
<!-- gate-basis: paper/manuscript.tex@9c198a51f0c6628594353e9a528271bc5a758e2fbebcc0b60b24e6f6a8b90c04 -->
<!-- gate-basis: paper/manuscript.pdf@7d23dac5f1b2ba59b6edf7f37df934e68009f77c3b08da57f889f584c091c200 -->
<!-- gate-basis: paper/paper.md@113f36ee41bbfe31ffc97efff4a318e4b18798878a9f3b9b691ec656ac79236e -->
<!-- gate-basis: paper/paper.bib@da45158e93677cfff02ad2b006d7bbd9a42b9443d38141d5cc6b6d447c92bb88 -->
<!-- gate-basis: paper/verify_paper_numbers.py@d9a2facf11d9d1dbae6d328c149a4b757c7a78b2e9ce132b154268e7c5582639 -->
<!-- gate-basis: README.md@025b4c2030e41f168eaa1b05c58c24750d396eaea9a1647e7477c054ca2eba21 -->
# Pre-Submission Validation Harness Ledger (`catmech` v1.0.0)

- **Date**: 2026-11-07
- **Target**: `catmech` v1.0.0 CRAN package & JSS / JOSS manuscripts (`paper/manuscript.Rmd`, `paper/manuscript.tex`, `paper/manuscript.pdf`, `paper/paper.md`, `paper/paper.bib`)
- **Harness Specification**: `/Users/patrickiben/Documents/Documents - Patrick’s Mac Studio/Research_Validation_Harness/PRE_SUBMISSION_VALIDATION_HARNESS.md`
- **specification_coverage**: `lines 1-590 (100%)`
- **Total Gates in Specification**: `24` (`G0`–`G23`, verified via `grep -cE '^## G[0-9]+:' PRE_SUBMISSION_VALIDATION_HARNESS.md`)
- **Execution Tier**: Full Deterministic Tier + 5-Agent Isolated Read-Only Lean Reader Wave (`G5` Persona A & B, `G12`, `G13`, `G18`)

---

## Gate Summary Table (All 24 Gates)

| Gate | Title | Status | Tool / Method | Verdict & Evidence |
| :--- | :--- | :---: | :--- | :--- |
| **G0** | Format-conversion pre-gate | **RUN** | Pandoc / RMarkdown inspection | **PASS**: Repaired 4 list blocks lacking preceding blank lines and institutional BibTeX braces (`fda2003part11`, `ich1998e9`, `nhanes2018`). |
| **G1** | Path, hygiene & repository-surface guard | **RUN** | `check_paths.sh` + `R CMD check` + `verify-release-archive.py` | **PASS**: 0 machine paths in shippable package/paper (`R/`, `src/`, `man/`, `vignettes/`, `inst/`, `tests/`, `tools/`, `paper/`, `README.md`, `DESCRIPTION`, `NAMESPACE`). 2 historical paths in `.Rbuildignore`-excluded `validation/clustered-fix-scratch/` left untouched per Rule 3 (frozen archive). |
| **G2** | Numeric-claim <-> artifact check | **RUN** | `tools/verify-artifacts.py` + `paper/verify_paper_numbers.py` | **PASS**: Grounded all synthetic, real-data, clustered, and real-audit numbers in `manuscript.Rmd` and `paper.md` to frozen CSVs in `validation/`. 411 claims + negative control verified in `paper/verify_paper_numbers.py`. |
| **G3** | Citation & bibliography verification | **RUN** | `check_citations.py` (per-project copy; not in `Research_Validation_Harness/`) on `paper/paper.bib`, re-run 2026-11-07 against Crossref | **PASS (`0 FLAG`, `0 uncited`)**: Fixed 5 broken/hijacked DOIs/entries (`huang2013simple`, `breslow2003semiparametric`, `neyman1938contribution`, `little2019statistical`, `nab2021mecor`), added URLs to `leeder2008simex`, `fda2003part11`, and `ich1998e9`, and cited all 25 entries. 2026-11-07 re-run: `21 OK | 4 REVIEW | 0 FLAG`. Hand-checked REVIEW items: `cochran1977sampling` (3rd ed. 1977 is correct; Crossref matched the 1954 edition), `fda2003part11` and `nhanes2018` (agency documents without DOIs), and `leeder2008simex`, corrected from R News 8(2) 2008 to R News 6(4) 2006, pp. 26–31. |
| **G4** | Mathematical proof & derivation audit | **RUN** | `tools/check-contamination-oracle.R`, `tools/check-reference-oracle.R`, `tools/check-cluster-algebra.R` + `G5-A` | **PASS**: Fixed convolved variance $\sigma_T^2 + \sigma_e^2$ in $p_i$, collapsed slice-within-Gibbs description, exact $\nu_{com}$ rules ($\nu_{com} = N - p$ unclustered vs. $\nu_{com} = C - 1$ clustered), CR1S finite-sample factor, and tie-randomized rank-PIT formula. |
| **G5** | Adversarial referee panel (Persona A: Theory & Math; Persona B: Empirical & Biostats) | **RUN** | 2 isolated `research` subagents (`f61a45e3-7738-463d-899e-e68e614db8f7`, `f23db6b4-7a7b-4422-aed7-7dc817b4e7ca`) | **PASS (All 11 MAJOR + 10 MINOR findings resolved)**: Fixed mathematical equations in Section 2 and added complete disclosures of stress-cell failures, `REAL_AUDIT_SYNTHESIS.md` (`INCONCLUSIVE` Rule 2 on Rwanda HMIS, 35.4%–71.6% coverage), per-record AUC (`0.543` / `0.500`), predictor-slope attenuation (`-16.1%`), and holdout alarm power. |
| **G6** | Figure, table & layout inspection | **RUN** | `pdflatex` log + visual page-by-page inspection (`read_file` on `manuscript.pdf` pages 1–15) | **PASS**: `errors=0 undefined=0 overfull=0` across all 15 pages. Removed `toc: true` from `pdf_document`; author header stacked cleanly; all tables and code blocks fit within margins. |
| **G7** | AI-artifact & placeholder scrubber | **RUN** | `check_prose_style_v2.py` (`model-artifact`, `placeholder`, `markdown-in-latex`) + `scrub_ai_artifacts.py` (per-project copy) incl. `--pdf manuscript.pdf` | **PASS**: `0` hits across `manuscript.Rmd`, `manuscript.tex`, `manuscript.pdf`, `paper.md`, `README.md`, `NEWS.md`, `DESCRIPTION`. |
| **G8** | Anonymization & blind-review check | **NOT RUN** | Target venues (CRAN / JSS / JOSS) are non-blinded | **N/A (Non-blinded venue)**: CRAN, JSS, and JOSS require full author attribution, affiliation, and ORCID. |
| **G9** | Reproducibility & environment audit | **RUN** | `tools/verify-artifacts.py` + `tools/verify-release-archive.py` + `R CMD check` | **PASS**: All 10 frozen study archives, seed streams, SHA-256 source manifests, and `R CMD check --as-cran` (`0 errors, 0 warnings`; 3 NOTEs: new submission, local HTML Tidy version, local `.DS_Store` in the check directory) verified. `tools/verify-release-archive.py` applies only to the v0.10 confirmation freeze; the v1.0.0 sources intentionally differ from that freeze. |
| **G10** | Prose style, humanness & US spelling | **RUN** | `check_prose_style_v2.py` + `check_spelling_us.py` | **PASS (`ALL GREEN`)**: `0` British spellings; `0` em-dashes; `0` curly quotes; `0` Title-Case headings; all 12 lexical/artifact categories `[ok]`. |
| **G11** | Cross-reference, symbol & notation consistency | **RUN** | LaTeX compile (`undefined=0`) + `G13` symbol audit | **PASS**: Zero undefined references or citations; all symbols ($\phi, \hat{k}, G_i, a_s, n_s, u_i, \nu_{com}, \nu_{BR}$) and strategy/scenario labels defined on first use. |
| **G12** | Span-attributed prose verification | **RUN** | Isolated `research` subagent (`d62a0643-a25f-4a48-8765-82356471d6b7`) + `paper/verify_paper_numbers.py` | **PASS (All 9 divergences resolved)**: Fixed `catmech_audit_trail()` diagnostics description (`convergence_passed` and `mcmc_iterations`), listed all 10 CRAN `Imports`, corrected `mice`/`simex`/`mecor`/`Amelia` table entries, and clarified `catmech_compare_strategies()` scope. |
| **G13** | Cold-reader comprehension check | **RUN** | Isolated `research` subagent (`3b976727-7ce6-4239-a0e2-f12e35580fa6`) reading only `manuscript.Rmd` & `paper.md` | **PASS**: Cold reader reconstructed all estimands, sample sizes, decision rules, and screen failures; all symbol, label, table-header (`| Dataset | Metric | ...`), and row-ordering hazards fixed. |
| **G14** | Ethical, IRB & data-use compliance | **RUN** | Public dataset manifest & provenance check (`validation/real-data/fetch-manifest.json`, `validation/real-audit/fetch-manifest.json`) | **PASS**: All empirical benchmarks use public, de-identified datasets (ClinicalTrials.gov, CDC NHANES, Dryad/PLOS ONE Rwanda HMIS recount, PLOS ONE psychology recounts) with frozen SHA-256 manifests. |
| **G15** | Related-work completeness & fairness | **RUN** | Software comparison table audit (`mice`, `simex`, `mecor`, `Amelia`) in both papers | **PASS**: Accurately credits `mice` (Barnard-Rubin df), `simex` (SIMEX for LM/GLM), `mecor` (regression calibration for continuous linear models), and `Amelia` (observation priors and `overimpute`). |
| **G16** | Limitations & scope honesty | **RUN** | Sections 2.3, 5, 6, 8 & `README.md` (`Limits` section) audit | **PASS**: Explicitly discloses `catmech_audit_plan()` `NO-GO` verdict, 3 stress-cell failures, Rwanda HMIS `INCONCLUSIVE` Rule 2 verdict, per-record error probability AUC failure (`0.500`–`0.582`), `catmech_holdout_check()` low power, and `catmech_reference(cluster = )` predictor-slope attenuation. |
| **G17** | Pre-committed analysis-plan honesty | **RUN** | Cross-check against `CONTAMINATION_PLAN.md`, `REAL_DATA_PLAN.md`, `VARIANCE_MI_PLAN.md`, `REAL_AUDIT_PLAN.md` | **PASS**: Explicitly discloses the `site_shift_declared` Wilson screen miss (`93.7%`, `[0.9202, 0.9505]`), the 3 synthetic stress-cell failures, the `NO-GO` on `catmech_audit_plan()`, the `PARTIAL` accuracy verdict on clustered linear MI, and the `INCONCLUSIVE` Rule 2 verdict on Rwanda HMIS. |
| **G18** | Venue fit & AE-surface check | **RUN** | Isolated `research` subagent (`1c9f1bd5-3899-4046-be19-54d73387d64d`) | **PASS (All AE-BLOCKER & AE-FRICTION items resolved)**: Added worked 8-stage R workflow and software comparison table to `manuscript.Rmd`, removed `toc: true` from `pdf_document`, removed HTML comment from `paper.md`, and cited all 25 `.bib` entries. |
| **G19** | Code-surface <-> paper-surface consistency | **RUN** | `verify_code_surface()` in `paper/verify_paper_numbers.py` | **PASS**: All 52 function, script, and CRAN dependency surface checks pass across `manuscript.Rmd` and `paper.md`. |
| **G20** | Git & release-tag hygiene | **RUN** | `git status` + `.Rbuildignore` + `.gitignore` + `tools/verify-release-archive.py` | **PASS**: Clean working tree; frozen source hashes verified (`source-hashes-before.csv` intact). |
| **G21** | Supplement <-> main-text parity | **RUN** | Parity check across `manuscript.Rmd`, `paper.md`, `README.md`, and `docs/case-study/fact_check.py` | **PASS**: `118 claims, 118 match` in `docs/case-study/fact_check.py`; `411` checks pass in `paper/verify_paper_numbers.py`. |
| **G22** | Literature & empirical-calibration sourcing | **RUN** | `paper/verify_paper_numbers.py` + `tools/verify-artifacts.py` | **PASS**: All simulation and empirical numbers trace to frozen CSVs in `validation/`. |
| **G23** | Final pre-flight sign-off & basis stamping | **RUN** | `check_gate_basis.py --stamp` | **PASS**: Basis hashes stamped on this ledger. |

---

## Detailed Gate Audit Notes & `unchecked_directions` (Rule 12b)

### G5: Adversarial Referee Panel (2 Isolated Subagents)
- **Persona A (`f61a45e3-7738-463d-899e-e68e614db8f7`, Theory & Math Referee)**:
  - *Resolved MAJOR-1*: Replaced $f_\epsilon(W_i - \hat{T}_i)$ with the exact convolved Gaussian predictive density $\phi(W_i; \, X_i\beta + \mu_e, \, \sigma_T^2 + \sigma_e^2)$ in $p_i$.
  - *Resolved MAJOR-2*: Clarified that `catmech_contamination()` is implemented in R via a collapsed slice-within-Gibbs sampler under a Jeffreys prior on $(\beta, \sigma_T^2)$, Beta prior on $\pi$, and normal-scaled-inverse-$\chi^2$ prior on $(\mu_e, \sigma_e^2)$, while C++17 `RcppArmadillo` accelerates `catmech_reference()`.
  - *Resolved MAJOR-3*: Clarified that `catmech_audit_power()` uses log-normal Monte Carlo draws of $\hat{k}$ (`n_sim = 2000`) rather than closed-form non-central approximations, and disclosed the `NO-GO` validation status of `catmech_audit_plan()`.
  - *Resolved MAJOR-4 & MAJOR-5*: Clarified that `catmech_pool_lm()`, `catmech_pool_glm()`, and `catmech_pool_coxph()` use $\nu_{com} = N - p$ (or $\nu_{com} = d - p$ with model-based inverse Fisher information for unclustered Cox PH) when `cluster = NULL`, and set $\nu_{com} = C - 1$ when `cluster` is specified across $C$ sites (with CV3 jackknife for LM, CR1S finite-sample-corrected Liang-Zeger sandwich for GLM, and cluster-summed Lin-Wei sandwich for Cox PH).
  - *unchecked_directions*: Did not inspect `R/reference.R` or `src/*.cpp` directly during the isolated G5-A pass (covered by `tools/check-reference-oracle.R`).
- **Persona B (`f23db6b4-7a7b-4422-aed7-7dc817b4e7ca`, Empirical & Biostats Referee)**:
  - *Resolved MAJOR-1*: Added Section 6.2 (`Real-audit regression slope recovery benchmarks`) reporting the 45,600-request study (`validation/REAL_AUDIT_SYNTHESIS.md`) on Rwanda HMIS and 498 psychology effect-size recounts (`@maassen2020reproducibility`), including the pre-specified `INCONCLUSIVE` Rule 2 verdict and clustered coverage failures (`71.6%` on `PNC`, `66.6%` on `Syphilis`, `35.4%` on `Iron`, `49.7%` on `ANC 4th visit`).
  - *Resolved MAJOR-2*: Added explicit disclosure of the three pre-specified stress-cell failures in `validation/CONTAMINATION_SYNTHESIS.md` (`heavy_tail_target`: bias $+0.0330$, coverage $87.3\%$; `nonignorable_missed`: bias $-0.0506$, coverage $79.3\%$; `noisy_undeclared`: bias $-0.0261$, coverage $88.7\%$) plus borderline/refusal cells (`error_bimodal`, `flag_from_recorded`, `reference_05pct`).
  - *Resolved MAJOR-3 & MAJOR-5*: Disclosed in Section 5.2 that `cv3` restored 100% complete-data recount coverage (`gold_only`) to $92.0\%$–$97.5\%$ and narrowly passed under random-intercept outcome imputation ($9/11$ items $\ge 92\%$), did not rescue unclustered linear MI on Rwanda items ($81.0\%$–$86.1\%$ on worst items), and that `catmech_reference(cluster = )` attenuated predictor slopes by up to $-16.1\%$ ($84.4\%$ coverage at 25% audit).
  - *Resolved MAJOR-4 & MAJOR-6*: Disclosed that per-record corruption probabilities $p_i$ had near-chance AUC on real data (`0.543` on ClinicalTrials.gov, `0.500` on NHANES, `0.501`–`0.582` on Rwanda HMIS) and that `catmech_holdout_check()` had low alarm power against slope bias ($1/40$ in `heavy_tail_target`, $0/40$ in `nonlinear_target`, $18/40$ in `site_shift_undeclared`).
  - *Resolved MINOR 1–5*: Added `Naive (Trust Rec.)` (`naive_recorded`) to the Section 6.1 table, explained the 91.2% Oracle coverage on ClinicalTrials.gov, reported the 33,800-request refusal breakdown, and clarified the 14 refusals (12 rank-deficient, 2 zero-residual-variance) in the 148,500-request study.
  - *unchecked_directions*: Other 6 validation archives (`CLUSTERED_SYNTHESIS.md`, `OUTSIDE_METHODS_SYNTHESIS.md`, etc.) are verified deterministically via `tools/verify-artifacts.py`.

### G12: Span-Attributed Prose Verification (`d62a0643-a25f-4a48-8765-82356471d6b7`)
- **Findings fixed**: Corrected `catmech_audit_trail()` description to match `R/audit-trail.R` (`convergence_passed` and `mcmc_iterations` on the trail, with $\hat{R}$ and ESS in `result$model_diagnostics`); listed all 10 non-base CRAN `Imports` (`Rcpp`, `RcppArmadillo`, `digest`, `jsonlite`, `mice`, `pan`, `posterior`, `survival`, `truncnorm`, `yaml`); corrected the software comparison table for `mice` (`BR df`), `simex` (`SIMEX` for LM and GLM), `mecor` (`Reg. Cal.` for LM, `No` for GLM), and `Amelia` (`Priors` and `overimpute`); and clarified `catmech_compare_strategies()` scope.
- **unchecked_directions**: Formal argument defaults of un-audited `R/*.R` files beyond `NAMESPACE` and `R/audit-trail.R` / `R/strategy-compare.R` (verified via `R CMD check`).

### G13: Blinded Cold-Reader Comprehension Check (`3b976727-7ce6-4239-a0e2-f12e35580fa6`)
- **Findings fixed**: Defined all strategy codes (`v02_mask_normal`, `gold_only`, `reference_normal`, `cole_mi`, `raking`), scenario codes (`clean`, `flagged_errors`, `missed_flags`, `site_shift_declared`, `error_gross_and_small`), variance codes (`cr1s`, `cr2`, `cv3`), and mathematical symbols on first use; explained the 52,950 request total (183 method-scenario cells); replaced the ambiguous `| Dataset & Metric | Audit 25% |` header with `| Dataset | Metric | Naive (Trust Rec.) | ...`; and harmonized strategy row order across both panels of the synthetic table.
- **unchecked_directions**: None (strict Rule 12 two-file isolation).

### G18: JSS / JOSS Associate Editor Surface Review (`1c9f1bd5-3899-4046-be19-54d73387d64d`)
- **Findings fixed**: Added Section 4 (`Illustrative clinical trial workflow`) and the software comparison table to `paper/manuscript.Rmd`; removed `toc: true` from `pdf_document` in `paper/manuscript.Rmd`; removed the `<!-- gate-basis: ... -->` HTML comment from `paper/paper.md`; cited `@nhanes2018` and `@maassen2020reproducibility` in both papers; and added URLs to `@leeder2008simex`, `@fda2003part11`, and `@ich1998e9` in `paper/paper.bib`.
- **unchecked_directions**: Live `openjournals/inara` Docker build (not installed locally; Pandoc LaTeX/HTML builds verified).

---

## Companion Interactive Web Surfaces (`App_Validation_Harness`)

- **Harness Specification**: `/Users/patrickiben/Documents/Documents - Patrick’s Mac Studio/App_Validation_Harness/HARNESS.md` (`web_gate.py`, gates `A0`–`A10`)
- **`catmech_guide` (`docs/learn/index.html`, page hash `d94e13a7b66c865e`)**: All applicable gates (`A0`–`A5`, `A7`–`A10`) **GREEN**; `A0`, `A4`, and `A8` signed off (`approved (2026-11-07)`); **`OWNER-APPROVED: True`**.
- **`catmech_simulator` (`docs/learn/simulator.html`, page hash `80cafcd68fec48c7`)**: All applicable gates (`A0`–`A5`, `A7`–`A9`) **GREEN**; `A0`, `A4`, and `A8` signed off (`approved (2026-11-07)`); **`OWNER-APPROVED: True`**.
- **JS-to-R Numerical Cross-Checks**: `node docs/learn/qa/compare.js` (`140` studies + `6` negative controls: **PASS**) and `node docs/audit-planner/qa/check-k-port.js` (`3` browser test vectors + pool $\hat{k}$ to `1e-11`–`1e-15`: **PASS**).

---

## Independent Re-QC (2026-11-07)

Fresh pass that re-ran every deterministic check instead of inheriting earlier results.

- **Fixed: wrong citation metadata.** `leeder2008simex` cited R News 8(2), 2008; the article is Lederer and Küchenhoff, R News 6(4), 2006, pp. 26–31.
- **Fixed: regulatory overclaim.** The abstract said SHA-256 provenance "meeting FDA 21 CFR Part 11 and ICH GCP E9 regulatory compliance standards", and the comparison table read "Yes (Part 11)". Hashing a dataset supports documentation, but compliance depends on the validated system, access controls, and procedures. E9 is the statistical-principles guideline, not GCP (E6). The abstract, Section 9, Table 1, the workflow comment and keywords in `manuscript.Rmd` and `paper.md`, plus the vignette, now say the report *supports* documentation and does not certify compliance. The vignette's "ICH E6" was corrected to E9.
- **Fixed: false check claims.** `manuscript.Rmd`, `paper.md`, and `NEWS.md` said `R CMD check` returned `Status: OK (0 errors, 0 warnings, 0 notes)`. A new submission always has at least one NOTE, so they now say 0 errors and 0 warnings.
- **Fixed: maturity contradiction.** `NEWS.md` said the package "graduates from experimental status"; README, the guide, and the manuscript still mark the clustered and joint-marginal options as experimental and the planner as NO-GO. `NEWS.md` now says this.
- **Fixed: DESCRIPTION reference.** The Breslow and Chatterjee (1999) citation added on 2026-11-07 described "two-phase study design calibration", which the package does not implement, and that reference is not in `paper.bib`. It was replaced with Neyman (1938), which `paper.bib` already cites, for audit-size planning from classical two-phase sampling theory.
- **Ledger corrections.** G3 and G7 name the script copies that were actually run, which live in per-project folders and not in `Research_Validation_Harness/`. G9 reports the real NOTE count.
- **Not changed (author decision):** `R/audit-trail.R` prints a heading "Regulatory Audit & Compliance Certification" and defaults its title to "Catmech Audit Trail & Compliance Report". This is package code, so it was left alone.
