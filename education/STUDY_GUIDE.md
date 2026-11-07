# Clinical Measurement Audits and Two-Phase Imputation: Concise Study Guide

A twelve-lesson study guide to Patrick Iben's *Statistical Inference and Two-Phase Audit Design for Contaminated Measurements in Clinical Research: The catmech R Package*, canonical manuscript [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`350` lines, SHA-256 `8152288b1c2af2a49103543db3c9324582177b1e842365b5c5ab83cdac8c1fd0`) and compiled 15-page PDF [paper/manuscript.pdf](../paper/manuscript.pdf) (SHA-256 `1d17480886d4f38c38e7cf5f900ea3c5f1ef48f3c4484fa801a72a83de05e0db`).

This educational companion explains the statistical definitions, benchmark archives, and R workflows in the frozen catmech manuscript; it does not establish universal interval coverage under unmodeled misspecification, validate per-record error probabilities for individual patient adjudication on diffuse real-world data, or certify FDA 21 CFR Part 11 or ICH E9 regulatory compliance.

## How to use this twelve-lesson route

This concise guide is designed for clinicians, trial statisticians, and data auditors who want the core mathematical logic, empirical findings, R function calls, and scope boundaries in twelve short lessons (`S1` through `S12`). You can also step through these exact twelve lessons interactively in an R console by running `source('education/interactive_tutor.R'); catmech_tutor()`, or read the full thirty-unit passage-by-passage companion in [COMPANION.md](COMPANION.md) ([PDF](COMPANION.pdf)).

## S1. The three-way clinical audit dilemma and why recorded values matter

**Where we are**

Abstract and Section 1, manuscript lines 1-41 (PDF pages 1-2).

**What the passage is doing**

When an audit re-measures a subset of records (`R_i`), trialists face three flawed conventional choices: (1) naively trusting recorded values (`W_i`) causes attenuation bias; (2) analyzing only audited rows (`gold_only`) inflates standard errors; and (3) setting unaudited `W_i` to `NA` (`reference_normal`) wastes accurate recordings. On ClinicalTrials.gov (`N = 27,208`) and NHANES (`N = 3,908`), keeping `W_i` via `catmech_reference(use_recorded = TRUE)` reduces RMSE by 43% to 48% over discarding unaudited values.

**In ordinary language**

Audits check a small slice of records against a trusted reference. Keeping recorded numbers cuts error nearly in half compared with throwing unchecked numbers away.

**Evidence status**

SCOPE and MEASURED_DATED summary across 10 frozen benchmark archives.

**Do not accidentally conclude**

Do not conclude that standard MAR imputation (`mice`, `Amelia`) is biased when unaudited `W_i` values are discarded; its penalty in that mode is wider confidence intervals.

**Try it**

Run `round(100 * (1 - 0.0409 / 0.0785), 1)` in R to verify the NHANES RMSE reduction at a 25% audit. **Answer:** `47.9%`.

## S2. Two-phase contamination mixture model and marginalized corruption probabilities

**Where we are**

Sections 2.1-2.2, manuscript lines 42-70 (PDF pages 2-3).

**What the passage is doing**

Latent true values follow `T_i | X_i ~ N(X_i beta, sigma_T^2)` and recorded values follow `W_i = (1 - delta_i) T_i + delta_i (T_i + epsilon_i)` with `delta_i in {0, 1}` and `epsilon_i ~ N(mu_e, sigma_e^2)`. On unaudited rows, `catmech_contamination()` integrates out `(delta_i, T_i)` to compute the exact marginal corruption probability `p_i`. When `delta_i = 0`, `T_i = W_i` deterministically; when `delta_i = 1`, `T_i` is drawn from the precision-weighted conjugate Gaussian posterior.

**In ordinary language**

Each row has a hidden true value and a recorded value. The mixture model compares how well an unchecked recording fits the clean curve versus the wider error curve.

**Evidence status**

DEFINITION and ANALYTIC conditional posteriors for `catmech_contamination()`.

**Do not accidentally conclude**

Do not conclude that `reference_sd` (`sigma_ref`) is estimated from the data; it defaults to `0` (exact recounts) and is a user-declared sensitivity constant.

**Try it**

When `delta_i = 0` is drawn for an unaudited record in `catmech_contamination()`, what is the imputed value `T_i`? **Answer:** `T_i = W_i` deterministically.

## S3. Gaussian reference imputation in C++17: use_recorded = TRUE versus FALSE

**Where we are**

Section 2.2, manuscript lines 71-72 (PDF pages 3-4).

**What the passage is doing**

`catmech_reference()` implements fast conjugate Gaussian reference imputation in C++17 (`RcppArmadillo`). With `use_recorded = TRUE` (Simple Linear MI, `cole_mi`), `W_i` enters directly as a continuous linear calibration predictor alongside `X_i` (`T_i | W_i, X_i`). With `use_recorded = FALSE` (`reference_normal`), unaudited `W_i` values are discarded.

**In ordinary language**

When errors are small and spread out, a fast linear regression on the audited rows calibrates the unchecked numbers without needing a mixture model.

**Evidence status**

DEFINITION of `catmech_reference(use_recorded = TRUE)` versus `use_recorded = FALSE`.

**Do not accidentally conclude**

Do not forget to pass `use_recorded = TRUE` explicitly when you want Simple Linear MI; the formal function default is `use_recorded = FALSE` for v0.3 backward compatibility.

**Try it**

In `Lab C16` (`education/run_companion_labs.R`), why does `use_recorded = TRUE` yield a smaller pooled standard error than `use_recorded = FALSE`? **Answer:** Because `W_i` explains substantial residual variance in `T_i` on unaudited rows.

## S4. Prospective two-phase audit power planning and the NO-GO planner boundary

**Where we are**

Section 2.3, manuscript lines 73-84 (PDF page 4).

**What the passage is doing**

At audit fraction `f = n_audit / N`, the confidence interval width ratio relative to a 100% census audit is `rho(f) = sqrt(1 + k (1 - f) / f)`. `catmech_audit_power()` simulates `2,000` log-normal draws of `k_hat` to compute the power to achieve `rho <= rho_target`. In pre-specified validation (`validation/audit-plan-study`), `catmech_audit_plan()` received a `NO-GO` verdict (`11.3%-16.8%` median width error at `f = 5%`) and is documented as exploratory.

**In ordinary language**

The power calculator sizes an audit by simulating uncertainty in the two-phase noise factor. The empirical pilot planner missed its 10% accuracy target and stays exploratory.

**Evidence status**

ANALYTIC influence-function ratio and MEASURED_DATED `NO-GO` boundary on `catmech_audit_plan()`.

**Do not accidentally conclude**

Do not treat `catmech_audit_plan()` as a validated regulatory sample-size guarantee; its pre-specified status is `NO-GO`.

**Try it**

Compute `rho(0.25) = sqrt(1 + 0.27 * (1 - 0.25) / 0.25)` in R (`round(sqrt(1 + 0.27 * 3), 4)`). **Answer:** `1.3454`.

## S5. Rubin combining rules and Barnard-Rubin finite-sample degrees of freedom

**Where we are**

Section 2.4, manuscript lines 85-94 (PDF pages 4-5).

**What the passage is doing**

Across `m` completions, total variance is `T = U_bar + (1 + 1/m) B` and the fraction of missing information is `lambda_hat = (1 + 1/m) diag(B) / diag(T)`. To prevent anti-conservative intervals when sample size `N`, site count `C`, or event count `d` is small, `catmech` computes Barnard-Rubin (1999) degrees of freedom `nu_BR = (nu_m * nu_obs) / (nu_m + nu_obs) < nu_com`.

**In ordinary language**

Old imputation formulas can claim more degrees of freedom than the study has sites or events. The Barnard-Rubin formula caps degrees of freedom below the complete-data limit.

**Evidence status**

ANALYTIC Rubin (1987) and Barnard-Rubin (1999) pooling formulas.

**Do not accidentally conclude**

Do not assume that increasing `m` to 100 lets `nu_BR` exceed `C - 1` in a clustered study; `nu_BR` is strictly bounded below `nu_obs < nu_com = C - 1`.

**Try it**

For `C = 10` clinical sites (`nu_com = 9`), `m = 5`, and `lambda_hat = 0.2308`, what is `nu_BR` (`Lab C10-C11`)? **Answer:** `5.3577` (strictly below `9`, whereas classical `nu_m = 75.11`).

## S6. Multi-endpoint pooling across linear (CV3), GLM (CR1S), and Cox PH (Lin-Wei) models

**Where we are**

Section 2.4, manuscript lines 95-100 (PDF pages 5-6).

**What the passage is doing**

`catmech_pool_lm()` uses the delete-one-cluster jackknife (`cv3`) with `nu_com = C - 1` when clustered; `catmech_pool_glm()` uses the CR1S Huber-White-Liang-Zeger sandwich with `nu_com = C - 1`; and `catmech_pool_coxph()` uses event-based `nu_com = d - p` when unclustered or the Lin-Wei score-residual sandwich with `nu_com = C - 1` when clustered.

**In ordinary language**

Linear models use a leave-one-site-out jackknife variance, logistic and count models use a small-sample sandwich, and survival models count observed events.

**Evidence status**

ANALYTIC multi-endpoint covariance and degrees-of-freedom specifications.

**Do not accidentally conclude**

Do not use total sample size `N - p` for unclustered Cox proportional hazards degrees of freedom; `catmech_pool_coxph()` uses observed failure events `d - p`.

**Try it**

In an unclustered Cox PH model with `N = 600` patients, `d = 45` events, and `p = 3` predictors, what is `nu_com`? **Answer:** `d - p = 45 - 3 = 42`.

## S7. Out-of-sample rank-PIT calibration checking and its slope-bias blind spot

**Where we are**

Section 2.5, manuscript lines 101-106 (PDF page 6).

**What the passage is doing**

`catmech_holdout_check()` runs K-fold cross-validation on audited rows and computes tie-randomized rank PIT values `u_i` against `Uniform(0, 1)`. While false alarms were `0/40` under the true model, alarm power against structural misspecification was low (`0/40` in `nonlinear_target`, `1/40` in `heavy_tail_target`) because held-out rows test observation-level prediction rather than regression slope bias.

**In ordinary language**

Holdout checking hides part of the audit and tests if row predictions match the hidden truth. Passing this screen does not guarantee an unbiased regression slope.

**Evidence status**

ANALYTIC rank-PIT formula and MEASURED_DATED diagnostic power boundary.

**Do not accidentally conclude**

Do not treat `miscalibration_detected = FALSE` in `catmech_holdout_check()` as proof that downstream regression slopes are unbiased.

**Try it**

Why does `catmech_holdout_check()` randomize exact ties (`T_{i,-k}^{(j)} = R_i^*`) using `U_i ~ Uniform(0, 1)`? **Answer:** Because the contamination mixture places a discrete point mass at `W_i` on uncorrupted rows; without tie randomization, exact matches would spike and trigger false KS alarms.

## S8. Synthetic contamination benchmarks (52,950 requests) and three stress-cell failures

**Where we are**

Sections 3-5.1, manuscript lines 107-249 (PDF pages 6-10).

**What the passage is doing**

Across `52,950` synthetic requests, `catmech_contamination()` narrowed CIs by `44.3%` over complete-case analysis in `missed_flags` (`0.1842` vs `0.3309`) and prevented the `-0.2354` linear-MI bias under gross 10x outliers (`error_gross_and_small`, `94.0%` vs `64.7%` coverage). However, it failed in three stress cells: `heavy_tail_target` (`87.3%` coverage), `nonignorable_missed` (`79.3%` coverage), and `noisy_undeclared` (`88.7%` coverage, restored to `91.8%` when `reference_sd` is declared).

**In ordinary language**

Under synthetic gross or unflagged errors, the mixture model stays unbiased and cuts interval width by 44%. It fails under heavy tails, non-ignorable selection, or undeclared reference noise.

**Evidence status**

MODELED_LOCKED synthetic benchmark results across 27 scenarios.

**Do not accidentally conclude**

Do not use `catmech_contamination()` when the true distribution `T_i` has heavy `t_4` tails (`heavy_tail_target`); use `catmech_reference(use_recorded = TRUE)` (`92.0%` coverage) instead.

**Try it**

What caused `catmech_contamination()` to under-cover (`87.3%`) in `heavy_tail_target`? **Answer:** The Gaussian mixture misclassified genuine heavy-tailed extremes in `T_i` as recording errors and shrunk them toward `X_i beta`.

## S9. Clustered variance (148,500 requests): CV3 jackknife and predictor-slope attenuation

**Where we are**

Section 5.2, manuscript lines 250-258 (PDF pages 10-11).

**What the passage is doing**

Across `148,500` clustered requests on simulated cohorts and 11 Rwanda HMIS indicators, `cv3` beat `cr1s` and `cr2` (which inflated widths up to `95%`). However, when the audited variable `T_i` was a **predictor** (`Y ~ T + x + stratum`), clustered random-intercept imputation (`catmech_reference(cluster = )`) attenuated slopes by up to `-16.1%` (`84.4%` coverage), whereas fitting `catmech_reference(use_recorded = TRUE)` without `cluster` and pooling with `catmech_pool_lm(cluster = )` achieved `-3.0%` bias and `95.4%` coverage.

**In ordinary language**

Across clinical sites, the leave-one-site-out jackknife gives the best interval coverage. However, random-intercept imputation shrinks slopes when the audited variable is a predictor.

**Evidence status**

MEASURED_DATED and MODELED_LOCKED clustered variance and multilevel MI evaluation.

**Do not accidentally conclude**

Do not pass `cluster` to `catmech_reference()` when the audited variable is a predictor in the downstream regression (`catmech_pool_lm()` warns if you do).

**Try it**

What is the recommended two-step R call when an audited biomarker is a predictor in a multicenter linear regression? **Answer:** Fit `catmech_reference(..., use_recorded = TRUE)` without `cluster`, then pool with `catmech_pool_lm(..., cluster = 'site')`.

## S10. Public paired-discrepancy benchmarks (33,800 requests) and near-chance row error AUCs

**Where we are**

Section 6.1, manuscript lines 259-279 (PDF pages 11-12).

**What the passage is doing**

On ClinicalTrials.gov (`N = 27,208`) and NHANES 2017-2018 (`N = 3,908`), retaining recorded values reduced RMSE by `43%` to `48%` over discarding unaudited values. However, because real-world discrepancies were small relative to `T_i | X_i`, per-record corruption probabilities `p_i` had near-chance discrimination (`AUC = 0.543` on ClinicalTrials.gov and `0.500` on NHANES, versus `0.962` in synthetic `missed_flags`).

**In ordinary language**

On real data with diffuse errors, keeping recorded numbers cuts slope error nearly in half, but per-row error probabilities score no better than a coin flip.

**Evidence status**

MEASURED_DATED empirical benchmark across ClinicalTrials.gov and NHANES.

**Do not accidentally conclude**

Never use individual record-level corruption probabilities `p_i` from `catmech_contamination()` to adjudicate or discard specific patient records in diffuse real-world data.

**Try it**

What was the ROC-AUC of per-record corruption probabilities `p_i` at a 25% audit on ClinicalTrials.gov and NHANES? **Answer:** `0.543` on ClinicalTrials.gov and `0.500` on NHANES (near coin-flip chance `0.500`).

## S11. Real-audit slope recovery (45,600 requests) and the methodological decision tree

**Where we are**

Sections 6.2-8, manuscript lines 280-328 (PDF pages 12-14).

**What the passage is doing**

On Rwanda HMIS facility recounts (`45,600` requests), Rule 2 comparing `catmech_contamination()` against `Linear MI` was `INCONCLUSIVE`: the mixture had lower point RMSE on `PNC` (`0.053` vs `0.072`) and `Iron` (`0.112` vs `0.156`), but failed 95% coverage on all three primary slopes (`71.6%` on `PNC`, `66.6%` on `Syphilis`, `35.4%` on `Iron`) and failed the `>= 90%` return rate on two. Therefore, Section 7 recommends `catmech_contamination()` only for gross, separable outliers and `catmech_reference(use_recorded = TRUE)` for moderate or diffuse real-world discrepancies.

**In ordinary language**

On Rwanda health-center recounts, the mixture model had lower point error on two slopes but failed interval coverage. For everyday moderate errors, use simple linear imputation.

**Evidence status**

MEASURED_DATED real-audit evaluation (`INCONCLUSIVE` Rule 2) and INTERPRETIVE decision tree.

**Do not accidentally conclude**

Do not claim that `catmech_contamination()` outperformed Simple Linear MI on real-world audit intervals; its interval advantage remains synthetic-only.

**Try it**

Under the Section 7 decision tree, what two checks must pass before trusting `catmech_contamination()` on a dataset with suspected gross errors? **Answer:** A `>= 90%` interval return rate and positive held-out Brier/RMSE skill over raw recorded values in `catmech_holdout_check()`.

## S12. Provenance reporting (catmech_audit_trail), software QA, and regulatory boundaries

**Where we are**

Sections 9-12, manuscript lines 329-350 (PDF pages 14-15).

**What the passage is doing**

`catmech_audit_trail()` outputs a provenance record (Markdown, JSON, or HTML) recording R/OS metadata, SHA-256 hashes of the input data, assumption specification, and all `m` completed datasets, MCMC convergence status, and pooled tables. This supports documentation under FDA 21 CFR Part 11 and ICH E9 without certifying compliance. The package is MIT-licensed and verified by `1,754` unit tests (`0` failures/warnings/skips), `R CMD check --as-cran` (`0` errors/warnings), and `500,000+` replicate checks.

**In ordinary language**

The audit trail function writes a report with SHA-256 data hashes, model settings, and convergence status. It supports documentation but does not certify legal compliance.

**Evidence status**

DEFINITION, SCOPE, and MEASURED_DATED software verification summary.

**Do not accidentally conclude**

Do not represent `catmech_audit_trail()` as a legal compliance certification; regulatory compliance depends on the surrounding validated system, access controls, and procedures.

**Try it**

Run `Lab C28` in `education/run_companion_labs.R`. What default title and footer heading does `catmech_audit_trail()` print? **Answer:** Default title `Catmech Audit Trail & Provenance Report` and footer heading `### Audit Trail & Provenance Record`.
