---
title: 'catmech: An Auditable R Package for Measurement Verification, Clinical Data Audits, and Multi-Endpoint Model Pooling'
tags:
  - R
  - biostatistics
  - clinical trials
  - measurement error
  - multiple imputation
  - audit planning
  - survival analysis
  - data integrity
  - regulatory compliance
authors:
  - name: Patrick Iben
    orcid: 0009-0005-2443-7973
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 6 November 2026
bibliography: paper.bib
---

# Summary

In multicenter clinical trials, electronic health record (EHR) registries, and observational epidemiological cohorts, continuous biological and clinical measurements (such as circulating cardiac troponin, viral loads, blood pressure, serum creatinine, or tumor volume) are frequently compromised by measurement errors, transcription slips, assay lot shifts, and missing values. Traditional statistical packages handle missing data by assuming values are Missing at Random (MAR) [@rubin1987multiple; @little2019statistical] and applying standard multiple imputation algorithms [@vanbuuren2011mice; @honaker2011amelia; @white2011multiple]. However, in clinical data auditing and source-data verification (SDV), the central challenge is fundamentally different: recorded values are often present but untrusted, and an internal gold-standard reference measurement has been audited on only a fraction of records [@angelopoulos2023prediction; @robins1994estimation; @breslow2003semiparametric].

`catmech` is an R package providing a mathematically rigorous, auditable framework for clinical data verification, two-phase audit sample size determination, predictive calibration checking, downstream multi-endpoint regression pooling, and regulatory audit trail generation. By explicitly modeling the contamination mixture process rather than discarding un-audited records or naively trusting recorded numbers, `catmech` retains the precision of verified recordings while correcting for contamination bias. Furthermore, the package integrates cryptographic provenance tracking (SHA-256 data and completion hashes) designed to meet strict regulatory compliance standards (FDA 21 CFR Part 11 [@fda2003part11] and ICH GCP E9 [@ich1998e9]).

# Statement of need

When an internal audit reveals recording discrepancies in an observational or clinical trial cohort, researchers face three unsatisfactory conventional choices:

1. **Naively trusting the recorded value**: Contaminated measurements act as classic errors-in-variables, attenuating regression coefficients, obscuring exposure-response gradients, and inflating residual variance [@carroll2006measurement].
2. **Complete-case analysis (audited cases only)**: Restricting downstream analysis exclusively to audited records drastically reduces sample size, severely inflating standard errors and sacrificing statistical power [@neyman1938contribution].
3. **Discarding recorded values and applying standard multiple imputation**: Treating all un-audited recorded entries as completely missing wastes the substantial information contained in records that were recorded accurately, leading to needlessly wide confidence intervals [@rubin1996multiple].

`catmech` resolves this dilemma by formalizing the two-phase sampling relationship [@cochran1977sampling] between observed recorded values $W_i$, latent true values $T_i$, and audited reference measurements $R_i$. It addresses both pre-audit planning (calculating how many medical records must be audited to achieve a target precision ratio) and post-audit analysis (pooling regression estimates with calibrated finite-sample degrees of freedom).

Crucially, `catmech` refuses to conceal model failures. If MCMC chains fail Gelman-Rubin convergence screens [@gelman1992inference], or if cross-validated holdout checks detect predictive miscalibration on audited reference cases [@dawid1984present], downstream inference is visibly refused.

# Comparison with existing software

While several mature R packages address missing data or measurement error, none integrate prospective two-phase audit planning, reference-calibrated contamination mixture modeling, multi-endpoint pooling with small-sample degrees of freedom, and cryptographic regulatory compliance tracking:

| Feature / Capability | `catmech` | `mice` [@vanbuuren2011mice] | `simex` [@leeder2008simex] | `mecor` [@nab2021mecor] | `Amelia` [@honaker2011amelia] |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Contamination Mixture Modeling** | **Yes** (Gibbs) | No | No | No | No |
| **Partial Gold-Standard Audits** | **Yes** (Two-phase) | Partial (MAR) | No | Yes (Validation) | Partial (MAR) |
| **Prospective Audit Power Planning** | **Yes** (Closed + MC) | No | No | No | No |
| **Predictive Calibration Diagnostics** | **Yes** (Holdout PIT) | No | No | No | No |
| **Linear Model Pooling** | **Yes** (BR df + CV3) | Yes (Rubin df) | No | Yes | Yes (Rubin df) |
| **GLM (Logistic/Poisson) Pooling** | **Yes** (Robust Sandwich) | Yes (Wald/Rubin) | Yes (GLM) | Yes | Yes |
| **Cox PH Survival Pooling** | **Yes** (Lin-Wei Sandwich) | Partial (via coxph) | Partial | No | No |
| **SHA-256 Provenance & Audit Trail** | **Yes** (21 CFR Part 11) | No | No | No | No |

# Mathematical framework and methodology

## 1. Reference-calibrated contamination mixture model

Let $T_i$ represent the latent true measurement for patient $i \in \{1, \dots, N\}$, conditionally distributed given auxiliary predictors $X_i$ as:
$$T_i \mid X_i \sim \mathcal{N}(X_i \beta, \sigma_T^2)$$

Let $\delta_i \in \{0, 1\}$ denote an unobserved contamination indicator with patient-specific corruption probability $\pi_i = \mathbb{P}(\delta_i = 1 \mid S_i)$, which may vary across clinical sites or quality flags $S_i$. The recorded measurement $W_i$ is generated as:
$$W_i = (1 - \delta_i) T_i + \delta_i (T_i + \epsilon_i), \quad \epsilon_i \sim \mathcal{N}(\mu_e, \sigma_e^2)$$
where $\epsilon_i$ represents an additive recording shift. On the audited subset $\mathcal{A} \subset \{1, \dots, N\}$, an independent reference measurement $R_i$ is observed:
$$R_i = T_i + \eta_i, \quad \eta_i \sim \mathcal{N}(0, \sigma_{\text{ref}}^2)$$
where $\sigma_{\text{ref}}^2 = 0$ corresponds to an exact gold-standard recount.

The model is estimated via Markov Chain Monte Carlo (MCMC) using a conjugate Gibbs sampler implemented in C++ via `RcppArmadillo`. For un-audited patients ($i \notin \mathcal{A}$), the posterior probability that a recorded value is corrupted is evaluated analytically at each Gibbs step:
$$p_i = \mathbb{P}(\delta_i = 1 \mid W_i, X_i, \Theta) = \frac{\pi f_{\epsilon}(W_i - \hat{T}_i)}{\pi f_{\epsilon}(W_i - \hat{T}_i) + (1 - \pi) f_{T}(W_i \mid X_i)}$$
Imputations are drawn from the exact conditional posterior, preserving uncorrupted recorded values while sampling fresh draws for corrupted entries.

## 2. Prospective power and audit sample size planning

Before committing clinical resources to chart re-abstraction, trialists must determine the necessary audit sample size $n_{\text{audit}}$. Using influence function expansions for efficient two-phase estimators [@neyman1938contribution; @cochran1977sampling; @breslow2003semiparametric], the ratio of the confidence interval width at audit fraction $f = n_{\text{audit}} / N$ relative to an ideal full audit ($f = 1$) is:
$$\rho(f) = \frac{\text{Width}(f)}{\text{Width}(1)} = \sqrt{1 + k \frac{1 - f}{f}}, \quad k = \frac{\text{Var}(\psi_i - \mathbb{E}[\psi_i \mid X_i])}{\text{Var}(\psi_i)}$$
where $\psi_i$ is the influence curve of the parameter of interest.

Because $k$ is subject to pilot sampling variability, the realized precision ratio $\hat{\rho}(f)$ is stochastic. `catmech_audit_power()` computes the exact statistical power to achieve a target precision ratio $\rho \le \rho_{\text{target}}$ (e.g., $\rho_{\text{target}} = 1.25$):
$$\text{Power}(n_{\text{audit}}) = \mathbb{P}\left(\sqrt{1 + k \frac{1 - f}{f}} \le \rho_{\text{target}}\right)$$
The function supports both prospective calculation from prior noise parameters and empirical calculation from a pilot audit plan (`catmech_audit_plan()`), generating visual power curves via `plot()`.

## 3. Downstream multi-endpoint model pooling

Parameter estimates $\hat{\theta}_j$ and variance-covariance matrices $\hat{U}_j$ across $j \in \{1, \dots, m\}$ completed datasets are combined using Rubin's rules [@rubin1987multiple]:
$$\bar{\theta} = \frac{1}{m} \sum_{j=1}^m \hat{\theta}_j, \quad \bar{U} = \frac{1}{m} \sum_{j=1}^m \hat{U}_j, \quad B = \frac{1}{m - 1} \sum_{j=1}^m (\hat{\theta}_j - \bar{\theta})(\hat{\theta}_j - \bar{\theta})^\top$$
$$T = \bar{U} + \left(1 + \frac{1}{m}\right) B, \quad \hat{\lambda} = \frac{(1 + 1/m) \text{diag}(B)}{\text{diag}(T)}$$

To prevent undercoverage in finite clinical samples, `catmech` implements the Barnard & Rubin (1999) degrees of freedom adjustment [@barnard1999small]:
$$\nu_m = (m - 1) \hat{\lambda}^{-2}, \quad \nu_{obs} = \left(\frac{\nu_{com} + 1}{\nu_{com} + 3}\right) \nu_{com} (1 - \hat{\lambda}), \quad \nu_{BR} = \frac{\nu_m \nu_{obs}}{\nu_m + \nu_{obs}}$$
where $\nu_{com}$ is the complete-data degrees of freedom. Downstream model pooling covers three core clinical endpoint families:

* **Linear models (`catmech_pool_lm`)**: For continuous endpoints, with default delete-one-cluster jackknife (CV3) variance estimation.
* **Generalized linear models (`catmech_pool_glm`)**: For binary (logistic) and count (Poisson) endpoints, computing exponentiated odds ratios and incidence rate ratios with Huber-White-Liang-Zeger cluster-robust sandwich covariance [@liang1986longitudinal]:
  $$\hat{V}_{\text{sandwich}} = \hat{J}^{-1} \left(\sum_{c=1}^C s_c s_c^\top\right) \hat{J}^{-1}$$
* **Cox proportional hazards models (`catmech_pool_coxph`)**: For time-to-event survival outcomes, defining $\nu_{com} = d - p$ (where $d$ is the number of observed failure events and $p$ is the parameter count) and providing Lin & Wei (1989) cluster-robust sandwich standard errors [@lin1989robust].

## 4. Diagnostic screens and visualizations

`catmech` provides unified S3 `plot()` diagnostics:

* **Holdout cross-validation (`catmech_holdout_check()`)**: Partitions audited records into $K$ cross-validation folds and computes the randomized Probability Integral Transform (PIT) [@dawid1984present]:
  $$u_i = F(R_i \mid \mathcal{D}_{-k}) \sim \text{Uniform}(0, 1)$$
  A two-panel diagnostic displays the PIT histogram with a Kolmogorov-Smirnov uniformity test alongside reference vs. predictive calibration scatter.
* **Recording disagreement audit (`catmech_recording_audit()`)**: Evaluates stratum-specific error rates with exact 95% Jeffreys beta intervals ($\text{Beta}(k + 0.5, n - k + 0.5)$) against the global error rate.
* **Multi-scenario sensitivity envelopes (`catmech_sensitivity()`)**: Evaluates the stability of clinical conclusions across competing, non-identifiable assumptions (e.g., contamination mixture vs. pure reference modeling vs. selection shifts).

# Illustrative clinical audit workflow

The following complete workflow illustrates the `catmech` pipeline on a simulated multicenter clinical cohort:

```r
library(catmech)
library(survival)

# 1. Simulate a multicenter cohort with recording contamination and 25% reference audit
set.seed(42)
cohort <- simulate_catmech_lab(n = 1000, seed = 42)
cohort$site <- rep(seq_len(10), length.out = nrow(cohort))
audited <- sample.int(nrow(cohort), size = 250)
cohort$lab_ref <- ifelse(seq_len(nrow(cohort)) %in% audited, cohort$lab_true, NA_real_)
cohort$adverse_event <- as.integer(cohort$outcome > median(cohort$outcome))
cohort$time <- rexp(nrow(cohort), rate = exp(0.25 * (cohort$lab_true - 6)))
cohort$event <- rbinom(nrow(cohort), 1, 0.75)

# 2. Prospective audit sample size and statistical power planning
power_plan <- catmech_audit_power(
  N = 1000,
  target_ratio = 1.25,
  power_target = 0.80,
  error_rate = 0.20,
  correlation = 0.85,
  fractions = seq(0.10, 0.50, by = 0.05)
)
plot(power_plan)

# 3. Discrepancy screening across clinical sites
audit_screen <- catmech_recording_audit(
  data = cohort,
  value = "lab_value",
  reference = "lab_ref",
  by = "site",
  agreement_tol = 0.5
)
plot(audit_screen)

# 4. Reference-calibrated Gibbs contamination model
fit <- catmech_contamination(
  data = cohort,
  value = "lab_value",
  reference = "lab_ref",
  predictors = c("age", "sex", "inflammation", "outcome"),
  corruption_strata = "outlier_flag",
  m = 10,
  seed = 2026
)

# 5. Out-of-sample predictive calibration verification
calib_check <- catmech_holdout_check(
  fit = fit,
  folds = 5L,
  m = 20L,
  seed = 123L
)
plot(calib_check)

# 6. Multi-endpoint downstream regression pooling and forest plots
# Continuous linear model
pool_continuous <- catmech_pool_lm(
  result = fit,
  formula = outcome ~ lab_value + age + sex,
  cluster = "site"
)

# Binary logistic GLM (Odds Ratios)
pool_binary <- catmech_pool_glm(
  result = fit,
  formula = adverse_event ~ lab_value + age + sex,
  family = binomial(link = "logit"),
  cluster = "site",
  exponentiate = TRUE
)
plot(pool_binary, exponentiate = TRUE)

# Survival Cox proportional hazards (Hazard Ratios)
pool_survival <- catmech_pool_coxph(
  result = fit,
  formula = Surv(time, event) ~ lab_value + age + sex,
  cluster = "site",
  exponentiate = TRUE
)
plot(pool_survival, exponentiate = TRUE)

# 7. Side-by-side comparison of all four audit strategies
comp <- catmech_compare_strategies(
  data = cohort,
  value = "lab_value",
  reference = "lab_ref",
  formula = outcome ~ lab_value + age + sex,
  m = 10,
  seed = 2026
)
plot(comp, term = "lab_value")

# 8. Regulatory compliance audit report (FDA 21 CFR Part 11)
report <- catmech_audit_trail(
  object = fit,
  pooled = pool_survival,
  format = "markdown",
  path = "audit_trail.md"
)
```

The resulting pooled summary for the Cox proportional hazards model yields an auditable table reporting log-hazard point estimates, Lin-Wei cluster-robust standard errors, Barnard-Rubin finite-sample degrees of freedom, exponentiated hazard ratios (`estimate_exp`), and 95% confidence intervals:

```text
-- Catmech Pooled Cox PH Model --
Completions (m): 10 | Observations (n): 1000 | Events: 735

      term estimate std.error     df   lower  upper estimate_exp lower_exp upper_exp
 lab_value   0.1632    0.0314 6.5933  0.0881 0.2384       1.1773    1.0921    1.2692
       age   0.0056    0.0041 7.4010 -0.0038 0.0151       1.0056    0.9962    1.0152
      sexM  -0.0541    0.0433 7.3074 -0.1555 0.0474       0.9474    0.8560    1.0485
```

# Empirical validation and simulation benchmarks

## 1. Synthetic Monte Carlo studies

To evaluate the statistical properties of `catmech`, pre-specified Monte Carlo simulation studies encompassing over 500,000 requests across 10 validation archives were conducted. In the primary contamination study (`validation/contamination-study`, 52,950 planned requests: 52,900 returned intervals, 16 expected design refusals, 33 convergence refusals, and 1 singular raking fit), seven analytical strategies were compared for estimating a true slope $\beta = 0.800$ in cohorts of size $N = 600$ with a 25% reference audit ($n_{\text{audit}} = 150$) and $m = 20$ imputations:

| Scenario & Strategy | Est. $\hat{\beta}$ | Bias (MCSE) | Emp. SD | Mean SE | Coverage | Mean Width | RMSE |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Unflagged Errors (`missed_flags`, $N_{\text{sim}} = 1,000$)** | | | | | | | |
| Oracle (Latent $T$ Known) | 0.7978 | -0.0022 (0.0012) | 0.0394 | 0.0411 | 95.7% | 0.1614 | 0.0394 |
| Flag Masking (`v02_mask_normal`) | 0.7185 | -0.0815 (0.0015) | 0.0472 | 0.0444 | 55.9% | 0.1750 | 0.0942 |
| Complete-Case (25% Audited Only) | 0.7988 | -0.0012 (0.0026) | 0.0815 | 0.0837 | 95.7% | 0.3309 | 0.0814 |
| Discard Recorded MI (`reference_normal`) | 0.7896 | -0.0104 (0.0021) | 0.0676 | 0.0684 | 95.4% | 0.2769 | 0.0684 |
| Simple Linear MI (`cole_mi`) | 0.7806 | -0.0194 (0.0018) | 0.0557 | 0.0528 | 92.4% | 0.2101 | 0.0590 |
| Generalized Raking (`raking`) | 0.7853 | -0.0147 (0.0020) | 0.0631 | 0.0589 | 92.9% | 0.2331 | 0.0648 |
| **`catmech` Contamination Mixture** | **0.7952** | **-0.0048 (0.0014)** | **0.0446** | **0.0467** | **94.9%** | **0.1842** | **0.0449** |
| **Gross & Small Errors (`error_gross_and_small`, $N_{\text{sim}} = 150$)** | | | | | | | |
| Oracle (Latent $T$ Known) | 0.7970 | -0.0030 (0.0034) | 0.0411 | 0.0412 | 94.7% | 0.1617 | 0.0411 |
| Flag Masking (`v02_mask_normal`) | 0.0609 | -0.7391 (0.0055) | 0.0674 | 0.0155 | 0.7% | 0.0614 | 0.7422 |
| Simple Linear MI (`cole_mi`) | 0.5646 | -0.2354 (0.0261) | 0.3200 | 0.0617 | 64.7% | 0.2562 | 0.3964 |
| Generalized Raking (`raking`) | 0.5156 | -0.2844 (0.0281) | 0.3436 | 0.0693 | 53.7% | 0.2740 | 0.4451 |
| Complete-Case (25% Audited Only) | 0.7962 | -0.0038 (0.0064) | 0.0787 | 0.0845 | 95.3% | 0.3339 | 0.0785 |
| Discard Recorded MI (`reference_normal`) | 0.7890 | -0.0110 (0.0055) | 0.0671 | 0.0695 | 95.3% | 0.2816 | 0.0678 |
| **`catmech` Contamination Mixture** | **0.7883** | **-0.0117 (0.0035)** | **0.0432** | **0.0436** | **94.0%** | **0.1714** | **0.0447** |

These findings illustrate the primary advantage of `catmech` in settings with unflagged or gross recording errors:

* Trusting unflagged recordings (`v02_mask_normal`) suffers severe attenuation bias ($-0.0815$ in `missed_flags` and $-0.7391$ in `error_gross_and_small`), collapsing 95% coverage to 55.9% and 0.7%.
* Complete-case analysis (`gold_only`) and discarding recorded values (`reference_normal`) eliminate bias but forfeit precision (mean CI widths of 0.3309 and 0.2769 in `missed_flags`).
* Simple linear MI (`cole_mi`) and generalized raking (`raking`) break down when gross errors are mixed with small errors (`error_gross_and_small`), suffering severe bias ($-0.2354$ and $-0.2844$) and undercoverage (64.7% and 53.7%).
* The `catmech` contamination mixture model preserves nominal coverage (94.9% in `missed_flags` and 94.0% in `error_gross_and_small`) while narrowing confidence intervals by 33.5% relative to discarding recorded values and by 44.3% relative to complete-case analysis (`0.1842` vs. `0.3309`). Across all four primary 1,000-cohort cells (`clean`, `flagged_errors`, `missed_flags`, `site_shift_declared`), `catmech_contamination()` passed the $[0.925, 0.975]$ Wilson coverage screen in three cells (96.5%, 96.3%, and 94.9%) and narrowly missed in `site_shift_declared` at 93.7% (Wilson $[0.9202, 0.9505]$, bias $-0.0085$, MCSE $0.0016$), while narrowing intervals by 30% to 38% relative to `reference_normal`.

## 2. Real-world empirical benchmarks

The package was further evaluated on two large-scale public datasets containing real recording differences without synthetic corruption (33,800 requests executed across 1,000 resampled cohorts of size $N = 600$ at the 25% audit fraction and 200 cohorts at secondary fractions):

| Dataset & Metric | Audit 25% Sample | Discard Recorded | Contamination Model | Simple Linear MI | Full Reference (Oracle) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **ClinicalTrials.gov** ($N = 27,208$) | | | | | |
| Enrolled vs. Started discrepancy rate: 14.2% | Bias | -0.0003 | +0.0068 | +0.0001 | -0.0003 |
| Target slope: 0.0357 log months / participant | Coverage | 92.3% | 91.2% | 91.8% | 91.2% |
| | Mean Width | 0.205 | 0.106 | 0.102 | 0.096 |
| | RMSE | 0.0539 | 0.0308 | 0.0286 | 0.0275 |
| **NHANES 2017--2018** ($N = 3,908$) | | | | | |
| Self-reported vs. measured weight discrepancy: 67.3% | Bias | -0.0063 | -0.0066 | -0.0025 | +0.0001 |
| Target slope: 0.127 mmHg / kg | Coverage | 94.2% | 94.5% | 93.9% | 94.9% |
| | Mean Width | 0.297 | 0.155 | 0.154 | 0.140 |
| | RMSE | 0.0785 | 0.0408 | 0.0409 | 0.0357 |

On both real datasets, retaining the recorded values cut RMSE by 43% to 48% compared to discarding un-audited records. However, when recording errors are diffuse rather than gross and separable, simple linear imputation matched or slightly outperformed the mixture model (RMSE 0.0286 vs. 0.0308 on ClinicalTrials.gov) with fewer parametric assumptions.

In addition, facility-level data verification was benchmarked against the Rwanda Health Management Information System recount audit [@nshimyiryo2020health] across 11 clinical indicators and 148,500 requests (148,486 returned intervals and 14 rank-deficiency refusals), confirming that delete-one-cluster jackknife (CV3) variance estimation achieved nominal 92.0% to 97.5% recount coverage across clinical facilities where uncorrected cluster variances suffered undercoverage.

## 3. Simulation-to-practice divergence and methodological guidance

A central methodological insight established by these empirical benchmarks is the divergence between synthetic Monte Carlo evaluations and real-world clinical audits:

* **Synthetic gross errors favor mixture modeling**: When data errors represent discrete, bimodal corruption events (e.g., transcription slips, decimal shifts, or assay calibration jumps), the Bayesian contamination mixture model is essential to eliminate severe attenuation bias (bias of $-0.0117$ vs. $-0.2354$ for linear MI in `error_gross_and_small`) and narrow confidence intervals by 44.3% over complete-case analysis.
* **Real-world subtle errors favor linear multiple imputation**: In actual registry audits where discrepancies are diffuse and continuous (e.g., self-report rounding or mild procedural differences), simple linear imputation matches or exceeds mixture modeling in predictive accuracy (RMSE 0.0286 vs. 0.0308 on ClinicalTrials.gov) with fewer parametric assumptions and zero convergence overhead.

This empirical divergence highlights why `catmech` incorporates prospective audit planning (`catmech_audit_power()`), diagnostic discrepancy auditing (`catmech_recording_audit()`), predictive cross-validation (`catmech_holdout_check()`), and automated strategy comparison (`catmech_compare_strategies()`), ensuring that trialists let empirical data guide the choice of model rather than relying on unverified assumptions.

# Regulatory compliance and data governance

Regulatory submissions in biopharmaceutical and medical device clinical trials (under FDA 21 CFR Part 11 [@fda2003part11] and ICH GCP E9 [@ich1998e9]) require an unbroken, tamper-evident audit trail for any algorithmic data imputation or measurement correction. `catmech_audit_trail()` generates complete provenance records in Markdown, JSON, or standalone HTML. The report documents:

* Full execution environment metadata (R version, operating system platform, package version, UTC execution timestamp).
* Cryptographic SHA-256 checksums of the raw input dataset, the mathematical assumption specification, and each of the $m$ completed datasets.
* Model convergence statistics (MCMC iterations, Gelman-Rubin $\hat{R}$, bulk and tail effective sample sizes).
* Downstream pooled parameter tables and multi-scenario sensitivity envelopes.

# Quality assurance and verification

`catmech` has been developed under an automated software verification harness:

* **Unit tests**: 1,754 unit tests executed via `testthat`, achieving 0 failures, 0 warnings, and 0 skips.
* **CRAN conformance**: Fully compliant with CRAN package policies, passing `R CMD check --no-manual` with `Status: OK` (0 errors, 0 warnings, 0 notes).
* **Independent artifact verification**: Verified against an independent Python verification suite (`verify-artifacts.py`, `verify-release-archive.py`, and `verify_paper_numbers.py`) confirming exact correspondence across 10 frozen simulation and empirical benchmark studies and over 500,000 replicate requests.

# Availability and dependencies

`catmech` is open-source software licensed under the MIT license. It requires R ($\ge$ 4.1.0) and a C++17 compiler, utilizing `Rcpp`, `RcppArmadillo`, `digest`, `jsonlite`, `mice`, `pan`, and `survival`. Source code, documentation, and vignettes are available at [https://github.com/patrickiben/catmech](https://github.com/patrickiben/catmech).

# References

<!-- gate-basis: paper/paper.bib@e4efb6523c93e3ea4ab40847542c68a6cc6b0e8dce0d33ac370f3d2e1b577cad -->
