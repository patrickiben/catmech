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
  - regulatory documentation
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

`catmech` is an R package providing a mathematically rigorous, auditable framework for clinical data verification, two-phase audit sample size determination, predictive calibration checking, downstream multi-endpoint regression pooling, and provenance audit trail generation. By modeling either a reference-calibrated contamination mixture (`catmech_contamination()`) or a linear reference calibration model (`catmech_reference(use_recorded = TRUE)`) rather than discarding unaudited records or naively trusting recorded numbers, `catmech` retains the precision of verified recordings while correcting for contamination bias. The package also records SHA-256 hashes of the input data, assumptions, and completed datasets to support regulatory documentation (FDA 21 CFR Part 11 [@fda2003part11] and ICH E9 [@ich1998e9]); this supports, but does not by itself establish, compliance.

# Statement of need

When an internal audit reveals recording discrepancies in an observational or clinical trial cohort, researchers face three unsatisfactory conventional choices:

1. **Naively trusting the recorded value**: Contaminated measurements act as classic errors-in-variables, attenuating regression coefficients, obscuring exposure-response gradients, and inflating residual variance [@carroll2006measurement].
2. **Complete-case analysis (audited cases only)**: Restricting downstream analysis exclusively to audited records drastically reduces sample size, severely inflating standard errors and sacrificing statistical power [@neyman1938contribution].
3. **Discarding recorded values and applying standard multiple imputation**: Treating all unaudited recorded entries as completely missing wastes the substantial information contained in records that were recorded accurately, leading to needlessly wide confidence intervals [@rubin1996multiple].

`catmech` resolves this dilemma by formalizing the two-phase sampling relationship [@cochran1977sampling] between observed recorded values $W_i$, latent true values $T_i$, and audited reference measurements $R_i$. It addresses both pre-audit planning (calculating how many medical records must be audited to achieve a target precision ratio) and post-audit analysis (pooling regression estimates with calibrated finite-sample degrees of freedom).

In addition, `catmech` refuses to conceal model failures. If MCMC chains fail Gelman-Rubin convergence screens [@gelman1992inference] (`convergence_pass = FALSE`), downstream pooling via `catmech_pool_lm()`, `catmech_pool_glm()`, and `catmech_pool_coxph()` is visibly refused, while `catmech_holdout_check()` provides a standalone cross-validated predictive calibration screen on audited reference cases [@dawid1984present].

# Comparison with existing software

While several mature R packages address missing data or measurement error, none integrate prospective two-phase audit power planning, reference-calibrated contamination mixture modeling, multi-endpoint pooling with cluster-robust covariance and Barnard-Rubin small-sample degrees of freedom, and cryptographic provenance tracking:

| Capability | `catmech` | `mice` [@vanbuuren2011mice] | `simex` [@leeder2008simex] | `mecor` [@nab2021mecor] | `Amelia` [@honaker2011amelia] |
| :------------------------ | :-----------------: | :--------------: | :--------------: | :--------------: | :----------------: |
| **Contamination Mixture** | **Yes** (Gibbs) | No | No | No | No |
| **Gold-Standard Audits** | **Yes** (Two-phase) | Partial (MAR) | No (Replicates) | Yes (Reg. Cal.) | Partial (Priors) |
| **Audit Power Planning** | **Yes** (MC) | No | No | No | No |
| **Calibration Screen** | **Yes** (Rank PIT) | No | No | No | Partial (`overimpute`) |
| **Linear Model Pooling** | **Yes** (BR + CV3) | Yes (BR df) | Yes (SIMEX) | Yes (Reg. Cal.) | Yes (Rubin df) |
| **GLM (Logit/Poisson)** | **Yes** (Sandwich) | Yes (BR df) | Yes (SIMEX) | No | Yes (Rubin df) |
| **Cox PH Pooling** | **Yes** (Lin-Wei) | Partial (`coxph`) | No | No | No |
| **SHA-256 Audit Trail** | **Yes** | No | No | No | No |

# Mathematical framework and methodology

## 1. Reference-calibrated contamination mixture model and linear reference MI

Let $T_i$ represent the latent true measurement for patient $i \in \{1, \dots, N\}$ across $C$ clinical sites, conditionally distributed given auxiliary predictors $X_i$ as:
$$T_i \mid X_i \sim \mathcal{N}(X_i \beta, \sigma_T^2)$$

Let $\delta_i \in \{0, 1\}$ denote an unobserved contamination indicator with stratum-specific corruption probability $\pi_{s_i} = \mathbb{P}(\delta_i = 1 \mid S_i = s_i)$ (`corruption_strata`). Writing $(\pi, \mu_e, \sigma_e^2)$ for a single stratum, the recorded measurement $W_i$ is generated as:
$$W_i = (1 - \delta_i) T_i + \delta_i (T_i + \epsilon_i), \quad \epsilon_i \sim \mathcal{N}(\mu_e, \sigma_e^2)$$
where $\epsilon_i$ represents an additive recording shift (`error_strata`). On the audited subset $\mathcal{A} \subset \{1, \dots, N\}$ of size $n_{\text{audit}} = |\mathcal{A}|$, an independent reference measurement $R_i$ is observed:
$$R_i = T_i + \eta_i, \quad \eta_i \sim \mathcal{N}(0, \sigma_{\text{ref}}^2)$$
where $\sigma_{\text{ref}}$ (`reference_sd`, defaulting to $0$) is a user-declared sensitivity constant, with $\sigma_{\text{ref}}^2 = 0$ corresponding to an exact gold-standard recount.

In `catmech_contamination()`, the parameter vector $\Theta = (\beta, \sigma_T^2, \pi, \mu_e, \sigma_e^2)$ is estimated in R via a collapsed slice-within-Gibbs sampler under a noninformative Jeffreys prior $p(\beta, \sigma_T^2) \propto 1/\sigma_T^2$ (for $\sigma_{\text{ref}} = 0$), a Beta prior $\pi \sim \text{Beta}(a_\pi, b_\pi)$, and a normal-scaled-inverse-$\chi^2$ prior on $(\mu_e, \sigma_e^2)$. Integrating out latent $(\delta_i, T_i)$ on unaudited rows ($i \notin \mathcal{A}$), $\text{logit}(\pi)$, $\mu_e$, and $\frac{1}{2}\log(\sigma_e^2)$ are updated via stepping-out slice sampling (plus a scale-ridge slice step when fewer than 5 discrepancies are observed), and the conditional posterior probability that an unaudited recorded value is corrupted is evaluated analytically from the convolved predictive density $W_i \mid \delta_i = 1, X_i, \Theta \sim \mathcal{N}(X_i \beta + \mu_e, \, \sigma_T^2 + \sigma_e^2)$:
$$p_i = \mathbb{P}(\delta_i = 1 \mid W_i, X_i, \Theta) = \frac{\pi \phi(W_i; \, X_i \beta + \mu_e, \, \sigma_T^2 + \sigma_e^2)}{\pi \phi(W_i; \, X_i \beta + \mu_e, \, \sigma_T^2 + \sigma_e^2) + (1 - \pi) \phi(W_i; \, X_i \beta, \, \sigma_T^2)}$$
where $\phi(w; \, m, v)$ is the Gaussian density with mean $m$ and variance $v$. Conditional on $\delta_i = 0$, $T_i = W_i$ is retained deterministically; conditional on $\delta_i = 1$, $T_i$ is drawn from the conjugate Gaussian posterior combining $X_i \beta$ and $W_i - \mu_e$.

For diffuse discrepancies where a discrete mixture is unnecessary, `catmech_reference()` provides Gaussian reference multiple imputation accelerated in C++17 via `RcppArmadillo`. Setting `use_recorded = TRUE` (Simple Linear MI, `cole_mi`) includes $W_i$ as a continuous linear predictor alongside $X_i$ without a mixture indicator $\delta_i$, whereas `use_recorded = FALSE` (Discard Recorded MI, `reference_normal`) discards $W_i$ on unaudited rows.

## 2. Prospective power and audit sample size planning

Before committing clinical resources to chart re-abstraction, trialists must determine the necessary audit sample size $n_{\text{audit}}$. Using influence function expansions for efficient two-phase estimators [@neyman1938contribution; @cochran1977sampling; @breslow2003semiparametric], the ratio of the confidence interval width at audit fraction $f = n_{\text{audit}} / N$ relative to an ideal full audit ($f = 1$) is:
$$\rho(f) = \frac{\text{Width}(f)}{\text{Width}(1)} = \sqrt{1 + k \frac{1 - f}{f}}, \quad k = \frac{\text{Var}(\psi_i - \mathbb{E}[\psi_i \mid G_i])}{\text{Var}(\psi_i)}$$
where $\psi_i$ is the influence curve of the parameter of interest and $G_i = (W_i, Y_i, X_i)$ denotes the Phase I auxiliary variables.

Because the pilot noise estimator $\hat{k}$ is subject to sampling variability, the realized precision ratio $\hat{\rho}(f)$ is stochastic. `catmech_audit_power()` computes the statistical power to achieve a target precision ratio $\rho \le \rho_{\text{target}}$ (e.g., $\rho_{\text{target}} = 1.25$):
$$\text{Power}(n_{\text{audit}}) = \mathbb{P}\left(\sqrt{1 + \hat{k} \frac{1 - f}{f}} \le \rho_{\text{target}}\right)$$
via log-normal Monte Carlo draws (`n_sim = 2000`) of $\hat{k}$, either from prospective error-rate and correlation parameters or calibrated to the BCa bootstrap interval bounds from `catmech_audit_plan()`, generating visual power curves via `plot()`. Note that `catmech_audit_plan()` returned a `NO-GO` verdict in pre-specified validation testing (its median error in predicting interval width at a 5% audit was 11.3% to 16.8%, missing the 10% limit) and is therefore provided as an exploratory calculation rather than a validated planning tool.

## 3. Downstream multi-endpoint model pooling

Parameter estimates $\hat{\theta}_j$ and variance-covariance matrices $\hat{U}_j$ across $j \in \{1, \dots, m\}$ completed datasets are combined using Rubin's rules [@rubin1987multiple]:
$$\bar{\theta} = \frac{1}{m} \sum_{j=1}^m \hat{\theta}_j, \quad \bar{U} = \frac{1}{m} \sum_{j=1}^m \hat{U}_j, \quad B = \frac{1}{m - 1} \sum_{j=1}^m (\hat{\theta}_j - \bar{\theta})(\hat{\theta}_j - \bar{\theta})^\top$$
$$T = \bar{U} + \left(1 + \frac{1}{m}\right) B, \quad \hat{\lambda} = \frac{(1 + 1/m) \text{diag}(B)}{\text{diag}(T)}$$

To prevent undercoverage in finite clinical samples, `catmech` implements the Barnard & Rubin (1999) degrees of freedom adjustment [@barnard1999small]:
$$\nu_m = (m - 1) \hat{\lambda}^{-2}, \quad \nu_{obs} = \left(\frac{\nu_{com} + 1}{\nu_{com} + 3}\right) \nu_{com} (1 - \hat{\lambda}), \quad \nu_{BR} = \frac{\nu_m \nu_{obs}}{\nu_m + \nu_{obs}}$$
where $\nu_{com}$ is the complete-data degrees of freedom. Downstream model pooling covers three core clinical endpoint families:

* **Linear models (`catmech_pool_lm`)**: For continuous endpoints, using $\nu_{com} = N - p$ when unclustered (`cluster = NULL`) and applying closed-form delete-one-cluster jackknife (`cv3`) variance estimation with $\nu_{com} = C - 1$ when `cluster` is specified across $C$ clinical sites.
* **Generalized linear models (`catmech_pool_glm`)**: For binary (logistic) and count (Poisson) endpoints, using $\nu_{com} = N - p$ when unclustered and setting $\nu_{com} = C - 1$ with CR1S finite-sample-corrected Huber-White-Liang-Zeger cluster-robust sandwich covariance [@liang1986longitudinal] when `cluster` is specified across $C$ sites:
  $$\hat{V}_{\text{sandwich}} = \frac{C}{C - 1}\frac{N - 1}{N - p} \hat{J}^{-1} \left(\sum_{c=1}^C s_c s_c^\top\right) \hat{J}^{-1}$$
  where $\hat{J}^{-1}$ is the unclustered GLM covariance matrix and $s_c$ is the summed score residual vector in cluster $c$.
* **Cox proportional hazards models (`catmech_pool_coxph`)**: For time-to-event survival outcomes, using model-based inverse partial Fisher information with $\nu_{com} = d - p$ (where $d$ is the number of observed failure events and $p$ is the parameter count) when `cluster = NULL`, and computing cluster-summed Lin & Wei (1989) score-residual sandwich covariance [@lin1989robust] with $\nu_{com} = C - 1$ when `cluster` is specified across $C$ sites.

## 4. Diagnostic screens and visualizations

`catmech` provides unified S3 `plot()`, `tidy()`, and `glance()` methods:

* **Holdout cross-validation (`catmech_holdout_check()`)**: Partitions audited records into $K$ cross-validation folds, refits the workflow on the remaining $K - 1$ folds $\mathcal{D}_{-k}$, and computes tie-randomized empirical rank Probability Integral Transform (PIT) values $u_i$ across the $m$ completions [@dawid1984present]. It flags `miscalibration_detected` if a two-sided Kolmogorov-Smirnov test yields $p < 0.01$, if the 95% Wilson upper bound on empirical predictive coverage falls below the nominal level, or if the corruption probability calibration $z$-statistic exceeds $z_{0.995}$.
* **Recording disagreement audit (`catmech_recording_audit()`)**: Evaluates stratum-specific error rates for $a_s$ disagreements among $n_s$ audited records in stratum $s$ with exact 95% Jeffreys beta intervals ($\text{Beta}(a_s + 0.5, n_s - a_s + 0.5)$) against the global error rate.
* **Multi-scenario sensitivity envelopes (`catmech_sensitivity()`)**: Evaluates the stability of clinical conclusions across competing, non-identifiable assumptions (e.g., contamination mixture vs. linear reference modeling vs. selection shifts).

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

# 7. Side-by-side comparison of four audit strategies (validation_subsample,
#    naive_trust_recorded, linear_mi_recorded, and contamination_model)
comp <- catmech_compare_strategies(
  data = cohort,
  value = "lab_value",
  reference = "lab_ref",
  formula = outcome ~ lab_value + age + sex,
  m = 10,
  seed = 2026
)
plot(comp, term = "lab_value")

# 8. SHA-256 provenance report (supports FDA 21 CFR Part 11 documentation)
report <- catmech_audit_trail(
  object = fit,
  pooled = pool_survival,
  format = "markdown",
  path = "audit_trail.md"
)
```

The resulting pooled summary for the Cox proportional hazards model yields an auditable table reporting log-hazard point estimates, clustered Lin-Wei sandwich standard errors, Barnard-Rubin finite-sample degrees of freedom (with $\nu_{com} = C - 1 = 9$ for $C = 10$ sites), exponentiated hazard ratios (`estimate_exp`), and 95% confidence intervals:

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

To evaluate the statistical properties of `catmech`, pre-specified Monte Carlo simulation studies encompassing over 500,000 requests across 10 validation archives were conducted. In the primary contamination study (`validation/contamination-study`, 52,950 planned requests across 183 method-scenario cells spanning 4 primary 1,000-cohort cells and 23 150-cohort stress cells: 52,900 returned intervals, 16 expected design refusals at 5% and 10% audit fractions, 33 convergence refusals, and 1 singular raking fit), seven core analytical strategies were compared for estimating a true slope $\beta = 0.800$ in cohorts of size $N = 600$ with a 25% reference audit ($n_{\text{audit}} = 150$) and $m = 20$ imputations: **Oracle** (true $T_i$ known), **Flag Masking** (`v02_mask_normal`, setting flagged rows $S_i = 1$ to `NA` for normal MI while trusting unflagged rows $S_i = 0$ as truth), **Complete-Case** (`gold_only`, 25% audited rows only), **Discard Recorded MI** (`reference_normal`, `catmech_reference(use_recorded = FALSE)`), **Simple Linear MI** (`cole_mi`, `catmech_reference(use_recorded = TRUE)`), **Generalized Raking** (`raking`), and **`catmech` Contamination Mixture** (`catmech_contamination()`):

| Scenario & Strategy | Est. $\hat{\beta}$ | Bias (MCSE) | Emp. SD | Mean SE | Coverage | Mean Width | RMSE |
| :-------------------------- | :------: | :-------------------: | :------: | :------: | :------: | :------: | :------: |
| **Unflagged Errors (`missed_flags`, $N_{\text{sim}} = 1,000$)** | | | | | | | |
| Oracle (Latent $T$ Known) | 0.7978 | -0.0022 (0.0012) | 0.0394 | 0.0411 | 95.7% | 0.1614 | 0.0394 |
| Flag Masking (`v02_mask_normal`) | 0.7185 | -0.0815 (0.0015) | 0.0472 | 0.0444 | 55.9% | 0.1750 | 0.0942 |
| Complete-Case (`gold_only`) | 0.7988 | -0.0012 (0.0026) | 0.0815 | 0.0837 | 95.7% | 0.3309 | 0.0814 |
| Discard Recorded MI (`reference_normal`) | 0.7896 | -0.0104 (0.0021) | 0.0676 | 0.0684 | 95.4% | 0.2769 | 0.0684 |
| Simple Linear MI (`cole_mi`) | 0.7806 | -0.0194 (0.0018) | 0.0557 | 0.0528 | 92.4% | 0.2101 | 0.0590 |
| Generalized Raking (`raking`) | 0.7853 | -0.0147 (0.0020) | 0.0631 | 0.0589 | 92.9% | 0.2331 | 0.0648 |
| **`catmech` Contamination Mixture** | **0.7952** | **-0.0048 (0.0014)** | **0.0446** | **0.0467** | **94.9%** | **0.1842** | **0.0449** |
| **Gross & Small Errors (`error_gross_and_small`, $N_{\text{sim}} = 150$)** | | | | | | | |
| Oracle (Latent $T$ Known) | 0.7970 | -0.0030 (0.0034) | 0.0411 | 0.0412 | 94.7% | 0.1617 | 0.0411 |
| Flag Masking (`v02_mask_normal`) | 0.0609 | -0.7391 (0.0055) | 0.0674 | 0.0155 | 0.7% | 0.0614 | 0.7422 |
| Complete-Case (`gold_only`) | 0.7962 | -0.0038 (0.0064) | 0.0787 | 0.0845 | 95.3% | 0.3339 | 0.0785 |
| Discard Recorded MI (`reference_normal`) | 0.7890 | -0.0110 (0.0055) | 0.0671 | 0.0695 | 95.3% | 0.2816 | 0.0678 |
| Simple Linear MI (`cole_mi`) | 0.5646 | -0.2354 (0.0261) | 0.3200 | 0.0617 | 64.7% | 0.2562 | 0.3964 |
| Generalized Raking (`raking`) | 0.5156 | -0.2844 (0.0281) | 0.3436 | 0.0693 | 53.7% | 0.2740 | 0.4451 |
| **`catmech` Contamination Mixture** | **0.7883** | **-0.0117 (0.0035)** | **0.0432** | **0.0436** | **94.0%** | **0.1714** | **0.0447** |

These results establish both the strengths and the pre-specified failure boundaries of `catmech`:

* Trusting unflagged recordings (`v02_mask_normal`) suffers severe attenuation bias ($-0.0815$ in `missed_flags` and $-0.7391$ in `error_gross_and_small`), collapsing 95% coverage to 55.9% and 0.7%.
* Complete-case analysis (`gold_only`) and discarding recorded values (`reference_normal`) eliminate bias but forfeit precision (mean CI widths of 0.3309 and 0.2769 in `missed_flags`).
* Simple linear MI (`cole_mi`) and generalized raking (`raking`) break down when gross errors are mixed with small errors (`error_gross_and_small`), suffering severe bias ($-0.2354$ and $-0.2844$) and undercoverage (64.7% and 53.7%).
* The `catmech` contamination mixture model preserves nominal coverage (94.9% in `missed_flags` and 94.0% in `error_gross_and_small`) while narrowing confidence intervals by 33.5% relative to discarding recorded values and by 44.3% relative to complete-case analysis in `missed_flags` (`0.1842` vs. `0.3309`). Across all four primary 1,000-cohort cells (`clean`, `flagged_errors`, `missed_flags`, `site_shift_declared`), `catmech_contamination()` passed the $[0.925, 0.975]$ Wilson coverage screen in three cells (96.5%, 96.3%, and 94.9%) and narrowly missed in `site_shift_declared` at 93.7% (Wilson $[0.9202, 0.9505]$, bias $-0.0085$, MCSE $0.0016$), while narrowing intervals by 30% to 38% relative to `reference_normal`.
* Under the pre-specified stress-cell failure threshold (coverage $< 90.0\%$ or $|\text{bias}| > 0.030$), `catmech_contamination()` failed under heavy-tailed target distributions (`heavy_tail_target`: bias $+0.0330$, coverage 87.3%, where genuine extremes are shrunk as errors while `reference_normal` achieved 93.3% and `cole_mi` 92.0%), non-ignorable audit selection on latent truth with errors present (`nonignorable_missed`: bias $-0.0506$, coverage 79.3%, vs. 87.3% for `raking`), and undeclared reference measurement noise (`noisy_undeclared`: bias $-0.0261$, coverage 88.7%, repaired to $-0.0071$ bias and 91.8% coverage when `reference_sd` is declared in `noisy_declared`), was borderline under bimodal two-sided errors (`error_bimodal`: bias $-0.0293$, coverage 90.0%), and experienced 17/150 (11.3%) refusals at a sparse 5% audit fraction (`reference_05pct`, 91.7% coverage) and in `nonignorable_clean` (92.5% coverage). In the 8-scenario holdout check evaluation (40 cohorts each), `catmech_holdout_check()` produced 0/40 false alarms in `clean` and `missed_flags`, but had low power against slope bias (flagging 18/40 in `site_shift_undeclared`, 1/40 in `heavy_tail_target`, and 0/40 in `nonlinear_target`), so a clean hold-out check is a weak screen rather than validation.

## 2. Real-world empirical benchmarks

The package was evaluated on two large-scale public paired-discrepancy datasets without synthetic corruption (`validation/REAL_DATA_SYNTHESIS.md`, 33,800 planned requests across 1,000 resampled cohorts of size $N = 600$ at the 25% audit fraction and 200 cohorts at 5%, 10%, 50%, and targeted fractions: 33,613 returned intervals, 76 unsupported-design refusals, 83 MCMC convergence refusals, and 28 singular raking errors) on ClinicalTrials.gov ($N = 27,208$) and NHANES 2017--2018 [@nhanes2018] ($N = 3,908$):

| Dataset | Metric | Naive (Trust Rec.) | Discard Rec. | Contamination | Linear MI | Full Ref. (Oracle) |
| :------------------------ | :-----------: | :-------------: | :-----------: | :-----------------: | :---------: | :--------------: |
| **ClinicalTrials.gov** ($N = 27,208$) | Bias | +0.0043 | -0.0003 | +0.0068 | +0.0001 | -0.0003 |
| Discrepancy rate: 14.2% | Coverage | 91.4% | 92.3% | 91.2% | 91.8% | 91.2% |
| Target: 0.0357 log mo / log N | Mean Width | 0.097 | 0.205 | 0.106 | 0.102 | 0.096 |
| | RMSE | 0.0282 | 0.0539 | 0.0308 | 0.0286 | 0.0275 |
| **NHANES 2017--2018** ($N = 3,908$) | Bias | -0.0180 | -0.0063 | -0.0066 | -0.0025 | +0.0001 |
| Discrepancy rate ($>1$ kg): 67.3% | Coverage | 92.7% | 94.2% | 94.5% | 93.9% | 94.9% |
| Target: 0.127 mmHg / kg | Mean Width | 0.149 | 0.297 | 0.155 | 0.154 | 0.140 |
| | RMSE | 0.0422 | 0.0785 | 0.0408 | 0.0409 | 0.0357 |

On both datasets, retaining the recorded values cut RMSE by 43% to 48% compared to discarding unaudited records (`Discard Rec.`). On ClinicalTrials.gov, the Oracle itself covered 91.2% due to heteroskedastic residuals (reported SE 0.0244 vs. empirical SD 0.0275), and `Naive (Trust Rec.)` had slightly lower RMSE (0.0282, bias $+0.0043$) than both `Linear MI` (0.0286, bias $+0.0001$) and `Contamination` (0.0308, bias $+0.0068$; at a 5% audit, `Contamination` returned 164/200 intervals with bias $+0.0159$ vs. 200/200 and $+0.0062$ for `Linear MI`). On NHANES at 25%, both `Linear MI` and `Contamination` corrected the naive bias ($-0.0180$). Because discrepancies in both datasets were diffuse relative to conditional target spread, per-record corruption probabilities $p_i$ had near-chance discrimination (AUC 0.543 on ClinicalTrials.gov and 0.500 on NHANES, with posterior mean RMSE worse than raw recorded values: 0.277 vs. 0.252 and 4.81 vs. 4.73 kg) and must not be used to label individual records.

In addition, two frozen studies benchmarked clustered facility recounts from the Rwanda Health Management Information System (HMIS) [@nshimyiryo2020health] and psychology effect-size recounts [@maassen2020reproducibility]:

1. **Clustered variance & multilevel MI (`validation/VARIANCE_MI_SYNTHESIS.md`, 148,500 requests: 148,486 returned intervals and 14 refusals [12 rank-deficient audited designs on `iron` and 2 zero-residual-variance fits on `stillbirths`])**: Standard cluster sandwich SEs (`cr1s`) and bias-reduced linearization with Satterthwaite df (`cr2`, which inflated widths by up to 95% and passed only 23/44 simulated cells) failed eligibility. Delete-one-cluster jackknife (`cv3`) restored 100% complete-data recount coverage (`gold_only`) to 92.0%--97.5% across the 11 Rwanda HMIS indicators and narrowly passed under random-intercept outcome imputation (9/11 Rwanda indicators $\ge 92\%$, with `syphilis` at 90.0% and `iron` at 91.9%, and 36/44 simulated cells, reducing mean distance from 95% coverage from 3.18 to 1.47 percentage points). However, `cv3` did not rescue unclustered linear MI on the Rwanda indicators (coverage stayed at 81.0% to 86.1% on the worst items, leaving the `PARTIAL` accuracy verdict unchanged), and when the audited variable was a **predictor** rather than the outcome, clustered random-intercept MI (`catmech_reference(cluster = )`) attenuated the predictor slope toward zero by up to $-16.1\%$ (84.4% coverage at 25% audit under `cv3`, vs. $-3.0\%$ bias and 95.4% coverage when fitting `catmech_reference(use_recorded = TRUE)` without `cluster` and pooling with `catmech_pool_lm(cluster = )`).
2. **Real-audit slope recovery (`validation/REAL_AUDIT_SYNTHESIS.md`, 45,600 requests)**: Under pre-specified Rule 2 on Rwanda HMIS at a 25% audit, `catmech_contamination()` was **`INCONCLUSIVE`**, failing coverage on all three primary slopes (`PNC`: 71.6% coverage, 468/500 [93.6%] returned, RMSE 0.053 vs. 0.072 for `Linear MI`; `Syphilis`: 66.6% coverage, 404/500 [80.8%] returned, RMSE 0.195 vs. 0.182; `Iron`: 35.4% coverage, 424/500 [84.8%] returned, RMSE 0.112 vs. 0.156) and on the size-dependent `ANC 4th visit` stress slope (49.7% coverage, RMSE 0.124 vs. 0.085), being closest to 95% coverage in only 1 of 20 cells (vs. 10 of 20 for `Linear MI`, worst case 79.5%, with per-record error AUC 0.501 to 0.582). Under Rule 1, a 25% audit with `Linear MI` beat raw recorded values clearly on only 1 of 4 slopes (`ANC 4th visits`, RMSE ratio $0.458$), tied on `PNC` ($1.009$ on 2,000 fresh seeds), and cost more error on `Syphilis` ($2.393\times$) and `Iron` ($1.125\times$). Under Rule 3 (498 psychology effect-size recounts [@maassen2020reproducibility] with near-zero reported bias), `Linear MI` at 25% inflated recovery RMSE by $2.705\times$ (SMD) and $1.578\times$ ($r$), whereas `catmech_contamination()` stayed near recorded values ($1.306\times$ and $0.950\times$).

## 3. Simulation-to-practice divergence and methodological guidance

A central methodological insight established by these empirical benchmarks is the divergence between synthetic Monte Carlo evaluations and real-world clinical audits:

* **Synthetic gross errors favor mixture modeling**: When data errors represent discrete, bimodal corruption events (e.g., transcription slips, decimal shifts, or assay calibration jumps) under Gaussian target residuals and ignorable audit selection, the Bayesian contamination mixture model eliminates severe attenuation bias (bias of $-0.0117$ vs. $-0.2354$ for linear MI in `error_gross_and_small`) and narrows confidence intervals by 44.3% over complete-case analysis in `missed_flags` (while remaining synthetic-only because real-world audits lacked separable gross errors).
* **Real-world subtle or moderate errors favor linear multiple imputation**: In actual registry and facility audits where discrepancies are diffuse and continuous (e.g., self-report rounding, protocol differences, or moderate recount discrepancies), simple linear imputation (`catmech_reference(use_recorded = TRUE)`) matches or exceeds mixture modeling in predictive accuracy and interval coverage with fewer parametric assumptions and zero MCMC convergence overhead.

This empirical divergence highlights why `catmech` incorporates prospective audit planning (`catmech_audit_power()`), diagnostic discrepancy auditing (`catmech_recording_audit()`), predictive cross-validation (`catmech_holdout_check()`), and automated strategy comparison (`catmech_compare_strategies()`), ensuring that trialists let empirical data guide the choice of model rather than relying on unverified assumptions.

# Regulatory documentation and data governance

Clinical trial analyses are expected to be documented and reproducible: FDA 21 CFR Part 11 [@fda2003part11] governs electronic records and audit trails, and ICH E9 [@ich1998e9] asks that data handling and statistical methods be prespecified and documented. `catmech_audit_trail()` supports that documentation; it does not certify compliance, which depends on the validated system, access controls, and procedures around the analysis. It generates provenance records in Markdown, JSON, or standalone HTML. The report documents:

* Full execution environment metadata (R version, operating system platform, package version, UTC execution timestamp).
* Cryptographic SHA-256 checksums of the raw input dataset, the mathematical assumption specification, and each of the $m$ completed datasets.
* Model convergence diagnostics recorded on the audit trail (`convergence_passed` pass/fail status and `mcmc_iterations` count, with full Gelman-Rubin $\hat{R}$ and bulk/tail effective sample size vectors stored in `result$model_diagnostics`).
* Downstream pooled parameter tables and multi-scenario sensitivity envelopes.

# Quality assurance and verification

`catmech` has been developed under an automated software verification harness:

* **Unit tests**: 1,754 unit tests executed via `testthat`, achieving 0 failures, 0 warnings, and 0 skips.
* **CRAN conformance**: `R CMD check --as-cran` reports 0 errors and 0 warnings.
* **Independent artifact verification**: Verified against an independent Python verification suite (`verify-artifacts.py`, `verify-release-archive.py`, and `verify_paper_numbers.py`) confirming exact correspondence across 10 frozen simulation and empirical benchmark studies and over 500,000 replicate requests.

# Availability and dependencies

`catmech` is open-source software licensed under the MIT license. It requires R ($\ge$ 4.1.0) and a C++17 compiler, utilizing the non-base CRAN packages `Rcpp`, `RcppArmadillo`, `digest`, `jsonlite`, `mice`, `pan`, `posterior`, `survival`, `truncnorm`, and `yaml` (alongside base R packages `graphics`, `grDevices`, `parallel`, `stats`, and `utils`). Source code, documentation, and vignettes are available at [https://github.com/patrickiben/catmech](https://github.com/patrickiben/catmech).

# References
