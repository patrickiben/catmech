# Auditable Clinical Data Verification and Analysis with catmech

## Introduction

In multi-center clinical trials, epidemiological studies, and electronic
health record (EHR) registries, continuous laboratory measurements (such
as troponin, serum creatinine, or viral load) frequently suffer from
data quality issues:

1.  **Missing values** due to uncollected or lost samples.
2.  **Gross recording errors and data contamination** resulting from
    transcription mistakes, assay batch errors, or EHR unit conversion
    discrepancies.
3.  **Internal reference audits** where a secondary gold-standard
    measurement is audited on a subset of patient records.

The **`catmech`** package provides an auditable framework for completing
contaminated measurements, checking predictive calibration, and pooling
downstream linear, generalized linear (GLM), and survival (Cox PH)
regressions under Rubin’s rules with Barnard–Rubin degrees of freedom.
It also has a rough, unvalidated aid for planning audit sizes. The
package has had no outside peer review or clinical validation; the
[README](https://github.com/patrickiben/catmech#readme) gives the
evidence and its limits.

------------------------------------------------------------------------

## 1. Simulating a Multi-Center Clinical Registry

We simulate a realistic 10-hospital clinical registry with $N = 250$
patients. Each patient has: \* An audited gold-standard reference
measurement (`lab_true`). \* A recorded measurement (`lab_value`)
corrupted in a fraction of patients. \* Patient covariates (`age`,
`sex`, `site`). \* Three downstream clinical endpoints: - Continuous
outcome (`cardiac_index`). - Binary 30-day major adverse cardiac event
(`mace_30d`). - Time-to-event overall survival (`time`, `death_event`).

``` r
library(catmech)
library(survival)

set.seed(2026)
n <- 250
site <- sample(1:10, n, replace = TRUE)
age <- round(rnorm(n, mean = 62, sd = 10))
sex <- sample(c("M", "F"), n, replace = TRUE, prob = c(0.55, 0.45))

# Latent true biomarker
lab_true <- rnorm(n, mean = 6.0 + 0.03 * (age - 62), sd = 1.2)

# Corrupted recorded measurement (12% gross contamination)
corrupted <- rbinom(n, 1, prob = 0.12) == 1
lab_recorded <- lab_true
lab_recorded[corrupted] <- lab_recorded[corrupted] + rnorm(sum(corrupted), mean = 4.5, sd = 2.0)

# Reference audit on a 25% random pilot
audited <- rbinom(n, 1, prob = 0.25) == 1
lab_reference <- ifelse(audited, lab_true, NA_real_)

# Clinical outcomes
cardiac_index <- 2.8 - 0.25 * lab_true + 0.01 * (age - 62) + rnorm(n, sd = 0.5)
pr_mace <- plogis(-2.0 + 0.45 * (lab_true - 6.0) + 0.02 * (age - 62))
mace_30d <- rbinom(n, 1, prob = pr_mace)

# Survival times
haz <- exp(0.35 * (lab_true - 6.0) + 0.02 * (age - 62)) * 0.08
surv_time <- rexp(n, rate = haz)
censor_time <- runif(n, 1, 10)
time <- pmin(surv_time, censor_time)
death_event <- as.numeric(surv_time <= censor_time)

df <- data.frame(
  patient_id = 1:n,
  site = site,
  age = age,
  sex = factor(sex),
  lab_value = lab_recorded,
  lab_ref = lab_reference,
  cardiac_index = cardiac_index,
  mace_30d = mace_30d,
  time = time,
  death_event = death_event
)
head(df[, c("patient_id", "site", "age", "lab_value", "lab_ref", "mace_30d")])
#>   patient_id site age lab_value lab_ref mace_30d
#> 1          1    9  43  6.111779      NA        0
#> 2          2    1  66  6.970042      NA        0
#> 3          3    6  70  7.333651      NA        1
#> 4          4    4  42  5.012186      NA        0
#> 5          5    5  60  4.937939      NA        0
#> 6          6   10  52  4.324780      NA        0
```

------------------------------------------------------------------------

## 2. Fitting the Reference-Calibrated Contamination Model

We fit
[`catmech_contamination()`](https://patrickiben.github.io/catmech/reference/catmech_contamination.md)
using the audited subset to model both the probability of corruption and
the measurement error mixture.

The biomarker is a *predictor* in every analysis below, so the
imputation model must include each analysis outcome. Leaving an outcome
out biases the biomarker’s coefficient toward zero; the pooling
functions warn when they see that. For the survival analysis, the
outcome enters as the event indicator and the Nelson–Aalen cumulative
hazard (White and Royston, 2009), which
[`mice::nelsonaalen()`](https://amices.org/mice/reference/nelsonaalen.html)
computes:

``` r
df$hazard <- mice::nelsonaalen(df, time, death_event)

fit_contam <- catmech_contamination(
  data = df,
  value = "lab_value",
  reference = "lab_ref",
  predictors = c("age", "cardiac_index", "mace_30d", "death_event", "hazard"),
  m = 4,
  seed = 42
)
print(fit_contam)
#> <catmech_result: reference_contamination>
#>   Target value:      lab_value
#>   Reference subset:  lab_ref (64/250 rows, 25.6%)
#>   Predictors:        age, cardiac_index, mace_30d, death_event, hazard
#>   Completions (m):   4
#>   Imputed rows:      186/250
#>   Model convergence: passed
```

The MCMC sampler generates $m = 4$ completions of the target biomarker,
reflecting parameter uncertainty. Four keeps this vignette fast; use 20
or more in practice.

------------------------------------------------------------------------

## 3. Diagnostic Checks & Base R Visualizations

### Predictive Holdout Calibration

[`catmech_holdout_check()`](https://patrickiben.github.io/catmech/reference/catmech_holdout_check.md)
tests whether the imputation model is well-calibrated using
cross-validation over the audited reference cases:

``` r
holdout <- catmech_holdout_check(
  fit_contam,
  folds = 3L,
  m = 20L,
  seed = 99
)
print(holdout)
#> -- Catmech Holdout Cross-Validation Check --
#> Status:              no_miscalibration_detected
#> Folds completed:     3 (holdout rows: 64, partitioned by: row)
#> Predictive coverage: 0.9844 (nominal: 0.9) (uniform PIT KS p-value: 0.1135)
#> Predictive RMSE:     model = 0.3086 vs recorded = 1.7378
#> Mean interval width: 0.5999
#> Flags:               none
plot(holdout)
```

![](catmech-clinical-audit_files/figure-html/holdout_check-1.png)

The left panel displays the uniform randomized Probability Integral
Transform (PIT) histogram with a Kolmogorov–Smirnov test, while the
right panel verifies mean predictive calibration.

------------------------------------------------------------------------

## 4. Rough Audit-Size Planning (Not Validated)

Before an audit,
[`catmech_audit_power()`](https://patrickiben.github.io/catmech/reference/catmech_audit_power.md)
gives a rough idea of its size. It is **not validated**, and the planner
it builds on,
[`catmech_audit_plan()`](https://patrickiben.github.io/catmech/reference/catmech_audit_plan.md),
returned NO-GO under a rule fixed in advance. Its default factor treats
the unexplained variation as if it sat in erroneous records that can be
identified, which is optimistic, so the printout also shows the audit
that linear imputation needs under the double-sampling regression
estimator at the same correlation:

``` r
pwr <- catmech_audit_power(
  N = 1000,
  target_ratio = 1.25,
  power_target = 0.80,
  error_rate = 0.12,
  correlation = 0.82
)
print(pwr)
#> -- Catmech Audit-Size Planning Aid (not validated) --
#> Cohort Size (N):    1000
#> Target Ratio:       <= 1.25 (within 25% of full audit width)
#> Target share of k:  80%
#> Two-Phase Factor k: 0.028 [0.010, 0.078] (prior_specification)
#> 
#>  fraction n_audit expected_ratio ratio_lower ratio_upper power
#>      0.02      20          1.587       1.216       2.194 0.046
#>      0.04      40          1.324       1.111       1.694 0.362
#>      0.06      60          1.222       1.074       1.490 0.680
#>      0.08      80          1.168       1.055       1.377 0.853
#>      0.10     100          1.134       1.043       1.304 0.938
#>      0.12     121          1.111       1.035       1.253 0.973
#>      0.14     140          1.094       1.030       1.216 0.986
#>      0.16     160          1.081       1.025       1.187 0.994
#>      0.18     180          1.070       1.022       1.164 0.998
#>      0.20     200          1.062       1.019       1.145 0.998
#>      0.22     220          1.055       1.017       1.130 0.999
#>      0.24     240          1.050       1.015       1.116 0.999
#>      0.26     260          1.045       1.014       1.105 1.000
#>      0.28     280          1.041       1.012       1.096 1.000
#>      0.30     301          1.037       1.011       1.087 1.000
#>      0.32     320          1.034       1.010       1.080 1.000
#>      0.34     340          1.031       1.009       1.073 1.000
#>      0.36     361          1.028       1.009       1.067 1.000
#>      0.38     380          1.026       1.008       1.062 1.000
#>      0.40     400          1.024       1.007       1.057 1.000
#> 
#> Recommended Audit Size: n = 80 records (8.0% audit fraction)
#> Meets the target in 85.3% of k draws (expected ratio: 1.168)
#> Linear imputation (double-sampling regression, k = 0.229): about 29.0% (n = 290) for ratio <= 1.25, ignoring uncertainty in k.
#> Power here is the share of plausible k values meeting the target, not the power of a test.
#> Not validated; see ?catmech_audit_power for the assumptions.
plot(pwr)
```

![](catmech-clinical-audit_files/figure-html/audit_power-1.png)

The curve shows, for each audit size, the share of plausible values of
the two-phase factor for which the expected interval width stays within
25% of a full audit’s. It guarantees nothing: treat it, and the
linear-imputation figure, as a starting point for planning.

------------------------------------------------------------------------

## 5. Downstream Multi-Endpoint Model Pooling

`catmech` supports downstream pooling across diverse outcome types:

### A. Linear Regression (`catmech_pool_lm`)

For continuous outcomes, pooling estimates slope parameters with Rubin’s
combining rules and cluster-robust standard errors:

``` r
res_lm <- catmech_pool_lm(
  fit_contam,
  cardiac_index ~ lab_value + age,
  cluster = "site"
)
print(res_lm)
#> -- Catmech Pooled Linear Model --
#> Completions (m): 4 | Observations (n): 250
#> 
#>         term estimate std.error     df   lower   upper
#>  (Intercept)   2.1856    0.4084 7.4740  1.2321  3.1390
#>    lab_value  -0.1739    0.0278 6.2173 -0.2414 -0.1065
#>          age   0.0024    0.0057 7.3694 -0.0110  0.0157
#>                                  interpretation
#>  conditional_contamination_model_cluster_robust
#>  conditional_contamination_model_cluster_robust
#>  conditional_contamination_model_cluster_robust
```

### B. Generalized Linear Models (`catmech_pool_glm`)

For binary outcomes (such as 30-day MACE),
[`catmech_pool_glm()`](https://patrickiben.github.io/catmech/reference/catmech_pool_glm.md)
fits logistic regression across all completions, reporting odds ratios
via `exponentiate = TRUE`:

``` r
res_glm <- catmech_pool_glm(
  fit_contam,
  mace_30d ~ lab_value + age,
  family = binomial(),
  exponentiate = TRUE
)
print(res_glm)
#> -- Catmech Pooled GLM: binomial (logit link) --
#> Completions (m): 4 | Observations (n): 250
#> 
#>         term estimate std.error       df    lower   upper estimate_exp
#>  (Intercept)  -7.1721    1.8284 146.5547 -10.7855 -3.5588       0.0008
#>    lab_value   0.2408    0.1842  50.3994  -0.1291  0.6107       1.2723
#>          age   0.0575    0.0228 244.8883   0.0126  0.1024       1.0592
#>  lower_exp upper_exp
#>     0.0000    0.0285
#>     0.8789    1.8417
#>     1.0127    1.1078
```

### C. Survival Proportional Hazards (`catmech_pool_coxph`)

For time-to-event outcomes,
[`catmech_pool_coxph()`](https://patrickiben.github.io/catmech/reference/catmech_pool_coxph.md)
fits Cox proportional hazards models using observed event counts for
Barnard–Rubin degrees of freedom:

``` r
res_cox <- catmech_pool_coxph(
  fit_contam,
  Surv(time, death_event) ~ lab_value + age,
  exponentiate = TRUE,
  cluster = "site"
)
print(res_cox)
#> -- Catmech Pooled Cox PH Model --
#> Completions (m): 4 | Observations (n): 250 | Events: 83
#> 
#>       term estimate std.error     df   lower  upper estimate_exp lower_exp
#>  lab_value   0.4160    0.1066 3.1900  0.0881 0.7440       1.5159    1.0921
#>        age   0.0214    0.0102 7.1523 -0.0027 0.0454       1.0216    0.9973
#>  upper_exp
#>     2.1043
#>     1.0464
```

------------------------------------------------------------------------

## 6. Sensitivity Envelope Across Modeling Assumptions

Investigators can prespecify alternative modeling assumptions (e.g.,
contamination mixture model vs. reference model) and evaluate whether
the conclusion is robust across all scenarios:

``` r
fit_ref <- catmech_reference(
  data = df,
  value = "lab_value",
  reference = "lab_ref",
  predictors = c("age", "cardiac_index", "mace_30d", "death_event", "hazard"),
  use_recorded = TRUE,
  m = 4,
  seed = 42
)

sens <- catmech_sensitivity(
  results = list(contamination = fit_contam, reference = fit_ref),
  formula = cardiac_index ~ lab_value + age,
  term = "lab_value",
  method = "pool_lm",
  reference = 0
)
print(sens)
#> -- Catmech Sensitivity Analysis --
#> Status:              direction_stable_over_supplied_scenarios (complete: TRUE)
#> Successful envelope: [-0.2307, -0.092]
#> Scenarios:           2 evaluated
#> 
#>       scenario status estimate std.error   lower   upper direction
#>  contamination     ok  -0.1739    0.0286 -0.2307 -0.1171     below
#>      reference     ok  -0.1606    0.0326 -0.2291 -0.0920     below
```

------------------------------------------------------------------------

## 7. Provenance Record

[`catmech_audit_trail()`](https://patrickiben.github.io/catmech/reference/catmech_audit_trail.md)
writes a record of one run: R and package versions, SHA-256 hashes of
the input data and every completion (taken over the serialized R data
frames, so they can be recomputed only in R) and of the assumption
specification’s canonical JSON text, the settings, the convergence
status and the pooled results. A rerun can be compared with it. It is
not tamper-evident, since anyone who edits the data can regenerate it,
and it can support documentation under FDA 21 CFR Part 11 or ICH E9
without making a system compliant:

``` r
trail <- catmech_audit_trail(
  object = fit_contam,
  pooled = res_cox,
  sensitivity = sens,
  format = "markdown"
)
cat(substr(as.character(trail), 1, 1200), "\n... [truncated for display]")
#> # Catmech Audit Trail & Provenance Report
#> 
#> **Generated:** 2026-11-08T16:29:01Z UTC | **Catmech Version:** 1.0.1
#> **R Version:** 4.6.2 | **Platform:** x86_64-pc-linux-gnu
#> 
#> ---
#> 
#> ## 1. Data and Specification Hashes (SHA-256)
#> 
#> - **Original Input Data Hash:** `957b56100d7d74bdb322e9bc67803d15c53aeff694d92de87dd4688608879b2f`
#> - **Assumption Specification Hash:** `67a39bb736d40bf6fbf1d0840fe417ba43247a444c6bf36b44dd1b3036019ce3`
#> - **Completed Dataset Checksums:**
#> 
#>   - Completion `m=1`: `cc84765505f9b1b3550973ba9fe392eb1b80bc623aaa095e6802dca20ae8d5a7`
#>   - Completion `m=2`: `51ba3565b86d4619a0ccae7bf978f3748b2c80ec23a72bc00563c26f26fc22ad`
#>   - Completion `m=3`: `3313ab9c8b71581798aad4e0988ba367e48a9cce93aa1e2369d7a67e465f6807`
#>   - Completion `m=4`: `fc964c70c9afb108b094f5dadd34cf9c248a6747deeef550682b3e668b7ee4ba`
#> 
#> ## 2. Cohort Specification & Imputation Profile
#> 
#> | Parameter | Value |
#> | :--- | :--- |
#> | Engine | `reference_contamination` |
#> | Target Value Column | `lab_value` |
#> | Reference Column | `lab_ref` |
#> | Predictors | age, cardiac_index, mace_30d, death_event, hazard |
#> | Cluster Variable | `none` |
#> | Total Observations | 250 |
#> | Completions (m) | 4 |
#> | Imputed Records | 186 (74.4%) |
#> | 
#> ... [truncated for display]
```

The record makes a rerun checkable; it does not by itself make an
analysis reproducible or compliant.

------------------------------------------------------------------------

## 8. Visualizing Pooled Models with Forest Plots

Clinical trial readouts and regulatory reports frequently require visual
representation of effect sizes and confidence intervals. `catmech`
provides native base R
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for all
pooled regression models:

``` r
# Forest plot for pooled logistic regression (Odds Ratios)
plot(res_glm, exponentiate = TRUE, main = "30-Day MACE Odds Ratios (95% CI)")
```

![](catmech-clinical-audit_files/figure-html/forest_plots-1.png)

``` r

# Forest plot for pooled Cox survival model (Hazard Ratios)
plot(res_cox, exponentiate = TRUE, main = "Overall Survival Hazard Ratios (95% CI)")
```

![](catmech-clinical-audit_files/figure-html/forest_plots-2.png)

------------------------------------------------------------------------

## 9. Side-by-Side Audit Strategy Comparison

To evaluate how much precision is recovered and whether point estimates
shift across different ways of handling un-audited records,
[`catmech_compare_strategies()`](https://patrickiben.github.io/catmech/reference/catmech_compare_strategies.md)
evaluates all four core workflows side by side. By default its
imputation uses every variable in the formula, the outcome included:

``` r
comp <- catmech_compare_strategies(
  data = df,
  value = "lab_value",
  reference = "lab_ref",
  formula = cardiac_index ~ lab_value + age,
  m = 4,
  seed = 2026
)
print(comp)
#> -- Catmech Audit Strategy Comparison --
#> Cohort Size (N): 250 | Audited Records (n): 64 (25.6% audit fraction)
#> Imputation predictors: cardiac_index, age
#> 
#>              strategy        term estimate std.error ci_width rel_efficiency
#>  validation_subsample (Intercept)   2.8243    0.5965   2.3854         1.0000
#>  validation_subsample   lab_value  -0.2098    0.0540   0.2158         1.0000
#>  validation_subsample         age  -0.0030    0.0082   0.0326         1.0000
#>  naive_trust_recorded (Intercept)   1.5703    0.2416   0.9517         6.0957
#>  naive_trust_recorded   lab_value  -0.0493    0.0171   0.0673         9.9663
#>  naive_trust_recorded         age   0.0006    0.0036   0.0142         5.1305
#>    linear_mi_recorded (Intercept)   2.4496    0.2891   1.1722         4.2570
#>    linear_mi_recorded   lab_value  -0.1976    0.0300   0.1231         3.2413
#>    linear_mi_recorded         age   0.0008    0.0034   0.0136         5.6694
#>   contamination_model (Intercept)   2.1943    0.2493   0.9821         5.7249
#>   contamination_model   lab_value  -0.1760    0.0270   0.1063         4.0050
#>   contamination_model         age   0.0024    0.0034   0.0134         5.7701

# Forest plot comparing all 4 strategies for the target measurement
plot(comp, term = "lab_value", main = "Audit Strategy Comparison: Biomarker Effect on Cardiac Index")
```

![](catmech-clinical-audit_files/figure-html/strategy_comparison-1.png)

------------------------------------------------------------------------

## 10. Practitioner Decision Flowchart

This table summarizes what the package’s own tests found. The
[README](https://github.com/patrickiben/catmech#readme) gives the
evidence, its limits and the cases where it was inconclusive.

| Scenario / Error Pattern                                              | Recommended Method                                                                                              | Rationale & Benchmark Findings                                                                                                                                                                                                                                                                      |
|-----------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Subtle or diffuse errors** (transcription slips, calibration drift) | `catmech_reference(use_recorded = TRUE)`                                                                        | On real recording differences (ClinicalTrials.gov and NHANES), linear MI matched or beat the contamination model, with fewer assumptions.                                                                                                                                                           |
| **Gross / separable spikes** (10x decimal errors, unit confusion)     | [`catmech_contamination()`](https://patrickiben.github.io/catmech/reference/catmech_contamination.md)           | The contamination mixture model explicitly estimates per-record error probability, achieving 94.0% coverage vs. 64.7% for linear MI under gross synthetic errors.                                                                                                                                   |
| **Multi-center / clustered sites**                                    | `catmech_reference(use_recorded = TRUE)` + `cluster = "site"` in the pooling function                           | The pooled variance is a delete-one-cluster jackknife (`cv3`) for linear models. `catmech_reference(cluster = )` is experimental, came out PARTIAL in its planned tests, and biased slopes when the audited value was a predictor, so use it only as a check when the audited value is the outcome. |
| **Audit size planning**                                               | `catmech_audit_power(N, target_ratio)`                                                                          | Rough and not validated; its default factor is optimistic, so also read its linear-imputation figure. [`catmech_audit_plan()`](https://patrickiben.github.io/catmech/reference/catmech_audit_plan.md) returned NO-GO.                                                                               |
| **Sensitivity check across strategies**                               | [`catmech_compare_strategies()`](https://patrickiben.github.io/catmech/reference/catmech_compare_strategies.md) | Compares validation subsample, naive analysis, linear MI, and contamination model side by side, imputing with the outcome included.                                                                                                                                                                 |
| **Provenance**                                                        | [`catmech_audit_trail()`](https://patrickiben.github.io/catmech/reference/catmech_audit_trail.md)               | Records SHA-256 hashes of the input data, settings and completions with the pooled results. Supports documentation; does not establish FDA 21 CFR Part 11 or ICH E9 compliance.                                                                                                                     |
