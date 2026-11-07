# Package index

## Audit Planning & Sample Size

Calculate required audit sample sizes and prospective statistical power.

- [`catmech_audit_power()`](https://patrickiben.github.io/catmech/reference/catmech_audit_power.md)
  : Prospective sample size and statistical power calculation for
  clinical data audits
- [`catmech_audit_plan()`](https://patrickiben.github.io/catmech/reference/catmech_audit_plan.md)
  : Plan how many records to audit from a random pilot
- [`plot(`*`<catmech_audit_power>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_audit_power.md)
  : Plot an audit power curve
- [`plot(`*`<catmech_audit_plan>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_audit_plan.md)
  : Plot a catmech audit plan

## Imputation & Error Models

Core algorithms for audit verification and multiple imputation.

- [`catmech_reference()`](https://patrickiben.github.io/catmech/reference/catmech_reference.md)
  : Complete a target using an independently measured reference subset
- [`catmech_contamination()`](https://patrickiben.github.io/catmech/reference/catmech_contamination.md)
  : Complete a target using a reference-calibrated recording-error model
- [`catmech_impute()`](https://patrickiben.github.io/catmech/reference/catmech_impute.md)
  : Generate auditable sensitivity completions
- [`catmech_pmm()`](https://patrickiben.github.io/catmech/reference/catmech_pmm.md)
  : Audited predictive mean matching after masking quality flags
- [`catmech_model()`](https://patrickiben.github.io/catmech/reference/catmech_model.md)
  : Gaussian completion with an explicit category observation model
- [`catmech_spec()`](https://patrickiben.github.io/catmech/reference/catmech_spec.md)
  : Declare assumptions for interval-constrained sensitivity completion
- [`read_catmech_spec()`](https://patrickiben.github.io/catmech/reference/read_catmech_spec.md)
  : Read a strict YAML assumption specification
- [`write_catmech()`](https://patrickiben.github.io/catmech/reference/write_catmech.md)
  : Write a self-contained result bundle
- [`catmech_derive_category()`](https://patrickiben.github.io/catmech/reference/catmech_derive_category.md)
  : Derive exact interval categories
- [`mice.impute.catmech_threshold_corruption()`](https://patrickiben.github.io/catmech/reference/mice.impute.catmech_threshold_corruption.md)
  : Category-constrained predictive mean matching with sensitivity
  replacement

## Downstream Regression Pooling

Rubin pooling for linear, generalized linear (GLM), and Cox proportional
hazards survival models.

- [`catmech_pool_lm()`](https://patrickiben.github.io/catmech/reference/catmech_pool_lm.md)
  : Pool a prespecified linear analysis of completed datasets
- [`catmech_pool_glm()`](https://patrickiben.github.io/catmech/reference/catmech_pool_glm.md)
  : Pool a prespecified generalized linear analysis of completed
  datasets
- [`catmech_pool_coxph()`](https://patrickiben.github.io/catmech/reference/catmech_pool_coxph.md)
  : Pool a prespecified Cox proportional hazards analysis of completed
  datasets
- [`catmech_compare_strategies()`](https://patrickiben.github.io/catmech/reference/catmech_compare_strategies.md)
  : Compare audit analysis strategies
- [`plot(`*`<catmech_pooled_glm>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_pooled_glm.md)
  : Plot pooled GLM coefficients (Forest Plot)
- [`plot(`*`<catmech_pooled_coxph>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_pooled_coxph.md)
  : Plot pooled Cox PH hazard ratios (Forest Plot)
- [`plot(`*`<catmech_pooled_lm>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_pooled_lm.md)
  : Plot pooled linear model coefficients (Forest Plot)
- [`plot(`*`<catmech_strategy_comparison>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_strategy_comparison.md)
  : Plot audit strategy comparison (Forest Plot)
- [`tidy()`](https://patrickiben.github.io/catmech/reference/tidy.md) :
  Turn a pooled catmech model into a tidy data frame
- [`glance()`](https://patrickiben.github.io/catmech/reference/glance.md)
  : Summarize a pooled catmech model with glance statistics

## Model Diagnostics & Calibration Checks

Evaluate holdout predictive accuracy and error distributions.

- [`catmech_holdout_check()`](https://patrickiben.github.io/catmech/reference/catmech_holdout_check.md)
  : Score an imputation model on reference rows it did not see
- [`catmech_recording_audit()`](https://patrickiben.github.io/catmech/reference/catmech_recording_audit.md)
  : Summarize what reference measurements reveal about recording errors
- [`catmech_sensitivity()`](https://patrickiben.github.io/catmech/reference/catmech_sensitivity.md)
  : Inspect a finite set of prespecified imputation assumptions
- [`catmech_marginal_summary()`](https://patrickiben.github.io/catmech/reference/catmech_marginal_summary.md)
  : Summarize population marginal coefficients from a joint covariate
  model
- [`catmech_joint_summary()`](https://patrickiben.github.io/catmech/reference/catmech_joint_summary.md)
  : Summarize the population coefficients of a fitted joint cluster
  model
- [`catmech_marginal_sensitivity_grid()`](https://patrickiben.github.io/catmech/reference/catmech_marginal_sensitivity_grid.md)
  : Canonical prior sensitivity grid for joint marginal cluster models
- [`plot(`*`<catmech_holdout_check>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_holdout_check.md)
  : Plot holdout check calibration diagnostics
- [`plot(`*`<catmech_recording_audit>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_recording_audit.md)
  : Plot a recording audit summary
- [`plot(`*`<catmech_marginal_sensitivity_grid>`*`)`](https://patrickiben.github.io/catmech/reference/plot.catmech_marginal_sensitivity_grid.md)
  : Plot a marginal prior sensitivity grid

## Regulatory Compliance & Provenance

Generate SHA-256 verified provenance reports and audit trails for
regulatory filings.

- [`catmech_audit_trail()`](https://patrickiben.github.io/catmech/reference/catmech_audit_trail.md)
  : Generate an auditable provenance and integrity trail for a catmech
  analysis
- [`catmech_hash()`](https://patrickiben.github.io/catmech/reference/catmech_hash.md)
  : Compute a canonical assumption hash

## Simulation Tools

Generate synthetic multi-center laboratory and trial datasets.

- [`simulate_catmech_lab()`](https://patrickiben.github.io/catmech/reference/simulate_catmech_lab.md)
  : Simulate measurements, quality flags, and interval categories

## S3 Print, Summary & Package Overview

Package overview and S3 formatting methods.

- [`catmech`](https://patrickiben.github.io/catmech/reference/catmech-package.md)
  [`catmech-package`](https://patrickiben.github.io/catmech/reference/catmech-package.md)
  : catmech: Auditable Model-Based Imputation and Category Sensitivity
- [`print(`*`<catmech_audit_plan>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_audit_plan.md)
  : Print a catmech audit size plan
- [`print(`*`<catmech_audit_power>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_audit_power.md)
  : Print an audit power calculation
- [`print(`*`<catmech_audit_trail>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_audit_trail.md)
  : Print an audit trail
- [`print(`*`<catmech_holdout_check>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_holdout_check.md)
  : Print a catmech holdout check result
- [`print(`*`<catmech_marginal_sensitivity_grid>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_marginal_sensitivity_grid.md)
  : Print a marginal sensitivity grid result
- [`print(`*`<catmech_model_diagnostics>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_model_diagnostics.md)
  : Print model diagnostics
- [`print(`*`<catmech_pooled_coxph>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_pooled_coxph.md)
  : Print a pooled Cox PH result
- [`print(`*`<catmech_pooled_glm>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_pooled_glm.md)
  : Print a pooled GLM result
- [`print(`*`<catmech_pooled_lm>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_pooled_lm.md)
  : Print a pooled linear model result
- [`print(`*`<catmech_recording_audit>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_recording_audit.md)
  : Print a catmech recording audit result
- [`print(`*`<catmech_result>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_result.md)
  : Print a catmech result
- [`print(`*`<catmech_sensitivity>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_sensitivity.md)
  : Print a catmech sensitivity analysis result
- [`print(`*`<catmech_strategy_comparison>`*`)`](https://patrickiben.github.io/catmech/reference/print.catmech_strategy_comparison.md)
  : Print strategy comparison results
- [`print(`*`<summary.catmech_audit_plan>`*`)`](https://patrickiben.github.io/catmech/reference/print.summary.catmech_audit_plan.md)
  : Print summary of a catmech audit plan
- [`print(`*`<summary.catmech_result>`*`)`](https://patrickiben.github.io/catmech/reference/print.summary.catmech_result.md)
  : Print a summary of a catmech result
- [`summary(`*`<catmech_audit_plan>`*`)`](https://patrickiben.github.io/catmech/reference/summary.catmech_audit_plan.md)
  : Summary of a catmech audit plan
- [`summary(`*`<catmech_result>`*`)`](https://patrickiben.github.io/catmech/reference/summary.catmech_result.md)
  : Summary of a catmech result
