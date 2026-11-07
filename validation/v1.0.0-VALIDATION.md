# Validation record: catmech 1.0.0

Date: 2026-11-05. Official production release graduating catmech from experimental
status to a frozen, auditable research and clinical toolkit. The [v0.10.1 record](v0.10.1-VALIDATION.md)
preserves the numerical convergence extension; all earlier version records (v0.1.0 through v0.10.0)
remain permanently frozen.

## Capabilities and interface freeze in v1.0.0

Version 1.0.0 expands downstream modeling, audit design, and regulatory compliance verification:
- **Multi-endpoint pooling & tidygenerics**: `catmech_pool_glm()` implements Barnard-Rubin pooling for generalized
  linear models (logistic, Poisson), and `catmech_pool_coxph()` supports time-to-event survival pooling
  using event-based residual degrees of freedom, with `tidy()`, `glance()`, and `catmech_compare_strategies()`.
- **Diagnostic & forest visualization**: base R S3 `plot()` methods provide visual inspection for holdout predictive
  calibration, recording audits, audit planning curves, prior sensitivity grids, pooled forest plots, and strategy comparisons.
- **Prospective audit planning**: `catmech_audit_power()` computes required audit sample sizes and statistical
  power empirically or under prospective assumptions.
- **Regulatory compliance trail**: `catmech_audit_trail()` generates SHA-256 verified provenance reports
  (Markdown, JSON, HTML) capturing inputs, seeds, options, convergence diagnostics, and pooled estimates.
- **End-to-end vignette**: `vignettes/catmech-clinical-audit.Rmd` demonstrates a full clinical audit workflow.

## Verification suite

The release was subjected to complete verification:
- **Unit test suite**: 1,754 testthat unit tests across 20 test files pass with 0 errors and 0 failures.
- **R CMD check**: Passes `R CMD check --no-manual catmech_1.0.0.tar.gz` with `Status: OK` (0 errors, 0 warnings, 0 notes).
- **Benchmark artifacts**: `tools/verify-artifacts.py` passes all 10 benchmark studies (synthetic contamination,
  real-data ClinicalTrials.gov and NHANES, Rwanda HMIS facility audit, clustered audits, variance checks,
  outside methods, and learning pages).
- **Fact-checking**: `docs/case-study/fact_check.py` recomputes 118 claims against frozen data sources with 100% agreement,
  0 broken links, and 1-page PDF compliance.

Core functions are auditable and production-ready for declared workflows with recorded values and audit subsets.
