---
title: 'catmech: An Auditable R Package for Measurement Verification, Clinical Data Audits, and Model Pooling'
tags:
  - R
  - biostatistics
  - clinical trials
  - measurement error
  - multiple imputation
  - audit planning
  - survival analysis
authors:
  - name: Patrick Iben
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 5 November 2026
bibliography: paper.bib
---

# Summary

In clinical trials, electronic health record (EHR) registries, and observational epidemiological cohorts, continuous biological and clinical measurements (such as cardiac biomarkers, viral loads, or serum creatinine) are frequently compromised by unmeasured recording errors, unit conversion bugs, and missingness. While standard multiple imputation packages (such as `mice` and `amelia`) impute completely missing covariates under missing-at-random (MAR) assumptions, they cannot address scenarios where recorded values are present but untrusted, or where an independent, gold-standard reference measurement has been audited on a subsample of records.

`catmech` is an R package for auditable measurement verification, clinical data audits, and multi-endpoint model pooling. It provides:
1. Reference-calibrated contamination mixture models that estimate per-record corruption probabilities and impute verified true measurements without discarding correct recorded values.
2. Cross-validated holdout checks and diagnostic visualizations (randomized Probability Integral Transform uniform histograms) to evaluate predictive calibration.
3. Prospective power and sample size calculators (`catmech_audit_power`) alongside empirical audit planning curves (`catmech_audit_plan`) based on two-phase sampling theory.
4. Downstream regression pooling for linear models (`catmech_pool_lm`), generalized linear models (`catmech_pool_glm`), and Cox proportional hazards survival models (`catmech_pool_coxph`) using Rubin's rules with Barnard–Rubin small-sample degrees of freedom.
5. Cryptographically verified audit trail reports (`catmech_audit_trail`) providing SHA-256 data and completion checksums designed to meet regulatory compliance standards (FDA 21 CFR Part 11 / ICH GCP E9).

# Statement of Need

A common practice in clinical research is source-data verification (SDV), where monitors audit a random or targeted subset of medical records against primary source documents. When recording errors are detected, researchers face a difficult dilemma:
- *Trusting the recorded value* introduces bias and distorts treatment effect estimates.
- *Discarding all un-audited recorded values* sacrifices substantial statistical precision, inflating confidence intervals.
- *Standard multiple imputation* ignores the valuable information present in un-audited measurements that were actually recorded correctly.

`catmech` addresses this gap by formalizing the two-phase sampling relationship between fallible recorded measurements $W_i$, latent true values $T_i$, and gold-standard audit references $R_i$. The package models the contamination process:
$$W_i = T_i + \delta_i \epsilon_i, \quad \delta_i \sim \text{Bernoulli}(\pi), \quad \epsilon_i \sim \mathcal{N}(\mu_e, \sigma_e^2)$$
where $\delta_i = 1$ indicates a corrupted entry and $\pi$ is the stratum-specific corruption rate.

Crucially, `catmech` embeds strict diagnostic screens: if the MCMC sampler fails convergence or if predictive cross-validation indicates miscalibration on held-out reference rows, pooled inference is explicitly refused. This guarantees that downstream clinical conclusions cannot be drawn from invalid or unconverged statistical completions.

# Key Features

* **Multi-Endpoint Pooling**: Supports continuous (Gaussian), binary (logistic), count (Poisson), and time-to-event (Cox proportional hazards) endpoints with cluster-robust sandwich covariance estimators (Huber–White–Liang–Zeger and Lin–Wei).
* **Two-Phase Sample Size Planning**: Prospective power calculations determine the minimum number of medical records $n_{\text{audit}}$ required to achieve a target confidence interval precision ratio $\rho \le \rho_{\text{target}}$.
* **Auditable Integrity Trail**: Computes SHA-256 hashes of original data, mathematical assumptions, and all completed datasets, formatting provenance reports in Markdown, JSON, or self-contained HTML.

# Mathematics and Software Architecture

The package is implemented in R and C++17 via `Rcpp` and `RcppArmadillo`. Sampler convergence is monitored using Gelman–Rubin $\hat{R}$, bulk effective sample size (ESS), and tail ESS. Small-sample degrees of freedom are calculated using the Barnard & Rubin (1999) adjustment:
$$\nu_m = (m - 1)\left(1 + \frac{\bar{U}}{(1 + m^{-1})B}\right)^2, \quad \nu_{BR} = \frac{\nu_m \nu_{obs}}{\nu_m + \nu_{obs}}$$
ensuring calibrated frequentist coverage across finite clinical trial samples.

# Acknowledgements

The package is available on GitHub and is tested against an extensive battery of unit tests and reproducible simulation studies.
