# catmech v1.0.0: One-Page Clinical Audit Decision & Boundary Card

**Canonical Manuscript:** `paper/manuscript.Rmd` (`SHA-256: 8152288b1c2af2a4...`) | **Companion Visual Case Study:** `docs/case-study/catmech-case-study.pdf`

This educational companion explains the statistical definitions, benchmark archives, and R workflows in the frozen catmech manuscript; it does not establish universal interval coverage under unmodeled misspecification, validate per-record error probabilities for individual patient adjudication on diffuse real-world data, or certify FDA 21 CFR Part 11 or ICH E9 regulatory compliance.

## 1. The Empirical Headline (500,000+ Verified Requests Across 10 Frozen Archives)

| Benchmark Archive | Audit Design | Key Empirical Finding | Scope / Operational Boundary |
| :--- | :--- | :--- | :--- |
| **ClinicalTrials.gov** (`N = 27,208`) & **NHANES** (`N = 3,908`) | `33,800` requests (`f = 25%` primary) | Retaining recorded values (`use_recorded = TRUE`) cuts RMSE by **43% to 48%** vs discarding unaudited values (`0.0286` vs `0.0539`; `0.0409` vs `0.0785`). | Per-record corruption probability `p_i` has near-chance AUC (**0.543** on CT.gov, **0.500** on NHANES). Never use `p_i` to adjudicate individual rows. |
| **Synthetic Contamination** (`validation/contamination-study`) | `52,950` requests (`27` scenarios) | `catmech_contamination()` narrows CIs by **44.3%** vs complete-case (`0.1842` vs `0.3309`) and prevents `-0.2354` linear-MI bias under 10x gross outliers. | Fails `<90%` coverage in 3 stress cells: `heavy_tail_target` (**87.3%**), `nonignorable_missed` (**79.3%**), `noisy_undeclared` (**88.7%**). |
| **Clustered Variance & MI** (`validation/variance-mi-study`) | `148,500` requests (`11` Rwanda HMIS + sim) | Delete-one-cluster jackknife (`cv3`) beats `cr1s`/`cr2` (`92.0%-97.5%` recount coverage). | When audited `T_i` is a **predictor**, fit `catmech_reference(use_recorded = TRUE)` *without* `cluster` and pool with `catmech_pool_lm(cluster = )`. |
| **Real-Audit Slope Recovery** (`Rwanda HMIS` & `Psychology`) | `45,600` requests (`Rules 1-3`) | `Linear MI` cuts RMSE to `0.458x` on size-dependent `ANC 4th visits`; Rule 2 (`Contamination` vs `Linear MI`) is **INCONCLUSIVE**. | `Contamination` had lower RMSE on 2 of 3 Rwanda slopes but **35.4% to 71.6%** clustered coverage and `6%` to `19%` MCMC refusal rates. |

## 2. Two-Step R Workflow & Decision Rule

1. **Screen & Compare:** Run `catmech_recording_audit(data, value, reference, by = 'site')` and `catmech_compare_strategies(data, value, reference, formula)`.
2. **Select Imputation & Pooling Engine:**
   - **Default (Subtle, diffuse, or moderate real-world discrepancies, or heavy-tailed target):** `fit <- catmech_reference(data, value, reference, predictors, use_recorded = TRUE, m = 20)`
   - **Gross, separable outliers only (with Gaussian target & ignorable audit selection):** `fit <- catmech_contamination(...)`, verifying `>= 90%` return rate and positive holdout Brier/RMSE skill in `catmech_holdout_check(fit)`.
   - **Pool with Barnard-Rubin df (`nu_BR < nu_com`) & Provenance:** `catmech_pool_lm(fit, formula, cluster)` (`cv3`), `catmech_pool_glm(...)` (`CR1S`), `catmech_pool_coxph(...)` (`d - p` or `Lin-Wei`), and `catmech_audit_trail(fit, pooled)`.
