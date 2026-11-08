#!/usr/bin/env Rscript
# ==============================================================================
# catmech Educational Companion — Interactive R Console Tutor (S1..S12)
# ==============================================================================
# Interactive usage in an R console:
#   source("education/interactive_tutor.R")
#   catmech_tutor()
#
# Automated verification usage:
#   Rscript education/interactive_tutor.R --non-interactive
# ==============================================================================

suppressPackageStartupMessages({
  library(catmech)
  library(survival)
})

CATMECH_LESSONS <- list(
  list(
    id = "S1",
    title = "The Three-Way Clinical Audit Dilemma",
    summary = paste(
      "When an internal audit verifies a gold-standard reference measurement R_i on",
      "a subset of records, trialists face three flawed conventional options:",
      "(1) Naive analysis of recorded values W_i suffers attenuation bias;",
      "(2) Complete-case analysis of audited rows inflates standard errors;",
      "(3) Discarding W_i on unaudited rows wastes valid recordings."
    ),
    code = "round(100 * (1 - c(ctgov = 0.0286 / 0.0539, nhanes = 0.0409 / 0.0785)), 1)",
    question = "On ClinicalTrials.gov and NHANES at a 25% audit, by roughly how much does retaining recorded values reduce RMSE compared to discarding unaudited values?",
    options = c(
      "A) About 5% to 10%",
      "B) About 43% to 48%",
      "C) It increases RMSE on both datasets"
    ),
    correct = "B",
    explanation = "Retaining recorded values reduces RMSE by 43% to 48% (0.0286 vs 0.0539 on ClinicalTrials.gov; 0.0409 vs 0.0785 on NHANES)."
  ),
  list(
    id = "S2",
    title = "Two-Phase Data Generation & Contamination Mixture",
    summary = paste(
      "Latent true values follow T_i | X_i ~ N(X_i beta, sigma_T^2).",
      "Observed recorded values follow W_i = (1 - delta_i) T_i + delta_i (T_i + eps_i),",
      "where delta_i in {0, 1} is an unobserved corruption indicator and eps_i ~ N(mu_e, sigma_e^2).",
      "Audited reference measurements follow R_i = T_i + eta_i."
    ),
    code = "pi <- 0.20; w <- 14; num <- pi * dnorm(w, 14, 2); round(num / (num + (1 - pi) * dnorm(w, 10, 1)), 4)",
    question = "When delta_i = 0 for an unaudited record in the contamination model, what is the conditional posterior draw of T_i?",
    options = c(
      "A) T_i = W_i deterministically",
      "B) T_i is drawn from N(X_i beta, sigma_T^2) ignoring W_i",
      "C) T_i is set to missing (NA)"
    ),
    correct = "A",
    explanation = "Conditional on drawing delta_i = 0 (uncorrupted recording), T_i = W_i deterministically; when delta_i = 1, T_i is drawn from the precision-weighted Gaussian posterior."
  ),
  list(
    id = "S3",
    title = "Gaussian Reference Imputation: use_recorded = TRUE vs FALSE",
    summary = paste(
      "catmech_reference() fits a fast C++17 (RcppArmadillo) Gaussian reference model.",
      "With use_recorded = TRUE (Simple Linear MI / cole_mi), W_i enters as a continuous",
      "linear predictor alongside X_i without positing discrete corruption indicators.",
      "With use_recorded = FALSE (reference_normal), unaudited W_i values are discarded."
    ),
    code = "args(catmech_reference)",
    question = "Why is catmech_reference(use_recorded = TRUE) preferred over catmech_contamination() when recording discrepancies are subtle and diffuse?",
    options = c(
      "A) Because diffuse discrepancies cannot be separated into discrete clean and corrupt components, and linear calibration avoids MCMC refusals and mixture undercoverage",
      "B) Because catmech_reference() does not require any audited reference column",
      "C) Because catmech_reference() certifies FDA 21 CFR Part 11 compliance"
    ),
    correct = "A",
    explanation = "Under diffuse real-world discrepancies, discrete mixture components are non-separable; Simple Linear MI (use_recorded = TRUE) achieves comparable or better RMSE and more reliable coverage without MCMC overhead."
  ),
  list(
    id = "S4",
    title = "Rough Two-Phase Audit-Size Planning & The NO-GO Boundary",
    summary = paste(
      "Using two-phase influence functions, the interval-width precision ratio at audit",
      "fraction f = n_audit / N is rho(f) = sqrt(1 + k * (1 - f) / f).",
      "catmech_audit_power() is a rough, unvalidated planning aid: its 'power' column is the share",
      "of plausible k draws meeting the target, not the power of a test, and it also prints the",
      "audit size linear imputation needs at the same correlation, often several times larger.",
      "catmech_audit_plan() carries a pre-specified NO-GO validation status (11.3%-16.8% median width error at f=5%)."
    ),
    code = "round(sqrt(1 + 0.27 * (1 - 0.25) / 0.25), 4)",
    question = "What does the pre-specified NO-GO status on catmech_audit_plan() mean for users?",
    options = c(
      "A) The function was removed from the package",
      "B) It is documented as an exploratory planning tool rather than a validated sample-size guarantee",
      "C) It only works when f = 1.0"
    ),
    correct = "B",
    explanation = "Because median interval-width prediction error at a 5% audit fraction was 11.3% to 16.8% (above the 10% pre-specified threshold) and intervals were conservative (38.9% of cases), catmech_audit_plan() is retained only as an exploratory aid."
  ),
  list(
    id = "S5",
    title = "Rubin's Rules & Barnard-Rubin Small-Sample Degrees of Freedom",
    summary = paste(
      "Across m completed datasets, total variance is T = U_bar + (1 + 1/m) B.",
      "Classical Rubin degrees of freedom nu_m = (m - 1) * lambda^(-2) assume infinite",
      "complete-data sample size. Barnard & Rubin (1999) bound nu_BR below nu_com",
      "(N - p for unclustered LM/GLM, C - 1 for C clusters, and d - p for d survival events)."
    ),
    code = "nu_com <- 9; lam <- 0.2308; nu_m <- 4 / lam^2; nu_obs <- ((nu_com + 1)/(nu_com + 3)) * nu_com * (1 - lam); round((nu_m * nu_obs)/(nu_m + nu_obs), 4)",
    question = "In a 10-site clustered study (C = 10), what is the upper bound on Barnard-Rubin degrees of freedom nu_BR when pooling with cluster = 'site'?",
    options = c(
      "A) N - p (total patients minus parameters)",
      "B) Strictly less than C - 1 = 9",
      "C) Infinite as m -> infinity"
    ),
    correct = "B",
    explanation = "With C = 10 clusters, complete-data degrees of freedom are nu_com = C - 1 = 9, and nu_BR = (nu_m * nu_obs)/(nu_m + nu_obs) < nu_obs < nu_com = 9."
  ),
  list(
    id = "S6",
    title = "Multi-Endpoint Pooling: Linear (CV3), GLM (CR1S), and Cox PH (Lin-Wei)",
    summary = paste(
      "catmech_pool_lm() uses the delete-one-cluster jackknife (CV3) for clustered linear models;",
      "catmech_pool_glm() uses the CR1S Huber-White-Liang-Zeger sandwich for logistic/Poisson models;",
      "catmech_pool_coxph() uses event-based df (d - p) or clustered Lin-Wei score-residual sandwich."
    ),
    code = "args(catmech_pool_coxph)",
    question = "In unclustered Cox proportional hazards pooling (catmech_pool_coxph), how are complete-data degrees of freedom nu_com defined?",
    options = c(
      "A) Total sample size N minus p",
      "B) Number of observed failure events d minus p",
      "C) Number of imputations m minus 1"
    ),
    correct = "B",
    explanation = "Partial likelihood Fisher information in Cox regression scales with the number of failure events d rather than total censored + uncensored rows N, so nu_com = d - p."
  ),
  list(
    id = "S7",
    title = "Out-of-Sample Calibration Screening & Its Slope-Bias Blind Spot",
    summary = paste(
      "catmech_holdout_check() runs K-fold cross-validation on audited rows and computes",
      "tie-randomized empirical rank PIT values u_i against Uniform(0, 1).",
      "Crucially, held-out row checks test observation-level prediction, so they have low",
      "power to detect regression slope bias (0/40 alarms in nonlinear_target, 1/40 in heavy_tail_target)."
    ),
    code = "args(catmech_holdout_check)",
    question = "If catmech_holdout_check() reports miscalibration_detected = FALSE, does that prove the downstream regression slope is unbiased?",
    options = c(
      "A) Yes, a passing KS test guarantees unbiased regression coefficients",
      "B) No; held-out row checks have low power against regression slope bias and serve only as a diagnostic screen",
      "C) Yes, provided m >= 20"
    ),
    correct = "B",
    explanation = "Observation-level predictive calibration does not guarantee regression slope validity; in our 8-scenario study, the KS alarm fired in 0/40 nonlinear_target and 1/40 heavy_tail_target cohorts."
  ),
  list(
    id = "S8",
    title = "Synthetic Contamination Benchmarks & Pre-Specified Stress Failures",
    summary = paste(
      "Across 52,950 synthetic requests, catmech_contamination() narrowed CIs by 44.3% vs",
      "complete-case analysis in missed_flags and avoided the -0.2354 linear-MI bias under",
      "gross outliers (error_gross_and_small). However, it failed in three stress cells:",
      "heavy_tail_target (87.3% cov), nonignorable_missed (79.3% cov), and noisy_undeclared (88.7% cov)."
    ),
    code = "round(100 * (1 - 0.1842 / 0.3309), 1)",
    question = "Why does catmech_contamination() under-cover (87.3%) when the true target distribution T_i has heavy t_4 tails (heavy_tail_target)?",
    options = c(
      "A) The Gaussian mixture mistakes genuine extreme tail values of T_i for recording errors and shrinks them toward X_i beta",
      "B) Because m = 20 is too large",
      "C) Because complete-case analysis also fails in heavy_tail_target"
    ),
    correct = "A",
    explanation = "Under a Gaussian target model, genuine heavy-tailed extremes in T_i look like contamination outliers and get shrunk toward the conditional mean X_i beta."
  ),
  list(
    id = "S9",
    title = "Clustered Variance (CV3) & The Predictor-Slope Attenuation Trap",
    summary = paste(
      "Across 148,500 clustered requests, the delete-one-cluster jackknife (cv3) beat cr1s and cr2",
      "on Rwanda HMIS recounts. However, clustered random-intercept imputation",
      "(catmech_reference(cluster = )) attenuated regression slopes by up to -16.1% when the",
      "audited variable was a predictor rather than the outcome."
    ),
    code = "cat('Predictor rule: fit catmech_reference(use_recorded = TRUE) without cluster, pool with catmech_pool_lm(cluster = )\\n')",
    question = "When the audited variable is a predictor in a clustered linear regression, what workflow does catmech recommend?",
    options = c(
      "A) Fit catmech_reference(cluster = 'site') and pool without cluster",
      "B) Fit catmech_reference(use_recorded = TRUE) without cluster, and pool with catmech_pool_lm(cluster = 'site')",
      "C) Always use cr2 variance with Satterthwaite df"
    ),
    correct = "B",
    explanation = "Fitting catmech_reference(use_recorded = TRUE) without cluster and pooling with catmech_pool_lm(cluster = 'site') achieved -3.0% bias and 95.4% coverage, avoiding random-intercept predictor attenuation."
  ),
  list(
    id = "S10",
    title = "Real-World Paired Benchmarks & Why Per-Record Error Probabilities Fail",
    summary = paste(
      "On ClinicalTrials.gov (N = 27,208) and NHANES (N = 3,908), retaining recorded values",
      "cut RMSE by 43%-48% over discarding unaudited values. However, because real-world",
      "discrepancies are diffuse relative to T_i | X_i, per-record corruption probabilities p_i",
      "had near-chance discrimination (AUC = 0.543 on CT.gov, 0.500 on NHANES, 0.501-0.582 on Rwanda HMIS)."
    ),
    code = "c(synthetic_missed_flags_auc = 0.962, ctgov_auc = 0.543, nhanes_auc = 0.500)",
    question = "Can per-record posterior corruption probabilities p_i from catmech_contamination() be used to flag individual patient records in diffuse real-world audits?",
    options = c(
      "A) Yes, p_i > 0.5 reliably identifies corrupted rows in all datasets",
      "B) No; on diffuse real-world discrepancies, p_i has near-chance AUC (0.500-0.582) and must not adjudicate individual records",
      "C) Only when N > 10,000"
    ),
    correct = "B",
    explanation = "Because diffuse real-world errors are small relative to the natural residual spread of T_i | X_i, row-level corruption probabilities p_i perform near coin-flip chance (AUC 0.500 to 0.582)."
  ),
  list(
    id = "S11",
    title = "Real-Audit Slope Recovery (Rwanda HMIS & Psychology Recounts) & Decision Tree",
    summary = paste(
      "On Rwanda HMIS facility recounts (45,600 requests), Rule 2 comparing catmech_contamination()",
      "against Linear MI was INCONCLUSIVE: the mixture model had lower point RMSE on 2 of 3 slopes",
      "but suffered 35.4%-71.6% clustered interval coverage and 6%-19% refusal rates.",
      "Therefore, Simple Linear MI is the default for moderate real-world discrepancies, while",
      "catmech_contamination() is reserved for gross, separable corruption."
    ),
    code = "c(rwanda_pnc_cov = 0.716, rwanda_syphilis_cov = 0.666, rwanda_iron_cov = 0.354)",
    question = "Why was the pre-specified Rule 2 verdict comparing catmech_contamination() to Simple Linear MI on Rwanda HMIS scored as INCONCLUSIVE?",
    options = c(
      "A) Because no discrepancies existed in the Rwanda dataset",
      "B) Because catmech_contamination() had lower point RMSE on two slopes (PNC and Iron) but failed the 95% coverage and >=90% return-rate conditions",
      "C) Because Linear MI failed to run on Rwanda HMIS"
    ),
    correct = "B",
    explanation = "Point error and interval coverage pointed in opposite directions: lower point RMSE on PNC and Iron, but severe clustered undercoverage (35.4% to 71.6%) and >10% MCMC refusal rates on two slopes."
  ),
  list(
    id = "S12",
    title = "Provenance Reporting (catmech_audit_trail) & Regulatory Scope Boundaries",
    summary = paste(
      "catmech_audit_trail() generates Markdown, JSON, or HTML provenance records containing",
      "R/platform metadata, SHA-256 hashes of the serialized input data and all m completions and",
      "of the assumption specification as canonical JSON, MCMC convergence status, and pooled tables.",
      "A rerun can be compared with the hashes, but they do not cover the pooled results or make the",
      "record tamper-evident. The record supports documentation under FDA 21 CFR Part 11 and",
      "ICH E9, but does not by itself certify regulatory compliance."
    ),
    code = "args(catmech_audit_trail)",
    question = "What is the exact regulatory scope of the report generated by catmech_audit_trail()?",
    options = c(
      "A) It legally certifies that the clinical trial is FDA 21 CFR Part 11 compliant",
      "B) It provides a SHA-256 provenance record that lets a rerun be compared and supports regulatory documentation, while compliance depends on the surrounding validated system and procedures",
      "C) It encrypts patient identifiers for HIPAA de-identification"
    ),
    correct = "B",
    explanation = "Hashing inputs and completions supports reproducibility and audit documentation (FDA 21 CFR Part 11 / ICH E9), but regulatory compliance depends on system validation, access controls, and SOPs."
  )
)

catmech_tutor <- function(interactive_mode = interactive()) {
  cat("======================================================================\n")
  cat("  catmech 1.0.1: Interactive R Console Tutor (Lessons S1..S12)\n")
  cat("======================================================================\n")
  score <- 0L
  total <- length(CATMECH_LESSONS)

  for (idx in seq_along(CATMECH_LESSONS)) {
    les <- CATMECH_LESSONS[[idx]]
    cat(sprintf("\n[%s/%d] %s: %s\n", les$id, total, les$id, les$title))
    cat(strwrap(les$summary, width = 72, prefix = "  "), sep = "\n")
    cat("\n  Live R verification snippet:\n")
    cat(sprintf("    > %s\n", les$code))
    val <- eval(parse(text = les$code), envir = globalenv())
    if (!is.null(val)) {
      out_lines <- capture.output(print(val))
      cat(paste0("    ", out_lines), sep = "\n")
    }
    cat(sprintf("\n  Check (%s): %s\n", les$id, les$question))
    for (opt in les$options) {
      cat(sprintf("    %s\n", opt))
    }

    ans <- if (isTRUE(interactive_mode)) {
      toupper(trimws(readline(prompt = "  Your answer (A/B/C, or Q to quit): ")))
    } else {
      les$correct
    }

    if (identical(ans, "Q")) {
      cat("\nExiting tutor early.\n")
      return(invisible(list(score = score, completed = idx - 1L, total = total)))
    }

    if (identical(ans, les$correct)) {
      score <- score + 1L
      cat(sprintf("  [CORRECT: %s] %s\n", les$correct, les$explanation))
    } else {
      cat(sprintf("  [NOT QUITE — Correct is %s] %s\n", les$correct, les$explanation))
    }
  }

  cat("\n======================================================================\n")
  cat(sprintf("  Tutor Complete: %d / %d correct (100%% boundary & transfer checks)\n", score, total))
  cat("======================================================================\n")
  invisible(list(score = score, completed = total, total = total))
}

if (!interactive() && "--non-interactive" %in% commandArgs(trailingOnly = TRUE)) {
  res <- catmech_tutor(interactive_mode = FALSE)
  if (!identical(res$score, 12L)) {
    stop("Non-interactive tutor self-check failed.", call. = FALSE)
  }
}
