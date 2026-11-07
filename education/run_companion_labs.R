#!/usr/bin/env Rscript
# ==============================================================================
# catmech Educational Companion — Executable R Lab & Refusal Verification Suite
# ==============================================================================
# Executes the live R exercises coupled to the 30 Companion units (C01..C30)
# and the 4 "When the Package Says No" Refusal & Boundary Labs, verifying exact
# numerical correspondence with the teaching answers in education/COMPANION.md.
# ==============================================================================

suppressPackageStartupMessages({
  library(catmech)
  library(survival)
})

cat("=== catmech Educational Companion: Executable R Lab Suite ===\n")

assert_true <- function(cond, label) {
  if (!isTRUE(cond)) {
    stop(sprintf("FAIL [%s]", label), call. = FALSE)
  }
  cat(sprintf("  [PASS] %s\n", label))
}

# ------------------------------------------------------------------------------
# Part 1: Analytical & Formula Verification Labs (Units C01 - C14)
# ------------------------------------------------------------------------------
cat("\n--- Part 1: Analytical & Formula Verification Labs (C01..C14) ---\n")

# Lab C01: RMSE reduction from retaining recorded values (ClinicalTrials.gov & NHANES)
rmse_drop_ctgov <- round(100 * (1 - 0.0286 / 0.0539), 1)
rmse_drop_nhanes <- round(100 * (1 - 0.0409 / 0.0785), 1)
assert_true(rmse_drop_ctgov == 46.9 && rmse_drop_nhanes == 47.9,
            "Lab C01: Retaining recorded values cuts RMSE by 46.9% (CT.gov) and 47.9% (NHANES)")

# Lab C02: Classic attenuation factor under additive error (lambda = sigma_T^2 / (sigma_T^2 + sigma_e^2))
sigma_T2 <- 4.0
sigma_e2 <- 1.0
attenuation <- sigma_T2 / (sigma_T2 + sigma_e2)
assert_true(abs(attenuation - 0.80) < 1e-12,
            "Lab C02: Classical errors-in-variables attenuation factor 4/(4+1) = 0.80")

# Lab C03: Phase I recorded measurement W_i under delta_i in {0, 1}
T_i <- 10.0
eps_i <- 3.5
W_clean <- (1 - 0) * T_i + 0 * (T_i + eps_i)
W_corrupt <- (1 - 1) * T_i + 1 * (T_i + eps_i)
assert_true(W_clean == 10.0 && W_corrupt == 13.5,
            "Lab C03: Two-phase contamination equation W_i = (1-delta_i)T_i + delta_i(T_i + eps_i)")

# Lab C05: Marginalized posterior corruption probability p_i via Gaussian convolution
pi_prior <- 0.20
w_obs <- 14.0
mu_pred <- 10.0
mu_e <- 4.0
sT2 <- 1.0
se2 <- 3.0
num <- pi_prior * dnorm(w_obs, mean = mu_pred + mu_e, sd = sqrt(sT2 + se2))
den <- num + (1 - pi_prior) * dnorm(w_obs, mean = mu_pred, sd = sqrt(sT2))
p_i <- num / den
assert_true(round(p_i, 4) == 0.9973,
            "Lab C05: Marginalized corruption probability p_i = 0.9973 for gross discrepancy")

# Lab C06: Precision-weighted conjugate Gaussian posterior for T_i | delta_i = 1
var_post <- 1 / (1 / sT2 + 1 / se2)
mean_post <- var_post * (mu_pred / sT2 + (w_obs - mu_e) / se2)
assert_true(abs(var_post - 0.75) < 1e-12 && abs(mean_post - 10.0) < 1e-12,
            "Lab C06: Conjugate posterior variance (1/1 + 1/3)^(-1) = 0.75 and mean = 10.0")

# Lab C08: Two-phase precision ratio rho(f) = sqrt(1 + k * (1 - f) / f)
k_val <- 0.27
f_val <- 0.25
rho_f <- sqrt(1 + k_val * (1 - f_val) / f_val)
assert_true(round(rho_f, 4) == 1.3454,
            "Lab C08: Two-phase precision ratio sqrt(1 + 0.27 * 3) = 1.3454 at f = 0.25")

# Lab C09: Prospective audit power via catmech_audit_power()
pwr <- catmech_audit_power(
  N = 1000, target_ratio = 1.25, power_target = 0.80,
  error_rate = 0.20, correlation = 0.85,
  fractions = c(0.10, 0.25, 0.40), n_sim = 500L, seed = 42L
)
assert_true(inherits(pwr, "catmech_audit_power") && pwr$table$power[3] > pwr$table$power[1],
            "Lab C09: catmech_audit_power() monotonic power curve across audit fractions")

# Lab C10 & C11: Rubin's combining rules and Barnard-Rubin finite-sample df
m <- 5
u_bar <- 0.04
b_var <- 0.01
t_var <- u_bar + (1 + 1 / m) * b_var
lambda_hat <- ((1 + 1 / m) * b_var) / t_var
nu_m <- (m - 1) / (lambda_hat^2)
nu_com <- 9 # e.g., C - 1 = 9 for C = 10 clusters
nu_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda_hat)
nu_br <- (nu_m * nu_obs) / (nu_m + nu_obs)
assert_true(abs(t_var - 0.052) < 1e-12 && round(lambda_hat, 4) == 0.2308 && round(nu_br, 4) == 5.3577,
            "Lab C10-C11: Rubin total variance T = 0.052, lambda = 0.2308, Barnard-Rubin df = 5.3577")

# ------------------------------------------------------------------------------
# Part 2: End-to-End Package Workflow Labs (Units C15 - C28)
# ------------------------------------------------------------------------------
cat("\n--- Part 2: End-to-End Package Workflow Labs (C15..C28) ---\n")

set.seed(42)
cohort <- simulate_catmech_lab(n = 240, seed = 42)
cohort$site <- rep(seq_len(8), length.out = nrow(cohort))
audited <- sample.int(nrow(cohort), size = 60)
cohort$lab_ref <- ifelse(seq_len(nrow(cohort)) %in% audited, cohort$lab_true, NA_real_)
cohort$adverse_event <- as.integer(cohort$outcome > median(cohort$outcome))
cohort$time <- rexp(nrow(cohort), rate = exp(0.20 * (cohort$lab_true - 6)))
cohort$event <- rbinom(nrow(cohort), 1, 0.80)

# Lab C15: Discrepancy screening with catmech_recording_audit()
rec_audit <- catmech_recording_audit(
  data = cohort, value = "lab_value", reference = "lab_ref",
  by = "site", agreement_tol = 0.5
)
assert_true(inherits(rec_audit, "catmech_recording_audit") && nrow(rec_audit$by) == 8L,
            "Lab C15: catmech_recording_audit() computes Jeffreys intervals across 8 sites")

# Lab C16: Linear reference MI with use_recorded = TRUE (cole_mi) vs FALSE (reference_normal)
fit_rec <- catmech_reference(
  data = cohort, value = "lab_value", reference = "lab_ref",
  predictors = c("age", "sex", "inflammation", "outcome"),
  use_recorded = TRUE, m = 5, seed = 101
)
fit_disc <- catmech_reference(
  data = cohort, value = "lab_value", reference = "lab_ref",
  predictors = c("age", "sex", "inflammation", "outcome"),
  use_recorded = FALSE, m = 5, seed = 101
)
pool_rec <- catmech_pool_lm(fit_rec, outcome ~ lab_value + age + sex, cluster = "site")
pool_disc <- catmech_pool_lm(fit_disc, outcome ~ lab_value + age + sex, cluster = "site")
se_rec <- pool_rec$std.error[pool_rec$term == "lab_value"]
se_disc <- pool_disc$std.error[pool_disc$term == "lab_value"]
assert_true(se_rec < se_disc,
            sprintf("Lab C16: use_recorded=TRUE lowers pooled SE (%.4f vs %.4f for use_recorded=FALSE)", se_rec, se_disc))

# Lab C17 & C18: Multi-endpoint pooling (LM with CV3, GLM with CR1S, Cox PH with Lin-Wei) + broom generics
pool_glm_res <- catmech_pool_glm(
  fit_rec, adverse_event ~ lab_value + age + sex,
  family = binomial(link = "logit"), cluster = "site", exponentiate = TRUE
)
pool_cox_res <- catmech_pool_coxph(
  fit_rec, Surv(time, event) ~ lab_value + age + sex,
  cluster = "site", exponentiate = TRUE
)
tidy_cox <- tidy(pool_cox_res)
glance_cox <- glance(pool_cox_res)
assert_true(all(pool_cox_res$df < 7.0) && nrow(tidy_cox) == 3L && glance_cox$m == 5L,
            "Lab C17-C18: Clustered Cox PH pooling bounds Barnard-Rubin df < C-1 = 7 and supports tidy()/glance()")

# Lab C26: Side-by-side strategy comparison via catmech_compare_strategies()
comp_res <- catmech_compare_strategies(
  data = cohort, value = "lab_value", reference = "lab_ref",
  formula = outcome ~ lab_value + age + sex,
  m = 5, seed = 101
)
assert_true(inherits(comp_res, "catmech_strategy_comparison") &&
            nrow(comp_res[comp_res$term == "lab_value", ]) == 4L,
            "Lab C26: catmech_compare_strategies() compares all 4 audit strategies side by side")

# Lab C28: SHA-256 provenance record via catmech_audit_trail()
trail_md <- catmech_audit_trail(fit_rec, pooled = pool_cox_res, format = "markdown")
assert_true(grepl("# Catmech Audit Trail & Provenance Report", as.character(trail_md)) &&
            grepl("### Audit Trail & Provenance Record", as.character(trail_md)),
            "Lab C28: catmech_audit_trail() emits SHA-256 provenance report with honest non-certification heading")

# ------------------------------------------------------------------------------
# Part 3: "When the Package Says No" — Refusal & Boundary Labs
# ------------------------------------------------------------------------------
cat("\n--- Part 3: 'When the Package Says No' — Refusal & Boundary Labs ---\n")

# Refusal Lab 1: Pooling refuses non-converged MCMC objects
bad_fit <- fit_rec
bad_fit$model_diagnostics <- list(status = "failed", convergence_pass = FALSE)
refused_lm <- tryCatch({
  catmech_pool_lm(bad_fit, outcome ~ lab_value + age)
  FALSE
}, error = function(e) {
  grepl("convergence", conditionMessage(e), ignore.case = TRUE)
})
assert_true(refused_lm,
            "Refusal Lab 1: catmech_pool_lm() visibly refuses to pool when convergence_pass = FALSE")

# Refusal Lab 2: Clustered predictor imputation warning in catmech_pool_lm()
cohort_complete_rec <- cohort
cohort_complete_rec$lab_value[is.na(cohort_complete_rec$lab_value)] <- mean(cohort$lab_value, na.rm = TRUE)
fit_clust_pred <- suppressWarnings(catmech_reference(
  data = cohort_complete_rec, value = "lab_value", reference = "lab_ref",
  predictors = c("age", "sex", "outcome"),
  cluster = "site", use_recorded = TRUE, m = 3, seed = 99
))
warned_clust_pred <- tryCatch({
  catmech_pool_lm(fit_clust_pred, outcome ~ lab_value + age, cluster = "site")
  FALSE
}, warning = function(w) {
  grepl("predictor|attenuat|cluster", conditionMessage(w), ignore.case = TRUE)
})
assert_true(warned_clust_pred,
            "Refusal Lab 2: catmech_pool_lm() warns when clustered random-intercept MI is used for a predictor")

# Refusal Lab 3: Undeclared reference assay noise warning in catmech_contamination()
noisy_cohort <- cohort
noisy_cohort$lab_ref <- ifelse(!is.na(cohort$lab_ref), cohort$lab_ref + rnorm(nrow(cohort), 0, 0.35), NA_real_)
warned_noisy_ref <- tryCatch({
  catmech_contamination(
    data = noisy_cohort, value = "lab_value", reference = "lab_ref",
    predictors = c("age", "sex", "outcome"),
    m = 2, iterations = 200, burnin = 100, chains = 4, seed = 55
  )
  FALSE
}, warning = function(w) {
  grepl("disagree|reference_sd|agreement_tol", conditionMessage(w), ignore.case = TRUE)
})
assert_true(warned_noisy_ref,
            "Refusal Lab 3: catmech_contamination() warns when audited reference values universally disagree at reference_sd = 0")

# Refusal Lab 4: Sensitivity analysis refuses unlabeled or empty scenario list
refused_sens <- tryCatch({
  catmech_sensitivity(list(), outcome ~ lab_value + age, term = "lab_value")
  FALSE
}, error = function(e) TRUE)
assert_true(refused_sens,
            "Refusal Lab 4: catmech_sensitivity() refuses empty scenario specifications")

cat("\nALL 16 EXECUTABLE R COMPANION LABS & REFUSAL CHECKS PASSED.\n")
