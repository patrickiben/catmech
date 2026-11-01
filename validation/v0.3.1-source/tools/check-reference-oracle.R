#!/usr/bin/env Rscript
# Independent checks of the reference-data contract and normal predictive law.
# Run from the repository root after the reference implementation is stable.
# This does not establish repeated-sampling coverage or clinical validity.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "validation/reference-oracle"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite a nonempty evidence directory: ", out)
stopifnot(file.exists("DESCRIPTION"), file.exists("R/reference.R"))
pkgload::load_all(".", quiet = TRUE)

source_files <- c(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
                  "tools/check-reference-oracle.R")
hash_files <- function(paths) vapply(paths, function(p)
  digest::digest(file = p, algo = "sha256"), character(1))
source_before <- hash_files(source_files)

g <- c(-1.4, -.3, -.8, .5, -.1, .2, .9, .4,
       .2, 1.1, .5, 1.5, 1.2, .8, 1.9, 1.4)
x_ref <- seq(-2, 2, length.out = length(g))
x_missing <- c(-1.5, -.3, .8, 1.7)
n_ref <- length(g)
missing_rows <- n_ref + seq_along(x_missing)
d <- data.frame(recorded = rep(777, n_ref + length(x_missing)),
                reference = c(g, rep(NA_real_, length(x_missing))),
                outcome = c(x_ref, x_missing),
                preserved = paste0("row", seq_len(n_ref + length(x_missing))))

# Direct linear algebra is independent of catmech's fit/draw helpers.
A <- cbind(1, x_ref)
A0 <- cbind(1, x_missing)
V <- solve(crossprod(A))
bhat <- V %*% crossprod(A, g)
nu <- length(g) - ncol(A)
SSE <- sum((g - as.vector(A %*% bhat))^2)
expected_mean <- as.vector(A0 %*% bhat)
expected_cov <- SSE / (nu - 2) * (diag(nrow(A0)) + A0 %*% V %*% t(A0))
stopifnot(nu > 4, SSE > 0)

m <- 12000L
fit <- catmech_reference(d, "recorded", "reference", "outcome",
                         method = "normal", m = m, seed = 9042026L)
draws <- t(vapply(fit$completed, function(z) z$recorded[missing_rows],
                  numeric(length(missing_rows))))
centered <- sweep(draws, 2, expected_mean, "-")
moment_rows <- list()
add_moment <- function(label, simulated, reference) {
  estimate <- mean(simulated)
  mcse <- stats::sd(simulated) / sqrt(length(simulated))
  tolerance <- 5 * mcse + 1e-12
  moment_rows[[length(moment_rows) + 1L]] <<- data.frame(
    quantity = label, analytic = reference, simulation = estimate,
    mcse = mcse, tolerance = tolerance,
    standardized_difference = (estimate - reference) / mcse,
    pass = abs(estimate - reference) <= tolerance)
}
for (j in seq_along(expected_mean))
  add_moment(paste0("mean_missing_", j), draws[, j], expected_mean[j])
for (j in seq_along(expected_mean)) for (k in seq_len(j))
  add_moment(paste0("covariance_missing_", j, "_", k),
             centered[, j] * centered[, k], expected_cov[j, k])
moments <- do.call(rbind, moment_rows)

controls <- list()
control <- function(name, pass, detail) {
  controls[[length(controls) + 1L]] <<- data.frame(
    control = name, pass = isTRUE(pass), detail = detail)
}
values <- function(z) lapply(z$completed, function(x) x$recorded)
small <- catmech_reference(d, "recorded", "reference", "outcome",
                           method = "normal", m = 5L, seed = 93L)
changed <- d
changed$recorded <- c(rep(-1e8, n_ref), NA_real_, 1e8, -1e8, NA_real_)
changed_fit <- catmech_reference(changed, "recorded", "reference", "outcome",
                                 method = "normal", m = 5L, seed = 93L)
control("unused_raw_measurement_invariance",
        identical(values(small), values(changed_fit)),
        "Changing recorded targets, including making them missing, cannot alter completions at a fixed seed when they are not selected predictors.")
control("reference_anchors_exact",
        all(vapply(fit$completed, function(z)
          identical(z$recorded[seq_len(n_ref)], g), logical(1))),
        "Observed reference measurements overwrite conflicting recorded values exactly.")
control("non_target_columns_preserved",
        all(vapply(small$completed, function(z)
          identical(z[setdiff(names(z), "recorded")],
                    d[setdiff(names(d), "recorded")]), logical(1))),
        "Reference values, selected predictors, and unrelated columns are unchanged.")
control("original_data_preserved", identical(small$original, d),
        "The original object remains the supplied recorded/reference data.")
control("mids_only_reference_observed",
        identical(small$mids$data$recorded, d$reference),
        "The mids input contains reference targets and NA elsewhere, not recorded pseudo-anchors.")

pmm <- catmech_reference(d, "recorded", "reference", "outcome",
                         method = "pmm", donors = 5L, m = 10L, seed = 94L)
pmm_changed <- catmech_reference(changed, "recorded", "reference", "outcome",
                                 method = "pmm", donors = 5L, m = 10L, seed = 94L)
control("pmm_unused_raw_measurement_invariance",
        identical(values(pmm), values(pmm_changed)),
        "PMM cannot use unvalidated recorded targets as donors or fitting responses.")
control("pmm_donors_are_reference_values",
        all(vapply(pmm$completed, function(z)
          all(z$recorded[missing_rows] %in% g), logical(1))),
        "Every unavailable PMM target is drawn from an actual reference value.")
duplicate_data <- data.frame(recorded = rep(999, 8),
  reference = c(0, 0, 1, 1, 2, 2, NA_real_, NA_real_),
  outcome = c(0, 3, 1, 4, 2, 5, 1.5, 3.5))
duplicate_fit <- catmech_reference(duplicate_data, "recorded", "reference", "outcome",
                                  method = "pmm", m = 4L, seed = 96L)
duplicate_audit <- duplicate_fit$audit[!duplicate_fit$audit$reference_observed, ]
control("duplicate_pmm_donor_identity_not_invented",
        all(is.na(duplicate_audit$donor_row)) &&
          all(duplicate_audit$donor_value_match_count == 2L) &&
          !any(duplicate_audit$donor_row_identified),
        "Duplicated reference values do not reveal which donor row mice selected; the audit admits that ambiguity.")

full <- d
full$reference[missing_rows] <- c(-.6, .1, .7, 1.3)
full_fit <- catmech_reference(full, "recorded", "reference", "outcome",
                             method = "normal", m = 3L, seed = 95L)
control("all_reference_deterministic",
        all(vapply(full_fit$completed, function(z)
          identical(z$recorded, full$reference), logical(1))),
        "With all rows referenced, every completed target is exactly the reference vector.")
gold <- full
gold$recorded <- gold$reference
gold_lm <- stats::lm(outcome ~ recorded, data = gold)
pooled <- catmech_pool_lm(full_fit, outcome ~ recorded)
gold_coeff <- summary(gold_lm)$coefficients
gold_ci <- stats::confint(gold_lm)
idx <- match(pooled$term, rownames(gold_coeff))
control("all_reference_matches_gold_lm",
        isTRUE(all.equal(pooled$estimate, unname(gold_coeff[idx, "Estimate"]),
                         tolerance = 1e-12)) &&
        isTRUE(all.equal(pooled$std.error, unname(gold_coeff[idx, "Std. Error"]),
                         tolerance = 1e-12)) &&
        isTRUE(all.equal(pooled$lower, unname(gold_ci[idx, 1]), tolerance = 1e-12)) &&
        isTRUE(all.equal(pooled$upper, unname(gold_ci[idx, 2]), tolerance = 1e-12)),
        "The explicitly all-reference case agrees with ordinary gold-data lm estimates, standard errors, and confidence limits.")

# A deliberately wrong implementation that retained raw values would be detected.
control("known_wrong_raw_retention_detected",
        !identical(d$recorded[seq_len(n_ref)], g) &&
          !identical(d$recorded[missing_rows],
                     small$completed[[1L]]$recorded[missing_rows]),
        "A retain-recorded negative control disagrees with both anchored and imputed targets.")
cross_check <- moments[moments$quantity == "covariance_missing_4_3", ]
variance_check <- moments[moments$quantity == "covariance_missing_4_4", ]
negative <- data.frame(
  quantity = c("fixed_parameter_zero_cross_covariance", "plug_in_residual_variance"),
  wrong_answer = c(0, SSE / nu),
  analytic = c(expected_cov[4, 3], expected_cov[4, 4]),
  tolerance = c(cross_check$tolerance, variance_check$tolerance))
negative$detected <- abs(negative$wrong_answer - negative$analytic) > negative$tolerance
control("known_wrong_parameter_uncertainty_detected", all(negative$detected),
        "Prespecified zero cross-record covariance and plug-in residual variance answers fail the independent moment tolerances.")
control_table <- do.call(rbind, controls)

source_after <- hash_files(source_files)
if (!identical(source_before, source_after))
  stop("Source changed during the oracle run; no evidence written.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(moments, file.path(out, "posterior-predictive-moments.csv"), row.names = FALSE)
utils::write.csv(control_table, file.path(out, "data-contract-controls.csv"), row.names = FALSE)
utils::write.csv(negative, file.path(out, "negative-controls.csv"), row.names = FALSE)
utils::write.csv(data.frame(file = source_files, sha256 = unname(source_before)),
                 file.path(out, "source-hashes.csv"), row.names = FALSE)
utils::write.csv(data.frame(row = seq_along(expected_mean), mean = expected_mean,
                           expected_cov), file.path(out, "analytic-target.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
passed <- all(moments$pass) && all(control_table$pass)
writeLines(c(
  "# Independent reference-model checks", "",
  paste0("Overall result: **", if (passed) "PASS" else "FAIL", "**."), "",
  paste0("The normal path used ", m, " independent completions of a fixed small dataset."),
  "Direct linear algebra supplies the multivariate Student posterior predictive mean and covariance.",
  "Each moment tolerance was prespecified as five estimated Monte Carlo standard errors plus 1e-12.",
  "Covariances are estimated by products centered at the analytic mean, including cross-record terms.",
  "No catmech fitting, posterior, or interval helper is used to calculate the analytic answer.", "",
  "Data-contract checks test unused-recorded-value invariance, exact references, preservation,",
  "PMM reference-donor membership, and the fully referenced ordinary-lm endpoint.",
  "Source hashes were identical before and after the run.", "",
  "These checks establish selected implementation properties. They do not establish repeated-sampling",
  "coverage, correct reference selection, true gold-standard status, or clinical validity."
), file.path(out, "README.md"))
print(moments)
print(control_table[, c("control", "pass")])
if (!passed) stop("One or more independent reference checks failed; evidence retained.")
