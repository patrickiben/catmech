contamination_fixture <- function(n = 240L, seed = 61L, every = 3L) {
  d <- simulate_catmech_lab(n, seed = seed, corruption_probability = .5)
  truth <- d$lab_true
  wrong <- d$corrupted
  d$gold <- NA_real_
  d$gold[seq(1L, n, by = every)] <- truth[seq(1L, n, by = every)]
  d[c("lab_true", "corrupted")] <- NULL
  list(data = d, truth = truth, wrong = wrong)
}
contamination_predictors <- c("age", "inflammation", "outcome")
contamination_fit <- function(d, ..., m = 4, seed = 5, predictors = contamination_predictors)
  catmech_contamination(d, "lab_value", "gold", predictors, corruption_strata = "outlier_flag", m = m, seed = seed, ...)

test_that("compound strata preserve ordinary labels and separate colliding tuples", {
  ordinary <- data.frame(flag = c("No", "Yes", "No", NA_character_),
                         site = c("B", "A", "A", "unused"))
  recorded <- c(TRUE, TRUE, TRUE, FALSE)
  old_keys <- paste(ordinary$flag, ordinary$site, sep = " | ")
  old_keys[!recorded] <- NA_character_
  old_levels <- sort(unique(old_keys[recorded]))
  expect_identical(catmech_contamination_strata(ordinary, c("flag", "site"), recorded, "test"),
    list(index = match(old_keys, old_levels), levels = old_levels))

  colliding <- data.frame(a = c("A | B", "A", "A | B"),
                          b = c("C", "B | C", "C"))
  mapped <- catmech_contamination_strata(colliding, c("a", "b"), rep(TRUE, 3), "test")
  expect_length(mapped$levels, 2L)
  expect_false(mapped$index[1] == mapped$index[2])
  expect_identical(mapped$index[1], mapped$index[3])
  expect_true(any(grepl('"a"="A | B", "b"="C"', mapped$levels, fixed = TRUE)))
  expect_true(any(grepl('"a"="A", "b"="B | C"', mapped$levels, fixed = TRUE)))
  permutation <- c(2L, 3L, 1L)
  reordered <- catmech_contamination_strata(colliding[permutation, ], c("a", "b"), rep(TRUE, 3), "test")
  expect_identical(reordered$levels, mapped$levels)
  expect_identical(reordered$index[order(permutation)], mapped$index)

  # Even a user label equal to a generated disambiguation label stays distinct.
  adversarial <- rbind(colliding[1:2, ],
    data.frame(a = "A", b = 'B | C ["a"="A | B", "b"="C"]'))
  unique_map <- catmech_contamination_strata(adversarial, c("a", "b"), rep(TRUE, 3), "test")
  expect_equal(length(unique(unique_map$index)), 3L)
  expect_false(anyDuplicated(unique_map$levels) > 0L)
})

test_that("compound label collisions cannot bypass calibration support", {
  d <- data.frame(w = c(1, 2, 3, 4), r = c(1, 2, NA_real_, NA_real_),
    a = c("A | B", "A | B", "A", "A"), b = c("C", "C", "B | C", "B | C"))
  expect_error(catmech_contamination(d, "w", "r", character(),
    corruption_strata = c("a", "b")), "No calibration rows support stratum")
  expect_error(catmech_contamination(d, "w", "r", character(),
    error_strata = c("a", "b")), "No calibration rows support stratum")
  d$r <- d$w
  complete <- catmech_contamination(d, "w", "r", character(),
    corruption_strata = c("a", "b"), m = 2L)
  expect_equal(length(unique(complete$contamination$corruption$corruption_stratum)), 2L)
  expect_identical(complete$completed[[1]]$w, d$r)
})

test_that("row terms equal brute-force integration over the target", {
  m <- 1.2; s2 <- .5^2; mu <- .7; t2 <- .3^2; w <- 2.1; r <- 1.6; rv <- .2^2
  exact <- catmech_contamination_row_terms(w, numeric(), m, s2, 0)
  expect_equal(exact$base, dnorm(w, m, .5, log = TRUE))
  joint <- integrate(function(t) dnorm(t, m, .5) * dnorm(w, t + mu, .3), -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(exp(exact$add) * dnorm(w, exact$mean1 + mu, sqrt(exact$var1 + t2)), joint, tolerance = 1e-7)
  noisy <- catmech_contamination_row_terms(c(w, w), c(r, NA), c(m, m), s2, rv)
  expect_equal(noisy$base[1], dnorm(w, m, .5, log = TRUE) + dnorm(r, w, .2, log = TRUE))
  joint <- integrate(function(t) dnorm(t, m, .5) * dnorm(w, t + mu, .3) * dnorm(r, t, .2), -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(exp(noisy$add[1]) * dnorm(w, noisy$mean1[1] + mu, sqrt(noisy$var1[1] + t2)), joint, tolerance = 1e-7)
  # A row without a reference uses the ordinary recorded-only terms.
  expect_equal(noisy$base[2], exact$base)
  expect_equal(noisy$mean1[2], exact$mean1)
  expect_equal(noisy$var1[2], exact$var1)
})

test_that("slice updates leave a known density invariant", {
  set.seed(14)
  x <- numeric(6000)
  current <- 3
  for (i in seq_along(x)) x[i] <- current <- catmech_slice(current, function(z) dnorm(z, 1, 2, log = TRUE), 1)
  kept <- x[-(1:200)]
  expect_lt(abs(mean(kept) - 1), .2)
  expect_lt(abs(var(kept) - 4), .5)
  expect_identical(catmech_slice(1, function(z) -Inf, 1), 1)
})

test_that("references are preserved and recorded values are kept or replaced whole", {
  f <- contamination_fixture()
  set.seed(21)
  before <- .Random.seed
  fit <- contamination_fit(f$data)
  expect_identical(.Random.seed, before)
  known <- !is.na(f$data$gold)
  recorded <- !is.na(f$data$lab_value)
  for (j in 1:4) {
    done <- fit$completed[[j]]
    expect_equal(done$lab_value[known], f$data$gold[known])
    expect_identical(done[setdiff(names(f$data), "lab_value")], f$data[setdiff(names(f$data), "lab_value")])
    expect_equal(mice::complete(fit$mids, j)$lab_value, done$lab_value)
    a <- fit$audit[fit$audit$imputation == j, ]
    retained <- a$reason == "recorded_retained"
    expect_true(all(done$lab_value[retained] == f$data$lab_value[retained]))
    replaced <- a$reason == "recorded_replaced"
    expect_true(all(done$lab_value[replaced] != f$data$lab_value[replaced]))
    expect_true(all(a$reason[known] == "gold_reference"))
    expect_true(all(a$reason[!known & !recorded] == "reference_missing"))
    expect_true(all(a$reconstructed == a$reason %in% c("recorded_replaced", "reference_missing")))
    expect_equal(fit$diagnostics$recorded_retained[j], sum(retained))
    expect_equal(fit$diagnostics$recorded_replaced[j], sum(replaced))
  }
  expect_identical(fit$original, f$data)
  expect_identical(fit$mids$data$lab_value, f$data$gold)
  expect_identical(fit$manifest$engine, "reference_contamination")
  expect_identical(fit$manifest$imputed_rows, which(!known))
  expect_equal(fit$manifest$assumption_hash, catmech_hash(fit$spec))
  expect_identical(fit$mids$call[[1]], as.name("catmech_contamination"))
  again <- contamination_fit(f$data)
  expect_identical(again$completed, fit$completed)
})

test_that("corruption probabilities separate real errors and use the flag", {
  f <- contamination_fixture(n = 400L)
  fit <- contamination_fit(f$data)
  table <- fit$contamination$corruption
  expect_equal(table$row, seq_len(400L))
  open <- table$evidence == "posterior_model"
  expect_gt(mean(table$posterior_corruption[open & f$wrong]), .95)
  expect_lt(mean(table$posterior_corruption[open & !f$wrong]), .05)
  expect_true(all(is.na(table$posterior_corruption[table$evidence == "no_recorded_value"])))
  disagree <- table$evidence == "reference_disagrees"
  expect_true(any(disagree))
  expect_true(all(table$posterior_corruption[disagree] == 1))
  expect_true(all(table$posterior_corruption[table$evidence == "reference_agrees"] == 0))
  p <- fit$contamination$parameters
  expect_gt(p$mean[p$parameter == "pi:Yes"], p$mean[p$parameter == "pi:No"])
  # The simulator's errors have mean 10.
  expect_lt(abs(p$mean[p$parameter == "mu:all"] - 10), 2)
  # Pooled analysis is labeled and recovers the simulator's slope of 0.8.
  pooled <- catmech_pool_lm(fit, outcome ~ lab_value + age + inflammation)
  row <- pooled[pooled$term == "lab_value", ]
  expect_identical(row$interpretation, "conditional_contamination_model")
  expect_true(row$lower < .8 && .8 < row$upper)
})

test_that("clean recordings stay almost untouched without a positive error signal", {
  d <- simulate_catmech_lab(240L, seed = 91, flag_rate = 0, corruption_probability = 0)
  d$gold <- NA_real_
  d$gold[seq(1L, 240L, 3L)] <- d$lab_true[seq(1L, 240L, 3L)]
  truth <- d$lab_true
  d[c("lab_true", "corrupted")] <- NULL
  fit <- catmech_contamination(d, "lab_value", "gold", contamination_predictors, m = 5, seed = 2)
  expect_identical(fit$contamination$recording_audit$status, "no_errors_observed")
  open <- !is.na(d$lab_value) & is.na(d$gold)
  expect_lt(mean(fit$contamination$corruption$posterior_corruption[open]), .05)
  expect_gt(mean(fit$completed[[1]]$lab_value[open] == d$lab_value[open]), .9)
})

test_that("complete references need no stochastic fit and match ordinary regression", {
  f <- contamination_fixture(every = 1L)
  fit <- contamination_fit(f$data)
  expect_length(fit$manifest$imputed_rows, 0L)
  expect_equal(fit$completed[[1]]$lab_value, f$data$gold)
  expect_identical(fit$model_diagnostics$status, "not_applicable")
  pooled <- catmech_pool_lm(fit, outcome ~ lab_value + age)
  expected <- lm(outcome ~ lab_value + age, data = fit$completed[[1]])
  expect_equal(unname(as.matrix(pooled[c("lower", "upper")])), unname(confint(expected)))
  expect_equal(pooled$df, rep(df.residual(expected), 3L))
})

test_that("unsupported or inconsistent requests are refused", {
  f <- contamination_fixture()
  d <- f$data
  expect_error(catmech_contamination(d, "lab_value", "absent", "age"), "existing column")
  expect_error(catmech_contamination(d, "lab_value", "gold", c("age", "gold")), "exclude")
  expect_error(contamination_fit(d, agreement_tol = .1, reference_sd = .1), "agreement_tol")
  expect_error(contamination_fit(d, agreement_tol = -1), "nonnegative")
  expect_error(contamination_fit(d, selection_delta = NA_real_), "selection_delta")
  expect_error(contamination_fit(d, prior = list(unknown = 1)), "named list")
  expect_error(contamination_fit(d, prior = list(nu0 = -1)), "positive")
  expect_error(contamination_fit(d, iterations = 50), "iterations")
  expect_error(contamination_fit(d, max_iterations = 100), "max_iterations")
  expect_error(catmech_contamination(d, "lab_value", "gold", "age", corruption_strata = "lab_value"), "target or reference")
  none <- d
  none$gold[] <- NA_real_
  expect_error(catmech_contamination(none, "lab_value", "gold", "age"), "independently observed")
  # References only where no value was recorded cannot calibrate recording errors.
  apart <- d
  apart$gold <- ifelse(is.na(d$lab_value), 1 + d$age / 50, NA_real_)
  expect_error(catmech_contamination(apart, "lab_value", "gold", "age"), "both recorded and reference")
  # An unaudited stratum has no evidence about its error rate.
  unaudited <- d
  unaudited$gold[d$outlier_flag == "Yes"] <- NA_real_
  expect_error(contamination_fit(unaudited), "No calibration rows support stratum: Yes")
  missing_flag <- d
  missing_flag$outlier_flag[which(!is.na(d$lab_value))[1]] <- NA
  expect_error(contamination_fit(missing_flag), "complete on recorded rows")
  # A predictor level seen only outside the reference sample is unsupported.
  sparse <- d
  sparse$only_unreferenced <- ifelse(is.na(d$gold) & seq_len(nrow(d)) %% 2 == 0, 1, 0)
  expect_error(catmech_contamination(sparse, "lab_value", "gold", "only_unreferenced"), "full column rank")
})

test_that("failed convergence refuses after the declared restarts", {
  f <- contamination_fixture()
  failed <- tryCatch(contamination_fit(f$data, burnin = 0, iterations = 100, max_iterations = 100), error = identity)
  expect_s3_class(failed, "catmech_convergence_error")
  expect_true(is.array(failed$model_diagnostics$traces))
  expect_identical(failed$model_diagnostics$engine, "gibbs_contamination")
  fit <- contamination_fit(f$data, iterations = 100, max_iterations = 1600)
  expect_gt(fit$manifest$iterations_used, 100L)
  expect_true(fit$model_diagnostics$convergence_pass)
})

test_that("agreement tolerance separates rounding from recording errors", {
  d <- simulate_catmech_lab(240L, seed = 33, flag_rate = 0, corruption_probability = 0)
  d$gold <- NA_real_
  d$gold[seq(1L, 240L, 3L)] <- d$lab_true[seq(1L, 240L, 3L)]
  d$lab_value <- round(d$lab_value, 1)
  d[c("lab_true", "corrupted")] <- NULL
  fit <- catmech_contamination(d, "lab_value", "gold", contamination_predictors, agreement_tol = .05, m = 4, seed = 3)
  expect_identical(fit$contamination$recording_audit$status, "no_errors_observed")
  expect_true(all(fit$contamination$corruption$evidence[!is.na(d$gold) & !is.na(d$lab_value)] == "reference_agrees"))
  # References keep their measured value even though the recorded value is rounded.
  expect_equal(fit$completed[[1]]$lab_value[!is.na(d$gold)], d$gold[!is.na(d$gold)])
  expect_warning(catmech_contamination(d, "lab_value", "gold", contamination_predictors, m = 4, seed = 3,
                                       max_iterations = 2400), "Nearly all calibration rows disagree")
})

test_that("a declared noisy reference is modeled rather than copied", {
  f <- contamination_fixture()
  d <- f$data
  known <- !is.na(d$gold)
  d$gold[known] <- d$gold[known] + rnorm(sum(known), 0, .2)
  fit <- contamination_fit(d, reference_sd = .2, m = 4, max_iterations = 2400)
  expect_false(any(fit$completed[[1]]$lab_value[known] == d$gold[known]))
  expect_identical(fit$completed[[1]]$gold, d$gold)
  expect_identical(fit$manifest$imputed_rows, seq_len(nrow(d)))
  expect_true(all(is.na(fit$mids$data$lab_value)))
  expect_equal(mice::complete(fit$mids, 2)$lab_value, fit$completed[[2]]$lab_value)
  expect_true(all(fit$contamination$corruption$evidence[known & !is.na(d$lab_value)] == "posterior_model"))
  expect_true("noisy_reference_model" %in% fit$audit$reason)
  expect_false("gold_reference" %in% fit$audit$reason)
  expect_true(all(catmech_pool_lm(fit, outcome ~ lab_value + age)$interpretation == "conditional_contamination_model"))
})

test_that("the selection shift moves only unreferenced imputations and shows in reports", {
  f <- contamination_fixture()
  fits <- list(minus = contamination_fit(f$data, selection_delta = -.5, seed = 8),
               zero = contamination_fit(f$data, seed = 8),
               plus = contamination_fit(f$data, selection_delta = .5, seed = 8))
  unavailable <- is.na(f$data$lab_value) & is.na(f$data$gold)
  means <- vapply(fits, function(z) mean(z$completed[[1]]$lab_value[unavailable]), numeric(1))
  expect_true(means[["minus"]] < means[["zero"]] && means[["zero"]] < means[["plus"]])
  known <- !is.na(f$data$gold)
  expect_equal(fits$plus$completed[[1]]$lab_value[known], f$data$gold[known])
  report <- catmech_sensitivity(fits, outcome ~ lab_value + age, "lab_value", expected_scenarios = names(fits))
  expect_true(report$complete)
  expect_true(all(report$table$interpretation == "conditional_contamination_model"))
  expect_false(identical(fits$plus$manifest$assumption_hash, fits$zero$manifest$assumption_hash))
})

test_that("results serialize with the existing writer", {
  f <- contamination_fixture(n = 120L)
  fit <- contamination_fit(f$data, seed = 4)
  out <- tempfile("contamination")
  write_catmech(fit, out)
  expect_true(all(file.exists(file.path(out, c("result.rds", "completed_long.csv", "audit.csv", "manifest.json", "assumptions.canonical.json")))))
  expect_identical(jsonlite::read_json(file.path(out, "manifest.json"))$engine, "reference_contamination")
})

test_that("the recording audit reports exact counts and clear statuses", {
  d <- data.frame(w = c(1, 2, 3, 4, 5, 6, NA, 8), r = c(1, 2.5, 3, NA, 5.5, 6, 7, NA),
                  site = c("a", "a", "a", "a", "b", "b", "b", "b"))
  a <- catmech_recording_audit(d, "w", "r", by = "site")
  expect_identical(a$status, "few_errors_observed")
  expect_equal(a$overall$rows, 8L)
  expect_equal(a$overall$reference_rows, 6L)
  expect_equal(a$overall$calibration_rows, 5L)
  expect_equal(a$overall$disagreements, 2L)
  expect_equal(a$overall$error_mean, -.5)
  expect_equal(a$overall$lower, qbeta(.025, 2.5, 3.5))
  site_a <- a$by[a$by$level == "a", ]
  expect_equal(site_a$calibration_rows, 3L)
  expect_equal(site_a$disagreements, 1L)
  expect_equal(a$heterogeneity$p_value, fisher.test(matrix(c(2, 1, 1, 1), 2))$p.value)
  expect_identical(catmech_recording_audit(transform(d, r = w), "w", "r")$status, "no_errors_observed")
  expect_identical(catmech_recording_audit(transform(d, r = w + 1), "w", "r")$status, "mostly_disagreeing")
  expect_identical(catmech_recording_audit(transform(d, r = ifelse(is.na(w), 1, NA)), "w", "r")$status, "no_calibration_rows")
  expect_identical(catmech_recording_audit(d, "w", "r", agreement_tol = .5)$status, "no_errors_observed")
  expect_error(catmech_recording_audit(d, "w", "r", by = "absent"), "existing columns")
  expect_error(catmech_recording_audit(d, "absent", "r"), "numeric recorded")
  # Selective audits are visible from recorded values alone.
  skew <- data.frame(w = 1:40)
  skew$r <- ifelse(skew$w > 30, skew$w, NA)
  selected <- catmech_recording_audit(skew, "w", "r")
  expect_gt(selected$selection$value, 1)
  expect_true(any(grepl("may not be simple random", selected$notes)))
})

test_that("the hold-out check scores unseen references and restores RNG", {
  f <- contamination_fixture(n = 300L)
  fit <- contamination_fit(f$data, m = 4)
  set.seed(77)
  before <- .Random.seed
  h <- catmech_holdout_check(fit, folds = 3, m = 20, seed = 6)
  expect_identical(.Random.seed, before)
  known <- which(!is.na(f$data$gold))
  expect_equal(sort(h$rows$row), known)
  expect_equal(h$summary$holdout_rows, length(known))
  expect_true(h$status %in% c("no_miscalibration_detected", "miscalibration_detected"))
  expect_true(all(h$rows$pit > 0 & h$rows$pit < 1))
  expect_lt(h$summary$brier_model, h$summary$brier_stratum_rate)
  expect_lt(h$summary$rmse_model, h$summary$rmse_recorded)
  normal <- catmech_reference(f$data, "lab_value", "gold", contamination_predictors, m = 4, seed = 2)
  n <- catmech_holdout_check(normal, folds = 3, m = 20, seed = 6)
  expect_true(all(is.na(n$rows$posterior_corruption)))
  expect_equal(n$summary$holdout_rows, length(known))
  # With exact references and a correct model the error count is close to expected.
  expect_lt(abs(h$summary$calibration_z), 4)
  pmm <- catmech_pmm(f$data, "lab_value", contamination_predictors, flag = "outlier_flag", m = 2)
  expect_error(catmech_holdout_check(pmm), "support")
  expect_error(catmech_holdout_check(fit, folds = 1), "folds")
  expect_error(catmech_holdout_check(fit, level = 1), "level")
})

test_that("the hold-out check reports folds that lose support", {
  f <- contamination_fixture(n = 150L)
  d <- f$data
  # One flagged reference row: hiding it leaves the flagged stratum unsupported.
  flagged_ref <- which(d$outlier_flag == "Yes" & !is.na(d$gold) & !is.na(d$lab_value))
  d$gold[flagged_ref[-1L]] <- NA_real_
  fit <- contamination_fit(d, m = 4)
  h <- catmech_holdout_check(fit, folds = 3, m = 20, seed = 2)
  expect_identical(h$status, "unresolved_failed_folds")
  expect_true(any(grepl("No calibration rows support", h$failures)))
  expect_false(is.null(h$summary))
})

test_that("bootstrap matching draws from references and adds between-imputation spread", {
  f <- contamination_fixture(n = 300L, every = 4L)
  d <- f$data
  predictors <- c("age", "outcome")
  set.seed(5)
  before <- .Random.seed
  abb <- catmech_reference(d, "lab_value", "gold", predictors, method = "pmm_abb", m = 20, seed = 9)
  expect_identical(.Random.seed, before)
  pmm <- catmech_reference(d, "lab_value", "gold", predictors, method = "pmm", m = 20, seed = 9)
  known <- !is.na(d$gold)
  for (j in c(1L, 20L)) {
    expect_true(all(abb$completed[[j]]$lab_value[!known] %in% d$gold[known]))
    expect_equal(abb$completed[[j]]$lab_value[known], d$gold[known])
  }
  expect_identical(abb$manifest$engine, "reference_pmm_abb")
  expect_identical(abb$spec$donor_resampling, "approximate_bayesian_bootstrap_each_imputation")
  expect_null(pmm$spec$donor_resampling)
  expect_identical(abb$manifest$claim_level, "conditional_reference_matching")
  pooled <- catmech_pool_lm(abb, outcome ~ lab_value + age)
  plain <- catmech_pool_lm(pmm, outcome ~ lab_value + age, allow_sensitivity = TRUE)
  expect_true(all(pooled$interpretation == "conditional_reference_matching"))
  expect_gt(pooled$between_variance[2], plain$between_variance[2])
  expect_true(all(abb$audit$donor_row[!abb$audit$reference_observed] %in% c(which(known), NA)))
  expect_true(all(abb$audit$eligible_donors[!abb$audit$reference_observed] == sum(known)))
  few <- d
  few$gold[-(1:3)] <- NA_real_
  expect_error(catmech_reference(few, "lab_value", "gold", "age", method = "pmm_abb"), "Fewer trusted/reference")
})

test_that("the hold-out check says when a model shows no skill", {
  # Recorded values are exact, so nothing can beat them.
  d <- simulate_catmech_lab(240L, seed = 91, flag_rate = 0, corruption_probability = 0)
  d$gold <- NA_real_
  d$gold[seq(1L, 240L, 3L)] <- d$lab_true[seq(1L, 240L, 3L)]
  d[c("lab_true", "corrupted")] <- NULL
  fit <- catmech_contamination(d, "lab_value", "gold", contamination_predictors, m = 4, seed = 2)
  h <- catmech_holdout_check(fit, folds = 3, m = 20, seed = 4)
  expect_true(any(grepl("no closer", h$advice)))
  expect_true(any(grepl("use_recorded", h$advice)))
  # Large, separable errors leave the model with real skill.
  f <- contamination_fixture(n = 300L)
  skilled <- catmech_holdout_check(contamination_fit(f$data), folds = 3, m = 20, seed = 6)
  expect_length(skilled$advice, 0L)
  # The linear workflow can be scored the same way.
  linear <- catmech_reference(f$data, "lab_value", "gold", contamination_predictors, use_recorded = TRUE, m = 4, seed = 5)
  scored <- catmech_holdout_check(linear, folds = 3, m = 20, seed = 7)
  expect_identical(scored$summary$engine, "reference_normal_recorded")
  expect_true(all(is.na(scored$rows$posterior_corruption)))
})

test_that("catmech_recording_audit and catmech_holdout_check S3 print methods format cleanly", {
  f <- contamination_fixture()
  audit <- catmech_recording_audit(f$data, "lab_value", "gold", by = "outlier_flag")
  expect_s3_class(audit, "catmech_recording_audit")
  a_txt <- utils::capture.output(a_val <- print(audit))
  expect_identical(a_val, audit)
  expect_true(any(grepl("Catmech Recording Audit", a_txt)))
  expect_true(any(grepl("Status:", a_txt)))
  expect_true(any(grepl("Disagreements:", a_txt)))

  linear <- catmech_reference(f$data, "lab_value", "gold", contamination_predictors, use_recorded = TRUE, m = 4, seed = 5)
  scored <- catmech_holdout_check(linear, folds = 3, m = 20, seed = 7)
  expect_s3_class(scored, "catmech_holdout_check")
  h_txt <- utils::capture.output(h_val <- print(scored))
  expect_identical(h_val, scored)
  expect_true(any(grepl("Catmech Holdout Cross-Validation Check", h_txt)))
  expect_true(any(grepl("Predictive coverage:", h_txt)))
  expect_true(any(grepl("Predictive RMSE:", h_txt)))

  # Test plot methods
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(audit))
  expect_silent(plot(scored))
})
