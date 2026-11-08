# Regression tests for the fixes from the 2026-11-07 review of 1.0.0
# (validation/REVIEW_FIXES_2026-11-07.md).

review_measurement_data <- function(n, seed, slope = 0.8, binary = FALSE) {
  set.seed(seed)
  t <- stats::rnorm(n)
  x <- 0.3 * t + stats::rnorm(n)
  w <- t + stats::rnorm(n, sd = 0.7)
  y <- if (binary) stats::rbinom(n, 1, stats::plogis(-0.3 + slope * t + 0.3 * x)) else
    1 + slope * t + 0.3 * x + stats::rnorm(n)
  audited <- seq_len(n) %in% sample.int(n, round(0.4 * n))
  data.frame(lab = w, ref = ifelse(audited, t, NA_real_), x = x, y = y,
             site = rep(sprintf("s%02d", seq_len(12L)), length.out = n))
}

collect_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = messages)
}

test_that("compare_strategies keeps the outcome among detected imputation predictors", {
  skip_on_cran() # slow; GitHub Actions runs it with NOT_CRAN=true
  d <- review_measurement_data(1500L, 4101L)
  run <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x, m = 5L, seed = 7L))
  comp <- run$value
  expect_setequal(attr(comp, "imputation_predictors"), c("x", "y"))
  expect_false(any(grepl("was not an imputation predictor", run$warnings)))
  expect_output(print(comp), "Imputation predictors: y, x")
  slope <- function(cmp, strategy) cmp$estimate[cmp$strategy == strategy & cmp$term == "lab"]
  # True slope 0.8. With the outcome the linear-imputation slope recovers it; the
  # old default (predictors without the outcome) attenuated it to about 0.63.
  expect_lt(abs(slope(comp, "linear_mi_recorded") - 0.8), 0.08)
  old <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x,
                                                     predictors = "x", m = 5L, seed = 7L))
  expect_lt(slope(old$value, "linear_mi_recorded"), 0.72)
  expect_true(any(grepl("outcome 'y' was not an imputation predictor", old$warnings)))

  # With the audited value as the outcome, the analysis predictors are the imputation predictors.
  as_outcome <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", lab ~ x, m = 3L, seed = 7L))$value
  expect_identical(attr(as_outcome, "imputation_predictors"), "x")
})

test_that("compare_strategies uses cluster only for variances, the same way in every strategy", {
  d <- review_measurement_data(360L, 4102L)
  run <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x,
                                                     m = 3L, seed = 9L, cluster = "site"))
  comp <- run$value
  expect_false(any(grepl("new intercept|Joint clustered", run$warnings)))
  expect_identical(attr(comp, "cluster"), "site")
  direct <- comp[comp$strategy %in% c("validation_subsample", "naive_trust_recorded"), ]
  expect_true(all(direct$df == 11))
  # The naive row equals a direct delete-one-cluster jackknife of the same fit.
  fit <- stats::lm(y ~ lab + x, data = d)
  v <- catmech:::catmech_cluster_jackknife(fit, d$site)
  naive <- comp[comp$strategy == "naive_trust_recorded", ]
  expect_equal(naive$std.error, unname(sqrt(diag(v))[naive$term]), tolerance = 1e-12)
  expect_error(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x, cluster = "nope"),
               "cluster must name a column")
})

test_that("compare_strategies explains a failed contamination fit", {
  d <- review_measurement_data(80L, 4103L)
  run <- collect_warnings(with_mocked_bindings(
    catmech_compare_strategies(d, "lab", "ref", y ~ lab + x, m = 2L, seed = 3L),
    catmech_contamination = function(...) stop("mock contamination failure")
  ))
  expect_true(any(grepl("contamination model failed.*mock contamination failure", run$warnings)))
  expect_true(all(is.na(run$value$estimate[run$value$strategy == "contamination_model"])))
})

test_that("pooling warns when the audited value is a predictor and the outcome was not imputed with", {
  d <- review_measurement_data(200L, 4104L, binary = TRUE)
  without <- catmech_reference(d, "lab", "ref", "x", m = 3L, seed = 1L, use_recorded = TRUE)
  with <- catmech_reference(d, "lab", "ref", c("x", "y"), m = 3L, seed = 1L, use_recorded = TRUE)
  expect_warning(catmech_pool_glm(without, y ~ lab + x), "outcome 'y' was not an imputation predictor")
  expect_warning(catmech_pool_lm(without, y ~ lab + x), "outcome 'y' was not an imputation predictor")
  expect_warning(catmech_pool_glm(without, (y > 0) ~ lab + x), "outcome 'y'")
  expect_no_warning(catmech_pool_glm(with, y ~ lab + x))
  expect_no_warning(catmech_pool_lm(with, y ~ lab + x))
  # The audited value as the outcome needs no outcome in the imputation model.
  expect_no_warning(catmech_pool_lm(without, lab ~ x))
})

test_that("Cox pooling warns without the event indicator in the imputation model", {
  set.seed(4105L)
  n <- 240L
  t <- stats::rnorm(n)
  w <- t + stats::rnorm(n, sd = 0.6)
  time <- stats::rexp(n, exp(0.5 * t) * 0.1)
  censor <- stats::runif(n, 1, 15)
  d <- data.frame(lab = w, ref = ifelse(seq_len(n) <= 100L, t, NA_real_),
                  age = stats::rnorm(n), time = pmin(time, censor),
                  status = as.numeric(time <= censor))
  d$hazard <- mice::nelsonaalen(d, time, status)
  without <- catmech_reference(d, "lab", "ref", "age", m = 3L, seed = 2L, use_recorded = TRUE)
  with <- catmech_reference(d, "lab", "ref", c("age", "status", "hazard"), m = 3L, seed = 2L,
                            use_recorded = TRUE)
  expect_warning(catmech_pool_coxph(without, survival::Surv(time, status) ~ lab + age),
                 "event indicator 'status' was not an imputation predictor")
  expect_no_warning(catmech_pool_coxph(with, survival::Surv(time, status) ~ lab + age))
})

test_that("the GLM sandwich does not depend on an estimated dispersion", {
  set.seed(4106L)
  g <- rep(seq_len(30L), each = 20L)
  x <- stats::rnorm(600L)
  mu <- exp(0.5 + 0.3 * x + stats::rnorm(30L, sd = 0.3)[g])
  y <- stats::rnbinom(600L, mu = mu, size = 1.5)
  pois <- stats::glm(y ~ x, family = stats::poisson())
  quasi <- stats::glm(y ~ x, family = stats::quasipoisson())
  expect_gt(summary(quasi)$dispersion, 2)
  v_pois <- catmech:::catmech_glm_cluster_vcov(pois, g)
  v_quasi <- catmech:::catmech_glm_cluster_vcov(quasi, g)
  expect_equal(v_quasi, v_pois, tolerance = 1e-10)
  # Against a direct CR1S sandwich built from the unscaled bread.
  X <- stats::model.matrix(pois)
  meat <- crossprod(rowsum(X * (pois$weights * stats::residuals(pois, "working")), g))
  bread <- summary(pois)$cov.unscaled
  manual <- bread %*% meat %*% bread * (30 / 29) * (599 / 598)
  expect_equal(unname(v_pois), unname(manual), tolerance = 1e-10)
  # For binomial and poisson the bread equals vcov(), so earlier results are unchanged.
  expect_equal(unname(summary(pois)$cov.unscaled), unname(stats::vcov(pois)), tolerance = 1e-12)
})

test_that("audit power reports the linear-imputation size without a plan", {
  pwr <- catmech_audit_power(N = 1000, target_ratio = 1.25, error_rate = 0.15,
                             correlation = 0.85, r_squared = 0.3, n_sim = 200L, seed = 1L)
  k <- (1 - 0.85^2) * (1 - 0.3)
  expect_equal(pwr$linear_mi$k, k)
  expect_equal(pwr$linear_mi$fraction, k / (1.25^2 - 1 + k))
  expect_identical(pwr$linear_mi$n_audit, as.integer(ceiling(1000 * k / (1.25^2 - 1 + k))))
  expect_identical(pwr$status, "not_validated")
  expect_gt(pwr$linear_mi$fraction, pwr$recommended$fraction)
  expect_output(print(pwr), "ignoring uncertainty in k")
  # The figure ignores uncertainty in k, so at very high error rates it can be smaller
  # than the main recommendation; the documentation says so rather than calling it an upper bound.
  high <- catmech_audit_power(N = 1000, error_rate = 0.9, n_sim = 200L, seed = 1L)
  expect_lt(high$linear_mi$fraction, high$recommended$fraction)
  expect_null(catmech_audit_power(plan = structure(list(
    k = data.frame(strategy = "recorded", estimate = 0.4, lower = 0.2, upper = 0.7),
    pilot = list(fit = list(original = data.frame(x = 1:500)))), class = "catmech_audit_plan"))$linear_mi)
})

test_that("the provenance record makes no tamper-evidence claim and escapes its HTML title", {
  set.seed(4107L)
  df <- data.frame(lab = stats::rnorm(30, 10, 2), age = stats::rnorm(30, 50, 10),
                   y = stats::rnorm(30, 5, 1))
  df$ref <- df$lab
  df$ref[11:30] <- NA
  res <- catmech_reference(df, "lab", "ref", c("age", "y"), m = 2L, seed = 1L)
  md <- catmech_audit_trail(res, pooled = catmech_pool_lm(res, y ~ lab + age))
  expect_match(md, "do not make this document tamper-evident", fixed = TRUE)
  expect_false(grepl("tied via SHA-256", md, fixed = TRUE))
  html <- catmech_audit_trail(res, format = "html", title = "Run <b>1</b> & checks")
  expect_match(html, "<title>Run &lt;b&gt;1&lt;/b&gt; &amp; checks</title>", fixed = TRUE)
})

test_that("cluster variances in compare_strategies drop coefficients the audited rows cannot estimate", {
  skip_on_cran() # slow; GitHub Actions runs it with NOT_CRAN=true
  set.seed(4108L)
  n <- 360L
  d <- data.frame(site = rep(sprintf("s%02d", seq_len(12L)), length.out = n), inpatient = stats::rbinom(n, 1, 0.5))
  t <- stats::rnorm(n)
  d$x <- stats::rnorm(n)
  d$lab <- t + stats::rnorm(n, sd = 0.7)
  d$y <- 1 + 0.8 * t + 0.3 * d$x + 0.5 * d$inpatient + stats::rnorm(n)
  d$yb <- as.integer(d$y > 1)
  audited <- sample(which(d$inpatient == 1L), 120L)   # inpatient is constant among audited rows
  d$ref <- ifelse(seq_len(n) %in% audited, t, NA_real_)
  for (family in list(stats::gaussian(), stats::binomial())) {
    gaussian <- identical(family$family, "gaussian")
    outcome <- if (gaussian) "y" else "yb"
    f <- stats::as.formula(paste(outcome, "~ inpatient + lab + x"))
    cmp <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", f, predictors = c("x", outcome),
                                                       family = family, m = 2L, seed = 5L, cluster = "site"))$value
    sub <- cmp[cmp$strategy == "validation_subsample", ]
    expect_setequal(sub$term, c("(Intercept)", "lab", "x"))
    # Against a direct computation on the estimable columns of the audited rows.
    a <- d[audited, ]
    a$lab <- a$ref
    reduced <- stats::as.formula(paste(outcome, "~ lab + x"))
    fit <- if (gaussian) stats::lm(reduced, data = a) else stats::glm(reduced, family = family, data = a)
    v <- if (gaussian) catmech:::catmech_cluster_jackknife(fit, a$site) else catmech:::catmech_glm_cluster_vcov(fit, a$site)
    expect_equal(sub$std.error, unname(sqrt(diag(v))[sub$term]), tolerance = 1e-8)
    expect_true(all(sub$df == 11))
  }
})

test_that("a working contamination fit gets the same cluster-robust variance as the others", {
  d <- review_measurement_data(240L, 4109L)
  # Stand in for a converged contamination fit, so the test does not hinge on MCMC.
  standin <- function(data, value, reference, predictors, m, seed, ...)
    catmech_reference(data, value, reference, predictors, m = m, seed = seed, use_recorded = TRUE)
  cmp <- with_mocked_bindings(
    collect_warnings(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x, m = 3L, seed = 4L, cluster = "site"))$value,
    catmech_contamination = standin)
  contamination <- cmp[cmp$strategy == "contamination_model", ]
  linear <- cmp[cmp$strategy == "linear_mi_recorded", ]
  expect_true(all(is.finite(contamination$estimate) & is.finite(contamination$std.error)))
  expect_equal(contamination$std.error, linear$std.error, tolerance = 1e-12)
  pooled <- catmech_pool_lm(catmech_reference(d, "lab", "ref", c("y", "x"), m = 3L, seed = 4L, use_recorded = TRUE),
                            y ~ lab + x, cluster = "site")
  expect_equal(linear$std.error, pooled$std.error, tolerance = 1e-12)
})

test_that("compare_strategies with a binomial family uses the CR1S variance with clusters minus one df", {
  skip_on_cran() # slow; GitHub Actions runs it with NOT_CRAN=true
  d <- review_measurement_data(360L, 4110L, binary = TRUE)
  cmp <- collect_warnings(catmech_compare_strategies(d, "lab", "ref", y ~ lab + x, family = stats::binomial(),
                                                     m = 2L, seed = 6L, cluster = "site"))$value
  sub <- cmp[cmp$strategy == "validation_subsample", ]
  expect_true(all(sub$df == 11))
  a <- d[!is.na(d$ref), ]
  a$lab <- a$ref
  fit <- stats::glm(y ~ lab + x, family = stats::binomial(), data = a)
  expect_equal(sub$std.error, unname(sqrt(diag(catmech:::catmech_glm_cluster_vcov(fit, a$site)))[sub$term]),
               tolerance = 1e-10)
})

test_that("GLM and Cox pooling give the clustered-imputation warning", {
  set.seed(4111L)
  n <- 240L
  t <- stats::rnorm(n)
  d <- data.frame(id = rep(sprintf("c%02d", seq_len(12L)), each = 20L), recorded = t + stats::rnorm(n, sd = 0.6))
  d$yb <- stats::rbinom(n, 1, stats::plogis(0.6 * t))
  time <- stats::rexp(n, exp(0.5 * t) * 0.1)
  censor <- stats::runif(n, 1, 15)
  d$time <- pmin(time, censor)
  d$status <- as.numeric(time <= censor)
  d$hazard <- mice::nelsonaalen(d, time, status)
  d$reference <- ifelse(d$id %in% sprintf("c%02d", 1:8) & stats::runif(n) < 0.5, t, NA_real_)
  fit <- suppressWarnings(catmech_reference(d, "recorded", "reference", c("yb", "status", "hazard"),
                                            m = 3L, seed = 3L, use_recorded = TRUE, cluster = "id"))
  expect_warning(catmech_pool_glm(fit, yb ~ recorded), "new intercept")
  expect_warning(catmech_pool_coxph(fit, survival::Surv(time, status) ~ recorded), "new intercept")
})

test_that("the outcome check skips unimputed fits, reads Surv() by argument and covers catmech_model()", {
  d <- review_measurement_data(150L, 4112L)
  d$ref <- d$lab + stats::rnorm(nrow(d), sd = 0.3)   # every row audited: nothing is imputed
  none <- catmech_reference(d, "lab", "ref", "x", m = 2L, seed = 1L, use_recorded = TRUE)
  expect_no_warning(catmech_pool_lm(none, y ~ lab + x))

  set.seed(4113L)
  n <- 200L
  t <- stats::rnorm(n)
  s <- data.frame(lab = t + stats::rnorm(n, sd = 0.5), ref = ifelse(seq_len(n) <= 90L, t, NA_real_),
                  age = stats::rnorm(n), time = stats::rexp(n, 0.1), status = stats::rbinom(n, 1, 0.6))
  fit <- catmech_reference(s, "lab", "ref", c("age", "time"), m = 2L, seed = 2L, use_recorded = TRUE)
  expect_warning(catmech_pool_coxph(fit, survival::Surv(event = status, time = time) ~ lab + age),
                 "event indicator 'status'")
  expect_warning(catmech_pool_coxph(fit, survival::Surv(time, status == 1) ~ lab + age),
                 "event indicator 'status'")

  m <- review_measurement_data(300L, 4114L)
  m$lab[sample.int(300L, 150L)] <- NA_real_
  model_fit <- catmech_model(m[c("lab", "x", "y")], "lab", "x", m = 2L, seed = 3L)
  expect_warning(catmech_pool_lm(model_fit, y ~ lab + x), "outcome 'y' was not an imputation predictor")
})
