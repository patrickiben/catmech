model_fixture <- function() {
  d <- simulate_catmech_lab(100L, seed = 9, corruption_probability = 0.5)
  d[c("lab_true", "corrupted")] <- NULL
  d
}

test_that("Gaussian default masks every flag and preserves all other input", {
  d <- model_fixture()
  set.seed(41)
  before <- .Random.seed
  fit <- catmech_model(d, "lab_value", c("age", "inflammation", "outcome"),
    flag = "outlier_flag", category = "nonexistent", m = 4L, seed = 2L)
  expect_identical(.Random.seed, before)
  expect_identical(fit$original, d)
  expect_identical(fit$manifest$engine, "gaussian_model")
  expect_identical(fit$spec$category_mode, "ignore")
  expect_true(fit$model_diagnostics$convergence_pass)
  expect_identical(fit$model_diagnostics$status, "not_applicable")
  trusted <- !is.na(d$lab_value) & as.character(d$outlier_flag) != "Yes"
  for (j in seq_along(fit$completed)) {
    expect_identical(fit$completed[[j]][setdiff(names(d), "lab_value")],
                     d[setdiff(names(d), "lab_value")])
    expect_equal(fit$completed[[j]]$lab_value[trusted], d$lab_value[trusted])
    expect_false(anyNA(fit$completed[[j]]$lab_value))
    expect_equal(mice::complete(fit$mids, j)$lab_value, fit$completed[[j]]$lab_value)
  }
  expect_true(all(fit$audit$reconstructed))
  again <- catmech_model(d, "lab_value", c("age", "inflammation", "outcome"),
    flag = "outlier_flag", m = 4L, seed = 2L)
  expect_identical(fit$completed, again$completed)
})

test_that("Gaussian regression draws match the analytic intercept posterior", {
  y <- c(1, 2, 4, 7)
  design <- catmech:::catmech_model_design(matrix(1, 4L, 1L))
  fit <- catmech:::catmech_model_fit(design, y)
  expect_equal(fit$beta, mean(y))
  expect_equal(fit$sse, sum((y - mean(y))^2))
  set.seed(81)
  expected_sigma <- sqrt(sum((y - mean(y))^2) / stats::rchisq(1, 3))
  expected_deviation <- expected_sigma * stats::rnorm(1) / sqrt(4)
  set.seed(81)
  draw <- catmech:::catmech_model_parameter_draw(design, fit)
  expect_equal(draw$sigma, expected_sigma)
  expect_equal(abs(draw$beta - mean(y)), abs(expected_deviation))
})

test_that("intercept-only adapter cannot overwrite a colliding target name", {
  d <- data.frame(.catmech_intercept = c(2, 3, 7, 9, NA_real_))
  fit <- catmech_model(d, ".catmech_intercept", character(0), m = 2, seed = 5)
  expect_identical(fit$original, d)
  expect_equal(mice::complete(fit$mids, 1)$.catmech_intercept,
               fit$completed[[1]]$.catmech_intercept)
})

test_that("noisy-bin sampling matches an analytic known-parameter law", {
  H <- matrix(c(.8, .2, .2, .8), 2L, 2L)
  set.seed(74)
  draws <- catmech:::catmech_model_latent(rep(0, 20000L), 1,
    rep(2L, 20000L), c(-Inf, 0, Inf), H)
  expected_mean <- .6 * sqrt(2 / pi)
  expect_equal(mean(draws > 0), .8, tolerance = .012)
  expect_equal(mean(draws), expected_mean, tolerance = .025)
  expect_equal(mean(draws^2), 1, tolerance = .04)
  set.seed(75)
  exact <- catmech:::catmech_model_latent(rep(0, 20000L), 1,
    rep(2L, 20000L), c(-Inf, 0, Inf), diag(2))
  expect_true(all(exact > 0))
  expect_equal(mean(exact), sqrt(2 / pi), tolerance = .02)
})

test_that("asymmetric category probabilities use true rows and observed columns", {
  H <- matrix(c(.9, .1, .35, .65), 2L, 2L, byrow = TRUE)
  set.seed(76)
  draws <- catmech:::catmech_model_latent(rep(0, 20000L), 1,
    rep(2L, 20000L), c(-Inf, 0, Inf), H)
  positive_probability <- .65 / (.1 + .65)
  expect_equal(mean(draws > 0), positive_probability, tolerance = .012)
  expect_equal(mean(draws), (2 * positive_probability - 1) * sqrt(2 / pi),
               tolerance = .025)
})

test_that("Gaussian interval probabilities remain stable in extreme tails", {
  actual <- catmech:::catmech_model_log_interval(c(-41, 40), c(-40, 41))
  expect_true(all(is.finite(actual)))
  expect_equal(actual[1], actual[2])
  expect_equal(catmech:::catmech_model_log_interval(-Inf, Inf), 0)
})

test_that("exact Gibbs completion respects known bins and keeps missing categories missing", {
  d <- model_fixture()
  d$lab_category[1L] <- NA
  fit <- catmech_model(d, "lab_value", c("age", "inflammation", "outcome"),
    flag = "outlier_flag", category = "lab_category", levels = c("Low", "Medium", "High"),
    thresholds = c(-Inf, 5.5, 7, Inf), category_mode = "exact", m = 4,
    burnin = 100, iterations = 400, chains = 4, seed = 14)
  expect_true(fit$model_diagnostics$convergence_pass)
  expect_equal(dim(fit$model_diagnostics$traces)[1:2], c(400L, 4L))
  expect_true(all(fit$diagnostics$threshold_violations == 0))
  expect_true(all(vapply(fit$completed, function(z) is.na(z$lab_category[1]), logical(1))))
})

test_that("soft categories preserve conflicting labels while exact mode refuses", {
  d <- model_fixture()
  trusted <- which(!is.na(d$lab_value) & d$outlier_flag == "No")[1L]
  d$lab_category[trusted] <- if (d$lab_category[trusted] == "Low") "High" else "Low"
  args <- list(data = d, value = "lab_value", predictors = c("age", "inflammation", "outcome"),
    flag = "outlier_flag", category = "lab_category", levels = c("Low", "Medium", "High"),
    thresholds = c(-Inf, 5.5, 7, Inf), m = 4, seed = 17,
    burnin = 100, iterations = 400, chains = 4)
  expect_error(do.call(catmech_model, c(args, list(category_mode = "exact"))), "zero probability")
  H <- matrix(.05, 3L, 3L, dimnames = list(args$levels, args$levels)); diag(H) <- .9
  fit <- do.call(catmech_model, c(args, list(category_mode = "noisy", confusion = H)))
  expect_true(fit$model_diagnostics$convergence_pass)
  expect_identical(fit$completed[[1]]$lab_category, d$lab_category)
  expect_equal(fit$completed[[1]]$lab_value[trusted], d$lab_value[trusted])
  expect_true(fit$diagnostics$category_disagreements[1] > 0)
})

test_that("uninformative confusion rows give exactly the direct engine", {
  d <- model_fixture()
  lev <- c("Low", "Medium", "High")
  H <- matrix(rep(c(.2, .5, .3), each = 3L), 3L, 3L, dimnames = list(lev, lev))
  direct <- catmech_model(d, "lab_value", c("age", "outcome"),
    flag = "outlier_flag", m = 3, seed = 1)
  noisy <- catmech_model(d, "lab_value", c("age", "outcome"),
    flag = "outlier_flag", m = 3, seed = 1, category = "lab_category", levels = lev,
    thresholds = c(-Inf, 5.5, 7, Inf), category_mode = "noisy", confusion = H)
  expect_identical(direct$completed, noisy$completed)
  expect_identical(noisy$model_diagnostics$engine, "independent_normal")
})

test_that("model rejects unsupported design, malformed likelihood and finite domain", {
  d <- model_fixture()
  d$alias <- 2 * d$age
  expect_error(catmech_model(d, "lab_value", c("age", "alias"), flag = "outlier_flag"), "full column rank")
  expect_error(catmech_model(d, "lab_value", "age", category = "lab_category",
    thresholds = c(0, 5.5, 7, Inf), levels = c("Low", "Medium", "High"),
    category_mode = "exact"), "outer endpoints")
  expect_error(catmech:::catmech_model_confusion(diag(3), c("Low", "Medium", "High")), "named square")
  expect_error(catmech_model(d, "lab_value", "age", confusion = diag(3)), "confusion is unused")
  H <- diag(3); dimnames(H) <- list(c("Low", "Medium", "High"), c("Low", "Medium", "High"))
  H[1, 1] <- .5
  expect_error(catmech:::catmech_model_confusion(H, rownames(H)), "summing to one")
  H[,] <- rep(c(.5, .5, 0), each = 3)
  expect_error(catmech_model(d, "lab_value", "age", flag = "outlier_flag",
    category = "lab_category", thresholds = c(-Inf, 5.5, 7, Inf),
    levels = rownames(H), category_mode = "noisy", confusion = H), "zero probability")
})

test_that("convergence checks reject deliberately separated and stuck chains", {
  set.seed(3)
  traces <- array(stats::rnorm(200 * 4), c(200L, 4L, 1L),
    dimnames = list(NULL, NULL, "parameter"))
  traces[, 1, 1] <- traces[, 1, 1] + 20
  result <- catmech:::catmech_model_chain_diagnostics(traces)
  expect_false(result$convergence_pass)
  expect_identical(result$status, "failed")
  expect_named(result$summary, c("parameter", "rhat", "ess_bulk", "ess_tail"))
  expect_equal(result$thresholds$min_bulk_ess, 400L)
  expect_equal(result$thresholds$min_tail_ess, 400L)
  traces[] <- 1
  expect_false(catmech:::catmech_model_chain_diagnostics(traces)$convergence_pass)
  # Even otherwise well-behaved short chains cannot satisfy the per-chain ESS floor.
  short <- array(stats::rnorm(20 * 4), c(20L, 4L, 1L),
    dimnames = list(NULL, NULL, "parameter"))
  expect_false(catmech:::catmech_model_chain_diagnostics(short)$convergence_pass)
})
